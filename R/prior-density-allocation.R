# Prior densities of variance-allocation products.
#
# A random-effect SD allocated from a scale prior T is T times a multiplier M
# that is independent of T: the inclusion gates and the Dirichlet allocation
# weights w ~ Dirichlet(a) over all K components of the allocation
# (R/JAGS-deterministic-nodes-random.R). M is a finite mixture of
# * point multipliers: 0 (a gate of the chain off, or no active component of
#   a total), 1 (every component of a total active, or a gate-only chain);
# * mapped Beta shares sqrt(k S), S ~ Beta(a, b): an SD component,
#   sqrt(k w_i) with w_i ~ Beta(a_i, a_- - a_i) (k = K for a mean-variance
#   allocation, 1 otherwise), and a total with the partial active set A,
#   sqrt(sum_{j in A} w_j) ~ sqrt(Beta(a_A, a_-A)).
# The density of T M is the mixture of T's leaves (the components of a
# mixture or spike-and-slab scale prior) times these components: atoms for
# point products (a zero factor or multiplier), the scaled scale prior for a
# point multiplier, and the scale-product quadrature
# f(y) = int Beta(s; a, b) f_T(y / sqrt(k s)) / sqrt(k s) ds for a share
# (R/priors-linear-density-combinations.R, the 'sqrt' multiplier map). A
# variance is the square of the SD, applied to every continuous leaf by the
# 'exp_lin' transformation node and to the atoms directly. The measure is
# recorded as the 'allocation_product' kind of 'adaptive_evaluation', so
# ordinates, region probabilities and plotted densities evaluate one route
# (R/prior-density-route.R).

# A point multiplier at 'location' with probability 'weight'.
.prior_allocation_point <- function(location, weight){

  list(type = "point", location = location, weight = weight)
}

# The mapped share sqrt(scale * S), S ~ Beta(alpha, beta), with probability
# 'weight'.
.prior_allocation_share <- function(alpha, beta, scale, weight){

  list(type = "share", alpha = alpha, beta = beta, scale = scale, weight = weight)
}

.prior_allocation_validate_multipliers <- function(multipliers){

  valid <- is.list(multipliers) && length(multipliers) > 0L &&
    all(vapply(multipliers, function(m){
      is.list(m) && is.numeric(m$weight) && length(m$weight) == 1L &&
        is.finite(m$weight) && m$weight >= 0 &&
        ((identical(m$type, "point") && is.numeric(m$location) &&
            length(m$location) == 1L && is.finite(m$location) && m$location >= 0) ||
           (identical(m$type, "share") &&
              all(vapply(m[c("alpha", "beta", "scale")], function(value){
                is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
              }, logical(1)))))
    }, logical(1)))
  if(!valid){
    stop("Allocation multipliers must be point or Beta-share components with nonnegative weights.",
         call. = FALSE)
  }
  total <- sum(vapply(multipliers, `[[`, numeric(1), "weight"))
  if(abs(total - 1) > sqrt(.Machine$double.eps)){
    stop("Allocation multiplier probabilities must sum to one.", call. = FALSE)
  }
  lapply(multipliers, function(m){
    m$weight <- m$weight / total
    m
  })
}

# The leaves of a scalar scale prior: its positive-probability components
# (recursively for nested mixtures), each a point or a simple continuous
# prior; NULL when a leaf is neither.
.prior_allocation_source_leaves <- function(prior, weight = 1){

  if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior)){
    probabilities <- .prior_density_ordinate_mixture_weights(prior)
    if(is.null(probabilities)){
      return(NULL)
    }
    leaves <- list()
    for(i in which(probabilities > 0)){
      component <- prior[[i]]
      attr(component, "prior_weights") <- NULL
      component_leaves <- .prior_allocation_source_leaves(component, weight * probabilities[[i]])
      if(is.null(component_leaves)){
        return(NULL)
      }
      leaves <- c(leaves, component_leaves)
    }
    return(leaves)
  }
  if(is.prior(prior) && is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    if(!is.numeric(location) || length(location) != 1L || !is.finite(location)){
      return(NULL)
    }
    return(list(list(weight = weight, prior = prior, point = location)))
  }
  if(!.prior_density_simple_continuous(prior)){
    return(NULL)
  }
  list(list(weight = weight, prior = prior, point = NULL))
}

# The (continuous or atomic) route of one leaf T_l * m of the product and its
# source-scale support hull. 'output' is NULL or the 'exp_lin' arguments of
# the output transformation (the square of a variance), composed into a
# point-factor share's own 'exp_lin' map and otherwise applied as a
# transformation node of the continuous leaf.
.prior_allocation_leaf_route <- function(leaf, multiplier, n_grid, output){

  atom <- function(location){
    if(!is.null(output)){
      location <- .density.prior_transformation_x(location, "exp_lin", output)
    }
    .prior_density_route_atom(
      locations   = location,
      probability = 1,
      provenance  = list(kind = "scalar_affine", offset = location, scale = 0)
    )
  }
  transformed <- function(route, hull){
    if(is.null(output)){
      return(route)
    }
    .prior_density_route_transform(route, "exp_lin", output, function() hull)
  }

  if(!is.null(leaf$point)){
    if(identical(multiplier$type, "point") || leaf$point == 0){
      return(atom(leaf$point * if(identical(multiplier$type, "point")) multiplier$location else 0))
    }
    # c sqrt(k S) = exp(log(c^2 k) / 2 + log(S) / 2), composed with the output
    share <- prior("beta", list(alpha = multiplier$alpha, beta = multiplier$beta))
    map <- list(a = log(leaf$point^2 * multiplier$scale) / 2, b = 1 / 2)
    if(!is.null(output)){
      map <- list(a = output$a + output$b * map$a, b = output$b * map$b)
    }
    return(.prior_density_route_transform(
      .prior_density_route_linear(list(source = share), c(source = 1), NULL, n_grid),
      "exp_lin", map, function() c(0, 1)
    ))
  }

  if(identical(multiplier$type, "point")){
    if(multiplier$location == 0){
      return(atom(0))
    }
    weights <- c(source = multiplier$location)
    return(transformed(
      .prior_density_route_linear(list(source = leaf$prior), weights, NULL, n_grid),
      .prior_linear_combination_support_hull(list(source = leaf$prior), weights)
    ))
  }

  spec <- .prior_scale_product_spec(
    offset     = 0,
    scale      = 1,
    factor     = leaf$prior,
    multiplier = prior("beta", list(alpha = multiplier$alpha, beta = multiplier$beta)),
    sources    = list(additive = character(), multiplied = "scale", multiplier = "share"),
    map        = list(type = "sqrt", scale = multiplier$scale)
  )
  transformed(
    list(type = "scale_product", spec = spec, n_grid = n_grid),
    .prior_scale_product_hull(spec)
  )
}

# Route of a recorded allocation product (the mixture over the scale prior's
# leaves and the multiplier components); an 'unknown' route when a leaf of
# the scale prior has no structural density.
.prior_density_route_allocation_product <- function(arguments){

  leaves <- .prior_allocation_source_leaves(arguments$source_prior)
  if(is.null(leaves)){
    return(.prior_density_route_unknown(
      reason     = "The scale prior of the allocation has no structural density route.",
      provenance = list(kind = "unsupported_provenance")
    ))
  }
  n_grid <- if(is.null(arguments$n_grid)) .prior_linear_density_default_grid() else arguments$n_grid
  output <- arguments$output_transformation_arguments
  components <- list()
  weights <- numeric()
  for(multiplier in arguments$multipliers){
    for(leaf in leaves){
      weight <- multiplier$weight * leaf$weight
      if(weight <= 0){
        next
      }
      components[[length(components) + 1L]] <- .prior_allocation_leaf_route(
        leaf, multiplier, n_grid, output
      )
      weights <- c(weights, weight)
    }
  }
  .prior_density_route_mixture(
    components       = components,
    weights          = weights,
    provenance_extra = list(context = "allocation_product")
  )
}

# Closed interval containing the support of T M on the SD scale.
.prior_allocation_product_hull <- function(arguments){

  leaves <- .prior_allocation_source_leaves(arguments$source_prior)
  if(is.null(leaves)){
    return(NULL)
  }
  hulls <- list()
  for(multiplier in arguments$multipliers){
    if(multiplier$weight <= 0){
      next
    }
    range_m <- if(identical(multiplier$type, "point")){
      rep(multiplier$location, 2L)
    }else{
      c(0, sqrt(multiplier$scale))
    }
    for(leaf in leaves){
      range_t <- if(is.null(leaf$point)){
        unlist(leaf$prior$truncation[c("lower", "upper")], use.names = FALSE)
      }else{
        rep(leaf$point, 2L)
      }
      products <- as.vector(outer(range_t, range_m, function(a, b){
        ifelse(a == 0 | b == 0, 0, a * b)
      }))
      hulls[[length(hulls) + 1L]] <- range(products)
    }
  }
  if(length(hulls) == 0L){
    return(NULL)
  }
  range(unlist(hulls))
}

# Numerical range of the continuous part of T M on the SD scale: the leaves'
# tail quantiles times the multipliers' ranges.
.prior_allocation_product_range <- function(arguments, tail_prob){

  leaves <- .prior_allocation_source_leaves(arguments$source_prior)
  lower <- Inf
  upper <- -Inf
  for(multiplier in arguments$multipliers){
    if(multiplier$weight <= 0){
      next
    }
    range_m <- if(identical(multiplier$type, "point")){
      rep(multiplier$location, 2L)
    }else{
      c(0, sqrt(multiplier$scale))
    }
    for(leaf in leaves){
      if(is.null(leaf$point)){
        range_t <- suppressWarnings(as.numeric(quant(leaf$prior, c(tail_prob / 2, 1 - tail_prob / 2))))
      }else if(identical(multiplier$type, "share") && leaf$point > 0){
        range_t <- rep(leaf$point, 2L)
      }else{
        next
      }
      products <- as.vector(outer(range_t, range_m, `*`))
      lower <- min(lower, products)
      upper <- max(upper, products)
    }
  }
  c(lower, upper)
}

# The prior density of T M ('prior_linear_density' with the atoms of the point
# products and the route density of the continuous part on at most 256 values
# of its numerical range; the grid is a display representation, since
# ordinates, heights and probabilities use the route). 'multipliers' are
# point and share components (.prior_allocation_point(),
# .prior_allocation_share()) with probabilities summing to one. NULL when the
# scale prior has a leaf without a structural density.
.prior_allocation_product_density <- function(source_prior, multipliers,
                                              n_grid = .prior_linear_density_default_grid(),
                                              tail_prob = .prior_linear_density_tail_prob(),
                                              square = FALSE){

  multipliers <- .prior_allocation_validate_multipliers(multipliers)
  multipliers <- multipliers[vapply(multipliers, `[[`, numeric(1), "weight") > 0]
  arguments <- list(
    source_prior                    = source_prior,
    multipliers                     = multipliers,
    n_grid                          = n_grid,
    tail_prob                       = tail_prob,
    output_transformation           = if(isTRUE(square)) "exp_lin",
    output_transformation_arguments = if(isTRUE(square)) list(a = 0, b = 2)
  )
  route <- .prior_density_route_allocation_product(arguments)
  if(identical(route$type, "unknown")){
    return(NULL)
  }

  atoms <- Filter(function(i) identical(route$components[[i]]$type, "atom"),
                  seq_along(route$components))
  points <- data.frame(
    x = vapply(route$components[atoms], `[[`, numeric(1), "locations"),
    p = route$weights[atoms]
  )
  continuous_mass <- max(0, 1 - sum(points$p))
  size <- min(max(16L, n_grid), 256L)
  densities <- list()
  dx <- NA_real_
  if(continuous_mass > 0){
    limits <- .prior_allocation_product_range(arguments, tail_prob)
    if(isTRUE(square)){
      limits <- limits^2
    }
    if(!all(is.finite(limits)) || limits[1L] >= limits[2L]){
      stop("The allocation product prior density has no finite numerical range.",
           call. = FALSE)
    }
    x <- seq(limits[1L], limits[2L], length.out = size)
    dx <- x[2L] - x[1L]
    y <- .prior_density_route_density(route, x)
    finite <- is.finite(y)
    densities[[1L]] <- list(x = x[finite], y = y[finite], mass = continuous_mass)
  }
  out <- .prior_linear_density_coalesce(
    densities = densities,
    points    = points,
    dx        = dx,
    n_grid    = size
  )
  attr(out, "adaptive_evaluation") <- list(
    kind      = "allocation_product",
    arguments = arguments
  )
  attr(out, "numerical_diagnostics") <- list(
    n_grid              = n_grid,
    display_grid        = size,
    tail_probability    = tail_prob,
    adaptive_evaluation = TRUE
  )
  out
}
