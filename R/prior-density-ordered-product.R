# The ordered-product leaf of prior-density routes (R/prior-density-route.R).
#
# A level, or a linear combination of levels, of ordered terms is
#
#   X = m + sum_j c_j(B) T_j,
#
# where T_j is the total of a slice of an ordered term (independent draws of the
# total prior per slice and term; one shared draw only where the total is a
# point), m a deterministic offset, and c_j(B) the multilinear function of the
# allocation shares B that the combination leaves in the slice
# (.prior_ordered_linear_slices()). The shares are the independent Beta
# coordinates of the Dirichlet allocations that stay random: an allocation
# whose coefficients take two values in the combination is one cumulative
# share, Beta(sum alpha_A, sum alpha_rest); otherwise the stick-breaking
# coordinates V_j ~ Beta(alpha_j, sum_{k > j} alpha_k) of the Dirichlet over the
# classes of increments that the combination does not distinguish (the
# aggregation property). An allocation record that several slices or terms use
# (an `id`) is one set of coordinates for all of them.
#
# Given B the law of X is closed form: with one random total T it is
# f_T((x - s(B)) / c(B)) / |c(B)|, where s(B) adds the point totals, for every
# simple continuous total with an exact density and distribution function; with
# several normal totals (and normal additive terms) it is the normal
# distribution N(s(B) + sum_j c_j mu_j, sum_j c_j^2 sigma_j^2 + additive
# variance). Ordinates and region probabilities integrate this conditional law
# over B by the tensor Gauss-Jacobi rule of the Beta weights (the weights absorb
# the endpoint singularities of the shares), doubling the rule until two rules
# agree at the relative criterion of the other quadrature leaves; a rule above
# the node cap, a strong singularity (a Beta exponent below 0.1), and a
# disagreement of the last rules reject the result with the integration
# diagnostics of the leaf.
#
# Spike-and-slab and mixture totals are expanded before the leaf is built
# (.prior_density_route_ordered_product()): a shared indicator is one
# component per term, never one per slice.

.prior_ordered_rule_cache <- new.env(parent = emptyenv())

# Largest tensor rule (nodes) the leaf evaluates.
.prior_ordered_product_node_cap <- function(){

  2^21
}

# Gauss-Jacobi rule for E[g(B)], B ~ Beta(a, b), by the Golub-Welsch algorithm:
# nodes u in (0, 1) and weights that sum to 1.
.prior_ordered_beta_rule <- function(n, a, b){

  key <- paste(n, sprintf("%a", a), sprintf("%a", b))
  rule <- .prior_ordered_rule_cache[[key]]
  if(!is.null(rule)){
    return(rule)
  }
  alpha <- b - 1
  beta  <- a - 1
  k <- seq_len(n) - 1
  s <- 2 * k + alpha + beta
  diagonal <- ifelse(s == 0, (beta - alpha) / (alpha + beta + 2),
                     (beta^2 - alpha^2) / (s * (s + 2)))
  jacobi <- diag(diagonal, n)
  if(n > 1L){
    j <- seq_len(n - 1L)
    s2 <- 2 * j + alpha + beta
    off <- 2 / s2 * sqrt(j * (j + alpha) * (j + beta) * (j + alpha + beta) / ((s2 + 1) * (s2 - 1)))
    # the first recurrence coefficient in the form without the 0/0 at alpha + beta = -1
    off[1L] <- 2 * sqrt((1 + alpha) * (1 + beta) / ((alpha + beta + 2)^2 * (alpha + beta + 3)))
    jacobi[cbind(j, j + 1L)] <- off
    jacobi[cbind(j + 1L, j)] <- off
  }
  decomposition <- eigen(jacobi, symmetric = TRUE)
  nodes <- sort((decomposition$values + 1) / 2, index.return = TRUE)
  weights <- decomposition$vectors[1L, nodes$ix]^2
  rule <- list(u = nodes$x, w = weights / sum(weights))
  assign(key, rule, envir = .prior_ordered_rule_cache)
  rule
}

# ---- structure: classes of increments and independent Beta coordinates -------

# The matrix of a tensor unfolded along the axis of the record 'key' (a row per
# increment).
.prior_ordered_product_unfold <- function(tensor, key){

  dimensions <- vapply(tensor$records, `[[`, integer(1), "dim")
  axis <- match(key, names(tensor$records))
  coefficients <- array(tensor$coefficients, dim = dimensions)
  matrix(aperm(coefficients, c(axis, setdiff(seq_along(dimensions), axis))), nrow = dimensions[[axis]])
}

.prior_ordered_product_row_labels <- function(unfolded){

  apply(unfolded, 1L, function(row) paste(sprintf("%a", row), collapse = "|"))
}

# The independent Beta coordinates of the random records of 'terms' (each with
# a reduced 'tensor'), and for every record the classes of its increments that
# no term distinguishes. A record with two classes is one cumulative share; a
# record with k > 2 classes has the k - 1 stick-breaking coordinates of the
# Dirichlet aggregated over its classes.
.prior_ordered_product_structure <- function(terms){

  keys <- unique(unlist(lapply(terms, function(term) names(term$tensor$records)), use.names = FALSE))
  records <- list()
  coordinates <- list()
  for(key in keys){
    labels <- NULL
    record <- NULL
    for(term in terms){
      if(!key %in% names(term$tensor$records)){
        next
      }
      record <- term$tensor$records[[key]]
      term_labels <- .prior_ordered_product_row_labels(.prior_ordered_product_unfold(term$tensor, key))
      labels <- if(is.null(labels)) term_labels else paste(labels, term_labels, sep = "/")
    }
    classes <- match(labels, unique(labels))
    n_classes <- max(classes)
    class_alpha <- vapply(seq_len(n_classes), function(k) sum(record$spec$alpha[classes == k]), numeric(1))
    first <- length(coordinates) + 1L
    if(n_classes == 2L){
      coordinates[[first]] <- list(record = key, a = class_alpha[[1L]], b = class_alpha[[2L]])
    }else{
      for(j in seq_len(n_classes - 1L)){
        coordinates[[length(coordinates) + 1L]] <- list(
          record = key, a = class_alpha[[j]], b = sum(class_alpha[(j + 1L):n_classes])
        )
      }
    }
    records[[key]] <- list(
      key            = key,
      n_classes      = n_classes,
      class_alpha    = class_alpha,
      representative = match(seq_len(n_classes), classes),
      coordinates    = seq(first, length(coordinates))
    )
  }
  list(records = records, coordinates = coordinates)
}

# The class tensor of a term: its coefficients on one increment per class of
# each of its axes (the coefficients are constant within a class).
.prior_ordered_product_class_tensor <- function(term, structure){

  keys <- names(term$tensor$records)
  if(length(keys) == 0L){
    return(list(axes = character(), coefficients = array(term$tensor$coefficients[[1L]], dim = 1L)))
  }
  dimensions <- vapply(term$tensor$records, `[[`, integer(1), "dim")
  coefficients <- array(term$tensor$coefficients, dim = dimensions)
  representatives <- lapply(keys, function(key) structure$records[[key]]$representative)
  list(axes = keys, coefficients = do.call(`[`, c(list(coefficients), representatives, list(drop = FALSE))))
}

# Allocation shares of every record at the nodes 'U' (a row per node, a column
# per coordinate): a column per class.
.prior_ordered_product_shares <- function(structure, U){

  lapply(structure$records, function(record){
    V <- U[, record$coordinates, drop = FALSE]
    if(record$n_classes == 2L){
      return(cbind(V[, 1L], 1 - V[, 1L]))
    }
    shares <- matrix(0, nrow(U), record$n_classes)
    remaining <- rep(1, nrow(U))
    for(j in seq_len(record$n_classes - 1L)){
      shares[, j] <- remaining * V[, j]
      remaining <- remaining * (1 - V[, j])
    }
    shares[, record$n_classes] <- remaining
    shares
  })
}

# The coefficient c_j(B) of every term at the nodes: the multilinear sum of the
# class coefficients times the shares of the classes.
.prior_ordered_product_coefficients <- function(terms, shares, n){

  lapply(terms, function(term){
    tensor <- term$class_tensor
    out <- numeric(n)
    for(index in which(tensor$coefficients != 0)){
      value <- rep(tensor$coefficients[[index]], n)
      if(length(tensor$axes) > 0L){
        position <- arrayInd(index, dim(tensor$coefficients))
        for(axis in seq_along(tensor$axes)){
          value <- value * shares[[tensor$axes[[axis]]]][, position[[axis]]]
        }
      }
      out <- out + value
    }
    out
  })
}

# Shift, scale and normal moments of the conditional law at the nodes 'U'.
.prior_ordered_product_state <- function(spec, U){

  n <- nrow(U)
  shares <- .prior_ordered_product_shares(spec$structure, U)
  coefficients <- .prior_ordered_product_coefficients(spec$terms, shares, n)
  shift <- rep(spec$offset, n)
  for(j in spec$point_terms){
    shift <- shift + spec$terms[[j]]$location * coefficients[[j]]
  }
  state <- list(shift = shift)
  if(identical(spec$mode, "scale")){
    state$scale <- coefficients[[spec$random_terms[[1L]]]]
  }else if(identical(spec$mode, "normal")){
    mean <- shift + spec$additive_mean
    variance <- rep(spec$additive_sd^2, n)
    for(j in spec$random_terms){
      mean <- mean + coefficients[[j]] * spec$terms[[j]]$mean
      variance <- variance + coefficients[[j]]^2 * spec$terms[[j]]$sd^2
    }
    state$mean <- mean
    state$sd <- sqrt(variance)
  }
  state
}

# ---- conditional law ---------------------------------------------------------

# Conditional density of X at the values 'x' given B at the nodes of 'state'
# (a row per node, a column per value).
.prior_ordered_product_conditional_density <- function(spec, state, x){

  n <- length(state$shift)
  values <- matrix(x, nrow = n, ncol = length(x), byrow = TRUE)
  if(identical(spec$mode, "normal")){
    out <- stats::dnorm((values - state$mean) / state$sd) / state$sd
    # every coefficient zero at a node (a point, a null set of the shares)
    out[!is.finite(out)] <- 0
    return(out)
  }
  term <- spec$terms[[spec$random_terms[[1L]]]]
  total <- term$total
  argument <- (values - state$shift) / state$scale
  if(!is.null(term$mean)){
    # a normal total: the density of the argument directly (no coefficient
    # c(B) = 0, where the conditional law is a point, has a density)
    out <- stats::dnorm(argument, term$mean, term$sd) / abs(state$scale)
    out[!is.finite(out)] <- 0
    return(out)
  }
  usable <- which(state$scale != 0 & is.finite(argument))
  out <- matrix(0, n, length(x))
  if(length(usable) > 0L){
    rows <- ((usable - 1L) %% n) + 1L
    log_density <- .prior_quadrature_log_density(total, argument[usable], spec$log_density) -
      log(abs(state$scale[rows]))
    out[usable] <- exp(log_density)
  }
  out
}

# Conditional probability of the continuous region 'intervals' given B at the
# nodes of 'state' (a vector over the nodes).
.prior_ordered_product_conditional_probability <- function(spec, state, intervals){

  n <- length(state$shift)
  probability <- numeric(n)
  for(i in seq_len(nrow(intervals))){
    lower <- intervals[i, 1L]
    upper <- intervals[i, 2L]
    if(identical(spec$mode, "normal")){
      probability <- probability +
        .prior_normal_interval_probability(lower, upper, state$mean, state$sd)
      next
    }
    total <- spec$terms[[spec$random_terms[[1L]]]]$total
    zero <- state$scale == 0
    scale <- ifelse(zero, 1, state$scale)
    a <- (lower - state$shift) / scale
    b <- (upper - state$shift) / scale
    flip <- scale < 0
    mass <- .prior_scalar_interval_probability(total, ifelse(flip, b, a), ifelse(flip, a, b))
    mass[zero] <- as.numeric(state$shift[zero] > lower & state$shift[zero] < upper)
    probability <- probability + mass
  }
  probability
}

# ---- tensor quadrature ---------------------------------------------------------

# The sequence of rule sizes (nodes per coordinate) of a leaf with 'p'
# coordinates: doubling from the first size while the tensor rule stays within
# the node cap.
.prior_ordered_product_sizes <- function(p){

  first <- if(p <= 2L) 16L else if(p == 3L) 12L else if(p == 4L) 8L else 6L
  sizes <- first * 2L^(0:6)
  sizes <- sizes[sizes <= 512L & as.numeric(sizes)^p <= .prior_ordered_product_node_cap()]
  as.integer(sizes)
}

# Nodes (a row per node, a column per coordinate) and weights of the tensor
# Gauss-Jacobi rule with 'n' nodes per coordinate.
.prior_ordered_product_nodes <- function(spec, n){

  coordinates <- spec$structure$coordinates
  p <- length(coordinates)
  total <- as.numeric(n)^p
  nodes <- matrix(0, total, p)
  weights <- rep(1, total)
  each <- 1
  for(k in seq_len(p)){
    rule <- .prior_ordered_beta_rule(n, coordinates[[k]]$a, coordinates[[k]]$b)
    index <- rep(rep(seq_len(n), each = each), length.out = total)
    nodes[, k] <- rule$u[index]
    weights <- weights * rule$w[index]
    each <- each * n
  }
  list(U = nodes, w = weights)
}

# The expectation over B of 'evaluate(state, columns)' (a matrix with a row per
# node and a column per selected target value) for 'm' target values, by the
# doubling tensor rule: a value is accepted when the two last rules agree at
# the relative criterion of the quadrature leaves (or 'tolerance'; a 'budget'
# of evaluations per rule ends the doubling of plotted values early). Returns, per value, the
# integral (NA unless accepted), the last estimate, the absolute change of the
# last two rules, and the total evaluations.
.prior_ordered_product_integrate <- function(spec, evaluate, m = 1L,
                                              tolerance = .prior_linear_density_refinement_tolerance()$relative,
                                              budget = Inf){

  p <- length(spec$structure$coordinates)
  sizes <- .prior_ordered_product_sizes(p)
  value <- rep(NA_real_, m)
  estimate <- rep(NA_real_, m)
  change <- rep(NA_real_, m)
  accepted <- rep(FALSE, m)
  previous <- rep(NA_real_, m)
  evaluations <- 0
  rules <- 0L
  message <- "OK"
  if(length(sizes) < 2L){
    message <- paste0("the node cap admits one rule only for ", p, " allocation coordinates")
  }else{
    for(size in sizes){
      active <- which(!accepted)
      if(length(active) == 0L){
        break
      }
      # a plotted density stops doubling when a rule would exceed its budget of
      # evaluations, once two rules can be compared
      if(as.numeric(size)^p * length(active) > budget && rules >= 2L){
        break
      }
      rules <- rules + 1L
      nodes <- .prior_ordered_product_nodes(spec, size)
      state <- .prior_ordered_product_state(spec, nodes$U)
      current <- rep(NA_real_, length(active))
      chunk <- max(1L, floor(4e6 / nrow(nodes$U)))
      for(start in seq.int(1L, length(active), by = chunk)){
        columns <- active[start:min(length(active), start + chunk - 1L)]
        values <- evaluate(state, columns)
        current[match(columns, active)] <- colSums(nodes$w * values)
      }
      evaluations <- evaluations + nrow(nodes$U) * length(active)
      difference <- abs(current - previous[active])
      estimate[active] <- current
      change[active] <- difference
      good <- is.finite(current) & current > 0 & is.finite(difference) &
        difference <= tolerance * current
      accepted[active[good]] <- TRUE
      value[active[good]] <- current[good]
      previous[active] <- current
    }
  }
  list(value = value, estimate = estimate, change = change, accepted = accepted,
       evaluations = evaluations, message = message, sizes = sizes)
}

# The integral record of one target value in the shape of the other quadrature
# leaves ('value' and the 'integration' diagnostics).
.prior_ordered_product_integral <- function(result, i, zero_message){

  accepted <- isTRUE(result$accepted[[i]])
  message <- result$message
  if(identical(message, "OK") && !accepted){
    message <- if(isTRUE(result$estimate[[i]] == 0)){
      zero_message
    }else if(!is.finite(result$estimate[[i]]) || !is.finite(result$change[[i]])){
      "the integral over the allocation shares is not finite"
    }else{
      paste0("relative error estimate ",
             format(signif(result$change[[i]] / abs(result$estimate[[i]]), 3)))
    }
  }
  integration <- list(
    kind = "scale_mixture", exact = FALSE,
    absolute_error = result$change[[i]],
    error_bound = .prior_linear_density_refinement_tolerance()$relative * abs(result$estimate[[i]]),
    evaluations = result$evaluations,
    budget = .prior_ordered_product_node_cap(),
    converged = accepted, message = message,
    refined = FALSE, rule_sizes = result$sizes
  )
  if(!accepted && is.finite(result$estimate[[i]]) && result$estimate[[i]] > 0){
    integration$estimate <- result$estimate[[i]]
  }
  list(value = if(accepted) result$value[[i]] else NA_real_, integration = integration)
}

# ---- support, offset and classification ----------------------------------------

# The corner nodes of the unit cube of the coordinates, where the multilinear
# coefficients and shifts attain their extremes.
.prior_ordered_product_corners <- function(structure){

  p <- length(structure$coordinates)
  as.matrix(expand.grid(rep(list(c(0, 1)), p), KEEP.OUT.ATTRS = FALSE))
}

# Range of the product of a coefficient vector with the bounds of the total.
.prior_ordered_product_bound_products <- function(coefficient, bounds){

  products <- as.vector(outer(coefficient, bounds, function(a, b) ifelse(a == 0, 0, a * b)))
  products[!is.nan(products)]
}

# The closed interval that contains the support of X.
.prior_ordered_product_hull <- function(spec){

  if(identical(spec$mode, "normal")){
    return(c(-Inf, Inf))
  }
  state <- .prior_ordered_product_state(spec, .prior_ordered_product_corners(spec$structure))
  total <- spec$terms[[spec$random_terms[[1L]]]]$total
  bounds <- unlist(total$truncation[c("lower", "upper")], use.names = FALSE)
  # X is multilinear in the shares for a fixed total: its extremes are at the
  # corners of the shares and at the bounds of the total
  range(unlist(lapply(seq_along(state$shift), function(i){
    state$shift[[i]] + .prior_ordered_product_bound_products(state$scale[[i]], bounds)
  })))
}

# The offset of the conditional law when it does not depend on B, else NULL.
.prior_ordered_product_constant_shift <- function(spec){

  if(!identical(spec$mode, "scale")){
    return(NULL)
  }
  shift <- .prior_ordered_product_state(spec, .prior_ordered_product_corners(spec$structure))$shift
  if(length(unique(shift)) == 1L) shift[[1L]] else NULL
}

# The factors Beta(a, b) whose product is |c(B)| / |kappa| when the single
# random term has one nonzero class coefficient (a monomial); NULL otherwise.
.prior_ordered_product_monomial <- function(spec){

  if(!identical(spec$mode, "scale")){
    return(NULL)
  }
  term <- spec$terms[[spec$random_terms[[1L]]]]
  tensor <- term$class_tensor
  nonzero <- which(tensor$coefficients != 0)
  if(length(nonzero) != 1L){
    return(NULL)
  }
  position <- if(length(tensor$axes) > 0L) arrayInd(nonzero, dim(tensor$coefficients)) else integer()
  factors <- list()
  for(axis in seq_along(tensor$axes)){
    record <- spec$structure$records[[tensor$axes[[axis]]]]
    coordinates <- spec$structure$coordinates[record$coordinates]
    class <- position[[axis]]
    if(record$n_classes == 2L){
      parameters <- c(coordinates[[1L]]$a, coordinates[[1L]]$b)
      factors[[length(factors) + 1L]] <- if(class == 1L) parameters else rev(parameters)
      next
    }
    # G_class = V_class * prod_{k < class} (1 - V_k)
    if(class < record$n_classes){
      factors[[length(factors) + 1L]] <- c(coordinates[[class]]$a, coordinates[[class]]$b)
    }
    for(k in seq_len(class - 1L)){
      factors[[length(factors) + 1L]] <- c(coordinates[[k]]$b, coordinates[[k]]$a)
    }
  }
  list(coefficient = tensor$coefficients[[nonzero]], factors = factors)
}

# ---- provenance ------------------------------------------------------------------

.prior_ordered_product_provenance <- function(spec){

  list(
    kind                = "scale_mixture",
    ordered_product     = TRUE,
    mode                = spec$mode,
    offset              = spec$offset,
    coordinates         = vapply(spec$structure$coordinates, function(coordinate){
      paste0(coordinate$record, ": Beta(", format(coordinate$a), ", ", format(coordinate$b), ")")
    }, character(1)),
    terms               = vapply(spec$terms, function(term){
      if(identical(term$kind, "point")) paste0("point ", format(term$location)) else
        .prior_density_ordinate_prior_provenance(term$total)$family
    }, character(1)),
    independent_sources = vapply(spec$terms, `[[`, character(1), "label"),
    support             = spec$hull
  )
}

.prior_ordered_product_route_provenance <- function(spec){

  provenance <- .prior_ordered_product_provenance(spec)
  provenance$offset_behavior <- "unknown"
  provenance
}

# ---- ordinate --------------------------------------------------------------------

.prior_ordered_product_ordinate <- function(spec, value, n_grid){

  provenance <- .prior_ordered_product_provenance(spec)
  hull <- spec$hull
  if(value < hull[1L] || value > hull[2L]){
    return(.prior_density_ordinate_result(
      value = value, behavior = "zero", log_density = -Inf, exact = TRUE,
      method = "scale_mixture", reason = "The requested value is outside the prior support.",
      provenance = provenance
    ))
  }
  shift <- .prior_ordered_product_constant_shift(spec)
  if(!is.null(shift) && value == shift){
    offset <- .prior_ordered_product_offset_ordinate(spec, value, provenance)
    if(!is.null(offset)){
      return(offset)
    }
  }
  if(value == hull[1L] || value == hull[2L]){
    return(.prior_density_ordinate_result(
      value = value, behavior = "unknown", log_density = NA_real_, exact = FALSE,
      method = "unsupported_provenance",
      reason = "The limit at a finite bound of a product support is not structurally classified.",
      provenance = provenance
    ))
  }
  provenance$structural_regularity <- "scale_mixture_inside_support"
  if(!is.null(shift) && !.prior_density_affine_full_precision(value, shift, 1)){
    return(.prior_density_ordinate_imprecise(
      value, "The distance of the requested value from the product's offset",
      "scale_mixture", provenance
    ))
  }
  strong <- .prior_ordered_product_strong_singularity(spec)
  if(!is.null(strong)){
    integral <- .prior_ordered_product_not_attempted(strong)
  }else{
    result <- .prior_ordered_product_integrate(spec, function(state, columns){
      .prior_ordered_product_conditional_density(spec, state, value)
    }, 1L)
    integral <- .prior_ordered_product_integral(result, 1L,
      "zero ordinate for a structurally positive density")
  }
  .prior_density_quadrature_ordinate(value, integral, "scale_mixture", provenance)
}

# The rejected integral of a quadrature that is not attempted ('message' is the
# reason): no value and no error bound.
.prior_ordered_product_not_attempted <- function(message){

  list(value = NA_real_, integration = list(
    kind = "scale_mixture", exact = FALSE, absolute_error = Inf, error_bound = NA_real_,
    evaluations = 0, budget = .prior_ordered_product_node_cap(), converged = FALSE,
    message = message, refined = FALSE
  ))
}

# A Beta coordinate with an exponent below the strong-singularity limit of the
# quadrature leaves: the tensor rule's error estimate is not reliable there.
.prior_ordered_product_strong_singularity <- function(spec){

  limit <- .prior_density_quadrature_strong_exponent()
  for(coordinate in spec$structure$coordinates){
    if(coordinate$a < limit || coordinate$b < limit){
      return(paste0("an allocation share has a Beta(", format(coordinate$a), ", ",
                    format(coordinate$b), ") coordinate with an exponent below ", format(limit),
                    ", a strong singularity of the density at its bound"))
    }
  }
  if(length(spec$structure$coordinates) > 0L && length(.prior_ordered_product_sizes(
    length(spec$structure$coordinates))) < 2L){
    return(paste0("the product involves ", length(spec$structure$coordinates),
                  " allocation coordinates, more than the node cap allows"))
  }
  NULL
}

# The ordinate at the offset of a single random total times a monomial of
# independent Beta variables: f_T(0) E[1 / |c|] with the closed form moments of
# the Beta variables. NULL where it is not classified (a general polynomial,
# several totals, a total whose support does not contain 0 as an interior
# point or whose density at 0 is not positive).
.prior_ordered_product_offset_ordinate <- function(spec, value, provenance){

  monomial <- .prior_ordered_product_monomial(spec)
  if(is.null(monomial)){
    return(NULL)
  }
  total <- spec$terms[[spec$random_terms[[1L]]]]$total
  bounds <- unlist(total$truncation[c("lower", "upper")], use.names = FALSE)
  if(!(bounds[[1L]] < 0 && bounds[[2L]] > 0)){
    return(NULL)
  }
  at_zero <- .prior_density_ordinate_primitive(total, 0)
  behavior <- .prior_density_ordinate_continuous_behavior(at_zero)
  provenance$structural_regularity <- "scale_mixture_offset"
  infinite_moment <- any(vapply(monomial$factors, function(f) f[[1L]] <= 1, logical(1)))
  if(identical(behavior, "infinite") || (identical(behavior, "regular") && infinite_moment)){
    return(.prior_density_ordinate_result(
      value = value, behavior = "infinite", log_density = Inf, exact = TRUE,
      method = "scale_mixture",
      reason = paste0(
        "A total multiplied by allocation shares whose density at zero is positive or ",
        "infinite has an infinite density at the requested value."
      ),
      provenance = provenance
    ))
  }
  if(!identical(behavior, "regular")){
    return(NULL)
  }
  log_moment <- sum(vapply(monomial$factors, function(f){
    lbeta(f[[1L]] - 1, f[[2L]]) - lbeta(f[[1L]], f[[2L]])
  }, numeric(1)))
  log_density <- at_zero$log_density - log(abs(monomial$coefficient)) + log_moment
  .prior_density_ordinate_result(
    value = value, behavior = "regular", log_density = log_density, exact = is.finite(log_density),
    method = "scale_mixture", provenance = provenance
  )
}

# ---- region probability -----------------------------------------------------------

.prior_ordered_product_region <- function(spec, region, n_grid){

  intervals <- region$intervals
  strong <- .prior_ordered_product_strong_singularity(spec)
  if(!is.null(strong)){
    return(.prior_ordered_product_not_attempted(strong))
  }
  result <- .prior_ordered_product_integrate(spec, function(state, columns){
    matrix(.prior_ordered_product_conditional_probability(spec, state, intervals),
           ncol = length(columns), nrow = length(state$shift))
  }, 1L)
  .prior_ordered_product_integral(result, 1L, "zero probability for a structurally positive region")
}

# ---- plotted density -----------------------------------------------------------------

# Continuous density at the values 'x' (display): the offset and the support
# bounds take their ordinates, other values the doubling rule accepted at a
# display criterion (a relative change of 1e-3 between the last two rules, which
# stop at a budget of 5e7 evaluations), else NA.
.prior_ordered_product_density <- function(spec, x){

  out <- rep(NA_real_, length(x))
  hull <- spec$hull
  finite <- is.finite(x)
  out[finite & (x < hull[1L] | x > hull[2L])] <- 0
  shift <- .prior_ordered_product_constant_shift(spec)
  at_shift <- if(is.null(shift)) rep(FALSE, length(x)) else x == shift
  special <- finite & (at_shift | x == hull[1L] | x == hull[2L])
  open <- which(finite & is.na(out) & !special)
  if(length(open) > 0L && is.null(.prior_ordered_product_strong_singularity(spec))){
    result <- .prior_ordered_product_integrate(spec, function(state, columns){
      .prior_ordered_product_conditional_density(spec, state, x[open][columns])
    }, length(open), tolerance = 1e-3, budget = 5e7)
    out[open[result$accepted]] <- result$estimate[result$accepted]
  }
  for(i in which(special)){
    out[[i]] <- .prior_density_ordinate_height_value(
      .prior_density_route_ordinate(list(type = "ordered_product", spec = spec, n_grid = NULL), x[[i]])
    )
  }
  out
}
