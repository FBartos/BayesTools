# Batched adaptive quadrature of prior densities at many values.
#
# A plotted density whose route has quadrature leaves (conditional-normal
# mixtures, scale mixtures and two-term convolutions) needs the ordinate's 1-D
# integral at every displayed value. Instead of one adaptive QUADPACK call per
# value, the integrals of all values are refined together: each value keeps the
# breakpoints of its ordinate (the same pieces), every round evaluates all
# pending intervals of all values in one call of an integrand vectorized over
# (node, value) pairs, and the intervals of values whose summed error estimate
# exceeds their tolerance are bisected. The rules and error estimates are
# QUADPACK's: qk21 on finite intervals, and qk15i on an interval with an
# infinite end, mapped to t in (0, 1] by x = a + d (1 - t) / t (d = +-1). The
# acceptance criterion, summed error <= relative * value, is applied per value
# with a positive finite value (a purely relative criterion, so a saturated
# first estimate of a far-tail value is refined as well); a value that does not
# meet it within the round cap is NA, and callers then evaluate its ordinate.
# Bisection without extrapolation converges too slowly next to an integrable
# singularity of the integrand (a term with an infinite density at a finite
# support bound, .prior_density_singular_bounds()) to reach 1e-8, so such
# leaves are batched only for display grids, with the ordinates' acceptance
# criterion of 1e-4 (.prior_density_route_quadrature_density()). Next to a
# strong singularity (.prior_density_strong_singularity()) the error estimate
# of the bisection is not reliable at that criterion either, so those leaves
# are not batched at all.

.prior_density_quadrature_tolerance <- function(){

  list(relative = 1e-8, max_rounds = 60L)
}

# Whether the declared density of a simple continuous prior is infinite at each
# of the finite values 'bounds' (its support bounds by default).
.prior_density_singular_bounds <- function(prior,
                                           bounds = unlist(prior$truncation[c("lower", "upper")],
                                                           use.names = FALSE)){

  vapply(bounds, function(bound){
    is.finite(bound) && isTRUE(is.infinite(suppressWarnings(exp(lpdf(prior, bound)))))
  }, logical(1))
}

# Smallest exponent p of an infinite density at a finite bound (a density
# ~ distance^(p - 1) there: the shape of a gamma density at 0, of a beta
# density at either bound, .prior_density_bound_exponent()) that is not a
# strong singularity.
.prior_density_quadrature_strong_exponent <- function(){

  0.1
}

# Whether a simple continuous prior has a strong singularity at one of the
# finite values 'bounds' (its support bounds by default): an infinite density
# whose exponent is below .prior_density_quadrature_strong_exponent() or not
# known (e.g. a Beta share of a Dirichlet(0.05) allocation). The batched
# bisection's error estimate is not reliable next to such a singularity, also
# where the other term's density cancels it at the bound (a zero bound of the
# share of a product): accepted at an estimate of 1e-4, grid values of a
# Beta(0.05, 0.15) share were 4.3e-4 off, and of a Beta(0.3, 0.05) multiplier
# 2.8e-4. Leaves with such a term therefore take per-value ordinates even in
# display grids (.prior_density_route_quadrature_density()).
.prior_density_strong_singularity <- function(prior,
                                              bounds = unlist(prior$truncation[c("lower", "upper")],
                                                              use.names = FALSE)){

  singular <- .prior_density_singular_bounds(prior, bounds)
  if(!any(singular)){
    return(FALSE)
  }
  exponents <- vapply(bounds[singular], function(bound){
    .prior_density_bound_exponent(prior, bound)
  }, numeric(1))
  any(is.na(exponents) | exponents < .prior_density_quadrature_strong_exponent())
}

# Gauss-Kronrod rules on [-1, 1] (QUADPACK qk21 and qk15i): nodes, Kronrod
# weights and the embedded Gauss weights (zero at Kronrod-only nodes).
.prior_density_quadrature_rules <- function(){

  mirror <- function(half){
    c(half[-length(half)], half[length(half)], rev(half[-length(half)]))
  }
  k21_nodes <- c(
    0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
    0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
    0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
    0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
    0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
    0
  )
  k21_kronrod <- c(
    0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077600525472847, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821
  )
  k21_gauss <- c(
    0, 0.066671344308688137593568809893332,
    0, 0.149451349150580593145776339657697,
    0, 0.219086362515982043995534934228163,
    0, 0.269266719309996355091226921569469,
    0, 0.295524224714752870173892994651338,
    0
  )
  k15_nodes <- c(
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0
  )
  k15_kronrod <- c(
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714
  )
  k15_gauss <- c(
    0, 0.129484966168869693270611432679082,
    0, 0.279705391489276667901467771423780,
    0, 0.381830050505118944950369775488975,
    0, 0.417959183673469387755102040816327
  )
  list(
    finite   = list(nodes   = c(-k21_nodes[-11L], k21_nodes[11L], rev(k21_nodes[-11L])),
                    kronrod = mirror(k21_kronrod),
                    gauss   = mirror(k21_gauss)),
    infinite = list(nodes   = c(-k15_nodes[-8L], k15_nodes[8L], rev(k15_nodes[-8L])),
                    kronrod = mirror(k15_kronrod),
                    gauss   = mirror(k15_gauss))
  )
}

# Kronrod estimates and QUADPACK error estimates of the intervals
# [lower, upper] (in t for an infinite end) of 'kind' 0 (finite), 1 (x = anchor
# + (1 - t) / t) or -1 (x = anchor - (1 - t) / t), with 'index' the value of
# each interval. 'integrand(nodes, index)' is vectorized over pairs; an
# integrand given as list(shared, value) is value(shared(nodes), index), where
# 'shared' returns a list of the terms that do not depend on the value. The
# values of a leaf share most of their intervals, so the nodes and the shared
# terms are computed once per distinct interval.
.prior_density_quadrature_evaluate <- function(integrand, lower, upper, kind,
                                               anchor, index, rules){

  value <- error <- numeric(length(lower))
  for(infinite in c(FALSE, TRUE)){
    selected <- which((kind != 0L) == infinite)
    if(length(selected) == 0L){
      next
    }
    rule <- if(infinite) rules$infinite else rules$finite
    k <- length(rule$nodes)
    interval <- complex(real = lower[selected], imaginary = upper[selected])
    end <- complex(real = kind[selected], imaginary = anchor[selected])
    key <- match(interval, interval) + (match(end, end) - 1) * length(selected)
    distinct <- which(!duplicated(key))
    map <- match(key, key[distinct])
    first <- selected[distinct]
    centre <- (lower[first] + upper[first]) / 2
    half <- (upper[first] - lower[first]) / 2
    nodes <- rep(centre, each = k) + rep(half, each = k) * rule$nodes
    source <- if(infinite){
      direction <- rep(kind[first], each = k)
      rep(anchor[first], each = k) + direction * (1 - nodes) / nodes
    }else{
      nodes
    }
    columns <- rep((map - 1L) * k, each = k) + seq_len(k)
    pair_index <- rep(index[selected], each = k)
    f <- if(is.function(integrand)){
      integrand(source[columns], pair_index)
    }else{
      integrand$value(lapply(integrand$shared(source), `[`, columns), pair_index)
    }
    if(infinite){
      f <- f / (nodes^2)[columns]
    }
    half <- half[map]
    f <- matrix(f, nrow = k)
    # an interval with a non-finite integrand value gets an infinite error
    finite <- is.finite(f)
    invalid <- colSums(!finite) > 0
    f[!finite] <- 0
    kronrod <- colSums(rule$kronrod * f)
    gauss <- colSums(rule$gauss * f)
    absolute <- colSums(rule$kronrod * abs(f))
    mean_f <- kronrod / 2
    spread <- colSums(rule$kronrod * abs(f - rep(mean_f, each = k)))
    estimate_error <- abs((kronrod - gauss) * half)
    absolute <- absolute * abs(half)
    spread <- spread * abs(half)
    scaled <- spread != 0 & estimate_error != 0
    estimate_error[scaled] <- spread[scaled] *
      pmin(1, (200 * estimate_error[scaled] / spread[scaled])^1.5)
    roundoff <- absolute > .Machine$double.xmin / (50 * .Machine$double.eps)
    estimate_error[roundoff] <- pmax(50 * .Machine$double.eps * absolute[roundoff],
                                     estimate_error[roundoff])
    estimate <- kronrod * half
    invalid <- invalid | !is.finite(estimate) | !is.finite(estimate_error)
    estimate_error[invalid] <- Inf
    value[selected] <- estimate
    error[selected] <- estimate_error
  }
  list(value = value, error = error)
}

# Integrals over the pieces between consecutive 'breakpoints[[j]]' (sorted;
# infinite ends allowed) of integrand(nodes, index) (or its split form, see
# .prior_density_quadrature_evaluate()) for each value j, refined
# in one batched loop. Returns the integrals, NA where the acceptance
# criterion is not met (including zero and non-finite integrals).
.prior_density_quadrature_batch <- function(integrand, breakpoints,
                                            tolerance = .prior_density_quadrature_tolerance()){

  n_values <- length(breakpoints)
  out <- rep(NA_real_, n_values)
  if(n_values == 0L){
    return(out)
  }
  rules <- .prior_density_quadrature_rules()

  pieces <- lengths(breakpoints) - 1L
  index <- rep(seq_len(n_values), pmax(pieces, 0L))
  a <- unlist(lapply(breakpoints, function(points) points[-length(points)]), use.names = FALSE)
  b <- unlist(lapply(breakpoints, function(points) points[-1L]), use.names = FALSE)
  if(length(a) == 0L){
    return(out)
  }
  # a doubly infinite piece is split at zero
  both <- is.infinite(a) & is.infinite(b)
  if(any(both)){
    index <- c(index[!both], index[both], index[both])
    a_new <- c(a[!both], a[both], rep(0, sum(both)))
    b <- c(b[!both], rep(0, sum(both)), b[both])
    a <- a_new
  }
  kind <- ifelse(is.finite(a) & is.finite(b), 0L, ifelse(is.finite(a), 1L, -1L))
  anchor <- ifelse(kind == 1L, a, ifelse(kind == -1L, b, 0))
  lower <- ifelse(kind == 0L, a, 0)
  upper <- ifelse(kind == 0L, b, 1)
  keep <- lower < upper
  lower <- lower[keep]; upper <- upper[keep]; kind <- kind[keep]
  anchor <- anchor[keep]; index <- index[keep]

  estimate <- .prior_density_quadrature_evaluate(
    integrand, lower, upper, kind, anchor, index, rules
  )
  value <- estimate$value
  error <- estimate$error
  active <- rep(TRUE, n_values)
  # per-value summaries of the intervals' 'x' (NA for a value without
  # intervals), each over the value's intervals in their order
  per_value <- function(x, summary){
    present <- which(tabulate(index, nbins = n_values) > 0L)
    codes <- integer(n_values)
    codes[present] <- seq_along(present)
    groups <- structure(codes[index], levels = as.character(seq_along(present)),
                        class = "factor")
    out <- rep(NA_real_, n_values)
    out[present] <- vapply(base::split(x, groups), summary, numeric(1))
    out
  }
  for(round in seq_len(tolerance$max_rounds + 1L)){
    total <- per_value(value, sum)
    total_error <- per_value(error, sum)
    bound <- tolerance$relative * abs(total)
    accepted <- active & is.finite(total) & total > 0 & is.finite(total_error) &
      total_error <= bound
    out[accepted] <- total[accepted]
    active <- active & !accepted
    if(!any(active) || round > tolerance$max_rounds){
      break
    }
    # only the intervals of pending values are refined further
    pending <- active[index]
    lower <- lower[pending]; upper <- upper[pending]; kind <- kind[pending]
    anchor <- anchor[pending]; index <- index[pending]
    value <- value[pending]; error <- error[pending]
    # bisect the intervals with the largest errors of each pending value;
    # an interval that can no longer be bisected ends that value's refinement
    largest <- per_value(error, max)
    split <- error > 0 & error >= largest[index] / 8
    middle <- lower[split] / 2 + upper[split] / 2
    stuck <- middle <= lower[split] | middle >= upper[split]
    if(any(stuck)){
      active[unique(index[split][stuck])] <- FALSE
      split[which(split)[stuck]] <- FALSE
      middle <- middle[!stuck]
    }
    if(!any(split)){
      break
    }
    split_at <- which(split)
    left <- list(lower = lower[split_at], upper = middle)
    right <- list(lower = middle, upper = upper[split_at])
    halves <- .prior_density_quadrature_evaluate(
      integrand,
      c(left$lower, right$lower), c(left$upper, right$upper),
      rep(kind[split_at], 2L), rep(anchor[split_at], 2L), rep(index[split_at], 2L),
      rules
    )
    n_split <- length(split_at)
    upper[split_at] <- middle
    value[split_at] <- halves$value[seq_len(n_split)]
    error[split_at] <- halves$error[seq_len(n_split)]
    lower <- c(lower, right$lower)
    upper <- c(upper, right$upper)
    kind <- c(kind, kind[split_at])
    anchor <- c(anchor, anchor[split_at])
    index <- c(index, index[split_at])
    value <- c(value, halves$value[n_split + seq_len(n_split)])
    error <- c(error, halves$error[n_split + seq_len(n_split)])
  }
  out
}
