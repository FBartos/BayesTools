# Integration of the ordered-product leaf with breakpoints for a total with a
# bounded support (R/prior-density-ordered-product.R).
#
# A single random total T with a finite support bound b has a conditional
# density f_T((x - s(B)) / c(B)) / |c(B)| with jumps, kinks, or integrable
# singularities in the shares B where the argument crosses a bound, and the
# probability of a region has such features where the argument of a region
# endpoint does; consecutive tensor rules miss them alike. The coefficient c(B)
# and the shift s(B) are multilinear in the shares, so for fixed other
# coordinates they are affine in one coordinate v, c(v) = c0 + (c1 - c0) v and
# s(v) likewise, and the points where the argument hits a bound,
#   x - s(v) = b c(v)  =>  v = (x - s0 - b c0) / ((s1 - s0) + b (c1 - c0)),
# and where c changes sign, v = -c0 / (c1 - c0), are closed form. The integral
# over the innermost coordinate is split there. The integral over an outer
# coordinate u is non-smooth where such a point of the inner coordinates reaches
# an end of [0, 1], that is where the same equations hold at a corner of the
# inner coordinates (an affine equation in u), so the outer coordinates are
# split at these corners as well. A piece that reaches an end of [0, 1] is a
# Gauss-Jacobi rule with the Beta weight of that end, so every end of a smooth
# piece is exact. A piece next to a breakpoint where the density of the total is
# infinite (a Beta total with a shape below 1, a gamma total with shape below 1)
# has an algebraic singularity dist^(e - 1) there, e the exponent of the
# density at the bound (known in closed form), which is not of the pure form
# that a Jacobi weight absorbs for the probabilities of regions and for the
# outer coordinates; the end is smoothed by the substitution tau = s^q with
# q e >= 6, a Gauss-Legendre rule in s. Totals whose density at a bound is
# infinite with an unknown or very small exponent are refused.
#
# The nodes next to a singular breakpoint that the floating-point resolution of
# its location does not resolve leave out a mass that is bounded by the
# envelope of the density of the total at their exact distances from the bound
# (.prior_ordered_product_unresolved_bias()); a value whose bound exceeds the
# criterion, or has no bound, is rejected, never set to 0.
#
# The rules are doubled, and a value is accepted only when the last two rules
# agree at the criterion, the difference at least halved from the one before,
# and the two rules before agreed at 100 times the criterion (the differences of
# an accepted sequence decrease); other values are rejected.

# Whether the total of the single random term has a finite support bound.
.prior_ordered_product_bounded_total <- function(spec){

  identical(spec$mode, "scale") &&
    any(is.finite(unlist(spec$terms[[spec$random_terms[[1L]]]]$total$truncation[c("lower", "upper")],
                         use.names = FALSE)))
}

# The most allocation coordinates of a product with a bounded total that the
# nested pieces integrate (their number grows by a factor of several per
# coordinate).
.prior_ordered_product_piece_coordinates <- function(){

  3L
}

# The smallest exponent of the density of a total at a bound that the pieces
# integrate (the density is ~ dist^(e - 1) there): the nodes of a piece next to
# a smaller exponent resolve distances of 1e-16 and neglect (1e-16)^e of the
# mass.
.prior_ordered_product_exponent_limit <- function(){

  0.5
}

# The exponent e of the density of a simple continuous prior at the finite
# bound 'bound' of its support (the density is ~ dist^(e - 1) there): the shape
# of the Beta and gamma families at their natural bounds, 1 for a finite
# positive density or a density that vanishes faster than any power, and NA
# for an infinite density of unknown exponent or an unclassified behavior.
.prior_ordered_product_bound_exponent <- function(prior, bound){

  parameters <- prior$parameters
  if(identical(prior$distribution, "beta")){
    if(bound == 0) return(parameters$alpha)
    if(bound == 1) return(parameters$beta)
  }
  if(identical(prior$distribution, "gamma") && bound == 0){
    return(parameters$shape)
  }
  behavior <- .prior_density_ordinate_continuous_behavior(.prior_density_ordinate_primitive(prior, bound))
  if(behavior %in% c("regular", "zero")) 1 else NA_real_
}

# The exponents of the finite bounds of the total of the single random term
# (a named numeric, by bound), NULL for a total without a finite bound.
.prior_ordered_product_bound_exponents <- function(total){

  bounds <- unlist(total$truncation[c("lower", "upper")], use.names = FALSE)
  bounds <- bounds[is.finite(bounds)]
  stats::setNames(vapply(bounds, function(bound){
    .prior_ordered_product_bound_exponent(total, bound)
  }, numeric(1)), bounds)
}

# Whether the density ~ dist^(e - 1) at a bound is not analytic there (an
# exponent below 1, an infinite density, or a non-integer exponent below 4, a
# vanishing density with an infinite derivative or higher derivative), which
# the substitution tau = s^q with q e >= 6 smooths.
.prior_ordered_product_needs_smoothing <- function(exponent){

  !is.na(exponent) & (exponent < 1 | (abs(exponent - round(exponent)) > 1e-12 & exponent < 4))
}

# Rows of a matrix sorted increasingly (columns compared elementwise), with the
# companion matrices in the list 'companions' permuted alike.
.prior_ordered_product_sort_rows <- function(x, companions){

  k <- ncol(x)
  for(pass in seq_len(k - 1L)){
    for(j in seq_len(k - pass)){
      swap <- x[, j] > x[, j + 1L]
      low <- pmin(x[, j], x[, j + 1L])
      x[, j + 1L] <- pmax(x[, j], x[, j + 1L])
      x[, j] <- low
      for(i in seq_along(companions)){
        kept <- companions[[i]][, j]
        companions[[i]][swap, j] <- companions[[i]][swap, j + 1L]
        companions[[i]][swap, j + 1L] <- kept[swap]
      }
    }
  }
  list(x = x, companions = companions)
}

# The breakpoints of a coordinate at one setting of the deeper coordinates: the
# candidate roots (a column each) of the equations of the bounds at the
# reference values of the groups and of the sign of the coefficient, from the
# states 'low' and 'high' of the coordinate at 0 and 1, with the exponent and
# the value of the bound of each root (1 and NA for the sign).
.prior_ordered_product_candidates <- function(low, high, references, bounds, exponents){

  dc <- high$scale - low$scale
  ds <- high$shift - low$shift
  values <- list()
  exponent <- numeric()
  bound <- numeric()
  for(r in seq_len(ncol(references))){
    for(k in seq_along(bounds)){
      values[[length(values) + 1L]] <-
        (references[, r] - low$shift - bounds[[k]] * low$scale) / (ds + bounds[[k]] * dc)
      exponent <- c(exponent, exponents[[k]])
      bound <- c(bound, bounds[[k]])
    }
  }
  values[[length(values) + 1L]] <- -low$scale / dc
  list(values = values, exponents = c(exponent, 1), bounds = c(bound, NA_real_))
}

# The sorted edges (from 0 to 1) of the pieces of a coordinate from the
# candidate roots, with the exponent and the value of the bound at each interior
# edge (NA at 0 and 1). A root outside the open interval (0, 1) is no
# breakpoint (the argument does not reach the bound there): it is moved to the
# end 1, where it bounds a piece of zero width, and carries no exponent.
.prior_ordered_product_finish_edges <- function(candidates){

  breaks <- do.call(cbind, candidates$values)
  exponents <- matrix(candidates$exponents, nrow = nrow(breaks), ncol = ncol(breaks), byrow = TRUE)
  bounds <- matrix(candidates$bounds, nrow = nrow(breaks), ncol = ncol(breaks), byrow = TRUE)
  outside <- !is.finite(breaks) | breaks <= 0 | breaks >= 1
  breaks[outside] <- 1
  exponents[outside] <- NA_real_
  bounds[outside] <- NA_real_
  sorted <- .prior_ordered_product_sort_rows(breaks, list(exponents, bounds))
  list(edges = cbind(0, sorted$x, 1),
       exponents = cbind(NA_real_, sorted$companions[[1L]], NA_real_),
       bounds = cbind(NA_real_, sorted$companions[[2L]], NA_real_))
}

# The nodes of the pieces of [0, 1] between the sorted 'edges' (a row of edges
# per integral, from 0 to 1) for a Beta(a, b) weight, with 'n' nodes per piece:
# the integral row ('row'), the node ('v') and the weight that multiplies the
# integrand there (the Beta density and the width of the piece included).
# 'exponents' are the exponents of the density of the total at the interior
# edges and 'edge_bounds' the bounds of its support that the edges belong to; a
# piece next to an edge with an exponent below 1 is smoothed there. A node of
# such a piece closer to the edge than a few units of the floating-point
# resolution of the edge is 'unresolved', and the exact (logarithmic) distance
# of the nodes of such a piece from its singular edges is returned for the
# bound of the mass they leave out (the log distance, the bound and the
# exponent of the left and of the right edge, NA without a singular edge).
.prior_ordered_product_split <- function(edges, exponents, edge_bounds, a, b, n){

  out <- list(row = integer(), v = numeric(), w = numeric(), unresolved = logical(),
              log_left = numeric(), log_right = numeric(), bound_left = numeric(), bound_right = numeric(),
              exponent_left = numeric(), exponent_right = numeric())
  log_beta <- lbeta(a, b)
  last <- ncol(edges)
  for(j in seq_len(last - 1L)){
    lower <- edges[, j]
    upper <- edges[, j + 1L]
    kept <- which(upper > lower)
    if(length(kept) == 0L){
      next
    }
    left_exponent <- if(j > 1L) exponents[kept, j] else rep(NA_real_, length(kept))
    right_exponent <- if(j + 1L < last) exponents[kept, j + 1L] else rep(NA_real_, length(kept))
    left_bound <- if(j > 1L) edge_bounds[kept, j] else rep(NA_real_, length(kept))
    right_bound <- if(j + 1L < last) edge_bounds[kept, j + 1L] else rep(NA_real_, length(kept))
    singular <- .prior_ordered_product_needs_smoothing(left_exponent) |
      .prior_ordered_product_needs_smoothing(right_exponent)
    reaches_end <- upper[kept] == 1
    # pieces without a singular breakpoint: Gauss-Jacobi with the Beta weight of an end
    for(right in c(TRUE, FALSE)){
      rows <- kept[!singular & reaches_end == right]
      if(length(rows) == 0L){
        next
      }
      left_parameter <- if(j == 1L) a else 1
      right_parameter <- if(right) b else 1
      rule <- .prior_ordered_beta_rule(n, left_parameter, right_parameter)
      width <- upper[rows] - lower[rows]
      v <- lower[rows] + outer(width, rule$u)
      weight <- (width * exp(lbeta(left_parameter, right_parameter) - log_beta)) *
        (if(j == 1L) width^(a - 1) else v^(a - 1)) *
        (if(right) width^(b - 1) else (1 - v)^(b - 1)) *
        matrix(rule$w, nrow = length(rows), ncol = n, byrow = TRUE)
      out$row <- c(out$row, rep(rows, times = n))
      out$v <- c(out$v, as.vector(v))
      out$w <- c(out$w, as.vector(weight))
      out$unresolved <- c(out$unresolved, rep(FALSE, length(rows) * n))
      for(field in c("log_left", "log_right", "bound_left", "bound_right", "exponent_left", "exponent_right")){
        out[[field]] <- c(out[[field]], rep(NA_real_, length(rows) * n))
      }
    }
    # pieces next to a singular breakpoint: tau = s^qL / (s^qL + (1 - s)^qR), Gauss-Legendre in s
    rows <- which(singular)
    if(length(rows) > 0L){
      rows_kept <- kept[rows]
      width <- upper[rows_kept] - lower[rows_kept]
      left_beta <- j == 1L
      right_beta <- upper[rows_kept] == 1
      smoothing <- function(exponent){
        ifelse(.prior_ordered_product_needs_smoothing(exponent), ceiling(6 / exponent), 1)
      }
      q_left <- if(left_beta) rep(smoothing(a), length(rows)) else smoothing(left_exponent[rows])
      q_right <- ifelse(right_beta, smoothing(b), smoothing(right_exponent[rows]))
      rule <- .prior_ordered_beta_rule(n, 1, 1)
      log_s <- matrix(log(rule$u), nrow = length(rows), ncol = n, byrow = TRUE)
      log_t <- matrix(log1p(-rule$u), nrow = length(rows), ncol = n, byrow = TRUE)
      left <- q_left * log_s
      right_term <- q_right * log_t
      log_d <- pmax(left, right_term) + log1p(exp(-abs(left - right_term)))
      log_tau <- left - log_d
      log_one_minus_tau <- right_term - log_d
      first <- log(q_left) + (q_left - 1) * log_s + right_term
      second <- log(q_right) + left + (q_right - 1) * log_t
      log_derivative <- pmax(first, second) + log1p(exp(-abs(first - second))) - 2 * log_d
      v <- lower[rows_kept] + width * exp(log_tau)
      log_v <- if(left_beta) log(width) + log_tau else log(v)
      log_one_minus_v <- ifelse(matrix(right_beta, nrow = length(rows), ncol = n),
                                log(width) + log_one_minus_tau, log1p(-v))
      weight <- width * matrix(rule$w, nrow = length(rows), ncol = n, byrow = TRUE) *
        exp((a - 1) * log_v + (b - 1) * log_one_minus_v - log_beta + log_derivative)
      # a node closer to a singular breakpoint than a few units of the resolution of
      # the edge is not resolved: its distance, and the density that is infinite or
      # vanishing there, are rounded (it may round onto the edge itself)
      resolution <- 8 * .Machine$double.eps
      singular_left <- !left_beta & .prior_ordered_product_needs_smoothing(left_exponent[rows])
      singular_right <- .prior_ordered_product_needs_smoothing(right_exponent[rows])
      near_left <- singular_left &
        width * exp(log_tau) <= resolution * pmax(abs(lower[rows_kept]), .Machine$double.xmin)
      near_right <- singular_right &
        width * exp(log_one_minus_tau) <= resolution * pmax(abs(upper[rows_kept]), .Machine$double.xmin)
      # the exact distances from the singular edges (the rounded nodes do not carry them)
      log_left <- log(width) + log_tau
      log_left[!singular_left, ] <- NA_real_
      log_right <- log(width) + log_one_minus_tau
      log_right[!singular_right, ] <- NA_real_
      per_node <- function(values) as.vector(matrix(values, nrow = length(rows), ncol = n))
      out$row <- c(out$row, rep(rows_kept, times = n))
      out$v <- c(out$v, as.vector(v))
      out$w <- c(out$w, as.vector(weight))
      out$unresolved <- c(out$unresolved, as.vector(near_left | near_right))
      out$log_left <- c(out$log_left, as.vector(log_left))
      out$log_right <- c(out$log_right, as.vector(log_right))
      out$bound_left <- c(out$bound_left, per_node(ifelse(singular_left, left_bound[rows], NA_real_)))
      out$bound_right <- c(out$bound_right, per_node(ifelse(singular_right, right_bound[rows], NA_real_)))
      out$exponent_left <- c(out$exponent_left, per_node(ifelse(singular_left, left_exponent[rows], NA_real_)))
      out$exponent_right <- c(out$exponent_right, per_node(ifelse(singular_right, right_exponent[rows], NA_real_)))
    }
  }
  out
}

# The sorted edges and their exponents of the coordinate after the fixed
# 'prefix' (a row per integral, with the group 'g' of target values in
# 'values'): the breakpoints at the 'deeper' coordinates' corners of the
# equations of the bounds and of the sign of the coefficient.
.prior_ordered_product_edges <- function(spec, prefix, g, values, bounds, exponents, deeper){

  rows <- nrow(prefix)
  corners <- as.matrix(expand.grid(rep(list(c(0, 1)), deeper), KEEP.OUT.ATTRS = FALSE))
  references <- values[g, , drop = FALSE]
  candidates <- list(values = list(), exponents = numeric(), bounds = numeric())
  for(corner in seq_len(nrow(corners))){
    tail <- matrix(corners[corner, ], nrow = rows, ncol = deeper, byrow = TRUE)
    found <- .prior_ordered_product_candidates(
      .prior_ordered_product_state(spec, cbind(prefix, 0, tail)),
      .prior_ordered_product_state(spec, cbind(prefix, 1, tail)),
      references, bounds, exponents
    )
    candidates$values <- c(candidates$values, found$values)
    candidates$exponents <- c(candidates$exponents, found$exponents)
    candidates$bounds <- c(candidates$bounds, found$bounds)
  }
  .prior_ordered_product_finish_edges(candidates)
}

# The expectation of 'integrand(group, c, s)' over the shares for the groups of
# target values in the rows of the matrix 'values' (NA entries ignored), by the
# nested pieces and the doubling rule, with the result of
# .prior_ordered_product_integrate(). With 'nested = FALSE' the outer
# coordinates are not split, and without 'stable' a value needs only the
# agreement of the last two rules (a plotted density).
.prior_ordered_product_piecewise <- function(spec, values, integrand,
                                              tolerance = .prior_linear_density_refinement_tolerance()$relative,
                                              budget = Inf, nested = TRUE, stable = TRUE,
                                              check = FALSE){

  coordinates <- spec$structure$coordinates
  p <- length(coordinates)
  m <- nrow(values)
  total <- spec$terms[[spec$random_terms[[1L]]]]$total
  exponents <- .prior_ordered_product_bound_exponents(total)
  bounds <- as.numeric(names(exponents))
  exponents <- unname(exponents)
  sizes <- .prior_ordered_product_sizes(p)
  value <- rep(NA_real_, m)
  estimate <- rep(NA_real_, m)
  change <- rep(NA_real_, m)
  accepted <- rep(FALSE, m)
  previous <- rep(NA_real_, m)
  previous_change <- rep(NA_real_, m)
  unresolved <- rep(FALSE, m)
  evaluations <- 0
  rules <- 0L
  message <- "OK"
  # the agreement of two earlier rules is required when three rules fit
  stable <- stable && length(sizes) >= 3L
  if(length(sizes) < 2L){
    message <- paste0("the node cap admits one rule only for ", p, " allocation coordinates")
  }else{
    for(size in sizes){
      active <- which(!accepted)
      if(length(active) == 0L){
        break
      }
      # about three pieces per outer coordinate and four for the innermost
      per_group <- as.numeric(size)^p * (if(nested) 3^(p - 1L) else 1) * 4
      if(per_group * length(active) > budget && rules >= 2L){
        break
      }
      rules <- rules + 1L
      current <- rep(NA_real_, length(active))
      chunk <- max(1L, floor(3e6 / per_group))
      bias <- rep(0, length(active))
      for(start in seq.int(1L, length(active), by = chunk)){
        groups <- active[start:min(length(active), start + chunk - 1L)]
        integrals <- .prior_ordered_product_nested(
          spec, values, bounds, exponents, integrand, groups, size, nested, check = check
        )
        current[match(groups, active)] <- integrals$value
        bias[match(groups, active)] <- integrals$bias
        evaluations <- evaluations + integrals$evaluations
      }
      difference <- abs(current - previous[active])
      estimate[active] <- current
      change[active] <- difference
      good <- is.finite(current) & current > 0 & is.finite(difference) &
        difference <= tolerance * current
      if(stable){
        # the differences of a convergent sequence decrease: at least halve (or
        # are at the rounding level) and follow an earlier difference that is small
        good <- good & is.finite(previous_change[active]) &
          previous_change[active] <= 100 * tolerance * current &
          (difference <= 0.5 * previous_change[active] | difference <= 1e-3 * tolerance * current)
      }
      # the unresolved nodes next to a singular bound could change the value by more than
      # the criterion: never certified (and never zeroed silently)
      representable <- is.finite(bias) & bias <= tolerance * abs(current)
      unresolved[active] <- !representable
      good <- good & representable
      accepted[active[good]] <- TRUE
      value[active[good]] <- current[good]
      previous_change[active] <- difference
      previous[active] <- current
    }
  }
  list(value = value, estimate = estimate, change = change, accepted = accepted,
       unresolved = unresolved, evaluations = evaluations, message = message, sizes = sizes)
}

# One rule size of the nested pieces for the 'groups': the integral of each
# group, the bound of the mass that its unresolved nodes leave out (with
# 'check'), and the evaluations of the integrand.
.prior_ordered_product_nested <- function(spec, values, bounds, exponents, integrand, groups, n, nested,
                                           check = FALSE){

  coordinates <- spec$structure$coordinates
  p <- length(coordinates)
  prefix <- matrix(0, length(groups), 0L)
  weight <- rep(1, length(groups))
  g <- groups
  flagged <- rep(0L, length(groups))
  for(k in seq_len(p - 1L)){
    if(nested){
      found <- .prior_ordered_product_edges(spec, prefix, g, values, bounds, exponents, p - k)
    }else{
      found <- list(edges = matrix(c(0, 1), nrow = nrow(prefix), ncol = 2L, byrow = TRUE),
                    exponents = matrix(NA_real_, nrow = nrow(prefix), ncol = 2L),
                    bounds = matrix(NA_real_, nrow = nrow(prefix), ncol = 2L))
    }
    split <- .prior_ordered_product_split(found$edges, found$exponents, found$bounds, coordinates[[k]]$a,
                                          coordinates[[k]]$b, n)
    prefix <- cbind(prefix[split$row, , drop = FALSE], split$v)
    weight <- weight[split$row] * split$w
    flagged <- flagged[split$row] + as.integer(split$unresolved)
    g <- g[split$row]
  }
  # the innermost coordinate: affine coefficients at its ends, breakpoints of
  # the bounds, and the integrand on the pieces
  inner <- coordinates[[p]]
  low <- .prior_ordered_product_state(spec, cbind(prefix, 0))
  high <- .prior_ordered_product_state(spec, cbind(prefix, 1))
  found <- .prior_ordered_product_finish_edges(.prior_ordered_product_candidates(
    low, high, values[g, , drop = FALSE], bounds, exponents
  ))
  split <- .prior_ordered_product_split(found$edges, found$exponents, found$bounds, inner$a, inner$b, n)
  dc <- high$scale - low$scale
  ds <- high$shift - low$shift
  scale <- low$scale[split$row] + dc[split$row] * split$v
  evaluated <- integrand(g[split$row], scale, low$shift[split$row] + ds[split$row] * split$v)
  # a node whose argument rounds onto a bound of the total has no density value
  near <- attr(evaluated, "unresolved", exact = TRUE)
  unresolved <- split$unresolved | (if(is.null(near)) FALSE else near)
  evaluated[!is.finite(evaluated)] <- 0
  inner_integral <- numeric(length(g))
  sums <- rowsum(split$w * evaluated, split$row)
  inner_integral[as.integer(rownames(sums))] <- sums[, 1L]
  totals <- rowsum(weight * inner_integral, g)
  value <- rep(NA_real_, length(groups))
  value[match(as.integer(rownames(totals)), groups)] <- totals[, 1L]
  bias <- rep(0, length(groups))
  if(check){
    bias <- .prior_ordered_product_unresolved_bias(
      spec, split, evaluated, unresolved, scale, dc[split$row], ds[split$row], weight, inner_integral,
      flagged, g, groups
    )
  }
  list(value = value, bias = bias, evaluations = length(split$row))
}

# The logarithm of an upper bound of the density of the total at a point inside
# its support at the distance 'distance' from its bound 'bound', where the
# density is ~ distance^(exponent - 1) (the families whose bound the pieces
# smooth: Beta at 0 and 1, gamma at 0); NA where no bound is established
# (another family, an exponent that is not the family's, a distance that is not
# small: the nodes it bounds are within a few units of the floating-point
# resolution of the bound, and a node farther from it belongs to another bound
# or to no singular breakpoint). The Beta density C t^(a - 1) (1 - t)^(b - 1)
# at the bound 0 has the distance d = t, and the factor (1 - t)^(b - 1) is
# monotone, at most 1 for b >= 1 and (1 - d)^(b - 1) for b < 1, so the density
# is at most C max(1, (1 - d)^(b - 1)) d^(a - 1) with C = 1 / (B(a, b) mass of
# the truncation); at the bound 1 the shapes are exchanged. The gamma density
# C t^(k - 1) exp(-rate t) is at most C t^(k - 1). The distance is enlarged
# for the monotone factor to cover the rounding of the distance, and the bound
# by a relative 1e-8 for the rounding of its evaluation.
.prior_ordered_product_density_envelope <- function(total, bound, exponent, distance){

  out <- rep(NA_real_, length(distance))
  parameters <- total$parameters
  log_mass <- tryCatch(
    if(.is_prior_default_range(total)) 0 else .prior_simple_log_C(total),
    error = function(error) NA_real_
  )
  if(!is.finite(log_mass)){
    return(out)
  }
  enlarged <- 2 * distance + 4 * .Machine$double.eps
  small <- is.finite(distance) & distance > 0 & enlarged < 1e-4
  if(identical(total$distribution, "beta")){
    at_zero <- !is.na(bound) & bound == 0
    at_one <- !is.na(bound) & bound == 1
    expected <- ifelse(at_zero, parameters$alpha, parameters$beta)
    other <- ifelse(at_zero, parameters$beta, parameters$alpha)
    valid <- small & (at_zero | at_one) & !is.na(exponent) & exponent == expected
    monotone <- pmax(0, (other - 1) * log1p(-pmin(enlarged, 0.25)))
    out[valid] <- (-lbeta(parameters$alpha, parameters$beta) - log_mass + monotone +
      (exponent - 1) * log(distance))[valid]
  }else if(identical(total$distribution, "gamma")){
    valid <- small & !is.na(bound) & bound == 0 & !is.na(exponent) & exponent == parameters$shape
    out[valid] <- (parameters$shape * log(parameters$rate) - lgamma(parameters$shape) - log_mass +
      (exponent - 1) * log(distance))[valid]
  }
  out + log1p(1e-8)
}

# A bound of the mass that the unresolved nodes (closer to a singular bound
# than a few units of the floating-point resolution of its location, or whose
# argument rounds onto the bound) leave out of an integral of a density, per
# group. The density at such a node is ~ distance^(e - 1) from the bound, which
# the rounded node cannot carry, and is bounded without sampling the integrand:
# by the envelope of its family at the exact distance of the node from the
# bound (.prior_ordered_product_density_envelope(); the distance in the
# argument is the distance in the coordinate times the slope |ds + bound dc| of
# the numerator of the argument over the coefficient there), the Beta density of
# the coordinate and the Jacobian of the smoothing being in the weight of the
# node. An unresolved node without such an envelope (a bound that is not a
# singular breakpoint of the piece, another family) makes the bound infinite.
# The unresolved nodes of the outer coordinates leave out the whole inner
# integral at their value, which is bounded by the computed inner integral plus
# the bound of its own unresolved nodes (the computed integral could be missing
# exactly that); a node is counted once for every outer coordinate that is
# unresolved. The caller rejects a value whose bound exceeds the criterion,
# never zeroes it silently (the representability limit of the double precision,
# e.g. a Beta(0.5, 0.5) total at a distance of 1e-10 from its bound in a
# product with uniform shares).
.prior_ordered_product_unresolved_bias <- function(spec, split, evaluated, unresolved, scale, dscale, dshift,
                                                    weight, inner_integral, flagged, g, groups){

  total <- spec$terms[[spec$random_terms[[1L]]]]$total
  n_rows <- length(g)
  lost <- numeric(n_rows)
  nodes <- which(unresolved)
  if(length(nodes) > 0L){
    log_distance <- matrix(NA_real_, length(nodes), 2L)
    bound <- exponent <- signed <- log_distance
    for(side in 1:2){
      suffix <- c("left", "right")[[side]]
      edge_bound <- split[[paste0("bound_", suffix)]][nodes]
      signed[, side] <- dshift[nodes] + edge_bound * dscale[nodes]
      log_distance[, side] <- log(abs(signed[, side])) + split[[paste0("log_", suffix)]][nodes] -
        log(abs(scale[nodes]))
      bound[, side] <- edge_bound
      exponent[, side] <- split[[paste0("exponent_", suffix)]][nodes]
    }
    log_distance[!is.finite(log_distance)] <- NA_real_
    # the bound that the node is closer to in the argument
    nearest <- ifelse(is.na(log_distance[, 1L]), 2L,
                      ifelse(is.na(log_distance[, 2L]), 1L,
                             ifelse(log_distance[, 1L] <= log_distance[, 2L], 1L, 2L)))
    pick <- cbind(seq_along(nodes), nearest)
    log_envelope <- .prior_ordered_product_density_envelope(
      total, bound[pick], exponent[pick], exp(log_distance[pick])
    ) - log(abs(scale[nodes]))
    # the argument t - bound = -(ds + bound dc) (v - v*) / c(v) has the sign of the side of the
    # breakpoint v* that the node lies on (right of the left edge of its piece, left of the right
    # edge): a node on the side of the bound outside the support has the exact density 0, so the
    # computed value alone is what it gets wrong
    known <- !is.na(bound[pick]) & !is.na(signed[pick]) & signed[pick] != 0
    argument_sign <- -sign(signed[pick]) * sign(scale[nodes]) * c(1, -1)[nearest]
    at_lower <- !is.na(bound[pick]) & bound[pick] == total$truncation$lower
    inside <- !known | ifelse(at_lower, argument_sign > 0, argument_sign < 0)
    # |computed - exact| is at most the larger of the two (both are non-negative)
    log_larger <- ifelse(inside, pmax(log_envelope, log(abs(evaluated[nodes]))), log(abs(evaluated[nodes])))
    log_larger[inside & is.na(log_envelope)] <- Inf
    # a node of zero weight adds nothing to the computed integral either
    node_lost <- numeric(length(nodes))
    positive <- split$w[nodes] > 0
    node_lost[positive] <- exp(log(split$w[nodes][positive]) + log_larger[positive])
    lost_by_row <- rowsum(node_lost, split$row[nodes])
    lost[as.integer(rownames(lost_by_row))] <- lost_by_row[, 1L]
  }
  # the inner nodes, and the outer unresolved nodes at the bound of their inner integral
  contribution <- ifelse(weight > 0, weight * lost, 0)
  outer <- flagged > 0L & weight > 0
  if(any(outer)){
    contribution[outer] <- contribution[outer] +
      weight[outer] * flagged[outer] * (inner_integral[outer] + lost[outer])
  }
  by_group <- rowsum(contribution, g)
  bias <- rep(0, length(groups))
  bias[match(as.integer(rownames(by_group)), groups)] <- by_group[, 1L]
  bias
}

# ---- entry points ------------------------------------------------------------------

# The integral of the density of X at the values 'x' (one group per value),
# with the rule doubling over the shares as .prior_ordered_product_integrate()
# or, for a total with a bounded support, with the nested pieces (only the
# innermost coordinate of a plotted density with three or more coordinates).
.prior_ordered_product_integrate_density <- function(spec, x,
                                                      tolerance = .prior_linear_density_refinement_tolerance()$relative,
                                                      budget = Inf, display = FALSE){

  if(.prior_ordered_product_bounded_total(spec)){
    nested <- !display || length(spec$structure$coordinates) <= 2L
    return(.prior_ordered_product_piecewise(
      spec, matrix(x, ncol = 1L),
      function(group, scale, shift){
        .prior_ordered_product_scale_density(spec, x[group], shift, scale)
      },
      tolerance = tolerance, budget = budget, nested = nested, stable = !display, check = TRUE
    ))
  }
  .prior_ordered_product_integrate(spec, function(state, columns){
    .prior_ordered_product_conditional_density(spec, state, x[columns])
  }, length(x), tolerance = tolerance, budget = budget)
}

# The probability of the region with the continuous 'intervals', as
# .prior_ordered_product_integrate_density().
.prior_ordered_product_integrate_region <- function(spec, intervals){

  if(.prior_ordered_product_bounded_total(spec)){
    endpoints <- unique(as.vector(intervals))
    endpoints <- endpoints[is.finite(endpoints)]
    return(.prior_ordered_product_piecewise(
      spec, matrix(if(length(endpoints)) endpoints else NA_real_, nrow = 1L),
      function(group, scale, shift){
        .prior_ordered_product_conditional_probability(
          spec, list(shift = shift, scale = scale), intervals
        )
      }
    ))
  }
  .prior_ordered_product_integrate(spec, function(state, columns){
    matrix(.prior_ordered_product_conditional_probability(spec, state, intervals),
           ncol = length(columns), nrow = length(state$shift))
  }, 1L)
}
