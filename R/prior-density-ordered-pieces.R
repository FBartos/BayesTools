# Integration of the ordered-product leaf with breakpoints for a total with a
# bounded support (R/prior-density-ordered-product.R).
#
# A single random total T with a finite support bound b has a conditional
# density f_T((x - s(B)) / c(B)) / |c(B)| with jumps or kinks in the shares B
# where the argument crosses a bound, and the probability of a region has kinks
# where the argument of a region endpoint does; consecutive tensor rules miss
# such features alike. The coefficient c(B) and the shift s(B) are multilinear
# in the shares, so for fixed other coordinates they are affine in one
# coordinate v, c(v) = c0 + (c1 - c0) v and s(v) likewise, and the points where
# the argument hits a bound,
#   x - s(v) = b c(v)  =>  v = (x - s0 - b c0) / ((s1 - s0) + b (c1 - c0)),
# and where c changes sign, v = -c0 / (c1 - c0), are closed form. The integral
# over the innermost coordinate is split there. The integral over an outer
# coordinate u is non-smooth where such a point of the inner coordinates reaches
# an end of [0, 1], that is where the same equations hold at a corner of the
# inner coordinates (an affine equation in u), so the outer coordinates are
# split at these corners as well. Pieces that reach an end of [0, 1] are
# Gauss-Jacobi rules with the Beta weight of that end, the other pieces
# Gauss-Jacobi rules without a weight and the Beta density as a factor, so every
# piece is a smooth integral. The rules are doubled and a value is accepted
# when the last two rules agree at the criterion and the two rules before them
# agreed at 100 times the criterion; other values are rejected.

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

# Rows of a matrix sorted increasingly (columns compared elementwise).
.prior_ordered_product_sort_rows <- function(x){

  k <- ncol(x)
  for(pass in seq_len(k - 1L)){
    for(j in seq_len(k - pass)){
      low <- pmin(x[, j], x[, j + 1L])
      x[, j + 1L] <- pmax(x[, j], x[, j + 1L])
      x[, j] <- low
    }
  }
  x
}

# The nodes of the pieces of [0, 1] between the sorted 'edges' (a row of edges
# per integral, from 0 to 1) for a Beta(a, b) weight, with 'n' nodes per piece:
# the integral row ('row'), the node ('v') and the weight that multiplies the
# integrand there (the Beta density and the width of the piece included).
.prior_ordered_product_split <- function(edges, a, b, n){

  out <- list(row = integer(), v = numeric(), w = numeric())
  log_beta <- lbeta(a, b)
  for(j in seq_len(ncol(edges) - 1L)){
    lower <- edges[, j]
    upper <- edges[, j + 1L]
    kept <- which(upper > lower)
    if(length(kept) == 0L){
      next
    }
    reaches_end <- upper[kept] == 1
    for(right in c(TRUE, FALSE)){
      rows <- kept[reaches_end == right]
      if(length(rows) == 0L){
        next
      }
      # the Beta weight of an end that the piece reaches; the factor of the
      # other end is a smooth function of the piece
      left_exponent <- if(j == 1L) a else 1
      right_exponent <- if(right) b else 1
      rule <- .prior_ordered_beta_rule(n, left_exponent, right_exponent)
      width <- upper[rows] - lower[rows]
      v <- lower[rows] + outer(width, rule$u)
      weight <- (width * exp(lbeta(left_exponent, right_exponent) - log_beta)) *
        (if(j == 1L) width^(a - 1) else v^(a - 1)) *
        (if(right) width^(b - 1) else (1 - v)^(b - 1)) *
        matrix(rule$w, nrow = length(rows), ncol = n, byrow = TRUE)
      out$row <- c(out$row, rep(rows, times = n))
      out$v <- c(out$v, as.vector(v))
      out$w <- c(out$w, as.vector(weight))
    }
  }
  out
}

# The sorted edges (from 0 to 1) of the pieces of the coordinate after the
# fixed 'prefix' (a row per integral, with the group 'g' of target values in
# 'values'): the breakpoints at the 'corners' (the values of the deeper
# coordinates, a row each) of the equations of the bounds and of the sign of the
# coefficient. Without 'corners' only the coordinate itself is split.
.prior_ordered_product_edges <- function(spec, prefix, g, values, bounds, deeper){

  rows <- nrow(prefix)
  corners <- as.matrix(expand.grid(rep(list(c(0, 1)), deeper), KEEP.OUT.ATTRS = FALSE))
  references <- values[g, , drop = FALSE]
  candidates <- list()
  for(corner in seq_len(nrow(corners))){
    tail <- matrix(corners[corner, ], nrow = rows, ncol = deeper, byrow = TRUE)
    low <- .prior_ordered_product_state(spec, cbind(prefix, 0, tail))
    high <- .prior_ordered_product_state(spec, cbind(prefix, 1, tail))
    dc <- high$scale - low$scale
    ds <- high$shift - low$shift
    for(r in seq_len(ncol(references))){
      for(bound in bounds){
        candidates[[length(candidates) + 1L]] <-
          (references[, r] - low$shift - bound * low$scale) / (ds + bound * dc)
      }
    }
    candidates[[length(candidates) + 1L]] <- -low$scale / dc
  }
  breaks <- do.call(cbind, candidates)
  breaks[!is.finite(breaks) | breaks <= 0 | breaks >= 1] <- 1
  cbind(0, .prior_ordered_product_sort_rows(breaks), 1)
}

# The expectation of 'integrand(group, c, s)' over the shares for the groups of
# target values in the rows of the matrix 'values' (NA entries ignored), by the
# nested pieces and the doubling rule, with the result of
# .prior_ordered_product_integrate(). With 'nested = FALSE' the outer
# coordinates are not split, and without 'stable' a value needs only the
# agreement of the last two rules (a plotted density).
.prior_ordered_product_piecewise <- function(spec, values, integrand,
                                              tolerance = .prior_linear_density_refinement_tolerance()$relative,
                                              budget = Inf, nested = TRUE, stable = TRUE){

  coordinates <- spec$structure$coordinates
  p <- length(coordinates)
  m <- nrow(values)
  total <- spec$terms[[spec$random_terms[[1L]]]]$total
  bounds <- unlist(total$truncation[c("lower", "upper")], use.names = FALSE)
  bounds <- bounds[is.finite(bounds)]
  sizes <- .prior_ordered_product_sizes(p)
  value <- rep(NA_real_, m)
  estimate <- rep(NA_real_, m)
  change <- rep(NA_real_, m)
  accepted <- rep(FALSE, m)
  previous <- rep(NA_real_, m)
  previous_change <- rep(NA_real_, m)
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
      for(start in seq.int(1L, length(active), by = chunk)){
        groups <- active[start:min(length(active), start + chunk - 1L)]
        integrals <- .prior_ordered_product_nested(
          spec, values, bounds, integrand, groups, size, nested
        )
        current[match(groups, active)] <- integrals$value
        evaluations <- evaluations + integrals$evaluations
      }
      difference <- abs(current - previous[active])
      estimate[active] <- current
      change[active] <- difference
      good <- is.finite(current) & current > 0 & is.finite(difference) &
        difference <= tolerance * current
      if(stable){
        good <- good & is.finite(previous_change[active]) &
          previous_change[active] <= 100 * tolerance * current
      }
      accepted[active[good]] <- TRUE
      value[active[good]] <- current[good]
      previous_change[active] <- difference
      previous[active] <- current
    }
  }
  list(value = value, estimate = estimate, change = change, accepted = accepted,
       evaluations = evaluations, message = message, sizes = sizes)
}

# One rule size of the nested pieces for the 'groups': the integral of each
# group and the evaluations of the integrand.
.prior_ordered_product_nested <- function(spec, values, bounds, integrand, groups, n, nested){

  coordinates <- spec$structure$coordinates
  p <- length(coordinates)
  prefix <- matrix(0, length(groups), 0L)
  weight <- rep(1, length(groups))
  g <- groups
  for(k in seq_len(p - 1L)){
    edges <- if(nested){
      .prior_ordered_product_edges(spec, prefix, g, values, bounds, p - k)
    }else{
      matrix(c(0, 1), nrow = nrow(prefix), ncol = 2L, byrow = TRUE)
    }
    split <- .prior_ordered_product_split(edges, coordinates[[k]]$a, coordinates[[k]]$b, n)
    prefix <- cbind(prefix[split$row, , drop = FALSE], split$v)
    weight <- weight[split$row] * split$w
    g <- g[split$row]
  }
  # the innermost coordinate: affine coefficients at its ends, breakpoints of
  # the bounds, and the integrand on the pieces
  inner <- coordinates[[p]]
  low <- .prior_ordered_product_state(spec, cbind(prefix, 0))
  high <- .prior_ordered_product_state(spec, cbind(prefix, 1))
  dc <- high$scale - low$scale
  ds <- high$shift - low$shift
  references <- values[g, , drop = FALSE]
  candidates <- list()
  for(r in seq_len(ncol(references))){
    for(bound in bounds){
      candidates[[length(candidates) + 1L]] <-
        (references[, r] - low$shift - bound * low$scale) / (ds + bound * dc)
    }
  }
  candidates[[length(candidates) + 1L]] <- -low$scale / dc
  breaks <- do.call(cbind, candidates)
  breaks[!is.finite(breaks) | breaks <= 0 | breaks >= 1] <- 1
  split <- .prior_ordered_product_split(
    cbind(0, .prior_ordered_product_sort_rows(breaks), 1), inner$a, inner$b, n
  )
  evaluated <- integrand(
    g[split$row],
    low$scale[split$row] + dc[split$row] * split$v,
    low$shift[split$row] + ds[split$row] * split$v
  )
  inner_integral <- numeric(length(g))
  sums <- rowsum(split$w * evaluated, split$row)
  inner_integral[as.integer(rownames(sums))] <- sums[, 1L]
  totals <- rowsum(weight * inner_integral, g)
  value <- rep(NA_real_, length(groups))
  value[match(as.integer(rownames(totals)), groups)] <- totals[, 1L]
  list(value = value, evaluations = length(split$row))
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
      tolerance = tolerance, budget = budget, nested = nested, stable = !display
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
