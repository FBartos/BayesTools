# Distributions of linear functionals of an R-squared budget prior.
#
# An R-squared budget gives the signal-to-noise ratio W = r / (1 - r),
# r ~ Beta(a, b), i.e. W ~ BetaPrime(a, b), and a Dirichlet allocation over
# n_units units, so that the share of one unit is phi ~ Beta(alpha, A - alpha),
# independent of W, when n_units > 1 (phi = 1 when n_units == 1). The signal of
# the unit is v = W * phi. A linear functional of the standardized coefficients
# of a term (coded columns d, scale c) has the law
#
#   exact   X = c * sqrt(v) * U1,  U1 one coordinate of a uniform direction on
#           the sphere of R^d (a fair sign for d = 1; for d >= 2 a symmetric
#           variable with U1^2 ~ Beta(1/2, (d - 1) / 2))
#   normal  X = c * sqrt(v) * Z,   Z ~ N(0, 1)
#   laplace X = c * sqrt(v) * Z,   Z Laplace with unit variance.
#
# Every law is symmetric about 0 and a product of independent factors, so that
# X = c * prod_k sqrt(S_k) * X0 with S_k ~ Beta(p_k, q_k) 'shrinkers' (the share
# of the unit, and the squared projection of the exact law with d >= 2) and a
# base X0 with a closed form (exact: sqrt(W) times a fair sign) or one 1-D
# integral over W (the kernel laws). The density and the tail probability of X
# are then nested one-dimensional integrals over the logarithms of the
# shrinkers (and of W), all computed deterministically by stats::integrate():
#
#   f_X(x) = int dbeta(s; p, q) f_X0(x / sqrt(s)) / sqrt(s) ds
#   P(X > x) = int dbeta(s; p, q) P(X0 > x / sqrt(s)) ds,     x > 0.
#
# The integrals run over l = log(s), which turns the integrands into smooth
# functions with O(1) features around the cut where x / sqrt(s) reaches the
# scale of X0 (a feature of s that is invisible to a quadrature in s itself for
# small x). Each integral is split at the locations of these features, and the
# integrand is evaluated on the log scale so that products of an overflowing
# weight and an underflowing tail stay finite.

.bt_R2_density_tolerance <- 1e-8

.bt_R2_density_stop_unavailable <- function(message){
  stop(structure(list(message = message, call = NULL),
    class = c("BayesTools_R2_density_unavailable", "BayesTools_R2_unavailable",
              "BayesTools_R2_condition", "error", "condition")))
}

# ---------------------------------------------------------------------------
# Specifications

.bt_R2_density_check_number <- function(value, name, count = FALSE){

  valid <- is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  if(valid && count){
    valid <- value >= 1 && value == round(value)
  }
  if(!valid){
    .bt_R2_stop(paste0("The R-squared functional specification needs '", name, "' to be ",
      if(count) "a positive integer" else "a single positive finite number", "."))
  }
  invisible(TRUE)
}

.bt_R2_functional_spec <- function(law, a, b, alpha = NA_real_, A = NA_real_, n_units = 1L,
                                   d = 1L, scale = 1){

  if(!is.character(law) || length(law) != 1L || !law %in% c("exact", "normal", "laplace")){
    .bt_R2_stop("The R-squared functional specification needs 'law' to be 'exact', 'normal' or 'laplace'.")
  }
  .bt_R2_density_check_number(a, "a")
  .bt_R2_density_check_number(b, "b")
  .bt_R2_density_check_number(scale, "scale")
  .bt_R2_density_check_number(n_units, "n_units", count = TRUE)
  .bt_R2_density_check_number(d, "d", count = TRUE)
  if(n_units > 1){
    .bt_R2_density_check_number(alpha, "alpha")
    .bt_R2_density_check_number(A, "A")
    if(alpha >= A){
      .bt_R2_stop(paste0("The R-squared functional specification needs 'alpha' to be smaller than ",
                         "the total concentration 'A'."))
    }
  }else{
    alpha <- NA_real_
    A     <- NA_real_
  }

  structure(
    list(law = law, a = a, b = b, alpha = alpha, A = A, n_units = as.integer(n_units),
         d = as.integer(d), scale = scale),
    class = "BayesTools_R2_functional")
}

# The Beta(p, q) factors S_k of the declared law, X = prod(sqrt(S_k)) * X0: the
# share of the unit (first, when there are several units) and the squared
# projection of the exact law with d >= 2.
.bt_R2_density_factors <- function(spec){

  c(if(spec$n_units > 1L) list(c(spec$alpha, spec$A - spec$alpha)),
    if(spec$law == "exact" && spec$d > 1L) list(c(1 / 2, (spec$d - 1) / 2)))
}

# The factors to integrate over and the base X0 (a fair sign times sqrt(W) for
# the exact law, sqrt(W) times the kernel otherwise). When a = A the share is
# absorbed exactly, BetaPrime(A, b) * Beta(alpha, A - alpha) = BetaPrime(alpha, b),
# so that the base takes the shape alpha and one level of integration is saved.
.bt_R2_density_structure <- function(spec){

  shrinkers <- .bt_R2_density_factors(spec)
  a         <- spec$a
  if(spec$n_units > 1L && spec$a == spec$A){
    a         <- spec$alpha
    shrinkers <- shrinkers[-1L]
  }
  list(base = list(type = spec$law, a = a, b = spec$b), shrinkers = shrinkers)
}

# ---------------------------------------------------------------------------
# Log density / log tail of the radial law at one x > 0

.bt_R2_density_softplus <- function(x){
  x * (x > 0) + log1p(exp(-abs(x)))
}

# The unit-variance kernel: log density and log upper tail at z >= 0 (z may be Inf).
.bt_R2_density_kernel <- function(type, z, what){

  if(type == "normal"){
    switch(what,
      pdf  = stats::dnorm(z, log = TRUE),
      tail = stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
  }else{
    switch(what,
      pdf  = log(sqrt(2) / 2) - sqrt(2) * z,
      tail = log(1 / 2) - sqrt(2) * z)
  }
}

# log P(B <= x), or log P(B > x), of B ~ Beta(p, q) at x = exp(log_x) <= 1 for a
# vector of log_x: pbeta() while x is a normal double; where x underflows, the
# lower tail x^p / (p B(p, q)) up to a relative error of order |q - 1| x (the
# next term of its series), used only where that error is below exp(-30);
# beyond that the probability is unavailable.
.bt_R2_density_log_pbeta <- function(log_x, p, q, lower.tail = TRUE){

  # pbeta() while x is a normal double; below, the leading term of the
  # series where its next term is negligible
  direct <- log_x >= log(.Machine$double.xmin)
  series <- !direct & (q == 1 | log(abs(q - 1)) + log_x < -30)
  if(any(!series & !direct)){
    .bt_R2_density_stop_unavailable(paste0(
      "The R-squared functional distribution needs a Beta(", format(p), ", ", format(q), ") probability below ",
      "the double range. Use less extreme shapes or evaluate the distribution by simulation."))
  }
  out   <- numeric(length(log_x))
  out[direct] <- stats::pbeta(exp(log_x[direct]), p, q, lower.tail = lower.tail, log.p = TRUE)
  # (the normalizer -log(p B(p, q)) = lgamma(p + q) - lgamma(q) - lgamma(p + 1)
  # is of order p; for a tiny p its first-order term p (digamma(q) -
  # digamma(1)), as p + q and p + 1 round; it is added to p log(x) before
  # anything larger; the complement by expm1())
  normalizer  <- if(p < 1e-8) p * (digamma(q) - digamma(1)) else (lgamma(p + q) - lgamma(q)) - lgamma(p + 1)
  lower       <- p * log_x[series] + normalizer
  out[series] <- if(lower.tail) lower else log(-expm1(lower))
  out
}

# The base of the exact law, sqrt(W) times a fair sign, at y >= 0 (vectorized):
# with r = y^2 / (1 + y^2), the density is y^(2a - 1) / {B(a, b) (1 + y^2)^(a + b)},
# P(|Y| <= y) = pbeta(r; a, b) and P(Y > y) = pbeta(1 - r; b, a) / 2 (written as
# 1 - pbeta(r; a, b) for y <= 1, where 1 - r is not accurate).
.bt_R2_density_exact_base <- function(base, y, what){

  out    <- rep(-Inf, length(y))
  finite <- is.finite(y) & y > 0
  ly     <- log(y[finite])
  far    <- y[finite] > 1
  # log r and log(1 - r), r = y^2 / (1 + y^2), without differences of the large
  # terms 2 log(y) and log(1 + y^2) for y > 1 (y^-2 from log(y), as y^2 can
  # overflow)
  inv2   <- exp(-2 * ly)
  lfar   <- log1p(inv2)
  log_r  <- ifelse(far, -lfar, 2 * ly - log1p(y[finite]^2))
  log_1r <- ifelse(far, -(2 * ly + lfar), -log1p(y[finite]^2))
  if(what == "pdf"){
    # y^(2a - 1) (1 + y^2)^-(a + b) = y^-(2b + 1) (1 + y^-2)^-(a + b) for y > 1,
    # where (a + b) log1p(y^-2) is (a + b) y^-2 once y^-2 nears the double
    # range (it can still be large for huge shapes)
    scaled <- ifelse(inv2 > 1e-200, (base$a + base$b) * lfar, exp(log(base$a + base$b) - 2 * ly))
    out[finite] <- ifelse(far, -(2 * base$b + 1) * ly - scaled,
                          (2 * base$a - 1) * ly + (base$a + base$b) * log_1r) - lbeta(base$a, base$b)
  }else if(what == "central"){
    out[finite] <- .bt_R2_density_log_pbeta(log_r, base$a, base$b)
    out[y == Inf] <- 0
  }else{
    out[finite] <- log(1 / 2) + ifelse(!far, .bt_R2_density_log_pbeta(log_r, base$a, base$b, lower.tail = FALSE),
                                       .bt_R2_density_log_pbeta(log_1r, base$b, base$a))
  }
  out
}

.bt_R2_density_logsumexp <- function(x){

  top <- max(x)
  if(is.finite(top)) top + log(sum(exp(x - top))) else top
}

# Merge candidate break points that lie within 1 of the previous kept one.
.bt_R2_density_breaks <- function(candidates, lower, upper){

  candidates <- sort.int(candidates[is.finite(candidates) & candidates > lower & candidates < upper])
  kept       <- numeric()
  for(candidate in candidates){
    if(length(kept) == 0L || candidate - kept[length(kept)] > 1){
      kept <- c(kept, candidate)
    }
  }
  kept
}

# The points at which a piece (lo, hi) of an integral is probed: its end
# points (one shifted by 1 when the other is infinite) and its midpoint.
.bt_R2_density_probes <- function(lo, hi){

  if(is.finite(lo) && is.finite(hi)){
    c(lo, (lo + hi) / 2, hi)
  }else if(is.finite(hi)){
    hi - c(0, 1)
  }else if(is.finite(lo)){
    lo + c(0, 1)
  }else{
    c(-1, 0, 1)
  }
}

# log of the integral of exp(log_f) over (lower, upper), split at 'breaks'. A
# piece whose largest probed value (the integrands are monotone towards the
# probed end points or smooth, see the break points of the callers) lies more
# than 50 below the largest of all pieces is dropped: it carries less than 1e-21
# of the integral (this removes the dead zones of the kernels, where the
# integrand underflows). Every other piece is a stats::integrate() call at the
# relative tolerance 'tolerance' on the integrand divided by its probed
# maximum, which also keeps the integrand away from underflow.
.bt_R2_density_integrate <- function(log_f, lower, upper, breaks, tolerance, level){

  nodes  <- c(lower, .bt_R2_density_breaks(breaks, lower, upper), upper)
  lo     <- nodes[-length(nodes)]
  hi     <- nodes[-1L]
  probes <- Map(.bt_R2_density_probes, lo, hi)
  values <- split(log_f(unlist(probes)), rep(seq_along(probes), lengths(probes)))
  if(anyNA(unlist(values)) || any(unlist(values) == Inf)){
    .bt_R2_density_stop_unavailable(paste0(
      "The R-squared functional density was rejected by diagnostics: the integrand over the ", level,
      " is not finite. Use less extreme shapes or evaluate the distribution by simulation."))
  }
  log_ref <- vapply(values, max, numeric(1))
  if(!is.finite(max(log_ref))){
    return(-Inf)
  }
  keep <- which(log_ref > max(log_ref) - 50)
  pieces <- vapply(keep, function(i){
    piece <- stats::integrate(function(t) exp(log_f(t) - log_ref[i]), lo[i], hi[i], rel.tol = tolerance,
                              abs.tol = 0, subdivisions = 200L, stop.on.error = FALSE)
    # the integrand reaches 1 at a probe of the piece, so a zero integral is
    # a feature that the quadrature missed, not a negligible piece
    if(piece$message != "OK" || !is.finite(piece$value) || piece$value <= 0){
      .bt_R2_density_stop_unavailable(paste0(
        "The R-squared functional density was rejected by diagnostics: the integration over the ",
        level, " reported '", piece$message, "'. Use less extreme shapes or evaluate the distribution by simulation."))
    }
    log_ref[i] + log(piece$value)
  }, numeric(1))
  .bt_R2_density_logsumexp(pieces)
}

# log of the density ('pdf'), the upper tail P(X > x) ('tail') or, for the exact
# law only (the unit signal's lower tail), the central probability P(|X| <= x)
# ('central') of prod(sqrt(S_k)) * X0 at one x > 0. The
# outermost shrinker S ~ Beta(p, q) is integrated over l = log(s); next to
# s = 1 a shape q < 1 is an algebraic singularity of the integrand in l, which
# the last piece (a break at -1 or later, so that it keeps s >= exp(-2) or
# more) removes exactly by the substitution 1 - s = w^(1 / q).
.bt_R2_density_eval <- function(base, shrinkers, x, what, tolerance){

  if(!is.finite(x)){
    return(if(what == "central") 0 else -Inf)
  }
  if(length(shrinkers) == 0L){
    return(.bt_R2_density_base(base, x, what, tolerance))
  }

  p        <- shrinkers[[1L]][1L]
  q        <- shrinkers[[1L]][2L]
  rest     <- shrinkers[-1L]
  shift    <- if(what == "pdf") -1 / 2 else 0
  log_norm <- lbeta(p, q)
  log_inner <- function(l){
    inner <- x * exp(-l / 2)
    if(length(rest) == 0L && base$type == "exact"){
      .bt_R2_density_exact_base(base, inner, what)
    }else{
      vapply(inner, function(z){
        .bt_R2_density_eval(base, rest, z, what, tolerance / 10)
      }, numeric(1))
    }
  }
  log_f  <- function(l) (p + shift) * l + (if(q == 1) 0 else (q - 1) * log(-expm1(l))) - log_norm + log_inner(l)
  breaks <- c(2 * log(x) - log(base$a / base$b), 2 * log(x), if(q > 1) log(p / (p + q - 1)))
  level  <- "log-share of a Beta factor"

  if(q >= 1){
    return(.bt_R2_density_integrate(log_f, -Inf, 0, breaks, tolerance, level))
  }
  breaks  <- .bt_R2_density_breaks(c(breaks, -1), -Inf, 0)
  cut     <- breaks[length(breaks)]
  log_f_w <- function(w){
    l <- log1p(-w^(1 / q))
    (p + shift - 1) * l - log(q) - log_norm + log_inner(l)
  }
  parts <- c(.bt_R2_density_integrate(log_f, -Inf, cut, breaks[-length(breaks)], tolerance, level),
             .bt_R2_density_integrate(log_f_w, 0, (-expm1(cut))^q, numeric(), tolerance, level))
  .bt_R2_density_logsumexp(parts)
}

# The base law at one y > 0.
.bt_R2_density_base <- function(base, y, what, tolerance){

  if(base$type == "exact"){
    return(.bt_R2_density_exact_base(base, y, what))
  }
  jacobian <- if(what == "pdf") -1 / 2 else 0
  log_f <- function(tau){
    .bt_R2_density_log_betaprime_at(tau, base$a, base$b) + jacobian * tau +
      .bt_R2_density_kernel(base$type, y * exp(-tau / 2), what)
  }
  # Beyond z_cut the kernel factor is below exp(-60) of its value at 0: the
  # integrand is dead for tau < 2 * log(y / z_cut).
  z_cut <- if(base$type == "normal") 11 else 42
  # tau = log(W) has its mode near log(a / b) and the standard deviation
  # sqrt(trigamma(a) + trigamma(b)); a concentrated total (large shapes) makes
  # that width much smaller than 1, the scale at which the break points
  # resolve features, so the integral runs over the standardized
  # t = (tau - centre) / width, with breaks at the bulk of the signal
  centre <- log(base$a / base$b)
  width  <- min(1, sqrt(trigamma(base$a) + trigamma(base$b)))
  log_f_t <- function(t) log_f(centre + width * t) + log(width)
  # (a ladder of breaks around the bulk keeps every piece next to it short)
  .bt_R2_density_integrate(log_f_t, -Inf, Inf,
                           c(0, c(-1, 1) %o% c(1.5, 3, 6, 12, 24), (c(2 * log(y), 2 * log(y / z_cut)) - centre) / width),
                           tolerance, "log signal-to-noise ratio")
}

# ---------------------------------------------------------------------------
# Two-level laws through the density of a log product
#
# In the two-level cases (the exact law of a block, d >= 2, of a model with
# several units; a kernel of a model with several units; a != A in both) the
# base depends on the two inner factors only through their product: with
# u = log(S1) + log(S2) (two Beta shrinkers for the exact law; the share and
# W for a kernel),
#
#   f_X(x) = int h(u) exp(-u / 2) k(x exp(-u / 2)) du,
#   P(X > x) = int h(u) K(x exp(-u / 2)) du,
#
# where k and K are the closed-form density and tail of sqrt(W) times a fair
# sign (exact law) or of the kernel, and h, the density of u, does not depend
# on x. h is a convolution of two log-Beta (or log-BetaPrime) densities,
# evaluated for a vector of u by double-exponential quadrature on pieces whose
# ends are its features (so that the endpoint singularities and the O(1)
# transitions are resolved), with the error estimated by halving the step. The
# outer integral is then one stats::integrate() call per x, instead of one per
# outer node. Values of u whose quadrature does not meet the tolerance are
# integrated by stats::integrate().

# Accepted difference between the log integrals of the steps h and 2 h: the
# double-exponential error is about squared by halving the step, so a
# difference of 1e-5 (the error of the step 2 h) leaves an error of order
# 1e-10 at the step h. Far below the double range of exp(), the difference is
# taken relative to the size of the log value, whose own rounding is larger.
.bt_R2_density_de_tolerance <- 1e-5

# Accepted share (on the log scale) of the integrand at the outermost nodes of
# a piece: a larger share means that the fixed limits of the rule truncate
# mass, whatever the step.
.bt_R2_density_de_truncation <- -30

# Nodes of a piece of length 'len' (a vector, one piece per row) whose left end
# is the origin of the distance d. With exponent s < 1 the distance is
# d = len * w^(1 / s), w in (0, 1), which turns an integrable singularity
# d^(s - 1) at the left end into a regular integrand in w; s = 1 is the plain
# tanh-sinh rule. Returns log(d), log(len - d) and the log weights (rows:
# pieces; columns: nodes), and the nodes of the step 2 h.
.bt_R2_density_de_piece <- function(len, s = 1, h = 1 / 8, limit = 4){

  t <- seq(-limit, limit, by = h)
  z <- pi * sinh(t)
  log_g <- stats::plogis(z, log.p = TRUE)
  log_g1 <- stats::plogis(-z, log.p = TRUE)
  log_len <- log(len)
  log_d_unit <- log_g / s
  log_rest_unit <- log(-expm1(log_d_unit))
  log_rest_unit[!is.finite(log_rest_unit)] <- log_g1[!is.finite(log_rest_unit)] - log(s)
  list(
    log_d      = outer(log_len, log_d_unit, "+"),
    log_rest   = outer(log_len, log_rest_unit, "+"),
    log_weight = outer(log_len, -log(s) + log_g / s + log_g1 + log(pi * cosh(t)) + log(h), "+"),
    coarse     = seq_along(t) %% 2L == 1L
  )
}

# Nodes of (start, Inf) (exp-sinh), one row per value of 'start'.
.bt_R2_density_de_tail <- function(start, h = 1 / 8, limit = 4.5){

  t <- seq(-limit, limit, by = h)
  z <- pi / 2 * sinh(t)
  list(
    e          = outer(start, exp(z), "+"),
    log_weight = matrix(z + log(pi / 2 * cosh(t)) + log(h), length(start), length(t), byrow = TRUE),
    coarse     = seq_along(t) %% 2L == 1L
  )
}

# log density of log(S) of S ~ Beta(p, q) at log(S) = -e (e > 0), from e, or
# from log(e) where e may underflow (log(1 - exp(-e)) = log(e) - e / 2 + O(e^2)).
.bt_R2_density_log_beta_at <- function(e, p, q, log_e = NULL){

  if(q == 1){
    return(-p * e - lbeta(p, q))
  }
  log_one_minus <- log(-expm1(-e))
  if(!is.null(log_e)){
    small <- log_e < -30
    log_one_minus[small] <- (log_e - e / 2)[small]
  }
  -p * e + (q - 1) * log_one_minus - lbeta(p, q)
}

# log density of tau = log(W), W ~ BetaPrime(a, b).
.bt_R2_density_log_betaprime_at <- function(tau, a, b){

  # centred at the mode tau0 = log(a / b): the log density at the mode plus
  # a (tau - tau0) - (a + b) (softplus(tau) - softplus(tau0)), whose
  # difference is log1p(p0 expm1(tau - tau0)) with p0 = a / (a + b) near the
  # mode, so that large shapes keep the relative precision of the density
  tau0 <- log(a / b)
  dt <- tau - tau0
  near <- abs(dt) <= 1
  # (keeps the shape of 'tau': vectors and matrices of nodes)
  difference <- dt
  difference[near] <- log1p(a / (a + b) * expm1(dt[near]))
  difference[!near] <- .bt_R2_density_softplus(tau[!near]) - .bt_R2_density_softplus(tau0)
  .bt_R2_density_log_betaprime_mode(a, b) + a * dt - (a + b) * difference
}

# The log density of log(W), W ~ BetaPrime(a, b), at its mode log(a / b):
# a log(a) + b log(b) - (a + b) log(a + b) - lbeta(a, b), which for large
# shapes is written with the Stirling corrections of the log-gamma functions
# (the direct form cancels terms of order a + b).
.bt_R2_density_log_betaprime_mode <- function(a, b){

  if(min(a, b) <= 10){
    # (log1p forms keep a small share of very unequal shapes; log(1 + x / y)
    # of a ratio beyond the double range is log(x) - log(y))
    log1p_ratio <- function(x, y){
      ratio <- x / y
      if(is.finite(ratio) && ratio < 1e300) log1p(ratio) else log(x) - log(y) + log1p(y / x)
    }
    return(-a * log1p_ratio(b, a) - b * log1p_ratio(a, b) - lbeta(a, b))
  }
  stirling <- function(x) 1 / (12 * x) - 1 / (360 * x^3) + 1 / (1260 * x^5)
  (log(a) + log(b) - log(a + b) - log(2 * pi)) / 2 - stirling(a) - stirling(b) + stirling(a + b)
}

# Row-wise log-sum-exp of a matrix (-Inf for rows without finite entries).
.bt_R2_density_logsumexp_rows <- function(m){

  m[!is.finite(m)] <- -Inf
  top <- .bt_R2_density_row_max(m)
  out <- rep(-Inf, nrow(m))
  ok <- is.finite(top)
  out[ok] <- top[ok] + log(rowSums(exp(m[ok, , drop = FALSE] - top[ok])))
  out
}

# Row maxima of a matrix without NA values.
.bt_R2_density_row_max <- function(m){
  m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
}

# Location (mean) and scale (SD) of log(S) of S ~ Beta(p, q), and of log(W) of
# W ~ BetaPrime(a, b).
.bt_R2_density_log_moments <- function(shapes, prime){

  if(prime){
    c(digamma(shapes[1L]) - digamma(shapes[2L]), sqrt(trigamma(shapes[1L]) + trigamma(shapes[2L])))
  }else{
    c(digamma(shapes[1L]) - digamma(sum(shapes)), sqrt(trigamma(shapes[1L]) - trigamma(sum(shapes))))
  }
}

# Break points of the distance e in (0, end) at the features of a factor
# located at e = centre with scale 'scale': the centre, and centre +- 1, 4 and
# 12 scales when the factor is narrower than the O(1) pieces (so that a
# concentrated factor is never missed between the nodes).
.bt_R2_density_feature_breaks <- function(centre, scale){

  steps <- if(scale < 1 / 2) c(0, -1, 1, -4, 4, -12, 12) else 0
  outer(centre, steps * scale, "+")
}

# log h(u) for a vector u: 'first' is the Beta shrinker c(p, q) of the share
# (or the squared projection), 'second' a Beta shrinker c(p, q) or, with
# 'prime', the BetaPrime shapes c(a, b) of W. The integral over the distance
# e from the end that the second factor (exact) or the share (kernels) makes
# singular is split at the features of both factors; its end pieces take the
# power substitution of their singular end, and its last piece is infinite for
# the kernels. Returns the log density and the indices whose quadrature fails
# its step-halving or truncation check.
.bt_R2_density_log_product <- function(u, first, second, prime){

  out <- rep(-Inf, length(u))
  inside <- which(is.finite(u) & (prime | u < 0))
  if(length(inside) == 0L){
    return(list(log = out, failed = integer()))
  }
  v <- u[inside]
  n <- length(v)
  first_moments <- .bt_R2_density_log_moments(first, FALSE)
  second_moments <- .bt_R2_density_log_moments(second, prime)
  if(prime){
    # e = tau - u: the share log(S) = -e, W = exp(u + e)
    end <- rep(Inf, n)
    centres <- cbind(.bt_R2_density_feature_breaks(rep(-first_moments[1L], n), first_moments[2L]),
                     .bt_R2_density_feature_breaks(second_moments[1L] - v, second_moments[2L]), 1)
    left_exponent <- min(first[2L], 1)
    log_f <- function(log_e, log_rest, e){
      .bt_R2_density_log_beta_at(e, first[1L], first[2L], log_e) +
        .bt_R2_density_log_betaprime_at(v + e, second[1L], second[2L])
    }
  }else{
    # e = l - u in (0, -u): log(S2) = -e, log(S1) = -(-u - e)
    end <- -v
    centres <- cbind(.bt_R2_density_feature_breaks(end + first_moments[1L], first_moments[2L]),
                     .bt_R2_density_feature_breaks(rep(-second_moments[1L], n), second_moments[2L]),
                     1, end - 1)
    left_exponent <- min(second[2L], 1)
    right_exponent <- min(first[2L], 1)
    log_f <- function(log_e, log_rest, e){
      .bt_R2_density_log_beta_at(e, second[1L], second[2L], log_e) +
        .bt_R2_density_log_beta_at(exp(log_rest), first[1L], first[2L], log_rest)
    }
  }
  # sorted interior break points of each row, padded by repeating the last
  # one (zero-length pieces carry no mass)
  centres[!(centres > 0 & centres < end)] <- NA
  breaks <- matrix(unlist(lapply(seq_len(n), function(i){
    row <- sort(unique(centres[i, !is.na(centres[i, ])]))
    if(length(row) == 0L) row <- if(prime) 1 else end[i] / 2
    c(row, rep(row[length(row)], ncol(centres) - length(row)))
  })), nrow = n, byrow = TRUE)

  fine <- coarse <- truncated <- matrix(-Inf, n, 0L)
  add <- function(values, weights, coarse_nodes){
    terms <- values + weights
    total <- .bt_R2_density_logsumexp_rows(terms)
    k <- ncol(terms)
    edge <- terms[, c(1:3, (k - 2):k), drop = FALSE]
    edge[!is.finite(edge)] <- -Inf
    edge <- .bt_R2_density_row_max(edge)
    fine <<- cbind(fine, total)
    coarse <<- cbind(coarse, .bt_R2_density_logsumexp_rows(terms[, coarse_nodes, drop = FALSE]) + log(2))
    truncated <<- cbind(truncated, ifelse(is.finite(total), edge - total, -Inf))
  }
  starts <- cbind(0, breaks)
  stops <- cbind(breaks, end)
  for(k in seq_len(ncol(starts))){
    lo <- starts[, k]
    hi <- stops[, k]
    if(k == ncol(starts) && prime){
      nodes <- .bt_R2_density_de_tail(lo)
      e <- nodes$e
      add(log_f(log(e), NULL, e), nodes$log_weight, nodes$coarse)
      next
    }
    len <- pmax(hi - lo, 0)
    if(k == 1L){
      # singular left end: distance from 0
      nodes <- .bt_R2_density_de_piece(len, left_exponent)
      log_e <- nodes$log_d
      e <- exp(log_e)
      log_rest <- log(end - e)
    }else if(k == ncol(starts)){
      # singular right end (exact): distance from the end
      nodes <- .bt_R2_density_de_piece(len, right_exponent)
      log_rest <- nodes$log_d
      e <- end - exp(log_rest)
      log_e <- log(e)
    }else{
      nodes <- .bt_R2_density_de_piece(len, 1)
      e <- lo + exp(nodes$log_d)
      log_e <- log(e)
      log_rest <- log(end - e)
    }
    add(log_f(log_e, log_rest, e), nodes$log_weight, nodes$coarse)
  }
  fine_total <- .bt_R2_density_logsumexp_rows(fine)
  coarse_total <- .bt_R2_density_logsumexp_rows(coarse)
  relevant <- fine > fine_total + .bt_R2_density_de_truncation
  truncation <- rowSums(relevant & truncated > .bt_R2_density_de_truncation) > 0
  out[inside] <- fine_total
  # Below the double range of exp() (log h < -745) a value carries no mass of
  # an integral of size one and needs only relative accuracy; the break points
  # bracket the modes of both factors, so that such a value is not a missed
  # peak (whose quadrature would not reach the double range).
  difference <- abs(fine_total - coarse_total)
  tolerance <- ifelse(fine_total < -745, 1e-3 * abs(fine_total), .bt_R2_density_de_tolerance)
  failed <- inside[!is.finite(fine_total) | truncation | difference > tolerance]
  list(log = out, failed = failed, difference = difference, truncation = truncation)
}

# log h(u) for a vector u, with stats::integrate() for the values whose
# double-exponential quadrature fails its checks (over the log distance from
# the singular end, so that no end of the range is singular).
.bt_R2_density_log_h <- function(u, first, second, prime, tolerance){

  quadrature <- .bt_R2_density_log_product(u, first, second, prime)
  out <- quadrature$log
  first_moments <- .bt_R2_density_log_moments(first, FALSE)
  second_moments <- .bt_R2_density_log_moments(second, prime)
  for(i in quadrature$failed){
    out[i] <- if(prime){
      log_f <- function(s){
        e <- exp(s)
        out <- s + .bt_R2_density_log_beta_at(e, first[1L], first[2L], s) +
          .bt_R2_density_log_betaprime_at(u[i] + e, second[1L], second[2L])
        # a share of exp(-e) with e beyond the double range has no mass
        out[!is.finite(e)] <- -Inf
        out
      }
      centres <- c(-first_moments[1L], second_moments[1L] - u[i])
      .bt_R2_density_integrate(log_f, -Inf, Inf, c(0, log(centres[centres > 0])), tolerance,
                               "log signal-to-noise ratio of the unit")
    }else{
      L <- -u[i]
      log_f <- function(s){
        log_g <- stats::plogis(s, log.p = TRUE)
        log_g1 <- stats::plogis(-s, log.p = TRUE)
        log(L) + log_g + log_g1 +
          .bt_R2_density_log_beta_at(L * exp(log_g1), first[1L], first[2L], log(L) + log_g1) +
          .bt_R2_density_log_beta_at(L * exp(log_g), second[1L], second[2L], log(L) + log_g)
      }
      centres <- c(L + first_moments[1L], -second_moments[1L])
      centres <- centres[centres > 0 & centres < L]
      .bt_R2_density_integrate(log_f, -Inf, Inf, c(0, stats::qlogis(centres / L)), tolerance,
                               "log share of the unit")
    }
  }
  out
}

# log density / tail / central probability of a two-level law at one x > 0
# through the density h of the log product (see above).
.bt_R2_density_eval_product <- function(base, shrinkers, x, what, tolerance){

  if(!is.finite(x)){
    return(if(what == "central") 0 else -Inf)
  }
  prime <- base$type != "exact"
  first <- shrinkers[[1L]]
  second <- if(prime) c(base$a, base$b) else shrinkers[[2L]]
  shift <- if(what == "pdf") -1 / 2 else 0
  log_kernel <- if(prime){
    function(z) .bt_R2_density_kernel(base$type, z, what)
  }else{
    unit <- list(type = "exact", a = base$a, b = base$b)
    function(z) .bt_R2_density_exact_base(unit, z, what)
  }
  # the exact base of a two-level law is sqrt(W) times a sign with the shapes
  # of the budget; the kernels' base is the kernel itself (W is in h)
  log_f <- function(u){
    .bt_R2_density_log_h(u, first, second, prime, tolerance / 10) + shift * u + log_kernel(x * exp(-u / 2))
  }
  centre <- log(first[1L] / sum(first)) + if(prime) log(second[1L] / second[2L]) else log(second[1L] / sum(second))
  breaks <- c(2 * log(x), centre, if(prime) 2 * log(x / if(base$type == "normal") 11 else 42) else 2 * log(x) - log(base$a / base$b))
  level <- "log signal of the unit"
  if(prime){
    return(.bt_R2_density_integrate(log_f, -Inf, Inf, breaks, tolerance, level))
  }
  # u < 0; next to 0 the density h behaves as (-u)^(Q - 1), Q = q1 + q2, which the
  # last piece removes for Q < 1 by the substitution -u = w^(1 / Q)
  Q <- first[2L] + second[2L]
  if(Q >= 1){
    return(.bt_R2_density_integrate(log_f, -Inf, 0, breaks, tolerance, level))
  }
  breaks <- .bt_R2_density_breaks(c(breaks, -1), -Inf, 0)
  cut <- breaks[length(breaks)]
  log_f_w <- function(w){
    u <- -w^(1 / Q)
    log_f(u) - log(Q) + (1 / Q - 1) * log(w)
  }
  parts <- c(.bt_R2_density_integrate(log_f, -Inf, cut, breaks[-length(breaks)], tolerance, level),
             .bt_R2_density_integrate(log_f_w, 0, (-cut)^Q, numeric(), tolerance, level))
  .bt_R2_density_logsumexp(parts)
}

# Whether a law has two levels of integration (see .bt_R2_density_eval_product()).
.bt_R2_density_two_level <- function(structure){

  length(structure$shrinkers) == 2L ||
    (length(structure$shrinkers) == 1L && structure$base$type != "exact")
}

# The spec-level evaluation at one y > 0 (scale 1). Two-level laws go through
# the density of the log product unless 'method' asks for the nested
# integrals (the independent reference of the tests).
.bt_R2_density_log <- function(spec, y, what, method = c("product", "nested")){

  method <- match.arg(method)
  structure <- .bt_R2_density_structure(spec)
  evaluate <- if(method == "product" && .bt_R2_density_two_level(structure)){
    .bt_R2_density_eval_product
  }else{
    .bt_R2_density_eval
  }
  tryCatch(
    evaluate(structure$base, structure$shrinkers, y, what, .bt_R2_density_tolerance),
    error = function(e){
      if(inherits(e, "BayesTools_R2_density_unavailable")){
        stop(e)
      }
      .bt_R2_density_stop_unavailable(paste0(
        "The R-squared functional density was rejected by diagnostics: the integration reported '",
        conditionMessage(e), "'. Use less extreme shapes or evaluate the distribution by simulation."))
    })
}

# ---------------------------------------------------------------------------
# Density at 0 (scale 1), from the shapes

# E[S^(-1/2)] of S ~ Beta(p, q), and of S = W ~ BetaPrime(p, q) (the 'kind').
.bt_R2_density_inverse_moment <- function(p, q, power, prime = FALSE){

  if(p <= power){
    return(Inf)
  }
  exp(lbeta(p - power, q + if(prime) power else 0) - lbeta(p, q))
}

# The density at 0 of a product of independent positive factors with Beta-type
# shapes p_k (left shape), critical constants crit_k (the density of the
# factor, in the 'power' scale, at 0+ when p_k equals the threshold) and
# inverse moments inv_k: infinite below the threshold or when several factors
# sit on it, 0 above it, and crit_k * prod(inv_j, j != k) when one factor sits
# on it.
.bt_R2_density_zero_rule <- function(p, crit, inv, threshold){

  p_min <- min(p)
  if(p_min < threshold){
    return(Inf)
  }
  if(p_min > threshold){
    return(0)
  }
  on_threshold <- which(p == threshold)
  if(length(on_threshold) > 1L){
    return(Inf)
  }
  crit[on_threshold] * prod(inv[-on_threshold])
}

.bt_R2_density_zero <- function(spec){

  shrinkers <- .bt_R2_density_factors(spec)
  inverse   <- c(.bt_R2_density_inverse_moment(spec$a, spec$b, 1 / 2, prime = TRUE),
                 vapply(shrinkers, function(s) .bt_R2_density_inverse_moment(s[1L], s[2L], 1 / 2), numeric(1)))
  if(spec$law != "exact"){
    kernel_zero <- if(spec$law == "normal") 1 / sqrt(2 * pi) else sqrt(2) / 2
    return(kernel_zero * prod(inverse))
  }
  shape    <- c(spec$a, vapply(shrinkers, `[`, numeric(1), 1L))
  critical <- c(1 / beta(1 / 2, spec$b), vapply(shrinkers, function(s) 1 / beta(1 / 2, s[2L]), numeric(1)))
  .bt_R2_density_zero_rule(shape, critical, inverse, 1 / 2)
}

# ---------------------------------------------------------------------------
# Functional distribution: pdf, cdf, quantile, rng, variance

.bt_R2_functional_pdf <- function(spec, x){

  exp(.bt_R2_functional_log_pdf(spec, x))
}

# The log density, finite also where the density underflows.
.bt_R2_functional_log_pdf <- function(spec, x){

  x   <- as.numeric(x)
  y   <- .bt_R2_density_standardize(abs(x), spec$scale)
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & is.infinite(x)] <- -Inf
  zero <- which(!is.na(x) & x == 0)
  if(length(zero) > 0L){
    out[zero] <- log(.bt_R2_density_zero(spec)) - log(spec$scale)
  }
  positive <- which(!is.na(x) & is.finite(x) & x != 0)
  out[positive] <- .bt_R2_density_unique_map(y[positive], function(z){
    .bt_R2_density_log(spec, z, "pdf")
  }) - log(spec$scale)
  out
}

# x / scale for the law of scale 1; a nonzero finite value whose quotient is
# not a positive finite double (underflow or overflow) is refused rather than
# taken for the origin or for infinity.
.bt_R2_density_standardize <- function(x, scale){

  y <- x / scale
  unrepresentable <- !is.na(x) & is.finite(x) & x != 0 & (y == 0 | !is.finite(y))
  if(any(unrepresentable)){
    .bt_R2_density_stop_unavailable(paste0(
      "The R-squared functional distribution is unavailable at values whose ratio to the scale of the ",
      "coefficient is not representable as a positive finite double (e.g., ", format(x[unrepresentable][1L]),
      " with scale ", format(scale), ")."))
  }
  y
}

# Evaluates 'f' once per distinct value of 'y' (values within a relative 1e-12
# are one value, e.g. the mirrored points of a display grid, whose rounding
# differs in the last bits), and maps the results back to 'y'.
.bt_R2_density_unique_map <- function(y, f){

  key <- signif(y, 12)
  distinct <- !duplicated(key)
  values <- vapply(y[distinct], f, numeric(1))
  values[match(key, key[distinct])]
}

.bt_R2_functional_cdf <- function(spec, q, lower.tail = TRUE){

  q     <- as.numeric(q)
  u     <- sign(q) * .bt_R2_density_standardize(abs(q), spec$scale)
  upper <- rep(NA_real_, length(u))
  upper[!is.na(u) & is.infinite(u)] <- 0
  upper[!is.na(u) & u == 0]         <- 1 / 2
  positive <- which(!is.na(u) & is.finite(u) & u != 0)
  upper[positive] <- .bt_R2_density_unique_map(abs(u[positive]), function(z){
    exp(.bt_R2_density_log(spec, z, "tail"))
  })
  # P(X > |u|) is the upper tail for u > 0 and the lower tail for u < 0
  below <- !is.na(u) & u < 0
  if(lower.tail){
    ifelse(below, upper, 1 - upper)
  }else{
    ifelse(below, 1 - upper, upper)
  }
}

# The positive x with P(X > x) = exp(log_tail), 0 < log_tail < log(1/2), of the
# scale-1 law; the bracket in log x is expanded from x = 1.
.bt_R2_density_root <- function(spec, log_tail){

  h <- function(u) .bt_R2_density_log(spec, exp(u), "tail") - log_tail
  h_one <- h(0)
  if(h_one == 0){
    return(1)
  }
  step <- if(h_one > 0) 2 else -2
  lower <- upper <- 0
  repeat{
    edge <- if(h_one > 0) upper + step else lower + step
    if(abs(edge) > 700){
      .bt_R2_density_stop_unavailable(paste0(
        "The R-squared functional quantile was rejected by diagnostics: the root of the tail ",
        "probability could not be bracketed within the double range. Use a less extreme probability."))
    }
    if(h_one > 0){
      lower <- upper
      upper <- edge
      if(h(upper) < 0) break
    }else{
      upper <- lower
      lower <- edge
      if(h(lower) > 0) break
    }
  }
  root <- tryCatch(stats::uniroot(h, c(lower, upper), tol = 1e-12)$root, error = function(e) e)
  if(inherits(root, "error")){
    .bt_R2_density_stop_unavailable(paste0(
      "The R-squared functional quantile was rejected by diagnostics: the root search reported '",
      conditionMessage(root), "'. Use a less extreme probability."))
  }
  exp(root)
}

.bt_R2_functional_quantile <- function(spec, p){

  p <- as.numeric(p)
  if(any(!is.na(p) & (p < 0 | p > 1))){
    .bt_R2_stop("The R-squared functional quantile needs probabilities in [0, 1].")
  }
  out <- rep(NA_real_, length(p))
  out[!is.na(p) & p == 0] <- -Inf
  out[!is.na(p) & p == 1] <- Inf
  out[!is.na(p) & p == 1 / 2] <- 0
  inner <- which(!is.na(p) & p > 0 & p < 1 & p != 1 / 2)
  out[inner] <- vapply(p[inner], function(prob){
    sign(prob - 1 / 2) * spec$scale * .bt_R2_density_root(spec, log(min(prob, 1 - prob)))
  }, numeric(1))
  out
}

.bt_R2_functional_rng <- function(spec, n){

  if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != round(n)){
    .bt_R2_stop("The R-squared functional generator needs 'n' to be a non-negative integer.")
  }
  r <- stats::rbeta(n, spec$a, spec$b)
  v <- r / (1 - r)
  if(spec$n_units > 1L){
    v <- v * stats::rbeta(n, spec$alpha, spec$A - spec$alpha)
  }
  sign <- 2 * stats::rbinom(n, 1L, 1 / 2) - 1
  u <- switch(spec$law,
    exact   = if(spec$d == 1L) sign else sign * sqrt(stats::rbeta(n, 1 / 2, (spec$d - 1) / 2)),
    normal  = stats::rnorm(n),
    laplace = sign * stats::rexp(n, rate = sqrt(2)))
  spec$scale * sqrt(v) * u
}

.bt_R2_functional_variance <- function(spec){

  signal     <- if(spec$b > 1) spec$a / (spec$b - 1) else Inf
  share      <- if(spec$n_units > 1L) spec$alpha / spec$A else 1
  projection <- if(spec$law == "exact") 1 / spec$d else 1
  spec$scale^2 * signal * share * projection
}

# ---------------------------------------------------------------------------
# Unit signal v = W * phi >= 0 (equal to the squared magnitude of the exact
# law with d = 1: v = X^2, so f_v(q) = f_X(sqrt(q)) / sqrt(q) and
# P(v <= q) = P(|X| <= sqrt(q))).

.bt_R2_signal_spec <- function(a, b, alpha = NA_real_, A = NA_real_, n_units = 1L){

  structure(
    list(functional = .bt_R2_functional_spec("exact", a, b, alpha, A, n_units, d = 1L, scale = 1)),
    class = "BayesTools_R2_signal")
}

.bt_R2_signal_pdf <- function(spec, v){

  v   <- as.numeric(v)
  out <- rep(NA_real_, length(v))
  out[!is.na(v) & (v < 0 | v == Inf)] <- 0
  zero <- which(!is.na(v) & v == 0)
  if(length(zero) > 0L){
    out[zero] <- .bt_R2_density_signal_zero(spec$functional)
  }
  positive <- which(!is.na(v) & is.finite(v) & v > 0)
  out[positive] <- vapply(v[positive], function(q){
    exp(.bt_R2_density_log(spec$functional, sqrt(q), "pdf") - log(q) / 2)
  }, numeric(1))
  out
}

.bt_R2_signal_cdf <- function(spec, v, lower.tail = TRUE){

  v   <- as.numeric(v)
  out <- rep(NA_real_, length(v))
  known <- !is.na(v)
  out[known & v <= 0] <- if(lower.tail) 0 else 1
  out[known & v == Inf] <- if(lower.tail) 1 else 0
  positive <- which(known & is.finite(v) & v > 0)
  out[positive] <- vapply(v[positive], function(q){
    if(lower.tail){
      exp(.bt_R2_density_log(spec$functional, sqrt(q), "central"))
    }else{
      2 * exp(.bt_R2_density_log(spec$functional, sqrt(q), "tail"))
    }
  }, numeric(1))
  out
}

# The density at 0 of v: its factors W ~ BetaPrime(a, b) and phi ~ Beta(alpha,
# A - alpha) have the critical shape 1, density b (resp. A - alpha) at 0 and
# E[1 / factor] = B(p - 1, q + prime) / B(p, q).
.bt_R2_density_signal_zero <- function(spec){

  shape    <- spec$a
  critical <- spec$b
  inverse  <- .bt_R2_density_inverse_moment(spec$a, spec$b, 1, prime = TRUE)
  if(spec$n_units > 1L){
    shape    <- c(shape, spec$alpha)
    critical <- c(critical, spec$A - spec$alpha)
    inverse  <- c(inverse, .bt_R2_density_inverse_moment(spec$alpha, spec$A - spec$alpha, 1))
  }
  .bt_R2_density_zero_rule(shape, critical, inverse, 1)
}
