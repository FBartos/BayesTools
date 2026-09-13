#' @title Creates a prior distribution for a weight function
#'
#' @description \code{prior_weightfunction} creates a prior distribution for
#' fitting a RoBMA selection model. The \code{side} and \code{steps} arguments
#' define the p-value bins, and the \code{weights} argument defines the prior on
#' the publication weights in those bins.
#'
#' @param side side geometry. Either \code{"one-sided"} or \code{"two-sided"}.
#' @param steps increasing p-value cut points between 0 and 1.
#' @param weights a weight-prior object created by \code{wf_cumulative()},
#' \code{wf_fixed()}, or \code{wf_independent()}.
#' @param reference reference bin. Currently only \code{"most_significant"} is
#' supported and fixes the most significant bin to \code{omega = 1}.
#' @param prior_weights prior odds associated with a given distribution.
#' @param model a fixed selection-model specification from
#' \code{selection_model()}. Its choices do not add weight-prior parameters or
#' change \code{prior_weights}.
#'
#' @examples
#' p1 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05, .10),
#'   weights = wf_cumulative(alpha = c(1, 1, 1))
#' )
#'
#' p2 <- prior_weightfunction(
#'   side = "one-sided",
#'   steps = c(.05),
#'   weights = wf_independent(prior("beta", list(1, 1)))
#' )
#'
#' @return \code{prior_weightfunction} returns an object of class 'prior'.
#'
#' @export prior_weightfunction
#' @seealso [plot.prior()]
prior_weightfunction <- function(side = "one-sided", steps = c(.025, .05),
                                 weights = wf_cumulative(),
                                 reference = "most_significant",
                                 prior_weights = 1,
                                 model = selection_model()){

  check_char(side, "side")
  check_real(steps, "steps", check_length = 0)
  check_char(reference, "reference", allow_values = "most_significant")
  .check_prior_weight(prior_weights)
  check_selection_model(model)

  side <- .weightfunction_normalize_side(side)
  steps <- .weightfunction_validate_steps(steps)
  bins <- .weightfunction_bins(steps)
  weights <- .weightfunction_validate_weights(weights, n_bins = nrow(bins), reference = reference)
  truncation <- .weightfunction_weights_truncation(weights)

  output <- list(
    distribution  = "weightfunction",
    side          = side,
    steps         = steps,
    bins          = bins,
    reference     = reference,
    weights       = weights,
    model         = model,
    parameters    = list(steps = steps),
    truncation    = truncation,
    prior_weights = prior_weights
  )

  class(output) <- c("prior", "prior.weightfunction")

  return(output)
}

#' @order 1
#' @title Specifies source conditioning in a selection model
#'
#' @description \code{selection_model()} declares which sources a selection
#' model conditions on or integrates out, and how estimate-level
#' publication weights combine within a selection group. These are fixed model
#' choices, not parameters with prior distributions.
#'
#' @param estimate_random_effects either \code{"integrate"} (the default) or
#' \code{"condition"} for the estimate-level random-effect term. In a random
#' formula, this is the sole declared term whose grouping map has one distinct
#' level per retained estimate. No such term is required; multiple qualifying
#' terms are an error when the model binds its data.
#' @param other_random_effects either \code{"condition"} (the default) or
#' \code{"integrate"} for all remaining declared random-effect terms.
#' @param known_sampling_variance either \code{"integrate"} (the default) or
#' \code{"condition"} for the complete sampling-error vector with known
#' covariance. This choice also applies to univariate models and independent
#' sampling errors. Conditioning retains the full unknown realized error;
#' integration averages its full distribution before selection normalization.
#' @param weight_rule either \code{"product"}, the product of estimate weights,
#' or \code{"best"}, the weight of the bin containing the smallest p-value.
#' The best-p-value rule does not select the largest publication weight.
#' @param group an optional data-column reference: a bare or backticked column
#' name, or a single character string. It is captured without evaluating the
#' column and bound by the consuming model when its data are available.
#' \code{NULL} requests automatic group resolution. Derived groups must first
#' be stored as data columns.
#' @param prior a single prior object or \code{NULL}. Mixtures must be examined
#' one branch at a time.
#' @param model a selection-model specification to validate.
#' @param name argument name used in validation messages.
#' @param x a selection-model specification.
#' @param ... additional arguments, currently unused.
#'
#' @details The three source choices are independent. Conditioning retains
#' unknown latent effects during selection normalization; it does not make them
#' observed or fixed. Integration averages the source before normalization.
#' Sampling versus marginalizing a JAGS node is a separate computational choice.
#' The consuming model binds source roles and publication groups to its data, so
#' a prior can be constructed before data exist.
#'
#' Known sampling covariance \eqn{V} describes the complete error vector
#' \eqn{e \sim N(0,V)}. With positive sampling standard deviations it can be
#' written as \eqn{V = S R S} and \eqn{e = S z}, where \eqn{S} contains the
#' standard deviations and \eqn{z \sim N(0,R)}. Conditioning retains all of
#' \eqn{e}; there is no automatically integrated residual sampling component.
#' Covariance factorizations are computational representations and do not
#' define additional selection sources. P-values retain the original marginal
#' sampling standard deviations under either choice.
#'
#' If every nonzero source is conditioned on, the candidate result is
#' deterministic. Positive selection weights then cancel in normalization,
#' leaving the unselected distribution and no likelihood information about
#' those weights. A retained context with zero acceptance probability has no
#' normalized selected distribution and must not be silently discarded.
#'
#' Estimate-level identification follows the declared grouping map, not numerical
#' independence, coefficient supports, or posterior variances. A one-to-one term
#' with correlated known group covariance remains estimate-level and retains
#' that full covariance. Other terms are not split into coefficient families.
#' See [random_effects_level_roles()] for the shared formula resolver.
#'
#' Ordinary forwarding through wrapper arguments preserves the original column
#' reference. Wrappers can also explicitly forward with \code{group = {{group}}}.
#' To pass a column name stored in a variable, use
#' \code{do.call(selection_model, list(group = column_name))}. A bare symbol
#' always denotes a column, even if a vector with that name exists elsewhere.
#'
#' @return \code{selection_model()} returns an object of class
#' \code{selection_model} containing only the four choices and a character
#' column reference or \code{NULL}. \code{selection_model_spec()} returns the
#' validated specification from a weightfunction prior or the selection child of
#' a composed \code{prior_bias()}; it returns \code{NULL} for a prior without a
#' selection component. \code{check_selection_model()} and \code{print()}
#' invisibly return their validated input.
#'
#' @examples
#' model <- selection_model(
#'   estimate_random_effects = "integrate",
#'   other_random_effects = "condition",
#'   known_sampling_variance = "integrate",
#'   weight_rule = "best",
#'   group = paper_id
#' )
#' p <- prior_weightfunction(steps = .05, model = model)
#' selection_model_spec(p)
#'
#' @export
selection_model <- function(estimate_random_effects = "integrate",
                            other_random_effects = "condition",
                            known_sampling_variance = "integrate",
                            weight_rule = "product", group = NULL){

  group <- .selection_model_group_reference(rlang::enquo(group))
  output <- structure(list(
    estimate_random_effects = estimate_random_effects,
    other_random_effects    = other_random_effects,
    known_sampling_variance = known_sampling_variance,
    weight_rule             = weight_rule,
    group                   = group
  ), class = "selection_model")
  check_selection_model(output)
  output
}

#' @rdname selection_model
#' @order 2
#' @export
selection_model_spec <- function(prior){

  if(is.null(prior)){
    return(NULL)
  }
  if(!is.prior(prior)){
    stop("'prior' must be a prior object or NULL.", call. = FALSE)
  }
  if(is.prior.mixture(prior)){
    stop("A selection-model specification is unavailable for a mixture as a whole. Inspect each prior branch with 'selection_model_spec()'.", call. = FALSE)
  }
  if(is_prior_bias(prior)){
    prior <- prior[["selection"]]
    if(is.null(prior)){
      return(NULL)
    }
    if(!is.prior.weightfunction(prior)){
      stop("The selection component must be a weightfunction prior.", call. = FALSE)
    }
  }
  if(!is.prior.weightfunction(prior)){
    return(NULL)
  }
  check_selection_model(prior[["model"]], name = "prior$model")
  prior[["model"]]
}

#' @rdname selection_model
#' @order 3
#' @export
check_selection_model <- function(model, name = "model"){

  fields <- c("estimate_random_effects", "other_random_effects",
              "known_sampling_variance", "weight_rule", "group")
  if(!is.list(model) || !inherits(model, "selection_model") ||
     !identical(names(model), fields)){
    stop(paste0("'", name, "' must be a specification from 'selection_model()'."), call. = FALSE)
  }
  for(field in fields[1:3]){
    check_char(model[[field]], field,
               allow_values = c("condition", "integrate"), allow_NA = FALSE)
  }
  check_char(model[["weight_rule"]], "weight_rule",
             allow_values = c("product", "best"), allow_NA = FALSE)
  if(!is.null(model[["group"]])){
    check_char(model[["group"]], "group", allow_NA = FALSE)
    if(!nzchar(model[["group"]])){
      stop("'group' must name a non-empty data column.", call. = FALSE)
    }
  }
  invisible(model)
}

.selection_model_group_reference <- function(captured){

  seen <- list()
  while(is.symbol(rlang::quo_get_expr(captured))){
    expression <- rlang::quo_get_expr(captured)
    environment <- rlang::quo_get_env(captured)
    frame <- which(vapply(sys.frames(), identical, logical(1), y = environment))
    if(length(frame) != 1L ||
       !as.character(expression) %in% names(formals(sys.function(frame)))){
      break
    }
    if(any(vapply(seen, identical, logical(1), y = captured))){
      stop("'group' contains a circular wrapper argument reference. Forward a data-column name.", call. = FALSE)
    }
    seen[[length(seen) + 1L]] <- captured
    captured <- rlang::eval_bare(
      rlang::call2("enquo", expression, .ns = "rlang"), environment
    )
  }
  reference <- rlang::quo_get_expr(captured)
  if(is.symbol(reference)){
    reference <- as.character(reference)
  }
  if(is.null(reference)){
    return(NULL)
  }
  if(!is.character(reference) || length(reference) != 1L ||
     is.na(reference) || !nzchar(reference)){
    stop("'group' must be a data-column name, a single character string, or NULL.", call. = FALSE)
  }
  unname(reference)
}

#' @rdname prior_weightfunction
#' @param alpha positive cumulative-Dirichlet concentration parameters. If
#' omitted, a flat Dirichlet prior is used with one concentration parameter per
#' bin.
#' @export
wf_cumulative <- function(alpha = NULL){

  if(!is.null(alpha)){
    check_real(alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = 0, allow_NA = FALSE)
    if(any(!is.finite(alpha))){
      stop("The 'alpha' argument must be finite.", call. = FALSE)
    }
  }

  out <- list(type = "cumulative", alpha = alpha)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.cumulative")
  return(out)
}

#' @rdname prior_weightfunction
#' @param omega fixed non-negative relative publication weights, one per bin.
#' The reference-bin weight must be exactly 1.
#' @export
wf_fixed <- function(omega){

  check_real(omega, "omega", lower = 0, check_length = 0, allow_NA = FALSE)
  if(any(!is.finite(omega))){
    stop("The 'omega' argument must be finite.", call. = FALSE)
  }

  out <- list(type = "fixed", omega = omega)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.fixed")
  return(out)
}

#' @rdname prior_weightfunction
#' @param prior prior distribution for each non-reference weight.
#' @param scale latent scale for independent weights. \code{"omega"} places the
#' prior directly on the non-negative publication weight. \code{"log_omega"}
#' places the prior on \code{log(omega)} and transforms with
#' \code{omega = exp(log_omega)}, allowing weights above one whenever the log
#' prior assigns mass above zero.
#' @export
wf_independent <- function(prior, scale = "omega"){

  .check_prior(prior)
  if(!is.prior.simple(prior) || is.prior.point(prior) || is.prior.discrete(prior)){
    stop("'prior' must be a continuous simple prior distribution.", call. = FALSE)
  }

  check_char(scale, "scale", allow_values = c("omega", "log_omega", "log"))
  if(scale == "log"){
    scale <- "log_omega"
  }

  if(scale == "omega"){
    if(prior$truncation[["lower"]] < 0){
      stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
    }
  }

  out <- list(type = "independent", scale = scale, prior = prior)
  class(out) <- c("weightfunction_weights", "weightfunction_weights.independent")
  return(out)
}

.weightfunction_normalize_side <- function(side){

  side_clean <- .prior_clean_input_name(side)
  if(side_clean %in% c("onesided", "one")){
    return("one-sided")
  }
  if(side_clean %in% c("twosided", "two")){
    return("two-sided")
  }

  stop("'side' must be either 'one-sided' or 'two-sided'.", call. = FALSE)
}
.weightfunction_validate_steps <- function(steps){

  check_real(steps, "steps", check_length = 0, allow_NA = FALSE)

  if(length(steps) == 0){
    stop("'steps' must contain at least one p-value cut point.", call. = FALSE)
  }
  if(any(steps >= 1) || any(steps <= 0)){
    stop("'steps' must be higher than 0 and lower than 1.", call. = FALSE)
  }
  if(anyDuplicated(steps)){
    stop("'steps' must not contain duplicate cut points.", call. = FALSE)
  }
  if(!all(steps == cummax(steps))){
    stop("'steps' must be monotonically increasing.", call. = FALSE)
  }

  steps
}
.weightfunction_bins <- function(steps){

  cuts <- c(0, steps, 1)
  data.frame(
    lower     = cuts[-length(cuts)],
    upper     = cuts[-1],
    reference = c(TRUE, rep(FALSE, length(cuts) - 2L))
  )
}
.weightfunction_validate_weights <- function(weights, n_bins, reference){

  if(!inherits(weights, "weightfunction_weights")){
    stop("'weights' must be created by wf_cumulative(), wf_fixed(), or wf_independent().", call. = FALSE)
  }

  if(weights$type == "cumulative"){
    if(is.null(weights$alpha)){
      weights$alpha <- rep(1, n_bins)
    }
    check_real(weights$alpha, "alpha", lower = 0, allow_bound = FALSE, check_length = n_bins, allow_NA = FALSE)
    if(any(!is.finite(weights$alpha))){
      stop("The 'alpha' argument must be finite.", call. = FALSE)
    }

  }else if(weights$type == "fixed"){
    check_real(weights$omega, "omega", lower = 0, check_length = n_bins, allow_NA = FALSE)
    if(any(!is.finite(weights$omega))){
      stop("The 'omega' argument must be finite.", call. = FALSE)
    }
    if(reference == "most_significant" && !isTRUE(all.equal(weights$omega[1], 1))){
      stop("The reference-bin fixed weight must be exactly 1.", call. = FALSE)
    }

  }else if(weights$type == "independent"){
    .check_prior(weights$prior)
    if(!is.prior.simple(weights$prior) || is.prior.point(weights$prior) || is.prior.discrete(weights$prior)){
      stop("'weights$prior' must be a continuous simple prior distribution.", call. = FALSE)
    }
    if(weights$scale == "omega"){
      if(weights$prior$truncation[["lower"]] < 0){
        stop("Independent omega-scale weight priors must have non-negative support.", call. = FALSE)
      }
    }else if(weights$scale != "log_omega"){
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

  }else{
    stop("Unsupported weightfunction weight prior type.", call. = FALSE)
  }

  weights
}
.weightfunction_weights_truncation <- function(weights){

  if(weights$type == "cumulative"){
    return(list(lower = 0, upper = 1))
  }

  if(weights$type == "fixed"){
    return(list(lower = 0, upper = max(1, weights$omega, na.rm = TRUE)))
  }

  if(weights$type == "independent"){
    if(weights$scale == "omega"){
      upper <- weights$prior$truncation[["upper"]]
    }else if(weights$scale == "log_omega"){
      upper <- weights$prior$truncation[["upper"]]
      upper <- if(is.infinite(upper)) Inf else exp(upper)
    }else{
      stop("Unsupported independent weight prior scale.", call. = FALSE)
    }

    return(list(lower = 0, upper = max(1, upper)))
  }

  stop("Unsupported weightfunction weight prior type.", call. = FALSE)
}
.weightfunction_n_bins <- function(prior){

  if(!is.prior.weightfunction(prior)){
    stop("'prior' must be a weightfunction prior.", call. = FALSE)
  }

  nrow(prior$bins)
}
.weightfunction_local_cuts <- function(prior){

  c(prior$bins$lower[1], prior$bins$upper)
}
.weightfunction_alpha_marginal <- function(alpha, index){

  if(index <= 1L){
    return(list(type = "point", location = 1))
  }

  list(
    type  = "beta",
    alpha = sum(alpha[index:length(alpha)]),
    beta  = sum(alpha[seq_len(index - 1L)])
  )
}
.weightfunction_rng <- function(prior, n){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    out <- matrix(rep(prior$weights$omega, each = n), nrow = n)

  }else if(prior$weights$type == "cumulative"){
    theta <- extraDistr::rdirichlet(n, alpha = prior$weights$alpha)
    out   <- t(apply(theta[,J:1, drop = FALSE], 1, cumsum))[,J:1, drop = FALSE]
    # Keep the reference bin exact; row-wise cumulative sums can leave 1 +/- eps.
    out[,1L] <- 1

  }else if(prior$weights$type == "independent"){
    out <- matrix(1, nrow = n, ncol = J)
    if(J > 1L){
      draws <- matrix(rng(prior$weights$prior, n * (J - 1L)), nrow = n, ncol = J - 1L)
      if(prior$weights$scale == "log_omega"){
        draws <- exp(draws)
      }
      out[,2:J] <- draws
    }
  }

  colnames(out) <- paste0("omega[", seq_len(J), "]")
  out
}
.weightfunction_marginal_components <- function(prior){

  J <- .weightfunction_n_bins(prior)

  if(prior$weights$type == "fixed"){
    return(lapply(prior$weights$omega, function(x){
      list(type = "point", location = x)
    }))
  }

  if(prior$weights$type == "cumulative"){
    return(lapply(seq_len(J), function(j){
      .weightfunction_alpha_marginal(prior$weights$alpha, j)
    }))
  }

  if(prior$weights$type == "independent"){
    return(lapply(seq_len(J), function(j){
      if(j == 1L){
        list(type = "point", location = 1)
      }else{
        list(type = "prior", prior = prior$weights$prior, scale = prior$weights$scale)
      }
    }))
  }
}
.prior_weightfunction_component_range <- function(component, quantiles = .005){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      if(component$scale == "omega"){
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) mquant(component$prior, quantiles) else lower
        upper <- if(is.infinite(upper)) mquant(component$prior, 1 - quantiles) else upper
        c(lower, upper)
      }else{
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) 0 else exp(lower)
        upper <- if(is.infinite(upper)) exp(mquant(component$prior, 1 - quantiles)) else exp(upper)
        c(lower, upper)
      }
    }
  )
}
.weightfunction_range <- function(prior, quantiles = .005){

  ranges <- do.call(rbind, lapply(
    .weightfunction_marginal_components(prior),
    .prior_weightfunction_component_range,
    quantiles = quantiles
  ))

  x_range <- range(c(0, 1, as.vector(ranges)), finite = TRUE)
  if(x_range[1] == x_range[2]){
    x_range <- range(c(0, 1, x_range), finite = TRUE)
  }

  x_range
}
.prior_weightfunction_component_cdf <- function(component, q){

  switch(
    component$type,
    "point" = ppoint(q, location = component$location),
    "beta"  = stats::pbeta(q, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mcdf(component$prior, q)
      }else{
        p <- rep(NA_real_, length(q))
        q_known <- !is.na(q)
        p[q_known & q <= 0] <- 0
        inside <- q_known & q > 0
        if(any(inside)){
          p[inside] <- mcdf(component$prior, log(q[inside]))
        }
        p
      }
    }
  )
}
.prior_weightfunction_component_ccdf <- function(component, q){

  switch(
    component$type,
    "point" = ppoint(q, location = component$location, lower.tail = FALSE),
    "beta"  = stats::pbeta(
      q,
      shape1     = component$alpha,
      shape2     = component$beta,
      lower.tail = FALSE
    ),
    "prior" = {
      if(component$scale == "omega"){
        mccdf(component$prior, q)
      }else{
        p <- rep(NA_real_, length(q))
        q_known <- !is.na(q)
        p[q_known & q <= 0] <- 1
        inside <- q_known & q > 0
        if(any(inside)){
          p[inside] <- mccdf(component$prior, log(q[inside]))
        }
        p
      }
    }
  )
}
.prior_weightfunction_component_lpdf <- function(component, x){

  switch(
    component$type,
    "point" = dpoint(x, location = component$location, log = TRUE),
    "beta"  = stats::dbeta(x, shape1 = component$alpha, shape2 = component$beta, log = TRUE),
    "prior" = {
      if(component$scale == "omega"){
        mlpdf(component$prior, x)
      }else{
        out <- rep(-Inf, length(x))
        inside <- x > 0
        out[inside] <- mlpdf(component$prior, log(x[inside])) - log(x[inside])
        out
      }
    }
  )
}
.prior_weightfunction_component_quant <- function(component, p){

  switch(
    component$type,
    "point" = ifelse(is.na(p), NA_real_, component$location),
    "beta"  = stats::qbeta(p, shape1 = component$alpha, shape2 = component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mquant(component$prior, p)
      }else{
        exp(mquant(component$prior, p))
      }
    }
  )
}
.prior_weightfunction_component_mean <- function(component){

  switch(
    component$type,
    "point" = component$location,
    "beta"  = component$alpha / (component$alpha + component$beta),
    "prior" = {
      if(component$scale == "omega"){
        mean(component$prior)
      }else{
        stats::integrate(
          f     = function(x, prior) {
            y <- exp(x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
      }
    }
  )
}
.prior_weightfunction_component_var <- function(component){

  switch(
    component$type,
    "point" = 0,
    "beta"  = (component$alpha * component$beta) /
      ((component$alpha + component$beta)^2 * (component$alpha + component$beta + 1)),
    "prior" = {
      if(component$scale == "omega"){
        var(component$prior)
      }else{
        m1 <- .prior_weightfunction_component_mean(component)
        m2 <- stats::integrate(
          f     = function(x, prior) {
            y <- exp(2 * x) * pdf(prior, x)
            y[!is.finite(y)] <- 0
            y
          },
          lower = component$prior$truncation[["lower"]],
          upper = component$prior$truncation[["upper"]],
          prior = component$prior
        )$value
        m2 - m1^2
      }
    }
  )
}
