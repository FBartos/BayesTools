#' Classify a prior-density ordinate
#'
#' @description
#' `prior_density_ordinate()` classifies the mathematical behavior of a prior at
#' one exact value. The classification uses the prior definition and
#' deterministic provenance recorded on linear-combination prior densities. It
#' does not use samples, kernel density estimates, numerical grids, or nearby
#' probe values to establish the behavior.
#'
#' @param x A BayesTools prior object or a `prior_linear_density` object produced
#'   by BayesTools' deterministic prior-density builders.
#' @param value One finite, non-missing numeric value at which to classify the
#'   prior ordinate.
#'
#' @return A named list of class `prior_density_ordinate` with fields:
#'
#' * `schema_version`: the stable result-schema version, currently `"1"`;
#' * `value`: the requested value;
#' * `behavior`: one of `"regular"`, `"zero"`, `"infinite"`,
#'   `"point_mass"`, `"undefined"`, or `"unknown"`;
#' * `log_density`: the continuous log-density ordinate, when available;
#' * `point_mass`: discrete probability at exactly `value`;
#' * `exact`: whether `behavior` follows from deterministic prior provenance;
#' * `method`: a machine-readable classification method;
#' * `reason`: `NULL` for an ordinary regular ordinate and a concise diagnostic
#'   otherwise;
#' * `provenance`: compact deterministic information used for classification.
#'
#' A `regular` classification is structural. Consequently, `log_density` may
#' be `-Inf` when an otherwise positive finite density underflows in ordinary
#' floating-point evaluation. General numerical convolutions, products, and
#' arbitrary user transformations are reported as `unknown` unless exact point
#' mass or singular-point metadata establishes the requested behavior.
#'
#' @examples
#' normal_prior <- prior("normal", list(mean = 0, sd = 1))
#' prior_density_ordinate(normal_prior, 0)
#'
#' gamma_prior <- prior("gamma", list(shape = 0.5, rate = 1))
#' prior_density_ordinate(gamma_prior, 0)
#'
#' mixture_prior <- prior_mixture(
#'   list(
#'     prior("point", list(location = 0), prior_weights = 1),
#'     prior("normal", list(mean = 0, sd = 1), prior_weights = 3)
#'   ),
#'   is_null = c(TRUE, FALSE)
#' )
#' prior_density_ordinate(mixture_prior, 0)
#'
#' @export
prior_density_ordinate <- function(x, value){

  check_real(value, "value", allow_NA = FALSE)
  if(!is.finite(value)){
    stop("The 'value' argument must be finite.", call. = FALSE)
  }
  value <- as.numeric(value)

  if(is.prior(x)){
    return(.prior_density_ordinate_prior(x, value))
  }
  if(inherits(x, "prior_linear_density")){
    return(.prior_density_ordinate_linear(x, value))
  }

  stop(
    "The 'x' argument must be a BayesTools prior or prior_linear_density object.",
    call. = FALSE
  )
}

.prior_density_ordinate_behaviors <- function(){
  c("regular", "zero", "infinite", "point_mass", "undefined", "unknown")
}

.prior_density_ordinate_methods <- function(){
  c(
    "primitive", "point", "finite_mixture", "scalar_affine",
    "linear_normal", "named_transform", "unsupported_provenance"
  )
}

.prior_density_ordinate_reason <- function(behavior){

  switch(
    behavior,
    "regular" = NULL,
    "zero" = "The continuous prior density is zero at the requested value.",
    "infinite" = paste0(
      "The continuous prior density tends to positive infinity at the ",
      "requested value."
    ),
    "point_mass" = paste0(
      "The prior assigns positive discrete probability to the requested value."
    ),
    "undefined" = paste0(
      "The prior-density ordinate is not mathematically defined at the ",
      "requested value."
    ),
    "unknown" = paste0(
      "The prior-density behavior is not available from supported ",
      "deterministic provenance."
    )
  )
}

.prior_density_ordinate_result <- function(value, behavior, log_density,
                                           point_mass = 0, exact = TRUE,
                                           method, reason = NULL,
                                           provenance = list(),
                                           continuous_behavior = behavior){

  if(!behavior %in% .prior_density_ordinate_behaviors()){
    stop("Unknown prior-density ordinate behavior.", call. = FALSE)
  }
  if(!is.character(method) || length(method) != 1L || is.na(method) ||
     !method %in% .prior_density_ordinate_methods()){
    stop("Unknown prior-density ordinate method.", call. = FALSE)
  }
  if(is.null(reason) && !identical(behavior, "regular")){
    reason <- .prior_density_ordinate_reason(behavior)
  }
  if(identical(behavior, "point_mass")){
    provenance$continuous_behavior <- continuous_behavior
  }

  out <- list(
    schema_version = "1",
    value       = as.numeric(value),
    behavior    = behavior,
    log_density = as.numeric(log_density),
    point_mass  = as.numeric(point_mass),
    exact       = isTRUE(exact),
    method      = as.character(method),
    reason      = reason,
    provenance  = provenance
  )
  class(out) <- c("prior_density_ordinate", "list")
  out
}

.prior_density_ordinate_continuous_behavior <- function(x){

  if(identical(x$behavior, "point_mass")){
    behavior <- x$provenance$continuous_behavior
    if(is.character(behavior) && length(behavior) == 1L &&
       behavior %in% .prior_density_ordinate_behaviors()){
      return(behavior)
    }
    return("unknown")
  }
  x$behavior
}

.prior_density_ordinate_compact <- function(x, max_length = 32L){

  if(is.null(x)){
    return(NULL)
  }
  if(is.function(x) || is.environment(x)){
    return("<omitted>")
  }
  if(is.expression(x) || is.call(x) || is.language(x)){
    return(paste(deparse(x, width.cutoff = 120L), collapse = ""))
  }
  if(is.atomic(x)){
    if(length(x) <= max_length){
      return(x)
    }
    return(list(
      length = length(x),
      first  = utils::head(x, 4L),
      last   = utils::tail(x, 4L)
    ))
  }
  if(is.list(x)){
    if(length(x) > max_length){
      return(list(
        length = length(x),
        first  = lapply(utils::head(x, 4L),
                        .prior_density_ordinate_compact,
                        max_length = max_length),
        last   = lapply(utils::tail(x, 4L),
                        .prior_density_ordinate_compact,
                        max_length = max_length)
      ))
    }
    out <- lapply(x, .prior_density_ordinate_compact, max_length = max_length)
    names(out) <- names(x)
    return(out)
  }
  as.character(class(x))
}

.prior_density_ordinate_prior_provenance <- function(prior){

  list(
    kind       = "primitive",
    family     = if(!is.null(prior$distribution)) prior$distribution else "unknown",
    parameters = .prior_density_ordinate_compact(prior$parameters),
    truncation = .prior_density_ordinate_compact(prior$truncation)
  )
}

.prior_density_ordinate_prior_definition_provenance <- function(prior){

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    weights <- .prior_density_ordinate_mixture_weights(prior)
    if(is.null(weights)){
      return(list(
        kind         = "unsupported_provenance",
        source_class = class(prior)
      ))
    }
    return(list(
      kind       = "finite_mixture",
      weights    = unname(weights),
      components = Map(function(component, weight){
        list(
          weight     = unname(weight),
          provenance =
            .prior_density_ordinate_prior_definition_provenance(component)
        )
      }, prior, weights)
    ))
  }
  if(is.prior.none(prior)){
    return(list(
      kind       = "primitive",
      family     = "none",
      location   = 0,
      truncation = c(lower = 0, upper = 0)
    ))
  }
  .prior_density_ordinate_prior_provenance(prior)
}

.prior_density_ordinate_unknown_prior <- function(prior, value, reason = NULL){

  if(is.null(reason)){
    reason <- .prior_density_ordinate_reason("unknown")
  }
  .prior_density_ordinate_result(
    value       = value,
    behavior    = "unknown",
    log_density = NA_real_,
    exact       = FALSE,
    method      = "unsupported_provenance",
    reason      = reason,
    provenance  = list(
      kind         = "unsupported_provenance",
      source_class = class(prior),
      family       = if(!is.null(prior$distribution)) prior$distribution else NULL
    )
  )
}

.prior_density_ordinate_prior <- function(prior, value){

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    return(.prior_density_ordinate_prior_mixture(prior, value))
  }
  if(is.prior.none(prior)){
    return(.prior_density_ordinate_atom_result(
      value      = value,
      locations  = 0,
      probability = 1,
      method     = "point",
      provenance = list(kind = "primitive", family = "none", location = 0)
    ))
  }
  if(is.prior.vector(prior) || is.prior.factor(prior) ||
     is.prior.weightfunction(prior) || is_prior_phacking(prior) ||
     is_prior_bias(prior) || is.prior.PET(prior) || is.prior.PEESE(prior)){
    return(.prior_density_ordinate_unknown_prior(prior, value))
  }
  if(is.prior.simple(prior)){
    return(.prior_density_ordinate_primitive(prior, value))
  }

  .prior_density_ordinate_unknown_prior(prior, value)
}

.prior_density_ordinate_parameters_numeric <- function(prior){

  parameters <- prior$parameters
  is.list(parameters) && length(parameters) > 0L &&
    all(vapply(parameters, function(parameter){
      is.numeric(parameter) && is.vector(parameter) && length(parameter) == 1L &&
        !is.na(parameter) && is.finite(parameter)
    }, logical(1)))
}

.prior_density_ordinate_atom_result <- function(value, locations, probability,
                                                method, provenance){

  keep <- is.finite(locations) & is.finite(probability) & probability > 0
  locations <- locations[keep]
  probability <- probability[keep]
  point_mass <- sum(probability[locations == value])

  if(point_mass > 0){
    return(.prior_density_ordinate_result(
      value               = value,
      behavior            = "point_mass",
      log_density         = -Inf,
      point_mass          = point_mass,
      exact               = TRUE,
      method              = method,
      provenance          = provenance,
      continuous_behavior = "zero"
    ))
  }

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "zero",
    log_density = -Inf,
    exact       = TRUE,
    method      = method,
    reason      = "The discrete prior has no mass at the requested value.",
    provenance  = provenance
  )
}

.prior_density_ordinate_primitive <- function(prior, value){

  provenance <- .prior_density_ordinate_prior_provenance(prior)
  family <- prior$distribution
  supported <- c(
    "normal", "lognormal", "t", "gamma", "invgamma", "beta", "exp",
    "uniform", "moment", "invmoment", "point", "bernoulli"
  )
  if(!is.character(family) || length(family) != 1L || !family %in% supported ||
     !.prior_density_ordinate_parameters_numeric(prior)){
    return(.prior_density_ordinate_unknown_prior(prior, value))
  }

  if(identical(family, "point")){
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = prior$parameters$location,
      probability = 1,
      method       = "point",
      provenance   = provenance
    ))
  }

  if(identical(family, "bernoulli")){
    discrete <- tryCatch(
      .prior_simple_truncated_discrete(prior),
      error = function(e) NULL
    )
    if(is.null(discrete)){
      return(.prior_density_ordinate_unknown_prior(prior, value))
    }
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = discrete$support,
      probability = discrete$prob,
      method       = "point",
      provenance   = provenance
    ))
  }

  truncation <- prior$truncation
  valid_truncation <- is.list(truncation) &&
    all(c("lower", "upper") %in% names(truncation)) &&
    all(vapply(truncation[c("lower", "upper")], function(bound){
      is.numeric(bound) && length(bound) == 1L && !is.na(bound)
    }, logical(1)))
  if(!valid_truncation){
    return(.prior_density_ordinate_unknown_prior(prior, value))
  }
  lower <- truncation$lower
  upper <- truncation$upper

  if(value < lower || value > upper){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "zero",
      log_density = -Inf,
      exact       = TRUE,
      method      = "primitive",
      reason      = "The requested value is outside the prior support.",
      provenance  = provenance
    ))
  }

  behavior <- switch(
    family,
    "normal" = "regular",
    "t" = "regular",
    "uniform" = "regular",
    "exp" = "regular",
    "lognormal" = if(value == 0) "zero" else "regular",
    "invgamma" = if(value == 0) "zero" else "regular",
    "gamma" = if(value == 0){
      if(prior$parameters$shape < 1) "infinite" else
        if(prior$parameters$shape == 1) "regular" else "zero"
    }else{
      "regular"
    },
    "beta" = if(value == 0){
      if(prior$parameters$alpha < 1) "infinite" else
        if(prior$parameters$alpha == 1) "regular" else "zero"
    }else if(value == 1){
      if(prior$parameters$beta < 1) "infinite" else
        if(prior$parameters$beta == 1) "regular" else "zero"
    }else{
      "regular"
    },
    "moment" = if(value == prior$parameters$location) "zero" else "regular",
    "invmoment" = if(value == prior$parameters$location) "zero" else "regular"
  )

  if(identical(behavior, "zero")){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = behavior,
      log_density = -Inf,
      exact       = TRUE,
      method      = "primitive",
      provenance  = provenance
    ))
  }
  if(identical(behavior, "infinite")){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = behavior,
      log_density = Inf,
      exact       = TRUE,
      method      = "primitive",
      provenance  = provenance
    ))
  }

  log_density <- tryCatch(
    .prior_simple_lpdf(prior, value),
    error = function(e) NA_real_
  )
  if(!is.numeric(log_density) || length(log_density) != 1L ||
     is.na(log_density) || (is.infinite(log_density) && log_density > 0)){
    log_density <- NA_real_
  }

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "regular",
    log_density = log_density,
    exact       = TRUE,
    method      = "primitive",
    reason      = if(identical(log_density, -Inf)){
      paste0(
        "The continuous prior density is structurally regular, but its ",
        "log-density underflows in ordinary floating-point arithmetic."
      )
    }else if(is.na(log_density)){
      paste0(
        "The continuous prior density is structurally regular, but its ",
        "log-density is unavailable in ordinary floating-point arithmetic."
      )
    }else{
      NULL
    },
    provenance  = provenance
  )
}

.prior_density_ordinate_mixture_weights <- function(prior){

  if(is.prior.spike_and_slab(prior)){
    inclusion <- tryCatch(
      mean(.get_spike_and_slab_inclusion(prior)),
      error = function(e) NA_real_
    )
    components <- attr(prior, "components", exact = TRUE)
    if(!is.numeric(inclusion) || length(inclusion) != 1L ||
       !is.finite(inclusion) || inclusion < 0 || inclusion > 1 ||
       !is.character(components) || length(components) != length(prior)){
      return(NULL)
    }
    weights <- ifelse(components == "alternative", inclusion, 1 - inclusion)
  }else{
    weights <- attr(prior, "prior_weights", exact = TRUE)
  }

  if(!is.numeric(weights) || length(weights) != length(prior) ||
     anyNA(weights) || any(!is.finite(weights)) || any(weights < 0) ||
     sum(weights) <= 0){
    return(NULL)
  }
  weights / sum(weights)
}

.prior_density_ordinate_log_sum <- function(log_density, weights){

  positive <- weights > 0
  log_density <- log_density[positive]
  weights <- weights[positive]
  if(length(log_density) == 0L){
    return(-Inf)
  }
  if(anyNA(log_density)){
    return(NA_real_)
  }
  terms <- log(weights) + log_density
  if(any(is.infinite(terms) & terms > 0)){
    return(Inf)
  }
  if(all(is.infinite(terms) & terms < 0)){
    return(-Inf)
  }
  maximum <- max(terms)
  maximum + log(sum(exp(terms - maximum)))
}

.prior_density_ordinate_combine <- function(results, weights, value,
                                            method = "finite_mixture",
                                            provenance_extra = list()){

  positive <- weights > 0
  results <- results[positive]
  weights <- weights[positive]
  weights <- weights / sum(weights)

  continuous <- vapply(
    results,
    .prior_density_ordinate_continuous_behavior,
    character(1)
  )
  continuous_behavior <- if(any(continuous == "undefined")){
    "undefined"
  }else if(any(continuous == "infinite")){
    "infinite"
  }else if(all(continuous %in% c("regular", "zero")) &&
           any(continuous == "regular")){
    "regular"
  }else if(all(continuous == "zero")){
    "zero"
  }else{
    "unknown"
  }

  point_mass <- sum(weights * vapply(results, `[[`, numeric(1), "point_mass"))
  log_density <- .prior_density_ordinate_log_sum(
    vapply(results, `[[`, numeric(1), "log_density"),
    weights
  )
  if(identical(continuous_behavior, "undefined")){
    log_density <- NA_real_
  }else if(identical(continuous_behavior, "infinite")){
    log_density <- Inf
  }else if(identical(continuous_behavior, "zero")){
    log_density <- -Inf
  }
  component_provenance <- Map(function(result, weight){
    list(
      weight              = unname(weight),
      behavior            = result$behavior,
      continuous_behavior = .prior_density_ordinate_continuous_behavior(result),
      point_mass           = result$point_mass,
      provenance           = result$provenance
    )
  }, results, weights)
  provenance <- c(
    list(
      kind       = "finite_mixture",
      weights    = unname(weights),
      components = component_provenance
    ),
    provenance_extra
  )

  if(point_mass > 0){
    return(.prior_density_ordinate_result(
      value               = value,
      behavior            = "point_mass",
      log_density         = log_density,
      point_mass          = point_mass,
      exact               = TRUE,
      method              = method,
      provenance          = provenance,
      continuous_behavior = continuous_behavior
    ))
  }

  component_exact <- vapply(results, `[[`, logical(1), "exact")
  exact <- if(identical(continuous_behavior, "undefined")){
    any(component_exact[continuous == "undefined"])
  }else if(identical(continuous_behavior, "infinite")){
    any(component_exact[continuous == "infinite"])
  }else if(continuous_behavior %in% c("regular", "zero")){
    all(component_exact)
  }else{
    FALSE
  }
  .prior_density_ordinate_result(
    value       = value,
    behavior    = continuous_behavior,
    log_density = log_density,
    exact       = exact,
    method      = method,
    provenance  = provenance
  )
}

.prior_density_ordinate_prior_mixture <- function(prior, value){

  weights <- .prior_density_ordinate_mixture_weights(prior)
  if(is.null(weights)){
    return(.prior_density_ordinate_unknown_prior(
      prior,
      value,
      "The finite-mixture weights are not available as exact numeric values."
    ))
  }
  results <- lapply(prior, .prior_density_ordinate_prior, value = value)
  .prior_density_ordinate_combine(results, weights, value)
}

.prior_density_ordinate_wrap <- function(source, value, log_jacobian,
                                         method, provenance){

  continuous_behavior <- .prior_density_ordinate_continuous_behavior(source)
  log_density <- source$log_density
  if(!is.na(log_density)){
    log_density <- log_density - log_jacobian
  }

  .prior_density_ordinate_result(
    value               = value,
    behavior            = source$behavior,
    log_density         = log_density,
    point_mass          = source$point_mass,
    exact               = source$exact,
    method              = method,
    reason              = source$reason,
    provenance          = provenance,
    continuous_behavior = continuous_behavior
  )
}

.prior_density_ordinate_prior_affine <- function(prior, value, offset, scale,
                                                 source_transform = NULL){

  provenance <- list(
    kind             = "scalar_affine",
    offset           = unname(offset),
    scale            = unname(scale),
    source_transform = source_transform,
    source           = if(is.prior(prior)){
      .prior_density_ordinate_prior_definition_provenance(prior)
    }else{
      list(kind = "unsupported_provenance")
    }
  )
  if(!is.finite(offset) || !is.finite(scale)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The scalar affine provenance is degenerate or non-finite.",
      provenance  = provenance
    ))
  }
  if(scale == 0){
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = offset,
      probability = 1,
      method       = "scalar_affine",
      provenance   = provenance
    ))
  }

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    weights <- .prior_density_ordinate_mixture_weights(prior)
    if(is.null(weights)){
      return(.prior_density_ordinate_unknown_prior(prior, value))
    }
    components <- lapply(prior, function(component){
      .prior_density_ordinate_prior_affine(
        component,
        value,
        offset,
        scale,
        source_transform
      )
    })
    mixture <- .prior_density_ordinate_combine(components, weights, value)
    continuous_behavior <-
      .prior_density_ordinate_continuous_behavior(mixture)
    mixture_provenance <- mixture$provenance
    mixture$method <- "scalar_affine"
    mixture$provenance <- provenance
    mixture$provenance$component_classification <-
      mixture_provenance$components
    if(identical(mixture$behavior, "point_mass")){
      mixture$provenance$continuous_behavior <- continuous_behavior
    }
    return(mixture)
  }

  if(is.prior.none(prior)){
    location <- 0
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = offset + scale * location,
      probability = 1,
      method       = "scalar_affine",
      provenance   = provenance
    ))
  }

  if(is.prior.simple(prior) &&
     prior$distribution %in% c("point", "bernoulli")){
    if(identical(prior$distribution, "point")){
      locations <- prior$parameters$location
      probability <- 1
    }else{
      discrete <- .prior_simple_truncated_discrete(prior)
      locations <- discrete$support
      probability <- discrete$prob
    }
    if(identical(source_transform, "log")){
      if(any(locations <= 0)){
        return(.prior_density_ordinate_result(
          value       = value,
          behavior    = "undefined",
          log_density = NA_real_,
          exact       = TRUE,
          method      = "named_transform",
          reason      = "The log source transformation is undefined for a prior atom.",
          provenance  = provenance
        ))
      }
      locations <- log(locations)
    }else if(!is.null(source_transform)){
      return(.prior_density_ordinate_unknown_prior(prior, value))
    }
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = offset + scale * locations,
      probability = probability,
      method       = "scalar_affine",
      provenance   = provenance
    ))
  }

  source_value <- (value - offset) / scale
  if(!is.finite(source_value)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = paste0(
        "The inverse affine value is not representable in ordinary ",
        "floating-point arithmetic."
      ),
      provenance  = provenance
    ))
  }

  if(is.null(source_transform)){
    source <- if(is.prior.simple(prior)){
      .prior_density_ordinate_primitive(prior, source_value)
    }else{
      .prior_density_ordinate_unknown_prior(prior, source_value)
    }
  }else if(identical(source_transform, "log")){
    source <- .prior_density_ordinate_log_source(prior, source_value)
  }else{
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The source transformation is not supported structurally.",
      provenance  = provenance
    ))
  }

  .prior_density_ordinate_wrap(
    source       = source,
    value        = value,
    log_jacobian = log(abs(scale)),
    method       = "scalar_affine",
    provenance   = provenance
  )
}

.prior_density_ordinate_log_source <- function(prior, value){

  provenance <- list(
    kind           = "named_transform",
    transformation = "log",
    arguments      = list(),
    source         = .prior_density_ordinate_prior_provenance(prior)
  )
  if(!is.prior.simple(prior) || is.prior.discrete(prior) ||
     is.prior.point(prior)){
    return(.prior_density_ordinate_unknown_prior(prior, value))
  }
  lower <- prior$truncation$lower
  upper <- prior$truncation$upper
  if(!is.numeric(lower) || length(lower) != 1L || is.na(lower) ||
     !is.numeric(upper) || length(upper) != 1L || is.na(upper)){
    return(.prior_density_ordinate_unknown_prior(prior, value))
  }
  if(lower < 0){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "undefined",
      log_density = NA_real_,
      exact       = TRUE,
      method      = "named_transform",
      reason      = "The log source transformation is not defined on the prior support.",
      provenance  = provenance
    ))
  }
  if((lower > 0 && value < log(lower)) ||
     (is.finite(upper) && value > log(upper))){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "zero",
      log_density = -Inf,
      exact       = TRUE,
      method      = "named_transform",
      reason      = "The requested value is outside the transformed prior support.",
      provenance  = provenance
    ))
  }

  original_value <- exp(value)
  if(!is.finite(original_value)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "regular",
      log_density = -Inf,
      exact       = TRUE,
      method      = "named_transform",
      reason      = paste0(
        "The continuous prior density is structurally regular, but its ",
        "log-density underflows in ordinary floating-point arithmetic."
      ),
      provenance  = provenance
    ))
  }
  source <- .prior_density_ordinate_primitive(prior, original_value)
  .prior_density_ordinate_wrap(
    source       = source,
    value        = value,
    log_jacobian = -value,
    method       = "named_transform",
    provenance   = provenance
  )
}

.prior_density_ordinate_point_group_location <- function(group, source_transforms){

  prior <- group$prior
  if(is.prior.none(prior)){
    return(0)
  }
  if(!is.prior.point(prior)){
    return(NULL)
  }
  location <- prior$parameters$location
  transforms <- source_transforms[names(group$weights)]
  transformed <- rep(location, length(group$weights))
  use_log <- !is.na(transforms) & transforms == "log"
  unsupported <- !is.na(transforms) & transforms != "log"
  if(any(unsupported) || (any(use_log) && location <= 0)){
    return(NA_real_)
  }
  transformed[use_log] <- log(location)
  sum(group$weights * transformed)
}

.prior_density_ordinate_linear_scalar <- function(prior_list, weights,
                                                  source_transforms, value){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }

  offset <- 0
  random_group <- NULL
  for(group in groups){
    point_location <- .prior_density_ordinate_point_group_location(
      group,
      source_transforms
    )
    if(length(point_location) == 1L){
      if(is.na(point_location)){
        return(NULL)
      }
      offset <- offset + point_location
      next
    }
    if(length(group$weights) != 1L || !is.null(random_group) ||
       is.prior.vector(group$prior) || is.prior.ordered(group$prior)){
      return(NULL)
    }
    random_group <- group
  }

  if(is.null(random_group)){
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = offset,
      probability = 1,
      method       = "scalar_affine",
      provenance   = list(
        kind    = "scalar_affine",
        offset  = offset,
        scale   = 0,
        weights = .prior_density_ordinate_compact(weights)
      )
    ))
  }

  parameter <- names(random_group$weights)[1L]
  .prior_density_ordinate_prior_affine(
    prior             = random_group$prior,
    value             = value,
    offset            = offset,
    scale             = unname(random_group$weights[[1L]]),
    source_transform  = .prior_linear_source_transform(source_transforms[parameter])
  )
}

.prior_density_ordinate_stable_norm <- function(x){

  maximum <- max(abs(x))
  if(maximum == 0){
    return(0)
  }
  maximum * sqrt(sum((x / maximum)^2))
}

.prior_density_ordinate_linear_normal <- function(prior_list, weights,
                                                  source_transforms, value){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }

  mean_terms <- numeric()
  scale_terms <- numeric()
  terms <- list()
  for(group in groups){
    prior <- group$prior
    transforms <- source_transforms[names(group$weights)]
    transforms[is.na(transforms)] <- "none"

    point_location <- .prior_density_ordinate_point_group_location(
      group,
      source_transforms
    )
    if(length(point_location) == 1L){
      if(is.na(point_location)){
        return(NULL)
      }
      mean_terms <- c(mean_terms, point_location)
      terms[[length(terms) + 1L]] <- list(
        parameter = group$parameter,
        family    = if(is.prior.none(prior)) "none" else prior$distribution,
        weights   = .prior_density_ordinate_compact(group$weights),
        source_transform = .prior_density_ordinate_compact(transforms)
      )
      next
    }

    if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior) ||
       is.prior.ordered(prior)){
      return(NULL)
    }
    full_support <- is.list(prior$truncation) &&
      identical(prior$truncation$lower, -Inf) &&
      identical(prior$truncation$upper, Inf)
    lognormal_support <- is.list(prior$truncation) &&
      identical(prior$truncation$lower, 0) &&
      identical(prior$truncation$upper, Inf)

    if(is.prior.vector(prior) && identical(prior$distribution, "mnormal") &&
       all(transforms == "none") && full_support){
      location <- prior$parameters$mean
      scale <- prior$parameters$sd
      family <- "mnormal"
    }else if(is.prior.simple(prior) && identical(prior$distribution, "normal") &&
             all(transforms == "none") && full_support){
      location <- prior$parameters$mean
      scale <- prior$parameters$sd
      family <- "normal"
    }else if(is.prior.simple(prior) && identical(prior$distribution, "lognormal") &&
             all(transforms == "log") && lognormal_support){
      location <- prior$parameters$meanlog
      scale <- prior$parameters$sdlog
      family <- "lognormal"
    }else{
      return(NULL)
    }

    mean_terms <- c(mean_terms, sum(group$weights) * location)
    scale_terms <- c(scale_terms, abs(group$weights) * scale)
    terms[[length(terms) + 1L]] <- list(
      parameter = group$parameter,
      family    = family,
      weights   = .prior_density_ordinate_compact(group$weights),
      source_transform = .prior_density_ordinate_compact(transforms)
    )
  }

  if(length(scale_terms) == 0L){
    return(NULL)
  }
  location <- sum(mean_terms)
  scale <- .prior_density_ordinate_stable_norm(scale_terms)
  provenance <- list(
    kind    = "linear_normal",
    mean    = location,
    sd      = scale,
    weights = .prior_density_ordinate_compact(weights),
    terms   = terms,
    support = c(lower = -Inf, upper = Inf)
  )
  if(!is.finite(location) || !is.finite(scale) || scale <= 0){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The analytic normal parameters are not representable.",
      provenance  = provenance
    ))
  }
  log_density <- stats::dnorm(
    value,
    mean = location,
    sd = scale,
    log = TRUE
  )

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "regular",
    log_density = log_density,
    exact       = TRUE,
    method      = "linear_normal",
    reason      = if(identical(log_density, -Inf)){
      paste0(
        "The continuous prior density is structurally regular, but its ",
        "log-density underflows in ordinary floating-point arithmetic."
      )
    }else{
      NULL
    },
    provenance  = provenance
  )
}

.prior_density_ordinate_deterministic_offset <- function(prior_list, weights,
                                                         source_transforms){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(0)
  }
  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  locations <- vapply(groups, function(group){
    location <- .prior_density_ordinate_point_group_location(
      group,
      source_transforms
    )
    if(length(location) != 1L || is.na(location)) NA_real_ else location
  }, numeric(1))
  if(anyNA(locations)){
    return(NULL)
  }
  sum(locations)
}

.prior_density_ordinate_additive_factor <- function(prior_list, weights,
                                                    source_transforms, value){

  weights <- weights[weights != 0]
  scalar <- .prior_density_ordinate_linear_scalar(
    prior_list,
    weights,
    source_transforms,
    value
  )
  if(!is.null(scalar)){
    return(scalar)
  }
  .prior_density_ordinate_linear_normal(
    prior_list,
    weights,
    source_transforms,
    value
  )
}

.prior_density_ordinate_product_singularity <- function(prior_list, split,
                                                        source_transforms,
                                                        value){

  if(length(split$product_groups) != 1L){
    return(NULL)
  }
  offset <- .prior_density_ordinate_deterministic_offset(
    prior_list,
    split$additive_weights,
    source_transforms
  )
  if(is.null(offset) || !is.finite(offset) || value != offset){
    return(NULL)
  }

  product_group <- split$product_groups[[1L]]
  factor <- .prior_density_ordinate_additive_factor(
    product_group$prior_list,
    product_group$weights,
    source_transforms,
    0
  )
  multiplier <- product_group$multiplier
  multiplier_prior <- prior_list[[multiplier]]
  if(is.null(factor) || !is.prior(multiplier_prior)){
    return(NULL)
  }
  multiplier_result <- .prior_density_ordinate_prior_affine(
    multiplier_prior,
    0,
    0,
    1,
    .prior_linear_source_transform(source_transforms[multiplier])
  )
  factor_behavior <- .prior_density_ordinate_continuous_behavior(factor)
  multiplier_behavior <-
    .prior_density_ordinate_continuous_behavior(multiplier_result)
  if(!identical(factor_behavior, "regular") ||
     !identical(multiplier_behavior, "regular")){
    return(NULL)
  }

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "infinite",
    log_density = Inf,
    exact       = TRUE,
    method      = "unsupported_provenance",
    reason      = paste0(
      "A product of two independent continuous factors with finite positive ",
      "densities at zero has a structural density singularity at the ",
      "requested value."
    ),
    provenance  = list(
      kind             = "product_singularity",
      singular_point   = offset,
      factor           = factor$provenance,
      multiplier       = multiplier_result$provenance,
      multiplier_name  = multiplier
    )
  )
}

.prior_density_ordinate_linear_base <- function(prior_list, weights,
                                                source_transforms, value){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = 0,
      probability = 1,
      method       = "scalar_affine",
      provenance   = list(kind = "scalar_affine", offset = 0, scale = 0)
    ))
  }
  if(is.null(source_transforms)){
    source_transforms <- rep(NA_character_, length(weights))
    names(source_transforms) <- names(weights)
  }else{
    source_transforms <- source_transforms[names(weights)]
  }
  unsupported_transform <- !is.na(source_transforms) & source_transforms != "log"
  if(any(unsupported_transform)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The source transformation is not supported structurally.",
      provenance  = list(
        kind              = "unsupported_provenance",
        weights           = .prior_density_ordinate_compact(weights),
        source_transforms = .prior_density_ordinate_compact(source_transforms)
      )
    ))
  }

  split <- tryCatch(
    .prior_linear_split_multiply_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(split)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      provenance  = list(kind = "unsupported_provenance")
    ))
  }
  if(length(split$product_groups) > 0L){
    singularity <- .prior_density_ordinate_product_singularity(
      prior_list,
      split,
      source_transforms,
      value
    )
    if(!is.null(singularity)){
      return(singularity)
    }
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "General products are not structurally classified.",
      provenance  = list(
        kind       = "general_product",
        weights    = .prior_density_ordinate_compact(weights),
        multipliers = names(split$product_groups)
      )
    ))
  }
  weights <- split$additive_weights
  weights <- weights[weights != 0]

  scalar <- .prior_density_ordinate_linear_scalar(
    prior_list,
    weights,
    source_transforms,
    value
  )
  if(!is.null(scalar)){
    scalar$provenance$weights <- .prior_density_ordinate_compact(weights)
    return(scalar)
  }

  normal <- .prior_density_ordinate_linear_normal(
    prior_list,
    weights,
    source_transforms,
    value
  )
  if(!is.null(normal)){
    return(normal)
  }

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "unknown",
    log_density = NA_real_,
    exact       = FALSE,
    method      = "unsupported_provenance",
    reason      = "General numerical convolutions are not structurally classified.",
    provenance  = list(
      kind              = "general_convolution",
      weights           = .prior_density_ordinate_compact(weights),
      source_families   = vapply(prior_list, function(prior){
        if(is.prior(prior) && !is.null(prior$distribution)) prior$distribution else "unknown"
      }, character(1)),
      source_transforms = .prior_density_ordinate_compact(source_transforms)
    )
  )
}

.prior_density_ordinate_transform_arguments <- function(transformation,
                                                        arguments){

  if(is.null(arguments)){
    arguments <- list()
  }
  if(!is.list(arguments)){
    return(NULL)
  }
  allowed <- if(transformation %in% c("lin", "exp_lin")) c("a", "b") else character()
  argument_names <- names(arguments)
  if((length(arguments) > 0L &&
      (is.null(argument_names) || anyNA(argument_names) ||
       any(!nzchar(argument_names)) || anyDuplicated(argument_names))) ||
     any(!argument_names %in% allowed)){
    return(NULL)
  }
  if(any(vapply(arguments, function(argument){
    !is.numeric(argument) || length(argument) != 1L ||
      is.na(argument) || !is.finite(argument)
  }, logical(1)))){
    return(NULL)
  }
  if(transformation %in% c("lin", "exp_lin")){
    if(is.null(arguments$a)) arguments$a <- 0
    if(is.null(arguments$b)) arguments$b <- 1
    arguments <- arguments[c("a", "b")]
  }
  arguments
}

.prior_density_ordinate_provenance_support <- function(provenance){

  if(!is.list(provenance)){
    return(NULL)
  }
  kind <- provenance$kind
  if(identical(kind, "primitive")){
    truncation <- provenance$truncation
    if(is.list(truncation)){
      lower <- truncation$lower
      upper <- truncation$upper
    }else{
      lower <- truncation[["lower"]]
      upper <- truncation[["upper"]]
    }
    if(is.numeric(lower) && length(lower) == 1L && !is.na(lower) &&
       is.numeric(upper) && length(upper) == 1L && !is.na(upper)){
      return(c(lower = lower, upper = upper))
    }
  }
  if(identical(kind, "linear_normal")){
    return(c(lower = -Inf, upper = Inf))
  }
  if(identical(kind, "scalar_affine") && is.list(provenance$source)){
    support <- .prior_density_ordinate_provenance_support(provenance$source)
    if(is.null(support)){
      return(NULL)
    }
    if(identical(provenance$source_transform, "log")){
      if(support[1L] < 0){
        return(NULL)
      }
      support <- c(
        if(support[1L] == 0) -Inf else log(support[1L]),
        if(is.infinite(support[2L])) Inf else log(support[2L])
      )
    }
    mapped <- provenance$offset + provenance$scale * support
    return(range(mapped))
  }
  if(identical(kind, "finite_mixture") && length(provenance$components) > 0L){
    components <- provenance$components[vapply(
      provenance$components,
      function(component){
        is.null(component$weight) ||
          (is.numeric(component$weight) && length(component$weight) == 1L &&
             is.finite(component$weight) && component$weight > 0)
      },
      logical(1)
    )]
    if(length(components) == 0L){
      return(NULL)
    }
    supports <- lapply(components, function(component){
      .prior_density_ordinate_provenance_support(component$provenance)
    })
    if(any(vapply(supports, is.null, logical(1)))){
      return(NULL)
    }
    return(c(
      lower = min(vapply(supports, `[[`, numeric(1), 1L)),
      upper = max(vapply(supports, `[[`, numeric(1), 2L))
    ))
  }
  NULL
}

.prior_density_ordinate_provenance_all_family <- function(provenance, family){

  if(!is.list(provenance)){
    return(FALSE)
  }
  if(identical(provenance$kind, "primitive")){
    return(identical(provenance$family, family))
  }
  if(identical(provenance$kind, "finite_mixture")){
    components <- provenance$components[vapply(
      provenance$components,
      function(component){
        is.null(component$weight) ||
          (is.numeric(component$weight) && length(component$weight) == 1L &&
             is.finite(component$weight) && component$weight > 0)
      },
      logical(1)
    )]
    return(length(components) > 0L && all(vapply(
      components,
      function(component){
        .prior_density_ordinate_provenance_all_family(
          component$provenance,
          family
        )
      },
      logical(1)
    )))
  }
  FALSE
}

.prior_density_ordinate_provenance_atoms <- function(provenance){

  if(!is.list(provenance)){
    return(NULL)
  }
  kind <- provenance$kind
  if(identical(kind, "linear_normal")){
    return(numeric())
  }
  if(identical(kind, "primitive")){
    family <- provenance$family
    if(identical(family, "none")){
      return(0)
    }
    if(identical(family, "point")){
      location <- provenance$parameters$location
      if(is.numeric(location) && length(location) == 1L &&
         is.finite(location)){
        return(unname(location))
      }
      return(NULL)
    }
    if(identical(family, "bernoulli")){
      probability <- provenance$parameters$probability
      if(!is.numeric(probability) || length(probability) != 1L ||
         !is.finite(probability) || probability < 0 || probability > 1){
        return(NULL)
      }
      locations <- c(0, 1)[c(1 - probability, probability) > 0]
      support <- .prior_density_ordinate_provenance_support(provenance)
      if(is.null(support)){
        return(NULL)
      }
      return(locations[locations >= support[1L] & locations <= support[2L]])
    }
    return(numeric())
  }
  if(identical(kind, "finite_mixture")){
    components <- provenance$components[vapply(
      provenance$components,
      function(component){
        is.null(component$weight) ||
          (is.numeric(component$weight) && length(component$weight) == 1L &&
             is.finite(component$weight) && component$weight > 0)
      },
      logical(1)
    )]
    atoms <- lapply(components, function(component){
      .prior_density_ordinate_provenance_atoms(component$provenance)
    })
    if(any(vapply(atoms, is.null, logical(1)))){
      return(NULL)
    }
    return(unlist(atoms, use.names = FALSE))
  }
  if(identical(kind, "scalar_affine") && is.list(provenance$source)){
    atoms <- .prior_density_ordinate_provenance_atoms(provenance$source)
    if(is.null(atoms)){
      return(NULL)
    }
    if(identical(provenance$source_transform, "log")){
      if(any(atoms <= 0)){
        return(NA_real_)
      }
      atoms <- log(atoms)
    }else if(!is.null(provenance$source_transform)){
      return(NULL)
    }
    return(provenance$offset + provenance$scale * atoms)
  }
  NULL
}

.prior_density_ordinate_tail_behavior <- function(provenance, direction,
                                                  transformation){

  support <- .prior_density_ordinate_provenance_support(provenance)
  if(!is.null(support)){
    bound <- if(direction < 0) support[1L] else support[2L]
    if(is.finite(bound)){
      return("zero")
    }
  }

  kind <- provenance$kind
  if(identical(kind, "linear_normal")){
    return("zero")
  }
  if(identical(kind, "scalar_affine") && is.list(provenance$source)){
    if(identical(provenance$source_transform, "log")){
      if(identical(transformation, "exp") && direction < 0){
        return(.prior_density_ordinate_exp_lin_boundary(
          provenance$source,
          provenance$scale
        ))
      }
      if(identical(transformation, "tanh") &&
         .prior_density_ordinate_provenance_all_family(
           provenance$source,
           "lognormal"
         )){
        return("zero")
      }
      return("unknown")
    }
    source_direction <- direction * sign(provenance$scale)
    return(.prior_density_ordinate_tail_behavior(
      provenance$source,
      source_direction,
      transformation
    ))
  }
  if(identical(kind, "finite_mixture")){
    components <- provenance$components[vapply(
      provenance$components,
      function(component){
        is.null(component$weight) ||
          (is.numeric(component$weight) && length(component$weight) == 1L &&
             is.finite(component$weight) && component$weight > 0)
      },
      logical(1)
    )]
    if(length(components) == 0L) return("unknown")
    behavior <- vapply(components, function(component){
      .prior_density_ordinate_tail_behavior(
        component$provenance,
        direction,
        transformation
      )
    }, character(1))
    if(any(behavior == "infinite")) return("infinite")
    if(any(behavior == "unknown")) return("unknown")
    if(any(behavior == "regular")) return("regular")
    return("zero")
  }
  if(identical(kind, "primitive")){
    family <- provenance$family
    if(family %in% c("normal", "moment")) return("zero")
    if(family %in% c("t", "invmoment")) return("infinite")
    if(identical(transformation, "tanh") &&
       identical(family, "lognormal") && direction > 0) return("infinite")
  }
  "unknown"
}

.prior_density_ordinate_exp_lin_boundary <- function(provenance, b){

  support <- .prior_density_ordinate_provenance_support(provenance)
  if(!is.null(support)){
    relevant <- if(b > 0) support[1L] else support[2L]
    if(is.finite(relevant) && relevant > 0){
      return("zero")
    }
  }
  if(identical(provenance$kind, "finite_mixture")){
    components <- provenance$components[vapply(
      provenance$components,
      function(component){
        is.null(component$weight) ||
          (is.numeric(component$weight) && length(component$weight) == 1L &&
             is.finite(component$weight) && component$weight > 0)
      },
      logical(1)
    )]
    if(length(components) == 0L) return("unknown")
    behavior <- vapply(components, function(component){
      .prior_density_ordinate_exp_lin_boundary(component$provenance, b)
    }, character(1))
    if(any(behavior == "infinite")) return("infinite")
    if(any(behavior == "unknown")) return("unknown")
    if(any(behavior == "regular")) return("regular")
    return("zero")
  }
  if(!identical(provenance$kind, "primitive")){
    return("unknown")
  }

  family <- provenance$family
  parameters <- provenance$parameters
  if(identical(family, "lognormal")){
    return("zero")
  }
  if(b < 0 && family %in% c("gamma", "exp")){
    return("zero")
  }
  if(b > 0 && identical(family, "invgamma")){
    return("zero")
  }
  exponent <- if(b > 0){
    if(identical(family, "gamma")) parameters$shape - b else
      if(identical(family, "exp")) 1 - b else
        if(identical(family, "beta")) parameters$alpha - b else
          if(identical(family, "uniform")) 1 - b else NA_real_
  }else{
    if(identical(family, "invgamma")) -(parameters$shape + b) else NA_real_
  }
  if(is.na(exponent)) return("unknown")
  if(exponent > 0) "zero" else if(exponent < 0) "infinite" else "regular"
}

.prior_density_ordinate_named_transform <- function(classifier, source_provenance,
                                                    transformation, arguments,
                                                    value){

  provenance <- list(
    kind           = "named_transform",
    transformation = if(is.character(transformation)) transformation else "custom",
    arguments      = .prior_density_ordinate_compact(arguments),
    source         = source_provenance
  )
  if(!is.character(transformation) || length(transformation) != 1L ||
     !transformation %in% c("lin", "exp", "exp_lin", "tanh")){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "Arbitrary user transformations are not structurally classified.",
      provenance  = provenance
    ))
  }
  arguments <- .prior_density_ordinate_transform_arguments(
    transformation,
    arguments
  )
  if(is.null(arguments)){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "undefined",
      log_density = NA_real_,
      exact       = TRUE,
      method      = "named_transform",
      reason      = "The named transformation arguments are invalid.",
      provenance  = provenance
    ))
  }
  provenance$arguments <- arguments

  if(identical(transformation, "lin")){
    if(arguments$b == 0){
      return(.prior_density_ordinate_atom_result(
        value       = value,
        locations   = arguments$a,
        probability = 1,
        method       = "named_transform",
        provenance   = provenance
      ))
    }
    source_value <- (value - arguments$a) / arguments$b
    if(!is.finite(source_value)){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "unknown",
        log_density = NA_real_,
        exact       = FALSE,
        method      = "unsupported_provenance",
        provenance  = provenance
      ))
    }
    source <- classifier(source_value)
    return(.prior_density_ordinate_wrap(
      source,
      value,
      log(abs(arguments$b)),
      "named_transform",
      provenance
    ))
  }

  if(identical(transformation, "exp")){
    if(value < 0){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "zero",
        log_density = -Inf,
        exact       = TRUE,
        method      = "named_transform",
        reason      = "The requested value is outside the transformed prior support.",
        provenance  = provenance
      ))
    }
    if(value == 0){
      behavior <- .prior_density_ordinate_tail_behavior(
        source_provenance,
        -1,
        "exp"
      )
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = behavior,
        log_density = if(behavior == "zero") -Inf else
          if(behavior == "infinite") Inf else NA_real_,
        exact       = !identical(behavior, "unknown"),
        method      = if(behavior == "unknown") "unsupported_provenance" else "named_transform",
        provenance  = provenance
      ))
    }
    source_value <- log(value)
    source <- classifier(source_value)
    return(.prior_density_ordinate_wrap(
      source,
      value,
      log(value),
      "named_transform",
      provenance
    ))
  }

  if(identical(transformation, "tanh")){
    if(abs(value) > 1){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "zero",
        log_density = -Inf,
        exact       = TRUE,
        method      = "named_transform",
        reason      = "The requested value is outside the transformed prior support.",
        provenance  = provenance
      ))
    }
    if(abs(value) == 1){
      behavior <- .prior_density_ordinate_tail_behavior(
        source_provenance,
        sign(value),
        "tanh"
      )
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = behavior,
        log_density = if(behavior == "zero") -Inf else
          if(behavior == "infinite") Inf else NA_real_,
        exact       = !identical(behavior, "unknown"),
        method      = if(behavior == "unknown") "unsupported_provenance" else "named_transform",
        provenance  = provenance
      ))
    }
    source_value <- atanh(value)
    source <- classifier(source_value)
    return(.prior_density_ordinate_wrap(
      source,
      value,
      log1p(-value^2),
      "named_transform",
      provenance
    ))
  }

  if(arguments$b == 0){
    location <- exp(arguments$a)
    if(!is.finite(location) || location == 0){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "unknown",
        log_density = NA_real_,
        exact       = FALSE,
        method      = "unsupported_provenance",
        reason      = paste0(
          "The constant transformed point is not representable in ordinary ",
          "floating-point arithmetic."
        ),
        provenance  = provenance
      ))
    }
    return(.prior_density_ordinate_atom_result(
      value       = value,
      locations   = location,
      probability = 1,
      method       = "named_transform",
      provenance   = provenance
    ))
  }
  support <- .prior_density_ordinate_provenance_support(source_provenance)
  if(is.null(support) || support[1L] < 0){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "undefined",
      log_density = NA_real_,
      exact       = TRUE,
      method      = "named_transform",
      reason      = paste0(
        "The exponential-linear transformation is not defined on the ",
        "established source support."
      ),
      provenance  = provenance
    ))
  }
  atoms <- .prior_density_ordinate_provenance_atoms(source_provenance)
  if(!is.null(atoms) && (anyNA(atoms) || any(atoms <= 0))){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "undefined",
      log_density = NA_real_,
      exact       = TRUE,
      method      = "named_transform",
      reason      = paste0(
        "The exponential-linear transformation is undefined for a positive-",
        "probability source atom."
      ),
      provenance  = provenance
    ))
  }
  if(value < 0){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "zero",
      log_density = -Inf,
      exact       = TRUE,
      method      = "named_transform",
      reason      = "The requested value is outside the transformed prior support.",
      provenance  = provenance
    ))
  }
  if(value == 0){
    behavior <- .prior_density_ordinate_exp_lin_boundary(
      source_provenance,
      arguments$b
    )
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = behavior,
      log_density = if(behavior == "zero") -Inf else
        if(behavior == "infinite") Inf else NA_real_,
      exact       = !identical(behavior, "unknown"),
      method      = if(behavior == "unknown") "unsupported_provenance" else "named_transform",
      provenance  = provenance
    ))
  }
  source_value <- exp((log(value) - arguments$a) / arguments$b)
  if(!is.finite(source_value) || source_value <= 0){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = "The inverse transformed value is not representable.",
      provenance  = provenance
    ))
  }
  source <- classifier(source_value)
  log_jacobian <- log(abs(arguments$b)) + arguments$a +
    (arguments$b - 1) * log(source_value)
  .prior_density_ordinate_wrap(
    source,
    value,
    log_jacobian,
    "named_transform",
    provenance
  )
}

.prior_density_ordinate_linear_arguments <- function(arguments, value){

  prior_list <- arguments$prior_list
  weights <- arguments$weights
  source_transforms <- arguments$source_transforms
  transformation <- arguments$output_transformation
  transformation_arguments <- arguments$output_transformation_arguments

  if(!is.list(prior_list) || !is.numeric(weights) || is.null(names(weights)) ||
     anyNA(weights) || any(!is.finite(weights))){
    return(.prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      provenance  = list(kind = "unsupported_provenance")
    ))
  }

  classifier <- function(source_value){
    .prior_density_ordinate_linear_base(
      prior_list,
      weights,
      source_transforms,
      source_value
    )
  }
  if(is.null(transformation)){
    return(classifier(value))
  }
  template <- classifier(0)
  .prior_density_ordinate_named_transform(
    classifier,
    template$provenance,
    transformation,
    transformation_arguments,
    value
  )
}

.prior_density_ordinate_context_classifier <- function(context, weights,
                                                       source_transforms,
                                                       transformation,
                                                       transformation_arguments,
                                                       value){

  if(inherits(context, "prior_density_context")){
    standardized <- tryCatch(
      .prior_density_context_standardized_weights(context, weights),
      error = function(e) NULL
    )
    if(is.null(standardized)){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "unknown",
        log_density = NA_real_,
        exact       = FALSE,
        method      = "unsupported_provenance",
        provenance  = list(kind = "density_context")
      ))
    }
    if(!is.null(source_transforms)){
      source_transforms <- source_transforms[names(standardized)]
    }
    result <- .prior_density_ordinate_linear_arguments(list(
      prior_list                      = context$prior_list,
      weights                         = standardized,
      source_transforms               = source_transforms,
      output_transformation           = transformation,
      output_transformation_arguments = transformation_arguments
    ), value)
    result$provenance$context <- list(
      kind                 = "prior_density_context",
      requested_weights    = .prior_density_ordinate_compact(weights),
      standardized_weights = .prior_density_ordinate_compact(standardized)
    )
    return(result)
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    component_classifier <- function(source_value){
      results <- lapply(seq_along(context$model_weights), function(model_i){
        model_prior_list <- lapply(context$prior_list, function(parameter_priors){
          if(is.prior(parameter_priors)) parameter_priors else parameter_priors[[model_i]]
        })
        for(parameter in names(model_prior_list)){
          if(is.null(model_prior_list[[parameter]])){
            model_prior_list[[parameter]] <- prior("point", list(location = 0))
          }
        }
        .prior_density_ordinate_linear_base(
          model_prior_list,
          weights,
          source_transforms,
          source_value
        )
      })
      .prior_density_ordinate_combine(
        results,
        context$model_weights,
        source_value,
        provenance_extra = list(context = "model_mixture")
      )
    }
    if(is.null(transformation)){
      return(component_classifier(value))
    }
    template <- component_classifier(0)
    return(.prior_density_ordinate_named_transform(
      component_classifier,
      template$provenance,
      transformation,
      transformation_arguments,
      value
    ))
  }

  if(inherits(context, "prior_density_conditional_context")){
    component_classifier <- function(source_value){
      results <- lapply(context$prior_lists, function(prior_list){
        if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
          component_context <- .prior_density_context(
            prior_list,
            context$column_names,
            context$formula_scale,
            context$n_grid,
            context$tail_prob
          )
          return(.prior_density_ordinate_context_classifier(
            component_context,
            weights,
            source_transforms,
            NULL,
            NULL,
            source_value
          ))
        }
        .prior_density_ordinate_linear_base(
          prior_list,
          weights,
          source_transforms,
          source_value
        )
      })
      .prior_density_ordinate_combine(
        results,
        context$model_weights,
        source_value,
        provenance_extra = list(context = "conditional_mixture")
      )
    }
    if(is.null(transformation)){
      return(component_classifier(value))
    }
    template <- component_classifier(0)
    return(.prior_density_ordinate_named_transform(
      component_classifier,
      template$provenance,
      transformation,
      transformation_arguments,
      value
    ))
  }

  .prior_density_ordinate_result(
    value       = value,
    behavior    = "unknown",
    log_density = NA_real_,
    exact       = FALSE,
    method      = "unsupported_provenance",
    provenance  = list(kind = "unknown_density_context")
  )
}

.prior_density_ordinate_from_adaptive <- function(adaptive, value){

  if(!is.list(adaptive) || !is.character(adaptive$kind) ||
     length(adaptive$kind) != 1L || !is.list(adaptive$arguments)){
    return(NULL)
  }
  arguments <- adaptive$arguments
  if(identical(adaptive$kind, "linear_combination")){
    return(.prior_density_ordinate_linear_arguments(arguments, value))
  }
  if(identical(adaptive$kind, "density_context")){
    return(.prior_density_ordinate_context_classifier(
      arguments$context,
      arguments$weights,
      arguments$source_transforms,
      arguments$output_transformation,
      arguments$output_transformation_arguments,
      value
    ))
  }
  if(identical(adaptive$kind, "density_context_rows")){
    weights <- arguments$weights
    if(!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights))){
      return(.prior_density_ordinate_result(
        value       = value,
        behavior    = "unknown",
        log_density = NA_real_,
        exact       = FALSE,
        method      = "unsupported_provenance",
        reason      = "The row-varying density weights are malformed.",
        provenance  = list(kind = "density_context_rows")
      ))
    }
    if(is.null(dim(weights))){
      weights <- matrix(weights, nrow = 1L, dimnames = list(NULL, names(weights)))
    }
    weights <- as.matrix(weights)
    if(nrow(weights) == 0L){
      return(.prior_density_ordinate_atom_result(
        value,
        0,
        1,
        "finite_mixture",
        list(kind = "density_context_rows", rows = 0L)
      ))
    }
    row_keys <- apply(weights, 1L, function(row){
      paste(sprintf("%a", row), collapse = "\r")
    })
    unique_keys <- unique(row_keys)
    row_counts <- tabulate(match(row_keys, unique_keys), nbins = length(unique_keys))
    row_indices <- match(unique_keys, row_keys)
    results <- lapply(row_indices, function(row_i){
      .prior_density_ordinate_context_classifier(
        arguments$context,
        weights[row_i, ],
        arguments$source_transforms,
        arguments$output_transformation,
        arguments$output_transformation_arguments,
        value
      )
    })
    combined <- .prior_density_ordinate_combine(
      results,
      row_counts,
      value
    )
    continuous_behavior <-
      .prior_density_ordinate_continuous_behavior(combined)
    combined$provenance <- list(
      kind              = "finite_mixture",
      context           = "density_context_rows",
      weight_dimensions = as.integer(dim(weights)),
      unique_rows       = length(unique_keys),
      total_rows        = nrow(weights),
      weights_hash      = .prior_density_ordinate_numeric_hash(weights),
      row_classifications = Map(function(result, count){
        list(
          count               = unname(count),
          behavior            = result$behavior,
          continuous_behavior =
            .prior_density_ordinate_continuous_behavior(result),
          point_mass           = result$point_mass,
          method               = result$method,
          source_kind          = result$provenance$kind
        )
      }, results, row_counts)
    )
    if(identical(combined$behavior, "point_mass")){
      combined$provenance$continuous_behavior <- continuous_behavior
    }
    return(combined)
  }
  NULL
}

.prior_density_ordinate_numeric_hash <- function(x){

  encoded <- paste(
    c(dim(x), sprintf("%a", as.numeric(x))),
    collapse = "\r"
  )
  bytes <- as.integer(charToRaw(encoded))
  hash <- 0
  for(byte in bytes){
    hash <- (hash * 131 + byte) %% 2147483647
  }
  sprintf("%08x", as.integer(hash))
}

.prior_density_ordinate_grid_log_density <- function(x, value){

  if(is.null(x$density) || !is.list(x$density) ||
     is.null(x$density$mass) || x$density$mass <= 0){
    return(-Inf)
  }
  height <- tryCatch(
    .prior_linear_density_grid_height(x, value),
    error = function(e) NA_real_
  )
  if(!is.numeric(height) || length(height) != 1L || is.na(height) || height < 0){
    return(NA_real_)
  }
  if(height == 0) -Inf else log(height)
}

.prior_density_ordinate_linear <- function(x, value){

  adaptive <- attr(x, "adaptive_evaluation", exact = TRUE)
  result <- .prior_density_ordinate_from_adaptive(adaptive, value)
  if(is.null(result)){
    result <- .prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = .prior_density_ordinate_grid_log_density(x, value),
      exact       = FALSE,
      method      = "unsupported_provenance",
      provenance  = list(
        kind         = "unsupported_provenance",
        source_class = class(x)
      )
    )
  }

  stored_point_mass <- .prior_linear_density_point_mass(x, value)
  continuous_behavior <- .prior_density_ordinate_continuous_behavior(result)
  result$point_mass <- stored_point_mass
  if(stored_point_mass > 0){
    result$behavior <- "point_mass"
    result$exact <- TRUE
    result$reason <- .prior_density_ordinate_reason("point_mass")
    result$provenance$continuous_behavior <- continuous_behavior
    result$provenance$point_locations <- if(!is.null(x$points)){
      .prior_density_ordinate_compact(x$points$x[x$points$p > 0])
    }else{
      numeric()
    }
  }else if(identical(result$behavior, "point_mass")){
    result$behavior <- continuous_behavior
    result$reason <- .prior_density_ordinate_reason(continuous_behavior)
  }

  if(result$behavior %in% c("unknown", "point_mass") &&
     is.na(result$log_density)){
    result$log_density <- .prior_density_ordinate_grid_log_density(x, value)
  }
  result$value <- value
  class(result) <- c("prior_density_ordinate", "list")
  result
}
