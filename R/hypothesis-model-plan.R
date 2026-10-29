# ============================================================================ #
# hypothesis-model-plan.R
# ============================================================================ #
#
# One plan per hypothesis statement. A plan records what the statement tests
# and what refuses it:
#   - kind: "scalar" (a scalar parameter or coefficient), "linear" (factor
#     levels of every contrast, ordered levels, and linear combinations of
#     the levels of one term: level contrasts and derived linear
#     expressions), "random" (random-effect quantities), "exp_affine"
#     (exponentiated affine original-scale coefficients), or "marginal_mean";
#   - targets: one per point value, with the target's weights on the fitted
#     coordinates, its prior density built by BayesTools, its declared
#     atoms, and the classification of the value by prior_ordinate_status();
#   - the refusals of the plan and of its density methods (method_refusals,
#     keyed by the name of a density backend or by its type).
# hypothesis_BF_model() and hypothesis_quantities() execute or render plans
# only: the eligibility of a statement has one source. Region statements are
# method-free: their plans refuse only statements on fixed quantities or
# unsupported targets. A plan is structural; whether a density backend can
# evaluate it is decided on every call (.hypothesis_model_plan_status()).
#
# ============================================================================ #


# A refusal: the reason and the condition classes of the error a refused
# statement stops with. BayesTools refusals keep their classes ('condition',
# e.g. "BayesTools_point_mass_at_null"). The refusals of the engine of a test
# that cannot be computed have the class "BayesTools_hypothesis_<type>" with
# the parent "BayesTools_hypothesis_unavailable": "fixed" (the quantity is
# fixed by the fitted model), "target" (no supported target or prior
# density), and "method" (the density method is unavailable for the target). A
# statement that has to be restated, or an argument that has to change, is no
# unavailable test: "statement" (the statement's form is unsupported) has the
# class "BayesTools_hypothesis_statement" only, and "ambiguous" (the
# statement names several parameters; 'namespace' selects one, and for
# marginal means an alias of several parameters is resolved by 'parameter')
# has it as its parent (.hypothesis_model_ambiguous_class()).
.hypothesis_model_refusal <- function(reason, type = NULL, condition = NULL) {

  class <- if (!is.null(condition)) {
    unique(condition)
  } else if (identical(type, "ambiguous")) {
    .hypothesis_model_ambiguous_class()
  } else if (identical(type, "statement")) {
    "BayesTools_hypothesis_statement"
  } else {
    c(paste0("BayesTools_hypothesis_", type), "BayesTools_hypothesis_unavailable")
  }

  list(reason = reason, class = class)
}


# A refusal returned by a provider or a density backend: a list with the
# 'reason', optionally the 'type' of the refusal (see
# .hypothesis_model_refusal()), additional condition classes 'class', and
# further 'fields' of the condition. Without a 'type' (and without a default
# one), 'class' holds the complete classes of the condition.
.hypothesis_model_as_refusal <- function(x, default_type = NULL) {

  if (is.null(x)) {
    return(NULL)
  }
  if (!is.list(x) || !is.character(x[["reason"]]) ||
      length(x[["reason"]]) != 1L || is.na(x[["reason"]])) {
    stop("A refusal must be a list with a 'reason'.", call. = FALSE)
  }
  type <- if (!is.null(x[["type"]])) x[["type"]] else default_type
  refusal <- if (!is.null(type)) {
    out <- .hypothesis_model_refusal(x[["reason"]], type)
    out[["class"]] <- c(out[["class"]], x[["class"]])
    out
  } else {
    list(reason = x[["reason"]], class = unique(x[["class"]]))
  }
  if (length(x[["fields"]]) > 0L) {
    refusal[["fields"]] <- x[["fields"]]
  }

  refusal
}


# The classes of an ambiguous statement, at every entry point: its own class
# and the parent "BayesTools_hypothesis_statement".
.hypothesis_model_ambiguous_class <- function() {

  c("BayesTools_hypothesis_ambiguous", "BayesTools_hypothesis_statement")
}


# A refusal of the ordinate backends for a generic 'cause' (a nonlinear
# original-scale coefficient, "original_scale"; a conditional random-effect
# statement, "conditional_random"; a random-effect quantity without a scalar
# coordinate, "random_target"): the refusal of 'type' with the reason, and the
# 'cause', for which a density backend may declare additional condition
# classes (.hypothesis_model_plan_status()), so that one cause has the same
# classes at every entry point of the backend.
.hypothesis_model_cause_refusal <- function(reason, cause, type = "method") {

  refusal <- .hypothesis_model_refusal(reason, type)
  refusal[["cause"]] <- cause

  refusal
}


# The error with which a helper of the planning reports an ordinate cause; the
# planning turns it into a refusal of the plan (the class
# "BayesTools_hypothesis_cause").
.hypothesis_model_stop_cause <- function(reason, cause) {

  refusal <- .hypothesis_model_cause_refusal(reason, cause)
  stop(structure(
    class = c(refusal[["class"]], "BayesTools_hypothesis_cause", "error",
              "condition"),
    list(message = reason, call = NULL, cause = cause)
  ))
}


# An unknown reference is a statement error with the classes with which
# BayesTools refuses it ("BayesTools_parameter_not_found" and its parent
# "BayesTools_parameter_resolution_error"), after "BayesTools_hypothesis_statement":
# the classes of an unknown name or level at every entry point.
.hypothesis_model_not_found_class <- function() {

  c("BayesTools_hypothesis_statement", "BayesTools_parameter_not_found",
    "BayesTools_parameter_resolution_error")
}


# The refusal of a statement that references an unknown level of a factor
# term (the resolver refuses unknown levels of fitted terms first; marginal
# means are refused by their plans).
.hypothesis_model_unknown_level_refusal <- function(level, parameter) {

  .hypothesis_model_refusal(
    paste0(
      "Hypothesis references unknown level '", level, "' for parameter '",
      parameter, "'."
    ),
    condition = .hypothesis_model_not_found_class()
  )
}


# The refusal of a condition: its message (or 'reason'), its classes, and
# its other fields (e.g. BayesTools' 'alias' and 'available'), with which
# .hypothesis_model_stop() raises it.
.hypothesis_model_refusal_from_condition <- function(condition, reason = NULL) {

  fields <- unclass(condition)
  list(
    reason = if (is.null(reason)) conditionMessage(condition) else reason,
    class  = setdiff(class(condition), c("error", "condition")),
    fields = fields[setdiff(names(fields), c("message", "call"))]
  )
}


# Raises a refusal (a list with 'reason', 'class', and optionally 'fields', or
# the refusal of a provider or backend with a 'type') as a classed error.
.hypothesis_model_stop <- function(refusal) {

  if (!is.null(refusal[["type"]])) {
    refusal <- .hypothesis_model_as_refusal(refusal)
  }
  stop(structure(
    class = c(refusal[["class"]], "error", "condition"),
    c(list(message = refusal[["reason"]], call = NULL), refusal[["fields"]])
  ))
}


#' @title Stop with a refusal of a hypothesis
#'
#' @description Raises the error with which the hypothesis engine
#' ([hypothesis_BF_model()]) refuses a statement or a test that it cannot
#' compute, so that other entry points of a package refuse the same cause with
#' the same condition classes.
#'
#' @param reason the message of the condition.
#' @param type the kind of refusal: \code{"target"} (no supported target or
#'   prior density), \code{"fixed"} (the quantity is fixed by the fitted
#'   model), or \code{"method"} (the density method is unavailable for the
#'   target) have the classes \code{BayesTools_hypothesis_<type>} and
#'   \code{BayesTools_hypothesis_unavailable}; \code{"statement"} (the
#'   statement has to be restated) has the class
#'   \code{BayesTools_hypothesis_statement}, and \code{"ambiguous"} (the
#'   statement names several parameters) the classes
#'   \code{BayesTools_hypothesis_ambiguous} and
#'   \code{BayesTools_hypothesis_statement}.
#' @param class additional condition classes appended to those of the type.
#' @param fields a named list of further fields of the condition.
#'
#' @return Does not return; it signals the error.
#'
#' @seealso [hypothesis_BF_model()]
#' @export
hypothesis_stop <- function(reason, type = c("target", "fixed", "method",
                                             "statement", "ambiguous"),
                            class = NULL, fields = NULL){

  check_char(reason, "reason", check_length = 1, allow_NA = FALSE)
  type <- match.arg(type)
  check_char(class, "class", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  if(!is.null(fields) && (!is.list(fields) || is.null(names(fields)) ||
                          any(!nzchar(names(fields))))){
    stop("'fields' must be a named list.", call. = FALSE)
  }
  refusal <- .hypothesis_model_refusal(reason, type)
  refusal[["class"]]  <- c(refusal[["class"]], class)
  refusal[["fields"]] <- fields

  .hypothesis_model_stop(refusal)
}


# The refusal of the density method 'backend' for a plan, or NULL when the
# method evaluates the statement. It is computed on every call from the
# registered backend and the normalized 'control', and never cached. The
# precedence is: a refusal of the whole plan (fixed quantities, unsupported
# targets and statements), then the refusals of the point targets (their
# values), then the refusal of a provider without the model the method needs
# (a marginal-means object without its source fit), then the capability of
# the backend for the model, then the quantity's own refusal of the method.
# The capability precedes the quantity's refusal, so that a request with both
# causes reports the capability at every entry point. Region statements need
# no density method.
.hypothesis_model_plan_status <- function(plan, provider, backend,
                                          control = NULL) {

  if (!is.null(plan[["refusal"]])) {
    return(plan[["refusal"]])
  }
  if (!plan[["point"]]) {
    return(NULL)
  }
  for (target in plan[["targets"]]) {
    refusal <- .hypothesis_model_plan_target_refusal(target)
    if (!is.null(refusal)) {
      return(refusal)
    }
  }
  if (!is.null(plan[["source_refusals"]][[backend[["type"]]]])) {
    return(plan[["source_refusals"]][[backend[["type"]]]])
  }
  if (is.function(backend[["available"]])) {
    capability_object <- provider[["capability_object"]]
    if (is.null(capability_object)) {
      capability_object <- provider[["object"]]
    }
    refusal <- .hypothesis_model_as_refusal(
      backend[["available"]](capability_object, plan, control),
      default_type = "method"
    )
    if (!is.null(refusal)) {
      return(refusal)
    }
  }
  refusal <- plan[["method_refusals"]][[backend[["name"]]]]
  if (is.null(refusal)) {
    refusal <- plan[["method_refusals"]][[backend[["type"]]]]
  }
  cause <- refusal[["cause"]]
  if (is.null(cause)) {
    cause <- refusal[["fields"]][["cause"]]
  }
  if (!is.null(cause)) {
    refusal[["class"]] <- c(
      refusal[["class"]],
      .hypothesis_model_backend_cause_classes(backend, cause)
    )
  }

  refusal
}


# Stops with the refusal of the density method 'backend' for the plan.
.hypothesis_model_plan_check <- function(plan, provider, backend,
                                         control = NULL) {

  refusal <- .hypothesis_model_plan_status(plan, provider, backend, control)
  if (!is.null(refusal)) {
    .hypothesis_model_stop(refusal)
  }

  invisible(plan)
}


# ---------------------------------------------------------------------------- #
# Plan caches
# ---------------------------------------------------------------------------- #

# The cache of the plans of a model object (a fit or a marginal-means object):
# mixed posteriors, marginal posteriors, random-effect selections, and the
# plans of the statements themselves are built once per object and kept for
# the session, whichever of hypothesis_BF_model() and hypothesis_quantities()
# asks for them again. An object is identified by the hash of its whole
# content (draws, data, priors and everything else the plans read), so another
# object never reads them; an object that cannot be hashed gets a cache of its
# own call only. The caches are environments of the package, so they never
# travel with a saved fit, leave with the session, and are dropped when the
# package is reloaded. Each is bounded by its entries and their size, and the
# package keeps the caches of a few objects, the oldest out.
.hypothesis_model_plan_registry <- new.env(parent = emptyenv())


# Objects whose caches are kept.
.hypothesis_model_plan_registry_limit <- function() 4L


# Entries, and bytes, kept per object.
.hypothesis_model_plan_cache_limit <- function() 256L


.hypothesis_model_plan_cache_max_bytes <- function() 256 * 1024^2


.hypothesis_model_plan_cache <- function(provider = NULL) {

  cache <- if (!is.null(provider)) .hypothesis_model_plan_object_cache(provider)

  if (is.null(cache)) new.env(parent = emptyenv()) else cache
}


.hypothesis_model_plan_object_cache <- function(provider) {

  key <- tryCatch(rlang::hash(provider[["object"]]), error = function(e) NULL)
  if (is.null(key)) {
    return(NULL)
  }
  registry <- .hypothesis_model_plan_registry
  if (!exists(key, envir = registry, inherits = FALSE)) {
    cache <- new.env(parent = emptyenv())
    assign(".bounded", TRUE, envir = cache)
    order <- c(setdiff(registry[[".order"]], key), key)
    assign(key, cache, envir = registry)
    if (length(order) > .hypothesis_model_plan_registry_limit()) {
      evicted <- order[seq_len(length(order) - .hypothesis_model_plan_registry_limit())]
      rm(list = evicted, envir = registry)
      order <- setdiff(order, evicted)
    }
    assign(".order", order, envir = registry)
  }

  get(key, envir = registry, inherits = FALSE)
}


# Drop the caches of every object, and those of the density backends. A test
# that replaces a function the plans are built with, between two plans of one
# statement, drops them too.
.hypothesis_model_plan_cache_clear <- function() {

  rm(list = ls(.hypothesis_model_plan_registry, all.names = TRUE),
     envir = .hypothesis_model_plan_registry)
  for (backend in .hypothesis_model_backends()) {
    .hypothesis_model_backend_clear_cache(backend)
  }

  invisible(NULL)
}


.hypothesis_model_plan_cached <- function(cache, key, compute) {

  if (is.null(cache)) {
    return(compute())
  }
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  value <- compute()
  .hypothesis_model_plan_cache_set(cache, key, value)

  value
}


# Keep `value` under `key`. A bounded cache drops its oldest entries beyond its
# entry and size bounds, and keeps no entry above the size bound: that value is
# computed again the next time it is asked for.
.hypothesis_model_plan_cache_set <- function(cache, key, value) {

  if (!isTRUE(cache[[".bounded"]])) {
    assign(key, value, envir = cache)
    return(invisible(TRUE))
  }
  size <- as.numeric(utils::object.size(value))
  if (size > .hypothesis_model_plan_cache_max_bytes()) {
    return(invisible(FALSE))
  }
  order <- c(setdiff(cache[[".order"]], key), key)
  sizes <- cache[[".sizes"]]
  sizes <- c(sizes[setdiff(names(sizes), key)], stats::setNames(size, key))
  assign(key, value, envir = cache)
  while (length(order) > .hypothesis_model_plan_cache_limit() ||
         sum(sizes[order]) > .hypothesis_model_plan_cache_max_bytes()) {
    rm(list = order[[1L]], envir = cache)
    order <- order[-1L]
  }
  assign(".order", order, envir = cache)
  assign(".sizes", sizes[order], envir = cache)

  invisible(TRUE)
}

# ---------------------------------------------------------------------------- #
# Plans
# ---------------------------------------------------------------------------- #

# The plans of the statements of a hypothesis, one per statement.
.hypothesis_model_plans <- function(provider, hypothesis, component = "auto",
                                    level_component = NULL,
                                    transform_scaled = FALSE,
                                    conditional = FALSE, metadata = NULL,
                                    n_samples = 10000L,
                                    cache = .hypothesis_model_plan_cache(provider)) {

  if (is.null(metadata)) {
    metadata <- provider[["metadata"]]()
  }
  ast        <- .hypothesis_model_ast(hypothesis, metadata[["catalog"]])
  statements <- hypothesis_render(ast)
  # A plan is a function of its statement and arguments (and of the object and
  # its metadata, which the cache is of).
  arguments  <- rlang::hash(list(
    component, level_component, transform_scaled, conditional,
    as.numeric(n_samples), metadata
  ))

  lapply(statements, function(statement) {
    .hypothesis_model_plan_cached(
      cache,
      paste("plan", arguments, rlang::hash(statement), sep = "\r"),
      function() {
        .hypothesis_model_plan(
          provider         = provider,
          statement        = statement,
          component        = component,
          level_component  = level_component,
          transform_scaled = transform_scaled,
          conditional      = conditional,
          metadata         = metadata,
          n_samples        = n_samples,
          cache            = cache
        )
      }
    )
  })
}


# The plan of one statement of a fitted model on the original ('transform_scaled')
# or fitted coefficient scale, conditional on the inclusion of the tested
# parameter ('conditional', mixture priors).
.hypothesis_model_plan <- function(provider, statement, component = "auto",
                                   level_component = NULL,
                                   transform_scaled = FALSE,
                                   conditional = FALSE, metadata = NULL,
                                   n_samples = 10000L,
                                   cache = .hypothesis_model_plan_cache()) {

  if (is.null(metadata)) {
    metadata <- provider[["metadata"]]()
  }
  statement <- .hypothesis_model_ast(statement, metadata[["catalog"]])
  if (length(statement[["statements"]]) != 1L) {
    stop("Internal error: a hypothesis plan covers one statement.",
         call. = FALSE)
  }
  selected <- .hypothesis_model_select_parameter(
    provider     = provider,
    hypothesis = statement,
    component  = component,
    metadata   = metadata,
    level_component = level_component
  )
  .hypothesis_model_check_supported_component(selected[["component"]], provider)
  parameter <- selected[["parameter"]]

  plan <- .hypothesis_model_plan_new(
    statement    = statement,
    hypothesis   = .hypothesis_model_rewrite(
      hypothesis = statement,
      aliases    = selected[["aliases"]],
      parameter  = parameter,
      resolved_roots = selected[["roots"]]
    ),
    parameter    = parameter,
    label        = .hypothesis_model_alias_label(selected[["aliases"]], parameter),
    component    = selected[["component"]],
    selected     = selected,
    transform_scaled = transform_scaled,
    conditional  = conditional,
    n_samples    = n_samples
  )
  entry <- selected[["entry"]]
  plan  <- if (identical(selected[["component"]], provider[["components"]][["random"]])) {
    .hypothesis_model_plan_random(plan, provider, cache)
  } else if (identical(entry[["role"]], "formula_coefficient_group")) {
    .hypothesis_model_plan_linear(plan, provider, metadata, cache)
  } else {
    .hypothesis_model_plan_scalar(plan, provider, cache)
  }

  .hypothesis_model_plan_finish(plan)
}


.hypothesis_model_plan_new <- function(statement, hypothesis, parameter, label,
                                       component, selected = NULL,
                                       transform_scaled = FALSE,
                                       conditional = FALSE,
                                       n_samples = 10000L) {

  sides <- unlist(lapply(hypothesis[["statements"]], function(x) {
    c(x[["left"]][["type"]], x[["right"]][["type"]])
  }), use.names = FALSE)
  refs  <- hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )

  list(
    kind            = NA_character_,
    route           = NA_character_,
    statement       = statement,
    hypothesis      = hypothesis,
    parameter       = parameter,
    label           = label,
    component       = component,
    selected        = selected,
    transform_scaled = transform_scaled,
    conditional     = conditional,
    n_samples       = n_samples,
    point           = any(sides %in% c("point", "not_point")),
    region          = any(sides == "region"),
    refs            = refs,
    targets         = list(),
    refusal         = NULL,
    method_refusals = list(),
    support         = NULL,
    group           = NA_character_
  )
}


# The label of a point side of the statement, as the user wrote it.
.hypothesis_model_plan_side_label <- function(plan, side) {

  statement <- plan[["statement"]][["statements"]][[1L]]
  label     <- statement[[side]][["label"]]
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    label <- statement[["source"]]
  }

  label
}


# A point target: the value, its label, the prior density, the declared
# atoms, the classification of the value, and the parameter spec that the
# provider declares for the ordinate backends (NULL when the method refusals
# of the plan cover the target).
.hypothesis_model_plan_target <- function(value, label, prior_density,
                                    spec = NULL, atoms = NULL, level = NA,
                                    refusal = NULL) {

  status <- if (is.null(refusal) && !is.null(prior_density)) {
    prior_ordinate_status(prior_density, value, labels = label)
  }

  list(
    value         = value,
    label         = label,
    level         = level,
    prior_density = prior_density,
    atoms         = atoms,
    status        = status,
    spec          = spec,
    refusal       = refusal
  )
}


# The refusal of a point target: its own refusal, or its prior ordinate
# classification when the value is not eligible (BayesTools classes).
.hypothesis_model_plan_target_refusal <- function(target) {

  if (!is.null(target[["refusal"]])) {
    return(target[["refusal"]])
  }
  if (is.null(target[["prior_density"]])) {
    return(.hypothesis_model_refusal(
      paste0(
        "The prior density of '", target[["label"]], "' is unavailable, so ",
        "the point hypothesis has no Savage-Dickey Bayes factor. Use a region ",
        "or directional hypothesis."
      ),
      "target"
    ))
  }
  status <- target[["status"]]
  if (is.null(status) || all(status[["eligible"]])) {
    return(NULL)
  }
  i         <- which(!status[["eligible"]])[[1L]]
  condition <- status[["condition"]][[i]]
  reason    <- status[["reason"]][[i]]
  if (identical(condition, "BayesTools_point_mass_at_null") &&
      !is.null(target[["point_mass_reason"]])) {
    reason <- target[["point_mass_reason"]]
  }

  refusal <- .hypothesis_model_refusal(
    reason    = reason,
    condition = c(condition, "BayesTools_hypothesis_ordinate")
  )
  refusal[["value_specific"]] <- TRUE

  refusal
}

# The generic refusal of the ordinate backends that a target's unsupported
# formula transform implies, unless the plan has one already.
.hypothesis_model_plan_finish <- function(plan) {

  unsupported <- Filter(function(target){
    identical(target[["spec"]][["type"]], "unsupported_formula_transform")
  }, plan[["targets"]])
  if(length(unsupported)){
    reason <- unsupported[[1L]][["spec"]][["reason"]]
    if(is.null(plan[["method_refusals"]][["ordinate"]])){
      plan[["method_refusals"]][["ordinate"]] <-
        .hypothesis_model_cause_refusal(reason, "original_scale")
    }
  }

  plan
}




# A quantity fixed by the fitted model (catalog status or fixed value).
.hypothesis_model_plan_fixed_refusal <- function(label) {

  .hypothesis_model_refusal(
    paste0(
      "The quantity '", label, "' is fixed by the fitted model; posterior ",
      "hypothesis tests are undefined."
    ),
    "fixed"
  )
}


.hypothesis_model_plan_quantity_fixed <- function(status, fixed_value) {

  identical(as.character(status), "structural") ||
    (length(fixed_value) == 1L && is.finite(fixed_value))
}


.hypothesis_model_plan_direct_refusal <- function(plan) {

  refs <- plan[["refs"]]
  if (nrow(refs) == 0L || all(refs[["direct"]])) {
    return(NULL)
  }

  .hypothesis_model_refusal(
    paste0(
      "Point-null hypotheses require a direct parameter or level reference; ",
      "unsupported point expression in: '",
      refs[["hypothesis"]][!refs[["direct"]]][[1L]], "'."
    ),
    "statement"
  )
}


# The draws of a fixed-effect target: the mixed posteriors of the tested
# parameter and of the parameters owning the coordinates that its prior
# density combines, on the fitted ('standardized') or original scale, and the
# marginal posterior of the tested parameter. 'prior_density' is set as the
# parameter's prior density before the marginal posterior is built.
.hypothesis_model_plan_draws <- function(plan, provider, cache,
                                   extra_parameters = character(),
                                   prior_density = NULL, key = "") {

  parameter        <- plan[["parameter"]]
  sample_parameter <- unique(c(
    provider[["sample_parameters"]](parameter),
    extra_parameters
  ))
  cache_key <- paste(
    "draws", parameter, paste(sort(sample_parameter), collapse = ","),
    plan[["transform_scaled"]], plan[["conditional"]], plan[["n_samples"]], key,
    sep = "\r"
  )

  .hypothesis_model_plan_cached(cache, cache_key, function() {
    samples <- provider[["mixed_posteriors"]](
      parameters       = sample_parameter,
      conditional      = if (plan[["conditional"]]) parameter,
      conditional_rule = "AND",
      transform_scaled = plan[["transform_scaled"]],
      n_prior_samples  = plan[["n_samples"]]
    )
    density <- if (is.function(prior_density)) {
      prior_density(samples)
    } else {
      prior_density
    }
    if (!is.null(density)) {
      prior_densities <- posterior_metadata(samples, "prior_densities")
      prior_densities[[parameter]] <- density
      posterior_metadata(samples, "prior_densities") <- prior_densities
    }
    density_parameter <- provider[["density_parameter"]](
      samples, parameter, sample_parameter
    )
    posterior <- marginal_posterior(
      samples       = samples,
      parameter     = density_parameter,
      prior_samples = TRUE,
      use_formula   = FALSE,
      n_samples     = plan[["n_samples"]]
    )
    if (!is.null(density)) {
      posterior_metadata(posterior, "prior_density") <- density
    }

    list(
      samples           = samples,
      posterior         = posterior,
      sample_parameter  = sample_parameter,
      density_parameter = density_parameter,
      prior_density     = density
    )
  })
}


# The route of an original-scale formula coefficient as BayesTools declares
# it in the fitted coefficient transform: the map type (identity, affine,
# exp_affine, or unsupported), the support of the map, and the target's
# weights on the fitted coordinates. An unsupported route has a reason and
# its cause: "metadata" (the transform lacks the certified metadata),
# "fixed" (the coefficient is structurally fixed), or "nonlinear" (a
# nonlinear joint transform). NULL for parameters that are no formula
# coefficients.
.hypothesis_model_formula_route <- function(provider, selected) {

  target_info <- .hypothesis_model_formula_coefficient_target(
    provider   = provider,
    selected = selected
  )
  if (is.null(target_info)) {
    return(NULL)
  }

  .hypothesis_model_formula_route_of(provider, target_info)
}


.hypothesis_model_formula_route_of <- function(provider, target_info) {

  transform <- target_info[["transform"]]
  target    <- target_info[["target"]]
  targets   <- transform[["targets"]]
  row       <- if (is.data.frame(targets) &&
                   all(c("target", "map_type", "support") %in% names(targets))) {
    match(target, targets[["target"]])
  } else {
    NA_integer_
  }
  if (!identical(transform[["target_scale"]], "original") || is.na(row)) {
    return(c(target_info, list(
      type   = "unsupported",
      reason = paste0(
        "The fitted coefficient transform for '", target,
        "' lacks the certified structural metadata required for hypothesis testing."
      ),
      cause  = "metadata"
    )))
  }
  map_type <- targets[["map_type"]][[row]]
  dynamic <- identical(map_type, "state_dependent")
  weights <- if(dynamic) stats::setNames(numeric(), character()) else
    stats::setNames(as.numeric(transform[["matrix"]][target_info[["target_i"]], , drop = FALSE]),
      colnames(transform[["matrix"]]))
  weights <- weights[weights != 0]
  state_sources <- transform[["state_dependencies"]][["source"]][
    transform[["state_dependencies"]][["target"]] == target]
  route <- c(target_info, list(type = map_type, weights = weights,
    state_sources = stats::setNames(rep(1, length(state_sources)), state_sources),
    support = targets[["support"]][[row]]))
  if(identical(targets[["structural_status"]][[row]], "structural")){
    route[["type"]] <- "unsupported"
    route[["reason"]] <- paste0("The fitted coefficient '", target, "' is structurally fixed and has no posterior hypothesis route.")
    route[["cause"]] <- "fixed"
  }else if(!map_type %in% c("identity", "affine", "exp_affine", "state_dependent")){
    route[["type"]] <- "unsupported"
    route[["reason"]] <- paste0("The fitted nonlinear joint coefficient transform for '", target,
      "' is not supported by ", provider[["caller"]], ".")
    route[["cause"]] <- "nonlinear"
  }

  route
}
#' @title Original-scale route of a formula coefficient
#'
#' @description Describes how a coefficient of a fitted formula is related to
#' the fitted coordinates on the original scale of the predictors, from the
#' fitted coefficient transform ([JAGS_formula_coefficient_transform()]): the
#' type of the map, its support, and the weights of the coefficient on the
#' fitted coordinates. Fit-level hypothesis tests use it to evaluate
#' hypotheses on original-scale coefficients ([hypothesis_BF()]).
#'
#' @param fit a model fitted with [JAGS_fit()].
#' @param parameter the formula parameter (for example \code{"mu"}).
#' @param target the name of the fitted coefficient (for example
#'   \code{"mu_x"}).
#' @param caller the name of the function reported in the reason of a
#'   coefficient with an unsupported map.
#' @param refit_class condition classes that precede
#'   \code{BayesTools_refit_required} when the fitted metadata are outdated.
#'
#' @return A list with the \code{formula_parameter}, the \code{target}, its
#'   position \code{target_i} in the \code{transform}, the \code{type} of the
#'   map (\code{"identity"}, \code{"affine"}, \code{"exp_affine"},
#'   \code{"state_dependent"}, or \code{"unsupported"}), the \code{weights}
#'   on the fitted coordinates, the \code{support} of the coefficient, and
#'   for an unsupported map the \code{reason} and its \code{cause}
#'   (\code{"metadata"}, \code{"fixed"}, or \code{"nonlinear"}).
#'
#' @seealso [JAGS_formula_coefficient_transform()], [hypothesis_BF()]
#' @export
JAGS_formula_coefficient_route <- function(fit, parameter, target,
                                           caller = "hypothesis_BF()",
                                           refit_class = NULL){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  check_char(parameter, "parameter", check_length = 1, allow_NA = FALSE)
  check_char(target, "target", check_length = 1, allow_NA = FALSE)
  check_char(caller, "caller", check_length = 1, allow_NA = FALSE)
  provider <- list(
    fit = fit, build = "BayesTools", refit_class = refit_class, caller = caller
  )
  target_info <- .hypothesis_model_formula_target_info(provider, parameter, target)

  .hypothesis_model_formula_route_of(provider, target_info)
}


# Point values outside the open support of an original-scale formula target.
.hypothesis_model_plan_support_refusal <- function(value, support, description) {

  if (is.null(support)) {
    return(NULL)
  }
  outside <- !is.finite(value) ||
    (is.finite(support[[1L]]) && value <= support[[1L]]) ||
    (is.finite(support[[2L]]) && value >= support[[2L]])
  if (!outside) {
    return(NULL)
  }

  refusal <- .hypothesis_model_refusal(
    paste0(
      "Point-null value ", value, " is outside or on the boundary of the ",
      "open support for ", description, "."
    ),
    "target"
  )
  refusal[["value_specific"]] <- TRUE

  refusal
}


# The reason of a prior point mass at a tested value that the provider words
# (for example by naming the inclusion Bayes factor of a model-averaged
# parameter), or NULL for the reason of the prior ordinate classification.
.hypothesis_model_plan_point_mass_reason <- function(provider) {

  if (!is.function(provider[["point_mass_reason"]])) {
    return(NULL)
  }

  provider[["point_mass_reason"]]()
}


# Whether BayesTools evaluates the region probabilities of a prior density
# (prior_density_has_provenance()); a density without that
# provenance, e.g. the plotting grid of a nested variance allocation, has none.
.hypothesis_model_plan_density_has_provenance <- function(prior_density) {

  !is.null(prior_density) && prior_density_has_provenance(prior_density)
}


# The exact support bounds BayesTools declares for posterior draws, or NULL.
.hypothesis_model_plan_support_bounds <- function(x) {

  support <- posterior_metadata(x, "support")
  if (is.list(support) && !inherits(support, "BayesTools_posterior_support") &&
      length(support) == 1L) {
    support <- support[[1L]]
  }
  if (!inherits(support, "BayesTools_posterior_support") ||
      !isTRUE(support[["exact"]])) {
    return(NULL)
  }

  as.numeric(support[["bounds"]])
}

.hypothesis_model_measure_check <- function(sample, measure) {

  unavailable <- posterior_metadata(sample, "measure_unavailable")
  if(is.null(unavailable)) return(invisible(TRUE))
  entries <- unavailable[unavailable[["measure"]] == measure, , drop = FALSE]
  if(!nrow(entries)) return(invisible(TRUE))
  posterior_measure_check(
    sample, measure, column = entries[["column"]])

}

.hypothesis_model_linear_weight_space <- function(provider, sample) {

  weights <- posterior_metadata(sample, "linear_weights")
  if(is.null(weights)) return(NULL)
  space <- posterior_metadata(sample, "linear_weight_space")
  if(!is.character(space) || length(space) != 1L || is.na(space) ||
     !space %in% c("coefficient", "formula_contribution")) .bt_stop_refit_required(
    "Stored linear weights have no supported coordinate space. Recreate the posterior with the current ", provider[["build"]], " build.",
    class = provider[["refit_class"]])
  context <- posterior_metadata(sample, "prior_context")
  context_schema <- if(inherits(context, "prior_density_model_mixture_context") ||
    inherits(context, "prior_density_conditional_context")) 2L else 1L
  if(length(context[["transforms"]]) && (!identical(context[["schema_version"]], context_schema) ||
     any(vapply(context[["transforms"]], function(transform){
       !identical(transform[["descriptor"]][["schema_version"]], 4L)
     }, logical(1))))) .bt_stop_refit_required(
       "Stored formula context has unsupported transformation metadata. Recreate the posterior with the current ", provider[["build"]], " build.",
       class = provider[["refit_class"]])
  if(identical(space, "formula_contribution") &&
     (length(context[["formula_scale"]]) || length(context[["transforms"]]))) .bt_stop_refit_required(
       "Formula contribution weights cannot be remapped through coefficient scaling. Recreate the posterior with the current ", provider[["build"]], " build.",
       class = provider[["refit_class"]])
  space
}

# The linear weights of a posterior as a fully named numeric vector (from the
# vector, one-row matrix, or one-row data frame of its metadata).
.hypothesis_model_linear_weights <- function(weights) {

  if (is.null(weights)) {
    return(NULL)
  }
  if (is.data.frame(weights)) {
    if (!all(vapply(weights, is.numeric, logical(1)))) {
      stop("Linear target weights must be numeric.", call. = FALSE)
    }
    weights <- as.matrix(weights)
  }
  if (is.matrix(weights)) {
    if (nrow(weights) != 1L) {
      stop("Linear target weights must contain exactly one row.", call. = FALSE)
    }
    weights <- weights[1, ]
  }
  if (!is.numeric(weights)) {
    stop("Linear target weights must be numeric.", call. = FALSE)
  }
  weight_names <- names(weights)
  weights <- as.numeric(weights)
  names(weights) <- weight_names
  if (is.null(names(weights)) || anyNA(names(weights)) ||
      any(names(weights) == "")) {
    stop("Linear target weights must be fully named.", call. = FALSE)
  }

  weights
}


# Only declared constant products have one raw fitted linear projection.
.hypothesis_model_fitted_linear_weights <- function(provider, sample, weights) {

  .hypothesis_model_measure_check(sample, "prior_density")
  space <- .hypothesis_model_linear_weight_space(provider, sample)
  if(!identical(space, "formula_contribution")) return(weights)
  context <- posterior_metadata(sample, "prior_context")
  priors <- context[["prior_list"]]
  out <- weights
  for(owner in names(priors)){
    columns <- names(out)[names(out) == owner | startsWith(names(out), paste0(owner, "["))]
    columns <- columns[out[columns] != 0]
    if(!length(columns)) next
    multiplier <- attr(priors[[owner]], "multiply_by", exact = TRUE)
    if(is.null(multiplier)) next
    if(!is.numeric(multiplier) || length(multiplier) != 1L || !is.finite(multiplier)){
      .hypothesis_model_stop_cause(paste0(
        provider[["ordinate_label"]],
        " are unavailable for stochastic formula contribution products. ",
        "Use density_method = 'KDE'."),
        "original_scale")
    }
    out[columns] <- out[columns] * multiplier
  }
  out[out != 0]
}


#' @title Metadata checks of the linear targets of posterior draws
#'
#' @description Checks of the metadata that the draws of a formula coefficient
#' or a factor level carry about the linear map of the fitted coordinates they
#' are a function of. The fit-level hypothesis tests ([hypothesis_BF()]) and
#' the density backends ([hypothesis_density_backend()]) use them to describe a
#' target.
#'
#' \code{posterior_formula_measure_check()} stops with the declared refusal
#' (\code{BayesTools_formula_measure_unavailable}) when the draws declare the
#' \code{measure} (\code{"prior_density"} or \code{"atoms"}) unavailable.
#' \code{posterior_linear_weight_space()} returns the coordinate space of the
#' linear weights of the draws (\code{"coefficient"} or
#' \code{"formula_contribution"}, \code{NULL} without weights) and stops with
#' \code{BayesTools_refit_required} for outdated metadata.
#' \code{posterior_fitted_linear_weights()} returns the weights on the fitted
#' coordinates, with the declared constant multipliers of formula
#' contributions folded in; stochastic multipliers have no raw fitted linear
#' projection and stop with the classes
#' \code{BayesTools_hypothesis_method}, \code{BayesTools_hypothesis_unavailable},
#' and \code{BayesTools_hypothesis_cause} (the field \code{cause} is
#' \code{"original_scale"}).
#'
#' @param sample posterior draws with metadata ([posterior_metadata()]).
#' @param measure the measure to check.
#' @param weights a named numeric vector of weights on the fitted coordinates.
#' @param ordinate_label how the message names the density methods that
#'   compute posterior ordinates.
#'
#' @seealso [posterior_measure_check()], [hypothesis_BF()]
#' @rdname posterior_linear_target_checks
#' @export
posterior_formula_measure_check <- function(sample, measure){

  check_char(measure, "measure", check_length = 1, allow_NA = FALSE)

  .hypothesis_model_measure_check(sample, measure)
}


#' @rdname posterior_linear_target_checks
#' @export
posterior_linear_weight_space <- function(sample){

  .hypothesis_model_linear_weight_space(
    list(build = "BayesTools", refit_class = NULL), sample
  )
}


#' @rdname posterior_linear_target_checks
#' @export
posterior_fitted_linear_weights <- function(sample, weights,
                                            ordinate_label = "Ordinate density methods"){

  check_char(ordinate_label, "ordinate_label", check_length = 1, allow_NA = FALSE)

  .hypothesis_model_fitted_linear_weights(
    list(build = "BayesTools", refit_class = NULL,
         ordinate_label = ordinate_label),
    sample, weights
  )
}


# ---------------------------------------------------------------------------- #
# Scalar and exp(affine) targets
# ---------------------------------------------------------------------------- #

.hypothesis_model_plan_scalar <- function(plan, provider, cache) {

  plan[["kind"]] <- "scalar"
  plan[["route"]] <- "scalar"
  plan[["group"]] <- paste("scalar", plan[["parameter"]], sep = "\r")
  tryCatch(.hypothesis_model_plan_scalar_available(plan, provider, cache),
    BayesTools_formula_measure_unavailable = function(condition){
      plan[["refusal"]] <- .hypothesis_model_refusal_from_condition(condition)
      plan
    })
}

.hypothesis_model_plan_scalar_available <- function(plan, provider, cache) {

  entry <- plan[["selected"]][["entry"]]
  plan[["kind"]]  <- "scalar"
  plan[["route"]] <- "scalar"
  plan[["group"]] <- paste("scalar", plan[["parameter"]], sep = "\r")
  route <- if(plan[["transform_scaled"]]) .hypothesis_model_formula_route(provider, plan[["selected"]])
  fixed <- if(is.null(route)) .hypothesis_model_plan_quantity_fixed(entry[["status"]], entry[["fixed_value"]]) else
    identical(route[["cause"]], "fixed")
  if(fixed){
    plan[["refusal"]] <- .hypothesis_model_plan_fixed_refusal(plan[["label"]])
    return(plan)
  }
  plan[["refusal"]] <- .hypothesis_model_plan_direct_refusal(plan)
  if (!is.null(plan[["refusal"]])) {
    return(plan)
  }

  if (identical(route[["type"]], "unsupported")) {
    plan[["refusal"]] <- .hypothesis_model_refusal(route[["reason"]], "target")
    return(plan)
  }
  if (!is.null(route) && identical(route[["transform"]][["output_transforms"]][[route[["target"]]]], "exp") &&
      !identical(route[["type"]], "identity")) {
    return(.hypothesis_model_plan_exp_affine(plan, provider, route, cache))
  }

  formula_density <- if (!is.null(route)) {
    function(samples) {
      JAGS_formula_prior_density(
        fit          = provider[["fit"]],
        parameter    = route[["formula_parameter"]],
        target       = route[["target"]],
        target_scale = "original",
        context      = posterior_metadata(samples, "prior_context")
      )
    }
  }
  draws <- .hypothesis_model_plan_draws(
    plan             = plan,
    provider           = provider,
    cache            = cache,
    extra_parameters = if (!is.null(route)) {
      .hypothesis_model_target_prior_parameters(provider,
        if(identical(route[["type"]], "state_dependent")) route[["state_sources"]] else route[["weights"]])
    } else {
      character()
    },
    prior_density    = formula_density,
    key              = "scalar"
  )
  prior_density <- posterior_metadata(draws[["posterior"]], "prior_density")
  spec <- if (is.null(route) || identical(route[["type"]], "identity")) {
    list(type = "primitive", prior_density = prior_density)
  } else {
    list(type = if(identical(route[["type"]], "state_dependent")) "unsupported_formula_transform" else "linear",
      weights = route[["weights"]], prior_density = prior_density,
      reason = if(identical(route[["type"]], "state_dependent"))
        paste0(provider[["ordinate_label"]], " are unavailable for state-dependent original coefficients. Use density_method = 'KDE'."))
  }
  plan[["draws"]]         <- draws
  plan[["route_info"]]    <- route
  plan[["prior_density"]] <- prior_density
  plan[["support"]]       <- if (is.null(route[["support"]])) {
    .hypothesis_model_plan_support_bounds(draws[["posterior"]])
  } else {
    route[["support"]]
  }
  plan[["weights"]]       <- if (is.null(route)) {
    stats::setNames(1, plan[["parameter"]])
  } else {
    route[["weights"]]
  }
  plan[["targets"]] <- .hypothesis_model_plan_point_targets(
    plan          = plan,
    provider        = provider,
    prior_density = prior_density,
    spec          = spec,
    atoms         = posterior_metadata(draws[["posterior"]], "atoms"),
    support       = route[["support"]],
    description   = paste0("transformed coefficient '", plan[["parameter"]], "'")
  )
  if(identical(route[["type"]], "state_dependent")){
    plan[["method_refusals"]] <- list(ordinate = .hypothesis_model_cause_refusal(
      paste0(provider[["ordinate_label"]], " are unavailable for a state-dependent original coefficient map. Use density_method = 'KDE'."),
      "original_scale"))
  }
  if (length(plan[["targets"]]) > 0L &&
      is.list(draws[["posterior"]])) {
    plan[["method_refusals"]] <- .hypothesis_model_plan_precomputed_refusals(paste0(
      provider[["ordinate_label"]], " point hypotheses for factor parameters must specify ",
      "a level, e.g. '", plan[["label"]], "[level] = ",
      plan[["targets"]][[1L]][["value"]], "'."
    ))
  }

  plan
}


# One point target per point side of the statement of a scalar quantity.
.hypothesis_model_plan_point_targets <- function(plan, provider, prior_density,
                                           spec = NULL, atoms = NULL,
                                           support = NULL, description = NULL) {

  refs <- plan[["refs"]]
  lapply(seq_len(nrow(refs)), function(i) {
    label  <- .hypothesis_model_plan_side_label(plan, refs[["side"]][[i]])
    target <- .hypothesis_model_plan_target(
      value         = refs[["value"]][[i]],
      label         = label,
      prior_density = prior_density,
      spec          = spec,
      atoms         = atoms,
      refusal       = if (!is.null(support)) {
        .hypothesis_model_plan_support_refusal(
          refs[["value"]][[i]], support, description
        )
      }
    )
    target[["point_mass_reason"]] <- .hypothesis_model_plan_point_mass_reason(provider)
    target
  })
}


.hypothesis_model_plan_precomputed_refusals <- function(reason, type = "method") {

  list(ordinate = .hypothesis_model_refusal(reason, type))
}


# An exponentiated affine original-scale coefficient: KDE point and region
# hypotheses on the coefficient's own scale for a certified atom-free,
# unconditional target, with the prior density of the fitted coefficient
# transform.
.hypothesis_model_plan_exp_affine <- function(plan, provider, route, cache) {

  plan[["kind"]] <- "exp_affine"
  plan[["route"]] <- "scalar"
  plan[["group"]] <- paste("scalar", plan[["parameter"]], sep = "\r")
  tryCatch(.hypothesis_model_plan_exp_affine_available(plan, provider, route, cache),
    BayesTools_formula_measure_unavailable = function(condition){
      plan[["refusal"]] <- .hypothesis_model_refusal_from_condition(condition)
      plan
    })
}

.hypothesis_model_plan_exp_affine_available <- function(plan, provider, route, cache) {

  plan[["kind"]]  <- "exp_affine"
  plan[["route"]] <- "scalar"
  target <- route[["target"]]
  if (plan[["conditional"]]) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Nonlinear fitted-scale KDE hypotheses are unavailable for conditional ",
        "product-space posteriors."
      ),
      "target"
    )
    return(plan)
  }
  draws <- .hypothesis_model_plan_exp_affine_draws(plan, provider, route, cache)
  plan[["refusal"]] <- .hypothesis_model_plan_exp_affine_certify(
    sample        = draws[["samples"]][[target]],
    prior_density = draws[["prior_density"]]
  )
  if (!is.null(plan[["refusal"]])) {
    return(plan)
  }

  plan[["draws"]]         <- draws
  plan[["route_info"]]    <- route
  plan[["prior_density"]] <- draws[["prior_density"]]
  plan[["support"]]       <- route[["support"]]
  plan[["weights"]]       <- route[["weights"]]
  plan[["targets"]]       <- .hypothesis_model_plan_point_targets(
    plan          = plan,
    provider        = provider,
    prior_density = draws[["prior_density"]],
    atoms         = posterior_metadata(draws[["posterior"]], "atoms"),
    support       = route[["support"]],
    description   = paste0("transformed coefficient '", target, "'")
  )
  reason <- paste0(
    "The requested nonlinear fitted-scale hypothesis for '",
    plan[["parameter"]], "' is supported only with density_method = 'KDE'. ",
    provider[["ordinate_label"]], " ordinates support only direct parameter or level point ",
    "hypotheses with an exact linear fitted-scale map."
  )
  # The ordinate backends are refused for the cause "original_scale"; the
  # normal approximation is no ordinate request.
  plan[["method_refusals"]] <- list(
    ordinate = .hypothesis_model_cause_refusal(reason, "original_scale"),
    normal   = .hypothesis_model_refusal(reason, "method")
  )

  plan
}


# The mixed posteriors of an exp(affine) target and of the parameters owning
# the fitted coordinates its prior combines, and the target's marginal
# posterior on its original scale (marginal_posterior()): its
# draws with the prior density of the fitted coefficient transform, the exact
# support of the map, and the declared atoms and conditioning of the mixed
# posterior.
.hypothesis_model_plan_exp_affine_draws <- function(plan, provider, route, cache) {

  parameter        <- plan[["parameter"]]
  sample_parameter <- unique(c(
    provider[["sample_parameters"]](parameter),
    .hypothesis_model_target_prior_parameters(provider,
        if(identical(route[["type"]], "state_dependent")) route[["state_sources"]] else route[["weights"]])
  ))
  .hypothesis_model_plan_cached(
    cache,
    paste("exp_affine", parameter, plan[["n_samples"]], sep = "\r"),
    function() {
      samples <- provider[["mixed_posteriors"]](
        parameters       = sample_parameter,
        conditional      = NULL,
        conditional_rule = "AND",
        transform_scaled = TRUE,
        n_prior_samples  = plan[["n_samples"]]
      )
      posterior <- marginal_posterior(
        samples       = samples,
        parameter     = route[["target"]],
        prior_samples = TRUE,
        use_formula   = FALSE,
        n_samples     = plan[["n_samples"]]
      )

      list(
        samples       = samples,
        posterior     = posterior,
        prior_density = posterior_metadata(posterior, "prior_density")
      )
    }
  )
}


# Structural certification of an exp(affine) target: the averaged,
# atom-free posterior of a scalar target whose prior density is atom-free.
.hypothesis_model_plan_exp_affine_certify <- function(sample, prior_density) {

  refusal <- function(reason) .hypothesis_model_refusal(reason, "target")
  if (is.null(sample) || !inherits(sample, "mixed_posteriors.simple")) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require a certified scalar ",
      "mixed posterior."
    )))
  }
  condition <- posterior_metadata(sample, "condition")
  if (!isTRUE(condition[["averaged"]])) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence for ",
      "an unconditional posterior."
    )))
  }
  if (!posterior_atoms_free(sample)) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the posterior is atom-free."
    )))
  }
  points <- prior_density[["points"]]
  if (!inherits(prior_density, "prior_density") || !is.data.frame(points) ||
      !all(c("x", "p") %in% names(points)) || nrow(points) > 0L) {
    return(refusal(paste0(
      "Nonlinear fitted-scale KDE hypotheses require structural evidence that ",
      "the prior is atom-free."
    )))
  }

  NULL
}


# ---------------------------------------------------------------------------- #
# Linear targets: factor levels and linear combinations of levels
# ---------------------------------------------------------------------------- #

# Levels of every contrast are linear in the fitted coordinates. A statement
# whose point sides reference levels directly ('g[a] = 0', 'g[a] = 0 vs
# g[b] = 0.1') and a region statement ('g[a] > g[b]') are evaluated on the
# marginal posterior of the factor term ("levels" route); a statement whose
# point sides are linear combinations of levels (level contrasts such as
# 'g[a] = g[b]' and derived linear expressions such as '2 * g[a] = 0.1')
# is evaluated on the scalar linear target of BayesTools
# hypothesis_linear_target() ("combination" route). On the original scale of
# formula coefficients, a target's prior density is
# JAGS_formula_prior_density() at the target's weights on the
# original coefficients; on the fitted scale it is the density BayesTools
# attaches to the level or linear target.
.hypothesis_model_plan_linear <- function(plan, provider, metadata, cache) {

  plan[["kind"]] <- "linear"
  entry     <- plan[["selected"]][["entry"]]
  parameter <- plan[["parameter"]]
  levels    <- .hypothesis_model_plan_term_levels(metadata, entry)
  occurrences <- hypothesis_symbols(
    plan[["hypothesis"]],
    occurrences = TRUE
  )
  referenced <- unique(occurrences[["level"]][occurrences[["parameter"]] == parameter])
  whole_term <- any(is.na(referenced))
  referenced <- referenced[!is.na(referenced)]
  unknown <- setdiff(referenced, levels[["level"]])
  if (length(unknown) > 0L) {
    plan[["refusal"]] <- .hypothesis_model_unknown_level_refusal(
      level     = unknown[[1L]],
      parameter = plan[["label"]]
    )
    return(plan)
  }
  # A statement tests a level the contrast fixes (the treatment reference
  # level) when it states a point on it directly, references the whole term
  # (a point or region event on every level jointly), or references no other
  # level; linear combinations and comparisons with other levels remain
  # defined.
  fixed_levels <- levels[["level"]][levels[["fixed"]]]
  refs         <- plan[["refs"]]
  tested       <- unique(c(
    refs[["level"]][refs[["direct"]] & !is.na(refs[["level"]])],
    if (whole_term) levels[["level"]],
    if (length(referenced) > 0L && all(referenced %in% fixed_levels)) referenced
  ))
  fixed <- intersect(tested, fixed_levels)
  if (length(fixed) > 0L) {
    plan[["refusal"]] <- .hypothesis_model_plan_fixed_refusal(
      paste0(plan[["label"]], "[", fixed[[1L]], "]")
    )
    return(plan)
  }

  scale <- .hypothesis_model_plan_linear_scale(
    plan        = plan,
    provider      = provider,
    entry       = entry,
    coordinates = unique(unlist(levels[["dependencies"]], use.names = FALSE))
  )
  draws <- .hypothesis_model_plan_draws(
    plan             = plan,
    provider           = provider,
    cache            = cache,
    extra_parameters = scale[["prior_parameters"]],
    key              = "linear"
  )
  posterior <- draws[["posterior"]]
  if (!is.list(posterior)) {
    stop("Internal error: factor levels have no factor marginal posterior.",
         call. = FALSE)
  }
  plan[["draws"]] <- draws
  plan[["scale"]] <- scale

  combination <- nrow(plan[["refs"]]) > 0L && any(!plan[["refs"]][["direct"]])
  if (combination) {
    return(.hypothesis_model_plan_combination(plan, provider))
  }

  plan[["route"]] <- "levels"
  plan[["group"]] <- .hypothesis_model_plan_levels_group(plan, parameter)
  refs <- plan[["refs"]]
  if (length(referenced) > 0L) {
    plan[["support"]] <- .hypothesis_model_plan_support_bounds(
      posterior[[referenced[[1L]]]]
    )
  }
  plan[["targets"]] <- unlist(lapply(seq_len(nrow(refs)), function(i) {
    ref_levels <- if (is.na(refs[["level"]][[i]])) {
      levels[["level"]]
    } else {
      refs[["level"]][[i]]
    }
    label <- .hypothesis_model_plan_side_label(plan, refs[["side"]][[i]])
    lapply(ref_levels, function(level) {
      .hypothesis_model_plan_level_target(
        plan   = plan,
        provider = provider,
        level  = level,
        value  = refs[["value"]][[i]],
        label  = if (is.na(refs[["level"]][[i]])) {
          paste0(plan[["label"]], "[", level, "] = ", refs[["value"]][[i]])
        } else {
          label
        }
      )
    })
  }), recursive = FALSE)
  if (whole_term && plan[["point"]]) {
    plan[["method_refusals"]] <- .hypothesis_model_plan_precomputed_refusals(paste0(
      provider[["ordinate_label"]], " point hypotheses for factor parameters must specify ",
      "a level, e.g. '", plan[["label"]], "[level] = ",
      refs[["value"]][[1L]], "'."
    ))
  }

  plan
}


# The levels of a factor term: their labels and whether the contrast fixes
# them (the treatment reference level), from the catalog.
.hypothesis_model_plan_term_levels <- function(metadata, entry) {

  quantities <- metadata[["catalog"]][["quantities"]]
  members    <- quantities[
    quantities[["quantity_id"]] %in%
      unlist(entry[["member_quantity_ids"]], use.names = FALSE),
    ,
    drop = FALSE
  ]

  out <- data.frame(
    level       = members[["component"]],
    quantity_id = members[["quantity_id"]],
    fixed       = vapply(seq_len(nrow(members)), function(i) {
      .hypothesis_model_plan_quantity_fixed(
        members[["status"]][[i]],
        members[["fixed_value"]][[i]]
      )
    }, logical(1)),
    stringsAsFactors = FALSE
  )
  out[["dependencies"]] <- lapply(members[["extraction_key"]], function(key) {
    as.character(key[["dependencies"]])
  })

  out
}


# The coefficient scale of the targets of a factor term: on the original
# scale of formula coefficients, the fitted coefficient transform (its
# matrix maps original-coefficient weights to fitted-coordinate weights) and
# the prior-list entries owning the fitted coordinates it combines, whose
# priors the prior densities combine.
.hypothesis_model_plan_linear_scale <- function(plan, provider, entry, coordinates) {

  formula_parameter <- entry[["formula_parameter"]]
  if (!plan[["transform_scaled"]] || !is.character(formula_parameter) ||
      length(formula_parameter) != 1L || is.na(formula_parameter) ||
      !nzchar(formula_parameter)) {
    return(list(original = FALSE, prior_parameters = character()))
  }
  transform <- JAGS_formula_coefficient_transform(
    fit          = provider[["fit"]],
    parameter    = formula_parameter,
    target_scale = "original"
  )
  term_coordinates <- intersect(coordinates, rownames(transform[["matrix"]]))
  static <- term_coordinates[transform[["targets"]][["map_type"]][
    match(term_coordinates, transform[["targets"]][["target"]])] != "state_dependent"]
  dependencies <- if(length(static)) colnames(transform[["matrix"]])[
    colSums(abs(transform[["matrix"]][static, , drop = FALSE])) > 0] else character()
  dependencies <- union(dependencies, transform[["state_dependencies"]][["source"]][
    transform[["state_dependencies"]][["target"]] %in% term_coordinates])

  list(
    original          = TRUE,
    formula_parameter = formula_parameter,
    transform         = transform,
    prior_parameters  = if (length(dependencies) > 0L) {
      .hypothesis_model_target_prior_parameters(
        provider,
        stats::setNames(rep(1, length(dependencies)), dependencies)
      )
    } else {
      character()
    }
  )
}


# A linear target's prior density and its weights on the fitted coordinates
# from its weights 'weights' on the coefficients of its posterior (original
# coefficients on the original scale, fitted coordinates otherwise).
.hypothesis_model_plan_linear_density <- function(plan, provider, weights, density, sample = NULL) {

  weights <- weights[weights != 0]
  if(!is.null(sample) && identical(.hypothesis_model_linear_weight_space(provider, sample), "formula_contribution")){
    folded <- tryCatch(.hypothesis_model_fitted_linear_weights(provider, sample, weights),
      BayesTools_hypothesis_cause = function(condition) condition)
    if(inherits(folded, "condition")) return(list(prior_density = density,
      weights = numeric(), static = FALSE, reason = conditionMessage(folded)))
    return(list(prior_density = density, weights = folded, static = TRUE))
  }
  scale   <- plan[["scale"]]
  if (!isTRUE(scale[["original"]])) {
    return(list(prior_density = density, weights = weights))
  }
  matrix <- scale[["transform"]][["matrix"]]
  selected <- matrix[names(weights), , drop = FALSE]
  static <- all(is.finite(selected))
  fitted <- if(static) drop(matrix(weights, nrow = 1L) %*% selected) else numeric()
  if(static) names(fitted) <- colnames(matrix)

  list(
    prior_density = JAGS_formula_prior_density(
      fit          = provider[["fit"]],
      parameter    = scale[["formula_parameter"]],
      weights      = weights,
      target_scale = "original",
      context      = posterior_metadata(
        plan[["draws"]][["samples"]],
        "prior_context"
      )
    ),
    weights = fitted[fitted != 0], static = static
  )
}


.hypothesis_model_plan_level_target <- function(plan, provider, level, value, label) {

  tryCatch(.hypothesis_model_plan_level_target_available(plan, provider, level, value, label),
    BayesTools_formula_measure_unavailable = function(condition){
      .hypothesis_model_plan_target(value = value, label = label, level = level, prior_density = NULL,
        refusal = .hypothesis_model_refusal_from_condition(condition))
    })
}

.hypothesis_model_plan_level_target_available <- function(plan, provider, level, value, label) {

  sample  <- plan[["draws"]][["posterior"]][[level]]
  .hypothesis_model_measure_check(sample, "prior_density")
  .hypothesis_model_measure_check(sample, "atoms")
  .hypothesis_model_linear_weight_space(provider, sample)
  weights <- .hypothesis_model_linear_weights(
    posterior_metadata(sample, "linear_weights")
  )
  if (length(weights) == 0L) {
    # Fitted metadata of an older build, not a refusal of the target.
    .bt_stop_refit_required(
      "The linear weights of factor level '", plan[["label"]], "[", level,
      "]' on the fitted coefficients are unavailable. Refit the model ",
      "with the current ", provider[["build"]], " build.",
      class = provider[["refit_class"]]
    )
  }
  density <- .hypothesis_model_plan_linear_density(
    plan    = plan,
    provider  = provider,
    weights = weights,
    density = posterior_metadata(sample, "prior_density"),
    sample = sample
  )
  target <- .hypothesis_model_plan_target(
    value         = value,
    label         = label,
    level         = level,
    prior_density = density[["prior_density"]],
    atoms         = posterior_metadata(sample, "atoms"),
    spec          = list(
      type          = if(isFALSE(density[["static"]])) "unsupported_formula_transform" else "linear",
      weights       = density[["weights"]],
      prior_density = density[["prior_density"]],
      reason        = if(isFALSE(density[["static"]])) if(is.null(density[["reason"]]))
        paste0(provider[["ordinate_label"]], " are unavailable for state-dependent original coefficients. Use density_method = 'KDE'.") else density[["reason"]]
    )
  )
  target[["point_mass_reason"]] <- .hypothesis_model_plan_point_mass_reason(provider)

  target
}


# A statement on one linear combination of the levels: the scalar linear
# target of BayesTools, which certifies the target's structural atom-freeness
# and joint prior context (classed refusals otherwise).
.hypothesis_model_plan_combination <- function(plan, provider) {

  tryCatch(.hypothesis_model_plan_combination_available(plan, provider),
    BayesTools_formula_measure_unavailable = function(condition){
      plan$route <- "combination"
      plan$group <- paste("combination", plan$parameter, sep = "\r")
      plan$refusal <- .hypothesis_model_refusal_from_condition(condition)
      plan
    })
}

# The group of a linear combination of factor levels: statements are
# evaluated together exactly when BayesTools accepts them as one linear
# target, i.e., when the level coefficients of their sides (scaled as
# BayesTools scales them, in the statement's own level expression) are
# identical. The fitted-coordinate weights of the target are not the key:
# different level expressions can give the same weights.
.hypothesis_model_plan_combination_group <- function(plan) {

  hypothesis <- plan[["hypothesis"]]
  if (!inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    hypothesis <- hypothesis_parse(hypothesis)
  }
  coefficients <- hypothesis_linear_target_side(
    hypothesis[["statements"]][[1L]][["left"]], plan[["parameter"]]
  )[["coefficients"]]

  paste(
    "combination", plan[["parameter"]],
    rlang::hash(list(names(coefficients), sprintf("%a", coefficients))),
    sep = "\r"
  )
}


.hypothesis_model_plan_combination_available <- function(plan, provider) {

  plan[["route"]] <- "combination"
  target <- tryCatch(
    hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = plan[["hypothesis"]],
      parameter  = plan[["parameter"]]
    ),
    error = .hypothesis_model_linear_target_condition
  )
  if (inherits(target, "condition")) {
    reason <- conditionMessage(target)
    if (inherits(target, "BayesTools_linear_target_unavailable") &&
        identical(target[["reason"]], "posterior_atoms") &&
        is.function(provider[["combination_hint"]])) {
      reason <- paste0(reason, provider[["combination_hint"]](plan))
    }
    plan[["refusal"]] <- .hypothesis_model_refusal_from_condition(target, reason)
    plan[["group"]] <- paste("combination", plan[["parameter"]], sep = "\r")
    return(plan)
  }
  weights <- .hypothesis_model_linear_weights(target[["weights"]])
  density <- .hypothesis_model_plan_linear_density(
    plan    = plan,
    provider  = provider,
    weights = weights,
    density = posterior_metadata(target[["posterior"]], "prior_density"),
    sample = target[["posterior"]]
  )
  posterior_metadata(target[["posterior"]], "prior_density") <-
    density[["prior_density"]]
  plan[["linear_target"]] <- target
  plan[["group"]] <- .hypothesis_model_plan_combination_group(plan)

  refs <- hypothesis_parse_point_reference(
    hypothesis     = target[["hypothesis"]],
    allow_compound = TRUE
  )
  spec <- list(
    type = if(isFALSE(density$static)) "unsupported_formula_transform" else "linear",
    weights = density$weights, prior_density = density$prior_density,
    reason = if(isFALSE(density$static)) if(is.null(density$reason))
      paste0(provider[["ordinate_label"]], " are unavailable for state-dependent original coefficients. Use density_method = 'KDE'.") else density$reason
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    .hypothesis_model_plan_target(
      value         = refs[["value"]][[i]],
      label         = .hypothesis_model_plan_side_label(plan, refs[["side"]][[i]]),
      prior_density = density[["prior_density"]],
      atoms         = posterior_metadata(target[["posterior"]], "atoms"),
      spec          = spec
    )
  })

  plan
}

.hypothesis_model_formula_coefficient_target <- function(
    provider, selected) {

  entry <- selected[["entry"]]
  if (selected[["component"]] %in% c(provider[["components"]][["random"]],
                                     names(provider[["components"]][["refused"]])) ||
      is.null(entry) ||
      identical(entry[["role"]], "formula_coefficient_group")) {
    return(NULL)
  }
  formula_parameter <- entry[["formula_parameter"]]
  if (is.null(formula_parameter) || length(formula_parameter) != 1L ||
      is.na(formula_parameter) || !nzchar(formula_parameter)) {
    return(NULL)
  }

  .hypothesis_model_formula_target_info(
    provider, formula_parameter, selected[["parameter"]]
  )
}


# The fitted coefficient transform of a formula parameter on the original
# scale, and the position of the coefficient 'target' in it.
.hypothesis_model_formula_target_info <- function(provider, formula_parameter,
                                                  target) {

  transform <- JAGS_formula_coefficient_transform(
    fit          = provider[["fit"]],
    parameter    = formula_parameter,
    target_scale = "original"
  )
  if (!inherits(transform, "BayesTools_formula_coefficient_transform") ||
      !identical(transform[["schema_version"]], 4L)) {
    .bt_stop_refit_required(
      "Formula coefficient transformation metadata are unsupported. Refit ",
      "the model with the current BayesTools version.",
      class = provider[["refit_class"]]
    )
  }
  target_i <- match(target, transform[["target_names"]])
  if (is.na(target_i)) {
    .bt_stop_refit_required(
      "Resolved formula coefficient '", target,
      "' is absent from the fitted coefficient transformation.",
      class = provider[["refit_class"]]
    )
  }

  return(list(
    formula_parameter = formula_parameter,
    target            = target,
    target_i          = target_i,
    transform         = transform
  ))
}


# The fitted prior-list entries (formula terms) that own the coordinates a
# formula coefficient target weights, from the coordinate table and the
# formula name maps.
.hypothesis_model_target_prior_parameters <- function(provider, weights) {

  coordinates <- parameter_coordinates(provider[["fit"]])
  rows <- match(names(weights), coordinates[["coordinate_name"]])
  if (anyNA(rows)) {
    .bt_stop_refit_required(
      "Fitted coordinate metadata of the hypothesis target are unavailable. ",
      "Refit the model with the current ", provider[["build"]], " build.",
      class = provider[["refit_class"]]
    )
  }
  coordinates <- coordinates[rows, , drop = FALSE]

  ordinary <- is.na(coordinates[["formula_parameter"]]) | !nzchar(coordinates[["formula_parameter"]])
  ordinary_owners <- intersect(coordinates[["monitor_name"]][ordinary], names(attr(provider[["fit"]], "prior_list", exact = TRUE)))
  formula_coordinates <- coordinates[!ordinary, , drop = FALSE]
  union(ordinary_owners, unique(unlist(lapply(
    unique(formula_coordinates[["formula_parameter"]]),
    function(formula_parameter) {
      name_map <- provider[["formula_name_map"]](formula_parameter)
      fixed    <- name_map[name_map[["kind"]] == "fixed", , drop = FALSE]
      terms    <- formula_coordinates[["term"]][
        formula_coordinates[["formula_parameter"]] == formula_parameter
      ]
      fixed[["jags_name"]][fixed[["term"]] %in% terms]
    }
  ), use.names = FALSE)))
}
