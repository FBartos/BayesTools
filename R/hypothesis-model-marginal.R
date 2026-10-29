# ============================================================================ #
# hypothesis-model-marginal.R
# ============================================================================ #
#
# Plans of statements on the stored marginal means of a model (a provider of
# kind "marginal_means"), and the helpers of their evaluation.
#
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Marginal means
# ---------------------------------------------------------------------------- #

# The plan of one statement on a marginal-means provider. Model-averaged point
# statements use the alternative-conditioned marginal means (their levels
# need not share one conditioning event, so a point statement tests one
# level); other statements use the averaged marginal means. A level's prior
# density is the one BayesTools attached to the marginal mean; linear
# combinations of levels are the scalar linear targets of
# hypothesis_linear_target(). Ordinate backends compute their ordinates from
# the stored source model (the capability object of the provider).
.hypothesis_model_plan_marginal_means <- function(provider, statement,
                                            parameter = NULL, cache = NULL) {

  .hypothesis_model_plan_cached(
    cache,
    paste("marginal_means", rlang::hash(list(statement, parameter)), sep = "\r"),
    function() {
      .hypothesis_model_plan_marginal_means_build(provider, statement, parameter)
    }
  )
}


.hypothesis_model_plan_marginal_means_build <- function(provider, statement,
                                                  parameter = NULL) {

  if (!inherits(statement, "BayesTools_hypothesis_ast")) {
    statement <- hypothesis_parse(statement)
  }
  selected <- provider[["select"]](statement, parameter)
  parameter <- selected[["parameter"]]
  plan <- .hypothesis_model_plan_new(
    statement  = statement,
    hypothesis = .hypothesis_model_rewrite(
      hypothesis = statement,
      aliases    = selected[["aliases"]],
      parameter  = parameter
    ),
    parameter  = parameter,
    label      = .hypothesis_model_alias_label(selected[["aliases"]], parameter),
    component  = "marginal_means",
    selected   = selected
  )
  plan[["kind"]]              <- "marginal_mean"
  plan[["route"]]             <- "levels"
  plan[["group"]]             <- .hypothesis_model_plan_levels_group(plan, parameter)
  model_averaged <- isTRUE(provider[["model_averaged"]])
  plan[["inference_type"]] <- if (model_averaged && plan[["point"]]) {
    "conditional"
  } else {
    "averaged"
  }
  if (model_averaged && plan[["point"]] && plan[["region"]]) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Model-averaged marginal-means hypotheses cannot mix point and ",
        "region events because they require different posterior ",
        "conditioning. Use a pure point-null or a pure region hypothesis."
      ),
      "statement"
    )
    return(.hypothesis_model_plan_finish(plan))
  }
  if (model_averaged && plan[["point"]] && length(unique(
    .hypothesis_model_marginal_ast_levels(plan[["hypothesis"]], parameter)
  )) > 1L) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Point hypotheses spanning multiple marginal-means levels are not ",
        "supported because the levels need not share one conditioning event."
      ),
      "statement"
    )
    return(.hypothesis_model_plan_finish(plan))
  }
  posterior <- provider[["inference"]][[plan[["inference_type"]]]][[parameter]]
  if (is.null(posterior)) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "'marginal_means' object does not contain ",
        plan[["inference_type"]], " marginal means."
      ),
      "target"
    )
    return(.hypothesis_model_plan_finish(plan))
  }
  plan[["draws"]] <- list(posterior = posterior)

  plan[["method_refusals"]] <- list(normal = .hypothesis_model_refusal(
    "'normal' density_method is not supported for marginal-means hypotheses.",
    "method"
  ))
  # Without its fitted source model, the ordinate backends are refused before
  # that model's capability is checked.
  if (is.null(provider[["capability_object"]])) {
    plan[["source_refusals"]] <- .hypothesis_model_plan_precomputed_refusals(
      provider[["source_missing_reason"]], "target"
    )
  }

  refs <- plan[["refs"]]
  plan <- if (nrow(refs) > 0L && any(!refs[["direct"]])) {
    .hypothesis_model_plan_marginal_combination(plan, provider)
  } else {
    .hypothesis_model_plan_marginal_levels(plan, provider)
  }

  .hypothesis_model_plan_finish(plan)
}


.hypothesis_model_plan_levels_group <- function(plan, parameter) {

  occurrences <- hypothesis_symbols(plan[["hypothesis"]], occurrences = TRUE)
  levels <- occurrences[["level"]][occurrences[["parameter"]] == parameter]
  mode <- if (anyNA(levels)) "whole" else "explicit"
  paste("levels", parameter, mode, sep = "\r")
}


# The refusal of a model-averaged marginal-means request as a whole, or NULL.
# Its statements are evaluated on one posterior conditioning (point
# statements on the alternative-conditioned marginal mean of one level,
# region statements on the averaged marginal means), so a request cannot
# combine point and region statements or point statements on several
# marginal means, whether the levels belong to one marginal-means parameter or
# to several. A certified linear-expression target without an explicit level
# is counted as the pair (parameter, "").
.hypothesis_model_plan_marginal_means_request_refusal <- function(plans, provider) {

  if (!isTRUE(provider[["model_averaged"]])) {
    return(NULL)
  }
  point <- vapply(plans, `[[`, logical(1), "point")
  if (any(point) && !all(point)) {
    return(.hypothesis_model_refusal(
      paste0(
        "A model-averaged marginal-means hypothesis request cannot mix ",
        "point-null and region statements because they require different ",
        "posterior conditioning."
      ),
      "statement"
    ))
  }
  levels <- unique(unlist(lapply(plans, function(plan) {
    vapply(plan[["targets"]], function(target) {
      level <- target[["level"]]
      if (is.null(level) || is.na(level)) level <- ""
      paste(plan[["parameter"]], level, sep = "\r")
    }, character(1L))
  })))
  if (all(point) && length(levels) > 1L) {
    return(.hypothesis_model_refusal(
      paste0(
        "Point hypotheses spanning multiple marginal-means levels are not ",
        "supported because the levels need not share one conditioning event."
      ),
      "statement"
    ))
  }

  NULL
}


# A marginal mean the fitted model fixes: its declared atoms carry the whole
# posterior mass.
.hypothesis_model_plan_draws_fixed <- function(sample) {

  atoms <- posterior_metadata(sample, "atoms")
  mass  <- atoms[["mass"]]
  locations <- atoms[["locations"]]
  isTRUE(atoms[["declared"]]) && !isFALSE(atoms[["joint_declared"]]) && length(mass) > 0L &&
    is.matrix(locations) && ncol(locations) == 1L && length(unique(as.numeric(locations))) == 1L &&
    abs(sum(mass) - 1) <= sqrt(.Machine[["double.eps"]])
}


.hypothesis_model_plan_marginal_label <- function(plan, level) {

  paste0(plan[["label"]], "[", level, "]")
}


.hypothesis_model_plan_marginal_levels <- function(plan, provider) {

  posterior <- plan[["draws"]][["posterior"]]
  parameter <- plan[["parameter"]]
  levels    <- names(posterior)
  occurrences <- hypothesis_symbols(
    plan[["hypothesis"]],
    occurrences = TRUE
  )
  referenced <- unique(occurrences[["level"]][occurrences[["parameter"]] == parameter])
  whole_term <- any(is.na(referenced))
  referenced <- referenced[!is.na(referenced)]
  unknown <- setdiff(referenced, levels)
  if (length(unknown) > 0L) {
    plan[["refusal"]] <- .hypothesis_model_unknown_level_refusal(
      level     = unknown[[1L]],
      parameter = parameter
    )
    return(plan)
  }
  sample_of <- function(level) {
    posterior[[level]]
  }
  level_fixed <- function(level) {
    .hypothesis_model_plan_draws_fixed(sample_of(level))
  }
  # As for the levels of a fitted factor (.hypothesis_model_plan_linear): a
  # statement tests a fixed level when it states a point on it directly,
  # references the whole term (an event on every level jointly), or
  # references no other level; comparisons of a fixed level with other levels
  # remain defined. Only direct references reach this route, so the level of
  # every point target below is among the tested ones: a fixed marginal mean
  # is refused here, before any target is built.
  refs   <- plan[["refs"]]
  tested <- unique(c(
    refs[["level"]][refs[["direct"]] & !is.na(refs[["level"]])],
    if (whole_term || length(referenced) == 0L) levels,
    if (length(referenced) > 0L && all(vapply(referenced, level_fixed, logical(1)))) referenced
  ))
  fixed  <- tested[vapply(tested, level_fixed, logical(1))]
  if (length(fixed) > 0L) {
    plan[["refusal"]] <- .hypothesis_model_plan_fixed_refusal(
      .hypothesis_model_plan_marginal_label(plan, fixed[[1L]])
    )
    return(plan)
  }

  # A factor marginal mean with one level resolves statements without a
  # level to that level.
  if (length(posterior) == 1L) {
    refs[["level"]][is.na(refs[["level"]])] <- levels
  }
  if (any(is.na(refs[["level"]]))) {
    plan[["method_refusals"]] <- c(
      plan[["method_refusals"]],
      .hypothesis_model_plan_precomputed_refusals(paste0(
        provider[["ordinate_label"]], " point hypotheses for marginal-means factor parameters ",
        "must specify a level, e.g. '", plan[["label"]], "[level] = ",
        refs[["value"]][is.na(refs[["level"]])][[1L]], "'."
      ))
    )
  }
  method_refusals <- plan[["method_refusals"]]
  plan[["targets"]] <- unlist(lapply(seq_len(nrow(refs)), function(i) {
    ref_levels <- if (is.na(refs[["level"]][[i]])) {
      levels
    } else {
      refs[["level"]][[i]]
    }
    label <- .hypothesis_model_plan_side_label(plan, refs[["side"]][[i]])
    lapply(ref_levels, function(level) {
      sample  <- sample_of(level)
      density <- posterior_metadata(sample, "prior_density")
      if(plan[["point"]]){
        weights <- .hypothesis_model_linear_weights(posterior_metadata(sample, "linear_weights"))
        if(length(weights)){
          available <- tryCatch(.hypothesis_model_fitted_linear_weights(provider, sample, weights),
            BayesTools_fitted_linear_weights_unavailable = function(condition) condition,
            BayesTools_formula_measure_unavailable = function(condition) condition)
          if(inherits(available, "condition")){
            # a stochastic multiplier is the method refusal for its cause (the
            # backend adds its classes); a declared measure refusal keeps its
            # classes and fields after those of the method refusal
            refusal <- if(inherits(available, "BayesTools_fitted_linear_weights_unavailable")){
              .hypothesis_model_weights_refusal(available)
            }else{
              refusal <- .hypothesis_model_refusal_from_condition(available)
              refusal[["class"]] <- unique(c(.hypothesis_model_refusal(conditionMessage(available), "method")[["class"]], refusal[["class"]]))
              refusal
            }
            if(is.null(method_refusals[["ordinate"]])) method_refusals[["ordinate"]] <<- refusal
          }
        }
      }
      .hypothesis_model_plan_target(
        value         = refs[["value"]][[i]],
        label         = label,
        level         = level,
        prior_density = density,
        atoms         = posterior_metadata(sample, "atoms"),
        refusal       = if (is.null(density)) {
          .hypothesis_model_refusal(
            paste0(
              "The prior density of marginal mean '",
              .hypothesis_model_plan_marginal_label(plan, level),
              "' is unavailable, so point hypotheses on it have no ",
              "Savage-Dickey Bayes factor. Use a region or directional ",
              "hypothesis."
            ),
            "target"
          )
        }
      )
    })
  }), recursive = FALSE)
  if (length(plan[["targets"]]) > 0L) {
    plan[["method_refusals"]] <- method_refusals
    plan[["support"]] <- .hypothesis_model_plan_support_bounds(
      sample_of(plan[["targets"]][[1L]][["level"]])
    )
  }

  plan
}


# The condition with which hypothesis_linear_target() refuses a
# statement, matched by its classes: a refit request of BayesTools
# ("BayesTools_refit_required") is no refusal of the statement and stops
# unconverted; BayesTools' refusals of the statement's references (its
# resolution errors, "BayesTools_parameter_resolution_error": an unknown
# level, "BayesTools_parameter_not_found", a statement without parameter
# symbols, "BayesTools_hypothesis_no_parameters", or a level of another
# component, "BayesTools_hypothesis_component_mismatch") are statement errors
# with BayesTools' classes and fields, as at the other entry points
# (.hypothesis_model_not_found_class()); other classed BayesTools conditions
# (linear targets BayesTools cannot certify,
# "BayesTools_linear_target_unavailable", or stale draw metadata) keep their
# classes. Declared measure refusals are handled narrowly; metadata and
# program errors propagate. Its unclassed grammar errors refuse a statement (e.g. a
# nonlinear expression of the levels), which is restated
# ("BayesTools_hypothesis_statement").
.hypothesis_model_linear_target_condition <- function(condition) {

  if (inherits(condition, "BayesTools_refit_required")) {
    stop(condition)
  }
  if (inherits(condition, "BayesTools_parameter_resolution_error")) {
    class(condition) <- unique(c("BayesTools_hypothesis_statement", class(condition)))
    return(condition)
  }
  if (inherits(condition, "BayesTools_formula_measure_unavailable") ||
      inherits(condition, "BayesTools_linear_target_unavailable")) {
    return(condition)
  }
  if (any(startsWith(class(condition), "BayesTools_"))) stop(condition)

  structure(
    class = c("BayesTools_hypothesis_statement", "error", "condition"),
    list(message = conditionMessage(condition), call = NULL)
  )
}


.hypothesis_model_plan_marginal_combination <- function(plan, provider) {

  tryCatch(.hypothesis_model_plan_marginal_combination_available(plan, provider),
    BayesTools_formula_measure_unavailable = function(condition){
      plan$route <- "combination"
      plan$group <- paste("combination", plan$parameter, sep = "\r")
      plan$refusal <- .hypothesis_model_refusal_from_condition(condition)
      plan
    })
}

.hypothesis_model_plan_marginal_combination_available <- function(plan, provider) {

  plan[["route"]] <- "combination"
  plan[["group"]] <- paste("combination", plan[["parameter"]], sep = "\r")
  target <- tryCatch(
    hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = plan[["hypothesis"]],
      parameter  = plan[["parameter"]]
    ),
    error = .hypothesis_model_linear_target_condition
  )
  if (inherits(target, "condition")) {
    plan[["refusal"]] <- .hypothesis_model_refusal_from_condition(target)
    return(plan)
  }
  weights <- .hypothesis_model_linear_weights(target[["weights"]])
  density <- posterior_metadata(target[["posterior"]], "prior_density")
  fitted <- tryCatch(.hypothesis_model_fitted_linear_weights(provider, target$posterior, weights),
    BayesTools_fitted_linear_weights_unavailable = function(condition) condition)
  spec <- if(inherits(fitted, "condition")) list(type = "unsupported_formula_transform",
    weights = numeric(), prior_density = density, reason = conditionMessage(fitted)) else
    list(type = "linear", weights = fitted, prior_density = density)
  plan[["linear_target"]] <- target
  plan[["group"]] <- .hypothesis_model_plan_combination_group(plan)
  refs <- hypothesis_parse_point_reference(
    hypothesis     = target[["hypothesis"]],
    allow_compound = TRUE
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    .hypothesis_model_plan_target(
      value         = refs[["value"]][[i]],
      label         = .hypothesis_model_plan_side_label(plan, refs[["side"]][[i]]),
      prior_density = density,
      atoms         = posterior_metadata(target[["posterior"]], "atoms"),
      spec          = spec
    )
  })

  plan
}
.hypothesis_model_marginal_ast_levels <- function(hypothesis, parameter) {

  find_levels <- function(node) {

    if (!is.list(node)) {
      return(character())
    }
    if (identical(node[["type"]], "level_reference") &&
        identical(node[["parameter"]], parameter)) {
      return(as.character(node[["level"]]))
    }

    return(unlist(lapply(node, find_levels), use.names = FALSE))
  }

  return(find_levels(hypothesis[["statements"]]))
}

.hypothesis_model_marginal_resolve_singleton_levels <- function(
    samples, point_refs) {

  missing_level <- is.na(point_refs[["level"]])
  if (length(samples) != 1L || !any(missing_level)) {
    return(point_refs)
  }
  level <- names(samples)
  if (length(level) != 1L || is.na(level) || !nzchar(level)) {
    return(point_refs)
  }
  point_refs[["level"]][missing_level] <- level

  return(point_refs)
}
