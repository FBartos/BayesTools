# ============================================================================ #
# hypothesis-model-random.R
# ============================================================================ #
#
# Random-effect routes: statements on the random-effect quantities of a model
# (SDs, variances, correlations, and allocation quantities), evaluated on the
# mixed posterior that the provider extracts.
#
# ============================================================================ #


# Statements on a random-effect quantity, evaluated on its BayesTools mixed
# posterior. Point statements on a variance are evaluated through its
# standard deviation (square-root values); the table then shows the variance
# statements and densities. Statements comparing a variance point with a
# region evaluate the point part through the standard deviation and the
# region part on the variance draws.
.hypothesis_model_plan_execute_random <- function(plans, hypothesis, provider,
                                                  backend, logBF, BF01, seed,
                                                  density_control, columns) {

  point_plans <- plans[vapply(plans, `[[`, logical(1), "point")]
  plan      <- if(length(point_plans)) point_plans[[1L]] else plans[[1L]]
  parameter <- plan[["parameter"]]
  samples   <- plan[["draws"]][["samples"]]
  selected  <- plan[["random"]]
  random    <- provider[["random"]]
  defined_footnote <- random[["defined_footnote"]](
    label             = selected[["spec"]][["label"]],
    samples           = samples[[parameter]],
    posterior_defined = random[["defined_share"]](
      samples[[parameter]],
      selected    = selected,
      conditional = plan[["conditional"]]
    )
  )
  method <- .hypothesis_model_plan_execution_method(plans, backend)

  if (isTRUE(plan[["prior_draws"]])) {
    # Without a prior density with deterministic provenance, region
    # hypotheses take the prior probabilities from prior draws (the plans
    # refused point hypotheses).
    return(.hypothesis_model_random_prior_draws(
      provider                    = provider,
      parameter                 = parameter,
      selected                  = selected,
      posterior                 = samples[[parameter]],
      hypothesis                = hypothesis,
      transform_scaled          = plan[["transform_scaled"]],
      logBF                     = logBF,
      BF01                      = BF01,
      seed                      = seed,
      n_samples                 = plan[["n_samples"]],
      columns                   = columns,
      density_method            = method,
      conditional               = plan[["conditional"]]
    ))
  }

  arguments <- list(
    plan            = plan,
    provider        = provider,
    backend         = backend,
    logBF           = logBF,
    BF01            = BF01,
    seed            = seed,
    method          = method,
    density_control = density_control,
    columns         = columns
  )
  values <- vapply(.hypothesis_model_plan_targets(plans), `[[`, numeric(1), "value")
  out <- if (is.null(plan[["evaluation"]][["parameter"]])) {
    do.call(.hypothesis_model_plan_random_evaluate, c(arguments, list(
      posterior  = samples[[parameter]],
      parameter  = parameter,
      hypothesis = hypothesis,
      values     = values
    )))
  } else if (plan[["region"]]) {
    do.call(.hypothesis_model_plan_random_sd_region, c(arguments, list(
      posterior  = samples[[parameter]],
      hypothesis = hypothesis,
      values     = values
    )))
  } else {
    do.call(.hypothesis_model_plan_random_sd_points, c(arguments, list(
      hypothesis = hypothesis,
      values     = values
    )))
  }
  if (!is.null(defined_footnote)) {
    attr(out, "footnotes") <- c(attr(out, "footnotes"), defined_footnote)
  }

  out
}


# Statements on the marginal posterior 'posterior' of the random-effect
# quantity 'parameter'. For ordinate backends, the ordinates of the point
# values 'values' of that quantity are attached first.
.hypothesis_model_plan_random_evaluate <- function(plan, provider, backend,
                                                   posterior, parameter,
                                                   hypothesis, values, logBF,
                                                   BF01, seed, method,
                                                   density_control, columns) {

  if (identical(method, "precomputed")) {
    target    <- plan[["density_target"]]
    attached  <- .hypothesis_model_plan_attach_scalar(
      provider          = provider,
      backend           = backend,
      plan              = plan,
      posterior         = posterior,
      parameter         = target[["parameter"]],
      parameter_label   = plan[["label"]],
      values            = values,
      spec              = target[["parameter_spec"]],
      conditional       = NULL,
      density_control   = density_control,
      display_transform = target[["display_transform"]]
    )
    posterior <- attached[["posterior"]]
    request   <- attached[["request"]]
  }

  out <- hypothesis_BF.default(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_model_append_ordinate_warnings(
      table     = out,
      posterior = posterior,
      backend   = backend,
      request   = request,
      parameter = parameter
    )
  }

  out
}


# Point statements on a variance, evaluated through its standard deviation:
# the table shows the variance statements, and the prior and posterior
# densities on the variance scale.
.hypothesis_model_plan_random_sd_points <- function(plan, hypothesis, values,
                                                    backend, ...) {

  evaluation <- plan[["evaluation"]]
  out <- tryCatch(.hypothesis_model_plan_random_evaluate(
    plan       = plan,
    backend    = backend,
    posterior  = evaluation[["samples"]][[evaluation[["parameter"]]]],
    parameter  = evaluation[["parameter"]],
    hypothesis = .hypothesis_model_plan_sd_hypothesis(hypothesis, evaluation[["parameter"]]),
    values     = sqrt(values),
    ...
  ), error = function(condition) {
    # The backend restates the failure of a request for the evaluated values
    # (the standard deviations) in terms of the requested ones (the variances).
    if (is.function(backend[["rescale_failure"]])) {
      restated <- backend[["rescale_failure"]](
        condition,
        requested      = values,
        evaluated      = sqrt(values),
        label          = plan[["label"]],
        density_method = backend[["name"]]
      )
      if (!is.null(restated)) condition <- restated
    }
    stop(condition)
  })

  .hypothesis_model_plan_sd_restore(out, hypothesis)
}


# Statements comparing a variance point with a region of the variance, e.g.
# 'tau2 = 0.09 vs tau2 > 0.09' (one point value per statement, 'values').
# Their Bayes factor is the point part, the point against the encompassing
# model ('tau2 = 0.09 vs tau2 != 0.09', evaluated through the standard
# deviation as a point statement), over the region part, the region against
# the encompassing model (evaluated on the variance draws); the inverse with
# the region on the left. BayesTools evaluates the statements on the variance
# with a unit point part (a posterior ordinate equal to the prior ordinate,
# carrying the Monte Carlo error of the point part, which BayesTools combines
# with that of the region part). The real point part is composed only after
# detecting the admitted numerical-zero point/region pair, as BayesTools does.
.hypothesis_model_plan_random_sd_region <- function(plan, provider, backend,
                                                    posterior, hypothesis,
                                                    values, logBF, BF01, seed,
                                                    method, density_control,
                                                    columns) {

  if (length(values) != length(hypothesis[["statements"]])) {
    stop("Internal error: variance point-region statements are misaligned.",
         call. = FALSE)
  }
  parameter <- plan[["parameter"]]
  points    <- unique(values)
  point     <- .hypothesis_model_plan_random_sd_points(
    plan            = plan,
    hypothesis      = hypothesis_rewrite(
      hypothesis_parse(sprintf(
        "theta = %.17g vs theta != %.17g", points, points
      )),
      c(theta = parameter)
    ),
    values          = points,
    provider        = provider,
    backend         = backend,
    logBF           = FALSE,
    BF01            = FALSE,
    seed            = seed,
    method          = method,
    density_control = density_control,
    columns         = "default"
  )
  unit_log_ordinate <- vapply(points, function(value) {
    prior_density_ordinate(
      plan[["prior_density"]],
      value
    )[["log_density"]]
  }, numeric(1))
  # The unit part carries the conditioning of the draws, so that BayesTools
  # matches it to conditional (product-space) draws.
  posterior_metadata(posterior, "posterior_ordinate") <- do.call(
    posterior_log_ordinate_attribute,
    c(
      list(
        value          = points,
        log_ordinate   = unit_log_ordinate,
        method         = "unit point part",
        density_method = backend[["name"]],
        diagnostics    = list(BF_error_percent = as.numeric(point[["BF_error"]]))
      ),
      .hypothesis_model_condition_metadata(posterior)
    )
  )

  out <- hypothesis_BF.default(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = "precomputed"
  )
  point_index <- match(values, points)
  point_log_BF <- attr(point, "raw_log_BF", exact = TRUE)[point_index]
  point_left  <- vapply(hypothesis[["statements"]], function(statement) {
    identical(statement[["left"]][["type"]], "point")
  }, logical(1))
  unit_log_BF <- attr(out, "raw_log_BF", exact = TRUE)
  # Here the accepted real point sign certifies its numerical zero ordinate;
  # the unit part's directional infinity certifies zero posterior region mass.
  # Finite-log natural underflow and pre-existing NA results do not match.
  undefined <- !is.na(point_log_BF) & point_log_BF == -Inf &
    !is.na(unit_log_BF) & unit_log_BF == ifelse(point_left, Inf, -Inf)
  raw_log_BF <- rep(NA_real_, length(unit_log_BF))
  raw_log_BF[!undefined] <- unit_log_BF[!undefined] +
    ifelse(point_left[!undefined], point_log_BF[!undefined], -point_log_BF[!undefined])
  attr(out, "raw_log_BF") <- raw_log_BF
  attr(out, "raw_BF") <- exp(raw_log_BF)
  numerical_diagnostics <- attr(point, "numerical_diagnostics", exact = TRUE)[point_index]
  attr(out, "numerical_diagnostics") <- numerical_diagnostics
  if ("BF" %in% names(out)) {
    out[["BF"]] <- format_BF_from_log(
      raw_log_BF, logBF = logBF, BF01 = BF01,
      diagnostics = numerical_diagnostics)
  }

  if(any(undefined)){
    warnings <- attr(out, "warnings", exact = TRUE)
    rows <- rownames(out)[undefined]
    matched <- names(warnings) %in% rows
    warnings[matched] <- trimws(gsub(
      "Posterior region mass is zero; reported BF is boundary-valued.", "",
      warnings[matched], fixed = TRUE))
    warnings <- warnings[nzchar(warnings)]
    attr(out, "warnings") <- c(warnings, stats::setNames(rep(
      "The posterior point ordinate and region mass are both numerically zero; transitive BF is undefined.",
      length(rows)), rows))
  }

  # The warnings and density diagnostics of the point parts, on the rows of
  # their statements.
  point_warnings <- attr(point, "warnings", exact = TRUE)
  point_rows     <- match(names(point_warnings), rownames(point))
  if (length(point_rows) != length(point_warnings)) {
    point_rows <- rep(NA_integer_, length(point_warnings))
  }
  row_warnings <- unlist(lapply(seq_along(values), function(i) {
    warnings <- as.character(point_warnings[point_rows %in% point_index[[i]]])
    stats::setNames(warnings, rep(rownames(out)[[i]], length(warnings)))
  }), use.names = TRUE)
  attr(out, "warnings") <- .hypothesis_model_unique_named_warnings(c(
    attr(out, "warnings", exact = TRUE),
    row_warnings,
    point_warnings[is.na(point_rows)]
  ))
  attr(out, "density_diagnostics") <- attr(point, "density_diagnostics", exact = TRUE)

  out
}


# Point statements on a variance written on its standard deviation 'sd':
# each point side's value is replaced by its square root.
.hypothesis_model_plan_sd_hypothesis <- function(hypothesis, sd) {

  statements <- vapply(hypothesis[["statements"]], function(statement) {

    sides <- vapply(c("left", "right"), function(side_name) {
      side <- statement[[side_name]]
      operator <- switch(
        side[["type"]],
        point     = "=",
        not_point = "!=",
        stop("Internal error: expected a point statement on a variance.",
             call. = FALSE)
      )
      paste("theta", operator, sprintf("%.17g", sqrt(side[["value"]])))
    }, character(1))
    if (isTRUE(statement[["explicit"]])) {
      paste(sides[[1L]], "vs", sides[[2L]])
    } else {
      sides[[1L]]
    }
  }, character(1))

  hypothesis_rewrite(
    hypothesis_parse(unname(statements)),
    c(theta = sd)
  )
}


# The table of variance point statements evaluated through the standard
# deviation: the variance statements, and the prior and posterior densities
# on the variance scale (divided by the derivative 2 * sd of the square).
.hypothesis_model_plan_sd_restore <- function(out, hypothesis) {

  values <- vapply(hypothesis[["statements"]], function(statement) {
    statement[["left"]][["value"]]
  }, numeric(1))
  if (nrow(out) != length(values)) {
    stop("Internal error: variance point rows are misaligned.", call. = FALSE)
  }
  for (column in intersect(c("prior", "posterior"), names(out))) {
    out[[column]] <- out[[column]] / (2 * sqrt(values))
  }
  diagnostics <- attr(out, "numerical_diagnostics", exact = TRUE)
  if(length(diagnostics) == length(values)){
    for(i in seq_along(values)){
      jacobian <- -log(2 * sqrt(values[[i]]))
      diagnostics[[i]]$log_prior_height <- diagnostics[[i]]$log_prior_height + jacobian
      diagnostics[[i]]$log_posterior_height <- diagnostics[[i]]$log_posterior_height + jacobian
      diagnostics[[i]]$log_jacobian <- diagnostics[[i]]$log_jacobian + jacobian
      diagnostics[[i]]$requested_value <- values[[i]]
      diagnostics[[i]]$evaluation_scale <- "sd"
      diagnostics[[i]]$requested_scale <- "variance"
    }
    attr(out, "numerical_diagnostics") <- diagnostics
    if("BF" %in% names(out)) attr(out[["BF"]], "numerical_diagnostics") <- diagnostics
  }
  attr(out, "hypothesis_ast") <- hypothesis

  out
}

# Region hypotheses on a random-effect quantity without a canonical prior
# density: prior probabilities from the quantity's prior draws.
.hypothesis_model_random_prior_draws <- function(
    provider, parameter, selected, posterior, hypothesis,
    transform_scaled, logBF, BF01, seed, n_samples,
    columns, density_method, conditional) {

  random <- provider[["random"]]
  prior  <- random[["prior_draws"]](
    parameter        = parameter,
    selected         = selected,
    transform_scaled = transform_scaled,
    n_samples        = n_samples,
    seed             = seed
  )
  prior_values  <- as.numeric(prior[["samples"]][, 1L])
  prior_defined <- random[["defined_draws"]](
    prior_values,
    prior[["samples"]],
    selected[["spec"]][["label"]]
  )

  out <- hypothesis_BF.default(
    posterior      = as.numeric(posterior),
    prior          = prior_values[prior_defined],
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = density_method
  )
  footnote <- random[["defined_footnote"]](
    label             = selected[["spec"]][["label"]],
    samples           = prior[["samples"]],
    posterior_defined = random[["defined_share"]](
      posterior, selected = selected, conditional = conditional
    ),
    prior_defined     = prior_defined
  )
  if (!is.null(footnote)) {
    attr(out, "footnotes") <- c(attr(out, "footnotes"), footnote)
  }

  out
}
# ---------------------------------------------------------------------------- #
# Random-effect quantities
# ---------------------------------------------------------------------------- #

# Random-effect quantities are tested on their BayesTools mixed posterior
# (catalog support, declared inclusion- and allocation-gate atoms, exact
# prior density). Point statements on a variance, and the point parts of
# statements comparing a variance point with a region, are evaluated through
# its standard deviation with the square display transform, so that a
# variance and its standard deviation give the same Bayes factor with every
# method.
.hypothesis_model_plan_random <- function(plan, provider, cache) {

  parameter <- plan[["parameter"]]
  plan[["kind"]]  <- "random"
  plan[["route"]] <- "random"
  plan[["group"]] <- paste("random", parameter, sep = "\r")
  selected <- .hypothesis_model_plan_cached(
    cache,
    paste("random_select", parameter, plan[["transform_scaled"]], sep = "\r"),
    function() {
      provider[["random"]][["select"]](
        parameter        = parameter,
        transform_scaled = plan[["transform_scaled"]]
      )
    }
  )
  label <- selected[["spec"]][["label"]]
  plan[["random"]]  <- selected
  plan[["support"]] <- provider[["random"]][["support"]](selected)
  if (identical(selected[["spec"]][["status"]], "structural")) {
    plan[["refusal"]] <- .hypothesis_model_plan_fixed_refusal(label)
    return(plan)
  }
  samples <- .hypothesis_model_plan_cached(
    cache,
    paste("random_samples", parameter, plan[["transform_scaled"]],
          plan[["conditional"]], sep = "\r"),
    function() {
      provider[["random"]][["mixed_posterior"]](
        parameter        = parameter,
        transform_scaled = plan[["transform_scaled"]],
        conditional               = plan[["conditional"]],
        selected                  = selected
      )
    }
  )
  prior_density <- posterior_metadata(samples[[parameter]], "prior_density")
  plan[["draws"]]         <- list(samples = samples, posterior = samples[[parameter]])
  plan[["prior_density"]] <- prior_density
  # Region probabilities need a prior density with deterministic provenance;
  # without one (no density, or a plotting grid such as the SD components of
  # nested allocations) they come from prior draws.
  plan[["prior_draws"]] <- !.hypothesis_model_plan_density_has_provenance(prior_density)
  if (plan[["prior_draws"]] && plan[["conditional"]]) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Conditional hypotheses are unavailable for random-effect quantity '",
        label, "' because its prior density is unavailable."
      ),
      "target"
    )
    return(plan)
  }
  if (!plan[["point"]]) {
    return(plan)
  }
  refs <- plan[["refs"]]
  if (any(!refs[["direct"]])) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Point-null tests for random-effect quantities require a direct ",
        "scalar parameter reference."
      ),
      "statement"
    )
    return(plan)
  }
  if (is.null(prior_density)) {
    plan[["refusal"]] <- .hypothesis_model_refusal(
      paste0(
        "Point-null Bayes factors are unavailable for random-effect quantity '",
        label, "' because its prior density is unavailable. Use a region or ",
        "directional hypothesis."
      ),
      "target"
    )
    return(plan)
  }

  evaluation <- .hypothesis_model_plan_random_evaluation(plan, provider, cache, selected)
  if (!is.null(evaluation[["refusal"]])) {
    plan[["refusal"]] <- evaluation[["refusal"]]
    return(plan)
  }
  plan[["evaluation"]] <- evaluation
  if (!is.null(evaluation[["parameter"]])) {
    plan[["group"]] <- paste(
      "random", parameter, if (plan[["region"]]) "sd_region" else "sd",
      sep = "\r"
    )
  }
  zero_boundary <- if (any(refs[["value"]] == 0)) {
    provider[["random"]][["zero_boundary"]](selected)
  }
  zero_alternative <- zero_boundary[["alternative"]]
  point_mass_reason <- paste0(
    "Point-null Bayes factors at %s are unavailable for random-effect ",
    "quantity '%s' because its prior has a point mass at %s. ",
    gsub("%", "%%", provider[["random_point_mass_advice"]], fixed = TRUE)
  )
  plan[["targets"]] <- lapply(seq_len(nrow(refs)), function(i) {
    value  <- refs[["value"]][[i]]
    label_value <- paste0(label, " = ", format(value, digits = 15L, trim = TRUE))
    # A prior point mass at 0 is refused as such; otherwise 0 is a
    # nonregular product boundary of an allocation-derived SD.
    point_mass <- !is.null(zero_boundary) && value == 0 && identical(
      prior_ordinate_status(prior_density, value, labels = label_value)[["condition"]],
      "BayesTools_point_mass_at_null"
    )
    target <- .hypothesis_model_plan_target(
      value         = value,
      label         = label_value,
      prior_density = prior_density,
      atoms         = posterior_metadata(samples[[parameter]], "atoms"),
      refusal       = if (!is.null(zero_boundary) && value == 0 && !point_mass) {
        c(
          .hypothesis_model_refusal(
            paste0(
              "Point-null Bayes factors are unavailable for allocation-derived ",
              "random-effect quantity '", label, "' at 0 because zero is a ",
              "nonregular product boundary of the common scale and allocation ",
              "weight. ",
              if (is.null(zero_alternative)) {
                "Use a region or directional hypothesis."
              } else {
                paste0("Test '", zero_alternative,
                  "' to compare omission of this component.")
              }
            ),
            "target"
          ),
          list(value_specific = TRUE)
        )
      }
    )
    value_text <- format(value, digits = 15L, trim = TRUE)
    target[["point_mass_reason"]] <- sprintf(point_mass_reason, value_text, label, value_text)
    target
  })
  precomputed <- .hypothesis_model_plan_random_precomputed(
    plan       = plan,
    provider     = provider,
    evaluation = evaluation,
    label      = label
  )
  plan[["method_refusals"]] <- precomputed[["refusals"]]
  plan[["density_target"]]  <- precomputed[["density_target"]]

  plan
}


# How point statements on a random-effect quantity are evaluated: variances
# with a standard-deviation counterpart in the catalog are evaluated on the
# standard deviation ('parameter', values mapped by the square root); other
# quantities on themselves (NULL 'parameter'). Of a statement comparing a
# variance point with a region, only the point part is evaluated on the
# standard deviation; the region part is evaluated on the variance draws.
.hypothesis_model_plan_random_evaluation <- function(plan, provider, cache, selected) {

  statement <- plan[["hypothesis"]][["statements"]][[1L]]
  sides <- vapply(list(statement[["left"]], statement[["right"]]),
    `[[`, character(1L), "type")
  values <- unique(plan[["refs"]][["value"]])
  unsupported <- length(values) > 1L ||
    ("not_point" %in% sides && "region" %in% sides) ||
    all(sides == "point")
  if (unsupported) {
    return(list(refusal = .hypothesis_model_refusal(
      "Unsupported hypothesis comparison. Compare one point value with its complement or a region hypothesis.",
      "statement"
    )))
  }
  pair <- provider[["random"]][["sd_pair"]](selected)
  # At a variance of 0 the square map is singular: the variance is evaluated
  # on its own draws there (its prior ordinate is classified on its own
  # scale).
  if (is.null(pair) || any(plan[["refs"]][["value"]] <= 0)) {
    return(list(parameter = NULL))
  }
  pair_samples <- .hypothesis_model_plan_cached(
    cache,
    paste("random_samples", pair, plan[["transform_scaled"]],
          plan[["conditional"]], sep = "\r"),
    function() {
      provider[["random"]][["mixed_posterior"]](
        parameter        = pair,
        transform_scaled = plan[["transform_scaled"]],
        conditional      = plan[["conditional"]],
        selected         = NULL
      )
    }
  )

  list(
    parameter = pair,
    samples   = pair_samples
  )
}


# Refusals of the ordinate backends for a random-effect point statement:
# conditional statements, quantities without a supported scalar
# random-component coordinate, and values where the public transformation of
# that coordinate is singular. The first two are generic causes
# ("conditional_random", "random_target") for which a backend may declare
# additional condition classes.
.hypothesis_model_plan_random_precomputed <- function(plan, provider, evaluation,
                                                      label) {

  refused <- function(reason, type = "method") {
    list(refusals = .hypothesis_model_plan_precomputed_refusals(reason, type))
  }
  refused_cause <- function(reason, cause, type = "method") {
    list(refusals = list(
      ordinate = .hypothesis_model_cause_refusal(reason, cause, type)
    ))
  }
  if (plan[["conditional"]]) {
    return(refused_cause(
      "Conditional random-effect hypotheses support 'density_method = \"KDE\"' only.",
      "conditional_random"
    ))
  }
  parameter <- if (is.null(evaluation[["parameter"]])) {
    plan[["parameter"]]
  } else {
    evaluation[["parameter"]]
  }
  target <- provider[["random"]][["density_target"]](
    parameter,
    transform_scaled = plan[["transform_scaled"]],
    operation        = "point hypotheses"
  )
  if (is.null(target[["parameter"]])) {
    # The refusal names the tested quantity (a variance rather than the
    # standard deviation it is evaluated through).
    if (!identical(parameter, plan[["parameter"]])) {
      own <- provider[["random"]][["density_target"]](
        plan[["parameter"]],
        transform_scaled = plan[["transform_scaled"]],
        operation        = "point hypotheses"
      )
      if (!is.null(own[["reason"]])) {
        target[["reason"]] <- own[["reason"]]
      }
    }
    # The refusal of the provider keeps its type (default "target": no
    # supported target) and its cause (default "random_target"), for which the
    # backend declares its additional classes.
    return(refused_cause(
      target[["reason"]],
      cause = .hypothesis_model_random_refusal_field(
        target, "cause", "random_target", c("original_scale", "random_target")
      ),
      type  = .hypothesis_model_random_refusal_field(
        target, "type", "target", c("target", "method", "fixed")
      )
    ))
  }
  # Values outside the support are refused by their prior ordinates.
  values <- .hypothesis_model_plan_random_source_values(plan, evaluation)
  values <- values[is.finite(values)]
  display_transform <- target[["display_transform"]]
  if (!is.null(display_transform) && length(values) > 0L) {
    source_values <- parameter_transform_inverse(values, display_transform)
    jacobian      <- parameter_transform_jacobian(
      source_values,
      display_transform
    )
    singular <- !is.finite(source_values) | !is.finite(jacobian) | jacobian <= 0
    if (any(singular)) {
      return(refused(paste0(
        "Point-null Bayes factors are unavailable for random-effect ",
        "quantity '", label, "' at ", values[which(singular)[[1L]]],
        " because its public transformation is singular at that support ",
        "boundary. Use the corresponding directly modeled scale or a ",
        "region hypothesis."
      ), "target"))
    }
  }

  list(refusals = list(), density_target = target)
}


# The 'field' of a refusal of random$density_target(): 'default' when absent,
# otherwise one of 'allowed'.
.hypothesis_model_random_refusal_field <- function(target, field, default,
                                                   allowed) {

  value <- target[[field]]
  if (is.null(value)) {
    return(default)
  }
  if (!is.character(value) || length(value) != 1L || !value %in% allowed) {
    stop(
      "Invalid hypothesis provider: the '", field, "' of a refusal of ",
      "'random$density_target()' must be one of ",
      paste0("'", allowed, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }

  value
}


# The tested values on the evaluated quantity's scale: the square roots of
# the values of a variance evaluated through its standard deviation.
.hypothesis_model_plan_random_source_values <- function(plan, evaluation) {

  values <- plan[["refs"]][["value"]]
  if (is.null(evaluation[["parameter"]])) {
    return(values)
  }

  ifelse(values >= 0, sqrt(pmax(values, 0)), NA_real_)
}
