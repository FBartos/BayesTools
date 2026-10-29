# ============================================================================ #
# hypothesis-model.R
# ============================================================================ #
#
# The hypothesis engine: hypothesis_BF_model() parses and plans the statements
# against the provider of a model, checks every plan against the requested
# density backend, evaluates the plans group by group, and binds the results.
# See .agents/instructions/hypothesis.md.
#
# ============================================================================ #


# The statements of the plans of one group as one hypothesis AST: the
# statements as written ('field' "statement") or rewritten onto the tested
# parameter ("hypothesis").
.hypothesis_model_plan_group_ast <- function(plans, field) {

  ast <- plans[[1L]][[field]]
  ast[["statements"]] <- unlist(lapply(plans, function(plan) {
    plan[[field]][["statements"]]
  }), recursive = FALSE)

  ast
}


# Evaluates the statements of one group of plans, which share one route and
# one tested target, and restores the labels the statements were written in.
.hypothesis_model_plan_execute <- function(plans, provider, backend,
                                           conditional_omitted, logBF, BF01,
                                           seed, density_control, columns) {

  plan       <- plans[[1L]]
  display    <- .hypothesis_model_plan_group_ast(plans, "statement")
  hypothesis <- .hypothesis_model_plan_group_ast(plans, "hypothesis")
  if (conditional_omitted) {
    .hypothesis_model_warn_full_ensemble(provider, plan[["parameter"]])
  }
  arguments <- list(
    plans           = plans,
    hypothesis      = hypothesis,
    provider        = provider,
    backend         = backend,
    logBF           = logBF,
    BF01            = BF01,
    seed            = seed,
    density_control = density_control,
    columns         = columns
  )
  out <- switch(
    plan[["route"]],
    scalar      = do.call(.hypothesis_model_plan_execute_scalar, arguments),
    levels      = do.call(.hypothesis_model_plan_execute_levels, arguments),
    combination = do.call(.hypothesis_model_plan_execute_combination, arguments),
    random      = do.call(.hypothesis_model_plan_execute_random, arguments),
    stop("Internal error: unknown hypothesis route.", call. = FALSE)
  )

  .hypothesis_model_restore_hypothesis_labels(
    out             = out,
    hypothesis      = display,
    parameter_label = plan[["label"]],
    producer_parameter = plan[["parameter"]],
    route           = plan[["route"]]
  )
}


# A coefficient test on the full ensemble of a prior with a null and an
# alternative component includes the models in which the coefficient is fixed
# by its null component (the provider decides which parameters those are).
.hypothesis_model_warn_full_ensemble <- function(provider, parameter) {

  if (!is.function(provider[["full_ensemble"]]) ||
      !isTRUE(provider[["full_ensemble"]](parameter))) {
    return(invisible(FALSE))
  }
  warning(
    "Model-averaged parameter test: this hypothesis uses the full ",
    "ensemble for ", parameter, ", including models where ",
    parameter, " is fixed by its null component. To test the ",
    "hypothesis only within models where ", parameter,
    " is active, use 'conditional = TRUE'.",
    call. = FALSE
  )

  invisible(TRUE)
}


# The method a group is evaluated with: ordinate backends only for point
# statements (region statements are method-free and use the draws).
.hypothesis_model_plan_execution_method <- function(plans, backend) {

  point <- any(vapply(plans, `[[`, logical(1), "point"))
  if (identical(backend[["type"]], "ordinate")) {
    return(if (point) "precomputed" else "KDE")
  }

  backend[["method"]]
}


.hypothesis_model_plan_targets <- function(plans) {

  unlist(lapply(plans, `[[`, "targets"), recursive = FALSE)
}


# The object that backends compute their ordinates from.
.hypothesis_model_capability_object <- function(provider) {

  if (!is.null(provider[["capability_object"]])) {
    provider[["capability_object"]]
  } else {
    provider[["object"]]
  }
}


# The request of a density backend for the point values 'values' of the
# posterior of a plan (see hypothesis_density_backend()). 'sample' are the
# draws that the ordinates belong to (the vector, or the draws of the level
# of a request for a factor level): their conditioning is the 'condition' of
# the request.
.hypothesis_model_request <- function(provider, backend, plan, kind, parameter,
                                      label, values, sample, spec = NULL,
                                      conditional = NULL, level = NA_character_,
                                      display_transform = NULL) {

  list(
    kind              = kind,
    route             = plan[["kind"]],
    object            = .hypothesis_model_capability_object(provider),
    model             = provider[["object"]],
    parameter         = parameter,
    label             = label,
    level             = level,
    values            = values,
    conditional       = conditional,
    condition         = .hypothesis_model_condition_metadata(sample),
    spec              = spec,
    display_transform = display_transform,
    density_method    = backend[["name"]],
    inference_type    = plan[["inference_type"]]
  )
}


# Attaches the ordinates of scalar point values to a posterior with the
# ordinate backend. Returns the posterior with the ordinates and the request
# that the backend received (which its warnings() and diagnostics() receive
# as well).
.hypothesis_model_plan_attach_scalar <- function(provider, backend, plan,
                                                 posterior, parameter,
                                                 parameter_label, values, spec,
                                                 conditional, density_control,
                                                 display_transform = NULL) {

  posterior <- .hypothesis_model_keep_requested_ordinates(
    posterior  = posterior,
    point_refs = data.frame(level = NA_character_, value = values)
  )
  request <- .hypothesis_model_request(
    provider, backend, plan, kind = "scalar", parameter = parameter,
    label = parameter_label, values = unique(values), sample = posterior,
    spec = spec, conditional = conditional, display_transform = display_transform
  )
  posterior <- .hypothesis_model_attach_request(
    posterior = posterior,
    backend   = backend,
    request   = request,
    control   = density_control,
    cache     = new.env(parent = emptyenv())
  )

  list(posterior = posterior, request = request)
}


.hypothesis_model_plan_execute_scalar <- function(plans, hypothesis, provider,
                                                  backend, logBF, BF01, seed,
                                                  density_control, columns) {

  plan      <- plans[[1L]]
  posterior <- plan[["draws"]][["posterior"]]
  method    <- .hypothesis_model_plan_execution_method(plans, backend)
  if (identical(method, "precomputed")) {
    targets  <- .hypothesis_model_plan_targets(plans)
    attached <- .hypothesis_model_plan_attach_scalar(
      provider        = provider,
      backend         = backend,
      plan            = plan,
      posterior       = posterior,
      parameter       = plan[["parameter"]],
      parameter_label = plan[["label"]],
      values          = vapply(targets, `[[`, numeric(1), "value"),
      spec            = targets[[1L]][["spec"]],
      conditional     = if (plan[["conditional"]]) plan[["parameter"]],
      density_control = density_control
    )
    posterior <- attached[["posterior"]]
    request   <- attached[["request"]]
  }

  out <- hypothesis_BF.default(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = plan[["parameter"]],
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
      request   = request
    )
  }

  out
}


# Statements on the marginal posterior of a factor term: each point value of
# a level gets the level's planned prior density and, for ordinate backends,
# its ordinate of the level's linear target.
.hypothesis_model_plan_execute_levels <- function(plans, hypothesis, provider,
                                                  backend, logBF, BF01, seed,
                                                  density_control, columns) {

  plan      <- plans[[1L]]
  posterior <- plan[["draws"]][["posterior"]]
  targets   <- .hypothesis_model_plan_targets(plans)
  for (target in targets) {
    posterior_metadata(posterior[[target[["level"]]]], "prior_density") <-
      target[["prior_density"]]
  }
  method <- .hypothesis_model_plan_execution_method(plans, backend)
  if (identical(method, "precomputed")) {
    refs <- unique(data.frame(
      level = vapply(targets, function(target) target[["level"]], character(1)),
      value = vapply(targets, `[[`, numeric(1), "value"),
      stringsAsFactors = FALSE
    ))
    posterior <- .hypothesis_model_keep_requested_ordinates(posterior, refs)
    cache     <- new.env(parent = emptyenv())
    requests  <- list()
    for (i in seq_len(nrow(refs))) {
      target <- targets[[which(
        vapply(targets, function(target) target[["level"]], character(1)) ==
          refs[["level"]][[i]]
      )[[1L]]]]
      request <- .hypothesis_model_request(
        provider, backend, plan, kind = "level",
        parameter = plan[["parameter"]], label = plan[["parameter"]],
        values = refs[["value"]][[i]],
        sample = posterior[[refs[["level"]][[i]]]], spec = target[["spec"]],
        conditional = if (plan[["conditional"]]) plan[["parameter"]],
        level = refs[["level"]][[i]]
      )
      posterior <- .hypothesis_model_attach_request(
        posterior = posterior,
        backend   = backend,
        request   = request,
        control   = density_control,
        cache     = cache
      )
      requests[[i]] <- request
    }
    # The requests of a group differ in their level, values, and spec only.
    request <- requests[[1L]]
  }

  out <- hypothesis_BF.default(
    posterior      = posterior,
    hypothesis     = hypothesis,
    parameter      = plan[["parameter"]],
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
      request   = request
    )
  }

  out
}


# Statements on one linear combination of the levels of a factor term,
# evaluated on the scalar linear target (level contrasts, derived linear
# expressions).
.hypothesis_model_plan_execute_combination <- function(plans, hypothesis,
                                                       provider, backend,
                                                       logBF, BF01, seed,
                                                       density_control,
                                                       columns) {

  plan   <- plans[[1L]]
  target <- if (length(plans) == 1L) {
    plan[["linear_target"]]
  } else {
    combined <- hypothesis_linear_target(
      posterior  = plan[["draws"]][["posterior"]],
      hypothesis = hypothesis,
      parameter  = plan[["parameter"]]
    )
    posterior_metadata(combined[["posterior"]], "prior_density") <-
      posterior_metadata(plan[["linear_target"]][["posterior"]], "prior_density")
    combined
  }
  method <- .hypothesis_model_plan_execution_method(plans, backend)
  if (identical(method, "precomputed")) {
    targets  <- .hypothesis_model_plan_targets(plans)
    attached <- .hypothesis_model_plan_attach_scalar(
      provider        = provider,
      backend         = backend,
      plan            = plan,
      posterior       = target[["posterior"]],
      parameter       = target[["parameter"]],
      parameter_label = plan[["label"]],
      values          = vapply(targets, `[[`, numeric(1), "value"),
      spec            = targets[[1L]][["spec"]],
      conditional     = .hypothesis_model_first_nonempty_condition(
        posterior_metadata(target[["posterior"]], "condition"),
        c("effective_conditional", "conditional")
      ),
      density_control = density_control
    )
    target[["posterior"]] <- attached[["posterior"]]
    request               <- attached[["request"]]
  }

  out <- hypothesis_BF.default(
    posterior      = target[["posterior"]],
    hypothesis     = target[["hypothesis"]],
    parameter      = target[["parameter"]],
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (identical(method, "precomputed")) {
    out <- .hypothesis_model_append_ordinate_warnings(
      table     = out,
      posterior = target[["posterior"]],
      backend   = backend,
      request   = request
    )
  }
  # One linear target: row i is statement i.
  .hypothesis_model_set_row_names(
    out       = out,
    row_names = .hypothesis_model_row_names(
      labels     = rep(plan[["parameter"]], nrow(out)),
      statements = seq_len(nrow(out))
    )
  )
}


# The first non-empty element of the 'condition' draw metadata among 'names'.
.hypothesis_model_first_nonempty_condition <- function(condition, names) {

  for (name in names) {
    value <- condition[[name]]
    if (!is.null(value) && length(value) > 0L) {
      return(value)
    }
  }

  return(NULL)
}


# The conditioning of the draws ('conditional', 'conditional_rule',
# 'condition_key', and the condition events) as the metadata of an ordinate
# that is to be matched to them.
.hypothesis_model_condition_metadata <- function(samples) {

  if (is.null(samples)) {
    return(list())
  }
  condition   <- posterior_metadata(samples, "condition")
  conditional <- .hypothesis_model_first_nonempty_condition(
    condition,
    c("effective_conditional", "conditional")
  )
  conditional_rule <- .hypothesis_model_first_nonempty_condition(
    condition,
    c("effective_conditional_rule", "conditional_rule")
  )
  condition_key <- condition[["condition_key"]]

  declared_empty <- isTRUE(condition[["averaged"]]) &&
    any(vapply(c("effective_conditional", "conditional"), function(field) {
      !is.null(condition[[field]]) && length(condition[[field]]) == 0L
    }, logical(1)))
  if (declared_empty && is.null(conditional)) conditional <- character()

  metadata <- list(
    conditional              = conditional,
    conditional_rule         = conditional_rule,
    condition_key            = condition_key,
    condition_event          = condition[["condition_event"]],
    resolved_condition_event = condition[["resolved_condition_event"]]
  )
  metadata[!vapply(metadata, is.null, logical(1))]
}


.hypothesis_model_result_statements <- function(result) {

  hypothesis <- attr(result, "hypothesis_ast", exact = TRUE)
  if (!inherits(hypothesis, "BayesTools_hypothesis_ast") ||
      !is.list(hypothesis[["statements"]]) ||
      length(hypothesis[["statements"]]) == 0L) {
    return(NULL)
  }

  hypothesis[["statements"]]
}


.hypothesis_model_bind_parameter_results <- function(
    results, groups, hypothesis, backend = NULL) {

  if (length(results) == 0L) {
    stop("Internal error: no hypothesis parameter-group results were produced.",
         call. = FALSE)
  }
  n_statements <- length(hypothesis[["statements"]])
  group_sizes  <- vapply(groups, length, integer(1))
  result_sizes <- vapply(results, nrow, integer(1))
  valid_results <- vapply(
    results,
    inherits,
    logical(1),
    what = "BayesTools_hypothesis_BF"
  )
  same_columns <- vapply(
    results,
    function(result) identical(names(result), names(results[[1L]])),
    logical(1)
  )
  result_statement_sizes <- vapply(results, function(result) {

    statements <- .hypothesis_model_result_statements(result)
    if (is.null(statements)) -1L else length(statements)
  }, integer(1))
  quantities_per_statement <- result_sizes / group_sizes
  if (!all(valid_results) ||
      any(group_sizes <= 0L) || any(!is.finite(quantities_per_statement)) ||
      any(quantities_per_statement < 1L | quantities_per_statement != floor(quantities_per_statement)) ||
      !identical(group_sizes, result_statement_sizes) ||
      sum(group_sizes) != n_statements || !all(same_columns)) {
    stop("Internal error: hypothesis parameter-group results are misaligned.",
         call. = FALSE)
  }

  out <- do.call(rbind, results)
  concatenated_rows <- unlist(Map(function(group, count) rep(group, each = count),
    groups, quantities_per_statement), use.names = FALSE)
  n_rows <- sum(result_sizes)
  restore_order     <- order(concatenated_rows)
  # Row names as in hypothesis_BF(): rows of several statements on
  # one quantity carry their statement number in the full hypothesis, not the
  # number within their parameter group.
  row_names <- .hypothesis_model_row_names(
    labels     = sub(
      " \\([0-9]+\\)$", "",
      unlist(lapply(results, rownames), use.names = FALSE)
    ),
    statements = concatenated_rows
  )

  raw_BF <- unlist(lapply(results, function(result) {
    attr(result, "raw_BF", exact = TRUE)
  }), use.names = FALSE)
  if (length(raw_BF) != n_rows) {
    stop("Internal error: hypothesis parameter-group metadata are misaligned.",
         call. = FALSE)
  }

  row_log_BF <- unlist(lapply(results, function(result){
    logs <- attr(result, "raw_log_BF", exact = TRUE)
    if(is.null(logs)) log(attr(result, "raw_BF", exact = TRUE)) else logs
  }), use.names = FALSE)
  row_diagnostics <- lapply(c("numerical_diagnostics", "prior_numerical_diagnostics"), function(attribute){
    unlist(lapply(results, function(result){
      diagnostics <- attr(result, attribute, exact = TRUE)
      if(length(diagnostics) == nrow(result)) diagnostics else rep(list(NULL), nrow(result))
    }), recursive = FALSE)
  })
  warnings <- list()
  row_offset <- 0L
  for (i in seq_along(results)) {
    result <- results[[i]]
    rows   <- row_offset + seq_len(nrow(result))
    result_warnings <- attr(result, "warnings", exact = TRUE)
    if (length(result_warnings) > 0L && !is.null(names(result_warnings))) {
      warning_rows <- match(names(result_warnings), rownames(result))
      matched      <- !is.na(warning_rows)
      names(result_warnings)[matched] <- row_names[rows[warning_rows[matched]]]
    }
    warnings[[i]] <- result_warnings
    row_offset    <- row_offset + nrow(result)
  }
  rownames(out) <- row_names

  bound_operator <- unlist(lapply(results, function(result) {
    operator <- attr(result[["BF"]], "bound_operator", exact = TRUE)
    if (is.null(operator)) {
      return(rep(NA_character_, nrow(result)))
    }
    rep_len(as.character(operator), nrow(result))
  }), use.names = FALSE)
  source_bound_operator <- if(isTRUE(attr(results[[1L]]$BF, "BF01"))){
    invert_BF_bound_operator(bound_operator)
  }else bound_operator

  out <- out[restore_order, , drop = FALSE]
  attr(out, "raw_BF")         <- raw_BF[restore_order]
  attr(out, "raw_log_BF") <- row_log_BF[restore_order]
  attr(out, "numerical_diagnostics") <- row_diagnostics[[1L]][restore_order]
  attr(out, "prior_numerical_diagnostics") <- row_diagnostics[[2L]][restore_order]
  out[["BF"]] <- format_BF_from_log(
    row_log_BF[restore_order], logBF = attr(results[[1L]]$BF, "logBF"),
    BF01 = attr(results[[1L]]$BF, "BF01"), bound_operator = source_bound_operator[restore_order],
    diagnostics = row_diagnostics[[1L]][restore_order])
  if ("BF_error" %in% names(out)) {
    attr(out[["BF_error"]], "name") <-
      attr(results[[1L]][["BF_error"]], "name", exact = TRUE)
  }
  attr(out, "hypothesis_ast") <- hypothesis
  warnings <- .hypothesis_model_unique_named_warnings(
    unlist(warnings, use.names = TRUE)
  )
  if (length(warnings) > 0L && !is.null(names(warnings))) {
    warning_rows <- match(names(warnings), rownames(out))
    warnings <- c(
      warnings[is.na(warning_rows)],
      warnings[order(warning_rows, na.last = NA)]
    )
  }
  attr(out, "warnings") <- warnings
  attr(out[["BF"]], "bound_operator") <- bound_operator[restore_order]

  footnotes <- unique(unlist(lapply(results, function(result) {
    attr(result, "footnotes", exact = TRUE)
  }), use.names = FALSE))
  attr(out, "footnotes") <- if (length(footnotes) > 0L) footnotes else NULL

  diagnostics <- lapply(results, function(result) {
    attr(result, "density_diagnostics", exact = TRUE)
  })
  diagnostics <- diagnostics[vapply(
    diagnostics,
    function(diagnostic) is.data.frame(diagnostic) && nrow(diagnostic) > 0L,
    logical(1)
  )]
  if (length(diagnostics) > 0L) {
    attr(out, "density_diagnostics") <- if (is.function(backend[["bind_diagnostics"]])) {
      backend[["bind_diagnostics"]](diagnostics)
    } else {
      do.call(rbind, diagnostics)
    }
  } else {
    attr(out, "density_diagnostics") <- NULL
  }
  attr(out, "rownames") <- FALSE

  return(out)
}

.hypothesis_model_restore_hypothesis_labels <- function(
    out, hypothesis, parameter_label = NULL, producer_parameter = NULL,
    route = "scalar") {

  if (!is.data.frame(out)) {
    return(out)
  }

  statements <- hypothesis[["statements"]]
  result_label <- function(statement, side_name) {

    implicit_equality <- !isTRUE(statement[["explicit"]]) &&
      identical(statement[["left"]][["type"]], "point")
    if (implicit_equality) {
      side_name <- switch(side_name, left = "right", right = "left")
    }

    statement[[side_name]][["label"]]
  }
  transformed_statements <- .hypothesis_model_result_statements(out)
  if (!is.list(transformed_statements) ||
      length(transformed_statements) != length(statements)) {
    stop("Internal error: transformed hypothesis metadata are misaligned.",
         call. = FALSE)
  }

  display_left       <- vapply(
    statements,
    result_label,
    side_name = "left",
    character(1)
  )
  display_right      <- vapply(
    statements,
    result_label,
    side_name = "right",
    character(1)
  )
  quantities_per_statement <- nrow(out) / length(statements)
  if (!is.finite(quantities_per_statement) || quantities_per_statement < 1L ||
      quantities_per_statement != floor(quantities_per_statement)) {
    stop("Internal error: transformed hypothesis result rows are misaligned.",
         call. = FALSE)
  }
  # BayesTools expands one common quantity list in statement-major order.
  # Equal statement text still denotes separate original statement positions.
  statement_i <- rep(seq_along(statements), each = quantities_per_statement)
  if ("Alternative" %in% names(out)) {
    out[["Alternative"]] <- display_left[statement_i]
  }
  if ("Null" %in% names(out)) {
    out[["Null"]] <- display_right[statement_i]
  }

  attr(out, "hypothesis_ast") <- hypothesis

  if (!is.null(parameter_label)) {
    old_rownames  <- rownames(out)
    suffix        <- rep("", length(old_rownames))
    if (identical(route, "levels")) {
      check_char(producer_parameter, "producer_parameter", check_length = 1L, allow_NA = FALSE)
      producer_rows <- sub(" \\([0-9]+\\)$", "", old_rownames)
      has_suffix <- startsWith(producer_rows, paste0(producer_parameter, "["))
      suffix[has_suffix] <- substring(producer_rows[has_suffix], nchar(producer_parameter) + 1L)
    }
    out <- .hypothesis_model_set_row_names(
      out       = out,
      row_names = .hypothesis_model_row_names(
        labels     = paste0(parameter_label, suffix),
        statements = statement_i
      )
    )
  }
  attr(out, "rownames") <- FALSE

  return(out)
}


# Row names of hypothesis tables, as in hypothesis_BF(): rows of
# several statements on one quantity carry the statement number ("mu (1)",
# "mu (2)"). "mu1" reads like another parameter, and "mu.1" collides with the
# names `[.data.frame` gives duplicated rows.
.hypothesis_model_row_names <- function(labels, statements) {

  repeated <- labels %in% labels[duplicated(labels)]
  labels[repeated] <- paste0(labels[repeated], " (", statements[repeated], ")")

  return(make.unique(labels, sep = " "))
}


# Rename table rows together with the table warnings keyed by them.
.hypothesis_model_set_row_names <- function(out, row_names) {

  old_rownames  <- rownames(out)
  rownames(out) <- row_names

  warnings <- attr(out, "warnings", exact = TRUE)
  if (!is.null(warnings) && !is.null(names(warnings))) {
    warning_rows <- match(names(warnings), old_rownames)
    matched      <- !is.na(warning_rows)
    names(warnings)[matched] <- row_names[warning_rows[matched]]
    attr(out, "warnings") <- warnings
  }

  return(out)
}
# Evaluates one group of marginal-means plans: statements on the marginal
# means of one parameter, or on one linear combination of its levels.
.hypothesis_model_plan_execute_marginal_means <- function(plans, provider,
                                                          backend, logBF, BF01,
                                                          seed, density_control,
                                                          columns) {

  plan        <- plans[[1L]]
  parameter   <- plan[["parameter"]]
  hypothesis  <- .hypothesis_model_plan_group_ast(plans, "hypothesis")
  method      <- .hypothesis_model_plan_execution_method(plans, backend)
  precomputed <- identical(method, "precomputed")

  if (identical(plan[["route"]], "combination")) {
    target <- if (length(plans) == 1L) {
      plan[["linear_target"]]
    } else {
      hypothesis_linear_target(
        posterior  = plan[["draws"]][["posterior"]],
        hypothesis = hypothesis,
        parameter  = parameter
      )
    }
    if (precomputed) {
      targets  <- .hypothesis_model_plan_targets(plans)
      attached <- .hypothesis_model_plan_attach_scalar(
        provider        = provider,
        backend         = backend,
        plan            = plan,
        posterior       = target[["posterior"]],
        parameter       = target[["parameter"]],
        parameter_label = parameter,
        values          = vapply(targets, `[[`, numeric(1), "value"),
        spec            = targets[[1L]][["spec"]],
        conditional     = .hypothesis_model_first_nonempty_condition(
          posterior_metadata(target[["posterior"]], "condition"),
          c("effective_conditional", "conditional")
        ),
        density_control = density_control
      )
      target[["posterior"]] <- attached[["posterior"]]
      request               <- attached[["request"]]
    }
    out <- hypothesis_BF.default(
      posterior      = target[["posterior"]],
      hypothesis     = target[["hypothesis"]],
      parameter      = target[["parameter"]],
      logBF          = logBF,
      BF01           = BF01,
      seed           = seed,
      columns        = columns,
      density_method = method
    )
    if (precomputed) {
      out <- .hypothesis_model_append_ordinate_warnings(
        table     = out,
        posterior = target[["posterior"]],
        backend   = backend,
        request   = request
      )
    }
    # One linear target: row i is statement i.
    return(.hypothesis_model_set_row_names(
      out       = out,
      row_names = .hypothesis_model_row_names(
        labels     = rep(parameter, nrow(out)),
        statements = seq_len(nrow(out))
      )
    ))
  }

  inference_type <- plan[["inference_type"]]
  if (precomputed) {
    attached <- .hypothesis_model_marginal_attach(
      provider        = provider,
      backend         = backend,
      plan            = plan,
      parameter       = parameter,
      hypothesis      = hypothesis,
      inference_type  = inference_type,
      density_control = density_control
    )
    provider <- attached[["provider"]]
    request  <- attached[["request"]]
  }
  inference <- provider[["inference"]]
  inference[["conditional"]] <- inference[[inference_type]]
  class(inference) <- unique(c(class(inference), "marginal_inference"))

  out <- hypothesis_BF.default(
    posterior      = inference,
    hypothesis     = hypothesis,
    parameter      = parameter,
    logBF          = logBF,
    BF01           = BF01,
    seed           = seed,
    columns        = columns,
    density_method = method
  )
  if (precomputed) {
    out <- .hypothesis_model_append_ordinate_warnings(
      table     = out,
      posterior = inference[["conditional"]][[parameter]],
      backend   = backend,
      request   = request
    )
  }

  out
}


# Attaches the ordinates of the point values of the marginal means of one
# parameter to the stored marginal inference of the provider. Returns the
# provider with the ordinates and the request that the backend received (the
# first one for a group with several levels; the requests of a group differ in
# their level, values, and spec only).
.hypothesis_model_marginal_attach <- function(provider, backend, plan,
                                              parameter, hypothesis,
                                              inference_type, density_control) {

  point_refs <- .hypothesis_model_point_refs(hypothesis, parameter)
  samples    <- provider[["inference"]][[inference_type]][[parameter]]
  if (nrow(point_refs) == 0L) {
    request <- .hypothesis_model_request(
      provider, backend, plan, kind = "level", parameter = parameter,
      label = parameter, values = numeric(), sample = NULL
    )
    return(list(provider = provider, request = request))
  }
  point_refs <- .hypothesis_model_marginal_resolve_singleton_levels(
    samples    = samples,
    point_refs = point_refs
  )

  provider[["inference"]][[inference_type]][[parameter]] <-
    .hypothesis_model_keep_requested_ordinates(
      posterior  = samples,
      point_refs = point_refs
    )
  cache    <- new.env(parent = emptyenv())
  requests <- list()
  for (i in seq_len(nrow(point_refs))) {
    samples <- provider[["inference"]][[inference_type]][[parameter]]
    if (is.null(samples)) {
      stop(
        "'marginal_means' object does not contain ", inference_type,
        " samples for '", parameter, "'.",
        call. = FALSE
      )
    }
    ref     <- point_refs[i, , drop = FALSE]
    request <- .hypothesis_model_request(
      provider, backend, plan, kind = "level", parameter = parameter,
      label = paste0(parameter, "[", ref[["level"]], "]"),
      values = ref[["value"]], sample = samples[[ref[["level"]]]],
      level = ref[["level"]]
    )
    provider[["inference"]][[inference_type]][[parameter]] <-
      .hypothesis_model_attach_request(
        posterior = samples,
        backend   = backend,
        request   = request,
        control   = density_control,
        cache     = cache
      )
    requests[[i]] <- request
  }

  list(provider = provider, request = requests[[1L]])
}


# ---------------------------------------------------------------------------- #
# Entry
# ---------------------------------------------------------------------------- #

#' @title Hypothesis Bayes factors of a fitted model
#'
#' @description \code{hypothesis_BF_model()} is the engine behind
#' [hypothesis_BF()] of a \code{BayesTools_fit}, exported for packages that
#' test hypotheses on the parameters of their own models. It parses the
#' statements with the parameter catalog of the model, resolves each to one
#' tested parameter, plans how it can be evaluated (the draws and the prior
#' density of the tested quantity, the classification of its point values),
#' checks every plan against the requested density backend, evaluates the
#' statements in groups that share a tested quantity, and binds the results
#' into one table in the order of the statements. What the engine needs to
#' know about the model comes from its provider ([hypothesis_provider()]); the
#' methods to evaluate a point hypothesis are the density backends
#' ([hypothesis_density_backend()]).
#'
#' @details Every statement is planned before any is evaluated: the plan
#' records the tested target (its weights on the fitted coefficients, its
#' prior density and declared atoms), the classification of each point value
#' by [prior_ordinate_status()], and what refuses it. A statement that the
#' requested density backend cannot evaluate stops with its refusal before
#' anything is computed. The refusals are classed conditions; callers match
#' the class and never the message:
#' \describe{
#'   \item{\code{BayesTools_hypothesis_statement}}{a statement that has to be
#'   restated, or an argument that has to change: point expressions that are
#'   no direct reference (for example \code{"2 * mu = 0"}), nonlinear
#'   expressions of factor levels, one quantity named by several names, and
#'   references that BayesTools refuses (\code{BayesTools_parameter_resolution_error},
#'   for example \code{BayesTools_parameter_not_found} for an unknown name or
#'   level, which keeps its fields \code{alias} and \code{available}).
#'   \code{BayesTools_hypothesis_ambiguous} (parent
#'   \code{BayesTools_hypothesis_statement}) marks a statement whose
#'   references still name several model parameters, and
#'   \code{BayesTools_hypothesis_component_mismatch} a statement whose
#'   parameters belong to another component than the requested one;
#'   \code{BayesTools_hypothesis_parameter_mismatch} a marginal-means
#'   statement that refers to another parameter than \code{parameter}.}
#'   \item{\code{BayesTools_hypothesis_unavailable}}{a test that cannot be
#'   computed, with the more specific \code{BayesTools_hypothesis_fixed} (the
#'   quantity is fixed by the fitted model), \code{BayesTools_hypothesis_target}
#'   (no supported target or prior density), or
#'   \code{BayesTools_hypothesis_method} (the density backend is unavailable
#'   for the target). A point hypothesis whose prior ordinate is not exact,
#'   finite, and regular stops with the class of [prior_ordinate_status()]
#'   (for example \code{BayesTools_point_mass_at_null}) and
#'   \code{BayesTools_hypothesis_ordinate}.}
#' }
#' A density backend owns its own conditions (a failed ordinate, a model it
#' does not support).
#'
#' Random numbers: the draws of every statement are extracted while it is
#' planned, without \code{seed}; \code{seed} only reaches the evaluation of
#' the Bayes factors of the draws ([hypothesis_BF()]) and the prior draws of
#' random-effect quantities that have no prior density.
#'
#' @param object a model object with a [hypothesis_provider()] method.
#' @param hypothesis character vector with hypothesis statements (see
#'   [hypothesis_BF()]). Names of the parameter catalog of the model need no
#'   backticks.
#' @param namespace the component of the model the statements refer to, in
#'   the vocabulary of its provider (for a \code{BayesTools_fit}: \code{"auto"},
#'   the formula parameters such as \code{"mu"}, \code{"model"}, or
#'   \code{"random"}); \code{NULL} selects the default.
#' @param component a catalog component (the level token of a factor term)
#'   passed to [hypothesis_resolve()]: a reference to a factor term without a
#'   level stands for this level, so that \code{"x_fac = 0.1"} with
#'   \code{component = "B"} is the statement \code{"x_fac[B] = 0.1"}. It
#'   combines with \code{namespace}, which selects the formula parameter (for
#'   example when the factor is shared by several formulas).
#' @param transform_scaled whether coefficients of scaled predictors, and
#'   random-effect SDs, variances, and correlations of such predictors, are
#'   tested on the original scale (\code{TRUE}) or on the fitted scale
#'   (\code{FALSE}, the default), as in [JAGS_estimates_table()].
#' @param conditional whether to test each quantity conditionally on its
#'   inclusion event: spike-and-slab and mixture coefficients through
#'   [as_mixed_posteriors()] (\code{conditional} parameter, rule \code{"AND"}),
#'   random-effect quantities through [parameter_mixed_posterior()]. Without
#'   it a point null on a quantity with a prior point mass at that value stops
#'   with \code{BayesTools_point_mass_at_null}; region hypotheses are tested.
#' @param logBF,BF01,seed,columns see [hypothesis_BF()].
#' @param density_method the name of a registered density backend
#'   ([hypothesis_density_backends()]); matching ignores case and accepts an
#'   unambiguous partial name. \code{NULL} selects the first backend of the
#'   provider (the stored method of a marginal-means provider).
#' @param density_control the settings of the density backend (a named list
#'   that the backend validates), or \code{NULL}.
#' @param n_samples the number of prior samples or grid points for the
#'   deterministic marginal prior densities and for the prior draws of
#'   random-effect quantities.
#' @param parameter for a marginal-means provider, the marginal-means
#'   parameter that the statements refer to.
#' @param conditional_omitted whether the caller omitted \code{conditional};
#'   the engine then warns when the full ensemble of a mixture prior is
#'   tested (\code{full_ensemble()} of the provider).
#'
#' @return A \code{BayesTools_hypothesis_BF} table, as returned by
#'   [hypothesis_BF()], with the rows labelled as the statements were written.
#'   When an ordinate backend computed ordinates, the table can carry the
#'   \code{density_diagnostics} attribute of the backend.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [hypothesis_BF()], [hypothesis_provider()],
#'   [hypothesis_density_backend()], [hypothesis_quantities()]
#' @export
hypothesis_BF_model <- function(object, hypothesis, namespace = NULL,
                                component = NULL, transform_scaled = FALSE,
                                conditional = FALSE, logBF = FALSE,
                                BF01 = FALSE, seed = NULL,
                                density_method = NULL, density_control = NULL,
                                n_samples = 10000, columns = "default",
                                parameter = NULL, conditional_omitted = FALSE){

  provider <- hypothesis_provider(object)
  .hypothesis_model_validate_provider(provider)
  check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  check_char(namespace, "namespace", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(component, "component", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(transform_scaled, "transform_scaled", allow_NA = FALSE)
  check_bool(conditional, "conditional", allow_NA = FALSE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01, "BF01", allow_NA = FALSE)
  check_real(seed, "seed", check_length = 1, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  check_int(n_samples, "n_samples", lower = 2)
  check_bool(conditional_omitted, "conditional_omitted", allow_NA = FALSE)

  if(identical(provider[["kind"]], "marginal_means")){
    return(.hypothesis_model_marginal_means(
      provider = provider, hypothesis = hypothesis, parameter = parameter,
      logBF = logBF, BF01 = BF01, seed = seed,
      density_method = density_method, density_control = density_control,
      columns = columns
    ))
  }
  if(!is.null(parameter)){
    stop("'parameter' is used by marginal-means providers only.", call. = FALSE)
  }

  density_method <- .hypothesis_model_density_method(
    density_method, .hypothesis_model_density_methods(provider)
  )
  backend <- .hypothesis_model_backend(density_method)
  if(identical(backend[["type"]], "ordinate") || !is.null(density_control)){
    density_control <- .hypothesis_model_backend_control(
      backend, density_control, provider
    )
  }
  metadata <- provider[["metadata"]]()

  # Every statement is planned before any is evaluated; a statement the
  # requested method cannot evaluate stops with its plan's refusal. The
  # planning and the evaluation have one catch point for BayesTools'
  # refusals of the statements' references, which stop as statement errors.
  .hypothesis_model_catch_statement_resolution({
    ast <- .hypothesis_model_ast(hypothesis, metadata[["catalog"]])
    plans <- .hypothesis_model_plans(
      provider         = provider,
      hypothesis       = ast,
      component        = namespace,
      level_component  = component,
      transform_scaled = transform_scaled,
      conditional      = conditional,
      metadata         = metadata,
      n_samples        = n_samples
    )
    for (plan in plans) {
      .hypothesis_model_plan_check(plan, provider, backend, density_control)
    }
    if (is.function(provider[["check_plans"]])) {
      provider[["check_plans"]](plans)
    }

    keys   <- vapply(plans, `[[`, character(1), "group")
    groups <- lapply(unique(keys), function(key) which(keys == key))
    results <- lapply(groups, function(rows) {
      .hypothesis_model_plan_execute(
        plans               = plans[rows],
        provider            = provider,
        backend             = backend,
        conditional_omitted = conditional_omitted,
        logBF               = logBF,
        BF01                = BF01,
        seed                = seed,
        density_control     = density_control,
        columns             = columns
      )
    })
    if (length(results) == 1L) {
      results[[1L]]
    } else {
      .hypothesis_model_bind_parameter_results(
        results    = results,
        groups     = groups,
        hypothesis = ast,
        backend    = backend
      )
    }
  })
}


# The control that a marginal-means request evaluates with: 'normalized' (the
# control of the call, or the stored one of the provider for an ordinate
# backend, normalized by the backend), 'check' (what the availability of the
# backend is checked with: the normalized control, or the control as given
# when it cannot be normalized), and the 'error' that stopped the
# normalization. Only ordinate backends and calls with a control normalize.
.hypothesis_model_marginal_means_control <- function(provider, backend,
                                                     density_method,
                                                     density_control) {

  ordinate <- identical(backend[["type"]], "ordinate")
  if (!ordinate && is.null(density_control)) {
    return(list(normalized = NULL, check = NULL, error = NULL))
  }
  given <- density_control
  if (is.null(given) && ordinate) {
    given <- provider[["density_control_default"]](density_method)
  }
  tryCatch(
    {
      normalized <- .hypothesis_model_backend_control(backend, given, provider)
      list(normalized = normalized, check = normalized, error = NULL)
    },
    error = function(error) {
      list(normalized = NULL, check = given, error = error)
    }
  )
}


# The hypotheses on the stored marginal means of a model: the statements are
# parsed without a catalog and each is planned on the marginal means of the
# parameter it names (or of 'parameter').
.hypothesis_model_marginal_means <- function(provider, hypothesis, parameter,
                                             logBF, BF01, seed,
                                             density_method, density_control,
                                             columns){

  .hypothesis_model_catch_statement_resolution({
    ast <- hypothesis_parse(hypothesis)
    if(is.null(density_method)){
      density_method <- provider[["density_method_default"]]()
    }
    density_method <- .hypothesis_model_density_method(
      density_method, .hypothesis_model_density_methods(provider)
    )
    backend <- .hypothesis_model_backend(density_method)
    cache   <- .hypothesis_model_plan_cache(provider)
    plans   <- lapply(hypothesis_render(ast), function(statement) {
      selected <- provider[["select"]](hypothesis_parse(statement), parameter)
      .hypothesis_model_plan_marginal_means(
        provider  = provider,
        statement = statement,
        parameter = selected[["parameter"]],
        cache     = cache
      )
    })
    # The control of the backend: the settings of the call, or for an ordinate
    # backend the stored settings of the provider, normalized by the backend as
    # on fitted models. The backend sees the normalized control in
    # available(). A control that cannot be normalized does not stop the
    # checks of the plans: they use the control as given, and the failure is
    # raised after them and the request validation.
    control <- .hypothesis_model_marginal_means_control(
      provider, backend, density_method, density_control
    )
    for (plan in plans) {
      .hypothesis_model_plan_check(plan, provider, backend, control[["check"]])
    }
    refusal <- .hypothesis_model_plan_marginal_means_request_refusal(plans, provider)
    if (!is.null(refusal)) {
      .hypothesis_model_stop(refusal)
    }
    if (is.function(provider[["check_plans"]])) {
      provider[["check_plans"]](plans)
    }
    point       <- any(vapply(plans, `[[`, logical(1), "point"))
    precomputed <- identical(backend[["type"]], "ordinate") && point
    if (precomputed || !is.null(density_control)) {
      if (!is.null(control[["error"]])) {
        stop(control[["error"]])
      }
    }
    if (precomputed) {
      density_control <- control[["normalized"]]
    }

    keys   <- vapply(plans, `[[`, character(1), "group")
    groups <- lapply(unique(keys), function(key) which(keys == key))
    results <- lapply(groups, function(rows) {
      .hypothesis_model_plan_execute_marginal_means(
        plans           = plans[rows],
        provider        = provider,
        backend         = backend,
        logBF           = logBF,
        BF01            = BF01,
        seed            = seed,
        density_control = density_control,
        columns         = columns
      )
    })
    if (length(results) == 1L) {
      results[[1L]]
    } else {
      .hypothesis_model_bind_parameter_results(
        results    = results,
        groups     = groups,
        hypothesis = .hypothesis_model_plan_group_ast(plans, "hypothesis"),
        backend    = backend
      )
    }
  })
}


#' @title Drop the caches of the hypothesis engine
#'
#' @description Drops the cached plans of [hypothesis_BF_model()] and
#' [hypothesis_quantities()] (draws, selections, and plans are kept per model
#' object for the session) and the caches of the registered density backends
#' ([hypothesis_density_backend()]). A package calls it when a function the
#' plans are built with has been replaced (for example after reloading).
#'
#' @return \code{NULL}, invisibly.
#' @export
hypothesis_cache_clear <- function(){

  .hypothesis_model_plan_cache_clear()
}


# A warning about arguments that a method does not use.
.hypothesis_model_warn_unused_dots <- function(dots, allowed, caller){

  if(length(dots) == 0L){
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  if(is.null(dot_names)){
    dot_names <- rep("", length(dots))
  }
  unused <- dot_names[!nzchar(dot_names) | !dot_names %in% allowed]
  if(length(unused) == 0L){
    return(invisible(TRUE))
  }
  unused[!nzchar(unused)] <- "<unnamed>"
  warning(
    "Unused argument", if(length(unused) > 1L) "s" else "",
    " in ", caller, ": ",
    paste0("'", unused, "'", collapse = ", "),
    call. = FALSE
  )

  invisible(FALSE)
}
