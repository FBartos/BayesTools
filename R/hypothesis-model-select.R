# ============================================================================ #
# hypothesis-model-select.R
# ============================================================================ #
#
# Parameter selection of the fit-level hypothesis engine: a statement is parsed
# against the catalog of a provider (hypothesis_provider()), resolved to one
# tested model parameter, and rewritten onto it. See
# .agents/instructions/hypothesis.md.
#
# ============================================================================ #


# The component of a request in the vocabulary of the provider: NULL selects
# the default, the synonyms of the provider map to their value, and the value
# must be one of the provider's.
.hypothesis_model_component_normalize <- function(component, components,
                                                  argument = "component") {

  if (is.null(component)) {
    component <- components[["default"]]
  }
  check_char(component, argument, check_length = 1, allow_NA = FALSE)
  component <- match.arg(component, components[["values"]])
  synonym   <- components[["synonyms"]][component]
  if (length(synonym) == 1L && !is.na(synonym)) {
    return(unname(synonym))
  }

  component
}


# A selection whose component differs from the requested 'component': the
# error of class "BayesTools_hypothesis_component_mismatch" at every entry
# point; 'class' prefixes the classes of the entry point.
.hypothesis_model_stop_component_mismatch <- function(message, class = NULL) {

  stop(structure(
    class = c(class, "BayesTools_hypothesis_component_mismatch", "error",
              "condition"),
    list(message = message, call = NULL)
  ))
}


.hypothesis_model_check_component_compatible <- function(component,
                                                         selected_component,
                                                         argument,
                                                         class = NULL) {

  if (!identical(component, "auto") &&
      !identical(component, selected_component)) {
    .hypothesis_model_stop_component_mismatch(
      paste0(
        "The '", argument, "' argument selects component = '",
        selected_component, "' but 'component' was set to '", component, "'."
      ),
      class = class
    )
  }

  invisible(TRUE)
}


# The rows of the metadata entries of the quantities 'quantity_ids' (their
# own entries and the entries grouping them).
.hypothesis_model_entries_for_quantities <- function(entries, quantity_ids) {

  direct <- entries[["quantity_id"]] %in% quantity_ids
  grouped <- vapply(
    entries[["member_quantity_ids"]],
    function(members) any(members %in% quantity_ids),
    logical(1)
  )

  entries[direct | grouped, , drop = FALSE]
}


.hypothesis_model_group_ambiguity <- function(entries, quantity_ids,
                                              component) {

  candidates <- .hypothesis_model_entries_for_quantities(
    entries      = entries,
    quantity_ids = quantity_ids
  )
  grouped <- candidates[
    candidates[["role"]] == "formula_coefficient_group",
    ,
    drop = FALSE
  ]
  retry_component <- if (identical(component, "auto") &&
                         nrow(candidates) == 1L && nrow(grouped) == 1L) {
    grouped[["component"]]
  } else {
    NULL
  }

  list(candidates = candidates, component = retry_component)
}


.hypothesis_model_select_parameter <- function(provider, hypothesis,
                                               component, metadata = NULL,
                                               level_component = NULL) {

  component          <- .hypothesis_model_component_normalize(
    component, provider[["components"]]
  )
  if (is.null(metadata)) {
    metadata <- provider[["metadata"]]()
  }
  ast                <- .hypothesis_model_ast(hypothesis, metadata[["catalog"]])
  occurrences        <- hypothesis_symbols(
    ast,
    occurrences = TRUE
  )
  resolver_component <- if (!is.null(level_component)) {
    level_component
  } else if (identical(component, "auto") ||
             any(!is.na(occurrences[["level"]]))) {
    NULL
  } else {
    component
  }
  # 'group_component' is the component with which an ambiguity between a
  # factor term and its contrast coefficients is resolved to the term
  # (.hypothesis_model_group_ambiguity()).
  resolve <- function(resolver_component, group_component = component) {
    tryCatch(
      hypothesis_resolve(
        ast       = ast,
        catalog   = metadata[["catalog"]],
        component = resolver_component,
        simplify_names = TRUE
      ),
      BayesTools_parameter_ambiguous = function(error) {
        ambiguity <- .hypothesis_model_group_ambiguity(
          entries      = metadata[["entries"]],
          quantity_ids = error[["candidates"]][["quantity_id"]],
          component    = group_component
        )
        if (!is.null(ambiguity[["component"]])) {
          return(hypothesis_resolve(
            ast       = ast,
            catalog   = metadata[["catalog"]],
            component = ambiguity[["component"]],
            simplify_names = TRUE
          ))
        }
        candidates <- ambiguity[["candidates"]]
        namespace  <- .hypothesis_model_component_namespace(
          candidates = candidates,
          quantities = error[["candidates"]],
          component  = component
        )
        if (!is.null(namespace)) {
          return(hypothesis_resolve(
            ast       = ast,
            catalog   = metadata[["catalog"]],
            namespace = namespace,
            component = resolver_component,
            simplify_names = TRUE
          ))
        }
        if (length(unique(candidates[["parameter"]])) > 1L) {
          .hypothesis_model_stop_multiple_parameters(candidates, provider)
        }
        # Candidates of no several model parameters (e.g. contrast
        # coefficients without parameter entries): BayesTools' ambiguity,
        # with the classes of an ambiguous statement first.
        class(error) <- unique(c(.hypothesis_model_ambiguous_class(), class(error)))
        stop(error)
      }
    )
  }
  # BayesTools' refusals of the statement's references are statement errors:
  # BayesTools' condition, with its classes and fields, after
  # "BayesTools_hypothesis_statement" (a statement without parameter symbols,
  # "BayesTools_hypothesis_no_parameters"; a level of another component,
  # "BayesTools_hypothesis_component_mismatch"; a contrast selector of a
  # level, "BayesTools_selector_unavailable"; an unknown reference). A
  # reference that is unknown within an explicit 'component' (whose aliases
  # it is not) but names quantities without the component is resolved
  # without it: a name of a quantity of the component (e.g. its display
  # label, such as '(mu) intercept' for the formula parameter 'mu') is
  # accepted; a name known only outside the component ('outside') is refused
  # by the checks below as a component mismatch (or as the target it is),
  # also next to references of the component. Every reference of the statement
  # counts, not only the unknown name that BayesTools reports first, so that a
  # display label of the component next to a name of another component is a
  # mismatch in either order. A reference unknown in every component is the
  # unknown name refused.
  outside_component <- function(resolution) {
    inside <- vapply(
      resolution[["occurrences"]][["quantity_id"]],
      function(quantity_id) {
        own <- .hypothesis_model_entries_for_quantities(
          entries      = metadata[["entries"]],
          quantity_ids = quantity_id
        )
        any(own[["component"]] == component)
      },
      logical(1)
    )
    !all(inside)
  }
  outside  <- FALSE
  resolved <- tryCatch(
    resolve(resolver_component),
    BayesTools_parameter_not_found = function(error) {
      if (!is.null(resolver_component) && is.null(level_component) &&
          NROW(metadata[["entries"]]) > 0L) {
        unrestricted <- tryCatch(
          resolve(NULL, group_component = "auto"),
          BayesTools_parameter_not_found = .hypothesis_model_stop_statement_condition
        )
        outside <<- outside_component(unrestricted)
        return(unrestricted)
      }
      .hypothesis_model_stop_statement_condition(error)
    },
    BayesTools_selector_unavailable          = .hypothesis_model_stop_statement_condition,
    BayesTools_hypothesis_no_parameters      = .hypothesis_model_stop_statement_condition,
    BayesTools_hypothesis_component_mismatch = .hypothesis_model_stop_statement_condition
  )
  # A parameter catalog without entries is metadata of an older build.
  if (NROW(metadata[["entries"]]) == 0L) {
    .bt_stop_refit_required(
      "Resolved hypothesis metadata are unavailable. Refit the model with ",
      "the current ", provider[["build"]], " build.",
      class = provider[["refit_class"]]
    )
  }
  quantity_ids <- unique(resolved[["occurrences"]][["quantity_id"]])
  if (length(quantity_ids) == 0L) {
    stop("Internal error: the hypothesis resolved to no catalog quantities.",
         call. = FALSE)
  }
  entries <- .hypothesis_model_entries_for_quantities(
    entries      = metadata[["entries"]],
    quantity_ids = quantity_ids
  )
  covered <- vapply(quantity_ids, function(quantity_id) {
    any(entries[["quantity_id"]] == quantity_id) ||
      any(vapply(
        entries[["member_quantity_ids"]],
        function(members) quantity_id %in% members,
        logical(1)
      ))
  }, logical(1))
  if (!all(covered)) {
    occurrences <- resolved[["occurrences"]]
    .hypothesis_model_stop(provider[["uncovered"]](
      resolution   = resolved,
      quantity_ids = quantity_ids[!covered],
      symbol       = occurrences[["symbol"]][
        match(quantity_ids[!covered][[1L]], occurrences[["quantity_id"]])
      ],
      metadata     = metadata
    ))
  }
  # A statement whose parameters do not belong to the requested component is
  # a component mismatch of the statement; so is a statement with a reference
  # known only outside the component. A statement with levels is resolved
  # without the component (see 'resolver_component'), so each of its
  # references is checked here the same way: a reference of another
  # component's parameter next to one of the component (e.g.
  # 'log_tau_g[a] > mu_g[a]' with component = "mods") is a mismatch, not a
  # reference dropped when the entries are narrowed to the component.
  if (!identical(component, "auto")) {
    if (is.null(resolver_component)) {
      outside <- outside_component(resolved)
    }
    compatible <- entries[["component"]] == component
    if (!any(compatible) || (outside && !all(compatible))) {
      selected_components <- unique(entries[["component"]][!compatible])
      if (length(selected_components) == 1L) {
        .hypothesis_model_check_component_compatible(
          component          = component,
          selected_component = selected_components,
          argument           = "hypothesis",
          class              = "BayesTools_hypothesis_statement"
        )
      }
      .hypothesis_model_stop_component_mismatch(
        paste0("The hypothesis does not resolve to component = '", component, "'."),
        class = "BayesTools_hypothesis_statement"
      )
    }
    entries <- entries[compatible, , drop = FALSE]
  }
  if (nrow(entries) > 1L) {
    .hypothesis_model_stop_multiple_parameters(entries, provider)
  }

  entry   <- entries[1L, , drop = FALSE]
  aliases <- as.list(rep(entry[["parameter"]], length(entry[["aliases"]][[1L]])))
  names(aliases) <- entry[["aliases"]][[1L]]
  aliases[[entry[["parameter"]]]] <- entry[["parameter"]]
  # The roots with which the statement references the selected parameter:
  # every resolved occurrence of its quantity or of its levels (its aliases,
  # and display labels such as '(mu) intercept', 'exp(intercept)' or the
  # root '(mu) g' of '(mu) g[a]'), which .hypothesis_model_rewrite() maps to
  # the parameter.
  occurrences <- resolved[["occurrences"]]
  own_ids     <- c(entry[["quantity_id"]],
                   unlist(entry[["member_quantity_ids"]], use.names = FALSE))
  roots       <- unique(occurrences[["parameter"]][
    occurrences[["quantity_id"]] %in% own_ids
  ])
  return(list(
    parameter  = entry[["parameter"]],
    aliases    = aliases,
    roots      = roots,
    component  = entry[["component"]],
    entry      = as.list(entry[1L, setdiff(names(entry), "aliases"), drop = FALSE]),
    resolution = resolved
  ))
}


# A statement whose references name several model parameters is ambiguous
# ("BayesTools_hypothesis_ambiguous" with the parent
# "BayesTools_hypothesis_statement", .hypothesis_model_ambiguous_class());
# 'component' selects one of them. The provider words the advice.
.hypothesis_model_stop_multiple_parameters <- function(entries, provider) {

  .hypothesis_model_stop(.hypothesis_model_refusal(
    paste0(
      "Hypothesis references multiple model parameters (",
      paste(unique(entries[["parameter"]]), collapse = ", "),
      "). ", provider[["components"]][["advice"]]
    ),
    "ambiguous"
  ))
}


# Level references are resolved without the requested component: the catalog
# component of a level quantity is its level. When a level reference names
# levels of several parameters (e.g. the alias of a factor term of both the
# location and the scale formula), an explicit component restricts them to
# its own parameter, as the resolver restricts unleveled references: the
# namespace of that parameter's candidate quantities, in which the
# hypothesis is resolved again. NULL when the component is "auto" or does not
# select one parameter with one namespace; stops when no candidate belongs to
# the component. 'candidates' are the catalog entries of the ambiguous
# quantities 'quantities'.
.hypothesis_model_component_namespace <- function(candidates, quantities,
                                                  component) {

  if (identical(component, "auto") ||
      length(unique(candidates[["parameter"]])) < 2L) {
    return(NULL)
  }
  selected <- candidates[candidates[["component"]] == component, , drop = FALSE]
  if (nrow(selected) == 0L) {
    .hypothesis_model_stop_component_mismatch(
      paste0("The hypothesis does not resolve to component = '", component, "'."),
      class = "BayesTools_hypothesis_statement"
    )
  }
  if (nrow(selected) != 1L) {
    return(NULL)
  }
  members   <- c(
    selected[["quantity_id"]],
    unlist(selected[["member_quantity_ids"]], use.names = FALSE)
  )
  namespace <- unique(quantities[["namespace"]][
    quantities[["quantity_id"]] %in% members
  ])
  if (length(namespace) != 1L) {
    return(NULL)
  }

  namespace
}


# Components the provider refuses as hypothesis targets (e.g. publication-bias
# parameters of a meta-analysis), with their reason.
.hypothesis_model_check_supported_component <- function(component, provider) {

  refused <- provider[["components"]][["refused"]]
  if (length(refused) > 0L && component %in% names(refused)) {
    .hypothesis_model_stop(.hypothesis_model_refusal(
      refused[[component]],
      "target"
    ))
  }

  invisible(TRUE)
}


.hypothesis_model_alias_label <- function(aliases, parameter) {

  keep <- vapply(aliases, identical, logical(1), y = parameter)
  candidates <- names(aliases)[keep]
  candidates <- candidates[nzchar(candidates)]
  candidates <- setdiff(candidates, parameter)
  if (length(candidates) > 0L) {
    return(candidates[1L])
  }

  return(parameter)
}


# The statement written against 'parameter': the roots that name it (its
# 'aliases', and the 'resolved_roots' of its quantities, which include
# display labels that are no alias of the parameter) are rewritten to it.
.hypothesis_model_rewrite <- function(hypothesis, aliases, parameter,
                                      resolved_roots = character()) {

  ast <- .hypothesis_model_ast(hypothesis)
  mapping <- unlist(aliases, use.names = TRUE)
  mapping <- mapping[mapping == parameter & names(mapping) != mapping]
  labels  <- setdiff(resolved_roots, c(names(mapping), parameter))
  if (length(labels) > 0L) {
    mapping <- c(mapping, stats::setNames(rep(parameter, length(labels)), labels))
  }
  roots   <- hypothesis_symbols(ast)
  mapping <- mapping[names(mapping) %in% roots]
  return(.hypothesis_model_rewrite_statement(ast, mapping))
}


# Rewrites the roots of one statement ('ast') by 'mapping'. A statement that
# names one quantity by several names (e.g. "intercept > mu_intercept", or
# "mu > 0 & intercept < 1" on the location intercept) would map them onto one
# symbol, which hypothesis_rewrite() refuses without a class (the symbols
# (root, level) that 'mapping' makes equal): the statement is restated with
# one name ("BayesTools_hypothesis_statement").
.hypothesis_model_rewrite_statement <- function(ast, mapping) {

  symbols <- unique(hypothesis_symbols(ast, occurrences = TRUE)[
    , c("parameter", "level"), drop = FALSE
  ])
  roots   <- symbols[["parameter"]]
  mapped  <- roots %in% names(mapping)
  roots[mapped] <- unname(mapping[roots[mapped]])
  keys    <- paste(roots, symbols[["level"]], sep = "\r")
  shared  <- keys %in% keys[duplicated(keys)]
  if (any(shared)) {
    names <- ifelse(
      is.na(symbols[["level"]][shared]),
      symbols[["parameter"]][shared],
      paste0(symbols[["parameter"]][shared], "[", symbols[["level"]][shared], "]")
    )
    .hypothesis_model_stop(.hypothesis_model_refusal(
      paste0(
        "Hypothesis names one quantity by several names (",
        paste0("'", names, "'", collapse = ", "), "). Use one of them ",
        "throughout the statement."
      ),
      "statement"
    ))
  }

  hypothesis_rewrite(ast, mapping)
}


# The level token of a request ('component' of hypothesis_BF_model(), the
# catalog component hypothesis_resolve() filters the references by) names the
# level of a factor term that the references to the term without a level
# stand for. The statement is restated with the level, so that it is planned,
# and refused, as the statement on that level: 'x_fac = 0.1' with the token
# "B" is 'x_fac[B] = 0.1'. 'selected' is the selection of
# .hypothesis_model_select_parameter() made with the token.
.hypothesis_model_insert_level <- function(hypothesis, selected, level) {

  entry <- selected[["entry"]]
  if (is.null(level) || !identical(entry[["role"]], "formula_coefficient_group")) {
    return(hypothesis)
  }
  occurrences <- selected[["resolution"]][["occurrences"]]
  members     <- unlist(entry[["member_quantity_ids"]], use.names = FALSE)
  roots <- unique(occurrences[["parameter"]][
    is.na(occurrences[["level"]]) & occurrences[["quantity_id"]] %in% members
  ])
  if (length(roots) == 0L) {
    return(hypothesis)
  }

  ast <- .hypothesis_model_ast(hypothesis)
  ast[["statements"]] <- lapply(ast[["statements"]], function(statement) {
    for (side_name in c("left", "right")) {
      statement[[side_name]][["expression"]] <- .hypothesis_model_insert_level_node(
        statement[[side_name]][["expression"]], roots, level
      )
      statement[[side_name]] <- .bt_hypothesis_refresh_side(statement[[side_name]])
    }
    statement[["source"]] <- .bt_hypothesis_render_statement(statement)
    statement
  })
  .bt_validate_hypothesis_ast(ast)

  ast
}


.hypothesis_model_insert_level_node <- function(node, roots, level) {

  if (identical(node[["type"]], "symbol")) {
    if (node[["name"]] %in% roots) {
      node <- list(type = "level_reference", source = node[["source"]],
                   parameter = node[["name"]], level = level)
    }
  } else if (node[["type"]] %in% c("parentheses", "negation")) {
    node[["expression"]] <- .hypothesis_model_insert_level_node(
      node[["expression"]], roots, level
    )
  } else if (node[["type"]] %in% c("boolean", "comparison")) {
    node[["left"]]  <- .hypothesis_model_insert_level_node(node[["left"]], roots, level)
    node[["right"]] <- .hypothesis_model_insert_level_node(node[["right"]], roots, level)
  } else if (node[["type"]] %in% c("arithmetic", "function")) {
    node[["arguments"]] <- lapply(
      node[["arguments"]], .hypothesis_model_insert_level_node,
      roots = roots, level = level
    )
  }
  node[["source"]] <- .bt_hypothesis_render_node(node)

  node
}


.hypothesis_model_symbol_roots <- function(hypothesis) {

  return(hypothesis_symbols(.hypothesis_model_ast(hypothesis)))
}


# Hypotheses are parsed against the parameter catalog of the provider, whose
# aliases include every rendered table label, level cell and contrast
# coefficient selector. A contrast-coefficient selector of a level coordinate
# (e.g. 'g{1}' of a treatment factor), which BayesTools refuses naming the
# level form, is a statement to restate on that level.
.hypothesis_model_ast <- function(hypothesis, catalog = NULL) {

  if (inherits(hypothesis, "BayesTools_hypothesis_ast")) {
    return(hypothesis)
  }
  return(tryCatch(
    hypothesis_parse(
      hypothesis     = hypothesis,
      catalog        = catalog,
      simplify_names = !is.null(catalog)
    ),
    BayesTools_selector_unavailable = .hypothesis_model_stop_statement_condition
  ))
}


# Re-raises a BayesTools refusal of a hypothesis statement (its message,
# classes, and fields) as a statement error, with
# "BayesTools_hypothesis_statement" before its classes; a condition that is a
# statement error already (e.g. an ambiguous statement,
# "BayesTools_hypothesis_ambiguous") stops unchanged.
.hypothesis_model_stop_statement_condition <- function(condition) {

  if (!inherits(condition, "BayesTools_hypothesis_statement")) {
    class(condition) <- c("BayesTools_hypothesis_statement", class(condition))
  }
  stop(condition)
}


# The one catch point of the hypothesis engine for BayesTools' refusals of the
# statements' references, around the planning and the evaluation of the
# statements (on fitted objects and on marginal means): every resolution
# error of BayesTools ("BayesTools_parameter_resolution_error", e.g. an
# ambiguous alias, "BayesTools_parameter_ambiguous", or an unknown name, also
# one that BayesTools raises only when it evaluates a planned statement, such
# as a level of a factor term next to the whole term) stops as a statement
# error, with its message, classes, and fields after
# "BayesTools_hypothesis_statement"
# (.hypothesis_model_stop_statement_condition()). The component checks of the
# planning run first: a statement whose mismatch with an explicit 'component'
# they detect stops with "BayesTools_hypothesis_component_mismatch".
.hypothesis_model_catch_statement_resolution <- function(expr) {

  tryCatch(
    expr,
    BayesTools_parameter_resolution_error = .hypothesis_model_stop_statement_condition
  )
}


# Point references (symbol, level, value) of a hypothesis for one parameter,
# refusing non-direct expressions.
.hypothesis_model_point_refs <- function(hypothesis, parameter) {

  refs <- hypothesis_parse_point_reference(
    hypothesis     = hypothesis,
    allow_compound = TRUE
  )
  if (nrow(refs) == 0L) {
    return(data.frame(
      symbol           = character(),
      level            = character(),
      value            = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  if (any(!refs[["direct"]])) {
    unsupported <- refs[["hypothesis"]][!refs[["direct"]]][1L]
    stop(
      "Point-null hypotheses require a direct parameter or level reference; ",
      "unsupported point expression in: '",
      unsupported, "'.",
      call. = FALSE
    )
  }

  refs <- refs[refs[["parameter"]] == parameter, , drop = FALSE]
  if (nrow(refs) == 0L) {
    return(data.frame(
      symbol           = character(),
      level            = character(),
      value            = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- unique(refs[, c("symbol", "level", "value"), drop = FALSE])
  rownames(out) <- NULL
  return(out)
}


#' @title Helpers of hypothesis tables
#'
#' @description Small helpers that the entry points of a package share with
#' [hypothesis_BF_model()]. \code{hypothesis_alias_label()} is the name under
#' which the rows of a hypothesis table show a tested parameter: the first
#' alias of the parameter other than its own name, or the parameter.
#' \code{hypothesis_unique_named_warnings()} drops repeated warnings of a
#' table (the same text for the same row).
#'
#' @param aliases a list named by aliases, each element the parameter the alias
#'   names.
#' @param parameter the name of the tested parameter.
#' @param warnings a character vector of warnings, named by the rows of a
#'   hypothesis table they belong to.
#'
#' @return \code{hypothesis_alias_label()} returns a string.
#'   \code{hypothesis_unique_named_warnings()} returns the unique, nonempty
#'   warnings, or \code{NULL} without any.
#'
#' @seealso [hypothesis_BF_model()]
#' @rdname hypothesis_table_helpers
#' @export
hypothesis_alias_label <- function(aliases, parameter){

  check_char(parameter, "parameter", check_length = 1, allow_NA = FALSE)
  if(!is.list(aliases)){
    stop("'aliases' must be a list.", call. = FALSE)
  }

  .hypothesis_model_alias_label(aliases, parameter)
}


#' @rdname hypothesis_table_helpers
#' @export
hypothesis_unique_named_warnings <- function(warnings){

  .hypothesis_model_unique_named_warnings(warnings)
}
