# ============================================================================ #
# hypothesis-model-ordinates.R
# ============================================================================ #
#
# Ordinate bookkeeping of the hypothesis engine: which stored posterior
# ordinates a request keeps or replaces, the attachment of the ordinates a
# density backend computes, and the lifting of the warnings and diagnostics
# of the attached ordinates to the rows of a hypothesis table. The backend
# computes and validates its ordinates (hypothesis_density_backend()); the
# engine only bookkeeps. See .agents/instructions/hypothesis.md.
#
# ============================================================================ #


# The single ordinates of a stored posterior-ordinate attribute.
.hypothesis_model_ordinate_entries <- function(posterior_ordinate) {

  if (is.null(posterior_ordinate)) {
    return(list())
  }

  .posterior_ordinate_entries(posterior_ordinate)
}


# One posterior-ordinate attribute holding the single ordinates 'entries'
# (NULL without any).
.hypothesis_model_ordinate_combine <- function(entries) {

  entries <- entries[!vapply(entries, is.null, logical(1))]
  if (length(entries) == 0L) {
    return(NULL)
  }
  out <- NULL
  for (entry in entries) {
    out <- posterior_ordinate_append(existing = out, ordinate = entry)
  }

  out
}


# Whether the ordinate 'entry' is at the point value 'value'.
.hypothesis_model_ordinate_value_matches <- function(entry, value) {

  entry_value <- suppressWarnings(as.numeric(entry[["value"]]))
  value       <- suppressWarnings(as.numeric(value))[1L]
  if (!is.finite(value)) {
    return(FALSE)
  }

  any(is.finite(entry_value) & entry_value == value)
}


# The stored ordinates of 'posterior_ordinate' at the point values 'values'.
.hypothesis_model_ordinate_keep_values <- function(posterior_ordinate, values) {

  entries <- .hypothesis_model_ordinate_entries(posterior_ordinate)
  if (length(entries) == 0L || length(values) == 0L) {
    return(NULL)
  }
  keep <- vapply(entries, function(entry) {
    any(vapply(values, function(value) {
      .hypothesis_model_ordinate_value_matches(entry, value)
    }, logical(1)))
  }, logical(1))

  .hypothesis_model_ordinate_combine(entries[keep])
}


# Prunes the stored ordinates of a posterior (a vector, or a list of level
# posteriors) to the requested point values 'point_refs' (columns 'level' and
# 'value').
.hypothesis_model_keep_requested_ordinates <- function(posterior, point_refs) {

  if (!is.list(posterior)) {
    values <- point_refs[["value"]][is.na(point_refs[["level"]])]
    posterior_metadata(posterior, "posterior_ordinate") <-
      .hypothesis_model_ordinate_keep_values(
        posterior_metadata(posterior, "posterior_ordinate"),
        values
      )
    return(posterior)
  }

  posterior_metadata(posterior, "posterior_ordinate") <- NULL
  for (level in names(posterior)) {
    values <- point_refs[["value"]][
      !is.na(point_refs[["level"]]) & point_refs[["level"]] == level
    ]
    posterior_metadata(posterior[[level]], "posterior_ordinate") <-
      .hypothesis_model_ordinate_keep_values(
        posterior_metadata(posterior[[level]], "posterior_ordinate"),
        values
      )
  }

  return(posterior)
}


# Attaches the ordinates that 'backend' computes for one request: the point
# values 'request$values' of a posterior (a vector) or of the level
# 'request$level' of a list of level posteriors. Stored ordinates at the
# requested values are kept only when the backend's reuse() accepts them (it
# receives the draws that carry the ordinate, to compare their conditioning
# and provenance with the request) and are removed otherwise (so that a
# replacement never collides with the unique-value rule of the ordinate
# attribute); the backend's attach() then computes the missing ones. The
# attached ordinates are checked with posterior_ordinate_supports_bf().
.hypothesis_model_attach_request <- function(posterior, backend, request,
                                             control, cache) {

  level     <- request[["level"]]
  has_level <- !is.null(level) && !is.na(level)
  sample    <- if (has_level) posterior[[level]] else posterior
  values    <- request[["values"]]

  entries   <- .hypothesis_model_ordinate_entries(
    posterior_metadata(sample, "posterior_ordinate")
  )
  requested <- vapply(entries, function(entry) {
    any(vapply(values, function(value) {
      .hypothesis_model_ordinate_value_matches(entry, value)
    }, logical(1)))
  }, logical(1))
  reusable  <- requested
  reusable[requested] <- vapply(entries[requested], function(entry) {
    is.function(backend[["reuse"]]) &&
      isTRUE(backend[["reuse"]](entry, request, control, sample))
  }, logical(1))
  posterior_metadata(sample, "posterior_ordinate") <-
    .hypothesis_model_ordinate_combine(entries[!requested | reusable])
  if (has_level) posterior[[level]] <- sample else posterior <- sample

  covered <- vapply(values, function(value) {
    any(vapply(entries[reusable], .hypothesis_model_ordinate_value_matches,
               logical(1), value = value))
  }, logical(1))
  if (!all(covered)) {
    request[["values"]] <- values[!covered]
    posterior <- backend[["attach"]](posterior, request, control, cache)
  }

  sample <- if (has_level) posterior[[level]] else posterior
  stored <- .hypothesis_model_ordinate_entries(
    posterior_metadata(sample, "posterior_ordinate")
  )
  for (value in values) {
    valid <- any(vapply(stored, function(entry) {
      .hypothesis_model_ordinate_value_matches(entry, value) &&
        posterior_ordinate_supports_bf(entry)
    }, logical(1)))
    if (!valid) {
      stop(structure(
        class = c("BayesTools_hypothesis_backend", "error", "condition"),
        list(
          message = paste0(
            "The density backend '", backend[["name"]], "' did not attach a ",
            "valid posterior ordinate for '", request[["label"]], "' at ",
            format(value, digits = 15L, trim = TRUE), "."
          ),
          call = NULL
        )
      ))
    }
  }

  posterior
}


# The warnings of the ordinates stored on 'posterior' (recursively over a list
# of level posteriors), one record per warning: the default of a backend
# without its own warnings(). A backend stores the text of a warning as
# 'warning' in the diagnostics of its ordinate.
.hypothesis_model_ordinate_warning_records <- function(posterior) {

  records <- .hypothesis_model_entry_warning_records(
    posterior_metadata(posterior, "posterior_ordinate")
  )
  if (is.list(posterior)) {
    child_records <- lapply(posterior, .hypothesis_model_ordinate_warning_records)
    child_records <- child_records[vapply(child_records, nrow, integer(1)) > 0L]
    if (length(child_records) > 0L) {
      records <- rbind(records, do.call(rbind, child_records))
    }
  }

  records <- records[!duplicated(records), , drop = FALSE]
  rownames(records) <- NULL

  records
}


.hypothesis_model_entry_warning_records <- function(posterior_ordinate) {

  entries <- .hypothesis_model_ordinate_entries(posterior_ordinate)
  records <- lapply(entries, function(entry) {

    warnings <- entry[["diagnostics"]][["warning"]]
    warnings <- unique(as.character(warnings))
    warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
    if (length(warnings) == 0L) {
      return(.hypothesis_model_empty_warning_records())
    }

    data.frame(
      warning   = warnings,
      parameter = .hypothesis_model_ordinate_character(entry, "parameter"),
      level     = .hypothesis_model_ordinate_character(entry, "level"),
      label     = .hypothesis_model_ordinate_label(entry),
      value     = .hypothesis_model_ordinate_value(entry),
      stringsAsFactors = FALSE
    )
  })
  records <- records[vapply(records, nrow, integer(1)) > 0L]
  if (length(records) == 0L) {
    return(.hypothesis_model_empty_warning_records())
  }

  do.call(rbind, records)
}


.hypothesis_model_empty_warning_records <- function() {

  data.frame(
    warning   = character(),
    parameter = character(),
    level     = character(),
    label     = character(),
    value     = numeric(),
    stringsAsFactors = FALSE
  )
}


.hypothesis_model_ordinate_character <- function(entry, name) {

  value <- entry[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_character_)
  }

  value <- as.character(value[[1L]])
  if (is.na(value) || !nzchar(value)) {
    return(NA_character_)
  }

  return(value)
}


.hypothesis_model_ordinate_label <- function(entry) {

  parameter <- .hypothesis_model_ordinate_character(entry, "parameter")
  level     <- .hypothesis_model_ordinate_character(entry, "level")
  if (!is.na(parameter) && nzchar(parameter) &&
      !is.na(level) && nzchar(level)) {
    return(paste0(parameter, "[", level, "]"))
  }
  if (!is.na(parameter) && nzchar(parameter)) {
    return(parameter)
  }

  return(NA_character_)
}


.hypothesis_model_ordinate_value <- function(entry) {

  for (name in c("value", "x", "null_hypothesis")) {
    value <- suppressWarnings(as.numeric(entry[[name]]))
    value <- value[is.finite(value)]
    if (length(value) > 0L) {
      return(value[[1L]])
    }
  }

  return(NA_real_)
}


# Lifts the diagnostics and the warnings of the ordinates stored on
# 'posterior' to the table 'table' of a hypothesis: the backend's public
# diagnostics become its 'density_diagnostics' attribute, and each warning is
# put on the rows of the statements it belongs to (the point side with its
# value and label). A backend without diagnostics() sets no attribute; one
# without warnings() has the warnings of its ordinate diagnostics lifted.
# 'request' is the request that the backend's attach() received for the group
# of statements (the first one of a group that is attached level by level).
.hypothesis_model_append_ordinate_warnings <- function(
    table, posterior, backend, request, parameter = NULL) {

  if (is.function(backend[["diagnostics"]])) {
    diagnostics <- backend[["diagnostics"]](posterior, request)
    if (is.data.frame(diagnostics) && nrow(diagnostics) > 0L) {
      attr(table, "density_diagnostics") <- diagnostics
    }
  }

  warning_records <- if (is.function(backend[["warnings"]])) {
    backend[["warnings"]](posterior, request)
  } else {
    .hypothesis_model_ordinate_warning_records(posterior)
  }
  if (nrow(warning_records) == 0L) {
    return(table)
  }
  if (!is.null(parameter)) {
    warning_records[["parameter"]] <- parameter
  }

  existing <- attr(table, "warnings", exact = TRUE)
  warnings <- unlist(lapply(seq_len(nrow(warning_records)), function(i) {

    record <- warning_records[i, , drop = FALSE]
    rows   <- .hypothesis_model_warning_rows(table, record)
    stats::setNames(rep(record[["warning"]], length(rows)), rows)
  }), use.names = TRUE)
  attr(table, "warnings") <- .hypothesis_model_unique_named_warnings(
    c(existing, warnings)
  )

  return(table)
}


.hypothesis_model_warning_rows <- function(table, record) {

  row_names <- rownames(table)
  if (is.null(row_names) || length(row_names) == 0L) {
    return(.hypothesis_model_warning_row(table))
  }

  candidates <- seq_len(nrow(table))

  point_rows <- .hypothesis_model_warning_point_rows(table, record)
  if (length(point_rows) > 0L) {
    candidates <- intersect(candidates, point_rows)
  } else if (!is.null(.hypothesis_model_result_statements(table))) {
    return(character())
  }

  label <- record[["label"]]
  if (!is.na(label) && nzchar(label)) {
    label_rows <- unique(c(
      which(row_names == label),
      which(sub(" \\([0-9]+\\)$", "", row_names) == label)
    ))
    narrowed <- intersect(candidates, label_rows)
    if (length(narrowed) == 0L) {
      narrowed <- intersect(candidates,
        .hypothesis_model_table_rows_with_labels(table, label))
    }
    if (length(narrowed) > 0L) {
      candidates <- narrowed
    }
  }

  if (length(candidates) == 0L) {
    return(.hypothesis_model_warning_row(table))
  }

  row_names[candidates]
}


.hypothesis_model_warning_point_rows <- function(table, record) {

  statements <- .hypothesis_model_result_statements(table)
  if (is.null(statements) || length(statements) == 0L ||
      nrow(table) %% length(statements) != 0L) {
    point_labels <- .hypothesis_model_warning_point_labels(table, record)
    return(.hypothesis_model_table_rows_with_labels(table, point_labels))
  }

  statement_matches <- vapply(statements, function(statement) {

    any(vapply(
      list(statement[["left"]], statement[["right"]]),
      .hypothesis_model_warning_side_matches,
      logical(1),
      record = record
    ))
  }, logical(1))
  if (!any(statement_matches)) {
    return(integer())
  }

  n_quantities <- nrow(table) %/% length(statements)
  rows <- unlist(lapply(which(statement_matches), function(i) {

    ((i - 1L) * n_quantities + 1L):(i * n_quantities)
  }), use.names = FALSE)

  rows
}


.hypothesis_model_warning_point_labels <- function(table, record) {

  statements <- .hypothesis_model_result_statements(table)
  if (is.null(statements) || length(statements) == 0L) {
    return(character())
  }

  labels <- unlist(lapply(statements, function(statement) {

    sides <- list(statement[["left"]], statement[["right"]])
    vapply(sides, function(side) {

      if (!.hypothesis_model_warning_side_matches(side, record)) {
        return(NA_character_)
      }
      side[["label"]]
    }, character(1))
  }), use.names = FALSE)
  labels <- unique(labels[!is.na(labels) & nzchar(labels)])

  return(labels)
}


.hypothesis_model_warning_side_matches <- function(side, record) {

  if (!is.list(side) || !identical(side[["type"]], "point")) {
    return(FALSE)
  }

  expression <- side[["expression"]]
  while (is.list(expression) &&
         identical(expression[["type"]], "parentheses")) {
    expression <- expression[["expression"]]
  }
  expr <- if (is.list(expression) &&
              identical(expression[["type"]], "symbol")) {
    expression[["name"]]
  } else if (is.list(expression) &&
             identical(expression[["type"]], "level_reference")) {
    paste0(expression[["parameter"]], "[", expression[["level"]], "]")
  } else {
    NULL
  }
  if (is.null(expr)) {
    label <- record[["label"]]
    # A compound ordinate can carry its expression as the declared label.
    if (!is.na(label) && nzchar(label)) {
      label_ast <- hypothesis_parse(paste0(label, " = 0"))
      expected <- label_ast[["statements"]][[1L]][["left"]][["expression"]]
      expr <- if (identical(expression, expected)) label else NULL
    }
  }
  if (is.null(expr) || length(expr) != 1L || is.na(expr) || !nzchar(expr)) {
    return(FALSE)
  }

  label     <- record[["label"]]
  parameter <- record[["parameter"]]
  terms     <- unique(c(label, parameter))
  terms     <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) > 0L && !expr %in% terms) {
    return(FALSE)
  }

  value <- suppressWarnings(as.numeric(record[["value"]]))
  if (!is.finite(value)) {
    return(TRUE)
  }
  side_value <- suppressWarnings(as.numeric(side[["value"]]))
  if (!is.finite(side_value)) {
    return(FALSE)
  }

  identical(value, side_value)
}


.hypothesis_model_table_rows_with_labels <- function(table, labels) {

  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0L) {
    return(integer())
  }

  columns <- intersect(c("Alternative", "Null"), colnames(table))
  if (length(columns) == 0L) {
    return(integer())
  }

  matches <- logical(nrow(table))
  for (column in columns) {
    values <- as.character(table[[column]])
    matches <- matches | values %in% labels
    for (i in which(!matches)) {
      parsed <- hypothesis_parse(values[[i]])
      references <- hypothesis_symbols(parsed, occurrences = TRUE)
      matches[[i]] <- matches[[i]] || any(references[["symbol"]] %in% labels)
    }
  }

  which(matches)
}


.hypothesis_model_warning_row <- function(table) {

  row_names <- rownames(table)
  if (!is.null(row_names) && length(row_names) > 0L) {
    return(row_names[[1L]])
  }

  return("density method")
}


.hypothesis_model_unique_named_warnings <- function(warnings) {

  if (length(warnings) == 0L) {
    return(NULL)
  }

  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (length(warnings) == 0L) {
    return(NULL)
  }

  warning_names <- names(warnings)
  if (is.null(warning_names)) {
    warning_names <- rep("", length(warnings))
  }
  key <- paste(warning_names, warnings, sep = "\r")

  warnings[!duplicated(key)]
}

