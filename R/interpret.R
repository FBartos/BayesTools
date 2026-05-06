#' @title Interpret ensemble inference and estimates
#'
#' @description Provides textual summary for posterior
#' distributions created by [mix_posteriors] and ensemble inference
#' created by [ensemble_inference].
#'
#' @param inference model inference created by [ensemble_inference]
#' @param samples posterior samples created by [mix_posteriors]
#' @param specification list of lists specifying the generated text.
#' Each inner list carries: (1) \code{inference} specifying the name of
#' in the \code{inference} entry and optionally \code{inference_name}
#' as a name to use in the text and \code{inference_BF_name} as a
#' symbol to be used instead of the default \code{"BF"}. Finite-sample
#' BF bounds can be supplied either as a \code{bound_operator} attribute
#' on the numeric BF or as \code{inference_BF_bound_operator} /
#' \code{BF_bound_operator} with value \code{"<"} or \code{">"}. (2) \code{samples}
#' specifying the name of in the \code{samples} entry and optionally
#' \code{samples_name} as a name to use in the text, \code{samples_units} as
#' a unit text to be appended after the estimate, and \code{samples_conditional}
#' specifying whether the estimate is conditional or model-averaged.
#' @param method character specifying name of the method to be
#' appended at the beginning of each sentence.
#'
#'
#' @return \code{interpret} returns character.
#'
#' @export interpret
#' @export interpret2
#' @name interpret
#'
#' @seealso [ensemble_inference] [mix_posteriors] [BayesTools_model_tables] [BayesTools_ensemble_tables]
NULL

#' @rdname interpret
interpret <- function(inference, samples, specification, method){

  # check input
  check_list(specification, "specification", check_length = 0)
  sapply(specification, function(s){
    check_char(s$inference, allow_NULL = TRUE)
    check_char(s$inference_name, allow_NULL = TRUE)
    check_char(s$inference_BF_name, allow_NULL = TRUE)
    check_char(s$inference_BF_bound_operator, allow_values = c("<", ">"), allow_NULL = TRUE)
    check_char(s$samples, allow_NULL = TRUE)
    check_char(s$samples_name, allow_NULL = TRUE)
    check_char(s$samples_units, allow_NULL = TRUE)
    check_bool(s$samples_conditional, allow_NULL = TRUE)
  })
  check_list(inference, "inference", check_names = sapply(specification, function(s) s$inference), all_objects = TRUE, allow_other = TRUE)
  check_list(samples, "samples", check_names = unlist(sapply(specification, function(s) s$samples)), all_objects = TRUE, allow_other = TRUE)
  check_char(method, allow_NULL = TRUE)


  output <- ""

  for(i in seq_along(specification)){
    output <- paste0(output, .interpret.specification(inference, samples, specification[[i]], method), if(i != length(specification)) " ")
  }

  return(output)
}

#' @rdname interpret
interpret2                <- function(specification, method = NULL){

  # check input
  check_list(specification, "specification", check_length = 0)
  sapply(specification, function(s){
    check_char(s$inference_name,       "inference_name",        allow_NULL = TRUE)
    check_char(s$inference_BF_name,    "inference_BF_name",     allow_NULL = TRUE)
    if(!is.null(s$inference_BF) && !is.numeric(s$inference_BF))
      stop("The 'inference_BF' argument must be a numeric vector.", call. = FALSE)
    check_real(if(is.null(s$inference_BF)) NULL else as.numeric(s$inference_BF), "inference_BF", allow_NULL = TRUE)
    check_char(s$inference_BF_bound_operator, "inference_BF_bound_operator", allow_values = c("<", ">"), allow_NULL = TRUE)
    check_char(s$estimate_name,        "estimate_name",         allow_NULL = TRUE)
    check_real(s$estimate_samples,     "estimate_samples",      allow_NULL = TRUE, check_length = 0)
    check_char(s$estimate_units,       "estimate_units",        allow_NULL = TRUE)
    check_bool(s$estimate_conditional, "estimate_conditional",  allow_NULL = TRUE)
  })
  check_char(method, allow_NULL = TRUE)


  output <- ""

  for(i in seq_along(specification)){
    output <- paste0(output, .interpret.specification2(specification[[i]], method), if(i != length(specification)) " ")
  }

  return(output)
}

.interpret.specification  <- function(inference, samples, specification, method){

  temp_inference <- inference[[specification[["inference"]]]]
  temp_BF <- temp_inference[["BF"]]
  temp_BF_bound_operator <- temp_inference[["BF_bound_operator"]]
  if(is.null(temp_BF_bound_operator)){
    temp_BF_bound_operator <- temp_inference[["inference_BF_bound_operator"]]
  }
  if(is.null(temp_BF_bound_operator)){
    temp_BF_bound_operator <- specification[["inference_BF_bound_operator"]]
  }
  text_BF <- .interpret.BF(temp_BF, if(!is.null(specification[["inference_name"]])) specification[["inference_name"]] else specification[["inference"]],
                           specification[["inference_BF_name"]], temp_BF_bound_operator)

  if(is.null(specification[["samples"]])){
    return(paste0(method, " found ", text_BF, "."))
  }

  temp_par <- samples[[specification[["samples"]]]]
  text_par <- .interpret.par(temp_par, if(!is.null(specification[["samples_name"]])) specification[["samples_name"]] else specification[["samples"]],
                             specification[["samples_units"]], specification[["samples_conditional"]])

  return(paste0(method, " found ", text_BF, ", ", text_par, "."))
}
.interpret.specification2 <- function(specification, method){

  text_BF <- .interpret.BF(specification[["inference_BF"]], specification[["inference_name"]], specification[["inference_BF_name"]],
                           specification[["inference_BF_bound_operator"]])

  if(is.null(specification[["estimate_samples"]])){
    return(paste0(method, " found ", text_BF, "."))
  }

  text_par <- .interpret.par(specification[["estimate_samples"]], specification[["estimate_name"]], specification[["estimate_units"]], specification[["estimate_conditional"]])

  return(paste0(method, " found ", text_BF, ", ", text_par, "."))
}
.interpret.BF             <- function(BF, name, BF_name, BF_bound_operator = NULL){

  if(!is.numeric(BF))
    stop("The 'inference_BF' argument must be a numeric vector.", call. = FALSE)
  bound_operator <- if(is.null(BF_bound_operator)) attr(BF, "bound_operator") else BF_bound_operator
  bound_operator <- .standardize_BF_bound_operator(bound_operator, length(BF))
  BF <- as.numeric(BF)
  check_real(BF, "inference_BF", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if(!is.finite(BF))
    stop("The 'inference_BF' argument must be finite.", call. = FALSE)

  bound_operator <- bound_operator[1]
  has_bound <- !is.na(bound_operator)
  bound_identifies_side <- !has_bound ||
    (bound_operator == ">" && BF >= 1) ||
    (bound_operator == "<" && BF <= 1)

  if(!bound_identifies_side){
    text <- paste0("a finite-sample bound for the evidence concerning the ", name)
  }else if(isTRUE(all.equal(BF, 1)) && !has_bound){
    text <- paste0("no evidence for or against the ", name)
  }else if(abs(log(BF)) > log(10)){
    text <- "strong evidence"
  }else if(abs(log(BF)) > log(3)){
    text <- "moderate evidence"
  }else{
    text <- "weak evidence"
  }

  if(bound_identifies_side && !(isTRUE(all.equal(BF, 1)) && !has_bound)){
    if(BF > 1 || (has_bound && isTRUE(all.equal(BF, 1)) && bound_operator == ">")){
      text <- paste0(text, " in favor of the ", name)
    }else{
      text <- paste0(text, " against the ", name)
    }
    if(has_bound){
      text <- paste0("at least ", text)
    }
  }

  BF <- format(round(BF, if(BF < 1) 3 else 2), nsmall = if(BF < 1) 3 else 2)

  if(is.null(BF_name)){
    text <- paste0(text, if(has_bound) ", BF " else ", BF = ", if(has_bound) paste(bound_operator, BF) else BF)
  }else{
    text <- paste0(text, ", ", BF_name, if(has_bound) " " else " = ", if(has_bound) paste(bound_operator, BF) else BF)
  }

  return(text)
}
.interpret.par            <- function(samples, name, unit, conditional){

  check_real(samples, "estimate_samples", allow_NA = FALSE, check_length = 0)
  if(length(samples) == 0)
    stop("The 'estimate_samples' argument cannot be empty.", call. = FALSE)
  if(any(!is.finite(samples)))
    stop("The 'estimate_samples' argument must be finite.", call. = FALSE)
  check_char(name, "estimate_name")

  est <- mean(samples)
  CI  <- unname(stats::quantile(samples, probs = c(0.025, 0.975)))

  est <- format(round(est, 3), nsmall = 3)
  CI  <- format(round(CI, 3),  nsmall = 3)

  text <- "with mean"

  if(is.null(conditional) || !conditional){
    text <- paste0(text, " model-averaged")
  }else{
    text <- paste0(text, " conditional")
  }

  if(!is.null(unit)){
    text <- paste0(text, " estimate ", name, " = ", est, " ", unit, ", 95% CI [", CI[1], ", ", CI[2], "]")
  }else{
    text <- paste0(text, " estimate ", name, " = ", est, ", 95% CI [", CI[1], ", ", CI[2], "]")
  }

  return(text)
}


#' @title Normalize structured interpretation sources
#'
#' @description \code{interpret_records} normalizes interpretation inputs from
#' one or more structured sources into traceable evidence, estimate, header,
#' note, prior, and footnote records. It is intended for downstream packages
#' that already own the domain-specific ordering and wording of an
#' interpretation, but want BayesTools to consistently extract table rows,
#' preserve Bayes-factor metadata, and carry traceability information.
#'
#' @param sources named list of interpretation sources. Each source can be a
#' \code{data.frame}/\code{BayesTools_table}, or a list with entries
#' \code{data}, optional \code{type}, and optional per-source \code{schema}.
#' Supported source types are \code{"table"}, \code{"record"}, and
#' \code{"records"}. Table schemas can define columns such as \code{row},
#' \code{BF}, \code{central}, \code{lower}, \code{upper}, and metadata such as
#' \code{BF_orientation}, \code{BF_scale}, \code{BF_name},
#' \code{BF_bound_operator}, \code{lower_prob}, \code{upper_prob},
#' \code{interval_level}, \code{units}, and \code{conditioning}.
#' @param plan ordered list specifying which records to create. Plan items can
#' have \code{kind = "evidence"}, \code{"estimate"}, \code{"pair"},
#' \code{"test_estimate"}, \code{"for_each"}, \code{"header"}, \code{"note"},
#' \code{"prior"}, or \code{"footnote"}. Pair items contain \code{evidence}
#' and/or \code{estimate} references. Dynamic \code{"for_each"} items generate
#' records from source row order.
#' @param output whether to return normalized \code{"records"} or a simple
#' generic \code{"text"} rendering. Downstream packages should usually request
#' \code{"records"} and render domain-specific text themselves.
#' @param missing how missing optional sources or rows should be handled.
#' @param method optional method name used only by the generic text renderer.
#' @param digits number of digits used only by the generic text renderer.
#' @param spec alias for \code{plan} used by \code{interpret_tables}.
#' @param ... additional arguments passed from \code{interpret_tables} to
#' \code{interpret_records}.
#'
#' @return \code{interpret_records} returns a data frame with class
#' \code{"BayesTools_interpret_records"} when \code{output = "records"}, or a
#' character vector when \code{output = "text"}. \code{interpret_tables} is a
#' convenience wrapper around \code{interpret_records}.
#'
#' @export
interpret_records <- function(sources, plan, output = c("records", "text"),
                              missing = c("error", "skip", "warn"),
                              method = NULL, digits = 3){

  output  <- match.arg(output)
  missing <- match.arg(missing)

  check_list(sources, "sources", check_length = 0)
  check_list(plan, "plan", check_length = 0)
  check_char(method, "method", allow_NULL = TRUE)
  check_int(digits, "digits", lower = 0)

  if(is.null(names(sources)) || any(!nzchar(names(sources)))){
    stop("The 'sources' argument must be a named list.", call. = FALSE)
  }

  sources <- .interpret_prepare_sources(sources)
  plan    <- .interpret_expand_plan(plan, sources, missing)

  record_list <- list()
  for(i in seq_along(plan)){
    item_records <- .interpret_plan_item_records(
      item          = plan[[i]],
      sources       = sources,
      missing       = missing,
      default_order = i
    )
    if(length(item_records) > 0){
      record_list <- c(record_list, item_records)
    }
  }

  records <- .interpret_bind_records(record_list)
  if(nrow(records) > 0){
    records <- records[order(records$order, seq_len(nrow(records))), , drop = FALSE]
    rownames(records) <- NULL
  }

  if(output == "text"){
    return(.interpret_records_text(records, method = method, digits = digits))
  }

  return(records)
}

#' @rdname interpret_records
#' @export
interpret_tables <- function(sources, spec, ...){
  interpret_records(sources = sources, plan = spec, ...)
}

.interpret_prepare_sources <- function(sources){

  out <- vector("list", length(sources))
  names(out) <- names(sources)

  for(source_name in names(sources)){
    source <- sources[[source_name]]

    if(inherits(source, "data.frame")){
      source <- list(type = "table", data = source, schema = list())
    }else if(is.list(source) && is.null(source[["data"]]) && !is.null(source[["kind"]])){
      source <- list(type = "record", data = source, schema = list())
    }else if(!is.list(source) || is.null(source[["data"]])){
      stop(paste0("The 'sources:", source_name, "' entry must be a table or a list with a 'data' entry."), call. = FALSE)
    }

    type <- source[["type"]]
    if(is.null(type)){
      type <- .interpret_infer_source_type(source[["data"]])
    }
    check_char(type, paste0("sources:", source_name, ":type"), allow_values = c("table", "record", "records"))

    schema <- source[["schema"]]
    if(is.null(schema)){
      schema <- list()
    }
    if(length(schema) > 0){
      check_list(schema, paste0("sources:", source_name, ":schema"), check_length = 0)
    }

    out[[source_name]] <- list(
      name     = source_name,
      type     = type,
      data     = source[["data"]],
      schema   = schema,
      metadata = source[["metadata"]]
    )
  }

  return(out)
}

.interpret_infer_source_type <- function(data){

  if(inherits(data, "data.frame")){
    return("table")
  }
  if(is.list(data) && !is.null(data[["kind"]])){
    return("record")
  }
  if(is.list(data)){
    return("records")
  }

  stop("Could not infer an interpretation source type.", call. = FALSE)
}

.interpret_expand_plan <- function(plan, sources, missing){

  expanded <- list()
  for(i in seq_along(plan)){
    item <- plan[[i]]
    check_list(item, paste0("plan:", i), check_length = 0)
    kind <- .interpret_plan_kind(item)

    if(kind == "for_each"){
      expanded <- c(expanded, .interpret_expand_for_each(item, sources, missing, i))
    }else{
      item[["kind"]] <- kind
      if(is.null(item[["order"]])){
        item[["order"]] <- i
      }
      expanded[[length(expanded) + 1L]] <- item
    }
  }

  return(expanded)
}

.interpret_plan_kind <- function(item){

  kind <- item[["kind"]]
  if(is.null(kind)){
    kind <- item[["type"]]
  }
  if(is.null(kind)){
    if(!is.null(item[["evidence"]]) && !is.null(item[["estimate"]])){
      kind <- "pair"
    }else if(!is.null(item[["text"]])){
      kind <- "note"
    }else{
      kind <- "evidence"
    }
  }
  check_char(kind, "plan:kind", allow_values = c(
    "evidence", "test", "estimate", "pair", "test_estimate",
    "for_each", "header", "note", "prior", "footnote", "record"
  ))

  if(kind == "test"){
    kind <- "evidence"
  }
  if(kind == "test_estimate"){
    kind <- "pair"
  }

  return(kind)
}

.interpret_expand_for_each <- function(item, sources, missing, index){

  source_name <- item[["source"]]
  if(is.null(source_name) && !is.null(item[["evidence"]])){
    source_name <- item[["evidence"]][["source"]]
  }
  if(is.null(source_name)){
    stop("A 'for_each' plan item must define 'source' or 'evidence$source'.", call. = FALSE)
  }

  source <- .interpret_get_source(sources, source_name, item, missing)
  if(is.null(source)){
    return(list())
  }

  rows <- item[["rows"]]
  if(is.null(rows) || identical(rows, "source_order")){
    rows <- .interpret_source_row_labels(source)
  }
  if(!is.null(item[["include"]])){
    rows <- rows[rows %in% item[["include"]]]
  }
  if(!is.null(item[["exclude"]])){
    rows <- rows[!rows %in% item[["exclude"]]]
  }
  if(!is.null(item[["row_regex"]])){
    rows <- rows[grepl(item[["row_regex"]], rows)]
  }

  template <- item[["template"]]
  if(is.null(template)){
    template <- if(!is.null(item[["pair_with"]]) || !is.null(item[["estimate"]])) "pair" else "evidence"
  }
  if(template == "test_estimate"){
    template <- "pair"
  }

  expanded <- list()
  base_order <- item[["order"]]
  if(is.null(base_order)){
    base_order <- index
  }

  for(j in seq_along(rows)){
    row <- rows[[j]]
    child <- item
    child[["kind"]] <- template
    child[["rows"]] <- NULL
    child[["include"]] <- NULL
    child[["exclude"]] <- NULL
    child[["row_regex"]] <- NULL
    child[["template"]] <- NULL
    child[["pair_with"]] <- NULL
    child[["row"]] <- row
    child[["order"]] <- base_order + j
    child[["item_id"]] <- .interpret_item_id(item, row)

    if(template == "pair"){
      evidence <- item[["evidence"]]
      if(is.null(evidence)){
        evidence <- list()
      }
      evidence[["source"]] <- .interpret_or(evidence[["source"]], source_name)
      evidence[["row"]]    <- .interpret_or(evidence[["row"]], row)

      estimate <- item[["estimate"]]
      if(is.null(estimate)){
        estimate <- list()
      }
      estimate[["source"]] <- .interpret_or(estimate[["source"]], item[["pair_with"]])
      estimate[["row"]]    <- .interpret_or(estimate[["row"]], row)

      child[["evidence"]] <- evidence
      child[["estimate"]] <- estimate
    }else{
      child[["source"]] <- source_name
    }

    expanded[[length(expanded) + 1L]] <- child
  }

  return(expanded)
}

.interpret_item_id <- function(item, row){

  prefix <- .interpret_or(item[["item_id"]], item[["id"]])
  row_id <- gsub("[^A-Za-z0-9_]+", "_", as.character(row))
  row_id <- gsub("^_+|_+$", "", row_id)
  if(is.null(prefix)){
    return(row_id)
  }

  paste(prefix, row_id, sep = ".")
}

.interpret_plan_item_records <- function(item, sources, missing, default_order){

  kind <- .interpret_plan_kind(item)
  item[["kind"]] <- kind
  if(is.null(item[["order"]])){
    item[["order"]] <- default_order
  }

  if(kind %in% c("header", "note", "prior", "footnote")){
    return(list(.interpret_complete_record(
      record = list(kind = kind, text = item[["text"]]),
      item   = item
    )))
  }

  if(kind == "record"){
    return(list(.interpret_complete_record(item, item)))
  }

  if(kind == "pair"){
    records <- list()
    if(!is.null(item[["evidence"]])){
      evidence <- .interpret_ref_record(
        ref     = .interpret_merge_ref(item[["evidence"]], item, "evidence", 0),
        kind    = "evidence",
        sources = sources,
        missing = missing
      )
      if(!is.null(evidence)){
        records[[length(records) + 1L]] <- evidence
      }
    }
    if(!is.null(item[["estimate"]])){
      estimate <- .interpret_ref_record(
        ref     = .interpret_merge_ref(item[["estimate"]], item, "estimate", 0.1),
        kind    = "estimate",
        sources = sources,
        missing = missing
      )
      if(!is.null(estimate)){
        records[[length(records) + 1L]] <- estimate
      }
    }
    return(records)
  }

  record <- .interpret_ref_record(
    ref     = .interpret_merge_ref(item, item, kind, 0),
    kind    = kind,
    sources = sources,
    missing = missing
  )
  if(is.null(record)){
    return(list())
  }

  return(list(record))
}

.interpret_merge_ref <- function(ref, item, kind, order_offset){

  if(is.null(ref)){
    ref <- list()
  }
  inherited <- c(
    "section", "item_id", "id", "label", "source", "row",
    "optional", "record_id"
  )
  for(field in inherited){
    if(is.null(ref[[field]]) && !is.null(item[[field]])){
      ref[[field]] <- item[[field]]
    }
  }

  ref[["kind"]] <- kind
  ref[["order"]] <- .interpret_or(ref[["order"]], item[["order"]] + order_offset)
  ref[["item_id"]] <- .interpret_or(ref[["item_id"]], ref[["id"]])

  return(ref)
}

.interpret_ref_record <- function(ref, kind, sources, missing){

  source_name <- ref[["source"]]
  if(is.null(source_name)){
    return(.interpret_handle_missing("Interpretation plan item is missing a source.", ref, missing))
  }

  source <- .interpret_get_source(sources, source_name, ref, missing)
  if(is.null(source)){
    return(NULL)
  }

  if(source[["type"]] == "table"){
    return(switch(
      kind,
      evidence = .interpret_table_evidence_record(source, ref, missing),
      estimate = .interpret_table_estimate_record(source, ref, missing),
      stop(paste0("Table sources cannot produce records of kind '", kind, "'."), call. = FALSE)
    ))
  }

  record <- .interpret_direct_source_record(source, ref, kind, missing)
  if(is.null(record)){
    return(NULL)
  }

  return(.interpret_complete_record(record, ref, source_name, ref[["row"]]))
}

.interpret_get_source <- function(sources, source_name, item, missing){

  if(!source_name %in% names(sources)){
    return(.interpret_handle_missing(paste0("Interpretation source '", source_name, "' was not found."), item, missing))
  }

  sources[[source_name]]
}

.interpret_direct_source_record <- function(source, ref, kind, missing){

  data <- source[["data"]]

  if(source[["type"]] == "record"){
    record <- data
  }else if(source[["type"]] == "records"){
    record <- .interpret_select_direct_record(data, ref, missing)
  }else{
    stop("Internal error: unsupported direct interpretation source.", call. = FALSE)
  }
  if(is.null(record)){
    return(NULL)
  }

  record <- .interpret_merge_lists(source[["schema"]], record)
  record <- .interpret_merge_lists(record, ref)
  record[["kind"]] <- .interpret_or(record[["kind"]], kind)

  if(record[["kind"]] == "evidence"){
    record <- .interpret_complete_BF_record(record)
  }

  return(record)
}

.interpret_select_direct_record <- function(data, ref, missing){

  row <- ref[["row"]]
  if(inherits(data, "data.frame")){
    if(is.null(row)){
      if(nrow(data) == 1L){
        return(as.list(data[1, , drop = FALSE]))
      }
      return(.interpret_handle_missing("Direct records source has multiple rows; specify 'row'.", ref, missing))
    }
    row_index <- .interpret_match_record_row(data, row)
    if(is.na(row_index)){
      return(.interpret_handle_missing(paste0("Record row '", row, "' was not found."), ref, missing))
    }
    return(as.list(data[row_index, , drop = FALSE]))
  }

  if(is.null(row)){
    if(length(data) == 1L && is.list(data[[1L]])){
      return(data[[1L]])
    }
    return(.interpret_handle_missing("Direct records source has multiple records; specify 'row'.", ref, missing))
  }

  if(!is.null(names(data)) && row %in% names(data)){
    return(data[[row]])
  }

  record_ids <- vapply(data, function(x) as.character(.interpret_or(x[["record_id"]], "")), character(1))
  row_index <- match(row, record_ids)
  if(is.na(row_index)){
    return(.interpret_handle_missing(paste0("Record row '", row, "' was not found."), ref, missing))
  }

  data[[row_index]]
}

.interpret_match_record_row <- function(data, row){

  row <- as.character(row)
  candidates <- list(rownames(data))
  if("record_id" %in% colnames(data)){
    candidates[[length(candidates) + 1L]] <- as.character(data[["record_id"]])
  }
  if("row" %in% colnames(data)){
    candidates[[length(candidates) + 1L]] <- as.character(data[["row"]])
  }

  for(candidate in candidates){
    row_index <- match(row, candidate)
    if(!is.na(row_index)){
      return(row_index)
    }
  }

  NA_integer_
}

.interpret_table_evidence_record <- function(source, ref, missing){

  data <- source[["data"]]
  schema <- source[["schema"]]
  row_index <- .interpret_table_row_index(source, ref[["row"]], ref, missing)
  if(is.null(row_index)){
    return(NULL)
  }
  row_label <- .interpret_source_row_labels(source)[row_index]

  BF_column <- .interpret_table_column(
    data       = data,
    schema     = schema,
    ref        = ref,
    field      = "BF",
    candidates = .interpret_BF_candidate_columns(data)
  )
  if(is.null(BF_column)){
    return(.interpret_handle_missing("Could not identify a Bayes-factor column.", ref, missing))
  }

  BF_vector <- data[[BF_column]]
  BF_value  <- as.numeric(BF_vector[row_index])

  BF_scale <- .interpret_table_scalar("BF_scale", data, row_index, schema, ref)
  if(is.null(BF_scale)){
    BF_scale <- if(isTRUE(attr(BF_vector, "logBF"))) "log" else "linear"
  }

  BF_orientation <- .interpret_table_scalar("BF_orientation", data, row_index, schema, ref)
  if(is.null(BF_orientation)){
    BF_orientation <- .interpret_infer_BF_orientation(data, BF_column, BF_vector)
  }

  BF_name <- .interpret_table_scalar("BF_name", data, row_index, schema, ref)
  if(is.null(BF_name)){
    BF_name <- .interpret_or(attr(BF_vector, "name"), BF_column)
  }

  BF_bound_operator <- .interpret_table_scalar("BF_bound_operator", data, row_index, schema, ref)
  if(is.null(BF_bound_operator)){
    BF_bound_operator <- attr(BF_vector[row_index], "bound_operator")
  }
  if(is.null(BF_bound_operator)){
    BF_bound_operator <- attr(BF_vector, "bound_operator")
    if(!is.null(BF_bound_operator)){
      BF_bound_operator <- BF_bound_operator[row_index]
    }
  }
  BF_bound_operator <- .standardize_BF_bound_operator(BF_bound_operator, 1)

  record <- list(
    kind              = "evidence",
    source            = source[["name"]],
    row               = row_label,
    label             = .interpret_or(.interpret_table_scalar("label", data, row_index, schema, ref), row_label),
    BF_value          = BF_value,
    BF_name           = BF_name,
    BF_scale          = BF_scale,
    BF_orientation    = BF_orientation,
    BF_bound_operator = BF_bound_operator,
    BF_target         = .interpret_table_scalar("BF_target", data, row_index, schema, ref),
    BF_complement     = .interpret_table_scalar("BF_complement", data, row_index, schema, ref)
  )
  record <- .interpret_merge_lists(record, ref)

  return(.interpret_complete_record(.interpret_complete_BF_record(record), ref, source[["name"]], row_label))
}

.interpret_table_estimate_record <- function(source, ref, missing){

  data <- source[["data"]]
  schema <- source[["schema"]]
  row_index <- .interpret_table_row_index(source, ref[["row"]], ref, missing)
  if(is.null(row_index)){
    return(NULL)
  }
  row_label <- .interpret_source_row_labels(source)[row_index]

  central_column <- .interpret_table_column(
    data       = data,
    schema     = schema,
    ref        = ref,
    field      = "central",
    candidates = c("Mean", "Median", "Mode", "mean", "median", "mode")
  )
  if(is.null(central_column)){
    return(.interpret_handle_missing("Could not identify a central estimate column.", ref, missing))
  }

  interval_columns <- .interpret_interval_columns(data, schema, ref, central_column)

  central_name <- .interpret_table_scalar("central_name", data, row_index, schema, ref)
  if(is.null(central_name)){
    central_name <- tolower(central_column)
  }

  record <- list(
    kind           = "estimate",
    source         = source[["name"]],
    row            = row_label,
    label          = .interpret_or(.interpret_table_scalar("label", data, row_index, schema, ref), row_label),
    parameter      = .interpret_or(.interpret_table_scalar("parameter", data, row_index, schema, ref), row_label),
    central_name   = central_name,
    central_value  = as.numeric(data[[central_column]][row_index]),
    lower_value    = if(is.null(interval_columns[["lower"]])) NA_real_ else as.numeric(data[[interval_columns[["lower"]]]][row_index]),
    upper_value    = if(is.null(interval_columns[["upper"]])) NA_real_ else as.numeric(data[[interval_columns[["upper"]]]][row_index]),
    lower_prob     = interval_columns[["lower_prob"]],
    upper_prob     = interval_columns[["upper_prob"]],
    interval_level = interval_columns[["interval_level"]],
    units          = .interpret_table_scalar("units", data, row_index, schema, ref),
    conditioning   = .interpret_table_scalar("conditioning", data, row_index, schema, ref)
  )
  record <- .interpret_merge_lists(record, ref)

  return(.interpret_complete_record(record, ref, source[["name"]], row_label))
}

.interpret_table_row_index <- function(source, row, item, missing){

  data <- source[["data"]]
  if(is.null(row)){
    if(nrow(data) == 1L){
      return(1L)
    }
    return(.interpret_handle_missing("Table source has multiple rows; specify 'row'.", item, missing))
  }

  if(is.numeric(row)){
    if(length(row) != 1L || is.na(row) || row < 1 || row > nrow(data)){
      return(.interpret_handle_missing(paste0("Table row index '", row, "' is out of range."), item, missing))
    }
    return(as.integer(row))
  }

  row <- as.character(row)
  candidates <- list(
    .interpret_source_row_labels(source),
    attr(data, "parameters")
  )

  schema_row <- source[["schema"]][["row"]]
  if(!is.null(schema_row) && schema_row %in% colnames(data)){
    candidates[[length(candidates) + 1L]] <- as.character(data[[schema_row]])
  }

  for(candidate in candidates){
    if(is.null(candidate)){
      next
    }
    row_index <- match(row, as.character(candidate))
    if(!is.na(row_index)){
      return(row_index)
    }
  }

  .interpret_handle_missing(paste0("Table row '", row, "' was not found."), item, missing)
}

.interpret_source_row_labels <- function(source){

  data <- source[["data"]]
  if(!inherits(data, "data.frame")){
    if(is.list(data) && !is.null(names(data))){
      return(names(data))
    }
    return(as.character(seq_along(data)))
  }

  schema_row <- source[["schema"]][["row"]]
  if(!is.null(schema_row) && schema_row %in% colnames(data)){
    return(as.character(data[[schema_row]]))
  }

  rows <- rownames(data)
  if(is.null(rows)){
    rows <- as.character(seq_len(nrow(data)))
  }

  rows
}

.interpret_table_column <- function(data, schema, ref, field, candidates){

  value <- .interpret_or(ref[[field]], schema[[field]])
  value <- .interpret_or(value, ref[[paste0(field, "_column")]])
  value <- .interpret_or(value, schema[[paste0(field, "_column")]])

  if(!is.null(value) && length(value) == 1L && value %in% colnames(data)){
    return(value)
  }

  for(candidate in candidates){
    if(candidate %in% colnames(data)){
      return(candidate)
    }
  }

  NULL
}

.interpret_table_scalar <- function(field, data, row_index, schema, ref){

  value <- .interpret_or(ref[[field]], schema[[field]])
  if(is.null(value)){
    return(NULL)
  }
  if(length(value) == 1L && is.character(value) && value %in% colnames(data)){
    return(data[[value]][row_index])
  }
  if(length(value) > 1L && length(value) >= row_index){
    return(value[row_index])
  }

  value
}

.interpret_BF_candidate_columns <- function(data){

  type <- attr(data, "type")
  if(!is.null(type)){
    type_match <- which(type %in% c("inclusion_BF", "BF"))
    if(length(type_match) > 0){
      return(colnames(data)[type_match])
    }
  }

  c("inclusion_BF", "BF", "Bayes factor", "bayes_factor")
}

.interpret_infer_BF_orientation <- function(data, BF_column, BF_vector){

  type <- attr(data, "type")
  BF_type <- NULL
  if(!is.null(type) && BF_column %in% colnames(data)){
    BF_type <- type[match(BF_column, colnames(data))]
  }
  BF01 <- isTRUE(attr(BF_vector, "BF01"))

  if(identical(BF_type, "inclusion_BF") || identical(BF_column, "inclusion_BF")){
    return(if(BF01) "exclusion_over_inclusion" else "inclusion_over_exclusion")
  }

  if(BF01){
    "null_over_alternative"
  }else{
    "alternative_over_null"
  }
}

.interpret_interval_columns <- function(data, schema, ref, central_column){

  lower <- .interpret_or(ref[["lower"]], schema[["lower"]])
  upper <- .interpret_or(ref[["upper"]], schema[["upper"]])
  lower_prob <- .interpret_or(ref[["lower_prob"]], schema[["lower_prob"]])
  upper_prob <- .interpret_or(ref[["upper_prob"]], schema[["upper_prob"]])

  if(is.null(lower) && !is.null(lower_prob)){
    lower <- .interpret_probability_column(data, lower_prob)
  }
  if(is.null(upper) && !is.null(upper_prob)){
    upper <- .interpret_probability_column(data, upper_prob)
  }

  if(is.null(lower) || is.null(upper)){
    probability_columns <- .interpret_probability_columns(data)
    probability_columns <- setdiff(probability_columns, central_column)
    if(length(probability_columns) >= 2L){
      probs <- as.numeric(probability_columns)
      lower <- probability_columns[which.min(probs)]
      upper <- probability_columns[which.max(probs)]
      lower_prob <- .interpret_or(lower_prob, as.numeric(lower))
      upper_prob <- .interpret_or(upper_prob, as.numeric(upper))
    }
  }

  if(!is.null(lower) && !lower %in% colnames(data)){
    lower <- NULL
  }
  if(!is.null(upper) && !upper %in% colnames(data)){
    upper <- NULL
  }

  if(is.null(lower_prob) && !is.null(lower)){
    lower_prob <- suppressWarnings(as.numeric(lower))
  }
  if(is.null(upper_prob) && !is.null(upper)){
    upper_prob <- suppressWarnings(as.numeric(upper))
  }

  interval_level <- .interpret_or(ref[["interval_level"]], schema[["interval_level"]])
  if(is.null(interval_level) && is.finite(lower_prob) && is.finite(upper_prob)){
    interval_level <- upper_prob - lower_prob
  }

  list(
    lower          = lower,
    upper          = upper,
    lower_prob     = if(is.null(lower_prob)) NA_real_ else as.numeric(lower_prob),
    upper_prob     = if(is.null(upper_prob)) NA_real_ else as.numeric(upper_prob),
    interval_level = if(is.null(interval_level)) NA_real_ else as.numeric(interval_level)
  )
}

.interpret_probability_column <- function(data, prob){

  prob <- as.numeric(prob)
  candidates <- c(as.character(prob), format(prob, trim = TRUE, scientific = FALSE))
  candidates <- unique(candidates)
  match <- candidates[candidates %in% colnames(data)]
  if(length(match) > 0){
    return(match[1])
  }

  column_probs <- suppressWarnings(as.numeric(colnames(data)))
  numeric_match <- which(!is.na(column_probs) & abs(column_probs - prob) < sqrt(.Machine$double.eps))
  if(length(numeric_match) > 0){
    return(colnames(data)[numeric_match[1L]])
  }

  NULL
}

.interpret_probability_columns <- function(data){

  probs <- suppressWarnings(as.numeric(colnames(data)))
  colnames(data)[!is.na(probs) & probs > 0 & probs < 1]
}

.interpret_complete_BF_record <- function(record){

  BF_value <- as.numeric(record[["BF_value"]])
  BF_scale <- .interpret_or(record[["BF_scale"]], "linear")
  BF_orientation <- .interpret_or(record[["BF_orientation"]], "target_over_complement")
  BF_bound_operator <- .standardize_BF_bound_operator(record[["BF_bound_operator"]], 1)

  check_char(BF_scale, "BF_scale", allow_values = c("linear", "log"))
  check_char(BF_orientation, "BF_orientation", allow_values = c(
    "inclusion_over_exclusion", "exclusion_over_inclusion",
    "alternative_over_null", "null_over_alternative",
    "target_over_complement", "complement_over_target"
  ))

  BF_linear <- if(BF_scale == "log") exp(BF_value) else BF_value
  canonical_bound <- BF_bound_operator
  if(BF_orientation %in% c("exclusion_over_inclusion", "null_over_alternative", "complement_over_target")){
    BF_linear <- 1 / BF_linear
    canonical_bound <- .invert_BF_bound_operator(canonical_bound)
  }

  record[["BF_value"]] <- BF_value
  record[["BF_scale"]] <- BF_scale
  record[["BF_orientation"]] <- BF_orientation
  record[["BF_bound_operator"]] <- BF_bound_operator
  record[["BF_canonical_value"]] <- BF_linear
  record[["BF_canonical_bound_operator"]] <- canonical_bound
  record[["BF_canonical_orientation"]] <- "target_over_complement"

  record
}

.interpret_complete_record <- function(record, item, source_name = NULL, row = NULL){

  empty <- .interpret_empty_record()
  record <- .interpret_merge_lists(empty, record)

  record[["section"]] <- .interpret_fill(record[["section"]], item[["section"]])
  record[["item_id"]] <- .interpret_fill(.interpret_fill(record[["item_id"]], item[["item_id"]]), item[["id"]])
  record[["source"]]  <- .interpret_fill(record[["source"]], source_name)
  record[["row"]]     <- .interpret_fill(record[["row"]], row)
  record[["order"]]   <- .interpret_fill(record[["order"]], item[["order"]])
  record[["record_id"]] <- .interpret_fill(record[["record_id"]], .interpret_record_id(record))

  record
}

.interpret_empty_record <- function(){

  list(
    record_id = NA_character_,
    kind = NA_character_,
    section = NA_character_,
    item_id = NA_character_,
    source = NA_character_,
    row = NA_character_,
    order = NA_real_,
    label = NA_character_,
    text = NA_character_,
    parameter = NA_character_,
    units = NA_character_,
    conditioning = NA_character_,
    central_name = NA_character_,
    central_value = NA_real_,
    lower_value = NA_real_,
    upper_value = NA_real_,
    lower_prob = NA_real_,
    upper_prob = NA_real_,
    interval_level = NA_real_,
    BF_value = NA_real_,
    BF_name = NA_character_,
    BF_scale = NA_character_,
    BF_orientation = NA_character_,
    BF_bound_operator = NA_character_,
    BF_canonical_value = NA_real_,
    BF_canonical_bound_operator = NA_character_,
    BF_canonical_orientation = NA_character_,
    BF_target = NA_character_,
    BF_complement = NA_character_
  )
}

.interpret_record_id <- function(record){

  pieces <- c(record[["section"]], record[["item_id"]], record[["kind"]])
  pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
  if(length(pieces) == 0){
    return(NA_character_)
  }

  paste(pieces, collapse = ".")
}

.interpret_bind_records <- function(records){

  columns <- names(.interpret_empty_record())
  if(length(records) == 0){
    out <- as.data.frame(.interpret_empty_record(), stringsAsFactors = FALSE)[0, , drop = FALSE]
    class(out) <- c("BayesTools_interpret_records", class(out))
    return(out)
  }

  rows <- lapply(records, function(record){
    record <- .interpret_merge_lists(.interpret_empty_record(), record)
    record <- record[columns]
    as.data.frame(record, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("BayesTools_interpret_records", class(out))

  out
}

.interpret_records_text <- function(records, method = NULL, digits = 3){

  if(nrow(records) == 0){
    return(character())
  }

  out <- character()
  for(i in seq_len(nrow(records))){
    record <- records[i, , drop = FALSE]
    kind <- record[["kind"]]
    text <- switch(
      kind,
      header = record[["text"]],
      note = record[["text"]],
      prior = record[["text"]],
      footnote = record[["text"]],
      evidence = .interpret_record_evidence_text(record),
      estimate = .interpret_record_estimate_text(record, digits),
      NA_character_
    )
    if(!is.na(text) && nzchar(text)){
      out <- c(out, text)
    }
  }

  if(!is.null(method) && length(out) > 0 && !any(records[["kind"]] == "header")){
    out[1] <- paste(method, out[1])
  }

  class(out) <- c("BayesTools_interpret_text", "character")
  out
}

.interpret_record_evidence_text <- function(record){

  BF <- record[["BF_canonical_value"]]
  if(!is.finite(BF)){
    return(paste0("evidence for ", record[["label"]], " unavailable"))
  }
  BF_name <- .interpret_canonical_BF_name(record)
  bound <- record[["BF_canonical_bound_operator"]]
  if(is.na(bound)){
    bound <- NULL
  }

  .interpret.BF(BF, record[["label"]], BF_name, bound)
}

.interpret_canonical_BF_name <- function(record){

  if(record[["BF_orientation"]] %in% c("inclusion_over_exclusion", "exclusion_over_inclusion")){
    return("Inclusion BF")
  }
  if(record[["BF_orientation"]] %in% c("alternative_over_null", "null_over_alternative")){
    return("BF")
  }

  BF_name <- record[["BF_name"]]
  BF_scale <- record[["BF_scale"]]
  if(is.na(BF_name) || !nzchar(BF_name) || is.na(BF_scale) || BF_scale != "linear"){
    return("BF")
  }

  BF_name
}

.interpret_record_estimate_text <- function(record, digits){

  central <- format(round(record[["central_value"]], digits), nsmall = digits)
  lower   <- format(round(record[["lower_value"]], digits), nsmall = digits)
  upper   <- format(round(record[["upper_value"]], digits), nsmall = digits)
  level   <- if(is.finite(record[["interval_level"]])) paste0(round(100 * record[["interval_level"]]), "%") else "uncertainty"

  text <- paste0("with ", record[["central_name"]], " estimate ", record[["parameter"]], " = ", central)
  if(!is.na(record[["units"]]) && nzchar(record[["units"]])){
    text <- paste(text, record[["units"]])
  }
  if(is.finite(record[["lower_value"]]) && is.finite(record[["upper_value"]])){
    text <- paste0(text, ", ", level, " interval [", lower, ", ", upper, "]")
  }
  if(!is.na(record[["conditioning"]]) && nzchar(record[["conditioning"]])){
    text <- paste0(text, " (", record[["conditioning"]], ")")
  }

  text
}

.interpret_merge_lists <- function(x, y){

  if(is.null(y)){
    return(x)
  }
  for(name in names(y)){
    x[[name]] <- y[[name]]
  }
  x
}

.interpret_or <- function(x, y){
  if(is.null(x)) y else x
}

.interpret_fill <- function(x, y){
  if(is.null(x) || length(x) == 0 || (length(x) == 1L && is.na(x))) y else x
}

.interpret_handle_missing <- function(message, item, missing){

  optional <- isTRUE(item[["optional"]])
  if(optional || missing == "skip"){
    return(NULL)
  }
  if(missing == "warn"){
    warning(message, call. = FALSE)
    return(NULL)
  }

  stop(message, call. = FALSE)
}
