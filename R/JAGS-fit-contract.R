# Versioned fitted-object compatibility contract.

.bt_parameter_encoding_version <- 1L
.bt_formula_name_map_version <- 1L
.bt_fit_contract_version <- 1L

.bt_fit_contract_components <- c(
  "name_encoding",
  "formula_name_map",
  "formula_design",
  "parameter_registry",
  "draw_geometry",
  "parameter_catalog"
)

#' JAGS formula parameter encoding and fitted-object contract
#'
#' @description
#' `JAGS_parameter_encode()` creates an injective JAGS-safe semantic identifier
#' from structured formula-coordinate fields. It does not replace the backend
#' column name stored in a fit; the persisted formula name map links both names.
#' `JAGS_parameter_decode()` reverses identifiers created by the current
#' encoding.
#'
#' `JAGS_formula_name_map()` returns the persisted mapping between opaque JAGS
#' base names and semantic formula coordinates. Downstream code should use this
#' map or [JAGS_parameter_registry()] rather than parsing a JAGS name.
#'
#' `JAGS_fit_contract()` returns the compatibility profile attached to a new
#' fit. `JAGS_validate_fit_contract()` checks only the components named in
#' `requires`, so methods unrelated to new metadata can retain their existing
#' legacy behavior.
#'
#' @param fields named list with scalar character fields `kind`,
#'   `formula_parameter`, `term`, and `role`.
#' @param name encoded scalar name returned by `JAGS_parameter_encode()`.
#' @param fit fitted object created by [JAGS_fit()].
#' @param parameter optional formula parameter selecting one name map.
#' @param requires character vector naming required contract components.
#'
#' @return The encoding helpers return a character scalar or decoded list. The
#' schema helpers return data frames. `JAGS_formula_name_map()` returns a
#' `BayesTools_formula_name_map` data frame (or a named list of them).
#' `JAGS_fit_contract()` returns a named list, and the validator returns it
#' invisibly.
#'
#' @export JAGS_parameter_encode
#' @export JAGS_parameter_decode
#' @export JAGS_parameter_encoding_schema
#' @export JAGS_formula_name_map
#' @export JAGS_fit_contract
#' @export JAGS_validate_fit_contract
#' @export JAGS_fit_contract_schema
#' @name JAGS_fit_contract
NULL

#' @rdname JAGS_fit_contract
JAGS_parameter_encode <- function(fields){

  .bt_validate_parameter_encoding_fields(fields)
  encoded <- vapply(fields, .bt_utf8_hex_encode, character(1))
  paste0("BT1_", paste(encoded, collapse = "_"))
}

#' @rdname JAGS_fit_contract
JAGS_parameter_decode <- function(name){

  check_char(name, "name", check_length = 1L, allow_NA = FALSE)
  if(!grepl("^BT[0-9]+_", name)){
    stop("'name' is not a BayesTools encoded parameter name.", call. = FALSE)
  }
  version_text <- sub("^BT([0-9]+)_.*$", "\\1", name)
  version <- suppressWarnings(as.integer(version_text))
  if(!identical(version, .bt_parameter_encoding_version)){
    stop(
      "Encoded parameter name uses unsupported encoding version '",
      version_text, "'.",
      call. = FALSE
    )
  }
  match <- regexec(
    "^BT1_([0-9A-F]*)_([0-9A-F]*)_([0-9A-F]*)_([0-9A-F]*)$",
    name
  )
  parts <- regmatches(name, match)[[1L]]
  if(length(parts) != 5L){
    stop("Encoded parameter name is malformed.", call. = FALSE)
  }
  decoded <- vapply(parts[-1L], .bt_utf8_hex_decode, character(1))
  fields <- stats::setNames(
    as.list(unname(decoded)),
    c("kind", "formula_parameter", "term", "role")
  )
  .bt_validate_parameter_encoding_fields(fields)
  c(list(encoding_version = .bt_parameter_encoding_version), fields)
}

#' @rdname JAGS_fit_contract
JAGS_parameter_encoding_schema <- function(){

  data.frame(
    field = c("encoding_version", "kind", "formula_parameter", "term", "role"),
    type = c("integer", rep("character", 4L)),
    description = c(
      "BayesTools parameter-encoding schema version.",
      "Coordinate kind, such as fixed, random, or formula_output.",
      "Formula output parameter owning the coordinate.",
      "Exact formula term, design column, or random-effect block.",
      "Generated semantic role within the coordinate kind."
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname JAGS_fit_contract
JAGS_formula_name_map <- function(fit, parameter = NULL){

  check_char(parameter, "parameter", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  JAGS_validate_fit_contract(
    fit,
    requires = c("name_encoding", "formula_name_map", "formula_design")
  )
  designs <- attr(fit, "formula_design", exact = TRUE)
  if(is.null(designs)){
    return(if(is.null(parameter)) list() else NULL)
  }
  maps <- lapply(designs, function(design){
    map <- design$name_map
    .bt_validate_formula_name_map(map)
    map
  })
  if(is.null(parameter)){
    return(maps)
  }
  if(!parameter %in% names(maps)){
    stop("Formula design for parameter '", parameter, "' was not found.",
         call. = FALSE)
  }
  maps[[parameter]]
}

#' @rdname JAGS_fit_contract
JAGS_fit_contract <- function(fit){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  contract <- attr(fit, "fit_contract", exact = TRUE)
  .bt_validate_fit_contract_object(contract)
  contract
}

#' @rdname JAGS_fit_contract
JAGS_validate_fit_contract <- function(fit, requires = character()){

  check_char(requires, "requires", check_length = FALSE, allow_NA = FALSE)
  unknown <- setdiff(requires, .bt_fit_contract_components)
  if(length(unknown) > 0L){
    stop(
      "Unknown fitted-contract component",
      if(length(unknown) > 1L) "s: " else ": ",
      paste0("'", unknown, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  contract <- JAGS_fit_contract(fit)
  supported <- .bt_supported_fit_contract()
  for(component in requires){
    observed <- contract[[paste0(component, "_version")]]
    expected <- supported[[component]]
    if(length(observed) != 1L || is.na(observed) ||
       !identical(observed, expected)){
      stop(
        "The fitted object has missing or unsupported '", component,
        "' metadata. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
  }
  invisible(contract)
}

#' @rdname JAGS_fit_contract
JAGS_fit_contract_schema <- function(){

  data.frame(
    field = c("schema_version", paste0(.bt_fit_contract_components, "_version")),
    type = rep("integer", length(.bt_fit_contract_components) + 1L),
    description = c(
      "Fitted-object compatibility-profile schema version.",
      paste("Schema version for", gsub("_", " ", .bt_fit_contract_components), "metadata.")
    ),
    stringsAsFactors = FALSE
  )
}

.bt_validate_parameter_encoding_fields <- function(fields){

  if(!is.list(fields) || !identical(
    names(fields),
    c("kind", "formula_parameter", "term", "role")
  )){
    stop(
      "'fields' must be a named list containing, in order, 'kind', ",
      "'formula_parameter', 'term', and 'role'.",
      call. = FALSE
    )
  }
  valid <- vapply(fields, function(field){
    is.character(field) && length(field) == 1L && !is.na(field)
  }, logical(1))
  if(!all(valid)){
    stop("Every parameter-encoding field must be one non-missing character value.",
         call. = FALSE)
  }
  if(!nzchar(fields$kind) || !nzchar(fields$formula_parameter) ||
     !nzchar(fields$role)){
    stop("The 'kind', 'formula_parameter', and 'role' fields must not be empty.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.bt_utf8_hex_encode <- function(x){

  bytes <- as.integer(charToRaw(enc2utf8(x)))
  paste(sprintf("%02X", bytes), collapse = "")
}

.bt_utf8_hex_decode <- function(x){

  if(nchar(x) %% 2L != 0L){
    stop("Encoded parameter name contains malformed UTF-8 bytes.", call. = FALSE)
  }
  if(!nzchar(x)){
    return("")
  }
  starts <- seq.int(1L, nchar(x), by = 2L)
  values <- strtoi(substring(x, starts, starts + 1L), base = 16L)
  if(anyNA(values)){
    stop("Encoded parameter name contains malformed hexadecimal bytes.",
         call. = FALSE)
  }
  out <- rawToChar(as.raw(values))
  Encoding(out) <- "UTF-8"
  if(!identical(enc2utf8(out), out)){
    stop("Encoded parameter name contains invalid UTF-8.", call. = FALSE)
  }
  out
}

.bt_formula_name_map_empty <- function(){

  out <- data.frame(
    encoded_name = character(),
    jags_name = character(),
    kind = character(),
    formula_parameter = character(),
    term = character(),
    role = character(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("BayesTools_formula_name_map", "data.frame")
  attr(out, "schema_version") <- .bt_formula_name_map_version
  out
}

.bt_formula_name_map <- function(jags_name, kind, formula_parameter, term, role){

  fields <- Map(function(kind_i, parameter_i, term_i, role_i){
    list(
      kind = kind_i,
      formula_parameter = parameter_i,
      term = term_i,
      role = role_i
    )
  }, kind, formula_parameter, term, role)
  encoded_name <- vapply(fields, JAGS_parameter_encode, character(1))
  out <- data.frame(
    encoded_name = encoded_name,
    jags_name = jags_name,
    kind = kind,
    formula_parameter = formula_parameter,
    term = term,
    role = role,
    stringsAsFactors = FALSE
  )
  class(out) <- c("BayesTools_formula_name_map", "data.frame")
  attr(out, "schema_version") <- .bt_formula_name_map_version
  .bt_validate_formula_name_map(out)
  out
}

.bt_validate_formula_name_map <- function(map){

  required <- c(
    "encoded_name", "jags_name", "kind", "formula_parameter", "term", "role"
  )
  valid <- inherits(map, "BayesTools_formula_name_map") &&
    is.data.frame(map) &&
    identical(attr(map, "schema_version", exact = TRUE),
              .bt_formula_name_map_version) &&
    all(required %in% names(map))
  if(!valid){
    stop(
      "Formula name-map metadata are missing or unsupported. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  if(anyNA(map[, required, drop = FALSE]) ||
     any(!nzchar(map$encoded_name)) || any(!nzchar(map$jags_name)) ||
     anyDuplicated(map$encoded_name) || anyDuplicated(map$jags_name)){
    stop("Formula name-map metadata contain missing or duplicate names. Refit the model with this version of BayesTools.",
         call. = FALSE)
  }
  if(nrow(map) > 0L){
    for(i in seq_len(nrow(map))){
      .bt_check_jags_node_name(map$jags_name[i], "name_map$jags_name")
      decoded <- JAGS_parameter_decode(map$encoded_name[i])
      expected <- as.list(map[i, c("kind", "formula_parameter", "term", "role")])
      if(!identical(unname(decoded[names(expected)]), unname(expected))){
        stop("Formula name-map encoded and semantic fields disagree. Refit the model with this version of BayesTools.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.bt_supported_fit_contract <- function(){

  list(
    name_encoding = .bt_parameter_encoding_version,
    formula_name_map = .bt_formula_name_map_version,
    formula_design = .bt_formula_design_schema_version(),
    parameter_registry = .bt_parameter_registry_version,
    draw_geometry = if(exists(".bt_draw_geometry_version", inherits = TRUE)){
      .bt_draw_geometry_version
    }else{
      NA_integer_
    },
    parameter_catalog = if(exists(".bt_parameter_catalog_version", inherits = TRUE)){
      .bt_parameter_catalog_version
    }else{
      NA_integer_
    }
  )
}

.bt_attach_fit_contract <- function(fit){

  if(inherits(fit, "error")){
    return(fit)
  }
  supported <- .bt_supported_fit_contract()
  attr(fit, "fit_contract") <- c(
    list(schema_version = .bt_fit_contract_version),
    stats::setNames(supported, paste0(names(supported), "_version"))
  )
  fit
}

.bt_validate_fit_contract_object <- function(contract){

  required <- JAGS_fit_contract_schema()$field
  valid <- is.list(contract) && identical(names(contract), required) &&
    identical(contract$schema_version, .bt_fit_contract_version)
  if(!valid){
    stop(
      "The fitted object does not contain a supported compatibility profile. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  values <- contract[required]
  scalar_integer <- vapply(values, function(value){
    is.integer(value) && length(value) == 1L
  }, logical(1))
  if(!all(scalar_integer)){
    stop(
      "The fitted-object compatibility profile is malformed. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
