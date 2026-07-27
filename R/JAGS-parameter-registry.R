# Canonical fitted-parameter registry.

.bt_parameter_registry_version <- 1L

.bt_parameter_registry_columns <- c(
  "canonical_name",
  "monitor_name",
  "formula_parameter",
  "role",
  "random_block",
  "random_name",
  "term",
  "column",
  "index",
  "dimensions",
  "fitted_scale",
  "monitor_status",
  "display_label",
  "random_grouping",
  "random_structure",
  "internal"
)

.bt_parameter_registry_empty <- function(){

  out <- data.frame(
    canonical_name = character(),
    monitor_name = character(),
    formula_parameter = character(),
    role = character(),
    random_block = character(),
    random_name = character(),
    term = character(),
    column = character(),
    index = character(),
    dimensions = character(),
    fitted_scale = character(),
    monitor_status = character(),
    display_label = character(),
    random_grouping = character(),
    random_structure = character(),
    internal = logical(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("BayesTools_parameter_registry", "data.frame")
  attr(out, "schema_version") <- .bt_parameter_registry_version
  out
}

#' Canonical JAGS parameter registry
#'
#' @description
#' `JAGS_parameter_registry()` returns the versioned parameter registry stored
#' by [JAGS_fit()]. The registry is the authoritative mapping between concrete
#' posterior columns and their formula, random-effect block, semantic role,
#' fitted coordinate scale, and display label. Downstream packages should use
#' this accessor instead of parsing JAGS parameter names.
#'
#' `JAGS_parameter_registry_schema()` documents the stable fields in the
#' current registry schema. Registry rows are unique by `canonical_name`.
#' Matrix/vector indices are stored as comma-separated text in `index`, and
#' their fitted dimensions are stored in `dimensions` (for example, `"3x2"`).
#' Empty strings mean that a field does not apply.
#'
#' Raw random-effect latent variables and realized group coefficients have
#' `internal = TRUE`. Their `fitted_scale` is respectively `"unit_latent"` and
#' `"fitted_standardized"`; they must not be presented as original-scale
#' coefficients.
#'
#' @param fit a fitted object created by [JAGS_fit()].
#'
#' @return `JAGS_parameter_registry()` returns a
#' `BayesTools_parameter_registry` data frame.
#' `JAGS_parameter_registry_schema()` returns a data frame describing each
#' registry field.
#'
#' @export
#' @name JAGS_parameter_registry
NULL

#' @rdname JAGS_parameter_registry
#' @export
JAGS_parameter_registry <- function(fit){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }

  registry <- attr(fit, "parameter_registry", exact = TRUE)
  if(is.null(registry)){
    stop(
      "The fitted object does not contain the canonical parameter registry. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }

  .bt_validate_parameter_registry(registry)
  registry
}

#' @rdname JAGS_parameter_registry
#' @export
JAGS_parameter_registry_schema <- function(){

  data.frame(
    field = .bt_parameter_registry_columns,
    type = c(
      rep("character", 15L),
      "logical"
    ),
    description = c(
      "Unique concrete posterior or structural parameter name.",
      "JAGS monitor node requested during fitting.",
      "Formula output parameter owning the coefficient or random block.",
      "Semantic role such as fixed coefficient, random SD, random correlation, latent variable, group coefficient, allocation, derived summary, or ordinary parameter.",
      "Canonical random-effect block identifier.",
      "Public random-effect name, distinct from the grouping label.",
      "Formula term represented by the parameter coordinate.",
      "Concrete fixed/random design-matrix column.",
      "Comma-separated concrete JAGS array indices.",
      "Fitted array dimensions joined by 'x'.",
      "Coordinate scale used by the fitted monitor.",
      "Whether the row is sampled, structural, or requested but unavailable.",
      "Default unambiguous user-facing label.",
      "Random grouping-variable label.",
      "Random covariance structure.",
      "Whether the coordinate is implementation-level and excluded from public coefficient tables."
    ),
    stringsAsFactors = FALSE
  )
}

.bt_validate_parameter_registry <- function(registry){

  if(!inherits(registry, "BayesTools_parameter_registry") ||
     !is.data.frame(registry)){
    stop(
      "The fitted parameter registry is malformed. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  version <- attr(registry, "schema_version", exact = TRUE)
  if(!identical(version, .bt_parameter_registry_version)){
    stop(
      "The fitted parameter registry uses unsupported schema version '",
      if(is.null(version)) "missing" else as.character(version),
      "'. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  missing <- setdiff(.bt_parameter_registry_columns, names(registry))
  if(length(missing) > 0L){
    stop(
      "The fitted parameter registry is missing required field",
      if(length(missing) > 1L) "s " else " ",
      paste0("'", missing, "'", collapse = ", "),
      ". Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(anyNA(registry$canonical_name) ||
     any(!nzchar(registry$canonical_name)) ||
     anyDuplicated(registry$canonical_name)){
    stop(
      "The fitted parameter registry must contain unique, non-missing canonical names. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(anyNA(registry$internal)){
    stop(
      "The fitted parameter registry contains an undefined 'internal' flag. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_parameter_registry_base <- function(x){

  sub("\\[.*$", "", x)
}

.bt_parameter_registry_index <- function(x){

  has_index <- grepl("\\[[^]]+\\]$", x)
  out <- rep("", length(x))
  out[has_index] <- sub("^.*\\[([^]]+)\\]$", "\\1", x[has_index])
  out
}

.bt_parameter_registry_dimensions <- function(base_name, columns){

  selected <- columns[.bt_parameter_registry_base(columns) == base_name]
  indices <- .bt_parameter_registry_index(selected)
  indices <- indices[nzchar(indices)]
  if(length(indices) == 0L){
    return("")
  }
  parsed <- strsplit(indices, ",", fixed = TRUE)
  n_dimensions <- unique(lengths(parsed))
  if(length(n_dimensions) != 1L){
    return("")
  }
  values <- suppressWarnings(lapply(parsed, as.integer))
  if(any(vapply(values, anyNA, logical(1)))){
    return("")
  }
  maxima <- vapply(seq_len(n_dimensions), function(i){
    max(vapply(values, `[[`, integer(1), i))
  }, integer(1))
  paste(maxima, collapse = "x")
}

.bt_parameter_registry_random_terms <- function(formula_design){

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(list())
  }
  .bt_random_effect_summary_random_terms(random_design)
}

.bt_parameter_registry_random_family <- function(random_term){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || is.na(stem) || !nzchar(stem)){
    return(data.frame(base = character(), role = character()))
  }
  data.frame(
    base = paste0(
      stem,
      c(
        "_xRE_Zx",
        "_xRE_GROUP_Zx",
        "_xRE_UNIT_COEFx",
        "_xRE_COEFx",
        "_xRE_CORx_R",
        "_xRE_CORx_L",
        "_rho",
        "_rho_z",
        "_rho_logit"
      )
    ),
    role = c(
      "random_latent",
      "random_latent",
      "random_latent",
      "random_group_coefficient",
      "random_correlation",
      "random_correlation",
      "random_correlation",
      "random_correlation",
      "random_correlation"
    ),
    stringsAsFactors = FALSE
  )
}

.bt_parameter_registry_prior_owner <- function(base_name, prior_list){

  if(length(prior_list) == 0L || is.null(names(prior_list))){
    return(NULL)
  }
  match <- which(names(prior_list) == base_name)
  if(length(match) == 0L){
    eta_names <- vapply(
      names(prior_list),
      .JAGS_prior_dirichlet_eta_name,
      character(1)
    )
    match <- which(eta_names == base_name)
  }
  if(length(match) != 1L){
    return(NULL)
  }
  prior_list[[match]]
}

.bt_parameter_registry_term_owner <- function(base_name, random_terms){

  if(length(random_terms) == 0L){
    return(NULL)
  }
  for(random_term in random_terms){
    family <- .bt_parameter_registry_random_family(random_term)
    match <- match(base_name, family$base)
    if(!is.na(match)){
      return(list(
        random_term = random_term,
        role = family$role[match]
      ))
    }
    sd_names <- unique(random_term$sd_parameter_names)
    sd_names <- sd_names[!is.na(sd_names)]
    if(base_name %in% sd_names){
      return(list(random_term = random_term, role = "random_sd"))
    }
  }
  NULL
}

.bt_parameter_registry_prior_role <- function(prior){

  if(is.null(prior)){
    return("parameter")
  }
  metadata <- .bt_random_effect_metadata(prior)
  if(nzchar(metadata$summary)){
    return("derived_summary")
  }
  if(isTRUE(metadata$allocation)){
    return("allocation")
  }
  if(isTRUE(metadata$raw_sd) || isTRUE(metadata$raw_sd_total)){
    return("random_sd")
  }
  if(isTRUE(metadata$raw_correlation)){
    return("random_correlation")
  }
  if(!is.null(attr(prior, "parameter", exact = TRUE))){
    return("fixed_coefficient")
  }
  "parameter"
}

.bt_parameter_registry_formula_parameter <- function(prior, random_term){

  if(!is.null(random_term) &&
     !is.null(random_term$parameter) &&
     length(random_term$parameter) == 1L){
    return(as.character(random_term$parameter))
  }
  if(!is.null(random_term) &&
     !is.null(random_term$parameter_stem)){
    stem <- as.character(random_term$parameter_stem)
    marker <- regexpr("__xREx__", stem, fixed = TRUE)
    if(marker[1L] > 1L){
      return(substr(stem, 1L, marker[1L] - 1L))
    }
  }
  parameter <- attr(prior, "parameter", exact = TRUE)
  if(is.null(parameter) || length(parameter) != 1L || is.na(parameter)){
    return("")
  }
  as.character(parameter)
}

.bt_parameter_registry_coordinate <- function(canonical_name, random_term,
                                              role){

  out <- list(term = "", column = "")
  if(is.null(random_term)){
    if(identical(role, "fixed_coefficient")){
      out$term <- .bt_parameter_registry_base(canonical_name)
      out$column <- canonical_name
    }
    return(out)
  }
  index <- .bt_parameter_registry_index(canonical_name)
  index <- if(nzchar(index)) suppressWarnings(as.integer(
    strsplit(index, ",", fixed = TRUE)[[1L]]
  )) else integer()

  column_names <- random_term$column_names
  if(!is.character(column_names)){
    column_names <- character()
  }
  if(role %in% c("random_latent", "random_group_coefficient",
                 "random_correlation") &&
     length(index) >= 2L && !is.na(index[2L]) &&
     index[2L] >= 1L && index[2L] <= length(column_names)){
    out$column <- column_names[index[2L]]
  }
  if(identical(role, "random_sd")){
    sd_names <- random_term$sd_parameter_names
    sd_match <- match(.bt_parameter_registry_base(canonical_name), sd_names)
    if(!is.na(sd_match) && sd_match <= length(column_names)){
      out$column <- column_names[sd_match]
    }
  }
  if(nzchar(out$column)){
    out$term <- out$column
  }
  out
}

.bt_parameter_registry_scale <- function(role, formula_parameter,
                                         formula_scale){

  if(identical(role, "random_latent")){
    return("unit_latent")
  }
  if(identical(role, "random_group_coefficient")){
    return("fitted_standardized")
  }
  if(role %in% c("random_sd", "random_correlation")){
    return("fitted_covariance")
  }
  if(identical(role, "allocation")){
    return("unitless")
  }
  if(nzchar(formula_parameter) &&
     !is.null(formula_scale) &&
     !is.null(formula_scale[[formula_parameter]])){
    return("fitted_standardized")
  }
  "fitted_original"
}

.bt_parameter_registry_display <- function(canonical_name, prior, random_term,
                                           role, formula_parameter){

  if(!is.null(prior)){
    label <- attr(prior, "random_summary_label", exact = TRUE)
    if(!is.null(label) && length(label) == 1L && !is.na(label) && nzchar(label)){
      prefix <- .bt_random_effect_summary_formula_prefix(
        formula_parameter,
        TRUE
      )
      return(paste0(prefix, label))
    }
  }
  if(is.null(random_term)){
    if(nzchar(formula_parameter)){
      prefix <- paste0(formula_parameter, "_")
      if(startsWith(canonical_name, prefix)){
        return(paste0(
          "(", formula_parameter, ") ",
          substring(canonical_name, nchar(prefix) + 1L)
        ))
      }
    }
    return(canonical_name)
  }

  names <- canonical_name
  raw_names <- canonical_name
  prefix <- .bt_random_effect_summary_formula_prefix(formula_parameter, TRUE)
  if(identical(role, "random_sd")){
    names <- .bt_random_effect_summary_raw_sd_display_names(
      names = names,
      raw_names = raw_names,
      prior_list = if(is.null(prior)) list() else stats::setNames(
        list(prior),
        .bt_parameter_registry_base(canonical_name)
      ),
      random_term = random_term,
      prefix = prefix
    )
  }else if(identical(role, "random_correlation")){
    names <- .bt_random_effect_summary_raw_rho_display_names(
      names = names,
      raw_names = raw_names,
      random_term = random_term,
      prefix = prefix
    )
    names <- .bt_random_effect_summary_raw_matrix_display_names(
      names = names,
      raw_names = raw_names,
      random_term = random_term,
      prefix = prefix
    )
  }else if(role %in% c("random_latent", "random_group_coefficient")){
    names <- .bt_random_effect_summary_raw_matrix_display_names(
      names = names,
      raw_names = raw_names,
      random_term = random_term,
      prefix = prefix
    )
  }
  names
}

.bt_build_parameter_registry <- function(columns, monitor_names = columns,
                                         prior_list = NULL,
                                         formula_design = NULL,
                                         formula_scale = NULL){

  if(is.null(columns)){
    columns <- character()
  }
  columns <- unique(as.character(columns))
  monitor_names <- unique(as.character(monitor_names))
  prior_list <- if(is.null(prior_list)) list() else prior_list

  structural <- character()
  if(length(prior_list) > 0L){
    for(parameter in names(prior_list)){
      prior <- prior_list[[parameter]]
      if(is.prior.point(prior)){
        structural <- c(
          structural,
          if(is.prior.factor(prior)){
            .JAGS_prior_factor_names(parameter, prior)
          }else{
            parameter
          }
        )
      }
    }
  }
  canonical_names <- unique(c(columns, setdiff(structural, columns)))
  if(length(canonical_names) == 0L){
    return(.bt_parameter_registry_empty())
  }

  random_terms <- .bt_parameter_registry_random_terms(formula_design)
  allocation_indicators <-
    .bt_random_variance_allocation_inclusion_indicator_names(formula_design)
  bases <- .bt_parameter_registry_base(canonical_names)
  registry <- .bt_parameter_registry_empty()
  registry <- registry[rep(NA_integer_, length(canonical_names)), , drop = FALSE]

  for(i in seq_along(canonical_names)){
    canonical_name <- canonical_names[i]
    base_name <- bases[i]
    prior <- .bt_parameter_registry_prior_owner(base_name, prior_list)
    owner <- .bt_parameter_registry_term_owner(base_name, random_terms)
    random_term <- if(is.null(owner)) NULL else owner$random_term
    role <- if(is.null(owner)){
      .bt_parameter_registry_prior_role(prior)
    }else{
      owner$role
    }
    if(canonical_name %in% allocation_indicators){
      role <- "allocation"
    }
    if(is.null(random_term) && !is.null(prior)){
      block <- .bt_random_effect_prior_effect(prior)
      if(nzchar(block)){
        matches <- vapply(random_terms, function(term){
          block %in% c(
            term$block_name,
            .bt_random_effect_public_name(term),
            term$group_label
          )
        }, logical(1))
        if(sum(matches) == 1L){
          random_term <- random_terms[[which(matches)]]
        }
      }
    }

    formula_parameter <- .bt_parameter_registry_formula_parameter(
      prior,
      random_term
    )
    coordinate <- .bt_parameter_registry_coordinate(
      canonical_name,
      random_term,
      role
    )
    if(identical(role, "fixed_coefficient") && nzchar(formula_parameter)){
      formula_prefix <- paste0(formula_parameter, "_")
      if(startsWith(coordinate$term, formula_prefix)){
        coordinate$term <- substring(
          coordinate$term,
          nchar(formula_prefix) + 1L
        )
      }
    }
    requested <- monitor_names[
      monitor_names == canonical_name |
        monitor_names == base_name
    ]
    monitor_name <- if(length(requested) > 0L) requested[1L] else base_name
    monitor_status <- if(canonical_name %in% columns){
      "sampled"
    }else{
      "structural"
    }
    random_block <- if(is.null(random_term)){
      if(is.null(prior)) "" else .bt_random_effect_prior_effect(prior)
    }else{
      as.character(random_term$block_name)
    }
    random_name <- if(is.null(random_term)){
      if(is.null(prior)) "" else .bt_random_effect_prior_name(prior)
    }else{
      .bt_random_effect_public_name(random_term)
    }
    grouping <- if(is.null(random_term)){
      if(is.null(prior)) "" else .bt_random_effect_prior_grouping(prior)
    }else{
      .bt_random_effect_summary_group_label(random_term)
    }
    structure <- if(is.null(random_term)){
      if(is.null(prior)) "" else .bt_random_effect_prior_structure(prior)
    }else{
      .bt_random_effect_summary_term_structure(random_term)
    }

    registry[i, ] <- list(
      canonical_name,
      monitor_name,
      formula_parameter,
      role,
      random_block,
      random_name,
      coordinate$term,
      coordinate$column,
      .bt_parameter_registry_index(canonical_name),
      .bt_parameter_registry_dimensions(base_name, columns),
      .bt_parameter_registry_scale(role, formula_parameter, formula_scale),
      monitor_status,
      .bt_parameter_registry_display(
        canonical_name,
        prior,
        random_term,
        role,
        formula_parameter
      ),
      grouping,
      structure,
      role %in% c(
        "random_latent",
        "random_group_coefficient",
        "random_correlation"
      )
    )
  }

  rownames(registry) <- NULL
  class(registry) <- c("BayesTools_parameter_registry", "data.frame")
  attr(registry, "schema_version") <- .bt_parameter_registry_version
  .bt_validate_parameter_registry(registry)
  registry
}

.bt_attach_parameter_registry <- function(fit, monitor_names = NULL){

  if(inherits(fit, "error")){
    return(fit)
  }
  samples <- .extract_posterior_samples(fit, as_list = FALSE)
  columns <- colnames(samples)
  if(is.null(monitor_names)){
    monitor_names <- if(is.list(fit) && !is.null(fit$monitor)){
      fit$monitor
    }else{
      columns
    }
  }
  attr(fit, "parameter_registry") <- .bt_build_parameter_registry(
    columns = columns,
    monitor_names = monitor_names,
    prior_list = attr(fit, "prior_list", exact = TRUE),
    formula_design = attr(fit, "formula_design", exact = TRUE),
    formula_scale = attr(fit, "formula_scale", exact = TRUE)
  )
  fit
}

.bt_parameter_registry_rows <- function(registry, canonical_names){

  .bt_validate_parameter_registry(registry)
  match <- match(canonical_names, registry$canonical_name)
  out <- registry[match, , drop = FALSE]
  rownames(out) <- NULL
  out
}
