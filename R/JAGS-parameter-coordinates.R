# Concrete fitted-coordinate table used by the parameter map.

.bt_parameter_coordinates_columns <- c(
  "coordinate_name",
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
  "fixed_value",
  "display_label",
  "random_grouping",
  "random_structure",
  "internal"
)

.bt_parameter_coordinates_empty <- function(){

  out <- data.frame(
    coordinate_name = character(),
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
    fixed_value = numeric(),
    display_label = character(),
    random_grouping = character(),
    random_structure = character(),
    internal = logical(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("BayesTools_parameter_coordinates", "data.frame")
  out
}

.bt_validate_parameter_coordinates <- function(coordinates){

  if(!inherits(coordinates, "BayesTools_parameter_coordinates") ||
     !is.data.frame(coordinates)){
    stop(
      "The fitted parameter-coordinate table is malformed. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  missing <- setdiff(.bt_parameter_coordinates_columns, names(coordinates))
  if(length(missing) > 0L){
    stop(
      "The fitted parameter-coordinate table is missing required field",
      if(length(missing) > 1L) "s " else " ",
      paste0("'", missing, "'", collapse = ", "),
      ". Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(anyNA(coordinates$coordinate_name) ||
     any(!nzchar(coordinates$coordinate_name)) ||
     anyDuplicated(coordinates$coordinate_name)){
    stop(
      "The fitted parameter-coordinate table must contain unique, non-missing coordinate names. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(anyNA(coordinates$internal)){
    stop(
      "The fitted parameter-coordinate table contains an undefined 'internal' flag. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(!is.numeric(coordinates$fixed_value) ||
     any(!coordinates$monitor_status %in% c("sampled", "structural", "unavailable")) ||
     any(!is.na(coordinates$fixed_value[coordinates$monitor_status != "structural"])) ||
     any(!is.finite(coordinates$fixed_value[coordinates$monitor_status == "structural"]))){
    stop(
      "The fitted parameter-coordinate table contains malformed structural fixed values. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_parameter_coordinates_base <- function(x){

  sub("\\[.*$", "", x)
}

.bt_parameter_coordinates_index <- function(x){

  has_index <- grepl("\\[[^]]+\\]$", x)
  out <- rep("", length(x))
  out[has_index] <- sub("^.*\\[([^]]+)\\]$", "\\1", x[has_index])
  out
}

.bt_parameter_coordinates_dimensions <- function(base_name, columns){

  selected <- columns[.bt_parameter_coordinates_base(columns) == base_name]
  indices <- .bt_parameter_coordinates_index(selected)
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

.bt_parameter_coordinates_random_terms <- function(formula_design){

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(list())
  }
  .bt_random_effect_summary_random_terms(random_design)
}

.bt_parameter_coordinates_name_map <- function(formula_design){

  if(is.null(formula_design) || length(formula_design) == 0L){
    return(.bt_formula_name_map_empty())
  }
  has_map <- vapply(formula_design, function(design){
    !is.null(design$name_map)
  }, logical(1))
  if(!any(has_map)){
    return(.bt_formula_name_map_empty())
  }
  maps <- lapply(formula_design[has_map], function(design){
    map <- design$name_map
    .bt_validate_formula_name_map(map)
    map
  })
  out <- do.call(rbind, unname(maps))
  rownames(out) <- NULL
  class(out) <- c("BayesTools_formula_name_map", "data.frame")
  attr(out, "schema_version") <- .bt_formula_name_map_version
  if(anyDuplicated(out$jags_name)){
    stop(
      "Formula name maps contain a duplicate JAGS base name. Refit the model after resolving the generated-name collision.",
      call. = FALSE
    )
  }
  .bt_validate_formula_name_map(out)
  out
}

.bt_parameter_coordinates_random_family <- function(random_term){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || is.na(stem) || !nzchar(stem)){
    return(data.frame(base = character(), role = character()))
  }
  family <- data.frame(
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

  correlation <- random_term$correlation
  if(is.list(correlation) && identical(correlation$type, "lkj")){
    coordinate_names <- unique(c(
      correlation$primitive_names,
      correlation$cpc_names
    ))
    coordinate_bases <- unique(.bt_parameter_coordinates_base(
      coordinate_names
    ))
    coordinate_bases <- coordinate_bases[
      !is.na(coordinate_bases) & nzchar(coordinate_bases)
    ]
    if(length(coordinate_bases) > 0L){
      family <- rbind(
        family,
        data.frame(
          base = coordinate_bases,
          role = rep("random_correlation_coordinate", length(coordinate_bases)),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  family
}

.bt_parameter_coordinates_prior_owner <- function(base_name, prior_list){

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

.bt_parameter_coordinates_point_values <- function(parameter, prior){

  if(is.null(prior) || !is.prior.point(prior)){
    return(stats::setNames(numeric(), character()))
  }
  parameter_names <- if(is.prior.vector(prior) || is.prior.factor(prior)){
    .JAGS_prior_factor_names(parameter, prior)
  }else{
    parameter
  }
  values <- rep(
    prior[["parameters"]][["location"]],
    length.out = length(parameter_names)
  )
  stats::setNames(as.numeric(values), parameter_names)
}

.bt_parameter_coordinates_term_owner <- function(base_name, random_terms){

  if(length(random_terms) == 0L){
    return(NULL)
  }

  # Resolve concrete random-effect parameters before generated auxiliary
  # suffixes. This preserves exact ownership when a legitimate SD parameter
  # itself ends in "_variable", "_indicator", or "_inclusion".
  for(random_term in random_terms){
    family <- .bt_parameter_coordinates_random_family(random_term)
    match <- match(base_name, family$base)
    if(!is.na(match)){
      return(list(
        random_term = random_term,
        role = family$role[match]
      ))
    }
    sd_names <- .bt_parameter_coordinates_random_sd_names(random_term)
    if(base_name %in% sd_names){
      return(list(random_term = random_term, role = "random_sd"))
    }
  }

  auxiliary_roles <- c(
    "_indicator" = "random_inclusion_indicator",
    "_inclusion" = "random_inclusion_probability",
    "_variable" = "random_sd_variable"
  )
  for(suffix in names(auxiliary_roles)){
    if(!endsWith(base_name, suffix)){
      next
    }
    owner_name <- substr(base_name, 1L, nchar(base_name) - nchar(suffix))
    for(random_term in random_terms){
      sd_names <- .bt_parameter_coordinates_random_sd_names(random_term)
      if(owner_name %in% sd_names){
        return(list(
          random_term = random_term,
          role = unname(auxiliary_roles[[suffix]])
        ))
      }
    }
  }

  NULL
}

.bt_parameter_coordinates_random_sd_names <- function(random_term){

  binding <- random_term$sd_binding
  if(!is.null(binding)){
    .bt_check_random_sd_binding(binding)
    if(.bt_random_sd_binding_has_external_source(binding) &&
       !isTRUE(binding$true_allocation)){
      return(character())
    }
  }

  sd_names <- unique(.bt_parameter_coordinates_base(
    random_term$sd_parameter_names
  ))
  sd_names[!is.na(sd_names) & nzchar(sd_names)]
}

.bt_parameter_coordinates_prior_role <- function(prior){

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
  if(isTRUE(metadata$raw_sd) || isTRUE(metadata$raw_allocation_sd)){
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

.bt_parameter_coordinates_formula_parameter <- function(prior, random_term,
                                                      name_map_row = NULL){

  if(!is.null(name_map_row) && nrow(name_map_row) == 1L){
    return(name_map_row$formula_parameter)
  }

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

.bt_parameter_coordinates_coordinate <- function(coordinate_name, random_term,
                                              role, name_map_row = NULL){

  out <- list(term = "", column = "")
  if(is.null(random_term)){
    if(identical(role, "fixed_coefficient")){
      out$term <- if(!is.null(name_map_row) && nrow(name_map_row) == 1L){
        name_map_row$term
      }else{
        .bt_parameter_coordinates_base(coordinate_name)
      }
      out$column <- coordinate_name
    }
    return(out)
  }
  index <- .bt_parameter_coordinates_index(coordinate_name)
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
  if(role %in% c("random_sd", "random_sd_variable")){
    sd_names <- random_term$sd_parameter_names
    sd_name <- coordinate_name
    if(identical(role, "random_sd_variable")){
      sd_name <- sub("_variable(?=\\[|$)", "", sd_name, perl = TRUE)
    }
    sd_match <- match(sd_name, sd_names)
    if(is.na(sd_match)){
      sd_base <- .bt_parameter_coordinates_base(sd_name)
      base_matches <- which(
        .bt_parameter_coordinates_base(sd_names) == sd_base
      )
      if(length(base_matches) == 1L){
        sd_match <- base_matches
      }
    }
    if(!is.na(sd_match) && sd_match <= length(column_names)){
      out$column <- column_names[sd_match]
    }
  }
  if(nzchar(out$column)){
    out$term <- out$column
  }
  out
}

.bt_parameter_coordinates_scale <- function(role, formula_parameter,
                                         formula_scale){

  if(identical(role, "random_latent")){
    return("unit_latent")
  }
  if(identical(role, "random_group_coefficient")){
    return("fitted_standardized")
  }
  if(role %in% c("random_sd", "random_sd_variable",
                 "random_correlation")){
    return("fitted_covariance")
  }
  if(role %in% c("allocation", "random_inclusion_indicator",
                 "random_inclusion_probability")){
    return("unitless")
  }
  if(identical(role, "random_correlation_coordinate")){
    return("unitless")
  }
  if(nzchar(formula_parameter) &&
     !is.null(formula_scale) &&
     !is.null(formula_scale[[formula_parameter]])){
    return("fitted_standardized")
  }
  "fitted_original"
}

.bt_parameter_coordinates_display <- function(coordinate_name, prior, random_term,
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
    return(format_parameter_names(
      coordinate_name,
      formula_parameters = if(nzchar(formula_parameter)){
        formula_parameter
      }else{
        NULL
      },
      formula_prefix = TRUE
    ))
  }

  names <- coordinate_name
  raw_names <- coordinate_name
  prefix <- .bt_random_effect_summary_formula_prefix(formula_parameter, TRUE)
  if(identical(role, "random_sd")){
    names <- .bt_random_effect_summary_raw_sd_display_names(
      names = names,
      raw_names = raw_names,
      prior_list = if(is.null(prior)) list() else stats::setNames(
        list(prior),
        .bt_parameter_coordinates_base(coordinate_name)
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

.bt_build_parameter_coordinates <- function(columns, monitor_names = columns,
                                         prior_list = NULL,
                                         formula_design = NULL,
                                         formula_scale = NULL,
                                         backend_anchor = NULL){

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
          if(is.prior.factor(prior) || is.prior.vector(prior)){
            .JAGS_prior_factor_names(parameter, prior)
          }else{
            parameter
          }
        )
      }
    }
  }
  coordinate_names <- unique(c(columns, setdiff(structural, columns)))
  if(length(coordinate_names) == 0L){
    return(.bt_parameter_coordinates_empty())
  }

  random_terms <- .bt_parameter_coordinates_random_terms(formula_design)
  name_map <- .bt_parameter_coordinates_name_map(formula_design)
  allocation_indicators <-
    .bt_random_variance_allocation_inclusion_indicator_names(formula_design)
  dirichlet_auxiliaries <- vapply(
    names(prior_list)[vapply(prior_list, is.prior.simplex, logical(1))],
    .JAGS_prior_dirichlet_eta_name,
    character(1)
  )
  bases <- .bt_parameter_coordinates_base(coordinate_names)
  coordinates <- .bt_parameter_coordinates_empty()
  coordinates <- coordinates[rep(NA_integer_, length(coordinate_names)), , drop = FALSE]

  for(i in seq_along(coordinate_names)){
    coordinate_name <- coordinate_names[i]
    base_name <- bases[i]
    name_map_row <- name_map[name_map$jags_name == base_name, , drop = FALSE]
    prior <- .bt_parameter_coordinates_prior_owner(base_name, prior_list)
    prior_metadata <- .bt_random_effect_metadata(prior)
    owner <- .bt_parameter_coordinates_term_owner(base_name, random_terms)
    random_term <- if(is.null(owner)) NULL else owner$random_term
    role <- if(is.null(owner)){
      .bt_parameter_coordinates_prior_role(prior)
    }else{
      owner$role
    }
    if(coordinate_name %in% allocation_indicators){
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

    formula_parameter <- .bt_parameter_coordinates_formula_parameter(
      prior,
      random_term,
      name_map_row
    )
    coordinate <- .bt_parameter_coordinates_coordinate(
      coordinate_name,
      random_term,
      role,
      name_map_row
    )
    if(identical(role, "fixed_coefficient") && nzchar(formula_parameter) &&
       nrow(name_map_row) == 0L){
      formula_prefix <- paste0(formula_parameter, "_")
      if(startsWith(coordinate$term, formula_prefix)){
        coordinate$term <- substring(
          coordinate$term,
          nchar(formula_prefix) + 1L
        )
      }
    }
    requested <- monitor_names[
      monitor_names == coordinate_name |
        monitor_names == base_name
    ]
    monitor_name <- if(length(requested) > 0L) requested[1L] else base_name
    monitor_status <- if(!is.null(prior) && is.prior.point(prior)){
      "structural"
    }else if(coordinate_name %in% columns){
      "sampled"
    }else{
      "structural"
    }
    fixed_value <- NA_real_
    if(identical(monitor_status, "structural") && !is.null(prior)){
      prior_values <- .bt_parameter_coordinates_point_values(base_name, prior)
      value_match <- match(coordinate_name, names(prior_values))
      if(!is.na(value_match)){
        fixed_value <- unname(prior_values[value_match])
      }
    }
    if(!is.null(backend_anchor) && identical(coordinate_name, backend_anchor)){
      role <- "backend_anchor"
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

    coordinates[i, ] <- list(
      coordinate_name,
      monitor_name,
      formula_parameter,
      role,
      random_block,
      random_name,
      coordinate$term,
      coordinate$column,
      .bt_parameter_coordinates_index(coordinate_name),
      .bt_parameter_coordinates_dimensions(base_name, columns),
      .bt_parameter_coordinates_scale(role, formula_parameter, formula_scale),
      monitor_status,
      fixed_value,
      .bt_parameter_coordinates_display(
        coordinate_name,
        prior,
        random_term,
        role,
        formula_parameter
      ),
      grouping,
      structure,
      role %in% c(
        "backend_anchor",
        "random_latent",
        "random_group_coefficient",
        "random_correlation",
        "random_correlation_coordinate",
        "random_inclusion_indicator",
        "random_inclusion_probability",
        "random_sd_variable"
      ) ||
        base_name %in% dirichlet_auxiliaries ||
        isTRUE(prior_metadata$allocation)
    )
  }

  rownames(coordinates) <- NULL
  class(coordinates) <- c("BayesTools_parameter_coordinates", "data.frame")
  .bt_validate_parameter_coordinates(coordinates)
  coordinates
}

