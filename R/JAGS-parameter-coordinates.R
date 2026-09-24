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
  "internal",
  "convergence_role"
)

# Convergence roles of fitted coordinates, set at fit time from declared
# metadata (see .bt_parameter_coordinates_convergence_roles()).
.bt_convergence_roles <- c(
  "sampled",
  "indicator",
  "structural",
  "derived",
  "auxiliary"
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
    convergence_role = character(),
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
  if(!identical(names(coordinates), .bt_parameter_coordinates_columns)){
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
    stop(
      "The fitted parameter-coordinate table is malformed. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  character_columns <- setdiff(
    .bt_parameter_coordinates_columns,
    c("fixed_value", "internal")
  )
  if(!all(vapply(coordinates[character_columns], is.character, logical(1))) ||
     !is.numeric(coordinates$fixed_value) ||
     !is.logical(coordinates$internal) ||
     anyNA(coordinates[character_columns]) ||
     anyNA(coordinates$internal) ||
     any(!nzchar(coordinates$coordinate_name)) ||
     anyDuplicated(coordinates$coordinate_name)){
    stop(
      "The fitted parameter-coordinate table must contain unique, non-missing coordinate names. ",
      "Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(any(!coordinates$monitor_status %in% c("sampled", "structural", "unavailable")) ||
     any(!is.na(coordinates$fixed_value[coordinates$monitor_status != "structural"])) ||
     any(!is.finite(coordinates$fixed_value[coordinates$monitor_status == "structural"]))){
    stop(
      "The fitted parameter-coordinate table contains malformed structural fixed values. Refit the model with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(any(!coordinates$convergence_role %in% .bt_convergence_roles)){
    stop(
      "The fitted parameter-coordinate table contains unknown convergence roles. Refit the model with the current BayesTools version.",
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

.bt_parameter_coordinates_dimensions <- function(columns){

  indices <- .bt_parameter_coordinates_index(columns)
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
  maps <- lapply(formula_design, function(design){
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
        "_xRE_MEANx",
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
      "random_mean_coordinate",
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
    auxiliary_owners <-
      .bt_parameter_coordinates_random_prior_auxiliary_owners(prior_list)
    if(base_name %in% names(auxiliary_owners)){
      match <- match(auxiliary_owners[[base_name]], names(prior_list))
    }
  }
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

.bt_parameter_coordinates_random_prior_auxiliary_owners <- function(prior_list){

  if(length(prior_list) == 0L || is.null(names(prior_list))){
    return(stats::setNames(character(), character()))
  }
  random_prior_names <- names(prior_list)[vapply(prior_list, function(prior){
    isTRUE(attr(prior, "random_allocation_sd", exact = TRUE))
  }, logical(1))]
  owners <- unlist(lapply(random_prior_names, function(parameter){
    stats::setNames(
      rep(parameter, 3L),
      paste0(parameter, c("_indicator", "_inclusion", "_variable"))
    )
  }), use.names = TRUE)

  owners
}

# Private implementation nodes of ordered-prior Dirichlet allocations and of
# weight functions, matching the nodes posterior extraction treats as
# auxiliary. They remain coordinate-only dependencies.
.bt_parameter_coordinates_private_auxiliaries <- function(prior_list){

  if(length(prior_list) == 0L){
    return(character())
  }
  unique(unlist(lapply(prior_list, function(prior){
    if(is.prior.ordered(prior)){
      return(vapply(
        .prior_ordered_dirichlet_records(prior),
        function(record) .JAGS_prior_dirichlet_eta_name(record$node),
        character(1)
      ))
    }
    if(is.prior.weightfunction(prior)){
      private <- .JAGS_monitor_private.weightfunction(prior)
      if(identical(prior$weights$type, "cumulative")){
        private <- c(private, "eta", "omega_ratio")
      }
      return(private)
    }
    character()
  }), use.names = FALSE))
}

# A point prior whose location is an expression is a deterministic function
# of other nodes, so its coordinate is derived rather than a structural
# constant and has no fixed value.
.bt_parameter_coordinates_structural_point <- function(prior){

  !is.null(prior) && is.prior.point(prior) && !.is_prior_expression(prior)
}

.bt_parameter_coordinates_point_values <- function(parameter, prior){

  if(!.bt_parameter_coordinates_structural_point(prior)){
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
  if(role %in% c("random_latent", "random_group_coefficient", "random_mean_coordinate",
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
  if(role %in% c("random_group_coefficient", "random_mean_coordinate")){
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

# Display name of one fixed factor coordinate: its level cell when the
# coordinate structurally is one, otherwise contrast coefficient `{j}`. The
# per-prior names are computed once per coordinate table through `cache`.
.bt_parameter_coordinates_factor_display_name <- function(coordinate_name,
                                                          prior, cache){

  if(is.null(prior)){
    return(coordinate_name)
  }
  base_name <- .bt_parameter_coordinates_base(coordinate_name)
  if(!exists(base_name, envir = cache, inherits = FALSE)){
    assign(
      base_name,
      .bt_factor_coordinate_display_names(base_name, prior),
      envir = cache
    )
  }
  display_names <- get(base_name, envir = cache, inherits = FALSE)
  if(is.null(display_names)){
    return(coordinate_name)
  }
  index <- .bt_parameter_coordinates_index(coordinate_name)
  coefficient <- if(nzchar(index)){
    suppressWarnings(as.integer(index))
  }else{
    1L
  }
  if(is.na(coefficient) || coefficient < 1L ||
     coefficient > length(display_names)){
    return(coordinate_name)
  }
  display_names[[coefficient]]
}

.bt_parameter_coordinates_display <- function(coordinate_name, prior, random_term,
                                           role, formula_parameter,
                                           factor_cache = new.env(parent = emptyenv())){

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
    if(role %in% c("fixed_coefficient", "parameter")){
      coordinate_name <- .bt_parameter_coordinates_factor_display_name(
        coordinate_name,
        prior,
        factor_cache
      )
    }
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
                                         backend_anchor = NULL,
                                         add_parameters = NULL,
                                         model_syntax = NULL,
                                         data_names = NULL){

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
      if(.bt_parameter_coordinates_structural_point(prior)){
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
  random_prior_auxiliaries <- names(
    .bt_parameter_coordinates_random_prior_auxiliary_owners(prior_list)
  )
  dirichlet_auxiliaries <- c(
    vapply(
      names(prior_list)[vapply(prior_list, is.prior.simplex, logical(1))],
      .JAGS_prior_dirichlet_eta_name,
      character(1)
    ),
    .bt_parameter_coordinates_private_auxiliaries(prior_list)
  )
  bases <- .bt_parameter_coordinates_base(coordinate_names)
  column_groups <- split(columns, bases[seq_along(columns)])
  base_dimensions <- vapply(
    column_groups, .bt_parameter_coordinates_dimensions, character(1)
  )
  dimensions <- unname(base_dimensions[match(bases, names(base_dimensions))])
  dimensions[is.na(dimensions)] <- ""
  coordinates <- .bt_parameter_coordinates_empty()
  coordinates <- coordinates[rep(NA_integer_, length(coordinate_names)), , drop = FALSE]
  factor_cache <- new.env(parent = emptyenv())

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
    random_prior_auxiliary <- base_name %in% random_prior_auxiliaries
    if(random_prior_auxiliary){
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
    monitor_status <- if(.bt_parameter_coordinates_structural_point(prior)){
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
      dimensions[i],
      .bt_parameter_coordinates_scale(role, formula_parameter, formula_scale),
      monitor_status,
      fixed_value,
      .bt_parameter_coordinates_display(
        coordinate_name,
        prior,
        random_term,
        role,
        formula_parameter,
        factor_cache
      ),
      grouping,
      structure,
      role %in% c(
        "backend_anchor",
        "allocation",
        "random_latent",
        "random_group_coefficient",
        "random_mean_coordinate",
        "random_correlation",
        "random_correlation_coordinate",
        "random_inclusion_indicator",
        "random_inclusion_probability",
        "random_sd_variable"
      ) ||
        base_name %in% dirichlet_auxiliaries ||
        isTRUE(prior_metadata$allocation),
      NA_character_
    )
  }

  coordinates$convergence_role <- .bt_parameter_coordinates_convergence_roles(
    coordinates = coordinates,
    columns = columns,
    prior_list = prior_list,
    formula_design = formula_design,
    backend_anchor = backend_anchor,
    add_parameters = add_parameters,
    model_syntax = model_syntax,
    data_names = data_names
  )
  rownames(coordinates) <- NULL
  class(coordinates) <- c("BayesTools_parameter_coordinates", "data.frame")
  .bt_validate_parameter_coordinates(coordinates)
  coordinates
}


# Convergence roles ------------------------------------------------------------
#
# Every fitted coordinate gets one convergence role at fit time, derived from
# declared metadata and the parsed model syntax, never from the draws:
#
# - "sampled": a stochastic coordinate; checked by default.
# - "indicator": a declared model indicator (mixture component, spike and
#   slab, variance-allocation inclusion); checked with 'check_indicators'.
# - "structural": a declared constant (a point prior, a one-state indicator,
#   reference and fixed publication-weight bins, a p-hacking kind shared by
#   every branch, an ordered point total and the coefficients it fixes, unit
#   correlation diagonals and Cholesky constants) or a monitored deterministic
#   node whose ancestors in the model syntax are all constants.
# - "derived": a deterministic function of sampled nodes that BayesTools
#   generates for formulas (correlation matrices, Cholesky factors, derived
#   latent effects and SDs), a point prior with an expression location, or a
#   mirrored two-sided publication-weight bin; checked only on request.
# - "auxiliary": implementation nodes and hyperparameters (inclusion
#   probabilities, the backend anchor, private weight-function nodes, Dirichlet
#   gamma draws, unreported p-hacking parameters); checked only on request.
# Mirrored bins and private implementation nodes are hidden from the summary
# tables and from the convergence targets, so they cannot be requested.
.bt_parameter_coordinates_convergence_roles <- function(coordinates, columns,
                                                        prior_list,
                                                        formula_design,
                                                        backend_anchor,
                                                        add_parameters,
                                                        model_syntax,
                                                        data_names){

  coordinate_names <- coordinates$coordinate_name
  bases <- .bt_parameter_coordinates_base(coordinate_names)
  roles <- rep(NA_character_, length(coordinate_names))
  assign_roles <- function(mask, values){
    mask <- mask & is.na(roles)
    mask[is.na(mask)] <- FALSE
    roles[mask] <<- rep_len(values, length(roles))[mask]
  }

  # The backend anchor only keeps an otherwise empty monitor non-empty.
  if(!is.null(backend_anchor)){
    assign_roles(coordinate_names == backend_anchor, "auxiliary")
  }

  # Implementation nodes that summaries do not show: private weight-function
  # nodes, mirrored two-sided bins, Dirichlet gamma draws, and unreported
  # p-hacking parameters.
  hidden <- setdiff(
    columns,
    .bt_convergence_visible_columns(columns, prior_list)$column
  )
  assign_roles(coordinate_names %in% hidden & bases == "omega", "derived")
  assign_roles(coordinate_names %in% hidden, "auxiliary")

  declarations <- .bt_convergence_role_declarations(prior_list)
  exact <- declarations[declarations$exact, , drop = FALSE]
  assign_roles(
    coordinate_names %in% exact$name,
    exact$role[match(coordinate_names, exact$name)]
  )

  # Correlation matrices of random-effect blocks have an exact unit diagonal;
  # their Cholesky factors have an exact unit first diagonal element and a
  # zero upper triangle.
  cell <- strsplit(coordinates$index, ",", fixed = TRUE)
  cell_row <- suppressWarnings(vapply(cell, function(index){
    if(length(index) == 2L) as.integer(index[[1L]]) else NA_integer_
  }, integer(1)))
  cell_column <- suppressWarnings(vapply(cell, function(index){
    if(length(index) == 2L) as.integer(index[[2L]]) else NA_integer_
  }, integer(1)))
  correlation <- coordinates$role == "random_correlation"
  assign_roles(
    correlation & endsWith(bases, "_xRE_CORx_R") & cell_row == cell_column,
    "structural"
  )
  assign_roles(
    correlation & endsWith(bases, "_xRE_CORx_L") &
      (cell_column > cell_row | (cell_row == 1L & cell_column == 1L)),
    "structural"
  )

  # Variance-allocation inclusion indicators.
  allocation_indicators <- .bt_convergence_role_allocation_indicators(
    formula_design,
    prior_list
  )
  assign_roles(
    coordinate_names %in% names(allocation_indicators),
    unname(allocation_indicators[coordinate_names])
  )

  # Base-name declarations, including the spike-and-slab prior of a
  # random-effect SD, whose point inclusion at 0 or 1 fixes the indicator.
  base_declarations <- declarations[!declarations$exact, , drop = FALSE]
  assign_roles(
    bases %in% base_declarations$name,
    base_declarations$role[match(bases, base_declarations$name)]
  )

  # Spike-and-slab auxiliaries of random-effect SDs that no prior declares.
  assign_roles(coordinates$role == "random_inclusion_indicator", "indicator")
  assign_roles(coordinates$role == "random_inclusion_probability", "auxiliary")
  assign_roles(coordinates$role == "random_sd_variable", "sampled")

  # Nodes monitored through 'add_parameters' (user supplied or generated for
  # formulas) are classified from the model syntax.
  monitored <- is.na(roles) &
    bases %in% .bt_parameter_coordinates_base(add_parameters)
  if(any(monitored)){
    node_roles <- .bt_convergence_role_monitored_nodes(
      nodes = unique(bases[monitored]),
      prior_list = prior_list,
      declarations = declarations,
      formula_design = formula_design,
      model_syntax = model_syntax,
      data_names = data_names
    )
    assign_roles(monitored, unname(node_roles[bases]))
  }

  assign_roles(rep(TRUE, length(roles)), "sampled")
  roles
}

# Posterior columns that summaries and convergence checks show, and their
# labels: the prior list's private implementation nodes are removed and
# publication-weight bins are named by their p-value interval. The optional
# 'remove_parameters' are removed as in the summary tables.
.bt_convergence_visible_columns <- function(columns, prior_list,
                                            remove_parameters = NULL){

  columns <- as.character(columns)
  if(length(columns) == 0L){
    return(data.frame(
      label = character(),
      column = character(),
      stringsAsFactors = FALSE
    ))
  }
  positions <- matrix(
    seq_along(columns),
    nrow = 1L,
    dimnames = list(NULL, columns)
  )
  visible <- .remove_auxiliary_parameters(
    positions,
    if(is.null(prior_list)) list() else prior_list,
    remove_parameters
  )$model_samples
  labels <- colnames(visible)
  if(is.null(labels)){
    labels <- character()
  }
  data.frame(
    label = labels,
    column = columns[as.integer(visible[1L, seq_along(labels)])],
    stringsAsFactors = FALSE
  )
}

# Convergence roles that the priors declare for their monitored nodes: exact
# coordinate names (such as a reference publication-weight bin) and base
# names (all elements of a node).
.bt_convergence_role_declarations <- function(prior_list){

  declarations <- .bt_convergence_role_declaration(character(), character())
  prior_names <- names(prior_list)
  if(length(prior_list) == 0L || is.null(prior_names)){
    return(declarations)
  }
  for(i in seq_along(prior_list)){
    parameter <- prior_names[[i]]
    if(is.na(parameter) || !nzchar(parameter)){
      next
    }
    declarations <- rbind(
      declarations,
      .bt_convergence_role_prior(parameter, prior_list[[i]])
    )
  }
  rownames(declarations) <- NULL
  declarations
}

.bt_convergence_role_declaration <- function(name, role, exact = FALSE){

  data.frame(
    name = as.character(name),
    exact = rep(exact, length(name)),
    role = rep(as.character(role), length.out = length(name)),
    stringsAsFactors = FALSE
  )
}

.bt_convergence_role_prior <- function(parameter, prior){

  declare <- .bt_convergence_role_declaration
  if(is.null(prior)){
    return(declare(character(), character()))
  }

  if(is.prior.weightfunction(prior)){
    # The first bin is the reference bin; fixed weights fix every bin.
    cuts <- weightfunctions_mapping(list(prior), cuts_only = TRUE)
    bins <- paste0("omega[", seq_len(max(length(cuts) - 1L, 1L)), "]")
    structural <- if(identical(prior$weights$type, "fixed")){
      bins
    }else{
      bins[[1L]]
    }
    return(rbind(
      declare(structural, "structural", exact = TRUE),
      declare("omega", "sampled")
    ))
  }

  if(is_prior_bias(prior) || is_prior_phacking(prior) ||
     inherits(prior, "prior.bias_mixture")){
    out <- rbind(
      declare(.bt_convergence_role_constant_bins(prior), "structural", exact = TRUE),
      declare(
        "phack_kind",
        if(.bt_convergence_role_constant_phacking_kind(prior)) "structural" else "sampled"
      ),
      declare("omega", "sampled")
    )
    if(inherits(prior, "prior.bias_mixture")){
      out <- rbind(out, declare(
        "bias_indicator",
        if(length(prior) == 1L) "structural" else "indicator"
      ))
    }
    return(out)
  }

  if(is.prior.PET(prior) || is.prior.PEESE(prior)){
    return(declare(
      unique(c(parameter, if(is.prior.PET(prior)) "PET" else "PEESE")),
      .bt_convergence_role_point(prior)
    ))
  }

  if(!is.null(attr(prior, "random_allocation_inclusion", exact = TRUE))){
    # A variance-allocation inclusion probability.
    return(declare(
      parameter,
      if(is.null(.bt_convergence_point_value(prior))) "auxiliary" else "structural"
    ))
  }

  if(is.prior.spike_and_slab(prior)){
    inclusion <- .bt_convergence_point_value(.get_spike_and_slab_inclusion(prior))
    return(rbind(
      declare(parameter, "sampled"),
      declare(
        paste0(parameter, "_indicator"),
        if(!is.null(inclusion) && inclusion %in% c(0, 1)) "structural" else "indicator"
      ),
      declare(
        paste0(parameter, "_inclusion"),
        if(is.null(inclusion)) "auxiliary" else "structural"
      ),
      .bt_convergence_role_prior(
        paste0(parameter, "_variable"),
        .get_spike_and_slab_variable(prior)
      )
    ))
  }

  if(is.prior.mixture(prior)){
    return(rbind(
      declare(
        parameter,
        if(length(prior) == 1L) .bt_convergence_role_point(prior[[1L]]) else "sampled"
      ),
      declare(
        paste0(parameter, "_indicator"),
        if(length(prior) == 1L) "structural" else "indicator"
      )
    ))
  }

  if(is.prior.factor(prior) && is.prior.ordered(prior)){
    # A point total fixes the level coefficients when it is zero or when every
    # allocation is fixed.
    total <- .bt_convergence_point_value(prior$total)
    metadata <- attr(prior, "ordered_metadata", exact = TRUE)
    fixed_allocation <- !is.null(metadata) &&
      all(vapply(metadata$allocations, function(record){
        identical(record$spec$type, "fixed")
      }, logical(1)))
    return(rbind(
      declare(
        parameter,
        if(!is.null(total) && (total == 0 || fixed_allocation)) "structural" else "sampled"
      ),
      .bt_convergence_role_prior(.prior_ordered_total_name(parameter), prior$total)
    ))
  }

  declare(parameter, .bt_convergence_role_point(prior))
}

# Role of a prior's own node: point priors are constants unless their location
# is an expression of other nodes.
.bt_convergence_role_point <- function(prior){

  if(!is.prior.point(prior)){
    return("sampled")
  }
  if(.is_prior_expression(prior)){
    return("derived")
  }
  "structural"
}

.bt_convergence_point_value <- function(prior){

  if(is.prior.mixture(prior) && length(prior) == 1L){
    return(.bt_convergence_point_value(prior[[1L]]))
  }
  if(!is.prior.point(prior) || .is_prior_expression(prior)){
    return(NULL)
  }

  location <- prior$parameters[["location"]]
  if(!is.numeric(location) || length(location) != 1L || !is.finite(location)){
    return(NULL)
  }

  location
}

# Composed bias priors, p-hacking priors, and publication-bias mixtures share
# one omega vector on the global one-sided cut grid of the selection backend.
# A global bin is constant by construction when every mixture branch fixes it
# to the same value: branches without a step selection contribute 1, the
# reference bin (local bin 1) of a selection is 1, and every bin of fixed
# weights is its declared weight. Mirrored two-sided bins map to their local
# bin through the component expansion. Returns the monitored coordinate names.
.bt_convergence_role_constant_bins <- function(prior){

  branch_info <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
  if(!any(has_selection) && !any(has_phacking)){
    return(character())
  }

  cuts <- if(any(has_selection)){
    weightfunctions_mapping(
      lapply(branch_info[has_selection], function(x) x$selection),
      cuts_only = TRUE,
      one_sided = TRUE
    )
  }else{
    c(0, 1)
  }
  n_bins <- length(cuts) - 1L

  branch_values <- vapply(branch_info, function(x){
    .bt_convergence_selection_bin_values(x$selection, cuts)
  }, numeric(n_bins))
  branch_values <- matrix(branch_values, nrow = n_bins)
  constant <- apply(branch_values, 1L, function(values){
    all(!is.na(values)) && all(values == values[[1L]])
  })

  bin_names <- if(n_bins == 1L){
    "omega"
  }else{
    paste0("omega[", seq_len(n_bins), "]")
  }

  bin_names[constant]
}

.bt_convergence_selection_bin_values <- function(selection, cuts){

  n_bins <- length(cuts) - 1L
  if(is.null(selection)){
    return(rep(1, n_bins))
  }

  expansion <- .weightfunction_mapping_expansion(selection, force_one_sided = TRUE)
  local_bins <- expansion$index[.weightfunction_global_bin_indices(cuts, expansion)]
  if(identical(selection$weights$type, "fixed")){
    return(as.numeric(selection$weights$omega[local_bins]))
  }

  ifelse(local_bins == 1L, 1, NA_real_)
}

# The monitored p-hacking kind is a declared constant of each mixture branch:
# the form code of a p-hacking branch and 0 for branches without p-hacking.
# It is structural when every branch declares the same code.
.bt_convergence_role_constant_phacking_kind <- function(prior){

  branch_info <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
  has_phacking <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
  if(!any(has_phacking)){
    return(FALSE)
  }
  kinds <- vapply(branch_info, function(x){
    if(is.null(x$phacking)) 0 else as.numeric(.phack_kind(x$phacking$form))
  }, numeric(1))

  all(kinds == kinds[[1L]])
}

# Variance-allocation inclusion indicators; an indicator whose inclusion
# probability is a point prior at 0 or 1 is a constant.
.bt_convergence_role_allocation_indicators <- function(formula_design,
                                                       prior_list){

  indicators <- unique(
    .bt_random_variance_allocation_inclusion_indicator_names(formula_design)
  )
  roles <- stats::setNames(rep("indicator", length(indicators)), indicators)
  for(prior in prior_list){
    indicator <- attr(prior, "random_allocation_indicator", exact = TRUE)
    if(!is.character(indicator) || length(indicator) != 1L ||
       !indicator %in% indicators){
      next
    }
    probability <- .bt_convergence_point_value(prior)
    if(!is.null(probability) && probability %in% c(0, 1)){
      roles[[indicator]] <- "structural"
    }
  }

  roles
}

# Roles of monitored 'add_parameters' nodes from the model syntax. A stochastic
# node is sampled, unless it is fully observed data. A deterministic node is
# structural when none of its ancestors is stochastic; otherwise it is derived
# when BayesTools generated it for a formula and sampled when the user
# requested it. A node the syntax does not define is sampled.
.bt_convergence_role_monitored_nodes <- function(nodes, prior_list,
                                                 declarations, formula_design,
                                                 model_syntax, data_names){

  graph <- .bt_jags_syntax_graph(model_syntax)
  name_map <- .bt_parameter_coordinates_name_map(formula_design)
  generated <- name_map$jags_name[name_map$kind != "formula_output"]
  prior_bases <- declarations[!declarations$exact, , drop = FALSE]
  constants <- unique(c(
    data_names,
    prior_bases$name[prior_bases$role == "structural"]
  ))

  roles <- vapply(nodes, function(node){
    if(node %in% graph$data){
      return("structural")
    }
    if(node %in% graph$stochastic){
      return(if(node %in% data_names) "structural" else "sampled")
    }
    if(!node %in% names(graph$deterministic)){
      return("sampled")
    }
    stochastic_parent <- .bt_jags_stochastic_ancestry(
      graph$deterministic[[node]],
      graph = graph,
      constants = constants
    )
    if(!any(stochastic_parent)){
      return("structural")
    }
    if(node %in% generated) "derived" else "sampled"
  }, character(1))

  stats::setNames(roles, nodes)
}

# Whether each node is stochastic or depends on a stochastic node. Data and
# the supplied constants (fully observed data names, point priors) are not;
# nodes the syntax does not define are stochastic unless declared constant,
# which is the conservative reading of prior nodes and unknown names.
.bt_jags_stochastic_ancestry <- function(nodes, graph, constants){

  memo <- new.env(parent = emptyenv())
  visit <- function(node){
    known <- memo[[node]]
    if(!is.null(known)){
      return(known)
    }
    # A definition cycle is invalid JAGS; treat it conservatively.
    assign(node, TRUE, envir = memo)
    result <- if(node %in% graph$data){
      FALSE
    }else if(node %in% graph$stochastic){
      !node %in% constants
    }else if(node %in% names(graph$deterministic)){
      any(vapply(graph$deterministic[[node]], visit, logical(1)))
    }else{
      !node %in% constants
    }
    assign(node, result, envir = memo)
    result
  }

  vapply(as.character(nodes), visit, logical(1), USE.NAMES = FALSE)
}

# Dependency graph of a JAGS model: the nodes defined stochastically ('~')
# outside a data block, the nodes defined in a data block, and for each node
# defined deterministically ('<-' or '=') the names its definitions read.
# Statements are delimited by the grammar, not by lines, so definitions that
# span several lines are read whole. Loop indices and function, distribution,
# and link-function names are not nodes; a link function on the left-hand
# side defines its argument.
.bt_jags_syntax_graph <- function(model_syntax){

  graph <- list(
    stochastic = character(),
    deterministic = list(),
    data = character()
  )
  tokens <- .bt_jags_syntax_tokens(model_syntax)
  n <- length(tokens)
  if(n == 0L){
    return(graph)
  }

  identifier <- grepl("^[A-Za-z]", tokens)
  next_token <- c(tokens[-1L], "")
  delta <- as.integer(tokens %in% c("(", "[")) -
    as.integer(tokens %in% c(")", "]"))
  level <- cumsum(delta) - delta
  call_name <- identifier & next_token == "("
  loop_variables <- unique(tokens[
    which(tokens == "for" & next_token == "(") + 2L
  ])

  in_data <- logical(n)
  for(start in which(tokens == "data" & next_token == "{" & level == 0L)){
    depth <- 0L
    for(k in seq.int(start + 1L, n)){
      if(tokens[[k]] == "{"){
        depth <- depth + 1L
      }else if(tokens[[k]] == "}"){
        depth <- depth - 1L
        if(depth == 0L){
          in_data[start:k] <- TRUE
          break
        }
      }
    }
  }

  relations <- which(tokens %in% c("<-", "~", "=") & level == 0L)
  if(length(relations) == 0L){
    return(graph)
  }
  lhs_start <- integer(length(relations))
  lhs_name <- rep(NA_character_, length(relations))
  for(r in seq_along(relations)){
    j <- relations[[r]] - 1L
    if(j >= 1L && tokens[[j]] == "]"){
      j <- .bt_jags_syntax_match_open(tokens, j) - 1L
    }
    if(j >= 1L && tokens[[j]] == ")"){
      open <- .bt_jags_syntax_match_open(tokens, j)
      lhs_start[[r]] <- max(open - 1L, 1L)
      if(open + 1L <= n && identifier[[open + 1L]]){
        lhs_name[[r]] <- tokens[[open + 1L]]
      }
    }else if(j >= 1L){
      lhs_start[[r]] <- j
      if(identifier[[j]]){
        lhs_name[[r]] <- tokens[[j]]
      }
    }else{
      lhs_start[[r]] <- relations[[r]]
    }
  }

  boundaries <- which(level == 0L & (
    tokens %in% c("{", "}", ";") | (tokens == "for" & next_token == "(")
  ))
  starts <- sort(unique(c(lhs_start, boundaries, n + 1L)))
  node_token <- identifier & !call_name &
    !tokens %in% c(loop_variables, "for", "in")

  for(r in seq_along(relations)){
    node <- lhs_name[[r]]
    if(is.na(node)){
      next
    }
    k <- relations[[r]]
    if(in_data[[k]]){
      graph$data <- union(graph$data, node)
      next
    }
    if(tokens[[k]] == "~"){
      graph$stochastic <- union(graph$stochastic, node)
      next
    }
    end <- starts[starts > k][[1L]] - 1L
    rhs <- if(end > k) seq.int(k + 1L, end) else integer()
    parents <- unique(tokens[rhs][node_token[rhs]])
    graph$deterministic[[node]] <- unique(c(graph$deterministic[[node]], parents))
  }

  graph
}

.bt_jags_syntax_tokens <- function(model_syntax){

  if(length(model_syntax) == 0L){
    return(character())
  }
  text <- paste(model_syntax, collapse = "\n")
  text <- gsub("#[^\n]*", "", text)
  pattern <- paste0(
    "[A-Za-z][A-Za-z0-9._]*|",
    "[0-9]*[.]?[0-9]+(?:[eE][-+]?[0-9]+)?|",
    "<-|[^[:space:]]"
  )
  tokens <- regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
  tokens[nzchar(tokens)]
}

.bt_jags_syntax_match_open <- function(tokens, close){

  closer <- tokens[[close]]
  opener <- if(identical(closer, ")")) "(" else "["
  depth <- 0L
  for(k in seq.int(close, 1L)){
    if(tokens[[k]] == closer){
      depth <- depth + 1L
    }else if(tokens[[k]] == opener){
      depth <- depth - 1L
      if(depth == 0L){
        return(k)
      }
    }
  }

  1L
}
