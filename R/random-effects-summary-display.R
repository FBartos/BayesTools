.bt_random_effect_summary_metadata_table <- function(parameter_names,
                                                     prior_list,
                                                     coordinates = NULL,
                                                     formula_design = NULL){

  if(is.null(coordinates)){
    coordinates <- .bt_build_parameter_coordinates(
      columns = parameter_names,
      prior_list = prior_list,
      formula_design = formula_design
    )
  }

  n_parameters <- length(parameter_names)
  out <- data.frame(
    "Random name" = rep("", n_parameters),
    "Random grouping" = rep("", n_parameters),
    "Random structure" = rep("", n_parameters),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if(n_parameters == 0L){
    return(out)
  }

  for(i in seq_along(parameter_names)){
    metadata <- .bt_random_effect_summary_metadata_for_parameter(
      parameter_name = parameter_names[i],
      prior_list = prior_list,
      coordinates = coordinates
    )
    out[i, ] <- metadata
  }

  out
}

.bt_random_effect_summary_add_metadata_columns <- function(table,
                                                          parameter_names,
                                                          prior_list,
                                                          coordinates = NULL,
                                                          formula_design = NULL){

  metadata <- .bt_random_effect_summary_metadata_table(
    parameter_names = parameter_names,
    prior_list = prior_list,
    coordinates = coordinates,
    formula_design = formula_design
  )
  for(i in rev(seq_len(ncol(metadata)))){
    table <- add_column(
      table = table,
      column_title = colnames(metadata)[i],
      column_values = metadata[[i]],
      column_position = 1,
      column_type = "string"
    )
  }
  attr(table, "random_effects_metadata") <- TRUE

  table
}

.bt_random_effect_summary_metadata_for_parameter <- function(parameter_name,
                                                            prior_list,
                                                            coordinates = NULL,
                                                            formula_design = NULL){

  prior <- .bt_random_effect_summary_prior_for_column(parameter_name, prior_list)
  if(!is.null(prior) && isTRUE(.bt_random_effect_metadata(prior)$any)){
    effect_label <- .bt_random_effect_prior_name(prior)
    effect <- .bt_random_effect_prior_effect(prior)
    if(!nzchar(effect_label)){
      effect_label <- effect
    }
    return(c(
      "Random name" = effect_label,
      "Random grouping" = .bt_random_effect_prior_grouping(prior),
      "Random structure" = .bt_random_effect_prior_structure(prior)
    ))
  }

  raw_metadata <- .bt_random_effect_summary_raw_metadata_for_parameter(
    parameter_name = parameter_name,
    coordinates = coordinates,
    formula_design = formula_design
  )
  c(
    "Random name" = raw_metadata$name,
    "Random grouping" = raw_metadata$grouping,
    "Random structure" = raw_metadata$structure
  )
}

.bt_random_effect_summary_prior_for_column <- function(parameter_name,
                                                       prior_list){

  if(length(prior_list) == 0L){
    return(NULL)
  }
  if(parameter_name %in% names(prior_list)){
    return(prior_list[[parameter_name]])
  }

  base_name <- sub("\\[[^\\]]+\\]$", "", parameter_name)
  if(base_name %in% names(prior_list)){
    return(prior_list[[base_name]])
  }

  NULL
}

.bt_random_effect_summary_raw_metadata_for_parameter <- function(parameter_name,
                                                                coordinates = NULL,
                                                                formula_design = NULL){

  out <- list(name = "", grouping = "", structure = "")
  if(is.null(coordinates)){
    coordinates <- .bt_build_parameter_coordinates(
      columns = parameter_name,
      formula_design = formula_design
    )
  }
  .bt_validate_parameter_coordinates(coordinates)
  row <- match(parameter_name, coordinates$coordinate_name)
  if(is.na(row)){
    return(out)
  }
  coordinate_row <- coordinates[row, , drop = FALSE]
  if(!startsWith(coordinate_row$role, "random_") &&
     !identical(coordinate_row$role, "allocation")){
    return(out)
  }

  list(
    name = coordinate_row$random_name,
    grouping = coordinate_row$random_grouping,
    structure = coordinate_row$random_structure
  )
}

.bt_random_effect_summary_random_terms <- function(random_design){

  random_terms <- list()
  for(design in random_design){
    random_terms <- c(random_terms, design$random_effects)
  }

  random_terms
}

.bt_random_effect_summary_filter_raw_columns <- function(model_samples,
                                                         coordinates = NULL,
                                                         formula_design = NULL,
                                                         remove_random_effects = NULL,
                                                         keep_random_effects = NULL,
                                                         remove_random_structures = NULL,
                                                         keep_random_structures = NULL){

  column_names <- colnames(model_samples)
  if(length(column_names) == 0L){
    return(model_samples)
  }

  if(is.null(coordinates)){
    coordinates <- .bt_build_parameter_coordinates(
      columns = column_names,
      formula_design = formula_design
    )
  }
  .bt_validate_parameter_coordinates(coordinates)
  coordinate_rows <- match(column_names, coordinates$coordinate_name)
  registered <- !is.na(coordinate_rows)
  row_role <- rep("", length(column_names))
  row_block <- rep("", length(column_names))
  row_name <- rep("", length(column_names))
  row_grouping <- rep("", length(column_names))
  row_structure <- rep("", length(column_names))
  row_role[registered] <- coordinates$role[coordinate_rows[registered]]
  row_block[registered] <- coordinates$random_block[coordinate_rows[registered]]
  row_name[registered] <- coordinates$random_name[coordinate_rows[registered]]
  row_grouping[registered] <- coordinates$random_grouping[coordinate_rows[registered]]
  row_structure[registered] <- coordinates$random_structure[coordinate_rows[registered]]
  random_columns <- registered &
    (startsWith(row_role, "random_") | row_role == "allocation")
  remove_columns <- rep(FALSE, length(column_names))
  blocks <- unique(row_block[random_columns])
  for(block in blocks){
    term_columns <- random_columns & row_block == block
    if(!any(term_columns)){
      next
    }
    if(!is.null(remove_random_effects)){
      term_matches_remove_effect <- block %in% remove_random_effects ||
        any(row_name[term_columns] %in% remove_random_effects) ||
        any(row_grouping[term_columns] %in% remove_random_effects)
      remove_columns <- remove_columns | (term_columns & term_matches_remove_effect)
    }
    if(!is.null(remove_random_structures)){
      term_matches_remove_structure <- any(
        row_structure[term_columns] %in% remove_random_structures
      )
      remove_columns <- remove_columns | (term_columns & term_matches_remove_structure)
    }
    if(!is.null(keep_random_effects) || !is.null(keep_random_structures)){
      term_matches <- TRUE
      if(!is.null(keep_random_effects)){
        term_matches <- term_matches &&
          (any(row_block[term_columns] %in% keep_random_effects) ||
             any(row_name[term_columns] %in% keep_random_effects) ||
             any(row_grouping[term_columns] %in% keep_random_effects))
      }
      if(!is.null(keep_random_structures)){
        term_matches <- term_matches &&
          any(row_structure[term_columns] %in% keep_random_structures)
      }
      remove_columns <- remove_columns | (term_columns & !term_matches)
    }
  }

  model_samples[, !remove_columns, drop = FALSE]
}

.bt_random_effect_summary_term_filter_matches <- function(random_term,
                                                         random_effects = NULL,
                                                         random_structures = NULL){

  matches <- TRUE
  if(!is.null(random_effects)){
    effect_names <- c(
      random_term$block_name,
      .bt_random_effect_public_name(random_term),
      random_term$group_label
    )
    matches <- matches && any(effect_names %in% random_effects)
  }
  if(!is.null(random_structures)){
    matches <- matches &&
      .bt_random_effect_summary_term_structure(random_term) %in% random_structures
  }

  matches
}

.bt_random_effect_summary_term_structure <- function(random_term){

  .bt_random_effect_structure(
    random_term,
    context = "Random-effect summary metadata"
  )
}

.bt_random_effect_summary_remove_raw <- function(model_samples, prior_list,
                                                 coordinates){

  raw_prior_names <- names(prior_list)[vapply(
    prior_list,
    .bt_random_effect_summary_is_raw_prior,
    logical(1)
  )]

  .bt_validate_parameter_coordinates(coordinates)
  coordinate_rows <- match(
    colnames(model_samples),
    coordinates$coordinate_name
  )
  registered <- !is.na(coordinate_rows)
  roles <- rep("", ncol(model_samples))
  roles[registered] <- coordinates$role[coordinate_rows[registered]]
  raw_cols <- registered &
    (startsWith(roles, "random_") | roles == "allocation")

  if(any(raw_cols)){
    model_samples <- model_samples[, !raw_cols, drop = FALSE]
  }
  if(length(raw_prior_names) > 0L){
    prior_list <- prior_list[!names(prior_list) %in% raw_prior_names]
  }

  list(model_samples = model_samples, prior_list = prior_list)
}

.bt_random_effect_summary_is_raw_prior <- function(prior){

  .bt_is_random_effect_prior(prior, include_summary = FALSE)
}

.bt_random_effect_summary_parameter_columns <- function(column_names,
                                                        parameter_names){

  if(length(parameter_names) == 0L || length(column_names) == 0L){
    return(rep(FALSE, length(column_names)))
  }

  out <- rep(FALSE, length(column_names))
  for(parameter_name in parameter_names){
    eta_name <- .JAGS_prior_dirichlet_eta_name(parameter_name)
    out <- out |
      column_names == parameter_name |
      startsWith(column_names, paste0(parameter_name, "[")) |
      column_names == eta_name |
      startsWith(column_names, paste0(eta_name, "["))
  }

  out
}

.bt_random_effect_summary_display_names <- function(names, raw_names,
                                                     prior_list,
                                                     formula_prefix = TRUE,
                                                     coordinates = NULL,
                                                     formula_design = NULL){

  if(length(raw_names) == 0L){
    return(names)
  }
  if(is.null(coordinates)){
    coordinates <- .bt_build_parameter_coordinates(
      columns = raw_names,
      prior_list = prior_list,
      formula_design = formula_design
    )
  }

  if(length(prior_list) > 0L){
    for(i in seq_along(raw_names)){
      prior <- prior_list[[raw_names[i]]]
      if(is.null(prior)){
        next
      }
      label <- attr(prior, "random_summary_label", exact = TRUE)
      if(is.null(label)){
        next
      }
      parameter <- attr(prior, "parameter")
      prefix <- .bt_random_effect_summary_formula_prefix(parameter, formula_prefix)
      names[i] <- paste0(prefix, label)
    }
  }

  .bt_validate_parameter_coordinates(coordinates)
  coordinate_rows <- match(raw_names, coordinates$coordinate_name)
  registered <- !is.na(coordinate_rows)
  registered[registered] <- nzchar(
    coordinates$random_block[coordinate_rows[registered]]
  )
  if(any(registered)){
    labels <- coordinates$display_label[coordinate_rows[registered]]
    if(!isTRUE(formula_prefix)){
      labels <- sub("^\\([^)]*\\) ", "", labels)
    }
    names[registered] <- labels
  }

  names
}

.bt_random_effect_summary_formula_prefix <- function(parameter, formula_prefix){

  if(isTRUE(formula_prefix) && !is.null(parameter) &&
     length(parameter) == 1L && nzchar(parameter)){
    return(paste0("(", parameter, ") "))
  }

  ""
}

.bt_parameter_catalog_random_summary_samples <- function(
    fit, model_samples, prior_list, coordinates,
    mode = c("standard", "full", "raw", "none"),
    simplify_names = FALSE){

  mode <- match.arg(mode)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
  if(identical(mode, "raw")){
    return(list(model_samples = model_samples, prior_list = prior_list))
  }
  cleaned <- .bt_random_effect_summary_remove_raw(
    model_samples = model_samples,
    prior_list = prior_list,
    coordinates = coordinates
  )
  if(identical(mode, "none")){
    return(cleaned)
  }

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  keep <- startsWith(quantities$role, "random_") &
    !quantities$internal & quantities$status != "unavailable"
  if(identical(mode, "standard")){
    keep <- keep & .bt_parameter_catalog_random_standard_quantities(quantities)
  }
  quantities <- quantities[keep, , drop = FALSE]
  quantities <- .bt_parameter_catalog_random_summary_order(
    quantities,
    mode = mode
  )
  if(nrow(quantities) == 0L){
    return(cleaned)
  }

  summary_columns <- vector("list", nrow(quantities))
  summary_priors <- vector("list", nrow(quantities))
  for(i in seq_len(nrow(quantities))){
    quantity <- quantities[i, , drop = FALSE]
    selection <- parameter_catalog_resolve(
      catalog,
      alias = quantity$canonical_name,
      namespace = quantity$namespace
    )
    values <- as.matrix(parameter_draws(fit, selection))
    if(ncol(values) != 1L || nrow(values) != nrow(model_samples)){
      stop(
        "A semantic random-effect summary could not be aligned with the fitted posterior draws.",
        call. = FALSE
      )
    }
    summary_columns[[i]] <- as.numeric(values[, 1L])
    summary_priors[[i]] <- .bt_parameter_catalog_random_summary_prior(
      fit,
      quantity,
      simplify_names = simplify_names
    )
  }
  names(summary_columns) <- quantities$canonical_name
  names(summary_priors) <- quantities$canonical_name
  summary_matrix <- do.call(cbind, summary_columns)
  colnames(summary_matrix) <- quantities$canonical_name

  list(
    model_samples = cbind(cleaned$model_samples, summary_matrix),
    prior_list = c(cleaned$prior_list, summary_priors)
  )
}

.bt_parameter_catalog_random_summary_order <- function(quantities, mode){

  if(nrow(quantities) < 2L || !identical(mode, "standard")){
    return(quantities)
  }
  quantity_order <- match(
    quantities$quantity,
    c(
      "sd_total", "sd_common", "sd", "sd_mult",
      "var_prop", "var_mult", "inclusion", "cor"
    )
  )
  quantities[order(quantity_order, na.last = TRUE), , drop = FALSE]
}

.bt_parameter_catalog_random_standard_quantities <- function(quantities){

  allocation_derived <- vapply(quantities$extraction_key, function(key){
    isTRUE(key$allocation_derived)
  }, logical(1))
  correlation <- quantities$quantity == "cor"
  inclusion   <- quantities$quantity == "inclusion"
  block_scale <- quantities$owner_type == "random_block" &
    quantities$quantity == "sd" &
    !allocation_derived
  allocation <- quantities$owner_type == "variance_allocation" &
    quantities$quantity %in% c(
      "sd_total", "sd_common", "var_prop", "var_mult", "inclusion"
    )

  correlation | inclusion | block_scale | allocation
}

.bt_parameter_catalog_random_summary_prior <- function(
    fit, quantity, simplify_names = FALSE){

  key <- quantity$extraction_key[[1L]]
  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  allocation <- if(!is.null(key$allocation_label)){
    .bt_parameter_catalog_find_allocation(fit, key, random_term)
  }else{
    NULL
  }
  component_label <- .bt_random_effect_semantic_quantity_name(
    quantity$quantity,
    quantity$arguments[[1L]]
  )
  summary_label <- if(simplify_names){
    quantity$display_label
  }else{
    quantity$canonical_name
  }
  prefix <- .bt_random_effect_summary_formula_prefix(
    quantity$formula_parameter,
    TRUE
  )
  if(nzchar(prefix) && startsWith(summary_label, prefix)){
    summary_label <- substring(summary_label, nchar(prefix) + 1L)
  }
  .bt_random_effect_summary_prior(
    parameter = quantity$formula_parameter,
    type = quantity$quantity,
    label = summary_label,
    component_label = component_label,
    block = if(is.null(random_term)) NULL else random_term$block_name,
    grouping = if(is.null(random_term)) NULL else random_term$group_label,
    structure = if(is.null(random_term)) NULL else
      .bt_random_effect_summary_term_structure(random_term),
    effect_label = if(is.null(random_term)) NULL else quantity$owner_name,
    allocation = if(is.null(allocation)) NULL else key$allocation_label,
    allocation_metadata = allocation,
    allocation_index = if(is.null(key$index)) NULL else key$index,
    component = if(nzchar(quantity$component)) quantity$component else NULL
  )
}

.bt_random_effect_summary_raw_display_names <- function(names, raw_names,
                                                        prior_list,
                                                        formula_prefix,
                                                        formula_design = NULL){

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(names)
  }

  for(design in random_design){
    parameter <- design$parameter
    prefix <- .bt_random_effect_summary_formula_prefix(parameter, formula_prefix)
    for(random_term in design$random_effects){
      names <- .bt_random_effect_summary_raw_sd_display_names(
        names = names,
        raw_names = raw_names,
        prior_list = prior_list,
        random_term = random_term,
        prefix = prefix
      )
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
    }
  }

  names
}

.bt_random_effect_summary_renamed_parameter_names <- function(parameter_names,
                                                              prior_list){

  if(length(parameter_names) == 0L || length(prior_list) == 0L){
    return(parameter_names)
  }

  dummy <- matrix(
    nrow = 0L,
    ncol = length(parameter_names),
    dimnames = list(NULL, parameter_names)
  )
  colnames(.rename_factor_levels(dummy, prior_list))
}

.bt_random_effect_summary_raw_sd_display_names <- function(names, raw_names,
                                                           prior_list,
                                                           random_term,
                                                           prefix){

  if(!is.null(random_term$sd_binding) &&
     .bt_random_sd_binding_has_external_source(random_term$sd_binding) &&
     !isTRUE(random_term$sd_binding$true_allocation)){
    return(names)
  }

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 0L){
    return(names)
  }

  display_sd_names <- .bt_random_effect_summary_renamed_parameter_names(
    parameter_names = sd_names,
    prior_list = prior_list
  )
  components <- .bt_random_effect_summary_sd_components(random_term, sd_names)
  labels <- paste0(prefix, vapply(
    components,
    .bt_random_effect_sd_summary_label,
    character(1),
    random_term = random_term
  ))

  for(i in seq_along(sd_names)){
    matches <- raw_names %in% c(sd_names[i], display_sd_names[i])
    names[matches] <- labels[i]
  }

  names
}

.bt_random_effect_summary_raw_rho_display_names <- function(names, raw_names,
                                                            random_term,
                                                            prefix){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || !nzchar(stem)){
    return(names)
  }

  owner <- .bt_random_effect_public_name(random_term)
  rho_names <- c(
    rho = paste0(stem, "_rho"),
    rho_z = paste0(stem, "_rho_z"),
    rho_logit = paste0(stem, "_rho_logit")
  )
  for(rho_label in names(rho_names)){
    names[raw_names == rho_names[[rho_label]]] <- .bt_random_effect_semantic_name(
      parameter = "",
      owner = owner,
      quantity = "cor",
      formula_prefix = FALSE
    )
    names[raw_names == rho_names[[rho_label]]] <- paste0(
      prefix,
      names[raw_names == rho_names[[rho_label]]]
    )
  }

  names
}

.bt_random_effect_summary_raw_matrix_display_names <- function(names,
                                                               raw_names,
                                                               random_term,
                                                               prefix){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || !nzchar(stem)){
    return(names)
  }

  components <- .bt_random_effect_summary_column_components(random_term)
  group <- .bt_random_effect_summary_group_label(random_term)
  group_levels <- random_term$group_levels

  names <- .bt_random_effect_summary_raw_correlation_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_CORx_R",
    label = "cor",
    components = components,
    group = group,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_correlation_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_CORx_L",
    label = "cor_chol",
    components = components,
    group = group,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_effect_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_Zx",
    label = "z",
    components = components,
    group = group,
    group_levels = group_levels,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_effect_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_COEFx",
    label = "coef",
    components = components,
    group = group,
    group_levels = group_levels,
    prefix = prefix
  )

  names
}

.bt_random_effect_summary_raw_correlation_matrix_names <- function(names,
                                                                   raw_names,
                                                                   stem,
                                                                   matrix,
                                                                   label,
                                                                   components,
                                                                   group,
                                                                   prefix){

  matrix_prefix <- paste0(stem, matrix)
  matches <- startsWith(raw_names, paste0(matrix_prefix, "["))
  for(i in which(matches)){
    index <- .bt_random_effect_summary_matrix_index(raw_names[i], matrix_prefix)
    if(anyNA(index) || any(index < 1L) || any(index > length(components))){
      next
    }
    names[i] <- paste0(
      prefix,
      label,
      "(",
      components[index[1L]],
      ",",
      components[index[2L]],
      " | ",
      group,
      ")"
    )
  }

  names
}

.bt_random_effect_summary_raw_effect_matrix_names <- function(names,
                                                              raw_names,
                                                              stem,
                                                              matrix,
                                                              label,
                                                              components,
                                                              group,
                                                              group_levels,
                                                              prefix){

  matrix_prefix <- paste0(stem, matrix)
  matches <- startsWith(raw_names, paste0(matrix_prefix, "["))
  for(i in which(matches)){
    index <- .bt_random_effect_summary_matrix_index(raw_names[i], matrix_prefix)
    if(anyNA(index) || index[2L] < 1L || index[2L] > length(components)){
      next
    }
    group_level <- as.character(index[1L])
    if(!is.null(group_levels) && index[1L] >= 1L &&
       index[1L] <= length(group_levels)){
      group_level <- as.character(group_levels[index[1L]])
    }
    names[i] <- paste0(
      prefix,
      label,
      "(",
      group,
      "[",
      group_level,
      "], ",
      components[index[2L]],
      ")"
    )
  }

  names
}

.bt_random_effect_summary_matrix_index <- function(x, prefix){

  rest <- substring(x, nchar(prefix) + 1L)
  if(!grepl("^\\[[0-9]+,[0-9]+\\]$", rest)){
    return(c(NA_integer_, NA_integer_))
  }

  as.integer(strsplit(substr(rest, 2L, nchar(rest) - 1L), ",", fixed = TRUE)[[1]])
}
