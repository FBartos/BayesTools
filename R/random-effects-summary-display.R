.bt_random_effect_summary_metadata_table <- function(parameter_names,
                                                     prior_list,
                                                     formula_design = NULL){

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
      formula_design = formula_design
    )
    out[i, ] <- metadata
  }

  out
}

.bt_random_effect_summary_add_metadata_columns <- function(table,
                                                          parameter_names,
                                                          prior_list,
                                                          formula_design = NULL){

  metadata <- .bt_random_effect_summary_metadata_table(
    parameter_names = parameter_names,
    prior_list = prior_list,
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
                                                                formula_design = NULL){

  out <- list(name = "", grouping = "", structure = "")
  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(out)
  }

  random_terms <- .bt_random_effect_summary_random_terms(random_design)
  for(random_term in random_terms){
    if(.bt_random_effect_summary_raw_parameter_matches(
      parameter_name = parameter_name,
      random_term = random_term,
      random_terms = random_terms
    )){
      return(list(
        name = .bt_random_effect_public_name(random_term),
        grouping = .bt_random_effect_summary_group_label(random_term),
        structure = .bt_random_effect_summary_term_structure(random_term)
      ))
    }
  }

  out
}

.bt_random_effect_summary_random_terms <- function(random_design){

  random_terms <- list()
  for(design in random_design){
    random_terms <- c(random_terms, design$random_effects)
  }

  random_terms
}

.bt_random_effect_summary_raw_parameter_owner_stem <- function(parameter_name,
                                                               random_terms){

  if(length(parameter_name) != 1L || is.na(parameter_name)){
    return(NA_character_)
  }

  stems <- vapply(random_terms, function(random_term){
    stem <- random_term$parameter_stem
    if(is.null(stem) || length(stem) != 1L || is.na(stem) || !nzchar(stem)){
      return(NA_character_)
    }
    as.character(stem)
  }, character(1))
  stems <- stems[!is.na(stems)]
  if(length(stems) == 0L){
    return(NA_character_)
  }

  matching_stems <- stems[
    startsWith(parameter_name, paste0(stems, "_"))
  ]
  if(length(matching_stems) == 0L){
    return(NA_character_)
  }

  matching_stems[which.max(nchar(matching_stems))]
}

.bt_random_effect_summary_raw_parameter_matches <- function(parameter_name,
                                                           random_term,
                                                           random_terms){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || is.na(stem) || !nzchar(stem)){
    return(FALSE)
  }

  identical(
    as.character(stem),
    .bt_random_effect_summary_raw_parameter_owner_stem(
      parameter_name = parameter_name,
      random_terms = random_terms
    )
  )
}

.bt_random_effect_summary_filter_raw_columns <- function(model_samples,
                                                         formula_design = NULL,
                                                         remove_random_effects = NULL,
                                                         keep_random_effects = NULL,
                                                         remove_random_structures = NULL,
                                                         keep_random_structures = NULL){

  column_names <- colnames(model_samples)
  if(length(column_names) == 0L){
    return(model_samples)
  }

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(model_samples)
  }

  random_terms <- .bt_random_effect_summary_random_terms(random_design)
  remove_columns <- rep(FALSE, length(column_names))
  for(random_term in random_terms){
    term_columns <- vapply(
      column_names,
      .bt_random_effect_summary_raw_parameter_matches,
      logical(1),
      random_term = random_term,
      random_terms = random_terms
    )
    if(!any(term_columns)){
      next
    }
    if(!is.null(remove_random_effects)){
      term_matches_remove_effect <- .bt_random_effect_summary_term_filter_matches(
        random_term = random_term,
        random_effects = remove_random_effects,
        random_structures = NULL
      )
      remove_columns <- remove_columns | (term_columns & term_matches_remove_effect)
    }
    if(!is.null(remove_random_structures)){
      term_matches_remove_structure <- .bt_random_effect_summary_term_filter_matches(
        random_term = random_term,
        random_effects = NULL,
        random_structures = remove_random_structures
      )
      remove_columns <- remove_columns | (term_columns & term_matches_remove_structure)
    }
    if(!is.null(keep_random_effects) || !is.null(keep_random_structures)){
      term_matches <- .bt_random_effect_summary_term_filter_matches(
        random_term = random_term,
        random_effects = keep_random_effects,
        random_structures = keep_random_structures
      )
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
                                                 random_design){

  raw_prior_names <- names(prior_list)[vapply(
    prior_list,
    .bt_random_effect_summary_is_raw_prior,
    logical(1)
  )]

  raw_prefixes <- character()
  for(design in random_design){
    for(random_term in design$random_effects){
      raw_prefixes <- c(
        raw_prefixes,
        paste0(random_term$parameter_stem, "_")
      )
    }
  }
  raw_prefixes <- unique(raw_prefixes)

  raw_cols <- .bt_random_effect_summary_parameter_columns(
    colnames(model_samples),
    raw_prior_names
  )
  allocation_indicator_names <- .bt_random_variance_allocation_inclusion_indicator_names(
    random_design
  )
  if(length(allocation_indicator_names) > 0L){
    raw_cols <- raw_cols | colnames(model_samples) %in% allocation_indicator_names
  }
  if(length(raw_prefixes) > 0L){
    raw_cols <- raw_cols | Reduce(
      "|",
      lapply(raw_prefixes, function(prefix) startsWith(colnames(model_samples), prefix))
    )
  }

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
                                                    formula_design = NULL){

  if(length(raw_names) == 0L){
    return(names)
  }

  if(length(prior_list) > 0L){
    for(i in seq_along(raw_names)){
      prior <- prior_list[[raw_names[i]]]
      if(is.null(prior)){
        next
      }
      label <- attr(prior, "random_summary_label")
      if(is.null(label)){
        next
      }
      parameter <- attr(prior, "parameter")
      prefix <- .bt_random_effect_summary_formula_prefix(parameter, formula_prefix)
      names[i] <- paste0(prefix, label)
    }
  }

  .bt_random_effect_summary_raw_display_names(
    names = names,
    raw_names = raw_names,
    prior_list = prior_list,
    formula_prefix = formula_prefix,
    formula_design = formula_design
  )
}

.bt_random_effect_summary_formula_prefix <- function(parameter, formula_prefix){

  if(isTRUE(formula_prefix) && !is.null(parameter) &&
     length(parameter) == 1L && nzchar(parameter)){
    return(paste0("(", parameter, ") "))
  }

  ""
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
  group <- .bt_random_effect_summary_group_label(random_term)
  labels <- paste0(prefix, vapply(
    components,
    .bt_random_effect_sd_summary_label,
    character(1),
    group = group,
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

  group <- .bt_random_effect_summary_group_label(random_term)
  rho_names <- c(
    rho = paste0(stem, "_rho"),
    rho_z = paste0(stem, "_rho_z"),
    rho_logit = paste0(stem, "_rho_logit")
  )
  for(rho_label in names(rho_names)){
    names[raw_names == rho_names[[rho_label]]] <- paste0(
      prefix,
      rho_label,
      "(",
      group,
      ")"
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
