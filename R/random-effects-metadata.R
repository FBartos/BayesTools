# Internal random-effect prior metadata helpers.

.bt_random_effect_metadata <- function(prior){

  summary_type <- attr(prior, "random_summary")
  if(is.null(summary_type)){
    summary_type <- ""
  }

  raw_sd <- isTRUE(attr(prior, "random_sd"))
  raw_sd_total <- isTRUE(attr(prior, "random_sd_total"))
  raw_correlation <- isTRUE(attr(prior, "random_correlation"))
  allocation <- !is.null(attr(prior, "random_allocation"))
  raw <- identical(summary_type, "") &&
    (raw_sd || raw_sd_total || raw_correlation || allocation)

  list(
    summary = summary_type,
    allocation = allocation,
    raw_sd = raw_sd,
    raw_sd_total = raw_sd_total,
    raw_correlation = raw_correlation,
    raw = raw,
    any = !identical(summary_type, "") || raw
  )
}

.bt_is_random_effect_prior <- function(prior, include_summary = TRUE){

  metadata <- .bt_random_effect_metadata(prior)
  if(isTRUE(include_summary)){
    return(isTRUE(metadata$any))
  }

  isTRUE(metadata$raw)
}

.bt_random_effect_set_prior_metadata <- function(prior, block, grouping,
                                                 type = c("sd", "correlation"),
                                                 structure = NULL,
                                                 name = NULL){

  type <- match.arg(type)
  if(identical(type, "sd")){
    attr(prior, "random_sd") <- TRUE
  }else if(identical(type, "correlation")){
    attr(prior, "random_correlation") <- TRUE
  }
  attr(prior, "random_factor") <- block
  attr(prior, "random_grouping_factor") <- grouping
  if(is.null(name)){
    name <- block
  }
  attr(prior, "random_name") <- name
  if(!is.null(structure)){
    attr(prior, "random_structure") <- structure
  }

  prior
}

.bt_random_effect_set_total_sd_metadata <- function(prior, allocation, terms){

  attr(prior, "random_sd_total") <- TRUE
  attr(prior, "random_allocation") <- allocation
  attr(prior, "random_allocation_terms") <- terms

  prior
}

.bt_random_effect_set_allocation_metadata <- function(prior, allocation,
                                                      terms, parent = NULL){

  attr(prior, "random_allocation") <- allocation
  attr(prior, "random_allocation_terms") <- terms
  if(!is.null(parent)){
    attr(prior, "random_allocation_parent") <- parent
  }

  prior
}

.bt_random_effect_prior_flags <- function(prior_list){

  if(length(prior_list) == 0L){
    return(data.frame(
      name = character(),
      summary = character(),
      effect_label = character(),
      effect = character(),
      grouping = character(),
      structure = character(),
      allocation = logical(),
      raw_sd = logical(),
      raw_sd_total = logical(),
      raw_correlation = logical(),
      raw = logical(),
      any = logical()
    ))
  }

  metadata <- lapply(prior_list, .bt_random_effect_metadata)
  prior_names <- names(prior_list)
  if(is.null(prior_names)){
    prior_names <- rep("", length(prior_list))
  }
  data.frame(
    name = prior_names,
    summary = vapply(metadata, `[[`, character(1), "summary"),
    effect_label = vapply(prior_list, .bt_random_effect_prior_name, character(1)),
    effect = vapply(prior_list, .bt_random_effect_prior_effect, character(1)),
    grouping = vapply(prior_list, .bt_random_effect_prior_grouping, character(1)),
    structure = vapply(prior_list, .bt_random_effect_prior_structure, character(1)),
    allocation = vapply(metadata, `[[`, logical(1), "allocation"),
    raw_sd = vapply(metadata, `[[`, logical(1), "raw_sd"),
    raw_sd_total = vapply(metadata, `[[`, logical(1), "raw_sd_total"),
    raw_correlation = vapply(metadata, `[[`, logical(1), "raw_correlation"),
    raw = vapply(metadata, `[[`, logical(1), "raw"),
    any = vapply(metadata, `[[`, logical(1), "any"),
    stringsAsFactors = FALSE
  )
}

.bt_random_effect_prior_attr <- function(prior, attribute){

  value <- attr(prior, attribute)
  if(is.null(value) || length(value) == 0L || is.na(value[1L])){
    return("")
  }

  as.character(value[1L])
}

.bt_random_effect_public_name <- function(random_term){

  if(isTRUE(random_term$has_explicit_name) &&
     !is.null(random_term$block_name) &&
     length(random_term$block_name) == 1L &&
     nzchar(random_term$block_name)){
    return(random_term$block_name)
  }

  if(!is.null(random_term$group_label) &&
     length(random_term$group_label) == 1L &&
     nzchar(random_term$group_label)){
    return(random_term$group_label)
  }

  random_term$block_name
}

.bt_random_effect_prior_name <- function(prior){

  .bt_random_effect_prior_attr(prior, "random_name")
}

.bt_random_effect_prior_effect <- function(prior){

  .bt_random_effect_prior_attr(prior, "random_factor")
}

.bt_random_effect_prior_grouping <- function(prior){

  .bt_random_effect_prior_attr(prior, "random_grouping_factor")
}

.bt_random_effect_prior_structure <- function(prior){

  .bt_random_effect_prior_attr(prior, "random_structure")
}

.bt_random_effect_filter_matches <- function(prior_list,
                                             random_effects = NULL,
                                             random_structures = NULL){

  flags <- .bt_random_effect_prior_flags(prior_list)
  if(nrow(flags) == 0L){
    return(logical())
  }

  matches <- flags$any
  if(!is.null(random_effects)){
    matches <- matches &
      (flags$effect %in% random_effects |
         flags$effect_label %in% random_effects |
         flags$grouping %in% random_effects)
  }
  if(!is.null(random_structures)){
    matches <- matches & flags$structure %in% random_structures
  }

  matches
}

.bt_random_effect_alias_matches <- function(prior_list, alias){

  flags <- .bt_random_effect_prior_flags(prior_list)
  if(nrow(flags) == 0L){
    return(character())
  }

  random_sd <- flags$raw_sd | flags$raw_sd_total
  random_rho <- flags$raw_correlation
  is_dirichlet_allocation <- flags$summary == "" &
    flags$allocation &
    vapply(prior_list, is.prior.simplex, logical(1))

  matches <- switch(
    alias,
    "random" = flags$any,
    "random_effects" = flags$any,
    "random_sd" = flags$summary %in% c("sd", "sd_total") | random_sd,
    "random_rho" = flags$summary == "rho" | random_rho,
    "random_cor" = flags$summary %in% c("rho", "cor") | random_rho,
    "random_correlation" = flags$summary %in% c("rho", "cor") | random_rho,
    "random_variance_fraction" = flags$summary == "var_frac" | is_dirichlet_allocation,
    "random_allocation" = flags$summary %in% c("sd_total", "var_frac", "sd_multiplier") |
      flags$allocation,
    "random_sd_multiplier" = flags$summary == "sd_multiplier",
    rep(FALSE, nrow(flags))
  )

  flags$name[matches]
}
