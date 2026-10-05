.posterior_density_for_method <- function(posterior_density, density_method){

  if(identical(density_method, "precomputed")){
    return(.posterior_density_from_attribute(posterior_density))
  }

  return(NULL)
}

.posterior_density_usable_for_null <- function(posterior_density,
                                               null_hypothesis){

  posterior_density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(posterior_density) ||
     !is.numeric(null_hypothesis) ||
     length(null_hypothesis) != 1L ||
     !is.finite(null_hypothesis)){
    return(FALSE)
  }

  support <- .posterior_support_from_attribute(posterior_density[["support"]])
  if(!is.null(support) && isTRUE(support$exact) &&
     .posterior_support_excludes_value(support, null_hypothesis)){
    return(TRUE)
  }

  height <- .posterior_density_height(posterior_density, null_hypothesis)
  is.finite(height) && height > 0
}

.posterior_density_aliases <- function(...){

  aliases <- unlist(list(...), use.names = FALSE)
  aliases <- as.character(aliases)
  aliases <- aliases[!is.na(aliases) & nzchar(aliases)]

  return(unique(aliases))
}

.posterior_density_sample_aliases <- function(samples, sample_name = NULL){

  return(.posterior_density_aliases(
    sample_name,
    attr(samples, "parameter", exact = TRUE),
    attr(samples, "level_name", exact = TRUE),
    attr(samples, "factor_cell_names", exact = TRUE)
  ))
}

.posterior_density_child_aliases <- function(parent, child, sample_name = NULL){

  parent_parameter <- attr(parent, "parameter", exact = TRUE)
  child_level <- attr(child, "level_name", exact = TRUE)
  if(is.null(child_level)){
    child_level <- attr(child, "level", exact = TRUE)
  }

  return(.posterior_density_aliases(
    sample_name,
    attr(child, "parameter", exact = TRUE),
    child_level,
    if(!is.null(parent_parameter) && !is.null(sample_name)){
      paste0(parent_parameter, "[", sample_name, "]")
    },
    if(!is.null(parent_parameter) && !is.null(child_level)){
      paste0(parent_parameter, "[", child_level, "]")
    },
    attr(child, "factor_cell_names", exact = TRUE)
  ))
}

# Metadata of a posterior density or ordinate attribute; containers carry none.
.posterior_density_metadata_values <- function(posterior_density, field){

  if(!inherits(posterior_density, c(
    "BayesTools_posterior_density",
    "BayesTools_posterior_ordinate",
    "BayesTools_posterior_ordinates"
  ))){
    return(character())
  }

  values <- as.character(unlist(posterior_density[[field]], use.names = FALSE))
  values <- values[!is.na(values) & nzchar(values)]

  return(unique(values))
}

.posterior_density_parameter_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(posterior_density, "parameter"))
}

.posterior_density_conditional_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(posterior_density, "conditional"))
}

.posterior_density_conditional_rule_metadata <- function(posterior_density){

  rule <- .posterior_density_metadata_values(posterior_density, "conditional_rule")
  if(length(rule) == 0L){
    return(NULL)
  }

  return(rule[[1]])
}

.posterior_density_condition_key_metadata <- function(posterior_density){

  key <- .posterior_density_metadata_values(posterior_density, "condition_key")
  if(length(key) == 0L){
    return(NULL)
  }

  key[[1]]
}

.posterior_density_normalize_condition <- function(conditional){

  conditional <- unlist(conditional, use.names = FALSE)
  conditional <- as.character(conditional)
  conditional <- conditional[!is.na(conditional) & nzchar(conditional)]

  return(unique(conditional))
}

.posterior_density_condition_matches <- function(posterior_density, conditional, conditional_rule,
                                                 condition_key = NULL){

  requested_conditional <- .posterior_density_normalize_condition(conditional)
  density_conditional   <- .posterior_density_normalize_condition(
    .posterior_density_conditional_metadata(posterior_density)
  )
  requested_key <- if(is.null(condition_key)){
    .condition_event_key(requested_conditional, conditional_rule)
  }else{
    as.character(condition_key)
  }
  density_key <- .posterior_density_condition_key_metadata(posterior_density)
  if(!is.null(density_key)){
    return(identical(as.character(density_key), requested_key))
  }

  if(length(requested_conditional) == 0L && length(density_conditional) == 0L){
    return(TRUE)
  }
  if(length(requested_conditional) == 0L || length(density_conditional) == 0L){
    return(FALSE)
  }

  density_rule <- .posterior_density_conditional_rule_metadata(posterior_density)
  if(length(requested_conditional) > 1L && is.null(density_rule)){
    return(FALSE)
  }
  if(!is.null(density_rule) && !identical(as.character(density_rule), as.character(conditional_rule))){
    return(FALSE)
  }
  if(is.null(density_rule)){
    density_rule <- conditional_rule
  }

  density_key <- .condition_event_key(density_conditional, density_rule)

  identical(requested_key, density_key)
}

.posterior_density_candidate_matches <- function(posterior_density, aliases, conditional, conditional_rule,
                                                 condition_key = NULL, allow_unlabeled = FALSE,
                                                 null_hypothesis = NULL){

  density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(density)){
    return(FALSE)
  }
  if(!is.null(null_hypothesis) &&
     !.posterior_density_usable_for_null(density, null_hypothesis)){
    return(FALSE)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_density)
  if(length(parameter_names) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(!allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_density,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key
  ))
}

.posterior_density_direct_candidate_matches <- function(posterior_density,
                                                        samples,
                                                        aliases = NULL,
                                                        allow_unlabeled = TRUE,
                                                        null_hypothesis = NULL){

  density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(density)){
    return(FALSE)
  }
  if(!is.null(null_hypothesis) &&
     !.posterior_density_usable_for_null(density, null_hypothesis)){
    return(FALSE)
  }
  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_density)
  if(length(parameter_names) > 0L && length(aliases) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(length(parameter_names) == 0L && !allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_density,
    conditional      = .bt_meta_condition(samples, "conditional"),
    conditional_rule = .bt_meta_condition(samples, "conditional_rule"),
    condition_key    = .bt_meta_condition(samples, "condition_key")
  ))
}

.posterior_density_direct_attribute <- function(samples, aliases = NULL,
                                                null_hypothesis = NULL){

  posterior_density <- .bt_meta_get(samples, "posterior_density")
  if(.posterior_density_direct_candidate_matches(
    posterior_density,
    samples         = samples,
    aliases         = aliases,
    null_hypothesis = null_hypothesis
  )){
    return(posterior_density)
  }

  return(NULL)
}

.posterior_density_direct_status <- function(samples, aliases = NULL,
                                             allow_unlabeled = TRUE){

  posterior_density <- .bt_meta_get(samples, "posterior_density")
  if(identical(.posterior_density_kind(posterior_density), "null")){
    return(list(present = FALSE, relevant = FALSE, valid = FALSE, value = NULL))
  }
  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_density)
  relevant <- TRUE
  if(length(parameter_names) > 0L && length(aliases) > 0L){
    relevant <- any(parameter_names %in% aliases)
  }else if(length(parameter_names) == 0L && !allow_unlabeled){
    relevant <- FALSE
  }
  if(isTRUE(relevant)){
    relevant <- .posterior_density_condition_matches(
      posterior_density,
      conditional      = .bt_meta_condition(samples, "conditional"),
      conditional_rule = .bt_meta_condition(samples, "conditional_rule"),
      condition_key    = .bt_meta_condition(samples, "condition_key")
    )
  }

  value <- if(isTRUE(relevant)){
    .posterior_density_from_attribute(posterior_density)
  }else{
    NULL
  }

  list(
    present  = TRUE,
    relevant = isTRUE(relevant),
    valid    = !is.null(value),
    value    = value
  )
}

.posterior_density_support_from_attribute <- function(posterior_density){

  density <- .posterior_density_from_attribute(posterior_density)
  if(is.null(density)){
    return(NULL)
  }

  .posterior_support_from_attribute(density[["support"]])
}

.posterior_density_fill_missing_support <- function(posterior_density,
                                                    support_source){

  if(!is.null(.posterior_density_support_from_attribute(posterior_density))){
    return(posterior_density)
  }

  source_support <- .posterior_density_support_from_attribute(support_source)
  if(is.null(source_support) || !is.list(posterior_density)){
    return(posterior_density)
  }

  posterior_density[["support"]] <- source_support
  posterior_density
}

# Filter validated original leaves before the exact-null parser drops metadata.
.posterior_ordinate_matching_attribute <- function(posterior_ordinate, aliases, conditional,
                                                   conditional_rule, condition_key = NULL,
                                                   allow_unlabeled = FALSE){

  if(!.posterior_ordinate_has_data(posterior_ordinate)) return(NULL)
  entries <- .posterior_ordinate_entries(posterior_ordinate)
  matched <- vapply(entries, function(entry){
    parameter_names <- .posterior_density_parameter_metadata(entry)
    if(length(parameter_names) > 0L && length(aliases) > 0L){
      if(!any(parameter_names %in% aliases)) return(FALSE)
    }else if(length(parameter_names) == 0L && !allow_unlabeled){
      return(FALSE)
    }
    .posterior_density_condition_matches(entry, conditional, conditional_rule, condition_key)
  }, logical(1))
  entries <- entries[matched]
  if(length(entries) == 0L) return(NULL)
  if(length(entries) == 1L) return(entries[[1L]])
  structure(list(status = "ok", ordinates = entries),
            class = c("BayesTools_posterior_ordinates", "list"))
}

.posterior_ordinate_candidate_matches <- function(posterior_ordinate, aliases, conditional, conditional_rule,
                                                  condition_key = NULL, allow_unlabeled = FALSE,
                                                  null_hypothesis = NULL){

  matched <- .posterior_ordinate_matching_attribute(
    posterior_ordinate, aliases, conditional, conditional_rule, condition_key, allow_unlabeled
  )
  !is.null(matched) && (is.null(null_hypothesis) ||
    !is.null(.posterior_ordinate_from_attribute(matched, null_hypothesis)))
}

.posterior_ordinate_direct_candidate_matches <- function(posterior_ordinate,
                                                         samples,
                                                         aliases = NULL,
                                                         allow_unlabeled = TRUE,
                                                         null_hypothesis = NULL){

  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }
  .posterior_ordinate_candidate_matches(
    posterior_ordinate,
    aliases          = aliases,
    conditional      = .bt_meta_condition(samples, "conditional"),
    conditional_rule = .bt_meta_condition(samples, "conditional_rule"),
    condition_key    = .bt_meta_condition(samples, "condition_key"),
    allow_unlabeled  = allow_unlabeled,
    null_hypothesis  = null_hypothesis
  )
}

.posterior_ordinate_direct_attribute <- function(samples, aliases = NULL,
                                                 null_hypothesis = NULL){

  posterior_ordinate <- .bt_meta_get(samples, "posterior_ordinate")
  if(is.null(aliases)) aliases <- .posterior_density_sample_aliases(samples)
  matched <- .posterior_ordinate_matching_attribute(
    posterior_ordinate,
    aliases          = aliases,
    conditional      = .bt_meta_condition(samples, "conditional"),
    conditional_rule = .bt_meta_condition(samples, "conditional_rule"),
    condition_key    = .bt_meta_condition(samples, "condition_key"),
    allow_unlabeled  = TRUE
  )
  if(!is.null(matched) && (is.null(null_hypothesis) ||
     !is.null(.posterior_ordinate_from_attribute(matched, null_hypothesis)))){
    return(matched)
  }

  return(NULL)
}

.posterior_ordinate_direct_status <- function(samples, aliases = NULL,
                                              null_hypothesis = NULL,
                                              allow_unlabeled = TRUE){

  posterior_ordinate <- .bt_meta_get(samples, "posterior_ordinate")
  if(identical(.posterior_ordinate_kind(posterior_ordinate), "null")){
    return(list(present = FALSE, relevant = FALSE, valid = FALSE, value = NULL))
  }
  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }

  matched <- .posterior_ordinate_matching_attribute(
    posterior_ordinate,
    aliases          = aliases,
    conditional      = .bt_meta_condition(samples, "conditional"),
    conditional_rule = .bt_meta_condition(samples, "conditional_rule"),
    condition_key    = .bt_meta_condition(samples, "condition_key"),
    allow_unlabeled  = allow_unlabeled
  )
  relevant <- !is.null(matched)

  value <- if(isTRUE(relevant) && !is.null(null_hypothesis)){
    .posterior_ordinate_from_attribute(matched, null_hypothesis)
  }else if(isTRUE(relevant)){
    matched
  }else{
    NULL
  }

  list(
    present  = TRUE,
    relevant = isTRUE(relevant),
    valid    = !is.null(value),
    value    = value
  )
}

# Whether the object is a valid posterior-ordinate attribute (containers are
# not); invalid attributes and unclassed metadata stop with an error.
.posterior_ordinate_has_data <- function(posterior_ordinate){

  kind <- .posterior_ordinate_kind(posterior_ordinate)
  if(identical(kind, "ordinate")){
    .posterior_ordinate_values(posterior_ordinate)
    return(TRUE)
  }
  if(identical(kind, "ordinates")){
    entries <- .posterior_ordinate_entries(posterior_ordinate)
    return(length(entries) > 0L &&
             all(vapply(entries, .posterior_ordinate_has_data, logical(1))))
  }

  return(FALSE)
}
