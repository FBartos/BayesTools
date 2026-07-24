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

  point_masses <- posterior_density[["point_masses"]]
  if(!is.null(point_masses) && nrow(point_masses) > 0L){
    point_tol <- sqrt(.Machine$double.eps) * max(1, abs(null_hypothesis))
    if(any(abs(point_masses[["x"]] - null_hypothesis) <= point_tol)){
      return(TRUE)
    }
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

.posterior_density_metadata_values <- function(posterior_density, fields){

  values <- NULL
  if(is.list(posterior_density)){
    for(field in fields){
      if(!is.null(posterior_density[[field]])){
        values <- c(values, unlist(posterior_density[[field]], use.names = FALSE))
      }
    }
    for(container_name in c("density", "ordinate", "ordinates")){
      if(!is.null(posterior_density[[container_name]]) &&
         is.list(posterior_density[[container_name]])){
        for(field in fields){
          if(!is.null(posterior_density[[container_name]][[field]])){
            values <- c(values, unlist(posterior_density[[container_name]][[field]], use.names = FALSE))
          }
        }
      }
    }
  }
  if(is.data.frame(posterior_density)){
    for(field in fields){
      if(field %in% colnames(posterior_density)){
        values <- c(values, unique(posterior_density[[field]]))
      }
    }
  }

  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]

  return(unique(values))
}

.posterior_density_parameter_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(
    posterior_density,
    c("parameter", "parameters", "parameter_name", "name")
  ))
}

.posterior_density_conditional_metadata <- function(posterior_density){

  return(.posterior_density_metadata_values(
    posterior_density,
    c("conditional", "condition", "conditioned_on", "conditioning")
  ))
}

.posterior_density_conditional_rule_metadata <- function(posterior_density){

  rule <- .posterior_density_metadata_values(
    posterior_density,
    c("conditional_rule", "condition_rule", "conditioning_rule")
  )
  if(length(rule) == 0L){
    return(NULL)
  }

  return(rule[[1]])
}

.posterior_density_condition_key_metadata <- function(posterior_density){

  key <- .posterior_density_metadata_values(
    posterior_density,
    c("condition_key", "conditional_key")
  )
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
    conditional      = attr(samples, "conditional", exact = TRUE),
    conditional_rule = attr(samples, "conditional_rule", exact = TRUE),
    condition_key    = attr(samples, "condition_key", exact = TRUE)
  ))
}

.posterior_density_direct_attribute <- function(samples, aliases = NULL,
                                                null_hypothesis = NULL){

  posterior_density <- attr(samples, "posterior_density", exact = TRUE)
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

  posterior_density <- attr(samples, "posterior_density", exact = TRUE)
  if(is.null(posterior_density)){
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
      conditional      = attr(samples, "conditional", exact = TRUE),
      conditional_rule = attr(samples, "conditional_rule", exact = TRUE),
      condition_key    = attr(samples, "condition_key", exact = TRUE)
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

.posterior_ordinate_candidate_matches <- function(posterior_ordinate, aliases, conditional, conditional_rule,
                                                  condition_key = NULL, allow_unlabeled = FALSE,
                                                  null_hypothesis = NULL){

  if(!.posterior_ordinate_has_data(posterior_ordinate)){
    return(FALSE)
  }
  if(!is.null(null_hypothesis) &&
     is.null(.posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis))){
    return(FALSE)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_ordinate)
  if(length(parameter_names) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(!allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_ordinate,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key
  ))
}

.posterior_ordinate_direct_candidate_matches <- function(posterior_ordinate,
                                                         samples,
                                                         aliases = NULL,
                                                         allow_unlabeled = TRUE,
                                                         null_hypothesis = NULL){

  if(!.posterior_ordinate_has_data(posterior_ordinate)){
    return(FALSE)
  }
  if(!is.null(null_hypothesis) &&
     is.null(.posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis))){
    return(FALSE)
  }
  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_ordinate)
  if(length(parameter_names) > 0L && length(aliases) > 0L){
    if(!any(parameter_names %in% aliases)){
      return(FALSE)
    }
  }else if(length(parameter_names) == 0L && !allow_unlabeled){
    return(FALSE)
  }

  return(.posterior_density_condition_matches(
    posterior_ordinate,
    conditional      = attr(samples, "conditional", exact = TRUE),
    conditional_rule = attr(samples, "conditional_rule", exact = TRUE),
    condition_key    = attr(samples, "condition_key", exact = TRUE)
  ))
}

.posterior_ordinate_direct_attribute <- function(samples, aliases = NULL,
                                                 null_hypothesis = NULL){

  posterior_ordinate <- attr(samples, "posterior_ordinate", exact = TRUE)
  if(.posterior_ordinate_direct_candidate_matches(
    posterior_ordinate,
    samples         = samples,
    aliases         = aliases,
    null_hypothesis = null_hypothesis
  )){
    return(posterior_ordinate)
  }

  return(NULL)
}

.posterior_ordinate_direct_status <- function(samples, aliases = NULL,
                                              null_hypothesis = NULL,
                                              allow_unlabeled = TRUE){

  posterior_ordinate <- attr(samples, "posterior_ordinate", exact = TRUE)
  if(is.null(posterior_ordinate)){
    return(list(present = FALSE, relevant = FALSE, valid = FALSE, value = NULL))
  }
  if(is.null(aliases)){
    aliases <- .posterior_density_sample_aliases(samples)
  }

  parameter_names <- .posterior_density_parameter_metadata(posterior_ordinate)
  relevant <- TRUE
  if(length(parameter_names) > 0L && length(aliases) > 0L){
    relevant <- any(parameter_names %in% aliases)
  }else if(length(parameter_names) == 0L && !allow_unlabeled){
    relevant <- FALSE
  }
  if(isTRUE(relevant)){
    relevant <- .posterior_density_condition_matches(
      posterior_ordinate,
      conditional      = attr(samples, "conditional", exact = TRUE),
      conditional_rule = attr(samples, "conditional_rule", exact = TRUE),
      condition_key    = attr(samples, "condition_key", exact = TRUE)
    )
  }

  value <- if(isTRUE(relevant) && !is.null(null_hypothesis)){
    .posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis)
  }else if(isTRUE(relevant) && .posterior_ordinate_has_data(posterior_ordinate)){
    posterior_ordinate
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

.posterior_ordinate_has_data <- function(posterior_ordinate){

  if(is.null(posterior_ordinate)){
    return(FALSE)
  }
  if(is.data.frame(posterior_ordinate)){
    return(
      any(c("x", "value", "null_hypothesis") %in% colnames(posterior_ordinate)) &&
        any(c("y", "ordinate", "height", "posterior_height") %in% colnames(posterior_ordinate))
    )
  }
  if(!is.list(posterior_ordinate)){
    return(FALSE)
  }
  if(!is.null(posterior_ordinate[["status"]]) &&
     !identical(posterior_ordinate[["status"]], "ok")){
    return(FALSE)
  }
  if((!is.null(posterior_ordinate[["x"]]) ||
      !is.null(posterior_ordinate[["value"]]) ||
      !is.null(posterior_ordinate[["null_hypothesis"]])) &&
     (!is.null(posterior_ordinate[["y"]]) ||
      !is.null(posterior_ordinate[["ordinate"]]) ||
      !is.null(posterior_ordinate[["height"]]) ||
      !is.null(posterior_ordinate[["posterior_height"]]))){
    return(TRUE)
  }
  if(!is.null(posterior_ordinate[["ordinate"]]) &&
     (is.list(posterior_ordinate[["ordinate"]]) ||
      is.data.frame(posterior_ordinate[["ordinate"]]))){
    return(.posterior_ordinate_has_data(posterior_ordinate[["ordinate"]]))
  }
  if(!is.null(posterior_ordinate[["ordinates"]]) &&
     (is.list(posterior_ordinate[["ordinates"]]) ||
      is.data.frame(posterior_ordinate[["ordinates"]]))){
    if(.posterior_ordinate_has_data(posterior_ordinate[["ordinates"]])){
      return(TRUE)
    }
    if(is.list(posterior_ordinate[["ordinates"]])){
      return(any(vapply(
        posterior_ordinate[["ordinates"]],
        .posterior_ordinate_has_data,
        logical(1)
      )))
    }
  }

  return(FALSE)
}
