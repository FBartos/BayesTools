.posterior_density_child_attributes <- function(samples, null_hypothesis = NULL){

  sample_names <- names(samples)
  out <- lapply(seq_along(samples), function(i) {
    sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
    .posterior_density_direct_attribute(
      samples[[i]],
      aliases         = .posterior_density_child_aliases(samples, samples[[i]], sample_name),
      null_hypothesis = null_hypothesis
    )
  })
  top_sources <- list(
    .bt_meta_get(samples, "posterior_density"),
    .bt_meta_get(samples, "posterior_densities")
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
      aliases <- .posterior_density_child_aliases(samples, samples[[i]], sample_name)
      conditional      <- .bt_meta_condition(samples[[i]], "conditional")
      conditional_rule <- .bt_meta_condition(samples[[i]], "conditional_rule")
      condition_key    <- .bt_meta_condition(samples[[i]], "condition_key")
      density <- .posterior_density_from_sources(
        sources          = top_sources,
        aliases          = aliases,
        conditional      = conditional,
        conditional_rule = conditional_rule,
        condition_key    = condition_key,
        allow_unlabeled  = FALSE,
        null_hypothesis  = null_hypothesis
      )
      if(is.null(out[[i]]) && !is.null(density)){
        out[[i]] <- density
      }else if(!is.null(out[[i]]) && !is.null(density)){
        out[[i]] <- .posterior_density_fill_missing_support(out[[i]], density)
      }
    }
  }
  for(top_level in top_sources){
    if(is.list(top_level) && is.null(names(top_level)) && length(top_level) == length(samples)){
      missing <- vapply(out, is.null, logical(1))
      if(any(missing)){
        sample_names <- names(samples)
        for(i in which(missing)){
          aliases <- .posterior_density_child_aliases(
            samples,
            samples[[i]],
            if(!is.null(sample_names)) sample_names[[i]] else NULL
          )
          conditional      <- .bt_meta_condition(samples[[i]], "conditional")
          conditional_rule <- .bt_meta_condition(samples[[i]], "conditional_rule")
          condition_key    <- .bt_meta_condition(samples[[i]], "condition_key")
          if(.posterior_density_candidate_matches(
            top_level[[i]],
            aliases          = aliases,
            conditional      = conditional,
            conditional_rule = conditional_rule,
            condition_key    = condition_key,
            allow_unlabeled  = TRUE,
            null_hypothesis  = null_hypothesis
          )){
            if(is.null(out[[i]])){
              out[[i]] <- top_level[[i]]
            }else{
              out[[i]] <- .posterior_density_fill_missing_support(
                out[[i]],
                top_level[[i]]
              )
            }
          }
        }
        names(out) <- sample_names
      }
    }
  }

  return(out)
}

.posterior_ordinate_for_method <- function(posterior_ordinate, null_hypothesis,
                                           density_method){

  if(identical(density_method, "precomputed")){
    return(.posterior_ordinate_from_attribute(posterior_ordinate, null_hypothesis))
  }

  return(NULL)
}

.posterior_ordinate_child_attributes <- function(samples, null_hypothesis = NULL){

  sample_names <- names(samples)
  out <- lapply(seq_along(samples), function(i) {
    sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
    .posterior_ordinate_direct_attribute(
      samples[[i]],
      aliases         = .posterior_density_child_aliases(samples, samples[[i]], sample_name),
      null_hypothesis = null_hypothesis
    )
  })
  top_sources <- list(
    .bt_meta_get(samples, "posterior_ordinate"),
    .bt_meta_get(samples, "posterior_ordinates")
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
      aliases <- .posterior_density_child_aliases(samples, samples[[i]], sample_name)
      conditional      <- .bt_meta_condition(samples[[i]], "conditional")
      conditional_rule <- .bt_meta_condition(samples[[i]], "conditional_rule")
      condition_key    <- .bt_meta_condition(samples[[i]], "condition_key")
      ordinate <- .posterior_ordinate_from_sources(
        sources          = top_sources,
        aliases          = aliases,
        conditional      = conditional,
        conditional_rule = conditional_rule,
        condition_key    = condition_key,
        allow_unlabeled  = FALSE,
        null_hypothesis  = null_hypothesis
      )
      if(is.null(out[[i]]) && !is.null(ordinate)){
        out[[i]] <- ordinate
      }
    }
  }
  for(top_level in top_sources){
    if(is.list(top_level) && is.null(names(top_level)) && length(top_level) == length(samples)){
      missing <- vapply(out, is.null, logical(1))
      if(any(missing)){
        sample_names <- names(samples)
        for(i in which(missing)){
          aliases <- .posterior_density_child_aliases(
            samples,
            samples[[i]],
            if(!is.null(sample_names)) sample_names[[i]] else NULL
          )
          conditional      <- .bt_meta_condition(samples[[i]], "conditional")
          conditional_rule <- .bt_meta_condition(samples[[i]], "conditional_rule")
          condition_key    <- .bt_meta_condition(samples[[i]], "condition_key")
          if(.posterior_ordinate_candidate_matches(
            top_level[[i]],
            aliases          = aliases,
            conditional      = conditional,
            conditional_rule = conditional_rule,
            condition_key    = condition_key,
            allow_unlabeled  = TRUE,
            null_hypothesis  = null_hypothesis
          )){
            out[[i]] <- .posterior_ordinate_matching_attribute(
              top_level[[i]], aliases, conditional, conditional_rule, condition_key,
              allow_unlabeled = TRUE
            )
          }
        }
        names(out) <- sample_names
      }
    }
  }

  return(out)
}

.posterior_precomputed_child <- function(parent, child, index, null_hypothesis,
                                         density_method){

  if(!identical(density_method, "precomputed") || is.null(parent) ||
     is.null(index)){
    return(child)
  }

  parent_names <- names(parent)
  if(is.character(index)){
    if(is.null(parent_names) || !index %in% parent_names){
      return(child)
    }
    child_index <- match(index, parent_names)
  }else{
    child_index <- as.integer(index)
    if(length(child_index) != 1L || is.na(child_index) ||
       child_index < 1L || child_index > length(parent)){
      return(child)
    }
  }

  child_name <- if(!is.null(parent_names)) parent_names[[child_index]] else NULL
  aliases <- .posterior_density_child_aliases(parent, child, child_name)

  child_density <- .posterior_density_for_method(
    .posterior_density_direct_attribute(
      child,
      aliases         = aliases,
      null_hypothesis = null_hypothesis
    ),
    density_method
  )
  if(is.null(child_density)){
    parent_densities <- .posterior_density_child_attributes(
      parent,
      null_hypothesis = null_hypothesis
    )
    if(length(parent_densities) >= child_index &&
       !is.null(parent_densities[[child_index]])){
      child <- .bt_meta_set(child, "posterior_density", parent_densities[[child_index]])
    }
  }

  child_ordinate <- .posterior_ordinate_for_method(
    .posterior_ordinate_direct_attribute(
      child,
      aliases         = aliases,
      null_hypothesis = null_hypothesis
    ),
    null_hypothesis,
    density_method
  )
  if(is.null(child_ordinate)){
    parent_ordinates <- .posterior_ordinate_child_attributes(
      parent,
      null_hypothesis = null_hypothesis
    )
    if(length(parent_ordinates) >= child_index &&
       !is.null(parent_ordinates[[child_index]])){
      child <- .bt_meta_set(child, "posterior_ordinate", parent_ordinates[[child_index]])
    }
  }

  child
}

.posterior_ordinate_from_attribute <- function(posterior_ordinate,
                                               null_hypothesis){

  kind <- .posterior_ordinate_kind(posterior_ordinate)
  if(identical(kind, "ordinates")){
    candidates <- lapply(
      .posterior_ordinate_entries(posterior_ordinate),
      .posterior_ordinate_from_attribute,
      null_hypothesis = null_hypothesis
    )
    matched <- !vapply(candidates, is.null, logical(1))
    if(sum(matched) != 1L){
      return(NULL)
    }

    return(candidates[[which(matched)]])
  }
  if(!identical(kind, "ordinate")){
    return(NULL)
  }

  values <- .posterior_ordinate_values(posterior_ordinate)
  index  <- which(values[["value"]] == null_hypothesis)
  if(length(index) != 1L){
    return(NULL)
  }

  return(list(
    x           = values[["value"]][index],
    y           = values[["ordinate"]][index],
    method      = posterior_ordinate[["method"]],
    diagnostics = .posterior_ordinate_subset_diagnostics(
      posterior_ordinate[["diagnostics"]],
      index
    )
  ))
}

# Validated value and ordinate vectors of one posterior-ordinate attribute.
.posterior_ordinate_values <- function(posterior_ordinate){

  value    <- posterior_ordinate[["value"]]
  ordinate <- posterior_ordinate[["ordinate"]]
  if(!identical(posterior_ordinate[["status"]], "ok") ||
     !is.numeric(value) || !is.numeric(ordinate) ||
     length(value) == 0L || length(value) != length(ordinate) ||
     any(!is.finite(value)) || any(!is.finite(ordinate)) ||
     any(ordinate <= 0) || anyDuplicated(value)){
    stop(
      "Posterior ordinate metadata is invalid: it needs status 'ok' and unique, ",
      "finite 'value' entries with finite, positive 'ordinate' heights.",
      call. = FALSE
    )
  }

  list(value = as.numeric(value), ordinate = as.numeric(ordinate))
}

.posterior_ordinate_subset_diagnostics <- function(diagnostics, index){

  if(is.null(diagnostics) || !is.list(diagnostics)){
    return(diagnostics)
  }

  if(is.data.frame(diagnostics)){
    diagnostics <- as.list(diagnostics[index, , drop = FALSE])
  }

  for(name in names(diagnostics)){
    value <- diagnostics[[name]]
    if(is.atomic(value) && length(value) > 1L && length(value) >= index){
      diagnostics[[name]] <- value[index]
    }
  }

  return(diagnostics)
}
