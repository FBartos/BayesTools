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
    attr(samples, "posterior_density", exact = TRUE),
    attr(samples, "posterior_densities", exact = TRUE)
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
      aliases <- .posterior_density_child_aliases(samples, samples[[i]], sample_name)
      conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
      conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
      condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
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
          conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
          conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
          condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
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
    attr(samples, "posterior_ordinate", exact = TRUE),
    attr(samples, "posterior_ordinates", exact = TRUE)
  )
  top_sources <- top_sources[!vapply(top_sources, is.null, logical(1))]
  if(length(top_sources) > 0L){
    sample_names <- names(samples)
    for(i in seq_along(samples)){
      sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
      aliases <- .posterior_density_child_aliases(samples, samples[[i]], sample_name)
      conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
      conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
      condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
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
          conditional      <- attr(samples[[i]], "conditional", exact = TRUE)
          conditional_rule <- attr(samples[[i]], "conditional_rule", exact = TRUE)
          condition_key    <- attr(samples[[i]], "condition_key", exact = TRUE)
          if(.posterior_ordinate_candidate_matches(
            top_level[[i]],
            aliases          = aliases,
            conditional      = conditional,
            conditional_rule = conditional_rule,
            condition_key    = condition_key,
            allow_unlabeled  = TRUE,
            null_hypothesis  = null_hypothesis
          )){
            out[[i]] <- top_level[[i]]
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
      attr(child, "posterior_density") <- parent_densities[[child_index]]
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
      attr(child, "posterior_ordinate") <- parent_ordinates[[child_index]]
    }
  }

  child
}

.posterior_ordinate_from_attribute <- function(posterior_ordinate,
                                               null_hypothesis){

  if(is.null(posterior_ordinate)){
    return(NULL)
  }

  if(is.list(posterior_ordinate) &&
     !is.null(posterior_ordinate[["status"]]) &&
     !identical(posterior_ordinate[["status"]], "ok")){
    return(NULL)
  }
  if(is.list(posterior_ordinate) &&
     !is.data.frame(posterior_ordinate[["ordinates"]]) &&
     is.list(posterior_ordinate[["ordinates"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["x"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["value"]]) &&
     is.null(posterior_ordinate[["ordinates"]][["null_hypothesis"]])){
    candidates <- lapply(
      posterior_ordinate[["ordinates"]],
      .posterior_ordinate_from_attribute,
      null_hypothesis = null_hypothesis
    )
    matched <- !vapply(candidates, is.null, logical(1))
    if(sum(matched) != 1L){
      return(NULL)
    }

    out <- candidates[[which(matched)]]
    if(is.null(out[["method"]])){
      if(!is.null(posterior_ordinate[["method"]])){
        out[["method"]] <- posterior_ordinate[["method"]]
      }else if(!is.null(posterior_ordinate[["estimator"]])){
        out[["method"]] <- posterior_ordinate[["estimator"]]
      }
    }
    if(!is.null(posterior_ordinate[["diagnostics"]])){
      diagnostics <- .posterior_ordinate_subset_diagnostics(
        posterior_ordinate[["diagnostics"]],
        which(matched)
      )
      out[["diagnostics"]] <- .posterior_ordinate_merge_diagnostics(
        diagnostics,
        out[["diagnostics"]]
      )
    }

    return(out)
  }

  source      <- posterior_ordinate
  method      <- NULL
  diagnostics <- NULL

  if(is.list(posterior_ordinate)){
    if(!is.null(posterior_ordinate[["method"]])){
      method <- posterior_ordinate[["method"]]
    }
    if(!is.null(posterior_ordinate[["estimator"]])){
      method <- posterior_ordinate[["estimator"]]
    }
    if(!is.null(posterior_ordinate[["diagnostics"]])){
      diagnostics <- posterior_ordinate[["diagnostics"]]
    }
    if(!is.null(posterior_ordinate[["ordinate"]]) &&
       (is.list(posterior_ordinate[["ordinate"]]) ||
        is.data.frame(posterior_ordinate[["ordinate"]]))){
      source <- posterior_ordinate[["ordinate"]]
    }
    if(!is.null(posterior_ordinate[["ordinates"]]) &&
       (is.list(posterior_ordinate[["ordinates"]]) ||
        is.data.frame(posterior_ordinate[["ordinates"]]))){
      source <- posterior_ordinate[["ordinates"]]
    }
  }

  source_diagnostics <- NULL
  if(is.data.frame(source)){
    x_name <- intersect(c("x", "value", "null_hypothesis"), colnames(source))[1]
    y_name <- intersect(c("y", "ordinate", "height", "posterior_height"), colnames(source))[1]
    if(is.na(x_name) || is.na(y_name)){
      return(NULL)
    }
    x <- source[[x_name]]
    y <- source[[y_name]]
    source_diagnostics <- source[, setdiff(colnames(source), c(x_name, y_name)), drop = FALSE]
  }else if(is.list(source)){
    if(is.null(method) && !is.null(source[["method"]])){
      method <- source[["method"]]
    }
    if(is.null(method) && !is.null(source[["estimator"]])){
      method <- source[["estimator"]]
    }
    if(!is.null(source[["x"]])){
      x <- source[["x"]]
    }else if(!is.null(source[["value"]])){
      x <- source[["value"]]
    }else if(!is.null(source[["null_hypothesis"]])){
      x <- source[["null_hypothesis"]]
    }else{
      return(NULL)
    }

    if(!is.null(source[["y"]])){
      y <- source[["y"]]
    }else if(!is.null(source[["ordinate"]])){
      y <- source[["ordinate"]]
    }else if(!is.null(source[["height"]])){
      y <- source[["height"]]
    }else if(!is.null(source[["posterior_height"]])){
      y <- source[["posterior_height"]]
    }else{
      return(NULL)
    }
    source_diagnostics <- source[setdiff(
      names(source),
      c("x", "value", "null_hypothesis", "y", "ordinate", "height",
        "posterior_height", "method", "estimator", "parameter", "parameters",
        "parameter_name", "name", "conditional", "condition", "conditioned_on",
        "conditioning", "conditional_rule", "condition_rule",
        "conditioning_rule", "condition_key", "conditional_key")
    )]
  }else{
    return(NULL)
  }

  x <- as.numeric(x)
  y <- as.numeric(y)
  if(length(x) != length(y)){
    return(NULL)
  }

  keep <- is.finite(x) & is.finite(y) & y > 0
  x <- x[keep]
  y <- y[keep]
  diagnostics <- .posterior_ordinate_subset_keep(diagnostics, keep)
  source_diagnostics <- .posterior_ordinate_subset_keep(source_diagnostics, keep)
  if(length(x) == 0L){
    return(NULL)
  }

  index <- which(x == null_hypothesis)
  if(length(index) != 1L){
    return(NULL)
  }
  if(is.null(diagnostics)){
    diagnostics <- source_diagnostics
  }
  diagnostics        <- .posterior_ordinate_subset_diagnostics(diagnostics, index)
  source_diagnostics <- .posterior_ordinate_subset_diagnostics(source_diagnostics, index)
  diagnostics        <- .posterior_ordinate_merge_diagnostics(diagnostics, source_diagnostics)

  return(list(
    x           = x[index],
    y           = y[index],
    method      = method,
    diagnostics = diagnostics
  ))
}

.posterior_ordinate_subset_keep <- function(diagnostics, keep){

  if(is.null(diagnostics)){
    return(NULL)
  }
  if(is.data.frame(diagnostics)){
    return(diagnostics[keep, , drop = FALSE])
  }
  if(!is.list(diagnostics)){
    return(diagnostics)
  }

  for(name in names(diagnostics)){
    value <- diagnostics[[name]]
    if(is.atomic(value) && length(value) == length(keep)){
      diagnostics[[name]] <- value[keep]
    }
  }

  return(diagnostics)
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

.posterior_ordinate_merge_diagnostics <- function(diagnostics, source_diagnostics){

  if(is.null(diagnostics)){
    return(source_diagnostics)
  }
  if(is.null(source_diagnostics)){
    return(diagnostics)
  }
  if(!is.list(diagnostics)){
    return(source_diagnostics)
  }
  if(!is.list(source_diagnostics)){
    return(diagnostics)
  }

  for(name in names(source_diagnostics)){
    diagnostics[[name]] <- source_diagnostics[[name]]
  }

  return(diagnostics)
}
