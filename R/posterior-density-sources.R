.posterior_density_from_sources <- function(sources, aliases, conditional = NULL, conditional_rule = "AND",
                                            condition_key = NULL, allow_unlabeled = FALSE,
                                            null_hypothesis = NULL){

  aliases <- .posterior_density_aliases(aliases)
  if(length(aliases) == 0L){
    return(NULL)
  }

  for(source in sources){
    out <- .posterior_density_from_source(
      source             = source,
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = FALSE,
      depth              = 0L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_density_from_source <- function(source, aliases, conditional, conditional_rule,
                                           condition_key, allow_unlabeled, selected_by_name, depth,
                                           null_hypothesis = NULL){

  if(is.null(source) || depth > 4L){
    return(NULL)
  }

  if(.posterior_density_candidate_matches(
    source,
    aliases          = aliases,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key,
    allow_unlabeled  = allow_unlabeled || selected_by_name,
    null_hypothesis  = null_hypothesis
  )){
    return(source)
  }

  if(!is.list(source)){
    return(NULL)
  }

  container_names <- c("posterior_density", "posterior_densities", "densities")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- .posterior_density_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L,
        null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    for(alias in aliases[aliases %in% source_names]){
      out <- .posterior_density_from_source(
        source             = source[[alias]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = TRUE,
        selected_by_name   = TRUE,
        depth              = depth + 1L,
        null_hypothesis    = null_hypothesis
      )
      if(!is.null(out)){
        return(out)
      }
    }
  }

  for(i in seq_along(source)){
    out <- .posterior_density_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_density_sources <- function(...){

  objects <- list(...)
  sources <- list()
  for(object in objects){
    if(is.null(object)){
      next
    }
    object_sources <- list(
      attr(object, "posterior_density", exact = TRUE),
      attr(object, "posterior_densities", exact = TRUE)
    )
    object_sources <- object_sources[!vapply(object_sources, is.null, logical(1))]
    sources <- c(sources, object_sources)
  }

  return(sources)
}

.posterior_density_attach <- function(samples, sources, parameter, conditional = NULL, conditional_rule = "AND",
                                      condition_key = NULL, allow_unlabeled = FALSE){

  if(length(sources) == 0L){
    return(samples)
  }

  if(is.matrix(samples) || is.data.frame(samples)){
    density_list <- attr(samples, "posterior_density", exact = TRUE)
    if(is.null(density_list) || !is.list(density_list)){
      density_list <- list()
    }

    density_aliases <- list()
    sample_names <- colnames(samples)
    if(!is.null(sample_names)){
      for(sample_name in sample_names){
        density_aliases[[sample_name]] <- .posterior_density_aliases(sample_name)
      }
    }

    level_names <- attr(samples, "level_names", exact = TRUE)
    if(is.list(level_names)){
      level_names <- .factor_cell_labels(level_names)
    }
    factor_cell_names <- attr(samples, "factor_cell_names", exact = TRUE)
    level_names <- .posterior_density_aliases(level_names, factor_cell_names)
    for(level_name in level_names){
      density_aliases[[level_name]] <- .posterior_density_aliases(
        level_name,
        paste0(parameter, "[", level_name, "]")
      )
    }

    if(length(density_aliases) == 0L && ncol(samples) == 1L){
      density_aliases[[parameter]] <- .posterior_density_aliases(parameter)
    }

    for(density_name in names(density_aliases)){
      density <- .posterior_density_from_sources(
        sources            = sources,
        aliases            = density_aliases[[density_name]],
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled && length(density_aliases) == 1L
      )
      if(!is.null(density)){
        density_list[[density_name]] <- density
      }
    }

    if(length(density_list) > 0L){
      attr(samples, "posterior_density") <- density_list
    }

    return(samples)
  }

  density <- .posterior_density_from_sources(
    sources            = sources,
    aliases            = .posterior_density_aliases(parameter, attr(samples, "parameter", exact = TRUE)),
    conditional        = conditional,
    conditional_rule   = conditional_rule,
    condition_key      = condition_key,
    allow_unlabeled    = allow_unlabeled
  )
  if(!is.null(density)){
    attr(samples, "posterior_density") <- density
  }

  return(samples)
}

.posterior_ordinate_sources <- function(...){

  objects <- list(...)
  sources <- list()
  for(object in objects){
    if(is.null(object)){
      next
    }
    object_sources <- list(
      attr(object, "posterior_ordinate", exact = TRUE),
      attr(object, "posterior_ordinates", exact = TRUE)
    )
    object_sources <- object_sources[!vapply(object_sources, is.null, logical(1))]
    sources <- c(sources, object_sources)
  }

  return(sources)
}

.posterior_ordinate_from_sources <- function(sources, aliases, conditional = NULL, conditional_rule = "AND",
                                             condition_key = NULL, allow_unlabeled = FALSE,
                                             null_hypothesis = NULL){

  aliases <- .posterior_density_aliases(aliases)
  if(length(aliases) == 0L){
    return(NULL)
  }

  if(is.null(null_hypothesis)){
    matches <- list()
    for(source in sources){
      matches <- c(matches, .posterior_ordinate_collect_from_source(
        source             = source,
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = FALSE,
        depth              = 0L
      ))
    }
    matches <- .posterior_ordinate_unique_sources(matches)
    if(length(matches) == 0L){
      return(NULL)
    }
    if(length(matches) == 1L){
      return(matches[[1]])
    }

    return(list(ordinates = matches))
  }

  for(source in sources){
    out <- .posterior_ordinate_from_source(
      source             = source,
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = FALSE,
      depth              = 0L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_ordinate_collect_from_source <- function(source, aliases, conditional, conditional_rule,
                                                    condition_key, allow_unlabeled,
                                                    selected_by_name, depth){

  if(is.null(source) || depth > 4L){
    return(list())
  }

  out <- list()
  if(.posterior_ordinate_candidate_matches(
    source,
    aliases          = aliases,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key,
    allow_unlabeled  = allow_unlabeled || selected_by_name
  )){
    out <- c(out, list(source))
  }

  if(!is.list(source)){
    return(out)
  }

  container_names <- c("posterior_ordinate", "posterior_ordinates", "ordinates")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- c(out, .posterior_ordinate_collect_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L
    ))
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    matched_names <- aliases[aliases %in% source_names]
    if(length(matched_names) > 0L){
      matched_count <- length(out)
      for(alias in matched_names){
        out <- c(out, .posterior_ordinate_collect_from_source(
          source             = source[[alias]],
          aliases            = aliases,
          conditional        = conditional,
          conditional_rule   = conditional_rule,
          condition_key      = condition_key,
          allow_unlabeled    = TRUE,
          selected_by_name   = TRUE,
          depth              = depth + 1L
        ))
      }

      if(length(out) > matched_count){
        return(out)
      }
    }
  }

  for(i in seq_along(source)){
    out <- c(out, .posterior_ordinate_collect_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L
    ))
  }

  return(out)
}

.posterior_ordinate_unique_sources <- function(sources){

  if(length(sources) <= 1L){
    return(sources)
  }

  keys <- vapply(sources, function(source) {
    paste(serialize(source, NULL, ascii = TRUE), collapse = "")
  }, character(1))

  return(sources[!duplicated(keys)])
}

.posterior_ordinate_from_source <- function(source, aliases, conditional, conditional_rule,
                                            condition_key, allow_unlabeled, selected_by_name, depth,
                                            null_hypothesis){

  if(is.null(source) || depth > 4L){
    return(NULL)
  }

  has_ordinate <- if(is.null(null_hypothesis)){
    !is.null(source)
  }else{
    !is.null(.posterior_ordinate_from_attribute(source, null_hypothesis))
  }
  if(has_ordinate &&
     .posterior_ordinate_candidate_matches(
       source,
     aliases          = aliases,
     conditional      = conditional,
     conditional_rule = conditional_rule,
     condition_key    = condition_key,
     allow_unlabeled  = allow_unlabeled || selected_by_name,
     null_hypothesis  = null_hypothesis
     )){
    return(source)
  }

  if(!is.list(source)){
    return(NULL)
  }

  container_names <- c("posterior_ordinate", "posterior_ordinates", "ordinates")
  for(container_name in container_names[container_names %in% names(source)]){
    out <- .posterior_ordinate_from_source(
      source             = source[[container_name]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled,
        selected_by_name   = selected_by_name,
        depth              = depth + 1L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  source_names <- names(source)
  if(!is.null(source_names)){
    for(alias in aliases[aliases %in% source_names]){
      out <- .posterior_ordinate_from_source(
        source             = source[[alias]],
        aliases            = aliases,
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = TRUE,
        selected_by_name   = TRUE,
        depth              = depth + 1L,
        null_hypothesis    = null_hypothesis
      )
      if(!is.null(out)){
        return(out)
      }
    }
  }

  for(i in seq_along(source)){
    out <- .posterior_ordinate_from_source(
      source             = source[[i]],
      aliases            = aliases,
      conditional        = conditional,
      conditional_rule   = conditional_rule,
      condition_key      = condition_key,
      allow_unlabeled    = allow_unlabeled,
      selected_by_name   = selected_by_name,
      depth              = depth + 1L,
      null_hypothesis    = null_hypothesis
    )
    if(!is.null(out)){
      return(out)
    }
  }

  return(NULL)
}

.posterior_ordinate_attach <- function(samples, sources, parameter, conditional = NULL, conditional_rule = "AND",
                                       condition_key = NULL, allow_unlabeled = FALSE){

  if(length(sources) == 0L){
    return(samples)
  }

  if(is.matrix(samples) || is.data.frame(samples)){
    ordinate_list <- attr(samples, "posterior_ordinate", exact = TRUE)
    if(is.null(ordinate_list) || !is.list(ordinate_list)){
      ordinate_list <- list()
    }

    ordinate_aliases <- list()
    sample_names <- colnames(samples)
    if(!is.null(sample_names)){
      for(sample_name in sample_names){
        ordinate_aliases[[sample_name]] <- .posterior_density_aliases(sample_name)
      }
    }

    level_names <- attr(samples, "level_names", exact = TRUE)
    if(is.list(level_names)){
      level_names <- .factor_cell_labels(level_names)
    }
    factor_cell_names <- attr(samples, "factor_cell_names", exact = TRUE)
    level_names <- .posterior_density_aliases(level_names, factor_cell_names)
    for(level_name in level_names){
      ordinate_aliases[[level_name]] <- .posterior_density_aliases(
        level_name,
        paste0(parameter, "[", level_name, "]")
      )
    }

    if(length(ordinate_aliases) == 0L && ncol(samples) == 1L){
      ordinate_aliases[[parameter]] <- .posterior_density_aliases(parameter)
    }

    for(ordinate_name in names(ordinate_aliases)){
      ordinate <- .posterior_ordinate_from_sources(
        sources            = sources,
        aliases            = ordinate_aliases[[ordinate_name]],
        conditional        = conditional,
        conditional_rule   = conditional_rule,
        condition_key      = condition_key,
        allow_unlabeled    = allow_unlabeled && length(ordinate_aliases) == 1L
      )
      if(!is.null(ordinate)){
        ordinate_list[[ordinate_name]] <- ordinate
      }
    }

    if(length(ordinate_list) > 0L){
      attr(samples, "posterior_ordinate") <- ordinate_list
    }

    return(samples)
  }

  ordinate <- .posterior_ordinate_from_sources(
    sources            = sources,
    aliases            = .posterior_density_aliases(parameter, attr(samples, "parameter", exact = TRUE)),
    conditional        = conditional,
    conditional_rule   = conditional_rule,
    condition_key      = condition_key,
    allow_unlabeled    = allow_unlabeled
  )
  if(!is.null(ordinate)){
    attr(samples, "posterior_ordinate") <- ordinate
  }

  return(samples)
}
