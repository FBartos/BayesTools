random_effects_vignette_cache_schema <- function(){
  c(
    stan_correlated = "stanreg",
    stan_correlated_scaled = "stanreg",
    fit_sleep = "BayesTools_fit",
    fit_sleep_independent = "BayesTools_fit",
    fit_cake_recipe = "BayesTools_fit",
    fit_cake_hcs = "BayesTools_fit",
    fit_cake_ar1 = "BayesTools_fit",
    fit_cake_har = "BayesTools_fit",
    fit_sleep_car = "BayesTools_fit",
    fit_nested = "BayesTools_fit",
    stan_crossed = "stanreg",
    fit_crossed_independent = "BayesTools_fit",
    fit_crossed_allocation = "BayesTools_fit"
  )
}

random_effects_vignette_cache_names <- function(){
  names(random_effects_vignette_cache_schema())
}

validate_random_effects_vignette_cache <- function(cache_file = file.path("models", "RandomEffects.RDS"),
                                                  load = FALSE){
  schema <- random_effects_vignette_cache_schema()
  status <- list(
    file_exists = file.exists(cache_file),
    read_error = NULL,
    is_list = FALSE,
    missing_objects = names(schema),
    invalid_objects = names(schema),
    valid = FALSE,
    cache = NULL
  )

  if(!status$file_exists){
    return(status)
  }

  cache <- tryCatch(readRDS(cache_file), error = function(e) e)
  if(inherits(cache, "error")){
    status$read_error <- conditionMessage(cache)
    return(status)
  }

  status$is_list <- is.list(cache)
  if(!status$is_list){
    return(status)
  }

  status$missing_objects <- setdiff(names(schema), names(cache))
  present <- intersect(names(schema), names(cache))
  status$invalid_objects <- present[
    !vapply(present, function(name) inherits(cache[[name]], schema[[name]]), logical(1))
  ]
  status$valid <- length(status$missing_objects) == 0L &&
    length(status$invalid_objects) == 0L
  if(isTRUE(load)){
    status$cache <- cache
  }

  status
}

format_random_effects_vignette_cache_error <- function(status){
  if(!isTRUE(status$file_exists)){
    return("Missing precomputed RandomEffects vignette cache.")
  }
  if(!is.null(status$read_error)){
    return(paste0("Could not read precomputed RandomEffects vignette cache: ", status$read_error))
  }
  if(!isTRUE(status$is_list)){
    return("Precomputed RandomEffects vignette cache must be a named list.")
  }
  if(length(status$missing_objects) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache is missing objects: ",
      paste(status$missing_objects, collapse = ", ")
    ))
  }
  if(length(status$invalid_objects) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache has unexpected object classes: ",
      paste(status$invalid_objects, collapse = ", ")
    ))
  }

  "Precomputed RandomEffects vignette cache is invalid."
}

stop_if_invalid_random_effects_vignette_cache <- function(cache_file = file.path("models", "RandomEffects.RDS")){
  status <- validate_random_effects_vignette_cache(cache_file)
  if(!isTRUE(status$valid)){
    stop(format_random_effects_vignette_cache_error(status), call. = FALSE)
  }
  invisible(status)
}
