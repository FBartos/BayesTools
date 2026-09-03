.format_column       <- function(x, type, n_models){
  if(is.null(x) || length(x) == 0){
    return(x)
  }else{
    return(switch(
      type,
      "integer"         = round(x),
      "prior"           = .center_priors(x),
      "string_left"     = .string_left(x),
      "string"          = x,
      "hypothesis_label" = x,
      "estimate"        = format(round(x, digits = 3), nsmall = 3),
      "prior_prob"      = format(round(x, digits = 3), nsmall = 3),
      "post_prob"       = format(round(x, digits = 3), nsmall = 3),
      "probability"     = format(round(x, digits = 3), nsmall = 3),
      "marglik"         = format(round(x, digits = 2), nsmall = 2),
      "BF"              = .format_BF_column(x),
      "inclusion_BF"    = .format_BF_column(x),
      "BF_error"        = format(round(x, digits = 3), nsmall = 3),
      "n_models"        = paste0(round(x), "/", n_models),
      "ESS"             = round(x),
      "R_hat"           = format(round(x, digits = 3), nsmall = 3),
      "MCMC_error"      = format(round(x, digits = 5), nsmall = 5),
      "MCMC_SD_error"   = format(round(x, digits = 3), nsmall = 3),
      "min_ESS"             = round(x),
      "max_R_hat"           = format(round(x, digits = 3), nsmall = 3),
      "max_MCMC_error"      = format(round(x, digits = 5), nsmall = 5),
      "max_MCMC_SD_error"   = format(round(x, digits = 3), nsmall = 3)
    ))
  }
}
.format_column_names <- function(x, type, values){
  if(is.null(x) || length(x) == 0){
    return(x)
  }else{
    return(switch(
      type,
      "integer"         = x,
      "prior"           = .string_center(paste0("Prior ", x), values),
      "string_left"     = .string_left(x, values),
      "string"          = x,
      "hypothesis_label" = paste0(x, ":"),
      "estimate"        = x,
      "probability"     = x,
      "prior_prob"      = "Prior prob.",
      "post_prob"       = "Post. prob.",
      "marglik"         = "log(marglik)",
      "BF"              = if(is.null(attr(values, "name"))) "BF"           else attr(values, "name"),
      "inclusion_BF"    = if(is.null(attr(values, "name"))) "Inclusion BF" else attr(values, "name"),
      "BF_error"        = if(is.null(attr(values, "name"))) x              else attr(values, "name"),
      "n_models"        = "Models",
      "ESS"             = "ESS",
      "R_hat"           = "R-hat",
      "MCMC_error"      = if(is.null(attr(values, "name"))) "error(MCMC)" else attr(values, "name"),
      "MCMC_SD_error"   = "error(MCMC)/SD",
      "min_ESS"             = "min(ESS)",
      "max_R_hat"           = "max(R-hat)",
      "max_MCMC_error"      = "max[error(MCMC)]",
      "max_MCMC_SD_error"   = "max[error(MCMC)/SD]",
    ))
  }
}
.center_priors <- function(x){

  if(any(grepl("~", x) | grepl("=", x))){

    position_tilda  <- regexpr("~", x)
    position_equal  <- regexpr("=", x)
    from_right      <- sapply(seq_along(x), function(i){
      if(position_tilda[i] != -1){
        return(nchar(x[i]) - position_tilda[i])
      }else if(position_equal[i] != -1){
        return(nchar(x[i]) - position_equal[i])
      }else{
        return(0)
      }
    })
    add_to_right  <- ifelse(from_right == 0, 0, max(from_right) - from_right)
    x <- paste0(x, sapply(seq_along(x), function(i)paste0(rep(" ", add_to_right[i]), collapse = "")))

  }

  return(x)
}
.string_left   <- function(x, reference = x){

  if(length(x) > 0){

    add_to_right <- max(nchar(reference)) - nchar(x)
    add_to_right <- ifelse(add_to_right < 0, 0, add_to_right)
    x <- paste0(x, sapply(seq_along(x), function(i)paste0(rep(" ", add_to_right[i]), collapse = "")))

  }

  return(x)
}
.string_center <- function(x, reference = x){

  if(length(x) > 0){

    add_to_sides <- max(nchar(reference)) - nchar(x)
    add_to_sides <- ifelse(add_to_sides < 0, 0, add_to_sides)
    x <- paste0( sapply(seq_along(x), function(i)paste0(rep(" ", round(add_to_sides[i]/2)), collapse = "")), x, sapply(seq_along(x), function(i)paste0(rep(" ", round(add_to_sides[i]/2)), collapse = "")))

  }

  return(x)
}

.check_table_types      <- function(x, name, allow_NULL = FALSE){
  check_char(x, name, allow_values = c(
    "integer", "prior", "string_left", "string", "hypothesis_label", "estimate",
    "probability", "prior_prob", "post_prob",
    "marglik", "BF", "inclusion_BF", "BF_error", "n_models",
    "ESS", "R_hat", "MCMC_error", "MCMC_SD_error", "min_ESS", "max_R_hat", "max_MCMC_error", "max_MCMC_SD_error"),
    allow_NULL = allow_NULL)
}
.JAGS_estimates_diagnostic_columns <- function(){
  c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
}
.JAGS_BF_diagnostic_columns <- function(){
  c("ESS", "MCMC_error", "BF_error_percent")
}
.JAGS_BF_diagnostic_column_types <- function(columns){
  unname(c(
    ESS              = "ESS",
    MCMC_error       = "MCMC_error",
    BF_error_percent = "BF_error"
  )[columns])
}
.normalize_diagnostic_columns <- function(columns, allowed_columns, name){

  if(is.null(columns) || length(columns) == 0){
    return(character())
  }

  if(is.logical(columns)){
    check_bool(columns, name, allow_NA = FALSE)
    if(columns){
      return(allowed_columns)
    }else{
      return(character())
    }
  }

  check_char(columns, name, check_length = 0, allow_NA = FALSE)

  shortcut_columns <- c("all", "none")
  if(any(columns %in% shortcut_columns)){
    if(length(columns) > 1){
      stop(paste0("The '", name, "' argument can use 'all' or 'none' only by itself."), call. = FALSE)
    }
    if(columns == "all"){
      return(allowed_columns)
    }else{
      return(character())
    }
  }

  check_char(columns, name, check_length = 0, allow_values = allowed_columns, allow_NA = FALSE)

  return(unique(columns))
}

.runjags_summary_fast   <- function(model_samples, n_samples, n_chains, conditional, probs = c(0.025, 0.975), remove_diagnostics = FALSE,
                                    diagnostic_columns = .JAGS_estimates_diagnostic_columns()){

  diagnostic_columns <- .normalize_diagnostic_columns(diagnostic_columns, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns")
  if(remove_diagnostics){
    diagnostic_columns <- character()
  }

  # compute quantiles dynamically
  quantile_cols <- lapply(probs, function(p) apply(model_samples, 2, stats::quantile, probs = p, na.rm = TRUE))
  names(quantile_cols) <- as.character(probs)

  # the chains needs to be kept merged for conditional summary (due to NAs in the chains)
  runjags_summary <- cbind.data.frame(
    "Mean"   = apply(model_samples, 2, mean,          na.rm = TRUE),
    "SD"     = apply(model_samples, 2, stats::sd,    na.rm = TRUE),
    as.data.frame(quantile_cols, check.names = FALSE)
  )

  # remove all but Mean for inclusions
  quantile_col_names <- as.character(probs)
  inclusion_rows <-
    grepl(" (inclusion", rownames(runjags_summary), fixed = TRUE) |
    grepl(": inclusion(", rownames(runjags_summary), fixed = TRUE)
  runjags_summary[inclusion_rows, c("SD", quantile_col_names)] <- NA

  # don't produce fit diagnostics for conditional samples (different chain lengths etc...) or if remove_diagnostics is TRUE
  if(conditional || length(diagnostic_columns) == 0){
    return(runjags_summary)
  }

  # Split back the chains for diagnostics. Quantities that are undefined for
  # some product-space branches retain NA diagnostics while their estimates are
  # summarized from the branches on which they are defined.
  complete_columns <- colSums(is.na(model_samples)) == 0L
  runjags_diagnostics <- as.data.frame(matrix(
    NA_real_,
    nrow = ncol(model_samples),
    ncol = 4L,
    dimnames = list(
      colnames(model_samples),
      c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
    )
  ))
  if(any(complete_columns)){
    diagnostic_samples <- model_samples[, complete_columns, drop = FALSE]
    model_samples_list <- split(
      as.data.frame(diagnostic_samples),
      rep(seq_len(n_chains), each = n_samples),
      drop = FALSE
    )
    model_samples_list <- coda::as.mcmc.list(lapply(
      model_samples_list,
      coda::as.mcmc
    ))
    mcmc_summary <- summary(model_samples_list, quantiles = NULL)$statistics

    # fix single parameter summaries
    if(is.null(dim(mcmc_summary))){
      mcmc_summary <- t(mcmc_summary)
    }

    complete_diagnostics <- cbind.data.frame(
      "MCMC_error"    = mcmc_summary[, "Time-series SE"],
      "MCMC_SD_error" = mcmc_summary[, "Time-series SE"] / mcmc_summary[, "SD"],
      "ESS"           = coda::effectiveSize(model_samples_list),
      "R_hat"         = if(n_chains > 1){
        coda::gelman.diag(
          model_samples_list,
          multivariate = FALSE,
          autoburnin = FALSE
        )$psrf[, 1]
      }else{
        NA
      }
    )
    runjags_diagnostics[complete_columns, ] <- complete_diagnostics
  }

  # remove incorrect NANs and NAs from the diagnostics
  runjags_diagnostics[is.nan(runjags_diagnostics[,"MCMC_SD_error"]),"MCMC_SD_error"] <- NA
  zero_ess <- !is.na(runjags_diagnostics[, "ESS"]) &
    runjags_diagnostics[, "ESS"] == 0
  runjags_diagnostics[zero_ess, "ESS"]                                               <- 0
  runjags_diagnostics[is.nan(runjags_diagnostics[,"R_hat"]),"R_hat"]                 <- NA

  # first omega parameter is always constant
  runjags_diagnostics[grepl("omega[0,", rownames(runjags_diagnostics), fixed = TRUE), ] <- NA
  runjags_summary <- cbind.data.frame(runjags_summary, runjags_diagnostics[, diagnostic_columns, drop = FALSE])

  return(runjags_summary)
}
.runjags_conditional_warning <- function(parameters, n_samples, warning_limit = 500){
  if(n_samples == 0){
    return(sprintf("Conditional summary for %1$s parameter could not be computed due to no posterior samples.", paste0(parameters, collapse = ", ")))
  }else if(n_samples <= warning_limit){
    return(sprintf("Conditional summary for %1$s is based on %2$i samples.", paste0(parameters, collapse = ", "), n_samples))
  }else{
    return()
  }
}

# Helper function to transform scaled samples in list format (for ensemble/marginal tables)
# Uses the combinatorial unscaling algorithm via the helper in JAGS-formula.R
.transform_scale_samples_list <- function(samples, formula_scale){

  if(is.null(formula_scale) || length(formula_scale) == 0){
    return(samples)
  }

  sample_names <- names(samples)
  nested_transformable <- vapply(samples, function(x){
    is.list(x) && length(x) > 0L &&
      all(vapply(x, function(value) is.numeric(value) || is.matrix(value), logical(1)))
  }, logical(1))
  transformable <- vapply(
    samples,
    function(x) is.numeric(x) || is.matrix(x),
    logical(1)
  ) | nested_transformable
  transformable_names <- sample_names[transformable]

  if(length(transformable_names) == 0){
    return(samples)
  }

  first_name <- transformable_names[1]
  n_samples <- if(nested_transformable[[first_name]]){
    length(samples[[first_name]][[1]])
  }else if(is.matrix(samples[[first_name]])){
    nrow(samples[[first_name]])
  }else{
    length(samples[[first_name]])
  }

  sample_columns <- lapply(transformable_names, function(name){
    if(nested_transformable[[name]]){
      temp_samples <- do.call(cbind, lapply(samples[[name]], as.numeric))
      colnames(temp_samples) <- paste0(name, "[", seq_along(samples[[name]]), "]")
      return(temp_samples)
    }else if(is.matrix(samples[[name]])){
      temp_samples <- samples[[name]]
      if(is.null(colnames(temp_samples))){
        colnames(temp_samples) <- if(ncol(temp_samples) == 1) name else paste0(name, "[", seq_len(ncol(temp_samples)), "]")
      }
      return(temp_samples)
    }else{
      temp_samples <- matrix(samples[[name]], ncol = 1)
      colnames(temp_samples) <- name
      return(temp_samples)
    }
  })
  names(sample_columns) <- transformable_names

  posterior_matrix <- do.call(cbind, sample_columns)
  if(nrow(posterior_matrix) != n_samples){
    stop("All posterior sample elements must contain the same number of samples.", call. = FALSE)
  }

  posterior_matrix <- .apply_unscale_transform(posterior_matrix, formula_scale)

  for(name in transformable_names){
    column_names <- colnames(sample_columns[[name]])

    if(nested_transformable[[name]]){
      for(i in seq_along(samples[[name]])){
        old_sample <- samples[[name]][[i]]
        transformed_sample <- posterior_matrix[, column_names[i]]
        if(is.matrix(old_sample)){
          old_sample[] <- transformed_sample
          samples[[name]][[i]] <- old_sample
        }else{
          old_attrs <- attributes(old_sample)
          samples[[name]][[i]] <- transformed_sample
          for(attr_name in setdiff(names(old_attrs), "names")){
            attr(samples[[name]][[i]], attr_name) <- old_attrs[[attr_name]]
          }
        }
      }
    }else if(is.matrix(samples[[name]])){
      samples[[name]][, seq_along(column_names)] <- posterior_matrix[, column_names, drop = FALSE]
    }else{
      old_attrs <- attributes(samples[[name]])
      samples[[name]] <- posterior_matrix[, column_names]
      for(attr_name in setdiff(names(old_attrs), "names")){
        attr(samples[[name]], attr_name) <- old_attrs[[attr_name]]
      }
    }
  }

  return(samples)
}
