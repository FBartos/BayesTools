.bt_apply_factor_prior_contrasts <- function(data,
                                             predictors_type,
                                             model_terms,
                                             model_terms_type,
                                             prior_list,
                                             context = "Factor predictor",
                                             validate_direct_factor_prior = TRUE){

  factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  if(length(factor_predictors) == 0L){
    return(data)
  }

  for(factor_name in factor_predictors){
    contrast_name <- NULL

    if(factor_name %in% names(prior_list)){
      contrast_name <- .factor_object_contrast_name(prior_list[[factor_name]])
      if(is.null(contrast_name) && isTRUE(validate_direct_factor_prior)){
        stop(paste0("Unsupported prior distribution defined for '", factor_name, "' factor variable. See '?prior_factor' for details."), call. = FALSE)
      }
    }

    if(is.null(contrast_name)){
      factor_terms <- model_terms[
        model_terms_type == "factor" &
          vapply(model_terms, function(term){
            factor_name %in% .bt_random_effect_term_components(term)
          }, logical(1))
      ]
      factor_term_contrasts <- vapply(factor_terms, function(term){
        if(term %in% names(prior_list)){
          contrast <- .factor_object_contrast_name(prior_list[[term]])
          if(is.null(contrast)) NA_character_ else contrast
        }else{
          NA_character_
        }
      }, character(1))
      factor_term_contrasts <- unique(factor_term_contrasts[!is.na(factor_term_contrasts)])

      if(length(factor_term_contrasts) == 1L){
        contrast_name <- factor_term_contrasts
      }else if(length(factor_term_contrasts) > 1L){
        stop(
          context, " '", factor_name,
          "' has conflicting contrast priors across formula terms.",
          call. = FALSE
        )
      }
    }

    if(!is.factor(data[[factor_name]])){
      data[[factor_name]] <- factor(data[[factor_name]])
    }
    if(!is.null(contrast_name) && .prior_ordered_is_contrast_name(contrast_name)){
      data[[factor_name]] <- ordered(data[[factor_name]], levels = levels(data[[factor_name]]))
    }

    if(!is.null(contrast_name)){
      stats::contrasts(data[[factor_name]]) <- contrast_name
    }else if(is.null(attr(data[[factor_name]], "contrasts"))){
      stats::contrasts(data[[factor_name]]) <- "contr.treatment"
    }
  }

  data
}

.bt_random_effect_apply_factor_prior_contrasts <- function(data,
                                                           predictors_type,
                                                           model_terms,
                                                           model_terms_type,
                                                           prior_list){

  .bt_apply_factor_prior_contrasts(
    data = data,
    predictors_type = predictors_type,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list,
    context = "Random-effect factor predictor",
    validate_direct_factor_prior = FALSE
  )
}

.bt_random_effect_term_components <- function(term){

  unlist(strsplit(term, "__xXx__|:", perl = TRUE), use.names = FALSE)
}

.bt_random_effect_factor_term_contrast <- function(model_term, predictors_type, data){

  components <- .bt_random_effect_term_components(model_term)
  factor_components <- components[
    components %in% names(predictors_type) &
      predictors_type[components] == "factor"
  ]
  if(length(factor_components) == 0L){
    return(NULL)
  }

  contrasts <- vapply(factor_components, function(factor_name){
    contrast_name <- attr(data[[factor_name]], "contrasts")
    if(is.null(contrast_name)) "contr.treatment" else contrast_name
  }, character(1))
  contrasts <- unique(contrasts)

  if(length(contrasts) == 1L){
    return(contrasts)
  }
  if(all(contrasts %in% c("contr.treatment", "contr.independent"))){
    return("contr.independent")
  }
  if(all(contrasts %in% c("contr.orthonormal", "contr.meandif"))){
    return("contr.orthonormal")
  }
  if(all(.prior_ordered_is_contrast_name(contrasts))){
    contrasts <- unique(contrasts)
    if(length(contrasts) == 1L){
      return(contrasts)
    }
  }

  stop(
    "Random-effect factor interaction '", model_term,
    "' uses mixed factor contrast families, which is not supported.",
    call. = FALSE
  )
}

.bt_random_term_structure <- function(random_term, prior_random = NULL){

  structure <- .bt_random_effect_structure(random_term)
  if(!is.null(prior_random)){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    requested <- block_prior$covariance$structure
    if(!is.null(requested)){
      requested <- tolower(.bt_random_covariance_normalize(requested))
      if(!identical(requested, structure)){
        stop(
          "Random-effect block '", random_term$block_name,
          "' uses covariance structure '", structure,
          "' in the formula but '", requested,
          "' in 'prior_random'. The formula owns the covariance structure.",
          call. = FALSE
        )
      }
    }
  }

  structure
}

.bt_random_effect_homogeneous_sd <- function(random_term, structure){

  switch(
    structure,
    id = TRUE,
    diag = identical(random_term$hom, TRUE),
    cs = TRUE,
    hcs = FALSE,
    ar1 = TRUE,
    car = TRUE,
    har = FALSE,
    us = FALSE,
    FALSE
  )
}

.bt_JAGS_structured_corr_cholesky <- function(node_prefix, prior_prefix, K,
                                              structure, block_prior,
                                              include_correlation = TRUE,
                                              require_rho = FALSE,
                                              distance_matrix = NULL){

  check_char(node_prefix, "node_prefix", allow_NA = FALSE)
  check_char(prior_prefix, "prior_prefix", allow_NA = FALSE)
  check_int(K, "K", lower = 1, allow_NA = FALSE)
  check_char(structure, "structure", allow_values = c("cs", "hcs", "ar1", "car", "har"), allow_NA = FALSE)
  check_bool(include_correlation, "include_correlation", allow_NA = FALSE)
  check_bool(require_rho, "require_rho", allow_NA = FALSE)
  if(identical(structure, "car")){
    distance_matrix <- .bt_random_effect_validate_car_distance_matrix(distance_matrix, K)
  }

  L_name <- paste0(node_prefix, "_xRE_CORx_L")
  R_name <- paste0(node_prefix, "_xRE_CORx_R")
  rho_name <- paste0(node_prefix, "_rho")

  syntax <- c(paste0("# Structured random-effect correlation: ", structure))
  prior_list <- list()
  monitor <- L_name

  if(K == 1L){
    if(!is.null(block_prior$covariance) && !is.null(block_prior$covariance$rho)){
      stop(
        "Single-column random-effect structure '", structure,
        "' has no correlation parameter; remove the 'rho' prior.",
        call. = FALSE
      )
    }
    syntax <- c(
      syntax,
      paste0(L_name, "[1,1] <- 1")
    )
    if(include_correlation){
      syntax <- c(syntax, paste0(R_name, "[1,1] <- 1"))
      monitor <- c(monitor, R_name)
    }
    return(list(
      syntax = paste0(paste(syntax, collapse = "\n"), "\n"),
      prior_list = prior_list,
      monitor = monitor,
      cholesky_name = L_name,
      correlation_name = if(include_correlation) R_name else NULL,
      rho_name = NULL,
      bridge = NULL
    ))
  }

  rho_info <- .bt_random_effect_structured_rho_prior(
    prior_prefix = prior_prefix,
    node_prefix = node_prefix,
    K = K,
    structure = structure,
    block_prior = block_prior,
    require_rho = require_rho
  )
  syntax <- c(syntax, rho_info$syntax)
  prior_list <- rho_info$prior_list
  monitor <- c(monitor, rho_info$monitor)

  for(row in seq_len(K)){
    for(column in seq_len(K)){
      target <- paste0(R_name, "[", row, ",", column, "]")
      expr <- if(row == column){
        "1"
      }else if(structure %in% c("cs", "hcs")){
        rho_name
      }else if(identical(structure, "car")){
        paste0("pow(", rho_name, ", ", .bt_JAGS_numeric_literal(distance_matrix[row, column]), ")")
      }else{
        paste0("pow(", rho_name, ", ", abs(row - column), ")")
      }
      syntax <- c(syntax, paste0(target, " <- ", expr))
    }
  }

  for(row in seq_len(K)){
    for(column in seq_len(K)){
      target <- paste0(L_name, "[", row, ",", column, "]")
      if(column > row){
        syntax <- c(syntax, paste0(target, " <- 0"))
      }else if(row == column){
        cholesky_sum <- .bt_JAGS_cholesky_crossprod_sum(L_name, row, row, column - 1L)
        syntax <- c(syntax, paste0(target, " <- sqrt(", R_name, "[", row, ",", row, "] - (", cholesky_sum, "))"))
      }else{
        cholesky_sum <- .bt_JAGS_cholesky_crossprod_sum(L_name, row, column, column - 1L)
        syntax <- c(syntax, paste0(
          target, " <- (", R_name, "[", row, ",", column, "] - (", cholesky_sum, ")) / ",
          L_name, "[", column, ",", column, "]"
        ))
      }
    }
  }

  if(include_correlation){
    monitor <- c(monitor, R_name)
  }

  list(
    syntax = paste0(paste(syntax, collapse = "\n"), "\n"),
    prior_list = prior_list,
    monitor = unique(monitor),
    cholesky_name = L_name,
    correlation_name = if(include_correlation) R_name else NULL,
    rho_name = rho_name,
    bridge = list(
      type = "rho",
      structure = structure,
      rho_name = rho_name,
      sample_name = rho_info$sample_name,
      prior_name = rho_info$prior_name,
      sample_fixed = rho_info$sample_fixed,
      rho_scale = rho_info$rho_scale,
      bounds = rho_info$bounds,
      distance_matrix = if(identical(structure, "car")) distance_matrix else NULL,
      cholesky_name = L_name,
      correlation_name = if(include_correlation) R_name else NULL
    )
  )
}

# Compile a scalar structured-correlation prior without materializing a dense
# Cholesky factor. Dense correlation entries are optional derived nodes.
.bt_JAGS_structured_corr_direct <- function(node_prefix, prior_prefix, K,
                                            structure, block_prior,
                                            include_correlation = FALSE,
                                            require_rho = FALSE,
                                            distance_matrix = NULL){

  check_char(node_prefix, "node_prefix", allow_NA = FALSE)
  check_char(prior_prefix, "prior_prefix", allow_NA = FALSE)
  check_int(K, "K", lower = 1, allow_NA = FALSE)
  check_char(
    structure,
    "structure",
    allow_values = c("cs", "hcs", "ar1", "car", "har"),
    allow_NA = FALSE
  )
  check_bool(include_correlation, "include_correlation", allow_NA = FALSE)
  check_bool(require_rho, "require_rho", allow_NA = FALSE)
  if(identical(structure, "car") && isTRUE(include_correlation)){
    distance_matrix <- .bt_random_effect_validate_car_distance_matrix(
      distance_matrix,
      K
    )
  }

  R_name   <- paste0(node_prefix, "_xRE_CORx_R")
  rho_name <- paste0(node_prefix, "_rho")
  syntax   <- paste0("# Direct structured random-effect correlation: ", structure)

  if(K == 1L){
    if(!is.null(block_prior$covariance) && !is.null(block_prior$covariance$rho)){
      stop(
        "Single-column random-effect structure '", structure,
        "' has no correlation parameter; remove the 'rho' prior.",
        call. = FALSE
      )
    }
    if(isTRUE(include_correlation)){
      syntax <- c(syntax, paste0(R_name, "[1,1] <- 1"))
    }
    return(list(
      syntax = paste0(paste(syntax, collapse = "\n"), "\n"),
      prior_list = list(),
      monitor = if(isTRUE(include_correlation)) R_name else character(),
      cholesky_name = NULL,
      correlation_name = if(isTRUE(include_correlation)) R_name else NULL,
      rho_name = NULL,
      bridge = NULL
    ))
  }

  rho_info <- .bt_random_effect_structured_rho_prior(
    prior_prefix = prior_prefix,
    node_prefix = node_prefix,
    K = K,
    structure = structure,
    block_prior = block_prior,
    require_rho = require_rho
  )
  syntax <- c(syntax, rho_info$syntax)

  if(isTRUE(include_correlation)){
    for(row in seq_len(K)){
      for(column in seq_len(K)){
        expression <- if(row == column){
          "1"
        }else if(structure %in% c("cs", "hcs")){
          rho_name
        }else if(identical(structure, "car")){
          paste0(
            "pow(", rho_name, ", ",
            .bt_JAGS_numeric_literal(distance_matrix[row, column]), ")"
          )
        }else{
          paste0("pow(", rho_name, ", ", abs(row - column), ")")
        }
        syntax <- c(
          syntax,
          paste0(R_name, "[", row, ",", column, "] <- ", expression)
        )
      }
    }
  }

  list(
    syntax = paste0(paste(syntax, collapse = "\n"), "\n"),
    prior_list = rho_info$prior_list,
    monitor = unique(c(
      rho_info$monitor,
      if(isTRUE(include_correlation)) R_name else character()
    )),
    cholesky_name = NULL,
    correlation_name = if(isTRUE(include_correlation)) R_name else NULL,
    rho_name = rho_name,
    bridge = list(
      type = "rho",
      structure = structure,
      rho_name = rho_name,
      sample_name = rho_info$sample_name,
      prior_name = rho_info$prior_name,
      sample_fixed = rho_info$sample_fixed,
      rho_scale = rho_info$rho_scale,
      bounds = rho_info$bounds,
      distance_matrix = if(identical(structure, "car")) distance_matrix else NULL,
      cholesky_name = NULL,
      correlation_name = if(isTRUE(include_correlation)) R_name else NULL
    )
  )
}

.bt_random_effect_validate_car_distance_matrix <- function(distance_matrix, K){

  if(is.null(distance_matrix)){
    stop("CAR random-effect structures require a distance matrix.", call. = FALSE)
  }
  if(!is.matrix(distance_matrix) || !is.numeric(distance_matrix) ||
     !all(dim(distance_matrix) == c(K, K))){
    stop("CAR distance matrix must be a numeric K by K matrix.", call. = FALSE)
  }
  if(any(is.na(distance_matrix)) || any(!is.finite(distance_matrix))){
    stop("CAR distance matrix must contain only finite values.", call. = FALSE)
  }
  if(any(distance_matrix < 0)){
    stop("CAR distance matrix cannot contain negative distances.", call. = FALSE)
  }
  if(!isTRUE(all.equal(distance_matrix, t(distance_matrix), tolerance = 1e-12))){
    stop("CAR distance matrix must be symmetric.", call. = FALSE)
  }
  if(any(abs(diag(distance_matrix)) > 1e-12)){
    stop("CAR distance matrix must have a zero diagonal.", call. = FALSE)
  }
  if(K > 1L && any(distance_matrix[row(distance_matrix) != col(distance_matrix)] <= 0)){
    stop("CAR distance matrix must have positive off-diagonal distances.", call. = FALSE)
  }

  distance_matrix
}

.bt_JAGS_numeric_literal <- function(x){

  if(length(x) != 1L || is.na(x) || !is.finite(x)){
    stop("JAGS numeric literal must be finite.", call. = FALSE)
  }

  format(x, scientific = FALSE, trim = TRUE, digits = 17)
}

.bt_JAGS_cholesky_crossprod_sum <- function(L_name, row, column, n_terms){

  if(n_terms < 1L){
    return("0")
  }

  paste0(
    L_name, "[", row, ",", seq_len(n_terms), "] * ",
    L_name, "[", column, ",", seq_len(n_terms), "]",
    collapse = " + "
  )
}

.bt_random_effect_structured_rho_prior <- function(prior_prefix, node_prefix, K,
                                                   structure, block_prior,
                                                   require_rho = FALSE){

  rho_prior <- block_prior$covariance$rho
  if(is.null(rho_prior)){
    if(isTRUE(require_rho)){
      stop(
        "Random-effect structure '", structure,
        "' requires a scalar correlation prior. Supply 'rho = prior(...)'.",
        call. = FALSE
      )
    }
    rho_prior <- prior("normal", list(0, 0.5))
  }

  rho_scale <- block_prior$covariance$rho_scale
  if(is.null(rho_scale)){
    rho_scale <- "fisher_z"
  }

  bounds <- .bt_random_effect_structured_rho_bounds(K = K, structure = structure)
  interior_bounds <- .bt_random_effect_representable_rho_bounds(
    bounds,
    structure
  )
  rho_name <- paste0(node_prefix, "_rho")
  syntax <- character()
  monitor <- character()
  sample_fixed <- NULL

  if(identical(rho_scale, "fisher_z")){
    lower <- if(bounds[["lower"]] <= -1) -Inf else atanh(bounds[["lower"]])
    upper <- if(bounds[["upper"]] >= 1) Inf else atanh(bounds[["upper"]])
    rho_prior <- .bt_random_effect_bound_scalar_prior(
      rho_prior,
      lower = lower,
      upper = upper,
      label = paste0(structure, " Fisher-z correlation prior"),
      warn = FALSE,
      lower_inclusive = .bt_random_effect_rho_lower_inclusive(structure)
    )
    prior_name <- paste0(prior_prefix, "_rho_z")
    sample_name <- paste0(node_prefix, "_rho_z")
    syntax <- c(syntax, paste0(
      rho_name, " <- max(",
      .bt_JAGS_numeric_literal(interior_bounds[["lower"]]), ", min(",
      .bt_JAGS_numeric_literal(interior_bounds[["upper"]]), ", 2 * ilogit(2 * ",
      sample_name, ") - 1))"
    ))
    monitor <- rho_name
  }else if(identical(rho_scale, "logit")){
    prior_name <- paste0(prior_prefix, "_rho_logit")
    sample_name <- paste0(node_prefix, "_rho_logit")
    syntax <- c(syntax, paste0(
      rho_name, " <- ",
      .bt_JAGS_numeric_literal(interior_bounds[["lower"]]), " + ",
      .bt_JAGS_numeric_literal(
        interior_bounds[["upper"]] - interior_bounds[["lower"]]
      ),
      " * ilogit(", sample_name, ")"
    ))
    monitor <- rho_name
  }else{
    rho_prior <- .bt_random_effect_bound_scalar_prior(
      rho_prior,
      lower = bounds[["lower"]],
      upper = bounds[["upper"]],
      label = paste0(structure, " raw correlation prior"),
      warn = TRUE,
      lower_inclusive = .bt_random_effect_rho_lower_inclusive(structure)
    )
    prior_name <- paste0(prior_prefix, "_rho")
    sample_name <- rho_name
  }
  if(is.prior.point(rho_prior)){
    sample_fixed <- rho_prior$parameters[["location"]]
  }

  prior_list <- stats::setNames(list(rho_prior), prior_name)

  list(
    prior_list = prior_list,
    syntax = syntax,
    monitor = monitor,
    prior_list_name = prior_name,
    prior_name = sample_name,
    sample_name = sample_name,
    sample_fixed = sample_fixed,
    rho_scale = rho_scale,
    bounds = bounds
  )
}

.bt_random_effect_structured_rho_bounds <- function(K, structure){

  if(structure %in% c("cs", "hcs")){
    lower <- -1 / (K - 1)
    upper <- 1
  }else if(structure %in% c("ar1", "har")){
    lower <- -1
    upper <- 1
  }else if(identical(structure, "car")){
    lower <- 0
    upper <- 1
  }else{
    stop("Unsupported structured random-effect correlation '", structure, "'.", call. = FALSE)
  }

  c(lower = lower, upper = upper)
}

.bt_random_effect_bound_scalar_prior <- function(x, lower, upper, label,
                                                warn = TRUE,
                                                lower_inclusive = FALSE){

  if(is.prior.none(x)){
    stop(label, " cannot use prior_none().", call. = FALSE)
  }

  if(is.prior.spike_and_slab(x) || is.prior.mixture(x)){
    for(i in seq_along(x)){
      if(is.prior.none(x[[i]])){
        next
      }
      x[[i]] <- .bt_random_effect_bound_scalar_prior(
        x[[i]],
        lower = lower,
        upper = upper,
        label = paste0(label, " component ", i),
        warn = warn,
        lower_inclusive = lower_inclusive
      )
    }
    return(x)
  }

  if(!is.prior.simple(x)){
    stop(label, " must be an ordinary scalar prior.", call. = FALSE)
  }

  if(is.prior.point(x)){
    location <- x$parameters[["location"]]
    lower_violation <- if(isTRUE(lower_inclusive)){
      location < lower
    }else{
      location <= lower
    }
    if(length(location) != 1L || is.na(location) || lower_violation || location >= upper){
      interval <- if(isTRUE(lower_inclusive)){
        paste0("[", lower, ", ", upper, ")")
      }else{
        paste0("(", lower, ", ", upper, ")")
      }
      stop(label, " point mass must lie inside ", interval, ".", call. = FALSE)
    }
    return(x)
  }

  new_lower <- max(x$truncation[["lower"]], lower)
  new_upper <- min(x$truncation[["upper"]], upper)
  if(new_lower >= new_upper){
    stop(label, " has no support inside (", lower, ", ", upper, ").", call. = FALSE)
  }
  if(isTRUE(warn) &&
     (!identical(new_lower, x$truncation[["lower"]]) || !identical(new_upper, x$truncation[["upper"]]))){
    warning(label, " was truncated to the valid correlation range.", immediate. = TRUE, call. = FALSE)
  }
  if(!identical(new_lower, x$truncation[["lower"]]) || !identical(new_upper, x$truncation[["upper"]])){
    x$truncation[["lower"]] <- new_lower
    x$truncation[["upper"]] <- new_upper
  }

  x
}

.bt_random_prior_terms_to_prior_list <- function(block_prior, model_terms,
                                                homogeneous_sd = FALSE){

  sd_prior <- block_prior$sd
  covariance_sd <- if(!is.null(block_prior$covariance)) block_prior$covariance$sd else NULL
  if(!is.null(sd_prior) && !is.null(covariance_sd)){
    stop(
      "Random-effect SD prior was supplied both as 'sd' and 'covariance = random_covariance(sd = ...)'. Supply it in only one place.",
      call. = FALSE
    )
  }
  if(is.null(sd_prior) && !is.null(covariance_sd)){
    sd_prior <- block_prior$covariance$sd
  }
  if(is.null(sd_prior)){
    stop("Random-effect SD prior is missing. Supply 'sd' in prior_random() or random_block().", call. = FALSE)
  }

  if(homogeneous_sd){
    out <- list(sd = sd_prior)
  }else{
    out <- rep(list(sd_prior), length(model_terms))
    names(out) <- model_terms
  }

  if(!is.null(block_prior$terms)){
    term_overrides <- block_prior$terms
    if(is.null(names(term_overrides)) || any(!nzchar(names(term_overrides)))){
      stop("Random-effect term overrides must be named.", call. = FALSE)
    }
    if(anyDuplicated(names(term_overrides))){
      stop("Random-effect term override names must be unique.", call. = FALSE)
    }
    unknown_terms <- setdiff(names(term_overrides), names(out))
    if(length(unknown_terms) > 0L){
      stop("Unknown random-effect term override(s): ", paste(unknown_terms, collapse = ", "), ".", call. = FALSE)
    }
    for(term in names(term_overrides)){
      override <- term_overrides[[term]]
      if(is.prior(override)){
        out[[term]] <- override
      }else if(inherits(override, "random_block") && !is.null(override$sd)){
        out[[term]] <- override$sd
      }else{
        stop("Random-effect term override '", term, "' must be a prior or random_block(sd = ...).", call. = FALSE)
      }
    }
  }

  out
}

.bt_random_effect_force_nonnegative_priors <- function(prior_list){

  for(i in seq_along(prior_list)){
    prior_list[[i]] <- .bt_random_effect_force_nonnegative_prior(
      prior = prior_list[[i]],
      name = names(prior_list)[i]
    )
  }

  prior_list
}

.bt_random_effect_force_nonnegative_prior <- function(prior, name){

  if(is.prior.ordered(prior)){
    prior$total <- .bt_random_effect_force_nonnegative_prior(
      prior = prior$total,
      name = paste0(name, "$total")
    )
    return(prior)
  }

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    for(j in seq_along(prior)){
      prior[[j]] <- .bt_random_effect_force_nonnegative_prior_component(
        prior = prior[[j]],
        label = paste0(j, "-th component in '", name, "'")
      )
    }
    return(prior)
  }

  .bt_random_effect_force_nonnegative_prior_component(
    prior = prior,
    label = paste0("'", name, "'")
  )
}

.bt_random_effect_force_nonnegative_prior_component <- function(prior, label){

  if(is.prior.none(prior)){
    stop(
      "Random-effect SD prior ", label, " cannot use prior_none().",
      call. = FALSE
    )
  }

  if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    if(any(is.na(location)) || any(location < 0)){
      stop(
        "Random-effect SD prior ", label,
        " point mass must be nonnegative.",
        call. = FALSE
      )
    }
    return(prior)
  }

  if(range(prior)[1] < 0){
    warning(
      paste0("The lower bound of the ", label, " prior distribution is below 0. Correcting to 0."),
      immediate. = TRUE,
      call. = FALSE
    )
    prior$truncation$lower <- 0
  }

  prior
}

