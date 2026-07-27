.bt_apply_factor_prior_contrasts <- function(data,
                                             predictors_type,
                                             model_terms,
                                             model_terms_type,
                                             prior_list,
                                             context = "Factor predictor",
                                             validate_direct_factor_prior = TRUE,
                                             contrast_overrides = NULL){

  factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  if(length(factor_predictors) == 0L){
    if(length(contrast_overrides) > 0L){
      stop(
        context,
        " contrast overrides reference predictors that are not factors: ",
        paste(names(contrast_overrides), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(data)
  }
  unknown_overrides <- setdiff(names(contrast_overrides), factor_predictors)
  if(length(unknown_overrides) > 0L){
    stop(
      context,
      " contrast overrides reference predictors outside this design: ",
      paste0("'", unknown_overrides, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for(factor_name in factor_predictors){
    .bt_validate_categorical_values(
      data[[factor_name]],
      factor_name,
      context = context
    )
    direct_contrast <- contrast_overrides[[factor_name]]

    if(factor_name %in% names(prior_list)){
      prior_contrast <- .factor_object_contrast_name(prior_list[[factor_name]])
      if(is.null(prior_contrast) && isTRUE(validate_direct_factor_prior)){
        stop(paste0("Unsupported prior distribution defined for '", factor_name, "' factor variable. See '?prior_factor' for details."), call. = FALSE)
      }
      direct_contrast <- c(direct_contrast, prior_contrast)
    }

    factor_terms <- model_terms[
      model_terms_type == "factor" &
        vapply(model_terms, function(term){
          factor_name %in% .bt_random_effect_term_components(term)
        }, logical(1))
    ]
    if(!is.null(direct_contrast)){
      factor_terms <- factor_terms[
        vapply(factor_terms, function(term){
          components <- .bt_random_effect_term_components(term)
          sum(components %in% factor_predictors) == 1L
        }, logical(1))
      ]
    }
    factor_term_contrasts <- vapply(factor_terms, function(term){
      if(term %in% names(prior_list)){
        contrast <- .factor_object_contrast_name(prior_list[[term]])
        if(is.null(contrast)) NA_character_ else contrast
      }else{
        NA_character_
      }
    }, character(1))
    contrast_names <- unique(c(
      direct_contrast,
      factor_term_contrasts[!is.na(factor_term_contrasts)]
    ))

    if(length(contrast_names) > 1L){
      stop(
        context, " '", factor_name,
        "' has conflicting contrast priors across formula terms.",
        call. = FALSE
      )
    }
    contrast_name <- if(length(contrast_names) == 1L) contrast_names else NULL

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
                                                           prior_list,
                                                           contrast_overrides = NULL){

  .bt_apply_factor_prior_contrasts(
    data = data,
    predictors_type = predictors_type,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list,
    context = "Random-effect factor predictor",
    validate_direct_factor_prior = FALSE,
    contrast_overrides = contrast_overrides
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

  "contr.mixed"
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
    stop(
      "Dense CAR Cholesky compilation is unsupported. Use the CAR Markov ",
      "compiler; use .bt_JAGS_structured_corr_direct() only when derived ",
      "dense correlation nodes are required.",
      call. = FALSE
    )
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
      distance_matrix = NULL,
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
          .bt_JAGS_car_correlation_expression(
            rho_name,
            distance_matrix[row, column]
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
  if(!isTRUE(all(distance_matrix == t(distance_matrix)))){
    stop("CAR distance matrix must be symmetric.", call. = FALSE)
  }
  if(any(diag(distance_matrix) != 0)){
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

  format(
    x,
    scientific = FALSE,
    trim = TRUE,
    digits = 17,
    decimal.mark = "."
  )
}

.bt_JAGS_car_correlation_expression <- function(rho_name, distance){

  paste0(
    "exp(", .bt_JAGS_numeric_literal(distance), " * log(", rho_name, "))"
  )
}

.bt_random_effect_car_rho_support_upper <- function(rho_prior, rho_scale,
                                                    bounds,
                                                    context = "CAR rho prior"){

  support <- .posterior_support_from_prior(rho_prior)
  if(is.null(support) || !is.numeric(support$bounds) ||
     length(support$bounds) != 2L ||
     any(is.na(support$bounds))){
    stop(
      context,
      " does not expose an upper support bound for JAGS innovation validation.",
      call. = FALSE
    )
  }

  sample_upper <- support$bounds[[2L]]
  interior <- .bt_random_effect_representable_rho_bounds(
    bounds,
    structure = "car"
  )
  if(identical(rho_scale, "fisher_z")){
    rho_upper <- tanh(sample_upper)
    rho_upper <- max(
      interior[["lower"]],
      min(interior[["upper"]], rho_upper)
    )
  }else if(identical(rho_scale, "logit")){
    rho_upper <- interior[["lower"]] +
      (interior[["upper"]] - interior[["lower"]]) *
      stats::plogis(sample_upper)
  }else if(identical(rho_scale, "rho")){
    rho_upper <- min(interior[["upper"]], sample_upper)
  }else{
    stop(context, " uses an unsupported rho scale.", call. = FALSE)
  }

  if(!is.numeric(rho_upper) || length(rho_upper) != 1L ||
     is.na(rho_upper) || !is.finite(rho_upper) ||
     rho_upper < bounds[["lower"]] || rho_upper >= bounds[["upper"]]){
    stop(
      context,
      " does not have a usable upper support bound after transformation.",
      call. = FALSE
    )
  }

  rho_upper
}

.bt_random_effect_car_centered_sd_support <- function(
    sd_prior, context = "Centered CAR SD prior"){

  support <- .posterior_support_from_prior(sd_prior)
  if(is.null(support) || !is.numeric(support$bounds) ||
     length(support$bounds) != 2L || any(is.na(support$bounds))){
    stop(
      context,
      " does not expose support bounds for JAGS precision validation.",
      call. = FALSE
    )
  }

  sd_support <- c(
    lower = support$bounds[[1L]],
    upper = support$bounds[[2L]]
  )
  if(sd_support[["lower"]] < 0 || sd_support[["upper"]] <= 0 ||
     sd_support[["lower"]] > sd_support[["upper"]]){
    stop(context, " does not have valid non-negative SD support.",
         call. = FALSE)
  }

  sd_support
}

.bt_random_effect_validate_car_centered_initial_precision <- function(
    sd_prior, block_name){

  centered_sd_support <- .bt_random_effect_car_centered_sd_support(
    sd_prior,
    context = paste0(
      "CAR random-effect block '", block_name, "' centered SD prior"
    )
  )
  initial_precision <- centered_sd_support ^ -2
  invalid_initial_precision <- !is.finite(initial_precision) |
    initial_precision <= 0
  if(any(invalid_initial_precision)){
    endpoint <- which(invalid_initial_precision)[1L]
    stop(
      "CAR random-effect block '", block_name,
      "' has an unrepresentable initial JAGS precision at the ",
      names(centered_sd_support)[endpoint],
      " centered SD support ",
      format(
        centered_sd_support[endpoint],
        digits = 17,
        scientific = TRUE
      ),
      ". The emitted initial precision pow(sd, -2) is non-finite or ",
      "non-positive. Centered CAR SD support endpoints must both be ",
      "representable; supports approaching zero or infinity are not ",
      "representable by the JAGS backend.",
      call. = FALSE
    )
  }

  centered_sd_support
}

.bt_random_effect_validate_car_jags_innovation_support <- function(
    coordinate_sets, rho_prior, rho_scale, bounds, block_name,
    centered = FALSE, sd_prior = NULL){

  if(is.numeric(coordinate_sets)){
    coordinate_sets <- list(coordinate_sets)
  }
  if(!is.list(coordinate_sets)){
    stop("CAR coordinate sets must be supplied as a list.", call. = FALSE)
  }
  check_bool(centered, "centered", allow_NA = FALSE)

  rho_upper <- .bt_random_effect_car_rho_support_upper(
    rho_prior = rho_prior,
    rho_scale = rho_scale,
    bounds = bounds,
    context = paste0("CAR random-effect block '", block_name, "' rho prior")
  )
  centered_sd_support <- if(isTRUE(centered)){
    .bt_random_effect_validate_car_centered_initial_precision(
      sd_prior = sd_prior,
      block_name = block_name
    )
  }else{
    NULL
  }

  for(set in seq_along(coordinate_sets)){
    coordinates <- coordinate_sets[[set]]
    if(!is.numeric(coordinates) || any(is.na(coordinates)) ||
       any(!is.finite(coordinates))){
      stop("CAR time coordinates must contain only finite values.",
           call. = FALSE)
    }
    if(length(coordinates) < 2L){
      next
    }
    gaps <- diff(coordinates)
    if(any(!is.finite(gaps)) || any(gaps <= 0)){
      stop("CAR time coordinates must be strictly increasing with finite gaps.",
           call. = FALSE)
    }

    log_phi <- gaps * log(rho_upper)
    innovation_var <- stats::pexp(-2 * log_phi, rate = 1)
    invalid_innovation <- is.na(log_phi) | is.nan(log_phi) |
      !is.finite(innovation_var) | innovation_var <= 0
    centered_precision <- if(!is.null(centered_sd_support)){
      precision <- vapply(
        centered_sd_support,
        function(sd_endpoint){
          sd_endpoint ^ -2 / innovation_var
        },
        numeric(length(innovation_var))
      )
      dim(precision) <- c(
        length(innovation_var),
        length(centered_sd_support)
      )
      dimnames(precision) <- list(NULL, names(centered_sd_support))
      precision
    }else{
      NULL
    }
    invalid_precision <- if(!is.null(centered_precision)){
      !is.finite(centered_precision) | centered_precision <= 0
    }else{
      matrix(
        FALSE,
        nrow = length(innovation_var),
        ncol = 0L
      )
    }
    if(any(invalid_innovation) || any(invalid_precision)){
      if(any(invalid_innovation)){
        transition <- which(invalid_innovation)[1L]
        reason <- "The stable innovation variance is zero or non-finite."
      }else{
        invalid_index <- which(invalid_precision, arr.ind = TRUE)[1L, ]
        transition <- invalid_index[[1L]]
        endpoint <- invalid_index[[2L]]
        reason <- paste0(
          "At the ", names(centered_sd_support)[endpoint],
          " centered SD support ",
          format(
            centered_sd_support[endpoint],
            digits = 17,
            scientific = TRUE
          ),
          ", the emitted conditional-normal precision ",
          "pow(sd, -2) / innovation_var is non-finite or non-positive."
        )
      }
      stop(
        "CAR random-effect block '", block_name,
        "' has an unrepresentable JAGS innovation between time coordinates ",
        format(coordinates[transition], digits = 17, scientific = TRUE),
        " and ",
        format(coordinates[transition + 1L], digits = 17, scientific = TRUE),
        " (gap = ",
        format(gaps[transition], digits = 17, scientific = TRUE),
        ") at the upper rho support ",
        format(rho_upper, digits = 17, scientific = TRUE),
        ". ", reason, " ",
        "Change the coordinate resolution or restrict the rho support explicitly; ",
        "BayesTools does not add epsilon, jitter coordinates, or rescale CAR time.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
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
      .bt_JAGS_numeric_literal(interior_bounds[["upper"]]), ", tanh(",
      sample_name, ")))"
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
