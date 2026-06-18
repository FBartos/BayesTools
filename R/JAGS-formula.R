#' @title Create JAGS formula syntax and data object
#'
#' @description Creates a JAGS formula syntax, prepares data input, and
#' returns modified prior list for further processing in the \code{JAGS_fit}
#' function.
#'
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored). If the formula contains \code{-1}, it will be
#' automatically converted to include an intercept with a spike(0) prior.
#' The formula can also have a \code{"log(intercept)"} attribute set to \code{TRUE}
#' to generate syntax of the form \code{log(intercept) + sum(beta_i * x_i)}, which
#' is useful for parameters that must be positive (e.g., standard deviation).
#' @param parameter name of the parameter to be created with the formula
#' @param data data.frame containing predictors included in the formula
#' @param prior_list named list of prior distribution of parameters specified within
#' the \code{formula}. When using \code{-1} in the formula, an "intercept" prior
#' can be explicitly specified; otherwise, \code{prior("spike", list(0))} is
#' automatically added. The list can also include two special entries:
#' @param formula_scale named list specifying whether to standardize continuous predictors.
#' If \code{NULL} (default), no standardization is applied. If a named list is provided,
#' continuous predictors with \code{TRUE} values will be standardized (mean-centered and
#' scaled by standard deviation). The intercept is never standardized.
#' @param prior_random optional `prior_random()` object defining random-effect
#' standard-deviation, covariance, monitoring, and prediction policies. Required
#' when \code{formula} contains random effects.
#' @param random_effects_compile optional `random_effects_compile()` object
#' specifying which resolved random-effect blocks should be compiled as sampled
#' random effects and which should be compiled as structural marginalized
#' blocks.
#' \describe{
#'   \item{\code{"__default_continuous"}}{A prior to use for any continuous predictors
#'     (including the intercept) that are not explicitly specified in the prior list.}
#'   \item{\code{"__default_factor"}}{A prior to use for any factor predictors
#'     (including interactions involving factors) that are not explicitly specified
#'     in the prior list.}
#' }
#' These default priors allow for more concise specification when many predictors
#' share the same prior distribution.
#'
#' @details When a formula with \code{-1} (no intercept) is specified, the
#' function automatically removes the \code{-1}, adds an intercept back to the
#' formula, and includes a spike(0) prior for the intercept to ensure equivalent
#' model behavior while maintaining consistent formula parsing.
#'
#' When using default priors (\code{"__default_continuous"} or \code{"__default_factor"}),
#' explicitly specified priors for individual terms take precedence over the defaults.
#' The defaults are only applied to terms that are not already in the prior list.
#'
#' Formula random effects require \code{prior_random}. Random-effect SD priors in
#' \code{prior_list} using \code{"term|group"} names are no longer supported.
#'
#' @examples
#' # simulate data
#' set.seed(1)
#' df <- data.frame(
#'   y      = rnorm(60),
#'   x_cont = rnorm(60),
#'   x_bin  = rbinom(60, 1, .5),
#'   x_fac3 = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
#'   x_fac4 = factor(rep(c("A", "B", "C", "D"), 15), levels = c("A", "B", "C", "D"))
#' )
#'
#' # specify priors with intercept
#' prior_list <- list(
#' "intercept"     = prior("normal", list(0, 1)),
#' "x_cont"        = prior("normal", list(0, .5)),
#' "x_fac3"        = prior_factor("normal",  list(0, 1),  contrast = "treatment"),
#' "x_fac4"        = prior_factor("mnormal", list(0, 1),  contrast = "orthonormal"),
#' "x_fac3:x_fac4" = prior_factor("mnormal", list(0, .5), contrast = "orthonormal")
#' )
#'
#' # create the formula object
#' formula_obj <- JAGS_formula(
#'   formula = ~ x_cont + x_fac3 * x_fac4,
#'   parameter = "mu", data = df, prior_list = prior_list)
#'
#' # using -1 notation (automatically adds spike(0) intercept)
#' prior_list_no_intercept <- list(
#'   "x_fac3" = prior_factor("normal", list(0, 1), contrast = "treatment")
#' )
#' formula_no_intercept <- JAGS_formula(
#'   formula = ~ x_fac3 - 1,
#'   parameter = "mu", data = df, prior_list = prior_list_no_intercept)
#' # Equivalent to specifying intercept = prior("spike", list(0))
#'
#' # using default priors for simpler specification
#' prior_list_defaults <- list(
#'   "__default_continuous" = prior("normal", list(0, 1)),
#'   "__default_factor"     = prior_factor("normal", list(0, 0.5), contrast = "treatment")
#' )
#' formula_defaults <- JAGS_formula(
#'   formula = ~ x_cont + x_fac3,
#'   parameter = "mu", data = df, prior_list = prior_list_defaults)
#' # intercept and x_cont get the default continuous prior
#' # x_fac3 gets the default factor prior
#'
#' @return \code{JAGS_formula} returns a list containing the formula JAGS syntax,
#' JAGS data object, modified prior_list, and (if standardization was applied) a
#' \code{formula_scale} list with standardization information for back-transformation.
#'
#' @seealso [JAGS_fit()]
#' @export
JAGS_formula <- function(formula, parameter, data, prior_list, formula_scale = NULL,
                         prior_random = NULL, random_effects_compile = NULL){

  if(!is.language(formula))
    stop("'formula' must be a formula")
  check_char(parameter, "parameter")
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .bt_check_prior_random(prior_random, allow_NULL = TRUE)
  .bt_check_random_effects_compile(random_effects_compile, allow_NULL = TRUE)
  # formula_scale can be TRUE/FALSE (apply to all continuous predictors) or a named list
  if(!is.null(formula_scale) && !is.logical(formula_scale) && !is.list(formula_scale)){
    stop("'formula_scale' must be NULL, TRUE, FALSE, or a named list.", call. = FALSE)
  }


  # remove the specified response
  formula <- .remove_response(formula)
  # store log(intercept) attribute (for models relying on mu = log(intercept) + sum(beta_i * x_i) trick
  # exp(mu) = intercept * exp(sum(beta_i * x_i)) (e.g., Poisson regression / regression with log link etc...)
  log_intercept  <- isTRUE(attr(formula, "log(intercept)"))
  # store expressions (included later as the literal character input)
  expressions    <- .extract_expressions(formula)
  # store random effects (included later via a formula interface)
  parsed_random_effects <- .bt_parse_random_effects(formula)$terms
  .bt_validate_random_effect_block_names(parsed_random_effects, prior_random)
  random_effects_compile <- .bt_resolve_random_effects_compile(
    random_effects = parsed_random_effects,
    random_effects_compile = random_effects_compile
  )
  random_effects_interface <- .bt_random_effects_interface(parsed_random_effects, prior_random)
  random_predictors_type <- .bt_random_effects_predictor_types(parsed_random_effects, data)
  # remove expressions and random effects from the formula
  formula <- .remove_expressions(formula)
  formula <- .remove_random_effects(formula)

  # handle -1 (no intercept) formulas: always add intercept back with spike(0) prior
  no_intercept_specified <- attr(stats::terms(formula), "intercept") == 0
  if(no_intercept_specified){
    # remove -1 from formula and add intercept back
    formula <- formula_add_intercept(formula)
    # add spike(0) prior for intercept if not already specified
    if(!"intercept" %in% names(prior_list)){
      prior_list[["intercept"]] <- prior("spike", list(0))
    }
  }

  # obtain predictors characteristics factors
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  if(any(!predictors %in% colnames(data)))
    stop(paste0("The ", paste0("'", predictors[!predictors %in% colnames(data)], "'", collapse = ", ")," predictor variable is missing in the data set."))
  predictors_type  <- sapply(predictors, function(predictor){
    if(is.factor(data[[predictor]]) | is.character(data[[predictor]])){
      return("factor")
    }else{
      return("continuous")
    }
  })
  if(length(predictors) == 0L){
    predictors_type <- stats::setNames(character(), character())
  }
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))
  model_terms_type <- sapply(model_terms, function(model_term){
    model_term <- strsplit(model_term, ":")[[1]]
    if(length(model_term) == 1 && model_term == "intercept"){
      return("continuous")
    }else if(any(predictors_type[model_term] == "factor")){
      return("factor")
    }else{
      return("continuous")
    }
  })

  scale_predictors_type <- .bt_merge_predictor_types(predictors_type, random_predictors_type)
  .bt_validate_formula_scale(formula_scale, scale_predictors_type)

  if(length(parsed_random_effects) > 0){
    if(any(.get_grouping_factor(names(prior_list)) != "")){
      stop(
        "Random-effect priors must be supplied through 'prior_random'; 'prior_list' names containing '|' are no longer supported.",
        call. = FALSE
      )
    }
  }

  # handle default priors: __default_factor and __default_continuous
  default_factor_prior     <- prior_list[["__default_factor"]]
  default_continuous_prior <- prior_list[["__default_continuous"]]
  has_defaults <- !is.null(default_factor_prior) || !is.null(default_continuous_prior)

  # remove default priors from prior_list before validation
  prior_list[["__default_factor"]]     <- NULL
  prior_list[["__default_continuous"]] <- NULL

  # fill in missing priors with defaults based on term type
  if(has_defaults){
    for(term in model_terms){
      if(!term %in% names(prior_list)){
        term_type <- model_terms_type[[term]]
        if(term_type == "factor" && !is.null(default_factor_prior)){
          prior_list[[term]] <- default_factor_prior
        }else if(term_type == "continuous" && !is.null(default_continuous_prior)){
          prior_list[[term]] <- default_continuous_prior
        }
      }
    }
  }

  # check that all predictors have a prior distribution
  check_list(prior_list, "prior_list", check_names = model_terms, allow_other = FALSE, all_objects = TRUE)

  # check the prior distribution for each predictor
  data <- .bt_apply_factor_prior_contrasts(
    data = data,
    predictors_type = predictors_type,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list
  )
  scale_info <- list()
  if(any(predictors_type == "continuous")){

    for(continuous in names(predictors_type[predictors_type == "continuous"])){

      # select the corresponding prior for the variable
      this_prior <- prior_list[[continuous]]

      if(is.prior.factor(this_prior)|| is.prior.discrete(this_prior) || is.prior.PET(this_prior) || is.prior.PEESE(this_prior) || is.prior.weightfunction(this_prior)){
        stop(paste0("Unsupported prior distribution defined for '", continuous, "' continuous variable. See '?prior' for details."))
      }
    }
  }

  random_effect_unscaled_data <- data

  # standardize continuous predictors if requested. This includes predictors
  # used only inside random-effect terms, excluding CAR time coordinates.
  if(!is.null(formula_scale) && any(scale_predictors_type == "continuous")){
    for(continuous in names(scale_predictors_type[scale_predictors_type == "continuous"])){
      if(!.bt_should_scale_predictor(formula_scale, continuous)){
        next
      }
      scale_info[[continuous]] <- list(
        mean = mean(data[, continuous], na.rm = TRUE),
        sd   = stats::sd(data[, continuous], na.rm = TRUE)
      )
      if(is.na(scale_info[[continuous]]$sd) || !is.finite(scale_info[[continuous]]$sd) || scale_info[[continuous]]$sd <= 0){
        stop(paste0("Cannot standardize predictor '", continuous, "' because its standard deviation must be positive and finite."), call. = FALSE)
      }
      data[, continuous] <- (data[, continuous] - scale_info[[continuous]]$mean) / scale_info[[continuous]]$sd
    }
  }

  # get the default design matrix
  model_frame  <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.pass),
    error = function(e){
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  if(anyNA(model_frame)){
    stop("Formula predictors contain missing values.", call. = FALSE)
  }
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)
  raw_column_names <- colnames(model_matrix)

  # check whether intercept is unique parameter
  if(sum(grepl("intercept", names(prior_list))) > 1)
    stop("only the intercept parameter can contain 'intercept' in its name.")
  # check whether any reserved term is in usage (note: __default_factor/__default_continuous are reserved but already removed from prior_list)
  .bt_validate_random_effect_reserved_name(
    colnames(data),
    context = "naming variables"
  )


  # replace interaction signs (due to JAGS incompatibility)
  colnames(model_matrix)  <- gsub(":", "__xXx__", colnames(model_matrix))
  column_names            <- colnames(model_matrix)
  names(prior_list)       <- gsub(":", "__xXx__", names(prior_list))
  names(model_terms_type) <- gsub(":", "__xXx__", names(model_terms_type))
  model_terms             <- gsub(":", "__xXx__", model_terms)

  # prepare syntax & data based on the formula
  formula_syntax <- NULL
  random_syntax  <- NULL
  JAGS_data      <- list()
  jags_data_names <- list()
  JAGS_data[[paste0("N_", parameter)]] <- nrow(data)

  # add intercept and prepare the indexing vector
  if(has_intercept){
    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    # use log(intercept) if the formula has the log(intercept) attribute
    if(log_intercept){
      formula_syntax <- c(formula_syntax, paste0("log(", parameter, "_intercept)"))
    }else{
      formula_syntax <- c(formula_syntax, paste0(parameter, "_intercept"))
    }
  }else{
    terms_indexes    <- attr(model_matrix, "assign")
  }

  # add remaining terms (omitting the intercept indexed as NA)
  for(i in unique(terms_indexes[terms_indexes > 0])){

    # extract the corresponding prior distribution for a given coefficient
    this_prior <- prior_list[[model_terms[i]]]

    # check whether the term is an interaction or not and save the corresponding attributes
    attr(this_prior, "interaction") <- grepl("__xXx__", model_terms[i])
    if(.is_prior_interaction(this_prior)){
      attr(this_prior, "interaction_terms") <- strsplit(model_terms[i], "__xXx__")[[1]]
    }


    if(model_terms_type[i] == "continuous"){

      # continuous variables or interactions of continuous variables are simple predictors
      data_name <- paste0(parameter, "_data_", model_terms[i])
      JAGS_data[[data_name]] <- model_matrix[,terms_indexes == i]
      jags_data_names[[model_terms[i]]] <- data_name

      formula_syntax <- c(formula_syntax, paste0(
        if(!is.null(attr(this_prior, "multiply_by"))) paste0(attr(this_prior, "multiply_by"), " * "),
        parameter, "_", model_terms[i],
        " * ",
        parameter, "_data_", model_terms[i], "[i]"
      ))

    }else if(model_terms_type[i] == "factor"){

      # factor variables or interactions with a factor requires factor style prior

      # add levels information attributes to factors
      if(is.prior.independent(this_prior)){
        attr(this_prior, "levels") <- sum(terms_indexes == i)
      }else{
        attr(this_prior, "levels") <- sum(terms_indexes == i) + 1
      }
      if(.is_prior_interaction(this_prior)){
        level_names <- list()
        for(sub_term in strsplit(model_terms[i], "__xXx__")[[1]]){
          if(predictors_type[sub_term] == "factor"){
            level_names[[sub_term]] <- levels(data[[sub_term]])
          }
        }
        attr(this_prior, "level_names") <- level_names
      }else{
        attr(this_prior, "level_names") <- levels(data[[model_terms[i]]])
      }
      attr(this_prior, "term_components") <- strsplit(model_terms[i], "__xXx__", fixed = TRUE)[[1]]
      attr(this_prior, "factor_terms") <- if(is.list(attr(this_prior, "level_names"))) {
        names(attr(this_prior, "level_names"))
      } else {
        model_terms[i]
      }
      attr(this_prior, "factor_contrasts") <- vapply(attr(this_prior, "factor_terms"), function(factor_term) {
        factor_contrast <- attr(data[[factor_term]], "contrasts")
        if(is.null(factor_contrast)){
          "contr.treatment"
        }else if(is.character(factor_contrast)){
          factor_contrast[1]
        }else{
          stop("Unsupported matrix-valued factor contrast metadata.", call. = FALSE)
        }
      }, character(1))
      factor_design_info <- .factor_term_design_from_formula(
        formula         = formula,
        data            = data,
        predictors      = predictors,
        predictors_type = predictors_type,
        term_index      = i,
        term_components = attr(this_prior, "term_components"),
        factor_terms    = attr(this_prior, "factor_terms"),
        has_intercept   = has_intercept
      )
      attr(this_prior, "factor_design")     <- factor_design_info[["design"]]
      attr(this_prior, "factor_cell_names") <- factor_design_info[["cell_names"]]

      data_name <- paste0(parameter, "_data_", model_terms[i])
      JAGS_data[[data_name]] <- model_matrix[,terms_indexes == i, drop = FALSE]
      jags_data_names[[model_terms[i]]] <- data_name
      formula_syntax <- c(formula_syntax, paste0(
        if(!is.null(attr(this_prior, "multiply_by"))) paste0(attr(this_prior, "multiply_by"), " * "),
        "inprod(",
        parameter, "_", model_terms[i],
        ", ",
        parameter, "_data_", model_terms[i], "[i,])"
      ))

    }else{
      stop("Unrecognized model term.")
    }

    # update the corresponding prior distribution back into the prior list
    # (and forward attributes to lower level components in the case of spike and slab and mixture priors)
    if(is.prior.spike_and_slab(this_prior) || is.prior.mixture(this_prior)){
      for(p in seq_along(this_prior)){
        attr(this_prior, "levels")            -> attr(this_prior[[p]], "levels")
        attr(this_prior, "level_names")       -> attr(this_prior[[p]], "level_names")
        attr(this_prior, "interaction")       -> attr(this_prior[[p]], "interaction")
        attr(this_prior, "interaction_terms") -> attr(this_prior[[p]], "interaction_terms")
        attr(this_prior, "term_components")   -> attr(this_prior[[p]], "term_components")
        attr(this_prior, "factor_terms")      -> attr(this_prior[[p]], "factor_terms")
        attr(this_prior, "factor_contrasts")  -> attr(this_prior[[p]], "factor_contrasts")
        attr(this_prior, "factor_design")     -> attr(this_prior[[p]], "factor_design")
        attr(this_prior, "factor_cell_names") -> attr(this_prior[[p]], "factor_cell_names")
      }
      this_prior -> prior_list[[model_terms[i]]]
    }else{
      this_prior -> prior_list[[model_terms[i]]]
    }

  }

  # add expressions input back to the formula
  for(i in seq_along(expressions)){
    formula_syntax <- c(formula_syntax, .clean_from_expression(expressions[[i]]))
  }

  # add random effects back to the formula
  random_scale_terms <- character()
  random_sd_leaves <- list()
  random_correlation_required <- character()
  add_parameters <- character()
  jags_modules <- character()
  required_packages <- character()
  random_sd_binding_context <- .bt_random_sd_binding_context(
    random_effects = parsed_random_effects,
    prior_random = prior_random,
    parameter = parameter
  )
  if(length(random_sd_binding_context$prior_list) > 0L){
    prior_list <- c(prior_list, random_sd_binding_context$prior_list)
  }
  if(length(random_sd_binding_context$syntax) > 0L){
    random_syntax <- c(random_syntax, random_sd_binding_context$syntax)
  }
  compiled_random_effects <- parsed_random_effects
  for(random_i in seq_along(parsed_random_effects)){
    random_effect_data <- data
    random_structure <- .bt_random_effect_structure(parsed_random_effects[[random_i]])
    if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
      random_effect_data <- random_effect_unscaled_data
    }
    compile_mode <- .bt_random_effects_compile_mode(
      random_effects_compile,
      parsed_random_effects[[random_i]]$block_name
    )
    temp_random   <- .JAGS_random_effect_formula(
      parsed_random_effects[[random_i]],
      parameter,
      random_effect_data,
      prior_random = prior_random,
      sd_binding_context = random_sd_binding_context,
      group_data = random_effect_unscaled_data,
      compile_mode = compile_mode
    )
    compiled_random_effects[[random_i]] <- temp_random[["random_effect"]]

    for(data_i in seq_along(temp_random[["data"]])){
      JAGS_data[[names(temp_random[["data"]])[data_i]]] <- temp_random[["data"]][[data_i]]
    }
    random_key <- paste0("__xREx__", attr(compiled_random_effects[[random_i]], "random_block"))
    jags_data_names[random_key] <- list(compiled_random_effects[[random_i]]$jags_data_names)
    random_sd_leaves[[random_key]] <- temp_random[["random_effect"]]$sd_leaves
    if(random_structure %in% c("us", "cs", "hcs", "ar1", "car", "har") &&
       is.numeric(compiled_random_effects[[random_i]]$n_columns) &&
       length(compiled_random_effects[[random_i]]$n_columns) == 1L &&
       !is.na(compiled_random_effects[[random_i]]$n_columns) &&
       compiled_random_effects[[random_i]]$n_columns > 1L){
      random_correlation_required <- c(random_correlation_required, random_key)
    }

    random_syntax  <- c(random_syntax,  temp_random[["random_syntax"]])
    formula_syntax <- c(formula_syntax, temp_random[["formula_term"]])
    prior_list     <- c(prior_list, temp_random[["prior_list"]])
    random_scale_terms <- c(random_scale_terms, temp_random[["random_scale_terms"]])
    add_parameters <- c(add_parameters, temp_random[["add_parameters"]])
    jags_modules <- c(jags_modules, temp_random[["jags_modules"]])
    required_packages <- c(required_packages, temp_random[["required_packages"]])
  }

  # finish the syntax
  formula_syntax <- paste0(
    "for(i in 1:N_", parameter, "){\n",
    "  ", parameter, "[i] = ", paste0(formula_syntax, collapse = " + "), "\n",
    "}\n")
  formula_syntax <- paste0(formula_syntax, paste0(random_syntax, collapse = "\n"), collapse = "\n")

  # add the parameter name as a prefix and attribute to each prior in the list
  names(prior_list) <- paste0(parameter, "_", names(prior_list))
  for(i in seq_along(prior_list)){
    attr(prior_list[[i]], "parameter") <- parameter
  }
  if(.JAGS_prior_list_uses_BayesTools_module(prior_list)){
    jags_modules <- c(jags_modules, "BayesTools")
    required_packages <- c(required_packages, "BayesTools")
  }

  # preserve log(intercept) attribute on output formula
  if(log_intercept){
    attr(formula, "log(intercept)") <- TRUE
  }

  output <- list(
    formula_syntax = formula_syntax,
    data           = JAGS_data,
    prior_list     = prior_list,
    formula        = formula,
    add_parameters = unique(add_parameters),
    jags_modules   = unique(jags_modules),
    required_packages = unique(required_packages)
  )

  # add scale information if standardization was applied
  if(exists("scale_info") && length(scale_info) > 0){
    # add parameter prefix to scale_info names for consistency
    names(scale_info) <- paste0(parameter, "_", names(scale_info))
    # store the parameter prefix as an attribute for later retrieval
    attr(scale_info, "parameter") <- parameter
    # store log_intercept attribute for proper unscaling transformation
    attr(scale_info, "log_intercept") <- log_intercept
    point_terms <- .formula_scale_point_terms(
      prior_list = prior_list,
      parameter = parameter,
      model_terms = model_terms
    )
    if(length(point_terms) > 0L){
      attr(scale_info, "point_terms") <- point_terms
    }
    if(length(random_scale_terms) > 0){
      names(random_scale_terms) <- paste0(parameter, "_", names(random_scale_terms))
      attr(scale_info, "random_effect_terms") <- random_scale_terms
      attr(scale_info, "random_effect_sd_leaves") <- random_sd_leaves
      if(length(random_correlation_required) > 0L){
        attr(scale_info, "random_effect_correlation_required") <- unique(random_correlation_required)
      }
    }
    output$formula_scale <- scale_info
  }

  design_formula <- formula
  attr(design_formula, "log(intercept)") <- NULL

  output$formula_design <- .JAGS_formula_design_object(
    parameter         = parameter,
    formula           = design_formula,
    log_intercept     = log_intercept,
    model_frame       = model_frame,
    source_data       = random_effect_unscaled_data,
    model_matrix      = model_matrix,
    raw_column_names  = raw_column_names,
    column_names      = column_names,
    predictors        = predictors,
    predictors_type   = predictors_type,
    model_terms       = model_terms,
    model_terms_type  = model_terms_type,
    prior_list        = prior_list,
    formula_scale     = output$formula_scale,
    expressions       = expressions,
    random_effects    = compiled_random_effects,
    random_effects_compile = random_effects_compile,
    jags_data_names   = jags_data_names,
    random_allocations = random_sd_binding_context$allocations,
    random_effects_interface = random_effects_interface
  )
  output$random_effects_interface <- random_effects_interface

  return(output)
}

.formula_scale_point_terms <- function(prior_list, parameter, model_terms){

  fixed_names <- paste0(parameter, "_", model_terms)
  fixed_names <- fixed_names[fixed_names %in% names(prior_list)]
  if(length(fixed_names) == 0L){
    return(stats::setNames(numeric(), character()))
  }

  point_terms <- numeric()
  for(fixed_name in fixed_names){
    this_prior <- prior_list[[fixed_name]]
    if(!is.prior.point(this_prior)){
      next
    }
    location <- this_prior$parameters[["location"]]
    if(!is.numeric(location) || length(location) != 1L ||
       is.na(location) || !is.finite(location)){
      next
    }
    point_terms[[fixed_name]] <- location
  }

  point_terms
}

.JAGS_formula_design_object <- function(parameter, formula, log_intercept,
                                        model_frame, source_data, model_matrix,
                                        raw_column_names, column_names,
                                        predictors, predictors_type,
                                        model_terms, model_terms_type,
                                        prior_list, formula_scale,
                                        expressions, random_effects,
                                        random_effects_compile = NULL,
                                        jags_data_names,
                                        random_allocations = list(),
                                        random_effects_interface = NULL){

  formula_terms <- stats::terms(formula)
  attr(formula_terms, ".Environment") <- emptyenv()
  formula_output <- formula
  environment(formula_output) <- emptyenv()
  model_frame_output <- model_frame
  attr(model_frame_output, "terms") <- formula_terms

  qr_info <- qr(model_matrix)
  aliased <- rep(FALSE, ncol(model_matrix))
  if(qr_info$rank < ncol(model_matrix)){
    aliased[qr_info$pivot[(qr_info$rank + 1L):ncol(model_matrix)]] <- TRUE
  }
  names(aliased) <- colnames(model_matrix)

  factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  xlevels <- lapply(factor_predictors, function(predictor){
    if(predictor %in% names(model_frame) && is.factor(model_frame[[predictor]])){
      levels(model_frame[[predictor]])
    }else{
      NULL
    }
  })
  names(xlevels) <- factor_predictors
  xlevels <- xlevels[!vapply(xlevels, is.null, logical(1))]

  if(is.null(random_effects_compile)){
    random_effects_compile <- .bt_random_effects_compile_resolved(
      sampled = vapply(random_effects, function(term) term$block_name, character(1)),
      marginalized = character(),
      mode = .bt_random_effects_compile_modes_from_terms(random_effects)
    )
  }

  out <- list(
    parameter          = parameter,
    formula            = formula_output,
    log_intercept      = isTRUE(log_intercept),
    model_frame        = model_frame_output,
    source_data        = source_data,
    model_matrix       = model_matrix,
    column_names       = column_names,
    raw_column_names   = raw_column_names,
    assign             = attr(model_matrix, "assign"),
    terms              = formula_terms,
    contrasts          = attr(model_matrix, "contrasts"),
    xlevels            = xlevels,
    predictors         = predictors,
    predictor_types    = predictors_type,
    model_terms        = model_terms,
    model_terms_type   = model_terms_type,
    prior_list         = prior_list,
    formula_scale      = formula_scale,
    rank               = qr_info$rank,
    qr_pivot           = qr_info$pivot,
    aliased            = aliased,
    transformed_terms  = expressions,
    random_effects     = random_effects,
    random_effects_compile = random_effects_compile,
    jags_data_names    = jags_data_names,
    random_allocations = random_allocations,
    random_effects_interface = random_effects_interface
  )
  class(out) <- c("BayesTools_formula_design", "list")

  return(out)
}

.bt_random_effects_predictor_types <- function(random_effects, data){

  if(length(random_effects) == 0L){
    return(stats::setNames(character(), character()))
  }

  scale_terms <- random_effects[
    !vapply(random_effects, function(random_term){
      .bt_random_effect_structure(random_term) %in% c("cs", "hcs", "ar1", "car", "har")
    }, logical(1))
  ]
  if(length(scale_terms) == 0L){
    return(stats::setNames(character(), character()))
  }

  predictors <- unique(unlist(lapply(scale_terms, function(random_term){
    formula_terms <- stats::terms(random_term$term_formula)
    as.character(attr(formula_terms, "variables"))[-1L]
  }), use.names = FALSE))
  if(length(predictors) == 0L){
    return(stats::setNames(character(), character()))
  }

  missing_predictors <- predictors[!predictors %in% colnames(data)]
  if(length(missing_predictors) > 0L){
    stop(
      paste0(
        "The ",
        paste0("'", missing_predictors, "'", collapse = ", "),
        " predictor variable is missing in the data set."
      ),
      call. = FALSE
    )
  }

  vapply(predictors, function(predictor){
    if(is.factor(data[, predictor]) || is.character(data[, predictor])){
      "factor"
    }else{
      "continuous"
    }
  }, character(1))
}

.bt_merge_predictor_types <- function(...){

  predictor_types <- list(...)
  predictor_types <- predictor_types[vapply(predictor_types, length, integer(1)) > 0L]
  if(length(predictor_types) == 0L){
    return(stats::setNames(character(), character()))
  }

  out <- predictor_types[[1L]]
  if(length(predictor_types) == 1L){
    return(out)
  }

  for(i in seq(2L, length(predictor_types))){
    current <- predictor_types[[i]]
    conflicts <- intersect(names(out), names(current))
    conflicts <- conflicts[out[conflicts] != current[conflicts]]
    if(length(conflicts) > 0L){
      stop(
        "Predictor type conflicts between fixed and random-effect formulas for: ",
        paste(conflicts, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    out[setdiff(names(current), names(out))] <- current[setdiff(names(current), names(out))]
  }

  out
}

.bt_should_scale_predictor <- function(formula_scale, predictor){

  if(is.logical(formula_scale) && length(formula_scale) == 1L){
    return(isTRUE(formula_scale))
  }
  if(is.list(formula_scale) && !is.null(formula_scale[[predictor]])){
    return(isTRUE(formula_scale[[predictor]]))
  }

  FALSE
}

.bt_validate_formula_scale <- function(formula_scale, predictor_types){

  if(is.null(formula_scale)){
    return(invisible(NULL))
  }

  if(is.logical(formula_scale)){
    check_bool(formula_scale, "formula_scale", allow_NA = FALSE)
    return(invisible(NULL))
  }

  if(!is.list(formula_scale)){
    stop("'formula_scale' must be NULL, a single TRUE or FALSE, or a named list.", call. = FALSE)
  }

  if(is.null(names(formula_scale)) || anyNA(names(formula_scale)) || any(names(formula_scale) == "")){
    stop("'formula_scale' must be a named list.", call. = FALSE)
  }

  if(any(duplicated(names(formula_scale)))){
    stop("'formula_scale' names must be unique.", call. = FALSE)
  }

  unknown_predictors <- setdiff(names(formula_scale), names(predictor_types))
  if(length(unknown_predictors) > 0L){
    stop(
      "The '",
      paste0(unknown_predictors, collapse = "', '"),
      "' entries in 'formula_scale' are not predictor variables in the formula.",
      call. = FALSE
    )
  }

  for(predictor in names(formula_scale)){
    check_bool(
      formula_scale[[predictor]],
      paste0("formula_scale[['", predictor, "']]"),
      allow_NA = FALSE
    )
  }

  scaled_predictors <- names(formula_scale)[vapply(formula_scale, isTRUE, logical(1))]
  noncontinuous_predictors <- scaled_predictors[predictor_types[scaled_predictors] != "continuous"]
  if(length(noncontinuous_predictors) > 0L){
    stop(
      "Only continuous predictors can be standardized; '",
      paste0(noncontinuous_predictors, collapse = "', '"),
      "' in 'formula_scale' is not continuous.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' @title Extract Fitted JAGS Formula Design Metadata
#'
#' @description Returns the fitted formula design metadata stored by
#' [JAGS_fit()]. The design contains the processed formula, fitted model frame,
#' exact model matrix used for JAGS data construction, JAGS-safe coefficient
#' names, contrast and factor-level metadata, rank diagnostics, prior metadata,
#' and formula-scale information.
#'
#' @param fit a fitted object returned by [JAGS_fit()].
#' @param parameter optional formula parameter name. If \code{NULL}, all stored
#' formula designs are returned.
#'
#' @return A named list of formula designs, or one formula design when
#' \code{parameter} is supplied. Returns \code{NULL} when no formula design
#' metadata is stored on \code{fit}.
#'
#' @seealso [JAGS_fit()] [JAGS_formula()]
#' @export
JAGS_formula_design <- function(fit, parameter = NULL){

  check_char(parameter, "parameter", allow_NULL = TRUE)

  formula_design <- attr(fit, "formula_design")
  if(is.null(formula_design)){
    return(NULL)
  }

  if(is.null(parameter)){
    return(formula_design)
  }

  if(!parameter %in% names(formula_design)){
    stop("Formula design for parameter '", parameter, "' was not found.", call. = FALSE)
  }

  return(formula_design[[parameter]])
}

.bt_random_effects_interface <- function(random_effects, prior_random = NULL){

  if(length(random_effects) == 0L){
    return("none")
  }
  if(!is.null(prior_random)){
    return("prior_random")
  }

  stop("Formula random effects require 'prior_random'.", call. = FALSE)
}

.JAGS_random_effect_formula <- function(formula, parameter, data,
                                        prior_random = NULL,
                                        sd_binding_context = NULL,
                                        group_data = data,
                                        compile_mode = c("sampled", "marginalized")){

  if(is.null(prior_random)){
    stop("Formula random effects require 'prior_random'.", call. = FALSE)
  }
  compile_mode <- match.arg(compile_mode)
  sampled_random_effect <- identical(compile_mode, "sampled")

  random_term <- .bt_as_random_effect_term(formula)
  .bt_validate_random_effect_term_supported(random_term)

  # extract the grouping factor information
  grouping_factor <- random_term$group_label
  grouping_values <- .bt_random_group_values(random_term, group_data)
  grouping_factor_levels <- levels(as.factor(grouping_values))
  grouping_mapping       <- as.numeric(factor(grouping_values, levels = grouping_factor_levels))

  formula <- random_term$term_formula
  random_structure <- .bt_random_term_structure(random_term, prior_random)
  structured_formula <- .bt_random_effect_normalize_structured_formula(
    formula = formula,
    data = data,
    structure = random_structure
  )
  formula <- structured_formula$formula
  data <- structured_formula$data
  random_term$term_formula <- formula
  random_term$structured_index <- structured_formula$index

  # obtain predictors characteristics factors (copy from formula)
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  if(any(!predictors %in% colnames(data)))
    stop(paste0("The ", paste0("'", predictors[!predictors %in% colnames(data)], "'", collapse = ", ")," predictor variable is missing in the data set."))
  predictors_type  <- sapply(predictors, function(predictor){
    if(is.factor(data[,predictor]) | is.character(data[,predictor])){
      return("factor")
    }else{
      return("continuous")
    }
  })
  if(length(predictors) == 0L){
    predictors_type <- stats::setNames(character(), character())
  }
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))
  model_terms_type <- sapply(model_terms, function(model_term){
    model_term <- strsplit(model_term, ":")[[1]]
    if(length(model_term) == 1 && model_term == "intercept"){
      return("continuous")
    }else if(any(predictors_type[model_term] == "factor")){
      return("factor")
    }else{
      return("continuous")
    }
  })

  homogeneous_sd <- .bt_random_effect_homogeneous_sd(random_term, random_structure)
  sd_binding <- .bt_random_sd_binding_for_block(
    sd_binding_context,
    random_term$block_name
  )

  block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
  if(is.null(sd_binding) && !is.null(block_prior$sd_source)){
    sd_binding <- .bt_random_sd_binding(
      source = block_prior$sd_source,
      application = "block",
      factors = list(),
      true_allocation = FALSE,
      allocations = list()
    )
  }
  bound_sd <- !is.null(sd_binding)
  .bt_validate_random_block_for_structure(
    block_prior,
    structure = random_structure,
    block_name = random_term$block_name
  )
  row_indexed_external_sd <- .bt_random_sd_binding_has_row_external_source(sd_binding)
  if(bound_sd){
    if(!is.null(block_prior$terms)){
      stop(
        "Random-effect block '", random_term$block_name,
        "' cannot use term-specific SD overrides while it is controlled by an SD binding.",
        call. = FALSE
      )
    }
    prior_list <- list()
    original_prior_names <- character()
  }else{
    prior_list <- .bt_random_prior_terms_to_prior_list(block_prior, model_terms, homogeneous_sd)
    original_prior_names <- names(prior_list)
  }
  monitor_policy <- block_prior$monitor
  if(isTRUE(sampled_random_effect) &&
     isTRUE(row_indexed_external_sd) && isTRUE(monitor_policy$coefficients)){
    stop(
      "Group-level coefficient monitoring is not supported for random-effect block '",
      random_term$block_name,
      "' with a row-indexed external SD source.",
      call. = FALSE
    )
  }
  if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
    monitor_policy$latent <- TRUE
  }

  if(!bound_sd){
    if(homogeneous_sd){
      check_list(prior_list, "prior_list", check_names = "sd", allow_other = TRUE, all_objects = TRUE)
    }else{
      check_list(prior_list, "prior_list", check_names = model_terms, allow_other = TRUE, all_objects = TRUE)
    }
    .bt_random_effect_check_structured_sd_priors(
      prior_list = prior_list,
      model_terms = model_terms,
      random_structure = random_structure,
      homogeneous_sd = homogeneous_sd
    )
    prior_list <- .bt_random_effect_force_nonnegative_priors(prior_list)
  }
  data <- .bt_random_effect_apply_factor_prior_contrasts(
    data = data,
    predictors_type = predictors_type,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list
  )

  # get the design matrix. For no-intercept random formulas, add an intercept
  # while constructing the matrix and drop it afterwards. This preserves
  # BayesTools factor-prior contrasts instead of forcing raw level indicators.
  random_design <- .bt_random_effect_design_matrix(
    formula,
    data,
    preserve_no_intercept_contrasts = !random_structure %in% c("cs", "hcs", "ar1", "car", "har"),
    structure = random_structure,
    block_name = random_term$block_name
  )
  model_frame <- random_design$model_frame
  model_matrix <- random_design$model_matrix
  car_metadata <- random_design$car
  raw_column_names <- colnames(model_matrix)
  random_factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  random_xlevels <- lapply(random_factor_predictors, function(predictor){
    if(predictor %in% names(model_frame) && is.factor(model_frame[[predictor]])){
      levels(model_frame[[predictor]])
    }else{
      NULL
    }
  })
  names(random_xlevels) <- random_factor_predictors
  random_xlevels <- random_xlevels[!vapply(random_xlevels, is.null, logical(1))]

  # check whether intercept is unique parameter
  if(sum(grepl("intercept", names(prior_list))) > 1)
    stop("only the intercept parameter can contain 'intercept' in its name.")
  # check whether any reserved term is in usage
  .bt_validate_random_effect_reserved_name(
    names(prior_list),
    context = "naming variables or prior distributions"
  )

  # replace interaction signs (due to JAGS incompatibility)
  colnames(model_matrix)  <- gsub(":", "__xXx__", colnames(model_matrix))
  column_names            <- colnames(model_matrix)
  names(prior_list)       <- gsub(":", "__xXx__", names(prior_list))
  names(model_terms_type) <- gsub(":", "__xXx__", names(model_terms_type))
  model_terms             <- gsub(":", "__xXx__", model_terms)

  # prepare syntax & data based on the formula
  parameter_suffix <- paste0("_xREx__", random_term$block_name) # priors should not be named with parameter name (done on exit from formula)
  parameter        <- paste0(parameter, "_", parameter_suffix) # variables should be named already here
  random_syntax    <- NULL
  JAGS_data        <- list()
  new_prior_list   <- list()
  random_scale_terms <- character()
  add_parameters <- character()
  jags_modules <- character()
  required_packages <- character()
  correlation_metadata <- NULL

  ### in essence, the following prepares constructors that:
  # 1) samples standardized random effects xRE_Zx[ids, predictors] from a multivariate normal distribution
  # 2) create a vector of by-parameter standard deviation of the random effects xRE_STDx[predictors]
  # 3) multiplies the standardized random effects by parameter-specific standard deviations to create xRE_COEFx[ids, predictors] matrix
  # 4) computes the per observation formula output based on indexing the by-id COEF and selecting the observation variables
  # 5) appends the per-observation output to the higher order formula (done in the formula call itself)

  n_id  <- length(grouping_factor_levels)
  n_par <- ncol(model_matrix)
  if(n_par < 1L){
    stop("Random-effect term '", random_term$block_name, "' does not generate any design columns.", call. = FALSE)
  }
  # step 1:
  if(isTRUE(sampled_random_effect)){
    random_syntax <- c(random_syntax, .add_JAGS_matrix(name = paste0(parameter, "_xRE_PRECx"), diag(1, n_par)))
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",n_id,"){\n",
      "   ",paste0(parameter, "_xRE_Zx"),"[i,1:", n_par ,"] ~ dmnorm(rep(0, ", n_par,"), ", paste0(parameter, "_xRE_PRECx"), ")\n",
      " }\n"
    ))
  }

  # step 2
  sd_spec <- .bt_random_effect_sd_spec(
    parameter = parameter,
    parameter_suffix = parameter_suffix,
    prior_list = prior_list,
    random_term = random_term,
    grouping_factor = grouping_factor,
    sd_binding = sd_binding,
    model_matrix = model_matrix,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    predictors_type = predictors_type,
    data = data,
    random_structure = random_structure,
    has_intercept = has_intercept,
    homogeneous_sd = homogeneous_sd
  )
  terms_indexes <- sd_spec$terms_indexes
  sd_parameter_names <- sd_spec$sd_parameter_names
  sd_leaves <- sd_spec$sd_leaves
  random_scale_terms <- c(random_scale_terms, sd_spec$random_scale_terms)
  add_parameters <- c(add_parameters, sd_spec$add_parameters)
  random_syntax <- c(random_syntax, sd_spec$syntax)
  new_prior_list <- c(new_prior_list, sd_spec$prior_list)
  sd_binding <- sd_spec$sd_binding
  row_indexed_external_sd <- .bt_random_sd_binding_has_row_external_source(sd_binding)

  # step 3
  if(random_structure == "us" && n_par == 1L){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    if(!is.null(block_prior$covariance$cor)){
      stop(
        "Single-column random-effect structure 'us' has no correlation parameter; remove the 'cor' prior.",
        call. = FALSE
      )
    }
  }
  if(random_structure %in% c("diag", "id") ||
     (random_structure == "us" && n_par == 1L)){
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(i in 1:",n_par,"){\n",
        "   ",paste0(parameter, "_xRE_UNIT_COEFx"),"[1:",n_id,",i] = ",paste0(parameter, "_xRE_Zx"),"[1:",n_id,",i]\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(i in 1:",n_par,"){\n",
        "   ",paste0(parameter, "_xRE_COEFx"),"[1:",n_id,",i] = ",paste0(parameter, "_xRE_Zx"),"[1:",n_id,",i] * ",paste0(parameter, "_xRE_STDx"),"[i]\n",
        " }\n"
      ))
    }
  }else if(random_structure == "us"){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    if(n_par > 1L && is.null(block_prior$covariance$cor)){
      stop(
        "Random-effect block '", random_term$block_name,
        "' with structure 'us' requires an LKJ correlation prior. ",
        "Supply 'cor = prior_lkj(eta = ...)'.",
        call. = FALSE
      )
    }
    lkj_prior <- .bt_random_block_lkj_prior(block_prior)
    lkj_module <- JAGS_lkj_corr_cholesky(
      name = paste0(parameter, "_xRE_CORx"),
      K = n_par,
      eta = lkj_prior$eta,
      include_correlation = isTRUE(monitor_policy$correlation) && isTRUE(lkj_prior$include_correlation),
      include_primitives = isTRUE(monitor_policy$lkj_primitives) || isTRUE(lkj_prior$include_primitives)
    )
    random_syntax <- c(random_syntax, lkj_module$syntax)
    add_parameters <- c(add_parameters, lkj_module$monitor, lkj_module$primitive_names)
    jags_modules <- c(jags_modules, lkj_module$jags_module)
    required_packages <- c(required_packages, lkj_module$required_packages)
    correlation_metadata <- list(
      type = "lkj",
      eta = lkj_prior$eta,
      primitive_names = lkj_module$primitive_names,
      primitive_bounds = lkj_module$primitive_bounds,
      cholesky_name = lkj_module$cholesky_name,
      correlation_name = lkj_module$correlation_name
    )
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_UNIT_COEFx"),"[g,i] = inprod(", lkj_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_COEFx"),"[g,i] = ",paste0(parameter, "_xRE_STDx"),"[i] * inprod(", lkj_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }
  }else if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    corr_module <- .bt_JAGS_structured_corr_cholesky(
      node_prefix = parameter,
      prior_prefix = parameter_suffix,
      K = n_par,
      structure = random_structure,
      block_prior = block_prior,
      include_correlation = isTRUE(monitor_policy$correlation),
      require_rho = TRUE,
      distance_matrix = if(identical(random_structure, "car")) car_metadata$distance_matrix else NULL
    )
    random_syntax <- c(random_syntax, corr_module$syntax)
    for(corr_prior_name in names(corr_module$prior_list)){
      corr_module$prior_list[[corr_prior_name]] <- .bt_random_effect_set_prior_metadata(
        corr_module$prior_list[[corr_prior_name]],
        block = random_term$block_name,
        grouping = grouping_factor,
        type = "correlation",
        structure = random_structure,
        name = .bt_random_effect_public_name(random_term)
      )
    }
    new_prior_list <- c(new_prior_list, corr_module$prior_list)
    add_parameters <- c(add_parameters, corr_module$monitor)
    correlation_metadata <- corr_module$bridge
    if(identical(random_structure, "car")){
      correlation_metadata$time_variable <- car_metadata$time_variable
      correlation_metadata$time_values <- car_metadata$time_values
    }
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_UNIT_COEFx"),"[g,i] = inprod(", corr_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_COEFx"),"[g,i] = ",paste0(parameter, "_xRE_STDx"),"[i] * inprod(", corr_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }
  }

  # step 4
  if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
    if(isTRUE(sd_binding$true_allocation)){
      allocation_target <- .bt_random_effect_allocation_target_metadata(
        sd_binding$allocations[[1L]],
        context = paste0(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "'"
        )
      )
      if(length(sd_binding$factors_by_column) > 0L &&
         !identical(allocation_target, "sd_component")){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' with row-indexed column factor chains require 'allocation$target' to be 'sd_component'.",
          call. = FALSE
        )
      }
    }
    if(isTRUE(sd_binding$true_allocation) &&
       length(sd_binding$allocations) > 0L &&
       identical(sd_binding$allocations[[1L]]$target, "sd_component")){
      .bt_check_random_sd_component_binding(
        binding = sd_binding,
        n_columns = n_par,
        context = paste0(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "'"
        )
      )
    }
    source_expression <- .bt_random_sd_binding_shared_source_expression(
      sd_binding,
      row_index = "i"
    )
    if(identical(sd_binding$application, "column")){
      if(length(sd_binding$factors_by_column) == 0L){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' are missing canonical 'binding$factors_by_column'.",
          call. = FALSE
        )
      }
      if(length(sd_binding$factors_by_column) != n_par){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' do not match the number of random-effect columns.",
          call. = FALSE
        )
      }
      column_terms <- vapply(seq_len(n_par), function(column){
        factor_expression <- .bt_random_sd_binding_factors_expression(
          sd_binding$factors_by_column[[column]]
        )
        paste(
          c(
            if(!identical(factor_expression, "1")) factor_expression,
            paste0(parameter, "_xRE_UNIT_COEFx[", parameter, "_xRE_MAPx[i],", column, "]"),
            paste0(parameter, "_xRE_DATAx[i,", column, "]")
          ),
          collapse = " * "
        )
      }, character(1))
      row_contribution <- paste0(source_expression, " * (", paste(column_terms, collapse = " + "), ")")
    }else{
      factor_expression <- .bt_random_sd_binding_factors_expression(sd_binding$factors)
      unit_contribution <- paste0(
        "inprod(",
        paste0(parameter, "_xRE_UNIT_COEFx[", paste0(parameter, "_xRE_MAPx[i]"), ", 1:", n_par, "]"),
        ", ",
        paste0(parameter, "_xRE_DATAx[i,1:", n_par, "]"),
        ")"
      )
      row_contribution <- paste0(source_expression, " * ", factor_expression, " * ", unit_contribution)
    }
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",nrow(model_matrix),"){\n",
      "   ",parameter,"[i] = ", row_contribution, "\n",
      " }\n"
    ))
  }else if(isTRUE(sampled_random_effect)){
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",nrow(model_matrix),"){\n",
      "   ",parameter,"[i] = inprod(", paste0(parameter, "_xRE_COEFx[", paste0(parameter, "_xRE_MAPx[i]"),", 1:",n_par,"]"), ", ", paste0(parameter, "_xRE_DATAx[i,1:", n_par,"]"),")\n",
      " }\n"
    ))
  }

  # create the JAGS data list
  if(isTRUE(sampled_random_effect)){
    JAGS_data[[paste0(parameter, "_xRE_DATAx")]] <- model_matrix
    JAGS_data[[paste0(parameter, "_xRE_MAPx")]]  <- grouping_mapping
  }

  if(isTRUE(sampled_random_effect) && isTRUE(monitor_policy$latent)){
    add_parameters <- c(add_parameters, paste0(parameter, "_xRE_Zx"))
  }
  if(isTRUE(sampled_random_effect) &&
     isTRUE(monitor_policy$coefficients) && !isTRUE(row_indexed_external_sd)){
    add_parameters <- c(add_parameters, paste0(parameter, "_xRE_COEFx"))
  }

  random_term$model_matrix    <- model_matrix
  random_term$raw_column_names <- raw_column_names
  random_term$column_names     <- column_names
  random_term$contrasts        <- attr(model_matrix, "contrasts")
  random_term$xlevels          <- random_xlevels
  random_term$assign           <- attr(model_matrix, "assign")
  random_term$model_terms      <- model_terms
  random_term$model_terms_type <- model_terms_type
  random_term$group_levels     <- grouping_factor_levels
  random_term$group_map        <- grouping_mapping
  random_term$n_groups         <- n_id
  random_term$n_columns        <- n_par
  random_term$prior_terms      <- original_prior_names
  random_term$parameter_stem   <- parameter
  random_term$sd_parameter_names <- sd_parameter_names
  random_term$sd_leaves          <- sd_leaves
  random_jags_data_names <- names(JAGS_data)
  if(is.null(random_jags_data_names)){
    random_jags_data_names <- character()
  }
  random_term["jags_data_names"] <- list(random_jags_data_names)
  random_term$structure        <- random_structure
  random_term$homogeneous_sd   <- homogeneous_sd
  random_term$interface        <- "prior_random"
  random_term$sd_binding       <- sd_binding
  random_term$correlation      <- correlation_metadata
  random_term$car              <- car_metadata
  random_term$new_levels       <- block_prior$new_levels
  random_term$compile_mode     <- compile_mode
  attr(random_term, "random_block") <- random_term$block_name
  attr(random_term, "compile_mode") <- compile_mode

  return(list(
    random_syntax  = random_syntax,
    formula_term   = if(isTRUE(sampled_random_effect)) paste0(parameter,"[i]") else character(),
    data           = JAGS_data,
    prior_list     = new_prior_list,
    random_scale_terms = random_scale_terms,
    add_parameters = unique(add_parameters),
    jags_modules   = unique(jags_modules),
    required_packages = unique(required_packages),
    random_effect  = random_term,
    formula        = formula
  ))
}

.bt_random_effect_design_matrix <- function(formula, data,
                                            preserve_no_intercept_contrasts = TRUE,
                                            structure = NULL,
                                            car_time_values = NULL,
                                            block_name = NULL){

  if(identical(structure, "car")){
    return(.bt_random_effect_car_design_matrix(
      formula = formula,
      data = data,
      car_time_values = car_time_values
    ))
  }

  formula_terms <- stats::terms(formula)
  has_intercept <- attr(formula_terms, "intercept") == 1L
  matrix_formula <- formula
  if(!has_intercept && isTRUE(preserve_no_intercept_contrasts)){
    matrix_formula <- formula_add_intercept(formula)
  }

  model_frame <- stats::model.frame(matrix_formula, data = data)
  model_matrix <- stats::model.matrix(model_frame, formula = matrix_formula, data = data)
  if(nrow(model_matrix) != nrow(data)){
    label <- "Random-effect term"
    if(!is.null(block_name) && nzchar(block_name)){
      label <- paste0("Random-effect block '", block_name, "'")
    }
    stop(
      label,
      " contains missing predictor values; random-effect design matrices must have one row per data row.",
      call. = FALSE
    )
  }

  if(!has_intercept && isTRUE(preserve_no_intercept_contrasts)){
    intercept_column <- which(colnames(model_matrix) == "(Intercept)")
    if(length(intercept_column) == 1L){
      assign <- attr(model_matrix, "assign")
      model_matrix <- model_matrix[, -intercept_column, drop = FALSE]
      attr(model_matrix, "assign") <- assign[-intercept_column]
    }
  }

  list(
    model_frame = model_frame,
    model_matrix = model_matrix
  )
}

.bt_random_effect_normalize_structured_formula <- function(formula, data,
                                                           structure){

  out <- list(formula = formula, data = data, index = NULL)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(out)
  }

  if(identical(structure, "car")){
    out$formula <- .bt_random_effect_formula_without_intercept(formula)
    return(out)
  }

  index_variables <- .bt_random_effect_structured_index_variables(formula, structure)
  if(structure %in% c("ar1", "har") && length(index_variables) > 1L){
    stop(
      "The '", structure,
      "' random-effect covariance structure requires a single ordered index variable. ",
      "Create an explicit ordered index column before calling '", structure, "()'.",
      call. = FALSE
    )
  }

  missing_variables <- index_variables[!index_variables %in% names(data)]
  if(length(missing_variables) > 0L){
    stop(
      paste0(
        "The ",
        paste0("'", missing_variables, "'", collapse = ", "),
        " structured random-effect index variable is missing in the data set."
      ),
      call. = FALSE
    )
  }

  index_name <- .bt_random_effect_structured_index_name(index_variables)
  data[[index_name]] <- .bt_random_effect_structured_index_values(data, index_variables)
  out$data <- data
  out$formula <- stats::as.formula(
    call("~", call("-", as.name(index_name), 1)),
    env = environment(formula)
  )
  out$index <- list(
    variables = index_variables,
    name = index_name,
    label = paste(index_variables, collapse = ":"),
    structure = structure
  )

  out
}

.bt_random_effect_formula_without_intercept <- function(formula){

  rhs_index <- if(length(formula) == 3L) 3L else 2L
  stats::as.formula(
    call("~", call("-", formula[[rhs_index]], 1)),
    env = environment(formula)
  )
}

.bt_random_effect_structured_index_variables <- function(formula, structure){

  rhs_index <- if(length(formula) == 3L) 3L else 2L
  rhs <- formula[[rhs_index]]
  usage <- if(structure %in% c("cs", "hcs")){
    paste0("Use '", structure, "(index | group)' or '", structure, "(index1 + index2 | group)'.")
  }else{
    paste0("Use '", structure, "(index | group)'.")
  }
  if(.bt_random_effect_structured_index_has_intercept_control(rhs)){
    stop(
      "The '", structure,
      "' random-effect covariance structure uses index variables and does not support explicit ",
      "'1', '0', or '-1' terms. ", usage,
      call. = FALSE
    )
  }

  variables <- .bt_random_effect_structured_index_plus_terms(rhs, structure)
  if(length(variables) == 0L || any(!nzchar(variables))){
    stop(
      "The '", structure,
      "' random-effect covariance structure requires at least one index variable.",
      call. = FALSE
    )
  }
  if(anyDuplicated(variables)){
    stop(
      "Structured random-effect index variables must be unique.",
      call. = FALSE
    )
  }

  variables
}

.bt_random_effect_structured_index_plus_terms <- function(expr, structure){

  if(is.symbol(expr)){
    return(as.character(expr))
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    return(c(
      .bt_random_effect_structured_index_plus_terms(expr[[2L]], structure),
      .bt_random_effect_structured_index_plus_terms(expr[[3L]], structure)
    ))
  }

  stop(
    "The '", structure,
    "' random-effect covariance structure supports index variables",
    if(structure %in% c("cs", "hcs")) " separated by '+'." else ".",
    " ",
    if(structure %in% c("cs", "hcs")){
      paste0("Use '", structure, "(index | group)' or '", structure, "(index1 + index2 | group)'.")
    }else{
      paste0("Use '", structure, "(index | group)'.")
    },
    call. = FALSE
  )
}

.bt_random_effect_structured_index_has_intercept_control <- function(expr){

  if(is.numeric(expr) && length(expr) == 1L && expr %in% c(0, 1)){
    return(TRUE)
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L){
    return(.bt_random_effect_structured_index_has_intercept_control(expr[[2L]]))
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L){
    return(
      .bt_random_effect_structured_index_has_intercept_control(expr[[2L]]) ||
        .bt_random_effect_structured_index_has_intercept_control(expr[[3L]])
    )
  }
  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    return(
      .bt_random_effect_structured_index_has_intercept_control(expr[[2L]]) ||
        .bt_random_effect_structured_index_has_intercept_control(expr[[3L]])
    )
  }

  FALSE
}

.bt_random_effect_structured_index_name <- function(variables){

  if(length(variables) == 1L){
    return(variables)
  }

  .bt_random_effect_sanitize_name(paste(variables, collapse = "_"))
}

.bt_random_effect_structured_index_values <- function(data, variables){

  values <- lapply(variables, function(variable){
    .bt_random_effect_structured_index_component(data[[variable]])
  })

  if(length(values) == 1L){
    return(values[[1L]])
  }

  do.call(interaction, c(values, list(drop = TRUE, lex.order = TRUE)))
}

.bt_random_effect_structured_index_component <- function(x){

  if(is.factor(x)){
    return(x)
  }
  if(is.character(x)){
    return(factor(x, levels = sort(unique(x))))
  }
  if(is.numeric(x) || is.integer(x) || is.logical(x)){
    return(factor(x, levels = sort(unique(x))))
  }

  factor(x, levels = sort(unique(x)))
}

.bt_random_effect_car_design_matrix <- function(formula, data,
                                                car_time_values = NULL){

  formula_terms <- stats::terms(formula)
  has_intercept <- attr(formula_terms, "intercept") == 1L
  if(has_intercept){
    stop(
      "CAR random-effect terms require exactly one untransformed time variable. Use 'car(time | group)'.",
      call. = FALSE
    )
  }

  term_labels <- attr(formula_terms, "term.labels")
  predictors <- as.character(attr(formula_terms, "variables"))[-1L]
  if(length(term_labels) != 1L || length(predictors) != 1L ||
     !identical(term_labels, predictors)){
    stop(
      "CAR random-effect terms require exactly one untransformed time variable. Use 'car(time | group)'.",
      call. = FALSE
    )
  }

  time_name <- predictors
  if(!time_name %in% names(data)){
    stop("The '", time_name, "' CAR time variable is missing in the data set.", call. = FALSE)
  }

  observed_time <- .bt_random_effect_car_observed_time(data[[time_name]], time_name)
  if(is.null(car_time_values)){
    car_time_values <- sort(unique(observed_time))
  }else{
    car_time_values <- .bt_random_effect_car_reference_time_values(car_time_values)
    new_time <- observed_time[!observed_time %in% car_time_values]
    if(length(new_time) > 0L){
      stop(
        "New CAR time coordinate(s) for random-effect prediction are not supported: ",
        paste(sort(unique(new_time)), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  time_index <- match(observed_time, car_time_values)
  model_matrix <- matrix(0, nrow = length(observed_time), ncol = length(car_time_values))
  if(length(observed_time) > 0L){
    model_matrix[cbind(seq_along(observed_time), time_index)] <- 1
  }
  colnames(model_matrix) <- paste0(time_name, .bt_random_effect_car_time_suffix(car_time_values))
  attr(model_matrix, "assign") <- rep(1L, ncol(model_matrix))
  attr(model_matrix, "contrasts") <- NULL

  model_frame <- data[time_name]
  attr(model_frame, "terms") <- formula_terms

  list(
    model_frame = model_frame,
    model_matrix = model_matrix,
    car = list(
      time_variable = time_name,
      time_values = car_time_values,
      distance_matrix = abs(outer(car_time_values, car_time_values, "-"))
    )
  )
}

.bt_random_effect_car_observed_time <- function(x, time_name){

  if(is.numeric(x) || is.integer(x)){
    out <- as.numeric(x)
  }else if(is.ordered(x)){
    level_values <- suppressWarnings(as.numeric(as.character(levels(x))))
    if(any(is.na(level_values))){
      stop(
        "Ordered factor CAR time variable '", time_name,
        "' must have numeric level labels.",
        call. = FALSE
      )
    }
    out <- level_values[match(as.character(x), levels(x))]
  }else{
    stop(
      "CAR time variable '", time_name,
      "' must be numeric or an ordered factor with numeric level labels.",
      call. = FALSE
    )
  }

  if(any(is.na(out)) || any(!is.finite(out))){
    stop("CAR time variable '", time_name, "' must contain only finite values.", call. = FALSE)
  }

  out
}

.bt_random_effect_car_reference_time_values <- function(x){

  if(!is.numeric(x) && !is.integer(x)){
    stop("CAR reference time values must be numeric.", call. = FALSE)
  }
  x <- as.numeric(x)
  if(length(x) == 0L || any(is.na(x)) || any(!is.finite(x))){
    stop("CAR reference time values must be finite and non-empty.", call. = FALSE)
  }
  if(anyDuplicated(x)){
    stop("CAR reference time values must be unique.", call. = FALSE)
  }
  sort(x)
}

.bt_random_effect_car_time_suffix <- function(x){

  labels <- vapply(x, function(value){
    format(value, scientific = FALSE, trim = TRUE, digits = 15)
  }, character(1))
  labels <- gsub("-", "m", labels, fixed = TRUE)
  labels <- gsub(".", "p", labels, fixed = TRUE)
  labels <- gsub("[^A-Za-z0-9_]", "_", labels)
  labels <- gsub("_+", "_", labels)
  labels <- gsub("^_|_$", "", labels)
  paste0("_", labels)
}

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

    if(!is.factor(data[[factor_name]])){
      data[[factor_name]] <- factor(data[[factor_name]])
    }

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

  format(x, scientific = FALSE, trim = TRUE, digits = 15)
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
    syntax <- c(syntax, paste0(rho_name, " <- 2 * ilogit(2 * ", node_prefix, "_rho_z) - 1"))
    monitor <- rho_name
  }else if(identical(rho_scale, "logit")){
    prior_name <- paste0(prior_prefix, "_rho_logit")
    sample_name <- paste0(node_prefix, "_rho_logit")
    syntax <- c(syntax, paste0(
      rho_name, " <- ", .bt_JAGS_numeric_literal(bounds[["lower"]]), " + ",
      .bt_JAGS_numeric_literal(bounds[["upper"]] - bounds[["lower"]]),
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

# formula helper functions
.remove_response        <- function(formula){
  # removes response from the expression
  # (prevents crash on formula evaluations)
  if(attr(stats::terms(formula), "response")  == 1){
    formula[2] <- NULL
  }
  return(formula)
}
.has_expression         <- function(formula){
  # check if there is any expression in the formula
  return(.bt_contains_expression_call(.bt_formula_rhs(formula)))
}
.extract_expressions    <- function(formula){
  # extract all expressions from the formula

  return(.bt_extract_expression_bodies(.bt_formula_rhs(formula)))
}
.clean_from_expression  <- function(x){
  # expression to character

  return(sub("expression\\((.*)\\)", "\\1", x))
}
.remove_expressions     <- function(formula){
  # remove all expressions from the formula

  rhs_index <- .bt_formula_rhs_index(formula)
  rhs <- .bt_remove_expression_terms(formula[[rhs_index]])
  if(is.null(rhs)){
    rhs <- 1
  }
  formula[[rhs_index]] <- rhs

  return(formula)
}
.bt_formula_rhs_index <- function(formula){
  if(length(formula) == 3L) 3L else 2L
}
.bt_formula_rhs <- function(formula){
  formula[[.bt_formula_rhs_index(formula)]]
}
.bt_is_expression_call <- function(x){
  is.call(x) && identical(as.character(x[[1L]]), "expression")
}
.bt_contains_expression_call <- function(x){
  if(.bt_is_expression_call(x)){
    return(TRUE)
  }
  if(is.call(x) || is.pairlist(x)){
    return(any(vapply(as.list(x)[-1], .bt_contains_expression_call, logical(1))))
  }

  FALSE
}
.bt_extract_expression_bodies <- function(x){
  if(.bt_is_expression_call(x)){
    return(lapply(as.list(x)[-1], function(expression_body){
      paste0(deparse(expression_body), collapse = " ")
    }))
  }
  if(is.call(x) || is.pairlist(x)){
    return(unlist(lapply(as.list(x)[-1], .bt_extract_expression_bodies), recursive = FALSE))
  }

  list()
}
.bt_remove_expression_terms <- function(x){
  if(.bt_is_expression_call(x)){
    return(NULL)
  }

  if(is.call(x)){
    call_name <- as.character(x[[1L]])
    if(call_name == "+" && length(x) == 3L){
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- .bt_remove_expression_terms(x[[3L]])
      if(is.null(lhs)){
        return(rhs)
      }
      if(is.null(rhs)){
        return(lhs)
      }
      return(call("+", lhs, rhs))
    }
    if(call_name == "-" && length(x) == 3L){
      lhs <- .bt_remove_expression_terms(x[[2L]])
      rhs <- .bt_remove_expression_terms(x[[3L]])
      if(is.null(lhs) && is.null(rhs)){
        return(NULL)
      }
      if(is.null(rhs)){
        return(lhs)
      }
      if(is.null(lhs)){
        stop("expression() terms must be additive formula terms.", call. = FALSE)
      }
      return(call("-", lhs, rhs))
    }
    if(.bt_contains_expression_call(x)){
      stop("expression() terms must be additive formula terms.", call. = FALSE)
    }
  }

  x
}
.has_random_effects     <- function(formula){
  return(length(.bt_parse_random_effects(formula)$terms) > 0L)
}
.remove_random_effects  <- function(formula){
  return(.bt_fixed_formula(formula))
}
.get_grouping_factor    <- function(x){
  has_grouping            <- grepl("\\|", x)
  grouping                <- rep("", length(x))
  grouping[has_grouping]  <- trimws(sub(".*\\|\\s*", "", x[has_grouping]))
  return(grouping)
}
.remove_grouping_factor <- function(formula){
  return(trimws(sub("\\|.*$", "", formula)))
}
#' @title Add an Intercept to a Formula
#'
#' @description Converts a no-intercept formula to the corresponding formula
#' with an intercept while preserving the formula environment. Top-level
#' no-intercept encodings such as \code{- 1}, \code{+ 0}, and \code{0 +} are
#' removed without editing transformed calls such as \code{I(x - 1)} or
#' \code{offset(x - 1)}.
#'
#' @param formula a formula object.
#'
#' @return A formula object with an intercept.
#'
#' @export
formula_add_intercept <- function(formula){

  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula.", call. = FALSE)
  }

  if(attr(stats::terms(formula), "intercept") == 1L){
    return(formula)
  }

  formula_env   <- environment(formula)
  formula_attrs <- attributes(formula)
  rhs_index     <- if(length(formula) == 3L) 3L else 2L
  rhs           <- .formula_strip_no_intercept(formula[[rhs_index]])

  if(is.null(rhs)){
    rhs <- 1
  }

  out <- formula
  out[[rhs_index]] <- rhs
  environment(out) <- formula_env

  for(attribute in setdiff(names(formula_attrs), c("class", ".Environment", "names"))){
    attr(out, attribute) <- formula_attrs[[attribute]]
  }

  if(attr(stats::terms(out), "intercept") == 0L){
    out[[rhs_index]] <- call("+", 1, out[[rhs_index]])
    environment(out) <- formula_env
  }

  return(out)
}
.add_intercept_to_formula <- formula_add_intercept

.formula_strip_no_intercept <- function(expr){

  if(.formula_is_no_intercept_additive_term(expr)){
    return(NULL)
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L){
    lhs <- .formula_strip_no_intercept(expr[[2L]])
    rhs <- .formula_strip_no_intercept(expr[[3L]])

    if(is.null(lhs)){
      return(rhs)
    }
    if(is.null(rhs)){
      return(lhs)
    }
    return(call("+", lhs, rhs))
  }

  if(is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 3L &&
     .formula_is_numeric_constant(expr[[3L]], 1)){
    return(.formula_strip_no_intercept(expr[[2L]]))
  }

  return(expr)
}

.formula_is_no_intercept_additive_term <- function(expr){

  .formula_is_numeric_constant(expr, 0) || .formula_is_negative_one(expr)
}

.formula_is_negative_one <- function(expr){

  (is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), -1)) ||
    (is.call(expr) && identical(expr[[1L]], as.name("-")) && length(expr) == 2L &&
       .formula_is_numeric_constant(expr[[2L]], 1))
}

.formula_is_numeric_constant <- function(expr, value){

  if(is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L){
    return(.formula_is_numeric_constant(expr[[2L]], value))
  }

  is.numeric(expr) && length(expr) == 1L && identical(as.numeric(expr), as.numeric(value))
}

#' @title Evaluate JAGS formula using posterior samples
#'
#' @description Evaluates a JAGS formula on a posterior distribution obtained
#' from a fitted model. Formula random effects can be evaluated for existing
#' grouping levels when either standardized latent random effects and covariance
#' hyperparameters were monitored via \code{random_monitor(latent = TRUE)}, or
#' the group-level coefficients were monitored via
#' \code{random_monitor(coefficients = TRUE)}.
#' Row-indexed external random-effect SD sources, such as
#' \code{random_sd_source("tau", shape = "row")}, are evaluated from latent
#' random effects only. Model generation automatically monitors the required
#' latent effects for these blocks. The posterior samples must either contain
#' source columns named \code{tau[1]}, ..., \code{tau[nrow(data)]} aligned to
#' the supplied prediction rows, or the source must provide a
#' \code{parameter_source()} \code{values} function for reconstructing row-wise
#' source values from the posterior samples and supplied prediction data.
#'
#' @param fit model fitted with either \link[runjags]{runjags} posterior
#' samples obtained with \link[rjags]{rjags-package}
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored). If `NULL`, the fitted formula stored in
#' `formula_design` metadata is used. If the formula has a
#' \code{"log(intercept)"} attribute set to \code{TRUE}, the intercept values
#' will be log-transformed before computing the linear predictor.
#' @param parameter name of the parameter created with the formula
#' @param data data.frame containing predictors included in the formula. If
#' `NULL`, fitted source data from `formula_design` metadata are used.
#' @param prior_list named list of prior distribution of parameters specified
#' within the \code{formula}. If `NULL`, fitted priors from `formula_design`
#' metadata are used.
#' @param formula_target optional formula prediction target. `NULL` preserves
#' the historical safety behavior. `"fixed"` evaluates only the fixed formula
#' contribution. `"conditional"` evaluates fixed effects plus fitted or
#' explicitly generated random-effect contributions.
#' @param blocks optional random-effect block names used with
#' `formula_target = "conditional"`.
#' @param new_levels optional new-level policy used only with
#' `formula_target = "conditional"`. Use a `random_new_levels()` object or one
#' of `"error"`, `"zero"`, or `"sample"`.
#'
#'
#' @return \code{JAGS_evaluate_formula} returns a matrix of the evaluated posterior samples on
#' the supplied data.
#'
#' @seealso [JAGS_fit()] [JAGS_formula()]
#' @export
JAGS_evaluate_formula <- function(fit, formula = NULL, parameter,
                                  data = NULL, prior_list = NULL,
                                  formula_target = NULL, blocks = NULL,
                                  new_levels = NULL){

  check_char(parameter, "parameter")
  formula_target <- .bt_formula_prediction_target(
    formula_target,
    allow_marginal = FALSE,
    context = "JAGS_evaluate_formula()"
  )
  if(!is.null(blocks)){
    check_char(blocks, "blocks", check_length = 0, allow_NA = FALSE)
    if(anyDuplicated(blocks)){
      stop("'blocks' must be unique.", call. = FALSE)
    }
  }
  if(!is.null(blocks) && !identical(formula_target, "conditional")){
    stop("'blocks' can be used only with formula_target = \"conditional\".", call. = FALSE)
  }
  if(!is.null(new_levels) && !identical(formula_target, "conditional")){
    stop("'new_levels' can be used only with formula_target = \"conditional\".", call. = FALSE)
  }
  if(!is.null(new_levels)){
    new_levels <- .bt_random_new_levels_resolve(new_levels)
  }
  fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  resolved_inputs <- .bt_JAGS_evaluate_formula_resolve_inputs(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    fitted_design = fitted_design
  )
  formula <- resolved_inputs$formula
  data <- resolved_inputs$data
  prior_list <- resolved_inputs$prior_list
  fitted_design <- resolved_inputs$fitted_design

  if(!is.language(formula))
    stop("'formula' must be a formula")
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")

  # extract the posterior distribution
  posterior <- as.matrix(.fit_to_posterior(fit))

  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)
  formula_has_random <- .has_random_effects(formula)
  fitted_has_random <- !is.null(fitted_design) &&
    .bt_formula_design_has_any_random_effects(fitted_design)
  if(!is.null(blocks) && !fitted_has_random){
    stop(
      "The fitted formula for parameter '", parameter,
      "' does not include random-effect blocks.",
      call. = FALSE
    )
  }
  if(identical(formula_target, "fixed")){
    formula <- .remove_random_effects(formula)
  }else if(formula_has_random ||
           (identical(formula_target, "conditional") && fitted_has_random)){
    return(.bt_JAGS_evaluate_formula_with_random_effects(
      fit = fit,
      formula = formula,
      parameter = parameter,
      data = data,
      prior_list = prior_list,
      posterior = posterior,
      formula_target = formula_target,
      blocks = blocks,
      new_levels = new_levels
    ))
  }
  if(is.null(formula_target) &&
     !inherits(fitted_design, "try-error") &&
     .bt_formula_design_has_sampled_random_effects(fitted_design)){
    stop(
      "The fitted formula for parameter '", parameter,
      "' includes random effects. JAGS_evaluate_formula() cannot currently evaluate ",
      "random-effect fits without silently dropping group-level contributions.",
      call. = FALSE
    )
  }
  log_intercept <- isTRUE(attr(formula, "log(intercept)"))
  if(attr(stats::terms(formula), "intercept") == 0){
    formula <- formula_add_intercept(formula)
    if(log_intercept){
      attr(formula, "log(intercept)") <- TRUE
    }
  }

  # select priors corresponding to the prior distribution
  prior_parameter <- sapply(prior_list, function(p) if(is.null(attr(p, "parameter"))) "__none" else attr(p, "parameter"))
  if(!any(parameter %in% unique(prior_parameter)))
    stop("The specified parameter '", parameter, "' was not used in any of the prior distributions.")
  prior_list_formula <- prior_list[prior_parameter == parameter]
  names(prior_list_formula) <- format_parameter_names(names(prior_list_formula), formula_parameters = parameter, formula_prefix = FALSE)

  # extract the terms information from the formula
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))

  # check that all predictors have data and prior distribution
  if(!all(predictors %in% colnames(data)))
    stop(paste0("The ", paste0("'", predictors[!predictors %in% colnames(data)], "'", collapse = ", ")," predictor variable is missing in the data."))
  if(!all(model_terms %in% names(prior_list_formula)))
    stop(paste0("The prior distribution for the ", paste0("'", predictors[!model_terms %in% format_parameter_names(names(prior_list_formula), formula_parameters = parameter, formula_prefix = FALSE)], "'", collapse = ", ")," term is missing in the prior_list."))

  # obtain predictors characteristics -- based on prior distributions used to fit the original model
  # (i.e., do not truest the supplied data -- probably passed by the user)
  model_terms_type <- sapply(model_terms, function(model_term){
    if(model_term == "intercept"){
      return("continuous")
    }else if(is.prior.factor(prior_list_formula[[model_term]]) || inherits(prior_list_formula[[model_term]], "prior.factor_mixture") || inherits(prior_list_formula[[model_term]], "prior.factor_spike_and_slab")){
      return("factor")
    }else if(is.prior.simple(prior_list_formula[[model_term]]) || inherits(prior_list_formula[[model_term]], "prior.simple_mixture") || inherits(prior_list_formula[[model_term]], "prior.simple_spike_and_slab")){
      return("continuous")
    } else {
      stop(paste0("Unrecognized prior distribution for the '", model_term, "' term."))
    }
  })
  predictors_type <- model_terms_type[predictors]

  # check that passed data correspond to the specified priors (factor levels etc...) and set the proper contrasts
  if(any(predictors_type == "factor")){

    # check the proper data input for each factor prior
    for(factor in names(predictors_type[predictors_type == "factor"])){

      # select the corresponding prior in the variable
      this_prior <- prior_list_formula[[factor]]

      if(is.factor(data[,factor])){
        if(all(levels(data[,factor]) %in% .get_prior_factor_level_names(this_prior))){
          # either the formatting is correct, or the supplied levels are a subset of the original levels
          # reformat to check ordering and etc...
          data[,factor] <- factor(data[,factor], levels = .get_prior_factor_level_names(this_prior))
        }else{
          # there are some additional levels
          stop(paste0("Levels specified in the '", factor, "' factor variable do not match the levels used for model specification."))
        }
      }else if(all(unique(data[,factor]) %in% .get_prior_factor_level_names(this_prior))){
        # the variable was not passed as a factor but the values matches the factor levels
        data[,factor] <- factor(data[,factor], levels = .get_prior_factor_level_names(this_prior))
      }else{
        # there are some additional mismatching values
        stop(paste0("Levels specified in the '", factor, "' factor variable do not match the levels used for model specification."))
      }

      # set the contrast
      if(is.prior.orthonormal(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.orthonormal"
      }else if(is.prior.meandif(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.meandif"
      }else if(is.prior.independent(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.independent"
      }else if(is.prior.treatment(this_prior)){
        stats::contrasts(data[[factor]]) <- "contr.treatment"
      }
    }
  }
  if(any(predictors_type == "continuous")){

    # check the proper data input for each continuous prior
    for(continuous in names(predictors_type[predictors_type == "continuous"])){

      # select the corresponding prior in the variable
      this_prior <- prior_list_formula[[continuous]]

      if(is.prior.factor(this_prior)|| is.prior.discrete(this_prior) || is.prior.PET(this_prior) || is.prior.PEESE(this_prior) || is.prior.weightfunction(this_prior)){
        stop(paste0("Unsupported prior distribution defined for '", continuous, "' continuous variable. See '?prior' for details."))
      }
    }

    data <- .bt_apply_formula_scale_to_data(
      fit = fit,
      parameter = parameter,
      data = data,
      predictors_type = predictors_type
    )
  }

  # get the design matrix
  model_frame  <- tryCatch(
    stats::model.frame(formula, data = data, na.action = stats::na.pass),
    error = function(e){
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  if(anyNA(model_frame)){
    stop("Formula predictors contain missing values.", call. = FALSE)
  }
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)

  ### evaluate the design matrix on the samples -> output[data, posterior]
  if(has_intercept){

    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    # check for scaling factors
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term = "intercept", prior_list = prior_list_formula, posterior = posterior, nrow = nrow(data), ncol = nrow(posterior))

    # get intercept values and apply log() transformation if log(intercept) attribute is set
    if(is.prior.point(prior_list_formula[["intercept"]])){
      intercept_values <- rep(
        prior_list_formula[["intercept"]]$parameters[["location"]],
        nrow(posterior)
      )
    }else{
      intercept_values <- posterior[, JAGS_parameter_names("intercept", formula_parameter = parameter)]
    }
    if(log_intercept){
      intercept_values <- log(intercept_values)
    }
    output           <- temp_multiply_by * matrix(intercept_values, nrow = nrow(data), ncol = nrow(posterior), byrow = TRUE)

  }else{

    terms_indexes    <- attr(model_matrix, "assign")
    output           <- matrix(0, nrow = nrow(data), ncol = nrow(posterior))

  }

  # add remaining terms (omitting the intercept indexed as NA)
  for(i in unique(terms_indexes[terms_indexes > 0])){

    # subset the model matrix
    temp_data <- model_matrix[,terms_indexes == i,drop = FALSE]

    # get the posterior (unless point prior was used)
    if(is.prior.point(prior_list_formula[[model_terms[i]]])){
      temp_posterior <- matrix(
        prior_list_formula[[model_terms[i]]]$parameters[["location"]],
        nrow = nrow(posterior),
        ncol = if(model_terms_type[i] == "factor") .get_prior_factor_levels(prior_list_formula[[model_terms[i]]]) else 1
      )
    }else{
      temp_posterior <- posterior[,paste0(
        JAGS_parameter_names(model_terms[i], formula_parameter = parameter),
        if(model_terms_type[i] == "factor" && .get_prior_factor_levels(prior_list_formula[[model_terms[i]]]) > 1) paste0("[", 1:.get_prior_factor_levels(prior_list_formula[[model_terms[i]]]), "]"))
        ,drop = FALSE]
    }

    # check for scaling factors
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term = model_terms[i], prior_list = prior_list_formula, posterior = posterior, nrow = nrow(data), ncol = nrow(posterior))

    output <- output + temp_multiply_by * (temp_data %*% t(temp_posterior))

  }

  return(output)
}

.bt_formula_prediction_target <- function(formula_target,
                                          allow_marginal = TRUE,
                                          context = "Formula prediction"){

  if(is.null(formula_target)){
    return(NULL)
  }
  check_char(formula_target, "formula_target", check_length = 1,
             allow_NULL = FALSE, allow_NA = FALSE)
  allowed <- c("fixed", "conditional", if(isTRUE(allow_marginal)) "marginal")
  if(!formula_target %in% allowed){
    stop(
      context, " supports formula_target = ",
      paste0("'", allowed, "'", collapse = ", "),
      if(!isTRUE(allow_marginal)) ". Use JAGS_predict_formula() for formula_target = 'marginal'." else ".",
      call. = FALSE
    )
  }

  formula_target
}

.bt_JAGS_evaluate_formula_design <- function(fit, parameter){

  fitted_design <- try(JAGS_formula_design(fit, parameter), silent = TRUE)
  if(inherits(fitted_design, "try-error")){
    return(NULL)
  }

  fitted_design
}

.bt_JAGS_evaluate_formula_resolve_inputs <- function(fit, formula,
                                                     parameter, data,
                                                     prior_list,
                                                     fitted_design = NULL){

  if(is.null(fitted_design)){
    fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  }

  if(is.null(formula)){
    if(is.null(fitted_design)){
      stop("'formula' must be a formula.", call. = FALSE)
    }
    formula <- fitted_design$formula
    # Stored design formulas intentionally drop their original environment.
    # Use the user workspace as a pragmatic evaluation environment for replay.
    environment(formula) <- globalenv()
    if(isTRUE(fitted_design$log_intercept)){
      attr(formula, "log(intercept)") <- TRUE
    }
  }
  if(is.null(data)){
    if(is.null(fitted_design)){
      stop("'data' must be a data.frame.", call. = FALSE)
    }
    data <- fitted_design$source_data
    if(is.null(data)){
      data <- as.data.frame(fitted_design$model_frame)
    }
  }
  if(is.null(prior_list)){
    fit_prior_list <- attr(fit, "prior_list", exact = TRUE)
    if(is.list(fit_prior_list) && length(fit_prior_list) > 0L){
      prior_list <- fit_prior_list
    }else if(!is.null(fitted_design) && is.list(fitted_design$prior_list)){
      prior_list <- fitted_design$prior_list
    }else{
      stop("'prior_list' must be supplied.", call. = FALSE)
    }
  }

  list(
    formula = formula,
    data = data,
    prior_list = prior_list,
    fitted_design = fitted_design
  )
}

#' Formula prediction with explicit random-effect target
#'
#' @description
#' Returns formula-scale predictions from fitted BayesTools formula metadata.
#' Unlike [JAGS_evaluate_formula()], this helper can return structured
#' prediction output for marginal random-effect targets, including marginal
#' random-effect covariance or simulated random-effect draws.
#'
#' @param fit model fitted with [JAGS_fit()] or posterior samples carrying
#' `formula_design` metadata.
#' @param parameter formula parameter name.
#' @param formula optional formula. If `NULL`, the fitted formula for
#' `parameter` is used.
#' @param data optional prediction data. If `NULL`, fitted source data are used.
#' @param prior_list optional named prior list. If `NULL`, fitted priors are
#' used.
#' @param formula_target prediction target: `"conditional"`, `"fixed"`, or
#' `"marginal"`. The default includes fitted random-effect contributions.
#' @param blocks optional random-effect block names.
#' @param new_levels new-level policy for conditional or marginal random-effect
#' prediction. Use a `random_new_levels()` object or one of `"error"`, `"zero"`,
#' or `"sample"`.
#' @param marginal_method marginal representation used only with
#' `formula_target = "marginal"`: `"covariance"` returns fixed means plus
#' marginal random-effect covariance; `"sample"` returns materialized
#' random-effect draws.
#' @param seed optional random seed used when random-effect draws are simulated.
#' @param components whether to include component matrices where available.
#'
#' @return A list of class `BayesTools_formula_prediction` with fields `value`,
#' `mean`, `random`, `vcov`, `components`, and `metadata`.
#'
#' @seealso [JAGS_evaluate_formula()] [random_effects_marginal_vcov()]
#' @export
JAGS_predict_formula <- function(fit, parameter, formula = NULL, data = NULL,
                                 prior_list = NULL,
                                 formula_target = c("conditional", "fixed", "marginal"),
                                 blocks = NULL, new_levels = NULL,
                                 marginal_method = c("covariance", "sample"),
                                 seed = NULL, components = FALSE){

  check_char(parameter, "parameter")
  marginal_method_supplied <- !missing(marginal_method)
  formula_target <- match.arg(formula_target)
  marginal_method <- match.arg(marginal_method)
  if(!is.null(blocks)){
    check_char(blocks, "blocks", check_length = 0, allow_NA = FALSE)
    if(anyDuplicated(blocks)){
      stop("'blocks' must be unique.", call. = FALSE)
    }
  }
  if(!is.null(blocks) && identical(formula_target, "fixed")){
    stop(
      "'blocks' can be used only with formula_target = \"conditional\" or \"marginal\".",
      call. = FALSE
    )
  }
  if(!is.null(new_levels) && identical(formula_target, "fixed")){
    stop(
      "'new_levels' can be used only with formula_target = \"conditional\" or \"marginal\".",
      call. = FALSE
    )
  }
  if(isTRUE(marginal_method_supplied) &&
     !identical(formula_target, "marginal")){
    stop(
      "'marginal_method' can be used only with formula_target = \"marginal\".",
      call. = FALSE
    )
  }
  if(!is.null(new_levels)){
    new_levels <- .bt_random_new_levels_resolve(new_levels)
  }
  check_bool(components, "components", allow_NA = FALSE)
  if(!is.null(seed)){
    check_int(seed, "seed", lower = 0, check_length = 1, allow_NA = FALSE)
  }

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    formula_target = "fixed"
  )
  if(identical(formula_target, "fixed")){
    return(.bt_formula_prediction_object(
      value = fixed,
      mean = fixed,
      random = NULL,
      vcov = NULL,
      components = NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = NA_character_,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  if(identical(formula_target, "conditional")){
    seed_state <- .bt_formula_prediction_seed(seed)
    on.exit(.bt_formula_prediction_restore_seed(seed_state), add = TRUE)
    value <- JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = parameter,
      data = data,
      prior_list = prior_list,
      formula_target = "conditional",
      blocks = blocks,
      new_levels = new_levels
    )
    random <- value - fixed
    return(.bt_formula_prediction_object(
      value = value,
      mean = fixed,
      random = random,
      vcov = NULL,
      components = if(isTRUE(components)) list(fixed = fixed, random = random) else NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = NA_character_,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  posterior <- .bt_random_effect_marginal_covariance_posterior(
    fit = fit,
    posterior_samples = NULL
  )
  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = parameter
  )
  resolved_inputs <- .bt_JAGS_evaluate_formula_resolve_inputs(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    fitted_design = design
  )
  prior_list <- resolved_inputs$prior_list

  if(identical(marginal_method, "covariance")){
    vcov <- random_effects_marginal_vcov(
      fit = fit,
      parameter = design$parameter,
      data = data,
      posterior_samples = posterior,
      prior_list = prior_list,
      blocks = blocks,
      new_levels = new_levels
    )
    return(.bt_formula_prediction_object(
      value = fixed,
      mean = fixed,
      random = NULL,
      vcov = vcov,
      components = if(isTRUE(components)) list(fixed = fixed) else NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = marginal_method,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  seed_state <- .bt_formula_prediction_seed(seed)
  on.exit(.bt_formula_prediction_restore_seed(seed_state), add = TRUE)
  random <- .bt_random_effects_marginal_sample(
    design = design,
    posterior = posterior,
    prior_list = prior_list,
    data = data,
    blocks = blocks,
    new_levels = new_levels
  )
  value <- fixed + random
  .bt_formula_prediction_object(
    value = value,
    mean = fixed,
    random = random,
    vcov = NULL,
    components = if(isTRUE(components)) list(fixed = fixed, random = random) else NULL,
    metadata = .bt_formula_prediction_metadata(
      parameter = parameter,
      formula_target = formula_target,
      marginal_method = marginal_method,
      blocks = blocks,
      new_levels = new_levels
    )
  )
}

.bt_formula_prediction_object <- function(value, mean, random, vcov,
                                          components, metadata){

  out <- list(
    value = value,
    mean = mean,
    random = random,
    vcov = vcov,
    components = components,
    metadata = metadata
  )
  class(out) <- c("BayesTools_formula_prediction", "list")

  out
}

.bt_formula_prediction_metadata <- function(parameter, formula_target,
                                            marginal_method, blocks,
                                            new_levels){

  list(
    parameter = parameter,
    formula_target = formula_target,
    marginal_method = marginal_method,
    blocks = blocks,
    new_levels = new_levels
  )
}

.bt_formula_prediction_seed <- function(seed){

  if(is.null(seed)){
    return(NULL)
  }
  state <- list(
    exists = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    value = if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }else{
      NULL
    }
  )
  set.seed(seed)
  state
}

.bt_formula_prediction_restore_seed <- function(state){

  if(is.null(state)){
    return(invisible(NULL))
  }
  if(isTRUE(state$exists)){
    assign(".Random.seed", state$value, envir = .GlobalEnv)
  }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(NULL)
}

.bt_random_effects_marginal_sample <- function(design, posterior, prior_list,
                                               data = NULL, blocks = NULL,
                                               new_levels = NULL){

  selected <- .bt_random_effect_marginal_covariance_terms(
    design = design,
    blocks = blocks
  )
  output <- NULL
  for(random_term in selected$terms){
    block_new_levels <- .bt_random_effect_new_levels_policy(
      random_term = random_term,
      override = new_levels
    )
    block_data <- .bt_random_effect_marginal_covariance_block_data(
      design = design,
      random_term = random_term,
      data = data
    )
    new_row <- block_data$group_map > random_term$n_groups
    if(any(new_row) && !isTRUE(block_new_levels$allow)){
      new_groups <- block_data$group_levels[unique(block_data$group_map[new_row])]
      stop(
        "New random-effect level(s) for block '", random_term$block_name,
        "' require an explicit new-level policy: ",
        paste(new_groups, collapse = ", "),
        ". Use new_levels = \"zero\" or new_levels = \"sample\".",
        call. = FALSE
      )
    }
    rows <- seq_len(nrow(block_data$model_matrix))
    if(any(new_row) && identical(block_new_levels$method, "zero")){
      rows <- rows[!new_row]
    }
    contribution <- .bt_random_effect_group_contribution_sample(
      random_term = random_term,
      model_matrix = block_data$model_matrix,
      group_map = block_data$group_map,
      rows = rows,
      posterior = posterior,
      prior_list = prior_list,
      source_data = block_data$source_data
    )
    if(is.null(output)){
      output <- contribution
    }else{
      if(!identical(dim(output), dim(contribution))){
        stop(
          "Random-effect marginal sample block '", random_term$block_name,
          "' produced incompatible dimensions.",
          call. = FALSE
        )
      }
      output <- output + contribution
    }
  }

  output
}

.bt_apply_formula_scale_to_data <- function(fit, parameter, data,
                                            predictors_type){

  continuous_predictors <- names(predictors_type[predictors_type == "continuous"])

  formula_scale <- attr(fit, "formula_scale")
  if(is.null(formula_scale)){
    return(data)
  }

  param_scale <- formula_scale[[parameter]]
  if(is.null(param_scale)){
    return(data)
  }

  scaled_predictors <- sub(paste0("^", parameter, "_"), "", names(param_scale))
  continuous_predictors <- unique(c(continuous_predictors, scaled_predictors))
  if(length(continuous_predictors) == 0L){
    return(data)
  }

  for(continuous in continuous_predictors){
    if(!continuous %in% colnames(data)){
      next
    }
    scaled_name <- paste0(parameter, "_", continuous)
    if(scaled_name %in% names(param_scale)){
      scale_info <- param_scale[[scaled_name]]
      data[, continuous] <- (data[, continuous] - scale_info$mean) / scale_info$sd
    }
  }

  data
}

.bt_JAGS_evaluate_formula_with_random_effects <- function(fit, formula,
                                                          parameter, data,
                                                          prior_list,
                                                          posterior,
                                                          formula_target = NULL,
                                                          blocks = NULL,
                                                          new_levels = NULL){

  fitted_design <- .bt_JAGS_evaluate_formula_design(fit, parameter)
  if(is.null(fitted_design)){
    stop(
      "JAGS_evaluate_formula() needs fitted formula design metadata to evaluate random effects. ",
      "Use a fit produced by JAGS_fit() with formula_list.",
      call. = FALSE
    )
  }
  if(!.bt_formula_design_has_any_random_effects(fitted_design)){
    stop(
      "The supplied formula contains random effects, but the fitted formula for parameter '",
      parameter, "' does not.",
      call. = FALSE
    )
  }

  random_terms <- if(.has_random_effects(formula)){
    .bt_parse_random_effects(formula)$terms
  }else{
    list()
  }
  if(length(random_terms) > 0L){
    .bt_validate_random_effect_prediction_terms(
      requested = random_terms,
      fitted = .bt_formula_design_random_effects(fitted_design),
      parameter = parameter
    )
  }
  selected_random_effects <- .bt_JAGS_evaluate_formula_random_effect_terms(
    fitted_design = fitted_design,
    requested_terms = random_terms,
    formula_target = formula_target,
    blocks = blocks,
    parameter = parameter
  )

  fixed_formula <- .remove_random_effects(formula)
  output <- JAGS_evaluate_formula(
    fit = fit,
    formula = fixed_formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    formula_target = "fixed"
  )
  random_data <- .bt_apply_formula_scale_to_data(
    fit = fit,
    parameter = parameter,
    data = data,
    predictors_type = fitted_design$predictor_types
  )

  for(random_term in selected_random_effects){
    random_structure <- .bt_random_effect_structure(
      random_term,
      context = "Random-effect prediction metadata"
    )
    output <- output + .bt_JAGS_evaluate_random_effect_term(
      random_term = random_term,
      data = if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")) data else random_data,
      group_data = data,
      posterior = posterior,
      prior_list = prior_list,
      new_levels = new_levels
    )
  }

  output
}

.bt_JAGS_evaluate_formula_random_effect_terms <- function(fitted_design,
                                                          requested_terms,
                                                          formula_target,
                                                          blocks,
                                                          parameter){

  fitted_terms <- .bt_formula_design_random_effects(fitted_design)
  fitted_names <- vapply(fitted_terms, `[[`, character(1), "block_name")
  if(!is.null(blocks)){
    unknown <- setdiff(blocks, fitted_names)
    if(length(unknown) > 0L){
      stop(
        "Unknown random-effect block(s): ",
        paste(unknown, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if(!identical(formula_target, "conditional")){
    return(.bt_formula_design_sampled_random_effects(fitted_design))
  }

  requested_names <- vapply(requested_terms, `[[`, character(1), "block_name")
  selected_names <- if(!is.null(blocks)){
    blocks
  }else if(length(requested_names) > 0L){
    requested_names
  }else{
    fitted_names
  }
  selected_terms <- fitted_terms[match(selected_names, fitted_names)]
  modes <- vapply(selected_terms, .bt_random_effect_term_compile_mode, character(1))
  marginalized <- selected_names[modes == "marginalized"]
  if(length(marginalized) > 0L){
    stop(
      "JAGS_evaluate_formula() cannot use formula_target = \"conditional\" for ",
      "random-effect block(s) compiled as marginalized: ",
      paste(marginalized, collapse = ", "),
      ". Use formula_target = \"marginal\" with JAGS_predict_formula() or refit ",
      "with the block(s) sampled.",
      call. = FALSE
    )
  }

  selected_terms
}

.bt_validate_random_effect_prediction_terms <- function(requested, fitted,
                                                        parameter){

  requested_names <- vapply(requested, function(term) term$block_name, character(1))
  if(anyDuplicated(requested_names)){
    stop(
      "Random-effect block names in the supplied formula must be unique.",
      call. = FALSE
    )
  }
  fitted_names <- vapply(fitted, function(term) term$block_name, character(1))
  unknown_names <- setdiff(requested_names, fitted_names)
  if(length(unknown_names) > 0L){
    stop(
      "Random-effect block(s) in the supplied formula do not match the fitted formula for parameter '",
      parameter, "'; block(s) were not found in the fitted formula: ",
      paste(unknown_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for(block in requested_names){
    requested_term <- requested[[match(block, requested_names)]]
    fitted_term <- fitted[[match(block, fitted_names)]]
    requested_structure <- .bt_random_effect_structure(
      requested_term,
      context = "Random-effect prediction metadata"
    )
    fitted_structure <- .bt_random_effect_structure(
      fitted_term,
      context = "Random-effect prediction metadata"
    )
    .bt_validate_random_effect_term_supported(requested_term)
    requested_homogeneous <- .bt_random_effect_homogeneous_sd(
      requested_term,
      requested_structure
    )
    fitted_homogeneous <- .bt_random_effect_homogeneous_sd_metadata(
      fitted_term,
      context = "Random-effect prediction metadata"
    )
    if(!identical(requested_structure, fitted_structure) ||
       !identical(requested_homogeneous, fitted_homogeneous) ||
       !identical(.bt_deparse_expr(requested_term$expr), .bt_deparse_expr(fitted_term$expr)) ||
       !identical(requested_term$group_label, fitted_term$group_label)){
      stop(
        "Random-effect block '", block,
        "' in the supplied formula does not match the fitted formula.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_JAGS_evaluate_random_effect_term <- function(random_term, data, posterior,
                                                prior_list,
                                                group_data = data,
                                                new_levels = NULL){

  new_levels <- .bt_random_effect_new_levels_policy(random_term, new_levels)
  prediction <- .bt_random_effect_prediction_data(
    random_term,
    data,
    group_data = group_data,
    allow_new_groups = isTRUE(new_levels$allow)
  )
  model_matrix <- prediction$model_matrix
  group_map <- prediction$group_map
  n_rows <- nrow(model_matrix)
  n_draws <- nrow(posterior)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  fitted_n_groups <- length(random_term$group_levels)
  new_row <- group_map > fitted_n_groups

  if(any(!new_row)){
    existing_rows <- which(!new_row)
    output[existing_rows, ] <- .bt_JAGS_evaluate_random_effect_term_existing(
      random_term = random_term,
      model_matrix = model_matrix[existing_rows, , drop = FALSE],
      group_map = group_map[existing_rows],
      posterior = posterior,
      prior_list = prior_list,
      group_data = group_data[existing_rows, , drop = FALSE]
    )
  }
  if(any(new_row) && identical(new_levels$method, "sample")){
    output <- output + .bt_random_effect_new_level_contribution_sample(
      random_term = random_term,
      prediction = prediction,
      new_row = new_row,
      posterior = posterior,
      prior_list = prior_list,
      source_data = group_data
    )
  }

  output
}

.bt_JAGS_evaluate_random_effect_term_existing <- function(random_term,
                                                          model_matrix,
                                                          group_map,
                                                          posterior,
                                                          prior_list,
                                                          group_data){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  n_groups <- length(random_term$group_levels)

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(
      .bt_random_effect_row_indexed_contribution_from_latent(
        random_term = random_term,
        model_matrix = model_matrix,
        group_map = group_map,
        posterior = posterior,
        prior_list = prior_list,
        data = group_data,
        context = "Prediction"
      )
    )
  }

  coefficient_names <- .bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(all(as.vector(coefficient_names) %in% colnames(posterior))){
    return(.bt_random_effect_contribution_from_coefficients(
      model_matrix = model_matrix,
      group_map = group_map,
      coefficient_names = coefficient_names,
      posterior = posterior
    ))
  }

  latent_contribution <- .bt_try_random_effect_contribution_from_latent(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    posterior = posterior,
    prior_list = prior_list
  )
  if(!is.null(latent_contribution)){
    return(latent_contribution)
  }

  stop(
    "Random-effect coefficients for block '", random_term$block_name,
    "' cannot be reconstructed from the posterior samples. Refit with ",
    "random_monitor(latent = TRUE) or random_monitor(coefficients = TRUE) ",
    "for that random-effect block before using JAGS_evaluate_formula() ",
    "with random effects.",
    call. = FALSE
  )
}

.bt_random_effect_contribution_from_coefficients <- function(model_matrix,
                                                            group_map,
                                                            coefficient_names,
                                                            posterior){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(column in seq_len(n_columns)){
    coefficient_matrix <- posterior[, coefficient_names[, column], drop = FALSE]
    output <- output +
      t(coefficient_matrix[, group_map, drop = FALSE]) *
      matrix(model_matrix[, column], nrow = n_rows, ncol = n_draws)
  }

  output
}

.bt_random_effect_new_level_contribution_sample <- function(random_term,
                                                            prediction,
                                                            new_row,
                                                            posterior,
                                                            prior_list,
                                                            source_data){

  .bt_random_effect_group_contribution_sample(
    random_term = random_term,
    model_matrix = prediction$model_matrix,
    group_map = prediction$group_map,
    rows = which(new_row),
    posterior = posterior,
    prior_list = prior_list,
    source_data = source_data
  )
}

.bt_random_effect_group_contribution_sample <- function(random_term,
                                                        model_matrix,
                                                        group_map,
                                                        rows = seq_len(nrow(model_matrix)),
                                                        posterior,
                                                        prior_list,
                                                        source_data){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  if(length(rows) == 0L){
    return(output)
  }

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(.bt_random_effect_group_contribution_sample_row_indexed(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      rows = rows,
      posterior = posterior,
      prior_list = prior_list,
      source_data = source_data
    ))
  }

  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    .bt_random_effect_marginal_covariance_missing_sd_stop(
      random_term = random_term,
      n_columns = n_columns
    )
  }
  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    G <- matrix(
      correlation[draw, , ],
      nrow = n_columns,
      ncol = n_columns
    ) * tcrossprod(sd_draws[draw, ])
    effects <- .bt_random_effect_mvn_group_draws(
      covariance = G,
      n_groups = length(groups)
    )
    output[rows, draw] <- rowSums(
      model_matrix[rows, , drop = FALSE] *
        effects[group_index, , drop = FALSE]
    )
  }

  output
}

.bt_random_effect_group_contribution_sample_row_indexed <- function(
    random_term,
    model_matrix,
    group_map,
    rows,
    posterior,
    prior_list,
    source_data){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = n_rows,
    posterior = posterior,
    data = source_data,
    context = "Prediction"
  )
  column_allocation <- .bt_random_effect_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = n_columns
  )
  if(is.null(column_allocation)){
    allocation <- .bt_random_effect_row_indexed_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list
    )
  }else{
    allocation <- NULL
  }
  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  groups <- sort(unique(group_map[rows]))
  group_index <- match(group_map[rows], groups)
  for(draw in seq_len(n_draws)){
    effects <- .bt_random_effect_mvn_group_draws(
      covariance = matrix(correlation[draw, , ], nrow = n_columns, ncol = n_columns),
      n_groups = length(groups)
    )
    Z <- model_matrix[rows, , drop = FALSE] *
      matrix(source_draws[draw, rows], nrow = length(rows), ncol = n_columns)
    if(is.null(column_allocation)){
      Z <- Z * allocation[draw]
    }else{
      Z <- Z * matrix(
        column_allocation[draw, ],
        nrow = length(rows),
        ncol = n_columns,
        byrow = TRUE
      )
    }
    output[rows, draw] <- rowSums(Z * effects[group_index, , drop = FALSE])
  }

  output
}

.bt_random_effect_mvn_group_draws <- function(covariance, n_groups){

  if(n_groups == 0L){
    return(matrix(numeric(), nrow = 0L, ncol = ncol(covariance)))
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  transform <- decomposition$vectors %*%
    (sqrt(values) * t(decomposition$vectors))
  z <- matrix(
    stats::rnorm(n_groups * ncol(covariance)),
    nrow = n_groups,
    ncol = ncol(covariance)
  )

  z %*% transform
}

.bt_random_effect_prediction_data <- function(random_term, data,
                                              group_data = data,
                                              allow_new_groups = FALSE){

  prediction_data <- data
  random_structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect prediction metadata"
  )
  prediction_data <- .bt_random_effect_prediction_structured_index_data(
    random_term,
    prediction_data
  )
  factor_levels <- random_term$xlevels
  if(!is.null(factor_levels) && length(factor_levels) > 0L){
    for(factor_name in names(factor_levels)){
      if(!factor_name %in% names(prediction_data)){
        stop(
          "The '", factor_name,
          "' predictor needed for random-effect prediction is missing in the data.",
          call. = FALSE
        )
      }
      if(is.factor(prediction_data[[factor_name]])){
        observed_levels <- unique(as.character(prediction_data[[factor_name]]))
        if(!all(observed_levels %in% factor_levels[[factor_name]])){
          stop(
            "Levels specified in the '", factor_name,
            "' factor variable do not match the levels used for model specification.",
            call. = FALSE
          )
        }
        prediction_data[[factor_name]] <- factor(
          as.character(prediction_data[[factor_name]]),
          levels = factor_levels[[factor_name]]
        )
      }else if(all(unique(prediction_data[[factor_name]]) %in% factor_levels[[factor_name]])){
        prediction_data[[factor_name]] <- factor(
          prediction_data[[factor_name]],
          levels = factor_levels[[factor_name]]
        )
      }else{
        stop(
          "Levels specified in the '", factor_name,
          "' factor variable do not match the levels used for model specification.",
          call. = FALSE
        )
      }
    }
  }

  contrasts <- random_term$contrasts
  if(!is.null(contrasts) && length(contrasts) > 0L){
    for(factor_name in names(contrasts)){
      if(factor_name %in% names(prediction_data) && is.factor(prediction_data[[factor_name]])){
        stats::contrasts(prediction_data[[factor_name]]) <- contrasts[[factor_name]]
      }
    }
  }

  random_design <- .bt_random_effect_design_matrix(
    random_term$term_formula,
    prediction_data,
    preserve_no_intercept_contrasts = !random_structure %in% c("cs", "hcs", "ar1", "car", "har"),
    structure = random_structure,
    car_time_values = if(identical(random_structure, "car")) random_term$car$time_values else NULL,
    block_name = random_term$block_name
  )
  model_matrix <- random_design$model_matrix
  colnames(model_matrix) <- gsub(":", "__xXx__", colnames(model_matrix))

  if(!identical(colnames(model_matrix), random_term$column_names)){
    stop(
      "Random-effect design columns for block '", random_term$block_name,
      "' do not match the fitted formula.",
      call. = FALSE
    )
  }

  grouping_values <- .bt_random_group_values(random_term, group_data)
  if(length(grouping_values) != nrow(model_matrix)){
    stop(
      "Random-effect grouping data for block '", random_term$block_name,
      "' must have one value per prediction row.",
      call. = FALSE
    )
  }
  group_levels <- random_term$group_levels
  group_map <- match(as.character(grouping_values), group_levels)
  if(any(is.na(group_map))){
    new_groups <- unique(as.character(grouping_values)[is.na(group_map)])
    if(isTRUE(allow_new_groups)){
      group_levels <- c(group_levels, new_groups)
      group_map <- match(as.character(grouping_values), group_levels)
      return(list(
        model_matrix = model_matrix,
        group_map = group_map,
        group_levels = group_levels
      ))
    }
    stop(
      "New random-effect level(s) for block '", random_term$block_name,
      "' are not supported by JAGS_evaluate_formula(): ",
      paste(new_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    model_matrix = model_matrix,
    group_map = group_map,
    group_levels = group_levels
  )
}

.bt_random_effect_prediction_structured_index_data <- function(random_term,
                                                               data){

  index <- random_term$structured_index
  if(is.null(index)){
    return(data)
  }

  missing_variables <- index$variables[!index$variables %in% names(data)]
  if(length(missing_variables) > 0L){
    stop(
      "The ",
      paste0("'", missing_variables, "'", collapse = ", "),
      " structured random-effect index variable is missing in the data.",
      call. = FALSE
    )
  }

  data[[index$name]] <- .bt_random_effect_structured_index_values(
    data,
    index$variables
  )
  data
}

.get_parameter_scaling_factor_matrix <- function(term, prior_list, posterior, nrow, ncol){

  if(!is.null(attr(prior_list[[term]], "multiply_by"))){
    if(is.numeric(attr(prior_list[[term]], "multiply_by"))){
      temp_multiply_by <- matrix(attr(prior_list[[term]], "multiply_by"), nrow = nrow, ncol = ncol)
    }else{
      temp_multiply_by <- matrix(posterior[,JAGS_parameter_names(attr(prior_list[[term]], "multiply_by"))], nrow = nrow, ncol = ncol, byrow = TRUE)
    }
  }else{
    temp_multiply_by <- matrix(1, nrow = nrow, ncol = ncol)
  }

  return(temp_multiply_by)
}

.factor_level_list <- function(x){

  level_names <- attr(x, "level_names")
  if(is.null(level_names)){
    factor_terms <- attr(x, "factor_terms")
    if(!is.null(factor_terms) && length(factor_terms) > 1){
      return(NULL)
    }

    n_levels <- attr(x, "levels")
    if(is.null(n_levels)){
      return(NULL)
    }

    if(is.prior.factor(x)){
      level_names <- .get_prior_factor_level_names(x)
    }else if(isTRUE(attr(x, "independent"))){
      level_names <- seq_len(n_levels)
    }else{
      level_names <- seq_len(n_levels + 1)
    }
  }

  if(is.list(level_names)){
    factor_terms <- attr(x, "factor_terms")
    if(is.null(factor_terms)){
      factor_terms <- names(level_names)
    }
    if(is.null(factor_terms) || any(!nzchar(factor_terms))){
      factor_terms <- paste0("factor", seq_along(level_names))
    }
    level_names <- level_names[factor_terms]
    names(level_names) <- factor_terms
  }else{
    factor_terms <- attr(x, "factor_terms")
    if(is.null(factor_terms) || length(factor_terms) != 1){
      factor_terms <- ".factor"
    }
    level_names <- setNames(list(level_names), factor_terms)
  }

  level_names <- lapply(level_names, as.character)
  return(level_names)
}

.factor_cell_grid <- function(level_names){

  if(is.null(level_names) || length(level_names) == 0){
    return(data.frame())
  }

  expand.grid(level_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

.factor_cell_labels <- function(level_names){

  level_grid <- .factor_cell_grid(level_names)
  if(ncol(level_grid) == 0){
    return(character(0))
  }

  if(ncol(level_grid) == 1){
    return(as.character(level_grid[[1]]))
  }

  apply(level_grid, 1, function(level_row) {
    paste0(names(level_row), "=", unname(level_row), collapse = ", ")
  })
}

.factor_contrast_parameter_names <- function(parameter, level_names, cell_names){

  if(is.null(level_names)){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  level_grid <- .factor_cell_grid(level_names)
  if(nrow(level_grid) != length(cell_names)){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  if(length(level_names) == 1){
    return(paste0(parameter, "[dif: ", level_grid[[1]], "]"))
  }

  parameter_terms <- strsplit(parameter, "__xXx__", fixed = TRUE)[[1]]
  factor_terms <- names(level_names)
  factor_positions <- vapply(factor_terms, function(factor_term) {
    factor_position <- which(
      parameter_terms == factor_term |
        endsWith(parameter_terms, paste0("_", factor_term))
    )

    if(length(factor_position) == 1){
      return(factor_position)
    }

    return(NA_integer_)
  }, integer(1))

  if(!all(!is.na(factor_positions))){
    return(paste0(parameter, "[dif: ", cell_names, "]"))
  }

  vapply(seq_len(nrow(level_grid)), function(level_i) {
    formatted_terms <- parameter_terms

    for(factor_term in factor_terms){
      factor_position <- factor_positions[[factor_term]]
      prefix_length <- nchar(formatted_terms[[factor_position]]) - nchar(factor_term)
      formatted_terms[[factor_position]] <- paste0(
        substr(formatted_terms[[factor_position]], 1, prefix_length),
        factor_term,
        "[dif: ",
        level_grid[[factor_term]][level_i],
        "]"
      )
    }

    paste0(formatted_terms, collapse = "__xXx__")
  }, character(1))
}

.factor_object_contrast_name <- function(x){

  if(is.prior.independent(x) || isTRUE(attr(x, "independent"))){
    return("contr.independent")
  }else if(is.prior.treatment(x) || isTRUE(attr(x, "treatment"))){
    return("contr.treatment")
  }else if(is.prior.orthonormal(x) || isTRUE(attr(x, "orthonormal"))){
    return("contr.orthonormal")
  }else if(is.prior.meandif(x) || isTRUE(attr(x, "meandif"))){
    return("contr.meandif")
  }

  return(NULL)
}

.factor_contrast_matrix <- function(level_names, contrast){

  switch(
    contrast,
    "contr.treatment"   = stats::contr.treatment(level_names),
    "contr.independent" = contr.independent(level_names),
    "contr.orthonormal" = contr.orthonormal(level_names),
    "contr.meandif"     = contr.meandif(level_names),
    stop("Unsupported factor contrast '", contrast, "'.", call. = FALSE)
  )
}

.complete_factor_metadata_prior_list <- function(prior_list){

  if(is.null(prior_list) || length(prior_list) == 0){
    return(prior_list)
  }

  for(parameter in names(prior_list)){
    prior_list[[parameter]] <- .complete_factor_metadata(prior_list[[parameter]], parameter)
  }

  return(prior_list)
}

.copy_missing_factor_metadata <- function(target, source){

  metadata_names <- c(
    "levels",
    "level_names",
    "interaction",
    "interaction_terms",
    "term_components",
    "factor_terms",
    "factor_contrasts",
    "factor_design",
    "factor_cell_names"
  )

  for(metadata_name in metadata_names){
    if(is.null(attr(target, metadata_name, exact = TRUE)) &&
       !is.null(attr(source, metadata_name, exact = TRUE))){
      attr(target, metadata_name) <- attr(source, metadata_name, exact = TRUE)
    }
  }

  return(target)
}

.complete_factor_metadata <- function(x, parameter = NULL){

  if((is.prior.mixture(x) || is.prior.spike_and_slab(x)) && length(x) > 0){
    prior_components <- vapply(x, is.prior, logical(1))
    for(component_i in which(prior_components)){
      x[[component_i]] <- .complete_factor_metadata(x[[component_i]], parameter)
    }

    factor_components <- which(vapply(x, is.prior.factor, logical(1)))
    if(length(factor_components) > 0){
      x <- .copy_missing_factor_metadata(x, x[[factor_components[[1]]]])
    }
  }

  if(!is.prior.factor(x)){
    return(x)
  }

  level_names <- .factor_level_list(x)
  if(is.null(level_names) || length(level_names) != 1L){
    return(x)
  }

  factor_terms <- attr(x, "factor_terms", exact = TRUE)
  if(is.null(factor_terms) || length(factor_terms) != 1L ||
     anyNA(factor_terms) || !nzchar(factor_terms)){
    factor_terms <- if(!is.null(parameter) && nzchar(parameter)){
      parameter
    }else{
      ".factor"
    }
    attr(x, "factor_terms") <- factor_terms
  }

  contrast <- .factor_object_contrast_name(x)
  if(is.null(contrast)){
    return(x)
  }

  factor_contrasts <- attr(x, "factor_contrasts", exact = TRUE)
  if(is.null(factor_contrasts)){
    factor_contrasts <- stats::setNames(contrast, factor_terms)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
    if(any(is.na(factor_contrasts))){
      factor_contrasts[is.na(factor_contrasts)] <- contrast
    }
  }
  attr(x, "factor_contrasts") <- factor_contrasts

  if(is.null(attr(x, "factor_design", exact = TRUE)) ||
     is.null(attr(x, "factor_cell_names", exact = TRUE))){
    design_info <- .factor_term_design_from_metadata(x)
    attr(x, "factor_design")     <- design_info[["design"]]
    attr(x, "factor_cell_names") <- design_info[["cell_names"]]
  }

  return(x)
}

.add_factor_metadata_from_named_objects <- function(x, parameter, objects){

  level_names <- .factor_level_list(x)
  if(is.null(level_names)){
    return(x)
  }

  factor_terms <- names(level_names)
  if(is.null(attr(x, "factor_terms"))){
    attr(x, "factor_terms") <- factor_terms
  }

  factor_contrasts <- attr(x, "factor_contrasts")
  if(is.null(factor_contrasts)){
    return(x)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
  }

  attr(x, "factor_contrasts") <- factor_contrasts
  return(x)
}

.factor_term_design_from_formula <- function(formula, data, predictors, predictors_type, term_index, term_components, factor_terms, has_intercept){

  level_names <- lapply(factor_terms, function(factor_term) levels(data[[factor_term]]))
  names(level_names) <- factor_terms
  cell_grid <- .factor_cell_grid(level_names)

  grid_data <- data[rep(1, nrow(cell_grid)), predictors, drop = FALSE]
  rownames(grid_data) <- NULL

  for(predictor in predictors){
    if(predictors_type[[predictor]] == "factor"){
      predictor_values <- if(predictor %in% factor_terms){
        cell_grid[[predictor]]
      }else{
        rep(levels(data[[predictor]])[1], nrow(cell_grid))
      }
      grid_data[[predictor]] <- factor(predictor_values, levels = levels(data[[predictor]]))
      stats::contrasts(grid_data[[predictor]]) <- attr(data[[predictor]], "contrasts")
    }else{
      grid_data[[predictor]] <- if(predictor %in% term_components) 1 else 0
    }
  }

  grid_model_frame <- stats::model.frame(formula, data = grid_data)
  grid_model_matrix <- stats::model.matrix(grid_model_frame, formula = formula, data = grid_data)
  grid_terms_indexes <- attr(grid_model_matrix, "assign")
  if(has_intercept){
    grid_terms_indexes <- grid_terms_indexes + 1
    grid_terms_indexes[1] <- 0
  }

  design <- grid_model_matrix[, grid_terms_indexes == term_index, drop = FALSE]

  return(list(
    design     = unname(design),
    cell_grid  = cell_grid,
    cell_names = .factor_cell_labels(level_names),
    level_names = level_names
  ))
}

.factor_term_design_from_metadata <- function(x){

  factor_design <- attr(x, "factor_design")
  if(!is.null(factor_design)){
    factor_design <- as.matrix(factor_design)
    cell_names <- attr(x, "factor_cell_names")
    level_names <- .factor_level_list(x)
    if(is.null(cell_names)){
      if(is.null(level_names)){
        stop("Factor level names are missing and the factor contrast cannot be transformed.", call. = FALSE)
      }
      cell_names <- .factor_cell_labels(level_names)
    }
    return(list(
      design     = factor_design,
      cell_names = cell_names,
      level_names = level_names
    ))
  }

  level_names <- .factor_level_list(x)
  if(is.null(level_names)){
    stop("Factor level names are missing and the factor contrast cannot be transformed.", call. = FALSE)
  }

  factor_terms <- names(level_names)
  factor_contrasts <- attr(x, "factor_contrasts")
  if(is.null(factor_contrasts)){
    stop("Factor contrast metadata is missing and cannot be inferred.", call. = FALSE)
  }else{
    factor_contrasts <- as.character(factor_contrasts)
    if(is.null(names(factor_contrasts))){
      names(factor_contrasts) <- factor_terms[seq_along(factor_contrasts)]
    }
    factor_contrasts <- factor_contrasts[factor_terms]
  }

  if(any(is.na(factor_contrasts))){
    stop("Factor contrast metadata is incomplete and cannot be inferred.", call. = FALSE)
  }

  contrast_matrices <- lapply(factor_terms, function(factor_term) {
    .factor_contrast_matrix(level_names[[factor_term]], factor_contrasts[[factor_term]])
  })
  names(contrast_matrices) <- factor_terms

  level_grid <- expand.grid(lapply(level_names, seq_along), KEEP.OUT.ATTRS = FALSE)
  coef_grid <- expand.grid(lapply(contrast_matrices, function(contrast_matrix) seq_len(ncol(contrast_matrix))), KEEP.OUT.ATTRS = FALSE)

  design <- matrix(NA_real_, nrow = nrow(level_grid), ncol = nrow(coef_grid))
  for(row_i in seq_len(nrow(level_grid))){
    for(col_i in seq_len(nrow(coef_grid))){
      design[row_i, col_i] <- prod(vapply(factor_terms, function(factor_term) {
        contrast_matrices[[factor_term]][level_grid[[factor_term]][row_i], coef_grid[[factor_term]][col_i]]
      }, numeric(1)))
    }
  }

  return(list(
    design     = design,
    cell_names = .factor_cell_labels(level_names),
    level_names = level_names
  ))
}

.transform_factor_contrast_samples <- function(coefficient_samples, metadata, parameter, transformed_class){

  if(!is.matrix(coefficient_samples)){
    coefficient_samples <- matrix(coefficient_samples, ncol = 1)
  }

  design_info <- .factor_term_design_from_metadata(metadata)
  design <- design_info[["design"]]

  if(ncol(coefficient_samples) != ncol(design)){
    stop(
      "The factor contrast design for '", parameter, "' has ", ncol(design),
      " coefficient columns, but the samples contain ", ncol(coefficient_samples), ".",
      call. = FALSE
    )
  }

  transformed_samples <- coefficient_samples %*% t(design)
  colnames(transformed_samples) <- .factor_contrast_parameter_names(
    parameter = parameter,
    level_names = design_info[["level_names"]],
    cell_names = design_info[["cell_names"]]
  )

  old_attributes <- attributes(coefficient_samples)
  old_class <- class(coefficient_samples)
  old_attributes <- old_attributes[
    !names(old_attributes) %in% c(
      "dim", "dimnames", "names", "class", "level_names",
      "posterior_support"
    )
  ]
  attributes(transformed_samples) <- c(attributes(transformed_samples), old_attributes)
  attr(transformed_samples, "level_names")       <- design_info[["cell_names"]]
  attr(transformed_samples, "factor_cell_names") <- design_info[["cell_names"]]
  class(transformed_samples) <- unique(c(old_class, class(transformed_samples), transformed_class))

  return(transformed_samples)
}

#' @title Transform factor posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on meandif/orthonormal prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_meandif_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [transform_meandif_samples] [transform_meandif_samples] [transform_orthonormal_samples]
#'
#' @export
transform_factor_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  samples <- transform_meandif_samples(samples)
  samples <- transform_orthonormal_samples(samples)

  return(samples)
}

#' @title Transform meandif posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on meandif prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_meandif_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.meandif]
#'
#' @export
transform_meandif_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.meandif_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "meandif"))){

      meandif_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = meandif_samples,
        metadata            = meandif_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.meandif_transformed"
      )
    }
  }

  return(samples)
}

#' @title Transform orthonomal posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on orthonormal prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#' @return \code{transform_orthonormal_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.orthonormal]
#'
#' @export
transform_orthonormal_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.orthonormal_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "orthonormal"))){

      orthonormal_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = orthonormal_samples,
        metadata            = orthonormal_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.orthonormal_transformed"
      )
    }
  }

  return(samples)
}

# not part of transform factor samples (as it's usefull only for marginal effects)
transform_treatment_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.treatment_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && isTRUE(attr(samples[[i]], "treatment"))){

      treatment_samples <- .add_factor_metadata_from_named_objects(samples[[i]], names(samples)[i], samples)
      samples[[i]] <- .transform_factor_contrast_samples(
        coefficient_samples = treatment_samples,
        metadata            = treatment_samples,
        parameter           = names(samples)[i],
        transformed_class   = "mixed_posteriors.treatment_transformed"
      )
    }
  }

  return(samples)
}


# Helper: Parse a term name into its component variable names
# e.g., "mu_x1__xXx__x2" with prefix "mu" -> c("x1", "x2")
# e.g., "mu_intercept" -> character(0) (intercept has no components)
# e.g., "mu_x1" -> c("x1")
.parse_term_components <- function(term_name, prefix) {
  # Remove prefix
  term_part <- sub(paste0("^", prefix, "_"), "", term_name)

  # Check if it's the intercept
  if (term_part == "intercept") {
    return(character(0))
  }

  # Indexed factor terms attach the index to the full term name. Strip the
  # trailing index so scaled continuous components inside interactions are
  # still detected correctly (e.g., alloc__xXx__year[1] -> alloc, year).
  term_part <- sub("\\[[^]]+\\]$", "", term_part)

  # Split by interaction separator
  components <- strsplit(term_part, "__xXx__", fixed = TRUE)[[1]]
  return(components)
}


# Helper: Parse a term into scaled components and an unscaled identity string
# used to determine which coefficients can contribute to each other during
# de-standardization.
.parse_unscale_term_structure <- function(term_name, prefix, scaled_vars) {

  term_part <- sub(paste0("^", prefix, "_"), "", term_name)

  if (term_part == "intercept") {
    return(list(
      components  = character(0),
      scaled      = character(0),
      unscaled_id = ""
    ))
  }

  term_index <- ""
  if (grepl("\\[[^]]+\\]$", term_part)) {
    term_index <- sub("^.*(\\[[^]]+\\])$", "\\1", term_part)
    term_core  <- sub("\\[[^]]+\\]$", "", term_part)
  } else {
    term_core <- term_part
  }

  components <- strsplit(term_core, "__xXx__", fixed = TRUE)[[1]]
  is_scaled  <- components %in% scaled_vars

  unscaled_components <- components[!is_scaled]
  unscaled_id <- paste(unscaled_components, collapse = "__xXx__")

  if (nzchar(term_index) && length(unscaled_components) > 0) {
    unscaled_id <- paste0(unscaled_id, term_index)
  }

  return(list(
    components  = components,
    scaled      = components[is_scaled],
    unscaled_id = unscaled_id
  ))
}


# Helper: Check if set A is a subset of set B (including equality)
.is_subset <- function(A, B) {

  length(A) == 0 || all(A %in% B)
}


# Helper: Compare lower-order unscaled term identities. Two-level factor main
# effects are stored without [1], while their interactions are indexed.
.unscale_ids_match <- function(target_id, source_id) {

  identical(target_id, source_id) ||
    (nzchar(target_id) && !grepl("\\[[^]]+\\]$", target_id) && identical(paste0(target_id, "[1]"), source_id))
}


# Helper: Validate the nested formula_scale structure used for unscaling.
.check_formula_scale_info <- function(formula_scale, name = "formula_scale") {

  check_list(formula_scale, name)

  if(is.null(names(formula_scale)) || anyNA(names(formula_scale)) || any(names(formula_scale) == ""))
    stop(paste0("The '", name, "' argument must be a named nested list keyed by parameter name."), call. = FALSE)

  for(param_name in names(formula_scale)){
    param_scale <- formula_scale[[param_name]]
    check_list(param_scale, paste0(name, "[['", param_name, "]]"))

    if(length(param_scale) == 0)
      next

    if(is.null(names(param_scale)) || anyNA(names(param_scale)) || any(names(param_scale) == ""))
      stop(paste0("The '", name, "[['", param_name, "]]" ,"' entry must be a named list keyed by parameter term."), call. = FALSE)

    for(term_name in names(param_scale)){
      term_scale <- param_scale[[term_name]]

      check_list(
        term_scale,
        paste0(name, "[['", param_name, "]][['", term_name, "]]"),
        check_names = c("mean", "sd"),
        all_objects = TRUE,
        allow_other = TRUE
      )
      check_real(
        term_scale[["mean"]],
        paste0(name, "[['", param_name, "]][['", term_name, "]][['mean'] ]"),
        allow_NA = FALSE
      )
      check_real(
        term_scale[["sd"]],
        paste0(name, "[['", param_name, "]][['", term_name, "]][['sd'] ]"),
        lower = 0,
        allow_bound = FALSE,
        allow_NA = FALSE
      )
    }
  }

  invisible(NULL)
}


# Helper: warn when formula_scale entries do not map to posterior terms.
.warn_unused_formula_scale_terms <- function(term_names, formula_scale, prefix) {

  scaled_terms <- names(formula_scale)
  if(length(scaled_terms) == 0)
    return(invisible(NULL))

  scaled_vars <- sub(paste0("^", prefix, "_"), "", scaled_terms)
  term_components <- unique(unlist(lapply(term_names, .parse_term_components, prefix = prefix), use.names = FALSE))
  random_scaled_vars <- .formula_scale_random_scaled_vars(
    formula_scale,
    scaled_vars = scaled_vars
  )
  unused_terms <- scaled_terms[
    !scaled_vars %in% term_components &
      !scaled_vars %in% random_scaled_vars
  ]

  if(length(unused_terms) == 0)
    return(invisible(NULL))

  if(length(unused_terms) == length(scaled_terms)){
    warning(
      paste0(
        "Ignoring all formula_scale[['", prefix, "]] entries because none match posterior terms: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Samples for this parameter prefix are returned unchanged."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }else{
    warning(
      paste0(
        "Ignoring unused formula_scale[['", prefix, "]] entries: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Matched entries are still applied."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  invisible(NULL)
}

.formula_scale_random_scaled_vars <- function(formula_scale, scaled_vars = NULL){

  if(is.null(scaled_vars)){
    prefix <- attr(formula_scale, "parameter")
    if(!is.null(prefix) && length(prefix) == 1L && !is.na(prefix) && nzchar(prefix)){
      scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))
    }else{
      scaled_vars <- names(formula_scale)
    }
  }

  random_terms <- character()

  metadata <- attr(formula_scale, "random_effect_terms")
  if(!is.null(metadata)){
    metadata <- unname(as.character(metadata))
    metadata <- metadata[is.na(metadata) | metadata != "sd"]
    random_terms <- c(random_terms, metadata)
  }

  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves")
  if(!is.null(sd_leaves) && length(sd_leaves) > 0L){
    for(leaves in sd_leaves){
      random_terms <- c(
        random_terms,
        .formula_scale_random_leaf_terms(
          leaves = leaves,
          scaled_vars = scaled_vars
        )
      )
    }
  }

  random_terms <- random_terms[!is.na(random_terms) & nzchar(random_terms)]
  random_terms <- setdiff(random_terms, "intercept")
  if(length(random_terms) == 0L){
    return(character())
  }

  random_components <- unique(unlist(lapply(random_terms, .bt_random_effect_term_components), use.names = FALSE))
  random_components <- setdiff(random_components, "intercept")
  if(length(scaled_vars) > 0L){
    random_components <- random_components[random_components %in% scaled_vars]
  }
  unique(random_components)
}

.formula_scale_random_leaf_terms <- function(leaves, scaled_vars){

  if(.random_sd_structured_leaves(leaves)){
    return(character())
  }

  leaf_terms <- character()
  if(!is.null(leaves$leaf_terms)){
    leaf_terms <- unname(as.character(leaves$leaf_terms))
  }

  homogeneous_sd <- length(leaf_terms) > 0L &&
    all(!is.na(leaf_terms) & leaf_terms == "sd")
  if(!homogeneous_sd){
    if(!is.null(leaves$leaf_terms_by_column)){
      leaf_terms <- c(leaf_terms, unname(as.character(leaves$leaf_terms_by_column)))
    }
    return(leaf_terms)
  }

  if(is.null(leaves$column_names)){
    return(character())
  }

  column_terms <- vapply(
    leaves$column_names,
    .random_sd_term_from_column_name,
    character(1)
  )
  scaled_columns <- vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  )

  column_terms[scaled_columns]
}


# Helper: Build the transformation matrix for unscaling coefficients
#
# For each target term T and source term S, computes the coefficient M[T,S] such that:
#   coef_orig[T] = sum over S of M[T,S] * coef_z[S]
#
# The formula is based on expanding products of (x_i - mu_i)/sigma_i terms.
# For S to contribute to T:
#   1. T_unscaled == S_unscaled (unscaled components must match exactly)
#   2. T_scaled is a subset of S_scaled
#
# The contribution is: (-1)^|extra| * prod(mu_extra) / prod(sigma_S_scaled)
# where extra = S_scaled \ T_scaled
#
# @param term_names Character vector of all term names in the posterior
# @param formula_scale Named list with scaling info (mean, sd) for scaled predictors
# @param prefix The parameter prefix (e.g., "mu")
# @return A square transformation matrix
.build_unscale_matrix <- function(term_names, formula_scale, prefix) {

  n_terms <- length(term_names)
  M <- diag(n_terms)  # Start with identity matrix
  rownames(M) <- colnames(M) <- term_names

  # Extract the variable names that are scaled (without prefix)
  scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))

  # Parse all terms into their scaled components and unscaled identity.
  term_structure <- lapply(
    term_names,
    .parse_unscale_term_structure,
    prefix = prefix,
    scaled_vars = scaled_vars
  )
  names(term_structure) <- term_names

  term_components <- lapply(term_structure, `[[`, "components")
  term_scaled <- lapply(term_structure, `[[`, "scaled")
  term_unscaled <- vapply(term_structure, `[[`, character(1), "unscaled_id")

  # Warn about high-order interactions
  max_order <- max(sapply(term_components, length))
  if (max_order >= 5) {
    warning("Model contains ", max_order, "-way or higher interactions. ",
            "Unscaling transformation may be computationally intensive.",
            immediate. = TRUE)
  }

  # Build the transformation matrix
  for (t_idx in seq_along(term_names)) {
    T_name <- term_names[t_idx]
    T_scaled <- term_scaled[[T_name]]
    T_unscaled <- term_unscaled[[T_name]]

    for (s_idx in seq_along(term_names)) {
      S_name <- term_names[s_idx]
      S_scaled <- term_scaled[[S_name]]
      S_unscaled <- term_unscaled[[S_name]]

      # Check contribution conditions
      # 1. Unscaled parts must match exactly
      if (!.unscale_ids_match(T_unscaled, S_unscaled)) next

      # 2. T_scaled must be a subset of S_scaled
      if (!.is_subset(T_scaled, S_scaled)) next

      # 3. S must have at least one scaled component (otherwise no transformation needed)
      if (length(S_scaled) == 0) {
        # No scaling for this source term - keep identity (already set)
        next
      }

      # Compute the coefficient
      extra_scaled <- setdiff(S_scaled, T_scaled)

      # Sign: (-1)^|extra|
      sign <- (-1)^length(extra_scaled)

      # Product of means for extra scaled components
      if (length(extra_scaled) > 0) {
        extra_params <- paste0(prefix, "_", extra_scaled)
        mean_product <- prod(sapply(extra_params, function(p) formula_scale[[p]]$mean))
      } else {
        mean_product <- 1
      }

      # Product of SDs for all scaled components in S
      S_scaled_params <- paste0(prefix, "_", S_scaled)
      sd_product <- prod(sapply(S_scaled_params, function(p) formula_scale[[p]]$sd))

      # Contribution coefficient
      M[t_idx, s_idx] <- sign * mean_product / sd_product
    }
  }

  return(M)
}


# Helper: Apply unscaling transformation to a matrix of posterior samples
#
# @param posterior Matrix with samples in rows, parameters in columns
# Apply the unscaling transformation to posterior samples
#
# @param posterior Matrix of posterior samples with parameter names as column names
# @param formula_scale Nested list with scaling info keyed by parameter name:
#   list(mu = list(mu_x1 = list(mean, sd)), log_sigma = list(log_sigma_x = list(mean, sd)))
# @return Transformed posterior matrix
.apply_unscale_transform <- function(posterior, formula_scale) {

  if (is.null(formula_scale) || length(formula_scale) == 0) {
    return(posterior)
  }

  # Handle nested structure: iterate over each parameter
  matched_prefix <- FALSE
  for (param_name in names(formula_scale)) {
    param_scale <- formula_scale[[param_name]]
    affected_cols <- grep(paste0("^", param_name, "_"), colnames(posterior), value = TRUE)
    if(length(affected_cols) == 0)
      next

    matched_prefix <- TRUE
    posterior <- .apply_unscale_transform_single(posterior, param_scale, prefix = param_name)
  }

  if(!matched_prefix){
    warning(
      "Ignoring formula_scale because none of its parameter prefixes match the posterior columns.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  return(posterior)
}

# Helper: Apply unscaling for a single parameter's predictors
# @param posterior Matrix of posterior samples
# @param formula_scale Flat list of scaling info: list(mu_x1 = list(mean, sd), mu_x2 = list(mean, sd))
# @param prefix Parameter prefix (e.g., "mu")
# @return Transformed posterior matrix
.apply_unscale_transform_single <- function(posterior, formula_scale, prefix) {

  if (is.null(formula_scale) || length(formula_scale) == 0) {
    return(posterior)
  }

  # Check if this parameter uses log(intercept)
  log_intercept <- isTRUE(attr(formula_scale, "log_intercept"))
  intercept_col <- paste0(prefix, "_intercept")

  # Identify which columns are affected by the transformation
  affected_cols <- grep(paste0("^", prefix, "_"), colnames(posterior), value = TRUE)
  if (length(affected_cols) > 0) {
    posterior <- .materialize_formula_scale_point_terms(
      posterior = posterior,
      formula_scale = formula_scale,
      prefix = prefix
    )
    affected_cols <- grep(paste0("^", prefix, "_"), colnames(posterior), value = TRUE)
  }
  random_sd_cols  <- grep(paste0("^", prefix, "__xREx__"), affected_cols, value = TRUE)
  random_aux_cols <- grep(paste0("^", prefix, "__(xRE_ALLOCx|xRE_SUMMARY__)"), affected_cols, value = TRUE)
  fixed_cols      <- setdiff(affected_cols, c(random_sd_cols, random_aux_cols))

  if (length(affected_cols) == 0) {
    return(posterior)
  }

  if(length(fixed_cols) > 0){
    .warn_unused_formula_scale_terms(fixed_cols, formula_scale, prefix)
  }

  # For log(intercept): transform to log scale before unscaling, then exp() back
  # This works because: log_sigma = log(intercept) + beta * x_z
  # is equivalent to: log_sigma = log_int + beta * x_z (standard additive form)
  # where log_int = log(intercept)
  if (length(fixed_cols) > 0 && log_intercept && intercept_col %in% colnames(posterior)) {
    posterior[, intercept_col] <- log(posterior[, intercept_col])
  }

  # Build and apply standard transformation matrix
  if(length(fixed_cols) > 0){
    M <- .build_unscale_matrix(fixed_cols, formula_scale, prefix)
    posterior[, fixed_cols] <- posterior[, fixed_cols, drop = FALSE] %*% t(M)
  }

  # Transform intercept back from log scale
  if (length(fixed_cols) > 0 && log_intercept && intercept_col %in% colnames(posterior)) {
    posterior[, intercept_col] <- exp(posterior[, intercept_col])
  }

  posterior <- .apply_random_sd_unscale(posterior, random_sd_cols, formula_scale, prefix)

  return(posterior)
}

.materialize_formula_scale_point_terms <- function(posterior, formula_scale,
                                                   prefix){

  point_terms <- attr(formula_scale, "point_terms", exact = TRUE)
  if(is.null(point_terms) || length(point_terms) == 0L){
    return(posterior)
  }

  point_terms <- point_terms[startsWith(names(point_terms), paste0(prefix, "_"))]
  point_terms <- point_terms[!names(point_terms) %in% colnames(posterior)]
  if(length(point_terms) == 0L){
    return(posterior)
  }

  point_matrix <- matrix(
    rep(unname(point_terms), each = nrow(posterior)),
    nrow = nrow(posterior),
    dimnames = list(NULL, names(point_terms))
  )

  cbind(posterior, point_matrix)
}

.apply_random_sd_unscale <- function(posterior, random_sd_cols, formula_scale,
                                     prefix,
                                     correlation_required_groups = NULL){

  if(length(random_sd_cols) == 0){
    return(posterior)
  }

  column_groups <- .random_sd_column_unscale_groups(
    random_sd_cols = random_sd_cols,
    formula_scale = formula_scale,
    prefix = prefix
  )
  if(!is.null(column_groups)){
    return(.apply_random_sd_column_unscale(
      posterior = posterior,
      column_groups = column_groups,
      formula_scale = formula_scale,
      prefix = prefix,
      correlation_required_groups = correlation_required_groups
    ))
  }

  term_map <- .random_sd_term_map(random_sd_cols, formula_scale, prefix)
  if(length(term_map) == 0){
    return(posterior)
  }

  random_sd_cols <- names(term_map)
  group_keys <- vapply(random_sd_cols, .random_sd_group_key,
                       character(1), term_map = term_map, prefix = prefix)
  required_groups <- .random_sd_correlation_required_groups(
    formula_scale = formula_scale,
    prefix = prefix,
    correlation_required_groups = correlation_required_groups
  )

  for(group_key in unique(group_keys)){
    group_cols <- random_sd_cols[group_keys == group_key]
    group_terms <- unname(term_map[group_cols])

    if(any(duplicated(group_terms))){
      duplicated_terms <- unique(group_terms[duplicated(group_terms)])
      stop(
        "Random-effect SD columns for group '", group_key,
        "' cannot be unscaled because multiple columns map to the same term: ",
        paste0("'", duplicated_terms, "'", collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    pseudo_terms <- paste0(prefix, "_", group_terms)
    names(pseudo_terms) <- group_cols
    M <- .build_unscale_matrix(unname(pseudo_terms), formula_scale, prefix)

    source_sd <- posterior[, group_cols, drop = FALSE]
    source_cor <- .random_sd_correlation_draws(
      posterior = posterior,
      prefix = prefix,
      group_key = group_key,
      n_terms = length(group_cols),
      required = group_key %in% required_groups
    )
    transformed_sd <- matrix(NA_real_, nrow = nrow(source_sd), ncol = ncol(source_sd))
    colnames(transformed_sd) <- group_cols

    if(is.null(source_cor)){
      for(target_i in seq_along(group_cols)){
        transformed_var <- rowSums(t(t(source_sd^2) * (M[target_i, ]^2)))
        transformed_sd[, target_i] <- sqrt(transformed_var)
      }
    }else{
      transformed_cor <- array(NA_real_, dim = dim(source_cor))
      valid_cor_draw <- rep(FALSE, nrow(source_sd))
      for(draw_i in seq_len(nrow(source_sd))){
        source_cov <- diag(source_sd[draw_i, ], nrow = ncol(source_sd)) %*%
          source_cor[draw_i, , ] %*%
          diag(source_sd[draw_i, ], nrow = ncol(source_sd))
        transformed_cov <- M %*% source_cov %*% t(M)
        transformed_sd[draw_i, ] <- sqrt(diag(transformed_cov))
        if(any(!is.finite(transformed_sd[draw_i, ]) | transformed_sd[draw_i, ] <= 0)){
          next
        }
        transformed_cor[draw_i, , ] <- transformed_cov /
          tcrossprod(transformed_sd[draw_i, ])
        diag(transformed_cor[draw_i, , ]) <- 1
        valid_cor_draw[draw_i] <- all(is.finite(transformed_cor[draw_i, , ]))
      }
      posterior <- .random_sd_assign_transformed_correlation(
        posterior = posterior,
        prefix = prefix,
        group_key = group_key,
        correlation = transformed_cor,
        valid_draw = valid_cor_draw
      )
    }

    posterior[, group_cols] <- transformed_sd
  }

  posterior
}

.apply_random_sd_column_unscale <- function(posterior, column_groups,
                                            formula_scale, prefix,
                                            correlation_required_groups = NULL){

  if(length(column_groups) == 0L){
    return(posterior)
  }

  required_groups <- .random_sd_correlation_required_groups(
    formula_scale = formula_scale,
    prefix = prefix,
    correlation_required_groups = correlation_required_groups
  )

  for(group in column_groups){
    group_key <- group$group_key
    pseudo_terms <- paste0(prefix, "_", group$column_terms)
    M <- .build_unscale_matrix(pseudo_terms, formula_scale, prefix)

    source_sd <- posterior[, group$leaf_names_by_column, drop = FALSE]
    source_cor <- .random_sd_correlation_draws(
      posterior = posterior,
      prefix = prefix,
      group_key = group_key,
      n_terms = length(group$column_terms),
      required = group_key %in% required_groups
    )
    transformed_sd_by_column <- matrix(
      NA_real_,
      nrow = nrow(source_sd),
      ncol = ncol(source_sd),
      dimnames = list(NULL, group$leaf_names_by_column)
    )

    if(is.null(source_cor)){
      for(target_i in seq_along(group$column_terms)){
        transformed_var <- rowSums(t(t(source_sd^2) * (M[target_i, ]^2)))
        transformed_sd_by_column[, target_i] <- sqrt(transformed_var)
      }
    }else{
      transformed_cor <- array(NA_real_, dim = dim(source_cor))
      valid_cor_draw <- rep(FALSE, nrow(source_sd))
      for(draw_i in seq_len(nrow(source_sd))){
        source_cov <- diag(source_sd[draw_i, ], nrow = ncol(source_sd)) %*%
          source_cor[draw_i, , ] %*%
          diag(source_sd[draw_i, ], nrow = ncol(source_sd))
        transformed_cov <- M %*% source_cov %*% t(M)
        transformed_sd_by_column[draw_i, ] <- sqrt(diag(transformed_cov))
        if(any(!is.finite(transformed_sd_by_column[draw_i, ]) |
               transformed_sd_by_column[draw_i, ] <= 0)){
          next
        }
        transformed_cor[draw_i, , ] <- transformed_cov /
          tcrossprod(transformed_sd_by_column[draw_i, ])
        diag(transformed_cor[draw_i, , ]) <- 1
        valid_cor_draw[draw_i] <- all(is.finite(transformed_cor[draw_i, , ]))
      }
      posterior <- .random_sd_assign_transformed_correlation(
        posterior = posterior,
        prefix = prefix,
        group_key = group_key,
        correlation = transformed_cor,
        valid_draw = valid_cor_draw
      )
    }

    for(sd_col in group$leaf_names){
      column_index <- which(group$leaf_names_by_column == sd_col)
      if(length(column_index) == 1L){
        posterior[, sd_col] <- transformed_sd_by_column[, column_index]
        next
      }
      shared_values <- transformed_sd_by_column[, column_index, drop = FALSE]
      max_difference <- max(
        abs(shared_values - shared_values[, 1L]),
        na.rm = TRUE
      )
      if(is.finite(max_difference) && max_difference > sqrt(.Machine$double.eps)){
        stop(
          "Random-effect SD column '", sd_col,
          "' cannot be unscaled because its shared design columns transform to different SDs.",
          call. = FALSE
        )
      }
      posterior[, sd_col] <- shared_values[, 1L]
    }
  }

  posterior
}

.random_sd_assign_transformed_correlation <- function(posterior, prefix,
                                                      group_key,
                                                      correlation,
                                                      valid_draw = NULL){

  n_terms <- dim(correlation)[2L]
  if(is.null(valid_draw)){
    valid_draw <- rep(TRUE, dim(correlation)[1L])
  }
  valid_draw <- valid_draw & is.finite(vapply(
    seq_len(dim(correlation)[1L]),
    function(draw_i) sum(correlation[draw_i, , ]),
    numeric(1)
  ))

  R_names <- .random_sd_correlation_matrix_names(
    prefix = prefix,
    group_key = group_key,
    suffix = "_xRE_CORx_R",
    n_terms = n_terms
  )
  L_names <- .random_sd_correlation_matrix_names(
    prefix = prefix,
    group_key = group_key,
    suffix = "_xRE_CORx_L",
    n_terms = n_terms
  )

  has_R <- all(as.vector(R_names) %in% colnames(posterior))
  has_L <- all(as.vector(L_names) %in% colnames(posterior))
  if(!has_R && !has_L){
    return(posterior)
  }

  if(has_R){
    posterior[, as.vector(R_names)] <- NA_real_
  }
  if(has_L){
    posterior[, as.vector(L_names)] <- NA_real_
  }
  if(!any(valid_draw)){
    return(posterior)
  }

  L <- NULL
  if(has_L){
    L <- array(NA_real_, dim = dim(correlation))
    for(draw_i in which(valid_draw)){
      this_L <- try(t(chol(correlation[draw_i, , ])), silent = TRUE)
      if(inherits(this_L, "try-error")){
        valid_draw[draw_i] <- FALSE
        next
      }
      L[draw_i, , ] <- this_L
    }
    if(!any(valid_draw)){
      return(posterior)
    }
  }

  if(has_R){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[valid_draw, R_names[row, column]] <-
          correlation[valid_draw, row, column]
      }
    }
  }
  if(has_L){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[valid_draw, L_names[row, column]] <-
          L[valid_draw, row, column]
      }
    }
  }

  posterior
}

.random_sd_correlation_matrix_names <- function(prefix, group_key, suffix,
                                                n_terms){

  stem <- paste0(prefix, "__xREx__", group_key, suffix)
  outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "[", row, ",", column, "]")
    })
  )
}

.random_sd_term_map <- function(random_sd_cols, formula_scale, prefix){

  scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))
  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves")
  if(!is.null(sd_leaves) && length(sd_leaves) > 0){
    descriptor_terms <- do.call(
      c,
      unname(lapply(sd_leaves, .random_sd_leaf_term_map, scaled_vars = scaled_vars))
    )
    term_map <- descriptor_terms[random_sd_cols]
    names(term_map) <- random_sd_cols
    term_map <- term_map[!is.na(term_map)]
    if(length(term_map) > 0){
      return(term_map)
    }
  }

  metadata <- attr(formula_scale, "random_effect_terms")
  if(!is.null(metadata) && length(metadata) > 0){
    base_cols <- sub("\\[[^]]+\\]$", "", random_sd_cols)
    term_map <- metadata[base_cols]
    names(term_map) <- random_sd_cols
    term_map <- term_map[!is.na(term_map)]
    if(length(term_map) > 0){
      return(term_map)
    }
  }

  possible_terms <- c("intercept", scaled_vars)
  names(possible_terms) <- possible_terms

  term_map <- vapply(random_sd_cols, function(col){
    rest <- sub(paste0("^", prefix, "__xREx__"), "", sub("\\[[^]]+\\]$", "", col))
    candidates <- possible_terms[vapply(possible_terms, function(term){
      endsWith(rest, paste0("_", term))
    }, logical(1))]
    if(length(candidates) == 0){
      return(NA_character_)
    }
    candidates[which.max(nchar(candidates))]
  }, character(1))

  term_map[!is.na(term_map)]
}

.random_sd_column_unscale_groups <- function(random_sd_cols, formula_scale,
                                             prefix){

  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves", exact = TRUE)
  if(is.null(sd_leaves)){
    return(NULL)
  }

  scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))
  groups <- list()
  leaf_keys <- names(sd_leaves)
  if(is.null(leaf_keys)){
    leaf_keys <- rep("", length(sd_leaves))
  }
  for(i in seq_along(sd_leaves)){
    group <- .random_sd_column_unscale_group(
      leaves = sd_leaves[[i]],
      leaf_key = leaf_keys[[i]],
      random_sd_cols = random_sd_cols,
      scaled_vars = scaled_vars,
      prefix = prefix
    )
    if(is.null(group)){
      next
    }
    groups[[length(groups) + 1L]] <- group
  }

  groups
}

.random_sd_column_unscale_group <- function(leaves, leaf_key, random_sd_cols,
                                            scaled_vars, prefix){

  if(.random_sd_structured_leaves(leaves)){
    return(NULL)
  }

  if(is.null(leaves$leaf_names_by_column) ||
     is.null(leaves$leaf_terms_by_column) ||
     is.null(leaves$leaf_names)){
    return(NULL)
  }

  leaf_names_by_column <- as.character(leaves$leaf_names_by_column)
  column_terms <- as.character(leaves$leaf_terms_by_column)
  leaf_names <- as.character(leaves$leaf_names)

  if(length(leaf_names_by_column) == 0L ||
     length(column_terms) != length(leaf_names_by_column)){
    return(NULL)
  }

  available_leaf_names <- leaf_names[leaf_names %in% random_sd_cols]
  if(length(available_leaf_names) == 0L){
    return(NULL)
  }

  homogeneous <- !is.null(leaves$leaf_terms) &&
    identical(unique(unname(as.character(leaves$leaf_terms))), "sd")
  if(isTRUE(homogeneous) && !is.null(leaves$column_names)){
    column_terms <- vapply(
      leaves$column_names,
      .random_sd_term_from_column_name,
      character(1)
    )
  }
  scaled_columns <- vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  )

  if(!any(scaled_columns)){
    return(NULL)
  }

  if(isTRUE(homogeneous)){
    if(length(column_terms) == 1L){
      column_terms <- .random_sd_term_from_column_name(leaves$column_names[[1L]])
    }else{
      stop(
        "Cannot unscale homogeneous random-effect SD '",
        paste(available_leaf_names, collapse = "', '"),
        "' because its block contains scaled random-slope columns. ",
        "Use a heterogeneous random-effect SD structure or leave samples on the fitted scale.",
        call. = FALSE
      )
    }
  }

  missing_leaf_names <- setdiff(leaf_names, random_sd_cols)
  group_key <- .random_sd_group_key_from_leaf_key(
    leaf_key = leaf_key,
    leaf_names = leaf_names,
    column_terms = column_terms,
    prefix = prefix
  )
  if(length(missing_leaf_names) > 0L){
    stop(
      "Random-effect SD unscaling for block '", group_key,
      "' requires all SD columns from the block. Missing: ",
      paste0("'", missing_leaf_names, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    group_key = group_key,
    leaf_names = leaf_names,
    leaf_names_by_column = leaf_names_by_column,
    column_terms = column_terms
  )
}

.random_sd_structured_leaves <- function(leaves){

  random_structure <- leaves$random_structure
  !is.null(random_structure) &&
    length(random_structure) == 1L &&
    random_structure %in% c("cs", "hcs", "ar1", "car", "har")
}

.random_sd_group_key_from_leaf_key <- function(leaf_key, leaf_names,
                                               column_terms, prefix){

  if(!is.null(leaf_key) && length(leaf_key) == 1L &&
     !is.na(leaf_key) && nzchar(leaf_key)){
    return(sub("^__xREx__", "", leaf_key))
  }

  term_map <- stats::setNames(column_terms[seq_along(leaf_names)], leaf_names)
  .random_sd_group_key(
    col = leaf_names[[1L]],
    term_map = term_map,
    prefix = prefix
  )
}

.random_sd_leaf_term_map <- function(leaves, scaled_vars){

  if(is.null(leaves$leaf_terms)){
    return(stats::setNames(character(), character()))
  }

  if(.random_sd_structured_leaves(leaves)){
    return(stats::setNames(character(), character()))
  }

  out <- leaves$leaf_terms
  homogeneous <- identical(unique(unname(leaves$leaf_terms)), "sd")
  if(!isTRUE(homogeneous)){
    return(out)
  }

  if(is.null(leaves$leaf_names_by_column) || is.null(leaves$column_names)){
    return(out)
  }

  leaf_names <- unique(leaves$leaf_names_by_column)
  if(length(leaf_names) != 1L){
    return(out)
  }

  column_terms <- vapply(
    leaves$column_names,
    .random_sd_term_from_column_name,
    character(1)
  )
  scaled_columns <- vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  )

  if(length(column_terms) == 1L){
    if(isTRUE(scaled_columns)){
      out[[leaf_names]] <- column_terms
    }
    return(out)
  }

  if(any(scaled_columns)){
    stop(
      "Cannot unscale homogeneous random-effect SD '", leaf_names,
      "' because its block contains scaled random-slope columns. ",
      "Use a heterogeneous random-effect SD structure or leave samples on the fitted scale.",
      call. = FALSE
    )
  }

  out
}

.random_sd_term_from_column_name <- function(column_name){

  if(identical(column_name, "(Intercept)")){
    return("intercept")
  }

  column_name
}

.random_sd_term_uses_scaled_var <- function(term, scaled_vars){

  if(length(scaled_vars) == 0L || identical(term, "intercept")){
    return(FALSE)
  }

  term <- sub("\\[[^]]+\\]$", "", term)
  components <- .bt_random_effect_term_components(term)
  any(components %in% scaled_vars)
}

.random_sd_group_key <- function(col, term_map, prefix){

  base_col <- sub("\\[[^]]+\\]$", "", col)
  rest <- sub(paste0("^", prefix, "__xREx__"), "", base_col)
  term <- unname(term_map[[col]])
  term_core <- gsub("\\[[^]]+\\]", "", term)
  .random_sd_strip_term_suffix(rest, term_core)
}

.random_sd_strip_term_suffix <- function(x, term){

  suffix <- paste0("_", term)
  if(endsWith(x, suffix)){
    return(substr(x, 1L, nchar(x) - nchar(suffix)))
  }

  x
}

.random_sd_correlation_required_groups <- function(formula_scale, prefix,
                                                   correlation_required_groups = NULL){

  required <- character()
  metadata <- attr(formula_scale, "random_effect_correlation_required")
  if(!is.null(metadata)){
    if(is.logical(metadata)){
      if(!is.null(names(metadata))){
        metadata <- names(metadata)[metadata]
      }else{
        metadata <- character()
      }
    }
    required <- c(required, as.character(metadata))
  }
  if(!is.null(correlation_required_groups)){
    required <- c(required, as.character(correlation_required_groups))
  }
  required <- required[!is.na(required) & nzchar(required)]
  if(length(required) == 0L){
    return(character())
  }

  required <- sub(paste0("^", prefix, "__xREx__"), "", required)
  required <- sub("^__xREx__", "", required)
  unique(required)
}

.random_sd_correlation_draws <- function(posterior, prefix, group_key, n_terms,
                                         required = FALSE){

  if(n_terms < 2L){
    return(NULL)
  }

  stem <- paste0(prefix, "__xREx__", group_key)
  R_names <- outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  R_present <- as.vector(R_names) %in% colnames(posterior)

  L_names <- outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "_xRE_CORx_L[", row, ",", column, "]")
    })
  )
  L_present <- as.vector(L_names) %in% colnames(posterior)

  if(any(R_present) && !all(R_present)){
    stop(
      "Random-effect correlation samples are incomplete for block '",
      group_key,
      "'. Expected a complete random-effect correlation matrix.",
      call. = FALSE
    )
  }
  if(any(L_present) && !all(L_present)){
    stop(
      "Random-effect Cholesky samples are incomplete for block '",
      group_key,
      "'. Expected a complete random-effect Cholesky matrix.",
      call. = FALSE
    )
  }

  if(all(R_present)){
    out <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        out[, row, column] <- posterior[, R_names[row, column]]
      }
    }
    return(out)
  }

  if(!all(L_present)){
    if(isTRUE(required)){
      stop(
        "Random-effect SD unscaling requires random-effect correlation samples for block '",
        group_key,
        "'. Expected monitored correlation or Cholesky coordinates.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  L <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
  for(row in seq_len(n_terms)){
    for(column in seq_len(n_terms)){
      L[, row, column] <- posterior[, L_names[row, column]]
    }
  }

  out <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
  for(draw_i in seq_len(nrow(posterior))){
    out[draw_i, , ] <- L[draw_i, , ] %*% t(L[draw_i, , ])
  }

  out
}


#' @title Transform standardized posterior samples back to original scale
#'
#' @description Transforms posterior samples from standardized continuous
#' predictors back to the original scale. This function is used when predictors
#' were standardized during model fitting via the \code{formula_scale} parameter.
#'
#' @param fit a fitted model object with \code{formula_scale} attribute, or
#' a matrix of posterior samples
#' @param formula_scale nested list containing standardization information keyed by
#' parameter name. Each parameter entry contains scaling info (mean and sd) for
#' each standardized predictor, e.g., \code{list(mu = list(mu_x1 = list(mean = 0, sd = 1)))}.
#' If \code{fit} is provided and has a \code{formula_scale} attribute, this will be used automatically.
#'
#' @details The function transforms regression coefficients and intercepts
#' to account for predictor standardization using a combinatorial approach that
#' correctly handles interactions of any order.
#'
#' For a k-way interaction between standardized predictors, the expansion of
#' \eqn{\prod_{i} (x_i - \mu_i)/\sigma_i} contributes to all lower-order terms.
#' The contribution to a target term T from a source term S (where T is a subset
#' of S's scaled components) is:
#' \deqn{(-1)^{|extra|} \cdot \prod_{i \in extra} \mu_i / \prod_{i \in S_{scaled}} \sigma_i}
#' where \eqn{extra = S_{scaled} \setminus T_{scaled}}.
#'
#' @return \code{transform_scale_samples} returns posterior samples transformed
#' back to the original predictor scale.
#'
#' @seealso [JAGS_formula()] [JAGS_fit()]
#'
#' @export
transform_scale_samples <- function(fit, formula_scale = NULL){

  # extract formula_scale from fit if available
  if(is.null(formula_scale) && !is.null(attr(fit, "formula_scale"))){
    formula_scale <- attr(fit, "formula_scale")
  }

  if(is.null(formula_scale) || length(formula_scale) == 0){
    # no scaling information, return as is
    return(fit)
  }

  .check_formula_scale_info(formula_scale)

  # extract posterior samples
  if(inherits(fit, "runjags") || inherits(fit, "BayesTools_fit")){
    posterior <- as.matrix(.fit_to_posterior(fit))
  }else if(is.matrix(fit)){
    posterior <- fit
  }else{
    stop("'fit' must be a fitted model object or a matrix of posterior samples.")
  }

  # Apply the combinatorial unscaling transformation
  posterior <- .apply_unscale_transform(posterior, formula_scale)

  return(posterior)
}


#' @title Transform prior samples to original scale
#'
#' @description Generate prior samples and transform them using the same
#' matrix transformation as posterior samples. This is the correct approach for
#' visualizing priors on the original (unscaled) scale, especially for the intercept
#' which depends on contributions from multiple coefficient priors.
#'
#' @param fit a fitted model object with \code{prior_list} and optionally
#' \code{formula_scale} attributes
#' @param n_samples number of samples to generate (default: 10000)
#' @param seed random seed for reproducibility (optional)
#' @param formula_scale optional nested list containing standardization information.
#' If not provided, extracted from \code{fit} attribute.
#'
#' @details When models use auto-scaling (standardizing predictors), the posterior
#' samples are on the standardized scale. To correctly visualize priors on the
#' original scale, we cannot simply apply a linear transformation to individual
#' priors because the intercept on the original scale is a weighted sum of
#' multiple priors:
#'
#' \deqn{\beta_0^{orig} = \beta_0^* - \sum_i \frac{\mu_i}{\sigma_i} \beta_i^*}
#'

#' This function generates samples from ALL priors simultaneously and applies
#' the same matrix transformation used for posterior samples, which correctly
#' handles the intercept and all other parameters.
#'
#' @return A matrix of prior samples on the original (unscaled) scale, with
#' columns matching the structure of posterior samples.
#'
#' @seealso [transform_scale_samples()] [plot_posterior()]
#'
#' @examples
#' # With a fitted model that used formula_scale:
#' # prior_samples <- transform_prior_samples(fit, n_samples = 10000)
#' # This can then be used with density() or for custom plotting
#'
#' @export
transform_prior_samples <- function(fit, n_samples = 10000, seed = NULL, formula_scale = NULL){

  check_int(n_samples, "n_samples", lower = 1)
  check_int(seed, "seed", allow_NULL = TRUE)

  # Extract prior_list from fit

  prior_list <- attr(fit, "prior_list")

  if(is.null(prior_list)){
    stop("'fit' must have 'prior_list' attribute.")
  }

  # Extract formula_scale from fit if not provided
  if(is.null(formula_scale)){
    formula_scale <- attr(fit, "formula_scale")
  }

  # Get posterior column names for structure matching
  if(inherits(fit, "runjags") || inherits(fit, "BayesTools_fit")){
    posterior <- as.matrix(.fit_to_posterior(fit))
  }else{
    stop("'fit' must be a fitted model object.")
  }

  prior_samples <- .generate_transformed_prior_samples(
    prior_list   = prior_list,
    column_names = colnames(posterior),
    n_samples    = n_samples,
    seed         = seed,
    formula_scale = formula_scale
  )

  return(prior_samples)
}


# Helper: Generate prior samples and apply the same unscaling transform used
# for posterior samples.
#
# @param prior_list Named list of prior objects
# @param column_names Column names to match from the posterior structure
# @param n_samples Number of samples to generate
# @param seed Optional random seed
# @param formula_scale Optional nested scaling information
# @return Matrix with transformed prior samples
.generate_transformed_prior_samples <- function(prior_list, column_names, n_samples, seed = NULL, formula_scale = NULL){

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    .check_formula_scale_info(formula_scale)
  }

  prior_samples <- .generate_prior_sample_matrix(
    prior_list     = prior_list,
    n_samples      = n_samples,
    column_names   = column_names,
    seed           = seed
  )

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    prior_samples <- .apply_unscale_transform(prior_samples, formula_scale)
  }

  return(prior_samples)
}

.generate_factor_prior_sample_matrix <- function(prior, parameter, n_samples){

  K <- .get_prior_factor_levels(prior)
  if(is.null(K) || is.na(K)){
    stop("The number of factor coefficients for prior '", parameter, "' is unknown.", call. = FALSE)
  }

  if(is.prior.spike_and_slab(prior)){
    prior_variable  <- .get_spike_and_slab_variable(prior)
    prior_inclusion <- .get_spike_and_slab_inclusion(prior)
    samples <- .generate_factor_prior_sample_matrix(prior_variable, parameter, n_samples)
    inclusion <- stats::rbinom(n_samples, size = 1, prob = rng(prior_inclusion, n_samples))
    samples <- samples * inclusion
  }else if(is.prior.mixture(prior)){
    prior_weights <- attr(prior, "prior_weights")
    prior_weights <- prior_weights / sum(prior_weights)
    components <- sample(seq_along(prior_weights), size = n_samples, replace = TRUE, prob = prior_weights)
    samples <- matrix(NA_real_, nrow = n_samples, ncol = K)

    for(component in unique(components)){
      samples[components == component, ] <- .generate_factor_prior_sample_matrix(
        prior      = prior[[component]],
        parameter  = parameter,
        n_samples  = sum(components == component)
      )
    }
  }else if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    samples <- matrix(rep(location, length.out = K), nrow = n_samples, ncol = K, byrow = TRUE)
  }else if(is.prior.orthonormal(prior) || is.prior.meandif(prior)){
    prior$parameters[["K"]] <- K
    samples <- rng(prior, n_samples, transform_factor_samples = FALSE)
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }else if(is.prior.treatment(prior) || is.prior.independent(prior)){
    samples <- replicate(K, rng(prior, n_samples, transform_factor_samples = FALSE))
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }else{
    samples <- rng(prior, n_samples, transform_factor_samples = FALSE)
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }

  colnames(samples) <- .JAGS_prior_factor_names(parameter, prior)
  return(samples)
}


# Helper: Generate a matrix of prior samples matching posterior structure
#
# @param prior_list Named list of prior objects
# @param n_samples Number of samples to generate
# @param column_names Optional vector of column names to match (filters output)
# @param seed Optional random seed
# @return Matrix with prior samples (rows = samples, columns = parameters)
.generate_prior_sample_matrix <- function(prior_list, n_samples, column_names = NULL, seed = NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  # Determine which parameters to sample
  param_names <- names(prior_list)

  if(is.null(param_names) || length(param_names) == 0){
    stop("'prior_list' must be a named list of priors.")
  }

  # Initialize list to collect samples (handles varying column counts per prior)
  samples_list <- list()

  for(param_name in param_names){
    prior <- prior_list[[param_name]]

    if(is.null(prior)){
      # No prior for this parameter - use zeros
      samples_list[[param_name]] <- matrix(0, nrow = n_samples, ncol = 1)
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.none(prior)){
      # No effect prior - use zeros
      samples_list[[param_name]] <- matrix(0, nrow = n_samples, ncol = 1)
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.factor(prior) || inherits(prior, "prior.factor_mixture") || inherits(prior, "prior.factor_spike_and_slab")){
      samples_list[[param_name]] <- .generate_factor_prior_sample_matrix(prior, param_name, n_samples)

    }else if(is.prior.point(prior)){
      # Point prior - constant values
      samples_list[[param_name]] <- matrix(
        prior$parameters[["location"]],
        nrow = n_samples,
        ncol = 1
      )
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.simple(prior)){
      # Simple priors - single column
      samples_list[[param_name]] <- matrix(
        rng(prior, n_samples),
        nrow = n_samples,
        ncol = 1
      )
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.vector(prior)){
      # Vector priors return matrix from rng
      temp_samples <- rng(prior, n_samples)
      if(is.matrix(temp_samples)){
        n_cols <- ncol(temp_samples)
        col_names <- paste0(param_name, "[", 1:n_cols, "]")
        colnames(temp_samples) <- col_names
        samples_list[[param_name]] <- temp_samples
      }else{
        samples_list[[param_name]] <- matrix(temp_samples, nrow = n_samples, ncol = 1)
        colnames(samples_list[[param_name]]) <- param_name
      }

    }else{
      temp_samples <- tryCatch(
        rng(prior, n_samples),
        error = function(e){
          stop(
            "Could not generate samples for prior '", param_name, "': ",
            e$message,
            call. = FALSE
          )
        }
      )
      if(!is.numeric(temp_samples)){
        stop(
          "Could not generate samples for prior '", param_name,
          "': rng() did not return numeric samples.",
          call. = FALSE
        )
      }

      if(is.matrix(temp_samples)){
        n_cols <- ncol(temp_samples)
        col_names <- paste0(param_name, "[", 1:n_cols, "]")
        colnames(temp_samples) <- col_names
        samples_list[[param_name]] <- temp_samples
      }else{
        samples_list[[param_name]] <- matrix(temp_samples, nrow = n_samples, ncol = 1)
        colnames(samples_list[[param_name]]) <- param_name
      }
    }
  }

  # Combine all samples into one matrix
  samples <- do.call(cbind, samples_list)

  # Filter to match column_names if provided
  if(!is.null(column_names)){
    available_cols <- intersect(column_names, colnames(samples))
    if(length(available_cols) > 0){
      samples <- samples[, available_cols, drop = FALSE]
    }
  }

  return(samples)
}


#' @title BayesTools Contrast Matrices
#'
#' @description BayesTools provides several contrast matrix functions for Bayesian factor analysis.
#' These functions create different types of contrast matrices that can be used with factor
#' variables in Bayesian models.
#'
#' @details
#' The package includes the following contrast functions:
#' \describe{
#'   \item{\code{contr.orthonormal}}{Return a matrix of orthonormal contrasts.
#'     Code is based on \code{stanova::contr.bayes} and corresponding to description
#'     by \insertCite{rouder2012default;textual}{BayesTools}. Returns a matrix with n rows and
#'     k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.meandif}}{Return a matrix of mean difference contrasts.
#'     This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#'     distributions on difference between the gran mean and factor level are identical independent
#'     of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#'     the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#'     distribution on difference between each factor level and the grand mean -- this is approximately
#'     twice the scale of \code{contr.orthonormal}. Returns a matrix with n rows and k columns,
#'     with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.independent}}{Return a matrix of independent contrasts -- a level for each term.
#'     Returns a matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#'     if \code{contrasts = FALSE}.}
#' }
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' # Orthonormal contrasts
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' # Mean difference contrasts
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' # Independent contrasts
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @aliases contr.orthonormal contr.meandif contr.independent
#' @name contr.BayesTools
NULL

#' @rdname contr.BayesTools
#' @export
contr.orthonormal <- function(n, contrasts = TRUE){
  # based on: stanova::contr.bayes
  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(n)
  if(contrasts){
    a       <- n
    I_a     <- diag(a)
    J_a     <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont    <- eigen(Sigma_a)$vectors[, seq_len(a - 1), drop = FALSE]
  }

  return(cont)
}

#' @rdname contr.BayesTools
#' @export
contr.meandif <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(n)
  if(contrasts){
    a       <- n
    I_a     <- diag(a)
    J_a     <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont    <- eigen(Sigma_a)$vectors[, seq_len(a - 1), drop = FALSE]
    cont    <- cont / (sqrt(1 - 1/n))
  }

  return(cont)
}


#' @rdname contr.BayesTools
#' @export
contr.independent <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n >= 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(x = 1, nrow = n, ncol = n)

  return(cont)
}


#' @title JAGS parameter-name helpers
#'
#' @description Format formula-derived parameter names and extract exact indexed
#' JAGS parameter columns.
#'
#' @param parameters a vector of parameter names
#' @param formula_parameter a formula parameter prefix name
#' @param formula_parameters a vector of formula parameter prefix names
#' @param formula_random a vector of random effects grouping factors
#' @param formula_prefix whether the \code{formula_parameters} names should be
#' kept. Defaults to \code{TRUE}.
#' @param formula_scale optional nested list containing scaling info. When provided,
#' intercepts from parameters with \code{log_intercept = TRUE} attribute will be
#' renamed to \code{exp(intercept)}.
#' @param x character vector to escape for regular-expression matching.
#' @param columns character vector of column names.
#' @param parameter exact unindexed JAGS parameter name.
#' @param samples matrix, data frame, \code{mcmc}, or \code{mcmc.list} object.
#' @param drop_missing whether missing indexed columns should return
#' \code{NULL}. If \code{FALSE}, a zero-column matrix with the input row count
#' is returned.
#' @param row named row vector, list, or one-row data frame.
#'
#' @examples
#' format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"),
#'                        formula_parameters = "mu")
#'
#' @return \code{format_parameter_names()}, \code{JAGS_parameter_names()}, and
#' \code{JAGS_regex_escape()} return character vectors.
#' \code{JAGS_indexed_parameter_columns()} returns a logical vector.
#' \code{JAGS_indexed_parameter_matrix()} returns a matrix of exact indexed
#' parameter columns sorted by numeric index, or \code{NULL} when none are
#' present and \code{drop_missing = TRUE}.
#' \code{JAGS_indexed_parameter_vector()} returns indexed row values sorted by
#' numeric index.
#'
#' @export format_parameter_names
#' @export JAGS_parameter_names
#' @export JAGS_regex_escape
#' @export JAGS_indexed_parameter_columns
#' @export JAGS_indexed_parameter_matrix
#' @export JAGS_indexed_parameter_vector
#' @name parameter_names
NULL

#' @rdname parameter_names
format_parameter_names <- function(parameters, formula_parameters = NULL, formula_random = NULL, formula_prefix = TRUE, formula_scale = NULL){

  check_char(parameters, "parameters", check_length = FALSE)
  check_char(formula_random, "formula_random", check_length = FALSE, allow_NULL = TRUE)
  check_char(formula_parameters, "formula_parameters", check_length = FALSE, allow_NULL = TRUE)
  check_bool(formula_prefix, "formula_prefix")
  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)

  # rename intercept to exp(intercept) for parameters with log_intercept attribute
  if(!is.null(formula_scale)){
    for(param_name in names(formula_scale)){
      if(isTRUE(attr(formula_scale[[param_name]], "log_intercept"))){
        intercept_name <- paste0(param_name, "_intercept")
        if(intercept_name %in% parameters){
          parameters[parameters == intercept_name] <- paste0(param_name, "_exp(intercept)")
        }
      }
    }
  }

  for(i in seq_along(formula_parameters)){
    parameters[grep(paste0(formula_parameters[i], "_"), parameters)] <- gsub(
      paste0(formula_parameters[i], "_"),
      if(formula_prefix) paste0("(", formula_parameters[i], ") ") else "",
      parameters[grep(paste0(formula_parameters[i], "_"), parameters)])
  }

  for(i in seq_along(formula_random)){
    temp_which <- grepl(paste0("_xREx__", formula_random[i], "_"), parameters)
    temp_incl  <- grepl("(inclusion)", parameters)
    parameters[temp_which] <- gsub(
      paste0("_xREx__", formula_random[i], "_"),
      "",
      parameters[temp_which]
    )
    if(any(temp_which &  temp_incl)){
      parameters[temp_which &  temp_incl] <- paste0(gsub("(inclusion)", "", parameters[temp_which & temp_incl], fixed = TRUE), "|", formula_random[i], " (inclusion)")
    }
    if(any(temp_which & !temp_incl)){
      parameters[temp_which & !temp_incl] <- paste0("sd(", parameters[temp_which & !temp_incl], "|", formula_random[i], ")")
    }
  }

  parameters[grep("__xXx__", parameters)] <- gsub("__xXx__", ":", parameters[grep("__xXx__", parameters)])

  return(parameters)
}
#' @rdname parameter_names
JAGS_parameter_names   <- function(parameters, formula_parameter = NULL){

  check_char(parameters, "parameters", check_length = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = TRUE, allow_NULL = TRUE)

  if(!is.null(formula_parameter)){
    parameters <- paste0(formula_parameter, "_", parameters)
  }
  parameters <- gsub(":", "__xXx__", parameters)

  return(parameters)
}

#' @rdname parameter_names
JAGS_regex_escape <- function(x){

  if(!is.character(x) || !is.vector(x)){
    stop("The 'x' argument must be a character vector.", call. = FALSE)
  }
  if(anyNA(x)){
    stop("The 'x' argument cannot contain NA/NaN values.", call. = FALSE)
  }
  if(length(x) == 0L){
    return(character())
  }

  special <- c("\\", ".", "|", "(", ")", "[", "]", "{", "}", "^", "$",
               "*", "+", "?")
  out <- vapply(x, function(x_i){
    chars <- strsplit(x_i, "", fixed = TRUE)[[1L]]
    chars <- vapply(chars, function(ch){
      if(ch %in% special){
        return(paste0("\\", ch))
      }
      ch
    }, character(1))
    paste0(chars, collapse = "")
  }, character(1))

  return(unname(out))
}

#' @rdname parameter_names
JAGS_indexed_parameter_columns <- function(columns, parameter){

  .JAGS_check_columns(columns)
  .JAGS_check_indexed_parameter(parameter)

  return(!is.na(.JAGS_indexed_parameter_indices(columns, parameter)))
}

#' @rdname parameter_names
JAGS_indexed_parameter_matrix <- function(samples, parameter,
                                         drop_missing = TRUE){

  check_bool(drop_missing, "drop_missing", allow_NA = FALSE)
  .JAGS_check_indexed_parameter(parameter)

  if(inherits(samples, "mcmc.list") || inherits(samples, "mcmc")){
    samples <- as.matrix(samples)
  }
  if(!is.matrix(samples) && !is.data.frame(samples)){
    stop("'samples' must be a matrix, data frame, mcmc, or mcmc.list object.",
         call. = FALSE)
  }
  if(is.data.frame(samples)){
    samples <- as.matrix(samples)
  }

  columns <- colnames(samples)
  if(is.null(columns)){
    if(drop_missing){
      return(NULL)
    }
    return(samples[, integer(0), drop = FALSE])
  }

  selected <- .JAGS_indexed_parameter_sorted_columns(columns, parameter)
  if(length(selected) == 0L){
    if(drop_missing){
      return(NULL)
    }
    return(samples[, integer(0), drop = FALSE])
  }

  return(samples[, selected, drop = FALSE])
}

#' @rdname parameter_names
JAGS_indexed_parameter_vector <- function(row, parameter){

  .JAGS_check_indexed_parameter(parameter)

  if(is.data.frame(row)){
    if(nrow(row) != 1L){
      stop("'row' data frames must contain exactly one row.", call. = FALSE)
    }
    row <- unlist(row[1L, , drop = TRUE], use.names = TRUE)
  }else if(is.list(row) && !is.atomic(row)){
    row <- unlist(row, use.names = TRUE)
  }

  row_names <- names(row)
  if(is.null(row_names)){
    stop("'row' must be a named vector or one-row data frame.", call. = FALSE)
  }

  selected <- .JAGS_indexed_parameter_sorted_columns(row_names, parameter)
  if(length(selected) == 0L){
    return(row[integer(0)])
  }

  return(row[selected])
}

.JAGS_check_columns <- function(columns){

  if(!is.character(columns) || !is.vector(columns)){
    stop("The 'columns' argument must be a character vector.",
         call. = FALSE)
  }
  if(anyNA(columns)){
    stop("The 'columns' argument cannot contain NA/NaN values.",
         call. = FALSE)
  }

  return(invisible(TRUE))
}

.JAGS_check_indexed_parameter <- function(parameter){

  check_char(parameter, "parameter", check_length = 1, allow_NA = FALSE)
  if(!nzchar(parameter)){
    stop("The 'parameter' argument must not be empty.", call. = FALSE)
  }

  return(invisible(TRUE))
}

.JAGS_indexed_parameter_indices <- function(columns, parameter){

  pattern <- paste0("^", JAGS_regex_escape(parameter), "\\[([0-9]+)\\]$")
  matches <- regexec(pattern, columns, perl = TRUE)
  parts <- regmatches(columns, matches)
  out <- rep(NA_integer_, length(columns))
  has_match <- vapply(parts, length, integer(1)) == 2L
  out[has_match] <- as.integer(vapply(
    parts[has_match],
    `[[`,
    character(1),
    2L
  ))

  return(out)
}

.JAGS_indexed_parameter_sorted_columns <- function(columns, parameter){

  indices <- .JAGS_indexed_parameter_indices(columns, parameter)
  keep <- which(!is.na(indices))
  if(length(keep) == 0L){
    return(integer())
  }

  return(keep[order(indices[keep])])
}

.JAGS_prior_factor_names <- function(parameter, prior){

  levels <- .get_prior_factor_levels(prior)
  if((is.null(levels) || length(levels) == 0L || is.na(levels)) && "K" %in% names(prior[["parameters"]])){
    levels <- prior[["parameters"]][["K"]]
  }
  if(is.null(levels) || length(levels) == 0L || is.na(levels)){
    stop("Factor-prior dimensions must be available before constructing JAGS parameter names.", call. = FALSE)
  }

  if(levels == 1){
    par_names <- parameter
  }else{
    par_names <- paste0(parameter, "[", 1:levels, "]")
  }

  return(par_names)
}
