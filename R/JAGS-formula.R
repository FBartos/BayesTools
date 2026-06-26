#' @title Create JAGS formula syntax and data object
#'
#' @description Creates a JAGS formula syntax, prepares data input, and
#' returns modified prior list for further processing in the \code{JAGS_fit}
#' function.
#'
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored), or a `BayesTools_random_effects` object returned
#' by [random_effects_formula()]. If the formula contains \code{-1}, it will be
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
#' \describe{
#'   \item{\code{"__default_continuous"}}{A prior to use for any continuous predictors
#'     (including the intercept) that are not explicitly specified in the prior list.
#'     This can also be a zero-argument function returning a prior object; in that
#'     case it is evaluated only if a missing continuous term needs the default.}
#'   \item{\code{"__default_factor"}}{A prior to use for any factor predictors
#'     (including interactions involving factors) that are not explicitly specified
#'     in the prior list. This can also be a zero-argument function returning a prior
#'     object; in that case it is evaluated only if a missing factor term needs the
#'     default.}
#' }
#' These default priors allow for more concise specification when many predictors
#' share the same prior distribution.
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

  formula_input <- formula
  formula <- .bt_formula_random_formula(formula)
  if(!is.language(formula))
    stop("'formula' must be a formula", call. = FALSE)
  resolved_random_terms <- if(inherits(formula_input, "BayesTools_random_effects")){
    formula_input$terms
  }else{
    attr(formula, "random_terms", exact = TRUE)
  }
  check_char(parameter, "parameter")
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  .JAGS_formula_check_prior_list(prior_list)
  .bt_check_prior_random(prior_random, allow_NULL = TRUE)
  .bt_check_random_effects_compile(random_effects_compile, allow_NULL = TRUE)
  # formula_scale can be TRUE/FALSE (apply to all continuous predictors) or a named list
  if(!is.null(formula_scale) && !is.logical(formula_scale) && !is.list(formula_scale)){
    stop("'formula_scale' must be NULL, TRUE, FALSE, or a named list.", call. = FALSE)
  }


  # remove the specified response
  formula <- .remove_response(formula)
  formula <- .bt_formula_preserve_random_terms(formula, resolved_random_terms)
  # store log(intercept) attribute (for models relying on mu = log(intercept) + sum(beta_i * x_i) trick
  # exp(mu) = intercept * exp(sum(beta_i * x_i)) (e.g., Poisson regression / regression with log link etc...)
  log_intercept  <- isTRUE(attr(formula, "log(intercept)"))
  # store expressions (included later as the literal character input)
  expressions    <- .extract_expressions(formula)
  # store random effects (included later via a formula interface)
  parsed_random_effects <- .bt_formula_random_terms(formula)
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
    missing_terms <- model_terms[!model_terms %in% names(prior_list)]
    if(any(model_terms_type[missing_terms] == "continuous") && !is.null(default_continuous_prior)){
      default_continuous_prior <- .JAGS_formula_resolve_default_prior(
        default_prior = default_continuous_prior,
        default_name  = "__default_continuous"
      )
    }
    if(any(model_terms_type[missing_terms] == "factor") && !is.null(default_factor_prior)){
      default_factor_prior <- .JAGS_formula_resolve_default_prior(
        default_prior = default_factor_prior,
        default_name  = "__default_factor"
      )
    }

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
      if(is.prior.ordered(this_prior) && !.is_prior_interaction(this_prior)){
        attr(this_prior, "levels") <- length(levels(data[[model_terms[i]]]))
      }else if(is.prior.independent(this_prior)){
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
      if(is.prior.ordered(this_prior)){
        this_prior <- .bt_bind_ordered_prior_metadata(this_prior, paste0(parameter, "_", model_terms[i]))
      }

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
        attr(this_prior, "coefficient_dim")   -> attr(this_prior[[p]], "coefficient_dim")
        attr(this_prior, "ordered_metadata")  -> attr(this_prior[[p]], "ordered_metadata")
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
  if(length(random_sd_binding_context$add_parameters) > 0L){
    add_parameters <- c(add_parameters, random_sd_binding_context$add_parameters)
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
  .bt_validate_ordered_shared_allocations(prior_list)
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

