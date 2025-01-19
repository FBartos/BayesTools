#' @title Create JAGS formula syntax and data object
#'
#' @description Creates a JAGS formula syntax, prepares data input, and
#' returns modified prior list for further processing in the \code{JAGS_fit}
#' function
#'
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored)
#' @param parameter name of the parameter to be created with the formula
#' @param data data.frame containing predictors included in the formula
#' @param prior_list named list of prior distribution of parameters specified within
#' the \code{formula}
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
#' # specify priors
#' prior_list <- list(
#' "intercept"     = prior("normal", list(0, 1)),
#' "x_cont"        = prior("normal", list(0, .5)),
#' "x_fac3"        = prior_factor("normal",  list(0, 1),  contrast = "treatment"),
#' "x_fac4"        = prior_factor("mnormal", list(0, 1),  contrast = "orthonormal"),
#' "x_fac3:x_fac4" = prior_factor("mnormal", list(0, .5), contrast = "orthonormal")
#' )
#'
#' # create the formula object
#' formula <- JAGS_formula(
#'   formula = ~ x_cont + x_fac3 * x_fac4,
#'   parameter = "mu", data = df, prior_list = prior_list)
#'
#' @return \code{JAGS_formula} returns a list containing the formula JAGS syntax,
#' JAGS data object, and modified prior_list.
#'
#' @seealso [JAGS_fit()]
#' @export
JAGS_formula <- function(formula, parameter, data, prior_list){

  if(!is.language(formula))
    stop("'formula' must be a formula")
  check_char("parameter", parameter)
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # remove the specified response
  formula <- .remove_response(formula)
  # store expressions (included later as the literal character input)
  expressions    <- .extract_expressions(formula)
  # store random effects (included later via a formula interface)
  random_effects <- .extract_random_effects(formula)
  # remove expressions and random effects from the formula
  formula <- .remove_expressions(formula)
  formula <- .remove_random_effects(formula)

  # obtain predictors characteristics factors
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

  # separate prior lists: extract random effects priors
  prior_list_random_effects <- list()
  if(length(random_effects)){
    # store the random effects specific priors in the corresponding entry
    for(i in seq_along(random_effects)){
      prior_list_random_effects[[i]] <- prior_list[which(.get_grouping_factor(names(prior_list)) == attr(random_effects[[i]], "grouping_factor"))]
    }
    # remove the random effects specific priors from the prior list
    prior_list <- prior_list[.get_grouping_factor(names(prior_list)) == ""]
  }
  # check that all predictors have a prior distribution
  check_list(prior_list, "prior_list", check_names = model_terms, allow_other = FALSE, all_objects = TRUE)

  # check the prior distribution for each predictor
  # assign factor contrasts to the data based on prior distributions
  if(any(predictors_type == "factor")){

    for(factor in names(predictors_type[predictors_type == "factor"])){

      # select the corresponding prior for the variable
      this_prior <- prior_list[[factor]]

      if(is.prior.treatment(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.treatment"
      }else if(is.prior.independent(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.independent"
      }else if(is.prior.orthonormal(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.orthonormal"
      }else if(is.prior.meandif(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.meandif"
      }else{
        stop(paste0("Unsupported prior distribution defined for '", factor, "' factor variable. See '?prior_factor' for details."))
      }
    }
  }
  if(any(predictors_type == "continuous")){

    for(continuous in names(predictors_type[predictors_type == "continuous"])){

      # select the corresponding prior for the variable
      this_prior <- prior_list[[continuous]]

      if(is.prior.factor(this_prior)|| is.prior.discrete(this_prior) || is.prior.PET(this_prior) || is.prior.PEESE(this_prior) || is.prior.weightfunction(this_prior)){
        stop(paste0("Unsupported prior distribution defined for '", continuous, "' continuous variable. See '?prior' for details."))
      }
    }
  }

  # get the default design matrix
  model_frame  <- stats::model.frame(formula, data = data)
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)

  # check whether intercept is unique parameter
  if(sum(grepl("intercept", names(prior_list))) > 1)
    stop("only the intercept parameter can contain 'intercept' in its name.")
  # check whether the interaction replacement is in usage
  if(any(grepl("__xXx__", names(prior_list))))
    stop("'__xXx__' string is internally used by the BayesTools package and can't be used for naming variables.")
  if(any(grepl("__xREx__", names(prior_list))))
    stop("'__xREx__' string is internally used by the BayesTools package and can't be used for naming variables.")

  # replace interaction signs (due to JAGS incompatibility)
  colnames(model_matrix)  <- gsub(":", "__xXx__", colnames(model_matrix))
  names(prior_list)       <- gsub(":", "__xXx__", names(prior_list))
  names(model_terms_type) <- gsub(":", "__xXx__", names(model_terms_type))
  model_terms             <- gsub(":", "__xXx__", model_terms)

  # prepare syntax & data based on the formula
  formula_syntax <- NULL
  JAGS_data      <- list()
  JAGS_data[[paste0("N_", parameter)]] <- nrow(data)

  # add intercept and prepare the indexing vector
  if(has_intercept){
    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    formula_syntax   <- c(formula_syntax, paste0(parameter, "_intercept"))
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
      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i]

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
          if(model_terms_type[sub_term] == "factor"){
            level_names[[sub_term]] <- levels(data[,sub_term])
          }
        }
        attr(this_prior, "level_names") <- level_names
      }else{
        attr(this_prior, "level_names") <- levels(data[,model_terms[i]])
      }

      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i, drop = FALSE]
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
    if(is.prior.spike_and_slab(this_prior)){
      attr(this_prior, "levels")            -> attr(this_prior[["variable"]], "levels")
      attr(this_prior, "level_names")       -> attr(this_prior[["variable"]], "level_names")
      attr(this_prior, "interaction")       -> attr(this_prior[["variable"]], "interaction")
      attr(this_prior, "interaction_terms") -> attr(this_prior[["variable"]], "interaction_terms")
      this_prior -> prior_list[[model_terms[i]]]
    }else if(is.prior.mixture(this_prior)){
      for(p in seq_along(this_prior)){
        attr(this_prior, "levels")            -> attr(this_prior[[p]], "levels")
        attr(this_prior, "level_names")       -> attr(this_prior[[p]], "level_names")
        attr(this_prior, "interaction")       -> attr(this_prior[[p]], "interaction")
        attr(this_prior, "interaction_terms") -> attr(this_prior[[p]], "interaction_terms")
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
  for(i in seq_along(random_effects)){
    # TODO
    temp_random   <- .JAGS_random_effect_formula(random_effects[[i]], parameter, data, prior_list_random_effects[[i]])

    formula_syntax <- c(formula_syntax, temp_random$formula_syntax)
  }

  # finish the syntax
  formula_syntax <- paste0(
    "for(i in 1:N_", parameter, "){\n",
    "  ", parameter, "[i] = ", paste0(formula_syntax, collapse = " + "), "\n",
    "}\n")

  # add the parameter name as a prefix and attribute to each prior in the list
  names(prior_list) <- paste0(parameter, "_", names(prior_list))
  for(i in seq_along(prior_list)){
    attr(prior_list[[i]], "parameter") <- parameter
  }

  return(list(
    formula_syntax = formula_syntax,
    data           = JAGS_data,
    prior_list     = prior_list,
    formula        = formula
  ))
}

.JAGS_random_effect_formula <- function(formula, parameter, data, prior_list){

  # extract the grouping factor information
  grouping_factor        <- attr(formula, "grouping_factor")
  grouping_factor_levels <- levels(as.factor(data[[grouping_factor]]))

  # remove the grouping factor from the formula
  formula <- .remove_grouping_factor(formula)
  formula <- as.formula(paste("~", formula))

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

  # TODO: expand to factor random effects
  if(any(model_terms_type) != "continuous")
    stop("Random effects can only be continuous.")

  # check that all priors have a lower bound on 0 or their range is > 0, if not, throw a warning and correct
  for(i in seq_along(prior_list)){
    if(range(prior_list[[i]])[1] < 0){
      warning(paste0("The lower bound of the '", names(prior_list)[i], "' prior distribution is below 0. Correcting to 0."), immediate. = TRUE)
      prior_list[[i]] <- prior_list[[i]]$truncation$lower <- 0
    }
  }

  # drop the grouping factor from the prior names
  names(prior_list) <- .remove_grouping_factor(names(prior_list))
  # check that all terms have a prior
  check_list(prior_list, "prior_list", check_names = model_terms, allow_other = TRUE, all_objects = TRUE)


  # prepare syntax & data based on the formula
  parameter      <- paste0(parameter, "__xREx__", grouping_factor)
  formula_syntax <- NULL
  JAGS_data      <- list()
  new_prior_list <- list()
  JAGS_data[[paste0("N_", parameter)]] <- length(grouping_factor_levels)

  # TODO: HERE - add expressions and dummy random_effect priros (maybe define new class?)
  if(has_intercept){
    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    formula_syntax   <- c(formula_syntax, paste0(parameter, "_intercept"))
  }else{
    terms_indexes    <- attr(model_matrix, "assign")
  }

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
  return(any(grepl("expression\\(", deparse(formula))))
}
.extract_expressions    <- function(formula){
  # extract all expressions from the formula

  # Convert the formula to a character string
  formula_string <- deparse(formula)

  # Use a regex to find all instances of "expression(...)"
  matches     <- gregexpr("expression\\(.*?\\)", formula_string)
  expressions <- regmatches(formula_string, matches)[[1]]

  # Use a regex to remove "expression(" and the closing ")"
  expressions <- lapply(expressions, .clean_from_expression)

  return(expressions)
}
.clean_from_expression  <- function(x){
  # expression to character

  return(sub("expression\\((.*)\\)", "\\1", x))
}
.remove_expressions     <- function(formula){
  # remove all expressions from the formula

  # Convert the formula to a character string
  formula_string <- deparse(formula)

  # Use a regex to remove all instances of "+ expression(...)" or "expression(...) +", considering spaces and newlines
  formula_string_clean <- gsub("\\+\\s*expression\\(.*?\\)\\s*", "", formula_string)
  formula_string_clean <- gsub("\\s*expression\\(.*?\\)\\s*\\+", "", formula_string_clean)

  # Handle the case where the expression is the first term in the formula
  formula_string_clean <- gsub("^\\s*expression\\(.*?\\)\\s*", "", formula_string_clean)

  # Handle the case where the formula reduces to just "y ~ expression(...)"
  if(grepl("^\\s*[a-zA-Z0-9._]+\\s*~\\s*expression\\(.*?\\)\\s*$", formula_string_clean)){
    formula_string_clean <- gsub("expression\\(.*?\\)", "1", formula_string_clean)
  }

  # Reconvert the cleaned string back to a formula
  return(as.formula(formula_string_clean))
}
.has_random_effects     <- function(formula){
  # Convert the formula to a character string
  formula_str <- paste(deparse(formula), collapse = " ")

  # Regular expression to match `( ... | ... )` patterns
  has_random <- grepl("\\([^\\)]+\\|[^\\)]+\\)", formula_str)

  # Return TRUE if at least one match is found, otherwise FALSE
  return(has_random)
}
.extract_random_effects <- function(formula){
  # Convert the formula to a character string
  formula_str <- paste(deparse(formula), collapse = " ")

  # Regular expression to match `( ... | ... )` patterns
  random_effects <- gregexpr("\\([^\\)]+\\|[^\\)]+\\)", formula_str)

  # Extract matches
  matches <- regmatches(formula_str, random_effects)

  # Clean up the parentheses and remove unnecessary spaces
  clean_matches <- lapply(unlist(matches), function(x) gsub("^\\(|\\)$", "", x))  # Remove outer parentheses
  clean_matches <- lapply(clean_matches, trimws)  # Remove extra spaces from each match

  # Store the conditioning variable at an attribute
  for(i in seq_along(clean_matches)){
    attr(clean_matches[[i]], "grouping_factor") <- .get_grouping_factor(clean_matches[[i]])
  }

  # Return the cleaned random effects as a list
  return(as.list(clean_matches))
}
.remove_random_effects  <- function(formula){
  # Convert the formula to a character string
  formula_str <- paste(deparse(formula), collapse = " ")

  # Regular expression to match and remove `( ... | ... )` patterns
  cleaned_formula <- gsub("\\+?\\s*\\([^\\)]+\\|[^\\)]+\\)", "", formula_str)

  # Normalize spacing around '+' and remove any leading '+'
  cleaned_formula <- gsub("\\s*\\+\\s*", " + ", cleaned_formula)  # Normalize '+' spacing
  cleaned_formula <- gsub("^\\s*~\\s*\\+\\s*", "~ ", cleaned_formula)  # Remove leading '+'

  # Ensure no excessive spaces remain
  cleaned_formula <- gsub("\\s{2,}", " ", cleaned_formula)  # Replace multiple spaces with a single space
  cleaned_formula <- trimws(cleaned_formula)  # Trim leading/trailing whitespace

  # Ensure at least "1" remains if formula is empty after cleaning
  if (grepl("^\\s*~\\s*$", cleaned_formula)) {
    cleaned_formula <- "~ 1"
  }

  # Return as a formula
  return(as.formula(cleaned_formula))
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

#' @title Evaluate JAGS formula using posterior samples
#'
#' @description Evaluates a JAGS formula on a posterior distribution
#' obtained from a fitted model.
#'
#' @param fit model fitted with either \link[runjags]{runjags} posterior
#' samples obtained with \link[rjags]{rjags-package}
#' @param formula formula specifying the right hand side of the assignment (the
#' left hand side is ignored)
#' @param parameter name of the parameter created with the formula
#' @param data data.frame containing predictors included in the formula
#' @param prior_list named list of prior distribution of parameters specified within
#' the \code{formula}
#'
#'
#' @return \code{JAGS_evaluate_formula} returns a matrix of the evaluated posterior samples on
#' the supplied data.
#'
#' @seealso [JAGS_fit()] [JAGS_formula()]
#' @export
JAGS_evaluate_formula <- function(fit, formula, parameter, data, prior_list){

  if(!is.language(formula))
    stop("'formula' must be a formula")
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_char(parameter, "parameter")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")

  # extract the posterior distribution
  posterior <- as.matrix(.fit_to_posterior(fit))

  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)

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
        stats::contrasts(data[,factor]) <- "contr.orthonormal"
      }else if(is.prior.meandif(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.meandif"
      }else if(is.prior.independent(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.independent"
      }else if(is.prior.treatment(this_prior)){
        stats::contrasts(data[,factor]) <- "contr.treatment"
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
  }

  # get the design matrix
  model_frame  <- stats::model.frame(formula, data = data)
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)

  ### evaluate the design matrix on the samples -> output[data, posterior]
  if(has_intercept){

    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    # check for scaling factors
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term = "intercept", prior_list = prior_list_formula, posterior = posterior, nrow = nrow(data), ncol = nrow(posterior))

    output           <- temp_multiply_by * matrix(posterior[,JAGS_parameter_names("intercept", formula_parameter = parameter)], nrow = nrow(data), ncol = nrow(posterior), byrow = TRUE)

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

#' @title Transform factor posterior samples into differences from the mean
#'
#' @description Transforms posterior samples from model-averaged posterior
#' distributions based on meandif/orthonormal prior distributions into differences from
#' the mean.
#'
#' @param samples (a list) of mixed posterior distributions created with
#' \code{mix_posteriors} function
#'
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
#'
#' @return \code{transform_meandif_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.meandif]
#'
#' @export
transform_meandif_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.meandif_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && attr(samples[[i]], "meandif")){

      meandif_samples     <- samples[[i]]
      transformed_samples <- meandif_samples %*% t(contr.meandif(1:(attr(samples[[i]], "levels")+1)))

      if(attr(samples[[i]], "interaction")){
        if(length(attr(samples[[i]], "level_names")) == 1){
          transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names")[[1]],"]")
        }else{
          stop("meandif de-transformation for interaction of multiple factors is not implemented.")
        }
      }else{
        transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names"),"]")
      }

      colnames(transformed_samples)   <- transformed_names
      attributes(transformed_samples) <- c(attributes(transformed_samples), attributes(meandif_samples)[!names(attributes(meandif_samples)) %in% names(attributes(transformed_samples))])
      class(transformed_samples)      <- c(class(transformed_samples), "mixed_posteriors.meandif_transformed")

      samples[[i]] <- transformed_samples
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
#'
#' @return \code{transform_orthonormal_samples} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors] [contr.orthonormal]
#'
#' @export
transform_orthonormal_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.orthonormal_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && attr(samples[[i]], "orthonormal")){

      orthonormal_samples <- samples[[i]]
      transformed_samples <- orthonormal_samples %*% t(contr.orthonormal(1:(attr(samples[[i]], "levels")+1)))

      if(attr(samples[[i]], "interaction")){
        if(length(attr(samples[[i]], "level_names")) == 1){
          transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names")[[1]],"]")
        }else{
          stop("orthonormal de-transformation for interaction of multiple factors is not implemented.")
        }
      }else{
        transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names"),"]")
      }

      colnames(transformed_samples)   <- transformed_names
      attributes(transformed_samples) <- c(attributes(transformed_samples), attributes(orthonormal_samples)[!names(attributes(orthonormal_samples)) %in% names(attributes(transformed_samples))])
      class(transformed_samples)      <- c(class(transformed_samples), "mixed_posteriors.orthonormal_transformed")

      samples[[i]] <- transformed_samples
    }
  }

  return(samples)
}

# not part of transform factor samples (as it's usefull only for marginal effects)
transform_treatment_samples <- function(samples){

  check_list(samples, "samples", allow_NULL = TRUE)

  for(i in seq_along(samples)){
    if(!inherits(samples[[i]],"mixed_posteriors.treatment_transformed") && inherits(samples[[i]], "mixed_posteriors.factor") && attr(samples[[i]], "treatment")){

      treatment_samples   <- samples[[i]]
      transformed_samples <- treatment_samples %*% t(stats::contr.treatment(1:(attr(samples[[i]], "levels")+1)))

      if(attr(samples[[i]], "interaction")){
        if(length(attr(samples[[i]], "level_names")) == 1){
          transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names")[[1]],"]")
        }else{
          stop("orthonormal de-transformation for interaction of multiple factors is not implemented.")
        }
      }else{
        transformed_names <- paste0(names(samples)[i], " [dif: ", attr(samples[[i]], "level_names"),"]")
      }

      colnames(transformed_samples)   <- transformed_names
      attributes(transformed_samples) <- c(attributes(transformed_samples), attributes(treatment_samples)[!names(attributes(treatment_samples)) %in% names(attributes(transformed_samples))])
      class(transformed_samples)      <- c(class(transformed_samples), "mixed_posteriors.treatment_transformed")

      samples[[i]] <- transformed_samples
    }
  }

  return(samples)
}


#' @title Orthornomal contrast matrix
#'
#' @description Return a matrix of orthornomal contrasts.
#' Code is based on \code{stanova::contr.bayes} and corresponding to description
#' by \insertCite{rouder2012default;textual}{BayesTools}
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
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

#' @title Mean difference contrast matrix
#'
#' @description Return a matrix of mean difference contrasts.
#' This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#' distributions on difference between the gran mean and factor level are identical independent
#' of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#' the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#' distribution on difference between each factor level and the grand mean -- this is approximately
#' twice the scale of \code{contr.orthonormal}.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
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


#' @title Independent contrast matrix
#'
#' @description Return a matrix of independent contrasts -- a level for each term.
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#' if \code{contrasts = FALSE}.
#'
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


#' @title Clean parameter names from JAGS
#'
#' @description Removes additional formatting from parameter names outputted from
#' JAGS.
#'
#' @param parameters a vector of parameter names
#' @param formula_parameter a formula parameter prefix name
#' @param formula_parameters a vector of formula parameter prefix names
#' @param formula_prefix whether the \code{formula_parameters} names should be
#' kept. Defaults to \code{TRUE}.
#'
#' @examples
#' format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"),
#'                        formula_parameters = "mu")
#'
#' @return A character vector with reformatted parameter names.
#'
#' @export format_parameter_names
#' @export JAGS_parameter_names
#' @name parameter_names
NULL

#' @rdname parameter_names
format_parameter_names <- function(parameters, formula_parameters = NULL, formula_prefix = TRUE){

  check_char(parameters, "parameters", check_length = FALSE)
  check_char(formula_parameters, "formula_parameters", check_length = FALSE, allow_NULL = TRUE)
  check_bool(formula_prefix, "formula_prefix")

  for(i in seq_along(formula_parameters)){
    parameters[grep(paste0(formula_parameters[i], "_"), parameters)] <- gsub(
      paste0(formula_parameters[i], "_"),
      if(formula_prefix) paste0("(", formula_parameters[i], ") ") else "",
      parameters[grep(paste0(formula_parameters[i], "_"), parameters)])
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

.JAGS_prior_factor_names <- function(parameter, prior){

  if(!.is_prior_interaction(prior)){
    if(.get_prior_factor_levels(prior) == 1){
      par_names <- parameter
    }else{
      par_names <- paste0(parameter,"[",1:.get_prior_factor_levels(prior),"]")
    }
  }else if(length(attr(prior, "levels")) == 1){
    par_names <-  paste0(parameter,"[",1:.get_prior_factor_levels(prior),"]")
  }

  return(par_names)
}
