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


  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)

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


  # check that all predictors have a prior distribution
  check_list(prior_list, "prior_list", check_names = model_terms, allow_other = FALSE, all_objects = TRUE)

  # assign factor contrasts to the data based on prior distributions
  if(any(predictors_type == "factor")){
    for(factor in names(predictors_type[predictors_type == "factor"])){
      if(is.prior.dummy(prior_list[[factor]])){
        stats::contrasts(data[,factor]) <- "contr.treatment"
      }else if(is.prior.orthonormal(prior_list[[factor]])){
        stats::contrasts(data[,factor]) <- "contr.orthonormal"
      }else if(is.prior.point(prior_list[[factor]])){
        class(prior_list[[factor]]) <-  c(class(prior_list[[factor]]), "prior.factor")
      }else{
        stop(paste0("Unsupported prior distribution defined for '", factor, "' factor variable. See '?prior_factor' for details."))
      }
    }
  }
  if(any(predictors_type == "continuous")){
    for(continuous in names(predictors_type[predictors_type == "continuous"])){
      if(!is.prior.simple(prior_list[[continuous]]) || is.prior.factor(prior_list[[continuous]])){
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
  # replace interaction signs (due to JAGS incompatibility)
  colnames(model_matrix)  <- gsub(":", "__xXx__", colnames(model_matrix))
  names(prior_list)       <- gsub(":", "__xXx__", names(prior_list))
  names(model_terms_type) <- gsub(":", "__xXx__", names(model_terms_type))
  model_terms             <- gsub(":", "__xXx__", model_terms)

  # prepare syntax & data based on the formula
  formula_syntax <- NULL
  prior_syntax   <- NULL
  JAGS_data      <- list()
  JAGS_data[[paste0("N_", parameter)]] <- nrow(data)

  # add intercept and prepare the indexing vector
  if(has_intercept){
    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- 0

    formula_syntax   <- c(formula_syntax, paste0(parameter, "_intercept"))
    prior_syntax     <- c(prior_syntax, .JAGS_prior.simple(prior_list[["intercept"]], paste0(parameter, "_intercept")))
  }else{
    terms_indexes    <- attr(model_matrix, "assign")
  }

  # add remaining terms (omitting the intercept indexed as NA)
  for(i in unique(terms_indexes[terms_indexes > 0])){

    attr(prior_list[[model_terms[i]]], "interaction") <- grepl("__xXx__", model_terms[i])

    if(model_terms_type[i] == "continuous"){
      # continuous variables or interactions of continuous variables are simple predictors

      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i]
      formula_syntax <- c(formula_syntax, paste0(
        if(!is.null(attr(prior_list[[model_terms[i]]], "multiply_by"))) paste0(attr(prior_list[[model_terms[i]]], "multiply_by"), " * "),
        parameter, "_", model_terms[i],
        " * ",
        parameter, "_data_", model_terms[i], "[i]"
      ))

    }else if(model_terms_type[i] == "factor"){
      # factor variables or interactions with a factor requires factor style prior

      # add levels information attributes to factors
      attr(prior_list[[model_terms[i]]], "levels") <- sum(terms_indexes == i) + 1
      if(attr(prior_list[[model_terms[i]]], "interaction")){
        level_names <- list()
        for(sub_term in strsplit(model_terms[i], "__xXx__")[[1]]){
          if(model_terms_type[sub_term] == "factor"){
            level_names[[sub_term]] <- levels(data[,sub_term])
          }
        }
        attr(prior_list[[model_terms[i]]], "level_names") <- level_names
      }else{
        attr(prior_list[[model_terms[i]]], "level_names") <- levels(data[,model_terms[i]])
      }

      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i, drop = FALSE]
      formula_syntax <- c(formula_syntax, paste0(
        if(!is.null(attr(prior_list[[model_terms[i]]], "multiply_by"))) paste0(attr(prior_list[[model_terms[i]]], "multiply_by"), " * "),
        "inprod(",
        parameter, "_", model_terms[i],
        ", ",
        parameter, "_data_", model_terms[i], "[i,])"
      ))

    }else{
      stop("Unrecognized model term.")
    }


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

.remove_response               <- function(formula){
  if(attr(stats::terms(formula), "response")  == 1){
    formula[2] <- NULL
  }
  return(formula)
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
    if(inherits(samples[[i]], "mixed_posteriors.factor") && attr(samples[[i]], "orthonormal")){

      orthonormal_samples <- samples[[i]]
      transformed_samples <- orthonormal_samples %*% t(contr.orthonormal(1:attr(samples[[i]], "levels")))

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


#' @title Orthornomal contrast natrix
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


#' @title Clean parameter names from JAGS
#'
#' @description Removes additional formatting from parameter names outputted from
#' JAGS.
#'
#' @param parameters a vector of parameter names
#' @param formula_parameters a vector of formula parameter prefix names
#' @param formula_prefix whether the \code{formula_parameters} names should be
#' kept. Defaults to \code{TRUE}.
#'
#' @examples
#' format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"), formula_parameters = "mu")
#'
#' @return A character vector with reformatted parameter names.
#'
#' @export
format_parameter_names <- function(parameters, formula_parameters = NULL, formula_prefix = TRUE){

  for(i in seq_along(formula_parameters)){
    parameters[grep(paste0(formula_parameters[i], "_"), parameters)] <- gsub(
      paste0(formula_parameters[i], "_"),
      if(formula_prefix) paste0("(", formula_parameters[i], ") ") else "",
      parameters[grep(paste0(formula_parameters[i], "_"), parameters)])
  }
  parameters[grep("__xXx__", parameters)] <- gsub("__xXx__", ":", parameters[grep("__xXx__", parameters)])

  return(parameters)
}
