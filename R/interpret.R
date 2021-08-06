#' @title Interpret ensemble inference and estimates
#'
#' @description Provides textual summary for posterior
#' distributions created by [mix_posteriors] and ensemble inference
#' created by [ensemble_inference].
#'
#' @param inference model inference created by [ensemble_inference]
#' @param samples posterior samples created by [mix_posteriors]
#' @param specification list of lists specifying the generated text.
#' Each inner list carries: (1) \code{inference} specifying the name of
#' in the \code{inference} entry and optionally \code{inference_name}
#' as a name to use in the text and \code{inference_BF_name} as a
#' symbol to be used instead of the default \code{"BF"}, (2) \code{samples}
#' specifying the name of in the \code{samples} entry and optionally
#' \code{samples_name} as a name to use in the text, \code{samples_units} as
#' a unit text to be appended after the estimate, and \code{samples_conditional}
#' specifying whether the estimate is conditional or model-averaged.
#' @param method character specifying name of the method to be
#' appended at the beginning of each sentence.
#'
#'
#' @return \code{interpret} returns character.
#'
#' @seealso [ensemble_inference] [mix_posteriors] [BayesTools_model_tables] [BayesTools_ensemble_tables]
#' @export
interpret <- function(inference, samples, specification, method){

  # check input
  check_list(specification, "specification", check_length = 0)
  sapply(specification, function(s){
    check_char(s$inference, allow_NULL = TRUE)
    check_char(s$inference_name, allow_NULL = TRUE)
    check_char(s$inference_BF_name, allow_NULL = TRUE)
    check_char(s$samples, allow_NULL = TRUE)
    check_char(s$samples_name, allow_NULL = TRUE)
    check_char(s$samples_units, allow_NULL = TRUE)
    check_bool(s$samples_conditional, allow_NULL = TRUE)
  })
  check_list(inference, "inference", check_names = sapply(specification, function(s) s$inference), all_objects = TRUE, allow_other = TRUE)
  check_list(samples, "samples", check_names = unlist(sapply(specification, function(s) s$samples)), all_objects = TRUE, allow_other = TRUE)
  check_char(method, allow_NULL = TRUE)


  output <- ""

  for(i in seq_along(specification)){
    output <- paste0(output, .interpret.specification(inference, samples, specification[[i]], method), if(i != length(specification)) " ")
  }

  return(output)
}

.interpret.specification <- function(inference, samples, specification, method){

  temp_BF <- inference[[specification[["inference"]]]][["BF"]]
  text_BF <- .interpret.BF(temp_BF, if(!is.null(specification[["inference_name"]])) specification[["inference_name"]] else specification[["inference"]],
                           specification[["inference_BF_name"]])

  if(is.null(specification[["samples"]])){
    return(paste0(method, " found ", text_BF, "."))
  }

  temp_par <- samples[[specification[["samples"]]]]
  text_par <- .interpret.par(temp_par, if(!is.null(specification[["samples_name"]])) specification[["samples_name"]] else specification[["samples"]],
                             specification[["samples_units"]], specification[["samples_conditional"]])

  return(paste0(method, " found ", text_BF, ", ", text_par, "."))
}
.interpret.BF            <- function(BF, name, BF_name){

  if(abs(log(BF)) > log(10)){
    text <- "strong evidence"
  }else if(abs(log(BF)) > log(3)){
    text <- "moderate evidence"
  }else{
    text <- "weak evidence"
  }

  if(BF > 1){
    text <- paste0(text, " in favor of the ", name)
  }else{
    text <- paste0(text, " against the ", name)
  }

  BF <- format(round(BF, if(BF < 1) 3 else 2), nsmall = if(BF < 1) 3 else 2)

  if(is.null(BF_name)){
    text <- paste0(text, ", BF = ", BF)
  }else{
    text <- paste0(text, ", ", BF_name, " = ", BF)
  }

  return(text)
}
.interpret.par           <- function(samples, name, unit, conditional){

  est <- mean(samples)
  CI  <- unname(stats::quantile(samples, probs = c(0.025, 0.975)))

  est <- format(round(est, 3), nsmall = 3)
  CI  <- format(round(CI, 3),  nsmall = 3)

  text <- "with"

  if(is.null(conditional) || !conditional){
    text <- paste0(text, " model-averaged")
  }else{
    text <- paste0(text, " conditional")
  }

  if(!is.null(unit)){
    text <- paste0(text, " mean estimate ", name, " = ", est, " ", unit, ", 95% CI [", CI[1], ", ", CI[2], "]")
  }else{
    text <- paste0(text, " mean estimate ", name, " = ", est, ", 95% CI [", CI[1], ", ", CI[2], "]")
  }

  return(text)
}
