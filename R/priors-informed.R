#' @title Creates an informed prior distribution based on research
#'
#' @description \code{prior_informed} creates an informed prior distribution based on past
#' research. The prior can be visualized by the \code{plot} function.
#'
#' @param name name of the prior distribution. There are many options based on prior psychological
#' or medical research.
#' For psychology, the possible options are
#' \describe{
#'   \item{\code{"van Erp"}}{for an informed prior distribution for the heterogeneity parameter tau
#'   of meta-analytic effect size estimates based on standardized mean differences
#'   \insertCite{erp2017estimates}{BayesTools},}
#'   \item{\code{"Oosterwijk"}}{for an informed prior distribution for the effect sizes expected in
#'   social psychology based on prior elicitation with dr. Oosterwijk
#'   \insertCite{gronau2017bayesian}{BayesTools}.}
#' }
#' For medicine, the possible options are based on \insertCite{bartos2021bayesian;textual}{BayesTools}
#' who developed empirical prior distributions for the effect size and heterogeneity parameters of the
#' continuous standardized outcomes based on the Cochrane database of systematic reviews.
#' Use \code{"Cochrane"} for a prior distribution based on the whole database or call
#' \code{print(prior_informed_medicine_names)} to inspect the names of
#' all 46 subfields and set the appropriate \code{parameter} and \code{type}.
#' @param parameter parameter name describing what prior distribution is supposed to be produced in cases
#' where the \code{name} corresponds to multiple prior distributions. Relevant only for the empirical medical
#' prior distributions.
#' @param type prior type describing what prior distribution is supposed to be produced in cases
#' where the \code{name} and \code{parameter} correspond to multiple prior distributions. Relevant only for
#' the empirical medical prior distributions.
#'
#' @examples
#' # prior distribution representing expected effect sizes in social psychology
#' # based on prior elicitation with dr. Oosterwijk
#' p1 <- prior_informed("Oosterwijk")
#'
#' # the prior distribution can be visualized using the plot function
#' # (see ?plot.prior for all options)
#' plot(p1)
#'
#'
#' # empirical prior distribution for the standardized mean differences from the oral health
#' # medical subfield based on meta-analytic effect size estimates from the
#' # Cochrane database of systematic reviews
#' p2 <- prior_informed("Oral Health", parameter = "effect", type = "smd")
#' print(p2)
#'
#' @return \code{prior_informed} returns an object of class 'prior'.
#'
#' @references
#' \insertAllCited{}
#' @seealso [prior()], [prior_informed_medicine_names]
#' @export
prior_informed <- function(name, parameter = NULL, type = "smd"){

  check_char(name,      "name")
  check_char(parameter, "parameter", allow_NULL = TRUE)
  check_char(type,      "type",      allow_NULL = TRUE)

  name <- .prior_clean_input_name(name)

  if(name %in% .prior_clean_input_name(prior_informed_medicine_names)){
    # check for implemented metrics
    check_char(type, "type", allow_values = c("smd"))

    p <- switch(
      parameter,
      "effect"        = .prior_informed.medicine_effect(name, type),
      "heterogeneity" = .prior_informed.medicine_heterogeneity(name, type),
      stop("unknown 'parameter' argument for an informed prior distribution from medicine.")
    )

    return(p)

  }else if(name %in% .prior_clean_input_name(.prior_informed_psychology_names)){

    p <- .prior_informed.psychology(name)

    return(p)

  }else{
    stop("unknown prior 'name'. See '?prior_informed' for help.")
  }
}


.prior_informed.medicine_effect        <- function(name, type){

  if(type == "smd"){
    # medical priors for continuous outcomes based on Bartoš et al. 2021
    # (can't use swucg if the argument is to be modified)
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.38, df = 5))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.38, df = 6))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.55, df = 4))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.37, df = 5))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.51, df = 5))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.55, df = 5))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.47, df = 5))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.18, df = 5))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.39, df = 5))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.39, df = 5))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.40, df = 6))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.38, df = 5))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.60, df = 4))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.48, df = 3))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.33, df = 6))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.59, df = 2))
    }else if(name == .prior_clean_input_name("Inflammatory Bowel Disease")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.35, df = 5))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.54, df = 5))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.49, df = 5))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.45, df = 6))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.51, df = 5))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Public Health")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.29, df = 4))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.48, df = 5))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.48, df = 5))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.44, df = 4))
    }else if(name == .prior_clean_input_name("Upper GI and Pancreatic Diseases")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.44, df = 5))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.46, df = 5))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.56, df = 5))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution = "t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else{
    stop(paste0("Type '", type, "' is not recognized."))
  }
}
.prior_informed.medicine_heterogeneity <- function(name, type){

  if(type == "smd"){
    # medical priors for continuous outcomes based on Bartoš et al. 2021
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.73, scale = 0.46))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 2.02, scale = 0.28))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.62, scale = 0.64))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.75, scale = 0.57))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.52, scale = 0.28))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.64, scale = 0.56))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.62, scale = 0.45))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.56, scale = 0.14))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.70, scale = 0.45))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.71, scale = 0.44))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.43, scale = 0.12))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.89, scale = 0.28))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.71, scale = 0.35))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.62, scale = 0.29))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.85, scale = 0.48))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.86, scale = 0.41))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.67, scale = 0.46))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.78, scale = 0.46))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.83, scale = 0.47))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.56, scale = 0.58))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.73, scale = 0.44))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 2.01, scale = 0.38))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.64, scale = 0.36))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.28, scale = 0.44))
    }else if(name == .prior_clean_input_name("Inflammatory Bowel Disease")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.76, scale = 0.39))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.80, scale = 0.34))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.72, scale = 0.53))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.71, scale = 0.37))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.72, scale = 0.51))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.88, scale = 0.33))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.87, scale = 0.38))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.68, scale = 0.38))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.79, scale = 0.28))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.69, scale = 0.42))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.86, scale = 0.32))
    }else if(name == .prior_clean_input_name("Public Health")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.76, scale = 0.23))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.60, scale = 0.27))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.70, scale = 0.59))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.64, scale = 0.51))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.71, scale = 0.40))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.73, scale = 0.42))
    }else if(name == .prior_clean_input_name("Upper GI and Pancreatic Diseases")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.76, scale = 0.38))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.73, scale = 0.45))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.66, scale = 0.50))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.76, scale = 0.39))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.54, scale = 0.41))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution = "invgamma", parameters = list(shape = 1.71, scale = 0.40))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else{
    stop(paste0("Type '", type, "' is not recognized."))
  }
}
.prior_informed.psychology             <- function(name){

  if(name == .prior_clean_input_name("van Erp")){
    p <- prior(distribution = "invgamma", parameters = list(shape = 1, scale = 0.15))
  }else if(name == .prior_clean_input_name("Oosterwijk")){
    p <- prior(distribution = "t", parameters = list(location = 0.35, scale = 0.102, df = 3), truncation = list(lower = 0))
  }else{
    stop("unknown prior 'name' for an informed prior distribution from psychology.")
  }

  return(p)
}

#' @title Names of medical subfields from the Cochrane database of systematic reviews
#'
#' @description Contain names identifying the individual subfields from the Cochrane database
#' of systematic reviews. The individual elements correspond to valid \code{name} arguments
#' for the [prior_informed()] function.
#'
#' @examples
#' print(prior_informed_medicine_names)
#'
#' @return returns a character vector with names of medical subfields from Cochrane database of
#' systematic reviews.
#'
#' @seealso [prior_informed()]
#' @export
prior_informed_medicine_names <- c(
  "Acute Respiratory Infections",
  "Airways",
  "Anaesthesia",
  "Back and Neck",
  "Bone, Joint and Muscle Trauma",
  "Colorectal",
  "Common Mental Disorders",
  "Consumers and Communication",
  "Cystic Fibrosis and Genetic Disorders",
  "Dementia and Cognitive Improvement",
  "Developmental, Psychosocial and Learning Problems",
  "Drugs and Alcohol",
  "Effective Practice and Organisation of Care",
  "Emergency and Critical Care",
  "ENT",
  "Eyes and Vision",
  "Gynaecological, Neuro-oncology and Orphan Cancer",
  "Gynaecology and Fertility",
  "Heart",
  "Hepato-Biliary",
  "HIV/AIDS",
  "Hypertension",
  "Incontinence",
  "Infectious Diseases",
  "Inflammatory Bowel Disease",
  "Injuries",
  "Kidney and Transplant",
  "Metabolic and Endocrine Disorders",
  "Methodology",
  "Movement Disorders",
  "Musculoskeletal",
  "Neonatal",
  "Oral Health",
  "Pain, Palliative and Supportive Care",
  "Pregnancy and Childbirth",
  "Public Health",
  "Schizophrenia",
  "Sexually Transmitted Infections",
  "Skin",
  "Stroke",
  "Tobacco Addiction",
  "Upper GI and Pancreatic Diseases",
  "Urology",
  "Vascular",
  "Work",
  "Wounds",
  "Cochrane"
)

.prior_informed_psychology_names <- c(
  "van Erp",
  "Oosterwijk"
)
