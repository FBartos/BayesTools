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
#' and  \insertCite{bartos2023empirical;textual}{BayesTools}
#' who developed empirical prior distributions for the effect size and heterogeneity parameters of the
#' continuous outcomes (standardized mean differences), dichotomous outcomes (logOR, logRR, and risk differences),
#' and time to event outcomes (logHR) based on the Cochrane database of systematic reviews.
#' Use \code{"Cochrane"} for a prior distribution based on the whole database or call
#' \code{print(prior_informed_medicine_names)} to inspect the names of
#' all 46 subfields and set the appropriate \code{parameter} and \code{type}.
#' @param parameter parameter name describing what prior distribution is supposed to be produced in cases
#' where the \code{name} corresponds to multiple prior distributions. Relevant only for the empirical medical
#' prior distributions.
#' @param type prior type describing what prior distribution is supposed to be produced in cases
#' where the \code{name} and \code{parameter} correspond to multiple prior distributions. Relevant only for
#' the empirical medical prior distributions with the following options
#' \describe{
#'   \item{\code{"smd"}}{for standardized mean differences}
#'   \item{\code{"logOR"}}{for log odds ratios}
#'   \item{\code{"logRR"}}{for log risk ratios}
#'   \item{\code{"RD"}}{for risk differences}
#'   \item{\code{"logHR"}}{for hazard ratios}
#' }
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
#' p2 <- prior_informed("Oral Health", parameter ="effect", type ="smd")
#' print(p2)
#'
#' @return \code{prior_informed} returns an object of class 'prior'.
#'
#' @references
#' \insertAllCited{}
#' @seealso [prior()], [prior_informed_medicine_names]
#' @export
prior_informed <- function(name, parameter = NULL, type ="smd"){

  check_char(name,"name")
  check_char(parameter,"parameter", allow_NULL = TRUE)
  check_char(type,"type",      allow_NULL = TRUE)

  name <- .prior_clean_input_name(name)

  if(name %in% .prior_clean_input_name(prior_informed_medicine_names)){
    # check for implemented metrics
    type <- tolower(type)
    check_char(type,"type", allow_values = c("smd","logor","logrr","loghr","rd"))

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

  if(type =="smd"){

    # medical priors for continuous outcomes based on Bartoš et al. 2021
    # (can't use switch if the argument is to be modified)
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 5))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 6))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.55, df = 4))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.37, df = 5))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.51, df = 5))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.55, df = 5))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.47, df = 5))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.18, df = 5))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 5))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 5))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 6))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 5))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.60, df = 4))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.48, df = 3))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 6))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.59, df = 2))
    }else if(name == .prior_clean_input_name("Inflammatory Bowel Disease")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 5))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.35, df = 5))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.54, df = 5))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.49, df = 5))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.45, df = 6))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.51, df = 5))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Public Health")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 5))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.29, df = 4))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.48, df = 5))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.48, df = 5))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 4))
    }else if(name == .prior_clean_input_name("Upper GI and Pancreatic Diseases")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 5))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.46, df = 5))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.42, df = 5))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.56, df = 5))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="logor"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.48, df = 3))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.37, df = 2))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.79, df = 6))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.62, df = 4))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 2))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 3))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.47, df = 4))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.65, df = 4))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.54, df = 4))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.48, df = 4))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 3))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.49, df = 4))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.83, df = 5))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 4))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.51, df = 4))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.35, df = 3))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.81, df = 4))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.88, df = 6))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.77, df = 5))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.46, df = 5))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.63, df = 5))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.53, df = 4))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 2))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.57, df = 4))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.20, df = 2))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.95, df = 4))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 3))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.32, df = 4))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.28, df = 5))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.75, df = 3))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.66, df = 3))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.60, df = 4))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.53, df = 4))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.61, df = 5))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.29, df = 2))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.60, df = 5))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.73, df = 5))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.76, df = 4))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.59, df = 4))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.29, df = 3))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.70, df = 5))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 1.13, df = 4))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 1.00, df = 7))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 2))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.58, df = 4))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.61, df = 2))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.81, df = 2))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.22, df = 2))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.49, df = 5))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.82, df = 5))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.68, df = 5))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.59, df = 4))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.60, df = 4))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.58, df = 4))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="logrr"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.27, df = 3))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.25, df = 3))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.45, df = 5))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.34, df = 3))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.26, df = 2))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.21, df = 3))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.28, df = 3))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 4))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.29, df = 3))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.13, df = 2))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.21, df = 2))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.29, df = 3))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 3))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.22, df = 4))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.30, df = 4))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.19, df = 2))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 4))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.58, df = 5))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.34, df = 4))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.34, df = 4))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.40, df = 4))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.35, df = 4))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.25, df = 2))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.31, df = 3))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.17, df = 3))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 3))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.23, df = 3))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.21, df = 3))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.25, df = 4))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 3))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.41, df = 3))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.27, df = 4))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.35, df = 4))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.38, df = 4))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.13, df = 1))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.27, df = 5))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.46, df = 4))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 3))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 3))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.18, df = 3))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.41, df = 4))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.55, df = 4))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.59, df = 5))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.24, df = 2))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.33, df = 3))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.39, df = 3))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.41, df = 1))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.13, df = 2))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.36, df = 5))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.44, df = 3))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.43, df = 5))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.23, df = 3))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.37, df = 4))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.32, df = 3))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

  }else if(type =="rd"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.012, df = 1))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.005, df = 1))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.020, df = 1))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.056, df = 2))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.014, df = 1))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.049, df = 2))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.021, df = 1))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.010, df = 1))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.051, df = 1))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.063, df = 2))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.015, df = 1))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.017, df = 1))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.082, df = 2))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.025, df = 1))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.056, df = 2))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.019, df = 1))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.042, df = 1))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.047, df = 2))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.041, df = 1))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.014, df = 1))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.047, df = 2))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.024, df = 1))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.021, df = 1))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.047, df = 1))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.005, df = 1))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.083, df = 1))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.014, df = 1))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.018, df = 2))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.006, df = 2))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.060, df = 2))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.019, df = 1))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.010, df = 1))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.028, df = 1))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.057, df = 2))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.003, df = 1))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.097, df = 2))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.055, df = 2))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.067, df = 2))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.010, df = 1))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.022, df = 1))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.034, df = 1))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.030, df = 1))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.111, df = 2))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.012, df = 1))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.028, df = 1))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.035, df = 1))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.033, df = 1))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.014, df = 1))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.042, df = 2))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.041, df = 1))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.012, df = 1))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.062, df = 2))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.012, df = 1))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.034, df = 1))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="loghr"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="t", parameters = list(location = 0, scale = 0.13, df = 2))
    }else{
      stop("There are no subfield specific prior distributions for logHR.")
    }

    return(p)
  }
}
.prior_informed.medicine_heterogeneity <- function(name, type){

  if(type =="smd"){

    # medical priors for continuous outcomes based on Bartoš et al. 2021
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.46))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.02, scale = 0.28))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.62, scale = 0.64))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.75, scale = 0.57))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.52, scale = 0.28))
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.56))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.62, scale = 0.45))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.56, scale = 0.14))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.70, scale = 0.45))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.44))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.43, scale = 0.12))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.89, scale = 0.28))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.35))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.62, scale = 0.29))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.85, scale = 0.48))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.86, scale = 0.41))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.67, scale = 0.46))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.78, scale = 0.46))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.83, scale = 0.47))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.56, scale = 0.58))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.44))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.01, scale = 0.38))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.36))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.28, scale = 0.44))
    }else if(name == .prior_clean_input_name("Inflammatory Bowel Disease")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.39))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.80, scale = 0.34))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.72, scale = 0.53))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.37))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.72, scale = 0.51))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.88, scale = 0.33))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.87, scale = 0.38))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.68, scale = 0.38))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.79, scale = 0.28))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.69, scale = 0.42))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.86, scale = 0.32))
    }else if(name == .prior_clean_input_name("Public Health")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.23))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.60, scale = 0.27))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.70, scale = 0.59))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.51))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.40))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.42))
    }else if(name == .prior_clean_input_name("Upper GI and Pancreatic Diseases")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.38))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.45))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.66, scale = 0.50))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.39))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.54, scale = 0.41))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.40))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="logor"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.67, scale = 0.45))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.35, scale = 0.27))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.12, scale = 0.86))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.84, scale = 0.68))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.01, scale = 0.36))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.59, scale = 0.48))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      stop("no informed subfield prior distribution for 'Childhood Cancer' category.")
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.60))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.63, scale = 0.45))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.37, scale = 0.86))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.82, scale = 0.58))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.28, scale = 0.25))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.83, scale = 0.82))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.59, scale = 0.42))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.90, scale = 0.68))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.46, scale = 0.34))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.71))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.43))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.09, scale = 0.94))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.21, scale = 0.71))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.94, scale = 0.62))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.80, scale = 0.56))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.24, scale = 0.28))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.91, scale = 0.66))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.29))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.83))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.58, scale = 0.40))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.36))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.25, scale = 0.10))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.07, scale = 1.09))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.08, scale = 0.86))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.52, scale = 0.49))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.68, scale = 0.44))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.04, scale = 0.68))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 0.92, scale = 0.11))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.04, scale = 0.50))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.14, scale = 0.64))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.09, scale = 0.71))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.76, scale = 0.59))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.80, scale = 0.42))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.74, scale = 0.49))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.85, scale = 0.70))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.50, scale = 0.40))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.73, scale = 0.47))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.92, scale = 0.69))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.72, scale = 0.47))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.64, scale = 0.49))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.36, scale = 0.21))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.51, scale = 0.63))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.72, scale = 0.50))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.70, scale = 0.45))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.79, scale = 0.73))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.16, scale = 0.86))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.77, scale = 0.55))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="logrr"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.58, scale = 0.25))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.36, scale = 0.15))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.21, scale = 0.21))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.52, scale = 0.22))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 0.97, scale = 0.14))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.49, scale = 0.24))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      stop("no informed subfield prior distribution for 'Childhood Cancer' category.")
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.42, scale = 0.27))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.43, scale = 0.20))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.48, scale = 0.14))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.71, scale = 0.24))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.53, scale = 0.19))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 0.97, scale = 0.11))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.58, scale = 0.20))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.40, scale = 0.24))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.65, scale = 0.25))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.50, scale = 0.24))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.74, scale = 0.27))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.60, scale = 0.26))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.59, scale = 0.31))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.69, scale = 0.29))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.43, scale = 0.25))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.29, scale = 0.16))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.22, scale = 0.11))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.96, scale = 0.24))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.62, scale = 0.25))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.44, scale = 0.19))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.66, scale = 0.23))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.48, scale = 0.11))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.53, scale = 0.33))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.30, scale = 0.27))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.46, scale = 0.19))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.51, scale = 0.23))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.49, scale = 0.26))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.55, scale = 0.19))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.74, scale = 0.23))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.54, scale = 0.24))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.58, scale = 0.35))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.54, scale = 0.31))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.89, scale = 0.30))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.44, scale = 0.15))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.50, scale = 0.37))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.98, scale = 0.39))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.40, scale = 0.21))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.47, scale = 0.26))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.43, scale = 0.21))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.49, scale = 0.23))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.52, scale = 0.15))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.02, scale = 0.36))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.37, scale = 0.18))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.26, scale = 0.19))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.50, scale = 0.26))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.63, scale = 0.38))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 1.51, scale = 0.23))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

  }else if(type =="rd"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Acute Respiratory Infections")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.102), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Airways")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.099), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Anaesthesia")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.108), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Back and Neck")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.102), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Bone, Joint and Muscle Trauma")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.119), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Breast Cancer")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.118), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Childhood Cancer")){
      stop("no informed subfield prior distribution for 'Childhood Cancer' category.")
    }else if(name == .prior_clean_input_name("Colorectal")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.096), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Common Mental Disorders")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.108), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Consumers and Communication")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.103), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Cystic Fibrosis and Genetic Disorders")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.092), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Dementia and Cognitive Improvement")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.092), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Developmental, Psychosocial and Learning Problems")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.132), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Drugs and Alcohol")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.101), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Effective Practice and Organisation of Care")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.106), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Emergency and Critical Care")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.068), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("ENT")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.123), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Epilepsy")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.086), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Eyes and Vision")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.112), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Fertility Regulation")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.082), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Gut")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.090), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Gynaecological, Neuro-oncology and Orphan Cancer")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.103), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Gynaecology and Fertility")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.086), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Haematology")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.089), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Heart")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.070), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Heart; Vascular")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.116), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Hepato-Biliary")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.093), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("HIV/AIDS")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.097), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Hypertension")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.059), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Incontinence")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.125), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Infectious Diseases")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.098), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Injuries")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.137), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Kidney and Transplant")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.089), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Lung Cancer")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.116), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Metabolic and Endocrine Disorders")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.091), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Methodology")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.102), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Movement Disorders")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.071), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Multiple Sclerosis and Rare Diseases of the CNS")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.101), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Musculoskeletal")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.109), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Neonatal")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.088), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Neuromuscular")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.075), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Oral Health")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.121), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Pain, Palliative and Supportive Care")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.106), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Pregnancy and Childbirth")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.082), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Schizophrenia")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.102), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Sexually Transmitted Infections")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.117), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Skin")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.123), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Stroke")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.060), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Tobacco Addiction")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.065), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Urology")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.102), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Vascular")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.119), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Work")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.123), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Wounds")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.088), truncation = list(lower = 0))
    }else if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="normal", parameters = list(mean = 0, sd = 0.099), truncation = list(lower = 0))
    }else{
      stop("unknown subfield 'name' argument for an informed prior distribution from medicine.")
    }

    return(p)

  }else if(type =="loghr"){

    # medical priors for dichotomous outcomes based on Bartoš et al. 2023
    if(name == .prior_clean_input_name("Cochrane")){
      p <- prior(distribution ="invgamma", parameters = list(shape = 2.42, scale = 0.30))
    }else{
      stop("There are no subfield specific prior distributions for logHR.")
    }

    return(p)
  }

}
.prior_informed.psychology             <- function(name){

  if(name == .prior_clean_input_name("van Erp")){
    p <- prior(distribution ="invgamma", parameters = list(shape = 1, scale = 0.15))
  }else if(name == .prior_clean_input_name("Oosterwijk")){
    p <- prior(distribution ="t", parameters = list(location = 0.35, scale = 0.102, df = 3), truncation = list(lower = 0))
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
  "Breast Cancer",
  "Childhood Cancer",
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
  "Epilepsy",
  "Eyes and Vision",
  "Fertility Regulation",
  "Gut",
  "Gynaecological, Neuro-oncology and Orphan Cancer",
  "Gynaecology and Fertility",
  "Haematology",
  "Heart",
  "Heart; Vascular",
  "Hepato-Biliary",
  "HIV/AIDS",
  "Hypertension",
  "Incontinence",
  "Infectious Diseases",
  "Inflammatory Bowel Disease",
  "Injuries",
  "Kidney and Transplant",
  "Lung Cancer",
  "Metabolic and Endocrine Disorders",
  "Methodology",
  "Movement Disorders",
  "Multiple Sclerosis and Rare Diseases of the CNS",
  "Musculoskeletal",
  "Neonatal",
  "Neuromuscular",
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
