#' @title Add 'JAGS' prior
#'
#' @description Adds priors to a 'JAGS' syntax.
#'
#' @param syntax JAGS model syntax
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#'
#' @details Syntax containing moment or inverse-moment priors uses the
#' BayesTools JAGS module. \code{JAGS_fit()} loads this module automatically;
#' users running syntax returned by \code{JAGS_add_priors()} directly in JAGS
#' need to load the module first.
#'
#' @return \code{JAGS_add_priors} returns a JAGS syntax.
#'
#' @export
JAGS_add_priors           <- function(syntax, prior_list){

  syntax <- .check_JAGS_syntax(syntax)

  # return the original syntax in case that no prior was specified
  if(length(prior_list) == 0){
    return(syntax)
  }

  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  .check_prior_list_unique_names(prior_list)

  # identify parts of the syntax
  opening_bracket <- regexpr("{", syntax, fixed = TRUE)[1]
  syntax_start    <- substr(syntax, 1, opening_bracket)
  syntax_end      <- substr(syntax, opening_bracket + 1, nchar(syntax))

  # create the priors relevant syntax
  syntax_priors <- .JAGS_add_priors.fun(prior_list)

  # merge everything back together
  syntax <- paste0(syntax_start, "\n", syntax_priors, "\n", syntax_end)

  return(syntax)
}

.JAGS_add_priors.fun       <- function(prior_list){

  syntax_priors <- ""
  syntax_attributes <- NULL

  for(i in seq_along(prior_list)){

    if(is.prior.weightfunction(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.weightfunction(prior_list[[i]]))

    }else if(is_prior_phacking(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.phacking(prior_list[[i]]))

    }else if(is_prior_bias(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.bias(prior_list[[i]]))

    }else if(is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.PP(prior_list[[i]]))

    }else if(is.prior.spike_and_slab(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.spike_and_slab(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.mixture(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.mixture(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.factor(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.factor(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.vector(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.vector(prior_list[[i]], names(prior_list)[i]))

    }else if(is.prior.simple(prior_list[[i]])){

      syntax_priors <- paste(syntax_priors, .JAGS_prior.simple(prior_list[[i]], names(prior_list)[i]))

    }
  }

  return(syntax_priors)
}
.JAGS_prior.simple         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.simple(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  # parse expressions to test
  if(.is_prior_expression(prior)){
    prior <- .prior_expression_to_character(prior)
  }

  # distribution
  syntax <- switch(
    prior[["distribution"]],
    "point"     = paste0(parameter_name," = ",prior$parameter[["location"]]),
    "normal"    = paste0(parameter_name," ~ dnorm(",prior$parameter[["mean"]],",", .JAGS_parameter_to_precision(prior$parameter[["sd"]]),")"),
    "lognormal" = paste0(parameter_name," ~ dlnorm(",prior$parameter[["meanlog"]],",", .JAGS_parameter_to_precision(prior$parameter[["sdlog"]]),")"),
    "t"         = paste0(parameter_name," ~ dt(",prior$parameter[["location"]],",", .JAGS_parameter_to_precision(prior$parameter[["scale"]]),",", prior$parameter[["df"]],")"),
    "gamma"     = paste0(parameter_name," ~ dgamma(",prior$parameter[["shape"]],",",prior$parameter[["rate"]],")"),
    "invgamma"  = paste0(parameter_name," ~ dbt_invgamma(",prior$parameter[["shape"]],",",prior$parameter[["scale"]],")"),
    "exp"       = paste0(parameter_name," ~ dexp(",prior$parameter[["rate"]],")"),
    "beta"      = paste0(parameter_name," ~ dbeta(",prior$parameter[["alpha"]],",",prior$parameter[["beta"]],")"),
    "bernoulli" = paste0(parameter_name," ~ dbern(",prior$parameter[["probability"]],")"),
    "uniform"   = paste0(parameter_name," ~ dunif(",prior$parameter[["a"]],",",prior$parameter[["b"]],")"),
    "moment"    = paste0(parameter_name," ~ dbt_moment(",prior$parameter[["location"]],",",prior$parameter[["tau"]],",",prior$parameter[["order"]],")"),
    "invmoment" = paste0(parameter_name," ~ dbt_invmoment(",prior$parameter[["location"]],",",prior$parameter[["tau"]],",",prior$parameter[["order"]],",",prior$parameter[["df"]],")")
  )

  # add truncation
  if(!.is_prior_default_range(prior)){
    syntax <- paste0(syntax, "T(",
                     ifelse(is.infinite(prior$truncation[["lower"]]),"",prior$truncation[["lower"]]),
                     ",",
                     ifelse(is.infinite(prior$truncation[["upper"]]),"",prior$truncation[["upper"]]),
                     ")")
  }

  # finish the line
  syntax <- paste0(syntax, "\n")

  return(syntax)
}
.JAGS_prior_dirichlet_eta_name <- function(parameter_name){
  paste0("prior_par_eta_", parameter_name)
}
.JAGS_prior.vector         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.vector(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(prior$parameters[["K"]], "K", lower = 1)
  if(prior[["distribution"]] != "mpoint")
    .check_vector_truncation_unsupported(prior$truncation)

  # parse expressions to test
  if(.is_prior_expression(prior)){
    prior <- .prior_expression_to_character(prior)
  }

  if(prior[["distribution"]] == "dirichlet"){
    eta_name <- .JAGS_prior_dirichlet_eta_name(parameter_name)
    syntax <- paste0(vapply(seq_len(prior$parameters[["K"]]), function(i){
      paste0(eta_name, "[", i, "] ~ dgamma(", prior$parameters[["alpha"]][i], ", 1)\n")
    }, character(1)), collapse = "")
    for(i in seq_len(prior$parameters[["K"]])){
      syntax <- paste0(
        syntax,
        parameter_name, "[", i, "] <- ", eta_name, "[", i, "] / sum(", eta_name, "[1:", prior$parameters[["K"]], "])\n"
      )
    }

  }else if(prior[["distribution"]] %in% c("mnormal", "mt")){
    # create the location/means vector the sigma matrix

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]]
    )
    par2 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["sd"]],
      "mt"      = prior$parameter[["scale"]]
    )

    par2 <- .JAGS_parameter_to_precision(par2)

    # TODO: beautify this code by specific JAGS distributions?
    if(prior[["distribution"]] == "mt"){
      # using the chisq * covariance parametrization since the mt fails with 1 df
      # (using a common df parameter as in Rouder et al. 2012)
      syntax <- paste0("prior_par1_", parameter_name, " = rep(0,", prior$parameter[["K"]], ")\n")
      syntax <- paste0(syntax, "prior_par_s_", parameter_name, " ~ dgamma(", prior$parameter[["df"]]/2, ", ", prior$parameter[["df"]]/2,")\n")
      syntax <- paste0(
        syntax,
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  prior_par2_", parameter_name, "[i,i] <- ", par2, "\n",
        "  for(j in 1:(i-1)){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "  for (j in (i+1):", prior$parameters[["K"]], "){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "}\n",
        "prior_par_z_", parameter_name, " ~ dmnorm(prior_par1_", parameter_name, ",prior_par2_", parameter_name, ")\n",
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  ", parameter_name, "[i] <- prior_par_z_", parameter_name, "[i]/sqrt(prior_par_s_", parameter_name, ") + ", par1, " \n",
        "}\n")
    }else if(prior[["distribution"]] == "mnormal"){
      syntax <- paste0("prior_par1_", parameter_name, " = rep(", par1, ",", prior$parameter[["K"]], ")\n")
      syntax <- paste0(
        syntax,
        "for(i in 1:", prior$parameters[["K"]], "){\n",
        "  prior_par2_", parameter_name, "[i,i] <- ", par2, "\n",
        "  for(j in 1:(i-1)){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "  for (j in (i+1):", prior$parameters[["K"]], "){\n",
        "    prior_par2_", parameter_name, "[i,j] <- 0\n",
        "  }\n",
        "}\n")
      syntax <- paste0(syntax, parameter_name," ~ dmnorm(prior_par1_", parameter_name, ",prior_par2_", parameter_name, ")\n")
    }

  }else if(prior[["distribution"]] == "mpoint"){

    syntax <- paste0(
      "for(i in 1:", prior$parameters[["K"]], "){\n",
      "  ", parameter_name, "[i] = ", prior$parameter[["location"]], " \n",
      "}\n")

  }


  return(syntax)
}
.JAGS_prior.factor         <- function(prior, parameter_name){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.factor(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")
  check_int(.get_prior_factor_levels(prior), "levels", lower = 1)

  if(is.prior.treatment(prior) | is.prior.independent(prior)){

    syntax <- paste0(
      "for(i in 1:", .get_prior_factor_levels(prior), "){\n",
      "  ", .JAGS_prior.simple(prior, paste0(parameter_name, "[i]")),
      "}\n")

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    prior$parameters[["K"]] <- .get_prior_factor_levels(prior)

    syntax <- .JAGS_prior.vector(prior, parameter_name)

  }

  return(syntax)
}
.JAGS_prior.PP             <- function(prior){

  .check_prior(prior, allow_expressions = TRUE)
  if(!is.prior.PET(prior) & !is.prior.PEESE(prior))
    stop("improper prior provided")

  if(is.prior.PET(prior)){
    syntax <- .JAGS_prior.simple(prior, "PET")
  }else if(is.prior.PEESE(prior)){
    syntax <- .JAGS_prior.simple(prior, "PEESE")
  }

  return(syntax)
}
.JAGS_prior.weightfunction <- function(prior){

  .check_prior(prior)
  if(!is.prior.weightfunction(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}
.JAGS_prior.phacking      <- function(prior){

  .check_prior(prior)
  if(!is_prior_phacking(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}
.JAGS_prior.bias          <- function(prior){

  .check_prior(prior)
  if(!is_prior_bias(prior))
    stop("improper prior provided")

  spec <- selection_backend_spec(prior)
  return(.JAGS_selection_backend_syntax(spec))
}

.JAGS_selection_backend_syntax <- function(spec){

  code <- c(spec$prior_code, spec$transform_code)
  code <- code[nzchar(code)]
  if(length(code) == 0L){
    return("")
  }
  code <- sub("[\r\n]+$", "", code)

  paste0(paste0(code, collapse = "\n"), "\n")
}

.JAGS_weightfunction_component_syntax <- function(prior, component_id = NULL, global_cuts = NULL, force_one_sided = FALSE){

  J <- .weightfunction_n_bins(prior)
  syntax <- character()

  expansion     <- .weightfunction_mapping_expansion(prior, force_one_sided)
  all_cuts      <- if(is.null(global_cuts)) expansion$cuts else global_cuts
  needs_mapping <- !identical(all_cuts, .weightfunction_local_cuts(prior)) ||
    !identical(expansion$index, seq_len(J))

  omega_local <- if(is.null(component_id) && !needs_mapping) "omega" else if(is.null(component_id)) "omega_local" else paste0("omega_local_component_", component_id)
  omega_target <- if(is.null(component_id)) "omega" else paste0("omega_component_", component_id)

  if(prior$weights$type == "cumulative"){
    eta_name <- if(is.null(component_id)) "eta" else paste0("eta_component_", component_id)
    std_eta_name <- if(is.null(component_id)) "std_eta" else paste0("std_eta_component_", component_id)

    for(i in seq_len(J)){
      syntax <- paste0(syntax, eta_name, "[", i, "] ~ dgamma(", prior$weights$alpha[i], ", 1)\n")
    }
    syntax <- paste0(syntax,
      "for(j in 1:", J, "){\n",
      "  ", std_eta_name, "[j] <- ", eta_name, "[j] / sum(", eta_name, ")\n",
      "}\n",
      omega_local, "[1] <- 1\n"
    )
    if(J > 1L){
      syntax <- paste0(syntax,
        "for(j in 2:", J, "){\n",
        "  ", omega_local, "[j] <- sum(", std_eta_name, "[j:", J, "])\n",
        "}\n"
      )
    }

  }else if(prior$weights$type == "fixed"){
    for(i in seq_len(J)){
      syntax <- paste0(syntax, omega_local, "[", i, "] <- ", prior$weights$omega[i], "\n")
    }

  }else if(prior$weights$type == "independent"){
    syntax <- paste0(syntax, omega_local, "[1] <- 1\n")
    if(J > 1L){
      for(i in 2:J){
        if(prior$weights$scale == "omega"){
          syntax <- paste0(syntax, .JAGS_prior.simple(prior$weights$prior, paste0(omega_local, "[", i, "]")))
        }else if(prior$weights$scale == "log_omega"){
          log_omega_name <- if(is.null(component_id)) "log_omega" else paste0("log_omega_component_", component_id)
          syntax <- paste0(
            syntax,
            .JAGS_prior.simple(prior$weights$prior, paste0(log_omega_name, "[", i, "]")),
            omega_local, "[", i, "] <- exp(", log_omega_name, "[", i, "])\n"
          )
        }
      }
    }
  }

  if(!is.null(component_id) || needs_mapping){
    global_bin_indices <- .weightfunction_global_bin_indices(all_cuts, expansion)
    for(i in seq_len(length(all_cuts) - 1L)){
      ind <- global_bin_indices[i]
      syntax <- paste0(syntax, omega_target, "[", i, "] <- ", omega_local, "[", expansion$index[ind], "]\n")
    }
  }

  syntax
}
.JAGS_weightfunction_none_component_syntax <- function(component_id, n_bins){

  syntax <- character()
  omega_target <- if(is.null(component_id)) "omega" else paste0("omega_component_", component_id)
  for(i in seq_len(n_bins)){
    syntax <- paste0(syntax, omega_target, "[", i, "] <- 1\n")
  }
  syntax
}
.JAGS_prior.spike_and_slab <- function(prior, parameter_name){

  .check_prior(prior)
  if(!is.prior.spike_and_slab(prior))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  prior_variable_list  <- list(.get_spike_and_slab_variable(prior))
  prior_inclusion_list <- list(.get_spike_and_slab_inclusion(prior))
  names(prior_variable_list)  <- paste0(parameter_name, "_variable")
  names(prior_inclusion_list) <- paste0(parameter_name, "_inclusion")

  syntax <- paste0(
    .JAGS_add_priors.fun(prior_variable_list),
    .JAGS_add_priors.fun(prior_inclusion_list),
    parameter_name, "_indicator ~ dbern(",   paste0(parameter_name, "_inclusion"), ")\n",
    parameter_name, " = ",  parameter_name, "_variable * ", parameter_name, "_indicator\n"
  )

  return(syntax)
}
.JAGS_prior.mixture        <- function(prior_list, parameter_name){

  .check_prior_list(prior_list, allow_expressions = TRUE)
  if(!is.prior.mixture(prior_list))
    stop("improper prior provided")
  check_char(parameter_name, "parameter_name")

  if(inherits(prior_list, "prior.bias_mixture")){

    # dispatch between publication bias prior mixture and a standard prior mixture
    is_PET            <- sapply(prior_list, is.prior.PET)
    is_PEESE          <- sapply(prior_list, is.prior.PEESE)
    is_weightfunction <- sapply(prior_list, is.prior.weightfunction)
    is_phacking       <- sapply(prior_list, is_prior_phacking)
    is_bias           <- sapply(prior_list, is_prior_bias)
    is_none           <- sapply(prior_list, is.prior.none)
    branch_info       <- lapply(prior_list, .selection_branch_info)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

    # if any prior is bias related, the whole component must be dispatching publication bias
    if(any(!(is_PET | is_PEESE | is_weightfunction | is_phacking | is_bias | is_none)))
      stop("Mixture of publication bias and standard priors is not supported.")

    prior_weights <- attr(prior_list, "prior_weights")
    if(any(has_selection) || any(has_phacking)){
      spec <- selection_backend_spec(prior_list)
      syntax <- .JAGS_selection_backend_syntax(spec)
    }else{
      syntax <- paste0(" bias_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))\n")
    }

    if(any(is_PET)){
      if(sum(is_PET) > 1) stop("Only one PET style publication bias adjustment is allowed.")

      named_prior_PET <- prior_list[[which(is_PET)]]
      class(named_prior_PET) <- class(named_prior_PET)[!class(named_prior_PET) %in% "prior.PET"]
      named_prior_PET <- list("PET_1" = named_prior_PET)

      syntax <- paste0(
        syntax,
        .JAGS_add_priors.fun(named_prior_PET),
        " PET <- PET_1 * equals(bias_indicator, ", which(is_PET), ")\n"
      )
    }
    if(any(is_PEESE)){
      if(sum(is_PEESE) > 1) stop("Only one PEESE style publication bias adjustment is allowed.")

      named_prior_PEESE <- prior_list[[which(is_PEESE)]]
      class(named_prior_PEESE) <- class(named_prior_PEESE)[!class(named_prior_PEESE) %in% "prior.PEESE"]
      named_prior_PEESE <- list("PEESE_1" = named_prior_PEESE)

      syntax <- paste0(
        syntax,
        .JAGS_add_priors.fun(named_prior_PEESE),
        " PEESE <- PEESE_1 * equals(bias_indicator, ", which(is_PEESE), ")\n"
      )
    }

  }else{

    prior_weights    <- attr(prior_list, "prior_weights")
    prior_components <- as.list(prior_list)
    class(prior_components) <- "list"
    names(prior_components) <- paste0(parameter_name, "_component_", seq_along(prior_components))

    syntax <- paste0(
      " ", parameter_name, "_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))\n",
      sapply(.JAGS_add_priors.fun(prior_components), paste, collapse = "\n"),
      " ", parameter_name, " = ",  paste0(names(prior_components), " * ", paste0("(", parameter_name, "_indicator == ", seq_along(prior_components), ")"), collapse = " + "), "\n"
    )
  }

  return(syntax)
}


.add_JAGS_vector   <- function(name, vector){

  if(!is.vector(vector))
    stop("vector must be a vector")
  check_char(name, "name")

  syntax <- paste0(" ", name, " = c(", paste0(vector, collapse = ", "), ")\n")

  return(syntax)
}
.add_JAGS_matrix   <- function(name, matrix){

  if(!is.matrix(matrix))
    stop("matrix must be a matrix")
  check_char(name, "name")

  syntax <- ""

  # this unfortunatelly cannot be defined on row/column basis
  # I tried simplifying this before but only possible initialization is elementwise
  for(i in 1:nrow(matrix)){
   syntax <- paste0(
     syntax, " ",
     paste0(name,"[", i, ",", seq_len(ncol(matrix)), "] = ", matrix[i,], collapse = "; "), "\n"
   )
  }

  return(syntax)
}
.check_JAGS_syntax <- function(syntax){

  check_char(syntax, "syntax", allow_NULL = TRUE)
  if(is.null(syntax)){
    syntax <- "model{}"
  }
  if(!grepl("model", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("{", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  if(!grepl("}", syntax, fixed = TRUE))
    stop("syntax must be a JAGS model syntax")
  return(syntax)
}
.JAGS_parameter_to_precision <- function(parameter){

  if(is.character(parameter)){
    parameter <- paste0("1/pow(", parameter, ", 2)")
  }else{
    parameter <- 1/parameter^2
  }

  return(parameter)
}
