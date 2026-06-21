#### marginal distribution functions ####
# weightfunctions require just marginals for plotting and etc
# (also, the joint are not implemented in general, since they differ between JAGS/Stan implementations)
#' @rdname prior_functions
mcdf.prior   <- function(x, q, ...){

  prior <- x

  .check_q(as.vector(q))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mcdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mcdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_cdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    p <- do.call(cbind, lapply(components, .prior_weightfunction_component_cdf, q = q))

  }else if(is.prior.vector(prior) && !is.prior.factor(prior)){

    p <- .prior_vector_marginal_cdf(prior, q, lower.tail = TRUE)

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2),
      "mpoint"  = ppoint(q, location = 0)
    )

  }

  return(p)
}
#' @rdname prior_functions
mccdf.prior  <- function(x, q, ...){

  prior <- x

  .check_q(as.vector(q))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mccdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mccdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    p <- .prior_simple_ccdf(prior, q)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    p <- do.call(cbind, lapply(components, .prior_weightfunction_component_ccdf, q = q))

  }else if(is.prior.vector(prior) && !is.prior.factor(prior)){

    p <- .prior_vector_marginal_cdf(prior, q, lower.tail = FALSE)

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    p <- switch(
      prior[["distribution"]],
      "mnormal" = stats::pnorm(q, mean = 0, sd = par2, lower.tail = FALSE),
      "mt"      = extraDistr::plst(q, df = prior$parameters[["df"]], mu = 0, sigma = par2, lower.tail = FALSE),
      "mpoint"  = ppoint(q, location = 0, lower.tail = FALSE),
    )

  }

  return(p)
}
#' @rdname prior_functions
mlpdf.prior  <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(as.vector(x))
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No mlpdf are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No mlpdf are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    log_lik <- .prior_simple_lpdf(prior, x)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    log_lik <- do.call(cbind, lapply(components, .prior_weightfunction_component_lpdf, x = x))

  }else if(is.prior.vector(prior) && !is.prior.factor(prior)){

    log_lik <- .prior_vector_marginal_lpdf(prior, x)

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    log_lik <- switch(
      prior[["distribution"]],
      "mnormal" = stats::dnorm(x, mean = 0, sd = par2, log = TRUE),
      "mt"      = extraDistr::dlst(x, df = prior$parameters[["df"]], mu = 0, sigma = par2, log = TRUE),
      "mpoint"  = dpoint(x, location = 0, log = TRUE),
    )

  }

  return(log_lik)
}
#' @rdname prior_functions
mpdf.prior   <- function(x, y, ...){

  prior <- x
  x     <- y

  .check_x(as.vector(x))
  .check_prior(prior)

  if(is.prior.simple(prior)){
    lik <- .prior_simple_pdf(prior, x)
  }else{
    log_lik <- mlpdf.prior(prior, x)
    lik     <- exp(log_lik)
  }

  return(lik)
}
#' @rdname prior_functions
mquant.prior <- function(x, p, ...){

  prior <- x

  .check_p(p, log.p = FALSE)
  .check_prior(prior)

  if(is.prior.spike_and_slab(prior)){

    stop("No quantile functions are implemented for spike and slab priors.")

  }else if(is.prior.mixture(prior)){

    stop("No quantile functions are implemented for prior mixtures.")

  }else if(is.prior.simple(prior)){

    q <- .prior_simple_quant(prior, p)

  }else if(is.prior.weightfunction(prior)){

    components <- .weightfunction_marginal_components(prior)
    q <- do.call(cbind, lapply(components, .prior_weightfunction_component_quant, p = p))

  }else if(is.prior.vector(prior) && !is.prior.factor(prior)){

    q <- .prior_vector_marginal_quant(prior, p)

  }else if(is.prior.orthonormal(prior) | is.prior.meandif(prior)){

    par1 <- switch(
      prior[["distribution"]],
      "mnormal" = prior$parameter[["mean"]],
      "mt"      = prior$parameter[["location"]],
      "mpoint"  = prior$parameter[["location"]]
    )

    if(length(par1) != 1){
      stop("unsported distribution specification in 'rng' -- non-symmetric")
    }

    if(par1 != 0){
      stop("the orthonormal/meandif prior distribution must be centered")
    }

    if(is.na(prior$parameters[["K"]]) && !is.null(attr(prior, "levels"))){
      prior$parameters[["K"]] <- .get_prior_factor_levels(prior)
    }else if(is.na(prior$parameters[["K"]])){
      prior$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    if(prior[["distribution"]] %in% c("mnormal", "mt")){
      if(is.prior.orthonormal(prior)){
        par2 <- sqrt(sum( (contr.orthonormal(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }else if(is.prior.meandif(prior)){
        par2 <- sqrt(sum( (contr.meandif(1:(prior$parameters[["K"]] + 1))[1,] * switch(
          prior[["distribution"]],
          "mnormal" = prior$parameters[["sd"]],
          "mt"      = prior$parameters[["scale"]]) )^2 ))
      }
    }

    q <- switch(
      prior[["distribution"]],
      "mnormal" = stats::qnorm(p, mean = 0, sd = par2),
      "mt"      = extraDistr::qlst(p, df = prior$parameters[["df"]], mu = 0, sigma = par2),
      "mpoint"  = qpoint(p, location = 0)
    )

  }

  return(q)
}
