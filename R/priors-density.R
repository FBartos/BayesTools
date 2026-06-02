#' @title Prior density
#'
#' @description Computes density of a prior
#' distribution across a range of values.
#'
#' @param x a prior
#' @param x_seq sequence of x coordinates
#' @param x_range vector of length two with
#' lower and upper range for the support
#' (used if \code{x_seq} is unspecified)
#' @param x_range_quant quantile used for
#' automatically obtaining \code{x_range}
#' if both \code{x_range} and \code{x_seq}
#' are unspecified. Defaults to \code{0.005}
#' for all but Cauchy, Student-t, Gamma, inverse-gamma,
#' moment, and inverse-moment distributions that use
#' \code{0.010}.
#' @param n_points number of equally spaced points
#' in the \code{x_range} if \code{x_seq} is unspecified
#' @param n_samples number of samples from the prior
#' distribution if the density cannot be obtained
#' analytically (or if samples are forced with
#' \code{force_samples = TRUE})
#' @param force_samples should prior be sampled instead
#' of obtaining analytic solution whenever possible
#' @param individual should individual densities be returned
#' (e.g., in case of weightfunction)
#' @param transformation transformation to be applied
#' to the prior distribution. Either a character
#' specifying one of the prepared transformations:
#' \describe{
#'   \item{lin}{linear transformation in form of \code{a + b*x}}
#'   \item{tanh}{hyperbolic tangent transformation}
#'   \item{exp}{exponential transformation}
#' }, or a list containing the transformation function \code{fun},
#' inverse transformation function \code{inv}, and derivative of the
#' transformation \code{jac}, evaluated on the original support. See examples
#' for details.
#' @param transformation_arguments a list with named arguments for
#' the \code{transformation}
#' @param transformation_settings boolean indicating whether the
#' settings the \code{x_seq} or \code{x_range} was specified on
#' the transformed support
#' @param truncate_end whether the density should be set to zero
#' for endpoints of truncated distributions
#' @param ... additional arguments
#'
#' @return \code{density.prior} returns an object of class 'density'.
#'
#' @details Sample-based density estimates for continuous priors with finite
#' support use boundary-reflected kernel density estimates. The plotting range
#' controls the evaluation grid, while reflection is based on the prior's true
#' truncation bounds.
#'
#' @importFrom stats density
#' @seealso [prior()]
#' @rdname density.prior
#' @export
density.prior <- function(x,
                          x_seq = NULL, x_range = NULL, x_range_quant = NULL, n_points = 1000,
                          n_samples = 10000, force_samples = FALSE, individual = FALSE,
                          transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE, truncate_end = TRUE, ...){

  # input check
  .check_prior(x, "x")
  check_real(x_seq, "x_seq", check_length = 0, allow_NULL = TRUE)
  check_real(x_range, "x_range", check_length = 2, allow_NULL = TRUE)
  if(!is.null(x_range) && x_range[1] > x_range[2])
    stop("The lower range limit must be lower than the upper range limit.")
  check_real(x_range_quant, "x_range_quant", lower = 0, upper = 1, allow_NULL = TRUE)
  check_int(n_points, "n_points",  lower = 2)
  check_int(n_samples, "n_samples", lower = 1)
  check_bool(force_samples, "force_samples")
  check_bool(individual, "individual")
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)
  check_bool(truncate_end, "truncate_end")


  ### setting the range
  # get plotting range if not specified
  if(is.null(x_range)){
    if(!is.null(x_seq)){
      x_range <- range(x_seq)
    }else{
      if(is_prior_phacking(x) || is_prior_bias(x)){
        .selection_prior_stop_unsupported_generic("density", x)
      }else if(!individual & (is.prior.PET(x) | is.prior.PEESE(x))){
        x_range <- c(0, 1)
      }else if(!individual & is.prior.weightfunction(x)){
        x_range <- c(0, 1)
      }else if(is.prior.spike_and_slab(x)){
        x_range <- range(c(range(.get_spike_and_slab_variable(x), if(is.null(x_range_quant)) .range.prior_quantile_default(.get_spike_and_slab_variable(x)) else x_range_quant), 0))
      }else if(is.prior.discrete(x)){
        x_range <- c(x[["truncation"]][["lower"]], x[["truncation"]][["upper"]])
      }else{
        x_range <- range(x, if(is.null(x_range_quant)) .range.prior_quantile_default(x) else x_range_quant)
      }
    }
  }

  # get the x_seq for plotting
  if(is.null(x_seq)){
    if(is.prior.discrete(x)){
      x_seq <- seq(x_range[1], x_range[2], by = 1)
    }else{
      x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
    }
  }

  # specify it on the transformed range if requested
  if(transformation_settings & !is.null(transformation)){
    x_seq   <- .density.prior_transformation_inv_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_inv_x(x_range, transformation, transformation_arguments)
  }


  # use the corresponding density subfunction
  if(is.prior.weightfunction(x)){
    out <- .density.prior.weightfunction(x, x_seq, x_range, n_points, n_samples, force_samples, individual)
  }else if(is_prior_phacking(x) || is_prior_bias(x)){
    .selection_prior_stop_unsupported_generic("density", x)
  }else if(is.prior.PET(x) | is.prior.PEESE(x)){
    out <- .density.prior.PETPEESE(x, x_seq, x_range, n_points, n_samples, force_samples, individual, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.spike_and_slab(x)){
    out <- .density.prior.spike_and_slab(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.point(x)){
    out <- .density.prior.point(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)
  }else if(is.prior.orthonormal(x) | is.prior.meandif(x)){
    out <- .density.prior.orthonormal_or_meandif(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }else if(is.prior.simple(x)){
    out <- .density.prior.simple(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  }

  if(!is.null(transformation)){
    attr(out, "transformation") <- transformation
  }

  return(out)
}

.density.prior.simple                 <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  boundary_reflection <- FALSE

  # get the samples to estimate density / obtain the density directly
  if(force_samples | .density.prior_need_samples(x)){
    x_sam <- rng(x, n_samples)
    if(is.prior.discrete(x)){
      x_seq <- unique(round(x_seq))
      x_den <- vapply(x_seq, function(x_i) mean(x_sam == x_i), numeric(1))
    }else{
      x_density <- .density_kde_boundary(
        x      = x_sam,
        n      = n_points,
        from   = x_range[1],
        to     = x_range[2],
        bounds = c(x$truncation[["lower"]], x$truncation[["upper"]])
      )
      x_seq               <- x_density$x
      x_den               <- x_density$y
      boundary_reflection <- isTRUE(attr(x_density, "boundary_reflection"))
    }
  }else{

    if(is.prior.discrete(x)){
      x_seq <- unique(round(x_seq))
    }

    x_den <- mpdf(x, x_seq)
    x_sam <- NULL
  }


  # set the endpoints to zero if they correspond to truncation
  if(truncate_end){
    if(isTRUE(all.equal(x$truncation[["lower"]], x_seq[1])) | x$truncation[["lower"]] >= x_seq[1]){
      x_den <- c(0, x_den)
      x_seq <- c(x_seq[1], x_seq)
    }
    if(isTRUE(all.equal(x$truncation[["upper"]], x_seq[length(x_seq)])) | x$truncation[["upper"]] <= x_seq[length(x_seq)]){
      x_den <- c(x_den, 0)
      x_seq <- c(x_seq, x_seq[length(x_seq)])
    }
  }



  # transform the output, if requested
  if(!is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
    x_den   <- .density.prior_transformation_y(x_seq, x_den, transformation, transformation_arguments)
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.simple")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))
  if(boundary_reflection){
    attr(out, "boundary_reflection") <- TRUE
  }

  return(out)
}
.density.prior.point                  <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments){

  # return the samples if requested
  if(force_samples){
    x_sam <- rng(x, n_samples)
  }else{
    x_sam <- NULL
  }

  x_seq <- x$parameters[["location"]]
  x_den <- 1


  # transform the output, if requested
  if(!is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.point", if(is.prior.orthonormal(x)) "density.prior.orthonormal" else if(is.prior.meandif(x)) "density.prior.meandif")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))

  return(out)
}
.density.prior.weightfunction         <- function(x, x_seq, x_range, n_points, n_samples, force_samples, individual){

  # create either distribution for the individual weights or the whole weightfunction
  if(individual){

    out <- list()
    out_types <- .density.prior_type(x)
    components <- .weightfunction_marginal_components(x)
    component_bounds <- lapply(components, .density.prior_weightfunction_component_bounds)
    density_boundary_reflection <- rep(FALSE, length(components))

    if(force_samples | .density.prior_need_samples(x)){
      x_sam <- rng(x, n_samples)
      density_ind <- which(out_types != "point")
      densities <- lapply(density_ind, function(i){
        .density_kde_boundary(
          x      = x_sam[,i],
          n      = n_points,
          from   = x_range[1],
          to     = x_range[2],
          bounds = component_bounds[[i]]
        )
      })
      if(length(densities) > 0L){
        x_seq <- densities[[1]]$x
        x_den <- matrix(0, nrow = length(x_seq), ncol = ncol(x_sam))
        for(j in seq_along(density_ind)){
          x_den[,density_ind[j]] <- densities[[j]]$y
          density_boundary_reflection[density_ind[j]] <- isTRUE(attr(densities[[j]], "boundary_reflection"))
        }
      }else{
        x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
        x_den <- matrix(0, nrow = length(x_seq), ncol = ncol(x_sam))
      }
    }else{
      x_den <- mpdf(x, x_seq)
      x_sam <- NULL
    }

    for(i in 1:ncol(x_den)){

      temp_samples <- if(is.null(x_sam)) NULL else x_sam[,i]
      temp_x_seq <- x_seq
      temp_y_den <- x_den[,i]
      temp_bounds <- component_bounds[[i]]

      if(out_types[i] != "point"){
        at_lower_bound <- is.finite(temp_bounds[1]) &&
          (isTRUE(all.equal(temp_bounds[1], temp_x_seq[1])) | temp_bounds[1] >= temp_x_seq[1])
        at_upper_bound <- is.finite(temp_bounds[2]) &&
          (isTRUE(all.equal(temp_bounds[2], temp_x_seq[length(temp_x_seq)])) | temp_bounds[2] <= temp_x_seq[length(temp_x_seq)])

        if(at_lower_bound){
          temp_y_den <- c(0, temp_y_den)
          temp_x_seq <- c(temp_x_seq[1], temp_x_seq)
        }
        if(at_upper_bound){
          temp_y_den <- c(temp_y_den, 0)
          temp_x_seq <- c(temp_x_seq, temp_x_seq[length(temp_x_seq)])
        }
      }

      # create the output object
      if(out_types[i] == "point"){
        temp_out <- list(
          call    = call("density", print(x, silent = TRUE)),
          bw      = NULL,
          n       = n_points,
          x       = components[[i]]$location,
          y       = 1,
          samples = temp_samples
        )
      }else{
        temp_out <- list(
          call    = call("density", print(x, silent = TRUE)),
          bw      = NULL,
          n       = n_points,
          x       = temp_x_seq,
          y       = temp_y_den,
          samples = temp_samples
        )
      }


      class(temp_out) <- c("density", "density.prior", paste0("density.prior.",out_types[i]))
      attr(temp_out, "x_range") <- x_range
      attr(temp_out, "y_range") <- if(out_types[i] == "point") c(0, 1) else c(0, max(temp_y_den))
      attr(temp_out, "steps")   <- c(x$bins$lower[i], x$bins$upper[i])
      if(density_boundary_reflection[i]){
        attr(temp_out, "boundary_reflection") <- TRUE
      }

      out[[i]] <- temp_out
    }

  }else{

    # weightfunction specific stuff
    x_seq     <- .weightfunction_local_cuts(x)
    x_seq_rep <- c(1, sort(rep(2:(length(x_seq)-1), 2)) ,length(x_seq))
    x_val_rep <- sort(rep(1:(length(x_seq)-1), 2))
    if(force_samples | .density.prior_need_samples(x)){
      x_sam  <- rng(x, n_samples)
      x_lCI  <- apply(x_sam, 2, stats::quantile, probs = .025)
      x_uCI  <- apply(x_sam, 2, stats::quantile, probs = .975)
      x_mean <- apply(x_sam, 2, mean)
    }else{
      x_sam  <- NULL
      x_lCI  <- mquant(x, .025)
      x_uCI  <- mquant(x, .975)
      x_mean <- mean(x)
    }

    out <- list(
      call    = call("density", print(x, silent = TRUE)),
      bw      = NULL,
      n       = n_points,
      x       = x_seq[x_seq_rep],
      y       = x_mean[x_val_rep],
      y_lCI   = x_lCI[x_val_rep],
      y_uCI   = x_uCI[x_val_rep],
      samples = x_sam
    )


    class(out) <- c("density", "density.prior", "density.prior.weightfunction")
    attr(out, "x_range") <- c(0, 1)
    attr(out, "y_range") <- c(0, max(1, x_mean, x_lCI, x_uCI, na.rm = TRUE))
  }

  return(out)
}
.density.prior.PETPEESE               <- function(x, x_seq, x_range, n_points, n_samples, force_samples, individual, transformation, transformation_arguments, truncate_end){

  # create either distribution for the parameter or the PET/PEESE function
  if(individual){

    if(is.prior.point(x)){
      out <- .density.prior.point(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)
    }else if(is.prior.simple(x)){
      out <- .density.prior.simple(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
    }

  }else{

    # get the samples to estimate density / obtain the density directly
    if(force_samples | .density.prior_need_samples(x)){
      x_sam  <- rng(x, n_samples)
      x_med  <- stats::quantile(x_sam, .500)
      x_lCI  <- stats::quantile(x_sam, .025)
      x_uCI  <- stats::quantile(x_sam, .975)
    }else{
      x_med  <- quant(x, .500)
      x_lCI  <- quant(x, .025)
      x_uCI  <- quant(x, .975)
      x_sam  <- NULL
    }


    # transform the output, if requested
    if(!is.null(transformation)){
      x_med   <- .density.prior_transformation_x(x_med,   transformation, transformation_arguments)
      x_lCI   <- .density.prior_transformation_x(x_lCI,   transformation, transformation_arguments)
      x_uCI   <- .density.prior_transformation_x(x_uCI,   transformation, transformation_arguments)
      x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
      if(!is.null(x_sam)){
        x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
      }
    }


    # compute the PET/PEESE
    if(is.prior.PET(x)){
      y_med   = x_med  * x_seq
      y_lCI   = x_lCI  * x_seq
      y_uCI   = x_uCI  * x_seq
    }else if(is.prior.PEESE(x)){
      y_med   = x_med  * x_seq^2
      y_lCI   = x_lCI  * x_seq^2
      y_uCI   = x_uCI  * x_seq^2
    }


    out <- list(
      call    = call("density", print(x, silent = TRUE)),
      bw      = NULL,
      n       = n_points,
      x       = x_seq,
      y       = y_med,
      y_lCI   = y_lCI,
      y_uCI   = y_uCI,
      samples = x_sam
    )


    class(out) <- c("density", "density.prior", if(is.prior.PET(x)) "density.prior.PET" else if(is.prior.PEESE(x)) "density.prior.PEESE")
    attr(out, "x_range") <- range(x_seq)
    attr(out, "y_range") <- range(y_med)
  }

  return(out)
}
.density.prior.orthonormal_or_meandif <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  boundary_reflection <- FALSE

  # get the samples to estimate density / obtain the density directly
  if(force_samples | .density.prior_need_samples(x)){

    if(is.na(x$parameters[["K"]]) && !is.null(attr(x, "levels"))){
      x$parameters[["K"]] <- .get_prior_factor_levels(x)
    }else if(is.na(x$parameters[["K"]])){
      x$parameters[["K"]] <- 1
      warning("number of factor levels / dimensionality of the prior distribution was not specified -- assuming two factor levels")
    }

    x_sam <- rng(x, n_samples)
    x_sam <- as.vector(x_sam)
    x_density <- .density_kde_boundary(
      x      = x_sam,
      n      = n_points,
      from   = x_range[1],
      to     = x_range[2],
      bounds = c(x$truncation[["lower"]], x$truncation[["upper"]])
    )
    x_seq               <- x_density$x
    x_den               <- x_density$y
    boundary_reflection <- isTRUE(attr(x_density, "boundary_reflection"))

  }else{
    x_den <- mpdf(x, x_seq)
    x_sam <- NULL
  }


  # set the endpoints to zero if they correspond to truncation
  if(truncate_end){
    if(isTRUE(all.equal(x$truncation[["lower"]], x_seq[1])) | x$truncation[["lower"]] >= x_seq[1]){
      x_den <- c(0, x_den)
      x_seq <- c(x_seq[1], x_seq)
    }
    if(isTRUE(all.equal(x$truncation[["upper"]], x_seq[length(x_seq)])) | x$truncation[["upper"]] <= x_seq[length(x_seq)]){
      x_den <- c(x_den, 0)
      x_seq <- c(x_seq, x_seq[length(x_seq)])
    }
  }


  # transform the output, if requested
  if(!is.null(transformation)){
    message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the orthonormal/meandif contrasts to the differences from the mean.")
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
    if(!is.null(x_sam)){
      x_sam <- .density.prior_transformation_x(x_sam,   transformation, transformation_arguments)
    }
    x_den   <- .density.prior_transformation_y(x_seq, x_den, transformation, transformation_arguments)
  }


  # create the output object
  out <- list(
    call    = call("density", print(x, silent = TRUE)),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_den,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", if(is.prior.orthonormal(x)) "density.prior.orthonormal" else if(is.prior.meandif(x)) "density.prior.meandif")
  attr(out, "x_range") <- x_range
  attr(out, "y_range") <- c(0, max(x_den))
  if(boundary_reflection){
    attr(out, "boundary_reflection") <- TRUE
  }

  return(out)
}
.density.prior.spike_and_slab         <- function(x, x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end){

  density_variable  <- .density.prior.simple(.get_spike_and_slab_variable(x), x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments, truncate_end)
  density_inclusion <- .density.prior.point(prior(distribution = "spike", parameters = list(location = 0)), x_seq, x_range, n_points, n_samples, force_samples, transformation, transformation_arguments)

  inclusion_prob <- mean(.get_spike_and_slab_inclusion(x))
  density_variable$y  <- density_variable[["y"]]  * inclusion_prob
  density_inclusion$y <- density_inclusion[["y"]] * (1 - inclusion_prob)

  attr(density_variable,  "y_range") <- attr(density_variable, "y_range")  * inclusion_prob
  attr(density_inclusion, "y_range") <- attr(density_inclusion, "y_range") * (1 - inclusion_prob)

  # create the output object
  out <- list(
    call      = call("density", print(x, silent = TRUE)),
    variable  = density_variable,
    inclusion = density_inclusion
  )


  class(out) <- c("density", "density.prior.spike_and_slab")
  attr(out, "x_range") <- range(c(attr(density_variable, "x_range"), attr(density_inclusion, "x_range")))
  attr(out, "y_range_variable")  <- attr(density_variable,  "y_range")
  attr(out, "y_range_inclusion") <- attr(density_inclusion, "y_range")

  return(out)
}


#' @title Prior range
#'
#' @description Computes range of a prior
#' distribution (if the prior distribution is
#' unbounded range from \code{quantiles} to
#' \code{1 -quantiles}) is returned.
#'
#' @param x a prior
#' @param quantiles quantile to be returned in
#' case of unbounded distribution.
#' @param ... additional arguments
#' @param na.rm unused
#'
#' @return \code{range.prior} returns a numeric vector of
#' length with a plotting range of a prior distribution.
#'
#' @seealso [prior()]
#' @rdname range.prior
#' @export
range.prior  <- function(x, quantiles = NULL, ..., na.rm = FALSE){

  .check_prior(x)
  if(!is.null(quantiles)){
    check_real(quantiles, "quantiles", upper = 0.5, allow_bound = FALSE)
  }else{
    quantiles <- .range.prior_quantile_default(x)
  }


  x_range <- c(NA, NA)

  if(is.prior.weightfunction(x)){
    return(.weightfunction_range(x, quantiles))
  }
  if(is_prior_phacking(x) || is_prior_bias(x)){
    .selection_prior_stop_unsupported_generic("range", x)
  }

  if(is.infinite(x[["truncation"]][["lower"]])){
    x_range[1] <- min(mquant(x, quantiles), na.rm = TRUE)
  }else{
    x_range[1] <- x[["truncation"]][["lower"]]
  }

  if(is.infinite(x[["truncation"]][["upper"]])){
    x_range[2] <- max(mquant(x, 1 - quantiles), na.rm = TRUE)
  }else{
    x_range[2] <- x[["truncation"]][["upper"]]
  }

  return(x_range)
}



# helper functions
.density_kde_boundary       <- function(x, n, from = NULL, to = NULL, bounds = c(-Inf, Inf), na.rm = FALSE, ...){

  if(!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds)){
    stop("'bounds' must be a numeric vector of length 2.", call. = FALSE)
  }
  if(bounds[1] > bounds[2]){
    stop("The lower boundary must be lower than or equal to the upper boundary.", call. = FALSE)
  }

  # Estimate the bandwidth on the original sample. The reflected sample is only
  # used to correct leakage at the true support boundaries.
  dots <- list(...)
  density_args <- c(list(x = x, n = n, na.rm = na.rm), dots)
  if(!is.null(from)){
    density_args$from <- from
  }
  if(!is.null(to)){
    density_args$to <- to
  }

  density_base <- do.call(stats::density, density_args)
  attr(density_base, "boundary_reflection") <- FALSE

  if(!any(is.finite(bounds))){
    return(density_base)
  }

  # 'from' and 'to' define only the evaluation grid. 'bounds' defines the true
  # support and therefore where reflection is applied.
  x_ref <- x
  weights_ref <- dots[["weights"]]
  if(na.rm){
    keep <- !is.na(x_ref)
    if(!is.null(weights_ref)){
      keep <- keep & !is.na(weights_ref)
    }
    x_ref <- x_ref[keep]
    if(!is.null(weights_ref)){
      weights_ref <- weights_ref[keep]
    }
  }

  if(is.null(weights_ref)){
    weights_ref <- rep(1 / length(x_ref), length(x_ref))
  }

  x_reflected <- x_ref
  weights_reflected <- weights_ref
  if(is.finite(bounds[1])){
    x_reflected <- c(x_reflected, 2 * bounds[1] - x_ref)
    weights_reflected <- c(weights_reflected, weights_ref)
  }
  if(is.finite(bounds[2])){
    x_reflected <- c(x_reflected, 2 * bounds[2] - x_ref)
    weights_reflected <- c(weights_reflected, weights_ref)
  }

  density_reflected_args <- dots
  density_reflected_args[c("adjust", "bw", "subdensity", "warnWbw", "weights", "width")] <- NULL
  density_reflected_args <- c(
    list(
      x          = x_reflected,
      bw         = density_base$bw,
      n          = n,
      from       = if(is.null(from)) density_base$x[1] else from,
      to         = if(is.null(to)) density_base$x[length(density_base$x)] else to,
      weights    = weights_reflected,
      subdensity = TRUE,
      warnWbw    = FALSE,
      na.rm      = FALSE
    ),
    density_reflected_args
  )

  density_reflected <- do.call(stats::density, density_reflected_args)
  if(is.finite(bounds[1])){
    density_reflected$y[density_reflected$x < bounds[1]] <- 0
  }
  if(is.finite(bounds[2])){
    density_reflected$y[density_reflected$x > bounds[2]] <- 0
  }
  attr(density_reflected, "boundary_reflection") <- TRUE

  return(density_reflected)
}
.density.prior_weightfunction_component_bounds <- function(component){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      lower <- component$prior$truncation[["lower"]]
      upper <- component$prior$truncation[["upper"]]
      if(component$scale == "omega"){
        c(lower, upper)
      }else{
        c(exp(lower), exp(upper))
      }
    },
    "one_minus_product_beta" = c(0, 1)
  )
}
.density.prior_weightfunction_components_bounds <- function(components){

  bounds <- do.call(rbind, lapply(components, .density.prior_weightfunction_component_bounds))
  c(min(bounds[,1], na.rm = TRUE), max(bounds[,2], na.rm = TRUE))
}
.density.prior_need_samples   <- function(prior){

  return(FALSE)
}
.density.prior_type           <- function(prior){
  if(is.prior.point(prior)){
    return("point")
  }else if(is.prior.simple(prior)){
    return("simple")
  }else if(is.prior.weightfunction(prior)){
    components <- .weightfunction_marginal_components(prior)
    return(vapply(components, function(component){
      if(component$type == "point") "point" else "simple"
    }, character(1)))
  }else if(is_prior_phacking(prior) || is_prior_bias(prior)){
    .selection_prior_stop_unsupported_generic("density", prior)
  }else if(is.prior.orthonormal(prior)){
    return("orthonormal")
  }else if(is.prior.meandif(prior)){
    return("meandif")
  }
}
.range.prior_quantile_default <- function(prior){

  switch(
    prior[["distribution"]],
    "normal"    = .005,
    "lognormal" = .005,
    "t"         = .010,
    "gamma"     = .010,
    "invgamma"  = .010,
    "moment"    = .010,
    "invmoment" = .010,
    "beta"      = .005,
    "exp"       = .005,
    "uniform"   = .005,
    "point"     = .005,
    "weightfunction" = .005,
    "mnormal"    = .005,
    "mt"         = .010,
    "dirichlet"  = .005
  )

}

# transformation functions
.density.prior_transformation_x         <- function(x, transformation, transformation_arguments = NULL){

  arg <- list(x = x)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  do.call(.density.prior_transformation_functions(transformation)$fun, arg)
}
.density.prior_transformation_inv_x     <- function(x, transformation, transformation_arguments = NULL){

  arg <- list(x = x)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  do.call(.density.prior_transformation_functions(transformation)$inv, arg)
}
.density.prior_transformation_y         <- function(x, y, transformation, transformation_arguments = NULL){

  x_inv <- .density.prior_transformation_inv_x(x, transformation, transformation_arguments)
  arg <- list(x = x_inv)
  for(i in seq_along(transformation_arguments)){
    arg[[names(transformation_arguments)[i]]] <- transformation_arguments[[i]]
  }

  y / abs(do.call(.density.prior_transformation_functions(transformation)$jac, arg))
}
.density.prior_transformation_functions <- function(transformation){

  if(is.character(transformation) & length(transformation) == 1){

    return(switch(
      transformation,
      "lin" = list(
        fun = function(x, a = 0, b = 1)a + b * x,
        inv = function(x, a = 0, b = 1)(x - a) / b,
        jac = function(x, a = 0, b = 1)b
      ),
      "exp_lin" = list(
        # Exponential-linear transformation: exp(a + b * log(x))
        # Used for log-intercept unscaling where: intercept_orig = exp(log(intercept_z) * b + a)
        # When a = 0 and b = 1, this is identity: exp(log(x)) = x
        fun = function(x, a = 0, b = 1) exp(a + b * log(x)),
        inv = function(x, a = 0, b = 1) exp((log(x) - a) / b),
        jac = function(x, a = 0, b = 1) b * exp(a + b * log(x)) / x
      ),
      "tanh" = list(
        fun = tanh,
        inv = atanh,
        jac = function(x)1 - tanh(x)^2
      ),
      "exp"  = list(
        fun = exp,
        inv = log,
        jac = exp
      )
    ))

  }else if(is.list(transformation) & length(transformation) == 3 & all(names(transformation) %in% c("fun", "inv", "jac"))){

    return(transformation)

  }else{

    stop("Transformation must be either a character vector of length 1 corresponding to one of known transformations ('lin' = linear, 'exp_lin' = exponential-linear for log-intercept, 'tanh' = hyperbolic tangent, 'exp' = exponential) or a list of three functions (fun = transformation function, inv = inverse transformation, jac = derivative of the transformation).")

  }

}
