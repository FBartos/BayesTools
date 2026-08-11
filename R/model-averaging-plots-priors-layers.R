#' @title Add list of prior objects to a plot
#'
#' @param scale_y2 scaling factor for a secondary axis
#' @param ... additional arguments
#' @inheritParams plot_prior_list
#' @inheritParams density.prior
#'
#' @details For base plots, the default \code{scale_y2 = NULL} reuses the
#' probability mapping established by the initial mixed-distribution plot.
#' This keeps separately added prior and posterior point masses comparable.
#'
#' @return \code{lines_prior_list} returns \code{NULL}.
#'
#' @seealso [plot_prior_list()] [geom_prior_list()]
#' @rdname lines_prior_list
#' @export
lines_prior_list <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             individual = FALSE, show_figures = if(individual) 1 else NULL,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, prior_list_mu = NULL, effect_direction = "positive", ...){

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0, allow_NULL = TRUE)
  check_char(effect_direction, "effect_direction", allow_values = c("positive", "negative"))


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else{
    prior_type <- "simple"
  }

  if(prior_type == "PETPEESE"){
    check_list(prior_list_mu, "prior_list_mu", check_length = length(prior_list))
    if(is.prior(prior_list_mu) | !all(sapply(prior_list_mu, is.prior)))
      stop("'prior_list_mu' must be a list of priors (priors for the mu parameter are required for plotting PET-PEESE).")
  }else{
    if(!is.null(prior_list_mu))
      stop("'prior_list_mu' is required only for PET-PEESE plots.")
  }


  # get the plotting range
  if(is.null(xlim) & is.null(x_seq)){
    if(prior_type %in% c("weightfunction", "PETPEESE") & !individual){
      xlim      <- c(0, 1)
    }else if(prior_type == "simple"){
      xlim   <- do.call(rbind, lapply(prior_list, range, quantiles = x_range_quant))
      xlim   <- range(pretty(range(as.vector(xlim))))
    }
  }


  # get the plotting data
  if(prior_type == "weightfunction"){

    # use analytical marginal summaries of the mapped cumulative Dirichlet weights
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    .lines.prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use deterministic linear-combination summaries when supported, with a sampling fallback
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu,
                                                effect_direction = effect_direction)
    .lines.prior.PETPEESE(plot_data = plot_data, ...)

  }else if(prior_type == "simple"){

    # solve analytically
    plot_data <- .plot_data_prior_list.simple(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                              n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                              transformation = transformation, transformation_arguments = transformation_arguments,
                                              transformation_settings = transformation_settings)

    scale_y2_state <- NULL
    if(is.null(scale_y2)){
      scale_y2_state <- .plot_scale_y2_state_current()
      if(!is.null(scale_y2_state)){
        scale_y2 <- scale_y2_state[["scale_y2"]]
      }else{
        scale_y2 <- .get_scale_y2(plot_data, ...)
      }
    }
    if(!is.null(scale_y2_state)){
      .plot_point_mass_warn_outside(plot_data, scale_y2_state[["ylim2"]])
    }
    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        .lines.prior.simple(plot_data[[i]], ...)
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        .lines.prior.point(plot_data[[i]], scale_y2 = scale_y2, ...)
      }
    }

  }

  return(invisible())
}

#' @title Add list of prior objects to a plot
#'
#' @inheritParams lines_prior_list
#' @inheritParams density.prior
#' @inheritParams plot_prior_list
#'
#' @return \code{geom_prior_list} returns an object of class 'ggplot'.
#'
#' @seealso [plot_prior_list()] [lines_prior_list()]
#' @rdname geom_prior_list
#' @export
geom_prior_list  <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             individual = FALSE, show_figures = if(individual) 1 else NULL,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, prior_list_mu = NULL, effect_direction = "positive", ...){

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0, allow_NULL = TRUE)
  check_char(effect_direction, "effect_direction", allow_values = c("positive", "negative"))


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else{
    prior_type <- "simple"
  }

  if(prior_type == "PETPEESE"){
    check_list(prior_list_mu, "prior_list_mu", check_length = length(prior_list))
    if(is.prior(prior_list_mu) | !all(sapply(prior_list_mu, is.prior)))
      stop("'prior_list_mu' must be a list of priors (priors for the mu parameter are required for plotting PET-PEESE).")
  }else{
    if(!is.null(prior_list_mu))
      stop("'prior_list_mu' is required only for PET-PEESE plots.")
  }

  # get the plotting range
  if(is.null(xlim) & is.null(x_seq)){
    if(prior_type %in% c("weightfunction", "PETPEESE") & !individual){
      xlim      <- c(0, 1)
    }else if(prior_type == "simple"){
      xlim   <- do.call(rbind, lapply(prior_list, range, quantiles = x_range_quant))
      xlim   <- range(pretty(range(as.vector(xlim))))
    }
  }


  # get the plotting data
  if(prior_type == "weightfunction"){

    # use analytical marginal summaries of the mapped cumulative Dirichlet weights
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    geom <- .geom_prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use deterministic linear-combination summaries when supported, with a sampling fallback
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu,
                                                effect_direction = effect_direction)
    geom <- .geom_prior.PETPEESE(plot_data = plot_data, ...)

  }else if(prior_type == "simple"){

    # solve analytically
    plot_data <- .plot_data_prior_list.simple(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                              n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                              transformation = transformation, transformation_arguments = transformation_arguments,
                                              transformation_settings = transformation_settings)
    geom     <- NULL
    if(is.null(scale_y2)){
      scale_y2 <- .get_scale_y2(plot_data, ...)
    }
    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        geom <- c(geom, list(.geom_prior.simple(plot_data[[i]], ...)))
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        geom <- c(geom, list(.geom_prior.point(plot_data[[i]], scale_y2 = scale_y2, ...)))
      }
    }

  }

  return(geom)
}
