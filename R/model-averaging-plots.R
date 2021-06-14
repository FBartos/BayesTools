#' @title Plot a list of prior distributions
#'
#' @param prior_list list of prior distributions
#' @param ... additional arguments
#' @inheritParams density.prior
#' @inheritParams plot.prior
#'
#' @seealso [prior()] [lines_prior_list()]  [geom_prior_list()]
#' @export
plot_prior_list <- function(prior_list, plot_type = "base",
                            x_seq = NULL, xlim = NULL, x_range_quant = NULL, n_points = 500,
                            n_samples = 10000, force_samples = FALSE,
                            transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                            rescale_x = FALSE, par_name = NULL, ...){

  # TODO: add plots for individual parameters for weightfunction and PET-PEESE
  individual = FALSE
  show_figures = if(individual) 1 else NULL

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  # check that there is no mixing of PET-PEESE and weightfunctions
  if(any(sapply(prior_list, is.prior.weightfunction)) & (any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))))
    stop("weightfunction and PET-PEESE priors cannot be mixed within a 'prior_list'.")


  # join the same priors
  prior_list <- .simplify_prior_list(prior_list)


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else{
    prior_type <- "simple"
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

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    plot <- .plot.prior.weightfunction(prior_list, plot_type = plot_type, plot_data = plot_data, rescale_x = rescale_x, par_name = par_name, ...)

  }else if(prior_type == "PETPEESE"){

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings)
    plot <- .plot.prior.PETPEESE(prior_list, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)

  }else if(prior_type == "simple"){

    # solve analytically
    plot_data <- .plot_data_prior_list.simple(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                              n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                              transformation = transformation, transformation_arguments = transformation_arguments,
                                              transformation_settings = transformation_settings)
    plot <- .plot_prior_list.both(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)

  }


  if(plot_type == "ggplot"){
    return(plot)
  }else{
    return(invisible())
  }
}

.plot_prior_list.both             <- function(plot_data, plot_type, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- range(as.vector(sapply(plot_data, attr, which = "x_range")))

  main      <- ""
  xlab      <- if(!is.null(par_name)){parse(text = bquote(.(par_name)))} else ""

  if(any(sapply(plot_data, inherits, what = "density.prior.simple")) & any(sapply(plot_data, inherits, what = "density.prior.point"))){
    type  <- "both"
    ylab  <- "Density"
    ylab2 <- "Probability"
    ylim  <- range(as.vector(sapply(plot_data[sapply(plot_data, inherits, what = "density.prior.simple")], attr, which = "y_range")))
    ylim2 <- range(as.vector(sapply(plot_data[sapply(plot_data, inherits, what = "density.prior.point")],  attr, which = "y_range")))
  }else if(any(sapply(plot_data, inherits, what = "density.prior.simple"))){
    type  <- "simple"
    ylab  <- "Density"
    ylim  <- range(as.vector(sapply(plot_data[sapply(plot_data, inherits, what = "density.prior.simple")], attr, which = "y_range")))
    ylab2 <- NULL
    ylim2 <- NULL
  }else if(any(sapply(plot_data, inherits, what = "density.prior.point"))){
    type  <- "point"
    ylab  <- "Probability"
    ylim  <- range(as.vector(sapply(plot_data[sapply(plot_data, inherits, what = "density.prior.point")],  attr, which = "y_range")))
    ylab2 <- NULL
    ylim2 <- NULL
  }


  # add it to the user input if desired
  if(is.null(dots[["main"]]))  dots$main  <-  main
  if(is.null(dots[["xlab"]]))  dots$xlab  <-  xlab
  if(is.null(dots[["ylab"]]))  dots$ylab  <-  ylab
  if(is.null(dots[["ylab2"]])) dots$ylab2 <-  ylab2
  if(is.null(dots[["xlim"]]))  dots$xlim  <-  xlim
  if(is.null(dots[["ylim"]]))  dots$ylim  <-  ylim
  if(is.null(dots[["ylim2"]])) dots$ylim2 <-  ylim2


  if(plot_type == "base"){

    .plot.prior_empty(type, dots)
    scale_y2 <- .get_scale_y2(plot_data, dots)
    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        do.call(.lines.prior.simple, c(list(plot_data[[i]]), dots))
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        do.call(.lines.prior.point,  c(list(plot_data[[i]]), scale_y2 = scale_y2, dots))
      }
    }
    plot <- list(scale_y2 = scale_y2)

  }else if(plot_type == "ggplot"){

    plot     <- .ggplot.prior_empty(type, dots)
    scale_y2 <- .get_scale_y2(plot_data, dots)
    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        plot <- plot + do.call(.geom_prior.simple, c(list(plot_data[[i]]), dots))
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        plot <- plot + do.call(.geom_prior.point,  c(list(plot_data[[i]]), scale_y2 = scale_y2, dots))
      }
    }

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

.plot_data_prior_list.weightfunction <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples){

  prior_odds  <- sapply(prior_list, function(p)p$prior_odds)
  mixing_prop <- prior_odds / sum(prior_odds)

  # get the samples
  samples_list <- list()
  for(i in seq_along(prior_list)){
    if(is.prior.weightfunction(prior_list[[i]])){
      samples_list[[i]] <- rng(prior_list[[i]], round(n_samples * mixing_prop[i]))
    }else{
      samples_list[[i]] <- list()
    }

  }

  # merge the samples
  omega_mapping <- weightfunctions_mapping(prior_list[sapply(prior_list, is.prior.weightfunction)])
  omega_cuts    <- weightfunctions_mapping(prior_list[sapply(prior_list, is.prior.weightfunction)], cuts_only = TRUE)

  # join samples
  samples    <- matrix(nrow = 0, ncol = length(omega_cuts) - 1)
  for(i in seq_along(samples_list)){
    if(is.prior.weightfunction(prior_list[[i]])){
      samples <- rbind(samples, samples_list[[i]][,omega_mapping[[i]]])
    }else{
      samples <- rbind(samples, matrix(1, ncol = length(omega_cuts) - 1, nrow = round(n_samples * mixing_prop[i])))
    }
  }

  x_lCI  <- apply(samples, 2, stats::quantile, probs = .025)
  x_uCI  <- apply(samples, 2, stats::quantile, probs = .975)
  x_mean <- apply(samples, 2, mean)

  x_seq     <- omega_cuts
  x_seq_rep <- c(1, sort(rep(2:(length(x_seq)-1), 2)) ,length(x_seq))
  x_val_rep <- sort(rep(1:(length(x_seq)-1), 2))


  out <- list(
    call    = call("density", "weightfunction list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq[x_seq_rep],
    y       = x_mean[x_val_rep],
    y_lCI   = x_lCI[x_val_rep],
    y_uCI   = x_uCI[x_val_rep],
    samples = samples
  )


  class(out) <- c("density", "density.prior", "density.prior.weightfunction")
  attr(out, "x_range") <- c(0, 1)
  attr(out, "y_range") <- c(0, 1)

  return(out)
}
.plot_data_prior_list.PETPEESE       <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples,
                                                 transformation, transformation_arguments, transformation_settings){
  if(is.null(x_seq)){
    x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  # specify it on the transformed range if requested
  if(transformation_settings & !is.null(transformation)){
    x_seq   <- .density.prior_transformation_inv_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_inv_x(x_range, transformation, transformation_arguments)
  }

  prior_odds  <- sapply(prior_list, function(p)p$prior_odds)
  mixing_prop <- prior_odds / sum(prior_odds)

  # get the samples
  samples_list <- list()
  for(i in seq_along(prior_list)){
    if(is.prior.PET(prior_list[[i]])){
      samples_list[[i]] <- cbind(rng(prior_list[[i]], round(n_samples * mixing_prop[i])), rep(0, length = round(n_samples * mixing_prop[i])))
    }else if(is.prior.PEESE(prior_list[[i]])){
      samples_list[[i]] <- cbind(rep(0, length = round(n_samples * mixing_prop[i])), rng(prior_list[[i]], round(n_samples * mixing_prop[i])))
    }else{
      samples_list[[i]] <- matrix(0, nrow = round(n_samples * mixing_prop[i]), ncol = 2)
    }
  }
  samples <- do.call(rbind, samples_list)


  # transform the parameter if requested
  if(!is.null(transformation)){
    samples   <- .density.prior_transformation_x(samples,   transformation, transformation_arguments)
  }


  # compute PET-PEESE
  x_sam  <- matrix(samples[,1], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq, nrow = length(samples), ncol = length(x_seq), byrow = TRUE) + matrix(samples[,2], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq^2, nrow = length(samples), ncol = length(x_seq), byrow = TRUE)
  x_med  <- apply(x_sam, 2, stats::quantile, prob = .500)
  x_lCI  <- apply(x_sam, 2, stats::quantile, prob = .025)
  x_uCI  <- apply(x_sam, 2, stats::quantile, prob = .975)


  out <- list(
    call    = call("density", "PET-PEESE list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_med,
    y_lCI   = x_lCI,
    y_uCI   = x_uCI,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.PETPEESE")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- range(x_med)

  return(out)
}
.plot_data_prior_list.simple         <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples, force_samples, individual,
                                                 transformation, transformation_arguments, transformation_settings){

  prior_odds  <- sapply(prior_list, function(p)p$prior_odds)
  mixing_prop <- prior_odds / sum(prior_odds)

  plot_data <- list()
  for(i in seq_along(prior_list)){
    plot_data[[i]] <- density(prior_list[[i]], x_seq = x_seq, x_range = x_range, x_range_quant = x_range_quant,
                              n_points = n_points, n_samples = round(n_samples * mixing_prop[i]), force_samples = force_samples,
                              transformation = transformation, transformation_arguments = transformation_arguments,
                              transformation_settings = transformation_settings, individual = individual, truncate_end = FALSE)
  }

  x_sam    <- NULL
  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL

  for(i in seq_along(plot_data)){

    if(force_samples){
      x_sam <- plot_data$samples
    }

    # align points and densities
    if(inherits(plot_data[[i]], "density.prior.point")){
      x_points <- c(x_points, plot_data[[i]]$x[plot_data[[i]]$y != 0])
      y_points <- c(y_points, mixing_prop[i])
    }else if(inherits(plot_data[[i]], "density.prior.simple")){
      x_den <- rbind(x_den, plot_data[[i]]$x)
      y_den <- rbind(y_den, plot_data[[i]]$y * mixing_prop[i])
    }
  }

  # deal with continuous densities
  if(!is.null(y_den)){
    y_den <- apply(y_den, 2, sum)
    if(any(sapply(1:nrow(x_den), function(i) isFALSE(all.equal(x_den[1,], x_den[i,])))))
      stop("non-matching x-coordinates")
    x_den <- x_den[1,]

    # set the endpoints to zero if they correspond to truncation
    prior_list_simple <- prior_list[!sapply(prior_list, is.prior.point)]
    prior_list_simple_lower <- min(sapply(prior_list_simple, function(p) p$truncation[["lower"]]))
    prior_list_simple_upper <- max(sapply(prior_list_simple, function(p) p$truncation[["upper"]]))
    if(!is.null(transformation)){
      prior_list_simple_lower   <- .density.prior_transformation_x(prior_list_simple_lower, transformation, transformation_arguments)
      prior_list_simple_upper   <- .density.prior_transformation_x(prior_list_simple_upper, transformation, transformation_arguments)
    }
    if(isTRUE(all.equal(prior_list_simple_lower, x_den[1])) | prior_list_simple_lower >= x_den[1]){
      y_den <- c(0, y_den)
      x_den <- c(x_den[1], x_den)
    }
    if(isTRUE(all.equal(prior_list_simple_upper, x_den[length(x_den)])) | prior_list_simple_upper <= x_den[length(x_den)]){
      y_den <- c(y_den, 0)
      x_den <- c(x_den, x_den[length(x_den)])
    }
  }


  # create the output object
  out <- list()

  # add continuous densities
  if(!is.null(y_den)){
    out_den    <- list(
      call    = call("density", "list priors"),
      bw      = NULL,
      n       = n_points,
      x       = x_den,
      y       = y_den,
      samples = x_sam
    )

    class(out_den) <- c("density", "density.prior", "density.prior.simple")
    attr(out_den, "x_range") <- x_range
    attr(out_den, "y_range") <- c(0, max(y_den))

    out[["density"]] <- out_den
  }

  # add spikes
  if(!is.null(y_points)){
    for(i in seq_along(y_points)){
      temp_points <- list(
        call    = call("density", paste0("point", i)),
        bw      = NULL,
        n       = n_points,
        x       = x_points[i],
        y       = y_points[i],
        samples = NULL
      )

      class(temp_points) <- c("density", "density.prior", "density.prior.point")
      attr(temp_points, "x_range") <- x_range
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}
.simplify_prior_list <- function(prior_list){

  # return the input with fewer than 2 priors
  if(length(prior_list) < 2){
    return(prior_list)
  }

  new_prior_list <- prior_list
  for(i in seq_along(new_prior_list)){
    new_prior_list[[i]][["prior_odds"]] <- NULL
  }

  are_equal <- do.call(rbind, lapply(new_prior_list, function(p)sapply(new_prior_list, identical, y = p)))
  are_equal <- are_equal[!duplicated(are_equal) & apply(are_equal, 1, sum) > 1,,drop = FALSE]

  # return the input with no matches
  if(nrow(are_equal) == 0){
    return(prior_list)
  }

  # find the duplicates and collect prior odds
  prior_odds <- unname(sapply(prior_list, function(p)p$prior_odds))
  to_remove  <- NULL
  for(i in 1:nrow(are_equal)){
    this_ind    <- c(1:ncol(are_equal))[are_equal[i,]]
    this_unique <- this_ind[1]
    prior_odds[this_unique] <- sum(prior_odds[this_ind])
    to_remove   <- c(to_remove, this_ind[-1])
  }

  # return prior odds
  for(i in seq_along(prior_list)){
    prior_list[[i]][["prior_odds"]] <- prior_odds[i]
  }

  # remove the duplicates
  prior_list[to_remove] <- NULL

  return(prior_list)
}


#' @title Add list of prior objects to a plot
#'
#' @param scale_y2 scaling factor for a secondary axis
#' @param ... additional arguments
#' @inheritParams plot_prior_list
#' @inheritParams density.prior
#'
#' @seealso [plot_prior_list()] [geom_prior_list()]
#' @rdname lines_prior_list
#' @export
lines_prior_list <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, ...){

  # TODO: add plots for individual parameters for weightfunction and PET-PEESE
  individual = FALSE
  show_parameter = if(individual) 1 else NULL

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0, allow_NULL = TRUE)


  # join the same priors
  prior_list <- .simplify_prior_list(prior_list)


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else{
    prior_type <- "simple"
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

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    .lines.prior.weightfunction(prior_list, plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings)
    .lines.prior.PETPEESE(prior_list, plot_data = plot_data, ...)

  }else if(prior_type == "simple"){

    # solve analytically
    plot_data <- .plot_data_prior_list.simple(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                              n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                              transformation = transformation, transformation_arguments = transformation_arguments,
                                              transformation_settings = transformation_settings)

    if(is.null(scale_y2)){
      scale_y2 <- .get_scale_y2(plot_data, ...)
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
#' @seealso [plot_prior_list()] [lines_prior_list()]
#' @rdname geom_prior_list
#' @export
geom_prior_list  <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, ...){

  # TODO: add plots for individual parameters for weightfunction and PET-PEESE
  individual = FALSE
  show_parameter = if(individual) 1 else NULL

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0, allow_NULL = TRUE)


  # join the same priors
  prior_list <- .simplify_prior_list(prior_list)


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else{
    prior_type <- "simple"
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

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    geom <- .geom_prior.weightfunction(prior_list, plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings)
    geom <- .geom_prior.PETPEESE(prior_list, plot_data = plot_data, ...)

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


#' @title Plot samples from the mixed posterior distributions
#'
#' @param samples samples from a posterior distribution for a
#' parameter generated by [mix_posteriors].
#' @param parameter parameter name to be plotted. Use \code{"PETPEESE"}
#' for PET-PEESE plot with parameters \code{"PET"} and \code{"PEESE"},
#' and \code{"weightfunction"} for plotting a weightfunction with
#' parameters \code{"omega"}.
#' @param prior whether prior distribution should be added to the figure
#' @param dots_prior additional arguments for the prior distribution plot
#' @param ... additional arguments
#' @inheritParams density.prior
#' @inheritParams plot.prior
#'
#' @seealso [prior()] [lines_prior_list()]  [geom_prior_list()]
#' @export
plot_posterior <- function(samples, parameter, plot_type = "base", prior = FALSE,
                           x_seq = NULL, xlim = NULL, x_range_quant = NULL, n_points = 1000,
                           n_samples = 10000, force_samples = FALSE,
                           transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                           rescale_x = FALSE, par_name = NULL, dots_prior = list(), ...){

  # TODO: add plots for individual parameters for weightfunction and PET-PEESE
  individual = FALSE
  show_figures = if(individual) 1 else NULL

  # check input
  check_list(samples, "prior_list")
  if(any(!sapply(samples, inherits, what = "mixed_posteriors")))
    stop("'samples' must be a be an object generated by 'mix_posteriors' function.")
  check_char(parameter, "parameter")
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  if(!is.null(transformation)){
    if(is.character(transformation)){
      check_char(transformation, "transformation")
    }else if(is.list(transformation)){
      check_list(check_list, "check_list", check_length = 3, check_names = c("fun", "inv", "jac"), all_objects = TRUE)
    }else{
      stop("Uknown format of the 'transformation' argument.")
    }
  }
  check_list(transformation_arguments, "transformation_arguments", allow_NULL = TRUE)
  check_bool(transformation_settings, "transformation_settings")

  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter <- "omega"
  }else if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) == "petpeese"){
    parameter <- "PETPEESE"
  }

  # get the plotting range
  if(is.null(xlim) & is.null(x_seq)){
    if(parameter %in% c("omega", "PETPEESE") & !individual){
      xlim      <- c(0, 1)
    }else{
      # use the data range otherwise
      xlim   <- NULL
    }
  }


  if(parameter == "omega"){

    plot_data <- .plot_data_samples.weightfunction(samples, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant, n_points = n_points)

    # add priors, if requested
    if(prior){
      prior_list      <- attr(samples[[parameter]], "prior_list")
      plot_data_prior <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                              n_points = n_points, n_samples = n_samples)

      # transplant common xlim and ylim
      plot_data_joined <- list(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      ylim <- range(as.vector(sapply(plot_data_joined, attr, which = "y_range")))
      attr(plot_data_prior, "x_range") <- xlim
      attr(plot_data_prior, "y_range") <- ylim
      dots_prior <- .transfer_dots(dots_prior, ...)

      plot <- do.call(.plot.prior.weightfunction, c(x = list(prior_list), plot_data = list(plot_data_prior), rescale_x = rescale_x, plot_type = plot_type, par_name = par_name, dots_prior))

      if(plot_type == "ggplot"){
        plot <- plot + .geom_prior.weightfunction(plot_data, rescale_x = rescale_x, ...)
      }else{
        .lines.prior.weightfunction(plot_data, rescale_x = rescale_x, ...)
      }

    }else{

      # plot just posterior otherwise
      plot <- .plot.prior.weightfunction(NULL, plot_data = plot_data, plot_type = plot_type, rescale_x = rescale_x, par_name = par_name, ...)

    }

  }else if(parameter == "PETPEESE"){

    plot_data <- .plot_data_samples.PETPEESE(samples, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant, n_points = n_points,
                                             transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)

    # add priors, if requested
    if(prior){

      # TODO: a bit of a hack - removing priors that were added as a fill for sampling
      if(!is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
        prior_list_PET   <- attr(samples[["PET"]],   "prior_list")
        prior_list_PEESE <- attr(samples[["PEESE"]], "prior_list")
        prior_fill       <- seq_along(prior_list_PET)[!sapply(prior_list_PET, is.prior.PET) & !sapply(prior_list_PEESE, is.prior.PEESE)]
        prior_list       <- c(prior_list_PET[sapply(prior_list_PET, is.prior.PET)], prior_list_PEESE[sapply(prior_list_PEESE, is.prior.PEESE)],
                              prior_list_PET[prior_fill])
      }else if(is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
        prior_list <- attr(samples[["PEESE"]], "prior_list")
      }else if(!is.null(samples[["PET"]]) & is.null(samples[["PEESE"]])){
        prior_list <- attr(samples[["PET"]], "prior_list")
      }
      #

      plot_data_prior <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                  n_points = n_points, n_samples = n_samples,
                                                  transformation = transformation, transformation_arguments = transformation_arguments,
                                                  transformation_settings = transformation_settings)

      # transplant common xlim and ylim
      plot_data_joined <- list(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      ylim <- range(as.vector(sapply(plot_data_joined, attr, which = "y_range")))
      attr(plot_data_prior, "x_range") <- xlim
      attr(plot_data_prior, "y_range") <- ylim
      dots_prior <- .transfer_dots(dots_prior, ...)

      plot <- do.call(.plot.prior.PETPEESE, c(x = list(prior_list), plot_data = list(plot_data_prior), plot_type = plot_type, par_name = par_name, dots_prior))

      if(plot_type == "ggplot"){
        plot <- plot + .geom_prior.PETPEESE(plot_data, ...)
      }else{
        .lines.prior.PETPEESE(plot_data, ...)
      }

    }else{

      # plot just posterior otherwise
      plot <- .plot.prior.PETPEESE(NULL, plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)

    }


  }else{

    plot_data <- .plot_data_samples.simple(samples, parameter = parameter, n_points = n_points,
                                           transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)

    # add priors, if requested
    if(prior){
      plot_data_prior <- .plot_data_prior_list.simple(attr(samples[[parameter]], "prior_list"), x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings)

      # transplant common xlim and ylim
      plot_data_joined <- c(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      attr(plot_data_prior[[1]], "x_range") <- xlim

      if(any(sapply(plot_data_joined, inherits, what = "density.prior.simple")) & any(sapply(plot_data_joined, inherits, what = "density.prior.point"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.simple")], attr, which = "y_range")))
        ylim2 <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.point")],  attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.simple"))]], "y_range") <- ylim
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.point"))]],  "y_range") <- ylim2
      }else if(any(sapply(plot_data_joined, inherits, what = "density.prior.simple"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.simple")], attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.simple"))]], "y_range") <- ylim
      }else if(any(sapply(plot_data_joined, inherits, what = "density.prior.point"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.point")],  attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.point"))]], "y_range") <- ylim
      }

      scale_y2   <- .get_scale_y2(plot_data_prior, ...)
      dots_prior <- .transfer_dots(dots_prior, ...)
      plot       <- do.call(.plot_prior_list.both, c(plot_data = list(plot_data_prior), plot_type = plot_type, par_name = par_name, dots_prior))

      if(plot_type == "ggplot"){
        for(i in seq_along(plot_data)){
          if(inherits(plot_data[[i]], what = "density.prior.simple")){
            plot <- plot + .geom_prior.simple(plot_data[[i]], ...)
          }else if(inherits(plot_data[[i]], what = "density.prior.point")){
            plot <- plot + .geom_prior.point(plot_data[[i]], scale_y2 = scale_y2, ...)
          }
        }
      }else{
        for(i in seq_along(plot_data)){
          if(inherits(plot_data[[i]], what = "density.prior.simple")){
            .lines.prior.simple(plot_data[[i]], ...)
          }else if(inherits(plot_data[[i]], what = "density.prior.point")){
            .lines.prior.point(plot_data[[i]], scale_y2 = scale_y2, ...)
          }
        }
      }

    }else{

      # plot just posterior otherwise
      plot <- .plot_prior_list.both(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)

    }


  }


  if(plot_type == "ggplot"){
    return(plot)
  }else{
    return(invisible())
  }
}

.plot_data_samples.simple         <- function(samples, parameter, n_points, transformation, transformation_arguments, transformation_settings){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL

  # extract the relevant data
  samples    <- samples[[parameter]]
  prior_list <- attr(samples, "prior_list")
  prior_list <- prior_list

  # deal with spikes
  if(any(sapply(prior_list, is.prior.point))){

    priors_point_map <- cbind(
      sapply(prior_list[sapply(prior_list, is.prior.point)], function(p)p$parameters[["location"]]),
      which(sapply(prior_list, is.prior.point))
    )

    # this should hopefully prevent computer precision issues
    samples_points <- priors_point_map[,1][attr(samples, "models_ind")[attr(samples, "models_ind") %in% priors_point_map[,2]]]
    points_table   <- table(samples_points) / length(samples)

    x_points <- as.numeric(names(points_table))
    y_points <- as.numeric(points_table)

  }

  # deal with the densities
  if(any(!sapply(prior_list, is.prior.point))){

    samples_density   <- samples[attr(samples, "models_ind") %in% which(!sapply(prior_list, is.prior.point))]

    if(length(samples_density) > 0){

      args <- list(x = samples_density, n = n_points)

      # set the endpoints for possible truncation
      prior_list_simple <- prior_list[!sapply(prior_list, is.prior.point)]
      prior_list_simple_lower <- min(sapply(prior_list_simple, function(p) p$truncation[["lower"]]))
      prior_list_simple_upper <- max(sapply(prior_list_simple, function(p) p$truncation[["upper"]]))
      if(!is.infinite(prior_list_simple_lower)){
        args <- c(args, from = prior_list_simple_lower)
      }
      if(!is.infinite(prior_list_simple_upper)){
        args <- c(args, to = prior_list_simple_upper)
      }

      # get the density estimate
      density_continuous <- do.call(stats::density, args)
      x_den    <- density_continuous$x
      y_den    <- density_continuous$y * (length(samples_density) / length(samples))

      # check for truncation
      if(isTRUE(all.equal(prior_list_simple_lower, x_den[1])) | prior_list_simple_lower >= x_den[1]){
        y_den <- c(0, y_den)
        x_den <- c(x_den[1], x_den)
      }
      if(isTRUE(all.equal(prior_list_simple_upper, x_den[length(x_den)])) | prior_list_simple_upper <= x_den[length(x_den)]){
        y_den <- c(y_den, 0)
        x_den <- c(x_den, x_den[length(x_den)])
      }

      # apply transformations
      if(!is.null(transformation)){
        x_den   <- .density.prior_transformation_x(x_den,   transformation, transformation_arguments)
        y_den   <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
        samples_density <- .density.prior_transformation_x(samples_density,   transformation, transformation_arguments)
      }

    }
  }


  # create the output object
  out <- list()

  # add continuous densities
  if(!is.null(y_den)){
    out_den    <- list(
      call    = call("density", "mixed samples"),
      bw      = NULL,
      n       = n_points,
      x       = x_den,
      y       = y_den,
      samples = samples_density
    )

    class(out_den) <- c("density", "density.prior", "density.prior.simple")
    attr(out_den, "x_range") <- range(x_den)
    attr(out_den, "y_range") <- c(0, max(y_den))

    out[["density"]] <- out_den
  }

  # add spikes
  if(!is.null(y_points)){
    for(i in seq_along(y_points)){
      temp_points <- list(
        call    = call("density", paste0("point", i)),
        bw      = NULL,
        n       = n_points,
        x       = x_points[i],
        y       = y_points[i],
        samples = NULL
      )

      class(temp_points) <- c("density", "density.prior", "density.prior.point")
      attr(temp_points, "x_range") <- range(x_points[i])
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}
.plot_data_samples.PETPEESE       <- function(samples, x_seq, x_range, x_range_quant, n_points, transformation, transformation_arguments, transformation_settings){

  check_list(samples, "samples")
  if(is.null(samples[["PET"]]) & is.null(samples[["PEESE"]]))
    stop("At least one 'PET' or 'PEESE' model needs to be specified.")

  # get the samples
  if(!is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
    if(!all(attr(samples[["PET"]], "models_ind") == attr(samples[["PEESE"]], "models_ind")))
      stop("non-matching dimensions")
    samples <- cbind(samples[["PET"]], samples[["PEESE"]])
  }else if(!is.null(samples[["PET"]])){
    samples <- cbind(rep(0, length(samples[["PEESE"]])), samples[["PEESE"]])
  }else if(is.null(samples[["PEESE"]])){
    samples <- cbind(samples[["PET"]], rep(0, length(samples[["PET"]])))
  }

  # get the plotting range
  if(is.null(x_range)){
    x_range <- c(0, 1)
  }
  if(is.null(x_seq)){
    x_seq   <- seq(x_range[1], x_range[2], length.out = n_points)
  }


  # transform the parameter if requested
  if(!is.null(transformation)){
    samples <- .density.prior_transformation_x(samples,   transformation, transformation_arguments)
  }


  # compute PET-PEESE
  x_sam  <- matrix(samples[,1], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq, nrow = length(samples), ncol = length(x_seq), byrow = TRUE) + matrix(samples[,2], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq^2, nrow = length(samples), ncol = length(x_seq), byrow = TRUE)
  x_med  <- apply(x_sam, 2, stats::quantile, prob = .500)
  x_lCI  <- apply(x_sam, 2, stats::quantile, prob = .025)
  x_uCI  <- apply(x_sam, 2, stats::quantile, prob = .975)


  out <- list(
    call    = call("density", "PET-PEESE list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = x_med,
    y_lCI   = x_lCI,
    y_uCI   = x_uCI,
    samples = x_sam
  )


  class(out) <- c("density", "density.prior", "density.prior.PETPEESE")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- range(x_med)

  return(out)
}
.plot_data_samples.weightfunction <- function(samples, x_seq, x_range, x_range_quant, n_points){

  check_list(samples, "samples", check_names = "omega", allow_other = TRUE)
  samples    <- samples[["omega"]]
  prior_list <- attr(samples, "prior_list")

  # get the plotting range
  if(is.null(x_range)){
    x_range <- c(0, 1)
  }
  if(is.null(x_seq)){
    x_seq   <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  # merge the samples
  omega_mapping <- weightfunctions_mapping(prior_list[sapply(prior_list, is.prior.weightfunction)])
  omega_cuts    <- weightfunctions_mapping(prior_list[sapply(prior_list, is.prior.weightfunction)], cuts_only = TRUE)

  x_lCI  <- apply(samples, 2, stats::quantile, probs = .025)
  x_uCI  <- apply(samples, 2, stats::quantile, probs = .975)
  x_mean <- apply(samples, 2, mean)

  x_seq     <- omega_cuts
  x_seq_rep <- c(1, sort(rep(2:(length(x_seq)-1), 2)) ,length(x_seq))
  x_val_rep <- sort(rep(1:(length(x_seq)-1), 2))


  out <- list(
    call    = call("density", "weightfunction list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq[x_seq_rep],
    y       = x_mean[x_val_rep],
    y_lCI   = x_lCI[x_val_rep],
    y_uCI   = x_uCI[x_val_rep],
    samples = samples
  )


  class(out) <- c("density", "density.prior", "density.prior.weightfunction")
  attr(out, "x_range") <- c(0, 1)
  attr(out, "y_range") <- c(0, 1)

  return(out)
}
