#' @title Plot a list of prior distributions
#'
#' @param prior_list list of prior distributions
#' @param prior_list_mu list of priors for the mu parameter
#' required when plotting PET-PEESE
#' @param ... additional arguments
#' @inheritParams density.prior
#' @inheritParams plot.prior
#'
#' @return \code{plot_prior_list} returns either \code{NULL} or
#' an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
#'
#' @seealso [prior()] [lines_prior_list()]  [geom_prior_list()]
#' @export
plot_prior_list <- function(prior_list, plot_type = "base",
                            x_seq = NULL, xlim = NULL, x_range_quant = NULL, n_points = 500,
                            n_samples = 10000, force_samples = FALSE,
                            transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                            rescale_x = FALSE, par_name = NULL, prior_list_mu = NULL, ...){

  # TODO: add plots for individual parameters for weightfunction and PET-PEESE
  individual = FALSE
  show_figures = if(individual) 1 else NULL

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  # check that there is no mixing of PET-PEESE and weightfunctions
  if(any(sapply(prior_list, is.prior.weightfunction)) & (any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))))
    stop("weightfunction and PET-PEESE priors cannot be mixed within a 'prior_list'.")


  # get the plotting type
  if(any(sapply(prior_list, is.prior.weightfunction))){
    prior_type <- "weightfunction"
  }else if(any(sapply(prior_list, is.prior.PET)) | any(sapply(prior_list, is.prior.PEESE))){
    prior_type <- "PETPEESE"
  }else if(any(sapply(prior_list, is.prior.orthonormal))){
    prior_type <- "orthonormal"
  }else if(any(sapply(prior_list, is.prior.meandif))){
    prior_type <- "meandif"
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
    }else if(prior_type %in% c("simple", "orthonormal", "meandif")){
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
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu)
    plot <- .plot.prior.PETPEESE(prior_list, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)

  }else if(prior_type %in% c("simple", "orthonormal", "meandif")){

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

.plot_prior_list.both             <- function(plot_data, plot_type, par_name = NULL, scale_y2 = NULL, add = FALSE, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- range(as.vector(sapply(plot_data, attr, which = "x_range")))

  main      <- ""
  xlab      <- if(!is.null(par_name)) par_name else ""

  if(is.null(scale_y2)) scale_y2  <- .get_scale_y2(plot_data, dots)

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

    if(!add){
      .plot.prior_empty(type, dots)
    }

    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        args           <- dots
        args$plot_data <- plot_data[[i]]
        do.call(.lines.prior.simple, args)
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        args           <- dots
        args$scale_y2  <- scale_y2
        args$plot_data <- plot_data[[i]]
        do.call(.lines.prior.point, args)
      }
    }
    plot <- list(scale_y2 = scale_y2)

  }else if(plot_type == "ggplot"){

    plot <- list()

    for(i in seq_along(plot_data)){
      if(inherits(plot_data[[i]], what = "density.prior.simple")){
        args           <- dots
        args$plot_data <- plot_data[[i]]
        plot           <- c(plot, do.call(.geom_prior.simple, args))
      }else if(inherits(plot_data[[i]], what = "density.prior.point")){
        args           <- dots
        args$scale_y2  <- scale_y2
        args$plot_data <- plot_data[[i]]
        plot           <- c(plot, do.call(.geom_prior.point, args))
      }
    }

    if(!add){
      plot <- .ggplot.prior_empty(type, dots) + plot
    }

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible(plot))
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot_prior_list.factor           <- function(plot_data, plot_type, par_name = NULL, scale_y2 = NULL, add = FALSE, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- range(as.vector(sapply(plot_data, attr, which = "x_range")))

  main      <- ""
  xlab      <- if(!is.null(par_name)) par_name else ""

  if(is.null(scale_y2)) scale_y2  <- .get_scale_y2(plot_data, dots)

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

  # split on points and factors
  plot_data_points  <- plot_data[sapply(plot_data, inherits, what = "density.prior.point")]
  plot_data_factors <- plot_data[sapply(plot_data, inherits, what = "density.prior.factor")]

  # prepare factor naming & formatting
  level_names <- sapply(plot_data_factors, attr, which = "level_name")
  if(any(grepl("[dif:", level_names, fixed = TRUE))){
    level_names <- substr(level_names, regexpr("[dif:", level_names, fixed = TRUE)[[1]] + 5, regexpr("]", level_names, fixed = TRUE) - 1)
  }else if(any(grepl("[", level_names, fixed = TRUE))){
    level_names <- substr(level_names, regexpr("[", level_names, fixed = TRUE)[[1]] + 1, regexpr("]", level_names, fixed = TRUE) - 1)
  }


  # prepare legend information
  if(!is.null(dots[["legend"]]) && !dots[["legend"]]){
    if(!is.null(dots[["col"]]))      dots[["col"]]      <- rep(dots[["col"]][1], length(level_names))
    if(!is.null(dots[["lty"]]))      dots[["lty"]]      <- rep(dots[["lty"]][1], length(level_names))
    if(!is.null(dots[["linetype"]])) dots[["linetype"]] <- rep(dots[["linetype"]][1], length(level_names))
  }else{
    if(is.null(dots[["col"]]) & (is.null(dots[["lty"]]) | is.null(dots[["linetype"]]))){
      dots$col <- grDevices::palette.colors(n = length(level_names) + 1)[-1]
    }
    if(length(dots[["col"]]) == 1)      dots[["col"]]      <- rep(dots[["col"]],      length(level_names))
    if(length(dots[["lty"]]) == 1)      dots[["lty"]]      <- rep(dots[["lty"]],      length(level_names))
    if(length(dots[["linetype"]]) == 1) dots[["linetype"]] <- rep(dots[["linetype"]], length(level_names))

    if(is.null(dots[["legend"]]))       dots[["legend"]]   <- TRUE
  }


  if(plot_type == "base"){

    if(!add){
      .plot.prior_empty(type, dots)
    }

    # plot points
    for(i in seq_along(plot_data_points)){
      args           <- dots
      args$scale_y2  <- scale_y2
      args$plot_data <- plot_data_points[[i]]
      args$col       <- if(unique(length(dots[["col"]])) > 1) .plot.prior_settings()[["col"]]
      args$lty       <- if(unique(length(dots[["lty"]])) > 1) .plot.prior_settings()[["lty"]]
      do.call(.lines.prior.point, args)
    }

    # plot factor levels
    for(i in seq_along(plot_data_factors)){
      args           <- dots
      args$plot_data <- plot_data_factors[[i]]
      args$level     <- i
      do.call(.lines.prior.factor, args)
    }

    if(dots[["legend"]]){
      graphics::legend(
        if(is.null(dots[["legend_position"]])) "topright" else dots[["legend_position"]],
        legend = level_names,
        col    = if(!is.null(dots[["col"]])) dots[["col"]] else rep(.plot.prior_settings()[["col"]], length(level_names)),
        lty    = if(!is.null(dots[["lty"]])) dots[["lty"]] else rep(.plot.prior_settings()[["lty"]], length(level_names)),
        lwd    = if(!is.null(dots[["lwd"]])) dots[["lwd"]] else rep(.plot.prior_settings()[["lwd"]], length(level_names)),
        bty    = "n")
    }

    plot <- list(scale_y2 = scale_y2)

  }else if(plot_type == "ggplot"){

    plot <- list()

    # plot points
    for(i in seq_along(plot_data_points)){
      args           <- dots
      args$scale_y2  <- scale_y2
      args$plot_data <- plot_data[[i]]
      args$col       <- if(unique(length(dots[["col"]])) > 1)      .plot.prior_settings()[["col"]]
      args$lty       <- if(unique(length(dots[["linetype"]])) > 1) .plot.prior_settings()[["linetype"]]
      plot           <- c(plot, do.call(.geom_prior.point, args))
    }

    # plot factor levels
    plot_data_factors <- data.frame(
      x     = do.call(c, lapply(plot_data_factors, function(x) x$x)),
      y     = do.call(c, lapply(plot_data_factors, function(x) x$y)),
      level = do.call(c, lapply(seq_along(plot_data_factors), function(i) rep(level_names[i], length(plot_data_factors[[i]]$x))))
    )

    args             <- dots
    args$level_names <- level_names
    args$plot_data   <- plot_data_factors

    plot <- c(plot, do.call(.geom_prior.factors, args))


    if(dots[["legend"]]){
      plot <- c(plot, list(ggplot2::theme(
        legend.title    = ggplot2::element_blank(),
        legend.position = if(is.null(dots[["legend_position"]])) "right" else dots[["legend_position"]])))
    }

    if(!add){
      plot <- .ggplot.prior_empty(type, dots) + plot
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

  # join the same priors
  prior_list <- .simplify_prior_list(prior_list)

  prior_weights  <- sapply(prior_list, function(p)p$prior_weights)
  mixing_prop    <- prior_weights / sum(prior_weights)

  prior_list     <- prior_list[round(n_samples * mixing_prop) > 0]
  mixing_prop    <- mixing_prop[round(n_samples * mixing_prop) > 0]

  # replace non-weighfunctions from prior mixture feneration
  if(any(!c(sapply(prior_list, is.prior.weightfunction) | sapply(prior_list, is.prior.none)))){
    for(i in seq_along(prior_list)){
      if(!(is.prior.weightfunction(prior_list[[i]]) | is.prior.none(prior_list[[i]]))){
        prior_list[[i]] <- prior_none(prior_weights = prior_weights[i])
      }
    }
  }

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
  omega_mapping <- weightfunctions_mapping(prior_list)
  omega_cuts    <- weightfunctions_mapping(prior_list, cuts_only = TRUE)

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
                                                 transformation, transformation_arguments, transformation_settings, prior_list_mu){

  # TODO: add dependency on the mu parameter as well
  if(is.null(x_seq)){
    x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  # specify it on the transformed range if requested
  if(transformation_settings & !is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
  }

  prior_weights  <- sapply(prior_list, function(p)p$prior_weights)
  mixing_prop    <- prior_weights / sum(prior_weights)

  prior_list     <- prior_list[round(n_samples * mixing_prop) > 0]
  prior_list_mu  <- prior_list_mu[round(n_samples * mixing_prop) > 0]
  mixing_prop    <- mixing_prop[round(n_samples * mixing_prop) > 0]

  # get the samples
  samples_list <- list()
  for(i in seq_along(prior_list)){
    if(is.prior.PET(prior_list[[i]])){
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), rng(prior_list[[i]], round(n_samples * mixing_prop[i])), rep(0, length = round(n_samples * mixing_prop[i])))
    }else if(is.prior.PEESE(prior_list[[i]])){
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), rep(0, length = round(n_samples * mixing_prop[i])), rng(prior_list[[i]], round(n_samples * mixing_prop[i])))
    }else{
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), matrix(0, nrow = round(n_samples * mixing_prop[i]), ncol = 2))
    }
  }
  samples <- do.call(rbind, samples_list)

  # compute PET-PEESE (mu + PET*se + PEESE*se^2)
  x_sam  <- matrix(samples[,1], nrow = length(samples), ncol = length(x_seq)) +
    matrix(samples[,2], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq,   nrow = length(samples), ncol = length(x_seq), byrow = TRUE) +
    matrix(samples[,3], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq^2, nrow = length(samples), ncol = length(x_seq), byrow = TRUE)

  # transform the PEESE parameter if requested
  if(!is.null(transformation)){
    x_sam <- .density.prior_transformation_x(x_sam, transformation, transformation_arguments)
  }

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

  # dispatching for spike and slab priors
  if(length(prior_list) == 1 && is.prior.spike_and_slab(prior_list[[1]])){

    prior_inclusion   <- .get_spike_and_slab_inclusion(prior_list[[1]])
    prior_variable    <- .get_spike_and_slab_variable(prior_list[[1]])

    if(mean(prior_inclusion) < 1 && mean(prior_inclusion) > 0){
      # create a dummy list for the simple mixture
      prior_null                        <- prior("spike", list(0), prior_weights = 1-mean(prior_inclusion))
      prior_variable[["prior_weights"]] <- mean(prior_inclusion)

      prior_list <- list(
        prior_variable,
        prior_null
      )
    }else if(mean(prior_inclusion) >= 1){
      prior_list <- list(prior_variable)
    }else if(mean(prior_inclusion) <= 0){
      prior_list <- list(prior("spike", list(0)))
    }
  }

  # join the same priors
  prior_list <- .simplify_prior_list(prior_list)

  # get common range to ascertain that all priors are aligned
  if(is.null(x_range)){
    if(!is.null(x_seq)){
      x_range <- range(x_seq)
    }else{
      x_range <- range(as.vector(do.call(rbind, lapply(prior_list, function(p) range(p, if(is.null(x_range_quant)) .range.prior_quantile_default(p) else x_range_quant)))))
    }
  }

  prior_weights  <- sapply(prior_list, function(p)p$prior_weights)
  mixing_prop    <- prior_weights / sum(prior_weights)

  prior_list  <- prior_list[round(n_samples * mixing_prop) > 1]
  mixing_prop <- mixing_prop[round(n_samples * mixing_prop) > 0]

  plot_data <- list()
  for(i in seq_along(prior_list)){
    plot_data[[i]] <- density(prior_list[[i]], x_seq = x_seq, x_range = x_range, x_range_quant = x_range_quant,
                              n_points = n_points, n_samples = round(n_samples * mixing_prop[i]), force_samples = force_samples,
                              transformation = transformation, transformation_arguments = transformation_arguments,
                              transformation_settings = transformation_settings, individual = individual, truncate_end = FALSE)
  }

  # the complete samples are added to each output object
  x_sam    <- NULL
  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL

  for(i in seq_along(plot_data)){

    if(force_samples){
      x_sam <- c(x_sam, plot_data[[i]]$samples)
    }

    # align points and densities
    if(inherits(plot_data[[i]], "density.prior.point")){
      x_points <- c(x_points, plot_data[[i]]$x[plot_data[[i]]$y != 0])
      y_points <- c(y_points, mixing_prop[i])
    }else if(inherits(plot_data[[i]], "density.prior.simple") | inherits(plot_data[[i]], "density.prior.orthonormal") | inherits(plot_data[[i]], "density.prior.meandif")){
      x_den <- rbind(x_den, plot_data[[i]]$x)
      y_den <- rbind(y_den, plot_data[[i]]$y * mixing_prop[i])
    }
  }

  # deal with continuous densities
  if(!is.null(y_den)){
    y_den <- apply(y_den, 2, sum)
    if(any(sapply(1:nrow(x_den), function(i) !isTRUE(all.equal(x_den[1,], x_den[i,])))))
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
        samples = x_sam
      )

      class(temp_points) <- c("density", "density.prior", "density.prior.point")
      attr(temp_points, "x_range") <- range(x_points)
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))

      out[[paste0("points",i)]] <- temp_points
    }
  }

  return(out)
}

.simplify_prior_list <- function(prior_list){

  # as_mixed_priors stores non-listed priors
  # (needs to be kept for marginal_ etc)
  if(is.prior.mixture(prior_list)){
    class(prior_list) <- NULL
  } else if(is.prior(prior_list)){
    prior_list <- list(prior_list)
  }


  # return the input with fewer than 2 priors
  if(length(prior_list) < 2){
    return(prior_list)
  }

  new_prior_list <- prior_list
  for(i in seq_along(new_prior_list)){
    new_prior_list[[i]][["prior_weights"]] <- NULL
  }

  # remove all attributes but names and class
  for(i in seq_along(new_prior_list)){
    attributes(new_prior_list[[i]])[!names(attributes(new_prior_list[[i]])) %in% c("names", "class")] <- NULL
  }

  # remove identical priors
  are_equal <- do.call(rbind, lapply(new_prior_list, function(p) sapply(new_prior_list, identical, y = p)))
  are_equal <- are_equal[!duplicated(are_equal) & apply(are_equal, 1, sum) > 1,,drop = FALSE]

  # return the input with no matches
  if(nrow(are_equal) == 0){
    return(prior_list)
  }

  # find the duplicates and collect prior odds
  prior_weights <- unname(sapply(prior_list, function(p)p$prior_weights))
  to_remove  <- NULL
  for(i in 1:nrow(are_equal)){
    this_ind    <- c(1:ncol(are_equal))[are_equal[i,]]
    this_unique <- this_ind[1]
    prior_weights[this_unique] <- sum(prior_weights[this_ind])
    to_remove   <- c(to_remove, this_ind[-1])
  }

  # return prior odds
  for(i in seq_along(prior_list)){
    prior_list[[i]][["prior_weights"]] <- prior_weights[i]
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
#' @return \code{lines_prior_list} returns \code{NULL}.
#'
#' @seealso [plot_prior_list()] [geom_prior_list()]
#' @rdname lines_prior_list
#' @export
lines_prior_list <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, prior_list_mu = NULL, ...){

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

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    .lines.prior.weightfunction(prior_list, plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu)
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
#' @return \code{geom_prior_list} returns an object of class 'ggplot'.
#'
#' @seealso [plot_prior_list()] [lines_prior_list()]
#' @rdname geom_prior_list
#' @export
geom_prior_list  <- function(prior_list, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 500,
                             n_samples = 10000, force_samples = FALSE,
                             transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                             rescale_x = FALSE, scale_y2 = NULL, prior_list_mu = NULL, ...){

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

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    geom <- .geom_prior.weightfunction(prior_list, plot_data = plot_data, rescale_x = rescale_x, ...)

  }else if(prior_type == "PETPEESE"){

    # use samples (not sure how to provide analytic solution for this yes)
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu)
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
#' @return \code{plot_posterior} returns either \code{NULL} or
#' an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
#'
#' @seealso [prior()] [lines_prior_list()]  [geom_prior_list()]
#' @export
plot_posterior <- function(samples, parameter, plot_type = "base", prior = FALSE,
                           n_points = 1000, n_samples = 10000, force_samples = FALSE,
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
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)

  # deal with bad parameter names for PET-PEESE, weightfunction
  if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) %in% c("weightfunction", "weigthfunction", "omega")){
    parameter <- "omega"
  }else if(tolower(gsub("-", "", gsub("_", "", gsub(".", "", parameter, fixed = TRUE),fixed = TRUE), fixed = TRUE)) == "petpeese"){
    parameter <- "PETPEESE"
  }

  # get the plotting range
  dots <- list(...)
  xlim <- dots[["xlim"]]
  if(is.null(xlim)){
    if(parameter %in% c("omega", "PETPEESE") & !individual){
      xlim      <- c(0, 1)
    }else{
      # use the data range otherwise
      xlim   <- NULL
    }
  }


  if(parameter == "omega"){

    plot_data <- .plot_data_samples.weightfunction(samples, x_seq = NULL, x_range = xlim, x_range_quant = NULL, n_points = n_points)

    # add priors, if requested
    if(prior){

      # extract the correct weightfunction samples
      if(!is.null(samples[[parameter]])){
        prior_list <- attr(samples[[parameter]], "prior_list")
      }else if(!is.null(samples[["bias"]])){
        prior_list <- attr(samples[["bias"]], "prior_list")
      }else{
        stop("No 'omega' or 'bias' samples found.")
      }

      prior_list      <- .simplify_prior_list(prior_list)
      plot_data_prior <- .plot_data_prior_list.weightfunction(prior_list, x_seq = NULL, x_range = xlim, x_range_quant = NULL,
                                                              n_points = n_points, n_samples = n_samples)

      # transplant common xlim and ylim
      plot_data_joined <- list(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      ylim <- range(as.vector(sapply(plot_data_joined, attr, which = "y_range")))
      attr(plot_data_prior, "x_range") <- xlim
      attr(plot_data_prior, "y_range") <- ylim
      dots_prior <- .transfer_dots(dots_prior, ...)

      args           <- dots_prior
      args$x         <- prior_list
      args$plot_data <- plot_data_prior
      args$rescale_x <- rescale_x
      args$plot_type <- plot_type
      args$par_name  <- par_name
      plot           <- do.call(.plot.prior.weightfunction, args)

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

    plot_data <- .plot_data_samples.PETPEESE(samples, x_seq = NULL, x_range = xlim, x_range_quant = NULL, n_points = n_points,
                                             transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)

    # add priors, if requested
    if(prior){

      if(is.null(samples[["mu"]]))
        stop("'mu' samples are required for plotting PET-PEESE.")
      prior_list_mu   <- attr(samples[["mu"]],   "prior_list")

      # TODO: a bit of a hack - removing priors that were added as a fill for sampling
      if(!is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
        prior_list_PET   <- attr(samples[["PET"]],   "prior_list")
        prior_list_PEESE <- attr(samples[["PEESE"]], "prior_list")
        prior_fill       <- seq_along(prior_list_PET)[!sapply(prior_list_PET, is.prior.PET) & !sapply(prior_list_PEESE, is.prior.PEESE)]
        prior_list       <- c(prior_list_PET[sapply(prior_list_PET, is.prior.PET)], prior_list_PEESE[sapply(prior_list_PEESE, is.prior.PEESE)],
                              prior_list_PET[prior_fill])
        prior_list_mu    <- prior_list_mu[c(c(1:length(prior_list_mu))[sapply(prior_list_PET, is.prior.PET)], c(1:length(prior_list_mu))[sapply(prior_list_PEESE, is.prior.PEESE)], c(1:length(prior_list_mu))[prior_fill])]
      }else if(is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
        prior_list <- attr(samples[["PEESE"]], "prior_list")
      }else if(!is.null(samples[["PET"]]) & is.null(samples[["PEESE"]])){
        prior_list <- attr(samples[["PET"]], "prior_list")
      }else{
        stop("Either PET or PEESE samples need to be provided.")
      }

      # cannot simplify prior_list - it would break the dependency with mu

      plot_data_prior <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = NULL, x_range = xlim, x_range_quant = NULL,
                                                  n_points = n_points, n_samples = n_samples,
                                                  transformation = transformation, transformation_arguments = transformation_arguments,
                                                  transformation_settings = transformation_settings, prior_list_mu = prior_list_mu)

      # transplant common xlim and ylim
      plot_data_joined <- list(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      ylim <- range(as.vector(sapply(plot_data_joined, attr, which = "y_range")))
      # make sure y-range does not collapse
      if(all(ylim < .01)){
        ylim <- c(0, 1)
      }
      attr(plot_data_prior, "x_range") <- xlim
      attr(plot_data_prior, "y_range") <- ylim
      dots_prior <- .transfer_dots(dots_prior, ...)

      args           <- dots_prior
      args$x         <- prior_list
      args$plot_data <- plot_data_prior
      args$plot_type <- plot_type
      args$par_name  <- par_name
      plot           <- do.call(.plot.prior.PETPEESE, args)

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

    prior_list  <- attr(samples[[parameter]], "prior_list")
    prior_list  <- .simplify_prior_list(prior_list)

    if(any(sapply(prior_list, is.prior.factor))){
      plot_data <- .plot_data_samples.factor(samples, parameter = parameter, n_points = n_points,
                                             transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)
    }else{
      plot_data <- .plot_data_samples.simple(samples, parameter = parameter, n_points = n_points,
                                             transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)
    }



    # add priors, if requested
    if(prior){

      plot_data_prior <- .plot_data_prior_list.simple(prior_list, x_seq = NULL, x_range = xlim, x_range_quant = NULL,
                                                n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings)

      # transplant common xlim and ylim
      plot_data_joined <- c(plot_data_prior, plot_data)

      xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
      attr(plot_data_prior[[1]], "x_range") <- xlim

      if(any(sapply(plot_data_prior, inherits, what = "density.prior.simple")) & any(sapply(plot_data_prior, inherits, what = "density.prior.point"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.simple")], attr, which = "y_range")))
        ylim2 <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.point")],  attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.simple"))]], "y_range") <- ylim
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.point"))]],  "y_range") <- ylim2
      }else if(any(sapply(plot_data_prior, inherits, what = "density.prior.simple"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.simple")], attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.simple"))]], "y_range") <- ylim
      }else if(any(sapply(plot_data_prior, inherits, what = "density.prior.point"))){
        ylim  <- range(as.vector(sapply(plot_data_joined[sapply(plot_data_joined, inherits, what = "density.prior.point")],  attr, which = "y_range")))
        attr(plot_data_prior[[which.max(sapply(plot_data_prior, inherits, what = "density.prior.point"))]], "y_range") <- ylim
      }

      scale_y2   <- .get_scale_y2(plot_data_prior, ...)
      dots_prior <- .transfer_dots(dots_prior, ...)


      # set the y/x ranges
      for(i in seq_along(plot_data)){
        if(inherits(plot_data[[i]], what = "density.prior.point")){
          attr(plot_data[[i]], which = "y_range") <- if(any(sapply(plot_data_prior, inherits, what = "density.prior.simple")) & any(sapply(plot_data_prior, inherits, what = "density.prior.point"))) ylim2 else ylim
        }else{
          attr(plot_data[[i]], which = "y_range") <- ylim
          attr(plot_data[[i]], which = "x_range") <- xlim
        }
      }

      # plot prior
      args_prior           <- dots_prior
      args_prior$plot_data <- plot_data_prior
      args_prior$plot_type <- plot_type
      args_prior$par_name  <- par_name
      args_prior$scale_y2  <- scale_y2

      plot <- do.call(.plot_prior_list.both, args_prior)


      # plot posterior
      args           <- list(...)
      args$plot_data <- plot_data
      args$plot_type <- plot_type
      args$par_name  <- par_name
      args$scale_y2  <- scale_y2
      args$add       <- TRUE

      if(plot_type == "base"){
        if(any(sapply(prior_list, is.prior.factor))){
          plot <- do.call(.plot_prior_list.factor, args)
        }else{
          plot <- do.call(.plot_prior_list.both, args)
        }
      }else if(plot_type == "ggplot"){
        if(any(sapply(prior_list, is.prior.factor))){
          plot <- plot + do.call(.plot_prior_list.factor, args)
        }else{
          plot <- plot + do.call(.plot_prior_list.both, args)
        }
      }


    }else{

      # plot just posterior otherwise
      if(any(sapply(prior_list, is.prior.factor))){
        plot <- .plot_prior_list.factor(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)
      }else{
        plot <- .plot_prior_list.both(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)
      }

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
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  # deal with spikes
  if(any(sapply(prior_list, is.prior.point))){

    # aggregate samples across spikes
    spikes_simplified <- .simplify_spike_samples(samples, prior_list)

    if(nrow(spikes_simplified) > 0){
      x_points <- spikes_simplified[,"location"]
      y_points <- spikes_simplified[,"probability"]
    }else{
      x_points <- NULL
      y_points <- NULL
    }

    # apply transformations
    if(!is.null(transformation)){
      x_points <- .density.prior_transformation_x(x_points, transformation, transformation_arguments)
    }
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
      if(!is.infinite(prior_list_simple_lower)){ # adding a small number for possible transformations (0 -> -Inf)
        args <- c(args, from = prior_list_simple_lower + if(!is.null(transformation)) 1e-5 else 0)
      }
      if(!is.infinite(prior_list_simple_upper)){
        args <- c(args, to = prior_list_simple_upper   - if(!is.null(transformation)) 1e-5 else 0)
      }

      # get the density estimate
      density_continuous <- do.call(stats::density, args)
      x_den    <- density_continuous$x
      y_den    <- density_continuous$y * (length(samples_density) / length(samples))

      # check for truncation
      if(!is.null(transformation)){
        if(isTRUE(all.equal(prior_list_simple_lower + 1e-5, x_den[1])) | prior_list_simple_lower + 1e-5 >= x_den[1]){
          y_den <- c(0, y_den)
          x_den <- c(x_den[1], x_den)
        }
        if(isTRUE(all.equal(prior_list_simple_upper - 1e-5, x_den[length(x_den)])) | prior_list_simple_upper + 1e-5 <= x_den[length(x_den)]){
          y_den <- c(y_den, 0)
          x_den <- c(x_den, x_den[length(x_den)])
        }
      }else{
        if(isTRUE(all.equal(prior_list_simple_lower, x_den[1])) | prior_list_simple_lower >= x_den[1]){
          y_den <- c(0, y_den)
          x_den <- c(x_den[1], x_den)
        }
        if(isTRUE(all.equal(prior_list_simple_upper, x_den[length(x_den)])) | prior_list_simple_upper <= x_den[length(x_den)]){
          y_den <- c(y_den, 0)
          x_den <- c(x_den, x_den[length(x_den)])
        }
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
  if(is.null(samples[["mu"]]))
    stop("'mu' samples need to be present.")

  # get the samples
  if(!is.null(samples[["PET"]]) & !is.null(samples[["PEESE"]])){
    if(!all(attr(samples[["PET"]], "models_ind") == attr(samples[["PEESE"]], "models_ind")))
      stop("non-matching dimensions")
    samples <- cbind(samples[["mu"]], samples[["PET"]], samples[["PEESE"]])
  }else if(is.null(samples[["PET"]])){
    samples <- cbind(samples[["mu"]], rep(0, length(samples[["PEESE"]])), samples[["PEESE"]])
  }else if(is.null(samples[["PEESE"]])){
    samples <- cbind(samples[["mu"]], samples[["PET"]], rep(0, length(samples[["PET"]])))
  }

  # get the plotting range
  if(is.null(x_range)){
    x_range <- c(0, 1)
  }
  if(is.null(x_seq)){
    x_seq   <- seq(x_range[1], x_range[2], length.out = n_points)
  }


  # compute PET-PEESE (mu + PET*se + PEESE*se^2)
  x_sam  <- matrix(samples[,1], nrow = length(samples), ncol = length(x_seq)) +
    matrix(samples[,2], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq, nrow = length(samples), ncol = length(x_seq), byrow = TRUE) +
    matrix(samples[,3], nrow = length(samples), ncol = length(x_seq)) * matrix(x_seq^2, nrow = length(samples), ncol = length(x_seq), byrow = TRUE)

  # transform the parameter if requested
  if(!is.null(transformation)){
    x_sam <- .density.prior_transformation_x(x_sam, transformation, transformation_arguments)
  }

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
  if(!is.null(samples[["omega"]])){
    samples <- samples[["omega"]]
  }else if(!is.null(samples[["bias"]])){
    samples <- samples[["bias"]]
  }else{
    stop("No 'omega' or 'bias' samples found.")
  }

  prior_list <- attr(samples, "prior_list")
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

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
.plot_data_samples.factor         <- function(samples, parameter, n_points, transformation, transformation_arguments, transformation_settings){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL

  # transform & extract the relevant data
  prior_list <- attr(samples[[parameter]], "prior_list")
  if (!(is.prior.mixture(prior_list) || is.prior.spike_and_slab(prior_list)) && is.prior(prior_list))
    prior_list <- list(prior_list)

  if(any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))){
    samples  <- transform_factor_samples(samples)
    if(!is.null(transformation)){
      message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the meandif/orthonormal contrasts to the differences from the mean.")
    }
  }

  samples    <- samples[[parameter]]

  # create the output object
  out <- list()

  # deal with spikes
  if(any(sapply(prior_list, is.prior.point))){

    # aggregate samples across spikes
    spikes_simplified <- .simplify_spike_samples(samples, prior_list)

    if(nrow(spikes_simplified) > 0){
      x_points <- spikes_simplified[,"location"]
      y_points <- spikes_simplified[,"probability"]
    }else{
      x_points <- NULL
      y_points <- NULL
    }

    # apply transformations
    if(!is.null(transformation)){
      x_points <- .density.prior_transformation_x(x_points, transformation, transformation_arguments)
    }

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

  # deal with the densities
  if(any(!sapply(prior_list, is.prior.point))){

    samples_density <- samples[attr(samples, "models_ind") %in% which(!sapply(prior_list, is.prior.point)),,drop=FALSE]

    if(nrow(samples_density) > 0){
      for(i in 1:ncol(samples_density)){

        args <- list(x = samples_density[,i], n = n_points)

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
        y_den    <- density_continuous$y * (length(samples_density[,i]) / length(samples))

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
          samples_density[,i] <- .density.prior_transformation_x(samples_density[,i],   transformation, transformation_arguments)
        }

        out_den    <- list(
          call    = call("density", "mixed samples"),
          bw      = NULL,
          n       = n_points,
          x       = x_den,
          y       = y_den,
          samples = samples_density
        )

        class(out_den) <- c("density", "density.prior", "density.prior.factor", "density.prior.simple")
        attr(out_den, "x_range")    <- range(x_den)
        attr(out_den, "y_range")    <- c(0, max(y_den))
        attr(out_den, "level")      <- i
        attr(out_den, "level_name") <- colnames(samples_density)[i]

        out[[paste0("density", i)]] <- out_den

      }
    }
  }



  return(out)
}


#' @title Plot estimates from models
#'
#' @details Plots prior and posterior estimates of the same parameter
#' across multiple models (prior distributions with orthonormal/meandif contrast
#' are always plotted as differences from the grand mean).
#'
#' @param parameter parameter name to be plotted. Does not support
#' PET-PEESE and weightfunction.
#' @param inference object created by [ensemble_inference] function
#' @param conditional whether conditional models should be displayed
#' @param order list specifying ordering of the models. The first
#' element describes whether the ordering should be \code{"increasing"}
#' or \code{"decreasing"} and the second element describes whether
#' the ordering should be based \code{"model"} order, \code{"estimate"}
#' size, posterior \code{"probability"}, or the inclusion \code{"BF"}.
#' @param ... additional arguments. E.g.:
#' \describe{
#'   \item{\code{"show_updating"}}{whether Bayes factors and change from
#'   prior to posterior odds should be shown on the secondary y-axis}
#'   \item{\code{"show_estimates"}}{whether posterior estimates and 95% CI
#'   should be shown on the secondary y-axis}
#'   \item{\code{"y_axis2"}}{whether the secondary y-axis should be shown}
#' }
#' @inheritParams ensemble_inference
#' @inheritParams plot.prior
#' @inheritParams plot_posterior
#' @inheritParams format_parameter_names
#'
#' @return \code{plot_models} returns either \code{NULL} or
#' an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
#'
#' @seealso [prior()] [lines_prior_list()]  [geom_prior_list()]
#' @export
plot_models <- function(model_list, samples, inference, parameter, plot_type = "base", prior = FALSE, conditional = FALSE,
                        order = NULL,
                        transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                        par_name = NULL, formula_prefix = TRUE, ...){

  # check input
  check_list(model_list, "model_list")
  check_char(parameter, "parameter")
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = "fit_summary", all_objects = TRUE, allow_other = TRUE))
  if(!all(unlist(sapply(model_list, function(m) sapply(attr(m[["fit"]], "prior_list"), function(p) is.prior(p))))))
    stop("model_list:priors must contain 'BayesTools' priors")
  if(!all(sapply(model_list, function(m)inherits(m[["fit_summary"]], what = "BayesTools_runjags_summary"))))
    stop("model_list:fit_summary must contain 'BayesTools' fit_summary")
  check_list(samples, "samples")
  if(any(!sapply(samples, inherits, what = "mixed_posteriors")))
    stop("'samples' must be a be an object generated by 'mix_posteriors' function.")
  check_list(order, "order", allow_NULL = TRUE)
  if(!is.null(order)){
    check_list(order, "order", check_length = 2)
    check_char(order[[1]], allow_values = c("decreasing", "increasing"))
    check_char(order[[2]], allow_values = c("model", "estimate", "probability", "BF"))
  }

  # extract the objects
  models_summary   <- lapply(model_list, function(m)m[["fit_summary"]])
  models_inference <- lapply(model_list, function(m)m[["inference"]])
  total_inference  <- inference[[parameter]]
  total_samples    <- samples[[parameter]]
  prior_list       <- attr(total_samples, "prior_list")

  # deal with factors
  if(inherits(total_samples, "mixed_posteriors.factor")){

    # make sure that the orthonormal/meandif priors are transformed to differences from the mean
    if(attr(total_samples, "orthonormal") | attr(total_samples, "meandif")){

      # transform the samples
      if(ncol(total_samples) == attr(total_samples, "levels")){
        total_samples <- transform_factor_samples(samples[parameter])[[parameter]]
      }
      # transform the model summaries
      if(!any(sapply(models_summary, function(m) all(colnames(total_samples) %in% attr(m, "parameters"))))){
        models_summary <- lapply(model_list, function(m) runjags_estimates_table(m[["fit"]], transform_factors = TRUE))
      }

    }

    parameter <- colnames(total_samples)
  }

  # prepare nice parameter names
  if(is.null(par_name)){
    par_name  <- format_parameter_names(parameter, attr(total_samples, "formula_parameter"), formula_prefix = formula_prefix)
  }

  plot <- list()
  for(i in seq_along(parameter)){
    plot[[i]] <- .plot_models.simple(models_summary = models_summary, models_inference = models_inference, total_inference = total_inference, total_samples = total_samples,
                                     prior_list = prior_list, parameter = parameter[i], par_name = par_name[i],
                                     plot_type = plot_type, prior = prior, conditional = conditional, order = order,
                                     transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings, ...)
  }


  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    if(length(parameter) == 1){
      plot <- plot[[i]]
    }
    return(plot)
  }
}

.plot_models_data_prior     <- function(prior_list, models_inference){
  return(data.frame(
    model      = 1:length(prior_list),
    y          = sapply(prior_list, mean),
    y_lCI      = sapply(prior_list, mquant, .025),
    y_uCI      = sapply(prior_list, mquant, .975),
    prior_prob = sapply(models_inference, function(m)m[["prior_prob"]]),
    post_prob  = sapply(models_inference, function(m)m[["post_prob"]]),
    BF         = sapply(models_inference, function(m)m[["inclusion_BF"]])
  ))
}
.plot_models_data_posterior <- function(models_summary, parameter, prior_list, models_inference){
  return(data.frame(
    model      = 1:length(models_summary),
    y          = sapply(1:length(models_summary), function(i){
      if(any(attr(models_summary[[i]], "parameters") == parameter)){
        return(models_summary[[i]][attr(models_summary[[i]], "parameters") == parameter, "Mean"])
      }else if(is.prior.point(prior_list[[i]])){
        return(prior_list[[i]]$parameters[["location"]])
      }else{
        stop(paste0("Posterior distribution summary for '", parameter, "' is not available."))
      }
    }),
    y_lCI      = sapply(1:length(models_summary), function(i){
      if(any(attr(models_summary[[i]], "parameters") == parameter)){
        return(models_summary[[i]][attr(models_summary[[i]], "parameters") == parameter, "lCI"])
      }else if(is.prior.point(prior_list[[i]])){
        return(prior_list[[i]]$parameters[["location"]])
      }else{
        stop(paste0("Posterior distribution for '", parameter, "' is not available."))
      }
    }),
    y_uCI      = sapply(1:length(models_summary), function(i){
      if(any(attr(models_summary[[i]], "parameters") == parameter)){
        return(models_summary[[i]][attr(models_summary[[i]], "parameters") == parameter, "uCI"])
      }else if(is.prior.point(prior_list[[i]])){
        return(prior_list[[i]]$parameters[["location"]])
      }else{
        stop(paste0("Posterior distribution for '", parameter, "' is not available."))
      }
    }),
    prior_prob = sapply(models_inference, function(m) m[["prior_prob"]]),
    post_prob  = sapply(models_inference, function(m) m[["post_prob"]]),
    BF         = sapply(models_inference, function(m) m[["inclusion_BF"]])
  ))
}
.plot_models.simple         <- function(models_summary, models_inference, total_inference, total_samples,
                                        prior_list, parameter, par_name,
                                        plot_type, prior, conditional, order,
                                        transformation, transformation_arguments, transformation_settings, ...){

  # create a table with results
  prior_data     <- .plot_models_data_prior(prior_list, models_inference)
  posterior_data <- .plot_models_data_posterior(models_summary, parameter, prior_list, models_inference)

  # remove null models if requested (assuming that the overall estimate is already supplied accordingly)
  if(conditional){
    prior_data     <- prior_data[!attr(total_inference, "is_null"),]
    posterior_data <- posterior_data[!attr(total_inference, "is_null"),]
  }

  # apply ordering
  if(!is.null(order)){
    prior_data <- switch(
      order[[2]],
      "model"       = prior_data[order(posterior_data$model,      decreasing = order[[1]] == "decreasing"),],
      "estimate"    = prior_data[order(posterior_data$y,          decreasing = order[[1]] == "decreasing"),],
      "probability" = prior_data[order(posterior_data$post_prob,  decreasing = order[[1]] == "decreasing"),],
      "BF"          = prior_data[order(posterior_data$BF,         decreasing = order[[1]] == "decreasing"),]
    )
    posterior_data <- switch(
      order[[2]],
      "model"       = posterior_data[order(posterior_data$model,      decreasing = order[[1]] == "decreasing"),],
      "estimate"    = posterior_data[order(posterior_data$y,          decreasing = order[[1]] == "decreasing"),],
      "probability" = posterior_data[order(posterior_data$post_prob,  decreasing = order[[1]] == "decreasing"),],
      "BF"          = posterior_data[order(posterior_data$BF,         decreasing = order[[1]] == "decreasing"),]
    )
  }

  # compute overall estimate
  overal_mean <- mean(total_samples)
  overal_lCI  <- unname(stats::quantile(total_samples, .025))
  overal_uCI  <- unname(stats::quantile(total_samples, .975))
  vertical_0  <- 0

  # apply transformations
  if(!is.null(transformation)){
    prior_data[,c("y", "y_lCI","y_uCI")]     <- .density.prior_transformation_x(prior_data[,c("y", "y_lCI","y_uCI")],   transformation, transformation_arguments)
    posterior_data[,c("y", "y_lCI","y_uCI")] <- .density.prior_transformation_x(posterior_data[,c("y", "y_lCI","y_uCI")],   transformation, transformation_arguments)
    overal_mean <- .density.prior_transformation_x(overal_mean, transformation, transformation_arguments)
    overal_lCI  <- .density.prior_transformation_x(overal_lCI,  transformation, transformation_arguments)
    overal_uCI  <- .density.prior_transformation_x(overal_uCI , transformation, transformation_arguments)
    vertical_0  <- .density.prior_transformation_x(vertical_0 , transformation, transformation_arguments)
  }

  # add names and locations
  prior_data$x     <- 3:(nrow(prior_data) + 2)     + .25
  posterior_data$x <- 3:(nrow(posterior_data) + 2)

  posterior_data$y_labels2  <- paste0(
    format(round(posterior_data$y, 2), nsmall = 2),
    " [", format(round(posterior_data$y_lCI, 2), nsmall = 2), ", ",
    format(round(posterior_data$y_uCI, 2), nsmall = 2), "]")
  prior_data$y_labels2      <- paste0(
    "BF = ", format(round(prior_data$BF, 2), nsmall = 2),
    " [", format(round(prior_data$prior_prob, 2), nsmall = 2), " -> ", format(round(prior_data$post_prob, 2), nsmall = 2),"]")
  overal_y_label  <- "Model-Averaged"
  overal_y_label2 <- paste0(
    format(round(overal_mean, 2), nsmall = 2), " [", format(round(overal_lCI, 2), nsmall = 2), ", ",
    format(round(overal_uCI, 2), nsmall = 2), "]")


  # set the plotting values
  dots      <- list(...)

  show_updating  <- is.null(dots[["show_updating"]])  || dots[["show_updating"]]
  show_estimates <- is.null(dots[["show_estimates"]]) || dots[["show_estimates"]]

  y_at      <- c(1,               posterior_data$x)
  y_labels  <- c(overal_y_label,  paste0("Model ", posterior_data$model))
  y_at2     <- c(1,               posterior_data$x[show_estimates],         prior_data$x[show_updating])
  y_labels2 <- c(overal_y_label2, posterior_data$y_labels2[show_estimates], prior_data$y_labels2[show_updating])
  ylim      <- c(0, max(posterior_data$x) + 1)


  if(!is.null(dots[["xlim"]])){
    xlim     <- dots[["xlim"]]
    x_labels <- pretty(xlim)
  }else{
    x_labels  <- pretty(range(c(posterior_data$y_lCI, posterior_data$y_uCI, overal_lCI, overal_uCI, if(prior) c(prior_data$y_lCI, prior_data$y_uCI))))
    xlim      <- range(x_labels)
  }
  if(!is.null(dots[["xlab"]])){
    xlab <- dots[["xlab"]]
  }else{
    xlab <- if(is.null(par_name)) parameter else bquote(.(par_name))
  }
  if(!is.null(dots[["cex"]])){
    cex  <- dots[["cex"]]
  }else{
    cex  <- 4
  }
  if(!is.null(dots[["col2"]])){
    col2  <- dots[["col2"]]
  }else{
    col2  <- "grey80"
  }
  if(!is.null(dots[["col"]])){
    col   <- dots[["col"]]
  }else{
    col   <- "black"
  }


  ### do the plotting
  if(plot_type == "base"){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))

    # set up margins
    if(length(dots[["mar"]]) == 0){
      graphics::par(mar = c(4, max(nchar(y_labels)) * 2/3, 0, max(nchar(y_labels2)) * 1/2))
    }else{
      graphics::par(mar = dots[["mar"]])
    }

    graphics::plot(NA, bty = "n", las = 1, xlab = xlab, ylab = "", main = "", yaxt = "n", ylim = ylim, xlim = xlim)
    graphics::axis(2, at = y_at,  labels = y_labels,  las = 1, col = NA)
    if(is.null(dots[["y_axis2"]]) || (!is.null(dots[["y_axis2"]]) && dots[["y_axis2"]])){
      graphics::axis(4, at = y_at2, labels = y_labels2, las = 1, col = NA, hadj = 0)
    }
    graphics::abline(v = vertical_0, lty = 3)

    graphics::arrows(
      x0 = posterior_data$y_lCI[posterior_data$y_uCI - posterior_data$y_lCI > 0],
      x1 = posterior_data$y_uCI[posterior_data$y_uCI - posterior_data$y_lCI > 0],
      y0 = posterior_data$x[posterior_data$y_uCI - posterior_data$y_lCI > 0],
      code = 3, angle = 90, length = 0.1, col = col)
    graphics::points(posterior_data$y, posterior_data$x, pch = 15, cex = cex*prior_data$post_prob, col = col)

    if(prior){
      graphics::arrows(
        x0 = prior_data$y_lCI[prior_data$y_uCI - prior_data$y_lCI > 0],
        x1 = prior_data$y_uCI[prior_data$y_uCI - prior_data$y_lCI > 0],
        y0 = prior_data$x[prior_data$y_uCI - prior_data$y_lCI > 0],
        code = 3, angle = 90, length = 0.1, col = col2)
      graphics::points(prior_data$y, prior_data$x, pch = 15, cex = cex*prior_data$prior_prob, col = col2)
    }


    graphics::polygon(
      x = c(overal_lCI, overal_mean , overal_uCI, overal_mean),
      y = c(1, 1.25, 1, 0.75),
      col = col
    )

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()

    # add the studies
    plot <- plot +ggplot2::geom_errorbarh(
      mapping = ggplot2::aes(
        xmin   = posterior_data[posterior_data$y_uCI - posterior_data$y_lCI > 0,]$y_lCI,
        xmax   = posterior_data[posterior_data$y_uCI - posterior_data$y_lCI > 0,]$y_uCI,
        y      = posterior_data[posterior_data$y_uCI - posterior_data$y_lCI > 0,]$x),
      color   = col,
      height  = .25)
    plot <- plot +ggplot2::geom_point(
      mapping = ggplot2::aes(
        x    = posterior_data$y,
        y    = posterior_data$x,
        size = cex * posterior_data$post_prob),
      col   = col,
      shape = 15)

    if(prior){
      plot <- plot +ggplot2::geom_errorbarh(
        mapping = ggplot2::aes(
          xmin   = prior_data[prior_data$y_uCI - prior_data$y_lCI > 0,]$y_lCI,
          xmax   = prior_data[prior_data$y_uCI - prior_data$y_lCI > 0,]$y_uCI,
          y      = prior_data[prior_data$y_uCI - prior_data$y_lCI > 0,]$x),
        color   = col2,
        height  = .25)
      plot <- plot +ggplot2::geom_point(
        mapping = ggplot2::aes(
          x    = prior_data$y,
          y    = prior_data$x,
          size = cex * prior_data$prior_prob),
        shape = 15,
        col   = col2)
    }

    # add the overall estimate
    plot <- plot + ggplot2::geom_polygon(
      mapping = ggplot2::aes(
        x = c(overal_lCI, overal_mean , overal_uCI, overal_mean),
        y = c(1, 1.25, 1, 0.75)),
      fill = col)

    # add the vertical line
    plot <- plot + ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = c(vertical_0, vertical_0),
        y = ylim),
      linetype = "dotted")

    # add all the other stuff
    if(is.null(dots[["y_axis2"]]) || (!is.null(dots[["y_axis2"]]) && dots[["y_axis2"]])){
      plot <- plot + ggplot2::scale_y_continuous(
        name = "", breaks = y_at, labels = y_labels, limits = ylim,
        sec.axis = ggplot2::sec_axis( ~ ., breaks = y_at2, labels = y_labels2))
      attr(plot, "sec_axis") <- TRUE
    }else{
      plot <- plot + ggplot2::scale_y_continuous(
        name = "", breaks = y_at, labels = y_labels, limits = ylim)
      attr(plot, "sec_axis") <- FALSE
    }
    plot <- plot + ggplot2::scale_x_continuous(
      name = xlab, breaks = x_labels, labels = x_labels, limits = xlim)
    plot <- plot + ggplot2::theme(
      axis.title.y      = ggplot2::element_blank(),
      axis.line.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_text(hjust = 0, color = "black"),
      axis.text.y.right = ggplot2::element_text(hjust = 1, color = "black"),
      legend.position   = "none")

  }

  return(plot)
}

.simplify_spike_samples <- function(samples, prior_list){

  # Check if we're dealing with spike_and_slab (0-based indexing) or mixture (1-based indexing)
  is_spike_and_slab <- is.prior.spike_and_slab(prior_list)
  
  # aggregate for each spike
  priors_point_map <- data.frame(do.call(rbind, lapply(seq_along(prior_list), function(i) {
    if(is.prior.point(prior_list[[i]])){
      if(is_spike_and_slab) {
        # For spike_and_slab: dbern() generates 0 (null) and 1 (alternative)
        # Components are stored as ["alternative", "null"], so:
        # - point component at index 2 (null) maps to JAGS index 0
        # - any other point component would be unusual but map accordingly
        component_name <- attr(prior_list[[i]], "component")
        model_index <- if(component_name == "null") 0 else 1
      } else {
        # For mixture: dcat() generates 1, 2, 3... so index i maps to JAGS index i
        model_index <- i
      }
      c("location" = prior_list[[i]]$parameters[["location"]], "frequency" = sum(attr(samples, "models_ind") == model_index))
    }
  })))


  # return the input with fewer than 2 inputs
  if(nrow(priors_point_map) < 2){
    spike_probability = data.frame(cbind(
      "location"    = priors_point_map[, "location"],
      "probability" = priors_point_map[, "frequency"] / if(!is.matrix(samples)) length(samples) else nrow(samples) ))
    spike_probability <- spike_probability[priors_point_map[, "frequency"] != 0, ]
    return(spike_probability)
  }

  # find unique spikes
  unique_map <- cbind("location" = unique(priors_point_map[, "location"]), "frequency" = 0)

  # collect them
  for(i in 1:nrow(unique_map)){
    unique_map[i, "frequency"] <- sum(priors_point_map[sapply(priors_point_map[, "location"], function(l) isTRUE(all.equal(l, unname(unique_map[i, "location"])))), "frequency"])
  }

  spike_probability = data.frame(cbind(
    "location"    = unique_map[, "location"],
    "probability" = unique_map[, "frequency"] / if(!is.matrix(samples)) length(samples) else nrow(samples) ))
  spike_probability <- spike_probability[unique_map[, "frequency"] != 0, ]

  return(spike_probability)
}


#' @title Plot samples from the marginal posterior distributions
#'
#' @param samples samples from a posterior distribution for a
#' parameter generated by [marginal_inference].
#' @param parameter parameter name to be plotted.
#' @param ... additional arguments
#' @inheritParams plot_posterior
#' @inheritParams density.prior
#' @inheritParams plot.prior
#'
#' @return \code{plot_marginal} returns either \code{NULL} or
#' an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
#'
#' @seealso [prior()] [marginal_inference()]  [plot_posterior()]
#' @export
plot_marginal <- function(samples, parameter, plot_type = "base", prior = FALSE,
                          n_points = 1000,
                          transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                          rescale_x = FALSE, par_name = NULL, dots_prior = list(), ...){

  # check input
  if(any(!sapply(samples, inherits, what = "marginal_posterior")))
    stop("'samples' must be a be an object generated by 'marginal_posterior' function.")
  check_char(parameter, "parameter")
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)


  # get the plotting range
  dots <- list(...)
  xlim <- dots[["xlim"]]


  plot_data <- .plot_data_marginal_samples(samples, parameter = parameter, prior = prior, n_points = n_points,
                                           transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings)



  # add priors, if requested
  if(prior){

    plot_data_prior <- lapply(plot_data, attr, which = "prior")

    # transplant common xlim and ylim
    plot_data_joined <- c(plot_data, plot_data_prior)

    xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "x_range")))
    attr(plot_data_prior[[1]], "x_range") <- xlim

    xlim <- range(as.vector(sapply(plot_data_joined, attr, which = "y_range")))
    attr(plot_data_prior[[1]], "y_range") <- xlim

    dots_prior <- .transfer_dots(dots_prior, ...)


    # plot prior
    args_prior           <- dots_prior
    args_prior$plot_data <- plot_data_prior
    args_prior$plot_type <- plot_type
    args_prior$par_name  <- par_name
    args_prior$hardcode  <- TRUE

    plot <- do.call(.plot_prior_list.factor, args_prior)


    # plot posterior
    args           <- list(...)
    args$plot_data <- plot_data
    args$plot_type <- plot_type
    args$par_name  <- par_name
    args$add       <- TRUE

    if(plot_type == "base"){
      plot <- do.call(.plot_prior_list.factor, args)
    }else if(plot_type == "ggplot"){
      plot <- plot + do.call(.plot_prior_list.factor, args)
    }

  }else{

    # plot just posterior otherwise
    plot <- .plot_prior_list.factor(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)

  }


  if(plot_type == "ggplot"){
    return(plot)
  }else{
    return(invisible())
  }

}

.plot_data_marginal_samples     <- function(samples, parameter, prior, n_points, transformation, transformation_arguments, transformation_settings){

  check_list(samples, "samples", check_names = parameter, allow_other = TRUE)

  x_points <- NULL
  y_points <- NULL
  x_den    <- NULL
  y_den    <- NULL

  # extract the relevant information
  if(is.list(samples[[parameter]]) && length(samples[[parameter]]) > 1){
    posterior_samples <- do.call(cbind, samples[[parameter]])
    if(prior){
      prior_samples <- do.call(cbind, lapply(samples[[parameter]], attr, which = "prior_samples"))
    }
  }else{
    posterior_samples  <- matrix(samples[[parameter]][[1]], ncol = 1)
    colnames(posterior_samples) <- names(samples[[parameter]])
    if(prior){
      prior_samples <- matrix(attr(samples[[parameter]][[1]], "prior_samples"), ncol = 1)
      if(is.null(prior_samples))
        stop("'samples' did not contain prior samples")
      colnames(prior_samples) <- names(samples[[parameter]])
    }
  }


  # create the output object
  out <- list()


  # deal with the densities
  for(i in 1:ncol(posterior_samples)){

    out_den <- .plot_data_marginal_samples.den(posterior_samples[,i], n_points, transformation, transformation_arguments, transformation_settings)
    attr(out_den, "level")      <- i
    attr(out_den, "level_name") <- colnames(posterior_samples)[i]

    if(prior){
      out_den.prior <- .plot_data_marginal_samples.den(prior_samples[,i], n_points, transformation, transformation_arguments, transformation_settings)
      attr(out_den.prior, "level")      <- i
      attr(out_den.prior, "level_name") <- colnames(prior_samples)[i]

      attr(out_den, "prior") <- out_den.prior
    }

    out[[paste0("density", i)]] <- out_den

  }

  return(out)
}
.plot_data_marginal_samples.den <- function(x, n_points, transformation, transformation_arguments, transformation_settings){

  args <- list(x = x, n = n_points)

  # get the density estimate
  density_continuous <- do.call(stats::density, args)
  x_den    <- density_continuous$x
  y_den    <- density_continuous$y


  # apply transformations
  if(!is.null(transformation)){
    x_den <- .density.prior_transformation_x(x_den, transformation, transformation_arguments)
    y_den <- .density.prior_transformation_y(x_den, y_den, transformation, transformation_arguments)
    x     <- .density.prior_transformation_x(x, transformation, transformation_arguments)
  }

  out_den <- list(
    call    = call("density", "mixed samples"),
    bw      = NULL,
    n       = n_points,
    x       = x_den,
    y       = y_den,
    samples = x
  )

  class(out_den) <- c("density", "density.prior", "density.prior.factor", "density.prior.simple")
  attr(out_den, "x_range")    <- range(x_den)
  attr(out_den, "y_range")    <- c(0, max(y_den))

  return(out_den)
}
