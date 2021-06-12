#' @title Plots a prior object
#'
#' @param x a prior
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param xlim x plotting range
#' @param par_name a type of parameter for which the prior is
#' specified. Only relevant if the prior corresponds to a mu
#' parameter that needs to be transformed.
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when priors for the omega
#' parameter are plotted and \code{individual = TRUE}.
#' @param ... additional arguments
#' @inheritParams density.prior
#'
#' @seealso [prior()] [lines.prior()]  [geom_prior()]
#' @rdname plot.prior
#' @export
plot.prior <- function(x, plot_type = "base",
                       x_seq = NULL, xlim = NULL, x_range_quant = NULL, n_points = 1000,
                       n_samples = 10000, force_samples = FALSE,
                       transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                       show_figures = if(individual) -1 else NULL, individual = FALSE, rescale_x = FALSE, par_name = NULL, ...){

  # check input (most arguments are checked within density)
  .check_prior(x)
  check_char(plot_type, "plot_type")
  if(!plot_type %in% c("base", "ggplot"))
    stop("The passed 'plot_type' is not supported for plotting. See '?plot.RoBMA' for more details.")
  if(plot_type == "ggplot" && !try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)


  # get the plotting data
  if(is.null(xlim) & is.null(x_seq)){
    if((is.prior.PET(x) | is.prior.PEESE(x) | is.prior.weightfunction(x)) & !individual){
      xlim   <- c(0, 1)
    }else{
      xlim   <- range(x, quantiles = x_range_quant)
      xlim   <- range(pretty(xlim))
    }
  }
  plot_data <- density(x = x, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                       n_points = n_points, n_samples = n_samples, force_samples = force_samples,
                       transformation = transformation, transformation_arguments = transformation_arguments,
                       transformation_settings = transformation_settings, individual = individual)


  # plot a weightfunction
  if(is.prior.weightfunction(x)){

    if(!individual){
      plot <- .plot.prior.weightfunction(x = x, plot_type = plot_type, plot_data = plot_data, rescale_x = rescale_x, par_name = par_name, ...)
      if(plot_type == "ggplot"){
        return(plot)
      }else{
        return(invisible())
      }
    }

    # plot individual weights
    if(individual){
      # deal with the multiple figures
      if(is.null(show_figures)){
        plots_ind <- c(1:length(plot_data))
      }else{
        plots_ind <- c(1:length(plot_data))[show_figures]
      }

      # reverse the order
      plot_data <- plot_data[length(plot_data):1]

      plots <- list()
      for(figure in plots_ind){
        if(inherits(plot_data[[figure]], "density.prior.simple")){
          plots[[figure]] <- .plot.prior.simple(x = x, plot_type = plot_type, plot_data = plot_data[[figure]], par_name = par_name, ...)
        }else if(inherits(plot_data[[figure]], "density.prior.point")){
          plots[[figure]] <- .plot.prior.point( x = x, plot_type = plot_type, plot_data = plot_data[[figure]], par_name = par_name, ...)
        }
      }
      if(plot_type == "ggplot"){
        if(length(plots_ind) == 1){
          plots <- plots[[plots_ind]]
        }
      }

      if(plot_type == "ggplot"){
        return(plots)
      }else{
        return(invisible())
      }
    }

  }


  # plot PET-PEESE
  if((is.prior.PET(x) | is.prior.PEESE(x)) & !individual){
    plots <- .plot.prior.PETPEESE(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }


  # default prior plots
  if(is.prior.simple(x)){
    if(inherits(plot_data, "density.prior.simple")){
      plots <- .plot.prior.simple(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      plots <- .plot.prior.point( x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    }

    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }
}

.plot.prior.point          <- function(x, plot_type, plot_data, par_name, ...){
  # get some common settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]

  cex       <- if(!is.null(dots[["cex"]]))      dots[["cex"]]      else .plot.prior_settings()[["cex"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]
  xlim      <- attr(plot_data, "x_range")
  ylim      <- c(0, 1)

  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else if (!is.null(attr(plot_data, "steps"))) print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names) else ""
  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else if (!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"]) else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names)))
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Density"

  if(plot_type == "base"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, axes = FALSE,
                   cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col = col, col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    graphics::axis(1, col = col.axis, cex = cex.axis)
    graphics::axis(2, at = c(0, 1), labels = c(0, expression(infinity)), col = col.axis, cex = cex.axis, las = 1)
    .lines.prior.point(plot_data, 1, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = c(0, 1), labels = c(0, expression(infinity)), limits = ylim)
    plot <- plot + .geom_prior.point(plot_data, 1, ...)

    plot <- plot + ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10 * cex.axis, color = col.axis),
      axis.title = ggplot2::element_text(size = 10 * cex.lab,  color = col.lab),
      title      = ggplot2::element_text(size = 10 * cex.main, color = col.main))
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.simple         <- function(x, plot_type, plot_data, par_name, ...){

  # get some common settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]

  cex       <- if(!is.null(dots[["cex"]]))      dots[["cex"]]      else .plot.prior_settings()[["cex"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]
  xlim      <- attr(plot_data, "x_range")
  ylim      <- if(!is.null(dots[["ylim"]])) dots[["ylim"]] else attr(plot_data, "y_range")

  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else if (!is.null(attr(plot_data, "steps"))) print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names) else ""
  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else if (!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"])  else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names)))
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Density"

  if(plot_type == "base"){

    graphics::plot(NA, type = "n", bty  = "n",
                   las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, lwd = lwd, lty = lty,
                   cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col = col, col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    .lines.prior.simple(plot_data, ...)
    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim)
    plot <- plot + .geom_prior.simple(plot_data, ...)

    plot <- plot + ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10 * cex.axis, color = col.axis),
      axis.title = ggplot2::element_text(size = 10 * cex.lab,  color = col.lab),
      title      = ggplot2::element_text(size = 10 * cex.main, color = col.main))
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.weightfunction <- function(x, plot_type, plot_data, rescale_x, par_name, ...){

  # get some common settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]

  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else bquote(italic(p)*"-value")
  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names)))
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Probability"
  cex       <- if(!is.null(dots[["cex"]]))      dots[["cex"]]      else .plot.prior_settings()[["cex"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]
  xlim      <- attr(plot_data, "x_range")
  ylim      <- if(!is.null(dots[["ylim"]])) dots[["ylim"]] else c(0, 1)


  # weightfunction specific stuff (required for axes)
  x_cuts <- plot_data$x

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(unique(plot_data$x)))
    x_at <- x_at[c(1, sort(rep(2:(length(x_at)-1), 2)), length(x_at))]
  }else{
    x_at <- x_cuts
  }


  if(plot_type == "base"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, axes = FALSE,
                   cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col = col, col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    graphics::axis(1, unique(x_at), unique(x_cuts), col = col.axis, cex = cex.axis)
    graphics::axis(2, pretty(ylim), col = col.axis, cex = cex.axis, las = 1)
    .lines.prior.weightfunction(plot_data, rescale_x, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = x_at, labels = x_cuts, limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim)
    plot <- plot + .geom_prior.weightfunction(plot_data, rescale_x, ...)

    plot <- plot + ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10 * cex.axis, color = col.axis),
      axis.title = ggplot2::element_text(size = 10 * cex.lab,  color = col.lab),
      title      = ggplot2::element_text(size = 10 * cex.main, color = col.main))
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.PETPEESE       <- function(x, plot_type, plot_data, par_name, ...){

  # get some common settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]

  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else "Standard error"
  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names)))
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Effect size"
  cex       <- if(!is.null(dots[["cex"]]))      dots[["cex"]]      else .plot.prior_settings()[["cex"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]
  xlim      <- attr(plot_data, "x_range")
  ylim      <- if(!is.null(dots[["ylim"]])) dots[["ylim"]] else c(0, max(plot_data$y))


  if(plot_type == "base"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim,
                   cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col = col, col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    .lines.prior.PETPEESE(plot_data, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = xlim, oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim, oob = scales::oob_keep)
    plot <- plot + .geom_prior.PETPEESE(plot_data, ...)

    plot <- plot + ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10 * cex.axis, color = col.axis),
      axis.title = ggplot2::element_text(size = 10 * cex.lab,  color = col.lab),
      title      = ggplot2::element_text(size = 10 * cex.main, color = col.main))
  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

.plot.prior_settings <- function(){
  return(list(
  cex       = 1,
  cex.axis  = 1,
  cex.lab   = 1,
  cex.main  = 1,
  col       = "black",
  col.axis  = "black",
  col.lab   = "black",
  col.main  = "black",
  col.fill  = "grey80",
  lwd       = 1,
  lty       = 1
  ))
}





#' @title Add prior object to a plot
#'
#' @param x a prior
#' @param xlim plotting range of the prior
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param point_top y range of the spike prior distribution
#' @param show_parameter which parameter should be returned in case of
#' multiple parameters per prior. Useful when priors for the omega
#' parameter are plotted and \code{individual = TRUE}.
#' @param ... additional arguments
#' @inheritParams density.prior
#'
#' @seealso [plot.prior()] [geom_prior()]
#' @rdname lines.prior
#' @export
lines.prior <- function(x, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 1000,
                        n_samples = 10000, force_samples = FALSE,
                        transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                        show_parameter = if(individual) 1 else NULL, individual = FALSE, rescale_x = FALSE, point_top = 1, ...){

  # check input (most arguments are checked within density)
  .check_prior(x)
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(point_top, "point_top", lower = 0)


  # get the plotting data
  if(is.null(xlim) & is.null(x_seq)){
    xlim   <- range(x, quantiles = x_range_quant)
    xlim   <- range(pretty(xlim))
  }
  plot_data <- density(x = x, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                       n_points = n_points, n_samples = n_samples, force_samples = force_samples,
                       transformation = transformation, transformation_arguments = transformation_arguments,
                       transformation_settings = transformation_settings, individual = individual)


  # plot a weightfunction
  if(is.prior.weightfunction(x) & !individual){
    .lines.prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)
    return(invisible())
  }else if(is.prior.weightfunction(x) & individual){
    if(inherits(plot_data[[show_parameter]], "density.prior.simple")){
      .lines.prior.simple(plot_data, ...)
    }else if(inherits(plot_data[[show_parameter]], "density.prior.point")){
      .lines.prior.point(plot_data, point_top = point_top, ...)
    }
    return(invisible())
  }


  # plot PET-PEESE
  if((is.prior.PET(x) | is.prior.PEESE(x)) | !individual){
    if(inherits(plot_data, "density.prior.PET") | inherits(plot_data, "density.prior.PEESE")){
      .lines.prior.PETPEESE(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.simple")){
      .lines.prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      .lines.prior.point(plot_data, point_top = point_top, ...)
    }
    return(invisible())
  }


  # default prior plots
  if(is.prior.simple(x)){
    if(inherits(plot_data, "density.prior.simple")){
      .lines.prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      .lines.prior.point(plot_data, point_top = point_top, ...)
    }
    return(invisible())
  }

  return(invisible())
}

#' @title Add prior object to a ggplot
#'
#' @inheritParams lines.prior
#' @inheritParams density.prior
#'
#' @seealso [plot.prior()] [lines.prior()]
#' @rdname geom_prior
#' @export
geom_prior  <- function(x, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 1000,
                        n_samples = 10000, force_samples = FALSE,
                        transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                        show_parameter = if(individual) 1 else NULL, individual = FALSE, rescale_x = FALSE, point_top = 1, ...){

  # check input (most arguments are checked within density)
  .check_prior(x)
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(point_top, "point_top", lower = 0)
  if(!try(requireNamespace("ggplot2")))
    stop("ggplot2 package needs to be installed. Run 'install.packages('ggplot2')'")


  # get the plotting data
  if(is.null(xlim) & is.null(x_seq)){
    xlim   <- range(x, quantiles = x_range_quant)
    xlim   <- range(pretty(xlim))
  }
  plot_data <- density(x = x, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                       n_points = n_points, n_samples = n_samples, force_samples = force_samples,
                       transformation = transformation, transformation_arguments = transformation_arguments,
                       transformation_settings = transformation_settings, individual = individual)


  # plot a weightfunction
  if(is.prior.weightfunction(x) & !individual){
    geom <- .geom_prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)
  }else if(is.prior.weightfunction(x[[show_parameter]]) & individual){
    if(inherits(plot_data[[show_parameter]], "density.prior.simple")){
      geom <- .geom_prior.simple(plot_data, ...)
    }else if(inherits(plot_data[[show_parameter]], "density.prior.point")){
      geom <- .geom_prior.point(plot_data, point_top = point_top, ...)
    }
    return(geom)
  }


  # plot PET-PEESE
  if((is.prior.PET(x) | is.prior.PEESE(x)) | !individual){
    if(inherits(plot_data, "density.prior.PET") | inherits(plot_data, "density.prior.PEESE")){
      geom <- .geom_prior.PETPEESE(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.simple")){
      .geom_prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      geom <- .geom_prior.point(plot_data, point_top = point_top, ...)
    }
    return(geom)
  }


  # default prior plots
  if(is.prior.simple(x)){
    if(inherits(plot_data, "density.prior.simple")){
      geom <- .geom_prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      geom <- .geom_prior.point(plot_data, point_top = point_top, ...)
    }
    return(geom)
  }

  return(invisible())
}


# base plot prior plot elements
.lines.prior.simple          <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines.prior.point           <- function(plot_data, point_top, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::arrows(x0 = unique(plot_data$x[is.infinite(plot_data$y)]), y0 = 0, y1 = point_top, lwd = 2*lwd, lty = lty, col = col)

  return(invisible())
}
.lines.prior.weightfunction  <- function(plot_data, rescale_x, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.fill  <- if(!is.null(dots[["col.fill"]])) dots[["col.fill"]] else .plot.prior_settings()[["col.fill"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]

  # weightfunction specific stuff
  x_cuts <- plot_data$x
  x_mean <- plot_data$y
  x_lCI  <- plot_data$y.lCI
  x_uCI  <- plot_data$y.uCI

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(unique(plot_data$x)))
    x_at <- x_at[c(1, sort(rep(2:(length(x_at)-1), 2)), length(x_at))]
  }else{
    x_at <- x_cuts
  }


  graphics::polygon(
    x = c(x_at,  rev(x_at)),
    y = c(x_lCI, rev(x_uCI)),
    col = col.fill, border = NA
  )
  graphics::lines(x_at, x_mean, lwd = lwd, lty = lty, col = col)


  return(invisible())
}
.lines.prior.PETPEESE        <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.fill  <- if(!is.null(dots[["col.fill"]])) dots[["col.fill"]] else .plot.prior_settings()[["col.fill"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::polygon(
    x = c(plot_data$x,  rev(plot_data$x)),
    y = c(plot_data$y.lCI, rev(plot_data$y.uCI)),
    col = col.fill, border = NA
  )
  graphics::lines(plot_data$x, plot_data$y, lwd = lwd, lty = lty, col = col)


  return(invisible())
}

# ggplot prior plot elements
.geom_prior.simple           <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  geom <- ggplot2::geom_line(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes_string(
      x = "x",
      y = "y"),
    size = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_prior.point            <- function(plot_data, point_top, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]


  geom <- ggplot2::geom_segment(
    data    = data.frame(
      x    = unique(plot_data$x[is.infinite(plot_data$y)]),
      xend = unique(plot_data$x[is.infinite(plot_data$y)]),
      y    = 0,
      yend = point_top),
    mapping = ggplot2::aes_string(
      x    = "x",
      xend = "xend",
      y    = "y",
      yend = "yend"),
    arrow   = ggplot2::arrow(length = ggplot2::unit(0.5, "cm")),
    size = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_prior.weightfunction   <- function(plot_data, rescale_x, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.fill  <- if(!is.null(dots[["col.fill"]])) dots[["col.fill"]] else .plot.prior_settings()[["col.fill"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  # weightfunction specific stuff
  x_cuts <- plot_data$x
  x_mean <- plot_data$y
  x_lCI  <- plot_data$y.lCI
  x_uCI  <- plot_data$y.uCI

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(unique(plot_data$x)))
    x_at <- x_at[c(1, sort(rep(2:(length(x_at)-1), 2)), length(x_at))]
  }else{
    x_at <- x_cuts
  }


  geom <- list(
    ggplot2::geom_polygon(
      data    = data.frame(
        x = c(x_at,  rev(x_at)),
        y = c(x_lCI, rev(x_uCI))),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      fill    = col.fill
    ),
    ggplot2::geom_line(
      data    = data.frame(
        x = x_at,
        y = x_mean),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      size = lwd, linetype = lty, color = col)
  )

  return(geom)
}
.geom_prior.PETPEESE         <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  col.fill  <- if(!is.null(dots[["col.fill"]])) dots[["col.fill"]] else .plot.prior_settings()[["col.fill"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]


  geom <-  list(
    ggplot2::geom_polygon(
      data    = data.frame(
        x = c(plot_data$x,  rev(plot_data$x)),
        y = c(plot_data$y.lCI, rev(plot_data$y.uCI))),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      fill    = col.fill
    ),
    ggplot2::geom_line(
      data    = data.frame(
        x = plot_data$x,
        y = plot_data$y),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      size = lwd, linetype = lty, color = col)
  )

  return(geom)
}
