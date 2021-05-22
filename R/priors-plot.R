#' @title Plots a prior object
#'
#' @param x a prior
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot2"} for plotting. The later
#' requires \pkg{ggplot2} package to be installed.
#' @param par_name a type of parameter for which the prior is
#' specified. Only relevant if the prior corresponds to a mu
#' parameter that needs to be transformed.
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when priors for the omega
#' parameter are plotted and \code{weights = TRUE}.
#' @param ... additional arguments
#' @inheritParams density.prior
#'
#' @seealso [prior()]
#' @rdname plot.prior
#' @export
plot.prior <- function(x, plot_type = "base",
                       x_seq = NULL, x_range = NULL, x_range_quant = NULL, n_points = 1000,
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
  if(is.null(x_range) & is.null(x_seq)){
    x_range   <- range(x, quantiles = x_range_quant)
    x_range   <- range(pretty(x_range))
  }
  plot_data <- density(x = x, x_seq = x_seq, x_range = x_range, x_range_quant = x_range_quant,
                       n_points = n_points, n_samples = n_samples, force_samples = force_samples,
                       transformation = transformation, transformation_arguments = transformation_arguments,
                       transformation_settings = transformation_settings, individual = individual)


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

  # plot a weightfunction
  if(is.prior.weightfunction(x) & !individual){
    plot <- .plot.prior.weigthfunction(x = x, plot_type = plot_type, plot_data = plot_data, rescale_x = rescale_x, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plot)
    }else{
      return(invisible())
    }
  }

  # plot individual weights
  if(is.prior.weightfunction(x) & individual){
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

    graphics::arrows(x0 = unique(plot_data$x[is.infinite(plot_data$y)]), y0 = 0, y1 = 1,
                     lwd = 2*lwd, lty = lty, col = col)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = c(0, 1), labels = c(0, expression(infinity)), limits = ylim)

    plot <- plot + ggplot2::geom_segment(
      data    = data.frame(
        x    = unique(plot_data$x[is.infinite(plot_data$y)]),
        xend = unique(plot_data$x[is.infinite(plot_data$y)]),
        y    = 0,
        yend = 1),
      mapping = ggplot2::aes_string(
        x    = "x",
        xend = "xend",
        y    = "y",
        yend = "yend"),
      arrow   = ggplot2::arrow(length = ggplot2::unit(0.5, "cm")),
      size = lwd, linetype = lty, color = col)

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

    graphics::plot(x = plot_data$x, y = plot_data$y, type = "l", bty  = "n",
                   las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, lwd = lwd, lty = lty,
                   cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col = col, col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim)
    plot <- plot + ggplot2::geom_line(
      data    = data.frame(
        x = plot_data$x,
        y = plot_data$y),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      size = lwd, linetype = lty, color = col)

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
.plot.prior.weigthfunction <- function(x, plot_type, plot_data, rescale_x, par_name, ...){


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

  # weightfunction specific stuff

  x_cuts <- plot_data$x
  x_mean <- plot_data$y
  x_lCI  <- plot_data$y.lCI
  x_uCI  <- plot_data$y.uCI

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(x$parameters[["steps"]]) + 2)
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

    graphics::polygon(
      x = c(x_at,  rev(x_at)),
      y = c(x_lCI, rev(x_uCI)),
      col = "grey80", border = NA
    )
    graphics::lines(x_at, x_mean, lwd = lwd, lty = lty, col = col)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = x_at, labels = x_cuts, limits = xlim)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim)

    plot <- plot + ggplot2::geom_polygon(
      data    = data.frame(
        x = c(x_at,  rev(x_at)),
        y = c(x_lCI, rev(x_uCI))),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      fill    = "grey80"
    )
    plot <- plot + ggplot2::geom_line(
      data    = data.frame(
        x = x_at,
        y = x_mean),
      mapping = ggplot2::aes_string(
        x = "x",
        y = "y"),
      size = lwd, linetype = lty, color = col)

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
  lwd       = 1,
  lty       = 1
  ))
}
