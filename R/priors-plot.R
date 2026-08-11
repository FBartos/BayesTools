#' @title Plots a prior object
#'
#' @param x a prior
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting.
#' @param xlim x plotting range
#' @param par_name a type of parameter for which the prior is
#' specified. Only relevant if the prior corresponds to a mu
#' parameter that needs to be transformed.
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param show_figures which figures should be returned in case of
#' multiple plots are generated. Useful when priors for the omega
#' parameter are plotted and \code{individual = TRUE}.
#' @param ... additional graphical arguments. For mixed continuous and point
#' distributions, \code{ylim} controls the density axis, \code{ylim2} controls
#' the probability-mass axis, and \code{ylab2} controls its label.
#' @inheritParams density.prior
#'
#' @examples
#' # create some prior distributions
#' p0 <- prior(distribution = "point",  parameters = list(location = 0))
#' p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
#' p2 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1), truncation = list(0, Inf))
#'
#' # a default plot
#' plot(p0)
#'
#' # manipulate line thickness and color, change the parameter name
#' plot(p1, lwd = 2, col = "blue", par_name = bquote(mu))
#'
#' # use ggplot
#' plot(p2, plot_type = "ggplot")
#'
#' # utilize the ggplot prior geom
#' plot(p2, plot_type = "ggplot", xlim = c(-2, 2)) + geom_prior(p1, col = "red", lty = 2)
#'
#' # apply transformation
#' plot(p1, transformation = "exp")
#'
#' @return \code{plot.prior} returns either \code{NULL} or
#' an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
#' Dirichlet simplex priors are plotted as one beta marginal per coordinate;
#' the ggplot method returns a list unless a single figure is selected. For an
#' ordered level with mixed probability measure, the continuous density and
#' exact probability-mass arrows are drawn together without rescaling either
#' component.
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
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_figures, "show_figures", allow_NULL = TRUE)

  if(is.prior.mixture(x)){
    class(x) <- NULL
    return(plot_prior_list(x, plot_type = plot_type, x_seq = x_seq, xlim = xlim, x_range_quant = x_range_quant, n_points = n_points,
                           n_samples = n_samples, force_samples = force_samples, transformation = transformation,
                           transformation_arguments = transformation_arguments, transformation_settings = transformation_settings,
                           show_figures = show_figures, individual = individual, rescale_x = rescale_x, par_name = par_name, ...))
  }


  # get the plotting data
  if(is.null(xlim) & is.null(x_seq)){
    if((is.prior.PET(x) | is.prior.PEESE(x) | is.prior.weightfunction(x)) & !individual){
      xlim   <- c(0, 1)
    }else if(is.prior.spike_and_slab(x)){
      xlim   <- range(c(0,range(.get_spike_and_slab_variable(x), quantiles = x_range_quant)))
      xlim   <- range(pretty(xlim))
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

  # spike and slab prior plots
  if(is.prior.spike_and_slab(x)){
    plots <- .plot.prior.spike_and_slab(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # point prior plots
  if(is.prior.point(x)){
    plots <- .plot.prior.point( x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # ordered factor prior plots
  if(is.prior.ordered(x)){
    plots <- .plot.prior.simplex(x = x, plot_type = plot_type, plot_data = plot_data, show_figures = show_figures, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # plot orthonormal and meandif priors
  if(is.prior.orthonormal(x) | is.prior.meandif(x)){
    plots <- .plot.prior.orthonormal_or_meandif(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # simplex prior plots
  if(is.prior.simplex(x)){
    plots <- .plot.prior.simplex(x = x, plot_type = plot_type, plot_data = plot_data, show_figures = show_figures, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # discrete prior plots
  if(is.prior.discrete(x)){
    plots <- .plot.prior.discrete(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }

  # default prior plots
  if(is.prior.simple(x)){
    plots <- .plot.prior.simple(x = x, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)
    if(plot_type == "ggplot"){
      return(plots)
    }else{
      return(invisible())
    }
  }
}

.plot.prior.simplex       <- function(x, plot_type, plot_data, show_figures = NULL, par_name = NULL, ...){

  if(is.null(show_figures)){
    plots_ind <- seq_along(plot_data)
  }else{
    plots_ind <- seq_along(plot_data)[show_figures]
  }

  plots <- list()
  for(figure in plots_ind){
    component_name <- if(is.null(par_name)){
      names(plot_data)[figure]
    }else{
      paste0(par_name, "[", figure, "]")
    }
    if(inherits(plot_data[[figure]], "density.prior.mixed_measure")){
      plots[[figure]] <- .plot_prior_mixed_measure(
        x = x,
        plot_type = plot_type,
        plot_data = plot_data[[figure]],
        par_name = component_name,
        ...
      )
    }else if(inherits(plot_data[[figure]], "density.prior.point")){
      plots[[figure]] <- .plot.prior.point(
        x = x,
        plot_type = plot_type,
        plot_data = plot_data[[figure]],
        par_name = component_name,
        ...
      )
    }else{
      plots[[figure]] <- .plot.prior.simple(
        x = x,
        plot_type = plot_type,
        plot_data = plot_data[[figure]],
        par_name = component_name,
        ...
      )
    }
  }

  if(plot_type == "ggplot" && length(plots_ind) == 1L){
    plots <- plots[[plots_ind]]
  }

  return(plots)
}

.plot_prior_mixed_measure  <- function(x, plot_type, plot_data,
                                      par_name = NULL, ...){

  dots <- list(...)
  xlim <- attr(plot_data, "x_range")
  ylim <- attr(plot_data, "y_range")

  short_name <- if(is.null(dots[["short_name"]])) FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main <- ""
  xlab <- bquote(
    .(if(!is.null(par_name)) bquote(.(par_name)~"~"))~.(prior_label)
  )
  ylab <- "Density / probability mass"

  if(is.null(dots[["main"]])) dots$main <- main
  if(is.null(dots[["xlab"]])) dots$xlab <- xlab
  if(is.null(dots[["ylab"]])) dots$ylab <- ylab
  if(is.null(dots[["xlim"]])) dots$xlim <- xlim
  if(is.null(dots[["ylim"]])) dots$ylim <- ylim

  if(plot_type == "base"){
    .plot.prior_empty("simple", dots)
    .lines_prior_mixed_measure(plot_data, ...)
    return(invisible())
  }

  plot <- .ggplot.prior_empty("simple", dots)
  geoms <- .geom_prior_mixed_measure(plot_data, ...)
  if(length(geoms) > 0L){
    plot <- plot + geoms
  }
  plot
}

.plot.prior.point          <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- attr(plot_data, "x_range")
  ylim      <- c(0, 1)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main      <- if(!is.null(attr(plot_data, "steps"))) prior_label else ""
  xlab      <- if(!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"]) else bquote(.(if(is.prior.orthonormal(x) | is.prior.meandif(x))"dif")*.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label))
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Probability"

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim


  if(plot_type == "base"){

    .plot.prior_empty("point", dots)
    .lines.prior.point(plot_data, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- .ggplot.prior_empty("point", dots)
    plot <- plot + .geom_prior.point(plot_data, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.simple         <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- attr(plot_data, "x_range")
  ylim      <- attr(plot_data, "y_range")

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main      <- if(!is.null(attr(plot_data, "steps"))) prior_label else ""
  xlab      <- if(!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"])  else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label))
  ylab      <- "Density"

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim


  if(plot_type == "base"){

    .plot.prior_empty("simple", dots)
    .lines.prior.simple(plot_data, ...)
    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- .ggplot.prior_empty("simple", dots)
    plot <- plot + .geom_prior.simple(plot_data, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.discrete       <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- attr(plot_data, "x_range")
  ylim      <- attr(plot_data, "y_range")

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main      <- if(!is.null(attr(plot_data, "steps"))) prior_label else ""
  xlab      <- if(!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"])  else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label))
  ylab      <- "Probability"

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim


  if(plot_type == "base"){

    .plot.prior_empty("simple", dots)
    .lines.prior.discrete(plot_data, ...)
    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- .ggplot.prior_empty("simple", dots)
    plot <- plot + .geom_prior.discrete(plot_data, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.weightfunction <- function(x, plot_type, plot_data, rescale_x, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else bquote(italic(p)*"-value")
  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else if(is.prior(x)) bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label)) else if(!is.null(par_name)) bquote(.(par_name)) else "Selection Models"
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Probability"

  xlim      <- attr(plot_data, "x_range")
  ylim      <- if(!is.null(dots[["ylim"]])) dots[["ylim"]] else attr(plot_data, "y_range")

  # weightfunction specific stuff (required for axes)
  x_cuts <- plot_data$x

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(unique(plot_data$x)))
    x_at <- x_at[c(1, sort(rep(2:(length(x_at)-1), 2)), length(x_at))]
  }else{
    x_at <- x_cuts
  }

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim
  dots$x_at     <- unique(x_at)
  dots$x_labels <- unique(x_cuts)


  if(plot_type == "base"){

    .plot.prior_empty("weightfunction", dots)
    .lines.prior.weightfunction(plot_data, rescale_x, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- .ggplot.prior_empty("weightfunction", dots)
    plot <- plot + .geom_prior.weightfunction(plot_data, rescale_x, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.PETPEESE       <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  xlab      <- if(!is.null(dots[["xlab"]])) dots[["xlab"]] else "Standard error"
  main      <- if(!is.null(dots[["main"]])) dots[["main"]] else if(is.prior(x)) bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label)) else if(!is.null(par_name)) bquote(.(par_name)) else "PET-PEESE"
  ylab      <- if(!is.null(dots[["ylab"]])) dots[["ylab"]] else "Effect size"

  xlim      <- attr(plot_data, "x_range")
  ylim      <- if(!is.null(dots[["ylim"]])) dots[["ylim"]] else attr(plot_data, "y_range")

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim


  if(plot_type == "base"){

    .plot.prior_empty("PETPEESE", dots)
    .lines.prior.PETPEESE(plot_data, ...)

    plot <- NULL

  }else if(plot_type == "ggplot"){


    plot <- .ggplot.prior_empty("PETPEESE", dots)
    plot <- plot + .geom_prior.PETPEESE(plot_data, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.orthonormal_or_meandif <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- attr(plot_data, "x_range")
  ylim      <- attr(plot_data, "y_range")

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main      <- if(!is.null(attr(plot_data, "steps"))) prior_label else ""
  xlab      <- if(!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"])  else bquote("dif"*.(if(!is.null(par_name)){" "*bquote(.(par_name)~"~")})~.(prior_label))
  ylab      <- "Density"

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim


  if(plot_type == "base"){

    .plot.prior_empty("simple", dots)
    .lines.prior.orthonormal_or_meandif(plot_data, ...)
    plot <- NULL

  }else if(plot_type == "ggplot"){

    plot <- .ggplot.prior_empty("simple", dots)
    plot <- plot + .geom_prior.orthonormal_or_meandif(plot_data, ...)

  }

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}
.plot.prior.spike_and_slab <- function(x, plot_type, plot_data, par_name = NULL, ...){

  # get default plot settings
  dots      <- list(...)

  xlim      <- attr(plot_data, "x_range")
  ylim      <- attr(plot_data, "y_range")

  short_name      <- if(is.null(dots[["short_name"]]))      FALSE else dots[["short_name"]]
  parameter_names <- if(is.null(dots[["parameter_names"]])) FALSE else dots[["parameter_names"]]
  prior_label     <- .plot.prior_label(x, plot_data, short_name, parameter_names)

  main      <- if(!is.null(attr(plot_data, "steps"))) prior_label else ""
  xlab      <- if(!is.null(attr(plot_data, "steps"))) bquote(omega["["*.(attr(plot_data, "steps")[1])*","~.(attr(plot_data, "steps")[2])*"]"])  else bquote(.(if(!is.null(par_name)){bquote(.(par_name)~"~")})~.(prior_label))
  ylab      <- "Density"

  # add it to the user input if desired
  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab
  if(is.null(dots[["xlim"]])) dots$xlim <-  xlim
  if(is.null(dots[["ylim"]])) dots$ylim <-  ylim

  args_prior           <- dots
  args_prior$plot_data <- plot_data
  args_prior$plot_type <- plot_type

  plot <- do.call(.plot_prior_list.both, args_prior)

  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    return(plot)
  }
}

.plot.prior_empty    <- function(type, dots = list(), ...){

  dots      <- c(dots, list(...))

  main      <- if(!is.null(dots[["main"]]))     dots[["main"]]     else ""
  xlab      <- if(!is.null(dots[["xlab"]]))     dots[["xlab"]]     else ""
  ylab      <- if(!is.null(dots[["ylab"]]))     dots[["ylab"]]     else ""
  xlim      <- if(!is.null(dots[["xlim"]]))     dots[["xlim"]]     else c(0, 1)
  ylim      <- if(!is.null(dots[["ylim"]]))     dots[["ylim"]]     else c(0, 1)
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  x_at      <- if(!is.null(dots[["x_at"]]))     dots[["x_at"]]     else NULL
  x_labels  <- if(!is.null(dots[["x_labels"]])) dots[["x_labels"]] else NULL

  if(type == "point"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, axes = FALSE,
                   cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    graphics::axis(1, col = col.axis, cex = cex.axis)
    graphics::axis(2, at = ylim, labels = ylim, col = col.axis, cex = cex.axis, las = 1)

  }else if(type == "simple"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim,
                   cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col.axis = col.axis, col.lab = col.lab, col.main = col.main)

  }else if(type == "weightfunction"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim, axes = FALSE,
                   cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    graphics::axis(1, at = x_at, labels = x_labels, col = col.axis, cex = cex.axis)
    graphics::axis(2, at = pretty(ylim), col = col.axis, cex = cex.axis, las = 1)

  }else if(type == "PETPEESE"){

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = ylim,
                   cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col.axis = col.axis, col.lab = col.lab, col.main = col.main)

  }else if(type == "both"){

    ylim2    <- if(!is.null(dots[["ylim2"]]))     dots[["ylim2"]]     else ylim
    ylab2    <- if(!is.null(dots[["ylab2"]]))     dots[["ylab2"]]     else ""
    scale_y2 <- if(!is.null(dots[[".scale_y2_resolved"]])){
      dots[[".scale_y2_resolved"]]
    }else{
      .plot_scale_y2_from_limits(ylim, ylim2, dots[["scale_y2"]])
    }

    graphics::plot(NA, type = "n", bty  = "n", las = 1, xlab = xlab, ylab = ylab, main = main,
                   xlim = xlim, ylim = range(c(pretty(ylim), pretty(ylim2) * scale_y2)), axes = FALSE,
                   cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
                   col.axis = col.axis, col.lab = col.lab, col.main = col.main)
    graphics::axis(1, at = x_at, labels = x_labels, col = col.axis, cex = cex.axis)
    graphics::axis(2, at = pretty(ylim),             labels = pretty(ylim),  col = col.axis, cex = cex.axis, las = 1)
    graphics::axis(4, at = pretty(ylim2) * scale_y2, labels = pretty(ylim2), col = col.axis, cex = cex.axis, las = 1)
    graphics::mtext(ylab2, side = 4, line = 3)

    return(invisible(list(scale_y2 = scale_y2)))
  }

  return(invisible())
}
.ggplot.prior_empty  <- function(type, dots = list(), ...){

  dots      <- c(dots, list(...))

  main      <- if(!is.null(dots[["main"]]))     dots[["main"]]     else ""
  xlab      <- if(!is.null(dots[["xlab"]]))     dots[["xlab"]]     else ""
  ylab      <- if(!is.null(dots[["ylab"]]))     dots[["ylab"]]     else ""
  xlim      <- if(!is.null(dots[["xlim"]]))     dots[["xlim"]]     else c(0, 1)
  ylim      <- if(!is.null(dots[["ylim"]]))     dots[["ylim"]]     else c(0, 1)
  col.main  <- if(!is.null(dots[["col.main"]])) dots[["col.main"]] else .plot.prior_settings()[["col.main"]]
  cex.axis  <- if(!is.null(dots[["cex.axis"]])) dots[["cex.axis"]] else .plot.prior_settings()[["cex.axis"]]
  cex.lab   <- if(!is.null(dots[["cex.lab"]]))  dots[["cex.lab"]]  else .plot.prior_settings()[["cex.lab"]]
  cex.main  <- if(!is.null(dots[["cex.main"]])) dots[["cex.main"]] else .plot.prior_settings()[["cex.main"]]
  col.axis  <- if(!is.null(dots[["col.axis"]])) dots[["col.axis"]] else .plot.prior_settings()[["col.axis"]]
  col.lab   <- if(!is.null(dots[["col.lab"]]))  dots[["col.lab"]]  else .plot.prior_settings()[["col.lab"]]
  x_at      <- if(!is.null(dots[["x_at"]]))     dots[["x_at"]]     else NULL
  x_labels  <- if(!is.null(dots[["x_labels"]])) dots[["x_labels"]] else NULL

  if(type == "point"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = range(pretty(xlim)), oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = range(pretty(ylim)), oob = scales::oob_keep)

    attr(plot, "sec_axis") <- FALSE

  }else if(type == "simple"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = range(pretty(xlim)), oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = range(pretty(ylim)), oob = scales::oob_keep)

    attr(plot, "sec_axis") <- FALSE

  }else if(type == "weightfunction"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = x_at, labels = x_labels, limits = xlim, oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = ylim,            oob = scales::oob_keep)

    attr(plot, "sec_axis") <- FALSE

  }else if(type == "PETPEESE"){

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = range(pretty(xlim)), oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = range(pretty(ylim)), oob = scales::oob_keep)

    attr(plot, "sec_axis") <- FALSE

  }else if(type == "both"){

    ylim2    <- if(!is.null(dots[["ylim2"]]))     dots[["ylim2"]]     else ylim
    ylab2    <- if(!is.null(dots[["ylab2"]]))     dots[["ylab2"]]     else ""
    scale_y2 <- if(!is.null(dots[[".scale_y2_resolved"]])){
      dots[[".scale_y2_resolved"]]
    }else{
      .plot_scale_y2_from_limits(ylim, ylim2, dots[["scale_y2"]])
    }

    plot <- ggplot2::ggplot()
    plot <- plot + ggplot2::ggtitle(main)
    plot <- plot + ggplot2::scale_x_continuous(name = xlab, breaks = pretty(xlim), limits = range(pretty(xlim)), oob = scales::oob_keep)
    plot <- plot + ggplot2::scale_y_continuous(name = ylab, breaks = pretty(ylim), limits = range(c(pretty(ylim), pretty(ylim2) * scale_y2)), oob = scales::oob_keep, sec.axis = ggplot2::sec_axis(~ ., name = ylab2, breaks = pretty(ylim2) * scale_y2, labels = pretty(ylim2)))

    attr(plot, "scale_y2") <- scale_y2
    attr(plot, "sec_axis") <- TRUE

  }

  plot <- plot + ggplot2::theme(
    axis.text  = ggplot2::element_text(size = 10 * cex.axis, color = col.axis),
    axis.title = ggplot2::element_text(size = 10 * cex.lab,  color = col.lab),
    title      = ggplot2::element_text(size = 10 * cex.main, color = col.main))

  return(plot)
}

.plot.prior_label <- function(x, plot_data, short_name = FALSE, parameter_names = FALSE){

  if(!is.prior(x)){
    return(NULL)
  }

  label <- print(x, plot = TRUE, short_name = short_name, parameter_names = parameter_names)
  transformation <- attr(plot_data, "transformation")

  if(is.null(transformation)){
    return(label)
  }

  paste0(.plot.prior_transformation_label(transformation), "(", print(x, silent = TRUE), ")")
}
.plot.prior_transformation_label <- function(transformation){

  if(is.character(transformation) && length(transformation) == 1){
    return(switch(
      transformation,
      "lin"     = "linear",
      "exp_lin" = "exp-linear",
      "tanh"    = "tanh",
      "exp"     = "exp",
      "transformed"
    ))
  }

  "transformed"
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
  lty       = 1,
  scale_y2  = 1.10,
  width     = 0.20
  ))
}
.plot_scale_y2_from_limits <- function(ylim, ylim2, scale_y2 = NULL){

  if(is.null(scale_y2)){
    scale_y2 <- .plot.prior_settings()[["scale_y2"]]
  }

  return(scale_y2 * max(pretty(ylim)) / max(pretty(ylim2)))
}
.plot_scale_y2_resolve <- function(plot_data, dots = list()){

  is_simple <- sapply(plot_data, inherits, what = "density.prior.simple")
  is_point  <- sapply(plot_data, inherits, what = "density.prior.point")

  if(any(is_simple) && (any(is_point) || !is.null(dots[["ylim2"]]))){
    ylim <- if(!is.null(dots[["ylim"]])){
      dots[["ylim"]]
    }else{
      range(as.vector(sapply(plot_data[is_simple], attr, which = "y_range")))
    }
    ylim2 <- if(!is.null(dots[["ylim2"]])){
      dots[["ylim2"]]
    }else{
      range(as.vector(sapply(plot_data[is_point], attr, which = "y_range")))
    }
    return(.plot_scale_y2_from_limits(ylim, ylim2, dots[["scale_y2"]]))
  }

  return(1)
}
.get_scale_y2        <- function(plot_data, ...){

  return(.plot_scale_y2_resolve(plot_data, list(...)))
}

.plot_scale_y2_remember <- function(scale_y2, ylim2 = NULL){

  device <- as.integer(grDevices::dev.cur())
  if(device == 1L){
    return(invisible(NULL))
  }

  states <- .BayesTools_private$plot_scale_y2_states
  if(is.null(states)){
    states <- list()
  }
  key <- as.character(device)

  if(is.null(scale_y2)){
    states[[key]] <- NULL
  }else{
    states[[key]] <- list(
      scale_y2 = scale_y2,
      ylim2    = ylim2,
      usr      = unname(graphics::par("usr"))
    )
  }
  .BayesTools_private$plot_scale_y2_states <- states

  return(invisible(NULL))
}

.plot_scale_y2_state_current <- function(){

  device <- as.integer(grDevices::dev.cur())
  if(device == 1L){
    return(NULL)
  }

  states <- .BayesTools_private$plot_scale_y2_states
  state  <- states[[as.character(device)]]
  if(is.null(state)){
    return(NULL)
  }

  usr <- tryCatch(
    unname(graphics::par("usr")),
    error = function(error) NULL
  )
  if(!isTRUE(all.equal(usr, state[["usr"]], tolerance = 1e-12))){
    .plot_scale_y2_remember(NULL)
    return(NULL)
  }

  return(state)
}
.plot_point_mass_warn_outside <- function(plot_data, ylim2){

  if(is.null(ylim2)){
    return(invisible(NULL))
  }

  if(inherits(plot_data, "density.prior.point")){
    plot_data <- list(plot_data)
  }

  is_point <- sapply(plot_data, inherits, what = "density.prior.point")
  if(!any(is_point)){
    return(invisible(NULL))
  }

  probabilities <- unlist(lapply(plot_data[is_point], function(x) x[["y"]]))
  if(any(probabilities < min(ylim2) | probabilities > max(ylim2))){
    warning(
      "Point-mass probabilities outside the active secondary-axis limits ",
      "will be clipped. Redraw the initial plot with a wider 'ylim2'.",
      call. = FALSE
    )
  }

  return(invisible(NULL))
}
.transfer_dots       <- function(dots, ...){

  dots_main <- list(...)

  dots$main      <- if(!is.null(dots_main[["main"]]))     dots_main[["main"]]
  dots$xlab      <- if(!is.null(dots_main[["xlab"]]))     dots_main[["xlab"]]
  dots$ylab      <- if(!is.null(dots_main[["ylab"]]))     dots_main[["ylab"]]
  dots$ylab2     <- if(!is.null(dots_main[["ylab2"]]))    dots_main[["ylab2"]]
  dots$xlim      <- if(!is.null(dots_main[["xlim"]]))     dots_main[["xlim"]]
  dots$ylim      <- if(!is.null(dots_main[["ylim"]]))     dots_main[["ylim"]]
  dots$ylim2     <- if(!is.null(dots_main[["ylim2"]]))    dots_main[["ylim2"]]
  dots$col.main  <- if(!is.null(dots_main[["col.main"]])) dots_main[["col.main"]]
  dots$cex.axis  <- if(!is.null(dots_main[["cex.axis"]])) dots_main[["cex.axis"]]
  dots$cex.lab   <- if(!is.null(dots_main[["cex.lab"]]))  dots_main[["cex.lab"]]
  dots$cex.main  <- if(!is.null(dots_main[["cex.main"]])) dots_main[["cex.main"]]
  dots$col.axis  <- if(!is.null(dots_main[["col.axis"]])) dots_main[["col.axis"]]
  dots$col.lab   <- if(!is.null(dots_main[["col.lab"]]))  dots_main[["col.lab"]]
  dots$x_at      <- if(!is.null(dots_main[["x_at"]]))     dots_main[["x_at"]]
  dots$x_labels  <- if(!is.null(dots_main[["x_labels"]])) dots_main[["x_labels"]]

  return(dots)
}
