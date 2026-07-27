#' @title Add prior object to a plot
#'
#' @param x a prior
#' @param xlim plotting range of the prior
#' @param rescale_x allows to rescale x-axis in case a
#' weightfunction is plotted.
#' @param show_parameter which parameter should be returned in case of
#' multiple parameters per prior. Useful when priors for the omega
#' parameter are plotted and \code{individual = TRUE}.
#' @param scale_y2 scaling factor for a secondary axis
#' @param ... additional arguments
#' @inheritParams density.prior
#'
#' @return \code{lines.prior} returns \code{NULL}.
#'
#' @seealso [plot.prior()] [geom_prior()]
#' @rdname lines.prior
#' @export
lines.prior <- function(x, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 1000,
                        n_samples = 10000, force_samples = FALSE,
                        transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                        show_parameter = if(individual) 1 else NULL, individual = FALSE, rescale_x = FALSE, scale_y2 = 1, ...){

  # check input (most arguments are checked within density)
  .check_prior(x)
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0)

  if(is.prior.mixture(x)){
    return(lines_prior_list(x, xlim = xlim, x_seq = x_seq, x_range_quant = x_range_quant, n_points = n_points,
                            n_samples = n_samples, force_samples = force_samples,
                            transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings,
                            rescale_x = rescale_x, scale_y2 = scale_y2, ...))
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
  if(is.prior.weightfunction(x) & !individual){
    .lines.prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)
    return(invisible())
  }else if(is.prior.weightfunction(x) & individual){
    if(inherits(plot_data[[show_parameter]], "density.prior.simple")){
      .lines.prior.simple(plot_data, ...)
    }else if(inherits(plot_data[[show_parameter]], "density.prior.point")){
      .lines.prior.point(plot_data, scale_y2 = scale_y2, ...)
    }
    return(invisible())
  }

  # plot PET-PEESE
  if((is.prior.PET(x) | is.prior.PEESE(x))){
    if(!individual){
      .lines.prior.PETPEESE(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.simple")){
      .lines.prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      .lines.prior.point(plot_data, scale_y2 = scale_y2, ...)
    }
    return(invisible())
  }

  # plot spike and slab prior
  if(is.prior.spike_and_slab(x)){
    .lines.prior.spike_and_slab(plot_data, ...)
    return(invisible())
  }

  # point prior plots
  if(is.prior.point(x)){
    .lines.prior.point(plot_data, scale_y2 = scale_y2, ...)
    return(invisible())
  }

  if(is.prior.ordered(x)){
    selected <- if(is.null(show_parameter)){
      seq_along(plot_data)
    }else{
      show_parameter
    }
    for(i in selected){
      if(inherits(plot_data[[i]], "density.prior.mixed_measure")){
        .lines_prior_mixed_measure(
          plot_data[[i]],
          scale_y2 = scale_y2,
          ...
        )
      }else if(inherits(plot_data[[i]], "density.prior.point")){
        .lines.prior.point(plot_data[[i]], scale_y2 = scale_y2, ...)
      }else{
        .lines.prior.simple(plot_data[[i]], ...)
      }
    }
    return(invisible())
  }

  # plot orthonormal and meandif plots
  if(is.prior.orthonormal(x) | is.prior.meandif(x)){
    .lines.prior.orthonormal_or_meandif(plot_data, ...)
    return(invisible())
  }

  # plot discrete prior
  if(is.prior.discrete(x)){
    .lines.prior.discrete(plot_data, ...)
    return(invisible())
  }

  # default prior plots
  if(is.prior.simple(x)){
    .lines.prior.simple(plot_data, ...)
    return(invisible())
  }

  return(invisible())
}

#' @title Add prior object to a ggplot
#'
#' @inheritParams lines.prior
#' @inheritParams density.prior
#'
#' @return \code{geom_prior_list} returns an object of class 'ggplot'.
#'
#' @seealso [plot.prior()] [lines.prior()]
#' @rdname geom_prior
#' @export
geom_prior  <- function(x, xlim = NULL, x_seq = NULL, x_range_quant = NULL, n_points = 1000,
                        n_samples = 10000, force_samples = FALSE,
                        transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                        show_parameter = if(individual) 1 else NULL, individual = FALSE, rescale_x = FALSE, scale_y2 = 1, ...){

  # check input (most arguments are checked within density)
  .check_prior(x)
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_int(show_parameter, "show_parameter", allow_NULL = TRUE)
  check_real(scale_y2, "scale_y2", lower = 0)

  if(is.prior.mixture(x)){
    class(x) <- NULL
    return(geom_prior_list(x, xlim = xlim, x_seq = x_seq, x_range_quant = x_range_quant, n_points = n_points,
                            n_samples = n_samples, force_samples = force_samples,
                            transformation = transformation, transformation_arguments = transformation_arguments, transformation_settings = transformation_settings,
                            rescale_x = rescale_x, scale_y2 = scale_y2, ...))
  }

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
      geom <- .geom_prior.weightfunction(plot_data = plot_data, rescale_x = rescale_x, ...)
    }else if(inherits(plot_data[[show_parameter]], "density.prior.simple")){
      geom <- .geom_prior.simple(plot_data, ...)
    }else if(inherits(plot_data[[show_parameter]], "density.prior.point")){
      geom <- .geom_prior.point(plot_data, ...)
    }
    return(geom)
  }


  # plot PET-PEESE
  if((is.prior.PET(x) | is.prior.PEESE(x))){
    if(!individual){
      geom <- .geom_prior.PETPEESE(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.simple")){
      .geom_prior.simple(plot_data, ...)
    }else if(inherits(plot_data, "density.prior.point")){
      geom <- .geom_prior.point(plot_data, ...)
    }
    return(geom)
  }


  # plot spike and slab prior
  if(is.prior.spike_and_slab(x)){
    geom <- .geom_prior.spike_and_slab(plot_data, ...)
    return(geom)
  }


  # plot point prior
  if(is.prior.point(x)){
    geom <- .geom_prior.point(plot_data, ...)
    return(geom)
  }

  if(is.prior.ordered(x)){
    selected <- if(is.null(show_parameter)){
      seq_along(plot_data)
    }else{
      show_parameter
    }
    geom <- list()
    for(i in selected){
      component_geom <- if(inherits(
        plot_data[[i]],
        "density.prior.mixed_measure"
      )){
        .geom_prior_mixed_measure(
          plot_data[[i]],
          scale_y2 = scale_y2,
          ...
        )
      }else if(inherits(plot_data[[i]], "density.prior.point")){
        list(.geom_prior.point(
          plot_data[[i]],
          scale_y2 = scale_y2,
          ...
        ))
      }else{
        list(.geom_prior.simple(plot_data[[i]], ...))
      }
      geom <- c(geom, component_geom)
    }
    geom <- geom[!vapply(geom, is.null, logical(1))]
    return(geom)
  }


  # plot orthonormal and meandif prior
  if(is.prior.orthonormal(x) | is.prior.meandif(x)){
    geom <- .geom_prior.orthonormal_or_meandif(plot_data, ...)
    return(geom)
  }


  # plot discrete prior
  if(is.prior.discrete(x)){
    geom <- .geom_prior.discrete(plot_data, ...)
    return(geom)
  }


  # default prior plots
  if(is.prior.simple(x)){
    geom <- .geom_prior.simple(plot_data, ...)
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
.lines.prior.discrete        <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))   dots[["col"]]   else .plot.prior_settings()[["col"]]
  width     <- if(!is.null(dots[["width"]])) dots[["width"]] else .plot.prior_settings()[["width"]]
  lwd       <- if(!is.null(dots[["lwd"]]))   dots[["lwd"]]   else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))   dots[["lty"]]   else .plot.prior_settings()[["lty"]]

  graphics::rect(
    xleft   = plot_data$x - width/2,
    ybottom = 0,
    xright  = plot_data$x + width/2,
    ytop    = plot_data$y,
    lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines.prior.point           <- function(plot_data, scale_y2 = 1, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]

  if(!all(plot_data$y == 0)){
    graphics::arrows(
      x0 = plot_data$x[plot_data$y != 0],
      y0 = 0,
      y1 = plot_data$y[plot_data$y != 0] * scale_y2,
      lwd = 2*lwd, lty = lty, col = col)
  }

  return(invisible())
}
.lines_prior_mixed_measure   <- function(plot_data, scale_y2 = 1, ...){

  if(!is.null(plot_data$continuous)){
    continuous <- list(
      x = plot_data$continuous$x,
      y = plot_data$continuous$density
    )
    .lines.prior.simple(continuous, ...)
  }
  if(!is.null(plot_data$atoms) && nrow(plot_data$atoms) > 0L){
    atoms <- list(
      x = plot_data$atoms$location,
      y = plot_data$atoms$mass
    )
    .lines.prior.point(atoms, scale_y2 = scale_y2, ...)
  }

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
  x_lCI  <- plot_data$y_lCI
  x_uCI  <- plot_data$y_uCI

  if(rescale_x){
    x_at <- seq(0, 1, length.out = length(unique(plot_data$x)))
    x_at <- x_at[c(1, sort(rep(2:(length(x_at)-1), 2)), length(x_at))]
  }else{
    x_at <- x_cuts
  }


  graphics::polygon(
    x   = c(x_at,  rev(x_at)),
    y   = c(x_lCI, rev(x_uCI)),
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
    x   = c(plot_data$x,     rev(plot_data$x)),
    y   = c(plot_data$y_lCI, rev(plot_data$y_uCI)),
    col = col.fill, border = NA
  )
  graphics::lines(plot_data$x, plot_data$y, lwd = lwd, lty = lty, col = col)


  return(invisible())
}
.lines.prior.orthonormal_or_meandif <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["lwd"]]))      dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))      dots[["lty"]]      else .plot.prior_settings()[["lty"]]


  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines.prior.factor          <- function(plot_data, ...){

  dots <- list(...)
  col  <- if(!is.null(dots[["col"]][dots[["level"]]])) dots[["col"]][dots[["level"]]] else .plot.prior_settings()[["col"]]
  lty  <- if(!is.null(dots[["lty"]][dots[["level"]]])) dots[["lty"]][dots[["level"]]] else .plot.prior_settings()[["lty"]]
  lwd  <- if(!is.null(dots[["lwd"]]))                  dots[["lwd"]]                  else .plot.prior_settings()[["lwd"]]

  graphics::lines(x = plot_data$x, y = plot_data$y, type = "l", lwd = lwd, lty = lty, col = col)

  return(invisible())
}
.lines.prior.spike_and_slab  <- function(plot_data, ...){

  .lines.prior.simple(plot_data[["variable"]], ...)
  .lines.prior.point(plot_data[["inclusion"]], ...)

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
    mapping = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["y"]]),
    linewidth = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_prior.discrete         <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))   dots[["col"]]   else .plot.prior_settings()[["col"]]
  width     <- if(!is.null(dots[["width"]])) dots[["width"]] else .plot.prior_settings()[["width"]]
  lwd       <- if(!is.null(dots[["lwd"]]))   dots[["lwd"]]   else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["lty"]]))   dots[["lty"]]   else .plot.prior_settings()[["lty"]]

  geom <-   geom <- ggplot2::geom_bar(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes(
      x      = .data[["x"]],
      weight = .data[["y"]]),
    linewidth = lwd, linetype = lty, color = col, fill = col, width = width)

  return(geom)
}
.geom_prior.point            <- function(plot_data, scale_y2 = 1, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  if(!all(plot_data$y == 0)){
    geom <- ggplot2::geom_segment(
      data    = data.frame(
        x    = unique(plot_data$x[plot_data$y != 0]),
        xend = unique(plot_data$x[plot_data$y != 0]),
        y    = 0,
        yend = plot_data$y[plot_data$y != 0] * scale_y2),
      mapping = ggplot2::aes(
        x    = .data[["x"]],
        xend = .data[["xend"]],
        y    = .data[["y"]],
        yend = .data[["yend"]]),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.5, "cm")),
      linewidth = 2*lwd, linetype = lty, color = col)
  }else{
    geom <- NULL
  }

  return(geom)
}
.geom_prior_mixed_measure    <- function(plot_data, scale_y2 = 1, ...){

  geom <- list()
  if(!is.null(plot_data$continuous)){
    continuous <- list(
      x = plot_data$continuous$x,
      y = plot_data$continuous$density
    )
    geom[[length(geom) + 1L]] <- .geom_prior.simple(continuous, ...)
  }
  if(!is.null(plot_data$atoms) && nrow(plot_data$atoms) > 0L){
    atoms <- list(
      x = plot_data$atoms$location,
      y = plot_data$atoms$mass
    )
    geom[[length(geom) + 1L]] <- .geom_prior.point(
      atoms,
      scale_y2 = scale_y2,
      ...
    )
  }

  geom[!vapply(geom, is.null, logical(1))]
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
  x_lCI  <- plot_data$y_lCI
  x_uCI  <- plot_data$y_uCI

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
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]]),
      fill    = col.fill
    ),
    ggplot2::geom_line(
      data    = data.frame(
        x = x_at,
        y = x_mean),
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]]),
      linewidth = lwd, linetype = lty, color = col)
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
        y = c(plot_data$y_lCI, rev(plot_data$y_uCI))),
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]]),
      fill    = col.fill
    ),
    ggplot2::geom_line(
      data    = data.frame(
        x = plot_data$x,
        y = plot_data$y),
      mapping = ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]]),
      linewidth = lwd, linetype = lty, color = col)
  )

  return(geom)
}
.geom_prior.orthonormal_or_meandif <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]]))      dots[["col"]]      else .plot.prior_settings()[["col"]]
  lwd       <- if(!is.null(dots[["size"]]))     dots[["size"]]     else  if(!is.null(dots[["lwd"]])) dots[["lwd"]] else .plot.prior_settings()[["lwd"]]
  lty       <- if(!is.null(dots[["linetype"]])) dots[["linetype"]] else  if(!is.null(dots[["lty"]])) dots[["lty"]] else .plot.prior_settings()[["lty"]]

  geom <- ggplot2::geom_line(
    data    = data.frame(
      x = plot_data$x,
      y = plot_data$y),
    mapping = ggplot2::aes(
      x = .data[["x"]],
      y = .data[["y"]]),
    linewidth = lwd, linetype = lty, color = col)

  return(geom)
}
.geom_prior.factors          <- function(plot_data, ...){

  # this function notably differs from the .line_prior.factor counterpart
  # - it's so much more difficult to draw custom legend in ggplot2 ... :(

  dots <- list(...)
  col  <- if(!is.null(dots[["col"]]))      dots[["col"]]      else rep(.plot.prior_settings()[["col"]], length(dots[["level_names"]]))
  lty  <- if(!is.null(dots[["linetype"]])) dots[["linetype"]]
  else  if(!is.null(dots[["lty"]]))        dots[["lty"]]      else rep(.plot.prior_settings()[["lty"]], length(dots[["level_names"]]))
  lwd  <- if(!is.null(dots[["size"]]))     dots[["size"]]
  else  if(!is.null(dots[["lwd"]]))        dots[["lwd"]]      else .plot.prior_settings()[["lwd"]]
  legend_title <- if(!is.null(dots[["legend_title"]])) dots[["legend_title"]] else NULL

  names(col) <- dots[["level_names"]]
  names(lty) <- dots[["level_names"]]

  if(!is.null(dots[["hardcode"]]) && dots[["hardcode"]]){
    geom <- lapply(unique(plot_data$level), function(lvl){
      ggplot2::geom_line(
        data    = plot_data[plot_data$level == lvl,],
        mapping = ggplot2::aes(
          x        = .data[["x"]],
          y        = .data[["y"]]),
        linewidth = 1, show.legend = dots[["legend"]], color = col[lvl], linetype = lty[lvl])
    })
  }else{
    geom <- list(
      ggplot2::geom_line(
        data    = plot_data,
        mapping = ggplot2::aes(
          x        = .data[["x"]],
          y        = .data[["y"]],
          color    = .data[["level"]],
          linetype = .data[["level"]],
          group    = .data[["level"]]),
        linewidth = 1, show.legend = dots[["legend"]]),
      ggplot2::scale_linetype_manual(
        name   = legend_title,
        values = lty,
        breaks = dots[["level_names"]],
        labels = dots[["level_names"]]),
      ggplot2::scale_color_manual(
        name   = legend_title,
        values = col,
        breaks = dots[["level_names"]],
        labels = dots[["level_names"]]))
  }


  return(geom)
}
.geom_prior.spike_and_slab   <- function(plot_data, ...){

  geom <- list(
    .geom_prior.simple(plot_data[["variable"]], ...),
    .geom_prior.point(plot_data[["inclusion"]], ...)
  )

  return(geom)
}
