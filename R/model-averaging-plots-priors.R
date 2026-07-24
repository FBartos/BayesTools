#' @title Plot a list of prior distributions
#'
#' @param prior_list list of prior distributions
#' @param prior_list_mu list of priors for the mu parameter
#' required when plotting PET-PEESE
#' @param effect_direction direction of the effect for PET-PEESE
#' regression. Use \code{"positive"} (default) for
#' \code{mu + PET*se + PEESE*se^2} or \code{"negative"} for
#' \code{mu - PET*se - PEESE*se^2}.
#' @param legend whether factor legends should be drawn.
#' @param legend_title optional title for factor legends.
#' @param legend_labels optional labels for factor legend levels.
#' @param legend_position optional legend position for factor legends.
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
                            individual = FALSE, show_figures = if(individual) 1 else NULL,
                            transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                            rescale_x = FALSE, par_name = NULL, prior_list_mu = NULL, effect_direction = "positive",
                            legend = TRUE, legend_title = NULL, legend_labels = NULL, legend_position = NULL, ...){

  # check input (most arguments are checked within density)
  check_list(prior_list, "prior_list")
  if(is.prior(prior_list) | !all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  check_bool(individual, "individual")
  check_bool(rescale_x, "rescale_x")
  check_bool(legend, "legend", allow_NA = FALSE)
  check_int(show_figures, "show_figures", allow_NULL = TRUE)
  check_char(effect_direction, "effect_direction", allow_values = c("positive", "negative"))
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

  if(prior_type == "PETPEESE" && !individual){
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
  if(prior_type == "weightfunction" && !individual){
    # special dispatching for visualizing the whole weightfunction

    # use analytical marginal summaries of the mapped cumulative Dirichlet weights
    plot_data <- .plot_data_prior_list.weightfunction(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                      n_points = n_points, n_samples = n_samples)
    plot <- .plot.prior.weightfunction(prior_list, plot_type = plot_type, plot_data = plot_data, rescale_x = rescale_x, par_name = par_name, ...)

  }else if(prior_type == "PETPEESE" && !individual){
    # special dispatching for visualizing the PET-PEESE regression

    # use deterministic linear-combination summaries when supported, with a sampling fallback
    plot_data <- .plot_data_prior_list.PETPEESE(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                                n_points = n_points, n_samples = n_samples,
                                                transformation = transformation, transformation_arguments = transformation_arguments,
                                                transformation_settings = transformation_settings, prior_list_mu = prior_list_mu,
                                                effect_direction = effect_direction)
    plot <- .plot.prior.PETPEESE(prior_list, plot_type = plot_type, plot_data = plot_data, par_name = par_name, ...)

  }else if(prior_type %in% c("simple", "orthonormal", "meandif")){
    # regular prior distributions (or individual plots for parameters from weightfunctions/PET-PEESE)

    # solve analytically
    plot_data <- .plot_data_prior_list.simple(prior_list, x_seq = x_seq, x_range = xlim, x_range_quant = x_range_quant,
                                              n_points = n_points, n_samples = n_samples, force_samples = force_samples, individual = individual,
                                              transformation = transformation, transformation_arguments = transformation_arguments,
                                              transformation_settings = transformation_settings)
    if(any(sapply(plot_data, inherits, what = "density.prior.factor"))){
      plot <- .plot_prior_list.factor(
        plot_data = plot_data, plot_type = plot_type, par_name = par_name,
        legend = legend, legend_title = legend_title, legend_labels = legend_labels,
        legend_position = legend_position, ...
      )
    }else{
      plot <- .plot_prior_list.both(plot_data = plot_data, plot_type = plot_type, par_name = par_name, ...)
    }

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
.plot_prior_factor_component_level_name <- function(component){

  level_name <- attr(component, "level_name")
  if(!is.null(level_name) && length(level_name) == 1L && !is.na(level_name)){
    return(as.character(level_name))
  }

  if(inherits(component, "density.prior.factor")){
    level <- attr(component, "level")
    if(!is.null(level) && length(level) == 1L && !is.na(level)){
      return(as.character(level))
    }
  }

  NA_character_
}
.plot_prior_factor_format_level_names <- function(level_names){

  if(length(level_names) == 0L){
    return(character())
  }

  level_labels <- level_names
  dif_matches  <- gregexpr("[dif:", level_labels, fixed = TRUE)
  has_dif       <- vapply(dif_matches, function(x) x[1] != -1L, logical(1))
  multi_dif     <- vapply(dif_matches, function(x) x[1] != -1L && length(x) > 1, logical(1))
  if(any(multi_dif)){
    level_labels[multi_dif] <- gsub("__xXx__", ":", level_labels[multi_dif], fixed = TRUE)
  }
  if(any(has_dif & !multi_dif)){
    single_dif <- has_dif & !multi_dif
    level_labels[single_dif] <- substr(
      level_labels[single_dif],
      regexpr("[dif:", level_labels[single_dif], fixed = TRUE) + 5,
      regexpr("]", level_labels[single_dif], fixed = TRUE) - 1
    )
  }
  no_dif <- !has_dif
  if(any(no_dif & grepl("[", level_labels, fixed = TRUE))){
    bracket_names <- no_dif & grepl("[", level_labels, fixed = TRUE)
    level_labels[bracket_names] <- substr(
      level_labels[bracket_names],
      regexpr("[", level_labels[bracket_names], fixed = TRUE) + 1,
      regexpr("]", level_labels[bracket_names], fixed = TRUE) - 1
    )
  }

  level_labels
}
.plot_prior_factor_normalize_data <- function(plot_data){

  is_point  <- vapply(plot_data, inherits, logical(1), what = "density.prior.point")
  is_factor <- vapply(plot_data, inherits, logical(1), what = "density.prior.factor")

  component_levels <- vapply(plot_data, .plot_prior_factor_component_level_name, character(1))
  level_names_raw  <- unique(component_levels[is_factor & !is.na(component_levels)])
  level_names      <- .plot_prior_factor_format_level_names(level_names_raw)

  plot_data <- lapply(seq_along(plot_data), function(i){
    component <- plot_data[[i]]
    level_id  <- match(component_levels[i], level_names_raw)

    if(is.na(level_id)){
      level_id <- NA_integer_
    }

    attr(component, "plot_component") <- if(is_point[i]) "point" else if(is_factor[i]) "density" else "other"
    attr(component, "component_id")    <- i
    attr(component, "level_id")       <- level_id
    attr(component, "level_label")    <- if(!is.na(level_id)) level_names[level_id] else NA_character_

    component
  })
  names(plot_data) <- names(is_point)

  list(
    plot_data       = plot_data,
    points          = plot_data[is_point],
    densities       = plot_data[is_factor & !is_point],
    level_names_raw = level_names_raw,
    level_names     = level_names
  )
}
.plot_prior_factor_style_value <- function(values, level_id, component_id, n_levels, default = NULL){

  if(is.null(values)){
    return(default)
  }

  if(length(values) == 1L){
    value <- values[[1L]]
  }else if(length(values) == n_levels && length(level_id) == 1L && !is.na(level_id) && is.finite(level_id) && level_id >= 1L && length(values) >= level_id){
    value <- values[[level_id]]
  }else if(length(component_id) == 1L && !is.na(component_id) && is.finite(component_id) && component_id >= 1L && length(values) >= component_id){
    value <- values[[component_id]]
  }else if(length(level_id) == 1L && !is.na(level_id) && is.finite(level_id) && level_id >= 1L && length(values) >= level_id){
    value <- values[[level_id]]
  }else{
    return(default)
  }

  if(is.null(value) || length(value) == 0L || all(is.na(value))){
    return(default)
  }

  value
}
.plot_prior_factor_level_style_values <- function(values, plot_data, n_levels, default = NULL){

  out <- rep(default, n_levels)
  if(n_levels == 0L){
    return(out)
  }

  for(level_id in seq_len(n_levels)){
    component_ind <- which(vapply(plot_data, function(component){
      identical(attr(component, "level_id"), level_id)
    }, logical(1)))[1]

    component_id <- if(!is.na(component_ind)) attr(plot_data[[component_ind]], "component_id") else NA_integer_
    out[level_id] <- .plot_prior_factor_style_value(
      values       = values,
      level_id     = level_id,
      component_id = component_id,
      n_levels     = n_levels,
      default      = default
    )
  }

  out
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

  # normalize factor component metadata before rendering
  plot_data_normalized <- .plot_prior_factor_normalize_data(plot_data)
  level_names          <- plot_data_normalized[["level_names"]]
  if(!is.null(dots[["legend_labels"]])){
    if(length(dots[["legend_labels"]]) != length(level_names)){
      stop("'legend_labels' must have the same length as the number of factor levels.", call. = FALSE)
    }
    level_names <- as.character(dots[["legend_labels"]])
    plot_data_normalized[["level_names"]] <- level_names
    plot_data_normalized[["plot_data"]] <- lapply(plot_data_normalized[["plot_data"]], function(component){
      level_id <- attr(component, "level_id")
      if(length(level_id) == 1L && !is.na(level_id)){
        attr(component, "level_label") <- level_names[level_id]
      }
      component
    })
    plot_data_normalized[["points"]] <- plot_data_normalized[["plot_data"]][
      vapply(plot_data_normalized[["plot_data"]], inherits, logical(1), what = "density.prior.point")
    ]
    plot_data_normalized[["densities"]] <- plot_data_normalized[["plot_data"]][
      vapply(plot_data_normalized[["plot_data"]], function(component){
        inherits(component, "density.prior.factor") && !inherits(component, "density.prior.point")
      }, logical(1))
    ]
  }
  plot_data_points     <- plot_data_normalized[["points"]]
  plot_data_factors    <- plot_data_normalized[["densities"]]
  style_components     <- if(length(plot_data_factors) > 0L) plot_data_factors else plot_data_points


  # prepare legend information
  if(is.null(dots[["legend"]])){
    dots[["legend"]] <- TRUE
  }else{
    check_bool(dots[["legend"]], "legend", allow_NA = FALSE)
  }
  draw_legend <- isTRUE(dots[["legend"]])
  if(identical(dots[["legend_position"]], "none")){
    draw_legend <- FALSE
  }
  dots[["legend"]] <- draw_legend

  if(is.null(dots[["col"]]) & (is.null(dots[["lty"]]) | is.null(dots[["linetype"]]))){
    dots$col <- grDevices::palette.colors(n = length(level_names) + 1)[-1]
  }
  if(length(dots[["col"]]) == 1)      dots[["col"]]      <- rep(dots[["col"]],      length(level_names))
  if(length(dots[["lty"]]) == 1)      dots[["lty"]]      <- rep(dots[["lty"]],      length(level_names))
  if(length(dots[["linetype"]]) == 1) dots[["linetype"]] <- rep(dots[["linetype"]], length(level_names))

  level_col <- .plot_prior_factor_level_style_values(
    dots[["col"]],
    style_components,
    length(level_names),
    .plot.prior_settings()[["col"]]
  )
  level_lty <- .plot_prior_factor_level_style_values(
    dots[["lty"]],
    style_components,
    length(level_names),
    .plot.prior_settings()[["lty"]]
  )
  level_linetype <- .plot_prior_factor_level_style_values(
    dots[["linetype"]],
    style_components,
    length(level_names),
    .plot.prior_settings()[["lty"]]
  )


  if(plot_type == "base"){

    if(!add){
      .plot.prior_empty(type, dots)
    }

    # plot points
    for(i in seq_along(plot_data_points)){
      args           <- dots
      args$scale_y2  <- scale_y2
      args$plot_data <- plot_data_points[[i]]
      point_level    <- attr(plot_data_points[[i]], "level_id")
      point_component <- attr(plot_data_points[[i]], "component_id")
      args$col       <- .plot_prior_factor_style_value(
        values       = dots[["col"]],
        level_id     = point_level,
        component_id = point_component,
        n_levels     = length(level_names),
        default      = if(length(dots[["col"]]) > 1) .plot.prior_settings()[["col"]]
      )
      args$lty <- .plot_prior_factor_style_value(
        values       = dots[["lty"]],
        level_id     = point_level,
        component_id = point_component,
        n_levels     = length(level_names),
        default      = if(length(dots[["lty"]]) > 1) .plot.prior_settings()[["lty"]]
      )
      do.call(.lines.prior.point, args)
    }

    # plot factor levels
    for(i in seq_along(plot_data_factors)){
      args           <- dots
      args$plot_data <- plot_data_factors[[i]]
      args$level     <- attr(plot_data_factors[[i]], "level_id")
      if(is.na(args$level)) args$level <- i
      args$col       <- level_col
      args$lty       <- level_lty
      do.call(.lines.prior.factor, args)
    }

    if(draw_legend && length(level_names) > 0){
      graphics::legend(
        if(is.null(dots[["legend_position"]])) "topright" else dots[["legend_position"]],
        legend = level_names,
        col    = level_col,
        lty    = level_lty,
        lwd    = if(!is.null(dots[["lwd"]])) dots[["lwd"]] else rep(.plot.prior_settings()[["lwd"]], length(level_names)),
        title  = dots[["legend_title"]],
        bty    = "n")
    }

    plot <- list(scale_y2 = scale_y2)

  }else if(plot_type == "ggplot"){

    plot <- list()

    # plot points
    for(i in seq_along(plot_data_points)){
      args           <- dots
      args$scale_y2  <- scale_y2
      args$plot_data <- plot_data_points[[i]]
      point_level    <- attr(plot_data_points[[i]], "level_id")
      point_component <- attr(plot_data_points[[i]], "component_id")
      args$col       <- .plot_prior_factor_style_value(
        values       = dots[["col"]],
        level_id     = point_level,
        component_id = point_component,
        n_levels     = length(level_names),
        default      = if(length(dots[["col"]]) > 1) .plot.prior_settings()[["col"]]
      )
      args$lty <- .plot_prior_factor_style_value(
        values       = if(!is.null(dots[["linetype"]])) dots[["linetype"]] else dots[["lty"]],
        level_id     = point_level,
        component_id = point_component,
        n_levels     = length(level_names),
        default      = .plot.prior_settings()[["lty"]]
      )
      plot           <- c(plot, do.call(.geom_prior.point, args))
    }

    # plot factor levels
    if(length(plot_data_factors) > 0){
      plot_data_factors <- data.frame(
        x     = do.call(c, lapply(plot_data_factors, function(x) x$x)),
        y     = do.call(c, lapply(plot_data_factors, function(x) x$y)),
        level = do.call(c, lapply(seq_along(plot_data_factors), function(i) {
          level <- attr(plot_data_factors[[i]], "level_label")
          if(length(level) != 1L || is.na(level)){
            level <- as.character(i)
          }
          rep(level, length(plot_data_factors[[i]]$x))
        }))
      )

      args             <- dots
      args$level_names <- level_names
      args$plot_data   <- plot_data_factors
      args$col         <- level_col
      args$lty         <- if(!is.null(dots[["linetype"]])) level_linetype else level_lty
      args$linetype    <- if(!is.null(dots[["linetype"]])) level_linetype else NULL

      plot <- c(plot, do.call(.geom_prior.factors, args))
    }


    if(draw_legend && length(level_names) > 0){
      plot <- c(plot, list(ggplot2::theme(
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
