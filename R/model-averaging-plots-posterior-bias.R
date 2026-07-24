.plot_data_samples.PETPEESE       <- function(samples, x_seq, x_range, x_range_quant, n_points, transformation, transformation_arguments, transformation_settings, effect_direction = "positive"){

  check_list(samples, "samples")
  if (is.null(samples[["mu"]]) && is.null(samples[["mu_intercept"]]))
    stop("'mu' or 'mu_intercept' samples need to be present.")

  if (!is.null(samples[["bias"]])) {

    if(length(c("PET", "PEESE") %in% samples[["bias"]]) == 0)
      stop("At least one 'PET' or 'PEESE' model needs to be specified.")

    # create mu-PET-PEESE samples matrix
    new_samples <- matrix(if(!is.null(samples[["mu"]])) samples[["mu"]] else samples[["mu_intercept"]], ncol = 1)
    for (par in c("PET", "PEESE")) {
      if (is.element(par, colnames(samples[["bias"]]))) {
        new_samples <- cbind(new_samples, samples[["bias"]][,par])
      } else {
        new_samples <- cbind(new_samples, 0)
      }
    }

  } else {

    if(is.null(samples[["PET"]]) & is.null(samples[["PEESE"]]))
      stop("At least one 'PET' or 'PEESE' model needs to be specified.")

    # create mu-PET-PEESE samples matrix
    new_samples <- matrix(if(!is.null(samples[["mu"]])) samples[["mu"]] else samples[["mu_intercept"]], ncol = 1)
    for (par in c("PET", "PEESE")) {
      if (!is.null(samples[[par]])) {
        new_samples <- cbind(new_samples, samples[[par]])
      } else {
        new_samples <- cbind(new_samples, 0)
      }
    }
  }

  # get the plotting range
  if(is.null(x_range)){
    x_range <- c(0, 1)
  }
  if(is.null(x_seq)){
    x_seq   <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  summary <- .petpeese_line_summary_from_samples(
    samples                  = new_samples,
    x_seq                    = x_seq,
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    effect_direction         = effect_direction
  )


  out <- list(
    call    = call("density", "PET-PEESE list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = summary$median,
    y_lCI   = summary$lCI,
    y_uCI   = summary$uCI,
    samples = summary$samples
  )


  class(out) <- c("density", "density.prior", "density.prior.PETPEESE")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- range(out$y, out$y_lCI, out$y_uCI)

  return(out)
}
.petpeese_line_summary_from_samples <- function(samples, x_seq, transformation, transformation_arguments,
                                                effect_direction = "positive"){

  samples <- as.matrix(samples)
  if(ncol(samples) != 3){
    stop("'samples' must contain mu, PET, and PEESE columns.", call. = FALSE)
  }

  direction_sign <- if(effect_direction == "negative") -1 else 1
  n_samples <- nrow(samples)

  x_sam <- matrix(samples[,1], nrow = n_samples, ncol = length(x_seq)) +
    direction_sign * matrix(samples[,2], nrow = n_samples, ncol = length(x_seq)) *
      matrix(x_seq, nrow = n_samples, ncol = length(x_seq), byrow = TRUE) +
    direction_sign * matrix(samples[,3], nrow = n_samples, ncol = length(x_seq)) *
      matrix(x_seq^2, nrow = n_samples, ncol = length(x_seq), byrow = TRUE)

  if(!is.null(transformation)){
    x_sam <- .density.prior_transformation_x(x_sam, transformation, transformation_arguments)
  }

  quantiles <- apply(x_sam, 2, stats::quantile, probs = c(.500, .025, .975), names = FALSE)
  quantiles <- matrix(quantiles, nrow = 3)

  list(
    median  = quantiles[1,],
    lCI     = quantiles[2,],
    uCI     = quantiles[3,],
    samples = x_sam
  )
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

  context    <- .weightfunction_prior_list_context(prior_list)
  omega_cuts <- context$omega_cuts

  omega_columns <- grepl("^omega\\[", colnames(samples))
  if(any(omega_columns)){
    samples <- samples[, omega_columns, drop = FALSE]
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
  attr(out, "y_range") <- c(0, max(1, x_mean, x_lCI, x_uCI, na.rm = TRUE))

  return(out)
}

.weightfunction_plot_data_pvalues <- function(data, show_data){

  if(!show_data){
    return(numeric())
  }
  if(is.null(data)){
    stop("'data' must be supplied when 'show_data = TRUE'.", call. = FALSE)
  }

  if(is.data.frame(data)){
    if(!"p" %in% names(data)){
      stop("'data' must be a numeric vector of p-values or a data frame with a 'p' column.", call. = FALSE)
    }
    data <- data[["p"]]
  }

  check_real(data, "data", lower = 0, upper = 1, check_length = 0, allow_NA = FALSE)
  data <- data[is.finite(data)]

  return(data)
}

.weightfunction_plot_data_x <- function(p, plot_data, rescale_x){

  if(length(p) == 0L || !rescale_x){
    return(p)
  }

  cuts <- unique(plot_data$x)
  scaled_cuts <- seq(0, 1, length.out = length(cuts))

  stats::approx(
    x    = cuts,
    y    = scaled_cuts,
    xout = p,
    rule = 2,
    ties = "ordered"
  )$y
}

.weightfunction_data_col <- function(dots_data){

  col <- if(!is.null(dots_data[["col"]])){
    dots_data[["col"]]
  }else if(!is.null(dots_data[["color"]])){
    dots_data[["color"]]
  }else{
    "black"
  }

  if(!is.null(dots_data[["alpha"]])){
    col <- grDevices::adjustcolor(col, alpha.f = dots_data[["alpha"]])
  }

  col
}

.weightfunction_data_lwd <- function(dots_data){

  if(!is.null(dots_data[["lwd"]])){
    dots_data[["lwd"]]
  }else if(!is.null(dots_data[["linewidth"]])){
    dots_data[["linewidth"]]
  }else if(!is.null(dots_data[["size"]])){
    dots_data[["size"]]
  }else{
    .5
  }
}

.weightfunction_data_side <- function(dots_data, ggplot = FALSE){

  side <- if(!is.null(dots_data[["side"]])){
    dots_data[["side"]]
  }else if(!is.null(dots_data[["rug_side"]])){
    dots_data[["rug_side"]]
  }else{
    if(ggplot) "b" else 1
  }

  if(ggplot){
    side <- as.character(side)
    if(side %in% c("1", "bottom", "b")){
      return("b")
    }
    if(side %in% c("3", "top", "t")){
      return("t")
    }
    return(side)
  }

  if(is.character(side)){
    if(side %in% c("bottom", "b")){
      return(1)
    }
    if(side %in% c("top", "t")){
      return(3)
    }
  }

  side
}

.weightfunction_data_height <- function(dots_data){

  if(!is.null(dots_data[["height"]])){
    dots_data[["height"]]
  }else if(!is.null(dots_data[["rug_height"]])){
    dots_data[["rug_height"]]
  }else if(!is.null(dots_data[["ticksize"]])){
    dots_data[["ticksize"]]
  }else{
    .03
  }
}

.lines.weightfunction_data <- function(p, plot_data, rescale_x, dots_data = list()){

  if(length(p) == 0L){
    return(invisible())
  }

  p <- .weightfunction_plot_data_x(p, plot_data, rescale_x)

  graphics::rug(
    p,
    side     = .weightfunction_data_side(dots_data, ggplot = FALSE),
    ticksize = .weightfunction_data_height(dots_data),
    col      = .weightfunction_data_col(dots_data),
    lwd      = .weightfunction_data_lwd(dots_data)
  )

  return(invisible())
}

.geom.weightfunction_data <- function(p, plot_data, rescale_x, dots_data = list()){

  if(length(p) == 0L){
    return(NULL)
  }

  p <- .weightfunction_plot_data_x(p, plot_data, rescale_x)

  ggplot2::geom_rug(
    data        = data.frame(p = p),
    mapping     = ggplot2::aes(x = .data[["p"]]),
    inherit.aes = FALSE,
    sides       = .weightfunction_data_side(dots_data, ggplot = TRUE),
    length      = grid::unit(.weightfunction_data_height(dots_data), "npc"),
    color       = .weightfunction_data_col(dots_data),
    linewidth   = .weightfunction_data_lwd(dots_data)
  )
}
