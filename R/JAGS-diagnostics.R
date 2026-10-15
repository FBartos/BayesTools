#' @title Plot diagnostics of a 'JAGS' model
#'
#' @description Creates density plots, trace plots, and autocorrelation plots
#' for a given parameter of a JAGS model.
#'
#' @param fit a JAGS model fitted via [JAGS_fit()] or a derived-draw view from
#' [JAGS_with_draws()]. Diagnostics describe the supplied view draws.
#' @param parameter parameter to be plotted
#' @param type what type of model diagnostic should be plotted. The available
#' options are \code{"density"}, \code{"trace"}, and \code{"autocorrelation"}
#' @param plot_type whether to use a base plot \code{"base"}
#' or ggplot2 \code{"ggplot"} for plotting.
#' @param xlim x plotting range
#' @param ylim y plotting range
#' @param lags number of lags to be shown for the autocorrelation plot.
#' Defaults to \code{30}.
#' @details The trace shows the retained MCMC iteration of every draw on its
#' x axis (thinned chains keep their original positions).
#' Autocorrelation uses the retained original iteration positions.
#' Requested lags without finite observation pairs, beyond the retained chain,
#' or with nonfinite autocorrelation are retained as \code{NA}. One warning per
#' parameter and chain has classes \code{BayesTools_autocorrelation_unavailable}
#' and \code{BayesTools_plot_condition}, with \code{parameter}, \code{chain},
#' \code{unavailable_lags}, \code{pair_count} and \code{reason} fields.
#' Computed x/y objects retain full-length \code{pair_count} and
#' \code{unavailability_reason} attributes; both renderers display finite bars
#' only. No iteration compression, replacement estimate or confidence band is
#' introduced.
#'
#' Weight bins that the prior fixes (for example by [wf_fixed()]) are constants
#' of the prior. The density and autocorrelation diagnostics omit them, and
#' when nothing else remains the call returns \code{NULL} with a message, as for
#' a spike prior; the trace is produced as for any other parameter.
#' @param ... additional arguments
#' @inheritParams density.prior
#' @inheritParams plot.prior
#' @inheritParams print.prior
#' @inheritParams BayesTools_model_tables
#'
#' @seealso [JAGS_fit()] [JAGS_check_convergence()]
#'
#' @return \code{diagnostics} returns either \code{NULL} if \code{plot_type = "base"}
#' or an object/list of objects (depending on the number of parameters to be plotted)
#' of class 'ggplot' if \code{plot_type = "ggplot"}.
#'
#' @name JAGS_diagnostics
#' @export JAGS_diagnostics

#' @rdname JAGS_diagnostics
JAGS_diagnostics                 <- function(fit, parameter, type, plot_type = "base",
                                             xlim = NULL, ylim = NULL, lags = 30, n_points = 1000,
                                             transformations = NULL, transform_factors = FALSE,
                                             short_name = FALSE, parameter_names = FALSE, formula_prefix = TRUE, ...){

  if("transform_orthonormal" %in% ...names()){
    stop("The 'transform_orthonormal' argument was removed; use 'transform_factors = TRUE'.",
         call. = FALSE)
  }

  # check fits
  if(!.bt_is_jags_analysis_fit(fit))
    stop("'fit' must be a runjags fit or a derived-draw view.", call. = FALSE)
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit.", call. = FALSE)

  check_char(type, "type", allow_values = c("density", "trace", "autocorrelation"))
  check_char(plot_type, "plot_type", allow_values = c("base", "ggplot"))
  if(type == "autocorrelation"){
    check_int(
      lags,
      "lags",
      lower = 0,
      allow_NA = FALSE
    )
  }
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(parameter, "parameter", allow_values = c(names(prior_list), if(any(names(prior_list) == "bias")) c("PET", "PEESE", "omega", "alpha", "pi_null")))
  if(!is.null(transformations))
    check_char(names(transformations), "names(transformations)", allow_values = parameter)
  if(parameter == "omega" && !is.null(prior_list[["bias"]]) && !.selection_prior_has_selection(prior_list[["bias"]]))
    stop("'omega' diagnostics require at least one weightfunction component in the 'bias' prior.", call. = FALSE)
  if(parameter %in% c("PET", "PEESE") && is.null(prior_list[[parameter]]) && !is.null(prior_list[["bias"]])){
    bias_components <- if(is.prior.mixture(prior_list[["bias"]])) unclass(prior_list[["bias"]]) else list(prior_list[["bias"]])
    has_branch <- if(parameter == "PET") is.prior.PET else is.prior.PEESE
    if(!any(vapply(bias_components, has_branch, logical(1))))
      stop("'", parameter, "' diagnostics require at least one ", parameter, " component in the 'bias' prior.", call. = FALSE)
  }
  if(parameter %in% c("alpha", "pi_null") && !is.null(prior_list[["bias"]]) && !.selection_prior_has_phacking(prior_list[["bias"]]))
    stop("'alpha' and 'pi_null' diagnostics require at least one p-hacking component in the 'bias' prior.", call. = FALSE)
  .bt_require_fit_contract(fit)

  # do not produce diagnostics for a spike prior
  if(is.prior.point(prior_list[[parameter]])){
    message("No diagnostic plots are produced for a spike prior distribution")
    return(NULL)
  }

  # prepare the plot data
  plot_data <- .diagnostics_plot_data(fit = fit, parameter = parameter, prior_list = prior_list, transformations = transformations, transform_factors = transform_factors)

  # columns declared constant by the prior (fixed weightfunction bins) have no
  # density or autocorrelation: with nothing else to show, behave as for a spike
  structural <- attr(plot_data, "structural", exact = TRUE)
  if(type != "trace" && length(structural) > 0L && all(structural)){
    message("No diagnostic plots are produced for a spike prior distribution")
    return(NULL)
  }

  plot_data <- switch(
    type,
    "density"         = .diagnostics_plot_data_density(plot_data, n_points, xlim),
    "trace"           = .diagnostics_plot_data_trace(plot_data, n_points, ylim),
    "autocorrelation" = .diagnostics_plot_data_autocorrelation(plot_data, n_points, lags)
  )


  # prepare nice parameter names
  display_prior <- if(parameter %in% c("PET", "PEESE", "omega", "alpha", "pi_null") && !is.null(prior_list[["bias"]])){
    prior_list[["bias"]]
  }else{
    prior_list[[parameter]]
  }
  parameter_name <- .bt_label(
    attr(plot_data, "label_parts", exact = TRUE),
    style          = "table",
    formula_prefix = formula_prefix
  )

  # get default plot settings
  dots      <- list(...)
  main      <- if(parameter %in% c("alpha", "pi_null")){
    parameter
  }else{
    print(display_prior, plot = TRUE, short_name = short_name, parameter_names = parameter_names)
  }
  xlab      <- switch(
    type,
    "density"         = parameter_name,
    "trace"           = "Iteration",
    "autocorrelation" = "Lag"
  )
  ylab      <- switch(
    type,
    "density"         = "Density",
    "trace"           = parameter_name,
    "autocorrelation" = paste0("Autocorrelation(", parameter_name, ")")
  )


  if(is.null(dots[["main"]])) dots$main <-  main
  if(is.null(dots[["xlab"]])) dots$xlab <-  xlab
  if(is.null(dots[["ylab"]])) dots$ylab <-  ylab

  if(is.null(dots[["col"]]))  dots$col  <-  "black"
  if(is.null(dots[["lty"]]))  dots$lty  <-  1
  if(is.null(dots[["lwd"]]))  dots$lwd  <-  1

  if(length(dots[["col"]]) == 1)  dots$col  <- rep(dots$col, attr(plot_data, "chains"))
  if(length(dots[["lty"]]) == 1)  dots$lty  <- rep(dots$lty, attr(plot_data, "chains"))
  if(length(dots[["lwd"]]) == 1)  dots$lwd  <- rep(dots$lwd, attr(plot_data, "chains"))

  plots <- list()

  for(i in seq_along(plot_data)){

    temp_dots <- dots

    temp_dots[["xlab"]] <- if(length(temp_dots[["xlab"]]) > 1) temp_dots[["xlab"]][i] else temp_dots[["xlab"]]
    temp_dots[["ylab"]] <- if(length(temp_dots[["ylab"]]) > 1) temp_dots[["ylab"]][i] else temp_dots[["ylab"]]

    if(is.null(temp_dots[["xlim"]])) temp_dots$xlim <-  attr(plot_data[[i]], "x_range")
    if(is.null(temp_dots[["ylim"]])) temp_dots$ylim <-  attr(plot_data[[i]], "y_range", exact = TRUE)


    if(plot_type == "ggplot"){
      plots[[i]] <- .ggplot.prior_empty("simple", temp_dots)
    }else{
      .plot.prior_empty("simple", temp_dots)
    }

    for(j in seq_along(plot_data[[i]])){

      temp_args <- list()
      temp_args[["plot_data"]] <- plot_data[[i]][[j]]
      temp_args[["col"]]       <- dots[["col"]][j]
      temp_args[["lwd"]]       <- dots[["lwd"]][j]
      temp_args[["lty"]]       <- dots[["lty"]][j]

      if(plot_type == "ggplot"){
        plots[[i]] <- plots[[i]] + do.call(switch(
          type,
          "density"         = .geom_diagnostics.density,
          "trace"           = .geom_diagnostics.trace,
          "autocorrelation" = .geom_diagnostics.autocorrelation
        ), temp_args)
      }else{
        do.call(switch(
          type,
          "density"         = .lines_diagnostics.density,
          "trace"           = .lines_diagnostics.trace,
          "autocorrelation" = .lines_diagnostics.autocorrelation
        ), temp_args)
      }

    }
  }


  # return the plots
  if(plot_type == "base"){
    return(invisible())
  }else if(plot_type == "ggplot"){
    if(length(plots) == 1){
      plots <- plots[[1]]
    }
    return(plots)
  }


}

.diagnostics_plot_data                 <- function(fit, parameter, prior_list, transformations, transform_factors){

  check_list(transformations, "transformations", allow_NULL = TRUE)
  if(!is.null(transformations) && any(!sapply(transformations, function(trans)is.function(trans[["fun"]]))))
    stop("'transformations' must be list of functions in the 'fun' element.")

  model_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  samples_chain <- lapply(seq_along(model_samples_list), function(i) {
    return(rep(i, nrow(model_samples_list[[i]])))
  })
  # the retained MCMC iteration of every draw (thinned and started chains
  # keep their original positions)
  samples_iter  <- lapply(seq_along(model_samples_list), function(i) {
    mcpar <- attr(model_samples_list[[i]], "mcpar", exact = TRUE)
    n     <- nrow(model_samples_list[[i]])
    if(is.numeric(mcpar) && length(mcpar) == 3L && all(is.finite(mcpar))){
      return(seq(from = mcpar[1L], by = mcpar[3L], length.out = n))
    }
    return(seq_len(n))
  })
  model_samples <- do.call(rbind, model_samples_list)

  # the transformation of the requested parameter applies to the columns the
  # call dispatches it to when the requested name is none of them (the omega
  # bins of a weightfunction, the PET or PEESE column of a bias prior)
  requested_parameter  <- parameter
  dispatched_to_columns <- FALSE

  # extract the relevant parameters
  if(is.prior.spike_and_slab(prior_list[[parameter]])){

    if(sum(model_samples[,colnames(model_samples) == paste0(parameter, "_indicator")] != 0) < 10)
      stop("The parameter with a spike and slab prior did not result in enough samples under the slab for producing a diagnostic figure.")

    # change the samples between conditional: the excluded rows of every
    # coefficient column of the parameter (the indexed columns of a factor)
    model_samples[
      model_samples[,colnames(model_samples) == paste0(parameter, "_indicator")] == 0,
      colnames(model_samples) == parameter | startsWith(colnames(model_samples), paste0(parameter, "["))] <- NA

    # modify the parameter list
    prior_list[[parameter]] <- .get_spike_and_slab_variable(prior_list[[parameter]])
  }

  # the prior giving the support of the diagnostics: the whole mixture of a
  # mixture prior (below, a factor mixture is reduced to a factor component,
  # which only resolves the dimensions and the coding of the columns)
  prior_support <- prior_list[[parameter]]

  # the draws of the point components of a mixture that also has a continuous
  # component are atoms, not part of its density: the kernel estimate drops
  # them like the spike of a spike-and-slab prior
  atom_rows <- NULL
  if(is.prior.mixture(prior_support) && !inherits(prior_support, "prior.bias_mixture")){
    is_atom <- vapply(prior_support, is.prior.point, logical(1))
    indicator_name <- paste0(parameter, "_indicator")
    if(any(is_atom) && !all(is_atom) && indicator_name %in% colnames(model_samples)){
      atom_rows <- model_samples[, indicator_name] %in% which(is_atom)
    }
  }

  if(is.prior.factor(prior_list[[parameter]]) || (is.prior.mixture(prior_list[[parameter]]) && any(sapply(prior_list[[parameter]], is.prior.factor)))){

    if(is.prior.mixture(prior_list[[parameter]])){
      prior_list[[parameter]] <- prior_list[[parameter]][[which(sapply(prior_list[[parameter]], is.prior.factor))[1]]]
    }

    if(.get_prior_factor_levels(prior_list[[parameter]]) > 1){
      model_samples <- model_samples[,paste0(parameter, "[", 1:.get_prior_factor_levels(prior_list[[parameter]]), "]"),drop = FALSE]
    }else{
      model_samples <- model_samples[,parameter,drop = FALSE]
    }
  }else if(is.prior.weightfunction(prior_list[[parameter]])){
    omega_cuts <- weightfunctions_mapping(list(prior_list[[parameter]]), cuts_only = TRUE)
    model_samples <- model_samples[,paste0("omega", "[", 2:(length(omega_cuts) - 1), "]"),drop = FALSE]
    dispatched_to_columns <- TRUE
  }else if(parameter == "omega" && !is.null(prior_list[["bias"]])){
    selection_priors <- .selection_prior_selection_priors(prior_list[["bias"]])
    model_samples <- model_samples[,paste0("omega", "[", 2:(length(weightfunctions_mapping(selection_priors, cuts_only = TRUE, one_sided = TRUE)) - 1), "]"),drop = FALSE]
    dispatched_to_columns <- TRUE
  }else if(is_prior_phacking(prior_list[[parameter]])){
    model_samples <- model_samples[,intersect(.phacking_report_parameter(prior_list[[parameter]]), colnames(model_samples)),drop = FALSE]
  }else if(is_prior_bias(prior_list[[parameter]]) || inherits(prior_list[[parameter]], "prior.bias_mixture")){
    bias_parameters <- .selection_bias_parameter_names(prior_list[[parameter]], include_kind = FALSE)
    bias_parameters <- bias_parameters[bias_parameters != "omega"]
    if("omega" %in% .selection_bias_parameter_names(prior_list[[parameter]], include_kind = FALSE)){
      omega_cuts <- weightfunctions_mapping(.selection_prior_selection_priors(prior_list[[parameter]]), cuts_only = TRUE, one_sided = TRUE)
      bias_parameters <- c(paste0("omega", "[", 2:(length(omega_cuts) - 1), "]"), bias_parameters)
    }
    model_samples <- model_samples[,intersect(bias_parameters, colnames(model_samples)),drop = FALSE]
  }else if(parameter %in% c("alpha", "pi_null") && !is.null(prior_list[["bias"]])){
    model_samples <- model_samples[,parameter,drop = FALSE]
  }else if(is.prior.vector(prior_list[[parameter]])){
    if(prior_list[[parameter]]$parameters[["K"]] > 1){
      model_samples <- model_samples[,paste0(parameter, "[", 1:prior_list[[parameter]]$parameters[["K"]], "]"),drop = FALSE]
    }else{
      model_samples <- model_samples[,parameter,drop = FALSE]
    }
  }else{
    if(is.prior.PET(prior_list[[parameter]])){
      parameter <- "PET"
    }else if(is.prior.PEESE(prior_list[[parameter]])){
      parameter <- "PEESE"
    }
    dispatched_to_columns <- !identical(parameter, requested_parameter)
    model_samples <- model_samples[,parameter,drop = FALSE]
  }

  # deal with bias mixture dispatching: a component of the bias prior takes
  # its support from the whole bias prior (a mixture or a single prior); a
  # single bias prior is wrapped in a named list, so that the blocks below
  # see a list of priors and never its fields
  if(parameter %in% c("PET", "PEESE", "omega", "alpha", "pi_null") && !is.null(prior_list[["bias"]])){
    prior_support   <- prior_list[["bias"]]
    prior_list      <- if(is.prior.mixture(prior_support)) prior_support else list(bias = prior_support)
    parameter_names <- parameter
  }else{
    prior_list      <- prior_list[parameter]
    parameter_names <- parameter
  }


  # mostly adapted from runjags_estimates_table
  # apply transformations
  if(dispatched_to_columns && !is.null(transformations[[requested_parameter]])){
    transformations <- stats::setNames(
      rep(list(transformations[[requested_parameter]]), ncol(model_samples)),
      colnames(model_samples)
    )
  }
  model_samples <- .apply_parameter_transformations(model_samples, transformations, prior_list, transform_factors)

  # transform meandif and orthonormal factors to differences
  model_samples <- .transform_factor_contrasts(model_samples, prior_list, transform_factors, transformations,
                                               with_label_parts = TRUE)
  created_parts <- attr(model_samples, "label_parts", exact = TRUE)
  attr(model_samples, "label_parts") <- NULL

  # extract parameter names from column names after transformations
  parameter_names <- colnames(model_samples)

  # the columns the prior declares constant (the bins of fixed weights) and the
  # draws agree: the density and autocorrelation diagnostics skip them
  is_structural <- .diagnostics_structural_columns(model_samples, prior_support)

  # rename weightfunctions factor levels (special case that overrides)
  if(any(sapply(prior_list, is.prior.weightfunction)) && !is.prior.mixture(prior_list)){
    for(par in names(prior_list)[sapply(prior_list, is.prior.weightfunction)]){
      omega_cuts      <- weightfunctions_mapping(prior_list[par], cuts_only = TRUE)
      parameter_names <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
      parameter_names <- parameter_names[-1]
    }
  }

  if(is.prior.mixture(prior_list) && parameter == "omega"){
    omega_cuts      <- weightfunctions_mapping(.selection_prior_selection_priors(prior_list), cuts_only = TRUE, one_sided = TRUE)
    parameter_names <- sapply(2:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
  }

  if(any(sapply(prior_list, is_prior_bias))){
    bias_prior <- prior_list[[which(sapply(prior_list, is_prior_bias))[1]]]
    if(.selection_prior_has_selection(bias_prior)){
      omega_cuts      <- weightfunctions_mapping(.selection_prior_selection_priors(bias_prior), cuts_only = TRUE, one_sided = TRUE)
      omega_names_old <- paste0("omega[", 1:(length(omega_cuts)-1), "]")
      omega_names     <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
      parameter_names[parameter_names %in% omega_names_old] <- omega_names[match(parameter_names[parameter_names %in% omega_names_old], omega_names_old)]
    }
  }

  # label every column from its label parts; the columns are named by their
  # selectors
  column_parts <- .bt_estimates_column_parts(
    columns = parameter_names,
    fit     = fit,
    created = if(is.null(created_parts)) list() else created_parts
  )
  parameter_names <- .bt_label(column_parts, style = "selector")

  # attach the relevant attributes
  colnames(model_samples)          <- parameter_names
  attr(model_samples, "label_parts") <- column_parts
  attr(model_samples, "chain")     <- do.call(c, samples_chain)
  attr(model_samples, "iter")      <- do.call(c, samples_iter)
  attr(model_samples, "atom_rows") <- atom_rows
  attr(model_samples, "parameter") <- parameter
  attr(model_samples, "prior")     <- prior_support
  attr(model_samples, "structural") <- is_structural
  # Diagnostic transformations are validated against the selected parameter;
  # reflected KDE bounds are therefore no longer on the same support scale.
  attr(model_samples, "density_support_transformed") <- !is.null(transformations)

  return(model_samples)
}
.diagnostics_plot_data_density         <- function(plot_data, n_points, xlim){

  plot_data <- .diagnostics_drop_structural_columns(plot_data)
  chain <- attr(plot_data, "chain", exact = TRUE)
  prior <- attr(plot_data, "prior", exact = TRUE)
  atom_rows <- attr(plot_data, "atom_rows", exact = TRUE)
  if(!is.null(atom_rows)){
    plot_data[atom_rows, ] <- NA
  }
  .bt_diagnostics_validate_plot_data(
    plot_data,
    chain,
    diagnostic = "Density"
  )

  out   <- list()

  for(i in 1:ncol(plot_data)){

    bounds <- .diagnostics_density_bounds(
      prior          = prior,
      parameter      = attr(plot_data, "parameter", exact = TRUE),
      parameter_name = colnames(plot_data)[i]
    )
    bounds <- c(bounds[["lower"]], bounds[["upper"]])
    if(isTRUE(attr(plot_data, "density_support_transformed"))){
      bounds <- c(-Inf, Inf)
    }

    if(is.null(xlim)){
      x_range <- range(plot_data[,i], na.rm = TRUE)
    }else{
      x_range <- xlim
    }

    for(j in seq_along(unique(chain))){

      density_continuous <- .density_kde_boundary(
        x      = plot_data[chain == j,i],
        n      = n_points,
        from   = x_range[1],
        to     = x_range[2],
        bounds = bounds,
        na.rm  = TRUE
      )

      x_den <- density_continuous$x
      y_den <- density_continuous$y
      boundary_reflection <- isTRUE(attr(density_continuous, "boundary_reflection"))

      temp_density    <- list(
        call    = call("density"),
        bw      = density_continuous$bw,
        n       = n_points,
        x       = x_den,
        y       = y_den
      )

      class(temp_density) <- c("density")
      attr(temp_density, "x_range")        <- range(x_den)
      attr(temp_density, "y_range")        <- c(0, max(y_den))
      attr(temp_density, "chain")          <- j
      attr(temp_density, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
      attr(temp_density, "parameter_name") <- colnames(plot_data)[i]
      if(boundary_reflection){
        attr(temp_density, "boundary_reflection") <- TRUE
      }

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_density
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- x_range
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- c(0, max(sapply(out[[i]], function(x) attr(x, "y_range", exact = TRUE))))
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
  attr(out, "parameter_name") <- colnames(plot_data)
  attr(out, "label_parts")    <- attr(plot_data, "label_parts", exact = TRUE)

  return(out)
}

.diagnostics_density_bounds <- function(prior, parameter, parameter_name){

  if(is_prior_bias(prior) && identical(parameter, "bias")){
    component <- sub("\\[.*$", "", parameter_name)
    if(component %in% c("omega", "alpha", "pi_null")){
      return(.diagnostics_prior_bounds(prior, component))
    }
    return(list(lower = -Inf, upper = Inf))
  }

  return(.diagnostics_prior_bounds(prior, parameter))
}

.diagnostics_prior_bounds <- function(prior, parameter){

  is_selection_prior <- .diagnostics_is_selection_prior(prior)

  if(is.prior.weightfunction(prior)){
    return(.weightfunction_weights_truncation(prior$weights))
  }

  if(is_selection_prior && parameter == "omega"){
    selection_priors <- .selection_prior_selection_priors(prior)
    if(length(selection_priors) == 0L){
      return(list(lower = 1, upper = 1))
    }
    upper <- max(vapply(selection_priors, function(x){
      .weightfunction_weights_truncation(x$weights)$upper
    }, numeric(1)))
    return(list(lower = 0, upper = upper))
  }

  if(is_selection_prior && parameter == "alpha"){
    return(list(lower = 0, upper = 1))
  }

  if(is_prior_phacking(prior) && !parameter %in% c("alpha", "pi_null")){
    parameter <- .phacking_report_parameter(prior)
  }

  if(is_selection_prior && parameter == "pi_null"){
    phacking_priors <- .selection_prior_phacking_priors(prior)
    if(length(phacking_priors) > 0L){
      upper <- max(vapply(phacking_priors, function(x){
        phack_backend_constants(x$form, x$source, x$destination, target = x$target)$pi_null_per_alpha
      }, numeric(1)))
      return(list(lower = 0, upper = upper))
    }
    return(list(lower = 0, upper = Inf))
  }

  if(is_prior_phacking(prior) || is_prior_bias(prior)){
    return(list(lower = 0, upper = 1))
  }

  if(is.prior.simplex(prior)){
    return(list(lower = 0, upper = 1))
  }

  if(is.prior.mixture(prior)){
    # the support of a mixture is the union of the continuous supports of its
    # components: a point (or a none) is an atom and no part of the support, a
    # simple component brings its truncation, and a factor component that is
    # not simple (orthonormal, meandif, ordered) the real line; components of
    # other kinds (the publication bias components) are no support of the
    # parameter
    is_atom   <- vapply(prior, function(x) is.prior.point(x) || is.prior.none(x), logical(1))
    is_simple <- vapply(prior, is.prior.simple, logical(1)) & !is_atom
    is_line   <- vapply(prior, is.prior.factor, logical(1)) & !is_atom & !is_simple
    prior_lower <- c(vapply(prior[is_simple], function(x) x$truncation[["lower"]], numeric(1)), rep(-Inf, sum(is_line)))
    prior_upper <- c(vapply(prior[is_simple], function(x) x$truncation[["upper"]], numeric(1)), rep(Inf, sum(is_line)))
    if(length(prior_lower) > 0L){
      return(list(lower = min(prior_lower), upper = max(prior_upper)))
    }

    # without a continuous component the atoms are all there is
    locations <- vapply(prior[vapply(prior, is.prior.point, logical(1))], function(x) x$parameters[["location"]], numeric(1))
    if(length(locations) == 0L){
      return(list(lower = -Inf, upper = Inf))
    }
    return(list(lower = min(locations), upper = max(locations)))
  }

  if(is.prior.point(prior)){
    return(list(lower = prior$parameters[["location"]], upper = prior$parameters[["location"]]))
  }

  return(list(lower = prior$truncation[["lower"]], upper = prior$truncation[["upper"]]))
}

.diagnostics_is_selection_prior <- function(prior){

  if(is.prior.weightfunction(prior) || is_prior_phacking(prior) || is_prior_bias(prior)){
    return(TRUE)
  }

  if(is.prior.mixture(prior)){
    return(any(vapply(prior, function(x){
      is.prior.none(x) || is.prior.PET(x) || is.prior.PEESE(x) ||
        is.prior.weightfunction(x) || is_prior_phacking(x) || is_prior_bias(x)
    }, logical(1))))
  }

  FALSE
}
.diagnostics_plot_data_trace           <- function(plot_data, n_points, ylim){

  chain <- attr(plot_data, "chain", exact = TRUE)
  iter  <- attr(plot_data, "iter", exact = TRUE)

  out   <- list()

  for(i in 1:ncol(plot_data)){

    if(is.null(ylim)){
      y_range <- range(plot_data[,i], na.rm = TRUE)
    }else{
      y_range <- ylim
    }
    x_range <- range(iter)

    for(j in seq_along(unique(chain))){

      temp_x       <- iter[chain == j]
      temp_y       <- plot_data[chain == j,i]

      temp_chain <- list(
        x       = temp_x,
        y       = temp_y
      )

      class(temp_chain) <- c("BayesTools_chain")
      attr(temp_chain, "x_range")        <- x_range
      attr(temp_chain, "y_range")        <- y_range
      attr(temp_chain, "chain")          <- j
      attr(temp_chain, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
      attr(temp_chain, "parameter_name") <- colnames(plot_data)[i]

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_chain
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- x_range
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- y_range
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
  attr(out, "parameter_name") <- colnames(plot_data)
  attr(out, "label_parts")    <- attr(plot_data, "label_parts", exact = TRUE)

  return(out)
}
.diagnostics_plot_data_autocorrelation <- function(plot_data, n_points, lags){

  plot_data <- .diagnostics_drop_structural_columns(plot_data)
  chain <- attr(plot_data, "chain", exact = TRUE)
  .bt_diagnostics_validate_plot_data(
    plot_data,
    chain,
    diagnostic = "Autocorrelation"
  )

  out   <- list()

  for(i in 1:ncol(plot_data)){

    for(j in seq_along(unique(chain))){

      values <- plot_data[chain == j,i]
      autocorrelation <- stats::acf(values, lag.max = lags, plot = FALSE, na.action = stats::na.pass)
      temp_x <- 0:lags
      temp_y <- rep(NA_real_, length(temp_x))
      returned_lags <- as.integer(round(as.numeric(autocorrelation$lag)))
      temp_y[returned_lags + 1L] <- as.numeric(autocorrelation$acf[, , 1L])
      temp_y[!is.finite(temp_y)] <- NA_real_
      n <- length(values)
      finite <- is.finite(values)
      pair_count <- vapply(temp_x, function(lag){
        if(lag >= n) return(0L)
        positions <- seq_len(n - lag)
        sum(finite[positions] & finite[lag + positions])
      }, integer(1))
      unavailable <- is.na(temp_y)
      reason <- rep(NA_character_, length(temp_x))
      reason[unavailable] <- ifelse(temp_x[unavailable] >= n, "chain_too_short",
        ifelse(pair_count[unavailable] == 0L, "no_finite_pairs", "nonfinite_autocorrelation"))
      if(any(unavailable)){
        warning(warningCondition(
          paste0("Autocorrelation diagnostics for '", colnames(plot_data)[i],
            "' in chain ", j, " are unavailable at lag(s) ",
            paste(temp_x[unavailable], collapse = ", "),
            ": the retained chain is too short, has insufficient finite observation pairs, or produced a nonfinite autocorrelation. These lag values were set to NA; original iteration spacing was retained."),
          call = NULL,
          class = c("BayesTools_autocorrelation_unavailable", "BayesTools_plot_condition"),
          parameter = colnames(plot_data)[i], chain = j,
          unavailable_lags = temp_x[unavailable], pair_count = pair_count[unavailable],
          reason = reason[unavailable]
        ))
      }

      temp_autocor <- list(
        x       = temp_x,
        y       = temp_y
      )

      class(temp_autocor) <- c("BayesTools_autocorrelation")
      attr(temp_autocor, "x_range")        <- range(temp_x)
      attr(temp_autocor, "y_range")        <- range(c(0, temp_y[is.finite(temp_y)]))
      attr(temp_autocor, "pair_count")     <- pair_count
      attr(temp_autocor, "unavailability_reason") <- reason
      attr(temp_autocor, "chain")          <- j
      attr(temp_autocor, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
      attr(temp_autocor, "parameter_name") <- colnames(plot_data)[i]

      out[[colnames(plot_data)[[i]]]][[j]] <- temp_autocor
    }

    attr(out[[colnames(plot_data)[[i]]]], "x_range")        <- range(unlist(lapply(out[[i]], function(x) attr(x, "x_range"))))
    attr(out[[colnames(plot_data)[[i]]]], "y_range")        <- range(unlist(lapply(out[[i]], function(x) attr(x, "y_range", exact = TRUE))))
    attr(out[[colnames(plot_data)[[i]]]], "chains")         <- length(unique(chain))
    attr(out[[colnames(plot_data)[[i]]]], "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
    attr(out[[colnames(plot_data)[[i]]]], "parameter_name") <- colnames(plot_data)[i]
  }

  attr(out, "chains")         <- length(unique(chain))
  attr(out, "parameter")      <- attr(plot_data, "parameter", exact = TRUE)
  attr(out, "parameter_name") <- colnames(plot_data)
  attr(out, "label_parts")    <- attr(plot_data, "label_parts", exact = TRUE)

  return(out)
}

# The columns of the diagnostic samples that the prior declares constant (the
# bins of fixed weights, see .bt_convergence_role_constant_bins()) and that are
# constant in the draws. A declared constant that varies is a sample to
# diagnose, not a structural constant.
.diagnostics_structural_columns <- function(model_samples, prior){

  structural <- rep(FALSE, ncol(model_samples))
  if(!.diagnostics_is_selection_prior(prior)){
    return(structural)
  }

  declared <- which(colnames(model_samples) %in% .bt_convergence_role_constant_bins(prior))
  structural[declared] <- vapply(declared, function(i){
    all(is.finite(model_samples[, i])) && length(unique(model_samples[, i])) == 1L
  }, logical(1))

  structural
}

# The diagnostic samples without their structural columns, with the attributes
# of the samples; a sample without the flag is returned as is.
.diagnostics_drop_structural_columns <- function(plot_data){

  structural <- attr(plot_data, "structural", exact = TRUE)
  if(!any(structural)){
    return(plot_data)
  }

  out <- plot_data[, !structural, drop = FALSE]
  for(name in setdiff(names(attributes(plot_data)), c("dim", "dimnames", "label_parts", "structural"))){
    attr(out, name) <- attr(plot_data, name, exact = TRUE)
  }
  attr(out, "label_parts") <- attr(plot_data, "label_parts", exact = TRUE)[!structural]
  attr(out, "structural")  <- structural[!structural]

  out
}

# Diagnostics that cannot be drawn from the posterior samples: class
# BayesTools_plot_condition (the family of plot conditions).
.bt_diagnostics_unavailable <- function(...){

  stop(structure(
    list(message = paste0(...), call = NULL),
    class = c("BayesTools_plot_condition", "error", "condition")
  ))
}

.bt_diagnostics_validate_plot_data <- function(
    plot_data,
    chain,
    diagnostic){

  if(ncol(plot_data) == 0L || nrow(plot_data) == 0L){
    .bt_diagnostics_unavailable(
      diagnostic,
      " diagnostics require at least one parameter with posterior samples. ",
      "Check that 'parameter' names a monitored parameter of the fit."
    )
  }
  if(length(chain) != nrow(plot_data)){
    stop(
      "Diagnostic chain metadata must identify every posterior sample.",
      call. = FALSE
    )
  }
  for(parameter in colnames(plot_data)){
    for(chain_id in unique(chain)){
      values <- plot_data[chain == chain_id, parameter]
      finite_values <- values[is.finite(values)]
      if(length(finite_values) < 2L){
        .bt_diagnostics_unavailable(
          diagnostic,
          " diagnostics for '", parameter, "' in chain ", chain_id,
          " require at least two finite posterior samples. ",
          "Use type = \"trace\" to inspect the chain."
        )
      }
      if(length(unique(finite_values)) < 2L){
        .bt_diagnostics_unavailable(
          diagnostic,
          " diagnostics for '", parameter, "' in chain ", chain_id,
          " are unavailable because the posterior samples are constant. ",
          "Use type = \"trace\" to inspect the chain."
        )
      }
    }
  }
  invisible(TRUE)
}


.lines_diagnostics.density         <- function(plot_data, ...){
  .lines.prior.simple(plot_data, ...)
}
.lines_diagnostics.trace           <- .lines_diagnostics.density
.lines_diagnostics.autocorrelation <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]])) dots[["col"]] else .plot.prior_settings()[["col"]]
  finite <- is.finite(plot_data$y)

  graphics::rect(
    xleft   = plot_data$x[finite] + 0.075 - 0.5,
    ybottom = 0,
    xright  = plot_data$x[finite] + 0.925 - 0.5,
    ytop    = plot_data$y[finite],
    col     = col)

  return(invisible())
}

.geom_diagnostics.density         <- function(plot_data, ...){
  .geom_prior.simple(plot_data, ...)
}
.geom_diagnostics.trace           <- .geom_diagnostics.density
.geom_diagnostics.autocorrelation <- function(plot_data, ...){

  dots      <- list(...)
  col       <- if(!is.null(dots[["col"]])) dots[["col"]] else .plot.prior_settings()[["col"]]
  finite <- is.finite(plot_data$y)

  geom <- ggplot2::geom_bar(
    # Retain the unit-lag width when finite display lags are widely separated.
    width = 0.9,
    data    = data.frame(
      x = plot_data$x[finite],
      y = plot_data$y[finite]),
    mapping = ggplot2::aes(
      x      = .data[["x"]],
      weight = .data[["y"]]),
    color = col, fill = col)

  return(geom)
}
