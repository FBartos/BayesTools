.plot_data_prior_list.simple         <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples, force_samples, individual,
                                                 transformation, transformation_arguments, transformation_settings){

  if(is.prior.spike_and_slab(prior_list))
    prior_list <- list(prior_list)

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

  prior_weights  <- sapply(prior_list, .prior_model_weight)
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
    }else if(inherits(plot_data[[i]], "density.prior.simple") |
             inherits(plot_data[[i]], "density.prior.orthonormal") |
             inherits(plot_data[[i]], "density.prior.meandif") |
             inherits(plot_data[[i]], "density.prior.PET") |
             inherits(plot_data[[i]], "density.prior.PEESE")){
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
  prior_weights <- unname(sapply(prior_list, .prior_model_weight))
  to_remove  <- NULL
  for(i in 1:nrow(are_equal)){
    this_ind    <- c(1:ncol(are_equal))[are_equal[i,]]
    this_unique <- this_ind[1]
    prior_weights[this_unique] <- sum(prior_weights[this_ind])
    to_remove   <- c(to_remove, this_ind[-1])
  }

  # return prior odds
  for(i in seq_along(prior_list)){
    prior_list[[i]] <- .set_prior_model_weight(prior_list[[i]], prior_weights[i])
  }

  # remove the duplicates
  prior_list[to_remove] <- NULL

  return(prior_list)
}


