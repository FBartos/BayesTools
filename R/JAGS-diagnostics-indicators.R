# Indicator traces use the recorded component of a draw. In particular a
# component can carry an exact zero while being part of the alternative.
.diagnostics_indicator_data <- function(fit, parameter, prior, statistic, batch_size){

  indicator <- paste0(parameter, "_indicator")
  label_parts <- .bt_label_parts_inclusion_term(parameter, prior)
  if(.posterior_components_is_ordered_mixture(prior)){
    indicator <- paste0(.prior_ordered_total_name(parameter), "_indicator")
    label_parts <- .bt_ordered_total_label_parts(parameter, prior)[[indicator]]
    prior <- prior$total
  }
  # The compiler samples every publication-bias mixture as 'bias_indicator'.
  if(inherits(prior, "prior.bias_mixture")) indicator <- "bias_indicator"
  gate <- .bt_is_random_allocation_inclusion_prior(prior)
  if(gate){
    indicator <- attr(prior, "random_allocation_indicator", exact = TRUE)
    label_parts <- .bt_random_allocation_inference_parts(fit, indicator)
    groups <- c("null", "alternative")
  }else if(is.prior.mixture(prior) && length(prior) > 1L){
    groups <- attr(prior, "components", exact = TRUE)
  }else{
    .bt_diagnostics_unavailable("No sampled component indicator is declared for '", parameter, "': ",
                                "indicator traces need a mixture or spike-and-slab prior with at least ",
                                "two components, or a random-effect allocation gate.")
  }
  has_groups <- setequal(groups, c("null", "alternative"))
  if(statistic == "inclusion" && !has_groups){
    .bt_diagnostics_unavailable("Inclusion for '", parameter, "' requires every component to be declared ",
                                "null or alternative, with both present. ",
                                "Use 'statistic = \"component\"' to inspect individual components.")
  }
  components <- seq_along(groups)
  alternative <- which(groups == "alternative")
  if(statistic == "component"){
    if(length(groups) == 2L && has_groups) components <- alternative
    panel_names <- if(gate) "Included" else .diagnostics_component_names(prior)[components]
  }else{
    panel_names <- parameter
  }

  chains <- .extract_posterior_samples(fit, as_list = TRUE)
  n_draws <- vapply(chains, nrow, integer(1))
  min_size <- if(statistic == "switching") 2L else 1L
  if(!length(chains) || any(n_draws < min_size)){
    .bt_diagnostics_unavailable("The '", statistic, "' trace for '", parameter,
                                "' requires at least ", min_size, " retained draw(s) in each chain.")
  }
  shortest <- min(n_draws)
  default_width <- shortest / min(shortest, max(2L, round(sqrt(shortest))))
  rows <- vector("list", length(chains))
  midpoints <- vector("list", length(chains))
  iter_range <- c(Inf, -Inf)
  for(j in seq_along(chains)){
    samples <- chains[[j]]
    if(!indicator %in% colnames(samples)){
      .bt_diagnostics_unavailable("The indicator '", indicator, "' is unavailable in chain ", j, ".")
    }
    states <- samples[, indicator]
    if(gate){
      if(any(!states %in% 0:1)){
        .bt_diagnostics_unavailable("Allocation indicator '", indicator, "' must contain only 0 and 1.")
      }
      states <- as.integer(states) + 1L
    }else{
      states <- tryCatch(.bt_component_from_indicator(prior, states), error = function(e){
        .bt_diagnostics_unavailable("Invalid draws for '", indicator, "' in chain ", j, ". ", conditionMessage(e))
      })
    }
    n <- length(states)
    pos <- seq_len(n) - 1L
    if(is.null(batch_size)){
      batches <- min(max(1L, round(n / default_width)), floor(n / min_size))
      batch <- floor(pos * as.double(batches) / n)
    }else{
      batch <- pos %/% batch_size
    }
    draws_by_batch <- split(states, batch)
    iterations <- .diagnostics_iterations(samples)
    iter_range <- range(iter_range[is.finite(iter_range)], iterations)
    midpoints[[j]] <- as.numeric(tapply(iterations, batch, mean))
    rows[[j]] <- do.call(rbind, lapply(draws_by_batch, function(z){
      switch(statistic,
        inclusion = mean(z %in% alternative),
        component = tabulate(z, nbins = length(groups))[components] / length(z),
        switching = if(length(z) > 1L) mean(diff(z) != 0L) else NA_real_
      )
    }))
  }
  out <- do.call(rbind, rows)
  colnames(out) <- panel_names
  attr(out, "iter") <- unlist(midpoints, use.names = FALSE)
  attr(out, "iter_range") <- iter_range
  attr(out, "chain") <- rep(seq_along(chains), lengths(midpoints))
  attr(out, "parameter") <- parameter
  attr(out, "label_parts") <- label_parts
  out
}

.diagnostics_component_names <- function(prior){

  labels <- names(prior)
  if(is.null(labels)) labels <- rep("", length(prior))
  labels[is.na(labels)] <- ""
  groups <- attr(prior, "components", exact = TRUE)
  custom <- !is.na(groups) & !groups %in% c("null", "alternative") &
    !duplicated(groups) & !duplicated(groups, fromLast = TRUE)
  use_group <- !nzchar(labels) & custom
  labels[use_group] <- groups[use_group]
  missing <- !nzchar(labels)
  labels[missing] <- vapply(prior[missing], function(branch){
    print(branch, silent = TRUE)
  }, character(1))
  make.unique(labels)
}
