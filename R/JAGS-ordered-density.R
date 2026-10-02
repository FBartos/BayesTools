#' @title Batched numerical density kernel for ordered primitives
#' @description Compiles the active ordered totals and distinct allocation
#' coordinates into a function of a numeric draws matrix.
#' @param prior_list named list of original bound ordered priors.
#' @param allocation_chart optional named list keyed by allocation key. Each
#' entry has \code{kind = "group"} and a nonempty proper index subset \code{J},
#' or \code{kind = "pair"} and two distinct \code{indices}. Its gamma factors
#' are replaced by the corresponding Beta density of the group or pair share.
#' @return A function of \code{samples}, a numeric matrix with named coordinate
#' columns, returning one log density or log probability per row.
#' @details Total mixture components are localized independently for each row;
#' multi-slice spike totals share one inclusion state. Each total slice and
#' allocation key contributes once. Mixture-state probabilities, inactive
#' component sources, and nuisance-only chart factors are omitted because they
#' are constant along the conditional chart. This is a numerical conditional
#' kernel, not a bridge density. Expression totals raise
#' \code{BayesTools_ordered_expression_unavailable}, missing coordinates raise
#' \code{BayesTools_ordered_coordinates_unavailable}, and invalid indicators
#' raise \code{BayesTools_ordered_invalid_state}; all inherit
#' \code{BayesTools_ordered_unavailable}. Ordinary values outside numeric support
#' return \code{-Inf}. Refit when bound ordered metadata are unavailable.
#' @export
JAGS_ordered_density_kernel <- function(prior_list, allocation_chart = NULL){

  check_list(prior_list, "prior_list", check_length = 0)
  .check_prior_list_unique_names(prior_list)
  if(!all(vapply(prior_list, is.prior.ordered, logical(1)))){
    stop("'prior_list' must contain bound ordered priors only.", call. = FALSE)
  }
  .bt_validate_ordered_shared_allocations(prior_list)
  compiled <- .bt_ordered_compile_density(prior_list, allocation_chart)
  function(samples){
    samples <- .bt_deterministic_draws_matrix(samples)
    compiled$evaluate(samples)
  }
}

.bt_ordered_total_components <- function(spec, draws){

  total <- spec$total_prior
  if(!is.prior.mixture(total)){
    return(list(priors = list(total), component = rep(1L, nrow(draws))))
  }
  name <- paste0(.prior_ordered_total_name(spec$parameter), "_indicator")
  if(!name %in% colnames(draws)){
    .bt_ordered_stop(paste0("Ordered total indicator '", name,
      "' is unavailable in 'draws'. Include its monitored source coordinate."))
  }
  indicator <- draws[, name]
  supported <- if(is.prior.spike_and_slab(total)) c(0, 1) else seq_along(total)
  if(any(!is.finite(indicator)) || any(!indicator %in% supported)){
    .bt_ordered_stop(paste0("Ordered total indicator '", name,
      "' does not select a declared component."), "BayesTools_ordered_invalid_state")
  }
  component <- if(is.prior.spike_and_slab(total)){
    .bt_component_from_indicator(total, indicator)
  }else as.integer(indicator)
  list(priors = unclass(total), component = component)
}

.bt_ordered_total_has_expression <- function(total){

  .is_prior_expression(total) || (is.prior.mixture(total) && any(vapply(total,.bt_ordered_total_has_expression,logical(1))))
}

.bt_ordered_localize_total <- function(prior, total){

  prior$total <- total
  prior
}

.bt_ordered_compile_density <- function(prior_list, allocation_chart = NULL,
                                        emitted_allocations = character(),
                                        strict_bridge = FALSE){

  specs <- Map(.bt_ordered_spec, names(prior_list), prior_list)
  for(spec in specs){
    if(.bt_ordered_total_has_expression(spec$total_prior)){
      .bt_ordered_stop(paste0("The ordered density kernel for '", spec$parameter,
        "' is unavailable for expression totals. Use a supported scalar total prior."),
        "BayesTools_ordered_expression_unavailable")
    }
  }
  records <- list()
  for(spec in specs){
    for(record in spec$allocations){
      if(identical(record$spec$type, "dirichlet") && !record$key %in% c(names(records), emitted_allocations)){
        records[[record$key]] <- record
      }
    }
  }
  if(!is.null(allocation_chart)){
    if(!is.list(allocation_chart) || is.null(names(allocation_chart)) ||
       anyNA(names(allocation_chart)) || anyDuplicated(names(allocation_chart)) ||
       any(!names(allocation_chart) %in% names(records))){
      stop("'allocation_chart' must be keyed by declared sampled allocation keys.", call. = FALSE)
    }
    for(key in names(allocation_chart)){
      chart <- allocation_chart[[key]]
      D <- records[[key]]$dim
      if(!is.list(chart) || !is.character(chart$kind) || length(chart$kind)!=1L || is.na(chart$kind)){
        stop("'allocation_chart' must declare a proper group or two distinct pair indices.", call. = FALSE)
      }
      indices <- if(identical(chart$kind, "group")) chart$J else chart$indices
      if(!chart$kind %in% c("group", "pair") || !is.numeric(indices) ||
         anyNA(indices) || any(!indices %in% seq_len(D)) || anyDuplicated(indices) ||
         (identical(chart$kind, "group") && (length(indices) == 0L || length(indices) == D)) ||
         (identical(chart$kind, "pair") && length(indices) != 2L)){
        stop("'allocation_chart' must declare a proper group or two distinct pair indices.", call. = FALSE)
      }
    }
  }
  evaluate <- function(samples){
    result <- numeric(nrow(samples))
    for(spec in specs){
      states <- .bt_ordered_total_components(spec, samples)
      if(strict_bridge && !is.prior.point(spec$total_prior) && !all(spec$total_names %in% colnames(samples))){
        .bt_JAGS_marglik_missing_columns("'samples' does not contain all monitored ordered total prior parameters.")
      }
      total_values <- .bt_ordered_total_values(spec, .bt_deterministic_lookup(samples))
      for(k in unique(states$component)){
        prior <- states$priors[[k]]
        rows <- which(states$component == k)
        if(is.prior.point(prior)) next
        if(is.null(total_values)){
          .bt_ordered_stop(paste0("Ordered total coordinates for '", spec$parameter,
            "' are unavailable in 'samples'. Include the declared total source coordinates."))
        }
        for(slice in seq_along(spec$total_names)){
          result[rows] <- result[rows] + lpdf(prior, total_values[rows, slice])
        }
      }
    }
    for(key in names(records)){
      record <- records[[key]]
      if(!all(record$gamma_coordinates %in% colnames(samples))){
        if(strict_bridge) .bt_JAGS_marglik_missing_columns("'samples' does not contain all monitored ordered Dirichlet allocation parameters.")
        .bt_ordered_stop(paste0("Ordered gamma coordinates of '", record$node,
          "' are unavailable in 'samples'. Include all declared gamma source coordinates."))
      }
      eta <- samples[, record$gamma_coordinates, drop = FALSE]
      eta_sum <- rowSums(eta)
      invalid <- rowSums(!is.finite(eta))>0L | rowSums(if(strict_bridge) eta<=0 else eta<0,na.rm=TRUE)>0L |
        !is.finite(eta_sum) | eta_sum<=0
      chart <- allocation_chart[[key]]
      if(is.null(chart)){
        for(j in seq_len(record$dim)){
          result <- result + stats::dgamma(eta[, j], shape = record$spec$alpha[[j]], rate = 1, log = TRUE)
        }
      }else if(identical(chart$kind, "group")){
        J <- chart$J
        K <- setdiff(seq_len(record$dim), J)
        result <- result + stats::dbeta(rowSums(eta[, J, drop = FALSE]) / rowSums(eta),
          sum(record$spec$alpha[J]), sum(record$spec$alpha[K]), log = TRUE)
      }else{
        indices <- chart$indices
        pair_sum <- rowSums(eta[,indices,drop=FALSE])
        if(any(!invalid & pair_sum<=0)) .bt_ordered_stop("The ordered gamma-pair chart requires a positive pair sum.","BayesTools_ordered_invalid_state")
        result <- result + stats::dbeta(eta[, indices[[1L]]] / pair_sum,
          record$spec$alpha[[indices[[1L]]]], record$spec$alpha[[indices[[2L]]]], log = TRUE)
        for(j in setdiff(seq_len(record$dim), indices)){
          result <- result + stats::dgamma(eta[, j], shape = record$spec$alpha[[j]], rate = 1, log = TRUE)
        }
      }
      result[invalid] <- -Inf
    }
    result
  }
  list(evaluate = evaluate, allocation_keys = names(records))
}
