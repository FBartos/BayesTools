# The readiness check between burn-in and sampling of JAGS_fit(): the chains
# start sampling once they agree over the last block of burn-in. The block runs
# with the fit's monitors on; its draws are checked and dropped.

# The split-chain point estimate of the potential scale reduction factor of one
# quantity over a block: each chain's block is split into its first and last
# half (the middle draw of an odd block is left out), and the Gelman-Rubin point
# estimate is computed over the 2 x chains half-chains. Degenerate blocks are
# decided on the half-chains before coda, which returns NaN for some of them:
# NA (not assessed) when every half-chain has a value that is not finite or
# all values are equal, Inf when only some half-chains have a value that is
# not finite (the chains do not agree) or every half-chain is constant but
# they differ, and 1 when the half-chains vary and have equal means.
.bt_burnin_split_Rhat <- function(values, chain_id, n_chains){

  halves <- unlist(lapply(seq_len(n_chains), function(chain){
    draws <- values[chain_id == chain]
    half  <- length(draws) %/% 2L
    list(draws[seq_len(half)], draws[length(draws) - half + seq_len(half)])
  }), recursive = FALSE)
  finite <- vapply(halves, function(half) all(is.finite(half)), logical(1))
  if(!any(finite)){
    return(NA_real_)
  }
  if(!all(finite)){
    return(Inf)
  }
  if(length(unique(unlist(halves))) < 2L){
    return(NA_real_)
  }
  if(all(vapply(halves, function(half) all(half == half[[1L]]), logical(1)))){
    return(Inf)
  }
  # equal means as sums over lengths: state occupancies sum exactly in any
  # order, while mean() refines its result in extended precision where the
  # platform has it and otherwise (Apple silicon) depends on the order. The
  # values are divided by their largest magnitude first (1 for occupancies),
  # so that the sums of huge values stay finite and do not compare equal as
  # infinities.
  scale <- max(abs(unlist(halves, use.names = FALSE)))
  means <- vapply(halves, function(half) sum(half / scale) / length(half), numeric(1))
  if(length(unique(means)) == 1L){
    return(1)
  }
  estimate <- tryCatch(
    suppressWarnings(coda::gelman.diag(
      coda::mcmc.list(lapply(halves, coda::mcmc)),
      autoburnin = FALSE,
      multivariate = FALSE
    )[["psrf"]][1L, 1L]),
    error = function(e) NA_real_
  )
  if(is.nan(estimate)) NA_real_ else estimate
}

# The split-chain R-hat of every quantity that the fit's convergence check
# would assess (the same selection: 'monitor', 'check_indicators' and the
# convergence roles of the coordinates; indicators as state occupancies) in a
# block of draws (an mcmc.list of the chains). Quantities that cannot be
# assessed are left out.
.bt_burnin_block_Rhat <- function(block, coordinates, prior_list, monitor, check_indicators){

  block_fit <- structure(list(mcmc = block), class = "BayesTools_fit")
  attr(block_fit, "prior_list") <- prior_list
  prepared <- .bt_convergence_prepare(
    fit = block_fit,
    coordinates = coordinates,
    add_parameters = NULL,
    monitor = monitor
  )
  targets  <- .bt_burnin_requested_states(prepared, monitor, prior_list)
  metadata <- targets$metadata
  selected <- .bt_convergence_select(metadata, monitor, check_indicators)
  rows <- which(metadata$parameter %in% selected & metadata$role != "structural")
  chain_lengths <- vapply(prepared$mcmc_samples_list, NROW, integer(1))
  chain_id <- rep.int(seq_along(chain_lengths), chain_lengths)
  values <- vapply(rows, function(row){
    .bt_burnin_split_Rhat(targets$samples[[row]], chain_id, length(chain_lengths))
  }, numeric(1))
  names(values) <- metadata$parameter[rows]
  values[!is.na(values)]
}

# The convergence targets of a block name an indicator's states by the states
# that the block visits (with two visited states, the occupancy carries the
# indicator's own name). A state-specific 'monitor' element such as
# "mu_indicator (state 3)" names a categorical indicator and one of its
# declared states whatever the block visits: a requested state without its
# own target is added as its occupancy in the block, which is constant (not
# assessable) when the block never visits it. A selector that names no
# categorical indicator or a state it cannot take is left to fail in the
# selection, as in the convergence check.
.bt_burnin_requested_states <- function(prepared, monitor, prior_list){

  targets   <- prepared$targets
  requested <- setdiff(unique(monitor), targets$metadata$parameter)
  if(length(requested) == 0L){
    return(targets)
  }
  draws <- do.call(rbind, prepared$mcmc_samples_list)
  declared <- .bt_burnin_indicator_states(prior_list)
  for(indicator in names(declared)){
    states <- declared[[indicator]]
    labels <- paste0(indicator, " (state ", states, ")")
    if(!any(labels %in% requested) || !indicator %in% colnames(draws) ||
       !any(targets$metadata$source == indicator & targets$metadata$role == "indicator")){
      next
    }
    for(i in which(labels %in% requested)){
      targets$samples[[length(targets$samples) + 1L]] <- as.numeric(draws[, indicator] == states[[i]])
      targets$metadata <- rbind(targets$metadata, data.frame(
        parameter = labels[[i]], source = indicator, role = "indicator", stringsAsFactors = FALSE))
    }
  }
  targets
}

# The declared states of the categorical model indicators, by indicator node:
# the component indicator of a mixture of three or more components (dcat over
# the components; the indicator of a mixture total of an ordered prior
# included). Binary indicators have no state-specific targets.
.bt_burnin_indicator_states <- function(prior_list){

  states <- list()
  add <- function(parameter, prior){
    if(is.prior.factor(prior) && is.prior.ordered(prior)){
      return(add(.prior_ordered_total_name(parameter), prior$total))
    }
    if(!is.prior.mixture(prior) || is.prior.spike_and_slab(prior) || length(prior) < 3L){
      return(NULL)
    }
    indicator <- if(inherits(prior, "prior.bias_mixture")) "bias_indicator" else paste0(parameter, "_indicator")
    states[[indicator]] <<- seq_along(prior)
    NULL
  }
  parameters <- names(prior_list)
  for(i in seq_along(parameters)){
    if(!is.na(parameters[[i]]) && nzchar(parameters[[i]])){
      add(parameters[[i]], prior_list[[i]])
    }
  }
  states
}

# The check and its extensions. 'first' holds the last 'block' draws of the
# requested burn-in, taken every 'thin' iterations; 'run_block(draws)' runs a
# further block of 'draws' draws ('draws * thin' iterations) with the
# monitors on and returns them; 'assess(draws)' returns the split-chain R-hat
# per assessable quantity. Without 'autofit' the chains are not extended.
# With it, blocks of 'burnin_extend' draws are added (each checked alone)
# until the chains are ready, 'max_extensions' blocks were added, a further
# block would exceed 'extension_room' (the burn-in iterations that keep the
# iteration number of the last retained draw within the integer range), or
# 'time_exceeded()' (read before a block is started). Returns the record
# stored in attr(fit, "burnin_check"), which counts its 'block', the 'draws'
# of its blocks and the 'extension' in draws at the fit's thinning, as
# 'fit$burnin' and 'fit$sample' count them.
.bt_burnin_readiness <- function(first, run_block, assess, burnin_max_Rhat, block, burnin_extend, thin,
                                 max_extensions, autofit, time_exceeded, extension_room = Inf){

  blocks      <- list()
  draws       <- first
  block_draws <- as.integer(block)
  extensions  <- 0L
  repeat{
    values <- sort(assess(draws), decreasing = TRUE)
    blocks[[length(blocks) + 1L]] <- data.frame(
      block     = length(blocks),
      draws     = block_draws,
      Rhat      = if(length(values)) unname(values[[1L]]) else NA_real_,
      parameter = if(length(values)) names(values)[[1L]] else NA_character_,
      stringsAsFactors = FALSE
    )
    stop_reason <- if(length(values) == 0L){
      "not_assessable"
    }else if(values[[1L]] <= burnin_max_Rhat){
      "ready"
    }else if(!autofit){
      "no_autofit"
    }else if(extensions >= max_extensions){
      "max_extend"
    }else if((extensions + 1) * burnin_extend * thin > extension_room){
      "iteration_limit"
    }else if(time_exceeded()){
      "max_time"
    }
    if(!is.null(stop_reason)){
      break
    }
    draws       <- run_block(burnin_extend)
    block_draws <- as.integer(burnin_extend)
    extensions  <- extensions + 1L
  }

  failing <- values[values > burnin_max_Rhat]
  list(
    burnin_max_Rhat = burnin_max_Rhat,
    block           = blocks[[1L]]$draws,
    burnin_extend   = burnin_extend,
    thin            = thin,
    blocks          = do.call(rbind, blocks),
    # in the type of 'burnin_extend' (the room check keeps the product within
    # the integer range)
    extension       = extensions * burnin_extend,
    status          = switch(stop_reason, ready = "ready", not_assessable = "not_assessable", "not_ready"),
    stop_reason     = stop_reason,
    failing         = data.frame(parameter = names(failing), Rhat = unname(failing), stringsAsFactors = FALSE)
  )
}

# The warning for chains that were not ready when sampling started, or NULL.
.bt_burnin_warning <- function(record){

  if(is.null(record) || !identical(record$status, "not_ready")){
    return(NULL)
  }
  failing <- record$failing
  listed  <- utils::head(failing, 5L)
  # formatted and compared with a decimal point; the session's decimal mark
  # applies to the displayed text only. Printed numbers are read back with
  # .JAGS_read_number(), which rounds correctly on every platform
  # (as.numeric() can be off in the last digit where long double is double).
  limit   <- .bt_burnin_format_limit(record$burnin_max_Rhat)
  bound   <- max(record$burnin_max_Rhat, .JAGS_read_number(limit, fallback = TRUE))
  shown   <- .bt_burnin_decimal_mark(.bt_burnin_format_Rhat(listed$Rhat, bound))
  limit   <- .bt_burnin_decimal_mark(limit)
  values  <- paste0(shown, " for '", listed$parameter, "'", collapse = ", ")
  if(nrow(failing) > 5L){
    values <- paste0(values, " and ", nrow(failing) - 5L, " more")
  }
  extensions <- nrow(record$blocks) - 1L
  paste0(
    "Chains did not agree at the end of burn-in",
    if(extensions > 0) paste0(" after ", extensions, ngettext(extensions, " extension", " extensions")),
    ": R-hat ", values, " (limit ", limit, "). ",
    switch(record$stop_reason,
      no_autofit = "Increase 'adapt' and 'burnin'; if the chains mix slowly, increase 'thin'. Or fit with 'autofit = TRUE'.",
      max_extend = ,
      max_time   = "Increase 'adapt', 'burnin' or 'thin'.",
      iteration_limit = paste0("A further extension would number the draws beyond the largest ",
        "iteration number (", .Machine$integer.max, "); reduce 'adapt' or the sampling iterations ",
        "('sample' times 'thin') to leave room for it.")
    )
  )
}

# The threshold with two decimals when that is exact (2.00, 2.15), otherwise
# rounded down to three (1.925, 1.828 for 1.8286). The printed number, read
# back as a double, is never above the threshold: a given value just below
# 1.501 prints 1.500.
.bt_burnin_format_limit <- function(threshold){

  printed <- formatC(threshold, format = "f", digits = 2, decimal.mark = ".")
  if(.JAGS_read_number(printed, fallback = TRUE) == threshold){
    return(printed)
  }
  thousandths <- round(threshold * 1000)
  printed <- formatC(thousandths / 1000, format = "f", digits = 3, decimal.mark = ".")
  if(.JAGS_read_number(printed, fallback = TRUE) > threshold){
    printed <- formatC((thousandths - 1) / 1000, format = "f", digits = 3, decimal.mark = ".")
  }
  printed
}

# Values above 'bound' (the threshold, or the printed limit when that reads
# higher) with two decimals, or as many more as it takes to read above it
# (2.00001 above 2); where fifteen are not enough, the shortest
# representation that reads back as the value (2.0000000000000004).
.bt_burnin_format_Rhat <- function(values, bound){

  vapply(values, function(value){
    if(!is.finite(value)){
      return(formatC(value))
    }
    for(digits in 2:15){
      printed <- formatC(value, format = "f", digits = digits, decimal.mark = ".")
      if(.JAGS_read_number(printed, fallback = TRUE) > bound){
        return(printed)
      }
    }
    for(digits in 15:17){
      printed <- formatC(value, format = "g", digits = digits, decimal.mark = ".")
      if(.JAGS_read_number(printed, fallback = TRUE) == value){
        break
      }
    }
    printed
  }, character(1))
}

# Numbers formatted with a decimal point, shown with the session's decimal
# mark (getOption("OutDec"), as formatC() and print() show them).
.bt_burnin_decimal_mark <- function(printed){

  sub(".", getOption("OutDec"), printed, fixed = TRUE)
}

# Signal the warning of chains that were not ready (whatever 'silent' is) and
# return its message for the fit's warnings.
.bt_burnin_warn <- function(record){

  message <- .bt_burnin_warning(record)
  if(!is.null(message)){
    warning(warningCondition(message,
      class = c("BayesTools_burnin_not_converged", "BayesTools_fit_diagnostic")))
  }
  message
}
