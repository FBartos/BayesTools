# Fold model inputs into the positive integer key used by scenario_fit().
regularized_fit_version <- function(...) {
  hash <- digest::digest(list(...), algo = "sha256")
  strtoi(substr(hash, 1L, 7L), base = 16L) + 1L
}

# The likelihood is Gaussian with residual SD 1 in both regressions.
# BayesTools reconstructs mu and adds the fitted priors. RStan evaluates
# its own normalized joint log density, including the constraint Jacobian.
regularized_bridge_estimate <- function(jags, stan, stan_data, seed) {
  likelihood <- function(parameters, data) {
    sum(stats::dnorm(data$y, parameters$mu, sd = 1, log = TRUE))
  }
  bt <- BayesTools::JAGS_bridgesampling(jags, log_posterior = likelihood,
                                       data = list(y = stan_data$y), seed = seed)
  # A reloaded stanfit retains its draws and compiled program but loses its
  # live C++ handle. Construct that handle at the fitted data without sampling.
  target <- rstan::sampling(rstan::get_stanmodel(stan), data = stan_data,
                            chains = 0L, iter = 1L, warmup = 0L, seed = seed + 1L)
  rs <- withr::with_seed(seed + 1L, bridgesampling::bridge_sampler(
    stan, stanfit_model = target, method = "normal", use_neff = TRUE,
    maxiter = 10000, silent = TRUE
  ))
  if (length(bt$diagnostics$upstream_warnings)) {
    warning(paste(bt$diagnostics$upstream_warnings, collapse = "\n"), call. = FALSE)
  }
  if (!isTRUE(bt$success) || !isFALSE(bt$iteration_limit_exceeded) ||
      rs$niter > 10000L || !all(is.finite(c(bt$logml, rs$logml))) ||
      !isTRUE(is.finite(bt$repetitions$mcse)) || !isTRUE(is.finite(rs$mcse_logml))) {
    stop("Bridge sampling did not return two converged estimates with finite MCSEs.")
  }
  if (any(!is.finite(rs$q11)) || any(!is.finite(rs$q21))) {
    stop("The Stan log density was nonfinite at a posterior draw or a bridge proposal.")
  }
  list(BayesTools = bt, RStan = rs)
}

# Refitting either source model invalidates the corresponding evidence cache.
regularized_bridge <- function(name, jags, stan, stan_data) {
  paths <- vapply(paste0(c("jags_", "stan_"), name), function(id) {
    .scenario_fit_paths(id)$fit
  }, character(1))
  if (!all(file.exists(paths))) stop("Missing source fit cache for ", name, ".")
  seed <- regularized_fit_version(name)
  version <- regularized_fit_version(unname(tools::md5sum(paths)), stan_data,
                                     deparse(body(regularized_bridge_estimate)), seed)
  scenario_fit(paste0("bridge_", name), {
    regularized_bridge_estimate(jags, stan, stan_data, seed)
  }, cache_version = version)
}

# The two engines compare the same four priors against ridge in each design.
# For a difference, sum its two estimated log-evidence MCSEs. This conservative
# approximation allows for unknown correlation in fits that reused sampler seeds.
regularized_evidence_table <- function(bridges) {
  if (!identical(names(bridges)[1L], "Ridge")) stop("Use Ridge as the first model.")
  table <- do.call(rbind, lapply(names(bridges), function(name) {
    fit <- bridges[[name]]
    data.frame(prior = name, engine = c("BayesTools", "RStan"),
               log_ML = c(fit$BayesTools$logml, fit$RStan$logml),
               log_ML_MCSE = c(fit$BayesTools$repetitions$mcse, fit$RStan$mcse_logml))
  }))
  ridge <- table[table$prior == "Ridge", ]
  reference <- match(table$engine, ridge$engine)
  table$log_BF_vs_ridge <- table$log_ML - ridge$log_ML[reference]
  table$approx_log_BF_SE <- ifelse(table$prior == "Ridge", 0, table$log_ML_MCSE + ridge$log_ML_MCSE[reference])
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], round, digits = 3)
  table
}
