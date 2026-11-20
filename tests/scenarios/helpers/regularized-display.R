# Fold model inputs into the positive integer key used by scenario_fit().
regularized_fit_version <- function(...) {
  hash <- digest::digest(list(...), algo = "sha256")
  strtoi(substr(hash, 1L, 7L), base = 16L) + 1L
}

# Arrange the panels, title and legend used by the coefficient density plots.
regularized_plot_layout <- function(main) {
  graphics::layout(rbind(c(1, 1, 1), 2:4, 5:7), heights = c(0.30, 1, 1))
  graphics::par(cex = 0.95, mar = rep(0, 4), oma = c(2.8, 1.8, 0.4, 0.5))
  graphics::plot.new()
  graphics::text(0.5, 0.72, main, font = 2, cex = 1.25)
  colours <- c("grey70", "#0072B2", "#D55E00")
  graphics::legend(0.5, 0.25, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
                   legend = c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
                   col = colours, lty = c(1, 1, 2), lwd = 2)

  graphics::par(mar = c(2.5, 2.6, 2.2, 0.6), mgp = c(1.6, 0.45, 0), tcl = -0.25)
}

# The two engines compare the same four priors against ridge in each design.
# For a difference, sum its two estimated log-evidence MCSEs. This conservative
# approximation allows for unknown correlation in fits that reused sampler seeds.
regularized_evidence_table <- function(bridges) {
  if (!identical(names(bridges)[1L], "Ridge")) stop("Use Ridge as the first model.")
  table <- do.call(rbind, lapply(names(bridges), function(name) {
    fit <- bridges[[name]]
    bt <- fit$BayesTools
    rs <- fit$RStan
    if (length(bt$diagnostics$upstream_warnings)) {
      warning(paste(bt$diagnostics$upstream_warnings, collapse = "\n"), call. = FALSE)
    }
    if (!isTRUE(bt$success) || !isFALSE(bt$iteration_limit_exceeded) ||
        rs$niter > 10000L || !all(is.finite(c(bt$logml, rs$logml))) ||
        !isTRUE(is.finite(bt$repetitions$mcse)) || !isTRUE(is.finite(rs$mcse_logml)) ||
        any(!is.finite(rs$q11)) || any(!is.finite(rs$q21))) {
      stop("Unreliable bridge estimate for ", name, ". Inspect the two bridge objects.")
    }
    data.frame(prior = name, engine = c("BayesTools", "RStan"),
               log_ML = c(bt$logml, rs$logml),
               log_ML_MCSE = c(bt$repetitions$mcse, rs$mcse_logml))
  }))
  ridge <- table[table$prior == "Ridge", ]
  reference <- match(table$engine, ridge$engine)
  table$log_BF_vs_ridge <- table$log_ML - ridge$log_ML[reference]
  table$approx_log_BF_SE <- ifelse(table$prior == "Ridge", 0, table$log_ML_MCSE + ridge$log_ML_MCSE[reference])
  numeric_columns <- vapply(table, is.numeric, logical(1))
  table[numeric_columns] <- lapply(table[numeric_columns], round, digits = 3)
  table
}
