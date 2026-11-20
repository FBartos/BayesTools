# The catalog puts BayesTools draws and priors on the selected effect's scale.
# Supply the corresponding RStan draws in the same units and order.
scenario_shrinkage_plot <- function(fit, posteriors, selectors, main) {

  catalog <- parameter_catalog(fit)
  effects <- names(selectors)
  mixed <- lapply(selectors, function(selector) {
    parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, selector))
  })
  xlim <- c(-1.5, 1.5)
  n_points <- 301
  stan_curves <- lapply(effects, function(effect) {
    stats::density(posteriors$RStan[, effect], n = n_points)
  })
  # The native plot evaluates the prior on a regular grid. An arrow continues
  # its last finite ordinate towards the unbounded density at zero.
  grid_step <- diff(xlim) / (n_points - 1)
  priors <- lapply(mixed, posterior_metadata, field = "prior_density")
  unbounded <- vapply(priors, function(p) prior_density_ordinate(p, 0)$behavior == "infinite", logical(1))
  prior_height <- vapply(seq_along(priors), function(i) {
    x <- if (unbounded[i]) grid_step else 0
    exp(prior_density_ordinate(priors[[i]], x)$log_density)
  }, numeric(1))
  peak <- function(curve) max(curve$y)
  ymax <- ceiling(1.12 * max(prior_height, vapply(stan_curves, peak, numeric(1)),
                             vapply(mixed, function(q) peak(stats::density(q)), numeric(1))))

  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::layout(rbind(c(1, 1, 1), 2:4, 5:7), heights = c(0.30, 1, 1))
  graphics::par(cex = 0.95, mar = rep(0, 4), oma = c(2.8, 1.8, 0.4, 0.5))
  graphics::plot.new()
  graphics::text(0.5, 0.72, main, font = 2, cex = 1.25)
  colours <- c("grey70", "#0072B2", "#D55E00")
  graphics::legend(0.5, 0.25, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
                   legend = c("BayesTools prior", "BayesTools posterior", "RStan posterior"),
                   col = colours, lty = c(1, 1, 2), lwd = 2)

  graphics::par(mar = c(2.5, 2.6, 2.2, 0.6), mgp = c(1.6, 0.45, 0), tcl = -0.25)
  for (i in seq_along(effects)) {
    effect <- effects[i]
    graphics::plot(priors[[i]], main = effect, xlab = "", ylab = "",
                   xlim = xlim, ylim = c(0, ymax), n_points = n_points, col = colours[1], lwd = 2)
    if (unbounded[i]) {
      graphics::arrows(0, prior_height[i], 0, ymax, col = colours[1], lwd = 1.5, length = 0.065)
    }
    plot_posterior(setNames(mixed[i], effect), effect, add = TRUE,
                   xlim = xlim, n_points = n_points, col = colours[2], lwd = 2)
    graphics::lines(stan_curves[[i]], col = colours[3], lty = 2, lwd = 2)
  }
  graphics::mtext("Density", side = 2, outer = TRUE, line = 0.3)
  graphics::mtext("Standardized effect", side = 1, outer = TRUE, line = 0)
  if (any(unbounded)) {
    graphics::mtext("Arrows mark unbounded prior density at zero.", side = 1, outer = TRUE, line = 1.4, cex = 0.8)
  }
  invisible(NULL)
}
