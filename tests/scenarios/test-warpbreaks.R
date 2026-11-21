if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("warpbreaks")

# Do wool type and loom tension explain the number of yarn breaks? Compare
# four Gaussian models for log(breaks / 25), with an intercept in every
# model and independent 50% inclusion probabilities for wool and tension.
# Wool is coded as B versus A. Tension uses orthonormal sum-to-zero
# contrasts, so its level rows describe deviations from the average level.
testthat::test_that("Wool and loom tension model probabilities", {

  data <- datasets::warpbreaks
  y    <- log(data[["breaks"]] / 25)
  tension <- stats::contr.helmert(levels(data[["tension"]]))
  tension <- sweep(tension, 2L, sqrt(colSums(tension^2)), "/")
  X <- stats::model.matrix(~ wool + tension, data,
                          contrasts.arg = list(wool = "contr.treatment", tension = tension))

  models <- expand.grid(wool = c(FALSE, TRUE), tension = c(FALSE, TRUE))
  model_names <- c("Neither", "Wool", "Tension", "Wool and tension")

  # All active coefficients have independent N(0, sigma^2) priors given the
  # residual variance. A single indicator switches both tension coefficients.
  fit <- scenario_fit("product-space", {
    JAGS_fit(
      "model {
        wool_on ~ dbern(0.5)
        tension_on ~ dbern(0.5)
        active[1] <- 1
        active[2] <- wool_on
        active[3] <- tension_on
        active[4] <- tension_on
        for (j in 1:4) {
          beta[j] ~ dnorm(0, precision)
          effect[j] <- beta[j] * active[j]
        }
        for (i in 1:N) {
          y[i] ~ dnorm(inprod(X[i, ], effect), precision)
        }
      }",
      data = list(y = y, X = X, N = length(y)),
      prior_list = list(precision = prior("gamma", list(shape = 2, rate = 0.5))),
      add_parameters = c("wool_on", "tension_on", "effect"),
      chains = 4, adapt = 1000, burnin = 1000, sample = 20000, seed = 1
    )
  })

  # Load the MCMC methods when the fit is read from cache.
  loadNamespace("runjags")
  draws <- as.matrix(coda::as.mcmc.list(fit))
  coefficient_draws <- draws[, paste0("effect[", 2:4, "]")]
  model <- 1L + draws[, "wool_on"] + 2L * draws[, "tension_on"]
  fit_prob <- stats::setNames(tabulate(model, nbins = nrow(models)) / length(model), model_names)
  fit_mean <- rowsum(coefficient_draws, factor(model, levels = 1:4),
                     reorder = TRUE) / tabulate(model, nbins = nrow(models))

  # Independent conjugate reference. Integrating the coefficients and the
  # Gamma(2, 0.5) residual precision gives each model's marginal likelihood
  # up to a common constant, and the posterior coefficient means.
  exact <- lapply(seq_len(nrow(models)), function(i) {
    active <- c(TRUE, models[["wool"]][i], rep(models[["tension"]][i], 2L))
    design <- X[, active, drop = FALSE]
    precision <- crossprod(design) + diag(ncol(design))
    score <- crossprod(design, y)
    mean <- solve(precision, score)
    rate <- 0.5 + (sum(y^2) - drop(crossprod(mean, score))) / 2
    coefficients <- numeric(ncol(X))
    coefficients[active] <- mean
    list(log_ml = -sum(log(diag(chol(precision)))) - (2 + length(y) / 2) * log(rate),
         mean = coefficients[-1L])
  })
  exact_log_ml <- vapply(exact, function(x) x[["log_ml"]], numeric(1))
  exact_weight <- exp(exact_log_ml - max(exact_log_ml))
  exact_prob <- exact_weight / sum(exact_weight)
  exact_mean <- do.call(rbind, lapply(exact, function(x) x[["mean"]]))

  # The level contrasts use the same design as the fit and reference. Order
  # the table by posterior probability, matching the plot columns.
  fit_effects   <- cbind(fit_mean[, 1L], fit_mean[, 2:3] %*% t(tension))
  exact_effects <- cbind(exact_mean[, 1L], exact_mean[, 2:3] %*% t(tension))
  display_order <- order(fit_prob, decreasing = TRUE)
  scenario_text("model-probabilities", {
    data.frame(model = model_names[display_order],
               BayesTools = round(fit_prob[display_order], 3),
               exact = round(exact_prob[display_order], 3), row.names = NULL)
  })
  scenario_text("conditional-contrasts", {
    colnames(fit_effects) <- colnames(exact_effects) <- c("Wool B - A", paste("Tension", levels(data[["tension"]])))
    data.frame(
      source = rep(c("BayesTools", "Exact"), each = nrow(models)),
      model = rep(model_names[display_order], 2L),
      round(rbind(fit_effects[display_order, ], exact_effects[display_order, ]), 3),
      row.names = NULL, check.names = FALSE
    )
  })

  # Apply the fitted level design to every draw before taking within-model medians.
  contrast_draws <- cbind(coefficient_draws[, 1L], coefficient_draws[, 2:3] %*% t(tension))
  fit_medians <- t(vapply(seq_len(nrow(models)), function(i) {
    apply(contrast_draws[model == i, , drop = FALSE], 2L, stats::median)
  }, numeric(ncol(contrast_draws))))

  # Make the whole-factor view and the signed level view of the same four
  # models. A blank segment means the corresponding effect is switched off.
  row <- c("wool", "tension", paste0("tension_", levels(data[["tension"]])))
  label <- c("Wool B \u2212 A", "Tension", paste("Tension", levels(data[["tension"]])))
  group <- c("Material", rep("Loom setting", 4L))
  signed <- ifelse(fit_medians > 0, "positive", ifelse(fit_medians < 0, "negative", "included"))
  signed[!models[["wool"]], 1L] <- "excluded"
  signed[!models[["tension"]], 2:4] <- "excluded"
  states <- cbind(signed[, 1L], ifelse(models[["tension"]], "included", "excluded"), signed[, 2:4])
  cells <- data.frame(model = rep(model_names, each = length(row)),
                      row = rep(row, times = nrow(models)),
                      label = rep(label, times = nrow(models)),
                      group = rep(group, times = nrow(models)),
                      state = as.vector(t(states)))

  scenario_plot("term-inclusion", {
    plot_model_probabilities(fit_prob, cells[cells[["row"]] %in% c("wool", "tension"), ], legend = TRUE)
  })
  scenario_plot("signed-tension-levels", {
    plot_model_probabilities(fit_prob, cells[cells[["row"]] != "tension", ],
                             group_labels = TRUE, legend = TRUE, plot_type = "ggplot")
  })
})
