test_that("JAGS diagnostics work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(150),
    x_fac2t = factor(rep(c("A", "B"), 75), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 50), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(150, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3o == "A", 0.0, ifelse(data_formula$x_fac3o == "B", -0.2, 0.4)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 150
  )


  # create model with mix of a formula and free parameters ---
  formula_list <- list(
    mu    = ~ x_cont1 + x_fac2t + x_fac3o
  )
  formula_data_list <- list(
    mu    = data_formula
  )
  formula_prior_list <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac2t"         = prior_factor("normal", contrast = "treatment", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )
  prior_list <- list(
    sigma = prior("lognormal", list(0, 1)),
    omega = prior_weightfunction("onesided", list(c(0.05, 0.10), c(1,1,1))),
    PET   = prior_PET("gamma", list(2, 2)),
    fac2i = prior_factor("normal", contrast = "independent", list(0, 1/2))
  )
  attr(prior_list$fac2i, "levels") <- 2
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list)


  ### density plots
  expect_doppelganger("diagnostics-plot-density-1", function() JAGS_diagnostics_density(fit, parameter = "mu_x_cont1", formula_prefix = FALSE))
  expect_doppelganger("diagnostics-plot-density-2", function() JAGS_diagnostics_density(fit, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  expect_doppelganger("diagnostics-plot-density-3", function() JAGS_diagnostics_density(fit, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", formula_prefix = FALSE, ylab = "Smth"))
  expect_doppelganger("diagnostics-plot-density-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit, parameter = "mu_x_fac3o")
  })
  expect_doppelganger("diagnostics-plot-density-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_density(fit, parameter = "mu_x_fac3o", formula_prefix = FALSE, transform_orthonormal = TRUE)
  })
  expect_doppelganger("diagnostics-plot-density-6", function()JAGS_diagnostics_density(fit, parameter = "PET"))
  expect_doppelganger("diagnostics-plot-density-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit, parameter = "omega")
  })
  expect_doppelganger("diagnostics-plot-density-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit, parameter = "fac2i")
  })

  expect_doppelganger("diagnostics-ggplot-density-1", JAGS_diagnostics_density(fit, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_density(fit, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_orthonormal = TRUE)
  expect_doppelganger("diagnostics-ggplot-density-2.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-density-2.2",temp_plot[[2]])
  expect_doppelganger("diagnostics-ggplot-density-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_density(fit, plot_type = "ggplot", parameter = "omega")
  expect_doppelganger("diagnostics-ggplot-density-3.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-density-3.2",temp_plot[[2]])


  ### trace plots
  expect_doppelganger("diagnostics-plot-trace-1", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_cont1", formula_prefix = FALSE))
  expect_doppelganger("diagnostics-plot-trace-2", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  expect_doppelganger("diagnostics-plot-trace-3", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", formula_prefix = FALSE, ylab = "Smth"))
  expect_doppelganger("diagnostics-plot-trace-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit, parameter = "mu_x_fac3o", formula_prefix = FALSE)
  })
  expect_doppelganger("diagnostics-plot-trace-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_trace(fit, parameter = "mu_x_fac3o", transform_orthonormal = TRUE)
  })
  expect_doppelganger("diagnostics-plot-trace-6", function() JAGS_diagnostics_trace(fit, parameter = "PET"))
  expect_doppelganger("diagnostics-plot-trace-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit, parameter = "omega")
  })
  expect_doppelganger("diagnostics-plot-trace-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit, parameter = "fac2i")
  })

  expect_doppelganger("diagnostics-ggplot-trace-1", JAGS_diagnostics_trace(fit, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_trace(fit, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_orthonormal = TRUE)
  expect_doppelganger("diagnostics-ggplot-trace-2.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-trace-2.2",temp_plot[[2]])
  expect_doppelganger("diagnostics-ggplot-trace-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_trace(fit, plot_type = "ggplot", parameter = "omega")
  expect_doppelganger("diagnostics-ggplot-trace-3.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-trace-3.2",temp_plot[[2]])


  ### autocorrelation plots
  expect_doppelganger("diagnostics-plot-autocorrelation-1", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_cont1", formula_prefix = FALSE))
  expect_doppelganger("diagnostics-plot-autocorrelation-2", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  expect_doppelganger("diagnostics-plot-autocorrelation-3", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", ylab = "Smth"))
  expect_doppelganger("diagnostics-plot-autocorrelation-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_fac3o")
  })
  expect_doppelganger("diagnostics-plot-autocorrelation-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_fac3o", formula_prefix = FALSE, transform_orthonormal = TRUE)
  })
  expect_doppelganger("diagnostics-plot-autocorrelation-6", function() JAGS_diagnostics_autocorrelation(fit, parameter = "PET"))
  expect_doppelganger("diagnostics-plot-autocorrelation-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit, parameter = "omega")
  })
  expect_doppelganger("diagnostics-plot-autocorrelation-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit, parameter = "fac2i")
  })

  expect_doppelganger("diagnostics-ggplot-autocorrelation-1", JAGS_diagnostics_autocorrelation(fit, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_autocorrelation(fit, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_orthonormal = TRUE)
  expect_doppelganger("diagnostics-ggplot-autocorrelation-2.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-autocorrelation-2.2",temp_plot[[2]])
  expect_doppelganger("diagnostics-ggplot-autocorrelation-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_autocorrelation(fit, plot_type = "ggplot", parameter = "omega")
  expect_doppelganger("diagnostics-ggplot-autocorrelation-3.1",temp_plot[[1]])
  expect_doppelganger("diagnostics-ggplot-autocorrelation-3.2",temp_plot[[2]])
})

test_that("JAGS diagnostics work (spike and slab)", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont  = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, .20 * data_formula$x_cont + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )


  # create model with mix of a formula and free parameters ---
  formula_list <- list(
    mu    = ~ x_cont
  )
  formula_data_list <- list(
    mu    = data_formula
  )
  formula_prior_list <- list(
    mu    = list(
      "intercept"   = prior("normal", list(0, 5)),
      "x_cont"      = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1,1)))
    )
  )
  prior_list <- list(
    sigma = prior("lognormal", list(0, 1))
  )
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list)

  ### density plots
  expect_doppelganger("diagnostics-plot-spike-and-slab-1", function() JAGS_diagnostics_density(fit, parameter = "mu_x_cont"))
  expect_doppelganger("diagnostics-plot-spike-and-slab-2", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_cont"))
  expect_doppelganger("diagnostics-plot-spike-and-slab-3", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_cont"))
})
