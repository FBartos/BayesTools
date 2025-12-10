context("JAGS diagnostics")

# Load common test helpers
source(testthat::test_path("common-functions.R"))

test_that("JAGS diagnostics work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Load pre-fitted models
  fit_formula_mix <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix.RDS"))
  fit_formula_fac <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_fac.RDS"))
  fit_pet         <- readRDS(file.path(temp_fits_dir, "fit_pet.RDS"))
  fit_wf          <- readRDS(file.path(temp_fits_dir, "fit_wf_onesided.RDS"))
  fit_independent <- readRDS(file.path(temp_fits_dir, "fit_factor_independent.RDS"))
  fit_meandif     <- readRDS(file.path(temp_fits_dir, "fit_spike_factors_alt.RDS"))

  ### density plots
  vdiffr::expect_doppelganger("diagnostics-plot-density-1", function() JAGS_diagnostics_density(fit_formula_mix, parameter = "mu_x_cont1", formula_prefix = FALSE))
  vdiffr::expect_doppelganger("diagnostics-plot-density-2", function() JAGS_diagnostics_density(fit_formula_mix, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  vdiffr::expect_doppelganger("diagnostics-plot-density-3", function() JAGS_diagnostics_density(fit_formula_fac, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", formula_prefix = FALSE, ylab = "Smth"))
  vdiffr::expect_doppelganger("diagnostics-plot-density-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit_formula_mix, parameter = "mu_x_fac3o")
  })
  vdiffr::expect_doppelganger("diagnostics-plot-density-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_density(fit_formula_mix, parameter = "mu_x_fac3o", formula_prefix = FALSE, transform_factors = TRUE)
  })
  vdiffr::expect_doppelganger("diagnostics-plot-density-6", function() JAGS_diagnostics_density(fit_pet, parameter = "PET"))
  vdiffr::expect_doppelganger("diagnostics-plot-density-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit_wf, parameter = "omega")
  })
  vdiffr::expect_doppelganger("diagnostics-plot-density-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    # Using p1 from fit_independent (prior only model)
    JAGS_diagnostics_density(fit_independent, parameter = "p1")
  })

  vdiffr::expect_doppelganger("diagnostics-ggplot-density-1", JAGS_diagnostics_density(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_density(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_factors = TRUE)
  vdiffr::expect_doppelganger("diagnostics-ggplot-density-2.1",temp_plot[[1]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-density-2.2",temp_plot[[2]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-density-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_density(fit_wf, plot_type = "ggplot", parameter = "omega")
  vdiffr::expect_doppelganger("diagnostics-ggplot-density-3.1",temp_plot)


  ### trace plots
  vdiffr::expect_doppelganger("diagnostics-plot-trace-1", function() JAGS_diagnostics_trace(fit_formula_mix, parameter = "mu_x_cont1", formula_prefix = FALSE))
  vdiffr::expect_doppelganger("diagnostics-plot-trace-2", function() JAGS_diagnostics_trace(fit_formula_mix, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  vdiffr::expect_doppelganger("diagnostics-plot-trace-3", function() JAGS_diagnostics_trace(fit_formula_fac, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", formula_prefix = FALSE, ylab = "Smth"))
  vdiffr::expect_doppelganger("diagnostics-plot-trace-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit_formula_mix, parameter = "mu_x_fac3o", formula_prefix = FALSE)
  })
  vdiffr::expect_doppelganger("diagnostics-plot-trace-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_trace(fit_formula_mix, parameter = "mu_x_fac3o", transform_factors = TRUE)
  })
  vdiffr::expect_doppelganger("diagnostics-plot-trace-6", function() JAGS_diagnostics_trace(fit_pet, parameter = "PET"))
  vdiffr::expect_doppelganger("diagnostics-plot-trace-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit_wf, parameter = "omega")
  })
  vdiffr::expect_doppelganger("diagnostics-plot-trace-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit_independent, parameter = "p1")
  })

  vdiffr::expect_doppelganger("diagnostics-ggplot-trace-1", JAGS_diagnostics_trace(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_trace(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_factors = TRUE)
  vdiffr::expect_doppelganger("diagnostics-ggplot-trace-2.1",temp_plot[[1]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-trace-2.2",temp_plot[[2]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-trace-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_trace(fit_wf, plot_type = "ggplot", parameter = "omega")
  vdiffr::expect_doppelganger("diagnostics-ggplot-trace-3.1",temp_plot)


  ### autocorrelation plots
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-1", function() JAGS_diagnostics_autocorrelation(fit_formula_mix, parameter = "mu_x_cont1", formula_prefix = FALSE))
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-2", function() JAGS_diagnostics_autocorrelation(fit_formula_mix, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-3", function() JAGS_diagnostics_autocorrelation(fit_formula_fac, parameter = "mu_x_fac2t", main = "Treatment", xlab = "Values", ylab = "Smth"))
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-4", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit_formula_mix, parameter = "mu_x_fac3o")
  })
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-5", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_autocorrelation(fit_formula_mix, parameter = "mu_x_fac3o", formula_prefix = FALSE, transform_factors = TRUE)
  })
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-6", function() JAGS_diagnostics_autocorrelation(fit_pet, parameter = "PET"))
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit_wf, parameter = "omega")
  })
  vdiffr::expect_doppelganger("diagnostics-plot-autocorrelation-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit_independent, parameter = "p1")
  })

  vdiffr::expect_doppelganger("diagnostics-ggplot-autocorrelation-1", JAGS_diagnostics_autocorrelation(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"), formula_prefix = FALSE, transformations = list(mu_x_cont1 = list(fun = function(x) exp(x)))))
  temp_plot <- JAGS_diagnostics_autocorrelation(fit_formula_mix, plot_type = "ggplot", parameter = "mu_x_fac3o", transform_factors = TRUE)
  vdiffr::expect_doppelganger("diagnostics-ggplot-autocorrelation-2.1",temp_plot[[1]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-autocorrelation-2.2",temp_plot[[2]])
  vdiffr::expect_doppelganger("diagnostics-ggplot-autocorrelation-2.3",temp_plot[[3]])
  temp_plot <- JAGS_diagnostics_autocorrelation(fit_wf, plot_type = "ggplot", parameter = "omega")
  vdiffr::expect_doppelganger("diagnostics-ggplot-autocorrelation-3.1",temp_plot)
})

test_that("JAGS diagnostics work (spike and slab)", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Use fit_complex_mixed which has spike and slab on x_cont1
  fit <- readRDS(file.path(temp_fits_dir, "fit_complex_mixed.RDS"))

  ### density plots
  vdiffr::expect_doppelganger("diagnostics-plot-spike-and-slab-1", function() JAGS_diagnostics_density(fit, parameter = "mu_x_cont1"))
  vdiffr::expect_doppelganger("diagnostics-plot-spike-and-slab-2", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_cont1"))
  vdiffr::expect_doppelganger("diagnostics-plot-spike-and-slab-3", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_cont1"))
})

test_that("JAGS diagnostics work (mixture priors)", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Use fit_complex_mixed which has mixture on intercept and x_fac3t
  fit <- readRDS(file.path(temp_fits_dir, "fit_complex_mixed.RDS"))

  ### density plots
  # Using mu_intercept as the first mixture example (was mu_x_cont in original, but x_cont1 is spike/slab in this model)
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-1", function() JAGS_diagnostics_density(fit, parameter = "mu_intercept"))
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-2", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_intercept"))
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-3", function() JAGS_diagnostics_trace(fit, parameter = "mu_intercept"))

  # Using mu_x_fac3t as the second mixture example
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-4", function() JAGS_diagnostics_density(fit, parameter = "mu_x_fac3t"))
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-5", function() JAGS_diagnostics_autocorrelation(fit, parameter = "mu_x_fac3t"))
  vdiffr::expect_doppelganger("diagnostics-plot-mixture-6", function() JAGS_diagnostics_trace(fit, parameter = "mu_x_fac3t"))
})

test_that("JAGS diagnostics work (meandif and independent)", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  fit_independent <- readRDS(file.path(temp_fits_dir, "fit_factor_independent.RDS"))
  fit_meandif     <- readRDS(file.path(temp_fits_dir, "fit_spike_factors_alt.RDS"))

  ### density plots
  vdiffr::expect_doppelganger("diagnostics3-plot-density-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit_independent, parameter = "p1")
  })
  vdiffr::expect_doppelganger("diagnostics3-plot-density-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_density(fit_meandif, parameter = "mu_x_fac3md")
  })
  vdiffr::expect_doppelganger("diagnostics3-plot-density-3", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_density(fit_meandif, parameter = "mu_x_fac3md", formula_prefix = FALSE, transform_factors = TRUE)
  })


  temp_plot <- JAGS_diagnostics_density(fit_independent, plot_type = "ggplot", parameter = "p1")
  vdiffr::expect_doppelganger("diagnostics3-ggplot-density-1.1",temp_plot[[1]])
  vdiffr::expect_doppelganger("diagnostics3-ggplot-density-1.2",temp_plot[[2]])

  temp_plot <- JAGS_diagnostics_density(fit_meandif, plot_type = "ggplot", parameter = "mu_x_fac3md", transform_factors = TRUE)
  vdiffr::expect_doppelganger("diagnostics3-ggplot-density-2.1",temp_plot[[1]])
  vdiffr::expect_doppelganger("diagnostics3-ggplot-density-2.2",temp_plot[[2]])
  vdiffr::expect_doppelganger("diagnostics3-ggplot-density-2.3",temp_plot[[3]])


  ### trace plots
  vdiffr::expect_doppelganger("diagnostics3-plot-trace-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_trace(fit_independent, parameter = "p1")
  })
  vdiffr::expect_doppelganger("diagnostics3-plot-trace-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_trace(fit_meandif, parameter = "mu_x_fac3md", formula_prefix = FALSE, transform_factors = TRUE)
  })


  ### autocorrelation plots
  vdiffr::expect_doppelganger("diagnostics3-plot-autocorrelation-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    JAGS_diagnostics_autocorrelation(fit_independent, parameter = "p1")
  })
  vdiffr::expect_doppelganger("diagnostics3-plot-autocorrelation-2", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 3))
    JAGS_diagnostics_autocorrelation(fit_meandif, parameter = "mu_x_fac3md", formula_prefix = FALSE, transform_factors = TRUE)
  })
})

test_that("JAGS diagnostics work (spike priors)", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  fit <- readRDS(file.path(temp_fits_dir, "fit_spike_factors.RDS"))
  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_simple.RDS"))

  ### density plots
  vdiffr::expect_doppelganger("diagnostics4-ggplot-density-fit_simple",JAGS_diagnostics_density(fit_simple, parameter = "mu"))

  # fit_spike_factors has factor spikes
  expect_message(JAGS_diagnostics_density(fit, parameter = "mu_x_fac2i"),         "No diagnostic plots are produced for a spike prior distribution")
  expect_message(JAGS_diagnostics_trace(fit, parameter = "mu_x_fac3md"),          "No diagnostic plots are produced for a spike prior distribution")

})
