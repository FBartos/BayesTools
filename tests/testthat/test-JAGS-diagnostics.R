test_that("JAGS evaluate formula works", {

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

  runjags_estimates_table(fit)

  JAGS_diagnostics_density(fit, parameter = "mu_x_cont1")
  JAGS_diagnostics_density(fit, parameter = "mu_x_cont1", col = c("red", "green", "blue", "yellow"))
  debugonce(JAGS_diagnostics_density)


})
