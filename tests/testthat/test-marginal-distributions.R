context("Marginal distributions")
set.seed(1)

test_that("helper functions work", {

  # check the posterior distributions with weak priors against a maximum likelihood estimates with ML
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip("under development")

  # complex formula including scaling
  set.seed(1)
  df_all <- data.frame(
    x_cont1  = rnorm(180),
    x_fac2t  = factor(rep(c("A", "B"), 90), levels = c("A", "B")),
    x_fac3md = factor(rep(c("A", "B", "C"), 60), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(180, 0.1, 0.5) + 0.5 + 0.30 * df_all$x_cont1 +
    ifelse(df_all$x_fac3md == "A", 0.3, ifelse(df_all$x_fac3md == "B", -0.3, 0))

  prior_list_0 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior("spike", contrast = "treatment", list = (0)),
    "x_fac3md"         = prior_factor("spike", contrast = "meandif", list(0)),
    "x_cont1:x_fac3md" = prior_factor("spike", contrast = "meandif", list(0))
  )
  prior_list_1 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("normal",  contrast = "treatment", list(0, 0.25)),
    "x_fac3md"         = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)),
    "x_cont1:x_fac3md" = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))
  )
  prior_list <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 1))
  )
  model_syntax <- paste0(
    "model{",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )
  log_posterior <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = parameters[["sigma"]], log = TRUE)))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = list(y = df_all$y, N = nrow(df_all)),
    prior_list = prior_list,
    formula_list       = list(mu = ~ x_cont1),
    formula_prior_list = list(mu = prior_list_0),
    formula_data_list  = list(mu = df_all))
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = list(y = df_all$y, N = nrow(df_all)),
    prior_list = prior_list,
    formula_list       = list(mu = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md),
    formula_prior_list = list(mu = prior_list_1),
    formula_data_list  = list(mu = df_all))
  marglik0 <- JAGS_bridgesampling(
    fit                = fit0,
    log_posterior      = log_posterior,
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list,
    formula_list       = list(mu = ~ x_cont1),
    formula_prior_list = list(mu = prior_list_0),
    formula_data_list  = list(mu = df_all))
  marglik1 <- JAGS_bridgesampling(
    fit                = fit1,
    log_posterior      = log_posterior,
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list,
    formula_list       = list(mu = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md),
    formula_prior_list = list(mu = prior_list_1),
    formula_data_list  = list(mu = df_all))

  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    is_null_list = list(
      "sigma"                = c(FALSE, FALSE),
      "mu_intercept"         = c(FALSE, FALSE),
      "mu_x_cont1"           = c(FALSE, FALSE),
      "x_fac2t"              = c(TRUE, FALSE),
      "x_fac3md"             = c(TRUE, FALSE),
      "x_cont1_xXx_x_fac3md" = c(TRUE, FALSE)
    ),
    conditional  = FALSE)
  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"),
    is_null_list = list(
      "sigma"                = c(FALSE, FALSE),
      "mu_intercept"         = c(FALSE, FALSE),
      "mu_x_cont1"           = c(FALSE, FALSE),
      "x_fac2t"              = c(TRUE, FALSE),
      "x_fac3md"             = c(TRUE, FALSE),
      "x_cont1_xXx_x_fac3md" = c(TRUE, FALSE)
    ),
    seed         = 1,
    conditional  = FALSE
  )


  samples <- mixed_posteriors
  formula <- ~ x_cont1 + x_fac2t + x_cont1*x_fac3md
  formula_parameter <- "mu"
  parameter <- "x_fac2t"
  at <- list(
    x_cont1  = 0,
    x_fac3md = c(NA, "A")
  )
  transformation <- NULL; transformation_arguments <- NULL; transformation_settings <- FALSE

  JAGS_estimates_table(fit1, transform_factors = TRUE)
  ensemble_estimates_table(mixed_posteriors, parameters = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"), transform_factors = TRUE)


})
