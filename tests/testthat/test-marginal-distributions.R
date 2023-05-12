context("Marginal distributions")
set.seed(1)

test_that("helper functions work", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip("under development")

  ### complex formula including scaling ----
  set.seed(1)
  df_all <- data.frame(
    x_cont1  = rnorm(180),
    x_fac2t  = factor(rep(c("A", "B"), 90), levels = c("A", "B")),
    x_fac3md = factor(rep(c("A", "B", "C"), 60), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(180, 0.1, 0.5) + 0.5 + 0.20 * df_all$x_cont1 +
    ifelse(df_all$x_fac3md == "A", 0.15, ifelse(df_all$x_fac3md == "B", -0.15, 0))

  prior_list_0 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("spike", contrast = "treatment", list(0)),
    "x_fac3md"         = prior_factor("spike", contrast = "meandif",   list(0)),
    "x_cont1:x_fac3md" = prior_factor("spike", contrast = "meandif",   list(0))
  )
  prior_list_1 <- list(
    "intercept"        = prior("normal", list(0, 1)),
    "x_cont1"          = prior("normal", list(0, 1)),
    "x_fac2t"          = prior_factor("normal",  contrast = "treatment", list(0, 0.25)),
    "x_fac3md"         = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25)),
    "x_cont1:x_fac3md" = prior_factor("mnormal", contrast = "meandif",   list(0, 0.25))
  )
  prior_list <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(prior_list_0$x_cont1, "multiply_by") <- "sigma"
  attr(prior_list_1$x_cont1, "multiply_by") <- "sigma"
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
  model_formula <- list(mu = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md)

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = list(y = df_all$y, N = nrow(df_all)),
    prior_list = prior_list,
    formula_list       = model_formula,
    formula_prior_list = list(mu = prior_list_0),
    formula_data_list  = list(mu = df_all))
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = list(y = df_all$y, N = nrow(df_all)),
    prior_list = prior_list,
    formula_list       = model_formula,
    formula_prior_list = list(mu = prior_list_1),
    formula_data_list  = list(mu = df_all))
  marglik0 <- JAGS_bridgesampling(
    fit                = fit0,
    log_posterior      = log_posterior,
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list,
    formula_list       = model_formula,
    formula_prior_list = list(mu = prior_list_0),
    formula_data_list  = list(mu = df_all))
  marglik1 <- JAGS_bridgesampling(
    fit                = fit1,
    log_posterior      = log_posterior,
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list,
    formula_list       = model_formula,
    formula_prior_list = list(mu = prior_list_1),
    formula_data_list  = list(mu = df_all))

  # make the mixing equal
  marglik1$logml <- marglik0$logml

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

  # manual mixing
  posterior_manual0 <- suppressWarnings(coda::as.mcmc(fit0))
  posterior_manual1 <- suppressWarnings(coda::as.mcmc(fit1))

  ### simple: continuous parameter ----
  marg_post_sigma <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = TRUE)

  hist(marg_post_sigma, freq = FALSE, main = "marginal posterior sigma")
  lines(density(c(posterior_manual0[,"sigma"], posterior_manual1[,"sigma"])))

  hist(attr(marg_post_sigma, "prior_samples"), freq = FALSE, main = "marginal prior sigma", breaks = 20)
  lines(density(prior_list$sigma))


  ### simple: factor ----
  marg_post_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    prior_samples     = TRUE,
    use_formula       = FALSE)

  par(mfrow = c(1, 2))
  hist(marg_post_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")

  hist(marg_post_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_x_fac2t"], posterior_manual1[,"mu_x_fac2t"])))


  par(mfrow = c(1, 2))
  hist(attr(marg_post_x_fac2t[["A"]], "prior_samples"), freq = FALSE, main = "marg_post_x_fac2t = A", breaks = 20)

  hist(attr(marg_post_x_fac2t[["B"]], "prior_samples"), freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  curve(dnorm(x, 0, 0.25)/2, add = T)

  ### formula: intercept ----
  marg_post_int <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_intercept",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  hist(marg_post_int, freq = FALSE, main = "marginal posterior intercept")
  lines(density(c(posterior_manual0[,"mu_intercept"], posterior_manual1[,"mu_intercept"] )))

  hist(attr(marg_post_int, "prior_samples"), freq = FALSE, main = "marginal prior intercept")
  lines(prior_list_0$intercept)

  ### formula: continuous parameter WHAT DO WE EVEN WANT HERE??? ----
  marg_post_x_cont1 <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  hist(marg_post_x_cont1, freq = FALSE, main = "marginal posterior x_cont1")
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0[,"mu_x_cont1"], posterior_manual1[,"mu_intercept"] + posterior_manual1[,"mu_x_cont1"])))

  hist(attr(marg_post_x_cont1, "prior_samples"), freq = FALSE, main = "marginal prior sigma", breaks = 20)


  ### formula: treatment factor ----
  marg_post_x_fac2t <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  par(mfrow = c(1, 2))
  hist(marg_post_x_fac2t[["A"]], freq = FALSE, main = "marg_post_x_fac2t = A")
  lines(density(c(posterior_manual0[,"mu_intercept"], posterior_manual1[,"mu_intercept"])))

  hist(marg_post_x_fac2t[["B"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0[,"mu_x_fac2t"], posterior_manual1[,"mu_intercept"] + posterior_manual1[,"mu_x_fac2t"])))


  ### formula: meandif factor ----
  marg_post_x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md <- posterior_manual0[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.md <- posterior_manual1[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))

  par(mfrow = c(1, 3))
  hist(marg_post_x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1])))

  hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2])))

  hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3])))

  ### formula: meandif factor interaction WHAT DO WE EVEN WANT HERE??? ----
  marg_post_x_cont1__xXx__x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1__xXx__x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md <- posterior_manual0[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.md <- posterior_manual1[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))

  par(mfrow = c(1, 3))
  hist(marg_post_x_cont1__xXx__x_fac3md[["A"]], freq = FALSE, main = "marg_post_x_fac3md = A", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1])))

  hist(marg_post_x_cont1__xXx__x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2])))

  hist(marg_post_x_cont1__xXx__x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3])))

  ### formula: meandif factor + at specification ----
  marg_post_x_fac3md_AT <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_fac3md",
    at                = list(
      x_cont1 = 1,
      x_fac2t = c("A", "B")
    ),
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  posterior_manual0.md  <- posterior_manual0[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.md  <- posterior_manual1[,c("mu_x_fac3md[1]", "mu_x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual0.mdi <- posterior_manual0[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))
  posterior_manual1.mdi <- posterior_manual1[,c("mu_x_cont1__xXx__x_fac3md[1]", "mu_x_cont1__xXx__x_fac3md[2]")] %*% t(contr.meandif(1:3))


  par(mfrow = c(2, 3))
  hist(marg_post_x_fac3md_AT[["A"]][1,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,A", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1] + posterior_manual0[,"mu_x_cont1"] + posterior_manual0.mdi[,1],
                  posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1] + posterior_manual1[,"mu_x_cont1"] + posterior_manual1.mdi[,1])))

  hist(marg_post_x_fac3md_AT[["A"]][2,], freq = FALSE, main = "marg_post_x_fac3md = A | 1,B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,1] + posterior_manual0[,"mu_x_cont1"] + posterior_manual0[,"mu_x_fac2t"] + posterior_manual0.mdi[,2],
                  posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,1] + posterior_manual1[,"mu_x_cont1"] + posterior_manual1[,"mu_x_fac2t"] + posterior_manual1.mdi[,2])))


  hist(marg_post_x_fac3md[["B"]], freq = FALSE, main = "marg_post_x_fac3md = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,2], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,2])))

  hist(marg_post_x_fac3md[["C"]], freq = FALSE, main = "marg_post_x_fac2t = B", breaks = 20)
  lines(density(c(posterior_manual0[,"mu_intercept"] + posterior_manual0.md[,3], posterior_manual1[,"mu_intercept"] + posterior_manual1.md[,3])))




  post_x_cont1_x_fac3md <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "mu_x_cont1__xXx__x_fac3md",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    prior_samples     = TRUE)

  class(samples[["mu_x_fac2t"]]) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")
  post_x_fac2t <- marginal_posterior(
    samples           = samples,
    parameter         = "mu_x_fac2t",
    formula           = ~ x_cont1 + x_fac2t + x_cont1*x_fac3md,
    at                = list(
      x_cont1  = 0,
      x_fac3md = c(NA, "A")
    ),
    prior_samples     = TRUE)

  class(samples[["mu_x_cont1__xXx__x_fac3md"]]) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")
  post_mu_x_cont1__xXx__x_fac3md <- marginal_posterior(
    samples           = samples,
    parameter         = "mu_x_cont1__xXx__x_fac3md",
    prior_samples     = TRUE)

  post_sigma <- marginal_posterior(
    samples           = mixed_posteriors,
    parameter         = "sigma",
    prior_samples     = TRUE)


  hist(post_sigma)
  hist(attr(post_sigma, "prior_samples"))

  JAGS_estimates_table(fit1, transform_factors = TRUE)
  ensemble_estimates_table(mixed_posteriors, parameters = c("sigma", "mu_intercept", "mu_x_cont1", "mu_x_fac2t", "mu_x_fac3md", "mu_x_cont1__xXx__x_fac3md"), transform_factors = TRUE)


})

