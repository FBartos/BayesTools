context("Summary tables functions")

test_that("Summary tables functions work",{

  runjags::runjags.options(silent.jags = T, silent.runjags = T)
  set.seed(1)
  data <- list(
    x = rnorm(20, 0, 1),
    N = 20
  )
  priors_list0 <- list(
    m     =prior("normal", list(0, 1)),
    omega = prior_none()
  )
  priors_list1 <- list(
    m  = prior("normal", list(0, .5)),
    omega = prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  )
  priors_list2 <- list(
    m  = prior("normal", list(0, .3)),
    omega = prior_weightfunction("one.sided", list(c(0.05, 0.50), c(1, 1, 1)))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, 1)
      }
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  fit0 <- JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0)
  fit1 <- JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  fit2 <- JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  marglik0 <- JAGS_bridgesampling(fit0, data, priors_list0, log_posterior)
  marglik1 <- JAGS_bridgesampling(fit1, data, priors_list1, log_posterior)
  marglik2 <- JAGS_bridgesampling(fit2, data, priors_list2, log_posterior)
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_odds = 1, fit_summary = runjags_estimates_table(fit0, priors_list0)),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_odds = 1, fit_summary = runjags_estimates_table(fit1, priors_list1)),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_odds = 1, fit_summary = runjags_estimates_table(fit2, priors_list2))
  )
  models <- models_inference(models)
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1)

  ### checking summary functions
  # model summary
  model_summary <- model_summary_table(models[[2]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  "))
  expect_equal(model_summary[,2], c("2", "0.333", "-0.61", "0.325", "0.964"))
  expect_equal(model_summary[,4], c("Parameter prior distributions", "m ~ Normal(0, 0.5)", "omega[one-sided: .05] ~ CumDirichlet(1, 1)", "", ""))

  # runjags summary
  runjags_summary <- models[[2]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC error", "MCMC SD error", "ESS", "R-hat"))
  expect_equal(rownames(runjags_summary), c("m", "omega[0,0.05]", "omega[0.05,1]"))
  expect_equal(round(unname(unlist(runjags_summary[1,])), 4), round(c(0.155080816, 0.197817354, -0.247495448, 0.167295089, 0.496803251, 0.009208408, 0.047000000, 461.000000000, NA), 4))

  # ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("m", "omega[0,0.05]", "omega[0.05,0.5]", "omega[0.5,1]"))
  expect_equal(round(unname(unlist(estimates_table[1,])), 4), round(c(0.1522389, 0.1519897, -0.2204951, 0.4610624), 4))
  expect_equal(round(unname(unlist(estimates_table[3,])), 4), round(c(0.6794735, 0.7447313,  0.0643561, 1.0000000), 4))

  # ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("Models", "Prior Prob.", "Post. Prob.", "BF10"))
  expect_equal(rownames(inference_table), c("m", "omega"))
  expect_equal(unname(unlist(inference_table[1,])), c(3,   1,   1, Inf))
  expect_equal(round(unname(unlist(inference_table[2,])), 4), round(c(2.0000000, 0.6666667, 0.8001882, 2.0023549), 4))

  # ensemble summary
  summary_table <- ensemble_summary_table(models, c("m", "omega"))
  expect_equal(colnames(summary_table), c("Model", "m", "omega", "Prior prob.", "log(marglik)", "Post. prob.", "Inclusion BF"))
  expect_equal(unname(unlist(summary_table[1,])), c("1", "Normal(0, 1)", "", "0.333333333333333", "-1.10230423916272", "0.199811784123417", "0.49941196373288"))
  expect_equal(unname(unlist(summary_table[2,])), c("2", "Normal(0, 0.5)", "omega[one-sided: .05] ~ CumDirichlet(1, 1)", "0.333333333333333", "-0.614989663490358", "0.32528132384813", "0.96419836991414"))

  # ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "omega"))
  expect_equal(colnames(diagnostics_table), c("Model", "m", "omega", "max(MCMC)", "max(MCMC SD)", "min(ESS)", "max(R-hat)"))
  expect_equal(unname(unlist(diagnostics_table[1,])), c("1", "Normal(0, 1)", "", "0.0101903900487724", "0.048", "434", NA))
  expect_equal(unname(unlist(diagnostics_table[2,])), c("2", "Normal(0, 0.5)", "omega[one-sided: .05] ~ CumDirichlet(1, 1)", "0.0134821149740382", "0.047", "461", NA))

  ### test additional settings
  # transformations
  runjags_summary2t <- runjags_estimates_table(fit2, priors_list2, transformations = list("m" = list(fun = exp)))
  expect_equal(exp(models[[3]]$fit_summary[1,c("Mean","lCI","Median","uCI","MCMC error")]), runjags_summary2t[1,c("Mean","lCI","Median","uCI","MCMC error")], tolerance = 1e-5)
  expect_equal(colnames(models[[3]]$fit_summary), colnames(runjags_summary2t))
  expect_equal(rownames(models[[3]]$fit_summary), rownames(runjags_summary2t))

})
