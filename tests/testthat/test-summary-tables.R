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
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0)),
    list(fit = fit1, marglik = marglik1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1)),
    list(fit = fit2, marglik = marglik2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2))
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
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("m", "omega[0,0.05]", "omega[0.05,1]"))
  expect_equal(unname(unlist(runjags_summary[1,])), c(0.155080816, 0.197817354, -0.247495448, 0.167295089, 0.496803251, 0.009208408, 0.047000000, 461.000000000, NA), tolerance = 1e-4)

  # ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("m", "omega[0,0.05]", "omega[0.05,0.5]", "omega[0.5,1]"))
  expect_equal(unname(unlist(estimates_table[1,])), c(0.1522389, 0.1519897, -0.2204951, 0.4610624), tolerance = 1e-4)
  expect_equal(unname(unlist(estimates_table[3,])), c(0.6794735, 0.7447313,  0.0643561, 1.0000000), tolerance = 1e-4)

  # ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(rownames(inference_table), c("m", "omega"))
  expect_equal(unname(unlist(inference_table[1,])), c(3,   1,   1, Inf))
  expect_equal(unname(unlist(inference_table[2,])), c(2.0000000, 0.6666667, 0.8001882, 2.0023549), tolerance = 1e-4)

  # ensemble summary
  summary_table <- ensemble_summary_table(models, c("m", "omega"))
  expect_equal(colnames(summary_table), c("Model", "m", "omega", "prior_prob", "marglik", "post_prob", "inclusion_BF"))
  expect_equal(unname(as.vector(summary_table[,1])), c(1, 2, 3))
  expect_equal(unname(as.vector(summary_table[,2])), c("Normal(0, 1)", "Normal(0, 0.5)", "Normal(0, 0.3)"))
  expect_equal(unname(as.vector(summary_table[,3])), c("", "omega[one-sided: .05] ~ CumDirichlet(1, 1)", "omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)"))
  expect_equal(unname(as.vector(summary_table[,4])), c(0.3333333, 0.3333333, 0.3333333),    tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,5])), c(-1.1023042, -0.6149897, -0.2365613), tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,6])), c(0.1998118, 0.3252813, 0.4749069),    tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,7])), c(0.4994120, 0.9641984, 1.8088483),    tolerance = 1e-4)

  # ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "omega"))
  expect_equal(colnames(diagnostics_table), c("Model", "m", "omega", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))

  expect_equal(unname(as.vector(diagnostics_table[,1])), c(1, 2, 3))
  expect_equal(unname(as.vector(diagnostics_table[,2])), c("Normal(0, 1)", "Normal(0, 0.5)", "Normal(0, 0.3)"))
  expect_equal(unname(as.vector(diagnostics_table[,3])), c("", "omega[one-sided: .05] ~ CumDirichlet(1, 1)", "omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)"))
  expect_equal(unname(as.vector(diagnostics_table[,4])), c(0.01019039, 0.01348211, 0.01061287), tolerance = 1e-4)
  expect_equal(unname(as.vector(diagnostics_table[,5])), c(0.048, 0.047, 0.045), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,6])), c(434, 461, 500))
  expect_equal(unname(as.vector(diagnostics_table[,7])), c(NA, NA, NA))


  ### test additional settings
  # transformations
  runjags_summary2t <- runjags_estimates_table(fit2, transformations = list("m" = list(fun = exp)))
  expect_equal(exp(models[[3]]$fit_summary[1,c("lCI","Median","uCI","MCMC_error")]), runjags_summary2t[1,c("lCI","Median","uCI","MCMC_error")], tolerance = 1e-5)
  expect_equal(colnames(models[[3]]$fit_summary), colnames(runjags_summary2t))
  expect_equal(rownames(models[[3]]$fit_summary), rownames(runjags_summary2t))

  ### test an empty tables
  runjags_summary_empty <- runjags_estimates_empty_table()
  expect_equivalent(nrow(runjags_summary_empty), 0)
  expect_equal(colnames(runjags_summary_empty), colnames(runjags_summary))
  expect_equal(capture_output_lines(runjags_summary_empty, width = 150)[1], capture_output_lines(runjags_summary, width = 150)[1])

  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  expect_equivalent(nrow(ensemble_estimates_empty), 0)
  expect_equal(colnames(ensemble_estimates_empty), colnames(estimates_table))
  expect_equal(capture_output_lines(ensemble_estimates_empty, width = 150)[1], capture_output_lines(estimates_table, width = 150)[1])

  ensemble_inference_empty <- ensemble_inference_empty_table()
  expect_equivalent(nrow(ensemble_inference_empty), 0)
  expect_equal(colnames(ensemble_inference_empty), colnames(inference_table))
  expect_equal(capture_output_lines(ensemble_inference_empty, width = 150)[1], capture_output_lines(inference_table, width = 150)[1])

  ensemble_summary_table <- ensemble_summary_empty_table()
  expect_equivalent(nrow(ensemble_summary_table), 0)
  summary_table.trimmed <- remove_column(summary_table, 2)
  summary_table.trimmed <- remove_column(summary_table.trimmed, 2)
  expect_equal(colnames(ensemble_summary_table), colnames(summary_table.trimmed))
  expect_equal(capture_output_lines(ensemble_summary_table, width = 150)[1], capture_output_lines(summary_table.trimmed, width = 150)[1])

  ensemble_diagnostics_empty <- ensemble_diagnostics_empty_table()
  expect_equivalent(nrow(ensemble_diagnostics_empty), 0)
  diagnostics_table.trimmed <- remove_column(diagnostics_table, 2)
  diagnostics_table.trimmed <- remove_column(diagnostics_table.trimmed, 2)
  expect_equal(colnames(ensemble_diagnostics_empty), colnames(diagnostics_table.trimmed))
  expect_equal(capture_output_lines(ensemble_diagnostics_empty, width = 150)[1], capture_output_lines(diagnostics_table.trimmed, width = 150)[1])

  model_summary_empty <- model_summary_empty_table()
  expect_equivalent(nrow(model_summary_empty), 5)
  expect_equal(model_summary_empty[,1], model_summary[,1])
  expect_equal(model_summary_empty[1,4], model_summary[1,4])

  ### test print functions
  expect_equal(capture_output_lines(model_summary, print = TRUE, width = 150),
               c("                                                                            ",
                 " Model              2                          Parameter prior distributions",
                 " Prior prob.    0.333                                 m ~ Normal(0, 0.5)    ",
                 " log(marglik)   -0.61             omega[one-sided: .05] ~ CumDirichlet(1, 1)",
                 " Post. prob.    0.325                                                       ",
                 " Inclusion BF   0.964                                                       "
               ))
  expect_equal(capture_output_lines(runjags_summary,   print = TRUE, width = 150),
               c("               Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat",
                 "m             0.155 0.198 -0.247  0.167 0.497     0.00921          0.047 461    NA",
                 "omega[0,0.05] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA",
                 "omega[0.05,1] 0.509 0.301  0.028  0.508 0.983     0.01348          0.045 500    NA"

               ))
  expect_equal(capture_output_lines(estimates_table,   print = TRUE, width = 150),
               c("                 Mean Median  0.025  0.95",
                 "m               0.152  0.152 -0.220 0.461",
                 "omega[0,0.05]   1.000  1.000  1.000 1.000",
                 "omega[0.05,0.5] 0.679  0.745  0.064 1.000",
                 "omega[0.5,1]    0.529  0.483  0.023 1.000"

               ))
  expect_equal(capture_output_lines(inference_table,   print = TRUE, width = 150),
               c("      Models Prior prob. Post. prob. Inclusion BF",
                 "m        3/3       1.000       1.000          Inf",
                 "omega    2/3       0.667       0.800        2.002"

               ))
  expect_equal(capture_output_lines(summary_table,     print = TRUE, width = 150),
               c(" Model     Prior m                        Prior omega                    Prior prob. log(marglik) Post. prob. Inclusion BF",
                 "     1    Normal(0, 1)                                                         0.333        -1.10       0.200        0.499",
                 "     2  Normal(0, 0.5)     omega[one-sided: .05] ~ CumDirichlet(1, 1)          0.333        -0.61       0.325        0.964",
                 "     3  Normal(0, 0.3) omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)       0.333        -0.24       0.475        1.809"
               ))
  expect_equal(capture_output_lines(diagnostics_table, print = TRUE, width = 150),
               c(" Model     Prior m                        Prior omega                    max[error(MCMC)] max[error(MCMC)/SD] min(ESS) max(R-hat)",
                 "     1    Normal(0, 1)                                                            0.01019               0.048      434         NA",
                 "     2  Normal(0, 0.5)     omega[one-sided: .05] ~ CumDirichlet(1, 1)             0.01348               0.047      461         NA",
                 "     3  Normal(0, 0.3) omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)          0.01061               0.045      500         NA"
               ))


  ### test adding columns
  expect_error(add_column(runjags_summary, column_title = "New Title", column_values = c(0.2, 0.3, 0.4, 0.5)),
               "The 'column_values' must be a vector of the same length as has the table rows.")
  expect_error(add_column(runjags_summary, column_title = "New Title", column_values = c(0.2, 0.3, 0.4), column_type = "random text"),
               "The 'random text' values are not recognized by the 'column_type' argument.")
  expect_error(add_column(runjags_summary, column_title = "New Title", column_values = c(0.2, 0.3, 0.4), column_position = 55),
               "The 'column_position' must be equal or lower than ")
  expect_error(add_column(data.frame(a = 1:3, b = c("A", "B", "C")), column_title = "New Title", column_values = c(0.2, 0.3, 0.4)),
               "The 'table' must be of class 'BayesTools_table'.")

  expect_equal(capture_output_lines(
    add_column(runjags_summary, column_title = "New Title", column_values = c(0.2, 0.3, 0.4)), print = TRUE, width = 150),
               c("               Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD ESS R-hat New Title",
                 "m             0.155 0.198 -0.247  0.167 0.497     0.00921          0.047 461    NA     0.200",
                 "omega[0,0.05] 1.000 0.000  1.000  1.000 1.000          NA             NA  NA    NA     0.300",
                 "omega[0.05,1] 0.509 0.301  0.028  0.508 0.983     0.01348          0.045 500    NA     0.400"
               ))
  expect_equal(capture_output_lines(
    add_column(estimates_table, column_title = "Models", column_values = c(1:4), column_position = 1), print = TRUE, width = 150),
    c("                Models  Mean Median  0.025  0.95",
      "m                    1 0.152  0.152 -0.220 0.461",
      "omega[0,0.05]        2 1.000  1.000  1.000 1.000",
      "omega[0.05,0.5]      3 0.679  0.745  0.064 1.000",
      "omega[0.5,1]         4 0.529  0.483  0.023 1.000"
    ))
  expect_equal(capture_output_lines(
    add_column(inference_table, column_title = "BF2", column_values = inference_table[,4], column_position = 5, column_type = "BF"), print = TRUE, width = 150),
    c("      Models Prior prob. Post. prob. Inclusion BF Inclusion BF",
      "m        3/3       1.000       1.000          Inf          Inf",
      "omega    2/3       0.667       0.800        2.002        2.002"
    ))
  expect_equal(capture_output_lines(
    add_column(summary_table, column_title = "Distribution", column_values = c("A", "B", "C"), column_position = 2, column_type = "string"), print = TRUE, width = 150),
    c(" Model Distribution     Prior m                        Prior omega                    Prior prob. log(marglik) Post. prob. Inclusion BF",
      "     1            A    Normal(0, 1)                                                         0.333        -1.10       0.200        0.499",
      "     2            B  Normal(0, 0.5)     omega[one-sided: .05] ~ CumDirichlet(1, 1)          0.333        -0.61       0.325        0.964",
      "     3            C  Normal(0, 0.3) omega[one-sided: .5, .05] ~ CumDirichlet(1, 1, 1)       0.333        -0.24       0.475        1.809"
    ))


  ### test removing columns
  expect_error(remove_column(runjags_summary, column_position = 10),
               "The 'column_position' must be equal or lower than 9.")

  expect_equal(capture_output_lines(
    remove_column(inference_table, column_position = 1), print = TRUE, width = 150),
    c("      Prior prob. Post. prob. Inclusion BF",
      "m           1.000       1.000          Inf",
      "omega       0.667       0.800        2.002"
    ))


  ### test explanatory texts
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1)

  expect_equal(interpret(inference, mixed_posteriors, list(
    list(
      inference         = "m",
      samples           = "m",
      inference_name    = "effect",
      inference_BF_name = "BF_10",
      samples_name      = "y",
      samples_units     = NULL
    )
  ), "Test"), "Test found strong evidence in favor of the effect, BF_10 = Inf, with mean model-averaged estimate y = 0.152, 95% CI [-0.220,  0.525].")

  inference[["m"]][["BF"]] <- 1/5
  expect_equal(interpret(inference, mixed_posteriors, list(
    list(
      inference           = "m",
      samples             = "m",
      inference_name      = "effect",
      inference_BF_name   = "BF_10",
      samples_name        = "y",
      samples_units       = "mm",
      samples_conditional = TRUE
    ),
    list(
      inference           = "omega",
      inference_name      = "bias",
      inference_BF_name   = "BF_pb"
    )
  ), "Test2"), "Test2 found moderate evidence against the effect, BF_10 = 0.200, with mean conditional estimate y = 0.152 mm, 95% CI [-0.220,  0.525]. Test2 found weak evidence in favor of the bias, BF_pb = 2.00.")

})

# skip the rest as it takes too long
skip_on_cran()

test_that("Summary tables functions work (formulas + factors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(60, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_fac2t)
  formula_list1 <- list(mu = ~ x_cont1 + x_fac3t)
  formula_list2 <- list(mu = ~ x_fac3o)
  formula_list3 <- list(mu = ~ x_cont1 * x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac2t"         = prior_factor("normal",  contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list2 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )
  formula_prior_list3 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
      "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)
  fit2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2, seed = 3)
  fit3 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3, seed = 4)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)
  marglik2 <- JAGS_bridgesampling(
    fit2, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2)
  marglik3 <- JAGS_bridgesampling(
    fit3, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3)


  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1), prior_weights = 1),
    list(fit = fit2, marglik = marglik2, fit_summary = runjags_estimates_table(fit2), prior_weights = 1),
    list(fit = fit3, marglik = marglik3, fit_summary = runjags_estimates_table(fit3), prior_weights = 1)
  )
  models <- models_inference(models)


  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1:x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1:x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)



  ### checking summary functions
  # model summary
  model_summary <- model_summary_table(models[[4]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  ", "  "))
  expect_equal(model_summary[,4], c("Parameter prior distributions","(mu) intercept ~ Normal(0, 5)","(mu) x_cont1 ~ Normal(0, 1)","(mu) x_fac3o ~ orthonormal contrast: mNormal(0, 1)","(mu) x_cont1:x_fac3o ~ orthonormal contrast: mNormal(0, 1)","sigma ~ Lognormal(0, 1)"))

  model_summary2 <- model_summary_table(models[[4]], formula_prefix = FALSE)
  expect_equal(model_summary2[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  ", "  "))
  expect_equal(model_summary2[,4], c("Parameter prior distributions","intercept ~ Normal(0, 5)","x_cont1 ~ Normal(0, 1)","x_fac3o ~ orthonormal contrast: mNormal(0, 1)","x_cont1:x_fac3o ~ orthonormal contrast: mNormal(0, 1)","sigma ~ Lognormal(0, 1)"))


  # runjags summary
  runjags_summary <- models[[2]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("(mu) intercept",  "(mu) x_cont1", "(mu) x_fac3t[B]", "(mu) x_fac3t[C]", "sigma"))
  expect_equal(unname(unlist(runjags_summary[3,])), c(5.746362e-03,  2.808364e-01, -5.496105e-01,  1.058318e-02,  5.504860e-01,  4.142589e-03,  1.500000e-02,  4.596000e+03,  1.000580e+00), tolerance = 1e-3)

  runjags_summary2 <- runjags_estimates_table(fit1, formula_prefix = FALSE)
  expect_equal(colnames(runjags_summary2), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary2), c("intercept",  "x_cont1", "x_fac3t[B]", "x_fac3t[C]", "sigma"))
  expect_equal(unname(unlist(runjags_summary2[3,])), c(5.746362e-03,  2.808364e-01, -5.496105e-01,  1.058318e-02,  5.504860e-01,  4.142589e-03,  1.500000e-02,  4.596000e+03,  1.000580e+00), tolerance = 1e-3)

  runjags_summary <- models[[4]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("(mu) intercept", "(mu) x_cont1", "(mu) x_fac3o[1]", "(mu) x_fac3o[2]", "(mu) x_cont1:x_fac3o[1]", "(mu) x_cont1:x_fac3o[2]", "sigma" ))
  expect_equal(unname(unlist(runjags_summary[1,])), c(1.876569e-01, 1.210763e-01, -5.091384e-02, 1.878474e-01, 4.285015e-01, 9.894116e-04, 8.000000e-03, 1.497500e+04, 1.000068e+00), tolerance = 1e-3)


  # ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_cont1", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("(mu) x_cont1", "(mu) x_fac3t[B]", "(mu) x_fac3t[C]", "(mu) x_fac3o[1]", "(mu) x_fac3o[2]", "(mu) x_cont1:x_fac3o[1]", "(mu) x_cont1:x_fac3o[2]"))
  expect_equal(unname(unlist(estimates_table[1,])), c(0.1224567, 0.0000000, 0.0000000, 0.4794182), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[3,])), c( 0.0397569,  0.0000000, -0.2895047,  0.4087159), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[5,])), c(-0.004121766,  0.000000000, -0.215131954,  0.036829714), tolerance = 1e-3)

  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_cont1", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("(mu) x_cont1", "(mu) x_fac3o[1]", "(mu) x_fac3o[2]", "(mu) x_cont1:x_fac3o[1]", "(mu) x_cont1:x_fac3o[2]"))
  expect_equal(unname(unlist(estimates_table[1,])), c(0.1224567, 0.0000000, 0.0000000, 0.4794182), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[3,])), c(-0.004121766,  0.000000000, -0.215131954,  0.036829714), tolerance = 1e-3)

  # ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(rownames(inference_table), c("(mu) x_cont1", "(mu) x_fac2t", "(mu) x_fac3t", "(mu) x_fac3o", "(mu) x_cont1:x_fac3o"))
  expect_equal(unname(unlist(inference_table[,1])),    c(2, 1, 1, 2, 1))
  expect_equal(unname(unlist(inference_table[,2])),    c(0.50, 0.25, 0.25, 0.50, 0.25))
  expect_equal(unname(unlist(inference_table[,3])),    c(0.37435772, 0.52598137, 0.33962193, 0.13439670, 0.03473579), tolerance = 1e-3)
  expect_equal(unname(as.vector(inference_table[,4])), c(0.5983575, 3.3288651, 1.5428523, 0.1552636, 0.1079573), tolerance = 1e-3)

  # ensemble summary
  summary_table <- ensemble_summary_table(models, c("mu_x_cont1", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"))
  expect_equal(colnames(summary_table), c("Model", "(mu) x_cont1", "(mu) x_fac3o", "(mu) x_cont1:x_fac3o", "prior_prob", "marglik", "post_prob", "inclusion_BF"))
  expect_equal(unname(as.vector(summary_table[,1])), c(1, 2, 3, 4))
  expect_equal(unname(as.vector(summary_table[,2])), c("", "Normal(0, 1)", "", "Normal(0, 1)"))
  expect_equal(unname(as.vector(summary_table[,3])), c("", "", "orthonormal contrast: mNormal(0, 1)", "orthonormal contrast: mNormal(0, 1)"))
  expect_equal(unname(as.vector(summary_table[,4])), c("", "", "", "orthonormal contrast: mNormal(0, 1)"))
  expect_equal(unname(as.vector(summary_table[,5])), c(0.25, 0.25, 0.25, 0.25), tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,6])), c(-88.22395, -88.66138, -89.88744, -90.94144),     tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,7])), c(0.52598137, 0.33962193, 0.09966091, 0.03473579), tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,8])), c(3.3288651, 1.5428523, 0.3320779, 0.1079573),     tolerance = 1e-3)

  # ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("mu_x_cont1", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"))
  expect_equal(colnames(diagnostics_table), c("Model", "(mu) x_cont1", "(mu) x_fac3o", "(mu) x_cont1:x_fac3o", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))

  expect_equal(unname(as.vector(diagnostics_table[,1])), c(1, 2, 3, 4))
  expect_equal(unname(as.vector(diagnostics_table[,2])), c("", "Normal(0, 1)", "", "Normal(0, 1)"))
  expect_equal(unname(as.vector(diagnostics_table[,3])), c("", "", "orthonormal contrast: mNormal(0, 1)", "orthonormal contrast: mNormal(0, 1)"))
  expect_equal(unname(as.vector(diagnostics_table[,4])), c("", "", "", "orthonormal contrast: mNormal(0, 1)"))
  expect_equal(unname(as.vector(diagnostics_table[,5])), c(0.003223670, 0.004142589, 0.001676136, 0.001959310), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,6])), c(0.013, 0.017, 0.011, 0.011), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,7])), c(5559, 3526, 8660, 7969))
  expect_equal(unname(as.vector(diagnostics_table[,8])), c(1.001154, 1.000955, 1.000125, 1.000658), tolerance = 1e-3)


  ### test additional settings
  # transformations of orthonormal contrast to differences from the mean
  runjags_summary_t <- runjags_estimates_table(fit3, transform_factors = TRUE)
  expect_equal(colnames(runjags_summary_t), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary_t), c("(mu) intercept","(mu) x_cont1","(mu) x_fac3o [dif: A]","(mu) x_fac3o [dif: B]","(mu) x_fac3o [dif: C]", "(mu) x_cont1:x_fac3o [dif: A]", "(mu) x_cont1:x_fac3o [dif: B]", "(mu) x_cont1:x_fac3o [dif: C]", "sigma" ))
  expect_equal(capture_output_lines(runjags_summary_t, print = TRUE, width = 150),
               c("                                Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD   ESS R-hat",
                 "(mu) intercept                 0.188 0.121 -0.051  0.188 0.429     0.00099          0.008 14975 1.000",
                 "(mu) x_cont1                   0.324 0.140  0.047  0.324 0.597     0.00112          0.008 15680 1.000",
                 "(mu) x_fac3o [dif: A]         -0.010 0.168 -0.337 -0.011 0.321     0.00132          0.008 15278 1.000",
                 "(mu) x_fac3o [dif: B]         -0.064 0.170 -0.397 -0.064 0.270     0.00134          0.008 15081 1.000",
                 "(mu) x_fac3o [dif: C]          0.074 0.167 -0.251  0.072 0.404     0.00132          0.008 15630 1.000",
                 "(mu) x_cont1:x_fac3o [dif: A] -0.283 0.197 -0.668 -0.283 0.105     0.00156          0.008 15581 1.000",
                 "(mu) x_cont1:x_fac3o [dif: B]  0.164 0.194 -0.221  0.164 0.539     0.00153          0.008 14954 1.000",
                 "(mu) x_cont1:x_fac3o [dif: C]  0.119 0.202 -0.275  0.118 0.521     0.00160          0.008 15372 1.000",
                 "sigma                          0.925 0.090  0.770  0.918 1.119     0.00100          0.011  7969 1.001"
               ))


  estimates_table_t <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_cont1", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"), probs = c(.025, 0.95), transform_factors = TRUE)
  expect_equal(colnames(estimates_table_t), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table_t), c("(mu) x_cont1","(mu) x_fac3o [dif: A]", "(mu) x_fac3o [dif: B]", "(mu) x_fac3o [dif: C]", "(mu) x_cont1:x_fac3o [dif: A]", "(mu) x_cont1:x_fac3o [dif: B]", "(mu) x_cont1:x_fac3o [dif: C]"))
  expect_equal(capture_output_lines(estimates_table_t, print = TRUE, width = 150),
               c("                                Mean Median  0.025  0.95",
                 "(mu) x_cont1                   0.122  0.000  0.000 0.479",
                 "(mu) x_fac3o [dif: A]         -0.003  0.000 -0.176 0.030",
                 "(mu) x_fac3o [dif: B]         -0.003  0.000 -0.181 0.039",
                 "(mu) x_fac3o [dif: C]          0.007  0.000 -0.105 0.100",
                 "(mu) x_cont1:x_fac3o [dif: A] -0.010  0.000 -0.183 0.000",
                 "(mu) x_cont1:x_fac3o [dif: B]  0.006  0.000  0.000 0.000",
                 "(mu) x_cont1:x_fac3o [dif: C]  0.005  0.000  0.000 0.000"
               ))



  ### test print functions
  expect_equal(capture_output_lines(model_summary, print = TRUE, width = 150),
               c("                                                                                             ",
                 " Model               4                                          Parameter prior distributions",
                 " Prior prob.     0.250                   (mu) intercept ~ Normal(0, 5)                       ",
                 " log(marglik)   -90.94                     (mu) x_cont1 ~ Normal(0, 1)                       ",
                 " Post. prob.     0.035                     (mu) x_fac3o ~ orthonormal contrast: mNormal(0, 1)",
                 " Inclusion BF    0.108             (mu) x_cont1:x_fac3o ~ orthonormal contrast: mNormal(0, 1)",
                 "                                                  sigma ~ Lognormal(0, 1)                    "
               ))
  expect_equal(capture_output_lines(runjags_summary,   print = TRUE, width = 150),
               c("                          Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD   ESS R-hat",
                 "(mu) intercept           0.188 0.121 -0.051  0.188 0.429     0.00099          0.008 14975 1.000",
                 "(mu) x_cont1             0.324 0.140  0.047  0.324 0.597     0.00112          0.008 15680 1.000",
                 "(mu) x_fac3o[1]          0.097 0.207 -0.314  0.096 0.508     0.00166          0.008 15450 1.000",
                 "(mu) x_fac3o[2]         -0.012 0.205 -0.412 -0.013 0.393     0.00164          0.008 15720 1.000",
                 "(mu) x_cont1:x_fac3o[1] -0.032 0.243 -0.507 -0.033 0.448     0.00196          0.008 15383 1.000",
                 "(mu) x_cont1:x_fac3o[2] -0.347 0.242 -0.818 -0.347 0.128     0.00193          0.008 15659 1.000",
                 "sigma                    0.925 0.090  0.770  0.918 1.119     0.00100          0.011  7969 1.001"

               ))
  expect_equal(capture_output_lines(estimates_table,   print = TRUE, width = 150),
               c("                          Mean Median  0.025  0.95",
                 "(mu) x_cont1             0.122  0.000  0.000 0.479",
                 "(mu) x_fac3o[1]          0.007  0.000 -0.145 0.125",
                 "(mu) x_fac3o[2]         -0.004  0.000 -0.215 0.037",
                 "(mu) x_cont1:x_fac3o[1] -0.001  0.000  0.000 0.000",
                 "(mu) x_cont1:x_fac3o[2] -0.013  0.000 -0.224 0.000"

               ))
  expect_equal(capture_output_lines(inference_table,   print = TRUE, width = 150),
               c("                     Models Prior prob. Post. prob. Inclusion BF",
                 "(mu) x_cont1            2/4       0.500       0.374        0.598",
                 "(mu) x_fac2t            1/4       0.250       0.526        3.329",
                 "(mu) x_fac3t            1/4       0.250       0.340        1.543",
                 "(mu) x_fac3o            2/4       0.500       0.134        0.155",
                 "(mu) x_cont1:x_fac3o    1/4       0.250       0.035        0.108"

               ))
  expect_equal(capture_output_lines(summary_table,     print = TRUE, width = 150),
               c(" Model Prior (mu) x_cont1          Prior (mu) x_fac3o              Prior (mu) x_cont1:x_fac3o     Prior prob. log(marglik) Post. prob. Inclusion BF",
                 "     1                                                                                                  0.250       -88.22       0.526        3.329",
                 "     2       Normal(0, 1)                                                                               0.250       -88.66       0.340        1.543",
                 "     3                    orthonormal contrast: mNormal(0, 1)                                           0.250       -89.89       0.100        0.332",
                 "     4       Normal(0, 1) orthonormal contrast: mNormal(0, 1) orthonormal contrast: mNormal(0, 1)       0.250       -90.94       0.035        0.108"
               ))
  expect_equal(capture_output_lines(diagnostics_table, print = TRUE, width = 180),
               c(" Model Prior (mu) x_cont1          Prior (mu) x_fac3o              Prior (mu) x_cont1:x_fac3o     max[error(MCMC)] max[error(MCMC)/SD] min(ESS) max(R-hat)",
                 "     1                                                                                                     0.00322               0.013     5559      1.001",
                 "     2       Normal(0, 1)                                                                                  0.00414               0.017     3526      1.001",
                 "     3                    orthonormal contrast: mNormal(0, 1)                                              0.00168               0.011     8660      1.000",
                 "     4       Normal(0, 1) orthonormal contrast: mNormal(0, 1) orthonormal contrast: mNormal(0, 1)          0.00196               0.011     7969      1.001"
               ))

})

test_that("Summary tables functions work (indepdent factors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac2i = factor(rep(c("A", "B"), 30), levels = c("A", "B"))
  )
  data <- list(
    y = rnorm(60, ifelse(data_formula$x_fac2i == "A", 0.0, -0.2), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_fac2i - 1)
  formula_list1 <- list(mu = ~ x_fac2i - 1)

  formula_prior_list0 <- list(
    mu    = list(
      "x_fac2i" = prior_factor("spike", contrast = "independent", list(0))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "x_fac2i" = prior_factor("normal", contrast = "independent", list(0, 1/4))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1), prior_weights = 1)
  )
  models <- models_inference(models)


  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_fac2i"),
    is_null_list = list(
      "mu_x_fac2i" = c(TRUE,  FALSE)
    ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_fac2i"),
    is_null_list = list(
      "mu_x_fac2i" = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)



  ### checking summary functions
  # model summary
  model_summary <- model_summary_table(models[[2]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  "))
  expect_equal(model_summary[,4], c("Parameter prior distributions", "(mu) x_fac2i ~ independent contrast: Normal(0, 0.25)","sigma ~ Lognormal(0, 1)", "", ""))

  # runjags summary
  runjags_summary <- models[[2]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("(mu) x_fac2i[A]", "(mu) x_fac2i[B]", "sigma" ))
  expect_equal(unname(unlist(runjags_summary[1,])), c(1.734095e-01, 1.340447e-01, -9.293281e-02, 1.747751e-01, 4.347246e-01, 1.067352e-03, 8.000000e-03, 1.577200e+04, 1.000033e+00), tolerance = 1e-3)

  # ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_fac2i"), probs = c(.025, 0.95))
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("(mu) x_fac2i[A]", "(mu) x_fac2i[B]"))
  expect_equal(unname(unlist(estimates_table[1,])), c(0.10208451, 0.03621004, -0.06041045, 0.35346681), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[2,])), c(-0.09355933, -0.01700284, -0.38746858,  0.02836426), tolerance = 1e-3)

  # ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(rownames(inference_table), c("(mu) x_fac2i"))
  expect_equal(unname(unlist(inference_table[,1])),    1)
  expect_equal(unname(unlist(inference_table[,2])),    0.5)
  expect_equal(unname(unlist(inference_table[,3])),    0.5876797, tolerance = 1e-3)
  expect_equal(unname(as.vector(inference_table[,4])), 1.425299, tolerance = 1e-3)

  # ensemble summary
  summary_table <- ensemble_summary_table(models, c("mu_x_fac2i"))
  expect_equal(colnames(summary_table), c("Model", "(mu) x_fac2i", "prior_prob", "marglik", "post_prob", "inclusion_BF"))
  expect_equal(unname(as.vector(summary_table[,1])), c(1, 2))
  expect_equal(unname(as.vector(summary_table[,2])), c("","independent contrast: Normal(0, 0.25)"))
  expect_equal(unname(as.vector(summary_table[,3])), c(0.5, 0.5), tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,4])), c(-79.15494, -78.80056), tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,5])), c(0.4123203, 0.5876797), tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,6])), c(0.7016071, 1.4252991), tolerance = 1e-3)

  # ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("mu_x_fac2i"), remove_spike_0 = FALSE)
  expect_equal(colnames(diagnostics_table), c("Model", "(mu) x_fac2i", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))

  expect_equal(unname(as.vector(diagnostics_table[,1])), c(1, 2))
  expect_equal(unname(as.vector(diagnostics_table[,2])), c("independent contrast: Spike(0)","independent contrast: Normal(0, 0.25)"))
  expect_equal(unname(as.vector(diagnostics_table[,3])), c(0.0008277888, 0.0010673515), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,4])), c(0.010, 0.011), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,5])), c(9564, 8145))

})

test_that("Summary tables functions work (meandif factors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(2)

  data_formula <- data.frame(
    x_fac3 = factor(rep(c("A", "B", "C"), 60), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(180, ifelse(data_formula$x_fac3 == "A", -0.2, ifelse(data_formula$x_fac3 == "B", 0.0, 0.2)), 1),
    N = 180
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ 1)
  formula_list1 <- list(mu = ~ x_fac3)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_fac3"    = prior_factor("mnormal", contrast = "meandif", list(0, 1/5))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)

  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, fit_summary = runjags_estimates_table(fit1), prior_weights = 1)
  )
  models <- models_inference(models)


  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_fac3"),
    is_null_list = list(
      "mu_x_fac3" = c(TRUE,  FALSE)
    ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_fac3"),
    is_null_list = list(
      "mu_x_fac3" = c(TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)



  ### checking summary functions
  # model summary
  model_summary <- model_summary_table(models[[2]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  "))
  expect_equal(model_summary[,4], c("Parameter prior distributions", "(mu) intercept ~ Normal(0, 5)", "(mu) x_fac3 ~ mean difference contrast: mNormal(0, 0.2)","sigma ~ Lognormal(0, 1)", ""))

  # runjags summary
  runjags_summary <- models[[2]]$fit_summary
  expect_equal(colnames(runjags_summary), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(runjags_summary), c("(mu) intercept", "(mu) x_fac3[1]", "(mu) x_fac3[2]", "sigma"))
  expect_equal(unname(unlist(runjags_summary[1,])), c(2.616574e-02,8.256672e-02,-1.369357e-01,2.621934e-02,1.851191e-01,6.471943e-04,8.000000e-03,1.627600e+04,9.999001e-01), tolerance = 1e-3)

  # ensemble estimates
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("mu_x_fac3"), probs = c(.025, 0.95), transform_factors = TRUE)
  expect_equal(colnames(estimates_table), c("Mean", "Median", "0.025",  "0.95"))
  expect_equal(rownames(estimates_table), c("(mu) x_fac3 [dif: A]", "(mu) x_fac3 [dif: B]", "(mu) x_fac3 [dif: C]"))
  expect_equal(unname(unlist(estimates_table[1,])), c(-0.2074503, -0.2206674, -0.4204564,  0.0000000), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[2,])), c(0.023169431,  0.008606852, -0.163666847,  0.185934433), tolerance = 1e-3)
  expect_equal(unname(unlist(estimates_table[3,])), c(0.1842808, 0.1938991, 0.0000000, 0.3678031), tolerance = 1e-3)

  # ensemble inference
  inference_table <- ensemble_inference_table(inference, names(inference))
  expect_equal(colnames(inference_table), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(rownames(inference_table), c("(mu) x_fac3"))
  expect_equal(unname(unlist(inference_table[,1])),    1)
  expect_equal(unname(unlist(inference_table[,2])),    0.5)
  expect_equal(unname(unlist(inference_table[,3])),    0.8737537, tolerance = 1e-3)
  expect_equal(unname(as.vector(inference_table[,4])), 6.921025, tolerance = 1e-3)

  # ensemble summary
  summary_table <- ensemble_summary_table(models, c("mu_x_fac3"))
  expect_equal(colnames(summary_table), c("Model", "(mu) x_fac3", "prior_prob", "marglik", "post_prob", "inclusion_BF"))
  expect_equal(unname(as.vector(summary_table[,1])), c(1, 2))
  expect_equal(unname(as.vector(summary_table[,2])), c("","mean difference contrast: mNormal(0, 0.2)"))
  expect_equal(unname(as.vector(summary_table[,3])), c(0.5, 0.5), tolerance = 1e-4)
  expect_equal(unname(as.vector(summary_table[,4])), c(-282.5467, -280.6121), tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,5])), c(0.1262463, 0.8737537), tolerance = 1e-3)
  expect_equal(unname(as.vector(summary_table[,6])), c(0.1444873, 6.9210254), tolerance = 1e-3)

  # ensemble diagnostics
  diagnostics_table <- ensemble_diagnostics_table(models, c("mu_x_fac3"), remove_spike_0 = FALSE)
  expect_equal(colnames(diagnostics_table), c("Model", "(mu) x_fac3", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))
  expect_equal(unname(as.vector(diagnostics_table[,1])), c(1, 2))
  expect_equal(unname(as.vector(diagnostics_table[,2])), c("", "mean difference contrast: mNormal(0, 0.2)"))
  expect_equal(unname(as.vector(diagnostics_table[,3])), c(0.0006707336, 0.0007978420), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,4])), c(0.01, 0.01), tolerance = 1e-3)
  expect_equal(unname(as.vector(diagnostics_table[,5])), c(9676, 9871))

})

test_that("Summary tables functions work (spike and slab priors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(60, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_cont1 + x_fac2t + x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_cont1"   = prior_spike_and_slab(
        prior_parameter = prior("normal", list(0, 0.5)),
        prior_inclusion = prior("beta", list(1, 1))
      ),
      "x_fac2t"   = prior_spike_and_slab(
        prior_parameter = prior_factor("normal",  contrast = "treatment", list(0, 1)),
        prior_inclusion = prior("beta", list(1, 1))
      ),
      "x_fac3o"   = prior_spike_and_slab(
        prior_parameter = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
        prior_inclusion = prior("spike", list(.5))
      )
    )
  )


  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)

  marglik0 <- list(logml = 0)
  class(marglik0) <- "bridge"
  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, fit_summary = runjags_estimates_table(fit0), prior_weights = 1)
  )
  models <- models_inference(models)

  model_summary_table(models[[1]])
  ### checking summary functions
  # model summary
  model_summary <- model_summary_table(models[[1]])
  expect_equal(model_summary[,1], c("Model  ", "Prior prob.  ", "log(marglik)  ", "Post. prob.  ", "Inclusion BF  ", "  "))
  expect_equal(model_summary[,4], c("Parameter prior distributions", "(mu) intercept ~ Normal(0, 5)", "(mu) x_cont1 ~ Normal(0, 0.5) * Beta(1, 1)", "(mu) x_fac2t ~ treatment contrast: Normal(0, 1) * Beta(1, 1)", "(mu) x_fac3o ~ orthonormal contrast: mNormal(0, 1) * Spike(0.5)", "sigma ~ Lognormal(0, 1)"))

  model_estimates <- runjags_estimates_table(fit0)
  expect_equal(colnames(model_estimates), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(model_estimates), c("(mu) intercept", "(mu) x_cont1", "(mu) x_cont1 (inclusion)", "(mu) x_fac2t[B]", "(mu) x_fac2t (inclusion)", "(mu) x_fac3o[1]", "(mu) x_fac3o[2]", "(mu) x_fac3o (inclusion)", "sigma"))
  expect_equal(unname(unlist(model_estimates[7,])), c(-1.776437e-03, 4.269207e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.428435e-04, 8.000000e-03, 1.550600e+04, 1.001065e+00), tolerance = 1e-3)

  model_estimates <- runjags_estimates_table(fit0, transform_factors = TRUE, conditional = TRUE)
  expect_equal(colnames(model_estimates), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(model_estimates), c("(mu) intercept", "(mu) x_cont1", "(mu) x_cont1 (inclusion)", "(mu) x_fac2t[B]", "(mu) x_fac2t (inclusion)", "(mu) x_fac3o [dif: A]", "(mu) x_fac3o [dif: B]", "(mu) x_fac3o [dif: C]", "(mu) x_fac3o (inclusion)", "sigma"))
  expect_equal(unname(unlist(model_estimates[8,])), c(0.0626582174, 0.1661778973, -0.2621073424, 0.0591205499, 0.3954805352,  NA, NA, NA, NA), tolerance = 1e-3)

  model_estimates <- runjags_estimates_table(fit0, transform_factors = TRUE, conditional = TRUE, remove_inclusion = TRUE)
  expect_equal(colnames(model_estimates), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(model_estimates), c("(mu) intercept", "(mu) x_cont1", "(mu) x_fac2t[B]", "(mu) x_fac3o [dif: A]", "(mu) x_fac3o [dif: B]", "(mu) x_fac3o [dif: C]", "sigma"))
  expect_equal(unname(unlist(model_estimates[2,])), c(3.040927e-01, 1.355633e-01, 3.256895e-02, 3.058346e-01, 5.678668e-01, 2.298187e-03, 1.300000e-02, 5.720000e+03, 1.000421e+00), tolerance = 1e-3)

  model_inference <- runjags_inference_table(fit0)
  expect_equal(colnames(model_inference), c("Parameter", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(model_inference[,1], c("(mu) x_cont1", "(mu) x_fac2t", "(mu) x_fac3o"))
  expect_equal(model_inference[,2], c(0.5, 0.5, 0.5))
  expect_equal(model_inference[,3], c(0.7798125, 0.1864375, 0.0399375), tolerance = 1e-3)
  expect_equal(model_inference[,4], c(3.54158388, 0.22916187, 0.04159885), tolerance = 1e-3)

  runjags_inference_empty <- runjags_inference_empty_table()
  expect_equivalent(nrow(runjags_inference_empty), 0)
  expect_equal(colnames(runjags_inference_empty), colnames(model_inference))
  expect_equal(capture_output_lines(runjags_inference_empty, width = 150)[1], capture_output_lines(model_inference, width = 150)[1])

})

test_that("Summary tables functions work (stan)",{

  skip_on_cran()
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes

  # prefitted model with RoBTT
  if(!file.exists(file.path("../results/fits", "fit_RoBTT.RDS")))
    skip(message = "Only runs locally")

  fit <- readRDS(file = file.path("../results/fits", "fit_RoBTT.RDS"))

  set.seed(1)

  ### checking summary functions
  model_estimates <- stan_estimates_table(fit)
  expect_equal(colnames(model_estimates), c("Mean", "SD", "lCI", "Median", "uCI", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(rownames(model_estimates), c("mu", "sigma2", "pooled_sigma", "sigma_i[1]", "sigma_i[2]", "mu_i[1]", "mu_i[2]" ))
  expect_equal(unname(unlist(model_estimates[1,])), c(1.43876353, 0.37708461, 0.81080656, 1.42486330, 2.15911838, 0.06223762, 0.16504949, 36.70892380, 1.01241771), tolerance = 1e-3)

})
