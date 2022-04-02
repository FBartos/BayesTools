context("JAGS formula")

test_that("JAGS formula works", {

  # check the posterior distributions with weak priors against a maximum likelihood estimates with ML
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)
  df_all <- data.frame(
    x_cont1 = rnorm(60),
    x_cont2 = rnorm(60),
    x_bin   = rbinom(60, 1, .5),
    x_fac2o = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(60, 0.1, 0.5) + 0.30 * df_all$x_cont1 - 0.15 * df_all$x_cont1 * df_all$x_cont2 + 0.2 * df_all$x_bin +
    ifelse(df_all$x_fac3t == "A", 0.2, ifelse(df_all$x_fac3t == "B", -0.2, 0)) +
    ifelse(df_all$x_fac3o == "A", 0.2, ifelse(df_all$x_fac3o == "B", -0.2, 0))
  prior_list_all <- list(
    "intercept"       = prior("normal", list(0, 5)),
    "x_cont1"         = prior("normal", list(0, 1)),
    "x_cont2"         = prior("normal", list(0, 1)),
    "x_cont1:x_cont2" = prior("normal", list(0, 1)),
    "x_fac2o"         = prior_factor("mcauchy", contrast = "orthonormal", list(0, 1)),
    "x_fac2t"         = prior_factor("normal",  contrast = "treatment",   list(0, 1)),
    "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
    "x_fac3t"         = prior_factor("uniform", contrast = "treatment",   list(-2, 2)),
    "x_fac2t:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0,  2)),
    "x_fac2o:x_fac3t" = prior_factor("normal",  contrast = "treatment",   list(0,  2)),
    "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0,  2)),
    "x_cont1:x_fac3t" = prior_factor("normal",  contrast = "treatment",   list(0,  2))
  )
  prior_list2  <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 1))
  )
  model_syntax <- paste0(
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n"
  )


  # simple linear regression ----
  formula_1      <- JAGS_formula(~ x_cont1, parameter = "mu", data = df_all[,"x_cont1", drop = FALSE], prior_list = prior_list_all[c("intercept", "x_cont1")])
  prior_list_1   <- c(formula_1$prior_list, prior_list2)
  model_syntax_1 <- JAGS_add_priors(paste0("model{", formula_1$formula_syntax, model_syntax, "}"), prior_list_1)
  data_1         <- c(formula_1$data, N = nrow(df_all), y = list(df_all$y))

  model_1   <- rjags::jags.model(file = textConnection(model_syntax_1), inits = JAGS_get_inits(prior_list_1, chains = 2, seed = 1), data = data_1, n.chains = 2, quiet = TRUE)
  samples_1 <- rjags::coda.samples(model = model_1, variable.names = JAGS_to_monitor(prior_list_1), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_1 <- do.call(rbind, samples_1)

  lm_1 <- stats::lm(y ~ x_cont1, data = df_all)

  expect_doppelganger("JAGS-formula-lm-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_1[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_1)["(Intercept)"], sd = summary(lm_1)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_1[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")
    curve(dnorm(x, mean = coef(lm_1)["x_cont1"], sd = summary(lm_1)$coefficients["x_cont1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_1[,"sigma"], freq = FALSE, main = "sigma")
    abline(v = sigma(lm_1), lwd = 3)
  })


  # linear regression with two continuous predictors and their interaction ----
  formula_2      <- JAGS_formula(~ x_cont1 * x_cont2, parameter = "mu", data = df_all[,c("x_cont1", "x_cont2")], prior_list = prior_list_all[c("intercept", "x_cont1", "x_cont2", "x_cont1:x_cont2")])
  prior_list_2   <- c(formula_2$prior_list, prior_list2)
  model_syntax_2 <- JAGS_add_priors(paste0("model{", formula_2$formula_syntax, model_syntax, "}"), prior_list_2)
  data_2         <- c(formula_2$data, N = nrow(df_all), y = list(df_all$y))

  model_2   <- rjags::jags.model(file = textConnection(model_syntax_2), inits = JAGS_get_inits(prior_list_2, chains = 2, seed = 1), data = data_2, n.chains = 2, quiet = TRUE)
  samples_2 <- rjags::coda.samples(model = model_2, variable.names = JAGS_to_monitor(prior_list_2), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_2 <- do.call(rbind, samples_2)

  lm_2 <- stats::lm(y ~ x_cont1 * x_cont2, data = df_all)

  expect_doppelganger("JAGS-formula-lm-2", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_2[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")
    curve(dnorm(x, mean = coef(lm_2)["x_cont1"], sd = summary(lm_2)$coefficients["x_cont1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_2[,"mu_x_cont2"], freq = FALSE, main = "x_cont2")
    curve(dnorm(x, mean = coef(lm_2)["x_cont2"], sd = summary(lm_2)$coefficients["x_cont2", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_2[,"mu_x_cont1__xXx__x_cont2"], freq = FALSE, main = "x_cont1:x_cont2")
    curve(dnorm(x, mean = coef(lm_2)["x_cont1:x_cont2"], sd = summary(lm_2)$coefficients["x_cont1:x_cont2", "Std. Error"]), add = TRUE, lwd = 2)
  })


  # linear regression with a dummy factor (2 levels) ----
  formula_3      <- JAGS_formula(~ x_fac2t, parameter = "mu", data = df_all[,"x_fac2t",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac2t")])
  prior_list_3   <- c(formula_3$prior_list, prior_list2)
  model_syntax_3 <- JAGS_add_priors(paste0("model{", formula_3$formula_syntax, model_syntax, "}"), prior_list_3)
  data_3         <- c(formula_3$data, N = nrow(df_all), y = list(df_all$y))

  model_3   <- rjags::jags.model(file = textConnection(model_syntax_3), inits = JAGS_get_inits(prior_list_3, chains = 2, seed = 1), data = data_3, n.chains = 2, quiet = TRUE)
  samples_3 <- rjags::coda.samples(model = model_3, variable.names = JAGS_to_monitor(prior_list_3), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_3 <- do.call(rbind, samples_3)

  lm_3 <- stats::lm(y ~ x_fac2t, data = df_all)

  expect_doppelganger("JAGS-formula-lm-3", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_3[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_3)["(Intercept)"], sd = summary(lm_3)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_3[,"mu_x_fac2t"], freq = FALSE, main = "x_fac2t")
    curve(dnorm(x, mean = coef(lm_3)["x_fac2tB"], sd = summary(lm_3)$coefficients["x_fac2tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_3[,"sigma"], freq = FALSE, main = "sigma")
    abline(v = sigma(lm_3), lwd = 3)
  })


  # linear regression with an orthornormal factor (2 levels) ----
  formula_4      <- JAGS_formula(~ x_fac2o, parameter = "mu", data = df_all[,"x_fac2o",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac2o")])
  prior_list_4   <- c(formula_4$prior_list, prior_list2)
  model_syntax_4 <- JAGS_add_priors(paste0("model{", formula_4$formula_syntax, model_syntax, "}"), prior_list_4)
  data_4         <- c(formula_4$data, N = nrow(df_all), y = list(df_all$y))

  model_4   <- rjags::jags.model(file = textConnection(model_syntax_4), inits = JAGS_get_inits(prior_list_4, chains = 2, seed = 1), data = data_4, n.chains = 2, quiet = TRUE)
  samples_4 <- rjags::coda.samples(model = model_4, variable.names = JAGS_to_monitor(prior_list_4), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_4 <- do.call(rbind, samples_4)

  df_4 <- df_all
  contrasts(df_4$x_fac2o) <- contr.orthonormal(levels(df_4$x_fac2o))
  lm_4 <- stats::lm(y ~ x_fac2o, data = df_4)

  expect_doppelganger("JAGS-formula-lm-4", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_4[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_4)["(Intercept)"], sd = summary(lm_4)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_4[,"mu_x_fac2o"], freq = FALSE, main = "x_fac2o")
    curve(dnorm(x, mean = coef(lm_4)["x_fac2o1"], sd = summary(lm_4)$coefficients["x_fac2o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_4[,"sigma"], freq = FALSE, main = "sigma")
    abline(v = sigma(lm_4), lwd = 3)
  })


  # linear regression with a dummy factor (3 levels) ----
  formula_5      <- JAGS_formula(~ x_fac3t, parameter = "mu", data = df_all[,"x_fac3t",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac3t")])
  prior_list_5   <- c(formula_5$prior_list, prior_list2)
  model_syntax_5 <- JAGS_add_priors(paste0("model{", formula_5$formula_syntax, model_syntax, "}"), prior_list_5)
  data_5         <- c(formula_5$data, N = nrow(df_all), y = list(df_all$y))

  model_5   <- rjags::jags.model(file = textConnection(model_syntax_5), inits = JAGS_get_inits(prior_list_5, chains = 2, seed = 1), data = data_5, n.chains = 2, quiet = TRUE)
  samples_5 <- rjags::coda.samples(model = model_5, variable.names = JAGS_to_monitor(prior_list_5), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_5 <- do.call(rbind, samples_5)

  lm_5 <- stats::lm(y ~ x_fac3t, data = df_all)

  expect_doppelganger("JAGS-formula-lm-5", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_5[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_5)["(Intercept)"], sd = summary(lm_5)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_5[,"mu_x_fac3t[1]"], freq = FALSE, main = "x_fac3t[1]")
    curve(dnorm(x, mean = coef(lm_5)["x_fac3tB"], sd = summary(lm_5)$coefficients["x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_5[,"mu_x_fac3t[2]"], freq = FALSE, main = "x_fac3t[2]")
    curve(dnorm(x, mean = coef(lm_5)["x_fac3tC"], sd = summary(lm_5)$coefficients["x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)
  })


  # linear regression with an orthornormal factor (3 levels) ----
  formula_6      <- JAGS_formula(~ x_fac3o, parameter = "mu", data = df_all[,"x_fac3o",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac3o")])
  prior_list_6   <- c(formula_6$prior_list, prior_list2)
  model_syntax_6 <- JAGS_add_priors(paste0("model{", formula_6$formula_syntax, model_syntax, "}"), prior_list_6)
  data_6         <- c(formula_6$data, N = nrow(df_all), y = list(df_all$y))

  model_6   <- rjags::jags.model(file = textConnection(model_syntax_6), inits = JAGS_get_inits(prior_list_6, chains = 2, seed = 1), data = data_6, n.chains = 2, quiet = TRUE)
  samples_6 <- rjags::coda.samples(model = model_6, variable.names = JAGS_to_monitor(prior_list_6), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_6 <- do.call(rbind, samples_6)

  df_6 <- df_all
  contrasts(df_6$x_fac3o) <- contr.orthonormal(levels(df_6$x_fac3o))
  lm_6 <- stats::lm(y ~ x_fac3o, data = df_6)

  expect_doppelganger("JAGS-formula-lm-6", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_6[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_6)["(Intercept)"], sd = summary(lm_6)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_6[,"mu_x_fac3o[1]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_6)["x_fac3o1"], sd = summary(lm_6)$coefficients["x_fac3o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_6[,"mu_x_fac3o[2]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_6)["x_fac3o2"], sd = summary(lm_6)$coefficients["x_fac3o2", "Std. Error"]), add = TRUE, lwd = 2)
  })


  # linear regression with an orthornormal interaction between factors ----
  formula_7      <- JAGS_formula(~ x_fac2t * x_fac3o, parameter = "mu", data = df_all[,c("x_fac2t", "x_fac3o")], prior_list = prior_list_all[c("intercept", "x_fac2t", "x_fac3o", "x_fac2t:x_fac3o")])
  prior_list_7   <- c(formula_7$prior_list, prior_list2)
  model_syntax_7 <- JAGS_add_priors(paste0("model{", formula_7$formula_syntax, model_syntax, "}"), prior_list_7)
  data_7         <- c(formula_7$data, N = nrow(df_all), y = list(df_all$y))

  model_7   <- rjags::jags.model(file = textConnection(model_syntax_7), inits = JAGS_get_inits(prior_list_7, chains = 2, seed = 1), data = data_7, n.chains = 2, quiet = TRUE)
  samples_7 <- rjags::coda.samples(model = model_7, variable.names = JAGS_to_monitor(prior_list_7), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_7 <- do.call(rbind, samples_7)

  df_7 <- df_all
  contrasts(df_7$x_fac3o) <- contr.orthonormal(levels(df_7$x_fac3o))
  lm_7 <- stats::lm(y ~ x_fac2t * x_fac3o, data = df_7)

  expect_doppelganger("JAGS-formula-lm-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    hist(samples_7[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_7)["(Intercept)"], sd = summary(lm_7)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_7[,"mu_x_fac3o[1]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_7)["x_fac3o1"], sd = summary(lm_7)$coefficients["x_fac3o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_7[,"mu_x_fac3o[2]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_7)["x_fac3o2"], sd = summary(lm_7)$coefficients["x_fac3o2", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_7[,"mu_x_fac2t"], freq = FALSE, main = "x_fac2t")
    curve(dnorm(x, mean = coef(lm_7)["x_fac2tB"], sd = summary(lm_7)$coefficients["x_fac2tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_7[,"mu_x_fac2t__xXx__x_fac3o[1]"], freq = FALSE, main = "x_fac2t:x_fac3o")
    curve(dnorm(x, mean = coef(lm_7)["x_fac2tB:x_fac3o1"], sd = summary(lm_7)$coefficients["x_fac2tB:x_fac3o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_7[,"mu_x_fac2t__xXx__x_fac3o[2]"], freq = FALSE, main = "x_fac2t:x_fac3o")
    curve(dnorm(x, mean = coef(lm_7)["x_fac2tB:x_fac3o2"], sd = summary(lm_7)$coefficients["x_fac2tB:x_fac3o2", "Std. Error"]), add = TRUE, lwd = 2)

  })


  # linear regression with a treatment interaction between factors ----
  formula_8      <- JAGS_formula(~  x_fac2o * x_fac3t, parameter = "mu", data = df_all[,c("x_fac2o", "x_fac3t")], prior_list = prior_list_all[c("intercept", "x_fac2o", "x_fac3t", "x_fac2o:x_fac3t")])
  prior_list_8   <- c(formula_8$prior_list, prior_list2)
  model_syntax_8 <- JAGS_add_priors(paste0("model{", formula_8$formula_syntax, model_syntax, "}"), prior_list_8)
  data_8         <- c(formula_8$data, N = nrow(df_all), y = list(df_all$y))

  model_8   <- rjags::jags.model(file = textConnection(model_syntax_8), inits = JAGS_get_inits(prior_list_8, chains = 2, seed = 1), data = data_8, n.chains = 2, quiet = TRUE)
  samples_8 <- rjags::coda.samples(model = model_8, variable.names = JAGS_to_monitor(prior_list_8), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_8 <- do.call(rbind, samples_8)

  df_8 <- df_all
  contrasts(df_8$x_fac2o) <- contr.orthonormal(levels(df_8$x_fac2o))
  lm_8 <- stats::lm(y ~ x_fac2o * x_fac3t, data = df_8)

  expect_doppelganger("JAGS-formula-lm-8", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    hist(samples_8[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_8)["(Intercept)"], sd = summary(lm_8)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_8[,"mu_x_fac3t[1]"], freq = FALSE, main = "x_fac3t")
    curve(dnorm(x, mean = coef(lm_8)["x_fac3tB"], sd = summary(lm_8)$coefficients["x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_8[,"mu_x_fac3t[2]"], freq = FALSE, main = "x_fac3t")
    curve(dnorm(x, mean = coef(lm_8)["x_fac3tC"], sd = summary(lm_8)$coefficients["x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_8[,"mu_x_fac2o"], freq = FALSE, main = "x_fac2o")
    curve(dnorm(x, mean = coef(lm_8)["x_fac2o1"], sd = summary(lm_8)$coefficients["x_fac2o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_8[,"mu_x_fac2o__xXx__x_fac3t[1]"], freq = FALSE, main = "x_fac2o:fac3t")
    curve(dnorm(x, mean = coef(lm_8)["x_fac2o1:x_fac3tB"], sd = summary(lm_8)$coefficients["x_fac2o1:x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_8[,"mu_x_fac2o__xXx__x_fac3t[2]"], freq = FALSE, main = "x_fac2o:fac3t")
    curve(dnorm(x, mean = coef(lm_8)["x_fac2o1:x_fac3tC"], sd = summary(lm_8)$coefficients["x_fac2o1:x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)

  })
  # linear regression with an interaction between continous variable and orthornormal factor ----
  formula_9      <- JAGS_formula(~ x_cont1 * x_fac3o , parameter = "mu", data = df_all[,c("x_cont1", "x_fac3o")], prior_list = prior_list_all[c("intercept", "x_cont1", "x_fac3o", "x_cont1:x_fac3o")])
  prior_list_9   <- c(formula_9$prior_list, prior_list2)
  model_syntax_9 <- JAGS_add_priors(paste0("model{", formula_9$formula_syntax, model_syntax, "}"), prior_list_9)
  data_9         <- c(formula_9$data, N = nrow(df_all), y = list(df_all$y))

  model_9   <- rjags::jags.model(file = textConnection(model_syntax_9), inits = JAGS_get_inits(prior_list_9, chains = 2, seed = 1), data = data_9, n.chains = 2, quiet = TRUE)
  samples_9 <- rjags::coda.samples(model = model_9, variable.names = JAGS_to_monitor(prior_list_9), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_9 <- do.call(rbind, samples_9)

  df_9 <- df_all
  contrasts(df_9$x_fac3o) <- contr.orthonormal(levels(df_9$x_fac3o))
  lm_9 <- stats::lm(y ~ x_cont1 * x_fac3o, data = df_9)

  expect_doppelganger("JAGS-formula-lm-9", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    hist(samples_9[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_9)["(Intercept)"], sd = summary(lm_9)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_9[,"mu_x_fac3o[1]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_9)["x_fac3o1"], sd = summary(lm_9)$coefficients["x_fac3o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_9[,"mu_x_fac3o[2]"], freq = FALSE, main = "x_fac3o")
    curve(dnorm(x, mean = coef(lm_9)["x_fac3o2"], sd = summary(lm_9)$coefficients["x_fac3o2", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_9[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")
    curve(dnorm(x, mean = coef(lm_9)["x_cont1"], sd = summary(lm_9)$coefficients["x_cont1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_9[,"mu_x_cont1__xXx__x_fac3o[1]"], freq = FALSE, main = "x_cont1:x_fac3o")
    curve(dnorm(x, mean = coef(lm_9)["x_cont1:x_fac3o1"], sd = summary(lm_9)$coefficients["x_cont1:x_fac3o1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_9[,"mu_x_cont1__xXx__x_fac3o[2]"], freq = FALSE, main = "x_cont1:x_fac3o")
    curve(dnorm(x, mean = coef(lm_9)["x_cont1:x_fac3o2"], sd = summary(lm_9)$coefficients["x_cont1:x_fac3o2", "Std. Error"]), add = TRUE, lwd = 2)

  })


  # linear regression with an interaction between continous variable and orthornormal factor ----
  formula_10      <- JAGS_formula(~ x_cont1 * x_fac3t , parameter = "mu", data = df_all[,c("x_cont1", "x_fac3t")], prior_list = prior_list_all[c("intercept", "x_cont1", "x_fac3t", "x_cont1:x_fac3t")])
  prior_list_10   <- c(formula_10$prior_list, prior_list2)
  model_syntax_10 <- JAGS_add_priors(paste0("model{", formula_10$formula_syntax, model_syntax, "}"), prior_list_10)
  data_10         <- c(formula_10$data, N = nrow(df_all), y = list(df_all$y))

  model_10   <- rjags::jags.model(file = textConnection(model_syntax_10), inits = JAGS_get_inits(prior_list_10, chains = 2, seed = 1), data = data_10, n.chains = 2, quiet = TRUE)
  samples_10 <- rjags::coda.samples(model = model_10, variable.names = JAGS_to_monitor(prior_list_10), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_10 <- do.call(rbind, samples_10)

  lm_10 <- stats::lm(y ~ x_cont1 * x_fac3t, data = df_all)

  expect_doppelganger("JAGS-formula-lm-10", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    hist(samples_10[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_10)["(Intercept)"], sd = summary(lm_10)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_10[,"mu_x_fac3t[1]"], freq = FALSE, main = "x_fac3t")
    curve(dnorm(x, mean = coef(lm_10)["x_fac3tB"], sd = summary(lm_10)$coefficients["x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_10[,"mu_x_fac3t[2]"], freq = FALSE, main = "x_fac3t")
    curve(dnorm(x, mean = coef(lm_10)["x_fac3tC"], sd = summary(lm_10)$coefficients["x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_10[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")
    curve(dnorm(x, mean = coef(lm_10)["x_cont1"], sd = summary(lm_10)$coefficients["x_cont1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_10[,"mu_x_cont1__xXx__x_fac3t[1]"], freq = FALSE, main = "x_cont1:x_fac3t")
    curve(dnorm(x, mean = coef(lm_10)["x_cont1:x_fac3tB"], sd = summary(lm_10)$coefficients["x_cont1:x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_10[,"mu_x_cont1__xXx__x_fac3t[2]"], freq = FALSE, main = "x_cont1:x_fac3t")
    curve(dnorm(x, mean = coef(lm_10)["x_cont1:x_fac3tC"], sd = summary(lm_10)$coefficients["x_cont1:x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)

  })


  # scaling formula parameters by another parameter works ----
  prior_list_1s  <- prior_list_all[c("intercept", "x_cont1")]
  attr(prior_list_1s$x_cont1, "multiply_by") <- "sigma"
  formula_1s     <- JAGS_formula(~ x_cont1, parameter = "mu", data = df_all[,"x_cont1", drop = FALSE], prior_list = prior_list_1s)
  prior_list_1s  <- c(formula_1s$prior_list, prior_list2)
  model_syntax_1s<- JAGS_add_priors(paste0("model{", formula_1s$formula_syntax, model_syntax, "}"), prior_list_1s)
  data_1         <- c(formula_1$data, N = nrow(df_all), y = list(df_all$y))

  model_1s  <- rjags::jags.model(file = textConnection(model_syntax_1s), inits = JAGS_get_inits(prior_list_1s, chains = 2, seed = 1), data = data_1, n.chains = 2, quiet = TRUE)
  samples_1s <- rjags::coda.samples(model = model_1s, variable.names = JAGS_to_monitor(prior_list_1s), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_1s <- do.call(rbind, samples_1s)

  expect_equal(formula_1s$formula_syntax, "for(i in 1:N_mu){\n  mu[i] = mu_intercept + sigma * mu_x_cont1 * mu_data_x_cont1[i]\n}\n")

  lm_1s <- stats::lm(y ~ I(sd(y) * x_cont1), data = df_all)

  expect_doppelganger("JAGS-formula-lm-1s", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_1s[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_1s)["(Intercept)"], sd = summary(lm_1s)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_1s[,"mu_x_cont1"], freq = FALSE, main = "I(sd(y) * x_cont1)")
    curve(dnorm(x, mean = coef(lm_1s)["I(sd(y) * x_cont1)"], sd = summary(lm_1s)$coefficients["I(sd(y) * x_cont1)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_1s[,"sigma"], freq = FALSE, main = "sigma")
    abline(v = sigma(lm_1s), lwd = 3)
  })


  # input checks work
  expect_error(JAGS_formula(~ x_cont1 , parameter = "mu", data = df_all[,c("x_cont1"), drop = FALSE], prior_list = prior_list_all[c("x_cont1")]),
               "The 'intercept' objects are missing in the 'prior_list' argument.")
  expect_error(JAGS_formula(~ x_cont1 , parameter = "mu", data = df_all[,c("x_cont1"), drop = FALSE], prior_list = prior_list_all[c("intercept")]),
               "The 'x_cont1' objects are missing in the 'prior_list' argument.")
  expect_error(JAGS_formula(~ x_fac2t , parameter = "mu", data = df_all[,c("x_cont1"), drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac2t")]),
               "The 'x_fac2t' predictor variable is missing in the data set.")
  expect_error(JAGS_formula(~ x_fac2t , parameter = "mu", data = as.matrix(df_all), prior_list = prior_list_all[c("intercept", "x_fac2t")]),
               "'data' must be a data.frame")
  expect_error(JAGS_formula(~ x_fac2t , parameter = "mu", data = df_all, prior_list = list(
    "intercept" = prior("normal", list(0, 1)),
    "x_fac2t"   = prior("normal", list(0, 1))
  )), "Unsupported prior distribution defined for 'x_fac2t' factor variable")
  expect_error(JAGS_formula(~ x_cont1 , parameter = "mu", data = df_all, prior_list = list(
    "intercept" = prior("normal", list(0, 1)),
    "x_cont1"   = prior_factor("normal", list(0, 1), contrast = "treatment")
  )), "Unsupported prior distribution defined for 'x_cont1' continuous variable.")

})



