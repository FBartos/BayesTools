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
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3i = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3md= factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
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
    "x_fac3i"         = prior_factor("normal",  contrast = "independent", list(0, 1)),
    "x_fac3md"        = prior_factor("mnormal", contrast = "meandif",     list(0, 1)),
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-1", function(){

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

  vdiffr::expect_doppelganger("JAGS-formula-lm-2", function(){

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


  # linear regression with a treatment factor (2 levels) ----
  formula_3      <- JAGS_formula(~ x_fac2t, parameter = "mu", data = df_all[,"x_fac2t",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac2t")])
  prior_list_3   <- c(formula_3$prior_list, prior_list2)
  model_syntax_3 <- JAGS_add_priors(paste0("model{", formula_3$formula_syntax, model_syntax, "}"), prior_list_3)
  data_3         <- c(formula_3$data, N = nrow(df_all), y = list(df_all$y))

  model_3   <- rjags::jags.model(file = textConnection(model_syntax_3), inits = JAGS_get_inits(prior_list_3, chains = 2, seed = 1), data = data_3, n.chains = 2, quiet = TRUE)
  samples_3 <- rjags::coda.samples(model = model_3, variable.names = JAGS_to_monitor(prior_list_3), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_3 <- do.call(rbind, samples_3)

  lm_3 <- stats::lm(y ~ x_fac2t, data = df_all)

  vdiffr::expect_doppelganger("JAGS-formula-lm-3", function(){

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


  # linear regression with an orthonormal factor (2 levels) ----
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-4", function(){

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


  # linear regression with a treatment factor (3 levels) ----
  formula_5      <- JAGS_formula(~ x_fac3t, parameter = "mu", data = df_all[,"x_fac3t",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac3t")])
  prior_list_5   <- c(formula_5$prior_list, prior_list2)
  model_syntax_5 <- JAGS_add_priors(paste0("model{", formula_5$formula_syntax, model_syntax, "}"), prior_list_5)
  data_5         <- c(formula_5$data, N = nrow(df_all), y = list(df_all$y))

  model_5   <- rjags::jags.model(file = textConnection(model_syntax_5), inits = JAGS_get_inits(prior_list_5, chains = 2, seed = 1), data = data_5, n.chains = 2, quiet = TRUE)
  samples_5 <- rjags::coda.samples(model = model_5, variable.names = JAGS_to_monitor(prior_list_5), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_5 <- do.call(rbind, samples_5)

  lm_5 <- stats::lm(y ~ x_fac3t, data = df_all)

  vdiffr::expect_doppelganger("JAGS-formula-lm-5", function(){

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


  # linear regression with an orthonormal factor (3 levels) ----
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-6", function(){
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


  # linear regression with an orthonormal interaction between factors ----
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-7", function(){
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-8", function(){
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
  # linear regression with an interaction between continuous variable and orthonormal factor ----
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

  vdiffr::expect_doppelganger("JAGS-formula-lm-9", function(){
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


  # linear regression with an interaction between continuous variable and orthonormal factor ----
  formula_10      <- JAGS_formula(~ x_cont1 * x_fac3t , parameter = "mu", data = df_all[,c("x_cont1", "x_fac3t")], prior_list = prior_list_all[c("intercept", "x_cont1", "x_fac3t", "x_cont1:x_fac3t")])
  prior_list_10   <- c(formula_10$prior_list, prior_list2)
  model_syntax_10 <- JAGS_add_priors(paste0("model{", formula_10$formula_syntax, model_syntax, "}"), prior_list_10)
  data_10         <- c(formula_10$data, N = nrow(df_all), y = list(df_all$y))

  model_10   <- rjags::jags.model(file = textConnection(model_syntax_10), inits = JAGS_get_inits(prior_list_10, chains = 2, seed = 1), data = data_10, n.chains = 2, quiet = TRUE)
  samples_10 <- rjags::coda.samples(model = model_10, variable.names = JAGS_to_monitor(prior_list_10), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_10 <- do.call(rbind, samples_10)

  lm_10 <- stats::lm(y ~ x_cont1 * x_fac3t, data = df_all)

  vdiffr::expect_doppelganger("JAGS-formula-lm-10", function(){

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

  vdiffr::expect_doppelganger("JAGS-formula-lm-1s", function(){

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

  # linear regression with an independent factor (3 levels) ----
  formula_11      <- JAGS_formula(~ x_fac3i - 1, parameter = "mu", data = df_all[,"x_fac3i",drop = FALSE], prior_list = prior_list_all[c("x_fac3i")])
  prior_list_11   <- c(formula_11$prior_list, prior_list2)
  model_syntax_11 <- JAGS_add_priors(paste0("model{", formula_11$formula_syntax, model_syntax, "}"), prior_list_11)
  data_11         <- c(formula_11$data, N = nrow(df_all), y = list(df_all$y))

  model_11   <- rjags::jags.model(file = textConnection(model_syntax_11), inits = JAGS_get_inits(prior_list_11, chains = 2, seed = 1), data = data_11, n.chains = 2, quiet = TRUE)
  samples_11 <- rjags::coda.samples(model = model_11, variable.names = JAGS_to_monitor(prior_list_11), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_11 <- do.call(rbind, samples_11)

  lm_11 <- stats::lm(y ~ x_fac3i - 1, data = df_all)

  vdiffr::expect_doppelganger("JAGS-formula-lm-11", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_11[,"mu_x_fac3i[1]"], freq = FALSE, main = "x_fac3i[1]")
    curve(dnorm(x, mean = coef(lm_11)["x_fac3iA"], sd = summary(lm_11)$coefficients["x_fac3iA", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_11[,"mu_x_fac3i[2]"], freq = FALSE, main = "x_fac3i[2]")
    curve(dnorm(x, mean = coef(lm_11)["x_fac3iB"], sd = summary(lm_11)$coefficients["x_fac3iB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_11[,"mu_x_fac3i[3]"], freq = FALSE, main = "x_fac3i[3]")
    curve(dnorm(x, mean = coef(lm_11)["x_fac3iC"], sd = summary(lm_11)$coefficients["x_fac3iC", "Std. Error"]), add = TRUE, lwd = 2)
  })


  # linear regression with a meandif factor (3 levels) ----
  formula_12      <- JAGS_formula(~ x_fac3md, parameter = "mu", data = df_all[,"x_fac3md",drop = FALSE], prior_list = prior_list_all[c("intercept", "x_fac3md")])
  prior_list_12   <- c(formula_12$prior_list, prior_list2)
  model_syntax_12 <- JAGS_add_priors(paste0("model{", formula_12$formula_syntax, model_syntax, "}"), prior_list_12)
  data_12         <- c(formula_12$data, N = nrow(df_all), y = list(df_all$y))

  model_12   <- rjags::jags.model(file = textConnection(model_syntax_12), inits = JAGS_get_inits(prior_list_12, chains = 2, seed = 1), data = data_12, n.chains = 2, quiet = TRUE)
  samples_12 <- rjags::coda.samples(model = model_12, variable.names = JAGS_to_monitor(prior_list_12), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_12 <- do.call(rbind, samples_12)

  df_12 <- df_all
  contrasts(df_12$x_fac3md) <- contr.meandif(levels(df_12$x_fac3o))
  lm_12 <- stats::lm(y ~ x_fac3md, data = df_12)

  vdiffr::expect_doppelganger("JAGS-formula-lm-12", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_12[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_12)["(Intercept)"], sd = summary(lm_12)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_12[,"mu_x_fac3md[1]"], freq = FALSE, main = "x_fac3md")
    curve(dnorm(x, mean = coef(lm_12)["x_fac3md1"], sd = summary(lm_12)$coefficients["x_fac3md1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_12[,"mu_x_fac3md[2]"], freq = FALSE, main = "x_fac3md")
    curve(dnorm(x, mean = coef(lm_12)["x_fac3md2"], sd = summary(lm_12)$coefficients["x_fac3md2", "Std. Error"]), add = TRUE, lwd = 2)
  })


  # linear regression with a spike independent factor (3 levels) ----
  prior_list_13   <- list("x_fac3i" = prior_factor("spike",  contrast = "independent", list(1.5)))
  formula_13      <- JAGS_formula(~ x_fac3i - 1, parameter = "mu", data = df_all[,"x_fac3i",drop = FALSE], prior_list = prior_list_13)
  prior_list_13   <- c(formula_13$prior_list, prior_list2)
  model_syntax_13 <- JAGS_add_priors(paste0("model{", formula_13$formula_syntax, model_syntax, "}"), prior_list_13)
  data_13         <- c(formula_13$data, N = nrow(df_all), y = list(df_all$y))

  model_13   <- rjags::jags.model(file = textConnection(model_syntax_13), inits = JAGS_get_inits(prior_list_13, chains = 2, seed = 1), data = data_13, n.chains = 2, quiet = TRUE)
  samples_13 <- rjags::coda.samples(model = model_13, variable.names = JAGS_to_monitor(prior_list_13), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_13 <- do.call(rbind, samples_13)
  expect_equal(diag(3), contr.independent(1:3))

  vdiffr::expect_doppelganger("JAGS-formula-lm-13", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_13[,"mu_x_fac3i[1]"], freq = FALSE, main = "x_fac3i[1]")
    hist(samples_13[,"mu_x_fac3i[2]"], freq = FALSE, main = "x_fac3i[2]")
    hist(samples_13[,"mu_x_fac3i[3]"], freq = FALSE, main = "x_fac3i[3]")
  })


  # linear regression with a meandif spike factor (3 levels) ----
  prior_list_14   <- list("intercept" = prior_list_all$intercept, "x_fac3md" = prior_factor("spike",  contrast = "meandif", list(0)))
  formula_14      <- JAGS_formula(~ x_fac3md, parameter = "mu", data = df_all[,"x_fac3md",drop = FALSE], prior_list = prior_list_14)
  prior_list_14   <- c(formula_14$prior_list, prior_list2)
  model_syntax_14 <- JAGS_add_priors(paste0("model{", formula_14$formula_syntax, model_syntax, "}"), prior_list_14)
  data_14         <- c(formula_14$data, N = nrow(df_all), y = list(df_all$y))

  model_14   <- rjags::jags.model(file = textConnection(model_syntax_14), inits = JAGS_get_inits(prior_list_14, chains = 2, seed = 1), data = data_14, n.chains = 2, quiet = TRUE)
  samples_14 <- rjags::coda.samples(model = model_14, variable.names = JAGS_to_monitor(prior_list_14), n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples_14 <- do.call(rbind, samples_14)

  df_14 <- df_all
  contrasts(df_14$x_fac3md) <- contr.meandif(levels(df_14$x_fac3o))
  lm_14 <- stats::lm(y ~ 1, data = df_14)

  vdiffr::expect_doppelganger("JAGS-formula-lm-14", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples_14[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_14)["(Intercept)"], sd = summary(lm_14)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(samples_14[,"mu_x_fac3md[1]"], freq = FALSE, main = "x_fac3md")
    hist(samples_14[,"mu_x_fac3md[2]"], freq = FALSE, main = "x_fac3md")
  })

})

test_that("JAGS formula tools work", {

  # additional tools work
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont")),
    c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t:x_cont")
  )
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"), formula_parameters = "mu"),
    c("(mu) x_cont", "(mu) x_fac3t", "(mu) x_fac3t:x_cont")
  )
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"), formula_parameters = "mu", formula_prefix = FALSE),
    c("x_cont", "x_fac3t", "x_fac3t:x_cont")
  )

  expect_equal(
    JAGS_parameter_names(c("x_cont", "x_fac3t", "x_fac3t:x_cont")),
    c("x_cont", "x_fac3t", "x_fac3t__xXx__x_cont")
  )
  expect_equal(
    JAGS_parameter_names(c("x_cont", "x_fac3t", "x_fac3t:x_cont"), formula_parameter = "mu"),
    c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont")
  )

})

test_that("JAGS evaluate formula works", {

  # check the posterior distributions with weak priors against a maximum likelihood estimates with ML
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  # complex formula including scaling
  set.seed(1)
  df_all <- data.frame(
    x_cont1 = rnorm(60),
    x_cont2 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(60, 0.1, 0.5) + 0.30 * df_all$x_cont1 - 0.15 * df_all$x_cont1 * df_all$x_cont2 +
    ifelse(df_all$x_fac3o == "A", 0.2, ifelse(df_all$x_fac3o == "B", -0.2, 0))

  prior_list_all <- list(
    "intercept"       = prior("normal", list(0, 5)),
    "x_cont1"         = prior("normal", list(0, 1)),
    "x_cont2"         = prior("normal", list(0, 1)),
    "x_fac2t"         = prior_factor("normal",  contrast = "treatment",   list(0, 1)),
    "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
    "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0,  2))
  )
  prior_list2  <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 1))
  )
  model_syntax <- paste0(
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n"
  )

  formula      <- JAGS_formula(~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, parameter = "mu", data = df_all, prior_list = prior_list_all)
  prior_list   <- c(formula$prior_list, prior_list2)
  model_syntax <- JAGS_add_priors(paste0("model{", formula$formula_syntax, model_syntax, "}"), prior_list)
  data         <- c(formula$data, N = nrow(df_all), y = list(df_all$y))

  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = JAGS_get_inits(prior_list, chains = 1, seed = 1), data = data, n.chains = 1, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = JAGS_to_monitor(prior_list), n.iter = 10, quiet = TRUE, progress.bar = "none")

  new_data <-  data.frame(
    x_cont1 = c(0, 0, 1, 1),
    x_cont2 = c(0, 1, 0, 1),
    x_fac2t = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
    x_fac3o = factor(c("A", "B", "C", "A"), levels = c("A", "B", "C"))
  )

  # test the results against the lm function (by passing the ML estimates)
  contrasts(df_all$x_fac3o) <- contr.orthonormal(levels(df_all$x_fac3o))
  fit_lm <- stats::lm(y~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, data = df_all)

  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 1)
  colnames(samples_new) <- colnames(samples[[1]])
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  expect_equal(predict(fit_lm, newdata = new_data), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)[,1])

  # for a posterior samples matrix
  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 5, ncol = length(samples_new), byrow = TRUE)
  colnames(samples_new) <- colnames(samples[[1]])
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  expect_equal(matrix(predict(fit_lm, newdata = new_data), nrow = 4, ncol = 5), unname(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)))

  # check filling in missing or miss ordered factor levels
  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 1)
  colnames(samples_new) <- colnames(samples[[1]])
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  new_data2         <- new_data
  new_data2$x_fac2t <- factor(as.character(new_data2$x_fac2t), levels = c("B", "A"))
  expect_equal(predict(fit_lm, newdata = new_data), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data2, prior_list)[,1])

  new_data3         <- new_data
  new_data3$x_fac3o <- factor(c("A", "B", "A", "B"), levels = c("B", "A"))
  expect_equal(predict(fit_lm, newdata = new_data3), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data3, prior_list)[,1])

  new_data4         <- new_data
  new_data4$x_fac3o <- c("A", "B", "A", "B")
  expect_equal(predict(fit_lm, newdata = new_data3), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data4, prior_list)[,1])

  # check scaling works (by multiplying be zero)
  prior_list2 <- prior_list
  attr(prior_list2$mu_x_cont2, "multiply_by") <- 0
  attr(prior_list2$mu_x_fac2t, "multiply_by") <- 0

  samples_new2 <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new2 <- matrix(samples_new2, nrow = 1)
  colnames(samples_new2) <- colnames(samples[[1]])
  samples_new2[,"mu_x_cont2"]    <- 0
  samples_new2[,"mu_x_fac2t"] <- 0
  samples_new2 <- coda::as.mcmc.list(coda::as.mcmc(samples_new2))

  expect_equal(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list2)[,1],
               JAGS_evaluate_formula(samples_new2, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)[,1])

  # check scaling by another parameter works
  prior_list2 <- prior_list
  attr(prior_list2$mu_x_cont2, "multiply_by") <- "sigma"

  expect_equal(unname(unlist(JAGS_evaluate_formula(samples, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list2)[,1])),
               c(0.4436353, -0.0658681, 0.1870391, 0.8548012), tolerance = 1e-5)

  ### test input tests
  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data[,1:3], prior_list),
               "The 'x_fac3o' predictor variable is missing in the data.")
  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list[-1]),
               "The prior distribution for the 'x_fac2t' term is missing in the prior_list.")

  bad_data         <- new_data
  bad_data$x_fac2t <- factor(c("C", "B", "C", "B"), levels = c("B", "C"))

  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", bad_data, prior_list),
               "Levels specified in the 'x_fac2t' factor variable do not match the levels used for model specification.")

  bad_data2         <- new_data
  bad_data2$x_fac2t <- c("C", "B", "C", "B")

  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", bad_data2, prior_list),
               "Levels specified in the 'x_fac2t' factor variable do not match the levels used for model specification.")

  # evaluate formula with spike prior distributions ----
  set.seed(1)
  df_all <- data.frame(
    x_fac2i  = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o  = factor(sample(c("A", "B", "C"), 60, replace = TRUE), levels = c("A", "B", "C")),
    x_fac3t  = factor(sample(c("A", "B", "C"), 60, replace = TRUE), levels = c("A", "B", "C")),
    x_fac3md = factor(sample(c("A", "B", "C"), 60, replace = TRUE), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(60, 0.1, 0.5)

  prior_list_all <- list(
    "intercept" = prior("normal", list(0, 5)),
    "x_fac2i"   = prior_factor("spike", contrast = "independent", list(1)),
    "x_fac3o"   = prior_factor("spike", contrast = "orthonormal", list(0)),
    "x_fac3t"   = prior_factor("spike", contrast = "treatment",   list(2)),
    "x_fac3md"  = prior_factor("spike", contrast = "meandif",     list(0))
  )
  prior_list2  <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 1))
  )
  model_syntax <- paste0(
    "model{",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit1 <- JAGS_fit(
    model_syntax       = model_syntax,
    formula_list       = list(mu = ~ x_fac2i + x_fac3o + x_fac3t + x_fac3md),
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list2,
    formula_data_list  = list(mu = df_all),
    formula_prior_list = list(mu = prior_list_all))

  new_data <-  data.frame(
    x_fac2i  = factor(c("A", "B", "A"), levels = c("A", "B")),
    x_fac3o  = factor(c("A", "A", "B"), levels = c("A", "B", "C")),
    x_fac3t  = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    x_fac3md = factor(c("B", "B", "C"), levels = c("A", "B", "C"))
  )
  new_samples <- JAGS_evaluate_formula(fit1, ~ x_fac2i + x_fac3o + x_fac3t + x_fac3md, "mu", new_data, attr(fit1, "prior_list"))
  new_samples <- apply(new_samples, 1, mean)

  intercept_estimate <- JAGS_estimates_table(fit1)["(mu) intercept", "Mean"]

  expect_equivalent(intercept_estimate + 1, new_samples[1])
  expect_equivalent(intercept_estimate + 1 + 2, new_samples[2])
  expect_equivalent(intercept_estimate + 1 + 2, new_samples[3])


  # dealing with spike and slab and mixture priors
  prior_list_all2 <- list(
    "intercept" = prior_spike_and_slab(prior("normal", list(0, 5))),
    "x_fac2i"   = prior_mixture(list(
      prior("spike", list(1)),
      prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    ), is_null = c(T,  F)),
    "x_fac3o"   = prior_spike_and_slab(prior_factor("mnormal", contrast = "orthonormal", list(0, 1))),
    "x_fac3t"   = prior_mixture(list(
      prior_factor("normal", contrast = "treatment",   list(0, 1)),
      prior("spike", list(0))
    ), is_null = c(T,  F))
  )
  fit2 <- JAGS_fit(
    model_syntax       = model_syntax,
    formula_list       = list(mu = ~ x_fac2i + x_fac3o + x_fac3t),
    data               = list(y = df_all$y, N = nrow(df_all)),
    prior_list         = prior_list2,
    formula_data_list  = list(mu = df_all),
    formula_prior_list = list(mu = prior_list_all2), chains = 1, adapt = 100, burnin = 100, sample = 200)

  new_samples <- JAGS_evaluate_formula(fit2, ~ x_fac2i + x_fac3o + x_fac3t, "mu", new_data, attr(fit2, "prior_list"))
  expect_equivalent(dim(new_samples), c(3, 200))

})

test_that("Expression handling functions works", {

  f1 <- formula(y ~ 1)
  f2 <- formula(y ~ z)
  f3 <- formula(y ~ expression(x))
  f4 <- formula(y ~ z + expression(x))
  f5 <- formula(y ~ expression(x) + z)
  f6 <- formula(y ~ expression(x) + z + expression(b))

  expect_true(!.has_expression(f1))
  expect_true(!.has_expression(f2))
  expect_true(.has_expression(f3))
  expect_true(.has_expression(f4))
  expect_true(.has_expression(f5))
  expect_true(.has_expression(f6))

  expect_equal(.extract_expressions(f3), list(expression(x)))
  expect_equal(.extract_expressions(f4), list(expression(x)))
  expect_equal(.extract_expressions(f5), list(expression(x)))
  expect_equal(.extract_expressions(f6), list(expression(x), expression(b)))

  expect_equal(.remove_expressions(f1), formula(y ~ 1))
  expect_equal(.remove_expressions(f2), formula(y ~ z))
  expect_equal(.remove_expressions(f3), formula(y ~ 1))
  expect_equal(.remove_expressions(f4), formula(y ~ z))
  expect_equal(.remove_expressions(f5), formula(y ~ z))
  expect_equal(.remove_expressions(f6), formula(y ~ z))
})
