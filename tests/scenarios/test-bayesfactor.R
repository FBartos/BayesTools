if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("bayesfactor")
# testthat::test_file("tests/scenarios/test-bayesfactor.R")

### Description
# Bayes factors of the default (JZS) models of the BayesFactor package,
# reproduced with BayesTools priors: one- and two-sample t-tests (sleep), a
# two-factor ANOVA (ToothGrowth), a random factor (the sleep subjects), a
# simple regression (attitude), and joint-multinomial contingency tables
# (raceDolls, mtcars). BayesFactor puts a flat prior on the intercept, the
# Jeffreys prior 1 / sigma on the residual SD, Cauchy priors on standardized
# effects (scale sqrt(2) / 2 for t-tests; 1 / 2 for fixed and 1 for random
# factors in orthonormal sum-to-zero coding; sqrt(2) / 4 for a predictor
# standardized by its population SD), and Dirichlet(1) priors on cell
# probabilities. BayesTools uses a uniform(-100, 100) intercept (standardized
# outcomes) and a Gamma(0.001, 0.001) SD truncated to [0.01, 100]; the compared
# models share both, so their normalizing constants cancel and the Bayes
# factors approach the improper limit. Effects get orthonormal multivariate
# Cauchy factor priors and Cauchy slope priors multiplied by sigma
# (multiply_by), and cell probabilities Dirichlet priors. References:
# BayesFactor's exact t-test, single-factor, regression, and contingency Bayes
# factors; a deterministic quadrature of the JZS Bayes factor, which
# reproduces those exact values and replaces BayesFactor's Monte Carlo
# estimates for several factors; and, for the one-sample t-test, a quadrature
# under the BayesTools priors, which measures the effect of the proper SD
# prior. Not covered: proportionBF (a logistic prior on the log odds, which
# BayesTools does not offer), correlationBF, and regressionBF with several
# predictors (a Zellner-Siow prior joint across continuous terms).

testthat::test_that("BayesFactor default Bayes factors", {

  sigma_prior <- prior("gamma", list(shape = 0.001, rate = 0.001), truncation = list(lower = 0.01, upper = 100))
  flat_prior  <- prior("uniform", list(a = -100, b = 100))
  factor_prior <- function(r) {
    out <- prior_factor("mcauchy", list(location = 0, scale = r), contrast = "orthonormal")
    attr(out, "multiply_by") <- "sigma"
    out
  }
  normal_likelihood <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu[i], pow(sigma, -2))
    }
  }
  "
  normal_log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["y"]], parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }
  marglik <- function(fit, data, log_posterior = normal_log_posterior, prior_list = list(sigma = sigma_prior)) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data, prior_list = prior_list, seed = 1)
  }
  log_bf <- function(fit_1, fit_0) scenario_log_bf_bridge(fit_1[["marglik"]], fit_0[["marglik"]])
  bayesfactor_log_bf <- function(bf) BayesFactor::extractBF(bf, logbf = TRUE)[["bf"]]

  # The JZS log Bayes factor against the intercept-only model for a balanced
  # design whose centred term columns are mutually orthogonal with a common
  # squared norm lambda_t per term (orthonormal factor coding, a predictor
  # standardized by its population SD); term t has p_t columns with
  # coefficients ~ N(0, g_t sigma^2 I) and g_t ~ inverse-gamma(1/2, r_t^2 / 2).
  # Given g, BF10 = prod_t (1 + g_t lambda_t)^(-p_t / 2) (1 - w)^(-k) with
  # w = sum_t a_t g_t lambda_t / (1 + g_t lambda_t), a_t the share of the
  # centred sum of squares that term t explains, and k = (n - 1) / 2. Writing
  # (1 - w)^(-k) = Gamma(k)^-1 int z^(k - 1) exp(-z (1 - w)) dz separates the
  # terms given z, so the Bayes factor is a one-dimensional integral over z of
  # a product of one-dimensional integrals over each g_t; both are trapezoid
  # sums in log g and log z (halving the step changes the result below 1e-8)
  jzs_log_bf <- function(y, terms, rscale, step = 0.01) {
    n  <- length(y)
    yc <- y - mean(y)
    Xc <- scale(do.call(cbind, terms), scale = FALSE)
    XtX <- crossprod(Xc)
    index  <- rep(seq_along(terms), vapply(terms, ncol, integer(1)))
    column <- diag(XtX)
    lambda <- tapply(column, index, mean)
    if (any(abs(XtX[upper.tri(XtX)]) > 1e-8 * max(column)) || any(abs(column - lambda[index]) > 1e-8 * max(column))) {
      stop("The quadrature requires orthogonal term columns with a common norm per term.", call. = FALSE)
    }
    p <- tabulate(index)
    a <- tapply(drop(crossprod(Xc, yc))^2, index, sum) / lambda / sum(yc^2)
    k <- (n - 1) / 2
    softplus  <- function(x) ifelse(x > 30, x, log1p(exp(x)))
    logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))
    u <- seq(-30, 40, by = step)
    v <- seq(-10, log(10 * k / (1 - sum(a))) + 3, by = step)
    log_F <- vapply(seq_along(terms), function(t) {
      x    <- u + log(lambda[[t]])
      base <- -p[[t]] / 2 * softplus(x) + 0.5 * log(rscale[[t]]^2 / 2) - lgamma(0.5) - 0.5 * u -
        rscale[[t]]^2 / 2 * exp(-u) + log(step)
      h    <- stats::plogis(x)
      vapply(exp(v), function(z) logsumexp(base + z * a[[t]] * h), numeric(1))
    }, numeric(length(v)))
    -lgamma(k) + logsumexp(k * v - exp(v) + rowSums(matrix(log_F, nrow = length(v))) + log(step))
  }

  ### t-tests: sleep ----
  # one-sample: the paired differences; two-sample: the groups as independent
  # samples (the classic BayesFactor examples)
  d  <- as.numeric(with(sleep, extra[group == "2"] - extra[group == "1"]))
  sl <- sleep
  sl$y <- as.numeric(scale(sl$extra))
  one_sample_data <- list(d = d, N = length(d))
  one_sample_log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["d"]], parameters[["sigma"]] * parameters[["delta"]], parameters[["sigma"]], log = TRUE))
  }
  one_sample_null_log_posterior <- function(parameters, data) {
    sum(stats::dnorm(data[["d"]], 0, parameters[["sigma"]], log = TRUE))
  }
  one_sample_prior <- list(delta = prior("cauchy", list(location = 0, scale = sqrt(2) / 2)), sigma = sigma_prior)

  one_sample <- list()
  one_sample[["alternative"]] <- scenario_fit("fit_one_sample", {
    fit <- JAGS_fit("model{ for(i in 1:N){ d[i] ~ dnorm(sigma * delta, pow(sigma, -2)) } }", one_sample_data,
                    prior_list = one_sample_prior, chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, one_sample_data, one_sample_log_posterior, one_sample_prior))
  })
  one_sample[["null"]] <- scenario_fit("fit_one_sample_null", {
    fit <- JAGS_fit("model{ for(i in 1:N){ d[i] ~ dnorm(0, pow(sigma, -2)) } }", one_sample_data,
                    prior_list = list(sigma = sigma_prior), chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, one_sample_data, one_sample_null_log_posterior))
  })
  two_sample <- list()
  two_sample[["alternative"]] <- scenario_fit("fit_two_sample", {
    data_list <- list(y = sl$y, N = nrow(sl))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ group), formula_data_list = list(mu = sl),
                    formula_prior_list = list(mu = list(intercept = flat_prior, group = factor_prior(1 / 2))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  two_sample[["null"]] <- scenario_fit("fit_two_sample_null", {
    data_list <- list(y = sl$y, N = nrow(sl))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ 1), formula_data_list = list(mu = sl),
                    formula_prior_list = list(mu = list(intercept = flat_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })

  # one-sample test under the BayesTools priors: delta = r tan(phi) with phi
  # uniform on (-pi / 2, pi / 2) under the Cauchy prior, and log sigma on a
  # grid over the truncation range; trapezoid sums on both grids
  exact_one_sample <- scenario_fit("exact_one_sample", {
    n <- length(d)
    log_sigma <- seq(log(0.01), log(100), by = 0.002)
    sigma     <- exp(log_sigma)
    log_prior <- stats::dgamma(sigma, 0.001, 0.001, log = TRUE) + log_sigma -
      log(stats::pgamma(100, 0.001, 0.001) - stats::pgamma(0.01, 0.001, 0.001)) + log(0.002)
    phi   <- seq(-pi / 2, pi / 2, length.out = 8001)[2:8000]
    delta <- sqrt(2) / 2 * tan(phi)
    log_lik <- function(sigma, delta) {
      -n / 2 * log(2 * pi) - n * log(sigma) - (sum(d^2) - 2 * sigma * delta * sum(d) + n * sigma^2 * delta^2) / (2 * sigma^2)
    }
    logsumexp <- function(x) max(x) + log(sum(exp(x - max(x))))
    rows <- t(vapply(sigma, function(s) {
      l <- log_lik(s, delta) + log(phi[2] - phi[1]) - log(pi)
      c(all = logsumexp(l), positive = logsumexp(l[delta > 0]))
    }, numeric(2))) + log_prior
    log_ml_1 <- logsumexp(rows[, "all"])
    log_ml_0 <- logsumexp(log_lik(sigma, 0) + log_prior)
    c(two_sided = log_ml_1 - log_ml_0,
      # BF(delta > 0 vs 0) = BF10 P(delta > 0 | y) / P(delta > 0)
      positive  = logsumexp(rows[, "positive"]) - log_ml_0 + log(2))
  })

  # for the paired differences delta = 0 lies far in the posterior tail, where
  # the kernel density estimate of the Savage-Dickey ordinate is raised by
  # smoothing, so the one-sample Savage-Dickey Bayes factors fall below the
  # exact ones; the one-sided Bayes factor uses the same ordinate
  scenario_text("bayes-factor-t-tests", {
    one_sided <- BayesFactor::ttestBF(d, nullInterval = c(0, Inf))
    rbind(
      scenario_bf_table("delta != 0", "delta = 0", c(
        "BayesFactor (exact)"         = bayesfactor_log_bf(BayesFactor::ttestBF(d)),
        "quadrature (BayesTools priors)" = exact_one_sample[["two_sided"]],
        "bridge sampling"             = log_bf(one_sample[["alternative"]], one_sample[["null"]]),
        "Savage-Dickey"               = log(hypothesis_BF(one_sample[["alternative"]][["fit"]], "delta = 0", seed = 1)[["BF"]])),
        note = "one-sample (paired differences)"),
      scenario_bf_table("delta > 0", "delta = 0", c(
        "BayesFactor (exact)"         = bayesfactor_log_bf(one_sided[1]),
        "quadrature (BayesTools priors)" = exact_one_sample[["positive"]],
        "hypothesis_BF"               = log(hypothesis_BF(one_sample[["alternative"]][["fit"]], "delta > 0 vs delta = 0", seed = 1)[["BF"]])),
        note = "one-sample, one-sided"),
      scenario_bf_table("group effect", "no group effect", c(
        "BayesFactor (exact)"         = bayesfactor_log_bf(BayesFactor::ttestBF(formula = extra ~ group, data = sleep)),
        "JZS quadrature"              = jzs_log_bf(sl$y, list(contr.orthonormal(levels(sl$group))[as.integer(sl$group), , drop = FALSE]), 1 / 2),
        "bridge sampling"             = log_bf(two_sample[["alternative"]], two_sample[["null"]]),
        "Savage-Dickey"               = log(hypothesis_BF(two_sample[["alternative"]][["fit"]], "group[1] = 0", seed = 1)[["BF"]])),
        note = "two-sample (independent groups)")
    )
  })

  ### two-factor ANOVA: ToothGrowth ----
  tg <- ToothGrowth
  tg$dose_f <- factor(tg$dose)
  tg$y      <- as.numeric(scale(tg$len))
  supp_columns <- contr.orthonormal(levels(tg$supp))[as.integer(tg$supp), , drop = FALSE]
  dose_columns <- contr.orthonormal(levels(tg$dose_f))[as.integer(tg$dose_f), , drop = FALSE]
  interaction_columns <- supp_columns[, 1L] * dose_columns
  anova <- list()
  anova[["null"]] <- scenario_fit("fit_anova_null", {
    data_list <- list(y = tg$y, N = nrow(tg))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ 1), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = flat_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  anova[["supp"]] <- scenario_fit("fit_anova_supp", {
    data_list <- list(y = tg$y, N = nrow(tg))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ supp), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = flat_prior, supp = factor_prior(1 / 2))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  anova[["dose"]] <- scenario_fit("fit_anova_dose", {
    data_list <- list(y = tg$y, N = nrow(tg))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ dose_f), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = flat_prior, dose_f = factor_prior(1 / 2))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  anova[["add"]] <- scenario_fit("fit_anova_add", {
    data_list <- list(y = tg$y, N = nrow(tg))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ supp + dose_f), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = flat_prior, supp = factor_prior(1 / 2),
                                                        dose_f = factor_prior(1 / 2))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  anova[["full"]] <- scenario_fit("fit_anova_full", {
    data_list <- list(y = tg$y, N = nrow(tg))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ supp * dose_f), formula_data_list = list(mu = tg),
                    formula_prior_list = list(mu = list(intercept = flat_prior, supp = factor_prior(1 / 2),
                                                        dose_f = factor_prior(1 / 2), "supp:dose_f" = factor_prior(1 / 2))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  # spike-and-slab on the interaction: one fit includes or excludes it; the
  # multiplier belongs to the spike-and-slab prior (a multiply_by attribute of
  # its slab is ignored)
  fit_anova_spike <- scenario_fit("fit_anova_spike", {
    interaction_prior <- prior_spike_and_slab(prior_factor("mcauchy", list(location = 0, scale = 1 / 2), contrast = "orthonormal"))
    attr(interaction_prior, "multiply_by") <- "sigma"
    JAGS_fit(normal_likelihood, list(y = tg$y, N = nrow(tg)), prior_list = list(sigma = sigma_prior),
             formula_list = list(mu = ~ supp * dose_f), formula_data_list = list(mu = tg),
             formula_prior_list = list(mu = list(intercept = flat_prior, supp = factor_prior(1 / 2), dose_f = factor_prior(1 / 2),
                                                 "supp:dose_f" = interaction_prior)),
             chains = 4, adapt = 1000, burnin = 1000, sample = 20000, seed = 1)
  })
  anova_terms <- list(
    supp = list(supp_columns),
    dose = list(dose_columns),
    add  = list(supp_columns, dose_columns),
    full = list(supp_columns, dose_columns, interaction_columns)
  )
  # the BayesFactor names of the four models
  anova_models <- c(supp = "supp", dose = "dose_f", add = "supp + dose_f", full = "supp + dose_f + supp:dose_f")
  anova_exact <- scenario_fit("exact_anova", {
    c(null = 0, vapply(anova_terms, function(terms) jzs_log_bf(tg$y, terms, rep(1 / 2, length(terms))), numeric(1)))
  })
  anova_bayesfactor <- scenario_fit("bayesfactor_anova", {
    set.seed(1)
    bf <- BayesFactor::anovaBF(len ~ supp * dose_f, data = tg, progress = FALSE)
    stats::setNames(BayesFactor::extractBF(bf, logbf = TRUE)[, c("bf", "error")], c("log_bf", "error"))
  })

  scenario_text("bayes-factor-anova", {
    tables <- lapply(names(anova_models), function(model) {
      scenario_bf_table(anova_models[[model]], "intercept only", c(
        "JZS quadrature"        = anova_exact[[model]],
        "BayesFactor (MC)"      = anova_bayesfactor[anova_models[[model]], "log_bf"],
        "bridge sampling"       = log_bf(anova[[model]], anova[["null"]])),
        note = sprintf("BayesFactor error %.2f%%", 100 * anova_bayesfactor[anova_models[[model]], "error"]))
    })
    interaction <- scenario_bf_table("supp:dose_f", "no interaction", c(
      "JZS quadrature"   = anova_exact[["full"]] - anova_exact[["add"]],
      "BayesFactor (MC)" = anova_bayesfactor[anova_models[["full"]], "log_bf"] - anova_bayesfactor[anova_models[["add"]], "log_bf"],
      "bridge sampling"  = log_bf(anova[["full"]], anova[["add"]]),
      "spike and slab"   = log(as.numeric(JAGS_inference_table(fit_anova_spike)[["inclusion_BF"]]))),
      note = "with supp + dose_f")
    do.call(rbind, c(tables, list(interaction)))
  })
  scenario_text("model-probabilities-anova", {
    nested <- c("null", "supp", "dose", "add", "full")
    inference <- models_inference(lapply(nested, function(name) {
      list(model = anova[[name]][["fit"]], marglik = anova[[name]][["marglik"]], prior_weights = 1)
    }))
    probability <- function(log_bf) round(exp(log_bf - max(log_bf)) / sum(exp(log_bf - max(log_bf))), 4)
    data.frame(
      model                = nested,
      "JZS quadrature"     = probability(anova_exact[nested]),
      "BayesFactor (MC)"   = probability(c(0, anova_bayesfactor[anova_models, "log_bf"])),
      "BayesTools"         = round(vapply(inference, function(model) model[["inference"]][["post_prob"]], numeric(1)), 4),
      check.names = FALSE
    )
  })

  ### random factor: sleep subjects ----
  # BayesFactor's random factor gets a g prior on all its levels; with a flat
  # intercept only the sum-to-zero part is identified, which is an
  # orthonormal factor with prior scale 1
  random <- list()
  random[["alternative"]] <- scenario_fit("fit_random_group", {
    data_list <- list(y = sl$y, N = nrow(sl))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ group + ID), formula_data_list = list(mu = sl),
                    formula_prior_list = list(mu = list(intercept = flat_prior, group = factor_prior(1 / 2), ID = factor_prior(1))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  random[["null"]] <- scenario_fit("fit_random_null", {
    data_list <- list(y = sl$y, N = nrow(sl))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ ID), formula_data_list = list(mu = sl),
                    formula_prior_list = list(mu = list(intercept = flat_prior, ID = factor_prior(1))),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  scenario_text("bayes-factor-random-factor", {
    group_columns <- contr.orthonormal(levels(sl$group))[as.integer(sl$group), , drop = FALSE]
    id_columns    <- contr.orthonormal(levels(sl$ID))[as.integer(sl$ID), , drop = FALSE]
    set.seed(1)
    bf <- BayesFactor::lmBF(extra ~ group + ID, data = sleep, whichRandom = "ID", progress = FALSE) /
      BayesFactor::lmBF(extra ~ ID, data = sleep, whichRandom = "ID", progress = FALSE)
    scenario_bf_table("group + ID", "ID", c(
      "JZS quadrature"   = jzs_log_bf(sl$y, list(group_columns, id_columns), c(1 / 2, 1)) - jzs_log_bf(sl$y, list(id_columns), 1),
      "BayesFactor (MC)" = bayesfactor_log_bf(bf),
      "bridge sampling"  = log_bf(random[["alternative"]], random[["null"]])),
      note = sprintf("BayesFactor error %.2f%%", 100 * BayesFactor::extractBF(bf)[["error"]]))
  })

  ### simple regression: attitude ----
  at <- attitude
  at$y <- as.numeric(scale(at$rating))
  at$complaints_s <- (at$complaints - mean(at$complaints)) / sqrt(mean((at$complaints - mean(at$complaints))^2))
  regression <- list()
  regression[["alternative"]] <- scenario_fit("fit_regression", {
    data_list <- list(y = at$y, N = nrow(at))
    slope_prior <- prior("cauchy", list(location = 0, scale = sqrt(2) / 4))
    attr(slope_prior, "multiply_by") <- "sigma"
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ complaints_s), formula_data_list = list(mu = at),
                    formula_prior_list = list(mu = list(intercept = flat_prior, complaints_s = slope_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  regression[["null"]] <- scenario_fit("fit_regression_null", {
    data_list <- list(y = at$y, N = nrow(at))
    fit <- JAGS_fit(normal_likelihood, data_list, prior_list = list(sigma = sigma_prior),
                    formula_list = list(mu = ~ 1), formula_data_list = list(mu = at),
                    formula_prior_list = list(mu = list(intercept = flat_prior)),
                    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
    list(fit = fit, marglik = marglik(fit, data_list))
  })
  scenario_text("bayes-factor-regression", {
    scenario_bf_table("complaints", "intercept only", c(
      "BayesFactor (exact)" = bayesfactor_log_bf(BayesFactor::lmBF(rating ~ complaints, data = attitude, progress = FALSE)),
      "JZS quadrature"      = jzs_log_bf(at$y, list(cbind(at$complaints_s)), sqrt(2) / 4),
      "bridge sampling"     = log_bf(regression[["alternative"]], regression[["null"]])))
  })

  ### contingency tables: joint multinomial ----
  # H1: the cell probabilities ~ Dirichlet(1); H0: independence with row and
  # column probabilities ~ Dirichlet(1) (Gunel and Dickey with a = 1); the
  # marginal likelihoods are Dirichlet-multinomial in closed form
  data("raceDolls", package = "BayesFactor")
  tables <- list(raceDolls = as.matrix(raceDolls), "mtcars cyl x am" = unclass(table(mtcars$cyl, mtcars$am)))
  multinomial_log_posterior <- function(parameters, data) {
    probabilities <- if (is.null(parameters[["p"]])) c(outer(parameters[["r"]], parameters[["cc"]])) else parameters[["p"]]
    stats::dmultinom(data[["y"]], prob = probabilities, log = TRUE)
  }
  contingency <- lapply(stats::setNames(names(tables), names(tables)), function(name) {
    scenario_fit(paste0("fit_contingency_", gsub("[^a-z]+", "_", tolower(name))), {
      table     <- tables[[name]]
      data_list <- list(y = c(table), K = length(table), N = sum(table), I = nrow(table), J = ncol(table))
      prior_1   <- list(p = prior("dirichlet", list(alpha = rep(1, length(table)))))
      prior_0   <- list(r = prior("dirichlet", list(alpha = rep(1, nrow(table)))),
                        cc = prior("dirichlet", list(alpha = rep(1, ncol(table)))))
      fit_1 <- JAGS_fit("model{ y[1:K] ~ dmulti(p[1:K], N) }", data_list[c("y", "K", "N")], prior_list = prior_1,
                        chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
      fit_0 <- JAGS_fit("model{
        for(i in 1:I){ for(j in 1:J){ q[(j - 1) * I + i] <- r[i] * cc[j] } }
        y[1:K] ~ dmulti(q[1:K], N)
      }", data_list, prior_list = prior_0, chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1)
      list(alternative = marglik(fit_1, data_list, multinomial_log_posterior, prior_1),
           null        = marglik(fit_0, data_list, multinomial_log_posterior, prior_0))
    })
  })
  scenario_text("bayes-factor-contingency", {
    log_dirichlet <- function(a) sum(lgamma(a)) - lgamma(sum(a))
    do.call(rbind, lapply(names(tables), function(name) {
      table <- tables[[name]]
      log_multinomial <- lgamma(sum(table) + 1) - sum(lgamma(table + 1))
      exact_1 <- log_multinomial + log_dirichlet(c(table) + 1) - log_dirichlet(rep(1, length(table)))
      exact_0 <- log_multinomial + log_dirichlet(rowSums(table) + 1) - log_dirichlet(rep(1, nrow(table))) +
        log_dirichlet(colSums(table) + 1) - log_dirichlet(rep(1, ncol(table)))
      rbind(
        scenario_bf_table("dependence", "independence", c(
          "BayesFactor (exact)" = bayesfactor_log_bf(BayesFactor::contingencyTableBF(table, sampleType = "jointMulti")),
          "closed form"         = exact_1 - exact_0,
          "bridge sampling"     = scenario_log_bf_bridge(contingency[[name]][["alternative"]], contingency[[name]][["null"]])),
          note = name),
        data.frame(H1 = c("log marginal likelihood", ""), H0 = c("dependence", "independence"),
                   method = "closed form - bridge sampling",
                   log_BF10 = round(c(exact_1 - contingency[[name]][["alternative"]][["logml"]],
                                      exact_0 - contingency[[name]][["null"]][["logml"]]), 3),
                   note = name)
      )
    }))
  })

  ### summary figure ----
  # log BF differences from the deterministic reference: BayesTools bridge
  # sampling (filled) and BayesFactor's Monte Carlo estimates (open)
  scenario_plot("bayes-factor-differences", {
    reference <- c(
      "one-sample t"    = exact_one_sample[["two_sided"]],
      "two-sample t"    = jzs_log_bf(sl$y, list(contr.orthonormal(levels(sl$group))[as.integer(sl$group), , drop = FALSE]), 1 / 2),
      "ANOVA supp"      = anova_exact[["supp"]],
      "ANOVA dose"      = anova_exact[["dose"]],
      "ANOVA add"       = anova_exact[["add"]],
      "ANOVA full"      = anova_exact[["full"]],
      "regression"      = jzs_log_bf(at$y, list(cbind(at$complaints_s)), sqrt(2) / 4)
    )
    bridge <- c(
      log_bf(one_sample[["alternative"]], one_sample[["null"]]),
      log_bf(two_sample[["alternative"]], two_sample[["null"]]),
      vapply(c("supp", "dose", "add", "full"), function(model) log_bf(anova[[model]], anova[["null"]]), numeric(1)),
      log_bf(regression[["alternative"]], regression[["null"]])
    )
    bayesfactor <- c(NA, NA, anova_bayesfactor[anova_models, "log_bf"], NA)
    position <- rev(seq_along(reference))
    limit <- max(0.05, abs(c(bridge, bayesfactor) - reference), na.rm = TRUE) * 1.1
    graphics::par(mar = c(4, 8, 1, 1))
    graphics::plot(bridge - reference, position, xlim = c(-limit, limit), pch = 19, yaxt = "n",
                   xlab = "log BF - reference", ylab = "")
    graphics::points(bayesfactor - reference, position, pch = 1)
    graphics::axis(2, at = position, labels = names(reference), las = 1)
    graphics::abline(v = 0, lty = 2, col = "grey50")
    graphics::legend("topright", c("BayesTools bridge sampling", "BayesFactor (Monte Carlo)"), pch = c(19, 1), bty = "n")
  })
})
