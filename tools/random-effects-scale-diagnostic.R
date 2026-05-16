pkgload::load_all(".")

if(!requireNamespace("ggplot2", quietly = TRUE)){
  stop("The diagnostic figure requires ggplot2.", call. = FALSE)
}

library(BayesTools)
library(lme4)

set.seed(20260513)

output_file <- file.path("tools", "random-effects-scale-diagnostic.png")
n_prior <- 50000
parameter_levels <- c(
  "sd(intercept | Subject)",
  "sd(Days | Subject)",
  "cor(intercept, Days | Subject)"
)

data("sleepstudy", package = "lme4")
sleepstudy$Subject <- factor(sleepstudy$Subject)
sleepstudy$Reaction_z <- as.numeric(scale(sleepstudy$Reaction))

days_mean <- mean(sleepstudy$Days)
days_sd <- stats::sd(sleepstudy$Days)

unit_prior <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
unit_sd_prior <- prior(
  distribution = "normal",
  parameters = list(mean = 0, sd = 1),
  truncation = list(0, Inf)
)
raw_slope_prior <- prior(
  distribution = "normal",
  parameters = list(mean = 0, sd = 1 / days_sd)
)
raw_slope_sd_prior <- prior(
  distribution = "normal",
  parameters = list(mean = 0, sd = 1 / days_sd),
  truncation = list(0, Inf)
)

sleep_likelihood <-
"model{
  for(i in 1:N){
    Reaction_z[i] ~ dnorm(mu[i], pow(sigma, -2))
  }
}
"

sleep_formula <- list(mu = ~ Days + random(1 + Days | Subject, name = "Subject", covariance = "us"))
sleep_data <- list(Reaction_z = sleepstudy$Reaction_z, N = nrow(sleepstudy))
monitor_cor <- random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)

scaled_fixed_priors <- list(intercept = unit_prior, Days = unit_prior)
scaled_random_priors <- list(
  mu = prior_random(
    Subject = random_block(
      sd = unit_sd_prior,
      cor = prior_lkj(eta = 2),
      monitor = monitor_cor
    )
  )
)

raw_fixed_priors <- list(intercept = unit_prior, Days = raw_slope_prior)
raw_random_priors <- list(
  mu = prior_random(
    Subject = random_block(
      sd = unit_sd_prior,
      cor = prior_lkj(eta = 2),
      monitor = monitor_cor,
      terms = list(Days = raw_slope_sd_prior)
    )
  )
)

fit_scaled <- JAGS_fit(
  model_syntax = sleep_likelihood,
  data = sleep_data,
  prior_list = list(sigma = unit_sd_prior),
  formula_list = sleep_formula,
  formula_data_list = list(mu = sleepstudy),
  formula_prior_list = list(mu = scaled_fixed_priors),
  formula_scale_list = list(mu = list(Days = TRUE)),
  formula_random_prior_list = scaled_random_priors,
  chains = 2,
  adapt = 500,
  burnin = 500,
  sample = 1500,
  seed = 101
)

fit_raw <- JAGS_fit(
  model_syntax = sleep_likelihood,
  data = sleep_data,
  prior_list = list(sigma = unit_sd_prior),
  formula_list = sleep_formula,
  formula_data_list = list(mu = sleepstudy),
  formula_prior_list = list(mu = raw_fixed_priors),
  formula_random_prior_list = raw_random_priors,
  chains = 2,
  adapt = 500,
  burnin = 500,
  sample = 1500,
  seed = 102
)

posterior_random_draws <- function(fit, case, transform_scaled = FALSE){
  samples <- JAGS_estimates_table(
    fit,
    transform_scaled = transform_scaled,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  data.frame(
    parameter = rep(parameter_levels, each = nrow(samples)),
    value = c(
      samples[, "(mu) sd(intercept | Subject)"],
      samples[, "(mu) sd(Days | Subject)"],
      samples[, "(mu) cor(intercept,Days | Subject)"]
    ),
    case = case,
    distribution = "Posterior"
  )
}

lkj2_rho <- function(n, eta = 2){
  2 * stats::rbeta(n, eta, eta) - 1
}

random_prior_draws <- function(n, case, scaled_then_backtransform){
  if(scaled_then_backtransform){
    sd_intercept_scaled <- abs(stats::rnorm(n, 0, 1))
    sd_days_scaled <- abs(stats::rnorm(n, 0, 1))
    cor_scaled <- lkj2_rho(n)
    cov_scaled <- cor_scaled * sd_intercept_scaled * sd_days_scaled

    raw_intercept_var <- sd_intercept_scaled^2 +
      (days_mean / days_sd)^2 * sd_days_scaled^2 -
      2 * (days_mean / days_sd) * cov_scaled
    raw_days_var <- sd_days_scaled^2 / days_sd^2
    raw_cov <- cov_scaled / days_sd - days_mean * sd_days_scaled^2 / days_sd^2

    sd_intercept <- sqrt(raw_intercept_var)
    sd_days <- sqrt(raw_days_var)
    cor <- raw_cov / (sd_intercept * sd_days)
  }else{
    sd_intercept <- abs(stats::rnorm(n, 0, 1))
    sd_days <- abs(stats::rnorm(n, 0, 1 / days_sd))
    cor <- lkj2_rho(n)
  }

  data.frame(
    parameter = rep(parameter_levels, each = n),
    value = c(sd_intercept, sd_days, cor),
    case = case,
    distribution = "Prior"
  )
}

lme4_raw <- lme4::lmer(
  Reaction_z ~ Days + (1 + Days | Subject),
  data = sleepstudy,
  REML = FALSE
)
lme4_vc <- lme4::VarCorr(lme4_raw)[["Subject"]]
lme4_reference <- data.frame(
  parameter = parameter_levels,
  value = c(
    attr(lme4_vc, "stddev")[["(Intercept)"]],
    attr(lme4_vc, "stddev")[["Days"]],
    attr(lme4_vc, "correlation")[1, 2]
  )
)

plot_data <- rbind(
  random_prior_draws(n_prior, "Scaled fit, backtransformed", scaled_then_backtransform = TRUE),
  posterior_random_draws(fit_scaled, "Scaled fit, backtransformed", transform_scaled = TRUE),
  random_prior_draws(n_prior, "Raw fit, scale-matched independent", scaled_then_backtransform = FALSE),
  posterior_random_draws(fit_raw, "Raw fit, scale-matched independent", transform_scaled = FALSE)
)

plot_data$case <- factor(
  plot_data$case,
  levels = c("Scaled fit, backtransformed", "Raw fit, scale-matched independent")
)
plot_data$distribution <- factor(plot_data$distribution, levels = c("Prior", "Posterior"))
plot_data$parameter <- factor(
  plot_data$parameter,
  levels = parameter_levels
)
lme4_reference$parameter <- factor(lme4_reference$parameter, levels = parameter_levels)

figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = value, color = case, linetype = distribution)
) +
  ggplot2::geom_density(linewidth = 0.85, adjust = 1.15) +
  ggplot2::geom_vline(
    data = lme4_reference,
    ggplot2::aes(xintercept = value),
    inherit.aes = FALSE,
    color = "gray25",
    linetype = "dotted",
    linewidth = 0.65
  ) +
  ggplot2::facet_wrap(~ parameter, scales = "free", ncol = 3) +
  ggplot2::scale_color_manual(values = c("#1F77B4", "#D95F02")) +
  ggplot2::scale_linetype_manual(values = c("Prior" = "longdash", "Posterior" = "solid")) +
  ggplot2::labs(
    title = "Random-effect prior and posterior on the raw Days scale",
    subtitle = "Dotted vertical lines are lme4 maximum-likelihood estimates",
    x = NULL,
    y = "Density",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.4, "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title.position = "plot"
  )

ggplot2::ggsave(output_file, figure, width = 11, height = 4.5, dpi = 180)

message("Wrote ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
