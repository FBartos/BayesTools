skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

# Savage-Dickey ordinates away from the posterior mode. With a known residual
# SD of 1 and N(0, 1) priors on the intercept and the slope, the slope has the
# normal posterior N(m, s^2), so the point null x = c has the exact
# BF10 = phi(c; 0, 1) / phi(c; m, s^2). The null is placed z = 0, 1, 2, 3
# posterior SDs from m. The normal ordinate errs only through the Monte Carlo
# error of the posterior mean and SD, sqrt(z^2 + (z^2 - 1)^2 / 2) / sqrt(ESS)
# on the log scale; the kernel ordinate additionally through its relative SD
# sqrt(R(K) / (ESS h f)) at the posterior density f, which grows in the tail,
# and its smoothing bias h^2 (z^2 - 1) / (2 s^2) (Gaussian kernel, R(K) =
# 1 / (2 sqrt(pi)), normal-reference bandwidth h over 40000 draws). Both are
# held to four of their Monte Carlo SDs plus the bias; at z = 3 the kernel
# tolerance is about 0.5 on the log scale, which is why tail nulls call for the
# normal (or a precomputed) ordinate.

ordinate_data <- local({
  set.seed(5)
  data <- data.frame(x = stats::rnorm(40))
  data$y <- 0.1 + 0.3 * data$x + stats::rnorm(40)
  data
})
ordinate_posterior <- local({
  X <- cbind(1, ordinate_data$x)
  V <- solve(crossprod(X) + diag(2))
  c(mean = drop(V %*% crossprod(X, ordinate_data$y))[[2L]], sd = sqrt(V[2L, 2L]))
})

test_that("Savage-Dickey ordinates stay within their Monte Carlo error away from the mode", {

  unit <- prior("normal", list(mean = 0, sd = 1))
  fit  <- JAGS_fit(
    "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) } }",
    list(y = ordinate_data$y, N = nrow(ordinate_data)),
    formula_list = list(mu = ~ x), formula_data_list = list(mu = ordinate_data),
    formula_prior_list = list(mu = list(intercept = unit, x = unit)),
    chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
  )
  ess <- min(as.data.frame(JAGS_estimates_table(fit))[["ESS"]])
  m <- ordinate_posterior[["mean"]]
  s <- ordinate_posterior[["sd"]]
  nulls <- round(m - c(0, 1, 2, 3) * s, 4)
  z     <- (m - nulls) / s
  exact <- stats::dnorm(nulls, 0, 1, log = TRUE) - stats::dnorm(nulls, m, s, log = TRUE)
  statements <- paste0("x = ", format(nulls, nsmall = 4))

  normal_sd <- sqrt(z^2 + (z^2 - 1)^2 / 2) / sqrt(ess)
  h         <- 0.9 * s * 40000^(-1 / 5)
  kernel_sd <- sqrt(1 / (2 * sqrt(pi)) / (ess * h * stats::dnorm(nulls, m, s)))
  bias      <- abs(h^2 * (z^2 - 1) / (2 * s^2))

  normal <- log(as.data.frame(hypothesis_BF(fit, statements, seed = 1, density_method = "normal"))[["BF10"]])
  expect_true(all(abs(normal - exact) < 4 * normal_sd))

  kernel <- log(as.data.frame(hypothesis_BF(fit, statements, seed = 1))[["BF10"]])
  expect_true(all(abs(kernel - exact) < 4 * (normal_sd + kernel_sd) + bias))
})
