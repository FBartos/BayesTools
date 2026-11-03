skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Bridge sampling with Dirichlet priors against closed forms. For a joint
# multinomial contingency table y with cell probabilities p ~ Dirichlet(alpha),
# log p(y) = log(N! / prod y!) + log B(y + alpha) - log B(alpha) with the
# multivariate beta function B; under independence p_ij = r_i c_j with
# r ~ Dirichlet(alpha_r) and c ~ Dirichlet(alpha_c), the row and column sums
# give log p(y) = log(N! / prod y!) + log B(y_r + alpha_r) - log B(alpha_r) +
# log B(y_c + alpha_c) - log B(alpha_c). The bridge estimates of 4 x 5000
# draws of these two- to six-dimensional simplex targets are held to 0.02 on
# the log scale (about four times the largest error seen in development,
# 0.008 in a two-way table with 32 observations).

simplex_table <- matrix(c(12, 7, 9, 15, 4, 11), nrow = 2L)
simplex_log_beta <- function(alpha) sum(lgamma(alpha)) - lgamma(sum(alpha))
simplex_exact <- function(table, alpha, alpha_rows, alpha_columns) {
  multinomial <- lgamma(sum(table) + 1) - sum(lgamma(table + 1))
  c(dependence   = multinomial + simplex_log_beta(c(table) + alpha) - simplex_log_beta(alpha),
    independence = multinomial + simplex_log_beta(rowSums(table) + alpha_rows) - simplex_log_beta(alpha_rows) +
      simplex_log_beta(colSums(table) + alpha_columns) - simplex_log_beta(alpha_columns))
}
simplex_log_posterior <- function(parameters, data) {
  probabilities <- if (is.null(parameters[["p"]])) c(outer(parameters[["r"]], parameters[["cc"]])) else parameters[["p"]]
  stats::dmultinom(data[["y"]], prob = probabilities, log = TRUE)
}

test_that("bridge sampling recovers the closed-form marginal likelihoods of Dirichlet models", {

  data_list <- list(y = c(simplex_table), K = length(simplex_table), N = sum(simplex_table),
                    I = nrow(simplex_table), J = ncol(simplex_table))
  # uniform and non-uniform concentrations
  for (alpha in list(rep(1, 6), c(2, 1, 3, 1, 2, 0.5))) {
    alpha_rows    <- if (all(alpha == 1)) c(1, 1) else c(0.5, 2)
    alpha_columns <- if (all(alpha == 1)) c(1, 1, 1) else c(3, 1, 0.5)
    prior_1 <- list(p = prior("dirichlet", list(alpha = alpha)))
    prior_0 <- list(r = prior("dirichlet", list(alpha = alpha_rows)), cc = prior("dirichlet", list(alpha = alpha_columns)))
    fit_1 <- JAGS_fit("model{ y[1:K] ~ dmulti(p[1:K], N) }", data_list[c("y", "K", "N")], prior_list = prior_1,
                      chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
    fit_0 <- JAGS_fit("model{
      for(i in 1:I){ for(j in 1:J){ q[(j - 1) * I + i] <- r[i] * cc[j] } }
      y[1:K] ~ dmulti(q[1:K], N)
    }", data_list, prior_list = prior_0, chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
    logml <- c(
      dependence   = JAGS_bridgesampling(fit_1, log_posterior = simplex_log_posterior, data = data_list, prior_list = prior_1, seed = 1)$logml,
      independence = JAGS_bridgesampling(fit_0, log_posterior = simplex_log_posterior, data = data_list, prior_list = prior_0, seed = 1)$logml
    )
    expect_lt(max(abs(logml - simplex_exact(simplex_table, alpha, alpha_rows, alpha_columns))), 0.02)
  }
})
