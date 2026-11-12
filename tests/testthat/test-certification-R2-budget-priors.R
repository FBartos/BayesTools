skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Bayes factors of R-squared budget priors (.R2_prior()) against the
# intercept-only model in a balanced design with mutually orthogonal centered
# term blocks: continuous x and z and a three-level factor f (orthonormal
# contrasts). With a flat intercept and p(sigma) proportional to 1/sigma, the
# Bayes factors depend on the data only through n and the term shares
# R2_t = ||P_t y_c||^2 / SST (per coded column for the kernels), and reduce
# to expectations over the budget of closed-form kernels, with
# s ~ Gamma((n - 1) / 2):
#   exact law: E[exp(-n Q / 2) E_s prod_t 0F1(; d_t / 2; n s q_t R2_t / 2)],
#     where E_s 0F1(; c; s y) = 1F1(k; c; y), and two one-column terms combine
#     by cosh(a) cosh(b) = (cosh(a + b) + cosh(a - b)) / 2;
#   normal kernel: E[prod_j (1 + n v_j)^(-1/2) (1 - sum_j c_j R2_j)^(-k)],
#     c_j = n v_j / (1 + n v_j), v_j = W phi_j;
#   Laplace kernel: the normal kernel with v = W phi E, E ~ Exp(1).
# The expectations over r ~ Beta(a, b) (W = r / (1 - r)), the Dirichlet share
# and E are adaptive quadrature (relative tolerance 1e-10), with 1F1 summed
# as a positive series in log space; these routes are independent of the
# package code. The fits use an N(0, 100^2) intercept and sigma ~
# lognormal(0, 10) truncated to [0.01, 100], whose density in log(sigma) is
# flat to 1e-3 over the posterior range, so the flat-intercept, 1/sigma Bayes
# factors are the targets up to Monte Carlo error. The bridge log Bayes
# factors (3 x 5000 draws per fit) are held to 0.05; in development the
# errors against the project's deterministic engine were at most 0.02 across
# laws, conventions and a separate compulsory-covariate budget.

R2_cert_data <- local({
  set.seed(20261111)
  n_per <- 16L
  f <- factor(rep(c("a", "b", "c"), each = n_per))
  n <- length(f)
  Xf <- contr.orthonormal(levels(f))[as.integer(f), , drop = FALSE]
  orthogonal <- function(v, M){
    v <- stats::residuals(stats::lm(v ~ M))
    (v - mean(v)) / sqrt(mean((v - mean(v))^2))
  }
  x <- orthogonal(stats::rnorm(n), Xf)
  z <- orthogonal(stats::rnorm(n), cbind(Xf, x))
  y <- 0.4 * Xf[, 2] + 0.35 * x + 0.1 * z + stats::rnorm(n)
  data.frame(y = y, f = f, x = x, z = z)
})

R2_cert_shares <- local({
  yc <- R2_cert_data$y - mean(R2_cert_data$y)
  SST <- sum(yc^2)
  share <- function(v){
    v <- v - mean(v)
    sum(v * yc)^2 / sum(v^2) / SST
  }
  Xf <- contr.orthonormal(levels(R2_cert_data$f))[as.integer(R2_cert_data$f), , drop = FALSE]
  list(x = share(R2_cert_data$x), z = share(R2_cert_data$z), f = c(share(Xf[, 1]), share(Xf[, 2])))
})

R2_cert_log1F1 <- function(k, c, y){

  vapply(y, function(value){
    if(value == 0) return(0)
    m <- 0:(ceiling(4 * value + 4 * sqrt(value) * sqrt(k) + 60))
    terms <- lgamma(k + m) - lgamma(k) + lgamma(c) - lgamma(c + m) - lgamma(m + 1) + m * log(value)
    top <- max(terms)
    top + log(sum(exp(terms - top)))
  }, numeric(1))
}

R2_cert_integrate <- function(f, lower = 0, upper = 1){
  stats::integrate(f, lower, upper, rel.tol = 1e-10, subdivisions = 2000L)$value
}

R2_cert_reference <- function(case, n = nrow(R2_cert_data)){

  k <- (n - 1) / 2
  a <- case$a
  b <- case$b
  W_of <- function(r) r / (1 - r)
  log_bf <- switch(
    case$reference,
    exact_one = {
      R2 <- sum(R2_cert_shares[[case$terms]])
      d <- length(R2_cert_shares[[case$terms]])
      log(R2_cert_integrate(function(r){
        W <- W_of(r)
        stats::dbeta(r, a, b) * exp(-n * W / 2 + R2_cert_log1F1(k, d / 2, n * W * R2 / 2))
      }))
    },
    exact_two = {
      alpha <- case$alpha
      R2 <- c(R2_cert_shares$x, R2_cert_shares$z)
      log(R2_cert_integrate(function(r){
        vapply(r, function(r_i){
          W <- W_of(r_i)
          R2_cert_integrate(function(phi){
            q1 <- W * phi
            q2 <- W * (1 - phi)
            a1 <- sqrt(2 * n * q1 * R2[1])
            a2 <- sqrt(2 * n * q2 * R2[2])
            plus <- R2_cert_log1F1(k, 1 / 2, (a1 + a2)^2 / 4)
            minus <- R2_cert_log1F1(k, 1 / 2, (a1 - a2)^2 / 4)
            top <- pmax(plus, minus)
            stats::dbeta(phi, alpha[1], alpha[2]) *
              exp(-n * W / 2 + top + log((exp(plus - top) + exp(minus - top)) / 2))
          })
        }, numeric(1)) * stats::dbeta(r, a, b)
      }))
    },
    normal_two = {
      R2 <- R2_cert_shares$f
      log(R2_cert_integrate(function(r){
        vapply(r, function(r_i){
          W <- W_of(r_i)
          R2_cert_integrate(function(phi){
            v1 <- W * phi
            v2 <- W * (1 - phi)
            c1 <- n * v1 / (1 + n * v1)
            c2 <- n * v2 / (1 + n * v2)
            stats::dbeta(phi, case$alpha, case$alpha) *
              exp(-0.5 * log1p(n * v1) - 0.5 * log1p(n * v2) - k * log1p(-(c1 * R2[1] + c2 * R2[2])))
          })
        }, numeric(1)) * stats::dbeta(r, a, b)
      }))
    },
    laplace_one = {
      R2 <- R2_cert_shares[[case$terms]]
      log(R2_cert_integrate(function(r){
        vapply(r, function(r_i){
          W <- W_of(r_i)
          R2_cert_integrate(function(e){
            v <- W * e
            cv <- n * v / (1 + n * v)
            stats::dexp(e) * exp(-0.5 * log1p(n * v) - k * log1p(-cv * R2))
          }, 0, Inf)
        }, numeric(1)) * stats::dbeta(r, a, b)
      }))
    }
  )
  log_bf
}

R2_cert_fit <- function(terms, budget, seed){

  formula <- if(length(terms)) stats::reformulate(terms) else ~ 1
  formula_priors <- list(intercept = prior("normal", list(0, 100)))
  if(length(terms)) formula_priors[["R2"]] <- budget
  arguments <- list(
    model_syntax       = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data               = list(y = R2_cert_data$y, N = nrow(R2_cert_data)),
    prior_list         = list(sigma = prior("lognormal", list(0, 10), truncation = list(lower = 0.01, upper = 100))),
    formula_list       = list(mu = formula),
    formula_data_list  = list(mu = R2_cert_data),
    formula_prior_list = list(mu = formula_priors)
  )
  fit <- do.call(JAGS_fit, c(arguments, list(chains = 3, adapt = 1000, burnin = 2000, sample = 5000,
    seed = seed, silent = TRUE)))
  marglik <- JAGS_bridgesampling(fit,
    log_posterior = function(parameters, data){
      sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
    },
    data = arguments$data, prior_list = arguments$prior_list,
    formula_list = arguments$formula_list, formula_data_list = arguments$formula_data_list,
    formula_prior_list = arguments$formula_prior_list, repetitions = 3, seed = seed)
  list(fit = fit, marglik = marglik, prior_weights = 1)
}

R2_cert_budget <- function(alpha, law, a = 1, b = 2, convention = "common"){
  .R2_prior(terms = names(alpha), total = prior("beta", list(if(identical(convention, "linked")) sum(alpha) else a, b)),
    allocation = if(length(alpha) > 1L) prior("dirichlet", list(alpha = unname(alpha))),
    law = law, convention = convention, concentration = "inherited", sigma = "sigma", contrast = "orthonormal")
}

R2_cert_cases <- list(
  list(name = "exact x", terms = "x", reference = "exact_one", a = 1, b = 2,
       budget = R2_cert_budget(c(x = 1, z = 1), "exact")),
  list(name = "exact f", terms = "f", reference = "exact_one", a = 0.5, b = 1,
       budget = R2_cert_budget(c(f = 1), "exact", a = 0.5, b = 1)),
  list(name = "exact x + z", terms = c("x", "z"), reference = "exact_two", a = 1, b = 2, alpha = c(1, 0.5),
       budget = R2_cert_budget(c(x = 1, z = 0.5), "exact")),
  list(name = "linked x + z", terms = c("x", "z"), reference = "exact_two", a = 1.5, b = 2, alpha = c(1, 0.5),
       budget = R2_cert_budget(c(x = 1, z = 0.5), "exact", convention = "linked")),
  list(name = "normal f", terms = "f", reference = "normal_two", a = 1, b = 2, alpha = 1,
       budget = R2_cert_budget(c(f = 1), "normal")),
  list(name = "Laplace x", terms = "x", reference = "laplace_one", a = 1, b = 2,
       budget = R2_cert_budget(c(x = 1), "laplace"))
)

R2_cert_null <- R2_cert_fit(character(), NULL, seed = 1)
R2_cert_fits <- lapply(seq_along(R2_cert_cases), function(i){
  R2_cert_fit(R2_cert_cases[[i]]$terms, R2_cert_cases[[i]]$budget, seed = i + 1)
})
names(R2_cert_fits) <- vapply(R2_cert_cases, `[[`, character(1), "name")

test_that("bridge-sampling Bayes factors of R2 budget priors match independent quadrature", {

  problems <- expectation_problems({
    for(i in seq_along(R2_cert_cases)){
      case <- R2_cert_cases[[i]]
      exact <- R2_cert_reference(case)
      bridge <- R2_cert_fits[[i]]$marglik$logml - R2_cert_null$marglik$logml
      cat(sprintf("%-14s bridge log BF %8.4f  quadrature %8.4f  difference %7.4f
", case$name, bridge, exact, bridge - exact))
      expect_true(is.finite(bridge) && abs(bridge - exact) <= 0.05,
        info = sprintf("%s: bridge log BF %.4f, exact %.4f", case$name, bridge, exact))
    }
  })
  expect_identical(problems, character())
})

test_that("monitored R2 coefficients equal the registered node evaluator on the posterior draws", {

  for(name in c("exact x + z", "exact f", "normal f", "Laplace x")){
    fit <- R2_cert_fits[[name]]$fit
    draws <- .extract_posterior_samples(fit, as_list = FALSE)
    priors <- attr(fit, "prior_list")
    for(term in names(priors)[vapply(priors, is.prior.R2, logical(1))]){
      node <- .bt_deterministic_nodes_fit(fit)[[term]]
      expect_identical(node$family, "R2_coefficient")
      reduced <- draws[, setdiff(colnames(draws), node$coordinates), drop = FALSE]
      values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(reduced, priors))
      monitored <- matrix(as.numeric(draws[, node$coordinates, drop = FALSE]), nrow = nrow(draws))
      expect_equal(matrix(as.numeric(values), nrow = nrow(values)), monitored, tolerance = 1e-12)
    }
  }
})

test_that("predictions of R2 fits use the residual-scaled coefficients", {

  fit <- R2_cert_fits[["exact x + z"]]$fit
  draws <- .extract_posterior_samples(fit, as_list = FALSE)
  newdata <- R2_cert_data[1:4, ]
  prediction <- JAGS_predict_formula(fit, parameter = "mu", data = newdata)
  expected <- t(draws[, "mu_intercept"] + draws[, "sigma"] *
    (draws[, "mu_x"] %o% newdata$x + draws[, "mu_z"] %o% newdata$z))
  expect_equal(unname(prediction$mean), unname(expected), tolerance = 1e-10)
})

test_that("model averaging over budget submodels reproduces the inclusion Bayes factor", {

  budget <- R2_cert_budget(c(x = 1, z = 0.5), "exact")
  models <- list(
    R2_cert_null,
    R2_cert_fit("x", budget, seed = 21),
    R2_cert_fit("z", budget, seed = 22),
    R2_cert_fits[["exact x + z"]]
  )
  models <- models_inference(models)
  inference <- ensemble_inference(models, parameters = "mu_z",
    is_null_list = list(mu_z = c(TRUE, TRUE, FALSE, FALSE)), conditional = FALSE)
  reference <- c(
    0,
    R2_cert_reference(list(terms = "x", reference = "exact_one", a = 1, b = 2)),
    R2_cert_reference(list(terms = "z", reference = "exact_one", a = 1, b = 2)),
    R2_cert_reference(R2_cert_cases[[3]])
  )
  log_inclusion <- log(sum(exp(reference[3:4]))) - log(sum(exp(reference[1:2])))
  expect_lt(abs(log(inference$mu_z$BF) - log_inclusion), 0.1)
  mixed <- mix_posteriors(models, parameters = c("mu_z", "mu_R2"),
    is_null_list = list(mu_z = c(TRUE, TRUE, FALSE, FALSE), mu_R2 = c(TRUE, FALSE, FALSE, FALSE)), seed = 1)
  expect_s3_class(mixed$mu_z, "mixed_posteriors")
  expect_true(all(mixed$mu_R2 >= 0 & mixed$mu_R2 < 1))
})

test_that("one gated fit reproduces the model average of its submodels", {

  # inclusion gates with a shared Beta(1, 1) rate and a spike-and-slab total
  # (null weight 1/2): the prior over the configurations of x and z is
  # P(none) = 1/2, P(x) = P(z) = 1/8, P(x, z) = 1/4, and given each
  # configuration the budget is the prior of that model as a separate fit;
  # the reference posterior combines it with the independent quadrature
  # Bayes factors of the separate models
  gated <- .R2_prior(terms = c("x", "z"),
    total = prior_spike_and_slab(prior("beta", list(1, 2)), prior_inclusion = prior("spike", list(0.5))),
    allocation = prior("dirichlet", list(alpha = c(1, 0.5))), law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma", inclusion = prior("beta", list(1, 1)))
  fit <- JAGS_fit(
    model_syntax       = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data               = list(y = R2_cert_data$y, N = nrow(R2_cert_data)),
    prior_list         = list(sigma = prior("lognormal", list(0, 10), truncation = list(lower = 0.01, upper = 100))),
    formula_list       = list(mu = ~ x + z),
    formula_data_list  = list(mu = R2_cert_data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 100)), R2 = gated)),
    chains = 4, adapt = 1000, burnin = 2000, sample = 20000, seed = 31, silent = TRUE)
  draws <- .extract_posterior_samples(fit, as_list = FALSE)
  open <- draws[, "mu_R2_indicator"] * draws[, c("mu_R2_gate[1]", "mu_R2_gate[2]")]
  observed <- c(none = mean(rowSums(open) == 0), x = mean(open[, 1] == 1 & open[, 2] == 0),
                z = mean(open[, 1] == 0 & open[, 2] == 1), xz = mean(rowSums(open) == 2))

  log_BF <- c(
    none = 0,
    x  = R2_cert_reference(list(terms = "x", reference = "exact_one", a = 1, b = 2)),
    z  = R2_cert_reference(list(terms = "z", reference = "exact_one", a = 1, b = 2)),
    xz = R2_cert_reference(R2_cert_cases[[3]])
  )
  log_post <- log(c(1 / 2, 1 / 8, 1 / 8, 1 / 4)) + log_BF
  expected <- exp(log_post - max(log_post))
  expected <- expected / sum(expected)
  cat("\n", sprintf("%-5s single-fit %.4f  model average %.4f\n", names(expected), observed, expected), sep = "")
  expect_true(all(abs(observed - expected) < 0.02))
  inference <- runjags_inference_table(fit)
  expect_equal(inference["(mu) R2{total}", "post_prob"], 1 - observed[["none"]])
})
