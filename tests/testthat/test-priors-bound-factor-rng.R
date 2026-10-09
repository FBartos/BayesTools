skip_if_not_test_profile("unit")

test_that("bound treatment mixture draws the declared independent product law", {
  mixture <- prior_factor_levels(prior_mixture(list(
    prior_factor("normal", list(1, 2), contrast = "treatment", prior_weights = 2),
    prior("point", list(0)))), c("c", "a", "b"))
  n <- 400L
  set.seed(71)
  actual <- rng(mixture, n, transform_factor_samples = FALSE)
  set.seed(71)
  gates <- sample(1:2, n, replace = TRUE, prob = c(2, 1))
  reference <- matrix(0, n, 2L)
  # The declared stream is column-major within each selected component.
  for(component in unique(gates)){
    rows <- gates == component
    reference[rows, ] <- if(component == 1L){
      matrix(stats::rnorm(sum(rows) * 2L, 1, 2), sum(rows), 2L)
    }else matrix(0, sum(rows), 2L)
  }
  expect_identical(unname(actual), structure(reference, components = gates))
  expect_false(isTRUE(all.equal(reference, reference[, rep(1L, 2L)])))
  set.seed(71)
  transformed <- rng(mixture, n, transform_factor_samples = TRUE)
  expect_equal(unname(transformed), structure(reference %*% t(attr(mixture, "factor_design")), components = gates))
  set.seed(71)
  expect_identical(rng(mixture, n, sample_components = TRUE), gates)
})

test_that("bound joint factor mixtures preserve the common t scale and exact design image", {
  for(contrast in c("meandif", "orthonormal")){
    mixture <- prior_factor_levels(prior_mixture(list(
      prior_factor("mt", list(0, 1.5, 5), contrast = contrast), prior("point", list(0)))), c("c", "a", "b"))
    set.seed(82)
    actual <- rng(mixture, 400, transform_factor_samples = FALSE)
    set.seed(82)
    gates <- sample(1:2, 400, replace = TRUE, prob = c(1, 1))
    reference <- matrix(0, 400, 2L)
    for(component in unique(gates)){
      rows <- gates == component
      reference[rows, ] <- if(component == 1L){
        mvtnorm::rmvt(sum(rows), delta = rep(0, 2), sigma = diag(1.5^2, 2), df = 5, type = "shifted")
      }else matrix(0, sum(rows), 2L)
    }
    expect_equal(as.vector(actual), as.vector(reference), tolerance = 0)
    set.seed(82)
    transformed <- rng(mixture, 400)
    expect_equal(as.vector(transformed), as.vector(reference %*% t(attr(mixture, "factor_design"))), tolerance = 0)
  }
})

test_that("bound spike gates and interaction designs use exact raw and cell dimensions", {
  for(contrast in c("treatment", "independent")){
    slab <- prior_factor("normal", list(0, 1), contrast = contrast)
    spike <- prior_factor_levels(prior_spike_and_slab(slab, prior("point", list(.4))), c("c", "a", "b"))
    set.seed(93)
    raw <- rng(spike, 200, transform_factor_samples = FALSE)
    set.seed(93)
    cells <- rng(spike, 200)
    expect_identical(attr(raw, "inclusion"), attr(cells, "inclusion"))
    expect_equal(as.vector(cells), as.vector(raw %*% t(attr(spike, "factor_design"))), tolerance = 0)
    expect_true(all(raw[attr(raw, "inclusion") == 0, ] == 0))
  }
  data <- expand.grid(g = factor(c("c", "a", "b")), h = factor(c("u", "v")), x = c(1, 2))
  mixture <- prior_mixture(list(prior_factor("normal", list(0, 1), contrast = "treatment"), prior("point", list(0))))
  for(formula in list(~ g * h, ~ g + g:x)){
    priors <- list(intercept = prior("point", list(0)), g = mixture,
      h = mixture, `g:h` = mixture, `g:x` = mixture)
    compiled <- JAGS_formula(formula, "mu", data, priors[c("intercept", attr(stats::terms(formula), "term.labels"))])
    bound <- tail(compiled$prior_list, 1L)[[1L]]
    set.seed(94)
    raw <- rng(bound, 50, transform_factor_samples = FALSE)
    set.seed(94)
    cells <- rng(bound, 50)
    expect_identical(dim(raw), c(50L, ncol(attr(bound, "factor_design"))))
    expect_identical(dim(cells), c(50L, nrow(attr(bound, "factor_design"))))
    expect_equal(as.vector(cells), as.vector(raw %*% t(attr(bound, "factor_design"))), tolerance = 0)
  }
})

test_that("standalone and initialization contracts remain separate from bound containers", {
  scalar <- prior_factor("normal", list(1, 2), contrast = "treatment")
  set.seed(14)
  standalone <- rng(scalar, 20, transform_factor_samples = FALSE)
  set.seed(14)
  expect_identical(standalone, stats::rnorm(20, 1, 2))
  joint <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  set.seed(15)
  expect_warning(rng(joint, 5), "assuming two factor levels")
  mixture <- prior_factor_levels(prior_mixture(list(joint, prior("point", list(0)))), 3)
  broken <- mixture
  broken[[1L]]$parameters$K <- 3
  expect_error(rng(broken, 5), "disagree with an explicit 'K'")
  prior_list <- list(b = prior_factor_levels(joint, 3))
  expect_identical(JAGS_get_inits(prior_list, 2, seed = 16), JAGS_get_inits(prior_list, 2, seed = 16))
  legacy <- prior_mixture(list(joint, prior("point", list(0))))
  for(i in seq_along(legacy)) legacy[[i]]$parameters$K <- 2
  expect_identical(dim(rng(legacy, 5)), c(5L, 3L))
})
