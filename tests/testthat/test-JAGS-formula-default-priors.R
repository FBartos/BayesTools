skip_if_not_test_profile("unit")

formula_default_prior_data <- function(){
  data.frame(
    x_cont1 = seq_len(12),
    x_cont2 = seq(12, 1),
    x_fac3  = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    x_fac2  = factor(rep(c("x", "y"), 6), levels = c("x", "y"))
  )
}

formula_default_prior_counters <- function(){
  counters <- new.env(parent = emptyenv())
  counters$continuous <- 0L
  counters$factor     <- 0L
  counters
}

formula_lazy_default_priors <- function(counters,
                                        continuous_prior = prior("normal", list(0, 1)),
                                        factor_prior = prior_factor("normal", list(0, 1), contrast = "treatment")){
  list(
    "__default_continuous" = function(){
      counters$continuous <- counters$continuous + 1L
      continuous_prior
    },
    "__default_factor" = function(){
      counters$factor <- counters$factor + 1L
      factor_prior
    }
  )
}

test_that("fully specified formula priors do not evaluate lazy defaults", {

  counters <- formula_default_prior_counters()
  prior_list <- c(
    list(
      intercept = prior("normal", list(0, 1)),
      x_cont1   = prior("normal", list(0, 2)),
      x_fac3    = prior_factor("normal", list(0, 3), contrast = "treatment")
    ),
    formula_lazy_default_priors(counters)
  )

  result <- JAGS_formula(
    formula    = ~ x_cont1 + x_fac3,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 0L)
  expect_equal(counters$factor, 0L)
  expect_equal(result$prior_list$mu_x_cont1$parameters$sd, 2)
  expect_equal(result$prior_list$mu_x_fac3$parameters$sd, 3)
})

test_that("missing continuous terms evaluate only the continuous lazy default once", {

  counters <- formula_default_prior_counters()
  prior_list <- formula_lazy_default_priors(counters)

  result <- JAGS_formula(
    formula    = ~ x_cont1 + x_cont2,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 1L)
  expect_equal(counters$factor, 0L)
  expect_equal(result$prior_list$mu_intercept$parameters$sd, 1)
  expect_equal(result$prior_list$mu_x_cont1$parameters$sd, 1)
  expect_equal(result$prior_list$mu_x_cont2$parameters$sd, 1)
})

test_that("missing factor terms evaluate only the factor lazy default once", {

  counters <- formula_default_prior_counters()
  prior_list <- c(
    list(
      intercept = prior("normal", list(0, 1)),
      x_cont1   = prior("normal", list(0, 2))
    ),
    formula_lazy_default_priors(counters)
  )

  result <- JAGS_formula(
    formula    = ~ x_cont1 + x_fac3 + x_fac2,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 0L)
  expect_equal(counters$factor, 1L)
  expect_equal(result$prior_list$mu_x_fac3$parameters$sd, 1)
  expect_equal(result$prior_list$mu_x_fac2$parameters$sd, 1)
})

test_that("missing continuous and factor terms evaluate both lazy defaults once", {

  counters <- formula_default_prior_counters()
  prior_list <- c(
    list(intercept = prior("normal", list(0, 1))),
    formula_lazy_default_priors(
      counters,
      continuous_prior = prior("normal", list(0, 2)),
      factor_prior     = prior_factor("normal", list(0, 3), contrast = "treatment")
    )
  )

  result <- JAGS_formula(
    formula    = ~ x_cont1 + x_fac3,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 1L)
  expect_equal(counters$factor, 1L)
  expect_equal(result$prior_list$mu_x_cont1$parameters$sd, 2)
  expect_equal(result$prior_list$mu_x_fac3$parameters$sd, 3)
})

test_that("factor interactions use the factor lazy default once", {

  counters <- formula_default_prior_counters()
  prior_list <- c(
    list(
      intercept = prior("normal", list(0, 1)),
      x_cont1   = prior("normal", list(0, 2))
    ),
    formula_lazy_default_priors(
      counters,
      factor_prior = prior_factor("mnormal", list(0, 3), contrast = "orthonormal")
    )
  )

  result <- JAGS_formula(
    formula    = ~ x_cont1 * x_fac3,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 0L)
  expect_equal(counters$factor, 1L)
  expect_equal(result$prior_list$mu_x_fac3$parameters$sd, 3)
  expect_equal(result$prior_list[["mu_x_cont1__xXx__x_fac3"]]$parameters$sd, 3)
})

test_that("lazy defaults returning non-priors error clearly", {

  data <- formula_default_prior_data()

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_continuous" = function() 1
      )
    ),
    "prior_list\\[\\[\"__default_continuous\"\\]\\].*BayesTools prior object"
  )

  expect_error(
    JAGS_formula(
      formula    = ~ x_fac3,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_factor" = function() list()
      )
    ),
    "prior_list\\[\\[\"__default_factor\"\\]\\].*BayesTools prior object"
  )
})

test_that("only zero-argument functions are accepted as lazy defaults", {

  data <- formula_default_prior_data()

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_continuous" = function(x) prior("normal", list(0, 1))
      )
    ),
    "zero-argument function",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_continuous" = function(x = 1) prior("normal", list(0, x))
      )
    ),
    "zero-argument function",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_continuous" = function(...) prior("normal", list(0, 1))
      )
    ),
    "zero-argument function",
    fixed = TRUE
  )
})

test_that("non-default prior list entries cannot be lazy functions", {

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = formula_default_prior_data(),
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x_cont1   = function() prior("normal", list(0, 1))
      )
    ),
    "'prior_list' must be a list of priors.",
    fixed = TRUE
  )
})

test_that("no-intercept formulas keep spike intercepts without consuming continuous defaults", {

  counters <- formula_default_prior_counters()
  prior_list <- c(
    list(x_cont1 = prior("normal", list(0, 2))),
    formula_lazy_default_priors(counters)
  )

  result <- JAGS_formula(
    formula    = ~ x_cont1 - 1,
    parameter  = "mu",
    data       = formula_default_prior_data(),
    prior_list = prior_list
  )

  expect_equal(counters$continuous, 0L)
  expect_equal(counters$factor, 0L)
  expect_true(is.prior.point(result$prior_list$mu_intercept))
  expect_equal(result$prior_list$mu_x_cont1$parameters$sd, 2)
})

test_that("wrong prior classes returned by lazy defaults use existing validation", {

  data <- formula_default_prior_data()

  expect_error(
    JAGS_formula(
      formula    = ~ x_cont1,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_continuous" = function(){
          prior_factor("normal", list(0, 1), contrast = "treatment")
        }
      )
    ),
    "Unsupported prior distribution defined for 'x_cont1' continuous variable.",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula    = ~ x_fac3,
      parameter  = "mu",
      data       = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "__default_factor" = function(){
          prior("normal", list(0, 1))
        }
      )
    ),
    "Unsupported prior distribution defined for 'x_fac3' factor variable.",
    fixed = TRUE
  )
})

test_that("existing prior-object defaults retain eager default behavior", {

  result <- JAGS_formula(
    formula   = ~ x_cont1 + x_fac3,
    parameter = "mu",
    data      = formula_default_prior_data(),
    prior_list = list(
      "__default_continuous" = prior("normal", list(0, 2)),
      "__default_factor"     = prior_factor("normal", list(0, 3), contrast = "treatment")
    )
  )

  expect_equal(result$prior_list$mu_intercept$parameters$sd, 2)
  expect_equal(result$prior_list$mu_x_cont1$parameters$sd, 2)
  expect_equal(result$prior_list$mu_x_fac3$parameters$sd, 3)
})
