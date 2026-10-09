skip_if_not_test_profile("unit")

test_that("R133 ordinary scalar and vector bridge inputs distinguish missing names", {
  scalar <- list(mu = prior("normal", list(0, 1)))
  vector <- list(o = prior("mnormal", list(mean = 0, sd = 1, K = 2)))
  for(samples in list(c(other = 1), list(other = 1))){
    condition <- tryCatch(JAGS_marglik_parameters(samples, scalar), error = identity)
    expect_s3_class(condition, "BayesTools_missing_monitored_columns")
    expect_s3_class(condition, "BayesTools_marglik_input")
    if(inherits(condition, "error")){
      expect_identical(conditionMessage(condition), "'samples' does not contain all monitored scalar prior parameters.")
      expect_null(conditionCall(condition))
    }
  }
  for(samples in list(c("o[1]" = 1), list("o[1]" = 1))){
    condition <- tryCatch(JAGS_marglik_parameters(samples, vector), error = identity)
    expect_s3_class(condition, "BayesTools_missing_monitored_columns")
    expect_s3_class(condition, "BayesTools_marglik_input")
    if(inherits(condition, "error")){
      expect_identical(conditionMessage(condition), "'samples' does not contain all monitored vector prior parameters.")
      expect_null(conditionCall(condition))
    }
  }
  expect_identical(JAGS_marglik_parameters(c(mu = 2), scalar), list(mu = 2))
  expect_identical(JAGS_marglik_parameters(list(mu = 2), scalar), list(mu = 2))
  expect_identical(JAGS_marglik_parameters(c("o[2]" = 2, other = 9, "o[1]" = 1), vector),
    list(o = c("o[1]" = 1, "o[2]" = 2)))
  expect_identical(JAGS_marglik_parameters(list("o[2]" = 2, "o[1]" = 1), vector),
    list(o = list("o[1]" = 1, "o[2]" = 2)))
  expect_identical(JAGS_marglik_parameters(c(mu = NA_real_), scalar), list(mu = NA_real_))
  expect_identical(JAGS_marglik_parameters(c("o[1]" = NA_real_, "o[2]" = 2), vector),
    list(o = c("o[1]" = NA_real_, "o[2]" = 2)))
  expect_identical(JAGS_marglik_parameters(numeric(), list(mu = prior("point", list(2)))), list(mu = 2))
  expect_identical(JAGS_marglik_parameters(numeric(), list(o = prior("mpoint", list(location = 2, K = 2)))), list(o = c(2, 2)))
})
