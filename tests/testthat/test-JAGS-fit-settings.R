skip_if_not_test_profile("unit")

test_that("JAGS fit settings reject missing and non-finite controls", {
  valid <- list(
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    thin = 1,
    autofit = FALSE,
    parallel = FALSE,
    cores = 1,
    silent = TRUE,
    seed = 1
  )
  invalid <- list(
    chains = NA_real_,
    adapt = Inf,
    burnin = NA_real_,
    sample = Inf,
    thin = NA_real_,
    autofit = NA,
    parallel = NA,
    cores = Inf,
    silent = NA,
    seed = NaN
  )

  expect_silent(do.call(JAGS_check_and_list_fit_settings, valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_check_and_list_fit_settings, settings),
      paste0("'", name, "'")
    )
  }
})

test_that("JAGS autofit settings reject missing and non-finite controls", {
  valid <- list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = FALSE
  )
  invalid <- list(
    max_Rhat = Inf,
    min_ESS = Inf,
    max_error = NaN,
    max_SD_error = NA_real_,
    max_time = list(time = Inf, unit = "mins"),
    sample_extend = Inf,
    restarts = NA_real_,
    max_extend = Inf,
    check_indicators = NA
  )

  expect_silent(JAGS_check_and_list_autofit_settings(valid))
  for(name in names(invalid)){
    settings <- valid
    settings[[name]] <- invalid[[name]]
    expect_error(
      JAGS_check_and_list_autofit_settings(settings),
      paste0("'", if(name == "max_time") "max_time:time" else name, "'")
    )
  }

  invalid_unit <- valid
  invalid_unit$max_time$unit <- NA_character_
  expect_error(
    JAGS_check_and_list_autofit_settings(invalid_unit),
    "'max_time:unit'"
  )
})

test_that("JAGS_extend validates runtime controls before extension", {
  fit <- structure(list(), class = "BayesTools_fit")
  invalid <- list(
    parallel = NA,
    cores = Inf,
    silent = NA,
    seed = NaN
  )

  for(name in names(invalid)){
    arguments <- list(fit = fit)
    arguments[[name]] <- invalid[[name]]
    expect_error(
      do.call(JAGS_extend, arguments),
      paste0("'", name, "'")
    )
  }
})
