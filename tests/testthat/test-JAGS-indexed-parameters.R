skip_if_not_test_profile("unit")

test_that("JAGS indexed parameter helpers match exact sorted indices", {

  columns <- c(
    "omega[1]", "omega[10]", "omega[2]", "omega", "omega_extra[1]",
    "log_omega[1]", "someomega[1]", "omega[1,2]"
  )

  expect_equal(
    JAGS_indexed_parameter_columns(columns, "omega"),
    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )

  samples <- matrix(
    seq_len(16),
    nrow = 2,
    dimnames = list(NULL, columns)
  )
  omega <- JAGS_indexed_parameter_matrix(samples, "omega")

  expect_equal(colnames(omega), c("omega[1]", "omega[2]", "omega[10]"))
  expect_equal(omega[, "omega[2]"], samples[, "omega[2]"])
})


test_that("JAGS indexed parameter helpers escape regex metacharacters", {

  samples <- data.frame(
    "omega.beta[1]" = 1:2,
    "omega.beta[2]" = 3:4,
    "omegaXbeta[1]" = 11:12,
    "omega+beta[1]" = 5:6,
    "theta(1)[2]" = 7:8,
    "theta(1)[1]" = 9:10,
    check.names = FALSE
  )

  expect_equal(
    colnames(JAGS_indexed_parameter_matrix(samples, "omega.beta")),
    c("omega.beta[1]", "omega.beta[2]")
  )
  expect_true(is.matrix(JAGS_indexed_parameter_matrix(samples, "omega.beta")))
  expect_equal(
    colnames(JAGS_indexed_parameter_matrix(samples, "omega+beta")),
    "omega+beta[1]"
  )
  expect_equal(
    colnames(JAGS_indexed_parameter_matrix(samples, "theta(1)")),
    c("theta(1)[1]", "theta(1)[2]")
  )
  expect_equal(JAGS_regex_escape("omega.beta[1]+x"), "omega\\.beta\\[1\\]\\+x")
  expect_equal(JAGS_regex_escape(character()), character())
  expect_equal(JAGS_indexed_parameter_columns(character(), "omega"), logical())
  expect_error(JAGS_indexed_parameter_columns("[1]", ""), "must not be empty")
})


test_that("JAGS indexed parameter helpers handle missing columns", {

  samples <- matrix(
    seq_len(4),
    nrow = 2,
    dimnames = list(NULL, c("theta", "phi[1]"))
  )

  expect_null(JAGS_indexed_parameter_matrix(samples, "omega"))

  empty <- JAGS_indexed_parameter_matrix(
    samples,
    "omega",
    drop_missing = FALSE
  )
  expect_equal(dim(empty), c(2L, 0L))
})


test_that("JAGS indexed parameter vector extraction sorts values", {

  row <- c("omega[10]" = 10, "omega[1]" = 1, other = 99, "omega[2]" = 2)

  out <- JAGS_indexed_parameter_vector(row, "omega")
  expect_equal(names(out), c("omega[1]", "omega[2]", "omega[10]"))
  expect_equal(unname(out), c(1, 2, 10))

  row_df <- data.frame(
    "omega[2]" = 2,
    "omega[1]" = 1,
    check.names = FALSE
  )
  expect_equal(
    unname(JAGS_indexed_parameter_vector(row_df, "omega")),
    c(1, 2)
  )
})


test_that("JAGS indexed parameter matrix supports coda objects", {

  samples <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(NULL, c("omega[2]", "omega[1]", "theta", "omega[10]"))
  )
  mcmc <- coda::mcmc(samples)
  mcmc_list <- coda::mcmc.list(coda::mcmc(samples), coda::mcmc(samples))

  expect_equal(
    colnames(JAGS_indexed_parameter_matrix(mcmc, "omega")),
    c("omega[1]", "omega[2]", "omega[10]")
  )
  expect_equal(
    colnames(JAGS_indexed_parameter_matrix(mcmc_list, "omega")),
    c("omega[1]", "omega[2]", "omega[10]")
  )
})
