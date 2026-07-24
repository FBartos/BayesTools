skip_if_not_test_profile("unit")

test_that("native loader inventory matches registered call routines", {

  dll <- getLoadedDLLs()[["BayesTools"]]
  expect_false(is.null(dll))

  registered <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expect_setequal(.BayesTools_native_symbols(), names(registered))
  expect_true(all(vapply(
    .BayesTools_native_symbols(),
    is.loaded,
    logical(1),
    PACKAGE = "BayesTools"
  )))
})
