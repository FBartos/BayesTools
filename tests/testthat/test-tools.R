context("Tools")


test_that("Check booleans", {

  # these should be allowed
  expect_null(check_bool(NULL,  "", allow_NULL = TRUE))
  expect_null(check_bool(TRUE,  ""))
  expect_null(check_bool(FALSE, ""))
  expect_null(check_bool(as.logical(stats::rbinom(5, 1, .5)), "", check_length = 0))
  expect_null(check_bool(c(FALSE, FALSE), "", check_length = 2))

  # these should fail
  expect_error(
    check_bool(as.matrix(as.logical(rbinom(5, 1, .5))), "test object"),
    "The 'test object' argument must be a logical vector."
  )
  expect_error(
    check_bool("string", "test object"),
    "The 'test object' argument must be a logical vector."
  )
  expect_error(
    check_bool(1, "test object"),
    "The 'test object' argument must be a logical vector."
  )
  expect_error(
    check_bool(list(TRUE), "test object"),
    "The 'test object' argument must be a logical vector."
  )
  expect_error(
    check_bool(TRUE, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  expect_error(
    check_bool(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})

test_that("Check strings", {

  # these should be allowed
  expect_null(check_char(NULL,  "", allow_NULL = TRUE))
  expect_null(check_char("string",  ""))
  expect_null(check_char(c("string", "string1"),  "", check_length = 0))
  expect_null(check_char(as.character(stats::rbinom(5, 1, .5)), "", check_length = 5))
  expect_null(check_char(c("string", "string1"),  "", check_length = 0, allow_values = c("string", "string1")))


  # these should fail
  expect_error(
    check_char(as.matrix(as.logical(as.character(5, 1, .5))), "test object"),
    "The 'test object' argument must be a character vector."
  )
  expect_error(
    check_char(TRUE, "test object"),
    "The 'test object' argument must be a character vector."
  )
  expect_error(
    check_char(1, "test object"),
    "The 'test object' argument must be a character vector."
  )
  expect_error(
    check_char(list("string"), "test object"),
    "The 'test object' argument must be a character vector."
  )
  expect_error(
    check_char("string", "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  expect_error(
    check_char(c("string", "string1"),  "test object", check_length = 0, allow_values = c("string")),
    "The 'string1' values are not recognized by the 'test object' argument."
  )
  expect_error(
    check_char(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})

test_that("Check reals", {

  # these should be allowed
  expect_null(check_real(NULL,  "", allow_NULL = TRUE))
  expect_null(check_real(pi,  ""))
  expect_null(check_real(c(pi, 2),  "", check_length = 0))
  expect_null(check_real(stats::rnorm(4, 1, .5), "", check_length = 4))
  expect_null(check_real(stats::rgamma(1, 1, 1), "", lower = 0))
  expect_null(check_real(stats::rbeta(1, 1, 1), "",  upper = 1))
  expect_null(check_real(c(0, 1), "", lower = 0, upper = 1, check_length = 2))

  # these should fail
  expect_error(
    check_real(as.matrix(stats::rnorm(4, 1, .5)), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_real(TRUE, "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_real("string", "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_real(list(3.2), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_real(1, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  expect_error(
    check_real(stats::rgamma(1, 1, 1), "test object", upper = 0),
    "The 'test object' must be equal or lower than 0."
  )
  expect_error(
    check_real(stats::rbeta(1, 1, 1), "test object",  lower = 1),
    "The 'test object' must be equal or higher than 1."
  )
  expect_error(
    check_real(0, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be higher than 0."
  )
  expect_error(
    check_real(1, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be lower than 1."
  )
  expect_error(
    check_real(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})

test_that("Check integers", {

  # these should be allowed
  expect_null(check_int(NULL,  "", allow_NULL = TRUE))
  expect_null(check_int(0,  ""))
  expect_null(check_int(c(-1, 2),  "", check_length = 0))
  expect_null(check_int(stats::rpois(4, 1), "", check_length = 4))
  expect_null(check_int(c(0, 5, 10), "", lower = 0, check_length = 3))
  expect_null(check_int(c(0, -2, 2), "", upper = 2, check_length = 3))
  expect_null(check_int(c(-3, -1), "", lower = -3, upper = -1, check_length = 2))

  # these should fail
  expect_error(
    check_int(as.matrix(stats::rpois(4, 1)), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_int(TRUE, "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_int("string", "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_int(list(3.2), "test object"),
    "The 'test object' argument must be a numeric vector."
  )
  expect_error(
    check_int(1, "test object", check_length = 2),
    "The 'test object' argument must have length '2'."
  )
  expect_error(
    check_int(1, "test object", upper = 0),
    "The 'test object' must be equal or lower than 0."
  )
  expect_error(
    check_int(0, "test object",  lower = 1),
    "The 'test object' must be equal or higher than 1."
  )
  expect_error(
    check_int(0, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be higher than 0."
  )
  expect_error(
    check_int(1, "test object", lower = 0, upper = 1, allow_bound = FALSE),
    "The 'test object' must be lower than 1."
  )
  expect_error(
    check_int(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})

test_that("Check lists", {

  # these should be allowed
  expect_null(check_list(NULL,  "", allow_NULL = TRUE))
  expect_null(check_list(list(), "", allow_NULL = TRUE))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), ""))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_names = c("a", "b")))
  expect_null(check_list(list("a" = c("a", "b")), "", check_names = c("a", "b")))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_names = c("a", "b"), all_objects = TRUE))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1), "", check_length = 2))
  expect_null(check_list(list("a" = c("a", "b"), "b" = 1, "c" = c("a", "b")), "", check_names = c("a", "b"), all_objects = TRUE, allow_other = TRUE))

  # these should fail
  expect_error(
    check_list("string", "test object"),
    "The 'test object' argument must be a list."
  )
  expect_error(
    check_list(1, "test object"),
    "The 'test object' argument must be a list."
  )
  expect_error(
    check_list(TRUE, "test object"),
    "The 'test object' argument must be a list."
  )
  expect_error(
    check_list(list(), "test object", check_length = 2),
    "The 'test object' argument cannot be NULL."
  )
  expect_error(
    check_list(list("c" = c("a", "b")), "test object", check_names = c("a", "b")),
    "The 'c' objects are not recognized by the 'test object' argument."
  )
  expect_error(
    check_list(list("a" = c("a", "b")), "test object", check_names = c("a", "b"), all_objects = TRUE),
    "The 'b' objects are missing in the 'test object' argument."
  )
  expect_error(
    check_list(NULL, "test object"),
    "The 'test object' argument cannot be NULL."
  )
})

test_that("Other tools",{

  expect_warning(.depreciate.transform_orthonormal(TRUE, FALSE),
                 "'transform_orthonormal' argument will be depreciated in favor of 'transform_factors' argument.")


  expect_error(.extract_stan(NULL), "'fit' must be an rstan fit")


  expect_null(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "jac" = function(x) exp(x)
  ), NULL, FALSE))

  expect_error(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "err" = function(x) exp(x)
  ), NULL, FALSE), "The 'jac' objects are missing in the 'transformation' argument.")

  expect_error(.check_transformation_input(transformation = 1, NULL, FALSE), "Uknown format of the 'transformation' argument.")

})
