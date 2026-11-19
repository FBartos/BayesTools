skip_if_not_test_profile("unit")

# The factor slice kernel header (inst/include/BayesTools/factor_slice_kernel.h)
# that BayesTools' FactorSlice sampler and RoBMA's block sampler build on. A toy
# translation unit that includes only the header and the standard library
# (native/factor_slice_kernel_toy.cpp) is compiled with R's compiler flags,
# samples a bivariate normal and runs the resolvability check. The draws of FactorSlice itself are certified in
# the 'fit-assertions' case.

.factor_slice_kernel_toy <- function(){

  skip_on_cran()
  skip_if_not_installed("pkgbuild")
  skip_if_not(pkgbuild::has_build_tools(), "No compiler for the toy translation unit.")
  include_dir <- system.file("include", package = "BayesTools")
  expect_true(file.exists(file.path(include_dir, "BayesTools", "factor_slice_kernel.h")))

  build_dir <- withr::local_tempdir(.local_envir = parent.frame())
  file.copy(testthat::test_path("native", "factor_slice_kernel_toy.cpp"), build_dir)
  withr::local_dir(build_dir)
  withr::local_envvar(PKG_CPPFLAGS = paste0('-I"', normalizePath(include_dir, winslash = "/"), '"'))
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "factor_slice_kernel_toy.cpp"),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status", exact = TRUE)
  log <- paste(output, collapse = "\n")
  expect_identical(if(is.null(status)) 0L else as.integer(status), 0L, info = log)
  expect_false(any(grepl("warning:", output, fixed = TRUE)), info = log)

  library_file <- file.path(build_dir, paste0("factor_slice_kernel_toy", .Platform$dynlib.ext))
  dll <- dyn.load(library_file)
  withr::defer(dyn.unload(library_file), envir = parent.frame())
  list(
    sample = function(n_adapt, n_iter, seed, mean, sd, correlation){
      out <- .C(getNativeSymbolInfo("bt_factor_slice_kernel_toy", dll),
                n_adapt = as.integer(n_adapt), n_iter = as.integer(n_iter), seed = as.integer(seed),
                mean = as.double(mean), sd = as.double(sd), correlation = as.double(correlation),
                draws = double(2L * n_iter), exhausted = integer(1L))
      list(draws = matrix(out$draws, ncol = 2L), exhausted = out$exhausted)
    },
    resolvable = function(lower, upper, x, mean, sd){
      out <- .C(getNativeSymbolInfo("bt_factor_slice_resolvable_toy", dll),
                lower = as.double(lower), upper = as.double(upper), x = as.double(x),
                mean = as.double(mean), sd = as.double(sd), resolvable = integer(1L), max_difference = double(1L),
                NAOK = TRUE)
      list(resolvable = out$resolvable == 1L, max_difference = out$max_difference)
    },
    exp_threshold = function(){
      out <- .C(getNativeSymbolInfo("bt_factor_slice_exp_threshold_toy", dll),
                finite_at = integer(1L), infinite_above = integer(1L))
      c(finite_at = out$finite_at == 1L, infinite_above = out$infinite_above == 1L)
    }
  )
}

test_that("the kernel header compiles alone and its slice updates sample a bivariate normal", {

  toy <- .factor_slice_kernel_toy()$sample

  # N(m, S) with sd (1, 3) and correlation 0.9; each moment is compared with
  # its exact value in units of its Monte Carlo standard error (sd / sqrt(ESS))
  m <- c(1, -2)
  s <- c(1, 3)
  r <- 0.9
  quantities <- list(
    mean_1 = list(series = function(x) x[, 1], target = m[1]),
    mean_2 = list(series = function(x) x[, 2], target = m[2]),
    var_1  = list(series = function(x) (x[, 1] - m[1])^2, target = s[1]^2),
    var_2  = list(series = function(x) (x[, 2] - m[2])^2, target = s[2]^2),
    cov    = list(series = function(x) (x[, 1] - m[1]) * (x[, 2] - m[2]), target = r * s[1] * s[2])
  )
  for(n_adapt in c(0L, 500L)){
    run <- toy(n_adapt = n_adapt, n_iter = 20000L, seed = 7L, mean = m, sd = s, correlation = r)
    expect_identical(run$exhausted, 0L)
    checks <- do.call(rbind, lapply(names(quantities), function(name){
      values <- quantities[[name]]$series(run$draws)
      ess <- coda::effectiveSize(values)
      mcse <- stats::sd(values) / sqrt(ess)
      data.frame(quantity = name, estimate = mean(values), target = quantities[[name]]$target,
                 ess = unname(ess), z = unname((mean(values) - quantities[[name]]$target) / mcse))
    }))
    info <- paste(c(paste("adapt =", n_adapt), utils::capture.output(print(checks, digits = 4))), collapse = "\n")
    expect_true(all(abs(checks$z) <= 4), info = info)
    expect_true(all(checks$ess >= 500), info = info)
    # the same seed gives the same chain
    expect_identical(toy(n_adapt = n_adapt, n_iter = 20000L, seed = 7L, mean = m, sd = s, correlation = r), run)
  }
})

test_that("the resolvability check separates members their coordinate can and cannot resolve", {

  toy <- .factor_slice_kernel_toy()

  # the back-transforms keep the original quotient exactly while exp() is finite
  expect_identical(toy$exp_threshold(), c(finite_at = TRUE, infinite_above = TRUE))

  # a standard normal (or a normal with sd 1e-5 at 0.5) at its mode: the
  # coordinate reaches a value to a relative precision of about |u| eps of its
  # distance from the anchoring bound, so bounds 1e9 or 1e6 away resolve it and
  # bounds 5e14 or 1e16 away (5e9 for sd 1e-5) do not
  cases <- list(
    list(lower = -1e9,  upper = 1e9,  x = 0,    mean = 0,   sd = 1,    resolvable = TRUE),
    list(lower = -1e6,  upper = Inf,  x = 0,    mean = 0,   sd = 1,    resolvable = TRUE),
    list(lower = -Inf,  upper = 1e6,  x = 0,    mean = 0,   sd = 1,    resolvable = TRUE),
    list(lower = 0,     upper = 1,    x = 0.3,  mean = 0.3, sd = 0.1,  resolvable = TRUE),
    list(lower = 0,     upper = 1e308, x = 0.08, mean = 0,  sd = 0.1,  resolvable = TRUE),
    list(lower = -Inf,  upper = Inf,  x = 0,    mean = 0,   sd = 1e-300, resolvable = TRUE),
    list(lower = -5e14, upper = 5e14, x = 0,    mean = 0,   sd = 1,    resolvable = FALSE),
    list(lower = -5e14, upper = 5e14, x = 0.3,  mean = 0,   sd = 1,    resolvable = FALSE),
    list(lower = -1e16, upper = Inf,  x = 0,    mean = 0,   sd = 1,    resolvable = FALSE),
    list(lower = -Inf,  upper = 1e16, x = 0,    mean = 0,   sd = 1,    resolvable = FALSE),
    list(lower = -1e10, upper = 1e10, x = 0.5,  mean = 0.5, sd = 1e-5, resolvable = FALSE)
  )
  results <- vapply(cases, function(case){
    toy$resolvable(case$lower, case$upper, case$x, case$mean, case$sd)$resolvable
  }, logical(1))
  expect_identical(results, vapply(cases, `[[`, logical(1), "resolvable"))
})
