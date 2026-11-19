skip_if_not_test_profile("unit")

# The factor slice kernel header (inst/include/BayesTools/factor_slice_kernel.h)
# that BayesTools' FactorSlice sampler and RoBMA's block sampler build on. A toy
# translation unit that includes only the header and the standard library
# (native/factor_slice_kernel_toy.cpp) is compiled with R's compiler flags and
# samples a bivariate normal. The draws of FactorSlice itself are certified in
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
  function(n_adapt, n_iter, seed, mean, sd, correlation){
    out <- .C(getNativeSymbolInfo("bt_factor_slice_kernel_toy", dll),
              n_adapt = as.integer(n_adapt), n_iter = as.integer(n_iter), seed = as.integer(seed),
              mean = as.double(mean), sd = as.double(sd), correlation = as.double(correlation),
              draws = double(2L * n_iter), exhausted = integer(1L))
    list(draws = matrix(out$draws, ncol = 2L), exhausted = out$exhausted)
  }
}

test_that("the kernel header compiles alone and its slice updates sample a bivariate normal", {

  toy <- .factor_slice_kernel_toy()

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
