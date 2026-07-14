test_that("loaded fit backend fingerprint is stable and well formed", {

  expect_true(exists(
    "value",
    envir    = .fit_backend_fingerprint_cache,
    inherits = FALSE
  ))
  fingerprint <- fit_backend_fingerprint()

  expect_type(fingerprint, "character")
  expect_length(fingerprint, 1L)
  expect_match(fingerprint, "^[[:xdigit:]]{32}$")
  expect_identical(fit_backend_fingerprint(), fingerprint)
})


test_that("fit backend fingerprint freezes the loaded implementation", {

  .clear_fit_backend_fingerprint_cache()
  on.exit({
    .clear_fit_backend_fingerprint_cache()
    .freeze_fit_backend_fingerprint()
  }, add = TRUE)
  loaded_fingerprint <- paste(rep("a", 32L), collapse = "")
  edited_fingerprint <- paste(rep("b", 32L), collapse = "")

  .freeze_fit_backend_fingerprint(loaded_fingerprint)
  .freeze_fit_backend_fingerprint(edited_fingerprint)

  expect_identical(fit_backend_fingerprint(), loaded_fingerprint)
})


test_that("fit backend fingerprint tracks R and native implementations", {

  package_root <- tempfile("BayesTools-backend-")
  dir.create(file.path(package_root, "R"), recursive = TRUE)
  dir.create(file.path(package_root, "src"), recursive = TRUE)
  on.exit(unlink(package_root, recursive = TRUE), add = TRUE)

  writeLines(
    c("Package: BayesTools", "Version: 1.0.0"),
    file.path(package_root, "DESCRIPTION"),
    useBytes = TRUE
  )
  writeLines(
    "export(fit_backend_fingerprint)",
    file.path(package_root, "NAMESPACE"),
    useBytes = TRUE
  )
  r_file   <- file.path(package_root, "R", "backend.R")
  cpp_file <- file.path(package_root, "src", "backend.cc")
  dll_file <- file.path(package_root, "src", "BayesTools.dll")
  writeLines("backend <- function() 1", r_file, useBytes = TRUE)
  writeLines("int backend() { return 1; }", cpp_file, useBytes = TRUE)
  writeBin(charToRaw("loaded implementation"), dll_file)

  files <- .fit_backend_files(package_root, dll_file)
  expect_true(all(c(
    normalizePath(r_file, winslash = "/"),
    normalizePath(cpp_file, winslash = "/"),
    normalizePath(dll_file, winslash = "/")
  ) %in% files))

  before <- .compute_fit_backend_fingerprint(package_root, dll_file)
  writeLines("backend <- function() 2", r_file, useBytes = TRUE)
  after_r <- .compute_fit_backend_fingerprint(package_root, dll_file)
  expect_false(identical(after_r, before))

  writeBin(charToRaw("changed loaded implementation"), dll_file)
  after_dll <- .compute_fit_backend_fingerprint(package_root, dll_file)
  expect_false(identical(after_dll, after_r))
})


test_that("installed backend fingerprint uses installed artifacts", {

  package_root <- tempfile("BayesTools-installed-")
  dir.create(file.path(package_root, "R"), recursive = TRUE)
  dir.create(file.path(package_root, "libs", "x64"), recursive = TRUE)
  on.exit(unlink(package_root, recursive = TRUE), add = TRUE)

  writeLines(
    c("Package: BayesTools", "Version: 1.0.0"),
    file.path(package_root, "DESCRIPTION"),
    useBytes = TRUE
  )
  writeLines(
    "useDynLib(BayesTools)",
    file.path(package_root, "NAMESPACE"),
    useBytes = TRUE
  )
  rdb_file <- file.path(package_root, "R", "BayesTools.rdb")
  dll_file <- file.path(package_root, "libs", "x64", "BayesTools.dll")
  writeBin(charToRaw("compiled R database"), rdb_file)
  writeBin(charToRaw("native library"), dll_file)

  files <- .fit_backend_files(package_root, character())
  expect_true(all(c(
    normalizePath(rdb_file, winslash = "/"),
    normalizePath(dll_file, winslash = "/")
  ) %in% files))
  expect_match(
    .compute_fit_backend_fingerprint(package_root, character()),
    "^[[:xdigit:]]{32}$"
  )
})
