skip_if_not_test_profile("unit")

cache_helper <- testthat::test_path(
  "..", "..", "vignettes", "precomputed-vignette-cache.R"
)
skip_if_not(
  file.exists(cache_helper),
  paste0(
    "Repository vignette cache sources are not available in this ",
    "installed-package test context."
  )
)

source(cache_helper, local = TRUE)

.precomputed_test_hash <- function(value){
  paste(rep(value, 32L), collapse = "")
}

.precomputed_test_state <- function(){
  list(
    source_hashes = c(
      "DESCRIPTION" = .precomputed_test_hash("1"),
      "NAMESPACE" = .precomputed_test_hash("2"),
      "R/JAGS-fit.R" = .precomputed_test_hash("3"),
      "vignettes/precomputed-vignette-cache.R" =
        .precomputed_test_hash("4")
    ),
    r_version = "4.6.0",
    package_versions = c(
      BayesTools = "0.3.1.6",
      runjags = "2.2.2-5",
      rjags = "4-17",
      coda = "0.19-4.1",
      bridgesampling = "1.2-1"
    )
  )
}

.precomputed_test_objects <- function(vignette){
  schema <- precomputed_vignette_cache_schema(vignette)
  objects <- lapply(names(schema), function(name){
    structure(list(identifier = name), class = schema[[name]])
  })
  stats::setNames(objects, names(schema))
}

.precomputed_test_write_text <- function(path, lines, eol = "\n"){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  text <- paste(lines, collapse = eol)
  if(length(lines) > 0L){
    text <- paste0(text, eol)
  }
  writeBin(charToRaw(text), path)
}

.precomputed_test_project <- function(root, eol = "\n"){
  files <- list(
    "DESCRIPTION" = c("Package: BayesTools", "Version: 0.0.1"),
    "NAMESPACE" = "exportPattern(\"^[[:alpha:]]+\")",
    "R/a.R" = "a <- function() 1",
    "R/nested/b.r" = "b <- function() 2",
    "src/native.cpp" = "int value = 1;",
    "src/native.h" = "#define VALUE 1",
    "src/Makevars" = "generated and excluded from source builds",
    "src/Makevars.win" = "PKG_CPPFLAGS=-DVALUE=1",
    "src/native.o" = "not a source dependency",
    "vignettes/ComparisonR.Rmd" = "# Comparison",
    "vignettes/SpikeAndSlab.Rmd" = "# Spike and slab",
    "vignettes/precomputed-vignette-cache.R" = "# shared helper"
  )
  for(path in names(files)){
    .precomputed_test_write_text(
      file.path(root, path),
      files[[path]],
      eol = eol
    )
  }
  data_file <- file.path(root, "data", "kitchen_rolls.RData")
  dir.create(dirname(data_file), recursive = TRUE, showWarnings = FALSE)
  writeBin(as.raw(c(1L, 2L, 3L, 4L)), data_file)
  invisible(root)
}

test_that("supported vignette schemas are exact and ordered", {
  expect_identical(
    precomputed_vignette_cache_schema("ComparisonR"),
    c(
      ttest_model_H0 = "BayesTools_fit",
      ttest_model_Hp = "BayesTools_fit",
      marglik_model_H0 = "BayesTools_marglik",
      marglik_model_Hp = "BayesTools_marglik"
    )
  )
  expect_identical(
    precomputed_vignette_cache_schema("SpikeAndSlab"),
    c(
      M0 = "BayesTools_fit",
      M1 = "BayesTools_fit",
      marglik_model_H0 = "BayesTools_marglik",
      marglik_model_H1 = "BayesTools_marglik",
      MS = "BayesTools_fit"
    )
  )
  expect_error(
    precomputed_vignette_cache_schema("RandomEffects"),
    "must be one of: ComparisonR, SpikeAndSlab",
    fixed = TRUE
  )
})

test_that("cache envelope and manifest are exact and versioned", {
  cache_file <- file.path(withr::local_tempdir(), "ComparisonR.RDS")
  state <- .precomputed_test_state()
  objects <- .precomputed_test_objects("ComparisonR")

  expect_identical(
    write_precomputed_vignette_cache(
      objects,
      vignette = "ComparisonR",
      cache_file = cache_file,
      state = state
    ),
    cache_file
  )

  envelope <- readRDS(cache_file)
  expect_identical(names(envelope), c("manifest", "payloads"))
  expect_identical(
    names(envelope$manifest),
    c(
      "format", "version", "vignette", "object_schema",
      "source_hashes", "payload_hashes", "producer"
    )
  )
  expect_identical(
    envelope$manifest$format,
    "BayesTools.precomputed-vignette-cache"
  )
  expect_identical(envelope$manifest$version, 1L)
  expect_identical(envelope$manifest$vignette, "ComparisonR")
  expect_identical(
    envelope$manifest$object_schema,
    precomputed_vignette_cache_schema("ComparisonR")
  )
  expect_identical(
    names(envelope$manifest$producer),
    c("generated_at_utc", "r_version", "package_versions")
  )
  expect_identical(envelope$manifest$producer$r_version, state$r_version)
  expect_identical(
    envelope$manifest$producer$package_versions,
    state$package_versions
  )
  expect_true(all(vapply(envelope$payloads, is.raw, logical(1))))

  status <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file = cache_file,
    state = state,
    load = TRUE
  )
  expect_true(status$valid)
  expect_identical(status$objects, objects)
  expect_identical(
    load_precomputed_vignette_cache(
      "ComparisonR",
      cache_file = cache_file,
      state = state
    ),
    objects
  )
})

test_that("writing enforces exact names, order, and required classes", {
  cache_file <- file.path(withr::local_tempdir(), "SpikeAndSlab.RDS")
  state <- .precomputed_test_state()
  objects <- .precomputed_test_objects("SpikeAndSlab")

  expect_error(
    write_precomputed_vignette_cache(
      objects[-1L],
      "SpikeAndSlab",
      cache_file,
      state = state
    ),
    "objects are missing: M0",
    fixed = TRUE
  )

  extra <- c(objects, list(unexpected = structure(list(), class = "list")))
  expect_error(
    write_precomputed_vignette_cache(
      extra,
      "SpikeAndSlab",
      cache_file,
      state = state
    ),
    "objects are unexpected: unexpected",
    fixed = TRUE
  )

  reordered <- objects[rev(names(objects))]
  expect_error(
    write_precomputed_vignette_cache(
      reordered,
      "SpikeAndSlab",
      cache_file,
      state = state
    ),
    "objects are not in the required order",
    fixed = TRUE
  )

  wrong_class <- objects
  class(wrong_class$MS) <- "list"
  expect_error(
    write_precomputed_vignette_cache(
      wrong_class,
      "SpikeAndSlab",
      cache_file,
      state = state
    ),
    "objects have unexpected classes: MS",
    fixed = TRUE
  )
  expect_false(file.exists(cache_file))
})

test_that("validation reports missing, unreadable, and malformed caches clearly", {
  cache_file <- file.path(withr::local_tempdir(), "ComparisonR.RDS")
  state <- .precomputed_test_state()

  missing <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_false(missing$valid)
  expect_identical(missing$reason, "missing")
  expect_match(missing$error, "Missing precomputed ComparisonR")
  expect_error(
    load_precomputed_vignette_cache(
      "ComparisonR",
      cache_file,
      state = state
    ),
    "Missing precomputed ComparisonR",
    fixed = TRUE
  )

  writeBin(charToRaw("not an RDS file"), cache_file)
  unreadable <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_identical(unreadable$reason, "unreadable")
  expect_match(unreadable$error, "Could not read precomputed ComparisonR")

  saveRDS(list(manifest = list(), payloads = list(), extra = TRUE), cache_file)
  malformed <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_identical(malformed$reason, "invalid")
  expect_match(malformed$error, "unexpected envelope schema")
})

test_that("payload fingerprints detect cache corruption", {
  cache_file <- file.path(withr::local_tempdir(), "ComparisonR.RDS")
  state <- .precomputed_test_state()
  write_precomputed_vignette_cache(
    .precomputed_test_objects("ComparisonR"),
    "ComparisonR",
    cache_file,
    state = state
  )

  envelope <- readRDS(cache_file)
  envelope$payloads$ttest_model_H0[[length(
    envelope$payloads$ttest_model_H0
  )]] <- as.raw(0L)
  saveRDS(envelope, cache_file)

  status <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_false(status$valid)
  expect_identical(status$reason, "corrupt")
  expect_match(status$error, "modified payloads: ttest_model_H0")
  expect_error(
    load_precomputed_vignette_cache(
      "ComparisonR",
      cache_file,
      state = state
    ),
    "modified payloads: ttest_model_H0",
    fixed = TRUE
  )
})

test_that("validation rejects unsupported manifests and payload classes", {
  cache_file <- file.path(withr::local_tempdir(), "ComparisonR.RDS")
  state <- .precomputed_test_state()
  write_precomputed_vignette_cache(
    .precomputed_test_objects("ComparisonR"),
    "ComparisonR",
    cache_file,
    state = state
  )

  envelope <- readRDS(cache_file)
  unsupported <- envelope
  unsupported$manifest$version <- 2L
  saveRDS(unsupported, cache_file)
  version_status <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_identical(version_status$reason, "invalid")
  expect_match(version_status$error, "manifest version is unsupported")

  wrong_class <- structure(list(), class = "list")
  envelope$payloads$ttest_model_H0 <-
    serialize(wrong_class, NULL, version = 3)
  envelope$manifest$payload_hashes <-
    .precomputed_vignette_payload_hashes(envelope$payloads)
  saveRDS(envelope, cache_file)
  class_status <- validate_precomputed_vignette_cache(
    "ComparisonR",
    cache_file,
    state = state
  )
  expect_identical(class_status$reason, "corrupt")
  expect_match(
    class_status$error,
    "objects have unexpected classes: ttest_model_H0"
  )
})

test_that("validation distinguishes source, R, and package staleness", {
  cache_file <- file.path(withr::local_tempdir(), "SpikeAndSlab.RDS")
  state <- .precomputed_test_state()
  write_precomputed_vignette_cache(
    .precomputed_test_objects("SpikeAndSlab"),
    "SpikeAndSlab",
    cache_file,
    state = state
  )

  changed_source <- state
  changed_source$source_hashes[[1L]] <- .precomputed_test_hash("a")
  source_status <- validate_precomputed_vignette_cache(
    "SpikeAndSlab",
    cache_file,
    state = changed_source
  )
  expect_identical(source_status$reason, "stale")
  expect_identical(source_status$stale, "source files")

  changed_runtime <- state
  changed_runtime$r_version <- "4.6.1"
  changed_runtime$package_versions[["coda"]] <- "0.20-0"
  runtime_status <- validate_precomputed_vignette_cache(
    "SpikeAndSlab",
    cache_file,
    state = changed_runtime
  )
  expect_identical(runtime_status$reason, "stale")
  expect_identical(
    runtime_status$stale,
    c("R version", "package versions")
  )
})

test_that("source hashes are relative, complete, and line-ending stable", {
  lf_root <- file.path(withr::local_tempdir(), "lf")
  crlf_root <- file.path(withr::local_tempdir(), "crlf")
  .precomputed_test_project(lf_root, eol = "\n")
  .precomputed_test_project(crlf_root, eol = "\r\n")

  lf_hashes <-
    .precomputed_vignette_source_hashes("ComparisonR", lf_root)
  crlf_hashes <-
    .precomputed_vignette_source_hashes("ComparisonR", crlf_root)

  expect_identical(lf_hashes, crlf_hashes)
  expect_identical(
    names(lf_hashes),
    sort(c(
      "DESCRIPTION",
      "NAMESPACE",
      "data/kitchen_rolls.RData",
      "R/a.R",
      "R/nested/b.r",
      "src/Makevars.win",
      "src/native.cpp",
      "src/native.h",
      "vignettes/ComparisonR.Rmd",
      "vignettes/precomputed-vignette-cache.R"
    ))
  )
  expect_false("src/Makevars" %in% names(lf_hashes))
  expect_false("src/native.o" %in% names(lf_hashes))
  expect_true(all(grepl("^[[:xdigit:]]{32}$", lf_hashes)))
})

test_that("DESCRIPTION fingerprints ignore only build-time normalization", {
  description_root <- withr::local_tempdir()
  source_description <- file.path(description_root, "source")
  staged_description <- file.path(description_root, "staged")
  writeLines(c(
    "Package: BayesTools",
    "Version: 0.3.1.7",
    "Authors@R: person(\"A\", \"Developer\", role = c(\"aut\", \"cre\"))",
    "Depends:",
    "    R (>= 4.3.0),",
    "    stats"
  ), source_description, useBytes = TRUE)
  writeLines(c(
    "Package: BayesTools",
    "Version: 0.3.1.7",
    "Authors@R: person(\"A\", \"Developer\", role = c(\"aut\", \"cre\"))",
    "Depends: R (>= 4.3.0), stats",
    "Packaged: 2026-07-27 17:02:42 UTC; builder",
    "Author: A Developer [aut, cre]"
  ), staged_description, useBytes = TRUE)

  expect_identical(
    .precomputed_vignette_description_md5(source_description),
    .precomputed_vignette_description_md5(staged_description)
  )

  staged <- readLines(staged_description, warn = FALSE)
  staged[staged == "Version: 0.3.1.7"] <- "Version: 0.3.1.8"
  writeLines(staged, staged_description, useBytes = TRUE)
  expect_false(identical(
    .precomputed_vignette_description_md5(source_description),
    .precomputed_vignette_description_md5(staged_description)
  ))
})

test_that("current state records all required producer versions", {
  project_root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  state <- precomputed_vignette_cache_state(
    "ComparisonR",
    cache_file = file.path(
      project_root,
      "vignettes",
      "ComparisonR.RDS"
    ),
    project_root = project_root
  )

  expect_identical(
    names(state),
    c("source_hashes", "r_version", "package_versions")
  )
  expect_identical(
    names(state$package_versions),
    c("BayesTools", "runjags", "rjags", "coda", "bridgesampling")
  )
  expect_identical(state$r_version, as.character(getRversion()))
  expect_true(all(c(
    "DESCRIPTION",
    "NAMESPACE",
    "vignettes/ComparisonR.Rmd",
    "vignettes/precomputed-vignette-cache.R"
  ) %in% names(state$source_hashes)))
  expect_true(any(startsWith(names(state$source_hashes), "R/")))
  expect_true(any(startsWith(names(state$source_hashes), "src/")))
})

test_that("atomic replacement installs only a validated candidate", {
  cache_root <- withr::local_tempdir()
  cache_file <- file.path(cache_root, "SpikeAndSlab.RDS")
  state <- .precomputed_test_state()
  original <- .precomputed_test_objects("SpikeAndSlab")
  write_precomputed_vignette_cache(
    original,
    "SpikeAndSlab",
    cache_file,
    state = state
  )

  replacement <- original
  replacement$MS$identifier <- "replacement"
  write_precomputed_vignette_cache(
    replacement,
    "SpikeAndSlab",
    cache_file,
    state = state
  )

  loaded <- load_precomputed_vignette_cache(
    "SpikeAndSlab",
    cache_file,
    state = state
  )
  expect_identical(loaded, replacement)
  expect_identical(
    list.files(cache_root, all.files = TRUE, no.. = TRUE),
    "SpikeAndSlab.RDS"
  )
})
