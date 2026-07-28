skip_if_not_test_profile("unit")

cache_helper <- testthat::test_path(
  "..", "..", "vignettes", "random-effects-vignette-cache.R"
)
skip_if_not(
  file.exists(cache_helper),
  paste0(
    "Repository vignette cache sources are not available in this ",
    "installed-package test context."
  )
)

source(
  cache_helper,
  local = TRUE
)

.random_effects_test_hash <- function(value){
  paste(rep(value, 32L), collapse = "")
}

.random_effects_test_sha256 <- function(value){
  paste(rep(value, 64L), collapse = "")
}

.random_effects_test_dependencies <- function(){
  list(
    generator_sources = c(
      "vignettes/RandomEffects.Rmd" = .random_effects_test_hash("1"),
      "vignettes/random-effects-vignette-cache.R" =
        .random_effects_test_hash("2")
    ),
    bayestools_sources = c(
      "DESCRIPTION" = .random_effects_test_hash("3"),
      "NAMESPACE" = .random_effects_test_hash("4"),
      "R/JAGS-fit.R" = .random_effects_test_hash("5"),
      "R/JAGS-formula.R" = .random_effects_test_hash("6")
    ),
    datasets = c(
      "lme4::sleepstudy" = .random_effects_test_hash("7"),
      "lme4::cake" = .random_effects_test_hash("8"),
      "lme4::Pastes" = .random_effects_test_hash("9"),
      "lme4::Penicillin" = .random_effects_test_hash("a")
    ),
    package_versions = c(
      BayesTools = "0.3.1.6",
      lme4 = "1.1-37",
      Matrix = "1.7-4",
      reformulas = "0.4.1",
      rstanarm = "2.32.1",
      rstan = "2.32.7",
      StanHeaders = "2.32.10",
      runjags = "2.2.2-5",
      rjags = "4-17",
      coda = "0.19-4.1"
    )
  )
}

.random_effects_test_producer <- function(
    backend = .random_effects_test_hash("b"),
    platform = "x86_64-test-platform",
    r_version = "4.6.0",
    jags_version = "4.3.2"){
  list(
    generated_at_utc = "2026-07-27T08:15:00Z",
    platform = platform,
    fit_backend_fingerprint = backend,
    loaded_package_versions = c(
      base = "4.6.0",
      BayesTools = "0.3.1.6"
    ),
    runtime = list(
      r_version = r_version,
      jags_version = jags_version,
      rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
      contrasts = c("contr.treatment", "contr.poly")
    )
  )
}

.random_effects_test_models <- function(){
  schema <- random_effects_vignette_cache_schema()
  models <- lapply(names(schema), function(name){
    structure(list(identifier = name), class = unname(schema[[name]]))
  })
  stats::setNames(models, names(schema))
}

.random_effects_test_implementation <- function(
    dependencies = .random_effects_test_dependencies(),
    backend = .random_effects_test_hash("b")){
  list(
    fit_backend_fingerprint = backend,
    bayestools_source_fingerprint =
      .random_effects_vignette_object_sha256(
        dependencies$bayestools_sources
      )
  )
}

.random_effects_test_write <- function(
    cache_file,
    models = .random_effects_test_models(),
    dependencies = .random_effects_test_dependencies(),
    producer = .random_effects_test_producer()){
  implementation <- .random_effects_test_implementation(
    dependencies,
    backend = producer$fit_backend_fingerprint
  )
  generation <- begin_random_effects_vignette_cache_regeneration(
    envir = new.env(parent = emptyenv()),
    cache_file = cache_file,
    dependency_state = dependencies,
    implementation_state = implementation
  )
  write_random_effects_vignette_cache(
    models = models,
    cache_file = cache_file,
    generation = generation,
    dependency_state = dependencies,
    implementation_state = implementation,
    producer = producer
  )
}

.random_effects_test_chunk <- function(lines, label){
  start <- grep(
    paste0("^```\\{r[[:space:]]+", label, "(?:,|\\})"),
    lines,
    perl = TRUE
  )
  if(length(start) != 1L){
    stop("Could not find exactly one vignette chunk: ", label, call. = FALSE)
  }
  relative_end <- which(lines[seq.int(start + 1L, length(lines))] == "```")
  if(length(relative_end) == 0L){
    stop("Could not find the end of vignette chunk: ", label, call. = FALSE)
  }
  end <- start + relative_end[[1L]]
  paste(lines[seq.int(start, end)], collapse = "\n")
}

test_that("RandomEffects cache schema is exact and ordered", {
  expected_classes <- c(
    stan_correlated = "stanreg",
    stan_correlated_scaled = "stanreg",
    fit_sleep = "BayesTools_fit",
    fit_sleep_independent = "BayesTools_fit",
    fit_cake_recipe = "BayesTools_fit",
    fit_cake_hcs = "BayesTools_fit",
    fit_cake_ar1 = "BayesTools_fit",
    fit_cake_har = "BayesTools_fit",
    fit_sleep_car = "BayesTools_fit",
    fit_nested = "BayesTools_fit",
    stan_crossed = "stanreg",
    fit_crossed_independent = "BayesTools_fit",
    fit_crossed_allocation = "BayesTools_fit"
  )
  expected_engines <- c(
    stan_correlated = "rstanarm",
    stan_correlated_scaled = "rstanarm",
    fit_sleep = "BayesTools/JAGS",
    fit_sleep_independent = "BayesTools/JAGS",
    fit_cake_recipe = "BayesTools/JAGS",
    fit_cake_hcs = "BayesTools/JAGS",
    fit_cake_ar1 = "BayesTools/JAGS",
    fit_cake_har = "BayesTools/JAGS",
    fit_sleep_car = "BayesTools/JAGS",
    fit_nested = "BayesTools/JAGS",
    stan_crossed = "rstanarm",
    fit_crossed_independent = "BayesTools/JAGS",
    fit_crossed_allocation = "BayesTools/JAGS"
  )
  expected_seeds <- c(
    stan_correlated = 11L,
    stan_correlated_scaled = 12L,
    fit_sleep = 1L,
    fit_sleep_independent = 13L,
    fit_cake_recipe = 34L,
    fit_cake_hcs = 31L,
    fit_cake_ar1 = 32L,
    fit_cake_har = 33L,
    fit_sleep_car = 14L,
    fit_nested = 41L,
    stan_crossed = 21L,
    fit_crossed_independent = 2L,
    fit_crossed_allocation = 3L
  )

  expect_identical(random_effects_vignette_cache_schema(), expected_classes)
  expect_identical(
    random_effects_vignette_cache_engines(),
    expected_engines
  )
  expect_identical(random_effects_vignette_cache_seeds(), expected_seeds)
  expect_identical(
    random_effects_vignette_cache_names(),
    names(expected_classes)
  )
  expect_identical(
    random_effects_vignette_cache_model_schema(),
    list(
      classes = expected_classes,
      engines = expected_engines,
      seeds = expected_seeds
    )
  )
})

test_that("writer installs one strict manifest envelope", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  models <- .random_effects_test_models()

  written <- .random_effects_test_write(
    cache_file = cache_file,
    models = models,
    dependencies = dependencies
  )
  envelope <- readRDS(cache_file)

  expect_true(written$valid)
  expect_identical(names(envelope), c("manifest", "models"))
  expect_identical(names(envelope$manifest), c(
    "format",
    "manifest_version",
    "cache_schema_version",
    "model_schema",
    "model_hashes",
    "dependencies",
    "compatibility_fingerprint",
    "producer",
    "generation_fingerprint"
  ))
  expect_identical(
    envelope$manifest$format,
    "BayesTools.RandomEffects.vignette-cache"
  )
  expect_identical(envelope$manifest$manifest_version, 1L)
  expect_identical(envelope$manifest$cache_schema_version, 1L)
  expect_identical(
    envelope$manifest$model_schema,
    random_effects_vignette_cache_model_schema()
  )
  expect_identical(
    envelope$manifest$model_hashes,
    .random_effects_vignette_model_hashes(models)
  )
  expect_identical(
    names(envelope$manifest$model_hashes),
    random_effects_vignette_cache_names()
  )
  expect_true(all(grepl(
    "^[[:xdigit:]]{64}$",
    envelope$manifest$model_hashes
  )))
  expect_identical(envelope$manifest$dependencies, dependencies)
  expect_match(
    envelope$manifest$compatibility_fingerprint,
    "^[[:xdigit:]]{32}$"
  )
  expect_match(
    envelope$manifest$generation_fingerprint,
    "^[[:xdigit:]]{64}$"
  )

  unloaded <- validate_random_effects_vignette_cache(
    cache_file,
    load = FALSE,
    dependency_state = dependencies
  )
  loaded <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_true(unloaded$valid)
  expect_null(unloaded$cache)
  expect_true(loaded$valid)
  expect_identical(loaded$cache, models)
})

test_that("writer canonicalizes models before payload hashing", {
  cache_file <- file.path(withr::local_tempdir(), "models", "RandomEffects.RDS")
  dependencies <- .random_effects_test_dependencies()
  implementation <- .random_effects_test_implementation(dependencies)
  producer <- .random_effects_test_producer()
  models <- .random_effects_test_models()
  live_environment <- new.env(parent = emptyenv())
  live_environment$value <- 1
  attr(models[[1L]], "live_environment") <- live_environment
  generation <- begin_random_effects_vignette_cache_regeneration(
    envir = new.env(parent = emptyenv()),
    dependency_state = dependencies,
    implementation_state = implementation
  )

  helper_environment <- environment(write_random_effects_vignette_cache)
  original_atomic_write <- get(
    ".random_effects_vignette_atomic_write",
    envir = helper_environment,
    inherits = FALSE
  )
  on.exit(assign(
    ".random_effects_vignette_atomic_write",
    original_atomic_write,
    envir = helper_environment
  ), add = TRUE)

  captured_envelope <- NULL
  assign(
    ".random_effects_vignette_atomic_write",
    function(envelope, cache_file, dependency_state, ...){
      captured_envelope <<- envelope
      list(valid = TRUE)
    },
    envir = helper_environment
  )

  written <- write_random_effects_vignette_cache(
    models = models,
    cache_file = cache_file,
    generation = generation,
    dependency_state = dependencies,
    implementation_state = implementation,
    producer = producer
  )

  expect_true(written$valid)
  expect_false(identical(
    attr(captured_envelope$models[[1L]], "live_environment"),
    live_environment
  ))
  expect_identical(
    captured_envelope$manifest$model_hashes,
    .random_effects_vignette_model_hashes(captured_envelope$models)
  )

  dir.create(dirname(cache_file), recursive = TRUE)
  saveRDS(captured_envelope, cache_file, version = 3)
  expect_true(validate_random_effects_vignette_cache(
    cache_file,
    dependency_state = dependencies
  )$valid)
})

test_that("canonicalization rejects payloads that do not stabilize", {
  models <- .random_effects_test_models()
  helper_environment <- environment(
    .random_effects_vignette_canonicalize_models
  )
  original_roundtrip <- get(
    ".random_effects_vignette_serialization_roundtrip",
    envir = helper_environment,
    inherits = FALSE
  )
  on.exit(assign(
    ".random_effects_vignette_serialization_roundtrip",
    original_roundtrip,
    envir = helper_environment
  ), add = TRUE)

  roundtrip_count <- 0L
  assign(
    ".random_effects_vignette_serialization_roundtrip",
    function(value){
      roundtrip_count <<- roundtrip_count + 1L
      value[[1L]]$serialization_generation <- roundtrip_count
      value
    },
    envir = helper_environment
  )

  expect_error(
    .random_effects_vignette_canonicalize_models(models),
    "fitted model payloads do not have a stable serialized representation",
    fixed = TRUE
  )
})

test_that("missing and legacy caches are graceful but never loadable", {
  cache_root <- withr::local_tempdir()
  cache_file <- file.path(cache_root, "models", "RandomEffects.RDS")
  impossible_root <- file.path(cache_root, "does-not-exist")

  missing <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    project_root = impossible_root
  )
  expect_false(missing$file_exists)
  expect_false(missing$valid)
  expect_null(missing$cache)
  expect_match(
    format_random_effects_vignette_cache_error(missing),
    "^Missing precomputed"
  )

  dir.create(dirname(cache_file), recursive = TRUE)
  saveRDS(.random_effects_test_models(), cache_file)
  legacy <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    project_root = impossible_root
  )
  expect_true(legacy$file_exists)
  expect_true(legacy$legacy_cache)
  expect_false(legacy$valid)
  expect_null(legacy$cache)
  expect_match(
    format_random_effects_vignette_cache_error(legacy),
    "legacy format without a manifest; regenerate it"
  )
  expect_error(
    stop_if_invalid_random_effects_vignette_cache(
      cache_file,
      project_root = impossible_root
    ),
    "legacy format without a manifest; regenerate it",
    fixed = TRUE
  )
})

test_that("malformed and corrupt envelopes return cache NULL", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  dir.create(dirname(cache_file), recursive = TRUE)

  writeBin(charToRaw("not an RDS file"), cache_file)
  corrupt <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(corrupt$valid)
  expect_true(nzchar(corrupt$read_error))
  expect_null(corrupt$cache)

  saveRDS(list(models = .random_effects_test_models()), cache_file)
  malformed <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(malformed$envelope_valid)
  expect_false(malformed$valid)
  expect_null(malformed$cache)

  .random_effects_test_write(cache_file, dependencies = dependencies)
  envelope <- readRDS(cache_file)
  envelope$manifest$compatibility_fingerprint <-
    .random_effects_test_hash("f")
  saveRDS(envelope, cache_file)
  tampered <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(tampered$manifest_valid)
  expect_false(tampered$valid)
  expect_null(tampered$cache)

  envelope <- readRDS(cache_file)
  envelope$manifest$dependencies$generator_sources[[1L]] <- NA_character_
  saveRDS(envelope, cache_file)
  invalid_hash <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(invalid_hash$manifest_valid)
  expect_false(invalid_hash$valid)
  expect_null(invalid_hash$cache)
})

test_that("dependency changes make a well-formed cache stale", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  .random_effects_test_write(cache_file, dependencies = dependencies)

  current_dependencies <- dependencies
  current_dependencies$generator_sources[[1L]] <-
    .random_effects_test_hash("c")
  stale <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = current_dependencies
  )

  expect_true(stale$manifest_valid)
  expect_identical(stale$stale_dependencies, "generator_sources")
  expect_false(stale$valid)
  expect_null(stale$cache)
  expect_match(
    format_random_effects_vignette_cache_error(stale),
    "changed dependencies: generator_sources"
  )
})

test_that("model payloads require exact names, order, and classes", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  .random_effects_test_write(cache_file, dependencies = dependencies)
  valid_envelope <- readRDS(cache_file)
  valid_models <- valid_envelope$models

  mutations <- list(
    missing = valid_models[-1L],
    extra = c(
      valid_models,
      unexpected = list(structure(list(), class = "BayesTools_fit"))
    ),
    reordered = valid_models[rev(seq_along(valid_models))],
    unnamed = unname(valid_models),
    wrong_class = within(valid_models, {
      fit_sleep <- structure(list(), class = "unexpected_fit")
    })
  )

  for(mutation in names(mutations)){
    envelope <- valid_envelope
    envelope$models <- mutations[[mutation]]
    saveRDS(envelope, cache_file)
    status <- validate_random_effects_vignette_cache(
      cache_file,
      load = TRUE,
      dependency_state = dependencies
    )
    expect_false(status$valid, info = mutation)
    expect_null(status$cache, info = mutation)
  }

  envelope <- valid_envelope
  envelope$models <- 1
  saveRDS(envelope, cache_file)
  non_list <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(non_list$is_list)
  expect_false(non_list$valid)
  expect_null(non_list$cache)
})

test_that("same-class model replacement fails payload integrity", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  .random_effects_test_write(cache_file, dependencies = dependencies)
  envelope <- readRDS(cache_file)

  envelope$models$fit_sleep <- structure(
    list(identifier = "same class, different payload"),
    class = "BayesTools_fit"
  )
  saveRDS(envelope, cache_file)
  replaced <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )

  expect_true(replaced$manifest_valid)
  expect_true(replaced$is_list)
  expect_length(replaced$invalid_objects, 0L)
  expect_false(replaced$payload_valid)
  expect_identical(replaced$invalid_payloads, "fit_sleep")
  expect_false(replaced$valid)
  expect_null(replaced$cache)
  expect_match(
    format_random_effects_vignette_cache_error(replaced),
    "modified payloads: fit_sleep"
  )
})

test_that("backend and runtime are provenance, not host compatibility", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  .random_effects_test_write(cache_file, dependencies = dependencies)
  envelope <- readRDS(cache_file)
  original_compatibility <- envelope$manifest$compatibility_fingerprint

  envelope$manifest$producer <- .random_effects_test_producer(
    backend = .random_effects_test_hash("d"),
    platform = "aarch64-another-platform",
    r_version = "4.7.1",
    jags_version = "4.4.0"
  )
  saveRDS(envelope, cache_file)
  tampered <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_false(tampered$manifest_valid)
  expect_match(
    tampered$manifest_error,
    "generation fingerprint does not match"
  )
  expect_null(tampered$cache)

  envelope$manifest$generation_fingerprint <-
    .random_effects_vignette_generation_fingerprint(
      envelope$manifest$compatibility_fingerprint,
      envelope$manifest$model_hashes,
      envelope$manifest$producer
    )
  saveRDS(envelope, cache_file)

  status <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_true(status$valid)
  expect_identical(
    envelope$manifest$compatibility_fingerprint,
    original_compatibility
  )
  expect_identical(status$cache, envelope$models)
})

test_that("regeneration removes cached fits and guards dependency drift", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  implementation <- .random_effects_test_implementation(dependencies)
  models <- .random_effects_test_models()
  .random_effects_test_write(
    cache_file,
    models = models,
    dependencies = dependencies
  )
  original_md5 <- unname(tools::md5sum(cache_file))

  regeneration_environment <- list2env(
    c(models, untouched = list(TRUE)),
    envir = new.env(parent = emptyenv())
  )
  generation <- begin_random_effects_vignette_cache_regeneration(
    envir = regeneration_environment,
    cache_file = cache_file,
    dependency_state = dependencies,
    implementation_state = implementation
  )
  expect_false(any(vapply(
    random_effects_vignette_cache_names(),
    exists,
    logical(1),
    envir = regeneration_environment,
    inherits = FALSE
  )))
  expect_true(exists(
    "untouched",
    envir = regeneration_environment,
    inherits = FALSE
  ))
  expect_error(
    begin_random_effects_vignette_cache_regeneration(
      envir = list(),
      dependency_state = dependencies
    ),
    "'envir' must be an environment",
    fixed = TRUE
  )

  expect_error(
    write_random_effects_vignette_cache(
      models = models[-1L],
      cache_file = cache_file,
      generation = generation,
      dependency_state = dependencies,
      implementation_state = implementation,
      producer = .random_effects_test_producer()
    ),
    "missing objects"
  )
  expect_identical(unname(tools::md5sum(cache_file)), original_md5)

  changed_dependencies <- dependencies
  changed_dependencies$package_versions[["BayesTools"]] <- "0.3.1.7"
  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = generation,
      dependency_state = changed_dependencies,
      implementation_state =
        .random_effects_test_implementation(changed_dependencies),
      producer = .random_effects_test_producer()
    ),
    "dependencies changed during regeneration: package_versions"
  )
  expect_identical(unname(tools::md5sum(cache_file)), original_md5)

  invalid_generation <- generation
  invalid_generation$dependencies$datasets[[1L]] <- NA_character_
  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = invalid_generation,
      dependency_state = dependencies,
      implementation_state = implementation,
      producer = .random_effects_test_producer()
    ),
    "regeneration checkpoint is invalid"
  )
  expect_identical(unname(tools::md5sum(cache_file)), original_md5)
})

test_that("loaded BayesTools implementation must match the current source", {
  dependencies <- .random_effects_test_dependencies()
  project_root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  current <- .random_effects_vignette_current_implementation(
    dependencies = dependencies,
    cache_file = file.path(project_root, "vignettes", "RandomEffects.RDS"),
    project_root = project_root
  )
  expect_null(
    .random_effects_vignette_implementation_state_error(
      current,
      dependencies
    )
  )
  expect_identical(
    current$fit_backend_fingerprint,
    BayesTools::fit_backend_fingerprint()
  )

  unrelated_root <- withr::local_tempdir()
  writeLines(
    c("Package: BayesTools", "Version: 0.0.0"),
    file.path(unrelated_root, "DESCRIPTION"),
    useBytes = TRUE
  )
  expect_error(
    .random_effects_vignette_current_implementation(
      dependencies = dependencies,
      project_root = unrelated_root
    ),
    "loaded BayesTools namespace is not the current source tree"
  )

  clear_fingerprint <- getFromNamespace(
    ".clear_fit_backend_fingerprint_cache",
    "BayesTools"
  )
  freeze_fingerprint <- getFromNamespace(
    ".freeze_fit_backend_fingerprint",
    "BayesTools"
  )
  original_fingerprint <- BayesTools::fit_backend_fingerprint()
  on.exit({
    clear_fingerprint()
    freeze_fingerprint(original_fingerprint)
  }, add = TRUE)
  clear_fingerprint()
  freeze_fingerprint(.random_effects_test_hash("f"))
  expect_error(
    .random_effects_vignette_current_implementation(
      dependencies = dependencies,
      cache_file = file.path(project_root, "vignettes", "RandomEffects.RDS"),
      project_root = project_root
    ),
    "changed after the namespace was loaded"
  )
})

test_that("writer rejects stale implementation checkpoints and reuse", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  models <- .random_effects_test_models()
  implementation <- .random_effects_test_implementation(dependencies)
  generation <- begin_random_effects_vignette_cache_regeneration(
    envir = new.env(parent = emptyenv()),
    dependency_state = dependencies,
    implementation_state = implementation
  )

  changed_implementation <- implementation
  changed_implementation$fit_backend_fingerprint <-
    .random_effects_test_hash("c")
  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = generation,
      dependency_state = dependencies,
      implementation_state = changed_implementation,
      producer = .random_effects_test_producer(
        backend = changed_implementation$fit_backend_fingerprint
      )
    ),
    "verified BayesTools implementation changed"
  )
  expect_false(file.exists(cache_file))

  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = generation,
      dependency_state = dependencies,
      implementation_state = implementation,
      producer = .random_effects_test_producer(
        backend = .random_effects_test_hash("d")
      )
    ),
    "producer backend does not match"
  )
  expect_false(file.exists(cache_file))

  written <- write_random_effects_vignette_cache(
    models = models,
    cache_file = cache_file,
    generation = generation,
    dependency_state = dependencies,
    implementation_state = implementation,
    producer = .random_effects_test_producer()
  )
  expect_true(written$valid)
  expect_true(generation$checkpoint_state$used)
  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = generation,
      dependency_state = dependencies,
      implementation_state = implementation,
      producer = .random_effects_test_producer()
    ),
    "checkpoint was already used"
  )
})

test_that("atomic replacement restores the previous cache on validation failure", {
  cache_file <- file.path(
    withr::local_tempdir(),
    "models",
    "RandomEffects.RDS"
  )
  dependencies <- .random_effects_test_dependencies()
  models <- .random_effects_test_models()
  .random_effects_test_write(
    cache_file,
    models = models,
    dependencies = dependencies
  )
  original_md5 <- unname(tools::md5sum(cache_file))

  helper_environment <- environment(.random_effects_vignette_atomic_write)
  original_validator <- get(
    "validate_random_effects_vignette_cache",
    envir = helper_environment,
    inherits = FALSE
  )
  on.exit(assign(
    "validate_random_effects_vignette_cache",
    original_validator,
    envir = helper_environment
  ), add = TRUE)
  validation_count <- 0L
  failing_validator <- function(...){
    validation_count <<- validation_count + 1L
    status <- original_validator(...)
    if(validation_count == 2L){
      status$manifest_valid <- FALSE
      status$manifest_error <- "forced installed validation failure"
      status$valid <- FALSE
      status$cache <- NULL
    }
    status
  }
  assign(
    "validate_random_effects_vignette_cache",
    failing_validator,
    envir = helper_environment
  )

  replacement_producer <- .random_effects_test_producer(
    backend = .random_effects_test_hash("e")
  )
  replacement_implementation <- .random_effects_test_implementation(
    dependencies,
    backend = replacement_producer$fit_backend_fingerprint
  )
  generation <- begin_random_effects_vignette_cache_regeneration(
    envir = new.env(parent = emptyenv()),
    dependency_state = dependencies,
    implementation_state = replacement_implementation
  )
  expect_error(
    write_random_effects_vignette_cache(
      models = models,
      cache_file = cache_file,
      generation = generation,
      dependency_state = dependencies,
      implementation_state = replacement_implementation,
      producer = replacement_producer
    ),
    "Installed RandomEffects cache failed validation"
  )
  assign(
    "validate_random_effects_vignette_cache",
    original_validator,
    envir = helper_environment
  )

  expect_identical(unname(tools::md5sum(cache_file)), original_md5)
  restored <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_true(restored$valid)
  expect_identical(restored$cache, models)
  expect_length(
    list.files(
      dirname(cache_file),
      pattern = "^\\.RandomEffects\\.RDS\\.(new|backup|journal.*)$",
      all.files = TRUE
    ),
    0L
  )
})

test_that("transaction journal recovers interrupted Windows replacement", {
  cache_root <- withr::local_tempdir()
  cache_file <- file.path(cache_root, "models", "RandomEffects.RDS")
  candidate_file <- file.path(cache_root, "candidate", "RandomEffects.RDS")
  dependencies <- .random_effects_test_dependencies()
  old_models <- .random_effects_test_models()
  new_models <- old_models
  new_models$fit_sleep <- structure(
    list(identifier = "new candidate"),
    class = "BayesTools_fit"
  )
  .random_effects_test_write(
    cache_file,
    models = old_models,
    dependencies = dependencies
  )
  .random_effects_test_write(
    candidate_file,
    models = new_models,
    dependencies = dependencies
  )

  paths <- .random_effects_vignette_transaction_paths(cache_file)
  old_sha256 <- .random_effects_vignette_sha256_file(cache_file)
  new_sha256 <- .random_effects_vignette_sha256_file(candidate_file)
  expect_true(file.copy(cache_file, paths$backup))
  expect_true(file.copy(candidate_file, paths$new))
  journal <- list(
    format = "BayesTools.RandomEffects.vignette-cache-transaction",
    version = 1L,
    had_cache = TRUE,
    previous_sha256 = old_sha256,
    new_sha256 = new_sha256
  )
  saveRDS(journal, paths$journal, version = 3)
  .random_effects_vignette_remove_file(
    cache_file,
    "simulated interrupted canonical cache"
  )

  recovered <- recover_random_effects_vignette_cache(
    cache_file,
    .platform = "windows"
  )
  expect_identical(recovered, "previous-restored")
  expect_identical(
    .random_effects_vignette_sha256_file(cache_file),
    old_sha256
  )
  expect_true(validate_random_effects_vignette_cache(
    cache_file,
    dependency_state = dependencies
  )$valid)
  expect_false(any(vapply(paths, file.exists, logical(1))))

  expect_true(file.copy(cache_file, paths$backup))
  expect_true(file.copy(candidate_file, paths$new))
  saveRDS(journal, paths$journal, version = 3)
  .random_effects_vignette_remove_file(
    cache_file,
    "canonical cache before simulated completed install"
  )
  expect_true(file.copy(candidate_file, cache_file))

  completed <- recover_random_effects_vignette_cache(
    cache_file,
    .platform = "windows"
  )
  expect_identical(completed, "canonical-complete")
  expect_identical(
    .random_effects_vignette_sha256_file(cache_file),
    new_sha256
  )
  loaded <- validate_random_effects_vignette_cache(
    cache_file,
    load = TRUE,
    dependency_state = dependencies
  )
  expect_true(loaded$valid)
  expect_identical(
    loaded$cache$fit_sleep$identifier,
    "new candidate"
  )
  expect_false(any(vapply(paths, file.exists, logical(1))))
})

test_that("recovery fails closed when backup does not match its journal", {
  cache_root <- withr::local_tempdir()
  cache_file <- file.path(cache_root, "models", "RandomEffects.RDS")
  dependencies <- .random_effects_test_dependencies()
  .random_effects_test_write(cache_file, dependencies = dependencies)
  paths <- .random_effects_vignette_transaction_paths(cache_file)
  previous_sha256 <- .random_effects_vignette_sha256_file(cache_file)
  expect_true(file.copy(cache_file, paths$backup))
  writeBin(charToRaw("modified recovery backup"), paths$backup)
  saveRDS(
    list(
      format = "BayesTools.RandomEffects.vignette-cache-transaction",
      version = 1L,
      had_cache = TRUE,
      previous_sha256 = previous_sha256,
      new_sha256 = .random_effects_test_sha256("e")
    ),
    paths$journal,
    version = 3
  )
  .random_effects_vignette_remove_file(
    cache_file,
    "canonical cache for corrupt-backup recovery test"
  )

  expect_error(
    recover_random_effects_vignette_cache(
      cache_file,
      .platform = "windows"
    ),
    "backup does not match the transaction journal"
  )
  expect_false(file.exists(cache_file))
  expect_true(file.exists(paths$backup))
  expect_true(file.exists(paths$journal))
})

test_that("POSIX replacement keeps canonical cache until direct rename", {
  helper_environment <- environment(.random_effects_vignette_replace_file)
  original_remove <- get(
    ".random_effects_vignette_remove_file",
    envir = helper_environment,
    inherits = FALSE
  )
  original_rename <- get(
    ".random_effects_vignette_rename_file",
    envir = helper_environment,
    inherits = FALSE
  )
  on.exit({
    assign(
      ".random_effects_vignette_remove_file",
      original_remove,
      envir = helper_environment
    )
    assign(
      ".random_effects_vignette_rename_file",
      original_rename,
      envir = helper_environment
    )
  }, add = TRUE)
  removed <- character()
  renamed <- list()
  assign(
    ".random_effects_vignette_remove_file",
    function(path, description){
      removed <<- c(removed, path)
      invisible(TRUE)
    },
    envir = helper_environment
  )
  assign(
    ".random_effects_vignette_rename_file",
    function(from, to, description){
      renamed <<- c(renamed, list(c(from = from, to = to)))
      invisible(TRUE)
    },
    envir = helper_environment
  )

  replacement_root <- withr::local_tempdir()
  prepared_file <- file.path(replacement_root, "prepared-new")
  canonical_file <- file.path(replacement_root, "canonical-cache")
  writeLines("new", prepared_file, useBytes = TRUE)
  writeLines("old", canonical_file, useBytes = TRUE)
  .random_effects_vignette_replace_file(
    prepared_file,
    canonical_file,
    "unix",
    "test POSIX replacement"
  )
  expect_length(removed, 0L)
  expect_identical(
    renamed,
    list(c(from = prepared_file, to = canonical_file))
  )

  assign(
    ".random_effects_vignette_remove_file",
    original_remove,
    envir = helper_environment
  )
  assign(
    ".random_effects_vignette_rename_file",
    original_rename,
    envir = helper_environment
  )
  atomic_source <- paste(
    deparse(.random_effects_vignette_atomic_write),
    collapse = "\n"
  )
  expect_false(grepl(
    "rename_file\\([[:space:]]*cache_file",
    atomic_source
  ))
})

test_that("transaction filesystem failures are checked", {
  helper_environment <- environment(.random_effects_vignette_remove_file)
  local_binding <- function(name, value){
    test_environment <- parent.frame()
    existed <- exists(name, envir = helper_environment, inherits = FALSE)
    original <- if(existed){
      get(name, envir = helper_environment, inherits = FALSE)
    }else{
      NULL
    }
    assign(name, value, envir = helper_environment)
    restore <- function(){
      if(existed){
        assign(name, original, envir = helper_environment)
      }else if(exists(name, envir = helper_environment, inherits = FALSE)){
        rm(list = name, envir = helper_environment)
      }
    }
    withr::defer(restore(), envir = test_environment)
  }

  transaction_root <- withr::local_tempdir()
  removable <- file.path(transaction_root, "removable")
  writeLines("canonical", removable, useBytes = TRUE)
  local_binding("unlink", function(...) 1L)
  expect_error(
    .random_effects_vignette_remove_file(
      removable,
      "the simulated locked cache"
    ),
    "Could not remove the simulated locked cache"
  )
  expect_true(file.exists(removable))

  prepared <- file.path(transaction_root, "prepared")
  destination <- file.path(transaction_root, "destination")
  writeLines("candidate", prepared, useBytes = TRUE)
  local_binding("file.rename", function(...) FALSE)
  expect_error(
    .random_effects_vignette_rename_file(
      prepared,
      destination,
      "the simulated prepared cache"
    ),
    "Could not the simulated prepared cache"
  )
  expect_true(file.exists(prepared))
  expect_false(file.exists(destination))
})

test_that("text fingerprints normalize line endings and use relative labels", {
  source_root <- withr::local_tempdir()
  lf_file <- file.path(source_root, "lf.R")
  crlf_file <- file.path(source_root, "crlf.R")
  writeBin(charToRaw("value <- 1\nnext_value <- 2\n"), lf_file)
  writeBin(charToRaw("value <- 1\r\nnext_value <- 2\r\n"), crlf_file)

  expect_identical(
    .random_effects_vignette_text_md5(lf_file),
    .random_effects_vignette_text_md5(crlf_file)
  )

  project_root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  generator_sources <-
    .random_effects_vignette_generator_sources(project_root)
  package_sources <-
    .random_effects_vignette_bayestools_sources(project_root)
  source_labels <- c(names(generator_sources), names(package_sources))

  expect_identical(names(generator_sources), c(
    "vignettes/RandomEffects.Rmd",
    "vignettes/random-effects-vignette-cache.R"
  ))
  expect_true(all(grepl(
    "^[^/]+(?:/[^/]+)*$",
    source_labels,
    perl = TRUE
  )))
  expect_false(any(grepl("\\\\", source_labels)))
  expect_false(any(grepl("^[[:alpha:]]:/", source_labels)))
  expect_false(any(startsWith(source_labels, "/")))
  expect_false("src/Makevars" %in% names(package_sources))
  expect_true(all(c(
    "src/Makevars.in",
    "src/Makevars.ucrt",
    "src/Makevars.win",
    "src/Makevars.win.common"
  ) %in% names(package_sources)))

  buildignore <- readLines(
    file.path(project_root, ".Rbuildignore"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  expect_true("^models($|/)" %in% buildignore)
  expect_true(file.exists(file.path(
    project_root,
    "vignettes",
    "RandomEffects.RDS"
  )))
  expect_false(any(grepl(
    "RandomEffects\\.RDS",
    buildignore
  )))
  expect_false("^vignettes($|/)" %in% buildignore)
  expect_true("^src/Makevars$" %in% buildignore)
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
    .random_effects_vignette_description_md5(source_description),
    .random_effects_vignette_description_md5(staged_description)
  )

  staged <- readLines(staged_description, warn = FALSE)
  staged[staged == "Version: 0.3.1.7"] <- "Version: 0.3.1.8"
  writeLines(staged, staged_description, useBytes = TRUE)
  expect_false(identical(
    .random_effects_vignette_description_md5(source_description),
    .random_effects_vignette_description_md5(staged_description)
  ))
})

test_that("vignette statically uses guarded regeneration and exact seeds", {
  rmd_file <- testthat::test_path("..", "..", "vignettes", "RandomEffects.Rmd")
  lines <- readLines(rmd_file, warn = FALSE, encoding = "UTF-8")
  setup_chunk <- .random_effects_test_chunk(lines, "setup")
  begin_chunk <- .random_effects_test_chunk(
    lines,
    "begin-random-effects-cache-regeneration"
  )
  save_chunk <- .random_effects_test_chunk(
    lines,
    "save-precomputed-random-effects"
  )

  expect_match(setup_chunk, "Missing required vignette cache")
  expect_match(
    setup_chunk,
    "format_random_effects_vignette_cache_error\\("
  )
  expect_match(
    begin_chunk,
    "begin_random_effects_vignette_cache_regeneration\\("
  )
  expect_match(begin_chunk, "pkgload::load_all\\(\"\\.\\.\"")
  expect_match(begin_chunk, "envir = knitr::knit_global\\(\\)")
  expect_false(grepl("implementation_state[[:space:]]*=", begin_chunk))
  expect_match(save_chunk, "write_random_effects_vignette_cache\\(")
  expect_match(save_chunk, "inherits = FALSE")
  expect_match(
    save_chunk,
    "generation = random_effects_cache_generation"
  )
  expect_false(any(grepl("saveRDS\\(", lines)))

  chunk_labels <- c(
    stan_correlated = "fit-stan-correlated-regenerate",
    stan_correlated_scaled = "fit-stan-correlated-regenerate",
    fit_sleep = "fit-sleep-regenerate",
    fit_sleep_independent = "fit-sleep-independent-regenerate",
    fit_cake_recipe = "fit-cake-recipe-regenerate",
    fit_cake_hcs = "fit-cake-hcs-regenerate",
    fit_cake_ar1 = "fit-cake-ar1-regenerate",
    fit_cake_har = "fit-cake-har-regenerate",
    fit_sleep_car = "fit-sleep-car-regenerate",
    fit_nested = "fit-nested-regenerate",
    stan_crossed = "fit-stan-crossed-regenerate",
    fit_crossed_independent = "fit-crossed-independent-regenerate",
    fit_crossed_allocation = "fit-crossed-allocation-regenerate"
  )
  seeds <- random_effects_vignette_cache_seeds()
  engines <- random_effects_vignette_cache_engines()
  regeneration_starts <- integer()
  for(model_name in names(chunk_labels)){
    chunk <- .random_effects_test_chunk(lines, chunk_labels[[model_name]])
    regeneration_starts <- c(
      regeneration_starts,
      grep(
        paste0(
          "^```\\{r[[:space:]]+",
          chunk_labels[[model_name]],
          "(?:,|\\})"
        ),
        lines,
        perl = TRUE
      )
    )
    expect_match(chunk, "eval = FALSE", info = model_name)
    expect_match(
      chunk,
      paste0(model_name, "[[:space:]]*<-"),
      info = model_name
    )
    expect_match(
      chunk,
      paste0("seed[[:space:]]*=[[:space:]]*", seeds[[model_name]], "\\b"),
      info = model_name
    )
    engine_call <- if(engines[[model_name]] == "rstanarm"){
      "rstanarm::stan_lmer\\("
    }else{
      "JAGS_fit\\("
    }
    expect_match(chunk, engine_call, info = model_name)
  }
  begin_start <- grep(
    "^```\\{r[[:space:]]+begin-random-effects-cache-regeneration(?:,|\\})",
    lines,
    perl = TRUE
  )
  writer_start <- grep(
    "^```\\{r[[:space:]]+save-precomputed-random-effects(?:,|\\})",
    lines,
    perl = TRUE
  )
  expect_lt(begin_start, min(regeneration_starts))
  expect_gt(writer_start, max(regeneration_starts))
})
