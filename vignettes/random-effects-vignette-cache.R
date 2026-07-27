random_effects_vignette_cache_schema <- function(){
  c(
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
}

random_effects_vignette_cache_engines <- function(){
  c(
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
}

random_effects_vignette_cache_seeds <- function(){
  c(
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
}

random_effects_vignette_cache_model_schema <- function(){
  list(
    classes = random_effects_vignette_cache_schema(),
    engines = random_effects_vignette_cache_engines(),
    seeds = random_effects_vignette_cache_seeds()
  )
}

random_effects_vignette_cache_names <- function(){
  names(random_effects_vignette_cache_schema())
}

.random_effects_vignette_cache_format <- function(){
  "BayesTools.RandomEffects.vignette-cache"
}

.random_effects_vignette_manifest_version <- function(){
  1L
}

.random_effects_vignette_cache_schema_version <- function(){
  1L
}

.random_effects_vignette_md5_raw <- function(value){
  hash_file <- tempfile("RandomEffects-hash-")
  on.exit(unlink(hash_file), add = TRUE)
  writeBin(value, hash_file)
  unname(tools::md5sum(hash_file))
}

.random_effects_vignette_object_md5 <- function(value){
  .random_effects_vignette_md5_raw(serialize(value, NULL, version = 3, xdr = TRUE))
}

.random_effects_vignette_sha256_raw <- function(value){
  tools_namespace <- asNamespace("tools")
  if(exists("sha256sum", envir = tools_namespace, inherits = FALSE)){
    sha256sum <- get("sha256sum", envir = tools_namespace, inherits = FALSE)
    return(unname(sha256sum(bytes = value)))
  }
  if(!requireNamespace("digest", quietly = TRUE)){
    stop(
      "SHA-256 cache integrity requires R >= 4.5.0 or the 'digest' package.",
      call. = FALSE
    )
  }
  unname(digest::digest(value, algo = "sha256", serialize = FALSE))
}

.random_effects_vignette_sha256_file <- function(path){
  tools_namespace <- asNamespace("tools")
  if(exists("sha256sum", envir = tools_namespace, inherits = FALSE)){
    sha256sum <- get("sha256sum", envir = tools_namespace, inherits = FALSE)
    return(unname(sha256sum(files = path)))
  }
  if(!requireNamespace("digest", quietly = TRUE)){
    stop(
      "SHA-256 cache integrity requires R >= 4.5.0 or the 'digest' package.",
      call. = FALSE
    )
  }
  unname(digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ))
}

.random_effects_vignette_object_sha256 <- function(value){
  .random_effects_vignette_sha256_raw(
    serialize(value, NULL, version = 3, xdr = TRUE)
  )
}

.random_effects_vignette_model_hashes <- function(models){
  model_names <- random_effects_vignette_cache_names()
  hashes <- vapply(
    model_names,
    function(name){
      .random_effects_vignette_object_sha256(models[[name]])
    },
    character(1)
  )
  stats::setNames(unname(hashes), model_names)
}

.random_effects_vignette_text_md5 <- function(path){
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  text <- paste(lines, collapse = "\n")
  if(length(lines) > 0L){
    text <- paste0(text, "\n")
  }
  .random_effects_vignette_md5_raw(charToRaw(enc2utf8(text)))
}

.random_effects_vignette_project_root <- function(cache_file, project_root = NULL){
  if(is.null(project_root)){
    project_root <- file.path(dirname(cache_file), "..")
  }
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if(!file.exists(file.path(project_root, "DESCRIPTION"))){
    stop(
      "Could not resolve the BayesTools project root for the RandomEffects cache.",
      call. = FALSE
    )
  }
  project_root
}

.random_effects_vignette_relative_paths <- function(paths, project_root){
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  prefix <- paste0(project_root, "/")
  if(any(!startsWith(paths, prefix))){
    stop("RandomEffects dependency paths must be inside the project root.", call. = FALSE)
  }
  substring(paths, nchar(prefix) + 1L)
}

.random_effects_vignette_source_hashes <- function(paths, project_root){
  paths <- unique(paths)
  labels <- .random_effects_vignette_relative_paths(paths, project_root)
  source_order <- order(enc2utf8(labels), method = "radix")
  paths <- paths[source_order]
  labels <- labels[source_order]
  hashes <- vapply(paths, .random_effects_vignette_text_md5, character(1))
  stats::setNames(unname(hashes), labels)
}

.random_effects_vignette_generator_sources <- function(project_root){
  paths <- file.path(
    project_root,
    "vignettes",
    c("RandomEffects.Rmd", "random-effects-vignette-cache.R")
  )
  missing <- paths[!file.exists(paths)]
  if(length(missing) > 0L){
    stop(
      "Missing RandomEffects vignette generator source(s): ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }
  .random_effects_vignette_source_hashes(paths, project_root)
}

.random_effects_vignette_bayestools_sources <- function(project_root){
  r_dir <- file.path(project_root, "R")
  src_dir <- file.path(project_root, "src")
  r_files <- if(dir.exists(r_dir)){
    list.files(
      r_dir,
      pattern = "\\.[Rr]$",
      recursive = TRUE,
      full.names = TRUE
    )
  }else{
    character()
  }
  src_files <- if(dir.exists(src_dir)){
    list.files(src_dir, recursive = TRUE, full.names = TRUE)
  }else{
    character()
  }
  source_extensions <- c(
    "c", "cc", "cpp", "cxx", "f", "f77", "f90", "f95",
    "h", "hh", "hpp", "hxx", "inc"
  )
  src_files <- src_files[
    tolower(tools::file_ext(src_files)) %in% source_extensions |
      grepl("^Makevars\\.", basename(src_files))
  ]
  root_files <- file.path(
    project_root,
    c(
      "DESCRIPTION", "NAMESPACE", "configure", "configure.win",
      "cleanup", "cleanup.win"
    )
  )
  root_files <- root_files[file.exists(root_files)]
  paths <- sort(unique(c(r_files, src_files, root_files)))
  if(length(paths) == 0L){
    stop("No BayesTools source dependencies were found.", call. = FALSE)
  }
  .random_effects_vignette_source_hashes(paths, project_root)
}

.random_effects_vignette_dataset_hashes <- function(){
  dataset_names <- c("sleepstudy", "cake", "Pastes", "Penicillin")
  data_environment <- new.env(parent = emptyenv())
  for(dataset_name in dataset_names){
    suppressWarnings(utils::data(
      list = dataset_name,
      package = "lme4",
      envir = data_environment
    ))
  }
  missing <- dataset_names[
    !vapply(
      dataset_names,
      exists,
      logical(1),
      envir = data_environment,
      inherits = FALSE
    )
  ]
  if(length(missing) > 0L){
    stop(
      "Could not load RandomEffects generator dataset(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  hashes <- vapply(dataset_names, function(dataset_name){
    .random_effects_vignette_object_md5(get(
      dataset_name,
      envir = data_environment,
      inherits = FALSE
    ))
  }, character(1))
  names(hashes) <- paste0("lme4::", dataset_names)
  hashes
}

.random_effects_vignette_package_versions <- function(){
  packages <- c(
    "BayesTools", "lme4", "Matrix", "reformulas", "rstanarm",
    "rstan", "StanHeaders", "runjags", "rjags", "coda"
  )
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if(length(missing) > 0L){
    stop(
      "Missing RandomEffects generator package(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  versions <- vapply(packages, function(package){
    as.character(utils::packageVersion(package))
  }, character(1))
  versions
}

random_effects_vignette_dependency_state <- function(
    cache_file = file.path("models", "RandomEffects.RDS"),
    project_root = NULL){
  project_root <- .random_effects_vignette_project_root(cache_file, project_root)
  list(
    generator_sources = .random_effects_vignette_generator_sources(project_root),
    bayestools_sources = .random_effects_vignette_bayestools_sources(project_root),
    datasets = .random_effects_vignette_dataset_hashes(),
    package_versions = .random_effects_vignette_package_versions()
  )
}

.random_effects_vignette_named_hashes_valid <- function(value){
  is.character(value) &&
    length(value) > 0L &&
    !is.null(names(value)) &&
    !anyNA(names(value)) &&
    all(nzchar(names(value))) &&
    !anyDuplicated(names(value)) &&
    !anyNA(value) &&
    all(grepl("^[[:xdigit:]]{32}$", value))
}

.random_effects_vignette_named_versions_valid <- function(value){
  is.character(value) &&
    length(value) > 0L &&
    !is.null(names(value)) &&
    !anyNA(names(value)) &&
    all(nzchar(names(value))) &&
    !anyDuplicated(names(value)) &&
    !anyNA(value) &&
    all(nzchar(value))
}

.random_effects_vignette_dependency_state_error <- function(value){
  if(!is.list(value) || !identical(
    names(value),
    c(
      "generator_sources", "bayestools_sources", "datasets",
      "package_versions"
    )
  )){
    return("dependency state has an unexpected schema")
  }
  if(!.random_effects_vignette_named_hashes_valid(value$generator_sources)){
    return("generator source fingerprints are invalid")
  }
  if(!.random_effects_vignette_named_hashes_valid(value$bayestools_sources)){
    return("BayesTools source fingerprints are invalid")
  }
  if(!.random_effects_vignette_named_hashes_valid(value$datasets)){
    return("dataset fingerprints are invalid")
  }
  if(!.random_effects_vignette_named_versions_valid(value$package_versions)){
    return("package versions are invalid")
  }
  NULL
}

.random_effects_vignette_implementation_state_error <- function(
    value, dependencies){
  if(!is.list(value) || !identical(
    names(value),
    c("fit_backend_fingerprint", "bayestools_source_fingerprint")
  )){
    return("implementation state has an unexpected schema")
  }
  if(!is.character(value$fit_backend_fingerprint) ||
      length(value$fit_backend_fingerprint) != 1L ||
      is.na(value$fit_backend_fingerprint) ||
      !grepl("^[[:xdigit:]]{32}$", value$fit_backend_fingerprint)){
    return("loaded backend fingerprint is invalid")
  }
  if(!is.character(value$bayestools_source_fingerprint) ||
      length(value$bayestools_source_fingerprint) != 1L ||
      is.na(value$bayestools_source_fingerprint) ||
      !grepl(
        "^[[:xdigit:]]{64}$",
        value$bayestools_source_fingerprint
      )){
    return("BayesTools source contract fingerprint is invalid")
  }
  expected_source_fingerprint <-
    .random_effects_vignette_object_sha256(
      dependencies$bayestools_sources
    )
  if(!identical(
    value$bayestools_source_fingerprint,
    expected_source_fingerprint
  )){
    return("BayesTools source contract does not match current dependencies")
  }
  NULL
}

.random_effects_vignette_current_implementation <- function(
    dependencies,
    cache_file = file.path("models", "RandomEffects.RDS"),
    project_root = NULL){
  project_root <- .random_effects_vignette_project_root(
    cache_file,
    project_root
  )
  if(!requireNamespace("BayesTools", quietly = TRUE)){
    stop(
      "BayesTools must be loaded from the current project before regeneration.",
      call. = FALSE
    )
  }
  namespace_root <- tryCatch(
    getNamespaceInfo(asNamespace("BayesTools"), "path"),
    error = function(e) ""
  )
  namespace_root <- tryCatch(
    normalizePath(namespace_root, winslash = "/", mustWork = TRUE),
    error = function(e) ""
  )
  if(!identical(namespace_root, project_root)){
    stop(
      paste0(
        "The loaded BayesTools namespace is not the current source tree. ",
        "Run pkgload::load_all() from the BayesTools project and restart ",
        "RandomEffects cache regeneration."
      ),
      call. = FALSE
    )
  }

  loaded_fingerprint <- BayesTools::fit_backend_fingerprint()
  compute_fingerprint <- getFromNamespace(
    ".compute_fit_backend_fingerprint",
    "BayesTools"
  )
  loaded_dll <- getFromNamespace(".fit_backend_dll", "BayesTools")()
  current_fingerprint <- compute_fingerprint(
    package_root = project_root,
    loaded_dll = loaded_dll
  )
  valid_fingerprint <- function(value){
    is.character(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      grepl("^[[:xdigit:]]{32}$", value)
  }
  if(!valid_fingerprint(loaded_fingerprint) ||
      !valid_fingerprint(current_fingerprint)){
    stop(
      "Could not verify the loaded BayesTools implementation.",
      call. = FALSE
    )
  }
  if(!identical(loaded_fingerprint, current_fingerprint)){
    stop(
      paste0(
        "BayesTools source or native code changed after the namespace was ",
        "loaded. Reload the current project and restart RandomEffects cache ",
        "regeneration."
      ),
      call. = FALSE
    )
  }

  list(
    fit_backend_fingerprint = loaded_fingerprint,
    bayestools_source_fingerprint =
      .random_effects_vignette_object_sha256(
        dependencies$bayestools_sources
      )
  )
}

.random_effects_vignette_contract <- function(dependencies){
  list(
    format = .random_effects_vignette_cache_format(),
    manifest_version = .random_effects_vignette_manifest_version(),
    cache_schema_version = .random_effects_vignette_cache_schema_version(),
    model_schema = random_effects_vignette_cache_model_schema(),
    dependencies = dependencies
  )
}

.random_effects_vignette_compatibility_fingerprint <- function(dependencies){
  .random_effects_vignette_object_md5(
    .random_effects_vignette_contract(dependencies)
  )
}

.random_effects_vignette_generation_fingerprint <- function(
    compatibility_fingerprint, model_hashes, producer){
  .random_effects_vignette_object_sha256(list(
    compatibility_fingerprint = compatibility_fingerprint,
    model_hashes = model_hashes,
    producer = producer
  ))
}

.random_effects_vignette_producer <- function(){
  if(!requireNamespace("BayesTools", quietly = TRUE)){
    stop("BayesTools must be installed before regenerating its vignette cache.", call. = FALSE)
  }
  if(!requireNamespace("rjags", quietly = TRUE)){
    stop(
      "The 'rjags' package is required to record the JAGS runtime.",
      call. = FALSE
    )
  }
  backend_fingerprint <- BayesTools::fit_backend_fingerprint()
  if(!is.character(backend_fingerprint) ||
      length(backend_fingerprint) != 1L ||
      is.na(backend_fingerprint) ||
      !grepl("^[[:xdigit:]]{32}$", backend_fingerprint)){
    stop("Could not fingerprint the loaded BayesTools fitting backend.", call. = FALSE)
  }
  namespaces <- sort(loadedNamespaces())
  namespace_versions <- vapply(namespaces, function(namespace){
    tryCatch(
      as.character(utils::packageVersion(namespace)),
      error = function(e) NA_character_
    )
  }, character(1))
  namespace_versions <- namespace_versions[!is.na(namespace_versions)]
  list(
    generated_at_utc = format(
      Sys.time(),
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    platform = unname(R.version[["platform"]]),
    fit_backend_fingerprint = backend_fingerprint,
    loaded_package_versions = namespace_versions,
    runtime = list(
      r_version = as.character(getRversion()),
      jags_version = as.character(rjags::jags.version()),
      rng_kind = RNGkind(),
      contrasts = unname(getOption("contrasts"))
    )
  )
}

.random_effects_vignette_producer_error <- function(value){
  if(!is.list(value) || !identical(
    names(value),
    c(
      "generated_at_utc", "platform", "fit_backend_fingerprint",
      "loaded_package_versions", "runtime"
    )
  )){
    return("producer metadata has an unexpected schema")
  }
  scalar_character <- function(x){
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if(!scalar_character(value$generated_at_utc) ||
      !grepl(
        "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
        value$generated_at_utc
      )){
    return("producer timestamp is invalid")
  }
  if(!scalar_character(value$platform)){
    return("producer platform is invalid")
  }
  if(!scalar_character(value$fit_backend_fingerprint) ||
      !grepl("^[[:xdigit:]]{32}$", value$fit_backend_fingerprint)){
    return("producer backend fingerprint is invalid")
  }
  if(!.random_effects_vignette_named_versions_valid(
    value$loaded_package_versions
  )){
    return("producer package versions are invalid")
  }
  if(!is.list(value$runtime) || !identical(
    names(value$runtime),
    c("r_version", "jags_version", "rng_kind", "contrasts")
  )){
    return("producer runtime has an unexpected schema")
  }
  if(!scalar_character(value$runtime$r_version)){
    return("producer R runtime version is invalid")
  }
  if(!scalar_character(value$runtime$jags_version)){
    return("producer JAGS runtime version is invalid")
  }
  if(!is.character(value$runtime$rng_kind) ||
      length(value$runtime$rng_kind) == 0L ||
      anyNA(value$runtime$rng_kind)){
    return("producer RNG runtime state is invalid")
  }
  if(!is.character(value$runtime$contrasts) ||
      length(value$runtime$contrasts) == 0L ||
      anyNA(value$runtime$contrasts)){
    return("producer contrast runtime state is invalid")
  }
  NULL
}

.random_effects_vignette_manifest <- function(
    dependencies, model_hashes, producer){
  compatibility_fingerprint <-
    .random_effects_vignette_compatibility_fingerprint(dependencies)
  list(
    format = .random_effects_vignette_cache_format(),
    manifest_version = .random_effects_vignette_manifest_version(),
    cache_schema_version = .random_effects_vignette_cache_schema_version(),
    model_schema = random_effects_vignette_cache_model_schema(),
    model_hashes = model_hashes,
    dependencies = dependencies,
    compatibility_fingerprint = compatibility_fingerprint,
    producer = producer,
    generation_fingerprint =
      .random_effects_vignette_generation_fingerprint(
        compatibility_fingerprint,
        model_hashes,
        producer
      )
  )
}

.random_effects_vignette_model_hashes_error <- function(value){
  if(!is.character(value) ||
      !identical(names(value), random_effects_vignette_cache_names()) ||
      anyNA(value) ||
      !all(grepl("^[[:xdigit:]]{64}$", value))){
    return("model payload hashes are invalid")
  }
  NULL
}

.random_effects_vignette_manifest_error <- function(manifest){
  if(!is.list(manifest) || !identical(
    names(manifest),
    c(
      "format", "manifest_version", "cache_schema_version", "model_schema",
      "model_hashes", "dependencies", "compatibility_fingerprint",
      "producer", "generation_fingerprint"
    )
  )){
    return("manifest has an unexpected schema")
  }
  if(!identical(manifest$format, .random_effects_vignette_cache_format())){
    return("manifest format is unsupported")
  }
  if(!identical(
    manifest$manifest_version,
    .random_effects_vignette_manifest_version()
  )){
    return("manifest version is unsupported")
  }
  if(!identical(
    manifest$cache_schema_version,
    .random_effects_vignette_cache_schema_version()
  )){
    return("cache schema version is unsupported")
  }
  if(!identical(
    manifest$model_schema,
    random_effects_vignette_cache_model_schema()
  )){
    return("model schema does not match the current vignette contract")
  }
  model_hashes_error <-
    .random_effects_vignette_model_hashes_error(manifest$model_hashes)
  if(!is.null(model_hashes_error)){
    return(model_hashes_error)
  }
  dependency_error <-
    .random_effects_vignette_dependency_state_error(manifest$dependencies)
  if(!is.null(dependency_error)){
    return(dependency_error)
  }
  expected_compatibility <-
    .random_effects_vignette_compatibility_fingerprint(
      manifest$dependencies
    )
  if(!identical(
    manifest$compatibility_fingerprint,
    expected_compatibility
  )){
    return("compatibility fingerprint does not match the manifest inputs")
  }
  producer_error <- .random_effects_vignette_producer_error(manifest$producer)
  if(!is.null(producer_error)){
    return(producer_error)
  }
  expected_generation <- .random_effects_vignette_generation_fingerprint(
    manifest$compatibility_fingerprint,
    manifest$model_hashes,
    manifest$producer
  )
  if(!identical(manifest$generation_fingerprint, expected_generation)){
    return(
      "generation fingerprint does not match the payload and producer metadata"
    )
  }
  NULL
}

.random_effects_vignette_model_status <- function(models){
  schema <- random_effects_vignette_cache_schema()
  status <- list(
    is_list = is.list(models),
    is_named = FALSE,
    missing_objects = names(schema),
    extra_objects = character(),
    order_valid = FALSE,
    invalid_objects = names(schema)
  )
  if(!status$is_list){
    return(status)
  }
  model_names <- names(models)
  status$is_named <- !is.null(model_names) &&
    length(model_names) == length(models) &&
    !anyNA(model_names) &&
    all(nzchar(model_names)) &&
    !anyDuplicated(model_names)
  if(!status$is_named){
    return(status)
  }
  status$missing_objects <- setdiff(names(schema), model_names)
  status$extra_objects <- setdiff(model_names, names(schema))
  status$order_valid <- identical(model_names, names(schema))
  present <- intersect(names(schema), model_names)
  status$invalid_objects <- present[
    !vapply(present, function(name){
      inherits(models[[name]], schema[[name]])
    }, logical(1))
  ]
  status
}

validate_random_effects_vignette_cache <- function(
    cache_file = file.path("models", "RandomEffects.RDS"),
    load = FALSE,
    dependency_state = NULL,
    project_root = NULL){
  schema <- random_effects_vignette_cache_schema()
  status <- list(
    file_exists = file.exists(cache_file),
    read_error = NULL,
    legacy_cache = FALSE,
    envelope_valid = FALSE,
    manifest_valid = FALSE,
    manifest_error = NULL,
    dependency_error = NULL,
    stale_dependencies = character(),
    is_list = FALSE,
    is_named = FALSE,
    missing_objects = names(schema),
    extra_objects = character(),
    order_valid = FALSE,
    invalid_objects = names(schema),
    payload_hash_error = NULL,
    invalid_payloads = names(schema),
    payload_valid = FALSE,
    valid = FALSE,
    cache = NULL
  )
  if(!status$file_exists){
    return(status)
  }

  envelope <- tryCatch(
    suppressWarnings(readRDS(cache_file)),
    error = function(e) e
  )
  if(inherits(envelope, "error")){
    status$read_error <- conditionMessage(envelope)
    return(status)
  }
  if(is.list(envelope) &&
      is.null(envelope$manifest) &&
      length(intersect(names(schema), names(envelope))) > 0L){
    status$legacy_cache <- TRUE
    return(status)
  }
  status$envelope_valid <- is.list(envelope) &&
    identical(names(envelope), c("manifest", "models"))
  if(!status$envelope_valid){
    return(status)
  }

  status$manifest_error <-
    .random_effects_vignette_manifest_error(envelope$manifest)
  status$manifest_valid <- is.null(status$manifest_error)
  if(!status$manifest_valid){
    return(status)
  }

  if(is.null(dependency_state)){
    dependency_state <- tryCatch(
      random_effects_vignette_dependency_state(
        cache_file = cache_file,
        project_root = project_root
      ),
      error = function(e) e
    )
    if(inherits(dependency_state, "error")){
      status$dependency_error <- conditionMessage(dependency_state)
      return(status)
    }
  }
  dependency_error <-
    .random_effects_vignette_dependency_state_error(dependency_state)
  if(!is.null(dependency_error)){
    status$dependency_error <- dependency_error
    return(status)
  }
  status$stale_dependencies <- names(dependency_state)[
    !vapply(names(dependency_state), function(name){
      identical(
        dependency_state[[name]],
        envelope$manifest$dependencies[[name]]
      )
    }, logical(1))
  ]

  model_status <- .random_effects_vignette_model_status(envelope$models)
  for(name in names(model_status)){
    status[[name]] <- model_status[[name]]
  }
  model_structure_valid <- status$is_list &&
    status$is_named &&
    length(status$missing_objects) == 0L &&
    length(status$extra_objects) == 0L &&
    status$order_valid &&
    length(status$invalid_objects) == 0L
  if(model_structure_valid){
    current_model_hashes <- tryCatch(
      .random_effects_vignette_model_hashes(envelope$models),
      error = function(e) e
    )
    if(inherits(current_model_hashes, "error")){
      status$payload_hash_error <- conditionMessage(current_model_hashes)
    }else{
      hash_matches <- current_model_hashes == envelope$manifest$model_hashes
      hash_matches[is.na(hash_matches)] <- FALSE
      status$invalid_payloads <- names(current_model_hashes)[!hash_matches]
      status$payload_valid <- length(status$invalid_payloads) == 0L
    }
  }
  status$valid <- status$manifest_valid &&
    length(status$stale_dependencies) == 0L &&
    model_structure_valid &&
    status$payload_valid
  if(isTRUE(status$valid) && isTRUE(load)){
    status$cache <- envelope$models
  }
  status
}

format_random_effects_vignette_cache_error <- function(status){
  if(!isTRUE(status$file_exists)){
    return("Missing precomputed RandomEffects vignette cache.")
  }
  if(!is.null(status$read_error)){
    return(paste0(
      "Could not read precomputed RandomEffects vignette cache: ",
      status$read_error
    ))
  }
  if(isTRUE(status$legacy_cache)){
    return(
      paste0(
        "Precomputed RandomEffects vignette cache uses the legacy format ",
        "without a manifest; regenerate it."
      )
    )
  }
  if(!isTRUE(status$envelope_valid)){
    return(
      "Precomputed RandomEffects vignette cache must contain manifest and models."
    )
  }
  if(!isTRUE(status$manifest_valid)){
    return(paste0(
      "Precomputed RandomEffects vignette cache manifest is invalid: ",
      status$manifest_error,
      "."
    ))
  }
  if(!is.null(status$dependency_error)){
    return(paste0(
      "Could not fingerprint current RandomEffects dependencies: ",
      status$dependency_error
    ))
  }
  if(length(status$stale_dependencies) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache is stale; changed dependencies: ",
      paste(status$stale_dependencies, collapse = ", "),
      "."
    ))
  }
  if(!isTRUE(status$is_list) || !isTRUE(status$is_named)){
    return(
      "Precomputed RandomEffects vignette cache models must be a named list."
    )
  }
  if(length(status$missing_objects) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache is missing objects: ",
      paste(status$missing_objects, collapse = ", ")
    ))
  }
  if(length(status$extra_objects) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache has unexpected objects: ",
      paste(status$extra_objects, collapse = ", ")
    ))
  }
  if(!isTRUE(status$order_valid)){
    return(
      "Precomputed RandomEffects vignette cache objects are in the wrong order."
    )
  }
  if(length(status$invalid_objects) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache has unexpected object classes: ",
      paste(status$invalid_objects, collapse = ", ")
    ))
  }
  if(!is.null(status$payload_hash_error)){
    return(paste0(
      "Could not verify RandomEffects model payloads: ",
      status$payload_hash_error
    ))
  }
  if(length(status$invalid_payloads) > 0L){
    return(paste0(
      "Precomputed RandomEffects vignette cache has modified payloads: ",
      paste(status$invalid_payloads, collapse = ", ")
    ))
  }
  "Precomputed RandomEffects vignette cache is invalid."
}

stop_if_invalid_random_effects_vignette_cache <- function(
    cache_file = file.path("models", "RandomEffects.RDS"),
    dependency_state = NULL,
    project_root = NULL){
  status <- validate_random_effects_vignette_cache(
    cache_file = cache_file,
    dependency_state = dependency_state,
    project_root = project_root
  )
  if(!isTRUE(status$valid)){
    stop(format_random_effects_vignette_cache_error(status), call. = FALSE)
  }
  invisible(status)
}

begin_random_effects_vignette_cache_regeneration <- function(
    envir = parent.frame(),
    cache_file = file.path("models", "RandomEffects.RDS"),
    dependency_state = NULL,
    implementation_state = NULL,
    project_root = NULL){
  if(!is.environment(envir)){
    stop("'envir' must be an environment.", call. = FALSE)
  }
  cache_names <- random_effects_vignette_cache_names()
  existing <- intersect(cache_names, ls(envir = envir, all.names = TRUE))
  if(length(existing) > 0L){
    rm(list = existing, envir = envir)
  }
  if(is.null(dependency_state)){
    dependency_state <- random_effects_vignette_dependency_state(
      cache_file = cache_file,
      project_root = project_root
    )
  }
  dependency_error <-
    .random_effects_vignette_dependency_state_error(dependency_state)
  if(!is.null(dependency_error)){
    stop(
      "Cannot begin RandomEffects cache regeneration: ",
      dependency_error,
      ".",
      call. = FALSE
    )
  }
  if(is.null(implementation_state)){
    implementation_state <- .random_effects_vignette_current_implementation(
      dependencies = dependency_state,
      cache_file = cache_file,
      project_root = project_root
    )
  }
  implementation_error <-
    .random_effects_vignette_implementation_state_error(
      implementation_state,
      dependency_state
    )
  if(!is.null(implementation_error)){
    stop(
      "Cannot begin RandomEffects cache regeneration: ",
      implementation_error,
      ".",
      call. = FALSE
    )
  }
  checkpoint_state <- new.env(parent = emptyenv())
  checkpoint_state$used <- FALSE
  structure(
    list(
      format = "BayesTools.RandomEffects.vignette-cache-regeneration",
      dependencies = dependency_state,
      implementation = implementation_state,
      checkpoint_state = checkpoint_state
    ),
    class = "BayesTools_random_effects_vignette_cache_regeneration"
  )
}

.random_effects_vignette_dependency_changes <- function(before, after){
  names(before)[
    !vapply(names(before), function(name){
      identical(before[[name]], after[[name]])
    }, logical(1))
  ]
}

.random_effects_vignette_validate_models_or_stop <- function(models){
  model_status <- .random_effects_vignette_model_status(models)
  if(!model_status$is_list || !model_status$is_named){
    stop("RandomEffects cache models must be a named list.", call. = FALSE)
  }
  if(length(model_status$missing_objects) > 0L){
    stop(
      "Cannot write RandomEffects cache; missing objects: ",
      paste(model_status$missing_objects, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if(length(model_status$extra_objects) > 0L){
    stop(
      "Cannot write RandomEffects cache; unexpected objects: ",
      paste(model_status$extra_objects, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if(!model_status$order_valid){
    stop("Cannot write RandomEffects cache; objects are in the wrong order.", call. = FALSE)
  }
  if(length(model_status$invalid_objects) > 0L){
    stop(
      "Cannot write RandomEffects cache; unexpected object classes: ",
      paste(model_status$invalid_objects, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.random_effects_vignette_transaction_paths <- function(cache_file){
  cache_directory <- dirname(cache_file)
  cache_name <- basename(cache_file)
  list(
    new = file.path(cache_directory, paste0(".", cache_name, ".new")),
    backup = file.path(cache_directory, paste0(".", cache_name, ".backup")),
    journal = file.path(cache_directory, paste0(".", cache_name, ".journal")),
    journal_new =
      file.path(cache_directory, paste0(".", cache_name, ".journal.new"))
  )
}

.random_effects_vignette_remove_file <- function(path, description){
  if(!file.exists(path)){
    return(invisible(TRUE))
  }
  result <- unlink(path, recursive = FALSE, force = TRUE)
  if(!identical(result, 0L) || file.exists(path)){
    stop("Could not remove ", description, ": ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.random_effects_vignette_rename_file <- function(from, to, description){
  if(!file.exists(from) || !isTRUE(file.rename(from, to)) ||
      file.exists(from) || !file.exists(to)){
    stop("Could not ", description, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.random_effects_vignette_copy_file <- function(from, to, description){
  if(!file.exists(from) ||
      !isTRUE(file.copy(from, to, overwrite = FALSE, copy.mode = TRUE)) ||
      !file.exists(to)){
    stop("Could not ", description, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.random_effects_vignette_transaction_journal_error <- function(journal){
  if(!is.list(journal) || !identical(
    names(journal),
    c(
      "format", "version", "had_cache", "previous_sha256", "new_sha256"
    )
  )){
    return("transaction journal has an unexpected schema")
  }
  if(!identical(
    journal$format,
    "BayesTools.RandomEffects.vignette-cache-transaction"
  ) || !identical(journal$version, 1L)){
    return("transaction journal version is unsupported")
  }
  if(!is.logical(journal$had_cache) ||
      length(journal$had_cache) != 1L ||
      is.na(journal$had_cache)){
    return("transaction journal cache state is invalid")
  }
  valid_sha256 <- function(value){
    is.character(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      grepl("^[[:xdigit:]]{64}$", value)
  }
  if(isTRUE(journal$had_cache)){
    if(!valid_sha256(journal$previous_sha256)){
      return("transaction journal previous-cache hash is invalid")
    }
  }else if(!identical(journal$previous_sha256, NA_character_)){
    return("transaction journal unexpectedly describes a previous cache")
  }
  if(!valid_sha256(journal$new_sha256)){
    return("transaction journal new-cache hash is invalid")
  }
  NULL
}

.random_effects_vignette_cleanup_transaction <- function(paths){
  .random_effects_vignette_remove_file(
    paths$new,
    "RandomEffects transaction new-file sidecar"
  )
  .random_effects_vignette_remove_file(
    paths$backup,
    "RandomEffects transaction backup"
  )
  .random_effects_vignette_remove_file(
    paths$journal_new,
    "RandomEffects transaction journal temporary file"
  )
  .random_effects_vignette_remove_file(
    paths$journal,
    "RandomEffects transaction journal"
  )
  invisible(TRUE)
}

.random_effects_vignette_replace_file <- function(
    from, to, platform, description){
  if(!platform %in% c("unix", "windows")){
    stop("Unsupported cache-installation platform.", call. = FALSE)
  }
  if(platform == "windows" && file.exists(to)){
    .random_effects_vignette_remove_file(
      to,
      paste0(description, " destination")
    )
  }
  .random_effects_vignette_rename_file(from, to, description)
}

.random_effects_vignette_restore_previous <- function(
    cache_file, paths, journal, platform){
  if(isTRUE(journal$had_cache)){
    canonical_hash <- if(file.exists(cache_file)){
      .random_effects_vignette_sha256_file(cache_file)
    }else{
      NA_character_
    }
    if(!identical(canonical_hash, journal$previous_sha256)){
      if(!file.exists(paths$backup) ||
          !identical(
            .random_effects_vignette_sha256_file(paths$backup),
            journal$previous_sha256
          )){
        stop(
          paste0(
            "Cannot recover the previous RandomEffects cache: the backup ",
            "does not match the transaction journal."
          ),
          call. = FALSE
        )
      }
      .random_effects_vignette_replace_file(
        paths$backup,
        cache_file,
        platform,
        "restore the previous RandomEffects cache"
      )
    }
    if(!identical(
      .random_effects_vignette_sha256_file(cache_file),
      journal$previous_sha256
    )){
      stop(
        "Recovered RandomEffects cache does not match its journal hash.",
        call. = FALSE
      )
    }
  }else if(file.exists(cache_file)){
    .random_effects_vignette_remove_file(
      cache_file,
      "failed RandomEffects cache installation"
    )
  }
  .random_effects_vignette_cleanup_transaction(paths)
  invisible(TRUE)
}

recover_random_effects_vignette_cache <- function(
    cache_file = file.path("models", "RandomEffects.RDS"),
    .platform = .Platform$OS.type){
  paths <- .random_effects_vignette_transaction_paths(cache_file)
  if(!file.exists(paths$journal)){
    .random_effects_vignette_cleanup_transaction(paths)
    return(invisible("no-transaction"))
  }

  journal <- tryCatch(
    suppressWarnings(readRDS(paths$journal)),
    error = function(e) e
  )
  if(inherits(journal, "error")){
    stop(
      "Could not read the RandomEffects cache transaction journal: ",
      conditionMessage(journal),
      call. = FALSE
    )
  }
  journal_error <-
    .random_effects_vignette_transaction_journal_error(journal)
  if(!is.null(journal_error)){
    stop(
      "Cannot recover the RandomEffects cache: ",
      journal_error,
      ".",
      call. = FALSE
    )
  }

  canonical_hash <- if(file.exists(cache_file)){
    .random_effects_vignette_sha256_file(cache_file)
  }else{
    NA_character_
  }
  if(identical(canonical_hash, journal$new_sha256) ||
      (isTRUE(journal$had_cache) &&
       identical(canonical_hash, journal$previous_sha256))){
    .random_effects_vignette_cleanup_transaction(paths)
    return(invisible("canonical-complete"))
  }

  if(isTRUE(journal$had_cache)){
    .random_effects_vignette_restore_previous(
      cache_file,
      paths,
      journal,
      .platform
    )
    return(invisible("previous-restored"))
  }
  if(file.exists(paths$new) &&
      identical(
        .random_effects_vignette_sha256_file(paths$new),
        journal$new_sha256
      )){
    .random_effects_vignette_replace_file(
      paths$new,
      cache_file,
      .platform,
      "complete the interrupted RandomEffects cache installation"
    )
    if(!identical(
      .random_effects_vignette_sha256_file(cache_file),
      journal$new_sha256
    )){
      stop(
        "Recovered RandomEffects cache does not match its journal hash.",
        call. = FALSE
      )
    }
    .random_effects_vignette_cleanup_transaction(paths)
    return(invisible("new-installed"))
  }

  stop(
    paste0(
      "Cannot recover the interrupted RandomEffects cache installation: ",
      "no journal-matching cache file remains."
    ),
    call. = FALSE
  )
}

.random_effects_vignette_write_journal <- function(paths, journal){
  journal_error <- .random_effects_vignette_transaction_journal_error(journal)
  if(!is.null(journal_error)){
    stop(
      "Refusing to write an invalid RandomEffects transaction journal: ",
      journal_error,
      ".",
      call. = FALSE
    )
  }
  saveRDS(journal, paths$journal_new, version = 3)
  if(!file.exists(paths$journal_new)){
    stop("Could not write the RandomEffects transaction journal.", call. = FALSE)
  }
  written_journal <- tryCatch(
    suppressWarnings(readRDS(paths$journal_new)),
    error = function(e) e
  )
  if(inherits(written_journal, "error") ||
      !identical(written_journal, journal)){
    stop(
      "Could not verify the RandomEffects transaction journal.",
      call. = FALSE
    )
  }
  .random_effects_vignette_rename_file(
    paths$journal_new,
    paths$journal,
    "install the RandomEffects transaction journal"
  )
}

.random_effects_vignette_atomic_write <- function(
    envelope,
    cache_file,
    dependency_state,
    .platform = .Platform$OS.type){
  cache_directory <- dirname(cache_file)
  dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(cache_directory)){
    stop("Could not create the RandomEffects cache directory.", call. = FALSE)
  }
  recover_random_effects_vignette_cache(
    cache_file = cache_file,
    .platform = .platform
  )
  paths <- .random_effects_vignette_transaction_paths(cache_file)

  preparation_error <- tryCatch({
    saveRDS(envelope, paths$new, version = 3)
    if(!file.exists(paths$new)){
      stop("Could not write the regenerated RandomEffects cache.", call. = FALSE)
    }
    temporary_status <- validate_random_effects_vignette_cache(
      cache_file = paths$new,
      dependency_state = dependency_state
    )
    if(!isTRUE(temporary_status$valid)){
      stop(
        "Refusing to install invalid RandomEffects cache: ",
        format_random_effects_vignette_cache_error(temporary_status),
        call. = FALSE
      )
    }

    had_cache <- file.exists(cache_file)
    previous_sha256 <- NA_character_
    if(had_cache){
      previous_sha256 <-
        .random_effects_vignette_sha256_file(cache_file)
      .random_effects_vignette_copy_file(
        cache_file,
        paths$backup,
        "create the RandomEffects recovery backup"
      )
      if(!identical(
        .random_effects_vignette_sha256_file(paths$backup),
        previous_sha256
      )){
        stop(
          "RandomEffects recovery backup does not match the canonical cache.",
          call. = FALSE
        )
      }
    }
    journal <- list(
      format = "BayesTools.RandomEffects.vignette-cache-transaction",
      version = 1L,
      had_cache = had_cache,
      previous_sha256 = previous_sha256,
      new_sha256 = .random_effects_vignette_sha256_file(paths$new)
    )
    .random_effects_vignette_write_journal(paths, journal)
    journal
  }, error = function(e) e)

  if(inherits(preparation_error, "error")){
    cleanup_error <- tryCatch({
      if(!file.exists(paths$journal)){
        .random_effects_vignette_cleanup_transaction(paths)
      }
      NULL
    }, error = function(e) e)
    if(inherits(cleanup_error, "error")){
      stop(
        conditionMessage(preparation_error),
        " Cleanup also failed: ",
        conditionMessage(cleanup_error),
        call. = FALSE
      )
    }
    stop(conditionMessage(preparation_error), call. = FALSE)
  }
  journal <- preparation_error

  installation_result <- tryCatch({
    .random_effects_vignette_replace_file(
      paths$new,
      cache_file,
      .platform,
      "atomically install the regenerated RandomEffects cache"
    )
    if(!identical(
      .random_effects_vignette_sha256_file(cache_file),
      journal$new_sha256
    )){
      stop(
        "Installed RandomEffects cache does not match its journal hash.",
        call. = FALSE
      )
    }
    installed_status <- validate_random_effects_vignette_cache(
      cache_file = cache_file,
      dependency_state = dependency_state
    )
    if(!isTRUE(installed_status$valid)){
      stop(
        "Installed RandomEffects cache failed validation: ",
        format_random_effects_vignette_cache_error(installed_status),
        call. = FALSE
      )
    }
    installed_status
  }, error = function(e) e)

  if(inherits(installation_result, "error")){
    recovery_error <- tryCatch({
      .random_effects_vignette_restore_previous(
        cache_file,
        paths,
        journal,
        .platform
      )
      NULL
    }, error = function(e) e)
    if(inherits(recovery_error, "error")){
      stop(
        conditionMessage(installation_result),
        " Automatic recovery failed: ",
        conditionMessage(recovery_error),
        " Transaction sidecars were retained.",
        call. = FALSE
      )
    }
    stop(conditionMessage(installation_result), call. = FALSE)
  }

  .random_effects_vignette_cleanup_transaction(paths)
  invisible(installation_result)
}

write_random_effects_vignette_cache <- function(
    models,
    cache_file = file.path("models", "RandomEffects.RDS"),
    generation,
    dependency_state = NULL,
    implementation_state = NULL,
    producer = NULL,
    project_root = NULL){
  if(!inherits(
    generation,
    "BayesTools_random_effects_vignette_cache_regeneration"
  ) || !identical(
    generation$format,
    "BayesTools.RandomEffects.vignette-cache-regeneration"
  )){
    stop(
      "Begin RandomEffects cache regeneration before writing the cache.",
      call. = FALSE
    )
  }
  if(!identical(
    names(generation),
    c("format", "dependencies", "implementation", "checkpoint_state")
  ) || !is.environment(generation$checkpoint_state) ||
      !exists(
        "used",
        envir = generation$checkpoint_state,
        inherits = FALSE
      ) ||
      !is.logical(generation$checkpoint_state$used) ||
      length(generation$checkpoint_state$used) != 1L ||
      is.na(generation$checkpoint_state$used)){
    stop(
      "RandomEffects cache regeneration checkpoint is invalid.",
      call. = FALSE
    )
  }
  if(isTRUE(generation$checkpoint_state$used)){
    stop(
      "RandomEffects cache regeneration checkpoint was already used.",
      call. = FALSE
    )
  }
  generation_dependency_error <-
    .random_effects_vignette_dependency_state_error(generation$dependencies)
  if(!is.null(generation_dependency_error)){
    stop(
      "RandomEffects cache regeneration checkpoint is invalid: ",
      generation_dependency_error,
      ". Restart regeneration.",
      call. = FALSE
    )
  }
  generation_implementation_error <-
    .random_effects_vignette_implementation_state_error(
      generation$implementation,
      generation$dependencies
    )
  if(!is.null(generation_implementation_error)){
    stop(
      "RandomEffects cache regeneration checkpoint is invalid: ",
      generation_implementation_error,
      ". Restart regeneration.",
      call. = FALSE
    )
  }
  if(is.null(dependency_state)){
    dependency_state <- random_effects_vignette_dependency_state(
      cache_file = cache_file,
      project_root = project_root
    )
  }
  dependency_error <-
    .random_effects_vignette_dependency_state_error(dependency_state)
  if(!is.null(dependency_error)){
    stop(
      "Cannot write RandomEffects cache: ",
      dependency_error,
      ".",
      call. = FALSE
    )
  }
  changed <- .random_effects_vignette_dependency_changes(
    generation$dependencies,
    dependency_state
  )
  if(length(changed) > 0L){
    stop(
      "RandomEffects dependencies changed during regeneration: ",
      paste(changed, collapse = ", "),
      ". Restart regeneration.",
      call. = FALSE
    )
  }
  if(is.null(implementation_state)){
    implementation_state <- .random_effects_vignette_current_implementation(
      dependencies = dependency_state,
      cache_file = cache_file,
      project_root = project_root
    )
  }
  implementation_error <-
    .random_effects_vignette_implementation_state_error(
      implementation_state,
      dependency_state
    )
  if(!is.null(implementation_error)){
    stop(
      "Cannot write RandomEffects cache: ",
      implementation_error,
      ".",
      call. = FALSE
    )
  }
  if(!identical(generation$implementation, implementation_state)){
    stop(
      paste0(
        "The verified BayesTools implementation changed during RandomEffects ",
        "cache regeneration. Reload the project and restart regeneration."
      ),
      call. = FALSE
    )
  }
  .random_effects_vignette_validate_models_or_stop(models)
  if(is.null(producer)){
    producer <- .random_effects_vignette_producer()
  }
  producer_error <- .random_effects_vignette_producer_error(producer)
  if(!is.null(producer_error)){
    stop(
      "Cannot write RandomEffects cache: ",
      producer_error,
      ".",
      call. = FALSE
    )
  }
  if(!identical(
    producer$fit_backend_fingerprint,
    implementation_state$fit_backend_fingerprint
  )){
    stop(
      paste0(
        "Cannot write RandomEffects cache: producer backend does not match ",
        "the verified current implementation."
      ),
      call. = FALSE
    )
  }
  model_hashes <- .random_effects_vignette_model_hashes(models)
  manifest <- .random_effects_vignette_manifest(
    dependency_state,
    model_hashes,
    producer
  )
  envelope <- list(manifest = manifest, models = models)
  status <- .random_effects_vignette_atomic_write(
    envelope = envelope,
    cache_file = cache_file,
    dependency_state = dependency_state
  )
  generation$checkpoint_state$used <- TRUE
  invisible(status)
}
