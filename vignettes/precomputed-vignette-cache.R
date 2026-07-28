precomputed_vignette_cache_schema <- function(vignette){
  schemas <- list(
    ComparisonR = c(
      ttest_model_H0 = "BayesTools_fit",
      ttest_model_Hp = "BayesTools_fit",
      marglik_model_H0 = "BayesTools_marglik",
      marglik_model_Hp = "BayesTools_marglik"
    ),
    SpikeAndSlab = c(
      M0 = "BayesTools_fit",
      M1 = "BayesTools_fit",
      marglik_model_H0 = "BayesTools_marglik",
      marglik_model_H1 = "BayesTools_marglik",
      MS = "BayesTools_fit"
    )
  )
  if(!is.character(vignette) ||
      length(vignette) != 1L ||
      is.na(vignette) ||
      !vignette %in% names(schemas)){
    stop(
      "'vignette' must be one of: ComparisonR, SpikeAndSlab.",
      call. = FALSE
    )
  }
  schemas[[vignette]]
}

precomputed_vignette_cache_names <- function(vignette){
  names(precomputed_vignette_cache_schema(vignette))
}

.precomputed_vignette_cache_format <- function(){
  "BayesTools.precomputed-vignette-cache"
}

.precomputed_vignette_cache_version <- function(){
  1L
}

.precomputed_vignette_cache_packages <- function(){
  c("BayesTools", "runjags", "rjags", "coda", "bridgesampling")
}

.precomputed_vignette_cache_file <- function(vignette){
  file.path("vignettes", paste0(vignette, ".RDS"))
}

.precomputed_vignette_md5_raw <- function(value){
  hash_file <- tempfile("BayesTools-vignette-hash-")
  on.exit(unlink(hash_file), add = TRUE)
  writeBin(value, hash_file)
  unname(tools::md5sum(hash_file))
}

.precomputed_vignette_text_md5 <- function(path){
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  text <- paste(lines, collapse = "\n")
  if(length(lines) > 0L){
    text <- paste0(text, "\n")
  }
  .precomputed_vignette_md5_raw(charToRaw(enc2utf8(text)))
}

.precomputed_vignette_description_md5 <- function(path){
  description <- read.dcf(path, all = TRUE)
  description[c("Author", "Built", "Packaged")] <- NULL
  description <- description[order(names(description), method = "radix")]
  values <- vapply(description, function(value){
    gsub("[[:space:]]+", " ", trimws(value))
  }, character(1))
  text <- paste(names(values), values, sep = ": ", collapse = "\n")
  .precomputed_vignette_md5_raw(
    charToRaw(enc2utf8(paste0(text, "\n")))
  )
}

.precomputed_vignette_project_root <- function(
    cache_file, project_root = NULL){
  if(is.null(project_root)){
    project_root <- file.path(dirname(cache_file), "..")
  }
  project_root <- tryCatch(
    normalizePath(project_root, winslash = "/", mustWork = TRUE),
    error = function(e) ""
  )
  if(!nzchar(project_root) ||
      !file.exists(file.path(project_root, "DESCRIPTION")) ||
      !file.exists(file.path(project_root, "NAMESPACE"))){
    stop(
      "Could not resolve the BayesTools project root for the vignette cache.",
      call. = FALSE
    )
  }
  project_root
}

.precomputed_vignette_relative_paths <- function(paths, project_root){
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  prefix <- paste0(project_root, "/")
  if(any(!startsWith(paths, prefix))){
    stop(
      "Vignette cache dependencies must be inside the project root.",
      call. = FALSE
    )
  }
  substring(paths, nchar(prefix) + 1L)
}

.precomputed_vignette_source_paths <- function(vignette, project_root){
  vignette_paths <- file.path(
    project_root,
    "vignettes",
    c(
      paste0(vignette, ".Rmd"),
      "precomputed-vignette-cache.R"
    )
  )
  root_paths <- file.path(project_root, c("DESCRIPTION", "NAMESPACE"))

  data_dir <- file.path(project_root, "data")
  data_paths <- if(dir.exists(data_dir)){
    list.files(
      data_dir,
      pattern = "\\.(?:RData|rda|rds)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  }else{
    character()
  }

  r_dir <- file.path(project_root, "R")
  r_paths <- if(dir.exists(r_dir)){
    list.files(
      r_dir,
      pattern = "\\.[Rr]$",
      recursive = TRUE,
      full.names = TRUE
    )
  }else{
    character()
  }

  src_dir <- file.path(project_root, "src")
  src_paths <- if(dir.exists(src_dir)){
    list.files(src_dir, recursive = TRUE, full.names = TRUE)
  }else{
    character()
  }
  native_extensions <- c(
    "c", "cc", "cpp", "cxx", "f", "f77", "f90", "f95",
    "h", "hh", "hpp", "hxx", "inc"
  )
  src_paths <- src_paths[
    tolower(tools::file_ext(src_paths)) %in% native_extensions |
      grepl("^Makevars\\..+$", basename(src_paths))
  ]

  paths <- unique(c(
    vignette_paths,
    root_paths,
    data_paths,
    r_paths,
    src_paths
  ))
  missing <- paths[!file.exists(paths)]
  if(length(missing) > 0L){
    stop(
      "Missing vignette cache source dependency: ",
      paste(basename(missing), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  paths[file.info(paths)$isdir %in% FALSE]
}

.precomputed_vignette_source_hashes <- function(vignette, project_root){
  paths <- .precomputed_vignette_source_paths(vignette, project_root)
  labels <- .precomputed_vignette_relative_paths(paths, project_root)
  source_order <- order(enc2utf8(labels), method = "radix")
  paths <- paths[source_order]
  labels <- labels[source_order]
  hashes <- vapply(
    seq_along(paths),
    function(index){
      if(identical(labels[[index]], "DESCRIPTION")){
        return(.precomputed_vignette_description_md5(paths[[index]]))
      }
      path <- paths[[index]]
      if(tolower(tools::file_ext(path)) %in% c("rdata", "rda", "rds")){
        return(unname(tools::md5sum(path)))
      }
      .precomputed_vignette_text_md5(path)
    },
    character(1)
  )
  stats::setNames(unname(hashes), labels)
}

.precomputed_vignette_package_versions <- function(){
  packages <- .precomputed_vignette_cache_packages()
  versions <- vapply(packages, function(package){
    tryCatch(
      as.character(utils::packageVersion(package)),
      error = function(e) NA_character_
    )
  }, character(1))
  if(anyNA(versions)){
    stop(
      "Missing vignette cache producer package(s): ",
      paste(names(versions)[is.na(versions)], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  versions
}

precomputed_vignette_cache_state <- function(
    vignette,
    cache_file = .precomputed_vignette_cache_file(vignette),
    project_root = NULL){
  precomputed_vignette_cache_schema(vignette)
  project_root <- .precomputed_vignette_project_root(
    cache_file,
    project_root
  )
  list(
    source_hashes =
      .precomputed_vignette_source_hashes(vignette, project_root),
    r_version = as.character(getRversion()),
    package_versions = .precomputed_vignette_package_versions()
  )
}

.precomputed_vignette_named_hashes_valid <- function(value){
  is.character(value) &&
    length(value) > 0L &&
    !is.null(names(value)) &&
    !anyNA(names(value)) &&
    all(nzchar(names(value))) &&
    !anyDuplicated(names(value)) &&
    !anyNA(value) &&
    all(grepl("^[[:xdigit:]]{32}$", value))
}

.precomputed_vignette_scalar_character <- function(value){
  is.character(value) &&
    length(value) == 1L &&
    !is.na(value) &&
    nzchar(value)
}

.precomputed_vignette_state_error <- function(value){
  if(!is.list(value) ||
      !identical(
        names(value),
        c("source_hashes", "r_version", "package_versions")
      )){
    return("dependency state has an unexpected schema")
  }
  if(!.precomputed_vignette_named_hashes_valid(value$source_hashes)){
    return("source fingerprints are invalid")
  }
  if(!.precomputed_vignette_scalar_character(value$r_version)){
    return("R version is invalid")
  }
  if(!is.character(value$package_versions) ||
      !identical(
        names(value$package_versions),
        .precomputed_vignette_cache_packages()
      ) ||
      anyNA(value$package_versions) ||
      any(!nzchar(value$package_versions))){
    return("package versions are invalid")
  }
  NULL
}

.precomputed_vignette_object_status <- function(objects, vignette){
  schema <- precomputed_vignette_cache_schema(vignette)
  if(!is.list(objects) ||
      is.null(names(objects)) ||
      anyNA(names(objects)) ||
      any(!nzchar(names(objects))) ||
      anyDuplicated(names(objects))){
    return("objects must be a uniquely named list")
  }
  if(!identical(names(objects), names(schema))){
    missing <- setdiff(names(schema), names(objects))
    extra <- setdiff(names(objects), names(schema))
    if(length(missing) > 0L){
      return(paste0("objects are missing: ", paste(missing, collapse = ", ")))
    }
    if(length(extra) > 0L){
      return(paste0(
        "objects are unexpected: ",
        paste(extra, collapse = ", ")
      ))
    }
    return("objects are not in the required order")
  }
  invalid <- names(schema)[
    !vapply(names(schema), function(name){
      inherits(objects[[name]], schema[[name]])
    }, logical(1))
  ]
  if(length(invalid) > 0L){
    return(paste0(
      "objects have unexpected classes: ",
      paste(invalid, collapse = ", ")
    ))
  }
  NULL
}

.precomputed_vignette_serialize_objects <- function(objects, vignette){
  object_error <- .precomputed_vignette_object_status(objects, vignette)
  if(!is.null(object_error)){
    stop(
      "Cannot write ", vignette, " vignette cache: ", object_error, ".",
      call. = FALSE
    )
  }
  payloads <- lapply(objects, serialize, connection = NULL, version = 3)
  names(payloads) <- names(objects)
  payloads
}

.precomputed_vignette_payload_error <- function(payloads, vignette){
  schema <- precomputed_vignette_cache_schema(vignette)
  if(!is.list(payloads) || !identical(names(payloads), names(schema))){
    return("serialized payloads do not match the required object schema")
  }
  if(!all(vapply(payloads, is.raw, logical(1)))){
    return("serialized payloads must be raw vectors")
  }
  NULL
}

.precomputed_vignette_payload_hashes <- function(payloads){
  hashes <- vapply(
    payloads,
    .precomputed_vignette_md5_raw,
    character(1)
  )
  stats::setNames(unname(hashes), names(payloads))
}

.precomputed_vignette_unserialize_payloads <- function(payloads, vignette){
  payload_error <- .precomputed_vignette_payload_error(payloads, vignette)
  if(!is.null(payload_error)){
    stop(payload_error, call. = FALSE)
  }
  objects <- lapply(payloads, unserialize)
  names(objects) <- names(payloads)
  object_error <- .precomputed_vignette_object_status(objects, vignette)
  if(!is.null(object_error)){
    stop(object_error, call. = FALSE)
  }
  objects
}

.precomputed_vignette_producer <- function(state){
  state_error <- .precomputed_vignette_state_error(state)
  if(!is.null(state_error)){
    stop("Cannot record vignette cache producer: ", state_error, ".", call. = FALSE)
  }
  list(
    generated_at_utc = format(
      Sys.time(),
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    ),
    r_version = state$r_version,
    package_versions = state$package_versions
  )
}

.precomputed_vignette_producer_error <- function(value){
  if(!is.list(value) ||
      !identical(
        names(value),
        c("generated_at_utc", "r_version", "package_versions")
      )){
    return("producer metadata has an unexpected schema")
  }
  if(!.precomputed_vignette_scalar_character(value$generated_at_utc) ||
      !grepl(
        "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$",
        value$generated_at_utc
      )){
    return("producer timestamp is invalid")
  }
  state <- list(
    source_hashes = c(source = paste(rep("0", 32L), collapse = "")),
    r_version = value$r_version,
    package_versions = value$package_versions
  )
  state_error <- .precomputed_vignette_state_error(state)
  if(!is.null(state_error)){
    return(sub("^dependency state ", "producer ", state_error))
  }
  NULL
}

.precomputed_vignette_manifest_error <- function(manifest, vignette){
  schema <- precomputed_vignette_cache_schema(vignette)
  if(!is.list(manifest) ||
      !identical(
        names(manifest),
        c(
          "format", "version", "vignette", "object_schema",
          "source_hashes", "payload_hashes", "producer"
        )
      )){
    return("manifest has an unexpected schema")
  }
  if(!identical(manifest$format, .precomputed_vignette_cache_format())){
    return("manifest format is unsupported")
  }
  if(!identical(manifest$version, .precomputed_vignette_cache_version())){
    return("manifest version is unsupported")
  }
  if(!identical(manifest$vignette, vignette)){
    return("manifest belongs to a different vignette")
  }
  if(!identical(manifest$object_schema, schema)){
    return("object schema does not match the current vignette contract")
  }
  if(!.precomputed_vignette_named_hashes_valid(manifest$source_hashes)){
    return("source fingerprints are invalid")
  }
  if(!.precomputed_vignette_named_hashes_valid(manifest$payload_hashes) ||
      !identical(names(manifest$payload_hashes), names(schema))){
    return("payload fingerprints are invalid")
  }
  producer_error <- .precomputed_vignette_producer_error(manifest$producer)
  if(!is.null(producer_error)){
    return(producer_error)
  }
  NULL
}

.precomputed_vignette_status <- function(vignette, cache_file){
  list(
    valid = FALSE,
    vignette = vignette,
    cache_file = cache_file,
    reason = NULL,
    error = NULL,
    stale = character(),
    manifest = NULL,
    objects = NULL
  )
}

validate_precomputed_vignette_cache <- function(
    vignette,
    cache_file = .precomputed_vignette_cache_file(vignette),
    project_root = NULL,
    state = NULL,
    load = FALSE){
  precomputed_vignette_cache_schema(vignette)
  status <- .precomputed_vignette_status(vignette, cache_file)
  if(!file.exists(cache_file)){
    status$reason <- "missing"
    status$error <- paste0(
      "Missing precomputed ", vignette, " vignette cache at '",
      cache_file, "'."
    )
    return(status)
  }

  envelope <- tryCatch(
    suppressWarnings(readRDS(cache_file)),
    error = function(e) e
  )
  if(inherits(envelope, "error")){
    status$reason <- "unreadable"
    status$error <- paste0(
      "Could not read precomputed ", vignette, " vignette cache: ",
      conditionMessage(envelope), "."
    )
    return(status)
  }
  if(!is.list(envelope) ||
      !identical(names(envelope), c("manifest", "payloads"))){
    status$reason <- "invalid"
    status$error <- paste0(
      "Precomputed ", vignette,
      " vignette cache has an unexpected envelope schema."
    )
    return(status)
  }

  manifest_error <-
    .precomputed_vignette_manifest_error(envelope$manifest, vignette)
  if(!is.null(manifest_error)){
    status$reason <- "invalid"
    status$error <- paste0(
      "Precomputed ", vignette, " vignette cache is invalid: ",
      manifest_error, "."
    )
    return(status)
  }
  status$manifest <- envelope$manifest

  payload_error <-
    .precomputed_vignette_payload_error(envelope$payloads, vignette)
  if(!is.null(payload_error)){
    status$reason <- "invalid"
    status$error <- paste0(
      "Precomputed ", vignette, " vignette cache is invalid: ",
      payload_error, "."
    )
    return(status)
  }
  payload_hashes <-
    .precomputed_vignette_payload_hashes(envelope$payloads)
  changed_payloads <- names(payload_hashes)[
    payload_hashes != envelope$manifest$payload_hashes
  ]
  if(length(changed_payloads) > 0L){
    status$reason <- "corrupt"
    status$error <- paste0(
      "Precomputed ", vignette,
      " vignette cache has modified payloads: ",
      paste(changed_payloads, collapse = ", "), "."
    )
    return(status)
  }

  objects <- tryCatch(
    .precomputed_vignette_unserialize_payloads(
      envelope$payloads,
      vignette
    ),
    error = function(e) e
  )
  if(inherits(objects, "error")){
    status$reason <- "corrupt"
    status$error <- paste0(
      "Could not restore precomputed ", vignette,
      " vignette objects: ", conditionMessage(objects), "."
    )
    return(status)
  }

  if(is.null(state)){
    state <- tryCatch(
      precomputed_vignette_cache_state(
        vignette = vignette,
        cache_file = cache_file,
        project_root = project_root
      ),
      error = function(e) e
    )
    if(inherits(state, "error")){
      status$reason <- "dependency-error"
      status$error <- paste0(
        "Could not fingerprint current ", vignette,
        " vignette dependencies: ", conditionMessage(state)
      )
      return(status)
    }
  }
  state_error <- .precomputed_vignette_state_error(state)
  if(!is.null(state_error)){
    status$reason <- "dependency-error"
    status$error <- paste0(
      "Could not fingerprint current ", vignette,
      " vignette dependencies: ", state_error, "."
    )
    return(status)
  }

  if(!identical(state$source_hashes, envelope$manifest$source_hashes)){
    status$stale <- c(status$stale, "source files")
  }
  if(!identical(state$r_version, envelope$manifest$producer$r_version)){
    status$stale <- c(status$stale, "R version")
  }
  if(!identical(
    state$package_versions,
    envelope$manifest$producer$package_versions
  )){
    status$stale <- c(status$stale, "package versions")
  }
  if(length(status$stale) > 0L){
    status$reason <- "stale"
    status$error <- paste0(
      "Precomputed ", vignette,
      " vignette cache is stale; changed dependencies: ",
      paste(status$stale, collapse = ", "), "."
    )
    return(status)
  }

  status$valid <- TRUE
  if(isTRUE(load)){
    status$objects <- objects
  }
  status
}

load_precomputed_vignette_cache <- function(
    vignette,
    cache_file = .precomputed_vignette_cache_file(vignette),
    project_root = NULL,
    state = NULL){
  status <- validate_precomputed_vignette_cache(
    vignette = vignette,
    cache_file = cache_file,
    project_root = project_root,
    state = state,
    load = TRUE
  )
  if(!isTRUE(status$valid)){
    stop(status$error, call. = FALSE)
  }
  status$objects
}

.precomputed_vignette_install_candidate <- function(
    candidate, cache_file){
  if(!file.exists(cache_file)){
    if(!file.rename(candidate, cache_file)){
      stop("Could not install the new vignette cache.", call. = FALSE)
    }
    return(invisible(TRUE))
  }

  backup <- tempfile(
    paste0(".", basename(cache_file), "-backup-"),
    tmpdir = dirname(cache_file)
  )
  on.exit(unlink(backup), add = TRUE)
  if(!file.rename(cache_file, backup)){
    stop(
      "Could not prepare the existing vignette cache for replacement.",
      call. = FALSE
    )
  }
  installed <- file.rename(candidate, cache_file)
  if(!installed){
    restored <- file.rename(backup, cache_file)
    if(!restored){
      stop(
        "Could not install the new vignette cache or restore the previous cache.",
        call. = FALSE
      )
    }
    stop(
      "Could not install the new vignette cache; the previous cache was restored.",
      call. = FALSE
    )
  }
  unlink(backup)
  invisible(TRUE)
}

write_precomputed_vignette_cache <- function(
    objects,
    vignette,
    cache_file = .precomputed_vignette_cache_file(vignette),
    project_root = NULL,
    state = NULL){
  precomputed_vignette_cache_schema(vignette)
  if(is.null(state)){
    state <- precomputed_vignette_cache_state(
      vignette = vignette,
      cache_file = cache_file,
      project_root = project_root
    )
  }
  state_error <- .precomputed_vignette_state_error(state)
  if(!is.null(state_error)){
    stop(
      "Cannot write ", vignette, " vignette cache: ",
      state_error, ".",
      call. = FALSE
    )
  }

  payloads <- .precomputed_vignette_serialize_objects(objects, vignette)
  producer <- .precomputed_vignette_producer(state)
  envelope <- list(
    manifest = list(
      format = .precomputed_vignette_cache_format(),
      version = .precomputed_vignette_cache_version(),
      vignette = vignette,
      object_schema = precomputed_vignette_cache_schema(vignette),
      source_hashes = state$source_hashes,
      payload_hashes = .precomputed_vignette_payload_hashes(payloads),
      producer = producer
    ),
    payloads = payloads
  )

  cache_dir <- dirname(cache_file)
  if(!dir.exists(cache_dir) &&
      !dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)){
    stop(
      "Could not create vignette cache directory: ", cache_dir, ".",
      call. = FALSE
    )
  }
  candidate <- tempfile(
    paste0(".", basename(cache_file), "-"),
    tmpdir = cache_dir,
    fileext = ".tmp"
  )
  on.exit(unlink(candidate), add = TRUE)
  saveRDS(envelope, candidate, version = 3, compress = "xz")

  candidate_status <- validate_precomputed_vignette_cache(
    vignette = vignette,
    cache_file = candidate,
    state = state,
    load = TRUE
  )
  if(!isTRUE(candidate_status$valid)){
    stop(
      "Prepared ", vignette, " vignette cache failed validation: ",
      candidate_status$error,
      call. = FALSE
    )
  }
  .precomputed_vignette_install_candidate(candidate, cache_file)
  invisible(cache_file)
}
