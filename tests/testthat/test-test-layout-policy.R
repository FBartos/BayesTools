skip_if_not_test_profile("unit")

test_that("live fitting is confined to the centralized fit file", {

  test_directory <- testthat::test_path()
  test_files <- list.files(
    test_directory,
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  centralized_file <- "test-00-model-fits.R"

  declares_fit_profile <- vapply(test_files, function(path){
    lines <- readLines(path, warn = FALSE, n = 20L)
    any(grepl(
      "^\\s*skip_if_not_test_profile\\([\"']fit[\"']\\)",
      lines
    ))
  }, logical(1))
  expect_setequal(
    basename(test_files[declares_fit_profile]),
    centralized_file
  )

  backend_patterns <- c(
    "rjags::jags\\.model\\s*\\(",
    "rjags::coda\\.samples\\s*\\(",
    "runjags::run\\.jags\\s*\\("
  )
  noncentralized <- test_files[basename(test_files) != centralized_file]
  violations <- unlist(lapply(noncentralized, function(path){
    lines <- readLines(path, warn = FALSE)
    matched <- unique(unlist(lapply(backend_patterns, function(pattern){
      grep(pattern, lines, value = TRUE)
    })))
    if(length(matched) == 0L){
      return(character())
    }
    paste0(basename(path), ": ", trimws(matched))
  }), use.names = FALSE)

  expect_identical(violations, character())
})

test_that("vignettes use the BayesTools marginal-likelihood contract", {

  vignette_directory <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(
    dir.exists(vignette_directory),
    paste0(
      "Repository vignette sources are not available in this ",
      "installed-package test context."
    )
  )

  vignette_files <- list.files(
    vignette_directory,
    pattern = "\\.Rmd$",
    full.names = TRUE
  )
  violations <- unlist(lapply(vignette_files, function(path){
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    matches <- grep(
      "bridgesampling::bf\\s*\\(",
      lines,
      value = TRUE
    )
    if(length(matches) == 0L){
      return(character())
    }
    paste0(basename(path), ": ", trimws(matches))
  }), use.names = FALSE)

  expect_identical(violations, character())
})

test_that("CI does not bypass vignette construction or checking", {

  workflow_file <- testthat::test_path(
    "..",
    "..",
    ".github",
    "workflows",
    "R-CMD-check.yaml"
  )
  skip_if_not(
    file.exists(workflow_file),
    "Repository workflow sources are not available in this installed-package context."
  )

  workflow <- readLines(
    workflow_file,
    warn = FALSE,
    encoding = "UTF-8"
  )
  bypasses <- workflow[
    grepl("--no-build-vignettes|--ignore-vignettes", workflow)
  ]

  expect_identical(bypasses, character())
  expect_true(any(grepl(
    'pattern = "\\\\.Rmd$"',
    workflow,
    fixed = TRUE
  )))
})

test_that("vignettes never fit models during ordinary rendering", {

  vignette_directory <- testthat::test_path("..", "..", "vignettes")
  skip_if_not(
    dir.exists(vignette_directory),
    paste0(
      "Repository vignette sources are not available in this ",
      "installed-package test context."
    )
  )

  vignette_files <- list.files(
    vignette_directory,
    pattern = "\\.Rmd$",
    full.names = TRUE
  )
  violations <- character()
  for(path in vignette_files){
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    chunk_header <- NULL
    for(line_number in seq_along(lines)){
      line <- lines[[line_number]]
      if(grepl("^```\\{r(?:[ ,}]|$)", line)){
        chunk_header <- line
        next
      }
      if(identical(line, "```")){
        chunk_header <- NULL
        next
      }
      if(is.null(chunk_header)){
        next
      }
      fits_model <- grepl(
        "(JAGS_fit|JAGS_bridgesampling|rstanarm::stan_lmer)\\s*\\(",
        line
      )
      if(fits_model && !grepl("eval\\s*=\\s*FALSE", chunk_header)){
        violations <- c(
          violations,
          paste0(
            basename(path), ":", line_number, ": ",
            trimws(line)
          )
        )
      }
    }
  }

  cache_files <- file.path(
    vignette_directory,
    paste0(tools::file_path_sans_ext(basename(vignette_files)), ".RDS")
  )
  missing_caches <- basename(cache_files[!file.exists(cache_files)])

  expect_identical(violations, character())
  expect_identical(missing_caches, character())
})


test_that("every test file is registered in a profile and declares a skip", {

  test_directory <- testthat::test_path()
  test_files <- list.files(
    test_directory,
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  contexts <- sub("^test-", "", sub("\\.R$", "", basename(test_files)))
  registered <- unique(unlist(
    bayestools_test_profile_contexts,
    use.names = FALSE
  ))

  expect_identical(sort(setdiff(contexts, registered)), character())
  expect_identical(sort(setdiff(registered, contexts)), character())

  skip_pattern <- paste(
    "skip_if_not_test_profile\\(",
    "skip_if_not_visual_tests\\(",
    "skip_if_not_visual_fixture_tests\\(",
    "skip_if_not_heavy_tests\\(",
    sep = "|"
  )
  missing_skip <- vapply(test_files, function(path){
    lines <- readLines(path, warn = FALSE, n = 20L)
    !any(grepl(skip_pattern, lines))
  }, logical(1))
  expect_identical(basename(test_files[missing_skip]), character())
})


# ---------------------------------------------------------------------------- #
# Build configuration: JAGS 4.x only
# ---------------------------------------------------------------------------- #

.layout_repository_file <- function(...){

  testthat::test_path("..", "..", ...)
}

.layout_write_lf <- function(lines, path){

  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), path)
}

# A fake JAGS installation: the headers the build scripts inspect
# (module/Module.h and a version.h declaring JAGS_MAJOR) and library
# directories. Nothing is compiled or linked against it.
.layout_fake_jags_tree <- function(root, header_major, include_dir = file.path("include", "JAGS"),
                                   lib_dir = "lib", libraries = character()){

  headers <- file.path(root, include_dir)
  dir.create(file.path(headers, "module"), recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(headers, "module", "Module.h"))
  .layout_write_lf(
    c("#ifndef JAGS_VERSION_H_", "#define JAGS_VERSION_H_", "",
      paste("#define JAGS_MAJOR", header_major), "", "#endif"),
    file.path(headers, "version.h")
  )
  dir.create(file.path(root, lib_dir), recursive = TRUE, showWarnings = FALSE)
  if(length(libraries) > 0L){
    file.create(file.path(root, lib_dir, libraries))
  }

  normalizePath(root, winslash = "/", mustWork = TRUE)
}

.layout_jags_env_unset <- function(){

  c(JAGS_ROOT = NA, JAGS_PREFIX = NA, JAGS_INCLUDE_DIR = NA, JAGS_LIB_DIR = NA,
    JAGS_VERSION = NA, PKG_CONFIG_PATH = NA)
}

.layout_run_tool <- function(command, args, directory, env){

  withr::local_envvar(env)
  withr::local_dir(directory)
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status", exact = TRUE)

  list(
    status = if(is.null(status)) 0L else as.integer(status),
    output = paste(output, collapse = "\n")
  )
}

test_that("configure accepts only JAGS 4.x from every version source", {

  skip_on_cran()
  configure_file <- .layout_repository_file("configure")
  makevars_file  <- .layout_repository_file("src", "Makevars.in")
  skip_if_not(
    file.exists(configure_file) && file.exists(makevars_file),
    "Repository build sources are not available in this installed-package context."
  )
  sh <- unname(Sys.which("sh"))
  skip_if(!nzchar(sh), "A POSIX shell is required to run configure.")

  work <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  package_dir <- file.path(work, "pkg")
  dir.create(file.path(package_dir, "src"), recursive = TRUE)
  file.copy(configure_file, file.path(package_dir, "configure"))
  file.copy(makevars_file, file.path(package_dir, "src", "Makevars.in"))

  jags4    <- .layout_fake_jags_tree(file.path(work, "jags4"), 4L)
  jags5    <- .layout_fake_jags_tree(file.path(work, "jags5"), 5L)
  prefix5  <- .layout_fake_jags_tree(file.path(work, "JAGS-5.0.0"), 4L)
  prefix43 <- .layout_fake_jags_tree(file.path(work, "JAGS-4.3.2"), 4L)

  # pkg-config stubs placed first on PATH: one without JAGS metadata (so a
  # system jags.pc cannot leak into the other version sources) and one that
  # reports JAGS 5.0.0 with the JAGS 4 test headers.
  stub_none <- file.path(work, "stub-none")
  stub_jags5 <- file.path(work, "stub-jags5")
  dir.create(stub_none)
  dir.create(stub_jags5)
  .layout_write_lf(c("#!/bin/sh", "exit 1"), file.path(stub_none, "pkg-config"))
  .layout_write_lf(c(
    "#!/bin/sh",
    "case \"$1\" in",
    "  --exists) exit 0 ;;",
    paste0("  --cflags) echo \"-I", jags4, "/include/JAGS\" ;;"),
    paste0("  --libs) echo \"-L", jags4, "/lib -ljags\" ;;"),
    "  --modversion) echo \"5.0.0\" ;;",
    "esac"
  ), file.path(stub_jags5, "pkg-config"))
  Sys.chmod(file.path(c(stub_none, stub_jags5), "pkg-config"), mode = "0755")

  run_configure <- function(args = character(), env = character(), stub = stub_none){

    unlink(file.path(package_dir, "src", "Makevars"))
    path <- paste(stub, Sys.getenv("PATH"), sep = .Platform$path.sep)
    result <- .layout_run_tool(
      sh, c("./configure", args), package_dir,
      c(.layout_jags_env_unset(), env, PATH = path)
    )
    result$makevars <- file.exists(file.path(package_dir, "src", "Makevars"))
    result
  }
  expect_rejected <- function(result, reported){

    expect_false(identical(result$status, 0L), info = result$output)
    expect_false(result$makevars, info = result$output)
    expect_match(
      result$output,
      paste0("BayesTools requires JAGS 4.x (>= 4.3.0, < 5.0.0); ", reported),
      fixed = TRUE
    )
  }
  expect_accepted <- function(result, detected){

    expect_identical(result$status, 0L, info = result$output)
    expect_true(result$makevars, info = result$output)
    expect_match(result$output, paste("detected JAGS", detected), fixed = TRUE)
  }
  jags4_dirs <- c(
    JAGS_INCLUDE_DIR = file.path(jags4, "include", "JAGS"),
    JAGS_LIB_DIR     = file.path(jags4, "lib")
  )

  expect_rejected(run_configure(env = c(JAGS_VERSION = "5.0.0", jags4_dirs)), "JAGS_VERSION reported 5.0.0")
  expect_rejected(run_configure(env = c(JAGS_VERSION = "4.2.0", jags4_dirs)), "JAGS_VERSION reported 4.2.0")
  expect_accepted(run_configure(env = c(JAGS_VERSION = "4.3.2", jags4_dirs)), "4.3.2 with JAGS_VERSION")
  expect_rejected(
    run_configure(args = c(
      "--with-jags-version=5.0.0",
      paste0("--with-jags-includedir=", jags4_dirs[["JAGS_INCLUDE_DIR"]]),
      paste0("--with-jags-libdir=", jags4_dirs[["JAGS_LIB_DIR"]])
    )),
    "configure argument reported 5.0.0"
  )
  expect_rejected(run_configure(args = paste0("--with-jags-prefix=", prefix5)), "JAGS prefix reported 5.0.0")
  expect_accepted(run_configure(args = paste0("--with-jags-prefix=", prefix43)), "4.3.2 with JAGS prefix")
  expect_rejected(run_configure(stub = stub_jags5), "pkg-config reported 5.0.0")
  # Headers declaring JAGS 5 stop configure before the jags_version() probe.
  expect_rejected(run_configure(args = paste0("--with-jags-prefix=", jags5)), "the JAGS headers in")
})

test_that("Windows builds select only JAGS 4.x installations", {

  skip_on_cran()
  makevars_file <- .layout_repository_file("src", "Makevars.win.common")
  skip_if_not(
    file.exists(makevars_file),
    "Repository build sources are not available in this installed-package context."
  )
  make <- unname(Sys.which("make"))
  skip_if(!nzchar(make), "GNU make is required to evaluate src/Makevars.win.common.")

  work <- normalizePath(withr::local_tempdir(), winslash = "/", mustWork = TRUE)
  skip_if(grepl("[[:space:]]", work), "Candidate JAGS roots are space-separated make words.")
  file.copy(makevars_file, file.path(work, "Makevars.win.common"))
  .layout_write_lf(c(
    "SHLIB = probe.dll",
    "print-selection:",
    "\t@echo \"selected $(notdir $(JAGS_ROOT)) version $(JAGS_VERSION) link -l$(JAGS_LINK)\"",
    "include Makevars.win.common"
  ), file.path(work, "probe.mk"))

  install_root <- function(version, header_major = sub("\\..*$", "", version)){

    major <- sub("\\..*$", "", version)
    .layout_fake_jags_tree(
      file.path(work, "JAGS", paste0("JAGS-", version)), header_major,
      include_dir = "include", lib_dir = file.path("x64", "lib"),
      libraries = c(paste0("libjags-", major, ".dll.a"), "libjrmath-0.dll.a")
    )
  }
  roots <- vapply(c("4.2.0", "4.3.1", "4.3.2", "5.0.0"), install_root, character(1))

  run_make <- function(candidates, vars = character()){

    .layout_run_tool(
      make,
      c("-s", "-f", "probe.mk", "R_ARCH=/x64",
        shQuote(paste0("JAGS_INSTALL_ROOTS=", paste(roots[candidates], collapse = " "))),
        vapply(vars, shQuote, character(1)), "print-selection"),
      work,
      .layout_jags_env_unset()
    )
  }

  selected <- run_make(c("4.2.0", "4.3.1", "5.0.0"))
  expect_identical(selected$status, 0L, info = selected$output)
  expect_match(selected$output, "selected JAGS-4.3.1 version 4.3.1 link -ljags-4", fixed = TRUE)

  newest <- run_make(c("4.3.1", "4.3.2", "5.0.0"))
  expect_identical(newest$status, 0L, info = newest$output)
  expect_match(newest$output, "selected JAGS-4.3.2 version 4.3.2", fixed = TRUE)

  only_other <- run_make("5.0.0")
  expect_false(identical(only_other$status, 0L), info = only_other$output)
  expect_match(
    only_other$output,
    "BayesTools requires JAGS 4.x (>= 4.3.0, < 5.0.0), but only other JAGS versions were found: JAGS-5.0.0",
    fixed = TRUE
  )

  explicit_root <- run_make("4.3.1", paste0("JAGS_ROOT=", roots[["5.0.0"]]))
  expect_false(identical(explicit_root$status, 0L), info = explicit_root$output)
  expect_match(explicit_root$output, "declare JAGS_MAJOR 5", fixed = TRUE)

  explicit_version <- run_make("4.3.1", "JAGS_VERSION=5.0.0")
  expect_false(identical(explicit_version$status, 0L), info = explicit_version$output)
  expect_match(explicit_version$output, "JAGS_VERSION/JAGS_ROOT reported 5.0.0", fixed = TRUE)
})

# Steps of a GitHub Actions workflow: name, shell, and the lines of the run
# block (a block scalar ends at the first non-blank line that is not indented
# deeper than its `run:` key).
.layout_workflow_steps <- function(path){

  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  indentation <- function(x) nchar(x) - nchar(sub("^ +", "", x))
  starts <- grep("^\\s*- (name|uses):", lines)
  ends <- c(starts[-1L] - 1L, length(lines))

  lapply(seq_along(starts), function(i){

    step <- lines[starts[i]:ends[i]]
    field <- function(key){
      value <- step[grepl(paste0("^\\s*(- )?", key, ":"), step)]
      if(length(value) == 0L) NA_character_ else trimws(sub(paste0("^\\s*(- )?", key, ":"), "", value[1L]))
    }
    run <- character()
    run_line <- grep("^\\s*run:", step)
    if(length(run_line) > 0L){
      run_line <- run_line[1L]
      inline <- field("run")
      if(inline %in% c("|", ">")){
        key_indent <- indentation(step[run_line])
        for(line in step[-seq_len(run_line)]){
          if(nzchar(trimws(line)) && indentation(line) <= key_indent){
            break
          }
          run <- c(run, line)
        }
      }else{
        run <- inline
      }
    }

    list(name = field("name"), shell = field("shell"), run = run)
  })
}

.layout_workflow_files <- function(){

  workflow_dir <- .layout_repository_file(".github", "workflows")
  skip_if_not(
    dir.exists(workflow_dir),
    "Repository workflow sources are not available in this installed-package context."
  )

  list.files(workflow_dir, pattern = "\\.ya?ml$", full.names = TRUE)
}

test_that("macOS CI builds against a cached JAGS 4 installer, never Homebrew's jags", {

  for(path in .layout_workflow_files()){
    workflow <- readLines(path, warn = FALSE, encoding = "UTF-8")
    expect_false(any(grepl("brew install jags", workflow, fixed = TRUE)), info = basename(path))

    steps <- .layout_workflow_steps(path)
    step_names <- vapply(steps, `[[`, character(1), "name")
    install <- steps[step_names %in% "Install JAGS (macOS)"]
    if(length(install) == 0L){
      next
    }

    cache <- steps[step_names %in% "Cache JAGS installer (macOS)"]
    expect_length(cache, 1L)
    cache_lines <- workflow[grepl("^\\s*path:\\s*~/jags-installer\\s*$", workflow)]
    expect_length(cache_lines, 1L)
    install_run <- paste(install[[1L]]$run, collapse = "\n")
    expect_match(install_run, "~/jags-installer/JAGS-4\\.3\\.[0-9]+\\.pkg", info = basename(path))
    expect_match(install_run, "sudo installer -pkg", fixed = TRUE, info = basename(path))
  }
})

test_that("CI installs BayesTools in a process that has loaded no packages", {

  # A namespace loaded before devtools::install() in the same process (RoBMA
  # loads BayesTools) keeps the old BayesTools DLL open, which Windows then
  # refuses to overwrite. The loaders are matched as names because they are
  # also passed as functions, e.g. sapply(packages, requireNamespace).
  loads_packages <- "\\b(requireNamespace|loadNamespace|library|require)\\b"
  install_steps <- 0L
  for(path in .layout_workflow_files()){
    for(step in .layout_workflow_steps(path)){
      run <- paste(step$run, collapse = "\n")
      if(!grepl("devtools::install(", run, fixed = TRUE)){
        next
      }
      install_steps <- install_steps + 1L
      expect_false(grepl(loads_packages, run, perl = TRUE), info = paste0(basename(path), ": ", step$name))
    }
  }

  expect_gt(install_steps, 0L)
})

test_that("DESCRIPTION declares the supported JAGS range", {

  description_file <- .layout_repository_file("DESCRIPTION")
  skip_if_not(
    file.exists(description_file),
    "Repository DESCRIPTION is not available in this installed-package context."
  )

  requirements <- unname(read.dcf(description_file, fields = "SystemRequirements")[1, 1])
  expect_match(requirements, "JAGS (>= 4.3.0, < 5.0.0)", fixed = TRUE)
})
