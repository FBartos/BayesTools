skip_if_not_test_profile("meta")

test_that("live fitting is confined to the centralized fit file and the certification cases", {

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
    "\\.bt_rjags_fit\\s*\\("
  )
  # the certification cases fit live outside the cached fixtures (testing.md)
  noncentralized <- test_files[basename(test_files) != centralized_file &
                                 !startsWith(basename(test_files), "test-certification-")]
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

# Rendering fits only missing or stale caches (see .agents/instructions/vignettes.md):
# every fitting call belongs to a chunk with a `cache_fits` option, or is not
# evaluated.
test_that("vignettes fit models only in cached chunks", {

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
      if(fits_model && !grepl("eval\\s*=\\s*FALSE|cache_fits\\s*=", chunk_header)){
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

  expect_identical(violations, character())
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
    sep = "|"
  )
  missing_skip <- vapply(test_files, function(path){
    lines <- readLines(path, warn = FALSE, n = 20L)
    !any(grepl(skip_pattern, lines))
  }, logical(1))
  expect_identical(basename(test_files[missing_skip]), character())
})

test_that("every test file's profile gate matches its registry lanes", {

  # The runner selects files by registry lane while each file skips itself
  # unless its gate lane is active: a file registered under a lane its gate
  # does not declare runs nowhere in that lane, and a gated lane without
  # registration never selects the file.
  gate_names <- c(
    "skip_if_not_test_profile", "skip_if_not_visual_tests",
    "skip_if_not_visual_fixture_tests"
  )
  gate_lanes <- function(call){

    switch(
      as.character(call[[1L]]),
      skip_if_not_test_profile         = bayestools_normalize_test_profiles(eval(call[[2L]], baseenv())),
      skip_if_not_visual_tests         = "visual",
      skip_if_not_visual_fixture_tests = "visual-fixture"
    )
  }
  is_gate <- function(expression){

    is.call(expression) && is.name(expression[[1L]]) &&
      as.character(expression[[1L]]) %in% gate_names
  }
  gate_line <- paste0("^\\s*(", paste(gate_names, collapse = "|"), ")\\(")

  test_files <- list.files(
    testthat::test_path(),
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  problems <- character()
  for(path in test_files){
    context <- sub("^test-", "", sub("\\.R$", "", basename(path)))
    registry_lanes <- names(bayestools_test_profile_contexts)[vapply(
      bayestools_test_profile_contexts,
      function(contexts) context %in% contexts,
      logical(1)
    )]

    file_gate <- Filter(is_gate, as.list(parse(path, keep.source = FALSE)))
    if(length(file_gate) == 0L){
      problems <- c(problems, paste0(basename(path), ": no file-level profile gate"))
      next
    }
    file_lanes <- gate_lanes(file_gate[[1L]])
    if(!setequal(file_lanes, registry_lanes)){
      problems <- c(problems, sprintf(
        "%s: gate lanes {%s}, registry lanes {%s}",
        basename(path),
        paste(sort(file_lanes), collapse = ", "),
        paste(sort(registry_lanes), collapse = ", ")
      ))
    }

    # Gates inside the file can only narrow the file gate.
    lines <- readLines(path, warn = FALSE)
    for(line in lines[grepl(gate_line, lines)]){
      inner_lanes <- gate_lanes(parse(text = trimws(line), keep.source = FALSE)[[1L]])
      if(!all(inner_lanes %in% file_lanes)){
        problems <- c(problems, sprintf(
          "%s: inner gate {%s} outside the file gate {%s}",
          basename(path),
          paste(inner_lanes, collapse = ", "),
          paste(file_lanes, collapse = ", ")
        ))
      }
    }
  }

  expect_identical(problems, character(), info = paste(problems, collapse = "\n"))
})


# ---------------------------------------------------------------------------- #
# JAGS build configuration and CI
# ---------------------------------------------------------------------------- #

.layout_repository_file <- function(...){

  testthat::test_path("..", "..", ...)
}

test_that("certification cases are registered, outside `all`, and keyed by the code they reach", {

  case_files <- list.files(testthat::test_path(), pattern = "^test-certification-.*\\.R$")
  cases      <- names(bayestools_certification_cases())
  expect_setequal(sub("^test-certification-(.*)\\.R$", "\\1", case_files), cases)
  expect_setequal(bayestools_test_profile_contexts$certification, paste0("certification-", cases))
  expect_false("certification" %in% bayestools_normalize_test_profiles("all"))

  # the key follows the case file into the helpers and the package objects
  # its code reaches, and not into the rest of the package
  helpers <- environment(bayestools_certification_key)
  key <- bayestools_certification_key("model-averaging", helpers)
  expect_identical(setdiff(c("case_file", "description", "r_version", "jags_version", "helper_fn_expect_true",
                             "package_fn_models_inference", "package_fn_JAGS_bridgesampling", "package_fn_.onLoad"),
                           names(key)), character())
  expect_false("package_fn_hypothesis_BF.BayesTools_fit" %in% names(key))
  expect_true(any(startsWith(names(key), "package_src_")))
  # the same key wherever the helpers are sourced (devtools::test() sources
  # them into the attached package environment)
  fresh <- new.env(parent = asNamespace("BayesTools"))
  for (name in names(key)[startsWith(names(key), "helper_fn_")]) {
    assign(sub("^helper_fn_", "", name), get(sub("^helper_fn_", "", name), envir = helpers), envir = fresh)
  }
  expect_identical(bayestools_certification_key("model-averaging", fresh), key)

  # a case without JAGS that reaches native routines (prior densities) keys
  # the native sources but not the JAGS version
  density_key <- bayestools_certification_key("density-kernels", helpers)
  expect_true(any(startsWith(names(density_key), "package_src_")))
  expect_false("jags_version" %in% names(density_key))
})


test_that("certification cases key their reference baselines and model files", {

  # Reference-table helpers and independently compiled Stan programs both
  # depend on inputs outside the case's R code.
  for (case in names(bayestools_certification_cases())) {
    code <- readLines(bayestools_certification_file(case), warn = FALSE)
    uses_references <- any(grepl("test_reference_|rstan::stan_model\\(", code))
    directories     <- bayestools_certification_cases()[[case]]$references
    expect_identical(length(directories) > 0L, uses_references, info = case)
    expect_true(all(dir.exists(file.path(testthat::test_path(), "..", directories))), info = case)
  }

  # a changed, a removed, and a new baseline change the key; a candidate does not
  results <- withr::local_tempdir()
  tests   <- file.path(results, "tests")
  dir.create(file.path(tests, "testthat"), recursive = TRUE)
  dir.create(file.path(tests, "results", "case"), recursive = TRUE)
  file.copy(bayestools_certification_file("hypothesis-bayes-factors"), file.path(tests, "testthat", "test-certification-case.R"))
  file.copy(testthat::test_path("..", "..", "DESCRIPTION"), results)
  writeLines("reviewed", file.path(tests, "results", "case", "a.txt"))
  cases   <- list(case = list(description = "", native = FALSE, references = "results/case"))
  helpers <- environment(bayestools_certification_key)
  key_of  <- function() {
    bayestools_certification_key("case", helpers, test_dir = file.path(tests, "testthat"), cases = cases)
  }
  key <- key_of()
  expect_true("reference_results/case/a.txt" %in% names(key))
  writeLines("candidate", file.path(tests, "results", "case", "a.new.txt"))
  expect_identical(key_of(), key)
  writeLines("changed", file.path(tests, "results", "case", "a.txt"))
  expect_false(identical(key_of(), key))
  unlink(file.path(tests, "results", "case", "a.txt"))
  expect_identical(unname(key_of()[["reference_results/case"]]), "missing")
})


test_that("a certification marker is current only while it holds the key", {

  marker <- withr::local_tempfile(fileext = ".txt")
  key    <- c(case_file = "a", package_fn_f = "b", helper_fn_g = "c")
  expect_identical(bayestools_certification_state("case", key, marker)$state, "missing")

  bayestools_write_certification_marker("case", key, marker)
  state <- bayestools_certification_state("case", key, marker)
  expect_identical(state$state, "current")
  expect_false(is.na(state$completed_at))

  changed <- bayestools_certification_state("case", c(case_file = "a", package_fn_f = "B", package_fn_h = "d"), marker)
  expect_identical(changed$state, "stale")
  expect_identical(changed$changed, c("helper_fn_g", "package_fn_f", "package_fn_h"))
  expect_error(bayestools_write_certification_marker("case", c(case_file = "a\nb"), marker), "spans lines")
})


test_that("a certification rerun invalidates the case before its process starts", {

  runner_path <- testthat::test_path("..", "..", "tools", "test-certification.R")
  skip_if_not(file.exists(runner_path), "Repository runner source is unavailable in installed-package tests.")
  runner <- new.env()
  source(runner_path, local = runner)

  # a case process that dies at startup, before any bookkeeping of its own
  directory <- withr::local_tempdir()
  crash     <- file.path(directory, "crash.R")
  writeLines("quit(save = 'no', status = 3L)", crash)
  helpers <- list(
    bayestools_certification_marker_file = function(case) file.path(directory, paste0(case, ".txt")),
    bayestools_certification_timeout     = 60
  )
  bayestools_write_certification_marker("case", c(case_file = "a"), helpers$bayestools_certification_marker_file("case"))

  passed <- suppressMessages(runner$.run_certification_case("case", crash, helpers))
  expect_false(passed)
  expect_false(file.exists(helpers$bayestools_certification_marker_file("case")))
})


test_that("native modules use the JAGS headers and libraries selected by rjags", {

  makevars <- .layout_repository_file("src", c("Makevars", "Makevars.win"))
  skip_if_not(all(file.exists(makevars)), "The package build sources are unavailable.")
  for(path in makevars){
    flags <- paste(readLines(path), collapse = "\n")
    for(type in c("CPPFLAGS", "CXXFLAGS", "LIBS")){
      expect_match(flags, sprintf("rjags::pkg.compile.flags('%s')", type), fixed = TRUE)
    }
    expect_match(flags, if(basename(path) == "Makevars.win") "-ljrmath-0" else "-ljrmath",
                 fixed = TRUE)
  }
  compiled <- BayesTools:::.BayesTools_jags_build$JAGS_VERSION
  expect_equal(compiled[, 1L], rjags::jags.version()[, 1L])
  expect_gte(as.integer(compiled[, 1L]), 5L)
})

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

.layout_action_files <- function(){

  action_dir <- .layout_repository_file(".github", "actions")
  skip_if_not(
    dir.exists(action_dir),
    "Repository action sources are not available in this installed-package context."
  )

  list.files(action_dir, pattern = "^action\\.ya?ml$", full.names = TRUE, recursive = TRUE)
}

test_that("every CI workflow installs JAGS through the shared JAGS 5 action", {

  # The action (also used by the dependent packages) is the one place that
  # installs JAGS; a workflow that runs jobs must not install it by itself.
  setup_jags <- "uses: (\\./|FBartos/BayesTools/)\\.github/actions/setup-jags5(@|\\s*$)"
  for(path in .layout_workflow_files()){
    workflow <- readLines(path, warn = FALSE, encoding = "UTF-8")
    if(!any(grepl("^\\s*runs-on:", workflow))){
      next
    }
    expect_true(any(grepl(setup_jags, workflow)), info = basename(path))
    expect_false(any(grepl("jags-installer|apt-get install.*\\bjags\\b|JAGS-4\\.3\\.[0-9]+\\.exe", workflow, perl = TRUE)), info = basename(path))
  }
})

test_that("JAGS 5 CI installs and checks a matching runtime and R interface", {

  path <- .layout_repository_file(".github", "actions", "setup-jags5", "action.yml")
  skip_if_not(file.exists(path), "The CI actions are unavailable.")
  action <- paste(readLines(path), collapse = "\n")
  expect_match(action, "JAGS-5[.][0-9]+[.][0-9]+-universal[.]pkg")
  expect_match(action, 'type = "source"', fixed = TRUE)
  expect_match(action, "rjags::jags.version()", fixed = TRUE)
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

# Whether `expr` calls `name` outside every function definition, i.e. where a
# script evaluates it at top level.
.layout_calls_outside_functions <- function(expr, name){

  if(!is.call(expr)){
    return(FALSE)
  }
  if(identical(expr[[1L]], as.name("function"))){
    return(FALSE)
  }
  if(identical(expr[[1L]], as.name(name))){
    return(TRUE)
  }
  for(i in seq_along(expr)[-1L]){
    if(identical(expr[[i]], quote(expr = ))){
      next
    }
    if(.layout_calls_outside_functions(expr[[i]], name)){
      return(TRUE)
    }
  }

  FALSE
}

test_that("workflow R scripts never rely on top-level on.exit()", {

  # `Rscript {0}` evaluates a run block at top level, where on.exit() never
  # fires, so a cleanup registered that way is silently skipped.
  scripts <- 0L
  for(path in .layout_workflow_files()){
    for(step in .layout_workflow_steps(path)){
      if(!identical(step$shell, "Rscript {0}")){
        next
      }
      scripts <- scripts + 1L
      expressions <- parse(text = step$run, keep.source = FALSE)
      top_level_exit <- vapply(expressions, .layout_calls_outside_functions, logical(1), name = "on.exit")
      expect_false(any(top_level_exit), info = paste0(basename(path), ": ", step$name))
    }
  }

  expect_gt(scripts, 0L)
})

test_that("tools/test.R points the test helpers at the cache directory it creates", {

  script <- .layout_repository_file("tools", "test.R")
  skip_if_not(
    file.exists(script),
    "Repository tools are not available in this installed-package context."
  )

  expressions <- as.list(parse(script, keep.source = FALSE))
  set_env <- Filter(function(expression){
    is.call(expression) && identical(expression[[1L]], as.name("Sys.setenv"))
  }, expressions)
  expect_length(set_env, 1L)

  files_dir <- as.list(set_env[[1L]])[["BAYESTOOLS_TEST_FILES_DIR"]]
  expect_false(is.null(files_dir))
  expect_true("cache_dir" %in% all.names(files_dir))
})

test_that("DESCRIPTION declares the supported JAGS range", {

  description_file <- .layout_repository_file("DESCRIPTION")
  skip_if_not(
    file.exists(description_file),
    "Repository DESCRIPTION is not available in this installed-package context."
  )

  requirements <- unname(read.dcf(description_file, fields = "SystemRequirements")[1, 1])
  expect_match(requirements, "JAGS (>= 5.0.0)", fixed = TRUE)
  fields <- read.dcf(description_file, fields = c("Imports", "Suggests"))
  expect_match(fields[1L, "Imports"], "rjags (>= 5-5)", fixed = TRUE)
  expect_false(any(grepl("\\brunjags\\b", fields)))
})

test_that("expectation problems report every failed check without recording it", {

  # The helper lets a test check many elements and assert once per unit; a
  # failed element must stay visible, with its own message.
  problems <- expectation_problems({
    expect_true(TRUE)
    expect_identical(1L, 2L, info = "first failure")
    expect_equal(1, 1 + 1e-9, tolerance = 1e-12, info = "second failure")
    expect_equal(1, 1 + 1e-9, tolerance = 1e-6, info = "within tolerance")
    expect_error(stop("expected"), "expected")
    expect_true(FALSE, info = "third failure")
    "the value of the code is not returned"
  })

  expect_length(problems, 3L)
  expect_match(problems[[1L]], "first failure", fixed = TRUE)
  expect_match(problems[[2L]], "second failure", fixed = TRUE)
  expect_match(problems[[3L]], "third failure", fixed = TRUE)
  expect_identical(expectation_problems(expect_true(TRUE)), character())

  # only expectations are collected: an error of the code stops it
  expect_error(
    expectation_problems({
      expect_true(FALSE)
      stop("not an expectation")
    }),
    "not an expectation",
    fixed = TRUE
  )
})

test_that("expect_equal_each() applies the criterion of one number to every element", {

  # relative to the expected value, absolute where that is below the tolerance;
  # infinities and NA agree with themselves
  expect_equal_each(
    c(1 + 1e-9, 1000 * (1 + 1e-9), 1e-10, Inf, NA),
    c(1, 1000, 0, Inf, NA),
    tolerance = 1e-8
  )
  # a large expected value does not hide the difference of another element
  expect_equal(c(1 + 3e-8, 1e6), c(1, 1e6 + 1e-3))
  problems <- expectation_problems({
    expect_equal_each(c(1 + 3e-8, 1e6), c(1, 1e6 + 1e-3), info = "one element")
    expect_equal_each(c(1, NA), c(1, 2), info = "an NA")
    expect_equal_each(1:2, 1, info = "lengths")
    expect_equal_each(1e-9, 0, tolerance = 1e-10, info = "absolute")
    expect_equal_each(c(1, 2), c(1, 2), info = "agreeing elements")
  })
  expect_length(problems, 4L)
  expect_match(problems[[1L]], "1 of 2 elements differ", fixed = TRUE)
  expect_match(problems[[1L]], "one element", fixed = TRUE)
  expect_match(problems[[2L]], "[2] NA vs 2", fixed = TRUE)
  expect_match(problems[[3L]], "The numbers of elements differ: 2 and 1.", fixed = TRUE)
  expect_match(problems[[4L]], "absolute", fixed = TRUE)
})

test_that("passing comparisons are recorded at once and failing ones read as testthat's own", {

  # every passing comparison, also within a tolerance, is recorded as a success
  recorded <- 0L
  withCallingHandlers(
    {
      expect_equal(1, 1)
      expect_identical(letters, letters)
      expect_true(TRUE)
      expect_true(c(a = TRUE))
      expect_false(FALSE)
      expect_equal(1, 1 + 1e-9, tolerance = 1e-6)
      expect_equal(1L, 1)
    },
    expectation_success = function(condition){
      recorded <<- recorded + 1L
      invokeRestart("continue_test")
    }
  )
  expect_identical(recorded, 7L)

  # a comparison that does not hold fails with testthat's message, the labels of
  # the call, and the tolerance and info of the call
  problems <- expectation_problems({
    expect_equal(1, 2)
    expect_equal(c(a = 1), c(a = 1.1), tolerance = 1e-3, info = "own info")
    expect_identical(1L, 1)
    expect_true(c(TRUE, FALSE))
    expect_true(NA)
    expect_false(TRUE)
    expect_equal(1, 1 + 1e-9, tolerance = 1e-12)
    expect_equal(1, 1.0001, tolerance = 1e-8, label = "own label", expected.label = "own expected")
  })
  expect_length(problems, 8L)
  expect_match(problems[[1L]], "Expected 1 to equal 2.", fixed = TRUE)
  expect_match(problems[[2L]], "Expected `c(a = 1)` to equal `c(a = 1.1)`.", fixed = TRUE)
  expect_match(problems[[2L]], "own info", fixed = TRUE)
  expect_match(problems[[3L]], "Expected 1L to be identical to 1.", fixed = TRUE)
  expect_match(problems[[4L]], "Expected `c(TRUE, FALSE)` to be TRUE.", fixed = TRUE)
  expect_match(problems[[5L]], "Expected NA to be TRUE.", fixed = TRUE)
  expect_match(problems[[6L]], "Expected TRUE to be FALSE.", fixed = TRUE)
  expect_match(problems[[8L]], "Expected own label to equal own expected.", fixed = TRUE)

  # a call longer than a line is labelled as testthat labels it
  long_call <- expectation_problems(
    expect_equal(sum(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)), 0)
  )
  expect_identical(long_call, expectation_problems(
    testthat::expect_equal(sum(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)), 0)
  ))
  expect_match(long_call, "Expected `sum(...)` to equal 0.", fixed = TRUE)

  # testthat decides, with its own outcome, what identical() alone does not
  # settle: a comparison option (here one stricter than identical()), an object
  # of expect_true() / expect_false() that is not a logical vector, a tolerance
  # it refuses, and an option it warns about
  latin1 <- iconv("\u00e9", "UTF-8", "latin1")
  masked_problems <- expectation_problems({
    expect_equal("\u00e9", latin1, ignore_encoding = FALSE)
    expect_identical("\u00e9", latin1, ignore_encoding = FALSE)
    expect_true(emptyenv())
    expect_false(emptyenv())
  })
  testthat_problems <- expectation_problems({
    testthat::expect_equal("\u00e9", latin1, ignore_encoding = FALSE)
    testthat::expect_identical("\u00e9", latin1, ignore_encoding = FALSE)
    testthat::expect_true(emptyenv())
    testthat::expect_false(emptyenv())
  })
  expect_length(testthat_problems, 4L)
  expect_identical(masked_problems, testthat_problems)
  expect_error(expect_equal(1, 1, tolerance = -1), "tolerance")
  expect_warning(expect_equal(1, 1, tol = 1e-3), "deprecated")
})
