skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Static Fixture Catalog
# ============================================================================ #
#
# PURPOSE:
#   Validates source-derived fixture catalog accountability without requiring
#   generated JAGS fixture artifacts.
#
# TAGS: @fixture, @catalog, @unit
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("source-derived fixture catalog covers every generated fit", {
  source_rows <- .bayestools_save_fit_catalog_rows()
  catalog <- bayestools_expected_fit_catalog()
  expect_expected_fit_catalog_schema(catalog)

  expect_equal(nrow(source_rows), 105L)
  expect_equal(sum(source_rows$has_marglik), 22L)
  expect_equal(sum(source_rows$assertion_only), 32L)
  expect_equal(nrow(catalog), nrow(source_rows))
  expect_equal(catalog$model_name, source_rows$model_name)
  expect_equal(catalog$fit_file, paste0(catalog$model_name, ".RDS"))
  expect_equal(catalog$has_marglik, source_rows$has_marglik)
  expect_equal(catalog$assertion_only, source_rows$assertion_only)
  expect_equal(catalog$note, source_rows$note)
  expect_false(anyDuplicated(catalog$model_name) > 0L)

  optional_catalog <- bayestools_optional_fit_requirements()
  expect_setequal(
    catalog$model_name[!is.na(catalog$requires_package)],
    optional_catalog$model_name
  )
  required_catalog <- bayestools_required_fit_catalog(catalog)
  unavailable_optional <- optional_catalog$model_name[
    !vapply(optional_catalog$requires_package, requireNamespace, logical(1), quietly = TRUE)
  ]
  expect_setequal(
    setdiff(catalog$model_name, required_catalog$model_name),
    unavailable_optional
  )

  semantic_catalog <- bayestools_semantic_fit_catalog_overrides()
  expect_setequal(
    catalog$model_name[catalog$oracle_type != "registry-metadata"],
    semantic_catalog$model_name
  )
  expect_equal(catalog$model_name[catalog$oracle_type != "registry-metadata"], catalog$model_name)
  expect_true(all(vapply(catalog$expected_monitor, length, integer(1)) > 0L))
})

test_that("source-derived fixture catalog preserves registry schema flags", {
  catalog <- bayestools_expected_fit_catalog()
  flag_cols <- c("has_marglik", bayestools_registry_flag_columns())

  for (flag in flag_cols) {
    expect_type(catalog[[flag]], "logical")
    expect_false(anyNA(catalog[[flag]]), info = paste("Catalog flag:", flag))
  }

  expect_true(all(nzchar(catalog$model_name)))
  expect_true(all(nzchar(catalog$note)))
  expect_type(catalog$requires_package, "character")
  expect_equal(is.na(catalog$marglik_file), !catalog$has_marglik)
  expect_equal(
    catalog$marglik_file[catalog$has_marglik],
    paste0(catalog$model_name[catalog$has_marglik], ".RDS")
  )
})

test_that("reviewed fits have summary-table baselines and assertion-only fits have none", {
  reference_dir <- testthat::test_path("..", "results", "JAGS-summary-tables")
  skip_if_not(dir.exists(reference_dir), "Reference tables are not available in this installed-package test context.")

  catalog <- bayestools_expected_fit_catalog()
  baselines <- file.exists(file.path(
    reference_dir,
    paste0(catalog$model_name, "_runjags_estimates.txt")
  ))
  expect_identical(catalog$model_name[!catalog$assertion_only & !baselines], character())
  expect_identical(catalog$model_name[catalog$assertion_only & baselines], character())
})

test_that("legacy cache markers are treated as stale metadata", {
  marker_name <- paste0("legacy-marker-", Sys.getpid())
  marker_file <- .test_cache_indicator_file(marker_name)
  withr::defer(unlink(marker_file), testthat::teardown_env())

  writeLines(
    c(
      paste("name:", marker_name),
      "completed_at: 2026-05-03 23:05:40 CEST",
      paste("test_files_dir:", test_files_dir)
    ),
    marker_file
  )

  expect_false(
    .test_cache_metadata_current(
      marker_name,
      required_fits = c("fit_a", "fit_b"),
      required_margliks = "fit_a"
    )
  )
})

test_that("fixture cache availability policy is profile-aware", {
  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "unit")
  expect_false(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "fixture")
  expect_true(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "fit")
  expect_true(.test_cache_required_for_active_profile())

  withr::local_envvar(BAYESTOOLS_TEST_PROFILE = "visual-fixture")
  expect_true(.test_cache_required_for_active_profile())
})

test_that("model-fit cache marker hashes only fit-generation sources", {
  skip_if_not(.test_cache_package_sources_available(), "Repository R source files are not available in this installed-package test context.")

  source_files <- .test_cache_source_files("model-fit")

  # package R code enters the key as the functions the fit generators reach and
  # the DESCRIPTION without its version as a fingerprint, not as whole files
  expect_false(any(startsWith(names(source_files), "package_R_")))
  expect_false("description" %in% names(source_files))
  expect_true("package_src_r_lkj_cc" %in% names(source_files))
  expect_true("package_src_lkj_BTLKJCore_cc" %in% names(source_files))
  expect_true("package_src_functions_BTLKJCholesky_cc" %in% names(source_files))
  expect_true("package_src_distributions_DBTLKJCPC_cc" %in% names(source_files))
  # every native source and build description under src/ (build products and
  # a Makevars generated by configure are not sources)
  src_dir <- testthat::test_path("..", "..", "src")
  native_sources <- list.files(
    src_dir,
    pattern = "\\.(c|cc|cpp|h|hpp)$|^Makevars\\.",
    recursive = TRUE,
    full.names = TRUE
  )
  expect_true(length(native_sources) > 0L)
  expect_setequal(
    normalizePath(unname(source_files[startsWith(names(source_files), "package_src_")]), winslash = "/"),
    normalizePath(native_sources, winslash = "/")
  )
  expect_false("common_functions" %in% names(source_files))
  expect_false("expected_fit_catalog" %in% names(source_files))
  expect_true("test_00_model_fits" %in% names(source_files))
  expect_true(all(file.exists(source_files)), info = paste(source_files[!file.exists(source_files)], collapse = ", "))

  source_functions <- .test_cache_source_functions("model-fit")
  expect_true("test_helper_save_fit" %in% names(source_functions))
  expect_true("catalog_expected_fit" %in% names(source_functions))
  expect_true("catalog_semantic_fit_overrides" %in% names(source_functions))
  # the catalog helpers that build the declared rows (their expected monitors
  # enter the fixture metadata of every saved fit) are hashed with them
  helper_env <- environment(bayestools_semantic_fit_catalog_overrides)
  row_helpers <- Filter(
    function(name) exists(name, envir = helper_env, mode = "function", inherits = FALSE),
    unique(all.names(body(bayestools_semantic_fit_catalog_overrides)))
  )
  expect_true(".bayestools_lkj_monitors" %in% row_helpers)
  expect_identical(setdiff(row_helpers, source_functions), character())

  source_hashes <- .test_cache_source_hashes("model-fit")
  is_package_object <- grepl("^package_(fn|obj)_", names(source_hashes))
  expect_true(any(is_package_object))
  expect_setequal(
    names(source_hashes)[!is_package_object],
    c("description", names(source_files), names(source_functions))
  )
  expect_false(anyNA(source_hashes))
})

test_that("the model-fit key hashes the package functions that the fit generators reach", {
  package_hashes <- .test_cache_package_function_hashes("model-fit")
  reached <- sub("^package_(fn|obj)_", "", names(package_hashes))

  # the session memo of the hashes is not an option that JAGS_runtime_cluster()
  # forwards to its workers
  expect_false(any(startsWith(names(options()), "BayesTools.test_cache")))
  expect_false(anyNA(package_hashes))
  expect_false(anyDuplicated(names(package_hashes)) > 0L)
  # fitting, marginal likelihoods, prior constructors and their computed
  # dispatch, the parameter map and catalog that a fit stores, and a
  # namespace constant that they read
  generator_code <- c(
    "JAGS_fit", "JAGS_bridgesampling", "prior", "prior_random", ".prior_normal",
    ".prior_mnormal", "rng.prior", "lpdf.prior", ".bt_attach_parameter_map",
    ".bt_build_parameter_catalog", ".bt_fit_contract_version", ".onLoad"
  )
  expect_identical(setdiff(generator_code, reached), character())
  expect_true(startsWith(
    names(package_hashes)[reached == ".bt_fit_contract_version"],
    "package_obj_"
  ))
  # inference, plotting, summaries, and the exports that act on fitted objects
  post_fit_code <- c(
    "marginal_inference", "mix_posteriors", "ensemble_inference",
    "plot_posterior", "plot_prior_list", "runjags_estimates_table",
    "JAGS_estimates_table", "interpret", "hypothesis_BF", "inclusion_BF",
    "parameter_prior_density", "JAGS_extend"
  )
  namespace_names <- ls(asNamespace("BayesTools"), all.names = TRUE)
  expect_identical(setdiff(post_fit_code, namespace_names), character())
  expect_identical(intersect(post_fit_code, reached), character())
  # environments of the namespace are state of a session, not code
  environments <- namespace_names[vapply(
    namespace_names,
    function(name) is.environment(get(name, envir = asNamespace("BayesTools"))),
    logical(1)
  )]
  expect_identical(intersect(environments, reached), character())
})

test_that("the reachable package objects follow calls, names, dispatch, and computed names", {
  toy <- new.env(parent = baseenv())
  local({
    entry <- function(x, scale = default_scale) helper(x) * scale
    helper <- function(y) deeper(y)
    deeper <- function(z) z + offset
    offset <- 2L
    default_scale <- 3L
    by_string <- function() do.call("string_target", list())
    string_target <- function() 1
    passed <- function() lapply(1:2, passed_function)
    passed_function <- function(i) i
    field_reader <- function(x) x$field_name
    field_name <- function() 1
    tag_writer <- function() list(tag_name = 1)
    tag_name <- function() 1
    local_shadow <- function() lapply(1:2, function(local_name) local_name)
    local_name <- function() 1
    replacement <- function(x) {
      names(x) <- "a"
      x
    }
    `names<-` <- function(x, value) x
    state <- new.env()
    reads_state <- function() state
    unreached <- function() entry(1)
    generic <- function(x) UseMethod("generic")
    generic.alpha <- function(x) alpha_helper()
    generic.default <- function(x) NULL
    alpha_helper <- function() 1
    not_a_method <- function() 1
    builder <- function(kind) {
      if (!kind %in% c("a", "b")) stop("unknown")
      do.call(paste0(".make_", kind), list())
    }
    .make_a <- function() make_shared()
    .make_b <- function() 1
    .make_c <- function() 1
    make_shared <- function() 1
    unmatched <- function(kind) get(paste0(".other_", kind))
    .other_x <- function() 1
    .other_y <- function() 1
    call_density <- function(x) density(x)
    call_other <- function(x) other_generic(x)
    names_a_class <- function(x) inherits(x, "named_class")
    unnamed_class_user <- function(x) density(x)
    density.named_class <- function(x) named_method_helper()
    density.unnamed_class <- function(x) unnamed_method_helper()
    density.default <- function(x) default_method_helper()
    other_generic.named_class <- function(x) 1
    named_method_helper <- function() 1
    unnamed_method_helper <- function() 1
    default_method_helper <- function() 1
    Ops.group_class <- function(e1, e2) group_helper()
    group_helper <- function() 1
    group_class_maker <- function() structure(1, class = "group_class")
  }, envir = toy)
  registry <- rbind(
    c("density", "named_class", "density.named_class"),
    c("density", "unnamed_class", "density.unnamed_class"),
    c("density", "default", "density.default"),
    c("other_generic", "named_class", "other_generic.named_class"),
    c("Ops", "group_class", "Ops.group_class")
  )
  reach <- function(roots, root_strings = character()) {
    .test_cache_reached_package_objects(roots, root_strings, namespace = toy, s3_registry = registry)
  }

  # calls, argument defaults, and namespace constants; not what nothing reaches
  expect_identical(reach("entry"), c("deeper", "default_scale", "entry", "helper", "offset"))
  # a string constant names a function; a function can be passed by name
  expect_identical(reach("by_string"), c("by_string", "string_target"))
  expect_identical(reach("passed"), c("passed", "passed_function"))
  # the field of $, the tag of an argument, and a local variable are not reads
  expect_identical(reach("field_reader"), "field_reader")
  expect_identical(reach("tag_writer"), "tag_writer")
  expect_identical(reach("local_shadow"), "local_shadow")
  # a replacement function is reached by the assignment that calls it
  expect_identical(reach("replacement"), c("names<-", "replacement"))
  # environments are not code
  expect_identical(reach("reads_state"), "reads_state")
  # methods of a generic that a reached function dispatches with, but not
  # functions that only share a prefix with a name
  expect_identical(
    reach("generic"),
    c("alpha_helper", "generic", "generic.alpha", "generic.default")
  )
  # a computed function name reaches the completions by the strings of its
  # function, and all functions of the prefix when no completion exists
  expect_identical(reach("builder"), c(".make_a", ".make_b", "builder", "make_shared"))
  expect_identical(reach("unmatched"), c(".other_x", ".other_y", "unmatched"))
  # a registered method needs a generic that reached code uses and a class
  # that it names, or the default
  expect_identical(
    reach("call_density"),
    c("call_density", "default_method_helper", "density.default")
  )
  expect_identical(
    reach(c("call_density", "names_a_class")),
    c(
      "call_density", "default_method_helper", "density.default",
      "density.named_class", "named_method_helper", "names_a_class"
    )
  )
  expect_identical(
    reach("call_other"),
    "call_other"
  )
  expect_identical(
    reach("call_other", root_strings = "named_class"),
    c("call_other", "other_generic.named_class")
  )
  # a group generic is used by any arithmetic: only the class is required
  expect_identical(
    reach("group_class_maker"),
    c("Ops.group_class", "group_class_maker", "group_helper")
  )
})

test_that("code references separate calls, reads, locals, strings, and computed names", {
  references <- .test_cache_code_references(function(x, y = default_value) {
    value <- called_function(x, "a string")
    for (item in seq_along(x)) {
      value <- value + read_value
    }
    lapply(x, function(inner) inner + helper_reference)
    x$field
    UseMethod("a_generic")
    do.call(paste0("prefix_", "kind_"), list(y))
    BayesTools:::qualified_call(x)
  })

  expect_true(all(c(
    "called_function", "read_value", "helper_reference", "default_value",
    "qualified_call", "lapply", "UseMethod", "do.call"
  ) %in% references$used))
  # locals are not reads
  expect_false(any(c("x", "y", "value", "item", "inner", "field") %in% references$used))
  expect_true(all(c("a string", "a_generic", "prefix_", "kind_") %in% references$strings))
  expect_identical(references$generics, "a_generic")
  expect_setequal(references$prefixes, c("prefix_", "kind_"))
})

test_that("the hash of a package object ignores comments and layout but not the digits of a number", {
  commented <- eval(parse(
    text = "function(x) {\n  # a comment\n  x +   1.8378770664093453\n}",
    keep.source = TRUE
  ))
  plain <- eval(parse(text = "function(x) {x + 1.8378770664093453}", keep.source = FALSE))
  # the two literals differ only beyond the 15th significant digit
  last_digit <- eval(parse(text = "function(x) {x + 1.8378770664093455}", keep.source = FALSE))

  expect_identical(.test_cache_object_hash(commented), .test_cache_object_hash(plain))
  expect_false(identical(.test_cache_object_hash(plain), .test_cache_object_hash(last_digit)))
})

test_that("the DESCRIPTION fingerprint of the model-fit key ignores the version and build fields", {
  description <- tempfile("DESCRIPTION")
  withr::defer(unlink(description))
  write_description <- function(lines) writeLines(lines, description)

  write_description(c("Package: Toy", "Version: 0.0.1", "Title: A toy", "Imports: a, b"))
  baseline <- .test_cache_description_md5(description)
  write_description(c(
    "Package: Toy",
    "Version: 9.9.9.9",
    "Packaged: 2026-10-01 10:00:00 UTC; someone",
    "Built: R 4.6.0; ; 2026-10-01 10:00:00 UTC; windows",
    "Author: Someone",
    "Title: A toy",
    "Imports: a,",
    "    b"
  ))
  expect_identical(.test_cache_description_md5(description), baseline)
  write_description(c("Package: Toy", "Version: 0.0.1", "Title: A toy", "Imports: a, b, c"))
  expect_false(identical(.test_cache_description_md5(description), baseline))

  repository_description <- testthat::test_path("..", "..", "DESCRIPTION")
  skip_if_not(file.exists(repository_description), "The repository DESCRIPTION is not available in this installed-package test context.")
  bumped <- sub(
    "^Version: .*$", "Version: 99.0.0.1",
    readLines(repository_description, warn = FALSE)
  )
  write_description(bumped)
  expect_identical(
    .test_cache_description_md5(description),
    .test_cache_description_md5(repository_description)
  )
  expect_identical(
    unname(.test_cache_description_hashes("model-fit")[["description"]]),
    .test_cache_description_md5(repository_description)
  )
  expect_length(.test_cache_description_hashes("fixture-consumer"), 0L)
})

test_that("a cache marker is stale when a source left the key", {
  marker_name <- paste0("source-key-marker-", Sys.getpid())
  marker_file <- .test_cache_indicator_file(marker_name)
  withr::defer(unlink(marker_file), testthat::teardown_env())
  write_marker <- function(...) {
    writeLines(
      c(
        paste("name:", marker_name),
        paste("test_files_dir:", test_files_dir),
        "required_fits: ",
        "required_margliks: ",
        ...
      ),
      marker_file
    )
  }

  # a name without sources has an empty key
  write_marker()
  expect_true(.test_cache_metadata_current(marker_name))
  write_marker("source_md5_package_fn_removed: 0123456789abcdef0123456789abcdef")
  expect_false(.test_cache_metadata_current(marker_name))
})

test_that("downstream cache consumer source scopes are separate from model fitting", {
  skip_if_not(.test_cache_package_sources_available(), "Repository R source files are not available in this installed-package test context.")

  model_fit_sources <- .test_cache_source_files("model-fit")
  fixture_sources <- .test_cache_source_files("fixture-consumer")
  visual_fixture_sources <- .test_cache_source_files("visual-fixture-consumer")

  expect_true("common_functions" %in% names(fixture_sources))
  expect_true("expected_fit_catalog" %in% names(fixture_sources))
  expect_true("package_R_JAGS-formula-random" %in% names(fixture_sources))
  expect_true("package_R_random-effects-metadata" %in% names(fixture_sources))
  expect_true("package_R_random-effects-reconstruction" %in% names(fixture_sources))
  expect_true("package_R_random-effects-summary" %in% names(fixture_sources))
  expect_true("package_R_random-priors" %in% names(fixture_sources))
  expect_true("package_R_summary-tables-model" %in% names(fixture_sources))
  expect_true("package_R_model-averaging" %in% names(fixture_sources))
  expect_false("test_00_model_fits" %in% names(fixture_sources))

  expect_true("common_functions" %in% names(visual_fixture_sources))
  expect_true("expected_fit_catalog" %in% names(visual_fixture_sources))
  expect_true("package_R_JAGS-formula-random" %in% names(visual_fixture_sources))
  expect_true("package_R_random-effects-summary" %in% names(visual_fixture_sources))
  expect_true("package_R_model-averaging-plots-posterior" %in% names(visual_fixture_sources))
  expect_true("test_test_JAGS_diagnostic_plots" %in% names(visual_fixture_sources))
  expect_false("test_00_model_fits" %in% names(visual_fixture_sources))

  expect_false("common_functions" %in% names(model_fit_sources))
  expect_false("expected_fit_catalog" %in% names(model_fit_sources))
  expect_false("package_R_summary-tables-model" %in% names(model_fit_sources))
  expect_false("package_R_model-averaging-plots-posterior" %in% names(model_fit_sources))

  expect_false(anyNA(.test_cache_source_hashes("fixture-consumer")))
  expect_false(anyNA(.test_cache_source_hashes("visual-fixture-consumer")))
})
