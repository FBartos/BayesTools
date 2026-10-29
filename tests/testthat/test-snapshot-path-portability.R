skip_if_not_test_profile("unit")


test_that("the package build ships only the CRAN smoke test", {

  package_root <- testthat::test_path("..", "..")
  build_ignore_path <- file.path(package_root, ".Rbuildignore")
  skip_if_not(
    file.exists(build_ignore_path),
    "Repository build-ignore policy is unavailable in this installed-package test context."
  )
  build_ignore <- readLines(build_ignore_path, warn = FALSE)
  test_paths   <- file.path("tests", list.files(
    file.path(package_root, "tests"),
    recursive  = TRUE,
    all.files  = TRUE,
    include.dirs = TRUE
  ))
  # R CMD build matches .Rbuildignore as case-insensitive Perl regexes and
  # drops everything below an excluded directory
  matches_build_ignore <- function(path) {
    any(vapply(build_ignore, grepl, logical(1), x = path,
               ignore.case = TRUE, perl = TRUE))
  }
  excluded <- vapply(test_paths, function(path) {
    parts     <- strsplit(path, "/", fixed = TRUE)[[1L]]
    ancestors <- vapply(seq_along(parts), function(i) {
      paste(parts[seq_len(i)], collapse = "/")
    }, character(1))
    any(vapply(ancestors, matches_build_ignore, logical(1)))
  }, logical(1))
  expect_setequal(test_paths[!excluded], c(
    "tests/namespace-only-contrasts.R",
    "tests/testthat.R",
    "tests/testthat",
    "tests/testthat/test-cran-smoke.R",
    "tests/testthat/helper-test-profiles.R"
  ))
})


test_that("visual snapshot paths are unique on case-insensitive file systems", {

  snapshot_paths <- list.files(
    testthat::test_path("_snaps"),
    pattern = "\\.svg$",
    recursive = TRUE,
    full.names = FALSE
  )
  package_paths <- gsub("\\", "/", snapshot_paths, fixed = TRUE)

  collision_keys <- tolower(package_paths)
  collision_paths <- sort(unique(package_paths[
    duplicated(collision_keys) |
      duplicated(collision_keys, fromLast = TRUE)
  ]))
  expect_identical(collision_paths, character())
})


test_that("visual snapshot labels have collision-free portable names", {

  snapshot_slug <- function(label) {
    slug <- gsub("[^a-z0-9]+", "-", tolower(label))
    gsub("(^-+|-+$)", "", slug)
  }

  collect_doppelganger_labels <- function(node) {
    if(is.expression(node) || is.pairlist(node)){
      return(unlist(
        lapply(as.list(node), collect_doppelganger_labels),
        use.names = FALSE
      ))
    }
    if(!is.call(node)){
      return(character())
    }

    call_head <- node[[1L]]
    call_name <- if(is.symbol(call_head)){
      as.character(call_head)
    }else if(
      is.call(call_head) &&
      identical(as.character(call_head[[1L]]), "::")
    ){
      as.character(call_head[[3L]])
    }else{
      ""
    }

    labels <- character()
    if(
      identical(call_name, "expect_doppelganger") &&
      length(node) >= 2L &&
      is.character(node[[2L]]) &&
      length(node[[2L]]) == 1L
    ){
      labels <- node[[2L]]
    }

    child_nodes <- if(length(node) > 1L){
      as.list(node)[-1L]
    }else{
      list()
    }
    c(
      labels,
      unlist(
        lapply(child_nodes, collect_doppelganger_labels),
        use.names = FALSE
      )
    )
  }

  snapshot_root <- testthat::test_path("_snaps")
  snapshot_contexts <- basename(list.dirs(
    snapshot_root,
    recursive = FALSE,
    full.names = TRUE
  ))
  test_files <- testthat::test_path(
    paste0("test-", snapshot_contexts, ".R")
  )

  missing_test_files <- basename(test_files[!file.exists(test_files)])
  expect_identical(missing_test_files, character())

  collisions <- unlist(lapply(seq_along(test_files), function(i){
    labels <- collect_doppelganger_labels(parse(
      file = test_files[[i]],
      keep.source = FALSE
    ))
    slugs <- snapshot_slug(labels)
    duplicated_slugs <- unique(slugs[
      duplicated(slugs) |
        duplicated(slugs, fromLast = TRUE)
    ])
    if(length(duplicated_slugs) == 0L){
      return(character())
    }
    paste0(snapshot_contexts[[i]], ": ", duplicated_slugs)
  }), use.names = FALSE)

  expect_identical(collisions, character())
})
