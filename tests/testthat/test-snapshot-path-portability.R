skip_if_not_test_profile("unit")


test_that("visual snapshot paths fit the portable tar name field", {

  snapshot_root <- testthat::test_path("_snaps")
  snapshot_paths <- list.files(
    snapshot_root,
    pattern = "\\.svg$",
    recursive = TRUE,
    full.names = FALSE
  )
  package_paths <- file.path(
    "BayesTools",
    "tests",
    "testthat",
    "_snaps",
    snapshot_paths
  )
  package_paths <- gsub("\\", "/", package_paths, fixed = TRUE)
  path_bytes <- nchar(enc2utf8(package_paths), type = "bytes")

  violations <- sprintf(
    "%d bytes: %s",
    path_bytes[path_bytes >= 100L],
    package_paths[path_bytes >= 100L]
  )
  expect_identical(violations, character())

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
