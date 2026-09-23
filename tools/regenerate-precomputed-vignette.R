main <- function() {
  if (!file.exists("DESCRIPTION") || !dir.exists("vignettes")) {
    stop(
      "Run this script from the BayesTools project root.",
      call. = FALSE
    )
  }

  args <- commandArgs(trailingOnly = TRUE)
  supported <- c("ComparisonR", "SpikeAndSlab")
  if (length(args) != 1L || !args[[1L]] %in% supported) {
    stop(
      "Supply exactly one vignette name: ",
      paste(supported, collapse = " or "),
      ".",
      call. = FALSE
    )
  }
  vignette <- args[[1L]]
  source_file <- file.path("vignettes", paste0(vignette, ".Rmd"))
  regeneration_rmd <- tempfile(
    paste0("BayesTools-regenerate-", vignette, "-"),
    fileext = ".Rmd"
  )
  regeneration_script <- tempfile(
    paste0("BayesTools-regenerate-", vignette, "-"),
    fileext = ".R"
  )
  on.exit(unlink(c(regeneration_rmd, regeneration_script)), add = TRUE)

  lines <- readLines(source_file, warn = FALSE, encoding = "UTF-8")
  regeneration_headers <- grepl("purl = TRUE", lines, fixed = TRUE) &
    grepl("eval = FALSE", lines, fixed = TRUE)
  expected_headers <- c(ComparisonR = 4L, SpikeAndSlab = 5L)
  if (sum(regeneration_headers) != expected_headers[[vignette]]) {
    stop(
      "Expected ",
      expected_headers[[vignette]],
      " regeneration chunks in ",
      vignette,
      ", found ",
      sum(regeneration_headers),
      ".",
      call. = FALSE
    )
  }
  lines[regeneration_headers] <- sub(
    "eval = FALSE",
    "eval = TRUE",
    lines[regeneration_headers],
    fixed = TRUE
  )
  writeLines(lines, regeneration_rmd, useBytes = TRUE)

  knitr::purl(
    regeneration_rmd,
    output = regeneration_script,
    documentation = 0L,
    quiet = TRUE
  )
  message("Regenerating the ", vignette, " vignette cache.")
  plot_file <- tempfile(
    paste0("BayesTools-regenerate-", vignette, "-"),
    fileext = ".pdf"
  )
  on.exit(unlink(plot_file), add = TRUE)
  grDevices::pdf(plot_file)
  plot_device <- grDevices::dev.cur()
  on.exit(
    {
      if (plot_device %in% grDevices::dev.list()) {
        grDevices::dev.off(plot_device)
      }
    },
    add = TRUE
  )
  withr::with_dir(
    "vignettes",
    sys.source(
      regeneration_script,
      envir = new.env(parent = globalenv())
    )
  )
  if (plot_device %in% grDevices::dev.list()) {
    grDevices::dev.off(plot_device)
  }

  source(file.path("vignettes", "precomputed-vignette-cache.R"))
  objects <- load_precomputed_vignette_cache(vignette)
  message(
    "Validated ", vignette, " cache with ",
    length(objects), " objects."
  )
}

main()
