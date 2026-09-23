main <- function() {
  if (!file.exists("DESCRIPTION") || !dir.exists("vignettes")) {
    stop(
      "Run this script from the BayesTools project root.",
      call. = FALSE
    )
  }

  source_file <- normalizePath(
    file.path("vignettes", "RandomEffects.Rmd"),
    winslash = "/",
    mustWork = TRUE
  )
  temporary_file <- file.path(
    dirname(source_file),
    "RandomEffects-regenerate.Rmd"
  )
  if (file.exists(temporary_file)) {
    stop(
      "The temporary RandomEffects regeneration vignette already exists.",
      call. = FALSE
    )
  }
  on.exit(
    {
      if (file.exists(temporary_file)) {
        unlink(temporary_file)
      }
    },
    add = TRUE
  )

  lines <- readLines(source_file, warn = FALSE, encoding = "UTF-8")
  can_evaluate_line <- grepl("can_evaluate <- ", lines, fixed = TRUE)
  if (sum(can_evaluate_line) != 1L) {
    stop(
      "Could not identify the RandomEffects evaluation gate.",
      call. = FALSE
    )
  }
  lines[can_evaluate_line] <- "can_evaluate <- TRUE"

  cache_guards <- grepl(
    "^if\\((!file\\.exists\\(random_effects_cache_file\\)|length\\(missing_packages\\) == 0L)",
    lines
  )
  if (sum(cache_guards) != 2L) {
    stop(
      "Could not identify both RandomEffects cache validity guards.",
      call. = FALSE
    )
  }
  lines[cache_guards] <- sub(
    "if(",
    "if(FALSE && ",
    lines[cache_guards],
    fixed = TRUE
  )

  load_header <- grepl(
    "load-precomputed-random-effects,",
    lines,
    fixed = TRUE
  )
  if (sum(load_header) != 1L) {
    stop(
      "Could not identify the RandomEffects cache-loading chunk.",
      call. = FALSE
    )
  }
  lines[load_header] <- sub(
    "purl = FALSE}",
    "eval = FALSE, purl = FALSE}",
    lines[load_header],
    fixed = TRUE
  )

  regeneration_labels <- c(
    "begin-random-effects-cache-regeneration,",
    "fit-stan-correlated-regenerate,",
    "fit-sleep-regenerate,",
    "fit-sleep-independent-regenerate,",
    "fit-cake-recipe-regenerate,",
    "fit-cake-hcs-regenerate,",
    "fit-cake-ar1-regenerate,",
    "fit-cake-har-regenerate,",
    "fit-sleep-car-regenerate,",
    "fit-nested-regenerate,",
    "fit-stan-crossed-regenerate,",
    "fit-crossed-independent-regenerate,",
    "fit-crossed-allocation-regenerate,",
    "save-precomputed-random-effects,"
  )
  regeneration_headers <- Reduce(
    "|",
    lapply(regeneration_labels, function(label) {
      grepl(label, lines, fixed = TRUE)
    })
  )
  if (sum(regeneration_headers) != length(regeneration_labels)) {
    stop(
      "Expected ",
      length(regeneration_labels),
      " RandomEffects regeneration chunks, found ",
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

  reload_line <- grepl(
    'pkgload::load_all("..", quiet = TRUE)',
    lines,
    fixed = TRUE
  )
  if (sum(reload_line) != 1L) {
    stop(
      "Could not identify the RandomEffects source-reload call.",
      call. = FALSE
    )
  }
  lines[reload_line] <- "# The regeneration driver loaded the current source before knitting."

  writeLines(lines, temporary_file, useBytes = TRUE)
  pkgload::load_all(".", quiet = TRUE)
  message("Rendering the complete RandomEffects regeneration document.")
  rmarkdown::render(
    temporary_file,
    output_dir = tempdir(),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  source(file.path("vignettes", "random-effects-vignette-cache.R"))
  stop_if_invalid_random_effects_vignette_cache()
  message(
    "Validated RandomEffects cache with ",
    length(random_effects_vignette_cache_names()),
    " models."
  )
}

main()
