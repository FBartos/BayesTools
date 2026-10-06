skip_if_not_test_profile("unit")

test_that("structured JAGS parameter encoding is injective and reversible", {

  cases <- list(
    list(kind = "fixed", formula_parameter = "mu", term = "a_b:c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu_a", term = "b c", role = "coefficient"),
    list(kind = "fixed", formula_parameter = "mu", term = "a__xXx__b", role = "coefficient"),
    list(kind = "random", formula_parameter = "mu", term = "group[one]", role = "sd"),
    list(kind = "fixed", formula_parameter = "mu", term = intToUtf8(c(946L, 233L)), role = "coefficient")
  )

  encoded <- vapply(cases, .bt_parameter_encode, character(1))
  expect_length(unique(encoded), length(cases))
  expect_true(all(grepl("^[A-Za-z][A-Za-z0-9_.]*$", encoded)))
  for(i in seq_along(cases)){
    decoded <- .bt_parameter_decode(encoded[i])
    expect_identical(decoded[names(cases[[i]])], cases[[i]])
    expect_identical(decoded$encoding_version, 1L)
  }

  left <- .bt_parameter_encode(list(
    kind = "fixed", formula_parameter = "a_b", term = "c", role = "coefficient"
  ))
  right <- .bt_parameter_encode(list(
    kind = "fixed", formula_parameter = "a", term = "b_c", role = "coefficient"
  ))
  expect_false(identical(left, right))
  expect_error(.bt_parameter_decode("BT2_00_00_00_00"), "unsupported")
  expect_error(.bt_parameter_decode("BT1_0_00_00_00"), "malformed")
  for(invalid_utf8 in c("FF", "C0AF", "EDA080", "C3")){
    expect_error(
      .bt_parameter_decode(paste0("BT1_", invalid_utf8, "_6D75__636F6566")),
      "invalid UTF-8",
      fixed = TRUE
    )
  }
})

test_that("formula designs persist a validated semantic name map", {

  result <- JAGS_formula(
    formula = ~ x + f,
    parameter = "mu",
    data = data.frame(
      x = c(-1, 0, 1, 2),
      f = factor(c("a", "b", "a", "b"))
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  map <- result$formula_design$name_map

  expect_identical(result$formula_design$schema_version, 5L)
  expect_s3_class(map, "BayesTools_formula_name_map")
  expect_identical(attr(map, "schema_version"), 1L)
  expect_setequal(map$jags_name, c("mu", "mu_intercept", "mu_x", "mu_f"))
  expect_identical(map$term[map$jags_name == "mu_f"], "f")

  auxiliary <- JAGS_formula(
    formula = ~ x,
    parameter = "mu",
    data = data.frame(x = c(-1, 0, 1, 2)),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior_spike_and_slab(prior("normal", list(0, 1)))
    )
  )$formula_design$name_map
  indicator <- auxiliary[auxiliary$jags_name == "mu_x_indicator", ]
  expect_identical(indicator$kind, "fixed_auxiliary")
  expect_identical(indicator$term, "x")
  expect_identical(indicator$role, "_indicator")

  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit <- .bt_attach_fit_contract(fit)
  expect_identical(JAGS_formula_name_map(fit, "mu"), map)
  expect_identical(
    JAGS_fit_contract(fit)$formula_design_version,
    5L
  )
  expect_silent(JAGS_validate_fit_contract(
    fit,
    requires = c("name_encoding", "formula_name_map", "formula_design")
  ))

  broken <- fit
  attr(broken, "fit_contract")$formula_name_map_version <- 2L
  expect_error(
    JAGS_formula_name_map(broken, "mu"),
    paste0(
      "The fitted object has missing or unsupported 'formula_name_map' ",
      "metadata. Refit the model with this version of BayesTools."
    ),
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )

  missing_map <- list(mu = result$formula_design)
  missing_map$mu$name_map <- NULL
  attr(broken, "formula_design") <- missing_map
  broken <- .bt_attach_fit_contract(broken)
  expect_error(
    JAGS_formula_name_map(broken),
    "name-map metadata are missing",
    class = "BayesTools_refit_required"
  )
  expect_error(
    .bt_build_parameter_coordinates(
      columns = "mu_intercept", formula_design = missing_map
    ),
    "name-map metadata are missing",
    class = "BayesTools_refit_required"
  )
})

test_that("refit errors have the class BayesTools_refit_required", {

  refit_class <- c("BayesTools_refit_required", "error", "condition")

  # the helper builds the message as stop() does and raises without a call
  parts <- list("Metadata of '", c("mu", "_x"), "' are missing. ", NULL,
                "Refit the model ", 2L, ".")
  condition <- tryCatch(
    do.call(.bt_stop_refit_required, parts),
    error = identity
  )
  expect_identical(class(condition), refit_class)
  expect_identical(
    conditionMessage(condition),
    conditionMessage(tryCatch(
      do.call(stop, c(parts, call. = FALSE)),
      error = identity
    ))
  )
  expect_identical(
    conditionMessage(condition),
    "Metadata of 'mu_x' are missing. Refit the model 2."
  )
  expect_null(conditionCall(condition))
  condition <- tryCatch(
    .bt_stop_refit_required(
      "Refit with monitoring.",
      class = "BayesTools_refit_monitoring"
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_refit_monitoring", refit_class)
  )

  # the fit-contract gate refuses a BayesTools 0.3.0 fit (no parameter map,
  # fit contract, or draw geometry) at every public entry
  legacy <- structure(list(), class = c("runjags", "BayesTools_fit"))
  missing_map <- paste0(
    "The fitted object does not contain parameter-map metadata. ",
    "Refit the model with the current BayesTools version."
  )
  missing_contract <- paste0(
    "The fitted object does not contain a supported schema contract. ",
    "Refit the model with this version of BayesTools."
  )
  calls <- list(
    .bt_require_fit_contract = function() .bt_require_fit_contract(legacy),
    parameter_map            = function() parameter_map(legacy),
    parameter_catalog        = function() parameter_catalog(legacy),
    JAGS_check_convergence   = function() JAGS_check_convergence(legacy),
    as_mixed_posteriors      = function() as_mixed_posteriors(legacy, "mu"),
    JAGS_fit_contract        = function() JAGS_fit_contract(legacy),
    JAGS_draw_geometry       = function() JAGS_draw_geometry(legacy)
  )
  expected <- c(
    .bt_require_fit_contract = missing_map,
    parameter_map            = missing_map,
    parameter_catalog        = missing_contract,
    JAGS_check_convergence   = missing_map,
    as_mixed_posteriors      = missing_map,
    JAGS_fit_contract        = missing_contract,
    JAGS_draw_geometry       = missing_contract
  )
  for(name in names(calls)){
    condition <- tryCatch(calls[[name]](), error = identity)
    expect_identical(class(condition), refit_class, info = name)
    expect_identical(conditionMessage(condition), expected[[name]], info = name)
  }

  # a plain runjags object, and one unsupported component of a current
  # contract
  condition <- tryCatch(
    JAGS_check_convergence(structure(list(), class = "runjags")),
    error = identity
  )
  expect_identical(class(condition), refit_class)
  expect_identical(
    conditionMessage(condition),
    paste0(
      "'fit' must be a 'BayesTools_fit' created by JAGS_fit(). ",
      "Refit the model with this version of BayesTools."
    )
  )
  current <- .bt_attach_fit_contract(structure(list(), class = "BayesTools_fit"))
  attr(current, "fit_contract")$draw_geometry_version <- 0L
  expect_silent(JAGS_validate_fit_contract(current, requires = "parameter_map"))
  expect_error(
    JAGS_validate_fit_contract(current, requires = "draw_geometry"),
    paste0(
      "The fitted object has missing or unsupported 'draw_geometry' ",
      "metadata. Refit the model with this version of BayesTools."
    ),
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
})

# String constants of the package sources that contain the word "refit" in
# any case ("refitting" and "refitted" are other words): the file, line, and
# whether a call of .bt_stop_refit_required() encloses the string.
.refit_message_sites <- function(files){

  sites <- lapply(files, function(file){
    pd <- .test_source_parse_data(file)
    strings <- pd[pd$token == "STR_CONST" &
                    grepl("\\brefit\\b", pd$text, ignore.case = TRUE, perl = TRUE), ]
    if(nrow(strings) == 0L){
      return(NULL)
    }
    helper <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" &
                          pd$text == ".bt_stop_refit_required"]
    helper_calls <- pd$parent[match(helper, pd$id)]
    enclosed <- vapply(strings$id, function(id){
      while(length(id) == 1L && id > 0L){
        if(id %in% helper_calls){
          return(TRUE)
        }
        id <- pd$parent[pd$id == id]
      }
      FALSE
    }, logical(1))
    data.frame(
      file     = basename(file),
      line     = strings$line1,
      enclosed = enclosed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, sites)
}

test_that("every message that asks for a refit is raised as BayesTools_refit_required", {

  # the scanner finds refit messages of every call form and their enclosure
  planted <- tempfile(fileext = ".R")
  on.exit(unlink(planted), add = TRUE)
  writeLines(c(
    'f <- function(){',
    '  stop("Refit the model.", call. = FALSE)',
    '  .bt_stop_refit_required("Refit the model.")',
    '  .bt_stop_refit_required(paste0("a", " or refit the model."), class = "x")',
    '  warning(paste0("Please ", "REFIT."))',
    '  message("refitting and refitted are other words")',
    '}'
  ), planted)
  planted_sites <- .refit_message_sites(planted)
  expect_identical(planted_sites$line, c(2L, 3L, 4L, 5L))
  expect_identical(planted_sites$enclosed, c(FALSE, TRUE, TRUE, FALSE))

  package_r_dir <- file.path(testthat::test_path("..", ".."), "R")
  skip_if_not(dir.exists(package_r_dir), "package sources are unavailable")
  sites <- .refit_message_sites(
    list.files(package_r_dir, pattern = "[.][Rr]$", full.names = TRUE)
  )
  expect_gt(sum(sites$enclosed), 0L)
  labels <- paste0(sites$file, ":", sites$line)
  expect_identical(labels[!sites$enclosed], character())
})
