skip_if_not_test_profile("unit")

# attr() call sites with a literal attribute name in the package sources: the
# file, line, attribute name, whether 'exact = TRUE' is given, and whether the
# call is an assignment target.
.draws_metadata_attr_sites <- function(files){

  sites <- lapply(files, function(file){
    exprs <- parse(file, keep.source = TRUE)
    pd <- utils::getParseData(exprs, includeText = TRUE)
    attr_tokens <- pd[pd$text == "attr" &
                        pd$token %in% c("SYMBOL_FUNCTION_CALL", "SYMBOL"), ]
    rows <- lapply(seq_len(nrow(attr_tokens)), function(k){
      token <- attr_tokens[k, ]
      call_id <- pd$parent[pd$id == token$parent]
      kids <- pd[pd$parent == call_id, ]
      kids <- kids[order(kids$line1, kids$col1), ]
      literal <- NULL
      for(r in which(kids$token == "expr")){
        sub <- pd[pd$parent == kids$id[r], ]
        if(nrow(sub) == 1L && sub$token == "STR_CONST"){
          literal <- gsub('^"|"$', "", sub$text)
          break
        }
      }
      if(is.null(literal)){
        return(NULL)
      }
      exact_at <- which(kids$token == "SYMBOL_SUB" & kids$text == "exact")
      parent_kids <- pd[pd$parent == pd$parent[pd$id == call_id], ]
      parent_kids <- parent_kids[order(parent_kids$line1, parent_kids$col1), ]
      data.frame(
        file   = basename(file),
        line   = token$line1,
        name   = literal,
        exact  = length(exact_at) == 1L &&
          identical(kids$text[exact_at + 2L], "TRUE"),
        assign = token$token == "SYMBOL_FUNCTION_CALL" &&
          nrow(parent_kids) == 3L && (
            (parent_kids$id[1L] == call_id &&
               parent_kids$token[2L] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")) ||
              (parent_kids$id[3L] == call_id &&
                 parent_kids$token[2L] == "RIGHT_ASSIGN")
          ),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  do.call(rbind, sites)
}

.draws_metadata_source_files <- function(){

  package_r_dir <- file.path(testthat::test_path("..", ".."), "R")
  skip_if_not(dir.exists(package_r_dir), "package sources are unavailable")
  list.files(package_r_dir, pattern = "[.][Rr]$", full.names = TRUE)
}

test_that("draw metadata is read and written only through its accessors", {

  sites <- .draws_metadata_attr_sites(.draws_metadata_source_files())

  # every storage name of the draw metadata and the conditioning fields;
  # 'formula_scale' also names the fit attribute of the formula scaling
  metadata_names <- setdiff(c(
    unname(BayesTools:::.bt_meta_attribute_names),
    BayesTools:::.bt_meta_condition_names
  ), "formula_scale")
  # 'conditional' and 'formula_parameter' also name attributes of the
  # ensemble-inference objects of ensemble_inference(), read by its table
  allowed <- sites$file == "draws-metadata.R" |
    (sites$file %in% c("model-averaging.R", "summary-tables-ensemble.R") &
       sites$name %in% c("conditional", "formula_parameter"))
  outside <- sites$name %in% metadata_names & !allowed
  labels <- paste0(sites$file, ":", sites$line, " ", sites$name)
  expect_identical(labels[outside], character())
})

test_that("attributes are never read by partial name matching", {

  sites <- .draws_metadata_attr_sites(.draws_metadata_source_files())
  names <- unique(sites$name)
  is_prefix <- vapply(names, function(name){
    any(startsWith(names, name) & names != name)
  }, logical(1))
  risky <- !sites$assign & !sites$exact & sites$name %in% names[is_prefix]
  labels <- paste0(sites$file, ":", sites$line, " ", sites$name)
  expect_identical(labels[risky], character())
})

test_that("draw-metadata accessors are exact and validate the field", {

  x <- 1:3
  x <- BayesTools:::.bt_meta_set(x, "prior_context", list(context = TRUE))
  # the prior-density context never answers for the prior density
  expect_null(BayesTools:::.bt_meta_get(x, "prior_density"))
  expect_identical(
    BayesTools:::.bt_meta_get(x, "prior_context"),
    list(context = TRUE)
  )
  x <- BayesTools:::.bt_meta_set(x, "prior_context", NULL)
  expect_null(BayesTools:::.bt_meta_get(x, "prior_context"))

  expect_error(
    BayesTools:::.bt_meta_get(x, "prior"),
    "'prior' is not a draw-metadata field.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_meta_set(x, "posterior_supp", NULL),
    "'posterior_supp' is not a draw-metadata field.",
    fixed = TRUE
  )

  x <- BayesTools:::.bt_meta_set(
    x, "condition", list(conditional = "mu", conditional_rule = "AND")
  )
  expect_identical(BayesTools:::.bt_meta_condition(x, "conditional"), "mu")
  expect_null(BayesTools:::.bt_meta_condition(x, "condition_key"))
  expect_error(
    BayesTools:::.bt_meta_set(x, "condition", list(conditions = "mu")),
    "Draw conditioning metadata must be a named list of condition fields.",
    fixed = TRUE
  )
})

test_that("Savage-Dickey reports missing prior densities, not their context", {

  posterior <- stats::qnorm(stats::ppoints(200))
  posterior <- BayesTools:::.posterior_atoms_set(
    posterior,
    BayesTools:::.posterior_atoms_new(source = "test")
  )
  context <- BayesTools:::.prior_density_build_context(
    prior_list   = list(mu = prior("normal", list(0, 1))),
    column_names = "mu"
  )
  posterior <- BayesTools:::.bt_meta_set(posterior, "prior_context", context)
  class(posterior) <- c("marginal_posterior.simple", "marginal_posterior")

  expect_error(
    Savage_Dickey_BF(posterior),
    "there are no prior densities for the posterior distribution",
    fixed = TRUE
  )
})
