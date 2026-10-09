skip_if_not_test_profile("unit")

test_that("empty public measure refusal tables clear metadata after validation", {
  x <- c(1, 2); attr(x, "parameter") <- "theta"
  table <- data.frame(column = "theta", measure = "atoms", reason = "Known unavailable")
  posterior_metadata(x, "measure_unavailable") <- table
  expect_identical(posterior_metadata(x, "measure_unavailable"), table)
  for(value in list(table[FALSE, , drop = FALSE], NULL)){
    y <- x
    posterior_metadata(y, "measure_unavailable") <- value
    expect_null(posterior_metadata(y, "measure_unavailable"))
    expect_identical(as.numeric(y), c(1, 2))
  }
  original <- serialize(x, NULL)
  condition <- tryCatch({posterior_metadata(x, "measure_unavailable") <- "a"; NULL}, error = identity)
  expect_s3_class(condition, "error")
  expect_identical(conditionMessage(condition), "Draw metadata 'measure_unavailable' is invalid: it must be a plain column/measure/reason table with unique target measures.")
  expect_null(conditionCall(condition))
  expect_identical(serialize(x, NULL), original)
})

test_that("affine hypothesis row recipes follow every selected draw", {
  x <- .bt_meta_set(c(1, 2, 3), "hypothesis_evaluation", list(numerator = c(1, 2, 3), divisor = 1, weights = c(theta = 1), offset = 0))
  original <- .bt_meta_get(x, "hypothesis_evaluation")
  for(rows in list(c(3L, 1L), c(2L, 2L, 1L), integer())){
    y <- .bt_draws_subset_rows(x, rows)
    recipe <- .bt_meta_get(y, "hypothesis_evaluation")
    expect_identical(recipe$numerator, original$numerator[rows])
    expect_identical(as.numeric(y), recipe$numerator / recipe$divisor)
    for(field in c("divisor", "weights", "offset")) expect_identical(recipe[[field]], original[[field]])
  }
  expect_identical(.bt_meta_get(x, "hypothesis_evaluation"), original)
})
source(testthat::test_path("common-functions.R"))

test_that("the public ordered source accessor preserves producer validation", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, .5)))
  x <- fixture$samples$mu_f
  source <- posterior_metadata(x, "ordered_source")
  expect_identical(source, .bt_meta_get(x, "ordered_source"))
  expect_null(posterior_metadata(1:3, "ordered_source"))
  posterior_metadata(x, "ordered_source") <- NULL
  expect_null(posterior_metadata(x, "ordered_source"))
  posterior_metadata(x, "ordered_source") <- source
  expect_identical(posterior_metadata(x, "ordered_source"), source)
  expect_error(posterior_metadata(x, "ordered_source") <- list(),
    "Draw metadata 'ordered_source' is invalid", fixed = TRUE)
  misaligned <- .bt_ordered_source_subset(source, seq_len(nrow(x) - 1L))
  expect_error(posterior_metadata(x, "ordered_source") <- misaligned,
    "Draw metadata 'ordered_source' must have one row per draw.", fixed = TRUE)
  malformed <- source
  malformed$models[[1L]]$parameter <- "wrong"
  expect_error(posterior_metadata(x, "ordered_source") <- malformed,
    "its models must contain authoritative bound ordered specifications", fixed = TRUE)
  stale <- x
  stale[1L, 1L] <- stale[1L, 1L] + 1
  expect_error(posterior_metadata(stale, "ordered_source"), class = "BayesTools_stale_metadata")
  expect_error(posterior_metadata(stale, "ordered_source") <- source, class = "BayesTools_stale_metadata")
})

test_that("scalar marginal atom declarations validate positional schema and transform without joint inference", {
  empty <- .posterior_atoms_new(column_names="a")
  point <- .posterior_atoms_new(matrix(2,1,1),.3,column_names="b")
  unknown <- setNames(list(NULL,point),c("a","b"))
  atoms <- .posterior_atoms_new(column_names=c("a","b"),marginals=unknown)
  expect_identical(.posterior_atoms_from_attribute(atoms)$marginals,unknown)
  expect_null(.posterior_atoms_for_column(atoms,1L))
  expect_identical(.posterior_atoms_for_column(atoms,2L),point)
  after <- .posterior_atoms_new(column_names=c("a","b"),marginals=setNames(list(empty,NULL),c("a","b")))
  expect_identical(length(after$marginals),2L)
  expect_null(.posterior_atoms_for_column(after,2L))
  bad <- list(list(empty,point),setNames(list(empty),"a"),setNames(list(empty,point),c("b","a")),
    setNames(list(empty,point),c("a","a")))
  for(marginals in bad) expect_error(.posterior_atoms_new(column_names=c("a","b"),marginals=marginals),
    "must name every location-matrix column in order",fixed=TRUE)
  multi <- .posterior_atoms_new(column_names=c("a","b"))
  nested <- .posterior_atoms_new(column_names="a",marginals=list(a=empty))
  for(marginal in list(multi,nested)) expect_error(.posterior_atoms_new(column_names="a",marginals=list(a=marginal)),
    "single column without nested marginals",fixed=TRUE)
  draws <- structure(matrix(seq_len(8),4,dimnames=list(NULL,c("a","b"))),class=c("mixed_posteriors","matrix","array"))
  declared <- .posterior_atoms_set(draws,.posterior_atoms_new(column_names=c("a","b")))
  expect_true(posterior_atoms_free(declared))
  expect_false(posterior_atoms_free(.posterior_atoms_set(draws,atoms)))
  expect_false(posterior_atoms_free(.posterior_atoms_set(draws,after)))
  both_empty <- .posterior_atoms_new(column_names=c("a","b"),marginals=list(a=empty,b=.posterior_atoms_new(column_names="b")))
  expect_true(posterior_atoms_free(.posterior_atoms_set(draws,both_empty)))
  joint <- .posterior_atoms_new(matrix(c(1,2),1,2),.2,column_names=c("a","b"),marginals=unknown)
  expect_false(posterior_atoms_free(.posterior_atoms_set(draws,joint)))
  renamed <- .posterior_atoms_rename_columns(atoms,c("x","y"))
  expect_identical(names(renamed$marginals),c("x","y"))
  expect_identical(colnames(renamed$marginals$y$locations),"y")
  subset <- .bt_draws_subset_rows(.posterior_atoms_set(draws,atoms),c(1L,3L))
  expect_identical(.posterior_atoms_get(subset)$marginals,unknown)
  transformed <- .posterior_atoms_transform(atoms,"lin",list(a=1,b=-2))
  expect_null(transformed$marginals$a)
  expect_identical(as.numeric(transformed$marginals$b$locations),-3)
  expect_identical(transformed$marginals$b$mass,.3)
  design <- rbind(c(0,-2),c(1,1),c(0,0),c(1,0))
  linear <- .posterior_atoms_linear_transform(atoms,design)
  expect_identical(names(linear$marginals),paste0("value",1:4))
  expect_identical(as.numeric(linear$marginals[[1]]$locations),-4)
  expect_identical(linear$marginals[[1]]$mass,.3)
  expect_null(linear$marginals[[2]])
  expect_identical(as.numeric(linear$marginals[[3]]$locations),0)
  expect_identical(linear$marginals[[3]]$mass,1)
  expect_null(linear$marginals[[4]])
  expect_identical(nrow(linear$locations),0L)
  preserved <- .posterior_atoms_linear_transform(joint,matrix(c(1,1),1,2),"sum")
  expect_identical(as.numeric(preserved$locations),3)
  expect_identical(preserved$mass,.2)
  expect_null(preserved$marginals[[1]])
  unnamed <- .posterior_atoms_new(matrix(2,1,1),.3)
  expect_identical(as.numeric(.posterior_atoms_linear_transform(unnamed,matrix(-2,1,1))$locations),-4)
  expect_null(colnames(.posterior_atoms_linear_transform(unnamed,matrix(-2,1,1))$locations))
})

# The scans read the parse data of a source file, which
# .test_source_parse_data() (common-functions.R) parses once per test run, with
# the row indices of the children of each node where a scan walks the tree.
# Only terminal tokens carry their text: the scans read tokens, never the text
# of whole expressions.

# The rows of the children of a parse-data node, in source order.
.draws_metadata_children <- function(pd, id){

  rows <- attr(pd, "children")[[as.character(id)]]
  rows[order(pd$line1[rows], pd$col1[rows])]
}

# Whether a parse-data node is the literal TRUE: the token itself, or an
# expression whose only token it is.
.draws_metadata_is_true <- function(pd, id){

  if(identical(pd$text[match(id, pd$id)], "TRUE")){
    return(TRUE)
  }
  value <- .draws_metadata_children(pd, id)
  length(value) == 1L && identical(pd$text[value], "TRUE")
}

# attr() call sites with a literal attribute name in the package sources: the
# file, line, attribute name, whether 'exact = TRUE' is given, and whether the
# call is an assignment target.
.draws_metadata_attr_sites <- function(files){

  sites <- lapply(files, function(file){
    pd <- .test_source_parse_data(file, children = TRUE)
    attr_tokens <- which(pd$text == "attr" &
                           pd$token %in% c("SYMBOL_FUNCTION_CALL", "SYMBOL"))
    rows <- lapply(attr_tokens, function(token){
      call_id <- pd$parent[match(pd$parent[token], pd$id)]
      kids <- .draws_metadata_children(pd, call_id)
      literal <- NULL
      for(kid in kids[pd$token[kids] == "expr"]){
        sub <- .draws_metadata_children(pd, pd$id[kid])
        if(length(sub) == 1L && pd$token[sub] == "STR_CONST"){
          literal <- gsub('^"|"$', "", pd$text[sub])
          break
        }
      }
      if(is.null(literal)){
        return(NULL)
      }
      exact_at <- which(pd$token[kids] == "SYMBOL_SUB" & pd$text[kids] == "exact")
      parent_kids <- .draws_metadata_children(pd, pd$parent[match(call_id, pd$id)])
      data.frame(
        file   = basename(file),
        line   = pd$line1[token],
        name   = literal,
        exact  = length(exact_at) == 1L &&
          .draws_metadata_is_true(pd, pd$id[kids[exact_at + 2L]]),
        assign = pd$token[token] == "SYMBOL_FUNCTION_CALL" &&
          length(parent_kids) == 3L && (
            (pd$id[parent_kids[1L]] == call_id &&
               pd$token[parent_kids[2L]] %in% c("LEFT_ASSIGN", "EQ_ASSIGN")) ||
              (pd$id[parent_kids[3L]] == call_id &&
                 pd$token[parent_kids[2L]] == "RIGHT_ASSIGN")
          ),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  do.call(rbind, sites)
}

# The attr() sites of the package sources, scanned once for the tests of this
# file that read them.
.draws_metadata_package_attr_sites <- local({
  sites <- NULL
  function(){
    if(is.null(sites)){
      sites <<- .draws_metadata_attr_sites(.draws_metadata_source_files())
    }
    sites
  }
})

.draws_metadata_source_files <- function(){

  package_r_dir <- file.path(testthat::test_path("..", ".."), "R")
  skip_if_not(dir.exists(package_r_dir), "package sources are unavailable")
  list.files(package_r_dir, pattern = "[.][Rr]$", full.names = TRUE)
}

test_that("draw metadata is read and written only through its accessors", {

  sites <- .draws_metadata_package_attr_sites()

  # the metadata container and the former free metadata attributes;
  # 'formula_scale' also names the fit attribute of the formula scaling
  metadata_names <- setdiff(c(
    BayesTools:::.bt_meta_attribute,
    BayesTools:::.bt_meta_legacy_names
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

# Attribute names set through structure(x, <name> = ) or accessed as
# attributes(x)$<name> / attributes(x)[["<name>"]] in the package sources:
# the file, line, and attribute name.
.draws_metadata_named_sites <- function(files){

  sites <- lapply(files, function(file){
    pd <- .test_source_parse_data(file)
    call_of <- function(function_name){
      symbols <- pd$parent[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == function_name]
      pd$parent[match(symbols, pd$id)]
    }
    rows <- list()
    for(call_id in call_of("structure")){
      arguments <- pd[pd$parent == call_id & pd$token == "SYMBOL_SUB", ]
      if(nrow(arguments) > 0L){
        rows[[length(rows) + 1L]] <- data.frame(
          file = basename(file), line = arguments$line1,
          name = gsub("`", "", arguments$text), stringsAsFactors = FALSE
        )
      }
    }
    for(call_id in call_of("attributes")){
      outer <- pd[pd$parent == pd$parent[pd$id == call_id], ]
      outer <- outer[order(outer$line1, outer$col1), ]
      if(nrow(outer) < 3L || outer$id[1L] != call_id){
        next
      }
      name <- if(identical(outer$token[2L], "'$'")){
        outer$text[3L]
      }else if(identical(outer$token[2L], "LBB")){
        index <- pd[pd$parent == outer$id[3L], ]
        if(nrow(index) == 1L && identical(index$token, "STR_CONST")) index$text else NA_character_
      }else{
        NA_character_
      }
      if(!is.na(name)){
        rows[[length(rows) + 1L]] <- data.frame(
          file = basename(file), line = outer$line1[1L],
          name = gsub("^[\"'`]|[\"'`]$", "", name), stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  })
  do.call(rbind, sites)
}

test_that("draw metadata is never set by structure() or read from attributes()", {

  metadata_names <- setdiff(c(
    BayesTools:::.bt_meta_attribute,
    BayesTools:::.bt_meta_legacy_names
  ), "formula_scale")

  # the scanners find every access form they check
  planted <- tempfile(fileext = ".R")
  on.exit(unlink(planted), add = TRUE)
  writeLines(c(
    'a <- attr(x, "prior_density")',
    'b <- structure(x, posterior_atoms = atoms, class = "y")',
    'd <- attributes(x)$undefined_draws',
    'e <- attributes(x)[["models_ind"]]'
  ), planted)
  expect_identical(.draws_metadata_attr_sites(planted)$name, "prior_density")
  planted_sites <- .draws_metadata_named_sites(planted)
  expect_identical(
    planted_sites$name[planted_sites$name %in% metadata_names],
    c("posterior_atoms", "undefined_draws", "models_ind")
  )

  sites <- .draws_metadata_named_sites(.draws_metadata_source_files())
  outside <- sites$name %in% metadata_names & sites$file != "draws-metadata.R"
  labels <- paste0(sites$file, ":", sites$line, " ", sites$name)
  expect_identical(labels[outside], character())
})

test_that("attributes are never read by partial name matching", {

  # the scanner reads exact literal 'exact = TRUE' arguments and assignments
  planted <- tempfile(fileext = ".R")
  on.exit(unlink(planted), add = TRUE)
  writeLines(c(
    'a <- attr(x, "prior", exact = TRUE)',
    'b <- attr(x, "prior")',
    'attr(x, "prior") <- 1',
    'd <- attr(x, "prior", exact = (TRUE))'
  ), planted)
  planted_sites <- .draws_metadata_attr_sites(planted)
  expect_identical(planted_sites$exact, c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(planted_sites$assign, c(FALSE, FALSE, TRUE, FALSE))

  sites <- .draws_metadata_package_attr_sites()
  names <- unique(sites$name)
  is_prefix <- vapply(names, function(name){
    any(startsWith(names, name) & names != name)
  }, logical(1))
  risky <- !sites$assign & !sites$exact & sites$name %in% names[is_prefix]
  labels <- paste0(sites$file, ":", sites$line, " ", sites$name)
  expect_identical(labels[risky], character())
})

test_that("selected fitted terms use the exact order attribute", {
  fitted <- stats::terms(~ x + y + x:y)
  design <- list(terms = fitted, assign = 0:3, raw_column_names = c("intercept", "x", "y", "x:y"))
  selected <- .bt_formula_selected_terms(~ x + x:y, design)$terms
  expect_identical(attr(selected, "order", exact = TRUE), c(1L, 2L))
  attr(design$terms, "order") <- NULL
  attr(design$terms, "ordered_extra") <- c(7L, 8L, 9L)
  selected <- .bt_formula_selected_terms(~ x + x:y, design)$terms
  expect_null(attr(selected, "order", exact = TRUE))
})

test_that("draw-metadata accessors are exact and validate the field", {

  x <- 1:3
  context <- BayesTools:::.prior_density_build_context(
    prior_list   = list(mu = prior("normal", list(0, 1))),
    column_names = "mu"
  )
  x <- BayesTools:::.bt_meta_set(x, "prior_context", context)
  # the prior-density context never answers for the prior density
  expect_null(BayesTools:::.bt_meta_get(x, "prior_density"))
  expect_identical(BayesTools:::.bt_meta_get(x, "prior_context"), context)
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
    "Draw metadata 'condition' is invalid: it must be a named list of condition fields.",
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

test_that("draw metadata is one validated container", {

  x <- stats::rnorm(20)
  posterior_metadata(x, "atoms") <- posterior_atom_attribute()
  posterior_metadata(x, "support") <- posterior_support_attribute(c(-Inf, Inf))
  expect_identical(names(attributes(x)), "bayestools_meta")
  meta <- attr(x, "bayestools_meta", exact = TRUE)
  expect_s3_class(meta, "BayesTools_draw_metadata")
  # the fields, and the fingerprint of the values that they describe
  expect_identical(names(meta), c("atoms", "support", "fingerprint"))
  expect_s3_class(posterior_metadata(x, "atoms"), "BayesTools_posterior_atoms")

  # removing the last field removes the container
  posterior_metadata(x, "atoms") <- NULL
  posterior_metadata(x, "support") <- NULL
  expect_null(attributes(x))

  # values are validated when they are set
  expect_error(
    posterior_metadata(x, "atoms") <- list(declared = TRUE),
    "Posterior atom metadata must be created with 'posterior_atom_attribute()'.",
    fixed = TRUE
  )
  expect_error(
    posterior_metadata(x, "support") <- c(0, 1),
    "Posterior support metadata must be created with 'posterior_support_attribute()'.",
    fixed = TRUE
  )
  expect_error(
    posterior_metadata(x, "posterior_density") <- list(x = 1:3, y = 1:3),
    "Posterior density metadata must be created with 'posterior_density_attribute()'.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_meta_set(x, "condition", list(conditional_rule = "XOR")),
    "Draw metadata 'condition' is invalid: 'conditional_rule' must be 'AND' or 'OR'.",
    fixed = TRUE
  )
  # the public accessor covers the fields that downstream packages attach
  expect_error(posterior_metadata(x, "models_ind"), "'field'")
  expect_error(posterior_metadata(x, "component"), "'field'")

  # level weights of marginal posteriors and the prior densities of
  # as_mixed_posteriors() are read and set by downstream packages
  weights <- c(mu_intercept = 1, mu_x = .5)
  posterior_metadata(x, "linear_weights") <- weights
  expect_identical(posterior_metadata(x, "linear_weights"), weights)
  expect_error(
    posterior_metadata(x, "linear_weights") <- "mu_x",
    "Draw metadata 'linear_weights' is invalid: it must be a finite numeric vector or matrix.",
    fixed = TRUE
  )
  densities <- list(mu = BayesTools:::.prior_linear_combination_density(
    prior_list = list(mu = prior("normal", list(0, 1))),
    weights    = c(mu = 1)
  ))
  posterior_metadata(x, "prior_densities") <- densities
  expect_identical(posterior_metadata(x, "prior_densities"), densities)
  expect_error(
    posterior_metadata(x, "prior_densities") <- list(mu = 1),
    "Draw metadata 'prior_densities' is invalid: it must be a list of prior densities.",
    fixed = TRUE
  )
})

test_that("free metadata attributes of development versions are not read", {

  x <- stats::qnorm(stats::ppoints(100))
  attr(x, "posterior_atoms") <- posterior_atom_attribute()
  attr(x, "posterior_support") <- posterior_support_attribute(c(-Inf, Inf))
  expect_null(BayesTools:::.posterior_atoms_get(x))
  expect_null(BayesTools:::.posterior_support_get(x))
  expect_null(posterior_metadata(x, "atoms"))
})

test_that("posterior atoms come only from the atom metadata", {

  x <- c(rep(0, 20), stats::qnorm(stats::ppoints(80)))
  posterior_metadata(x, "prior_density") <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(source = prior("normal", list(0, 1))),
    weights    = c(source = 1)
  )
  posterior_metadata(x, "posterior_density") <- posterior_density_attribute(
    x            = seq(-3, 3, length.out = 61),
    y            = stats::dnorm(seq(-3, 3, length.out = 61)) * .8,
    method         = "test",
    density_method = "precomputed"
  )
  class(x) <- c("marginal_posterior.simple", "marginal_posterior")
  # a stored density (the continuous part) does not declare the posterior atoms
  expect_null(BayesTools:::.posterior_atoms_get(x))
  expect_error(Savage_Dickey_BF(x), "Posterior atom status is unknown", fixed = TRUE)

  posterior_metadata(x, "atoms") <- posterior_atom_attribute(data.frame(x = 0, mass = .2))
  expect_error(
    Savage_Dickey_BF(x),
    "declared point mass at the exact null",
    fixed = TRUE
  )
})

test_that("a posterior point mass at the null stops Savage-Dickey with its class", {

  x <- c(rep(0, 20), stats::qnorm(stats::ppoints(80)))
  posterior_metadata(x, "prior_density") <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(source = prior("normal", list(0, 1))),
    weights    = c(source = 1)
  )
  posterior_metadata(x, "atoms") <- posterior_atom_attribute(data.frame(x = 0, mass = .2))
  class(x) <- c("marginal_posterior.simple", "marginal_posterior")

  condition <- tryCatch(Savage_Dickey_BF(x), error = function(e) e)
  expect_s3_class(condition, "BayesTools_posterior_point_mass_at_null")
  expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
  expect_identical(
    conditionMessage(condition),
    paste0(
      "The posterior contains a declared point mass at the exact null ",
      "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
    )
  )
  # an atom elsewhere leaves the ratio at the null defined
  expect_true(is.finite(Savage_Dickey_BF(x, null_hypothesis = .5)))
})

test_that("producers store their draw metadata in the container only", {

  set.seed(1)
  posterior <- cbind(mu = c(rep(0, 10), stats::rnorm(30)),
                     mu_indicator = rep(c(0, 1), c(10, 30)))
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(
    mu = prior_spike_and_slab(prior("normal", list(0, 1)))
  )
  fit <- attach_test_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_parameter_map(fit, monitor_names = colnames(posterior))
  mixed <- as_mixed_posteriors(fit, "mu")

  element_attributes <- names(attributes(mixed$mu))
  expect_true("bayestools_meta" %in% element_attributes)
  expect_length(intersect(element_attributes, BayesTools:::.bt_meta_legacy_names), 0L)
  expect_length(intersect(names(attributes(mixed)), BayesTools:::.bt_meta_legacy_names), 0L)
  atoms <- posterior_metadata(mixed$mu, "atoms")
  expect_equal(as.numeric(atoms$locations[, 1L]), 0)
  expect_equal(atoms$mass, 10 / 40)

  marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  expect_length(intersect(names(attributes(marginal)), BayesTools:::.bt_meta_legacy_names), 0L)
  expect_s3_class(posterior_metadata(marginal, "prior_density"), "prior_density")
  expect_s3_class(posterior_metadata(marginal, "atoms"), "BayesTools_posterior_atoms")
})

.draws_metadata_mixed_for_test <- function(){

  set.seed(1)
  posterior <- cbind(
    mu           = c(rep(0, 100), stats::rnorm(300, .5, .2)),
    mu_indicator = rep(c(0, 1), c(100, 300)),
    sigma        = stats::rlnorm(400, 0, .2)
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(
    mu    = prior_spike_and_slab(prior("normal", list(0, 1))),
    sigma = prior("lognormal", list(0, 1))
  )
  fit <- attach_test_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_parameter_map(fit, monitor_names = colnames(posterior))
  as_mixed_posteriors(fit, c("mu", "sigma"))
}

test_that("producers declare averaged draws and posterior_atoms_free() reads the declared atoms", {

  set.seed(1)
  posterior <- cbind(
    mu           = c(rep(0, 100), stats::rnorm(300, .5, .2)),
    mu_indicator = rep(c(0, 1), c(100, 300)),
    sigma        = stats::rlnorm(400, 0, .2)
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(
    mu    = prior_spike_and_slab(prior("normal", list(0, 1))),
    sigma = prior("lognormal", list(0, 1))
  )
  fit <- attach_test_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_parameter_map(fit, monitor_names = colnames(posterior))

  averaged <- as_mixed_posteriors(fit, c("mu", "sigma"))
  conditional <- as_mixed_posteriors(fit, c("mu", "sigma"), conditional = "mu")
  expect_true(posterior_metadata(averaged, "condition")$averaged)
  expect_true(posterior_metadata(averaged$sigma, "condition")$averaged)
  expect_false(posterior_metadata(conditional, "condition")$averaged)
  expect_false(posterior_metadata(conditional$sigma, "condition")$averaged)
  expect_true(posterior_metadata(
    marginal_posterior(averaged, "sigma", prior_samples = TRUE), "condition"
  )$averaged)
  expect_false(posterior_metadata(
    marginal_posterior(conditional, "sigma", prior_samples = TRUE), "condition"
  )$averaged)
  expect_error(
    .bt_meta_set(averaged$sigma, "condition", list(averaged = NA)),
    "Draw metadata 'condition' is invalid: 'averaged' must be TRUE or FALSE.",
    fixed = TRUE
  )

  # level comparisons read 'averaged', not the condition keys
  level <- function(values, key){
    .bt_meta_update(
      structure(values, class = c("marginal_posterior.simple", "numeric")),
      condition = list(condition_key = key, averaged = TRUE)
    )
  }
  levels <- list(a = level(stats::rnorm(20), "first"), b = level(stats::rnorm(20), "second"))
  expect_true(BayesTools:::.hypothesis_validate_level_conditionals(levels, "mu", c("a", "b")))
  levels$b <- .bt_meta_set(levels$b, "condition", list(condition_key = "second", averaged = FALSE))
  expect_error(
    BayesTools:::.hypothesis_validate_level_conditionals(levels, "mu", c("a", "b")),
    "different conditional posterior subsets", fixed = TRUE
  )

  # the spike of 'mu' is a declared atom; 'sigma' is declared atom-free
  expect_false(posterior_atoms_free(averaged$mu))
  expect_true(posterior_atoms_free(averaged$sigma))
  expect_true(posterior_atoms_free(conditional$sigma))
  undeclared <- structure(stats::rnorm(10), class = c("marginal_posterior.simple", "marginal_posterior"))
  expect_false(posterior_atoms_free(undeclared))
  posterior_metadata(undeclared, "atoms") <- posterior_atom_attribute()
  expect_true(posterior_atoms_free(undeclared))
  expect_error(posterior_atoms_free(as.numeric(averaged$sigma)),
               "'posterior_atoms_free' requires BayesTools posterior draws, not plain numeric draws",
               fixed = TRUE)
  expect_error(posterior_atoms_free(averaged$sigma + 1), "plain numeric draws", fixed = TRUE)
})

test_that("arithmetic, math, and subsetting of draws return plain numerics", {

  mixed    <- .draws_metadata_mixed_for_test()
  mp_x     <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  mp_sigma <- marginal_posterior(mixed, "sigma", prior_samples = TRUE)

  plain <- list(
    mp_x * 2, 2 * mp_x, mp_x + 1, -mp_sigma, mp_sigma / 10, mp_x > 0,
    log(mp_sigma), exp(mp_x), round(mp_x, 2), cumsum(mp_x),
    mixed$mu * 2, log(mixed$sigma), c(mp_x, mp_x), as.numeric(mp_x),
    mp_x[1:10], mixed$mu[1:10]
  )
  for(value in plain){
    expect_null(attributes(value))
  }
  # draws of different classes combine without incompatible-method dispatch
  expect_silent(combined <- mixed$mu - mp_sigma)
  expect_null(attributes(combined))

  # matrices keep their dimensions only
  draws <- structure(
    matrix(1:4, 2, dimnames = list(NULL, c("a", "b"))),
    class = c("mixed_posteriors", "mixed_posteriors.vector")
  )
  draws <- BayesTools:::.bt_meta_set(draws, "atoms", posterior_atom_attribute())
  expect_identical(names(attributes(draws * 2)), c("dim", "dimnames"))
  expect_identical(names(attributes(sqrt(draws))), c("dim", "dimnames"))

  # transformed inference goes through marginal_posterior(transformation = ):
  # BF(2 mu = 0.5) equals BF(mu = 0.25) exactly (kernel bandwidths scale)
  mp_2x <- marginal_posterior(
    mixed, "mu", prior_samples = TRUE,
    transformation = "lin", transformation_arguments = list(a = 0, b = 2)
  )
  expect_equal(
    as.numeric(Savage_Dickey_BF(mp_2x, 0.5)),
    as.numeric(Savage_Dickey_BF(mp_x, 0.25)),
    tolerance = 1e-10
  )

  # consumers that need the metadata stop on plain numeric draws
  for(value in list(mp_x * 2, mp_sigma / 10, mp_x + 1, log(mp_sigma), mp_sigma[1:100])){
    expect_error(
      Savage_Dickey_BF(value, .5),
      "'Savage_Dickey_BF' requires an object of class 'marginal_posterior', not plain numeric draws",
      fixed = TRUE
    )
  }
  # the error names the route for transformed draws
  expect_error(
    Savage_Dickey_BF(mp_x * 2, .5),
    paste0(
      "'Savage_Dickey_BF' requires an object of class 'marginal_posterior', not ",
      "plain numeric draws: arithmetic, mathematical functions, and subsetting ",
      "of posterior draws return plain numeric draws without their supports, ",
      "atoms, and prior densities. Use posterior_transform() (or ",
      "marginal_posterior(transformation = )) for transformed posterior ",
      "distributions."
    ),
    fixed = TRUE
  )
  scaled <- mixed
  scaled$mu <- scaled$mu * 2
  expect_error(
    plot_posterior(scaled, "mu", plot_type = "ggplot"),
    "The posterior samples of 'mu' must be created by 'mix_posteriors' or 'as_mixed_posteriors', not plain numeric draws",
    fixed = TRUE
  )
  expect_error(
    marginal_posterior(scaled, "mu"),
    "The posterior samples of 'mu' must be created by 'mix_posteriors' or 'as_mixed_posteriors', not plain numeric draws",
    fixed = TRUE
  )
  # subset draws no longer produce an empty plot
  expect_error(
    BayesTools:::.plot_data_samples.simple(
      list(mu = mixed$mu[1:100]), "mu", 64, NULL, NULL, FALSE
    ),
    "Posterior atom status is unknown",
    fixed = TRUE
  )
})

test_that("posterior_transform() maps spike-and-slab draws with their metadata", {

  mixed <- .draws_metadata_mixed_for_test()
  mp    <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  # a precomputed posterior ordinate at the null 0.5 (any positive value)
  posterior_metadata(mp, "posterior_ordinate") <- posterior_ordinate_attribute(
    value = .5, ordinate = .8, method = "test", density_method = "precomputed"
  )
  tr <- posterior_transform(mp, "exp")

  expect_equal(as.numeric(tr), exp(as.numeric(mp)), tolerance = 0)
  # support: the real line maps to (0, Inf) and the spike at 0 to 1
  support <- posterior_metadata(tr, "support")
  expect_identical(support$bounds, c(0, Inf))
  expect_identical(support$points, 1)
  # atoms: the location maps, the mass (100 of 400 draws) is kept
  atoms <- posterior_metadata(tr, "atoms")
  expect_equal(as.numeric(atoms$locations), 1)
  expect_equal(atoms$mass, .25)
  # prior density: 0.5 point mass at 1 and 0.5 x lognormal(0, 1) (the slab
  # normal(0, 1) through exp; prior inclusion probability mean(beta(1, 1)))
  prior_density <- posterior_metadata(tr, "prior_density")
  for(value in c(.2, .5, 1.5, 3, 10)){
    ordinate <- prior_density_ordinate(prior_density, value)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), .5 * stats::dlnorm(value),
                 tolerance = 1e-12, info = value)
  }
  expect_identical(prior_density_ordinate(prior_density, 1)$behavior, "point_mass")
  # the Savage-Dickey ratio is invariant at the mapped null (prior and
  # posterior ordinates share the Jacobian exp(0.5))
  transformed_ordinate <- posterior_metadata(tr, "posterior_ordinate")
  expect_equal(transformed_ordinate$value, exp(.5))
  expect_equal(transformed_ordinate$ordinate, .8 / exp(.5))
  expect_equal(
    as.numeric(Savage_Dickey_BF(tr, exp(.5), density_method = "precomputed")),
    as.numeric(Savage_Dickey_BF(mp, .5, density_method = "precomputed")),
    tolerance = 1e-12
  )
  # conditioning is kept; the quantity records the transformation, keeps its
  # label, and is no longer the catalog quantity
  expect_identical(posterior_metadata(tr, "condition"), posterior_metadata(mp, "condition"))
  quantities <- posterior_metadata(tr, "quantities")
  expect_identical(quantities$label_parts[[1L]]$transformation, c("none", "exp"))
  expect_identical(quantities$quantity_id, "")
  expect_identical(quantities$dependencies[[1L]], character())
  expect_identical(parameter_labels(quantities), parameter_labels(posterior_metadata(mp, "quantities")))
  # the linear-combination weights no longer apply
  expect_null(posterior_metadata(tr, "linear_weights"))
  expect_identical(BayesTools:::.bt_meta_get(tr, "joint_prior_transformation"), "exp")

  # marginal_posterior(transformation = ) is this transformation
  via_marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE,
                                     transformation = "exp")
  expect_identical(unclass(via_marginal), unclass(posterior_transform(
    marginal_posterior(mixed, "mu", prior_samples = TRUE), "exp"
  )))
})

test_that("posterior_transform() swaps supports and reverses densities of decreasing maps", {

  mixed <- .draws_metadata_mixed_for_test()
  sigma <- marginal_posterior(mixed, "sigma", prior_samples = TRUE)
  posterior_metadata(sigma, "posterior_density") <- posterior_density_attribute(
    x = c(.5, 1, 2), y = c(.2, .6, .1), method = "test", density_method = "precomputed",
    support = posterior_support_attribute(c(0, Inf))
  )
  dec <- posterior_transform(sigma, "lin", list(a = 1, b = -2))

  expect_equal(as.numeric(dec), 1 - 2 * as.numeric(sigma), tolerance = 0)
  expect_identical(posterior_metadata(dec, "support")$bounds, c(-Inf, 1))
  expect_true(posterior_atoms_free(dec))
  # the stored density: locations 1 - 2 x in increasing order, heights / 2
  density <- posterior_metadata(dec, "posterior_density")
  expect_equal(density$x, c(-3, -1, 0))
  expect_equal(density$y, c(.1, .6, .2) / 2)
  expect_identical(density$support$bounds, c(-Inf, 1))
  # prior density of 1 - 2 sigma, sigma ~ lognormal(0, 1)
  prior_density <- posterior_metadata(dec, "prior_density")
  for(value in c(-5, -1, .5)){
    expect_equal(
      exp(prior_density_ordinate(prior_density, value)$log_density),
      stats::dlnorm((1 - value) / 2) / 2,
      tolerance = 1e-12, info = value
    )
  }
  # linear weights follow the affine map
  expect_equal(posterior_metadata(dec, "linear_weights")[["sigma"]], -2)
  expect_identical(BayesTools:::.bt_meta_get(dec, "linear_offset"), 1)

  # the spike of 'mu' maps to 1 - 2 * 0
  mu <- posterior_transform(marginal_posterior(mixed, "mu"), "lin", list(a = 1, b = -2))
  expect_equal(as.numeric(posterior_metadata(mu, "atoms")$locations), 1)
  expect_identical(posterior_metadata(mu, "support")$bounds, c(-Inf, Inf))
})

test_that("posterior_transform() refuses transformations that are not monotone and invertible", {

  mixed <- .draws_metadata_mixed_for_test()
  mp    <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  square <- list(fun = function(x) x^2, inv = sqrt, jac = function(x) 2 * x)
  expect_error(posterior_transform(mp, square),
               class = "BayesTools_nonmonotone_transformation")
  expect_error(
    posterior_transform(mp, square),
    paste0(
      "The transformation of the posterior draws is unavailable: it must be ",
      "strictly monotone and invertible, but its derivative 'jac' is not finite ",
      "and nonzero with one sign on the draws."
    ),
    fixed = TRUE
  )
  expect_error(posterior_transform(mp, "lin", list(a = 1, b = 0)),
               class = "BayesTools_nonmonotone_transformation")
  expect_error(posterior_transform(mp, "exp_lin", list(a = 0, b = 0)),
               class = "BayesTools_transformation")
  # exp_lin is defined for nonnegative draws only
  expect_error(posterior_transform(mp, "exp_lin", list(a = 0, b = 2)),
               class = "BayesTools_transformation_domain")
  expect_error(posterior_transform(as.numeric(mp), "exp"),
               "not plain numeric draws", fixed = TRUE)

  # transformed mixed posteriors are not marginalized with their untransformed prior
  expect_error(
    marginal_posterior(posterior_transform(mixed, "exp"), "mu"),
    "Pass the transformation to 'marginal_posterior(transformation = )' instead.",
    fixed = TRUE
  )

  # the draws record the applied transformations in their own metadata, so
  # mixed posteriors without a column table are refused as well
  theta <- stats::rnorm(400, .3, .2)
  class(theta) <- c("mixed_posteriors", "mixed_posteriors.simple", class(theta))
  theta <- .bt_meta_set(theta, "draw_index", seq_along(theta))
  theta <- .bt_meta_set(theta, "model_probabilities", list(
    prior = .model_probability_pair(1, 0, "prior", "raw"),
    posterior = .model_probability_pair(1, 0, "posterior", "raw")))
  theta <- .bt_draws_set_component(theta, source = "model", component = rep(1, length(theta)))
  attr(theta, "parameter")  <- "theta"
  attr(theta, "prior_list") <- prior("normal", list(0, 1))
  hand_built <- list(theta = theta)
  class(hand_built) <- c("mixed_posteriors", "list")
  expect_null(posterior_metadata(hand_built$theta, "quantities"))
  transformed <- posterior_transform(hand_built, "exp")
  expect_identical(posterior_metadata(transformed$theta, "output_transformations"), "exp")
  expect_identical(
    posterior_metadata(posterior_transform(transformed, "lin", list(a = 0, b = 2))$theta,
                       "output_transformations"),
    c("exp", "lin")
  )
  expect_error(
    marginal_posterior(transformed, "theta", prior_samples = TRUE),
    "Pass the transformation to 'marginal_posterior(transformation = )' instead.",
    fixed = TRUE
  )
  expect_null(posterior_metadata(marginal_posterior(hand_built, "theta"), "output_transformations"))
  expect_identical(
    posterior_metadata(marginal_posterior(hand_built, "theta", transformation = "exp"),
                       "output_transformations"),
    "exp"
  )
  expect_error(
    posterior_metadata(theta, "output_transformations") <- "log",
    "Draw metadata 'output_transformations' is invalid",
    fixed = TRUE
  )
})

test_that("draw components index the declared component list", {

  # spike-and-slab: the slab is component 1 and the spike component 2 of the
  # prior (the fitted inclusion indicator is 1 for the slab)
  mixed <- .draws_metadata_mixed_for_test()
  expect_identical(BayesTools:::.bt_meta_get(mixed$mu, "component"), rep(2:1, c(100L, 300L)))
  expect_identical(BayesTools:::.bt_meta_get(mixed$mu, "component_source"), "spike_and_slab")
  expect_identical(
    attr(attr(mixed$mu, "prior_list"), "components")[c(2L, 1L)],
    c("null", "alternative")
  )
  # parameters without a mixture prior have one component by explicit rule
  expect_null(BayesTools:::.bt_meta_get(mixed$sigma, "component"))
  expect_identical(BayesTools:::.bt_draws_component(mixed$sigma), rep(1L, 400L))
  expect_null(BayesTools:::.bt_meta_get(mixed$mu, "draw_index"))

  # the per-component Savage-Dickey keys use the same component indices
  marginal <- marginal_posterior(mixed, "mu", prior_samples = TRUE)
  components <- BayesTools:::.bt_meta_get(marginal, "components")
  expect_setequal(components$keys[, "mu"], c(1, 2))
  expect_identical(
    as.integer(components$keys[components$index, "mu"]),
    BayesTools:::.bt_meta_get(mixed$mu, "component")
  )

  # mixtures keep the fitted component index
  set.seed(2)
  posterior <- cbind(
    theta           = c(rep(0, 20), stats::rnorm(40)),
    theta_indicator = rep(c(1, 2), c(20L, 40L))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(theta = prior_mixture(
    list(prior("spike", list(0)), prior("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  ))
  fit <- attach_test_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_parameter_map(fit, monitor_names = colnames(posterior))
  theta <- as_mixed_posteriors(fit, "theta")$theta
  expect_identical(BayesTools:::.bt_meta_get(theta, "component"), rep(1:2, c(20L, 40L)))
  expect_identical(BayesTools:::.bt_meta_get(theta, "component_source"), "mixture")
  atoms <- posterior_metadata(theta, "atoms")
  expect_equal(as.numeric(atoms$locations[, 1L]), 0)
  expect_equal(atoms$mass, 1 / 3)

  # the component helpers translate fitted indicators and reject others
  spike_and_slab <- prior_spike_and_slab(prior("normal", list(0, 1)))
  expect_identical(
    BayesTools:::.bt_component_from_indicator(spike_and_slab, c(0, 1, 1)),
    c(2L, 1L, 1L)
  )
  expect_error(
    BayesTools:::.bt_component_from_indicator(spike_and_slab, c(0, 2)),
    "Spike-and-slab indicator draws must be 0 or 1.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_meta_set(1:3, "component", c(0, 1, 1)),
    "Draw metadata 'component' is invalid",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_meta_set(1:3, "component_source", "models"),
    "Draw metadata 'component_source' is invalid",
    fixed = TRUE
  )
})

test_that("metadata build failures propagate instead of switching estimators", {

  mixed <- .draws_metadata_mixed_for_test()
  # the per-component Savage-Dickey estimator needs the component of every
  # mixture term: without it the marginal stops instead of silently pooling
  broken <- mixed
  broken$mu <- BayesTools:::.bt_meta_set(broken$mu, "component", NULL)
  expect_error(
    marginal_posterior(broken, "mu", prior_samples = TRUE),
    "Mixture component indices are unavailable.",
    fixed = TRUE
  )

  # malformed atom declarations stop instead of reading as undeclared
  atoms <- posterior_atom_attribute(data.frame(x = 0, mass = .25))
  atoms$mass <- 2
  expect_error(
    BayesTools:::.posterior_atoms_from_attribute(atoms),
    "Posterior atom masses cannot sum to more than one.",
    fixed = TRUE
  )
  expect_null(BayesTools:::.posterior_atoms_from_attribute(NULL))

  # columns outside the prior-density context have no support by rule
  context <- BayesTools:::.prior_density_build_context(
    prior_list   = list(mu = prior("normal", list(0, 1), list(0, Inf))),
    column_names = c("mu", "PET")
  )
  draws <- matrix(abs(stats::rnorm(20)), ncol = 2, dimnames = list(NULL, c("mu", "PET")))
  draws <- BayesTools:::.posterior_support_set_from_prior_context(draws, context)
  support <- BayesTools:::.bt_meta_get(draws, "support")
  expect_identical(names(support), "mu")
  expect_equal(support$mu$bounds, c(0, Inf))
})

test_that("metadata of draws whose values changed stop instead of describing the old values", {

  mixed <- .draws_metadata_mixed_for_test()
  mp <- marginal_posterior(mixed, "sigma", prior_samples = TRUE)
  reference <- Savage_Dickey_BF(mp, 1)
  message <- paste0(
    "The metadata of these posterior draws are unavailable: their values ",
    "changed after the metadata were attached (for example by 'x[] <- ', ",
    "'x[i] <- ', or 'pmin()'), so the supports, atoms, and prior densities ",
    "no longer describe them. Use 'posterior_transform()' (or ",
    "'marginal_posterior(transformation = )') for transformed posterior ",
    "distributions."
  )

  # replacing values keeps the attributes; the review's probe returned the
  # Bayes factor of the untransformed support and prior (0.160)
  replaced <- mp
  replaced[] <- 2 * mp
  element <- mp
  element[3] <- 10
  stale <- list(
    replaced = replaced,
    element  = element,
    pmin     = pmin(mp, 1),
    pmax     = pmax(mp, 1)
  )
  for(name in names(stale)){
    expect_error(Savage_Dickey_BF(stale[[name]], 2), message, fixed = TRUE, info = name)
    expect_error(Savage_Dickey_BF(stale[[name]], 2), class = "BayesTools_stale_metadata")
    expect_error(posterior_metadata(stale[[name]], "support"), class = "BayesTools_metadata")
    # the metadata cannot be updated without rebuilding the container
    expect_error(posterior_metadata(stale[[name]], "atoms") <- NULL,
                 class = "BayesTools_stale_metadata")
    expect_error(BayesTools:::.bt_meta_set(stale[[name]], "support", NULL),
                 class = "BayesTools_stale_metadata")
  }
  factor_draws <- BayesTools:::.bt_meta_set(
    matrix(stats::rnorm(20), ncol = 2, dimnames = list(NULL, c("a", "b"))),
    "atoms", posterior_atom_attribute()
  )
  factor_draws[, "b"] <- 0
  expect_error(posterior_metadata(factor_draws, "atoms"), class = "BayesTools_stale_metadata")

  # draws with their values pass, whatever kept or restored the attributes
  restored <- mp
  restored[] <- as.numeric(mp)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(mp, path)
  unchanged <- list(mp, restored, pmin(mp, Inf), readRDS(path))
  for(value in unchanged){
    expect_identical(Savage_Dickey_BF(value, 1), reference)
  }

  # missing values count: defining an undefined draw changes the draws
  undefined <- BayesTools:::.bt_meta_set(c(1, NA, 3), "undefined_draws", "correlation")
  expect_identical(posterior_metadata(undefined, "undefined_draws"), "correlation")
  undefined[2] <- 2
  expect_error(posterior_metadata(undefined, "undefined_draws"), class = "BayesTools_stale_metadata")

  # the sums agree up to the rounding bound of summing in double precision
  values <- stats::rnorm(1000)
  current <- BayesTools:::.bt_meta_fingerprint(values)
  rounded <- current$value
  rounded[["sum"]] <- rounded[["sum"]] + 8 * .Machine$double.eps * current$scale[["sum"]]
  expect_true(BayesTools:::.bt_meta_fingerprint_matches(rounded, current))
  shifted <- current$value
  shifted[["weighted_sum"]] <- shifted[["weighted_sum"]] + 1e-9 * current$scale[["weighted_sum"]]
  expect_false(BayesTools:::.bt_meta_fingerprint_matches(shifted, current))
  # a permutation keeps the sum but not the position-weighted sum
  expect_false(BayesTools:::.bt_meta_fingerprint_matches(
    current$value, BayesTools:::.bt_meta_fingerprint(rev(values))
  ))
})

test_that("subsetting a list of mixed posteriors keeps the list's metadata", {

  set.seed(1)
  data <- data.frame(x = stats::rnorm(50, 3, 2))
  formula_result <- JAGS_formula(
    ~ x, parameter = "mu", data = data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE)
  )
  n <- 400
  posterior <- cbind(mu_intercept = stats::rnorm(n), mu_x = stats::rnorm(n),
                     sigma = stats::rlnorm(n))
  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = n
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- c(formula_result$prior_list, list(sigma = prior("lognormal", list(0, 1))))
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)
  samples <- as_mixed_posteriors(
    fit, c("mu_intercept", "mu_x", "sigma"), transform_scaled = TRUE, n_prior_samples = 500
  )

  subset <- samples[c("mu_x", "sigma")]
  expect_identical(class(subset), class(samples))
  expect_identical(names(subset), c("mu_x", "sigma"))
  expect_identical(subset$mu_x, samples$mu_x)
  for(field in c("prior_context", "formula_scale", "transform_scaled")){
    expect_identical(BayesTools:::.bt_meta_get(subset, field),
                     BayesTools:::.bt_meta_get(samples, field), info = field)
  }
  expect_true(isTRUE(BayesTools:::.bt_meta_get(subset, "transform_scaled")))
  # the per-parameter prior densities follow the kept parameters
  densities <- BayesTools:::.bt_meta_get(samples, "prior_densities")
  expect_true("mu_intercept" %in% names(densities))
  kept_densities <- BayesTools:::.bt_meta_get(subset, "prior_densities")
  expect_identical(names(kept_densities), intersect(names(densities), c("mu_x", "sigma")))
  expect_identical(class(kept_densities), class(densities))
  for(key in names(kept_densities)){
    expect_identical(kept_densities[[key]], densities[[key]], info = key)
  }
  expect_identical(attr(subset, "prior_list"), attr(samples, "prior_list"))
  # the prior overlay of the subset uses the transformed prior of the full list
  plot_layers <- function(samples){
    set.seed(1)
    plot <- plot_posterior(samples, "mu_x", plot_type = "ggplot", prior = TRUE)
    lapply(seq_along(plot$layers), function(i) ggplot2::layer_data(plot, i))
  }
  expect_identical(plot_layers(subset), plot_layers(samples))
  # an empty index keeps every element
  expect_identical(samples[], samples)
})

# Forces the R evaluator of draw fingerprints (or, with 'r = FALSE', the
# native pass) until the calling test ends.
.local_r_fingerprint <- function(r = TRUE, envir = parent.frame()){
  private <- BayesTools:::.BayesTools_private
  old <- private$draw_fingerprint_r
  private$draw_fingerprint_r <- r
  withr::defer(private$draw_fingerprint_r <- old, envir = envir)
  invisible(private)
}

test_that("draw fingerprints are computed in one native pass in a fixed order", {

  .local_r_fingerprint(r = FALSE)

  fingerprint <- function(x) .Call("BayesTools_draw_fingerprint", x, PACKAGE = "BayesTools")
  set.seed(3)
  values <- c(stats::rnorm(1001), NA, NaN, stats::rnorm(6))
  native <- fingerprint(values)
  n <- length(values)
  observed <- !is.na(values)
  position <- which(observed)
  expect_identical(native[1:2], c(n, 2))
  # value i is added to partial sum i mod 4 in increasing i, in double
  # precision, and the partial sums are combined pairwise
  lane <- (seq_len(n) - 1L) %% 4L
  partial <- vapply(0:3, function(l){
    Reduce(`+`, values[observed & lane == l], accumulate = FALSE)
  }, numeric(1))
  expect_identical(native[[3L]], (partial[[1L]] + partial[[2L]]) + (partial[[3L]] + partial[[4L]]))
  absolute <- vapply(0:3, function(l){
    Reduce(`+`, abs(values[observed & lane == l]), accumulate = FALSE)
  }, numeric(1))
  expect_identical(native[[5L]], (absolute[[1L]] + absolute[[2L]]) + (absolute[[3L]] + absolute[[4L]]))
  # the position-weighted sum within its rounding bound
  k <- ceiling(n / 4) + 2
  gamma <- k * .Machine$double.eps / 2
  expect_lte(abs(native[[4L]] - sum(values[observed] * position)),
             2 * gamma * sum(abs(values[observed] * position)))
  expect_equal(native[[6L]], sum(abs(values[observed] * position)), tolerance = 1e-12)

  # integer and logical draws, whose sums are exact
  expect_identical(fingerprint(c(3L, NA, -2L, 7L, 1L)), c(5, 1, 9, 30, 13, 42))
  expect_identical(fingerprint(c(TRUE, FALSE, NA, TRUE)), c(4, 1, 2, 5, 2, 5))
  expect_identical(fingerprint(numeric()), c(0, 0, 0, 0, 0, 0))
  expect_error(fingerprint("a"), "Draw fingerprints require numeric or logical values.", fixed = TRUE)

  # the container records the first four
  x <- .bt_meta_set(values, "undefined_draws", "correlation")
  expect_identical(unname(attr(x, "bayestools_meta")$fingerprint), native[1:4])
})

test_that("the R evaluator of draw fingerprints agrees with the native pass", {

  native_fingerprint <- function(x) .Call("BayesTools_draw_fingerprint", x, PACKAGE = "BayesTools")
  as_fingerprint <- function(out) list(
    value = c(length = out[[1L]], missing = out[[2L]], sum = out[[3L]], weighted_sum = out[[4L]]),
    scale = c(sum = out[[5L]], weighted_sum = out[[6L]])
  )
  set.seed(1)
  cases <- list(
    normal_1e6    = stats::rnorm(1e6),
    odd_length    = stats::rnorm(1e4 + 3),
    n1            = 0.3,
    n3            = c(1.5, -2, 3),
    n5            = stats::rnorm(5),
    with_na_nan   = c(stats::rnorm(100), NA, NaN, stats::rnorm(3)),
    with_inf      = c(stats::rnorm(10), Inf),
    inf_minus_inf = c(Inf, 1, -Inf),
    cancel        = c(1e16, 1, -1e16, 3.5, 2^-40),
    tiny          = stats::rnorm(1e4) * 1e-300,
    huge          = stats::rnorm(1e4) * 1e300,
    weighted_inf  = c(stats::rnorm(10), 1e308, 1e308 / 2),
    integer       = sample.int(100L, 1e5, TRUE),
    integer_na    = c(1L, NA, 3L),
    logical       = c(TRUE, NA, FALSE, TRUE),
    compact_int   = 1:100000,
    compact_real  = as.numeric(1:100000),
    matrix        = matrix(stats::rnorm(1e4), 100),
    all_na        = rep(NA_real_, 7),
    empty         = numeric()
  )
  # identical where the arithmetic is exact (integer-valued sums, lanes of at
  # most one value, no observed values) or where R's long double is double
  exact <- c("n1", "n3", "integer", "integer_na", "logical", "compact_int",
             "compact_real", "all_na", "empty")
  long_double_is_double <- is.null(.Machine$longdouble.digits) ||
    .Machine$longdouble.digits <= .Machine$double.digits
  for(name in names(cases)){
    native <- native_fingerprint(cases[[name]])
    r <- .bt_meta_fingerprint_r(cases[[name]])
    expect_identical(r[1:2], native[1:2], info = name)
    if(name %in% exact || long_double_is_double){
      expect_identical(r, native, info = name)
    }
    # a fingerprint stored by either evaluator matches the other one's
    expect_true(.bt_meta_fingerprint_matches(as_fingerprint(native)$value, as_fingerprint(r)), info = name)
    expect_true(.bt_meta_fingerprint_matches(as_fingerprint(r)$value, as_fingerprint(native)), info = name)
  }
  expect_error(.bt_meta_fingerprint_r("a"), "Draw fingerprints require numeric or logical values.", fixed = TRUE)
})

test_that("draw metadata work with the R evaluator of the fingerprint", {

  set.seed(2)
  native_draws <- .bt_meta_set(stats::rnorm(1e4), "atoms", posterior_atom_attribute())
  native_mixed <- .draws_metadata_mixed_for_test()
  reference <- Savage_Dickey_BF(marginal_posterior(native_mixed, "sigma", prior_samples = TRUE), 1)

  .local_r_fingerprint()
  expect_false(.bt_meta_fingerprint_native())
  draws <- stats::rnorm(100)
  posterior_metadata(draws, "atoms") <- posterior_atom_attribute()
  expect_s3_class(posterior_metadata(draws, "atoms"), "BayesTools_posterior_atoms")
  expect_identical(unname(attr(draws, "bayestools_meta")$fingerprint), .bt_meta_fingerprint_r(draws)[1:4])

  # metadata stored by the native pass are read with the R evaluator
  expect_s3_class(posterior_metadata(native_draws, "atoms"), "BayesTools_posterior_atoms")
  mixed <- .draws_metadata_mixed_for_test()
  mp <- marginal_posterior(mixed, "sigma", prior_samples = TRUE)
  expect_identical(Savage_Dickey_BF(mp, 1), reference)

  # and changed values stop
  replaced <- mp
  replaced[] <- 2 * mp
  element <- mp
  element[3] <- 10
  for(stale in list(replaced, element, pmin(mp, 1))){
    expect_error(Savage_Dickey_BF(stale, 2), class = "BayesTools_stale_metadata")
    expect_error(posterior_metadata(stale, "atoms") <- NULL, class = "BayesTools_stale_metadata")
  }
})

test_that("draw fingerprints use the R evaluator when the native routines are not loaded", {

  testthat::local_mocked_bindings(
    .BayesTools_native_routines_loaded = function(pkgname = "BayesTools") FALSE,
    .package = "BayesTools"
  )
  expect_false(.bt_meta_fingerprint_native())
  set.seed(3)
  draws <- stats::rnorm(100)
  posterior_metadata(draws, "support") <- posterior_support_attribute(c(-Inf, Inf))
  expect_identical(unname(attr(draws, "bayestools_meta")$fingerprint), .bt_meta_fingerprint_r(draws)[1:4])
  draws[1] <- draws[1] + 1
  expect_error(posterior_metadata(draws, "support"), class = "BayesTools_stale_metadata")
})

test_that("setting several draw-metadata fields checks the draws once", {

  checks <- 0L
  original <- BayesTools:::.bt_meta_check_current
  testthat::local_mocked_bindings(
    .bt_meta_check_current = function(meta, current){
      checks <<- checks + 1L
      original(meta, current)
    },
    .package = "BayesTools"
  )
  x <- .bt_meta_set(stats::rnorm(100), "atoms", posterior_atom_attribute())
  checks <- 0L
  x <- .bt_meta_update(
    x,
    support = posterior_support_attribute(c(-Inf, Inf)),
    undefined_draws = "correlation",
    atoms = NULL
  )
  expect_identical(checks, 1L)
  expect_identical(names(attr(x, "bayestools_meta")), c("support", "undefined_draws", "fingerprint"))
  expect_identical(posterior_metadata(x, "undefined_draws"), "correlation")
})

test_that("reading several draw-metadata fields fingerprints the draws once", {

  passes <- 0L
  original <- BayesTools:::.bt_meta_fingerprint
  testthat::local_mocked_bindings(
    .bt_meta_fingerprint = function(x){
      passes <<- passes + 1L
      original(x)
    },
    .package = "BayesTools"
  )
  support <- posterior_support_attribute(c(-Inf, Inf))
  condition <- list(conditional = "mu", conditional_rule = "AND", condition_key = "mu")
  set.seed(4)
  x <- .bt_meta_update(
    stats::rnorm(100),
    support       = support,
    prior_density = prior("normal", list(0, 1)),
    condition     = condition
  )

  passes <- 0L
  fields <- .bt_meta_get_fields(x, c("support", "atoms", "prior_density"))
  expect_identical(passes, 1L)
  expect_identical(fields, list(support = support, atoms = NULL, prior_density = prior("normal", list(0, 1))))

  # the continuous part of draws for Savage-Dickey Bayes factors carries their
  # support and prior density: one read and one update (and the components)
  passes <- 0L
  continuous <- .Savage_Dickey_BF.continuous_posterior(x, posterior_atom_attribute())
  expect_identical(passes, 3L)
  expect_identical(.bt_meta_get(continuous$samples, "support"), support)
  expect_identical(.bt_meta_get(continuous$samples, "prior_density"), prior("normal", list(0, 1)))

  # the conditioning metadata of the draws of a marginal posterior are read once
  passes <- 0L
  metadata <- .marginal_posterior_condition_metadata(list(mu = x), condition_source = x)
  expect_identical(passes, 1L)
  expect_identical(metadata, c(condition, list(condition_event = NULL)))

  # and so are the posterior density and ordinate sources of draws
  passes <- 0L
  expect_identical(.posterior_density_sources(x), list())
  expect_identical(.posterior_ordinate_sources(x), list())
  expect_identical(passes, 2L)

  x[] <- 2 * x
  expect_error(.bt_meta_get_fields(x, "support"), class = "BayesTools_stale_metadata")
  expect_error(.marginal_posterior_condition_metadata(list(mu = x), condition_source = x),
               class = "BayesTools_stale_metadata")
})
