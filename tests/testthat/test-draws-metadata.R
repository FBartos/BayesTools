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
  expect_identical(names(meta), c("atoms", "support"))
  expect_s3_class(posterior_metadata(x, "atoms"), "BayesTools_posterior_atoms")

  # removing the last field removes the container
  posterior_metadata(x, "atoms") <- NULL
  posterior_metadata(x, "support") <- NULL
  expect_null(attributes(x))

  # values are validated when they are set
  expect_error(
    posterior_metadata(x, "atoms") <- list(declared = TRUE),
    "Draw metadata 'atoms' is invalid: it must be created with 'posterior_atom_attribute()'.",
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
    density_method = "precomputed",
    point_masses   = data.frame(x = 0, mass = .2)
  )
  class(x) <- c("marginal_posterior.simple", "marginal_posterior")
  # the stored density's point masses do not declare the posterior atoms
  expect_null(BayesTools:::.posterior_atoms_get(x))
  expect_error(Savage_Dickey_BF(x), "Posterior atom status is unknown", fixed = TRUE)

  posterior_metadata(x, "atoms") <- posterior_atom_attribute(data.frame(x = 0, mass = .2))
  expect_error(
    Savage_Dickey_BF(x),
    "declared point mass at the exact null",
    fixed = TRUE
  )
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
