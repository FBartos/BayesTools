skip_if_not_test_profile("unit")

test_that("missing ordered recipes cannot default-bind former formula owners", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2L)), x = 1:6)
  info <- JAGS_formula(~f+x+f:x, "mu", data, list(intercept = prior("point", list(0)), x = prior("normal", list(0, 1)),
    f = prior_ordered(prior("normal", list(0, 1)), id = c(f = "shape")),
    "f:x" = prior_ordered(prior("normal", list(0, 1)), id = c(f = "shape"))))
  original <- info$prior_list
  draws <- .generate_prior_sample_matrix(original, 8L, seed = 17)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), original,
    formula_design = list(mu = info$formula_design))
  stripped <- original
  for(parameter in names(stripped)[vapply(stripped, is.prior.ordered, logical(1))]) attr(stripped[[parameter]], "ordered_metadata") <- NULL
  calls <- 0L; leaf_rng <- rng
  testthat::local_mocked_bindings(rng = function(...){calls <<- calls + 1L; leaf_rng(...)}, .package = "BayesTools")
  for(priors in list(stripped["mu_f"], stripped)){
    set.seed(113); before <- .Random.seed
    condition <- expect_error(.generate_prior_sample_matrix(priors, 4L, seed = 17), class = "BayesTools_ordered_metadata_unavailable")
    expect_s3_class(condition, "BayesTools_refit_required")
    expect_null(conditionCall(condition))
    expect_identical(.Random.seed, before)
  }
  for(design in list(list(mu = info$formula_design), NULL, "absent")){
    for(scale in list(NULL, list())){
      bad <- fit; attr(bad, "prior_list") <- stripped
      attr(bad, "formula_design") <- design; attr(bad, "formula_scale") <- scale
      set.seed(113); before <- .Random.seed
      condition <- expect_error(transform_prior_samples(bad, 4L, seed = 17, formula_scale = scale), class = "BayesTools_ordered_metadata_unavailable")
      expect_null(conditionCall(condition))
      expect_identical(.Random.seed, before)
    }
  }
  expect_identical(calls, 0L)
  expect_type(transform_prior_samples(fit, 4L, seed = 17, formula_scale = list()), "double")
  unbound <- prior_ordered(prior("point", list(10)), c(.25, .75)); attr(unbound, "levels") <- 3L
  expect_equal(unname(.generate_prior_sample_matrix(list(f = unbound), 4L, seed = 17)[, c("f[1]", "f[2]")]), matrix(rep(c(2.5, 7.5), each = 4L), 4L), tolerance = 0)
  container <- prior_mixture(list(prior_factor_levels(unbound, c("lo", "mid", "hi")), prior_factor_levels(unbound, c("lo", "mid", "hi"))))
  expect_type(.generate_prior_sample_matrix(list(f = container), 4L, seed = 17), "double")
})

test_that("ordered persisted recipes match declarations and their actual owner", {
  fixture <- ordered_plot_test_fixture(prior("point", list(10)), c(.25, .75),
    levels = c("lo", "mid", "hi"))
  original <- fixture$prior
  valid <- function(p, owner = "mu_f"){
    if("parameter" %in% names(formals(.bt_ordered_metadata_valid))) .bt_ordered_metadata_valid(p, owner) else .bt_ordered_metadata_valid(p)
  }
  expect_true(valid(original))
  corruptions <- list(
    duplicate_grid = function(m) {m$coefficient_grid$f <- c(1L, 1L); m},
    permuted_grid = function(m) {m$coefficient_grid$f <- c(2L, 1L); m},
    out_of_range = function(m) {m$coefficient_grid$f[2L] <- 3L; m},
    slice = function(m) {m$slice_index[2L] <- 2L; m},
    coverage = function(m) {m$allocations <- list(); m},
    node = function(m) {m$allocations[[1L]]$node <- NULL; m},
    id = function(m) {m$allocations[[1L]]$id <- "other"; m},
    contrast = function(m) {m$allocations[[1L]]$contrast <- NULL; m},
    key = function(m) {names(m$allocations) <- "other"; m},
    dimension = function(m) {m$allocations[[1L]]$dim <- 3L; m},
    allocation = function(m) {m$allocations[[1L]]$spec$weights <- c(.5, .5); m},
    provenance = function(m) {m$numeric_literals$total <- NULL; m},
    partition = function(m) {m$ordinary_terms <- "f"; m},
    owner = function(m) {m$parameter_name <- "other"; m})
  for(mutate in corruptions){
    bad <- original
    attr(bad, "ordered_metadata") <- mutate(attr(bad, "ordered_metadata"))
    expect_false(valid(bad))
    priors <- list(mu_f = bad)
    for(action in list(function() JAGS_add_priors("model{}", priors),
      function() JAGS_to_monitor(priors), function() JAGS_get_inits(priors, chains = 1L, seed = 17),
      function() .generate_prior_sample_matrix(priors, 4L, seed = 17))){
      condition <- tryCatch({action(); NULL}, error = identity)
      expect_s3_class(condition, "BayesTools_ordered_metadata_unavailable")
      expect_s3_class(condition, "BayesTools_refit_required")
      if(inherits(condition, "condition")) expect_null(conditionCall(condition))
    }
    fit <- fixture$fit
    attr(fit, "prior_list")$mu_f <- bad
    expect_s3_class(tryCatch({JAGS_ordered_parameter_spec(fit); NULL}, error = identity), "BayesTools_ordered_metadata_unavailable")
    expect_s3_class(tryCatch({JAGS_evaluate_deterministic(fit, fixture$draws, nodes = "mu_f"); NULL}, error = identity),
      "BayesTools_ordered_metadata_unavailable")
  }
  for(field in c("factor_terms", "factor_contrasts", "level_names", "factor_design")){
    bad <- original; attr(bad, field) <- NULL
    expect_false(valid(bad))
  }
  bad <- original; attr(bad, "coefficient_dim") <- 3L
  expect_false(valid(bad))
  expect_false(valid(original, "theta_f"))
  expect_error(JAGS_to_monitor(list(theta_f = original)), class = "BayesTools_ordered_metadata_unavailable")
  data <- data.frame(f = ordered(c("lo", "mid", "hi")))
  bind <- function(p) JAGS_formula(~f, "theta", data, list(intercept = prior("point", list(0)), f = p))$prior_list$theta_f
  rebound <- bind(original)
  expect_identical(rebound, bind(prior_ordered(prior("point", list(10)), c(.25, .75))))
  expect_identical(fixture$prior, original)
})

test_that("fixed ordered replay retains named empty matrix dimensions", {
  original <- JAGS_formula(~f, "mu", data.frame(f = ordered(c("lo", "mid", "hi"))),
    list(intercept = prior("point", list(0)), f = prior_ordered(prior("point", list(10)), c(.25, .75))))$prior_list$mu_f
  spec <- .bt_ordered_spec("mu_f", original)
  for(n in c(0L, 2L)){
    draws <- matrix(numeric(n), n, 1L, dimnames = list(NULL, "unused"))
    out <- .bt_deterministic_node_evaluate(.bt_dnode_ordered_coefficients(spec), .bt_deterministic_lookup(draws))
    expect_identical(dim(out), c(n, 2L))
    expect_identical(colnames(out), spec$coefficient_names)
    expect_equal(unname(out), matrix(rep(c(2.5, 7.5), each = n), n, 2L), tolerance = 0)
  }
})

test_that("ordered roots share only compatible emitted allocation owners", {
  data <- data.frame(shape = ordered(c("lo", "mid", "hi")), ordered_total = ordered(c("lo", "mid", "hi")))
  first <- JAGS_formula(~shape, "ordered_alloc", data, list(intercept = prior("point", list(0)),
    shape = prior_ordered(prior("normal", list(0, 1)), c(.4, .6))))$prior_list
  second <- JAGS_formula(~ordered_total, "mu", data, list(intercept = prior("point", list(0)),
    ordered_total = prior_ordered(prior("normal", list(0, 1)), id = "shape")))$prior_list
  collision <- c(first, second)
  for(action in list(function() JAGS_add_priors("model{}", collision), function() JAGS_to_monitor(collision),
    function() JAGS_get_inits(collision, chains = 1L, seed = 17),
    function() .generate_prior_sample_matrix(collision, 4L, seed = 17))){
    expect_error(action(), "ordered_alloc_shape_ordered_total", fixed = TRUE)
  }
  ordered <- second["mu_ordered_total"]
  ordinary <- list(prior("dirichlet", list(c(2, 3))))
  names(ordinary) <- "ordered_alloc_shape_ordered_total"
  expect_error(.generate_prior_sample_matrix(c(ordered, ordinary), 4L, seed = 17),
    "ordered_alloc_shape_ordered_total", fixed = TRUE)
  expect_type(JAGS_add_priors("model{}", first), "character")
  # A fixed split emits neither its normalized allocation nor Gamma roots.
  fixed <- first
  fixed$ordered_alloc_shape_ordered_alloc_shape_1 <- prior("normal", list(0, 1))
  fixed$prior_par_eta_ordered_alloc_shape_ordered_alloc_shape_1 <- prior("normal", list(0, 1))
  expect_type(JAGS_add_priors("model{}", fixed), "character")
  data$q <- data$shape
  auxiliary <- JAGS_formula(~q, "prior_par_s", data, list(intercept = prior("point", list(0)),
    q = prior_ordered(prior("normal", list(0, 1)), c(.4, .6))))$prior_list
  auxiliary$q_ordered_total <- prior("mt", list(location = 0, scale = 1, df = 3, K = 2))
  expect_error(JAGS_add_priors("model{}", auxiliary), "prior_par_s_q_ordered_total", fixed = TRUE)
})

test_that("ordered named maps require unique nonempty non-NA names", {
  total <- prior("point", list(10))
  for(map_names in list(c("f", "f"), c("f", NA_character_), c("f", ""))){
    expect_error(prior_ordered(total, id = setNames(c("a", "b"), map_names)), "unique, nonempty, non-NA", fixed = TRUE)
    expect_error(prior_ordered(total, allocation = setNames(list(c(.25, .75), c(.5, .5)), map_names)), "unique, nonempty, non-NA", fixed = TRUE)
  }
})

test_that("direct ordered container compiler refusals are explicit", {
  total <- prior("point", list(10))
  p <- prior_mixture(list(prior_factor_levels(prior_ordered(total, c(.25, .75)), c("lo", "mid", "hi")),
    prior_factor_levels(prior_ordered(total, c(.5, .5)), c("lo", "mid", "hi"))))
  for(action in list(function() JAGS_add_priors("model{}", list(mu_f = p)),
    function() JAGS_to_monitor(list(mu_f = p)), function() JAGS_get_inits(list(mu_f = p), chains = 1L, seed = 17))){
    condition <- expect_error(action(), class = "BayesTools_ordered_unavailable")
    expect_identical(conditionMessage(condition), "Mixtures of ordered prior containers cannot be bound to JAGS formulas. Put mixture behavior on 'total' instead.")
    expect_null(conditionCall(condition))
  }
})

test_that("joint ordered sharing uses first-owner leaf draws and preserves continued streams", {
  data <- data.frame(f = ordered(c("lo", "mid", "hi")))
  alpha <- c(.7, 3)
  make <- function(parameter, id) JAGS_formula(~f, parameter, data, list(intercept = prior("point", list(0)),
    f = prior_ordered(prior("normal", list(2, .7)), prior("dirichlet", list(alpha)), id = id)))$prior_list[[paste0(parameter, "_f")]]
  priors <- list(one_f = make("one", "shape"), two_f = make("two", "shape"), after = prior("normal", list(-1, 2)))
  n <- 32L
  set.seed(413)
  first_total <- stats::rnorm(n, 2, .7)
  first_share <- rng(prior("dirichlet", list(alpha)), n)
  second_total <- stats::rnorm(n, 2, .7)
  discarded_share <- rng(prior("dirichlet", list(alpha)), n)
  after <- stats::rnorm(n, -1, 2)
  joint <- .generate_prior_sample_matrix(priors, n, seed = 413)
  expect_identical(unname(joint[, c("one_f[1]", "one_f[2]")]), unname(first_total * first_share))
  expect_identical(unname(joint[, c("two_f[1]", "two_f[2]")]), unname(second_total * first_share))
  expect_identical(unname(joint[, "after"]), after)
  expect_false(identical(first_share, discarded_share))
  expect_identical(unname(joint[, c("ordered_alloc_shape_f[1]", "ordered_alloc_shape_f[2]")]), unname(first_share))
  unshared <- priors; unshared$two_f <- make("two", NULL)
  control <- .generate_prior_sample_matrix(unshared, n, seed = 413)
  expect_identical(control[, "after"], joint[, "after"])
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(joint[1:2, , drop = FALSE])), priors)
  public <- transform_prior_samples(fit, n_samples = n, seed = 413, formula_scale = list())
  columns <- c("one_f[1]", "one_f[2]", "two_f[1]", "two_f[2]", "after")
  expect_identical(public[, columns], joint[, columns])
  expect_identical(unname(control[, c("two_f[1]", "two_f[2]")]), unname(second_total * discarded_share))
})

test_that("ordered IDs reject malformed vectors and keep reusable factor maps", {
  total <- prior("point", list(1))
  for(id in list(character(), setNames(character(), character()))){
    expect_error(prior_ordered(total, id = id),
      "The 'id' argument must be NULL or contain at least one identifier.", fixed = TRUE)
  }
  expect_error(prior_ordered(total, id = c("a", "b")),
    "The unnamed 'id' argument must contain exactly one identifier.", fixed = TRUE)
  for(id in list(NULL, "shape", c(f = "shape"), c(f = "shape", g = "other"))){
    expect_s3_class(prior_ordered(total, id = id), "prior.ordered")
  }
  data <- data.frame(f = ordered(c("lo", "mid", "hi")))
  expect_type(JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(total, allocation = list(f = c(.4, .6), g = c(.2, .8)), id = c(g = "other")))), "list")
  expect_error(JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(total, allocation = list(g = c(.4, .6))))), "f", fixed = TRUE)
})

test_that("ordered indicator-coded terms refuse even zero totals", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2)), x = 1:6)
  for(total in c(0, 10)){
    condition <- expect_error(JAGS_formula(~f:x, "mu", data,
      list(intercept = prior("point", list(0)),
        "f:x" = prior_ordered(prior("point", list(total)), allocation = c(.2, .3, .5),
          contrast = "cumulative_levels"))), "codes 'f' by level indicators", fixed = TRUE,
      class = "BayesTools_ordered_unavailable")
    expect_identical(conditionMessage(condition), paste0(
      "The 'cumulative_levels' prior of the factor term 'f:x' is unavailable: the formula has no term 'x', ",
      "so 'f:x' codes 'f' by level indicators and has one coefficient per level instead of 'cumulative_levels' ",
      "contrast coefficients. Add 'x' to the formula to keep the 'cumulative_levels' contrast, or use ",
      "prior_factor(contrast = \"independent\") for one independent coefficient per level."))
    expect_null(conditionCall(condition))
    expect_false(inherits(condition, "BayesTools_ordered_coordinates_unavailable"))
  }
  data <- expand.grid(f = ordered(c("lo", "mid", "hi")), g = ordered(c("a", "b", "c")))
  p <- prior_ordered(prior("point", list(10)), allocation = list(f = c(.2, .3, .5), g = c(.1, .4, .5)),
    contrast = "cumulative_levels")
  condition <- expect_error(JAGS_formula(~f:g, "mu", data, list(intercept = prior("point", list(0)), "f:g" = p)),
    "codes 'f' and 'g' by level indicators", fixed = TRUE, class = "BayesTools_ordered_unavailable")
  expect_identical(conditionMessage(condition), paste0(
    "The 'cumulative_levels' prior of the factor term 'f:g' is unavailable: the formula has no term 'g' or 'f', ",
    "so 'f:g' codes 'f' and 'g' by level indicators and has one coefficient per level instead of 'cumulative_levels' ",
    "contrast coefficients. Add 'g' and 'f' to the formula to keep the 'cumulative_levels' contrast, or use ",
    "prior_factor(contrast = \"independent\") for one independent coefficient per level."))
  expect_null(conditionCall(condition))
  expect_false(inherits(condition, "BayesTools_ordered_coordinates_unavailable"))
  expect_type(JAGS_formula(~f+g+f:g, "mu", data, list(intercept = prior("point", list(0)), f = p, g = p, "f:g" = p)), "list")
})

test_that("fixed ordered replay declares only stochastic source dependencies", {
  data <- data.frame(f = ordered(c("lo", "mid", "hi")))
  info <- JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(prior("point", list(10)), allocation = c(.25, .75))))
  spec <- .bt_ordered_spec("mu_f", info$prior_list$mu_f)
  node <- .bt_dnode_ordered_coefficients(spec)
  expect_identical(node$dependencies, character())
  expect_length(spec$allocations[[1L]]$coordinates, 2L)
  empty <- matrix(numeric(), 2L, 0L, dimnames = list(NULL, character()))
  expect_equal(.bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(empty)),
    matrix(rep(c(2.5, 7.5), each = 2L), 2L, dimnames = list(NULL, spec$coefficient_names)), tolerance = 0)
  random <- JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(prior("normal", list(0, 1)))))
  random_spec <- .bt_ordered_spec("mu_f", random$prior_list$mu_f)
  expect_identical(.bt_dnode_ordered_coefficients(random_spec)$dependencies,
    c(random_spec$total_names, random_spec$allocations[[1L]]$coordinates))
})

test_that("ordered mixture replay keeps ordered ownership and conditions", {
  fixture <- ordered_plot_test_fixture(prior_mixture(list(prior("point", list(2)), prior("point", list(4)))),
    allocation = prior("dirichlet", list(c(2, 2, 2))))
  nodes <- JAGS_deterministic_nodes(fixture$fit)
  expect_identical(nodes$parameter[nodes$node == "mu_f_ordered_total"], "mu_f")
  for(value in c(3, 1.5, NA_real_, Inf)){
    bad <- fixture$draws; bad[, "mu_f_ordered_total_indicator"] <- value
    for(node in c("mu_f_ordered_total", "mu_f")){
      condition <- expect_error(JAGS_evaluate_deterministic(fixture$fit, bad, nodes = node),
        "Ordered total indicator 'mu_f_ordered_total_indicator' does not select a declared component.",
        fixed = TRUE, class = "BayesTools_ordered_invalid_state")
      expect_null(conditionCall(condition))
    }
  }
  missing <- fixture$draws[, setdiff(colnames(fixture$draws), c("mu_f_ordered_total_indicator", "mu_f_ordered_total")), drop = FALSE]
  condition <- expect_error(JAGS_evaluate_deterministic(fixture$fit, missing, nodes = "mu_f_ordered_total"),
    "Ordered deterministic node 'mu_f_ordered_total' is unavailable from 'draws'. Include its declared primitive source coordinates.",
    fixed = TRUE, class = "BayesTools_ordered_coordinates_unavailable")
  expect_null(conditionCall(condition))
  ordinary <- .bt_dnode_prior_mixture("theta", prior_mixture(list(prior("point", list(2)), prior("point", list(4)))))
  condition <- expect_error(.bt_deterministic_node_evaluate(ordinary,
    .bt_deterministic_lookup(matrix(3, 2L, 1L, dimnames = list(NULL, "theta_indicator")))),
    "Mixture indicator draws of 'theta' must index a mixture component.", fixed = TRUE)
  expect_identical(class(condition), c("simpleError", "error", "condition"))
})

test_that("expression-point ordered totals retain multi-slice snapshot monitors", {
  data <- expand.grid(f = ordered(c("lo", "mid", "hi")), g = factor(c("a", "b", "c")))
  info <- JAGS_formula(~f*g, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(prior("point", list(0)), allocation = c(.2, .8)),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"),
    "f:g" = prior_ordered(prior("point", list(location = expression(sigma))), allocation = c(.2, .8))))
  expect_true("mu_f__xXx__g_ordered_total" %in% JAGS_to_monitor(info$prior_list))
  spec <- .bt_ordered_spec("mu_f__xXx__g", info$prior_list$mu_f__xXx__g)
  snapshots <- matrix(c(2, 3, 4, 5), 2L, dimnames = list(NULL, spec$total_names))
  expect_identical(.bt_ordered_total_values(spec, .bt_deterministic_lookup(snapshots)), snapshots)
  expect_error(JAGS_ordered_density_kernel(info$prior_list["mu_f__xXx__g"]), class = "BayesTools_ordered_expression_unavailable")
  literal <- info$prior_list$mu_f__xXx__g
  literal$total <- prior("point", list(0))
  literal <- .bt_bind_ordered_prior_metadata(literal, "mu_f__xXx__g")
  expect_false("mu_f__xXx__g_ordered_total" %in% .JAGS_monitor.ordered(literal, "mu_f__xXx__g"))
})

test_that("ordered generated roots cannot alias ordinary declarations", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2L)), x = 1:6,
    f_ordered_alloc_f_1 = 1:6, f_ordered_total_indicator = 1:6)
  bind <- function(formula, p, ordinary) JAGS_formula(formula, "mu", data,
    c(list(intercept = prior("point", list(0)), f = p), ordinary))
  expect_error(bind(~f+f_ordered_alloc_f_1, prior_ordered(prior("point", list(1))),
    list(f_ordered_alloc_f_1 = prior("normal", list(0, 1)))), "mu_f_ordered_alloc_f_1", fixed = TRUE)
  spike <- prior_ordered(prior_spike_and_slab(prior("normal", list(0, 1)), prior("bernoulli", list(.5))))
  expect_error(bind(~f+f_ordered_total_indicator, spike,
    list(f_ordered_total_indicator = prior("normal", list(0, 1)))), "mu_f_ordered_total_indicator", fixed = TRUE)
  base <- bind(~f, spike, list())$prior_list
  for(root in c("mu_f_ordered_total_variable", "mu_f_ordered_total_inclusion", "prior_par_eta_mu_f_ordered_alloc_f_1")){
    expect_error(JAGS_add_priors("model{}", c(base, setNames(list(prior("normal", list(0, 1))), root))), root, fixed = TRUE)
  }
  mixture <- bind(~f, prior_ordered(prior_mixture(list(prior("normal", list(0, 1)), prior("point", list(0))))), list())$prior_list
  expect_error(JAGS_add_priors("model{}", c(mixture, list(mu_f_ordered_total_component_1 = prior("normal", list(0, 1))))),
    "mu_f_ordered_total_component_1", fixed = TRUE)
  expect_type(bind(~f+f_ordered_alloc_f_1,
    prior_ordered(prior("point", list(1)), allocation = c(.4, .6)),
    list(f_ordered_alloc_f_1 = prior("normal", list(0, 1)))), "list")
})

test_that("joint ordered prior draws honor shared allocation IDs without shifting other streams", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2L)), x = 1:6)
  make <- function(shared) {
    id <- if(shared) "shape" else NULL
    info <- JAGS_formula(~f+x+f:x, "mu", data, list(intercept = prior("point", list(0)),
      x = prior("normal", list(0, 1)), f = prior_ordered(prior("point", list(1)), id = id),
      "f:x" = prior_ordered(prior("point", list(1)), id = id)))
    columns <- c("mu_intercept", "mu_x", "mu_f[1]", "mu_f[2]", "mu_f__xXx__x[1]", "mu_f__xXx__x[2]",
      "mu_f_ordered_total", "mu_f__xXx__x_ordered_total")
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(matrix(0, 2L, length(columns),
      dimnames = list(NULL, columns)))), c(info$prior_list, list(after = prior("normal", list(0, 1)))),
      formula_design = list(mu = info$formula_design))
    list(fit = fit, priors = attr(fit, "prior_list"))
  }
  shared <- make(TRUE); unshared <- make(FALSE)
  set.seed(716); previous_seed <- .Random.seed
  joint <- transform_prior_samples(shared$fit, n_samples = 16L, seed = 17, formula_scale = list())
  expect_identical(.Random.seed, previous_seed)
  expect_identical(unname(joint[, c("mu_f[1]", "mu_f[2]")]), unname(joint[, c("mu_f__xXx__x[1]", "mu_f__xXx__x[2]")]))
  control <- transform_prior_samples(unshared$fit, n_samples = 16L, seed = 17, formula_scale = list())
  expect_false(identical(unname(control[, c("mu_f[1]", "mu_f[2]")]), unname(control[, c("mu_f__xXx__x[1]", "mu_f__xXx__x[2]")])))
  shared_raw <- .generate_prior_sample_matrix(shared$priors, 16L, seed = 17)
  unshared_raw <- .generate_prior_sample_matrix(unshared$priors, 16L, seed = 17)
  expect_identical(shared_raw[, c("mu_x", "after")], unshared_raw[, c("mu_x", "after")])
  expect_false(anyDuplicated(colnames(shared_raw)) > 0L)
  set.seed(17)
  first <- rng(shared$priors$mu_f, 16L, transform_factor_samples = FALSE)
  second <- rng(shared$priors$mu_f__xXx__x, 16L, transform_factor_samples = FALSE)
  expect_false(identical(first, second))
})

test_that("unbound ordered joint helpers retain their standalone streams and columns", {
  p <- prior_ordered(prior("normal", list(0, 1)), id = "shape")
  attr(p, "levels") <- 3L
  ordinary <- prior("normal", list(0, 1))
  actual <- .generate_prior_sample_matrix(list(a = p, b = p, ordinary = ordinary), 16L, seed = 17)
  set.seed(17)
  first <- rng(p, 16L, transform_factor_samples = FALSE)
  second <- rng(p, 16L, transform_factor_samples = FALSE)
  later <- rng(ordinary, 16L)
  expect_identical(as.numeric(actual[, c("a[1]", "a[2]")]), as.numeric(first))
  expect_identical(as.numeric(actual[, c("b[1]", "b[2]")]), as.numeric(second))
  expect_identical(as.numeric(actual[, "ordinary"]), as.numeric(later))
  expect_identical(colnames(actual), c("a[1]", "a[2]", "a_ordered_total", "b[1]", "b[2]", "b_ordered_total", "ordinary"))
})

test_that("formula binding refuses ordered container mixtures without changing their sampler", {
  levels <- c("lo", "mid", "hi")
  data <- data.frame(f = ordered(levels, levels = levels))
  bind <- function(p) JAGS_formula(~f, "mu", data,
    list(intercept = prior("point", list(0)), f = p))
  for(contrast in c("cumulative", "cumulative_levels")){
    allocation <- if(contrast == "cumulative") c(.4, .6) else c(.2, .3, .5)
    components <- lapply(c(2, 4), function(total) prior_factor_levels(
      prior_ordered(prior("point", list(total)), allocation, contrast = contrast), levels))
    mixture <- prior_mixture(components)
    set.seed(2)
    draws <- rng(mixture, 4L, transform_factor_samples = FALSE)
    expect_identical(dim(draws), c(4L, length(allocation)))
    expect_true(all(apply(draws, 1L, function(row) any(vapply(c(2, 4),
      function(total) identical(as.numeric(row), total * allocation), logical(1))))))
    condition <- expect_error(bind(mixture),
      "JAGS formula binding is unavailable for mixtures of ordered prior containers. Put mixture or spike-and-slab behavior on 'prior_ordered(total = )' instead.",
      fixed = TRUE, class = "BayesTools_ordered_unavailable")
    expect_s3_class(condition, "BayesTools_ordered_unavailable")
    expect_false(inherits(condition, "BayesTools_ordered_coordinates_unavailable"))
    expect_null(conditionCall(condition))
  }
  total_mixture <- prior_ordered(prior_mixture(list(prior("normal", list(0, 1)),
    prior("point", list(0)))), allocation = c(.4, .6))
  expect_type(bind(total_mixture), "list")
  for(p in list(prior_mixture(list(prior_factor("normal", list(0, 1), contrast = "independent"),
      prior("point", list(0)))),
    prior_mixture(list(prior_factor("normal", list(0, 1), contrast = "treatment"), prior_none())))){
    expect_type(bind(p), "list")
  }
  wf <- prior_mixture(list(prior_weightfunction("one-sided", c(.05), wf_cumulative(c(2, 3))),
    prior_weightfunction("one-sided", c(.05), wf_independent(prior("beta", list(1, 1)))),
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5))), prior_none()),
    is_null = c(FALSE, FALSE, FALSE, TRUE))
  expect_s3_class(wf, "prior.bias_mixture")
})

test_that("D6 ordered declarations reuse the Dirichlet minimum and binding revalidates it", {

  total <- prior("normal", list(0, 1))
  expect_error(prior_ordered(total, allocation = prior("dirichlet", list(alpha = c(.0099, 1)))), "The 'alpha' must be equal or higher than 0.01.", fixed = TRUE)
  p <- prior_ordered(total, allocation = prior("dirichlet", list(alpha = c(.01, 1))))
  expect_s3_class(p, "prior.ordered")
  expect_identical(BayesTools:::.prior_ordered_bind_allocation(list(type = "default_dirichlet"), 2, "f")$alpha, c(1, 1))
  expect_identical(BayesTools:::.prior_ordered_bind_allocation(list(type = "default_dirichlet"), 1, "f"), list(type = "fixed", weights = 1))

  # Internal persisted-spec revalidation; public construction is above.
  allocation <- list(type = "dirichlet", alpha = c(.0099, 1))
  expect_error(BayesTools:::.prior_ordered_bind_allocation(allocation, 2, "f"), "Dirichlet allocation concentrations must be finite and at least 0.01.", fixed = TRUE)
  expect_error(BayesTools:::.prior_ordered_bind_allocation(allocation, 3, "f"), "has length 2, but 3", fixed = TRUE)
  allocation$alpha <- c(.01, 1)
  expect_identical(BayesTools:::.prior_ordered_bind_allocation(allocation, 2, "f"), allocation)
})

test_that("ordered recipes replay primitive chains and batched kernels independently", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi", "top"), 2L),levels=c("lo","mid","hi","top")))
  info <- JAGS_formula(~ f, "mu", data, list(intercept = prior("point", list(0)),
    f = prior_ordered(prior("normal", list(0, 2)), allocation = prior("dirichlet", list(alpha = c(2, 3, 4))))))
  prior <- info$prior_list$mu_f
  spec <- .bt_ordered_spec("mu_f", prior)
  gamma_names <- spec$allocations[[1L]]$gamma_coordinates
  samples <- cbind(mu_f_ordered_total = c(2, -3), matrix(c(1, 3, 2, 2, 3, 1), 2))
  colnames(samples)[-1L] <- gamma_names
  nodes <- .bt_deterministic_nodes(info$prior_list)
  coefficient <- nodes$mu_f
  replay <- .bt_deterministic_node_evaluate(coefficient, .bt_deterministic_lookup(samples))
  expected <- matrix(c(2/6, -9/6, 4/6, -6/6, 6/6, -3/6), 2,
    dimnames = list(NULL, spec$coefficient_names))
  expect_equal(replay, expected, tolerance = 1e-15)
  normalized <- samples[, gamma_names, drop = FALSE] / rowSums(samples[, gamma_names, drop = FALSE])
  colnames(normalized) <- spec$allocations[[1L]]$coordinates
  expect_identical(.bt_deterministic_node_evaluate(coefficient,
    .bt_deterministic_lookup(cbind(samples[, 1, drop = FALSE], normalized))), replay)
  overflow <- samples
  overflow[, gamma_names] <- .Machine$double.xmax
  expect_error(.bt_deterministic_node_evaluate(coefficient,
    .bt_deterministic_lookup(overflow)), class = "BayesTools_ordered_invalid_state")
  kernel <- JAGS_ordered_density_kernel(info$prior_list["mu_f"])
  oracle <- dnorm(samples[, 1], 0, 2, log = TRUE) +
    dgamma(samples[, 2], 2, 1, log = TRUE) + dgamma(samples[, 3], 3, 1, log = TRUE) +
    dgamma(samples[, 4], 4, 1, log = TRUE)
  expect_equal(kernel(samples), oracle, tolerance = 1e-14)
  expect_equal(JAGS_marglik_priors_rows_evaluator(info$prior_list)(samples), oracle, tolerance = 1e-14)
  chart <- setNames(list(list(kind = "group", J = 1L)), spec$allocations[[1L]]$key)
  expect_equal(JAGS_ordered_density_kernel(info$prior_list["mu_f"], chart)(samples),
    dnorm(samples[, 1], 0, 2, log = TRUE) + dbeta(samples[, 2]/6, 2, 7, log = TRUE), tolerance = 1e-14)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(cbind(mu_intercept=0,samples, replay))),
    info$prior_list, formula_design = list(mu = info$formula_design))
  stale <- cbind(samples, replay * 99)
  stale[, gamma_names] <- stale[, gamma_names] * c(3, 1, 2)
  expect_equal(JAGS_evaluate_deterministic(fit, stale, nodes = "mu_f"),
    .bt_deterministic_node_evaluate(coefficient, .bt_deterministic_lookup(stale)), tolerance = 0)
  expect_equal(JAGS_ordered_parameter_spec(fit)$mu_f$coefficient_names, spec$coefficient_names)
  selected <- parameter_catalog_resolve(parameter_catalog(fit),"f[mid]",namespace="mu")
  supplied <- cbind(samples,replay*99)
  supplied[,gamma_names] <- supplied[,gamma_names]*c(3,1,2)
  expect_equal(as.numeric(as.matrix(parameter_draws(fit,selected,model_samples=supplied))),
    supplied[,1]*supplied[,gamma_names[1]]/rowSums(supplied[,gamma_names,drop=FALSE]),tolerance=1e-15)
  expect_error(parameter_draws(fit,selected,model_samples=supplied[,setdiff(colnames(supplied),gamma_names),drop=FALSE]),
    class="BayesTools_ordered_coordinates_unavailable")
  expect_equal(as.numeric(as.matrix(parameter_draws(fit,selected))),replay[,1],tolerance=1e-15)
  supplied_view <- .posterior_atoms_set(supplied,.posterior_atoms_new(column_names=colnames(supplied)))
  supplied_view <- posterior_transform(supplied_view,"exp")
  expect_identical(as.numeric(as.matrix(parameter_draws(fit,selected,model_samples=supplied_view))),
    as.numeric(supplied_view[,spec$coefficient_names[1]]))
  old_priors <- attr(fit, "prior_list", exact = TRUE)
  old_metadata <- attr(old_priors$mu_f, "ordered_metadata", exact = TRUE)
  old_metadata$numeric_literals <- NULL
  attr(old_priors$mu_f, "ordered_metadata") <- old_metadata
  attr(fit, "prior_list") <- old_priors
  expect_error(JAGS_ordered_parameter_spec(fit), class = "BayesTools_ordered_metadata_unavailable")
  for(read in list(function() JAGS_deterministic_nodes(fit),
      function() JAGS_deterministic_evaluator(fit, nodes = "mu_f"),
      function() JAGS_evaluate_deterministic(fit, supplied, nodes = "mu_f"))){
    condition <- expect_error(read(),
      "Fitted ordered numeric provenance is unavailable. Refit the model with this version of BayesTools.",
      fixed = TRUE, class = "BayesTools_ordered_metadata_unavailable")
    expect_s3_class(condition, "BayesTools_refit_required")
    expect_null(conditionCall(condition))
  }
  expect_error(parameter_draws(fit,selected,model_samples=supplied_view),class="BayesTools_ordered_metadata_unavailable")
  expect_identical(as.numeric(as_mixed_posteriors(fit,"mu_intercept")$mu_intercept),rep(0,2))
})

test_that("ordered kernels preserve alternating mixture states and pair chart Jacobians", {
  data <- data.frame(f=ordered(c("lo","mid","hi","top")))
  total <- prior_mixture(list(prior("normal",list(-1,2),prior_weights=.2),prior("gamma",list(3,2),prior_weights=.8)))
  info <- JAGS_formula(~f,"mu",data,list(intercept=prior("point",list(0)),
    f=prior_ordered(total,allocation=prior("dirichlet",list(c(2,3,4))))))
  spec <- .bt_ordered_spec("mu_f",info$prior_list$mu_f)
  record <- spec$allocations[[1L]]
  samples <- cbind(c(1,2,1,2),c(-2,99,1,99),c(99,.5,99,2),
    matrix(c(1,2,3,4,2,3,4,5,3,4,5,6),4))
  colnames(samples) <- c(spec$total_node$spec$indicator,
    unlist(lapply(spec$total_node$spec$components,`[[`,"coordinates")),record$gamma_coordinates)
  total_density <- c(dnorm(-2,-1,2,log=TRUE),dgamma(.5,3,2,log=TRUE),
    dnorm(1,-1,2,log=TRUE),dgamma(2,3,2,log=TRUE))
  eta <- samples[,record$gamma_coordinates,drop=FALSE]
  gamma_density <- dgamma(eta[,1],2,1,log=TRUE)+dgamma(eta[,2],3,1,log=TRUE)+dgamma(eta[,3],4,1,log=TRUE)
  expect_equal(JAGS_ordered_density_kernel(info$prior_list["mu_f"])(samples),total_density+gamma_density,tolerance=1e-14)
  chart <- setNames(list(list(kind="pair",indices=c(1L,3L))),record$key)
  pair_sum <- eta[,1]+eta[,3]
  share <- eta[,1]/pair_sum
  pair_kernel <- total_density+dbeta(share,2,4,log=TRUE)+dgamma(eta[,2],3,1,log=TRUE)
  expect_equal(JAGS_ordered_density_kernel(info$prior_list["mu_f"],chart)(samples),pair_kernel,tolerance=1e-14)
  # The eta-pair -> (h,p) change of variables has Jacobian h.
  expect_equal(gamma_density+log(pair_sum),
    dbeta(share,2,4,log=TRUE)+dgamma(pair_sum,6,1,log=TRUE)+dgamma(eta[,2],3,1,log=TRUE),tolerance=1e-14)
  malformed <- list(1,list(list(kind="group",J=1L)),setNames(chart,"unknown"),
    setNames(c(chart,chart),rep(record$key,2)),setNames(chart,NA_character_),
    setNames(list(NULL),record$key),setNames(list(list(kind=c("pair","group"),indices=1:2)),record$key),
    setNames(list(list(kind="bad",indices=1:2)),record$key),
    setNames(list(list(kind="group",J=integer())),record$key),
    setNames(list(list(kind="group",J=1:3)),record$key),
    setNames(list(list(kind="group",J=c(1,1))),record$key),
    setNames(list(list(kind="pair",indices=c(1,1.5))),record$key),
    setNames(list(list(kind="pair",indices=c(1,NA))),record$key),
    setNames(list(list(kind="pair",indices=1L)),record$key))
  for(candidate in malformed) expect_error(JAGS_ordered_density_kernel(info$prior_list["mu_f"],candidate),"allocation_chart",fixed=TRUE)
})

test_that("ordered active kernels localize each nested total event", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2)), g = factor(rep(c("a", "b"), 3)))
  total <- prior_spike_and_slab(prior("normal", list(0, 2)), prior("point", list(.5)))
  info <- JAGS_formula(~ f*g, "mu", data, list(intercept=prior("point", list(0)),
    f=prior_ordered(prior("normal", list(0,1))),
    g=prior_factor("normal", list(0,1), contrast="independent"), `f:g`=prior_ordered(total, allocation=c(.2,.8))))
  p <- info$prior_list[["mu_f:g"]]
  parameter <- names(info$prior_list)[vapply(info$prior_list, function(p) is.prior.ordered(p) && is.prior.spike_and_slab(p$total), logical(1))]
  p <- info$prior_list[[parameter]]
  spec <- .bt_ordered_spec(parameter, p)
  samples <- cbind(c(0, 1, 0, 1), c(1, 2, 3, 4), c(2, 3, 4, 5))
  colnames(samples) <- c(spec$total_node$spec$indicator, spec$total_node$spec$components[[1L]]$coordinates)
  expected <- c(0, dnorm(2,0,2,log=TRUE)+dnorm(3,0,2,log=TRUE), 0,
    dnorm(4,0,2,log=TRUE)+dnorm(5,0,2,log=TRUE))
  expect_equal(JAGS_ordered_density_kernel(setNames(list(p), parameter))(samples), expected, tolerance = 1e-14)
  totals <- .bt_ordered_total_values(spec, .bt_deterministic_lookup(samples))
  expect_equal(unname(totals), unname(samples[, -1L, drop=FALSE] * samples[, 1L]), tolerance=0)
  samples[1, 1] <- 2
  expect_error(JAGS_ordered_density_kernel(setNames(list(p), parameter))(samples), class="BayesTools_ordered_invalid_state")
})

test_that("ordered numeric literals round trip supplied doubles", {
  values <- c(1/6, 1/3, .6, 2.5, 1.2345678901234567)
  expect_identical(as.double(vapply(values, .prior_ordered_format_number, character(1))), values)
  p <- prior_ordered(prior("point", list(2.5)), allocation = c(1, 2, 3)/6)
  attr(p, "levels") <- 4L
  p <- .prior_ordered_default_bound(p, "mu_f")
  literal <- .prior_ordered_metadata(p)$numeric_literals
  expect_identical(as.double(literal$total$location), 2.5)
  expect_identical(as.double(literal$allocations[[1L]]), c(1, 2, 3)/6)
})

test_that("ordered source retention follows selected model and row without extra sampling", {
  data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 2)))
  make <- function(allocation = NULL, id = NULL){
    info <- JAGS_formula(~f, "mu", data, list(intercept=prior("point",list(0)),
      f=prior_ordered(prior("normal",list(0,1)), allocation=allocation, id=id)))
    p <- info$prior_list$mu_f
    spec <- .bt_ordered_spec("mu_f", p)
    theta <- c(0, 1, 3, 2, -1, 4)
    weights <- if(is.null(allocation)) cbind(c(1,2,3,4,5,6)/7, c(6,5,4,3,2,1)/7) else{
      matrix(rep(allocation, each=6), 6)
    }
    coefficients <- weights * theta
    colnames(coefficients) <- spec$coefficient_names
    draws <- cbind(mu_intercept=0, coefficients, mu_f_ordered_total=theta)
    if(is.null(allocation)){
      gamma <- weights * 14
      colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
      draws <- cbind(draws, gamma)
    }
    fit <- structure(list(mcmc=coda::mcmc.list(coda::mcmc(draws))), class=c("runjags","BayesTools_fit","list"))
    attr(fit,"prior_list") <- info$prior_list
    attr(fit,"formula_design") <- list(mu=info$formula_design)
    attr(fit,"parameter_map") <- .bt_build_parameter_map(colnames(draws), prior_list = info$prior_list,
      formula_design=list(mu=info$formula_design))
    fit <- .bt_attach_fit_contract(.bt_attach_draw_geometry(fit))
    fit <- attach_test_parameter_map(fit)
    list(fit=fit, prior=p, theta=theta, weights=weights)
  }
  models <- list(make(), make(c(.25,.75)))
  fits <- lapply(models, `[[`, "fit")
  priors <- lapply(models, `[[`, "prior")
  mixed <- .mix_posteriors.factor(fits,priors,"mu_f",c(.3,.7),seed=13,n_samples=40)
  source <- .bt_meta_get(mixed,"ordered_source")
  raw <- .bt_ordered_source_display(source)$samples
  expected <- t(vapply(seq_len(nrow(raw)),function(i){
    m <- models[[source$model[[i]]]]
    row <- source$draw_index[[i]]
    c(m$theta[[row]],m$weights[row,])
  }, numeric(3)))
  expect_identical(unname(raw), unname(expected))
  subset <- .bt_draws_subset_rows(mixed,c(2L,7L,14L))
  expect_identical(.bt_meta_get(subset,"ordered_source")$primitives,source$primitives[c(2L,7L,14L),,drop=FALSE])
  expect_identical(.bt_meta_get(subset,"draw_index"),source$draw_index[c(2L,7L,14L)])
  reversed <- .mix_posteriors.factor(rev(fits),rev(priors),"mu_f",c(.7,.3),seed=13,n_samples=40)
  expect_identical(colnames(.bt_ordered_source_display(.bt_meta_get(reversed,"ordered_source"))$samples),colnames(raw))
  absent <- .mix_posteriors.factor(c(fits[1],list(structure(list(),class="null_model"))),
    c(priors[1],list(prior("point",list(0)))),"mu_f",c(.5,.5),seed=13,n_samples=40)
  absent_source <- .bt_meta_get(absent,"ordered_source")
  absent_raw <- .bt_ordered_source_display(absent_source)
  expect_true(all(is.na(absent_raw$samples[absent_source$model==2L,-1L,drop=FALSE])))
  expect_true(all(absent_raw$samples[absent_source$model==2L,1L]==0))
  expect_identical(unname(absent_raw$undefined_draws),rep("ordered_parameterization",2L))
  zero_count <- .mix_posteriors.factor(fits,priors,"mu_f",c(0,1),seed=13,n_samples=20)
  expect_identical(colnames(.bt_ordered_source_display(.bt_meta_get(zero_count,"ordered_source"))$samples),colnames(raw))
  expect_error(.bt_meta_set(mixed,"draw_index",1:2),"one row per draw")
  table <- ensemble_estimates_table(list(mu_f=absent),"mu_f")
  expect_match(paste(attr(table,"footnotes"),collapse=" "),"ordered allocation parameterization")
})

test_that("ordered scalar measures use declared contractions and exact primitive identities", {
  data <- data.frame(f=ordered(rep(c("lo","mid","hi","top"),10),levels=c("lo","mid","hi","top")))
  make <- function(total_prior,allocation=NULL,theta=NULL,indicator=NULL,gamma=NULL){
    info <- JAGS_formula(~f,"mu",data,list(intercept=prior("point",list(0)),f=prior_ordered(total_prior,allocation=allocation)))
    spec <- .bt_ordered_spec("mu_f",info$prior_list$mu_f)
    if(is.null(theta)) theta <- if(is.prior.point(total_prior)) rep(total_prior$parameters$location,40) else seq_len(40)/10
    draws <- matrix(theta,ncol=1,dimnames=list(NULL,spec$total_names))
    if(!is.null(indicator)){
      variable <- theta
      draws[,1] <- theta * indicator
      draws <- cbind(draws,variable,indicator)
      colnames(draws)[2:3] <- c(spec$total_node$spec$components[[1L]]$coordinates,spec$total_node$spec$indicator)
    }
    if(is.null(allocation)){
      if(is.null(gamma)) gamma <- matrix(rep(c(1,2,3),each=40),40)
      if(!is.null(indicator)) gamma[indicator==1,] <- matrix(rep(c(3,2,1),each=sum(indicator==1)),ncol=3)
      colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
      draws <- cbind(draws,gamma)
    }
    coefficients <- .bt_deterministic_node_evaluate(.bt_dnode_ordered_coefficients(spec),.bt_deterministic_lookup(draws))
    chains <- coda::mcmc.list(coda::mcmc(cbind(mu_intercept=0,coefficients,draws)))
    fit <- structure(list(mcmc=chains,sample=40L,summary.pars=list(mutate=NULL),monitor=colnames(chains[[1L]])),
      class=c("runjags","BayesTools_fit","list"))
    attr(fit,"prior_list") <- info$prior_list
    attr(fit,"formula_design") <- list(mu=info$formula_design)
    fit <- attach_test_parameter_map(fit)
    list(fit=fit,spec=spec,draws=draws)
  }
  zero <- make(prior("normal",list(0,1)),c(0,.5,.5))
  raw <- as_mixed_posteriors(zero$fit,"mu_f")
  levels <- transform_factor_samples(raw)$mu_f
  atoms <- .posterior_atoms_get(levels)
  expect_identical(.posterior_atoms_for_column(atoms,2)$locations,matrix(0,1,1,dimnames=list(NULL,"mu_f[mid]")))
  expect_identical(.posterior_atoms_for_column(atoms,2)$mass,1)
  expect_false(posterior_atoms_free(levels))
  fixed <- make(prior("point",list(2.5)),c(1,2,3)/6)
  fixed_levels <- transform_factor_samples(as_mixed_posteriors(fixed$fit,"mu_f"))$mu_f
  fixed_atoms <- .posterior_atoms_get(fixed_levels)
  expect_identical(as.numeric(fixed_levels[1,]),c(0,2.5*(1/6),1.25,2.5))
  expect_identical(as.numeric(fixed_atoms$locations[1,]),as.numeric(fixed_levels[1,]))
  random <- make(prior("point",list(2.5)))
  random_levels <- transform_factor_samples(as_mixed_posteriors(random$fit,"mu_f"))$mu_f
  random_atoms <- .posterior_atoms_get(random_levels)
  expect_length(.posterior_atoms_for_column(random_atoms,2)$mass,0)
  expect_identical(.posterior_atoms_for_column(random_atoms,4)$mass,1)
  expect_true(all(as.numeric(random_levels[,4])==2.5))
  expect_identical(nrow(random_atoms$locations),0L)
  selection <- parameter_catalog_resolve(parameter_catalog(random$fit),"f[top]",namespace="mu")
  scalar <- parameter_mixed_posterior(random$fit,selection)
  expect_identical(.posterior_atoms_get(scalar)$mass,1)
  expect_identical(.posterior_atoms_get(scalar)$locations,matrix(2.5,1,1))
  expect_true(all(as.numeric(scalar)==2.5))
  bernoulli <- make(prior("bernoulli",list(.5)),c(.25,.25,.5),theta=rep(c(0,1),20))
  bernoulli_levels <- transform_factor_samples(as_mixed_posteriors(bernoulli$fit,"mu_f"))$mu_f
  expect_identical(.posterior_atoms_for_column(.posterior_atoms_get(bernoulli_levels),4)$mass,c(.5,.5))
  expect_error(as_mixed_posteriors(bernoulli$fit,"mu_f",conditional="mu_f"),"not a conditional parameter")
  spike <- make(prior_spike_and_slab(prior("normal",list(0,1)),prior("point",list(.5))),indicator=rep(c(0,1),20))
  conditional <- as_mixed_posteriors(spike$fit,"mu_f",conditional="mu_f")
  expect_identical(.bt_meta_get(conditional$mu_f,"draw_index"),seq(2L,40L,2L))
  conditional_levels <- transform_factor_samples(conditional)$mu_f
  expect_equal(mean(as.numeric(conditional_levels[,4])),2.1,tolerance=1e-15)
  total_selection <- parameter_catalog_resolve(parameter_catalog(spike$fit),"mu_f_ordered_total")
  total_conditional <- parameter_mixed_posterior(spike$fit,total_selection,conditional=TRUE)
  expect_equal(mean(as.numeric(total_conditional)),2.1,tolerance=1e-15)
  raw_table <- JAGS_estimates_table(spike$fit,conditional=TRUE,transform_factors=FALSE,remove_diagnostics=TRUE)
  level_table <- JAGS_estimates_table(spike$fit,conditional=TRUE,transform_factors=TRUE,remove_diagnostics=TRUE)
  expect_equal(raw_table["(mu) f_ordered_total","Mean"],2.1,tolerance=1e-15)
  expect_equal(raw_table["(mu) f_ordered_allocation[mid]","Mean"],.5,tolerance=1e-15)
  expect_equal(level_table["(mu) f[top]","Mean"],2.1,tolerance=1e-15)
  unconditional_levels <- transform_factor_samples(as_mixed_posteriors(spike$fit,"mu_f"))$mu_f
  subset <- .bt_draws_subset_rows(unconditional_levels,seq(2L,40L,2L))
  expect_length(.posterior_atoms_for_column(.posterior_atoms_get(subset),4)$mass,0L)
  transformed <- posterior_transform(random_levels,"lin",list(a=1,b=-2))
  transformed_atoms <- .posterior_atoms_get(transformed)
  expect_identical(.posterior_atoms_for_column(transformed_atoms,4)$locations,
    matrix(-4,1,1,dimnames=list(NULL,"mu_f[top]")))
  transformed_subset <- .bt_draws_subset_rows(transformed,c(1L,3L,5L))
  expect_identical(as.numeric(transformed_subset[,4]),rep(-4,3))
  expect_identical(.posterior_atoms_for_column(.posterior_atoms_get(transformed_subset),4)$mass,1)
  renamed <- random_levels
  colnames(renamed) <- paste0("renamed_",seq_len(4))
  expect_identical(names(.posterior_atoms_get(renamed)$marginals),colnames(renamed))
  expect_identical(.posterior_atoms_for_column(.posterior_atoms_get(renamed),"renamed_4")$mass,1)
  sampled <- JAGS_ordered_parameter_spec(spike$fit,draws=spike$draws)$mu_f$sampled_parameters
  expect_identical(colnames(sampled$values),c("mu_f_ordered_total","mu_f_ordered_allocation[mid]","mu_f_ordered_allocation[hi]","mu_f_ordered_allocation[top]"))
  expect_identical(sampled$allocation_names,colnames(sampled$values)[-1L])
  expect_identical(sampled$values[1,-1L],setNames(c(1,2,3)/6,sampled$allocation_names))
  first <- make(prior_spike_and_slab(prior("normal",list(0,1)),prior("point",list(.2))),
    indicator=c(rep(0,30),rep(1,10)))
  second <- make(prior_spike_and_slab(prior("normal",list(0,1)),prior("point",list(.8))),
    indicator=c(rep(0,10),rep(1,30)))
  models <- list(list(fit=first$fit,marglik=bridgesampling_object(log(.6)),prior_weights=1),
    list(fit=second$fit,marglik=bridgesampling_object(log(.4)),prior_weights=1))
  ensemble <- mix_posteriors(models,"mu_f",list(mu_f=c(FALSE,FALSE)),conditional=TRUE,seed=13,n_samples=100)
  source <- .bt_meta_get(ensemble$mu_f,"ordered_source")
  expect_equal(source$conditioning$post_probs,c(1/3,2/3),tolerance=1e-15)
  expect_equal(source$conditioning$prior_probs,c(.2,.8),tolerance=1e-15)
  expect_identical(source$conditioning$posterior_fractions,c(.25,.75))
  expect_true(all(source$draw_index[source$model==1L] %in% 31:40))
  expect_true(all(source$draw_index[source$model==2L] %in% 11:40))
  context <- .bt_meta_get(ensemble$mu_f,"prior_context")
  expect_equal(context$model_weights,c(.2,.8),tolerance=1e-15)
  expect_true(all(vapply(context$prior_lists,function(priors) is.prior.ordered(priors$mu_f) && !is.prior.mixture(priors$mu_f$total),logical(1))))
  expect_true(is.prior.spike_and_slab(source$models[[1L]]$prior$total))
  expect_length(.posterior_atoms_get(ensemble$mu_f)$mass,0L)
  expression_normal <- make(prior("normal",list(0,expression(tau))),theta=rep(0,40))
  expression_raw <- as_mixed_posteriors(expression_normal$fit,"mu_f")
  expect_true(posterior_atoms_free(expression_raw$mu_f))
  expression_levels <- transform_factor_samples(expression_raw)$mu_f
  expect_length(.posterior_atoms_for_column(.posterior_atoms_get(expression_levels),4)$mass,0L)
  expect_length(.plot_data_factor_column_atoms(expression_levels),4L)
  expression_varying <- make(prior("normal",list(0,expression(tau))),theta=seq_len(40)/10)
  expect_length(.plot_data_samples.factor(as_mixed_posteriors(expression_varying$fit,"mu_f"),"mu_f",128,NULL,NULL,NULL),3L)
  expression_marginal <- marginal_posterior(expression_raw,"mu_f",use_formula=FALSE,prior_samples=FALSE)
  expect_true(posterior_atoms_free(expression_marginal[[4L]]))
  expression_spike <- make(prior_spike_and_slab(prior("normal",list(0,expression(tau))),prior("point",list(.5))),
    theta=rep(2,40),indicator=rep(c(0,1),20))
  expression_spike_raw <- as_mixed_posteriors(expression_spike$fit,"mu_f")
  expression_spike_levels <- transform_factor_samples(expression_spike_raw)$mu_f
  expect_identical(.posterior_atoms_for_column(.posterior_atoms_get(expression_spike_levels),4)$mass,.5)
  expect_identical(as.numeric(expression_spike_levels[seq(1,40,2),4]),rep(0,20))
  expect_length(.plot_data_factor_column_atoms(expression_spike_levels),4L)
  expect_identical(.posterior_atoms_get(marginal_posterior(expression_spike_raw,"mu_f",use_formula=FALSE,prior_samples=FALSE)[[4L]])$mass,.5)
  expression_point <- make(prior("point",list(expression(tau))),theta=rep(2,40))
  expression_point_raw <- as_mixed_posteriors(expression_point$fit,"mu_f")
  expect_error(posterior_atoms_free(expression_point_raw$mu_f), class = "BayesTools_formula_measure_unavailable")
  point_levels <- transform_factor_samples(expression_point_raw)$mu_f
  expect_error(.posterior_atoms_get(point_levels), class = "BayesTools_formula_measure_unavailable")
  expect_identical(length(.posterior_atoms_get(point_levels, allow_partial = TRUE)$marginals),4L)
  expect_identical(vapply(.posterior_atoms_get(point_levels, allow_partial = TRUE)$marginals,is.null,logical(1)),
    setNames(c(FALSE,TRUE,TRUE,TRUE),colnames(point_levels)))
  point_marginal <- marginal_posterior(expression_point_raw,"mu_f",use_formula=FALSE,prior_samples=FALSE)
  expect_equal(.posterior_atoms_get(point_marginal[[1L]])$mass, 1, tolerance = 0)
  for(level in point_marginal[-1L]){
    class(level) <- unique(c(class(level), "marginal_posterior"))
    expect_error(Savage_Dickey_BF(level, silent = TRUE), class = "BayesTools_formula_measure_unavailable")
  }
  expect_error(.plot_data_factor_column_atoms(point_levels),class="BayesTools_ordered_expression_unavailable")
  expect_error(.plot_data_samples.factor(expression_point_raw,"mu_f",128,NULL,NULL,NULL),class="BayesTools_ordered_expression_unavailable")
  zero_expression_point <- make(prior("point",list(expression(tau))),c(0,.5,.5),theta=rep(2,40))
  zero_selection <- parameter_catalog_resolve(parameter_catalog(zero_expression_point$fit),"f[mid]",namespace="mu")
  expect_identical(.posterior_atoms_get(parameter_mixed_posterior(zero_expression_point$fit,zero_selection))$mass,1)
  expect_identical(.prior_ordered_linear_range(zero_expression_point$spec$prior,1,1L,1e-4),c(0,0))
  zero_distribution <- .prior_ordered_linear_distribution(zero_expression_point$spec$prior,1,1L,n_grid=32)
  expect_identical(zero_distribution$points,data.frame(x=0,p=1))
  expect_identical(attr(zero_distribution,"ordered_measure",exact=TRUE)$allocation,
    zero_expression_point$spec$allocations[[1L]]$spec)
  point_formula_marginal <- marginal_posterior(expression_point_raw,"mu_f",formula=~0+f,prior_samples=FALSE)
  expect_equal(.posterior_atoms_get(point_formula_marginal[[1L]])$mass, 1, tolerance = 0)
  for(level in point_formula_marginal[-1L]){
    class(level) <- unique(c(class(level), "marginal_posterior"))
    expect_error(Savage_Dickey_BF(level, silent = TRUE), class = "BayesTools_formula_measure_unavailable")
  }
  expression_selection <- parameter_catalog_resolve(parameter_catalog(expression_point$fit),"f[top]",namespace="mu")
  expect_identical(as.numeric(as.matrix(parameter_draws(expression_point$fit,expression_selection))),rep(2,40))
  condition <- tryCatch(parameter_mixed_posterior(expression_point$fit,expression_selection),error=identity)
  expect_s3_class(condition,"BayesTools_ordered_expression_unavailable")
  expect_identical(conditionMessage(condition),paste0(
    "Ordered scalar measure is unavailable for expression totals with unclassified stochastic ancestry. ",
    "Use a supported scalar total prior with numeric point locations, or inspect fitted snapshot values with 'parameter_draws()'."))
  expression_bernoulli <- make(prior("bernoulli",list(expression(prob))),c(.25,.25,.5),theta=rep(c(0,1),20))
  expression_bernoulli_levels <- transform_factor_samples(as_mixed_posteriors(expression_bernoulli$fit,"mu_f"))$mu_f
  expect_identical(.posterior_atoms_for_column(.posterior_atoms_get(expression_bernoulli_levels),4)$mass,c(.5,.5))
  primitive_selection <- parameter_catalog_resolve(parameter_catalog(random$fit),"f[top]",namespace="mu")
  supplied <- cbind(random$draws,matrix(99,40,3,dimnames=list(NULL,random$spec$coefficient_names)))
  supplied[,random$spec$total_names] <- seq_len(40)/5
  expect_identical(as.numeric(as.matrix(parameter_draws(random$fit,primitive_selection,model_samples=supplied))),rep(2.5,40))
  normal_selection <- parameter_catalog_resolve(parameter_catalog(zero$fit),"f[top]",namespace="mu")
  supplied_normal <- cbind(zero$draws,matrix(99,40,3,dimnames=list(NULL,zero$spec$coefficient_names)))
  supplied_normal[,zero$spec$total_names] <- seq_len(40)/5
  expect_identical(as.numeric(as.matrix(parameter_draws(zero$fit,normal_selection,model_samples=supplied_normal))),seq_len(40)/5)
  expect_error(parameter_draws(zero$fit,normal_selection,model_samples=supplied_normal[,zero$spec$coefficient_names,drop=FALSE]),
    class="BayesTools_ordered_coordinates_unavailable")

  # Semantic producers use the same primitive contraction even for continuous
  # intermediate levels and when a derived monitor is stale.
  continuous <- make(prior("normal",list(0,1)),
    gamma=matrix((seq_len(120)^.7 + .3)/7,40))
  original_monitors <- continuous$fit$mcmc
  stale_fit <- continuous$fit
  stale_fit$mcmc[[1L]][1L,continuous$spec$coefficient_names] <-
    stale_fit$mcmc[[1L]][1L,continuous$spec$coefficient_names] + 1
  continuous_raw <- as_mixed_posteriors(stale_fit,"mu_f")
  continuous_levels <- transform_factor_samples(continuous_raw)$mu_f
  continuous_marginal <- marginal_posterior(continuous_raw,"mu_f",use_formula=FALSE,prior_samples=FALSE)
  expected <- vapply(seq_len(3),function(i){
    weights <- setNames(as.numeric(seq_len(3)==i),continuous$spec$coefficient_names)
    JAGS_ordered_parameter_spec(stale_fit,weights=weights,draws=continuous$draws)$values
  },numeric(40))
  expect_identical(unname(.bt_draws_plain(continuous_raw$mu_f)),expected)
  level_names <- c("mid","hi","top")
  expected_levels <- vapply(level_names,function(level){
    selection <- parameter_catalog_resolve(parameter_catalog(stale_fit),paste0("f[",level,"]"),namespace="mu")
    key <- selection$quantities$extraction_key[[1L]]
    projection <- JAGS_ordered_parameter_spec(stale_fit,weights=setNames(key$weights,key$dependencies),draws=continuous$draws)
    expect_identical(as.numeric(as.matrix(parameter_draws(stale_fit,selection))),projection$values)
    expect_identical(as.numeric(parameter_mixed_posterior(stale_fit,selection)),projection$values)
    expect_identical(projection$state,rep("continuous",40))
    projection$values
  },numeric(40))
  expect_identical(unname(.bt_draws_plain(continuous_levels)[,-1L]),unname(expected_levels))
  expect_identical(vapply(continuous_marginal[level_names],as.numeric,numeric(40)),expected_levels)
  formula_marginal <- marginal_posterior(continuous_raw,"mu_f",formula=~0+f,prior_samples=FALSE)
  expect_identical(vapply(formula_marginal[level_names],as.numeric,numeric(40)),expected_levels)
  expect_identical(as.numeric(continuous_levels[,4L]),as.numeric(continuous$draws[,continuous$spec$total_names]))
  expect_true(posterior_atoms_free(continuous_raw$mu_f))
  expect_identical(continuous$fit$mcmc,original_monitors)
  expect_identical(stale_fit$mcmc[[1L]][1L,continuous$spec$coefficient_names],
    original_monitors[[1L]][1L,continuous$spec$coefficient_names] + 1)
  source <- .bt_meta_get(continuous_raw$mu_f,"ordered_source")
  expect_identical(.bt_meta_get(continuous_levels,"ordered_source")$primitives,source$primitives)
  expect_identical(.bt_meta_get(continuous_levels,"draw_index"),.bt_meta_get(continuous_raw$mu_f,"draw_index"))
  expect_identical(.bt_meta_get(continuous_levels,"condition"),.bt_meta_get(continuous_raw$mu_f,"condition"))

  scalar <- parameter_mixed_posterior(stale_fit,
    parameter_catalog_resolve(parameter_catalog(stale_fit),"f[hi]",namespace="mu"))
  scalar <- .posterior_atoms_set(scalar,.posterior_atoms_rename_columns(.posterior_atoms_get(scalar),"mu_f[hi]"))
  scalar_columns <- colnames(.posterior_atoms_get(scalar)$locations)
  for(transformation in c("lin","exp","tanh")){
    arguments <- if(transformation=="lin") list(a=1,b=-2) else NULL
    view <- posterior_transform(scalar,transformation,arguments)
    corrupted_view <- .bt_draws_transform_values(view,function(values) values + 1)
    restored <- .bt_ordered_source_semantics(corrupted_view,matrix(1,1,1),scalar_columns)
    expect_identical(as.numeric(restored),as.numeric(view))
    expect_identical(.bt_meta_get(restored,"ordered_source")$primitives,.bt_meta_get(view,"ordered_source")$primitives)
    expect_identical(.bt_meta_get(restored,"condition"),.bt_meta_get(view,"condition"))
  }
  custom <- posterior_transform(scalar,list(fun=function(x) x+2,inv=function(x) x-2,jac=function(x) rep(1,length(x))))
  custom_values <- as.numeric(custom)
  custom <- .bt_ordered_source_semantics(custom,matrix(1,1,1),scalar_columns)
  expect_identical(as.numeric(custom),custom_values)
  expect_true(all(.bt_ordered_source_project(.bt_meta_get(custom,"ordered_source"),1)$state=="unavailable"))
  unsupported_source <- .bt_meta_get(scalar,"ordered_source")
  unsupported_source$models[[2L]] <- .bt_ordered_source_model(
    prior_factor_levels(prior_factor("normal",list(0,1),contrast="treatment"),levels(data$f)),"mu_f")
  unsupported_source$model[21:40] <- 2L
  unsupported_source$model_probabilities <- c(.5,.5)
  probability_pair <- .model_probability_pair(c(.5,.5),log(c(.5,.5)),"posterior","raw")
  unsupported_source$model_log_probabilities <- probability_pair$logs
  unsupported_source$model_probability_declaration <- probability_pair$declaration
  unsupported_source$projection_context <- NULL
  unsupported <- .bt_meta_set(.bt_draws_transform_values(scalar,function(values) values+1),"ordered_source",unsupported_source)
  restored <- .bt_ordered_source_semantics(unsupported,matrix(1,1,1),scalar_columns)
  expect_identical(as.numeric(restored),c(as.numeric(scalar)[1:20],as.numeric(scalar)[21:40]+1))
  undefined <- .bt_draws_transform_values(scalar,function(values){ values[1L] <- NA_real_; values })
  undefined <- .bt_meta_set(undefined,"undefined_draws",setNames("ordered_parameterization",scalar_columns))
  restored <- .bt_ordered_source_semantics(undefined,matrix(1,1,1),scalar_columns)
  expect_identical(as.numeric(restored),as.numeric(undefined))
  expect_identical(.bt_meta_get(restored,"undefined_draws"),.bt_meta_get(undefined,"undefined_draws"))

  # A stored near-cancellation weight remains nonzero and continuous.
  weight <- (.1+.2)-.3
  expect_true(weight!=0)
  mid <- parameter_mixed_posterior(stale_fit,
    parameter_catalog_resolve(parameter_catalog(stale_fit),"f[mid]",namespace="mu"))
  mid <- .posterior_atoms_set(mid,.posterior_atoms_rename_columns(.posterior_atoms_get(mid),"mu_f[mid]"))
  near <- .bt_ordered_source_semantics(mid,matrix(weight,1,1),colnames(.posterior_atoms_get(mid)$locations))
  near_projection <- JAGS_ordered_parameter_spec(stale_fit,
    weights=setNames(c(weight,0,0),continuous$spec$coefficient_names),draws=continuous$draws)
  expect_identical(as.numeric(near),near_projection$values)
  expect_identical(near_projection$state,rep("continuous",40))
  expect_true(all(as.numeric(near)!=0))
})

test_that("mixture numeric literals preserve round-trip dcat weights only when requested", {
  weights <- c(1/7,6/7)
  total <- prior_mixture(list(prior("point",list(0),prior_weights=weights[1]),prior("normal",list(0,1),prior_weights=weights[2])))
  stored <- attr(total,"prior_weights",exact=TRUE)
  syntax <- .JAGS_prior.mixture(total,"theta",numeric_literals=TRUE)
  literal <- regmatches(syntax,regexec("dcat\\(c\\(([^)]+)\\)\\)",syntax))[[1L]][2L]
  expect_identical(as.numeric(strsplit(literal,",",fixed=TRUE)[[1L]]),stored)
  legacy <- .JAGS_prior.mixture(total,"theta",numeric_literals=FALSE)
  expect_match(legacy,paste0("dcat(c(",paste0(stored,collapse=", "),"))"),fixed=TRUE)
  expect_false(identical(legacy,syntax))
})

test_that("ordered projections localize ordinary point and continuous components", {
  data <- data.frame(f=ordered(c("lo","mid","hi"),levels=c("lo","mid","hi")))
  info <- JAGS_formula(~f,"mu",data,list(
    intercept=prior_mixture(list(prior("point",list(0)),prior("point",list(1)),prior("normal",list(0,1)))),
    f=prior_ordered(prior("point",list(2.5)))))
  p <- info$prior_list$mu_f
  spec <- .bt_ordered_spec("mu_f",p)
  draws <- cbind(mu_intercept=c(0,1,.4),mu_intercept_indicator=c(1,2,3),
    matrix(c(1,1,1,3,3,3),3))
  colnames(draws)[3:4] <- spec$allocations[[1L]]$gamma_coordinates
  weights <- c(mu_intercept=1,setNames(c(1,1),spec$coefficient_names))
  projection <- .bt_ordered_projection(list(mu_f=spec),weights,draws,info$prior_list)
  expect_identical(projection$atom,c(2.5,3.5,NA_real_))
  expect_identical(projection$state,c("point","point","continuous"))
  expect_identical(projection$values,c(2.5,3.5,2.9))
  adjacent <- info$prior_list
  adjacent$mu_intercept <- prior_mixture(list(prior("point",list(1)),prior("point",list(1+.Machine$double.eps))))
  adjacent_draws <- draws[1:2,,drop=FALSE]
  adjacent_draws[,"mu_intercept_indicator"] <- 1:2
  adjacent_draws[,"mu_intercept"] <- c(1,1+.Machine$double.eps)
  adjacent_projection <- .bt_ordered_projection(list(mu_f=spec),c(mu_intercept=1),adjacent_draws,adjacent)
  expect_identical(adjacent_projection$atom,c(1,1+.Machine$double.eps))
  expect_identical(adjacent_projection$values,adjacent_projection$atom)
  expression_prior <- p
  expression_prior$total <- prior("point",list(expression(tau)))
  expression_prior <- .bt_bind_ordered_prior_metadata(expression_prior, "mu_f")
  expression_spec <- .bt_ordered_spec("mu_f",expression_prior)
  expression_draws <- cbind(draws,mu_f_ordered_total=c(.2,.3,.4))
  guarded <- .bt_ordered_projection(list(mu_f=expression_spec),setNames(c(1,1),spec$coefficient_names),expression_draws)
  expect_identical(guarded$state,rep("unavailable",3))
  expect_s3_class(guarded$reason,"BayesTools_ordered_expression_unavailable")
  expect_error(JAGS_ordered_density_kernel(list(mu_f=expression_prior)),class="BayesTools_ordered_expression_unavailable")
})

test_that("ordered projections preserve slice events, shared unions and original-scale identities", {
  data <- data.frame(f=ordered(rep(c("lo","mid","hi","top"),2),levels=c("lo","mid","hi","top")),
    g=factor(rep(c("A","B","C","A"),2),levels=c("A","B","C")),x=seq_len(8))
  fixture <- function(formula,priors,scale=NULL){
    info <- JAGS_formula(formula,"mu",data,priors,formula_scale=scale)
    bound <- info$prior_list
    specs <- Map(.bt_ordered_spec,names(bound)[vapply(bound,is.prior.ordered,logical(1))],bound[vapply(bound,is.prior.ordered,logical(1))])
    draws <- matrix(numeric(),8,0,dimnames=list(NULL,character()))
    for(spec in specs){
      total <- matrix(seq_len(8)/10,8,length(spec$total_names))
      if(is.prior.point(spec$total_prior)) total[] <- spec$total_prior$parameters$location
      if(is.prior.spike_and_slab(spec$total_prior)){
        variable <- total + rep(seq_along(spec$total_names)-1L,each=8)
        indicator <- rep(c(0,1),4)
        colnames(variable) <- spec$total_node$spec$components[[1L]]$coordinates
        draws <- cbind(draws,variable,matrix(indicator,ncol=1,dimnames=list(NULL,spec$total_node$spec$indicator)))
        total <- variable * indicator
      }
      colnames(total) <- spec$total_names
      draws <- cbind(draws,total)
      for(record in spec$allocations){
        if(!identical(record$spec$type,"dirichlet") || all(record$gamma_coordinates %in% colnames(draws))) next
        gamma <- matrix(rep(seq_len(record$dim),each=8),8) * seq_len(8)
        colnames(gamma) <- record$gamma_coordinates
        draws <- cbind(draws,gamma)
      }
    }
    for(name in names(bound)){
      p <- bound[[name]]
      if(is.prior.ordered(p)){
        value <- .bt_deterministic_node_evaluate(.bt_dnode_ordered_coefficients(specs[[name]]),.bt_deterministic_lookup(draws))
      }else{
        columns <- if(is.prior.factor(p)) .JAGS_prior_factor_names(name,p) else name
        value <- matrix(if(is.prior.point(p)) p$parameters$location else .1,8,length(columns),dimnames=list(NULL,columns))
      }
      draws <- cbind(draws,value)
    }
    fit <- structure(list(mcmc=coda::mcmc.list(coda::mcmc(draws)),sample=8L,summary.pars=list(mutate=NULL),monitor=colnames(draws)),
      class=c("runjags","BayesTools_fit","list"))
    attr(fit,"prior_list") <- bound
    attr(fit,"formula_design") <- list(mu=info$formula_design)
    if(!is.null(scale)) attr(fit,"formula_scale") <- list(mu=info$formula_scale)
    fit <- attach_test_parameter_map(fit)
    list(fit=fit,info=info,specs=specs,draws=draws)
  }
  N01 <- prior("normal",list(0,1))
  priors <- list(intercept=prior("point",list(0)),f=prior_ordered(N01),
    g=prior_factor("normal",list(0,1),contrast="treatment"),
    `f:g`=prior_ordered(prior_spike_and_slab(N01,prior("point",list(.5)))))
  sliced <- fixture(~f*g,priors)
  name <- names(sliced$specs)[vapply(sliced$specs,function(spec) spec$metadata$theta_dim==2L,logical(1))]
  mixed <- as_mixed_posteriors(sliced$fit,name)[[name]]
  expect_identical(.posterior_atoms_get(mixed)$mass,.5)
  expect_true(all(.posterior_atoms_get(mixed)$locations==0))
  conditional <- as_mixed_posteriors(sliced$fit,name,conditional=name)[[name]]
  expect_identical(nrow(conditional),4L)
  expect_length(.posterior_atoms_get(conditional)$mass,0L)
  shared_priors <- priors
  shared_priors[["f:g"]] <- prior_ordered(N01,id="shape")
  shared <- fixture(~f*g,shared_priors)
  source_mix <- .mix_posteriors.factor(list(sliced$fit,shared$fit),
    list(sliced$info$prior_list[[name]],shared$info$prior_list[[name]]),name,c(.5,.5),seed=9,n_samples=20)
  display <- .bt_ordered_source_display(.bt_meta_get(source_mix,"ordered_source"))$samples
  expect_identical(ncol(display),8L)
  expect_identical(unname(display[,3:5,drop=FALSE]),unname(display[,6:8,drop=FALSE]))
  sampled <- JAGS_ordered_parameter_spec(shared$fit,name,draws=shared$draws)[[name]]$sampled_parameters
  expect_identical(length(sampled$allocation_names),3L)
  fixed_priors <- list(intercept=prior("point",list(1)),
    f=prior_ordered(prior("point",list(2.5)),allocation=c(1,2,3)/6),x=prior("point",list(0)),
    `f:x`=prior_ordered(prior("point",list(-2.5)),allocation=c(1,2,3)/6))
  scaled <- fixture(~f+x+f:x,fixed_priors,scale=list(x=TRUE))
  raw <- as_mixed_posteriors(scaled$fit,names(scaled$info$prior_list))
  original <- .transform_scale_samples_list(raw,attr(scaled$fit,"formula_scale",exact=TRUE))
  levels <- transform_factor_samples(original)$mu_f
  scale <- scaled$info$formula_scale[[1L]]
  expected <- 2.5 + 2.5 * scale$mean/scale$sd
  expect_true(all(as.numeric(levels[,4])==.posterior_atoms_for_column(.posterior_atoms_get(levels),4)$locations[1,1]))
  expect_equal(as.numeric(levels[1,4]),expected,tolerance=1e-14)
  public <- as_mixed_posteriors(scaled$fit,names(scaled$info$prior_list),transform_scaled=TRUE,n_prior_samples=128)
  public_levels <- transform_factor_samples(public)$mu_f
  expect_identical(as.numeric(public_levels[,4]),as.numeric(levels[,4]))
  expect_identical(.bt_meta_get(public$mu_f,"ordered_source")$primitives,.bt_meta_get(raw$mu_f,"ordered_source")$primitives)
  data$o <- data$f
  two <- fixture(~f*o,list(intercept=prior("point",list(0)),f=prior_ordered(N01),o=prior_ordered(N01),
    `f:o`=prior_ordered(prior("point",list(-2.5)),allocation=list(
      f=prior("dirichlet",list(alpha=c(2,3,4))),o=prior("dirichlet",list(alpha=c(3,2,1)))))))
  term <- names(two$specs)[vapply(two$specs,function(spec) length(spec$metadata$ordered_terms)==2L,logical(1))]
  spec <- two$specs[[term]]
  target <- setNames(as.numeric(spec$coefficient_grid$f==1L & spec$coefficient_grid$o<=2L),spec$coefficient_names)
  projection <- JAGS_ordered_parameter_spec(two$fit,term,weights=target,draws=two$draws)
  keys <- vapply(spec$allocations,`[[`,character(1),"factor")
  f_key <- names(keys)[keys=="f"]
  o_key <- names(keys)[keys=="o"]
  expect_identical(unname(projection$allocation_contractions[[f_key]][1,]),c(-1.25,0,0))
  expect_identical(unname(projection$allocation_contractions[[o_key]][1,]),c(-2.5*(1/6),-2.5*(1/6),0))
  expect_identical(projection$state,rep("continuous",8))
  full <- JAGS_ordered_parameter_spec(two$fit,term,weights=setNames(rep(1,9),spec$coefficient_names),draws=two$draws)
  expect_identical(full$atom,rep(-2.5,8))
  expect_identical(full$values,rep(-2.5,8))

  # Selected model rows and fitted scale maps agree with the public primitive
  # provider for every cell, including continuous intermediate interactions.
  canonical_cells <- function(fits,samples){
    parameter <- .bt_meta_get(samples,"ordered_source")$parameter
    values <- transform_factor_samples(setNames(list(samples),parameter))[[parameter]]
    source <- .bt_meta_get(values,"ordered_source")
    design <- source$projection_design
    expected <- matrix(NA_real_,nrow(values),ncol(values))
    for(model in unique(source$model)){
      rows <- which(source$model==model)
      for(column in seq_len(ncol(values))){
        projection <- JAGS_ordered_parameter_spec(fits[[model]],weights=design[column,],
          draws=as.matrix(fits[[model]]$mcmc))
        expected[rows,column] <- projection$values[source$draw_index[rows]]
      }
    }
    expect_identical(unname(.bt_draws_plain(values)),expected)
    expect_identical(source$primitives,.bt_meta_get(samples,"ordered_source")$primitives)
    expect_identical(source$model,.bt_meta_get(samples,"ordered_source")$model)
    expect_identical(source$draw_index,.bt_meta_get(samples,"ordered_source")$draw_index)
  }
  canonical_cells(list(sliced$fit),mixed)
  canonical_cells(list(sliced$fit),conditional)
  canonical_cells(list(sliced$fit,shared$fit),source_mix)
  canonical_cells(list(scaled$fit),public$mu_f)
  continuous_priors <- fixed_priors
  continuous_priors$f <- prior_ordered(N01,id="shared")
  continuous_priors[["f:x"]] <- prior_ordered(N01,id="shared")
  continuous_scaled <- fixture(~f+x+f:x,continuous_priors,scale=list(x=TRUE))
  continuous_original <- as_mixed_posteriors(continuous_scaled$fit,names(continuous_scaled$info$prior_list),
    transform_scaled=TRUE,n_prior_samples=128)
  canonical_cells(list(continuous_scaled$fit),continuous_original$mu_f)
  canonical_cells(list(two$fit),as_mixed_posteriors(two$fit,term)[[term]])
})

test_that("scaled ordered formula marginals reuse aligned fitted ordinary sources", {
  data <- data.frame(f=ordered(rep(c("early","middle","late","last"),2L),
    levels=c("early","middle","late","last")),x=c(10,12,17,21,11,13,18,22))
  info <- JAGS_formula(~f+x,"mu",data,list(intercept=prior("normal",list(0,1)),
    f=prior_ordered(prior("normal",list(0,1))),x=prior("normal",list(0,1))),
    formula_scale=list(x=TRUE))
  spec <- .bt_ordered_spec("mu_f",info$prior_list$mu_f)
  total <- c(2,-3,5,7)
  gamma <- rbind(c(1,2,3),c(4,1,2),c(3,5,1),c(2,3,7))
  beta <- cbind(mu_intercept=c(1,3,-2,4),mu_x=c(.5,-1,2,3))
  coefficients <- total * (gamma/rowSums(gamma))
  colnames(coefficients) <- spec$coefficient_names
  colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
  draws <- cbind(beta,coefficients,
    matrix(total,ncol=1L,dimnames=list(NULL,spec$total_names)),gamma)
  fit <- structure(list(mcmc=coda::mcmc.list(coda::mcmc(draws,start=11L,thin=2L)),
    sample=4L,summary.pars=list(mutate=NULL),monitor=colnames(draws)),
    class=c("runjags","BayesTools_fit","list"))
  attr(fit,"prior_list") <- info$prior_list
  attr(fit,"formula_design") <- list(mu=info$formula_design)
  attr(fit,"formula_scale") <- list(mu=info$formula_scale)
  fit <- attach_test_parameter_map(fit)
  fit_bytes <- serialize(fit,NULL)
  scale <- info$formula_scale$mu_x
  raw <- as_mixed_posteriors(fit,names(info$prior_list))
  original <- as_mixed_posteriors(fit,names(info$prior_list),transform_scaled=TRUE,n_prior_samples=128L)
  expected_original <- as.matrix(transform_scale_samples(fit))
  transform <- JAGS_formula_coefficient_transform(fit,"mu",target_scale="original")
  for(parameter in colnames(beta)){
    expect_identical(as.numeric(original[[parameter]]),as.numeric(expected_original[,parameter]))
    selected <- parameter_catalog_resolve(parameter_catalog(fit),parameter)
    expect_identical(as.numeric(as.matrix(parameter_draws(fit,selected))),as.numeric(beta[,parameter]))
  }
  problems <- expectation_problems({
    for(scaled in c(FALSE,TRUE)) for(x in c(12,19)){
      samples <- if(scaled) original else raw
      marginal <- marginal_posterior(samples,"mu_f",formula=~f+x,at=list(x=x),prior_samples=FALSE)
      fitted_x <- if(scaled) (x-scale$mean)/scale$sd else x
      for(level in seq_along(marginal)){
        share <- if(level==1L) numeric(4L) else if(level==4L) rep(1,4L) else{
          rowSums(gamma[,seq_len(level-1L),drop=FALSE])/rowSums(gamma)
        }
        independent <- beta[,"mu_intercept"] + fitted_x*beta[,"mu_x"] + total*share
        weights <- c(mu_intercept=1,mu_x=fitted_x,
          setNames(as.numeric(seq_len(3L)<level),spec$coefficient_names))
        # Both views describe this physical fitted row directly. Remapping
        # original coefficient weights through C changes its rounding.
        public <- JAGS_ordered_parameter_spec(fit,weights=weights,draws=draws)$values
        expect_equal(as.numeric(marginal[[level]]),independent,tolerance=5e-15)
        expect_identical(as.numeric(marginal[[level]]),public)
      }
    }
  })
  expect_identical(problems,character())
  expect_identical(serialize(fit,NULL),fit_bytes)
  source <- .bt_meta_get(original$mu_f,"ordered_source")
  expect_identical(source$projection_context$primitives[,colnames(beta),drop=FALSE],beta)
  expect_identical(source$draw_index,seq_len(4L))
  expect_identical(source$primitives,.bt_meta_get(raw$mu_f,"ordered_source")$primitives)

  # A clone with the fitted context removed cannot combine original ordinary
  # values with fitted formula weights. Pure ordered targets remain available.
  missing <- original
  for(parameter in names(missing)){
    own <- .bt_meta_get(missing[[parameter]],"ordered_source")
    own$projection_context <- NULL
    missing[[parameter]] <- .bt_meta_set(missing[[parameter]],"ordered_source",own)
  }
  expect_error(marginal_posterior(missing,"mu_f",formula=~f+x,at=list(x=12),prior_samples=FALSE),
    "Ordered formula projections are unavailable without retained fitted ordinary coefficients for scaled samples. Recreate mixed posteriors from the source fits with this version of BayesTools.",
    fixed=TRUE,class="BayesTools_ordered_coordinates_unavailable")
  pure_weights <- matrix(1,1,3,dimnames=list(NULL,spec$coefficient_names))
  pure <- .bt_ordered_formula_projections(missing,pure_weights)
  expect_identical(pure[[1L]]$values,total)
  ordinary_weights <- matrix(c(1,2),1,2,dimnames=list(NULL,colnames(beta)))
  ordinary <- .bt_ordered_formula_projections(original,ordinary_weights)
  expect_identical(ordinary[[1L]]$values,beta[,1L]+2*beta[,2L])
  expect_identical(ordinary[[1L]]$state,rep("continuous",4L))
  reordered <- original
  reordered$mu_x <- .bt_draws_subset_rows(reordered$mu_x,c(2L,1L,3L,4L))
  expect_error(.bt_ordered_formula_projections(reordered,ordinary_weights),
    "model and draw rows do not align",class="BayesTools_ordered_coordinates_unavailable")
  altered <- original
  changed <- .bt_meta_get(altered$mu_x,"ordered_source")
  changed$projection_context$primitives[1L,"mu_x"] <- 50
  altered$mu_x <- .bt_meta_set(altered$mu_x,"ordered_source",changed)
  expect_error(.bt_ordered_formula_projections(altered,ordinary_weights),
    "retained sources do not align",class="BayesTools_ordered_coordinates_unavailable")

  # The fitted context keeps the declared model union, including an absent
  # ordered term and a model from which no draws are selected.
  absent_info <- JAGS_formula(~x,"mu",data,list(intercept=prior("normal",list(0,1)),
    x=prior("normal",list(0,1))),formula_scale=list(x=TRUE))
  absent_fit <- fit
  absent_fit$mcmc <- coda::mcmc.list(coda::mcmc(beta,start=11L,thin=2L))
  absent_fit$monitor <- colnames(beta)
  attr(absent_fit,"prior_list") <- absent_info$prior_list
  attr(absent_fit,"formula_design") <- list(mu=absent_info$formula_design)
  attr(absent_fit,"formula_scale") <- list(mu=absent_info$formula_scale)
  absent_fit <- attach_test_parameter_map(absent_fit)
  models <- list(list(fit=fit,marglik=bridgesampling_object(log(.6)),prior_weights=1),
    list(fit=absent_fit,marglik=bridgesampling_object(log(.4)),prior_weights=1),
    list(fit=fit,marglik=bridgesampling_object(0),prior_weights=0))
  null <- setNames(rep(list(c(FALSE,FALSE,FALSE)),3L),names(info$prior_list))
  mixed <- mix_posteriors(models,names(info$prior_list),null,seed=9L,n_samples=20L)
  mixed <- .transform_scale_samples_list(mixed,attr(fit,"formula_scale",exact=TRUE))
  mixed_source <- .bt_meta_get(mixed$mu_f,"ordered_source")
  expect_length(mixed_source$models,3L)
  expect_false(any(mixed_source$model==3L))
  expect_true(all(c(1L,2L) %in% mixed_source$model))
  marginal <- marginal_posterior(mixed,"mu_f",formula=~f+x,at=list(x=12),prior_samples=FALSE)
  fitted_x <- (12-scale$mean)/scale$sd
  independent <- vapply(seq_along(marginal),function(level){
    index <- mixed_source$draw_index
    share <- if(level==1L) numeric(4L) else if(level==4L) rep(1,4L) else{
      rowSums(gamma[,seq_len(level-1L),drop=FALSE])/rowSums(gamma)
    }
    beta[index,1L] + fitted_x*beta[index,2L] +
      ifelse(mixed_source$model==1L,(total*share)[index],0)
  },numeric(20L))
  colnames(independent) <- names(marginal)
  expect_equal(vapply(marginal,as.numeric,numeric(20L)),independent,tolerance=5e-15)
  declared <- .bt_ordered_formula_projections(mixed,pure_weights)[[1L]]
  absent_rows <- mixed_source$model==2L
  expect_identical(declared$state,ifelse(absent_rows,"point","continuous"))
  expect_identical(declared$atom[absent_rows],rep(0,sum(absent_rows)))
  expect_identical(declared$values[absent_rows],rep(0,sum(absent_rows)))
  ordinary <- .bt_ordered_formula_projections(mixed,ordinary_weights)[[1L]]
  expect_identical(ordinary$state,rep("continuous",20L))
  expect_identical(ordinary$values,beta[mixed_source$draw_index,1L]+2*beta[mixed_source$draw_index,2L])
  subset <- mixed
  for(parameter in names(subset)) subset[[parameter]] <- .bt_draws_subset_rows(subset[[parameter]],c(7L,2L,10L))
  subset_marginal <- marginal_posterior(subset,"mu_f",formula=~f+x,at=list(x=12),prior_samples=FALSE)
  expect_identical(vapply(subset_marginal,as.numeric,numeric(3L)),
    vapply(marginal,as.numeric,numeric(20L))[c(7L,2L,10L),,drop=FALSE])
})

test_that("ordered formula projections use each contributing formula prefix", {
  data <- data.frame(f=ordered(rep(c("early","middle","late","last"),2L),
    levels=c("early","middle","late","last")),x=c(10,12,17,21,11,13,18,22),
    z=c(-5,-2,4,9,-3,0,6,12))
  normal <- function() prior("normal",list(0,1))
  mu <- JAGS_formula(~f+x,"mu",data,list(intercept=normal(),
    f=prior_ordered(normal()),x=normal()),formula_scale=list(x=TRUE))
  make_draws <- function(info, beta, total, gamma){
    parameter <- names(info$prior_list)[vapply(info$prior_list,is.prior.ordered,logical(1))]
    spec <- .bt_ordered_spec(parameter,info$prior_list[[parameter]])
    coefficients <- total*(gamma/rowSums(gamma))
    colnames(coefficients) <- spec$coefficient_names
    colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
    cbind(beta,coefficients,matrix(total,ncol=1L,dimnames=list(NULL,spec$total_names)),gamma)
  }
  mu_beta <- cbind(mu_intercept=c(1,3,-2,4),mu_x=c(.5,-1,2,3))
  mu_total <- c(2,-3,5,7)
  mu_gamma <- rbind(c(1,2,3),c(4,1,2),c(3,5,1),c(2,3,7))
  tau_beta <- cbind(tau_intercept=c(-7,5,11,-13),tau_z=c(-2,4,-1,.25))
  tau_total <- c(-11,13,-17,19)
  tau_gamma <- rbind(c(7,2,1),c(1,6,3),c(2,1,8),c(5,7,2))
  problems <- expectation_problems({
    for(tau_scaled in c(FALSE,TRUE)){
      tau <- if(tau_scaled){
        JAGS_formula(~f+z,"tau",data,list(intercept=normal(),
          f=prior_ordered(normal()),z=normal()),formula_scale=list(z=TRUE))
      }else JAGS_formula(~f,"tau",data,list(intercept=normal(),f=prior_ordered(normal())))
      tau_ordinary <- if(tau_scaled) tau_beta else tau_beta[,1L,drop=FALSE]
      draws <- cbind(make_draws(mu,mu_beta,mu_total,mu_gamma),
        make_draws(tau,tau_ordinary,tau_total,tau_gamma))
      fit <- structure(list(mcmc=coda::mcmc.list(coda::mcmc(draws,start=11L,thin=2L)),
        sample=4L,summary.pars=list(mutate=NULL),monitor=colnames(draws)),
        class=c("runjags","BayesTools_fit","list"))
      attr(fit,"prior_list") <- c(mu$prior_list,tau$prior_list)
      attr(fit,"formula_design") <- list(mu=mu$formula_design,tau=tau$formula_design)
      attr(fit,"formula_scale") <- if(tau_scaled) list(mu=mu$formula_scale,tau=tau$formula_scale) else list(mu=mu$formula_scale)
      fit <- attach_test_parameter_map(fit)
      fit_bytes <- serialize(fit,NULL)
      for(scaled in c(FALSE,TRUE)){
        samples <- as_mixed_posteriors(fit,names(attr(fit,"prior_list")),
          transform_scaled=scaled,n_prior_samples=128L)
        expect_identical(.bt_meta_get(samples$mu_f,"ordered_source")$draw_index,1:4)
        expect_identical(.bt_meta_get(samples$tau_f,"ordered_source")$draw_index,1:4)
        if(!tau_scaled){
          own <- .bt_meta_get(samples$tau_f,"ordered_source")
          identity <- diag(3L)
          dimnames(identity) <- list(colnames(samples$tau_f),own$models[[1L]]$coefficient_names)
          expect_identical(own$projection_design,identity)
          expect_identical(as.numeric(samples$tau_f),as.numeric(tau_total*(tau_gamma/rowSums(tau_gamma))))
          expect_true(!is.null(.posterior_atoms_get(samples$tau_f)$marginals))
        }
        expected <- list()
        for(prefix in c("mu","tau")){
          info <- if(prefix=="mu") mu else tau
          beta <- if(prefix=="mu") mu_beta else tau_ordinary
          total <- if(prefix=="mu") mu_total else tau_total
          gamma <- if(prefix=="mu") mu_gamma else tau_gamma
          value <- if(prefix=="mu") 12 else 6
          slope <- if(prefix=="mu") "mu_x" else "tau_z"
          predictor <- if(scaled && !is.null(info$formula_scale)){
            (value-info$formula_scale[[slope]]$mean)/info$formula_scale[[slope]]$sd
          }else value
          expected[[prefix]] <- vapply(1:4,function(level){
            share <- if(level==1L) numeric(4L) else if(level==4L) rep(1,4L) else{
              rowSums(gamma[,seq_len(level-1L),drop=FALSE])/rowSums(gamma)
            }
            beta[,1L] + (if(ncol(beta)==2L) predictor*beta[,2L] else 0) + total*share
          },numeric(4L))
          formula <- if(prefix=="mu") ~f+x else if(tau_scaled) ~f+z else ~f
          at <- if(prefix=="mu") list(x=12) else if(tau_scaled) list(z=6) else NULL
          result <- marginal_posterior(samples,paste0(prefix,"_f"),formula=formula,at=at,prior_samples=FALSE)
          expect_equal_each(as.numeric(vapply(result,as.numeric,numeric(4L))),as.numeric(expected[[prefix]]),tolerance=5e-15)
          subset <- samples[names(info$prior_list)]
          result <- marginal_posterior(subset,paste0(prefix,"_f"),formula=formula,at=at,prior_samples=FALSE)
          expect_equal_each(as.numeric(vapply(result,as.numeric,numeric(4L))),as.numeric(expected[[prefix]]),tolerance=5e-15)
        }
        # An existing combined target uses both prefixes, including independent
        # ordinary draws and allocations. Retained roots supply its last levels.
        columns <- unique(unlist(lapply(names(samples),function(parameter){
          .posterior_atoms_coefficient_columns(samples[[parameter]],parameter)
        }),use.names=FALSE))
        weights <- matrix(0,2L,length(columns),dimnames=list(NULL,columns))
        for(prefix in c("mu","tau")){
          info <- if(prefix=="mu") mu else tau
          spec <- .bt_ordered_spec(paste0(prefix,"_f"),info$prior_list[[paste0(prefix,"_f")]])
          weights[,paste0(prefix,"_intercept")] <- 1
          weights[1L,spec$coefficient_names[1:2]] <- 1
          weights[2L,spec$coefficient_names] <- 1
          if(prefix=="mu" || tau_scaled){
            slope <- if(prefix=="mu") "mu_x" else "tau_z"
            value <- if(prefix=="mu") 12 else 6
            weights[,slope] <- if(scaled) (value-info$formula_scale[[slope]]$mean)/info$formula_scale[[slope]]$sd else value
          }
        }
        combined <- .bt_ordered_formula_projections(samples,weights)
        expect_equal_each(as.numeric(do.call(cbind,lapply(combined,`[[`,"values"))),
          as.numeric(expected$mu[,3:4]+expected$tau[,3:4]),tolerance=5e-15)
        unrelated <- samples
        unrelated$mu_f <- posterior_transform(unrelated$mu_f,"exp")
        unrelated$mu_x <- .bt_draws_subset_rows(unrelated$mu_x,c(2L,1L,3L,4L))
        tau_weights <- weights
        mu_columns <- unique(unlist(lapply(names(mu$prior_list),function(parameter){
          .posterior_atoms_coefficient_columns(samples[[parameter]],parameter)
        }),use.names=FALSE))
        tau_weights[,mu_columns] <- 0
        if(!tau_scaled && scaled){
          # An unscaled formula emits NULL, even if its prefix is represented
          # explicitly in a caller's retained scale list.
          scales <- .bt_meta_get(unrelated,"formula_scale")
          unrelated <- .bt_meta_set(unrelated,"formula_scale",c(scales,list(tau=NULL)))
        }
        tau_projection <- .bt_ordered_formula_projections(unrelated,tau_weights)
        expect_equal_each(as.numeric(do.call(cbind,lapply(tau_projection,`[[`,"values"))),
          as.numeric(expected$tau[,3:4]),tolerance=5e-15)
        expect_identical(as.numeric(unrelated$mu_x),as.numeric(samples$mu_x)[c(2L,1L,3L,4L)])
        expect_identical(attr(tau_projection,"model",exact=TRUE),rep(1L,4L))
        if(scaled && tau_scaled){
          for(field in c("specs","priors")){
            inconsistent <- samples
            mu_context <- .bt_meta_get(samples$mu_f,"ordered_source")$projection_context
            for(parameter in names(tau$prior_list)){
              source <- .bt_meta_get(inconsistent[[parameter]],"ordered_source")
              if(field=="specs"){
                overlap <- mu_context$models[[1L]]$specs$mu_f
                overlap$total_names <- "contradictory_total"
                source$projection_context$models[[1L]]$specs$mu_f <- overlap
              }else{
                source$projection_context$models[[1L]]$priors$mu_intercept <- prior("normal",list(2,1))
              }
              inconsistent[[parameter]] <- .bt_meta_set(inconsistent[[parameter]],"ordered_source",source)
            }
            expect_error(.bt_ordered_formula_projections(inconsistent,weights),
              "overlapping sources do not agree",class="BayesTools_ordered_coordinates_unavailable")
          }
          # Every tau owner agrees with its own context, but that context now
          # contradicts an overlapping mu source. The union must refuse it.
          for(parameter in names(tau$prior_list)){
            source <- .bt_meta_get(samples[[parameter]],"ordered_source")
            source$projection_context$primitives <- cbind(source$projection_context$primitives,mu_intercept=mu_beta[,1L]+1)
            samples[[parameter]] <- .bt_meta_set(samples[[parameter]],"ordered_source",source)
          }
          expect_error(.bt_ordered_formula_projections(samples,weights),
            "overlapping sources do not agree",class="BayesTools_ordered_coordinates_unavailable")
        }
      }
      expect_identical(serialize(fit,NULL),fit_bytes)
    }
  })
  expect_identical(problems,character())
})

test_that("public no-intercept ordered reference has exact zero values and a unit atom", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))
  levels <- marginal_posterior(fixture$samples, "mu_f", formula = ~0 + f,
    prior_samples = FALSE)
  reference <- levels[[1L]]
  expect_identical(as.numeric(reference), rep(0, 120L))
  atoms <- .posterior_atoms_get(reference)
  expect_identical(as.numeric(atoms$locations), 0)
  expect_identical(atoms$mass, 1)
})

test_that("unscaled ordered producers always restore raw primitive semantics", {
  data <- data.frame(f=ordered(rep(c("early","middle","late","last"),2L),
    levels=c("early","middle","late","last")))
  normal <- function() prior("normal",list(0,1))
  cases <- list(
    point0=list(total_prior=prior("point",list(0)),total=rep(0,4L),allocation=NULL),
    fixed_point=list(total_prior=prior("point",list(2.5)),total=rep(2.5,4L),allocation=c(1,2,3)/6),
    fixed_normal_zero=list(total_prior=normal(),total=c(2,-3,5,7),allocation=c(0,.5,.5)),
    normal=list(total_prior=normal(),total=c(2,-3,5,7),allocation=NULL))
  problems <- expectation_problems({
    for(case in cases){
      info <- JAGS_formula(~f,"tau",data,list(intercept=prior("point",list(0)),
        f=prior_ordered(case$total_prior,allocation=case$allocation)))
      spec <- .bt_ordered_spec("tau_f",info$prior_list$tau_f)
      gamma <- rbind(c(1,2,3),c(4,1,2),c(3,5,1),c(2,3,7))
      shares <- if(is.null(case$allocation)) gamma/rowSums(gamma) else matrix(rep(case$allocation,each=4L),4L)
      expected <- case$total*shares
      colnames(expected) <- spec$coefficient_names
      identity <- diag(3L)
      primitives <- matrix(case$total,ncol=1L,dimnames=list(NULL,spec$total_names))
      if(is.null(case$allocation)){
        colnames(gamma) <- spec$allocations[[1L]]$gamma_coordinates
        primitives <- cbind(primitives,gamma)
      }
      draws <- cbind(tau_intercept=rep(0,4L),expected,primitives)
      fit <- structure(list(mcmc=coda::mcmc.list(coda::mcmc(draws,start=11L,thin=2L)),
        sample=4L,summary.pars=list(mutate=NULL),monitor=colnames(draws)),
        class=c("runjags","BayesTools_fit","list"))
      attr(fit,"prior_list") <- info$prior_list
      attr(fit,"formula_design") <- list(tau=info$formula_design)
      fit <- attach_test_parameter_map(fit)
      fit_bytes <- serialize(fit,NULL)
      for(stale in c(FALSE,TRUE)){
        source_fit <- fit
        if(stale){
          # Only this controlled scratch clone has stale backend increments;
          # its declared total and allocation primitives remain unchanged.
          stale_draws <- draws
          stale_draws[,spec$coefficient_names] <- stale_draws[,spec$coefficient_names,drop=FALSE]+1
          source_fit$mcmc <- coda::mcmc.list(coda::mcmc(stale_draws,start=11L,thin=2L))
        }
        source_bytes <- serialize(source_fit,NULL)
        outputs <- lapply(c(FALSE,TRUE),function(scaled){
          as_mixed_posteriors(source_fit,"tau_f",transform_scaled=scaled,n_prior_samples=64L)$tau_f
        })
        for(output in outputs){
          expect_equal_each(as.numeric(output),as.numeric(expected),tolerance=5e-15)
          source <- .bt_meta_get(output,"ordered_source")
          dimnames(identity) <- list(colnames(output),spec$coefficient_names)
          expect_identical(source$projection_design,identity)
          expect_identical(source$draw_index,1:4)
          expect_identical(source$model,rep(1L,4L))
          atoms <- .posterior_atoms_get(output)
          if(is.prior.point(case$total_prior)){
            expect_identical(unname(vapply(atoms$marginals,function(atom) atom$mass,numeric(1))),rep(1,3L))
            expect_equal_each(vapply(atoms$marginals,function(atom) atom$locations[1L,1L],numeric(1)),
              as.numeric(expected[1L,]),tolerance=5e-15)
          }else if(!is.null(case$allocation)){
            expect_identical(atoms$marginals[[1L]]$locations,matrix(0,1L,1L,dimnames=list(NULL,colnames(output)[1L])))
            expect_identical(atoms$marginals[[1L]]$mass,1)
          }
        }
        expect_identical(.bt_meta_get(outputs[[1L]],"ordered_source"),.bt_meta_get(outputs[[2L]],"ordered_source"))
        expect_identical(.posterior_atoms_get(outputs[[1L]]),.posterior_atoms_get(outputs[[2L]]))
        expect_identical(serialize(source_fit,NULL),source_bytes)
      }
      expect_identical(serialize(fit,NULL),fit_bytes)
    }
  })
  expect_identical(problems,character())
})

test_that("ordered formula weights are validated before source selection", {
  weights <- matrix(0,1L,3L,dimnames=list(NULL,c("a","b","c")))
  expect_null(.bt_ordered_formula_projections(list(),weights))
  for(value in c(NA_real_,NaN)){
    invalid <- weights
    invalid[,] <- value
    expect_error(.bt_ordered_formula_projections(list(),invalid),"cannot contain NA/NaN values",fixed=TRUE)
  }
  for(value in c(Inf,-Inf)){
    invalid <- weights
    invalid[,] <- value
    expect_error(.bt_ordered_formula_projections(list(),invalid),"must be finite named fitted-coordinate weights",fixed=TRUE)
  }
  colnames(weights) <- NULL
  expect_error(.bt_ordered_formula_projections(list(),weights),"must be finite named fitted-coordinate weights",fixed=TRUE)
  colnames(weights) <- c("a","a","c")
  expect_error(.bt_ordered_formula_projections(list(),weights),"must be finite named fitted-coordinate weights",fixed=TRUE)
})

test_that("prior_ordered() validates constructor inputs", {
  p <- prior_ordered(prior("normal", list(0, 1)))

  expect_true(is.prior(p))
  expect_true(is.prior.factor(p))
  expect_true(is.prior.ordered(p))
  expect_equal(p$contrast, "cumulative")
  expect_equal(p$allocation$type, "default_dirichlet")

  expect_error(
    prior_ordered(prior_factor("normal", list(0, 1), contrast = "treatment")),
    "scalar prior"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .2)),
    "sum to one"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), id = ""),
    "cannot contain empty strings"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), id = "  "),
    "cannot contain empty strings"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), allocation = prior("normal", list(0, 1))),
    "Dirichlet"
  )
  expect_error(
    prior_spike_and_slab(p),
    "inside prior_ordered"
  )
})

test_that("ordered cumulative contrasts encode level effects", {
  expect_equal(
    contr.ordered_cumulative(c("low", "mid", "high")),
    matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    contr.ordered_cumulative_levels(c("low", "mid", "high")),
    matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  )
})

test_that("formula binding converts explicit ordered priors with current factor rules", {
  df <- data.frame(
    y  = seq_len(6),
    f  = factor(rep(c("mid", "low", "high"), 2), levels = c("mid", "low", "high")),
    ch = rep(c("b", "a", "c"), 2)
  )

  formula_info <- JAGS_formula(
    y ~ f + ch,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f  = prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .8)),
      ch = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )

  expect_equal(formula_info$formula_design$xlevels$f, c("mid", "low", "high"))
  expect_equal(formula_info$formula_design$xlevels$ch, c("a", "b", "c"))
  expect_equal(formula_info$formula_design$contrasts$f, "contr.ordered_cumulative")
  expect_equal(formula_info$formula_design$contrasts$ch, "contr.ordered_cumulative")
  expect_equal(unname(formula_info$data$mu_data_f[1:3, ]), matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = TRUE))

  expect_error(
    JAGS_formula(
      y ~ f,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .3, .5))
      )
    ),
    "has length 3, but 2 value"
  )
})

test_that("cumulative_levels uses one allocation piece per level", {
  p <- prior_ordered(
    prior("point", list(location = 10)),
    allocation = c(.2, .3, .5),
    contrast = "cumulative_levels"
  )
  attr(p, "levels") <- 3
  attr(p, "level_names") <- c("low", "mid", "high")

  samples <- rng(p, 2)
  expect_equal(unname(samples[1, ]), c(2, 5, 10))

  df <- data.frame(
    y = seq_len(6),
    f = factor(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ 0 + f,
    "mu",
    data = df,
    prior_list = list(
      f = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = c(.2, .3, .5),
        contrast = "cumulative_levels"
      )
    )
  )

  expect_equal(unname(formula_info$data$mu_data_f[1:3, ]), contr.ordered_cumulative_levels(1:3))
})

test_that("JAGS syntax, inits, and monitors use latent total and allocation nodes", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_fixed <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
    )
  )
  syntax_fixed <- JAGS_add_priors("model{}", formula_fixed$prior_list)

  expect_match(syntax_fixed, "mu_f_ordered_total ~ dnorm\\(0,1\\)")
  expect_match(syntax_fixed, "mu_f\\[1\\] <- mu_f_ordered_total \\* 0.25")
  expect_false(grepl("prior_par_eta_mu_f_ordered_alloc", syntax_fixed))

  formula_dirichlet <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    )
  )
  syntax_dirichlet <- JAGS_add_priors("model{}", formula_dirichlet$prior_list)
  expect_match(syntax_dirichlet, "prior_par_eta_mu_f_ordered_alloc_f_1\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax_dirichlet, "mu_f\\[2\\] <- mu_f_ordered_total \\* mu_f_ordered_alloc_f_1\\[2\\]")

  monitors <- JAGS_to_monitor(formula_dirichlet$prior_list)
  expect_true(all(c("mu_f", "mu_f_ordered_total", "prior_par_eta_mu_f_ordered_alloc_f_1") %in% monitors))

  inits <- JAGS_get_inits(formula_dirichlet$prior_list, chains = 1, seed = 1)[[1]]
  expect_true("mu_f_ordered_total" %in% names(inits))
  expect_true("prior_par_eta_mu_f_ordered_alloc_f_1" %in% names(inits))
  expect_equal(length(inits$prior_par_eta_mu_f_ordered_alloc_f_1), 2L)

  formula_spike_slab <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(
        prior_spike_and_slab(prior("normal", list(0, 1))),
        allocation = c(.25, .75)
      )
    )
  )
  syntax_spike_slab <- JAGS_add_priors("model{}", formula_spike_slab$prior_list)
  expect_match(syntax_spike_slab, "mu_f_ordered_total_indicator ~ dbern")
  expect_match(syntax_spike_slab, "mu_f_ordered_total = mu_f_ordered_total_variable \\* mu_f_ordered_total_indicator")
})

test_that("ordered interactions expand through the formula binder", {
  df <- data.frame(
    y = seq_len(12),
    x = rep(c(-1, 1), 6),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    g = factor(rep(c("A", "B"), each = 6))
  )

  formula_info <- JAGS_formula(
    y ~ f + x + g + f:x + f:g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      x = prior("normal", list(0, 1)),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:x" = prior_ordered(prior("normal", list(0, 1))),
      "f:g" = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = list(f = prior("dirichlet", list(alpha = c(2, 3))))
      )
    )
  )

  expect_equal(dim(attr(formula_info$prior_list$mu_f__xXx__x, "factor_design")), c(3L, 2L))
  expect_equal(dim(attr(formula_info$prior_list$mu_f__xXx__g, "factor_design")), c(6L, 2L))

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_match(syntax, "mu_f__xXx__x\\[1\\] <- mu_f__xXx__x_ordered_total")
  expect_match(syntax, "prior_par_eta_mu_f__xXx__g_ordered_alloc_f_1\\[2\\] ~ dgamma\\(3, 1\\)")

  expect_error(
    JAGS_formula(
      y ~ f:g,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "f:g" = prior_ordered(prior("normal", list(0, 1)))
      )
    ),
    "by level indicators", fixed = TRUE
  )
})

test_that("ordered allocation id sharing emits one shared allocation", {
  df <- data.frame(
    y = seq_len(6),
    x = rep(c(-1, 1), 3),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )

  formula_info <- JAGS_formula(
    y ~ f + x + f:x,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape"),
      x = prior("normal", list(0, 1)),
      "f:x" = prior_ordered(prior("normal", list(0, 1)), id = "shape")
    )
  )

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_equal(lengths(regmatches(syntax, gregexpr("prior_par_eta_ordered_alloc_shape_f[1] ~", syntax, fixed = TRUE))), 1L)
  expect_match(syntax, "mu_f__xXx__x\\[1\\] <- mu_f__xXx__x_ordered_total \\* ordered_alloc_shape_f\\[1\\]")

  monitors <- JAGS_to_monitor(formula_info$prior_list)
  expect_equal(sum(monitors == "prior_par_eta_ordered_alloc_shape_f"), 1L)

  formula_info_other <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape")
    )
  )
  combined_priors <- c(formula_info$prior_list, formula_info_other$prior_list)
  combined_syntax <- JAGS_add_priors("model{}", combined_priors)
  expect_equal(
    lengths(regmatches(
      combined_syntax,
      gregexpr("prior_par_eta_ordered_alloc_shape_f[1] ~", combined_syntax, fixed = TRUE)
    )),
    1L
  )

  expect_error(
    JAGS_formula(
      y ~ f + x + f:x,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(
          prior("normal", list(0, 1)),
          allocation = c(.25, .75),
          id = "shape"
        ),
        x = prior("normal", list(0, 1)),
        "f:x" = prior_ordered(
          prior("normal", list(0, 1)),
          allocation = c(.5, .5),
          id = "shape"
        )
      )
    ),
    "incompatible allocation specifications"
  )

  formula_info_incompatible <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = prior("dirichlet", list(alpha = c(2, 3))),
        id = "shape"
      )
    )
  )
  incompatible_priors <- c(formula_info$prior_list, formula_info_incompatible$prior_list)
  expect_error(
    JAGS_add_priors("model{}", incompatible_priors),
    "incompatible allocation specifications"
  )
  expect_error(
    JAGS_get_inits(incompatible_priors, chains = 1, seed = 1),
    "incompatible allocation specifications"
  )
  expect_error(
    JAGS_to_monitor(incompatible_priors),
    "incompatible allocation specifications"
  )

  formula_info_collision_a <- JAGS_formula(
    y ~ f,
    "phi",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape-a")
    )
  )
  formula_info_collision_b <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape a")
    )
  )
  colliding_priors <- c(
    formula_info_collision_a$prior_list,
    formula_info_collision_b$prior_list
  )
  expect_error(
    JAGS_add_priors("model{}", colliding_priors),
    "JAGS node 'ordered_alloc_shape_a_f' conflicts"
  )
})

test_that("ordered hidden total nodes cannot collide with formula coefficients", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(
      rep(c("low", "mid", "high"), 2),
      levels = c("low", "mid", "high")
    ),
    f_ordered_total = rep(c(-1, 1), 3)
  )

  expect_error(
    JAGS_formula(
      y ~ f + f_ordered_total,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior("normal", list(0, 1))),
        f_ordered_total = prior("normal", list(0, 1))
      )
    ),
    "JAGS node 'mu_f_ordered_total' conflicts"
  )
})

test_that("multi-slice ordered expression totals omit initialization", {
  df <- expand.grid(
    f = ordered(
      c("low", "mid", "high"),
      levels = c("low", "mid", "high")
    ),
    g = factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  )
  formula_info <- JAGS_formula(
    ~ f * g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_ordered(
        prior("normal", list(0, expression(sigma)))
      )
    )
  )

  expect_equal(
    attr(formula_info$prior_list$mu_f__xXx__g, "ordered_metadata")$theta_dim,
    2L
  )
  inits <- JAGS_get_inits(
    formula_info$prior_list,
    chains = 1,
    seed = 1
  )[[1L]]
  expect_false("mu_f__xXx__g_ordered_total" %in% names(inits))
  expect_true(any(grepl(
    "prior_par_eta_mu_f__xXx__g_ordered_alloc",
    names(inits),
    fixed = TRUE
  )))
})

test_that("multi-slice spike-and-slab ordered totals skip expression inits", {
  df <- expand.grid(
    f = ordered(
      c("low", "mid", "high"),
      levels = c("low", "mid", "high")
    ),
    g = factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  )
  formula_info <- JAGS_formula(
    ~ f * g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_ordered(prior_spike_and_slab(
        prior("normal", list(0, expression(sigma))),
        prior_inclusion = prior("beta", list(1, 1))
      ))
    )
  )
  prior_list <- c(
    formula_info$prior_list,
    list(sigma = prior("gamma", list(2, 2)))
  )

  expect_equal(
    attr(prior_list$mu_f__xXx__g, "ordered_metadata")$theta_dim,
    2L
  )
  expect_match(
    JAGS_add_priors("model{}", prior_list),
    "mu_f__xXx__g_ordered_total_variable[2]",
    fixed = TRUE
  )
  inits <- JAGS_get_inits(prior_list, chains = 1, seed = 1)[[1L]]
  # The expression slab is initialized by JAGS from its parent 'sigma'; the
  # sampled inclusion probability keeps its own initial value.
  expect_false("mu_f__xXx__g_ordered_total_variable" %in% names(inits))
  expect_true("mu_f__xXx__g_ordered_total_inclusion" %in% names(inits))
  expect_true("sigma" %in% names(inits))
})

test_that("the slices of a spike-and-slab ordered total share one inclusion indicator", {

  # The fitted model draws one inclusion probability and one indicator for a
  # spike-and-slab total and multiplies the slab of every theta slice by it
  # (all slices are included or all are excluded).
  df <- expand.grid(
    f = ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")),
    g = factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  )
  formula_info <- JAGS_formula(
    ~ f * g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_ordered(prior_spike_and_slab(
        prior("normal", list(0, 1)),
        prior_inclusion = prior("beta", list(2, 3))
      ))
    )
  )
  interaction <- formula_info$prior_list$mu_f__xXx__g
  metadata <- attr(interaction, "ordered_metadata")
  expect_identical(metadata$theta_dim, 2L)
  expect_match(
    JAGS_add_priors("model{}", formula_info$prior_list),
    "mu_f__xXx__g_ordered_total[2] <- mu_f__xXx__g_ordered_total_variable[2] * mu_f__xXx__g_ordered_total_indicator",
    fixed = TRUE
  )

  n <- 100000
  set.seed(1)
  totals <- rng(interaction, n, quantity = "total")
  included <- totals != 0
  # the marginal inclusion probability E[Beta(2, 3)] = 0.4 and its binomial
  # standard error for n draws
  p <- 2 / 5
  se <- sqrt(p * (1 - p) / n)
  expect_identical(unname(included[, 1L]), unname(included[, 2L]))
  expect_lt(abs(mean(included[, 1L]) - p), 4 * se)
  # jointly, all slices are in with probability p and out with probability
  # 1 - p (independent slices would give p^2 = 0.16, 0.36, and a mixed
  # pattern with probability 2 p (1 - p) = 0.48)
  expect_lt(abs(mean(rowSums(included) == 2L) - p), 4 * se)
  expect_lt(abs(mean(rowSums(included) == 0L) - (1 - p)), 4 * se)
  expect_identical(sum(rowSums(included) == 1L), 0L)
  # one total component per draw
  component <- .bt_meta_get(totals, "ordered_total_component")
  expect_identical(length(component), as.integer(n))
  expect_identical(component == which(attr(interaction$total, "components") == "alternative"),
                   unname(included[, 1L]))

  # every coefficient of a draw is in or out with the shared indicator
  set.seed(2)
  coefficients <- rng(interaction, n, transform_factor_samples = FALSE)
  coefficient_included <- rowSums(unclass(coefficients) != 0)
  expect_true(all(coefficient_included %in% c(0L, metadata$coefficient_dim)))
  expect_lt(abs(mean(coefficient_included > 0L) - p), 4 * se)
})

test_that("ordered random slope contrasts are specified independently", {
  df <- data.frame(
    y = seq_len(12),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    id = factor(rep(seq_len(4), each = 3))
  )

  formula_info <- JAGS_formula(
    ~ 1 + f + (f || id),
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        contrasts = c(f = "cumulative")
      )
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(
    random_term$sd_leaves$leaf_names,
    c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_f[1]",
      "mu__xREx__id_f[2]"
    )
  )
  expect_equal(
    random_term$sd_leaves$leaf_names_by_column,
    c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_f[1]",
      "mu__xREx__id_f[2]"
    )
  )
  expect_equal(
    unname(formula_info$data$mu__xREx__id_xRE_DATAx[1:3, ]),
    matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  )
})

test_that("prior_ordered() can define ordered random slope SD components", {
  df <- data.frame(
    y = seq_len(12),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    id = factor(rep(seq_len(4), each = 3))
  )

  formula_info <- JAGS_formula(
    ~ 1 + (0 + f || id),
    "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior_ordered(
          prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
          allocation = c(.4, .6)
        )
      )
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(random_term$sd_parameter_names, c("mu__xREx__id_f[1]", "mu__xREx__id_f[2]"))
  expect_s3_class(formula_info$prior_list$mu__xREx__id_f, "prior.ordered")

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_match(syntax, "mu__xREx__id_f_ordered_total ~ dnorm\\(0,1\\)T\\(0,\\)")
  expect_match(syntax, "mu__xREx__id_f\\[1\\] <- mu__xREx__id_f_ordered_total \\* 0.4")
  expect_match(syntax, "mu__xREx__id_f\\[2\\] <- mu__xREx__id_f_ordered_total \\* 0.59999999999999998")
})

test_that("prior_ordered() indexes a sole two-level random slope SD", {
  df <- data.frame(
    y = seq_len(8),
    f = ordered(rep(c("low", "high"), 4), levels = c("low", "high")),
    id = factor(rep(seq_len(4), each = 2))
  )

  formula_info <- JAGS_formula(
    ~ 1 + (0 + f || id),
    "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior_ordered(
          prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
          allocation = 1
        )
      )
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(random_term$n_columns, 1L)
  expect_equal(random_term$sd_parameter_names, "mu__xREx__id_f")
  expect_s3_class(formula_info$prior_list$mu__xREx__id_f, "prior.ordered")

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_match(syntax, "mu__xREx__id_f_ordered_total ~ dnorm\\(0,1\\)T\\(0,\\)")
  expect_match(syntax, "mu__xREx__id_f <- mu__xREx__id_f_ordered_total \\* 1")
})

test_that("ordered posterior extraction transforms coefficients to public levels", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
    )
  )
  samples <- matrix(
    c(1, 2, 4, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )

  transformed <- BayesTools:::.transform_factor_contrasts(
    samples,
    formula_info$prior_list,
    transform_factors = TRUE
  )

  # the level effects, labelled by their cells, without the reference level
  # that the cumulative contrast fixes at zero
  expect_equal(colnames(transformed), c("mu_f[mid]", "mu_f[high]"))
  expect_equal(unname(transformed[1, ]), c(1, 3))
  expect_equal(unname(transformed[2, ]), c(4, 12))
})

test_that("prior sample generation uses stored ordered-mixture dimensions", {
  make_prior <- function(total, prior_weights){
    p <- prior_ordered(
      prior("point", list(location = total)),
      allocation = c(.2, .3, .5),
      contrast = "cumulative_levels",
      prior_weights = prior_weights
    )
    attr(p, "levels") <- 3
    attr(p, "level_names") <- c("low", "mid", "high")
    p
  }

  p1 <- make_prior(10, 1)
  p2 <- make_prior(20, 3)

  single <- BayesTools:::.generate_prior_sample_matrix(
    prior_list = list(mu_f = p1),
    n_samples  = 4,
    seed       = 1
  )
  expect_equal(
    unname(single[, paste0("mu_f[", 1:3, "]"), drop = FALSE]),
    matrix(c(2, 3, 5), nrow = 4, ncol = 3, byrow = TRUE)
  )
  # the coefficients and the fitted total node
  expect_equal(colnames(single), c(paste0("mu_f[", 1:3, "]"), "mu_f_ordered_total"))
  expect_equal(unname(single[, "mu_f_ordered_total"]), rep(10, 4))

  mixture_prior <- prior_mixture(list(p1, p2))
  mixed <- BayesTools:::.generate_prior_sample_matrix(
    prior_list = list(mu_f = mixture_prior),
    n_samples  = 32,
    seed       = 2
  )
  # the coefficients and the mixture's component indicator
  expect_equal(colnames(mixed), c(paste0("mu_f[", 1:3, "]"), "mu_f_indicator"))
  row_values <- apply(mixed[, paste0("mu_f[", 1:3, "]")], 1L, paste, collapse = ",")
  expect_setequal(
    unique(row_values),
    c("2,3,5", "4,6,10")
  )
  expect_true(all(row_values[mixed[, "mu_f_indicator"] == 1] == "2,3,5"))
  expect_true(all(row_values[mixed[, "mu_f_indicator"] == 2] == "4,6,10"))
})

test_that("fixed ordered allocations canonicalize only roundoff drift", {

  allocation <- c(.2, .3, .5 + .Machine$double.eps)
  spec <- BayesTools:::.prior_ordered_allocation_spec(allocation)

  expect_equal(sum(spec$weights), 1)
  expect_equal(
    spec$canonicalization$original_sum,
    sum(allocation)
  )
  expect_lte(
    spec$canonicalization$max_correction,
    spec$canonicalization$roundoff_bound
  )

  expect_error(
    BayesTools:::.prior_ordered_allocation_spec(c(.2, .3, .5 + 1e-10)),
    "exceeding the roundoff bound"
  )
})

test_that("public posterior mixing projects ordered primitive rows and preserves metadata", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
    )
  )
  ordered_prior <- formula_info$prior_list$mu_f
  posterior <- matrix(
    c(1, 2, 4, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )
  source_posterior <- cbind(posterior, mu_f_ordered_total = c(3, 12))

  fit <- coda::mcmc(source_posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  fit <- attach_test_parameter_map(fit)
  single <- as_mixed_posteriors(fit, parameters = "mu_f")

  expect_identical(unname(.bt_draws_plain(single$mu_f)),cbind(c(3,12)*.25,c(3,12)*.75))
  expect_true(isTRUE(attr(single$mu_f, "ordered")))
  expect_equal(
    attr(single$mu_f, "ordered_metadata"),
    attr(ordered_prior, "ordered_metadata")
  )

  single_levels <- transform_factor_samples(single)$mu_f
  expect_equal(
    unname(single_levels[, , drop = FALSE]),
    matrix(c(0, .75, 3, 0, 3, 12), nrow = 2, byrow = TRUE)
  )
  expect_equal(
    BayesTools:::.prior_factor_level_weight_matrix(single_levels, "mu_f"),
    structure(
      attr(ordered_prior, "factor_design"),
      dimnames = list(
        colnames(single_levels),
        c("mu_f[1]", "mu_f[2]")
      )
    )
  )

  make_model <- function(samples, prior){
    samples <- coda::mcmc(samples)
    model_fit <- structure(
      list(
        mcmc = coda::mcmc.list(samples),
        sample = nrow(samples),
        summary.pars = list(mutate = NULL),
        monitor = colnames(samples)
      ),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(model_fit, "prior_list") <- list(mu_f = prior)
    model_fit <- attach_test_parameter_map(model_fit)
    list(
      fit = model_fit,
      marglik = bridgesampling_object(0),
      prior_weights = 1
    )
  }

  posterior_1 <- matrix(
    c(1, 10, 2, 20, 3, 30),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )
  posterior_2 <- matrix(
    c(4, 40, 5, 50, 6, 60),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )
  ordered_prior_2 <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )$prior_list$mu_f

  mixed <- mix_posteriors(
    list(
      make_model(cbind(posterior_1, mu_f_ordered_total = c(11, 22, 33)), ordered_prior),
      make_model(cbind(posterior_2, mu_f_ordered_total = c(44, 55, 66)), ordered_prior_2)
    ),
    parameters = "mu_f",
    is_null_list = list(mu_f = c(FALSE, FALSE)),
    seed = 3,
    n_samples = 6
  )

  source_samples <- list(cbind(c(11,22,33)*.25,c(11,22,33)*.75),
    cbind(c(44,55,66)*.4,c(44,55,66)*.6))
  for(row_i in seq_len(nrow(mixed$mu_f))){
    model_i <- .bt_meta_get(mixed$mu_f, "component")[[row_i]]
    sample_i <- .bt_meta_get(mixed$mu_f, "draw_index")[[row_i]]
    expect_equal(
      unname(mixed$mu_f[row_i, ]),
      unname(source_samples[[model_i]][sample_i, ])
    )
  }
  expect_true(isTRUE(attr(mixed$mu_f, "ordered")))
  expect_false(is.null(attr(mixed$mu_f, "ordered_metadata")))
  # The first cumulative coordinate is level "mid"; the second is an
  # increment, contrast coefficient 2, never the position "[2]".
  expect_identical(colnames(mixed$mu_f), c("mu_f[mid]", "mu_f{2}"))
})

test_that("marginal posterior uses the stored full-rank ordered design", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ 0 + f,
    "mu",
    data = df,
    prior_list = list(
      f = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = c(.2, .3, .5),
        contrast = "cumulative_levels"
      )
    )
  )
  posterior <- matrix(
    c(2, 3, 5, 4, 6, 10),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, paste0("mu_f[", 1:3, "]"))
  )
  fit <- coda::mcmc(cbind(posterior, mu_f_ordered_total = c(10, 20)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  fit <- attach_test_parameter_map(fit)
  mixed <- as_mixed_posteriors(fit, parameters = "mu_f")

  formula_marginal <- marginal_posterior(
    mixed,
    parameter = "mu_f",
    formula = ~ 0 + f,
    prior_samples = FALSE
  )
  direct_marginal <- marginal_posterior(
    mixed,
    parameter = "mu_f",
    use_formula = FALSE,
    prior_samples = FALSE
  )

  expected <- list(
    low = c(2, 4),
    mid = c(5, 10),
    high = c(10, 20)
  )
  expect_equal(names(formula_marginal), names(expected))
  expect_equal(names(direct_marginal), names(expected))
  for(level in names(expected)){
    expect_equal(as.numeric(formula_marginal[[level]]), expected[[level]])
    expect_equal(as.numeric(direct_marginal[[level]]), expected[[level]])
  }
})

test_that("ordered levels combined with an intercept convolve on the common grid", {

  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  set.seed(8801)
  posterior <- cbind(
    mu_intercept = stats::rnorm(200, .2, .1),
    "mu_f[1]"    = stats::rnorm(200, .1, .05),
    "mu_f[2]"    = stats::rnorm(200, .15, .05)
  )
  level_densities <- function(allocation){
    formula_info <- JAGS_formula(
      y ~ f, "mu", data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior("normal", list(0, 1)), allocation = allocation)
      )
    )
    primitive_draws <- cbind(posterior, mu_f_ordered_total = seq_len(200)/100)
    if(is.null(allocation)) primitive_draws <- cbind(primitive_draws,
      "prior_par_eta_mu_f_ordered_alloc_f_1[1]" = 1,
      "prior_par_eta_mu_f_ordered_alloc_f_1[2]" = 1)
    fit <- coda::mcmc(primitive_draws)
    class(fit) <- c("mcmc", "BayesTools_fit")
    attr(fit, "prior_list") <- formula_info$prior_list
    fit <- attach_test_parameter_map(fit)
    mixed <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_f"))
    marginal <- marginal_posterior(mixed, parameter = "mu_f", formula = ~ f,
                                   prior_samples = TRUE, n_samples = 200)
    bf <- unlist(Savage_Dickey_BF(marginal, null_hypothesis = 0, silent = TRUE))
    expect_true(all(is.finite(bf) & bf > 0))
    lapply(marginal, .bt_meta_get, field = "prior_density")
  }

  # Each level is intercept + share * total. Heights are compared with the
  # analytic normal (fixed shares) or the one-dimensional mixture integral
  # over the Beta(1, 1) share; 1e-3 covers the omitted 1e-4 tails of both
  # sources and the product-grid quadrature (observed <= 4.2e-4).
  fixed <- level_densities(c(.4, .6))
  dirichlet <- level_densities(NULL)
  share_mixture <- function(value){
    stats::integrate(function(c) stats::dnorm(value, 0, sqrt(1 + c^2)),
                     0, 1, rel.tol = 1e-12)$value
  }
  for(value in c(0, 1)){
    expect_equal(.prior_linear_density_grid_height(fixed$mid, value),
                 stats::dnorm(value, 0, sqrt(1.16)), tolerance = 1e-3)
    expect_equal(.prior_linear_density_grid_height(fixed$high, value),
                 stats::dnorm(value, 0, sqrt(2)), tolerance = 1e-3)
    expect_equal(.prior_linear_density_grid_height(dirichlet$mid, value),
                 share_mixture(value), tolerance = 1e-3)
    expect_equal(.prior_linear_density_grid_height(dirichlet$high, value),
                 stats::dnorm(value, 0, sqrt(2)), tolerance = 1e-3)
    expect_equal(as.numeric(.prior_linear_density_height(dirichlet$high, value)),
                 stats::dnorm(value, 0, sqrt(2)), tolerance = 1e-3)
  }
})

test_that("ordered densities are direct when supported and bridge sampling stops for complex totals", {
  p <- prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
  attr(p, "levels") <- 3
  density_fixed <- density(p, n_points = 21)
  expect_s3_class(density_fixed, "density.prior.ordered")
  expect_equal(attr(density_fixed, "method"), "direct")
  expect_equal(density_fixed[[2]]$y[11], stats::dnorm(density_fixed[[2]]$x[11] / .25) / .25)

  p_dir <- prior_ordered(prior("normal", list(0, 1)), allocation = prior("dirichlet", list(alpha = c(1, 1))))
  attr(p_dir, "levels") <- 3
  density_dir <- density(p_dir, n_points = 11)
  expect_equal(attr(density_dir, "method"), "direct")

  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    )
  )
  samples <- matrix(
    1,
    nrow = 2,
    ncol = 6,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_f[1]",
      "mu_f[2]",
      "mu_f_ordered_total",
      "prior_par_eta_mu_f_ordered_alloc_f_1[1]",
      "prior_par_eta_mu_f_ordered_alloc_f_1[2]"
    ))
  )
  bridge_samples <- JAGS_bridgesampling_posterior(samples, formula_info$prior_list)
  expect_true("mu_f_ordered_total" %in% colnames(bridge_samples))
  expect_false("mu_f[1]" %in% colnames(bridge_samples))

  formula_info$prior_list$mu_f$total <- prior_spike_and_slab(prior("normal", list(0, 1)))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "only available when 'total' is a simple scalar prior"
  )

  formula_info$prior_list$mu_f$total <- prior("normal", list(0, expression(sigma)))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "does not support parameter expressions in 'total'"
  )

  formula_info$prior_list$mu_f$total <- prior("bernoulli", list(.5))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "requires a continuous or point-valued 'total' prior"
  )

  formula_info$prior_list$mu_f$total <- prior("point", list(1))
  expect_no_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list)
  )
})

test_that("ordered totals with boundary-infinite densities keep direct densities", {

  p <- prior_ordered(prior("gamma", list(shape = .5, rate = 1)))
  attr(p, "levels") <- 3
  densities <- density(p, n_points = 101)
  expect_equal(attr(densities, "method"), "direct")

  # Level 2 is total * c with c ~ Beta(1, 1): its density diverges at zero
  # with the total's density, and interior values are the product integral.
  middle <- densities[[2]]
  expect_identical(middle$y[middle$x == 0], Inf)
  interior <- which(middle$x > 0)[c(1, 25, 60)]
  reference <- vapply(middle$x[interior], function(value){
    stats::integrate(
      function(c) stats::dgamma(value / c, .5, 1) / c,
      lower = 0, upper = 1, rel.tol = 1e-10
    )$value
  }, numeric(1))
  expect_equal(middle$y[interior], reference, tolerance = 1e-6)

  highest <- densities[[3]]
  positive <- highest$x > 0
  expect_equal(highest$y[positive], stats::dgamma(highest$x[positive], .5, 1))

  bounded <- prior_ordered(prior("beta", list(.5, 1)))
  attr(bounded, "levels") <- 3
  expect_equal(attr(density(bounded, n_points = 51), "method"), "direct")

  mixed <- prior_ordered(prior_spike_and_slab(
    prior("gamma", list(shape = .5, rate = 1)),
    prior_inclusion = prior("spike", list(.5))
  ))
  attr(mixed, "levels") <- 3
  mixed_densities <- density(mixed, n_points = 101)
  expect_equal(attr(mixed_densities, "method"), "analytic_mixed_measure")
  expect_equal(mixed_densities[[3]]$atoms, data.frame(location = 0, mass = .5))
  expect_equal(attr(mixed_densities[[3]]$continuous, "mass"), .5)
})

test_that("ordered levels and products with a structural route are constructed from it", {

  # The first level of an ordered Cauchy total with a Dirichlet(2, 3)
  # allocation is C * S, S ~ Beta(2, 3), a scale mixture. Its capped product
  # grid (the total's 1e-4 tail range is +-3183) had no finite positive mass,
  # so the density could not be constructed; it is now built from its route.
  # References: f(0) = dcauchy(0) E[1 / S] = 4 dcauchy(0) (closed form), and
  # f(x) = int_0^1 dcauchy(x / s) dbeta(s, 2, 3) / s ds by integrate() at
  # rel.tol 1e-12. The same product as a 'multiply_by' term is the same
  # measure.
  ordered <- prior_ordered(prior("cauchy", list(0, 1)),
                           allocation = prior("dirichlet", list(alpha = c(2, 3))))
  attr(ordered, "levels") <- 3L
  ordered <- .prior_ordered_default_bound(ordered, "g")
  level <- .prior_linear_combination_density(list(g = ordered), c("g[1]" = 1))
  coefficient <- prior("cauchy", list(0, 1))
  attr(coefficient, "multiply_by") <- "s"
  product <- .prior_linear_combination_density(
    list(b = coefficient, s = prior("beta", list(2, 3))), c(b = 1)
  )
  reference <- function(value){
    stats::integrate(function(s) stats::dcauchy(value / s) * stats::dbeta(s, 2, 3) / s,
                     0, 1, rel.tol = 1e-12)$value
  }
  for(density in list(level, product)){
    expect_equal(density$density$mass, 1)
    at_zero <- prior_density_ordinate(density, 0)
    expect_identical(at_zero$behavior, "regular")
    expect_true(at_zero$exact)
    expect_identical(at_zero$method, "scale_mixture")
    expect_equal(exp(at_zero$log_density), 4 * stats::dcauchy(0), tolerance = 1e-12)
    for(value in c(-2, .4, 3)){
      expect_equal(as.numeric(.prior_linear_density_height(density, value)), reference(value),
                   tolerance = 1e-8)
    }
    curve <- .prior_linear_density_to_plot_data(density, x_range = c(-3, 3))$density
    expect_equal(curve$y[c(1, 75, 150)], vapply(curve$x[c(1, 75, 150)], reference, numeric(1)),
                 tolerance = 1e-8)
  }
})

test_that("mixed-measure ordered levels are evaluated on their structural route", {

  # A spike-and-slab Cauchy total (inclusion 1/2) with a Dirichlet(2, 3)
  # allocation: the first increment level is 1/2 point(0) + 1/2 C * S,
  # S ~ Beta(2, 3), and the last level 1/2 point(0) + 1/2 C. The capped level
  # grid cannot resolve the Cauchy scale (the plotted height at zero was
  # 0.0157 instead of 0.637); the route gives the exact heights. References:
  # integrate() of the scale mixture at rel.tol 1e-12 and dcauchy().
  total <- prior_spike_and_slab(prior("cauchy", list(0, 1)),
                                prior_inclusion = prior("spike", list(.5)))
  ordered <- prior_ordered(total, allocation = prior("dirichlet", list(alpha = c(2, 3))))
  attr(ordered, "levels") <- 3
  scale_mixture <- function(value){
    .5 * stats::integrate(function(s) stats::dcauchy(value / s) * stats::dbeta(s, 2, 3) / s,
                          0, 1, rel.tol = 1e-12)$value
  }

  densities <- density(ordered, n_points = 201, x_range = c(-3, 3))
  expect_equal(attr(densities, "method"), "analytic_mixed_measure")
  shared <- densities[[2]]
  expect_equal(shared$atoms, data.frame(location = 0, mass = .5))
  expect_equal(attr(shared$continuous, "mass"), .5)
  check <- c(1, 50, 101, 150, 201)
  expect_equal(
    shared$continuous$density[check],
    vapply(shared$continuous$x[check], scale_mixture, numeric(1)),
    tolerance = 1e-8
  )
  expect_identical(shared$diagnostics$grid$route, "mixture")
  full <- densities[[3]]
  expect_equal(full$continuous$density, .5 * stats::dcauchy(full$continuous$x), tolerance = 1e-12)

  # a transformation maps the display values and the heights (Jacobian)
  transformed <- density(ordered, n_points = 201, x_range = c(-3, 3), transformation = "exp")
  level <- transformed[[2]]
  expect_equal(level$atoms, data.frame(location = 1, mass = .5))
  values <- level$continuous$x[c(20, 101, 180)]
  expect_equal(
    level$continuous$density[c(20, 101, 180)],
    vapply(log(values), scale_mixture, numeric(1)) / values,
    tolerance = 1e-8
  )
})

test_that("prior curves without a structural route omit unresolved heavy-tailed product grids", {

  # A sum without a structural route is plotted from its numerical grid. The
  # grid of a product component with a heavy-tailed factor (a Cauchy ordered
  # total or 'multiply_by' coefficient, whose 1e-4 tail range spans about
  # 6000 scales on at most 1024 values) cannot resolve the product's scale:
  # the Riemann sum of the product density on it is about 80% below its mass,
  # and the plotted sums were 129-456% off. Such curves are omitted with a
  # classed warning; heights stay refused, and curves with a structural
  # route or a resolved product grid are drawn.
  bound <- function(total, alpha){
    p <- prior_ordered(total, allocation = prior("dirichlet", list(alpha = alpha)))
    attr(p, "levels") <- length(alpha) + 1L
    .prior_ordered_default_bound(p, "g")
  }
  plot_data <- function(density){
    warnings <- list()
    value <- withCallingHandlers(
      .prior_linear_density_to_plot_data(density, x_range = c(-4, 4), n_points = 101),
      warning = function(w){
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, warnings = warnings)
  }
  coefficient <- prior("cauchy", list(0, 1))
  attr(coefficient, "multiply_by") <- "s"
  unresolved <- list(
    .prior_linear_combination_density(
      list(a = prior("normal", list(0, 1)), g = bound(prior("cauchy", list(0, 1)), c(2, 3))),
      c(a = 1, "g[1]" = 1)
    ),
    .prior_linear_combination_density(
      list(a = prior("t", list(0, 1, 3)), b = coefficient, s = prior("beta", list(2, 3))),
      c(a = 1, b = 1)
    )
  )
  for(density in unresolved){
    expect_false(attr(density, "product_grid_resolution")$resolved)
    data <- plot_data(density)
    expect_null(data$value$density)
    expect_length(data$warnings, 1L)
    expect_s3_class(data$warnings[[1L]], "BayesTools_prior_curve_unavailable")
    expect_s3_class(data$warnings[[1L]], "BayesTools_plot_condition")
    expect_error(.prior_linear_density_height(density, 0), "numerical product grids are not used")
  }

  # a resolved product grid (t3 total) and a Cauchy level with its own route
  resolved <- .prior_linear_combination_density(
    list(a = prior("t", list(0, 1, 3)), g = bound(prior("t", list(0, 1, 3)), c(2, 3))),
    c(a = 1, "g[1]" = 1)
  )
  expect_true(attr(resolved, "product_grid_resolution")$resolved)
  level <- .prior_linear_combination_density(list(g = bound(prior("cauchy", list(0, 1)), c(2, 3))), c("g[1]" = 1))
  expect_false(attr(level, "product_grid_resolution")$resolved)
  for(density in list(resolved, level)){
    data <- plot_data(density)
    expect_length(data$warnings, 0L)
    expect_false(is.null(data$value$density))
  }

  # plot_marginal of y ~ f with a Cauchy ordered total: the level with the
  # unresolved product omits its prior curve, and the plot is drawn
  df <- data.frame(y = seq_len(6), f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high")))
  set.seed(8801)
  posterior <- cbind(mu_intercept = stats::rnorm(500, .2, .1),
                     "mu_f[1]" = stats::rnorm(500, .1, .05), "mu_f[2]" = stats::rnorm(500, .15, .05))
  formula_info <- JAGS_formula(y ~ f, "mu", data = df, prior_list = list(
    intercept = prior("normal", list(0, 1)), f = prior_ordered(prior("cauchy", list(0, 1)))))
  fit <- coda::mcmc(cbind(posterior, mu_f_ordered_total = seq_len(500)/100,
    "prior_par_eta_mu_f_ordered_alloc_f_1[1]" = 1,
    "prior_par_eta_mu_f_ordered_alloc_f_1[2]" = 1))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  fit <- attach_test_parameter_map(fit)
  mixed <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_f"))
  marginal <- marginal_posterior(mixed, parameter = "mu_f", formula = ~ f, prior_samples = TRUE, n_samples = 500)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(plot_marginal(list(mu_f = marginal), "mu_f", prior = TRUE),
                 class = "BayesTools_prior_curve_unavailable")
  expect_warning(plot <- plot_marginal(list(mu_f = marginal), "mu_f", prior = TRUE, plot_type = "ggplot"),
                 class = "BayesTools_prior_curve_unavailable")
  expect_s3_class(plot, "ggplot")
})

test_that("sampled ordered level densities reflect at the exact level support", {

  gamma_total <- prior_ordered(prior("gamma", list(2, 2)))
  attr(gamma_total, "levels") <- 3
  expect_equal(
    BayesTools:::.density.prior.ordered_level_bounds(gamma_total, 3L),
    list(c(0, 0), c(0, Inf), c(0, Inf))
  )
  fixed_normal <- prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
  attr(fixed_normal, "levels") <- 3
  expect_equal(
    BayesTools:::.density.prior.ordered_level_bounds(fixed_normal, 3L),
    list(c(0, 0), c(-Inf, Inf), c(-Inf, Inf))
  )
  bounded <- prior_ordered(prior("uniform", list(-1, 2)), contrast = "cumulative_levels")
  attr(bounded, "levels") <- 3
  expect_equal(
    BayesTools:::.density.prior.ordered_level_bounds(bounded, 3L),
    rep(list(c(-1, 2)), 3)
  )

  # Level 2 of exp(total * share) near its support bound exp(0) = 1. Before
  # reflection the KDE lost 38% of this bin; the remaining boundary bias of
  # the reflected KDE is about 4% here (the level density has a nonzero slope
  # at 0), and the Monte Carlo SE of the reference bin is about 0.7%.
  set.seed(4104)
  transformed <- density(gamma_total, transformation = "exp", n_samples = 1e5,
                         n_points = 2001)
  set.seed(4105)
  reference <- exp(rng(gamma_total, 2e5, transform_factor_samples = TRUE))
  level <- transformed[[2]]
  expect_true(isTRUE(attr(level, "boundary_reflection")))
  keep <- level$x >= 1 & level$x <= 1.05
  bin <- sum(diff(level$x[keep]) * (utils::head(level$y[keep], -1) + utils::tail(level$y[keep], -1)) / 2)
  expect_equal(bin, mean(reference[, 2] > 1 & reference[, 2] <= 1.05), tolerance = .1)
})

test_that("ordered ranges cover zero and scaled total densities", {
  p_positive <- prior_ordered(
    prior("normal", list(10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p_positive, "levels") <- 3

  expected_range <- c(0, stats::qnorm(.995, mean = 10, sd = 1))
  expect_equal(range(p_positive), expected_range)
  expect_equal(
    range(p_positive, quantiles = .1),
    c(0, stats::qnorm(.9, mean = 10, sd = 1))
  )

  density_positive <- density(p_positive, n_points = 101)
  grid_step <- diff(expected_range) / 100
  expect_equal(attr(density_positive, "x_range"), expected_range)
  expect_lte(
    abs(density_positive[[2]]$x[which.max(density_positive[[2]]$y)] - 2.5),
    grid_step
  )
  expect_gt(max(density_positive[[2]]$y), 1)

  p_negative <- prior_ordered(
    prior("normal", list(-10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  expect_equal(
    range(p_negative),
    c(stats::qnorm(.005, mean = -10, sd = 1), 0)
  )
})

test_that("ordered continuous-mixture totals have deterministic density ranges", {
  total <- prior_mixture(
    list(
      prior("normal", list(-20, 1)),
      prior("normal", list(10, 1))
    ),
    components = c("lower", "upper")
  )
  p <- prior_ordered(
    total,
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  expected_range <- c(
    stats::qnorm(.005, mean = -20, sd = 1),
    stats::qnorm(.995, mean = 10, sd = 1)
  )
  expect_equal(range(p), expected_range)

  set.seed(7403)
  density_mixture <- density(p, n_points = 101, n_samples = 2000)

  expect_s3_class(density_mixture, "density.prior.ordered")
  expect_equal(attr(density_mixture, "x_range"), expected_range)
  expect_s3_class(density_mixture[[1]], "density.prior.point")
  expect_true(all(vapply(density_mixture, function(component){
    all(is.finite(component$y))
  }, logical(1))))
  expect_equal(
    density_mixture[[2]]$samples * 4,
    density_mixture[[3]]$samples
  )
  expect_true(all(vapply(density_mixture[2:3], function(component){
    max(component$y) > 0
  }, logical(1))))
})

test_that("ordered sampled densities transform components and preserve cumulative atoms", {
  p <- prior_ordered(
    prior("normal", list(0, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  set.seed(7401)
  density_identity <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE
  )
  set.seed(7401)
  density_exp <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE,
    transformation = "exp"
  )

  expect_s3_class(density_exp[[1]], "density.prior.point")
  expect_equal(density_exp[[1]]$x, 1)
  expect_equal(unique(density_exp[[1]]$samples), 1)
  expect_equal(density_exp[[1]]$y, 1)
  expect_null(density_exp[[1]]$bw)

  expect_s3_class(density_exp[[2]], "density.prior.simple")
  expect_equal(density_exp[[2]]$x, exp(density_identity[[2]]$x))
  expect_equal(density_exp[[2]]$samples, exp(density_identity[[2]]$samples))
  expect_equal(
    density_exp[[2]]$y,
    density_identity[[2]]$y / density_exp[[2]]$x,
    tolerance = 1e-10
  )
  expect_true(all(density_exp[[2]]$x > 0))

  plot_exp <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE,
    transformation = "exp"
  )
  expect_s3_class(plot_exp, "ggplot")
  expect_s3_class(plot_exp$layers[[1]]$geom, "GeomSegment")
  expect_equal(NROW(plot_exp$layers[[1]]$data), 1)
  expect_equal(plot_exp$layers[[1]]$data$x, 1)
  expect_equal(plot_exp$layers[[1]]$data$xend, 1)
  expect_equal(plot_exp$layers[[1]]$data$yend, 1)
})

test_that("ordered cumulative atoms contribute to direct and transformed plot ranges", {
  p <- prior_ordered(
    prior("normal", list(10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  density_direct <- density(p, n_points = 31)
  expect_s3_class(density_direct[[1]], "density.prior.point")
  expect_equal(density_direct[[1]]$x, 0)
  expect_equal(min(attr(density_direct[[1]], "x_range")), 0)
  expect_equal(min(attr(density_direct, "x_range")), 0)

  plot_direct <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31
  )
  expect_s3_class(plot_direct$layers[[1]]$geom, "GeomSegment")
  expect_equal(plot_direct$layers[[1]]$data$x, 0)
  direct_limits <- plot_direct$scales$get_scales("x")$limits
  expect_lte(min(direct_limits), 0)
  expect_gte(max(direct_limits), 0)

  set.seed(7402)
  density_exp <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    transformation = "exp"
  )
  expect_s3_class(density_exp[[1]], "density.prior.point")
  expect_equal(density_exp[[1]]$x, 1)
  expect_equal(min(attr(density_exp[[1]], "x_range")), 1)
  expect_equal(min(attr(density_exp, "x_range")), 1)

  set.seed(7402)
  plot_exp <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31,
    n_samples = 2000,
    transformation = "exp"
  )
  expect_s3_class(plot_exp$layers[[1]]$geom, "GeomSegment")
  expect_equal(plot_exp$layers[[1]]$data$x, 1)
  exp_limits <- plot_exp$scales$get_scales("x")$limits
  expect_lte(min(exp_limits), 1)
  expect_gte(max(exp_limits), 1)
})

test_that("ordered mixed-measure densities preserve atoms and continuous mass", {
  total <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("spike", list(.5))
  )
  fixed <- prior_ordered(
    total,
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(fixed, "levels") <- 3

  density_fixed <- density(fixed, n_points = 201)
  expect_equal(attr(density_fixed, "method"), "analytic_mixed_measure")
  expect_equal(attr(density_fixed, "measure_schema_version"), 1L)
  expect_true(all(vapply(
    density_fixed,
    inherits,
    logical(1),
    what = "density.prior.mixed_measure"
  )))

  expect_equal(density_fixed[[1]]$atoms, data.frame(location = 0, mass = 1))
  expect_null(density_fixed[[1]]$continuous)
  for(i in 2:3){
    expect_equal(
      density_fixed[[i]]$atoms,
      data.frame(location = 0, mass = .5)
    )
    expect_equal(attr(density_fixed[[i]]$continuous, "mass"), .5)
    curve <- density_fixed[[i]]$continuous
    integral <- BayesTools:::.density.prior.ordered_curve_integral(
      curve$x, curve$density
    )
    scale <- c(.25, 1)[[i - 1L]]
    expected <- .5 * diff(stats::pnorm(range(curve$x) / scale))
    # The source grid drops 1e-4 in each normal tail; allow that
    # normalization error and trapezoid discretization, not display clipping.
    expect_lt(abs(integral - expected), 2e-4)
    expect_equal(density_fixed[[i]]$diagnostics$continuous_integral, integral)
    expect_equal(density_fixed[[i]]$diagnostics$atom_mass, .5)
    expect_equal(density_fixed[[i]]$diagnostics$continuous_mass, .5)
    expect_equal(
      density_fixed[[i]]$diagnostics$method,
      "analytic_components"
    )
  }

  density_exp <- density(
    fixed,
    n_points = 201,
    transformation = "exp"
  )
  expect_equal(density_exp[[2]]$atoms, data.frame(location = 1, mass = .5))
  expect_equal(attr(density_exp[[2]]$continuous, "mass"), .5)
  expect_equal(density_exp[[2]]$transformation$name, "exp")
  expect_true(all(density_exp[[2]]$continuous$x > 0))

  random_allocation <- prior_ordered(
    total,
    allocation = prior("dirichlet", list(alpha = c(2, 3))),
    contrast = "cumulative"
  )
  attr(random_allocation, "levels") <- 3
  density_random <- density(random_allocation, n_points = 201)
  expect_equal(density_random[[2]]$atoms$mass, .5)
  expect_equal(density_random[[3]]$atoms$mass, .5)
  expect_equal(attr(density_random[[2]]$continuous, "mass"), .5)
  expect_equal(attr(density_random[[3]]$continuous, "mass"), .5)

  plot_middle <- plot(
    fixed,
    plot_type = "ggplot",
    show_figures = 2,
    n_points = 201
  )
  geom_classes <- vapply(
    plot_middle$layers,
    function(layer) class(layer$geom)[1L],
    character(1)
  )
  expect_true("GeomLine" %in% geom_classes)
  expect_true("GeomSegment" %in% geom_classes)
  point_layer <- plot_middle$layers[[which(geom_classes == "GeomSegment")[1L]]]
  expect_equal(point_layer$data$x, 0)
  expect_equal(point_layer$data$yend, .5)

  layer_geoms <- geom_prior(fixed, show_parameter = 2, n_points = 201)
  expect_s3_class(layer_geoms, "BayesTools_prior_overlay")
  expect_equal(
    vapply(layer_geoms$geoms, function(layer) class(layer$geom)[1L], character(1)),
    c("GeomLine", "GeomSegment")
  )
})

test_that("clipping ordered mixed densities preserves their heights", {

  dist <- list(
    density = list(x = seq(-3, 3, length.out = 601),
                   y = stats::dnorm(seq(-3, 3, length.out = 601)), mass = .5),
    points = data.frame(x = 0, p = .5),
    n_grid = 601L
  )
  grid <- seq(-.5, .5, length.out = 101)
  clipped <- BayesTools:::.density.prior.ordered_regrid_mixed(dist, grid)
  expect_equal(clipped$density$y, stats::dnorm(grid), tolerance = 1e-14)
  expect_identical(clipped$density$mass, .5)
  expect_identical(clipped$points, dist$points)
  expect_lt(abs(
    attr(clipped, "ordered_grid_diagnostics")$captured_continuous_shape_integral -
    diff(stats::pnorm(c(-.5, .5)))
  ), 4e-6)

  total <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("spike", list(.5))
  )
  ordered <- prior_ordered(total, allocation = c(.25, .75),
                           contrast = "cumulative")
  attr(ordered, "levels") <- 3
  result <- density(ordered, x_range = c(-.5, .5), n_points = 1001)
  curve <- result[[3]]$continuous
  expect_lt(max(abs(curve$density - .5 * stats::dnorm(curve$x))), 1e-4)
  expect_lt(result[[3]]$diagnostics$continuous_integral, .2)
  expect_equal(result[[3]]$atoms$mass, .5)
  expect_equal(result[[3]]$diagnostics$continuous_mass, .5)
  transformed <- density(ordered, x_range = c(-.5, .5), n_points = 1001,
                         transformation = "exp")
  curve <- transformed[[3]]$continuous
  expect_lt(max(abs(curve$density - .5 * stats::dlnorm(curve$x))), 1e-4)
  expect_lt(transformed[[3]]$diagnostics$continuous_integral, .2)
  expect_equal(transformed[[3]]$atoms, data.frame(location = 1, mass = .5))
})

test_that("ordered mixed measures propagate through marginal inference", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(
      rep(c("low", "mid", "high"), 2),
      levels = c("low", "mid", "high")
    )
  )
  formula_info <- JAGS_formula(
    y ~ 0 + f,
    "mu",
    data = df,
    prior_list = list(
      f = prior_ordered(
        prior_spike_and_slab(
          prior("normal", list(0, 1)),
          prior_inclusion = prior("spike", list(.5))
        ),
        allocation = c(.25, .75),
        contrast = "cumulative"
      )
    )
  )
  posterior <- cbind(
    "mu_f[1]" = c(0, .1, .2, 0),
    "mu_f[2]" = c(0, .2, .3, 0),
    "mu_f_ordered_total" = c(0, .3, .5, 0),
    "mu_f_ordered_total_indicator" = c(0, 1, 1, 0)
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(fit, parameters = "mu_f")
  coefficient_atoms <- BayesTools:::.posterior_atoms_get(samples$mu_f)
  # The first cumulative coordinate is level "mid"; the second is the
  # increment, contrast coefficient 2.
  expect_equal(coefficient_atoms$locations, matrix(
    0,
    nrow = 1,
    ncol = 2,
    dimnames = list(NULL, c("mu_f[mid]", "mu_f{2}"))
  ))
  expect_equal(coefficient_atoms$mass, .5)
  # the component of the spike-and-slab total: its slab (1) or spike (2)
  expect_equal(
    .bt_meta_get(samples$mu_f, "ordered_total_component"),
    c(2L, 1L, 1L, 2L)
  )

  for(use_formula in c(FALSE, TRUE)){
    marginal <- marginal_posterior(
      samples,
      parameter = "mu_f",
      formula = if(use_formula) ~ 0 + f else NULL,
      use_formula = use_formula,
      prior_samples = TRUE,
      n_samples = 128
    )
    expected_mass <- c(low = 1, mid = .5, high = .5)
    for(level in names(expected_mass)){
      prior_density <- .bt_meta_get(marginal[[level]], "prior_density")
      expect_equal(
        BayesTools:::.prior_linear_density_point_mass(prior_density, 0),
        expected_mass[[level]],
        tolerance = 1e-12
      )
      if(expected_mass[[level]] < 1){
        expect_equal(
          prior_density$density$mass,
          1 - expected_mass[[level]],
          tolerance = 1e-12
        )
      }
      posterior_atoms <- BayesTools:::.posterior_atoms_get(
        marginal[[level]]
      )
      expect_equal(posterior_atoms$mass, expected_mass[[level]])
      expect_equal(unname(posterior_atoms$locations[, 1L]), 0)
    }
  }

  expect_error(
    BayesTools:::.as_mixed_posteriors.factor(
      posterior[, c("mu_f[1]", "mu_f[2]"), drop = FALSE],
      formula_info$prior_list$mu_f,
      "mu_f"
    ),
    "required total-prior indicator"
  )
})

test_that("ordered Dirichlet marginals reject non-level combinations", {
  p <- prior_ordered(
    prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("spike", list(.5))
    ),
    allocation = prior("dirichlet", list(alpha = c(2, 3))),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3
  p <- BayesTools:::.prior_ordered_default_bound(p, "mu_f")

  expect_error(
    BayesTools:::.prior_ordered_linear_distribution(
      p,
      weights = c(1, 2),
      indices = 1:2,
      dx = .01,
      n_grid = 128
    ),
    "not a level or a single allocation subset"
  )
})

test_that("ordered levels are exact products of the total and its allocation share", {

  # A level (or allocation subset) is total * w with the allocation share
  # w ~ Beta(a, b) (.prior_ordered_linear_share()): f(x) = int f_T(x / w) /
  # w f_w(w) dw, and at 0 f(0) = f_T(0) E[1 / w] = f_T(0) (a + b - 1) /
  # (a - 1) for a > 1, infinite for a <= 1. References: integrate() at
  # rel.tol 1e-12 with stats:: densities; for the Beta(1, 1) share
  # f(x) = int_|x|^Inf phi(t) / t dt. The capped product grid gave -2.5% at 0.
  split_reference <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12,
                       subdivisions = 5000L)$value
    }, numeric(1)))
  }
  bound <- function(total, alpha, name){
    p <- prior_ordered(total, allocation = prior("dirichlet", list(alpha = alpha)),
                       contrast = "cumulative")
    attr(p, "levels") <- length(alpha) + 1L
    .prior_ordered_default_bound(p, name)
  }
  height <- function(density, value) as.numeric(.prior_linear_density_height(density, value))

  normal_total <- bound(prior("normal", list(0, 1)), c(2, 3), "f")
  for(level in c("f[1]", "f[2]")){
    share <- .prior_ordered_linear_share(normal_total, stats::setNames(1, level),
                                         match(level, c("f[1]", "f[2]")))
    expect_identical(share$type, "beta")
    a <- share$alpha
    density <- .prior_linear_combination_density(list(f = normal_total), stats::setNames(1, level))
    zero <- prior_density_ordinate(density, 0)
    expect_identical(zero$behavior, "regular")
    expect_true(zero$exact)
    expect_equal(height(density, 0), stats::dnorm(0) * (a[1] + a[2] - 1) / (a[1] - 1), tolerance = 1e-12)
    for(value in c(.1, 1)){
      expect_equal(
        height(density, value),
        split_reference(function(w) stats::dnorm(value / w) / w * stats::dbeta(w, a[1], a[2]),
                        c(0, 1e-3, .01, .1, .5, 1)),
        tolerance = 1e-8
      )
    }
  }

  # Beta(1, 1) share of the first level: an infinite ordinate at 0
  flat <- bound(prior("normal", list(0, 1)), c(1, 1), "g")
  density <- .prior_linear_combination_density(list(g = flat), c("g[1]" = 1))
  singular <- prior_density_ordinate(density, 0)
  expect_identical(singular$behavior, "infinite")
  expect_true(singular$exact)
  expect_equal(height(density, .1), 0.94271965, tolerance = 1e-8)
  expect_equal(height(density, .1),
               split_reference(function(t) stats::dnorm(t) / t, c(.1, 1, 10, Inf)),
               tolerance = 1e-8)

  # a non-normal two-sided total is a scale product; with a Beta(1, b) share
  # (positive density at 0) its level density at 0 is infinite from both
  # sides, f(x) = int_0^1 f_T(x / w) / w dw ~ -f_T(0) log|x| (no jump), and
  # with a Beta(2, 3) share it is f_T(0) (2 + 3 - 1) / (2 - 1)
  for(total in list(prior("t", list(0, 1, 3)), prior("cauchy", list(0, 1)))){
    f_total <- function(x) mpdf(total, x)
    for(alpha in list(c(1, 1), c(1, 4))){
      density <- .prior_linear_combination_density(list(t = bound(total, alpha, "t")), c("t[1]" = 1))
      singular <- prior_density_ordinate(density, 0)
      expect_identical(singular$behavior, "infinite")
      expect_true(singular$exact)
      expect_identical(singular$method, "scale_mixture")
      expect_equal(height(density, .4),
                   split_reference(function(w) f_total(.4 / w) / w * stats::dbeta(w, alpha[1], alpha[2]),
                                   c(0, 1e-3, .01, .1, .4, 1)),
                   tolerance = 1e-8)
    }
  }
  t_total <- bound(prior("t", list(0, 1, 3)), c(2, 3), "t")
  density <- .prior_linear_combination_density(list(t = t_total), c("t[1]" = 1))
  expect_equal(height(density, 0), stats::dt(0, 3) * (2 + 3 - 1) / (2 - 1), tolerance = 1e-12)
  # the Savage-Dickey point hypothesis at the infinite ordinate is refused
  flat_t <- .prior_linear_combination_density(
    list(t = bound(prior("t", list(0, 1, 3)), c(1, 1), "t")), c("t[1]" = 1)
  )
  draws <- .bt_meta_update(
    structure(stats::qnorm(seq(.001, .999, length.out = 500), .2, .3), class = c("marginal_posterior.simple", "marginal_posterior", "numeric")),
    prior_density = flat_t,
    atoms = posterior_atom_attribute()
  )
  expect_error(hypothesis_BF(draws, hypothesis = "theta = 0", parameter = "theta"),
               class = "BayesTools_infinite_ordinate")

  # a half-normal total is a scale product: zero below 0, and at 0
  # f_T(0) E[1 / w] with f_T(0) = 2 phi(0)
  positive_total <- bound(prior("normal", list(0, 1), list(0, Inf)), c(2, 3), "h")
  density <- .prior_linear_combination_density(list(h = positive_total), c("h[1]" = 1))
  expect_identical(prior_density_ordinate(density, -.1)$behavior, "zero")
  expect_identical(prior_density_ordinate(density, .1)$method, "scale_mixture")
  expect_equal(height(density, 0), 2 * stats::dnorm(0) * (2 + 3 - 1) / (2 - 1), tolerance = 1e-12)
  expect_equal(height(density, .1),
               split_reference(function(w) 2 * stats::dnorm(.1 / w) / w * stats::dbeta(w, 2, 3),
                               c(0, 1e-3, .01, .1, .5, 1)),
               tolerance = 1e-8)
  side <- hypothesis_parse("theta < 0.5")$statements[[1L]]$left
  expect_equal(.hypothesis_prior_density_prob(density, side, "theta"),
               split_reference(function(w) (2 * stats::pnorm(.5 / w) - 1) * stats::dbeta(w, 2, 3),
                               c(0, .1, .5, 1)),
               tolerance = 1e-8)

  # with an intercept the level is a conditional-normal scale mixture
  density <- .prior_linear_combination_density(
    list(mu = prior("normal", list(0, 1)), f = normal_total), c(mu = 1, "f[1]" = 1)
  )
  expect_identical(prior_density_ordinate(density, 0)$method, "conditional_normal_mixture")
  for(value in c(0, .5)){
    expect_equal(height(density, value),
                 split_reference(function(w) stats::dnorm(value, 0, sqrt(1 + w^2)) * stats::dbeta(w, 2, 3),
                                 c(0, .1, .5, 1)),
                 tolerance = 1e-8)
  }
})

test_that("estimates tables show ordered totals and shares or the level effects", {

  # A fit with draws of the sampled parameters of its ordered priors (totals
  # and gamma allocation nodes, each node its share times a positive per-draw
  # scale that the normalization removes) and the increments they give; the
  # other coefficients are standard normal draws.
  data <- data.frame(
    f = ordered(rep(c("lo", "mid", "hi"), 6), levels = c("lo", "mid", "hi")),
    o = ordered(rep(c("p", "q", "r"), each = 6), levels = c("p", "q", "r")),
    g = factor(rep(c("A", "B", "C"), each = 6)),
    x = sin(seq_len(18))
  )
  ordered_table_fit <- function(formula, priors, formula_scale = NULL){
    result <- JAGS_formula(formula, "mu", data, priors)
    prior_list <- result$prior_list
    set.seed(1)
    n <- 40L
    columns <- list()
    shares  <- list()
    for(name in names(prior_list)){
      prior <- prior_list[[name]]
      if(!is.prior.ordered(prior)){
        coefficient_names <- if(BayesTools:::.bt_prior_is_factor_family(prior)){
          BayesTools:::.JAGS_prior_factor_names(name, prior)
        }else{
          name
        }
        columns[[name]] <- matrix(stats::rnorm(n * length(coefficient_names)), nrow = n,
                                  dimnames = list(NULL, coefficient_names))
        next
      }
      draws <- BayesTools:::.prior_ordered_draws(prior, n)
      coefficients <- draws$coefficients
      colnames(coefficients) <- BayesTools:::.JAGS_prior_factor_names(name, prior)
      total <- draws$theta
      colnames(total) <- BayesTools:::.prior_ordered_total_monitor_names(prior, name)
      columns[[name]] <- cbind(coefficients, total)
      for(record in BayesTools:::.prior_ordered_dirichlet_records(prior)){
        eta_name <- BayesTools:::.JAGS_prior_dirichlet_eta_name(record$node)
        if(eta_name %in% names(shares)){
          next
        }
        share <- draws$allocation_samples[[record$key]]
        shares[[eta_name]] <- share
        columns[[eta_name]] <- matrix(
          share * stats::rexp(n), nrow = n,
          dimnames = list(NULL, paste0(eta_name, "[", seq_len(record$dim), "]"))
        )
      }
    }
    samples <- do.call(cbind, unname(columns))
    fit <- structure(
      list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = n),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(fit, "prior_list") <- prior_list
    attr(fit, "formula_design") <- list(mu = result$formula_design)
    if(!is.null(formula_scale)){
      attr(fit, "formula_scale") <- list(mu = formula_scale)
    }
    list(fit = attach_test_parameter_map(fit), samples = samples, shares = shares)
  }
  allocation_rows <- function(case){
    table <- JAGS_estimates_table(case$fit, transform_factors = FALSE,
                                  remove_diagnostics = TRUE)
    rows <- grep("_ordered_allocation", rownames(table), value = TRUE)
    list(rows = rows, means = unname(table[rows, "Mean"]))
  }
  share_means <- function(case, nodes){
    unname(unlist(lapply(nodes, function(node) colMeans(case$shares[[node]]))))
  }
  N01 <- prior("normal", list(0, 1))

  # an ordered factor in an interaction with an ordinary factor: one total
  # and one allocation per slice (level of the ordinary factor's contrast)
  slices <- ordered_table_fit(~ f * g, list(
    intercept = N01, f = prior_ordered(N01),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"),
    "f:g" = prior_ordered(N01)
  ))
  out <- allocation_rows(slices)
  expect_identical(out$rows, c(
    "(mu) f_ordered_allocation[mid]", "(mu) f_ordered_allocation[hi]",
    "(mu) f:g_ordered_allocation[1][mid]", "(mu) f:g_ordered_allocation[1][hi]",
    "(mu) f:g_ordered_allocation[2][mid]", "(mu) f:g_ordered_allocation[2][hi]"
  ))
  expect_equal(out$means, share_means(slices, c(
    "prior_par_eta_mu_f_ordered_alloc_f_1",
    "prior_par_eta_mu_f__xXx__g_ordered_alloc_f_1",
    "prior_par_eta_mu_f__xXx__g_ordered_alloc_f_2"
  )), tolerance = 1e-12)
  # the transformed table shows the level effects: no zero reference cells
  # and no totals (the last level of each slice equals its total)
  transformed <- JAGS_estimates_table(slices$fit, transform_factors = TRUE,
                                      remove_diagnostics = TRUE)
  expect_identical(rownames(transformed), c(
    "(mu) intercept", "(mu) f[mid]", "(mu) f[hi]", "(mu) g[B]", "(mu) g[C]",
    "(mu) f[mid]:g[B]", "(mu) f[hi]:g[B]", "(mu) f[mid]:g[C]", "(mu) f[hi]:g[C]"
  ))
  expect_equal(
    unname(transformed[c("(mu) f[hi]:g[B]", "(mu) f[hi]:g[C]"), "Mean"]),
    unname(colMeans(slices$samples[, c("mu_f__xXx__g_ordered_total[1]",
                                       "mu_f__xXx__g_ordered_total[2]")])),
    tolerance = 1e-12
  )

  # two ordered factors in one term: one allocation per factor
  two <- ordered_table_fit(~ f * o, list(
    intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
    "f:o" = prior_ordered(N01)
  ))
  out <- allocation_rows(two)
  expect_identical(out$rows, c(
    "(mu) f_ordered_allocation[mid]", "(mu) f_ordered_allocation[hi]",
    "(mu) o_ordered_allocation[q]", "(mu) o_ordered_allocation[r]",
    "(mu) f:o_ordered_allocation_f[mid]", "(mu) f:o_ordered_allocation_f[hi]",
    "(mu) f:o_ordered_allocation_o[q]", "(mu) f:o_ordered_allocation_o[r]"
  ))
  expect_equal(out$means, share_means(two, c(
    "prior_par_eta_mu_f_ordered_alloc_f_1",
    "prior_par_eta_mu_o_ordered_alloc_o_1",
    "prior_par_eta_mu_f__xXx__o_ordered_alloc_f_1",
    "prior_par_eta_mu_f__xXx__o_ordered_alloc_o_1"
  )), tolerance = 1e-12)

  # a shared allocation is shown with every term that uses it
  shared_priors <- list(
    intercept = N01, f = prior_ordered(N01, id = "shape"), x = N01,
    "f:x" = prior_ordered(N01, id = "shape")
  )
  shared <- ordered_table_fit(~ f + x + f:x, shared_priors)
  out <- allocation_rows(shared)
  expect_identical(out$rows, c(
    "(mu) f_ordered_allocation[mid]", "(mu) f_ordered_allocation[hi]",
    "(mu) f:x_ordered_allocation[mid]", "(mu) f:x_ordered_allocation[hi]"
  ))
  expect_equal(out$means, rep(share_means(shared, "prior_par_eta_ordered_alloc_shape_f"), 2L),
               tolerance = 1e-12)

  # the totals and shares are the fitted parameters: when the original-scale
  # transformation changes the ordered terms (a standardized 'x' in 'f:x'),
  # the table notes that they remain on the fitted scale
  fitted_scale_note <- paste0(
    "Ordered-factor totals and allocations are summarized on the fitted ",
    "(standardized) scale. Use 'transform_factors = TRUE' for the ",
    "original-scale level effects."
  )
  expect_null(attr(JAGS_estimates_table(shared$fit, transform_scaled = TRUE), "footnotes"))
  scaled <- ordered_table_fit(~ f + x + f:x, shared_priors, formula_scale_for_test(
    ~ f + x + f:x, list(x = list(mean = 2, sd = 3)), data = data, prior_list = shared_priors
  ))
  expect_identical(
    attr(JAGS_estimates_table(scaled$fit, transform_scaled = TRUE), "footnotes"),
    fitted_scale_note
  )
  expect_null(attr(JAGS_estimates_table(scaled$fit, transform_scaled = FALSE), "footnotes"))
  expect_null(attr(JAGS_estimates_table(scaled$fit, transform_scaled = TRUE,
                                        transform_factors = TRUE), "footnotes"))
})
