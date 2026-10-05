skip_if_not_test_profile("unit")

.image_point_draws <- function(location, with_prior = TRUE){
  x <- structure(rep(location, 4), class = c("marginal_posterior.simple", "marginal_posterior"))
  posterior_metadata(x, "support") <- posterior_support_attribute(c(location, location), points = location, type = "points")
  posterior_metadata(x, "atoms") <- posterior_atom_attribute(point_masses = data.frame(x = location, mass = 1))
  if(with_prior) posterior_metadata(x, "prior_density") <- prior("point", list(location))
  x
}

test_that("R133 finite named images refuse false boundary atoms across public and point producers", {
  cases <- list(list(-1000, "exp", NULL), list(20, "tanh", NULL), list(-20, "tanh", NULL),
                list(1e-300, "exp_lin", list(b = 2)))
  for(case in cases){
    location <- case[[1L]]
    transformation <- case[[2L]]
    arguments <- case[[3L]]
    point <- .prior_linear_density_point(location)
    draws <- .image_point_draws(location)
    routes <- list(
      function() posterior_transform(draws, transformation, arguments),
      function() posterior_transform(.image_point_draws(location, FALSE), transformation, arguments),
      function() .prior_linear_density_transform(point, transformation, arguments),
      function() density(prior("point", list(location)), transformation = transformation, transformation_arguments = arguments),
      function() .prior_linear_density_to_plot_data(point, transformation = transformation, transformation_arguments = arguments),
      function() .plot_data_samples.simple(list(theta = draws), "theta", 64, transformation, arguments, FALSE),
      function() .plot_data_factor_points(location, .3, 64, transformation, arguments),
      function() .plot_data_marginal_samples.den(draws, 64, transformation, arguments, FALSE),
      function() .posterior_support_transform(posterior_metadata(draws, "support"), transformation, arguments))
    conditions <- lapply(routes, function(route) tryCatch(route(), error = identity))
    expect_true(all(vapply(conditions, inherits, logical(1), "BayesTools_transformation_image_unavailable")),
      info = paste(vapply(conditions, function(x) paste(class(x), if(inherits(x, "error")) conditionMessage(x) else "returned", collapse = "/"), character(1)), collapse = "\n"))
    for(condition in conditions[vapply(conditions, inherits, logical(1), "BayesTools_transformation_image_unavailable")]){
      expect_s3_class(condition, "BayesTools_transformation")
      expect_null(conditionCall(condition))
      expect_identical(condition$transformation, transformation)
      expect_true(all(condition$source_values == location))
      expected_image <- switch(transformation, exp = 0, tanh = sign(location), exp_lin = 0)
      expect_true(all(condition$images == expected_image))
      expect_identical(conditionMessage(condition), paste0("The '", transformation,
        "' transformation is numerically unavailable: finite source values do not have representable images under this transformation. Use a transformation whose images remain representable."))
    }
  }
})

test_that("R133 ordinary point images and true source zero preserve atoms and masses", {
  for(case in list(list(-1, "exp", NULL, exp(-1)), list(.5, "tanh", NULL, tanh(.5)),
                  list(0, "exp_lin", list(b = 2), 0))){
    x <- posterior_transform(.image_point_draws(case[[1L]]), case[[2L]], case[[3L]])
    expect_identical(as.numeric(x), rep(case[[4L]], 4))
    expect_identical(posterior_metadata(x, "atoms")$mass, 1)
    expect_identical(as.numeric(posterior_metadata(x, "atoms")$locations), case[[4L]])
    expect_identical(prior_density_ordinate(posterior_metadata(x, "prior_density"), case[[4L]])$behavior, "point_mass")
    expect_identical(prior_density_ordinate(posterior_metadata(x, "prior_density"), case[[4L]])$point_mass, 1)
  }
  point <- density(prior("point", list(0)), transformation = "exp_lin", transformation_arguments = list(a = 1, b = 0))
  expect_identical(point$x, exp(1))
  expect_identical(point$y, 1)
  expect_error(posterior_transform(.image_point_draws(0), "exp_lin", list(b = 0)), class = "BayesTools_nonmonotone_transformation")
  expect_error(posterior_transform(.image_point_draws(-1), "exp_lin"), class = "BayesTools_transformation_domain")
  expect_error(posterior_transform(.image_point_draws(1000), "exp"), class = "BayesTools_transformation_domain")
  custom <- list(fun = function(x) x, inv = function(x) x, jac = function(x) rep(0, length(x)))
  expect_error(posterior_transform(.image_point_draws(.5), custom), class = "BayesTools_nonmonotone_transformation")
  expect_error(.prior_linear_density_transform(.prior_linear_density_point(0), "exp_lin", list(a = -1000, b = 0)),
    "The constant prior-density transformation is not representable.", fixed = TRUE)
})

test_that("R133 matrix atoms scalar marginals and infinite support limits retain their meaning", {
  margin <- posterior_atom_attribute(point_masses = data.frame(x = .5, mass = .3))
  atoms <- .posterior_atoms_new(locations = matrix(c(.5, -1), 1, dimnames = list(NULL, c("a", "b"))),
    mass = .3, column_names = c("a", "b"),
    marginals = list(a = margin, b = posterior_atom_attribute(point_masses = data.frame(x = -1, mass = .3))))
  mapped <- .posterior_atoms_transform(atoms, "exp")
  expect_identical(dim(mapped$locations), c(1L, 2L))
  expect_identical(colnames(mapped$locations), c("a", "b"))
  expect_identical(as.numeric(mapped$locations), exp(c(.5, -1)))
  expect_identical(mapped$mass, .3)
  expect_identical(mapped$marginals$a$mass, .3)
  custom <- list(fun = function(x) x + 1, inv = function(x) x - 1, jac = function(x) rep(1, length(x)))
  expect_identical(as.numeric(.posterior_atoms_transform(atoms, custom)$locations), c(1.5, 0))
  support <- posterior_support_attribute(c(-Inf, Inf))
  expect_identical(.posterior_support_transform(support, "exp")$bounds, c(0, Inf))
  expect_identical(.posterior_support_transform(support, "tanh")$bounds, c(-1, 1))
  expect_null(.posterior_support_transform(support, custom))
  expect_error(.posterior_atoms_transform(.posterior_atoms_new(c(1e308, 1), c(.2, .3)), "lin", list(b = 2)),
    class = "BayesTools_transformation_image_unavailable")
  expect_identical(as.numeric(.posterior_atoms_transform(.posterior_atoms_new(1, .3), "lin", list(a = -1))$locations), 0)
})

test_that("R133 named replay allocation and joint point producers check declared images", {
  # Synthetic declared-source metadata, no fitted model or adaptive provenance.
  source <- function(location) list(model = c(1L, 2L), models = list(
    list(parameterization = "ordinary", prior = prior("point", list(location))),
    list(parameterization = "ordinary", prior = prior("normal", list(0, 1)))),
    view_transformations = list(list(transformation = "exp", arguments = NULL)))
  expect_error(.bt_ordered_source_project(source(-1000), 1), class = "BayesTools_transformation_image_unavailable")
  replay <- .bt_ordered_source_project(source(-1), 1)
  expect_identical(replay$atom, c(exp(-1), NA_real_))
  expect_identical(replay$state, c("point", "unavailable"))
  leaf <- list(prior = prior("point", list(1e-300)), point = 1e-300)
  expect_error(.prior_allocation_leaf_route(leaf, list(type = "point", location = 1), 64, list(a = 0, b = 2)),
    class = "BayesTools_transformation_image_unavailable")
  leaf$point <- 0
  expect_identical(.prior_allocation_leaf_route(leaf, list(type = "point", location = 1), 64, list(a = 0, b = 2))$locations, 0)
  plan <- list(components = matrix(1L, 1L, 1L, dimnames = list(NULL, "x")), probabilities = 1, model_mixture = FALSE)
  design <- matrix(1, 1, 1, dimnames = list("value", "x"))
  expect_error(.posterior_atoms_joint_linear(list(x = prior("point", list(-1000))), plan, design,
    output_transforms = c(value = "exp")), class = "BayesTools_transformation_image_unavailable")
  mapped <- .posterior_atoms_joint_linear(list(x = prior("point", list(-1))), plan, design,
    output_transforms = c(value = "exp"))
  expect_identical(as.numeric(mapped$locations), exp(-1))
  expect_identical(mapped$mass, 1)
})

test_that("R133 region indicators refuse positive-source underflow and preserve true-zero comparisons", {
  region <- function(inclusive) list(intervals = .prior_region_intervals(-Inf, 0),
    indicator = function(x) if(inclusive) x <= 0 else x < 0)
  evaluate <- function(location, inclusive) .prior_region_transformed("exp_lin", list(b = 2), region(inclusive),
    function(r) .prior_region_atoms(location, 1, r), function() c(0, Inf))
  expect_error(evaluate(1e-300, TRUE), class = "BayesTools_transformation_image_unavailable")
  expect_identical(evaluate(0, TRUE)$probability, 1)
  expect_identical(evaluate(0, FALSE)$probability, 0)
})

.image_factor_posterior <- function(){

  formula <- JAGS_formula(~ fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = list(intercept = prior("normal", list(0, 1)),
                      fac = prior_factor("normal", list(0, 1), contrast = "treatment")))
  fit <- coda::mcmc(cbind(mu_intercept = seq(-1, 1, length.out = 201),
    "mu_fac[1]" = seq(-1, 2, length.out = 201),
    "mu_fac[2]" = seq(-2, 1, length.out = 201)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  marginal_posterior(as_mixed_posteriors(fit, "mu_fac"), "mu_fac",
                    use_formula = FALSE, prior_samples = TRUE)
}

test_that("R133 D3 transports the specific image leaf while computing valid siblings", {
  # Synthetic compiled factor metadata and mocked prior transport, not fitted numerical evidence.
  posterior <- .image_factor_posterior()
  posterior$A <- NULL
  region_mass <- .hypothesis_region_mass
  testthat::local_mocked_bindings(.package = "BayesTools", .hypothesis_region_mass = function(quantity, side, prior){
    if(prior && identical(quantity$label, "mu_fac[B]")){
      region <- list(intervals = .prior_region_intervals(-Inf, 0), indicator = function(x) x <= 0)
      result <- .prior_region_transformed("exp_lin", list(b = 2), region,
        function(r) .prior_region_atoms(1e-300, 1, r), function() c(0, Inf))
      return(result$probability)
    }
    region_mass(quantity, side, prior)
  })
  result <- hypothesis_BF(posterior, hypothesis = "mu_fac > 0", parameter = "mu_fac", columns = "all")
  expect_identical(result$method[1L], "unavailable")
  expect_true(is.na(as.numeric(result$BF[1L])))
  expect_true(is.finite(as.numeric(result$BF[2L])))
  expect_match(attr(result, "warnings")[["mu_fac[B]"]], "finite source values", fixed = TRUE)
  expect_error(hypothesis_BF(.hypothesis_marginal_child(posterior$B), hypothesis = "mu_fac > 0", parameter = "mu_fac"),
    class = "BayesTools_transformation_image_unavailable")
})
