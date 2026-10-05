skip_if_not_test_profile("unit")

test_that("prior density ordinate validates input and has a stable schema", {

  normal_prior <- prior("normal", list(mean = 0, sd = 1))
  out <- prior_density_ordinate(normal_prior, 0)

  expect_s3_class(out, "prior_density_ordinate")
  expect_identical(
    names(unclass(out)),
    c(
      "schema_version", "value", "behavior", "log_density", "point_mass",
      "exact", "method", "reason", "provenance"
    )
  )
  expect_identical(out$schema_version, "1")
  expect_identical(out$value, 0)
  expect_identical(out$behavior, "regular")
  expect_equal(out$log_density, stats::dnorm(0, log = TRUE))
  expect_identical(out$point_mass, 0)
  expect_true(out$exact)
  expect_identical(out$method, "primitive")
  expect_null(out$reason)
  expect_true(out$method %in% c(
    "primitive", "point", "finite_mixture", "scalar_affine",
    "linear_normal", "named_transform", "unsupported_provenance"
  ))

  expect_error(
    prior_density_ordinate(normal_prior, NA_real_),
    "cannot contain NA/NaN",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(normal_prior, Inf),
    "must be finite",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(normal_prior, c(0, 1)),
    "must have length '1'",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(list(), 0),
    "must be a BayesTools prior or prior_linear_density object",
    fixed = TRUE
  )
})

test_that("primitive ordinates distinguish support and boundary limits", {

  truncated_normal <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  expect_identical(
    prior_density_ordinate(truncated_normal, -1)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(truncated_normal, 0)$behavior,
    "regular"
  )

  finite_boundary <- prior_density_ordinate(
    prior("exp", list(rate = 2)),
    0
  )
  expect_identical(finite_boundary$behavior, "regular")
  expect_equal(finite_boundary$log_density, log(2))

  zero_boundary <- prior_density_ordinate(
    prior("gamma", list(shape = 2, rate = 1)),
    0
  )
  expect_identical(zero_boundary$behavior, "zero")
  expect_identical(zero_boundary$log_density, -Inf)

  infinite_boundary <- prior_density_ordinate(
    prior("gamma", list(shape = 0.5, rate = 1)),
    0
  )
  expect_identical(infinite_boundary$behavior, "infinite")
  expect_identical(infinite_boundary$log_density, Inf)

  expect_identical(
    prior_density_ordinate(prior("lognormal", list(0, 1)), 0)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(prior("beta", list(0.5, 2)), 0)$behavior,
    "infinite"
  )
  expect_identical(
    prior_density_ordinate(prior("beta", list(2, 1)), 1)$behavior,
    "regular"
  )

  uniform_prior <- prior("uniform", list(a = -2, b = 3))
  expect_identical(
    prior_density_ordinate(uniform_prior, -2)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(uniform_prior, 3)$behavior,
    "regular"
  )
})

test_that("nonlocal primitive zeros and exact discrete masses are structural", {

  expect_identical(
    prior_density_ordinate(
      prior("moment", list(mode = 0.5, location = 0.25)),
      0.25
    )$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(
      prior("invmoment", list(mode = 0.5, df = 3, location = -0.25)),
      -0.25
    )$behavior,
    "zero"
  )

  point_at_zero <- prior_density_ordinate(prior("point", list(0)), 0)
  expect_identical(point_at_zero$behavior, "point_mass")
  expect_identical(point_at_zero$point_mass, 1)
  expect_identical(point_at_zero$log_density, -Inf)

  point_elsewhere <- prior_density_ordinate(prior("point", list(0)), 1)
  expect_identical(point_elsewhere$behavior, "zero")
  expect_identical(point_elsewhere$point_mass, 0)

  bernoulli <- prior_density_ordinate(
    prior("bernoulli", list(probability = 0.3)),
    1
  )
  expect_identical(bernoulli$behavior, "point_mass")
  expect_equal(bernoulli$point_mass, 0.3)

  none <- prior_density_ordinate(prior_none(), 0)
  expect_identical(none$behavior, "point_mass")
  expect_identical(none$point_mass, 1)
  expect_identical(none$method, "point")

  vector_prior <- prior("mnormal", list(mean = 0, sd = 1, K = 2))
  unsupported <- prior_density_ordinate(vector_prior, 0)
  expect_identical(unsupported$behavior, "unknown")
  expect_false(unsupported$exact)
})

test_that("finite mixtures combine continuous behavior and atoms exactly", {

  make_mixture <- function(priors){
    prior_mixture(
      priors,
      is_null = rep(FALSE, length(priors))
    )
  }

  regular <- make_mixture(list(
    prior("normal", list(0, 1), prior_weights = 1),
    prior("gamma", list(2, 1), prior_weights = 1)
  ))
  regular_out <- prior_density_ordinate(regular, 0)
  expect_identical(regular_out$behavior, "regular")
  expect_equal(regular_out$log_density, log(0.5) + stats::dnorm(0, log = TRUE))

  zero <- make_mixture(list(
    prior("gamma", list(2, 1), prior_weights = 1),
    prior("beta", list(2, 2), prior_weights = 1)
  ))
  expect_identical(prior_density_ordinate(zero, 0)$behavior, "zero")

  infinite <- make_mixture(list(
    prior("gamma", list(0.5, 1), prior_weights = 1),
    prior("normal", list(0, 1), prior_weights = 1)
  ))
  expect_identical(
    prior_density_ordinate(infinite, 0)$behavior,
    "infinite"
  )

  mixed_measure <- make_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("normal", list(0, 1), prior_weights = 3)
  ))
  mixed_out <- prior_density_ordinate(mixed_measure, 0)
  expect_identical(mixed_out$behavior, "point_mass")
  expect_equal(mixed_out$point_mass, 0.25)
  expect_equal(
    mixed_out$log_density,
    log(0.75) + stats::dnorm(0, log = TRUE)
  )
  expect_identical(
    mixed_out$provenance$continuous_behavior,
    "regular"
  )

  point_and_infinite <- make_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("gamma", list(0.5, 1), prior_weights = 1)
  ))
  point_and_infinite_out <- prior_density_ordinate(point_and_infinite, 0)
  expect_identical(point_and_infinite_out$behavior, "point_mass")
  expect_equal(point_and_infinite_out$point_mass, 0.5)
  expect_identical(
    point_and_infinite_out$provenance$continuous_behavior,
    "infinite"
  )

  zero_weight <- make_mixture(list(
    prior("normal", list(0, 1), prior_weights = 1),
    prior("gamma", list(0.5, 1), prior_weights = 1)
  ))
  attr(zero_weight, "prior_weights") <- c(1, 0)
  zero_weight_out <- prior_density_ordinate(zero_weight, 0)
  expect_identical(zero_weight_out$behavior, "regular")
  expect_equal(zero_weight_out$log_density, stats::dnorm(0, log = TRUE))

  spike_and_slab <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0.25))
  )
  spike_and_slab_out <- prior_density_ordinate(spike_and_slab, 0)
  expect_identical(spike_and_slab_out$behavior, "point_mass")
  expect_equal(spike_and_slab_out$point_mass, 0.75)
  expect_equal(
    spike_and_slab_out$log_density,
    log(0.25) + stats::dnorm(0, log = TRUE)
  )
})

test_that("finite-mixture structural precedence is deterministic", {

  result <- BayesTools:::.prior_density_ordinate_result
  undefined <- result(
    0, "undefined", NA_real_, method = "named_transform"
  )
  infinite <- result(
    0, "infinite", Inf, method = "primitive"
  )
  unknown <- result(
    0, "unknown", NA_real_, exact = FALSE,
    method = "unsupported_provenance"
  )

  expect_identical(
    BayesTools:::.prior_density_ordinate_combine(
      list(undefined, infinite),
      c(1, 1),
      0
    )$behavior,
    "undefined"
  )
  unknown_and_infinite <- BayesTools:::.prior_density_ordinate_combine(
    list(unknown, infinite),
    c(1, 1),
    0
  )
  expect_identical(unknown_and_infinite$behavior, "infinite")
  expect_identical(unknown_and_infinite$log_density, Inf)
  expect_true(unknown_and_infinite$exact)
})

test_that("nonzero scalar affine transformations retain behavior and Jacobian", {

  normal_prior <- prior("normal", list(mean = 1, sd = 2))
  positive <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = 2),
    n_grid = 128
  )
  negative <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = -2),
    n_grid = 128
  )

  positive_out <- prior_density_ordinate(positive, 2)
  negative_out <- prior_density_ordinate(negative, -2)
  expect_identical(positive_out$behavior, "regular")
  expect_identical(negative_out$behavior, "regular")
  expect_identical(positive_out$method, "scalar_affine")
  expect_identical(negative_out$method, "scalar_affine")
  expect_equal(
    positive_out$log_density,
    stats::dnorm(1, mean = 1, sd = 2, log = TRUE) - log(2)
  )
  expect_equal(negative_out$log_density, positive_out$log_density)

  point <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("point", list(location = 0.1))),
    weights = c(beta = 3),
    n_grid = 128
  )
  point_out <- prior_density_ordinate(point, point$points$x)
  expect_identical(point_out$behavior, "point_mass")
  expect_identical(point_out$point_mass, 1)

  zero <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = 0),
    n_grid = 128
  )
  zero_out <- prior_density_ordinate(zero, 0)
  expect_identical(zero_out$behavior, "point_mass")
  expect_identical(zero_out$point_mass, 1)

  affine_mixture <- prior_mixture(
    list(
      prior("point", list(0), prior_weights = 1),
      prior("normal", list(0, 1), prior_weights = 1)
    ),
    is_null = c(TRUE, FALSE)
  )
  affine_mixture_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = affine_mixture),
    weights = c(beta = 2),
    n_grid = 128
  )
  affine_mixture_out <- prior_density_ordinate(affine_mixture_density, 0)
  expect_identical(affine_mixture_out$behavior, "point_mass")
  expect_equal(affine_mixture_out$point_mass, 0.5)
  expect_identical(
    affine_mixture_out$provenance$continuous_behavior,
    "regular"
  )
})

test_that("scalar affine ordinates at a bound of the mapped support are the one-sided limit inside it", {

  # offset + w S with a scalar prior S and a point term o: the mapped bounds
  # are rounded and their inverse image (o + w - o) / w can round just outside
  # the support of S (above 1 for S ~ U(0, 1) with (w, o) = (.1, .2), (.3, .1),
  # (.3, .7)), where the ordinate was 0. The ordinate at a bound is the
  # one-sided limit inside the support, the density of S at its bound divided
  # by |w| (analytic references); outside the support it stays 0
  ulp <- function(x) 2^(floor(log2(abs(x))) - 52)
  affine <- function(source, w, o, source_transform = NULL){
    BayesTools:::.prior_linear_combination_density(
      prior_list        = list(s = source, p = prior("point", list(o))),
      weights           = c(s = w, p = 1),
      source_transforms = if(!is.null(source_transform)) c(s = source_transform, p = NA_character_),
      n_grid            = 64
    )
  }
  # the density of S at each finite source bound (the log source Z = log(S)
  # has the density f_S(s) s at log(s))
  normal_mass <- stats::pnorm(2) - stats::pnorm(-1)
  sources <- list(
    list(prior = prior("uniform", list(0, 1)), bounds = c(0, 1), limits = c(1, 1)),
    list(prior = prior("uniform", list(1, 3)), bounds = c(1, 3), limits = c(.5, .5)),
    list(prior = prior("beta", list(1, 3)), bounds = c(0, 1), limits = c(3, 0)),
    list(prior = prior("beta", list(3, 1)), bounds = c(0, 1), limits = c(0, 3)),
    list(prior = prior("normal", list(0, 1), list(-1, 2)), bounds = c(-1, 2),
         limits = stats::dnorm(c(-1, 2)) / normal_mass),
    list(prior = prior("exp", list(2)), bounds = 0, limits = 2)
  )
  expect_bound_limits <- function(source, w, o, log_source = FALSE){
    density <- affine(source$prior, w, o, if(log_source) "log")
    positive <- source$bounds > 0
    bounds <- if(log_source) log(source$bounds[positive]) else source$bounds
    limits <- source$limits[if(log_source) positive else TRUE]
    if(log_source) limits <- limits * source$bounds[positive]
    mapped <- o + w * bounds
    for(i in seq_along(bounds)){
      value <- mapped[i]
      label <- sprintf("%s, w = %g, o = %g, bound %g", source$prior$distribution, w, o, bounds[i])
      ordinate <- prior_density_ordinate(density, value)
      expect_identical(ordinate$method, "scalar_affine", info = label)
      expect_true(ordinate$exact, info = label)
      if(limits[i] > 0){
        expect_identical(ordinate$behavior, "regular", info = label)
        expect_equal(exp(ordinate$log_density), limits[i] / abs(w), tolerance = 1e-12, info = label)
      }else{
        expect_identical(ordinate$behavior, "zero", info = label)
        expect_identical(ordinate$log_density, -Inf, info = label)
      }
      # outside the support (4 ulps of the largest operand of the mapping, which
      # the rounding of the inverse image is in, and 1e-9 beyond the mapped
      # bound) the ordinate is 0; a one-sided bound has its outside on one side
      # only
      outward <- if(length(bounds) == 2L) {
        if(value == max(mapped)) 1 else -1
      }else{
        if(w > 0) -1 else 1
      }
      operands <- max(abs(c(value, o, w * bounds[i])))
      for(beyond in value + outward * c(1e-9 * max(1, abs(value)), 4 * ulp(operands))){
        outside <- prior_density_ordinate(density, beyond)
        expect_identical(outside$behavior, "zero", info = label)
        expect_identical(outside$log_density, -Inf, info = label)
        expect_true(outside$exact, info = label)
      }
    }
  }

  # S ~ U(0, 1): both bounds of 32 combinations (the three rounded ones among
  # them), and each source on a subset
  for(w in c(.1, .3, .7, 1.3, -.1, -.3, -.7, -1.3)){
    for(o in c(.1, .2, .7, -.35)){
      expect_bound_limits(sources[[1L]], w, o)
    }
  }
  for(source in sources[-1L]){
    for(w in c(.1, .3, 1.3, -.3, -.7)){
      for(o in c(.15, .2, .7)){
        expect_bound_limits(source, w, o)
      }
    }
  }
  # the three rounded combinations of the brief, by name
  for(combination in list(c(.1, .2), c(.3, .1), c(.3, .7))){
    density <- affine(prior("uniform", list(0, 1)), combination[1L], combination[2L])
    expect_equal(
      exp(prior_density_ordinate(density, combination[2L] + combination[1L])$log_density),
      1 / combination[1L], tolerance = 1e-12
    )
  }

  # a log source term: the bounds of S ~ U(5, 9) and U(1, 3) are rounded by
  # exp(log(s)) as well (it is below 5 and above 9 and 3)
  for(source in list(
    list(prior = prior("uniform", list(5, 9)), bounds = c(5, 9), limits = c(.25, .25)),
    sources[[2L]]
  )){
    for(w in c(.1, .3, .7, 1.3, -.1, -.3, -.7, -1.3)){
      for(o in c(.1, .2, .7, -.35)){
        expect_bound_limits(source, w, o, log_source = TRUE)
      }
    }
  }
  # no offset: the log source alone (weight 1) at the rounded bound
  log_alone <- BayesTools:::.prior_linear_combination_density(
    list(s = prior("uniform", list(1, 3))), c(s = 1), source_transforms = c(s = "log")
  )
  expect_equal(exp(prior_density_ordinate(log_alone, log(3))$log_density), 1.5, tolerance = 1e-12)

  # the adjacent double of a mapped bound is that bound, as for a named
  # linear transformation of the same route (the mapped bound is rounded, so
  # its neighbour cannot be told from it)
  source_route <- BayesTools:::.prior_density_route_linear(
    list(s = prior("uniform", list(0, 1))), c(s = 1), NULL
  )
  for(combination in list(c(.1, .2), c(.3, .1), c(.3, .7), c(.7, .2))){
    w <- combination[1L]
    o <- combination[2L]
    density <- affine(prior("uniform", list(0, 1)), w, o)
    lin <- BayesTools:::.prior_density_route_transform(
      source_route, "lin", list(a = o, b = w), function() c(o, o + w)
    )
    upper <- o + w
    for(k in c(-1, 0, 1)){
      value <- upper + k * ulp(upper)
      scalar <- prior_density_ordinate(density, value)
      named <- BayesTools:::.prior_density_route_ordinate(lin, value)
      label <- sprintf("w = %g, o = %g, %+d ulp", w, o, k)
      expect_identical(scalar$behavior, "regular", info = label)
      expect_identical(named$behavior, "regular", info = label)
      expect_equal(scalar$log_density, named$log_density, tolerance = 1e-12, info = label)
    }
  }

  # a mixture slab keeps its own bounds: each component is matched against
  # its own support
  mixture <- prior_mixture(
    list(
      prior("uniform", list(0, 1), prior_weights = 1),
      prior("uniform", list(0, 2), prior_weights = 1)
    ),
    is_null = c(FALSE, FALSE)
  )
  density <- affine(mixture, .3, .1)
  expect_equal(
    exp(prior_density_ordinate(density, .1 + .3 * 1)$log_density),
    (.5 * 1 + .5 * .5) / .3, tolerance = 1e-12
  )
  expect_equal(
    exp(prior_density_ordinate(density, .1 + .3 * 2)$log_density),
    (.5 * 0 + .5 * .5) / .3, tolerance = 1e-12
  )
})

test_that("normal linear combinations are classified analytically", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid = 128
  )
  out <- prior_density_ordinate(density, 0)
  expect_identical(out$behavior, "regular")
  expect_identical(out$method, "linear_normal")
  expect_equal(
    out$log_density,
    stats::dnorm(0, mean = -0.5, sd = sqrt(5), log = TRUE)
  )
  expect_identical(out$provenance$weights, c(x = 2, y = -0.5))

  vector_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      beta = prior("mnormal", list(mean = 1, sd = 2, K = 2))
    ),
    weights = c("beta[1]" = 1, "beta[2]" = -0.5),
    n_grid = 128
  )
  vector_out <- prior_density_ordinate(vector_density, 0)
  expect_identical(vector_out$behavior, "regular")
  expect_identical(vector_out$method, "linear_normal")
  expect_equal(
    vector_out$log_density,
    stats::dnorm(0, mean = 0.5, sd = sqrt(5), log = TRUE)
  )

  shifted <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      shift = prior("point", list(1.5)),
      x = prior("normal", list(1, 2)),
      y = prior("normal", list(-1, 1))
    ),
    weights = c(shift = 2, x = -0.5, y = 2),
    n_grid = 128
  )
  shifted_out <- prior_density_ordinate(shifted, 0.5)
  expect_identical(shifted_out$behavior, "regular")
  expect_identical(shifted_out$method, "linear_normal")
  expect_equal(
    shifted_out$log_density,
    stats::dnorm(0.5, mean = 0.5, sd = sqrt(5), log = TRUE)
  )
})

test_that("density contexts retain structural ordinate provenance", {

  context <- BayesTools:::.prior_density_context(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    column_names = c("x", "y"),
    n_grid = 128
  )
  density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(x = 2, y = -0.5)
  )
  out <- prior_density_ordinate(density, 0)
  expect_identical(out$behavior, "regular")
  expect_identical(out$method, "linear_normal")
  expect_identical(out$provenance$context$kind, "prior_density_context")
  expect_identical(
    out$provenance$context$standardized_weights,
    c(x = 2, y = -0.5)
  )

  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, 0, 2, 0),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("x", "y"))
    )
  )
  row_out <- prior_density_ordinate(row_density, 0)
  expect_identical(row_out$behavior, "regular")
  expect_identical(row_out$method, "finite_mixture")
  expect_identical(row_out$provenance$unique_rows, 2L)
  expect_identical(row_out$provenance$weight_dimensions, c(2L, 2L))
  expect_match(row_out$provenance$weights_hash, "^[0-9a-f]{8}$")
  expect_false("representative_weights" %in% names(row_out$provenance))

  model_context <- BayesTools:::.prior_density_model_mixture_context(
    prior_list = list(
      x = list(
        prior("point", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 3)
      )
    ),
    column_names = "x",
    n_grid = 128
  )
  model_density <- BayesTools:::.prior_density_from_context(
    model_context,
    weights = c(x = 1)
  )
  model_out <- prior_density_ordinate(model_density, 0)
  expect_identical(model_out$behavior, "point_mass")
  expect_equal(model_out$point_mass, 0.25)
  expect_equal(
    model_out$log_density,
    log(0.75) + stats::dnorm(0, log = TRUE)
  )
  expect_identical(model_out$provenance$context, "model_mixture")
})

test_that("named monotone transformations classify interiors and boundaries", {

  normal_prior <- list(beta = prior("normal", list(0, 1)))
  exponential <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp"
  )
  exp_interior <- prior_density_ordinate(exponential, 1)
  expect_identical(exp_interior$behavior, "regular")
  expect_identical(exp_interior$method, "named_transform")
  expect_equal(exp_interior$log_density, stats::dnorm(0, log = TRUE))
  expect_identical(
    prior_density_ordinate(exponential, 0)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(exponential, -1)$behavior,
    "zero"
  )

  hyperbolic <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "tanh"
  )
  expect_equal(
    prior_density_ordinate(hyperbolic, 0)$log_density,
    stats::dnorm(0, log = TRUE)
  )
  expect_identical(
    prior_density_ordinate(hyperbolic, 1)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(hyperbolic, 2)$behavior,
    "zero"
  )

  linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "lin",
    output_transformation_arguments = list(a = 1, b = -2)
  )
  expect_equal(
    prior_density_ordinate(linear, 1)$log_density,
    stats::dnorm(0, log = TRUE) - log(2)
  )

  positive_prior <- list(
    beta = prior(
      "lognormal",
      list(0, 1),
      truncation = list(lower = 0.01, upper = Inf)
    )
  )
  exponential_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = positive_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = 2)
  )
  expect_identical(
    prior_density_ordinate(exponential_linear, 1)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(exponential_linear, 0)$behavior,
    "zero"
  )

  constant_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "lin",
    output_transformation_arguments = list(a = 2, b = 0)
  )
  constant_linear_out <- prior_density_ordinate(constant_linear, 2)
  expect_identical(constant_linear_out$behavior, "point_mass")
  expect_identical(constant_linear_out$point_mass, 1)

  constant_exp_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 2, b = 0)
  )
  constant_exp_linear_out <- prior_density_ordinate(
    constant_exp_linear,
    exp(2)
  )
  expect_identical(constant_exp_linear_out$behavior, "point_mass")
  expect_identical(constant_exp_linear_out$point_mass, 1)
})

test_that("composed named-transform boundary limits use source provenance", {

  lognormal_power <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 2),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "exp"
  )
  expect_identical(
    prior_density_ordinate(lognormal_power, 0)$behavior,
    "zero"
  )

  gamma_power <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("gamma", list(1, 1))),
    weights = c(x = 2),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "exp"
  )
  expect_identical(
    prior_density_ordinate(gamma_power, 0)$behavior,
    "infinite"
  )

  inverse_gamma_power <- function(shape){
    BayesTools:::.prior_linear_combination_density(
      prior_list = list(x = prior(
        "invgamma",
        list(shape = shape, scale = 1),
        truncation = list(lower = 0.01, upper = Inf)
      )),
      weights = c(x = -2),
      source_transforms = c(x = "log"),
      n_grid = 128,
      output_transformation = "exp"
    )
  }
  boundary_behaviors <- vapply(c(0.5, 2, 3), function(shape){
    prior_density_ordinate(inverse_gamma_power(shape), 0)$behavior
  }, character(1))
  expect_identical(boundary_behaviors, c("infinite", "regular", "zero"))
  expect_identical(
    prior_density_ordinate(inverse_gamma_power(3), -1)$behavior,
    "zero"
  )

  shape <- 3
  scale <- 1
  lower <- 0.01
  power <- -2
  density <- inverse_gamma_power(shape)
  interior <- 0.25
  source_value <- interior^(1 / power)
  log_normalizer <- stats::pgamma(
    1 / lower,
    shape = shape,
    rate = scale,
    log.p = TRUE
  )
  expected_log_density <-
    shape * log(scale) - lgamma(shape) -
    (shape + 1) * log(source_value) - scale / source_value -
    log_normalizer + log(abs(1 / power)) +
    (1 / power - 1) * log(interior)
  expect_equal(
    prior_density_ordinate(density, interior)$log_density,
    expected_log_density,
    tolerance = 1e-12
  )
  expect_identical(
    prior_density_ordinate(density, lower^power * (1 - 1e-12))$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(density, lower^power + 1)$behavior,
    "zero"
  )

  bounded_power <- function(power){
    BayesTools:::.prior_linear_combination_density(
      prior_list = list(x = prior(
        "lognormal",
        list(0, 1),
        truncation = list(lower = 0.01, upper = 2)
      )),
      weights = c(x = power),
      source_transforms = c(x = "log"),
      n_grid = 128,
      output_transformation = "exp"
    )
  }
  increasing <- bounded_power(2)
  decreasing <- bounded_power(-2)
  expect_identical(
    vapply(c(0.01^2, 2^2), function(endpoint){
      prior_density_ordinate(increasing, endpoint)$behavior
    }, character(1)),
    c("regular", "regular")
  )
  expect_identical(
    vapply(c(2^-2, 0.01^-2), function(endpoint){
      prior_density_ordinate(decreasing, endpoint)$behavior
    }, character(1)),
    c("regular", "regular")
  )

  bounded_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior(
      "normal",
      list(0, 1),
      truncation = list(lower = 0.1, upper = 1)
    )),
    weights = c(x = 1),
    n_grid = 128,
    output_transformation = "lin",
    output_transformation_arguments = list(a = 0.1, b = 2)
  )
  expect_false(identical(0.1 + 2 * 0.1, 0.3))
  expect_identical(
    prior_density_ordinate(bounded_linear, 0.3)$behavior,
    "regular"
  )

  lognormal_tanh <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 1),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "tanh"
  )
  expect_identical(
    prior_density_ordinate(lognormal_tanh, 1)$behavior,
    "zero"
  )

  student_tanh <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("t", list(0, 1, 3))),
    weights = c(x = 1),
    n_grid = 128
  )
  student_adaptive <- attr(student_tanh, "adaptive_evaluation", exact = TRUE)
  student_adaptive$arguments$output_transformation <- "tanh"
  attr(student_tanh, "adaptive_evaluation") <- student_adaptive
  expect_identical(
    prior_density_ordinate(student_tanh, 1)$behavior,
    "infinite"
  )

  positive_t <- prior(
    "t",
    list(location = 0, scale = 1, df = 3),
    truncation = list(lower = 0, upper = Inf)
  )
  dependent_boundary <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = positive_t),
    weights = c(x = 1),
    n_grid = 128
  )
  dependent_adaptive <- attr(
    dependent_boundary,
    "adaptive_evaluation",
    exact = TRUE
  )
  dependent_adaptive$arguments$output_transformation <- "exp_lin"
  dependent_adaptive$arguments$output_transformation_arguments <-
    list(a = 0, b = 1)
  attr(dependent_boundary, "adaptive_evaluation") <- dependent_adaptive
  dependent_out <- prior_density_ordinate(dependent_boundary, 0)
  expect_identical(dependent_out$behavior, "unknown")
  expect_false(dependent_out$exact)
})

test_that("unsupported transformations, convolutions, and products stay unknown", {

  custom_identity <- list(
    fun = function(x) x,
    inv = function(x) x,
    jac = function(x) rep(1, length(x))
  )
  custom <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("normal", list(0, 1))),
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = custom_identity
  )
  custom_out <- prior_density_ordinate(custom, 0)
  expect_identical(custom_out$behavior, "unknown")
  expect_false(custom_out$exact)

  # three non-normal terms have no structural convolution route
  convolution <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("gamma", list(2, 1)),
      y = prior("gamma", list(2, 1)),
      z = prior("gamma", list(2, 1))
    ),
    weights = c(x = 1, y = 1, z = 1),
    n_grid = 128
  )
  convolution_out <- prior_density_ordinate(convolution, 1)
  expect_identical(convolution_out$behavior, "unknown")
  expect_false(convolution_out$exact)
  expect_true(is.finite(convolution_out$log_density))

  # a product with a non-normal additive term has no structural route
  product_priors <- list(
    alpha = prior("t", list(0, 1, 5)),
    beta = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1))
  )
  attr(product_priors$beta, "multiply_by") <- "sigma"
  general_product <- BayesTools:::.prior_linear_combination_density(
    prior_list = product_priors,
    weights = c(alpha = 1, beta = 1),
    n_grid = 128
  )
  expect_identical(prior_density_ordinate(general_product, 1)$behavior, "unknown")
  expect_false(prior_density_ordinate(general_product, 1)$exact)
})

test_that("two-term convolutions and Cauchy sums have structural ordinates", {

  combination <- function(priors, ...){
    weights <- c(...)
    names(weights) <- names(priors)
    .prior_linear_combination_density(priors, weights)
  }
  ordinate_height <- function(density, value){
    exp(prior_density_ordinate(density, value)$log_density)
  }
  probability <- function(density, hypothesis){
    side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
    .hypothesis_prior_density_prob(density, side, "theta")
  }
  split_integral <- function(f, points){
    sum(vapply(seq_len(length(points) - 1L), function(i){
      stats::integrate(f, points[i], points[i + 1L], rel.tol = 1e-12, abs.tol = 0,
                       subdivisions = 1000L)$value
    }, numeric(1)))
  }

  # Cauchy terms sum to a Cauchy term with the summed scales
  cauchy <- combination(list(a = prior("cauchy", list(0, 1)), b = prior("cauchy", list(0, .01))), 1, -1)
  for(value in c(0, .5, 3)){
    ordinate <- prior_density_ordinate(cauchy, value)
    expect_identical(ordinate$method, "scalar_affine")
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), stats::dcauchy(value, 0, 1.01), tolerance = 1e-14)
  }
  expect_equal(probability(cauchy, "theta > 1000"),
               stats::pcauchy(1000, 0, 1.01, lower.tail = FALSE), tolerance = 1e-12)
  transformed <- .prior_linear_combination_density(
    list(a = prior("cauchy", list(0, 1)), b = prior("cauchy", list(0, 1))), c(a = 1, b = 1),
    output_transformation = "exp"
  )
  expect_equal(ordinate_height(transformed, 5), stats::dcauchy(log(5), 0, 2) / 5, tolerance = 1e-14)

  # two t3 terms: quadrature against hand-split integrals, also in the far tail
  student <- combination(list(a = prior("t", list(0, 1, 3)), b = prior("t", list(0, 1, 3))), 1, 1)
  for(value in c(0, 2, 20)){
    ordinate <- prior_density_ordinate(student, value)
    expect_identical(ordinate$method, "convolution")
    expect_true(ordinate$exact)
    expect_true(ordinate$provenance$integration$converged)
    expect_equal(exp(ordinate$log_density), split_integral(
      function(t) stats::dt(t, 3) * stats::dt(value - t, 3), c(-Inf, -10, 0, value, value + 10, Inf)
    ), tolerance = 1e-10)
  }
  expect_equal(probability(student, "theta > 100"), split_integral(
    function(t) stats::dt(t, 3) * stats::pt(100 - t, 3, lower.tail = FALSE), c(-Inf, -10, 0, 100, 110, Inf)
  ), tolerance = 1e-10)

  # gamma(1/2, 1) - gamma(1/2, 1) has the density K0(|v|) / pi, infinite at 0
  # where the two singular bounds meet
  gamma_difference <- combination(list(a = prior("gamma", list(.5, 1)), b = prior("gamma", list(.5, 1))), 1, -1)
  expect_identical(prior_density_ordinate(gamma_difference, 0)$behavior, "infinite")
  for(value in c(-1, .01, 1)){
    expect_equal(ordinate_height(gamma_difference, value), besselK(abs(value), 0) / pi, tolerance = 1e-9)
  }
  expect_equal(probability(gamma_difference, "theta > 1"),
               split_integral(function(z) besselK(z, 0) / pi, c(1, 10, Inf)), tolerance = 1e-9)

  # meeting bounds: two arcsine (beta(1/2, 1/2)) terms are infinite where
  # bounds meet inside the support; their positive limits at the ends of the
  # support are not classified. Two uniform terms are triangular, zero at the
  # ends; two gamma(2, 1) terms are gamma(4, 1).
  arcsine <- combination(list(a = prior("beta", list(.5, .5)), b = prior("beta", list(.5, .5))), 1, 1)
  expect_identical(prior_density_ordinate(arcsine, 1)$behavior, "infinite")
  expect_identical(prior_density_ordinate(arcsine, 0)$behavior, "unknown")
  expect_identical(prior_density_ordinate(arcsine, 2)$behavior, "unknown")
  expect_equal(ordinate_height(arcsine, .5), stats::integrate(function(u){
    t <- .5 * sin(u)^2
    stats::dbeta(t, .5, .5) * stats::dbeta(.5 - t, .5, .5) * sin(u) * cos(u)
  }, 0, pi / 2, rel.tol = 1e-12)$value, tolerance = 1e-8)
  uniform <- combination(list(a = prior("uniform", list(0, 1)), b = prior("uniform", list(0, 1))), 1, 1)
  expect_identical(prior_density_ordinate(uniform, 0)$behavior, "zero")
  expect_identical(prior_density_ordinate(uniform, 2)$behavior, "zero")
  for(value in c(.5, 1, 1.5)){
    expect_equal(ordinate_height(uniform, value), 1 - abs(value - 1), tolerance = 1e-12)
  }
  expect_equal(probability(uniform, "theta > 1.9"), .005, tolerance = 1e-12)
  gamma_sum <- combination(list(a = prior("gamma", list(2, 1)), b = prior("gamma", list(2, 1))), 1, 1)
  expect_identical(prior_density_ordinate(gamma_sum, 0)$behavior, "zero")
  for(value in c(1, 5)){
    expect_equal(ordinate_height(gamma_sum, value), stats::dgamma(value, 4, 1), tolerance = 1e-12)
  }

  # the convolution is a component of mixtures: t3 + spike-and-slab(t3)
  slab <- prior_spike_and_slab(prior("t", list(0, 1, 3)), prior_inclusion = prior("spike", list(.5)))
  mixture <- combination(list(a = prior("t", list(0, 1, 3)), b = slab), 1, 1)
  ordinate <- prior_density_ordinate(mixture, .5)
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density), .5 * stats::dt(.5, 3) + .5 * split_integral(
    function(t) stats::dt(t, 3) * stats::dt(.5 - t, 3), c(-Inf, -10, 0, .5, 10.5, Inf)
  ), tolerance = 1e-10)
})

test_that("a truncated normal term plus normal terms has a closed-form ordinate", {

  # X = offset + w T + G with T a truncated normal and G the sum of the normal
  # terms. References: 40-digit mpmath quadratures of f_T(t) phi(x - offset -
  # w t; 0, sd(G)) over T's support, split into 200 pieces around the mass
  # (reported relative error estimates below 1e-40), independent of the
  # closed form; log densities agree within 1e-12 (a relative density error
  # of 1e-12), also far in both tails.
  half <- prior("normal", list(0.5, 0.35), list(0, Inf))
  cases <- list(
    list(priors = list(t = prior("normal", list(0, 1), list(0, Inf)), g = prior("normal", list(0, 1))),
         weights = c(t = 1, g = 1), value = .3, log_reference = -1.13272268820496627806),
    list(priors = list(t = half, g = prior("normal", list(.3, .3)), h = prior("normal", list(0, .4))),
         weights = c(t = -2, g = 1, h = 1), value = -.7, log_reference = -0.6957464505949913255644),
    list(priors = list(t = prior("normal", list(0, 1), list(-1, 2)), g = prior("normal", list(0, .1)),
                       d = prior("point", list(1))),
         weights = c(t = 1, g = 1, d = 1), value = 2.5, log_reference = -1.837608904894240190235),
    list(priors = list(t = prior("normal", list(1, 2), list(-Inf, 0)), g = prior("normal", list(-1, 3))),
         weights = c(t = 1, g = 1), value = -4, log_reference = -2.22940924376757506274),
    list(priors = list(t = prior("normal", list(0, 1), list(0, Inf)), g = prior("normal", list(0, .5))),
         weights = c(t = 1, g = 1), value = 12, log_reference = -57.93736312830183231025),
    list(priors = list(t = half, g = prior("normal", list(.2, .4))),
         weights = c(t = 1, g = 1), value = -9, log_reference = -269.2956623278153067464852)
  )
  for(case in cases){
    density <- .prior_linear_combination_density(case$priors, case$weights)
    ordinate <- prior_density_ordinate(density, case$value)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "truncated_normal_convolution")
    expect_identical(ordinate$provenance$kind, "truncated_normal_convolution")
    expect_null(ordinate$provenance$integration)
    expect_lt(abs(ordinate$log_density - case$log_reference), 1e-12)
    height <- as.numeric(.prior_linear_density_height(density, case$value))
    expect_lt(abs(log(height) - case$log_reference), 1e-12)
  }

  # the density integrates to one, the plotted density is the closed form, and
  # region probabilities (the Gaussian-convolution quadrature) match an
  # independent integral of T's density times the normal tail probability
  priors <- list(t = half, g = prior("normal", list(.2, .4)))
  density <- .prior_linear_combination_density(priors, c(t = 1, g = 1))
  route <- .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation"))
  expect_identical(route$type, "truncated_normal_convolution")
  closed <- function(x) .prior_density_route_density(route, x)
  expect_equal(stats::integrate(closed, -Inf, Inf, rel.tol = 1e-12)$value, 1, tolerance = 1e-10)
  values <- c(-2, -.3, 0, .7, 1.4, 3)
  expect_equal(closed(values), vapply(values, function(value){
    exp(prior_density_ordinate(density, value)$log_density)
  }, numeric(1)), tolerance = 1e-14)
  plotted <- .prior_linear_density_to_plot_data(density)$density
  expect_equal(plotted$y, closed(plotted$x), tolerance = 1e-14)
  probability <- .hypothesis_prior_density_prob(
    density, hypothesis_parse("theta > 0.3")$statements[[1L]]$left, "theta"
  )
  expect_equal(as.numeric(probability), stats::integrate(function(t){
    stats::dnorm(t, .5, .35) / stats::pnorm(.5 / .35) *
      stats::pnorm(.3 - t, .2, .4, lower.tail = FALSE)
  }, 0, Inf, rel.tol = 1e-12)$value, tolerance = 1e-9)

  # a spike-and-slab truncated normal is a mixture of the closed form and the
  # normal term alone
  slab <- prior_spike_and_slab(half, prior_inclusion = prior("spike", list(.4)))
  mixture <- .prior_linear_combination_density(list(t = slab, g = prior("normal", list(.2, .4))), c(t = 1, g = 1))
  ordinate <- prior_density_ordinate(mixture, .6)
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density),
               .4 * closed(.6) + .6 * stats::dnorm(.6, .2, .4), tolerance = 1e-14)
})

test_that("prior_density_has_provenance() signals densities without deterministic provenance", {

  region_probability <- function(density, hypothesis){
    .hypothesis_prior_density_prob(
      density, hypothesis_parse(hypothesis)$statements[[1L]]$left, "theta"
    )
  }
  expect_true(prior_density_has_provenance(prior("normal", list(0, 1))))
  normal_sum <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, 1)), b = prior("normal", list(1, 2))), c(a = 1, b = 1)
  )
  expect_true(prior_density_has_provenance(normal_sum))

  # combinations without a structural route have provenance: their ordinates
  # are unknown (method "unsupported_provenance", or "named_transform" under an
  # output transformation), while their region probabilities are refined grid
  # probabilities
  gammas <- .prior_linear_combination_density(
    list(x = prior("gamma", list(2, 1)), y = prior("gamma", list(2, 1)), z = prior("gamma", list(2, 1))),
    c(x = 1, y = 1, z = 1)
  )
  expect_true(prior_density_has_provenance(gammas))
  expect_identical(prior_density_ordinate(gammas, 6)$method, "unsupported_provenance")
  expect_equal(region_probability(gammas, "theta > 6"), stats::pgamma(6, 6, 1, lower.tail = FALSE),
               tolerance = 1e-3)
  log_intercept <- .prior_linear_combination_density(
    list(t = prior("normal", list(0, .35), list(0, Inf)), g = prior("t", list(0, .5, 3))),
    c(t = 1, g = -1), source_transforms = c(t = "log", g = NA),
    output_transformation = "exp"
  )
  expect_true(prior_density_has_provenance(log_intercept))
  ordinate <- prior_density_ordinate(log_intercept, .3)
  expect_identical(ordinate$behavior, "unknown")
  expect_identical(ordinate$method, "named_transform")

  # a density grid without its provenance record and a product without a
  # structural route are for plotting only: their heights stop
  grid <- normal_sum
  attr(grid, "adaptive_evaluation") <- NULL
  expect_false(prior_density_has_provenance(grid))
  expect_identical(prior_density_ordinate(grid, 0)$method, "unsupported_provenance")
  expect_error(.prior_linear_density_height(grid, 0), "has no deterministic provenance",
               fixed = TRUE)
  product_priors <- list(
    alpha = prior("t", list(0, 1, 5)),
    beta  = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1))
  )
  attr(product_priors$beta, "multiply_by") <- "sigma"
  product <- .prior_linear_combination_density(product_priors, c(alpha = 1, beta = 1), n_grid = 128)
  expect_false(prior_density_has_provenance(product))
  expect_error(.prior_linear_density_height(product, 1), "no structural density route",
               fixed = TRUE)

  # point masses alone are exact without a record
  points <- structure(list(density = NULL, points = data.frame(x = c(0, 1), p = c(.4, .6))),
                      class = c("prior_linear_density", "prior_density"))
  expect_true(prior_density_has_provenance(points))
  expect_error(prior_density_has_provenance(list()),
               "must be a BayesTools prior or prior_linear_density object", fixed = TRUE)
})

test_that("the exp of a log-source term plus a Gaussian part is a scale product", {

  # Y = exp(log(b0) + c b1) = b0 W, W = exp(c b1) ~ lognormal(0, |c| .5), with
  # the positive intercept b0 ~ N(0, sqrt(2) / 4)T(0, Inf) and the slope
  # b1 ~ N(0, .5) (the unscaled heterogeneity intercept of a log-intercept
  # scale regression, bangertdrowns2004). References: integrals over b1 of
  # f_b0(y e^(-c t)) e^(-c t) phi(t; 0, .5) at rel.tol 1e-13.
  intercept_sd <- sqrt(2) / 4
  slope <- -0.9894605
  exp_density <- function(intercept, slope_prior = prior("normal", list(0, .5)), weight = 1){
    .prior_linear_combination_density(
      list(b0 = intercept, b1 = slope_prior), c(b0 = weight, b1 = slope),
      source_transforms = c(b0 = "log", b1 = NA), output_transformation = "exp"
    )
  }
  reference <- function(y, log_f_intercept){
    stats::integrate(function(t){
      exp(log_f_intercept(y * exp(-slope * t)) - slope * t + stats::dnorm(t, 0, .5, log = TRUE))
    }, -Inf, Inf, rel.tol = 1e-13, abs.tol = 0)$value
  }
  f_half <- function(x) log(2) + stats::dnorm(x, 0, intercept_sd, log = TRUE)
  half <- prior("normal", list(0, intercept_sd), list(0, Inf))
  density <- exp_density(half)
  route <- .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation"))
  expect_identical(route$type, "scale_product")
  for(value in c(.05, .2, 1, 3)){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "scale_mixture")
    expect_equal(exp(ordinate$log_density), reference(value, f_half), tolerance = 1e-10)
  }
  # at 0, the image of log(b0) = -Inf: f_b0(0) E[1 / W] with W lognormal
  zero <- prior_density_ordinate(density, 0)
  expect_true(zero$exact)
  expect_equal(exp(zero$log_density), exp(f_half(0)) * exp((slope * .5)^2 / 2), tolerance = 1e-12)
  expect_identical(prior_density_ordinate(density, -.1)$behavior, "zero")
  # region probabilities against integrals of the probability of b0
  probability <- function(hypothesis){
    .hypothesis_prior_density_prob(
      density, hypothesis_parse(hypothesis)$statements[[1L]]$left, "theta"
    )
  }
  interval_reference <- function(lower, upper){
    stats::integrate(function(t){
      scale <- exp(-slope * t)
      upper_tail <- if(is.finite(upper)) stats::pnorm(upper * scale, 0, intercept_sd, lower.tail = FALSE) else 0
      2 * (stats::pnorm(lower * scale, 0, intercept_sd, lower.tail = FALSE) - upper_tail) *
        stats::dnorm(t, 0, .5)
    }, -Inf, Inf, rel.tol = 1e-13)$value
  }
  expect_equal(probability("theta > 0.1"), interval_reference(.1, Inf), tolerance = 1e-10)
  expect_equal(probability("theta > 0.2 & theta < 1"), interval_reference(.2, 1), tolerance = 1e-10)
  # the plotted density is the scale product's
  plotted <- .prior_linear_density_to_plot_data(density)$density
  checked <- plotted$x[plotted$x > 0][c(10, 60, 150)]
  expect_equal(plotted$y[match(checked, plotted$x)], vapply(checked, function(value){
    exp(prior_density_ordinate(density, value)$log_density)
  }, numeric(1)), tolerance = 1e-8)

  # a mixture prior of b0 is the mixture of the components' scale products
  gamma <- prior("gamma", list(2, 4))
  mixture <- exp_density(prior_mixture(list(half, gamma), is_null = c(FALSE, FALSE)))
  mixture_route <- .prior_density_route_from_adaptive(attr(mixture, "adaptive_evaluation"))
  expect_identical(mixture_route$type, "mixture")
  expect_identical(vapply(mixture_route$components, `[[`, character(1), "type"),
                   c("scale_product", "scale_product"))
  for(value in c(.2, 1)){
    ordinate <- prior_density_ordinate(mixture, value)
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density),
                 .5 * reference(value, f_half) +
                   .5 * reference(value, function(x) stats::dgamma(x, 2, 4, log = TRUE)),
                 tolerance = 1e-10)
  }

  # not applicable, unchanged: a non-Gaussian other term, and a log-source
  # term with another weight (b0^2 W, which the scale-product leaf does not
  # represent), remain general convolutions with unknown ordinates
  for(unchanged in list(exp_density(half, slope_prior = prior("t", list(0, .5, 3))),
                        exp_density(half, weight = 2))){
    ordinate <- prior_density_ordinate(unchanged, .5)
    expect_identical(ordinate$behavior, "unknown")
    expect_false(ordinate$exact)
    expect_identical(ordinate$reason, "General numerical convolutions are not structurally classified.")
    expect_identical(
      .prior_density_route_from_adaptive(attr(unchanged, "adaptive_evaluation"))$type,
      "transform"
    )
  }
})

test_that("the log-scale sum of a log-source term and a Gaussian part is the log image of the scale product", {

  # Z = log(b0) + c b1 = log(Y), Y = b0 exp(c b1) the scale product above: its
  # density is f_Z(z) = f_Y(e^z) e^z. References: the convolution integral over
  # b1 of the density of log(b0), f_b0(e^u) e^u at u = z - c t, times
  # phi(t; 0, .5), at rel.tol 1e-13; region probabilities integrate
  # P(log(b0) in (a - c t, b - c t)) phi(t; 0, .5).
  intercept_sd <- sqrt(2) / 4
  slope <- -0.9894605
  log_density <- function(intercept, slope_prior = prior("normal", list(0, .5)), weight = 1){
    .prior_linear_combination_density(
      list(b0 = intercept, b1 = slope_prior), c(b0 = weight, b1 = slope),
      source_transforms = c(b0 = "log", b1 = NA)
    )
  }
  reference <- function(z, log_f_intercept){
    stats::integrate(function(t){
      u <- z - slope * t
      exp(log_f_intercept(exp(u)) + u + stats::dnorm(t, 0, .5, log = TRUE))
    }, -Inf, Inf, rel.tol = 1e-13, abs.tol = 0)$value
  }
  f_half <- function(x) log(2) + stats::dnorm(x, 0, intercept_sd, log = TRUE)
  half <- prior("normal", list(0, intercept_sd), list(0, Inf))
  density <- log_density(half)
  route <- .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation"))
  expect_identical(route$type, "log_scale_product")
  for(value in c(-8, log(c(.05, .2, 1, 3)), 4)){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$behavior, "regular")
    expect_true(ordinate$exact)
    expect_identical(ordinate$method, "scale_mixture")
    expect_identical(ordinate$provenance$kind, "log_scale_product")
    expect_true(ordinate$provenance$integration_scale > 0)
    expect_equal(exp(ordinate$log_density), reference(value, f_half), tolerance = 1e-10)
    expect_equal(as.numeric(.prior_linear_density_height(density, value)),
                 reference(value, f_half), tolerance = 1e-10)
  }
  expect_true(all(prior_ordinate_status(density, c(-1, 0, 1))$eligible))
  # a value whose exponential is not representable has no ordinate value
  underflow <- prior_density_ordinate(density, -800)
  expect_identical(underflow$behavior, "regular")
  expect_false(underflow$exact)
  expect_match(underflow$reason, "not representable", fixed = TRUE)
  # neither has a value whose exponential is subnormal (below
  # .Machine$double.xmin, z < -708.4): for b0 ~ gamma(2, 4), whose density is
  # 16 x e^(-4x), f_Z(z) = 16 e^(2z) E[e^(-2 c b1)] (1 + O(e^z)), so
  # log f_Z(z) = log(16) + 2z + c^2 / 2 to double precision at z = -708,
  # while exp(-745) rounds to the smallest subnormal 4.9e-324 = e^(-744.44),
  # which shifted the log density by 0.56 with exact = TRUE
  subnormal_density <- log_density(prior("gamma", list(2, 4)))
  normal_edge <- prior_density_ordinate(subnormal_density, -708)
  expect_true(normal_edge$exact)
  expect_equal(normal_edge$log_density, log(16) - 2 * 708 + slope^2 / 2, tolerance = 1e-12)
  for(value in c(-740, -745)){
    subnormal <- prior_density_ordinate(subnormal_density, value)
    expect_identical(subnormal$behavior, "regular")
    expect_false(subnormal$exact)
    expect_match(subnormal$reason, "not representable", fixed = TRUE)
  }
  expect_true(is.na(.prior_density_route_density(
    .prior_density_route_from_adaptive(attr(subnormal_density, "adaptive_evaluation")), -745
  )))

  probability <- function(hypothesis){
    .hypothesis_prior_density_prob(
      density, hypothesis_parse(hypothesis)$statements[[1L]]$left, "theta"
    )
  }
  interval_reference <- function(lower, upper){
    stats::integrate(function(t){
      lower_tail <- stats::pnorm(exp(lower - slope * t), 0, intercept_sd, lower.tail = FALSE)
      upper_tail <- if(is.finite(upper)){
        stats::pnorm(exp(upper - slope * t), 0, intercept_sd, lower.tail = FALSE)
      }else{
        0
      }
      2 * (lower_tail - upper_tail) * stats::dnorm(t, 0, .5)
    }, -Inf, Inf, rel.tol = 1e-13)$value
  }
  expect_equal(probability("theta > -0.7"), interval_reference(-.7, Inf), tolerance = 1e-10)
  expect_equal(probability("theta > -1 & theta < 0.5"), interval_reference(-1, .5), tolerance = 1e-10)
  expect_equal(probability("theta < -1"), 1 - interval_reference(-1, Inf), tolerance = 1e-10)

  # the plotted density is the log image of the product's
  plotted <- .prior_linear_density_to_plot_data(density)$density
  checked <- plotted$x[c(10, 60, 150)]
  expect_equal(plotted$y[c(10, 60, 150)], vapply(checked, function(value){
    exp(prior_density_ordinate(density, value)$log_density)
  }, numeric(1)), tolerance = 1e-8)

  # a mixture prior of b0 is the mixture of the components' log images
  gamma <- prior("gamma", list(2, 4))
  mixture <- log_density(prior_mixture(list(half, gamma), is_null = c(FALSE, FALSE)))
  mixture_route <- .prior_density_route_from_adaptive(attr(mixture, "adaptive_evaluation"))
  expect_identical(vapply(mixture_route$components, `[[`, character(1), "type"),
                   c("log_scale_product", "log_scale_product"))
  for(value in c(-1, .5)){
    ordinate <- prior_density_ordinate(mixture, value)
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density),
                 .5 * reference(value, f_half) +
                   .5 * reference(value, function(x) stats::dgamma(x, 2, 4, log = TRUE)),
                 tolerance = 1e-10)
  }

  # not applicable, unchanged: a non-Gaussian other term, and a log-source
  # term with another weight, remain general convolutions
  for(unchanged in list(log_density(half, slope_prior = prior("t", list(0, .5, 3))),
                        log_density(half, weight = 2))){
    ordinate <- prior_density_ordinate(unchanged, -.5)
    expect_identical(ordinate$behavior, "unknown")
    expect_false(ordinate$exact)
    expect_identical(ordinate$reason, "General numerical convolutions are not structurally classified.")
    expect_identical(
      .prior_density_route_from_adaptive(attr(unchanged, "adaptive_evaluation"))$type,
      "unknown"
    )
  }
})

test_that("scale-product ordinates at subnormal distances from the offset have no value", {

  # Y = exp(log(b0) + .4 b1 + .3) = b0 W with b0 ~ gamma(2, 4) and the
  # lognormal W = exp(.3 + .4 b1), b1 ~ N(0, 1): f_Y(y) = 16 y E[W^-2]
  # (1 + O(y)), E[W^-2] = exp(-.6 + .32), so log f_Y(y) = log(16) + log(y) -
  # .28 to double precision for y <= 1e-300 (analytic reference). Below
  # .Machine$double.xmin the quadrature's factor argument y / W is rounded to
  # a multiple of the smallest subnormal: the log density was off by 1.35e-4
  # at 1e-320 and by 7.2e-2 at 4.9e-324 with exact = TRUE.
  density <- .prior_linear_combination_density(
    list(b0 = prior("gamma", list(2, 4)), b1 = prior("normal", list(0, 1)),
         k = prior("point", list(.3))),
    c(b0 = 1, b1 = .4, k = 1), source_transforms = c(b0 = "log", b1 = NA, k = NA),
    output_transformation = "exp"
  )
  route <- .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation"))
  expect_identical(route$type, "scale_product")
  small_log_density <- function(y) log(16) + log(y) - .28
  for(value in c(1e-300, 2.3e-308, .Machine$double.xmin)){
    ordinate <- prior_density_ordinate(density, value)
    expect_true(ordinate$exact)
    expect_lte(abs(ordinate$log_density - small_log_density(value)), 1e-12)
  }
  for(value in c(1e-310, 1e-320, 4.9e-324)){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$behavior, "regular")
    expect_false(ordinate$exact)
    expect_true(is.na(ordinate$log_density))
    expect_match(ordinate$reason, "not representable at full precision", fixed = TRUE)
  }
  status <- prior_ordinate_status(density, c(1e-300, 1e-320))
  expect_identical(status$eligible, c(TRUE, FALSE))
  expect_identical(status$condition[[2L]], "BayesTools_inexact_ordinate")
  # the plotted density takes the ordinate there (no value)
  plotted <- .prior_density_route_density(route, c(1e-300, 1e-320))
  expect_equal(plotted[[1L]], exp(small_log_density(1e-300)), tolerance = 1e-10)
  expect_true(is.na(plotted[[2L]]))
  # a mixture with such a component has no value either
  mixture <- .prior_linear_combination_density(
    list(b0 = prior_mixture(list(prior("gamma", list(2, 4)), prior("normal", list(0, .5), list(0, Inf))),
                            is_null = c(FALSE, FALSE)),
         b1 = prior("normal", list(0, 1)), k = prior("point", list(.3))),
    c(b0 = 1, b1 = .4, k = 1), source_transforms = c(b0 = "log", b1 = NA, k = NA),
    output_transformation = "exp"
  )
  expect_true(prior_density_ordinate(mixture, 1e-300)$exact)
  expect_false(prior_density_ordinate(mixture, 1e-320)$exact)

  # the rule is on the leaf's own distance (value - c) / w, whichever route
  # reaches it. Leaf: b ~ gamma(2, 4) times s ~ Beta(3, 2), with
  # f(x) = 16 x E[s^-2] (1 + O(x)), E[s^-2] = B(1, 2) / B(3, 2) = 6.
  leaf <- list(type = "scale_product", n_grid = 1024L, spec = .prior_scale_product_spec(
    offset = 0, scale = 1, factor = prior("gamma", list(2, 4)),
    multiplier = prior("beta", list(3, 2)), sources = list()
  ))
  leaf_log_density <- function(x) log(96) + log(x)
  # a lin transformation with slope 1e10 maps 1e-290 to 1e-300 and 1e-300 to
  # the subnormal 1e-310
  transformed <- .prior_density_route_transform(
    leaf, "lin", list(a = 0, b = 1e10), function() .prior_scale_product_hull(leaf$spec)
  )
  ordinate <- .prior_density_route_ordinate(transformed, 1e-290)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (leaf_log_density(1e-300) - log(1e10))), 1e-12)
  ordinate <- .prior_density_route_ordinate(transformed, 1e-300)
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "not representable at full precision", fixed = TRUE)
  # a product weight: the standardized distance 1e-300 / 1e10 is subnormal
  # (the log density was off by 8.4e-7 with exact = TRUE); at 1e-290 the
  # distance is normal but the density 9.6e-309 is subnormal, at 1e-280 both
  # are normal
  weighted <- leaf
  weighted$spec$scale <- 1e10
  ordinate <- .prior_density_route_ordinate(weighted, 1e-280)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (leaf_log_density(1e-290) - log(1e10))), 1e-12)
  expect_false(.prior_density_route_ordinate(weighted, 1e-290)$exact)
  expect_false(.prior_density_route_ordinate(weighted, 1e-300)$exact)
  # a subnormal value is refused also where its standardized distance is not
  expect_false(.prior_density_route_ordinate(leaf, 1e-310)$exact)
  small_weight <- leaf
  small_weight$spec$scale <- 1e-10
  expect_false(.prior_density_route_ordinate(small_weight, 1e-310)$exact)
  # the square-root share map of allocated SDs: gamma(2, 2) scale prior times
  # sqrt(2 s), s ~ Beta(3, 2), f(x) = 4 x E[1 / (2 s)] (1 + O(x)) with
  # E[1 / s] = B(2, 2) / B(3, 2) (the log density was off by 1.2e-4 at 1e-320)
  allocated <- list(type = "scale_product", n_grid = 1024L, spec = .prior_scale_product_spec(
    offset = 0, scale = 1, factor = prior("gamma", list(2, 2)),
    multiplier = prior("beta", list(3, 2)), sources = list(),
    map = list(type = "sqrt", scale = 2)
  ))
  ordinate <- .prior_density_route_ordinate(allocated, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density -
                   (log(4) + log(1e-300) + lbeta(2, 2) - lbeta(3, 2) - log(2))), 1e-12)
  expect_false(.prior_density_route_ordinate(allocated, 1e-320)$exact)
})

# Full-precision rule on every ordinate route: no ordinate is exact when its
# value is computed from a subnormal, zero (underflowed) or non-finite
# intermediate. References are analytic log densities evaluated without such
# intermediates (logs of the exact inputs), compared at 1e-12 in log, well
# above the double rounding of values near -1000.
refused_at_full_precision <- function(ordinate){
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "not representable at full precision", fixed = TRUE)
}

test_that("primitive and scalar affine ordinates of subnormal arguments have no value", {

  # gamma(3, 0.7): log f(x) = 3 log(0.7) - lgamma(3) + 2 log(x) - 0.7 x.
  # dgamma() rescales x, which rounds a subnormal: the log density was off by
  # 2.8e-4 at 1e-320 and 4.0e-2 at 3.5e-323 with exact = TRUE
  gamma_prior <- prior("gamma", list(3, .7))
  gamma_log_density <- function(x) 3 * log(.7) - lgamma(3) + 2 * log(x) - .7 * x
  ordinate <- prior_density_ordinate(gamma_prior, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - gamma_log_density(1e-300)), 1e-12)
  for(value in c(1e-310, 1e-320, 3.5e-323)){
    ordinate <- prior_density_ordinate(gamma_prior, value)
    expect_identical(ordinate$behavior, "regular")
    expect_true(is.na(ordinate$log_density))
    refused_at_full_precision(ordinate)
  }
  # the rule is on the value, for every family; 0 keeps its structural class
  refused_at_full_precision(prior_density_ordinate(prior("normal", list(0, 3)), 1e-320))
  expect_identical(prior_density_ordinate(gamma_prior, 0)$behavior, "zero")
  expect_true(prior_density_ordinate(gamma_prior, 0)$exact)

  # a gamma(2, 4) term with weight 3: f(y) = f_b(y / 3) / 3, log f(y) =
  # log(16) + log(y) - 2 log(3) - 4 y / 3; the rounded subnormal y / 3 was off
  # by 4.9e-4 at 1e-320 and by 0.15 at 3.5e-323 with exact = TRUE
  weighted <- .prior_linear_combination_density(list(b = prior("gamma", list(2, 4))), c(b = 3))
  weighted_log_density <- function(y) log(16) + log(y) - 2 * log(3) - 4 * y / 3
  ordinate <- prior_density_ordinate(weighted, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - weighted_log_density(1e-300)), 1e-12)
  for(value in c(1e-320, 3.5e-323)){
    ordinate <- prior_density_ordinate(weighted, value)
    refused_at_full_precision(ordinate)
    expect_match(ordinate$reason, "inverse affine value", fixed = TRUE)
  }
  expect_identical(prior_ordinate_status(weighted, c(1e-300, 1e-320, 3.5e-323))$eligible,
                   c(TRUE, FALSE, FALSE))
  # weight 1e30: the inverse value of 1e-300 underflows to 0, which was
  # classified as the gamma bound (a structural zero with exact = TRUE)
  large <- .prior_linear_combination_density(list(b = prior("gamma", list(2, 4))), c(b = 1e30))
  ordinate <- prior_density_ordinate(large, 1e-270)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density -
                   (log(16) + log(1e-270) - 2 * log(1e30))), 1e-12)
  ordinate <- prior_density_ordinate(large, 1e-300)
  expect_false(identical(ordinate$behavior, "zero"))
  refused_at_full_precision(ordinate)
})

test_that("log-source and Jacobian ordinates from subnormal or non-finite intermediates have no value", {

  log_source <- function(source_prior){
    .prior_linear_combination_density(list(b = source_prior), c(b = 1),
                                      source_transforms = c(b = "log"))
  }
  # Z = log(X), X ~ gamma(2, 4): log f_Z(z) = log(16) + 2 z - 4 e^z. e^z is
  # subnormal below z = -708.4 (off by 2.6e-3 at -740) and 0 below -745,
  # where the bound was classified (a structural zero, and an infinite
  # density for gamma(0.5, 1)) with exact = TRUE
  gamma_log <- log_source(prior("gamma", list(2, 4)))
  ordinate <- prior_density_ordinate(gamma_log, -700)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(16) - 1400 - 4 * exp(-700))), 1e-9)
  for(value in c(-709, -740, -800)){
    ordinate <- prior_density_ordinate(gamma_log, value)
    expect_identical(ordinate$behavior, "regular")
    refused_at_full_precision(ordinate)
  }
  ordinate <- prior_density_ordinate(log_source(prior("gamma", list(.5, 1))), -800)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)

  # X half-Cauchy: log f_Z(z) = log(2 / pi) + z - log1p(e^(2 z)). The t
  # density (evaluated before its log) underflows at x = e^700, so its log
  # density is -Inf; the Jacobian e^z made the ordinate at z = 700
  # exp(-700.45), reported as -Inf with exact = TRUE
  cauchy_log <- log_source(prior("t", list(0, 1, 1), list(0, Inf)))
  ordinate <- prior_density_ordinate(cauchy_log, 300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(2 / pi) - 300)), 1e-12)
  ordinate <- prior_density_ordinate(cauchy_log, 700)
  expect_identical(ordinate$behavior, "regular")
  expect_true(is.na(ordinate$log_density))
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "Jacobian", fixed = TRUE)
  # X half-t(0.5): e^720 overflows; log f_Z(720) = -361.1 was reported -Inf
  t_log <- log_source(prior("t", list(0, 1, .5), list(0, Inf)))
  refused_at_full_precision(prior_density_ordinate(t_log, 720))
})

test_that("named-transformation ordinates from subnormal intermediates have no value", {

  # transformations of the scalar route of b ~ gamma(2, 4), log f_b(x) =
  # log(16) + log(x) - 4 x
  source <- .prior_density_route_linear(list(b = prior("gamma", list(2, 4))), c(b = 1), NULL)
  hull <- function() c(0, Inf)
  # lin with slope 1e30: the inverse value of 1e-300 underflows to 0 (it was
  # classified as the bound, a structural zero with exact = TRUE)
  lin <- .prior_density_route_transform(source, "lin", list(a = 0, b = 1e30), hull)
  ordinate <- .prior_density_route_ordinate(lin, 1e-270)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(16) + log(1e-270) - 2 * log(1e30))), 1e-12)
  ordinate <- .prior_density_route_ordinate(lin, 1e-300)
  expect_false(identical(ordinate$behavior, "zero"))
  refused_at_full_precision(ordinate)
  # the smallest subnormal is not the endpoint 0 (it was classified as the
  # bound, a structural zero with exact = TRUE)
  identity <- .prior_density_route_transform(source, "lin", list(a = 0, b = 1), hull)
  ordinate <- .prior_density_route_ordinate(identity, 4.9e-324)
  expect_false(identical(ordinate$behavior, "zero"))
  refused_at_full_precision(ordinate)
  # exp_lin with b = 1/2 (y = sqrt(x), x = y^2): log f(y) = log f_b(y^2) +
  # log(2 y); at y = 1e-160, x = exp(2 log(y)) = 1e-320 was rounded (off by
  # 1.7e-5 with exact = TRUE)
  root <- .prior_density_route_transform(source, "exp_lin", list(a = 0, b = .5), hull)
  ordinate <- .prior_density_route_ordinate(root, 1e-150)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(16) + 3 * log(1e-150) - 4e-300 + log(2))), 1e-12)
  refused_at_full_precision(.prior_density_route_ordinate(root, 1e-160))
  # tanh: atanh(y) = y for tiny y (Jacobian 1)
  tanh <- .prior_density_route_transform(source, "tanh", list(), function() c(0, 1))
  ordinate <- .prior_density_route_ordinate(tanh, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(16) + log(1e-300))), 1e-12)
  refused_at_full_precision(.prior_density_route_ordinate(tanh, 1e-320))
  # exp keeps full precision (log(y) of a positive double): the exp of b ~
  # N(0, 1) at y = 1e-320 has log f(y) = log phi(log(y)) - log(y)
  normal <- .prior_density_route_linear(list(b = prior("normal", list(0, 1))), c(b = 1), NULL)
  exp_normal <- .prior_density_route_transform(normal, "exp", list(), function() c(0, Inf))
  ordinate <- .prior_density_route_ordinate(exp_normal, 1e-320)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density -
                   (stats::dnorm(log(1e-320), log = TRUE) - log(1e-320))), 1e-9)
})

test_that("quadrature ordinates at subnormal distances from their offset have no value", {

  # a pure scale mixture b s, b ~ N(0, 1), s ~ Beta(2, 2), is continuous at
  # its offset 0 with f(0) = phi(0) E[1 / s] = 3 phi(0); a Gaussian
  # convolution (additive SD 0.5) keeps full precision at subnormal values
  scale_mixture <- list(type = "conditional_normal", n_grid = 1024L, spec = list(
    additive_mean = 0, additive_sd = 0, product_mean = 0, product_sd = 1,
    multiplier = prior("beta", list(2, 2)), bounds = c(0, 1), sources = list()
  ))
  ordinate <- .prior_density_route_ordinate(scale_mixture, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - log(3 * stats::dnorm(0))), 1e-12)
  ordinate <- .prior_density_route_ordinate(scale_mixture, 1e-320)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)
  plotted <- .prior_density_route_density(scale_mixture, c(1e-300, 1e-320))
  # (the batched quadrature of plotted values accepts at a relative 1e-8)
  expect_equal(plotted[[1L]], 3 * stats::dnorm(0), tolerance = 1e-8)
  expect_true(is.na(plotted[[2L]]))
  convolution_normal <- scale_mixture
  convolution_normal$spec$additive_sd <- .5
  continuous_at_zero <- stats::integrate(function(s) stats::dnorm(0, 0, sqrt(.25 + s^2)) * 6 * s * (1 - s),
                                         0, 1, rel.tol = 1e-13)$value
  ordinate <- .prior_density_route_ordinate(convolution_normal, 1e-320)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - log(continuous_at_zero)), 1e-10)

  # two-term convolutions: A - B with A, B ~ Beta(0.5, 1) is infinite only at
  # its inner meeting point 0, f(x) = log(1 + sqrt(1 - x)) / 2 - log(x) / 4;
  # the smallest subnormal counted as that point (infinite with exact = TRUE)
  difference <- .prior_linear_combination_density(
    list(a = prior("beta", list(.5, 1)), b = prior("beta", list(.5, 1))), c(a = 1, b = -1)
  )
  expect_identical(prior_density_ordinate(difference, 0)$behavior, "infinite")
  ordinate <- prior_density_ordinate(difference, 4.9e-324)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)
  # gamma(2, 4) + t3 is continuous at 0 (reference: integrate() of the
  # convolution at 0, rel.tol 1e-13)
  sum_density <- .prior_linear_combination_density(
    list(a = prior("gamma", list(2, 4)), b = prior("t", list(0, 1, 3))), c(a = 1, b = 1)
  )
  at_zero <- stats::integrate(function(a) stats::dgamma(a, 2, 4) * stats::dt(-a, 3), 0, Inf,
                              rel.tol = 1e-13)$value
  ordinate <- prior_density_ordinate(sum_density, 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - log(at_zero)), 1e-10)
  refused_at_full_precision(prior_density_ordinate(sum_density, 1e-320))
})

test_that("primitive ordinates from subnormal internal intermediates have no value", {

  # dgamma() evaluates at value * rate and dlnorm() at log(value * sdlog):
  # gamma(2, 1e-20) at 1e-300 has the subnormal argument 1e-320 (its log
  # density 2 log(1e-20) + log(1e-300) - 1e-320 was off by 1.1e-5 with
  # exact = TRUE), gamma(2, 1e-10) at 1e-290 the normal argument 1e-300
  refused_at_full_precision(prior_density_ordinate(prior("gamma", list(2, 1e-20)), 1e-300))
  ordinate <- prior_density_ordinate(prior("gamma", list(2, 1e-10)), 1e-290)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (2 * log(1e-10) + log(1e-290) - 1e-300)), 1e-12)
  refused_at_full_precision(prior_density_ordinate(prior("lognormal", list(log(1e-300), 1e-20)), 1e-300))
  expect_true(prior_density_ordinate(prior("gamma", list(1, 1e-20)), 0)$exact)
  # with shape 1 the argument does not enter the log (exponential density
  # rate e^(-rate x), log(1e-10) - 1e-310 at 1e-300)
  ordinate <- prior_density_ordinate(prior("gamma", list(1, 1e-10)), 1e-300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - log(1e-10)), 1e-12)

  # extraDistr::dlst() evaluates the t density before its log, which is
  # rounded where the density is subnormal: the Cauchy log density
  # -log(pi) - 2 log(x) (up to x^-2) was off by 4.1e-4 at 1e160 with
  # exact = TRUE; an underflowed density keeps -Inf
  cauchy <- prior("t", list(0, 1, 1))
  ordinate <- prior_density_ordinate(cauchy, 1e150)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (-log(pi) - 2 * log(1e150))), 1e-12)
  ordinate <- prior_density_ordinate(cauchy, 1e160)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)
  expect_identical(prior_density_ordinate(cauchy, 1e200)$log_density, -Inf)
  # a Jacobian factor above 1 carried that error to ordinary ordinates: the
  # log image of a half-Cauchy term, log f_Z(z) = log(2 / pi) - z up to
  # e^(-2 z), was off by 4.95e-4 at z = 370, and the term with weight 1e-100
  # by 4.1e-4 at 1e60
  cauchy_log <- .prior_linear_combination_density(
    list(b = prior("t", list(0, 1, 1), list(0, Inf))), c(b = 1), source_transforms = c(b = "log")
  )
  ordinate <- prior_density_ordinate(cauchy_log, 340)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(2 / pi) - 340)), 1e-12)
  refused_at_full_precision(prior_density_ordinate(cauchy_log, 370))
  expect_false(prior_ordinate_status(cauchy_log, 370)$eligible)
  weighted <- .prior_linear_combination_density(
    list(b = prior("t", list(0, 1, 1), list(0, Inf))), c(b = 1e-100)
  )
  ordinate <- prior_density_ordinate(weighted, 1e50)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(2 / pi) - 2 * log(1e150) + log(1e100))), 1e-12)
  refused_at_full_precision(prior_density_ordinate(weighted, 1e60))
})

test_that("quadrature ordinates with a subnormal value have no value", {

  # the log image of X W, X half-Cauchy and W = exp(.3 + .4 b1) lognormal:
  # f_Z(z) = E[sech(z - G)] / pi with G ~ N(.3, .4), so log f_Z(z) =
  # log(2 / pi) - z + .3 + .4^2 / 2 up to e^(-2 z) (analytic reference). The
  # product's density at e^370 is about 4e-322, and the log of that
  # subnormal quadrature value was off by 1.2e-2 with exact = TRUE
  log_image <- .prior_linear_combination_density(
    list(b0 = prior("t", list(0, 1, 1), list(0, Inf)), b1 = prior("normal", list(0, 1)),
         k = prior("point", list(.3))),
    c(b0 = 1, b1 = .4, k = 1), source_transforms = c(b0 = "log", b1 = NA, k = NA)
  )
  route <- .prior_density_route_from_adaptive(attr(log_image, "adaptive_evaluation"))
  expect_identical(route$type, "log_scale_product")
  # (the ordinates' quadrature acceptance criterion is a relative 1e-4)
  ordinate <- prior_density_ordinate(log_image, 300)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - (log(2 / pi) - 300 + .38)), 1e-4)
  ordinate <- prior_density_ordinate(log_image, 370)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)

  # a pure scale mixture b s, b ~ N(0, 1), s ~ Beta(2, 2), in its far tail
  # (reference: integrate() of the Gaussian factor scaled by e^(x^2 / 2),
  # rel.tol 1e-13): the density is 2.3e-320 at 38 (its log was off by 3e-5)
  scale_mixture <- list(type = "conditional_normal", n_grid = 1024L, spec = list(
    additive_mean = 0, additive_sd = 0, product_mean = 0, product_sd = 1,
    multiplier = prior("beta", list(2, 2)), bounds = c(0, 1), sources = list()
  ))
  tail_log_density <- function(x){
    log(stats::integrate(function(s) exp(stats::dnorm(x / s, log = TRUE) + x^2 / 2) * 6 * (1 - s),
                         0, 1, rel.tol = 1e-13)$value) - x^2 / 2
  }
  ordinate <- .prior_density_route_ordinate(scale_mixture, 37)
  expect_true(ordinate$exact)
  expect_lte(abs(ordinate$log_density - tail_log_density(37)), 1e-4)
  ordinate <- .prior_density_route_ordinate(scale_mixture, 38)
  expect_identical(ordinate$behavior, "regular")
  refused_at_full_precision(ordinate)
  # plotted densities keep the value as a display estimate
  expect_lte(abs(log(ordinate$provenance$integration$estimate) - tail_log_density(38)), 1e-3)
})

test_that("products of normal terms are classified on the scale-mixture route", {

  # beta * sigma with beta, sigma ~ N(0, 1) has the density K0(|x|) / pi
  # (the product of two standard normals), which is infinite at 0 because
  # sigma's density is positive there.
  product_priors <- list(
    beta = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1))
  )
  attr(product_priors$beta, "multiply_by") <- "sigma"
  product <- BayesTools:::.prior_linear_combination_density(
    prior_list = product_priors,
    weights = c(beta = 1),
    n_grid = 128
  )
  regular <- prior_density_ordinate(product, 1)
  expect_identical(regular$behavior, "regular")
  expect_identical(regular$method, "conditional_normal_mixture")
  expect_true(regular$exact)
  expect_equal(exp(regular$log_density), besselK(1, 0) / pi, tolerance = 1e-8)
  singular <- prior_density_ordinate(product, 0)
  expect_identical(singular$behavior, "infinite")
  expect_true(singular$exact)
  expect_identical(singular$method, "conditional_normal_mixture")
  expect_identical(singular$provenance$kind, "product_singularity")

  transformed_product <- BayesTools:::.prior_linear_combination_density(
    prior_list = product_priors,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp"
  )
  transformed_singular <- prior_density_ordinate(transformed_product, 1)
  expect_identical(transformed_singular$behavior, "infinite")
  expect_identical(transformed_singular$method, "named_transform")
})

test_that("invalid named transformation provenance is undefined", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("normal", list(0, 1))),
    weights = c(beta = 1),
    n_grid = 128
  )
  adaptive <- attr(density, "adaptive_evaluation", exact = TRUE)
  adaptive$arguments$output_transformation <- "exp_lin"
  adaptive$arguments$output_transformation_arguments <- list(a = 0, b = 1)
  attr(density, "adaptive_evaluation") <- adaptive

  out <- prior_density_ordinate(density, 1)
  expect_identical(out$behavior, "undefined")
  expect_true(out$exact)
  expect_identical(out$log_density, NA_real_)

  invalid_arguments <- adaptive
  invalid_arguments$arguments$output_transformation <- "lin"
  invalid_arguments$arguments$output_transformation_arguments <-
    list(a = Inf, b = 1)
  attr(density, "adaptive_evaluation") <- invalid_arguments
  invalid_out <- prior_density_ordinate(density, 1)
  expect_identical(invalid_out$behavior, "undefined")

  atom <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("point", list(0))),
    weights = c(x = 1),
    n_grid = 128
  )
  atom_adaptive <- attr(atom, "adaptive_evaluation", exact = TRUE)
  atom_adaptive$arguments$output_transformation <- "exp_lin"
  atom_adaptive$arguments$output_transformation_arguments <-
    list(a = 0, b = 1)
  attr(atom, "adaptive_evaluation") <- atom_adaptive
  atom_out <- prior_density_ordinate(atom, 1)
  expect_identical(atom_out$behavior, "undefined")
})

test_that("ordinary density underflow does not imply a structural zero", {

  normal_prior <- prior("normal", list(0, 1))
  expect_identical(pdf(normal_prior, 40), 0)

  out <- prior_density_ordinate(normal_prior, 40)
  expect_identical(out$behavior, "regular")
  expect_true(out$exact)
  expect_equal(out$log_density, stats::dnorm(40, log = TRUE))

  extreme <- prior_density_ordinate(normal_prior, 1e200)
  expect_identical(extreme$behavior, "regular")
  expect_identical(extreme$log_density, -Inf)
  expect_match(extreme$reason, "structurally regular")

  log_source <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 1),
    source_transforms = c(x = "log"),
    n_grid = 128
  )
  # the log source at z = 1000 is regular, but e^1000 overflows, so its
  # density f_X(e^z) e^z is not evaluated: no value rather than an asserted
  # -Inf (which is wrong for heavy-tailed sources, see the full-precision
  # tests above)
  log_source_extreme <- prior_density_ordinate(log_source, 1000)
  expect_identical(log_source_extreme$behavior, "regular")
  expect_false(log_source_extreme$exact)
  expect_true(is.na(log_source_extreme$log_density))
  expect_match(log_source_extreme$reason, "not representable at full precision", fixed = TRUE)

  normal_sum <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(0, 1))
    ),
    weights = c(x = 1, y = 1),
    n_grid = 128
  )
  normal_sum_extreme <- prior_density_ordinate(normal_sum, 1e200)
  expect_identical(normal_sum_extreme$behavior, "regular")
  expect_identical(normal_sum_extreme$log_density, -Inf)
  expect_match(normal_sum_extreme$reason, "structurally regular")
})

test_that("legacy numerical metadata never establishes structural behavior", {

  stale <- structure(
    list(
      density = list(x = c(-1, 0, 1), y = c(1, 0, 1), mass = 1),
      points = data.frame(x = numeric(), p = numeric()),
      n_grid = 3L
    ),
    class = c("prior_linear_density", "prior_density")
  )
  attr(stale, "singular_density_points") <- 0
  attr(stale, "fft_clipping") <- list(clipped_value_count = 1L)
  attr(stale, "density_evaluator") <- function(value) rep(0, length(value))

  stale_zero <- prior_density_ordinate(stale, 0)
  expect_identical(stale_zero$behavior, "unknown")
  expect_false(stale_zero$exact)
  expect_identical(stale_zero$log_density, -Inf)

  stale_outside <- prior_density_ordinate(stale, 100)
  expect_identical(stale_outside$behavior, "unknown")

  structural <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("normal", list(0, 1))),
    weights = c(x = 1),
    n_grid = 128
  )
  structural$density$y[] <- 0
  expect_identical(
    prior_density_ordinate(structural, 0)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(structural, 100)$behavior,
    "regular"
  )
})

test_that("ordinate provenance is compact and contains no closures", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid = 128
  )
  out <- prior_density_ordinate(density, 0)

  contains_function <- function(x){
    if(is.function(x)) return(TRUE)
    if(!is.list(x)) return(FALSE)
    any(vapply(x, contains_function, logical(1)))
  }

  expect_false(contains_function(unclass(out)))
  expect_lt(length(serialize(unclass(out), NULL)), 20000)
  expect_false("adaptive_evaluation" %in% names(out$provenance))
  expect_false("density_evaluator" %in% names(out$provenance))

  context <- BayesTools:::.prior_density_context(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    column_names = c("x", "y"),
    n_grid = 128
  )
  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, 0, 2, 0, 3, 0),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("x", "y"))
    )
  )
  row_out <- prior_density_ordinate(row_density, 0)
  contains_matrix <- function(x){
    if(is.matrix(x)) return(TRUE)
    if(!is.list(x)) return(FALSE)
    any(vapply(x, contains_matrix, logical(1)))
  }
  forbidden_names <- function(x){
    if(!is.list(x)) return(character())
    c(
      intersect(
        names(x),
        c("grid", "density", "draws", "samples", "raw_weights",
          "representative_weights", "adaptive_evaluation",
          "density_evaluator")
      ),
      unlist(lapply(x, forbidden_names), use.names = FALSE)
    )
  }
  expect_false(contains_matrix(row_out$provenance))
  expect_length(forbidden_names(row_out$provenance), 0)
})

test_that("original-scale intercepts with mixture sources are classified per component", {

  # b0 - (m / s) b1 of a model-averaged meta-regression with a standardized
  # moderator (m = 1, s = 2): a mixture intercept (spike at .2 with weight 1,
  # N(0, s0) with weight 3) and a spike-and-slab slope (inclusion .3). The
  # measure is the product of the terms' components: an atom of .25 * .7 at
  # .2 and three normal components. References: the analytic normal mixture.
  s0 <- sqrt(.5)
  s1 <- sqrt(.125)
  scaled <- .5 * s1
  context <- .prior_density_context(
    prior_list = list(
      mu_intercept = prior_mixture(list(
        prior("point", list(.2), prior_weights = 1),
        prior("normal", list(0, s0), prior_weights = 3)
      ), is_null = c(TRUE, FALSE)),
      mu_x = prior_spike_and_slab(prior("normal", list(0, s1)),
                                  prior_inclusion = prior("spike", list(.3)))
    ),
    column_names = c("mu_intercept", "mu_x"),
    formula_scale = list(mu = formula_scale_for_test(~ x, list(x = list(mean = 1, sd = 2))))
  )
  density <- .prior_density_from_context(context, c(mu_intercept = 1))
  weights <- c(.25 * .3, .75 * .7, .75 * .3)
  means <- c(.2, 0, 0)
  sds <- c(scaled, s0, sqrt(s0^2 + scaled^2))
  continuous <- function(value) sum(weights * stats::dnorm(value, means, sds))
  above <- function(value) sum(weights * stats::pnorm(value, means, sds, lower.tail = FALSE))

  for(value in c(-.5, 0, .1, .3, 1)){
    ordinate <- prior_density_ordinate(density, value)
    expect_identical(ordinate$behavior, "regular", label = format(value))
    expect_identical(ordinate$method, "finite_mixture")
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density), continuous(value), tolerance = 1e-12)
    expect_equal(as.numeric(.prior_linear_density_height(density, value)),
                 continuous(value), tolerance = 1e-12)
  }
  # the spike-spike component is an atom; the continuous part stays regular
  atom <- prior_density_ordinate(density, .2)
  expect_identical(atom$behavior, "point_mass")
  expect_equal(atom$point_mass, .25 * .7, tolerance = 1e-15)
  expect_identical(atom$provenance$continuous_behavior, "regular")
  expect_equal(exp(atom$log_density), continuous(.2), tolerance = 1e-12)
  expect_error(.hypothesis_prior_density_height(density, .2),
               "There is a point mass in the prior at the exact null hypothesis value.",
               fixed = TRUE)
  support <- .posterior_support_from_prior_context_weights(context, c(mu_intercept = 1))
  expect_identical(support$points, .2)
  expect_identical(support$type, "mixed")

  probability <- function(hypothesis){
    side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
    .hypothesis_prior_density_prob(density, side, "theta")
  }
  expect_equal(probability("theta > 0.2"), above(.2), tolerance = 1e-14)
  expect_equal(probability("theta >= 0.2"), above(.2) + .25 * .7, tolerance = 1e-14)
  expect_equal(probability("theta < -0.3"), sum(weights * stats::pnorm(-.3, means, sds)),
               tolerance = 1e-14)
  expect_equal(probability("theta > 0 & theta < 0.5"),
               above(0) - above(.5) + .25 * .7, tolerance = 1e-14)

  # a Cauchy alternative gives a Gaussian-convolution component; reference by
  # integrating over the normal intercept (the package integrates over the
  # slope)
  context <- .prior_density_context(
    prior_list = list(
      mu_intercept = prior_mixture(list(
        prior("point", list(0), prior_weights = 1),
        prior("normal", list(0, s0), prior_weights = 1)
      ), is_null = c(TRUE, FALSE)),
      mu_x = prior_mixture(list(
        prior("point", list(0), prior_weights = 1),
        prior("cauchy", list(0, .5), prior_weights = 1)
      ), is_null = c(TRUE, FALSE))
    ),
    column_names = c("mu_intercept", "mu_x"),
    formula_scale = list(mu = formula_scale_for_test(~ x, list(x = list(mean = 1, sd = 2))))
  )
  density <- .prior_density_from_context(context, c(mu_intercept = 1))
  convolution <- function(value){
    pieces <- s0 * c(-40, -10, -3, -1, 0, 1, 3, 10, 40)
    sum(vapply(seq_len(length(pieces) - 1L), function(i){
      stats::integrate(function(g) stats::dnorm(g, 0, s0) * stats::dcauchy(value - g, 0, .25),
                       pieces[i], pieces[i + 1L], rel.tol = 1e-12)$value
    }, numeric(1)))
  }
  convolution_above <- function(value){
    pieces <- s0 * c(-40, -10, -3, -1, 0, 1, 3, 10, 40)
    sum(vapply(seq_len(length(pieces) - 1L), function(i){
      stats::integrate(function(g) stats::dnorm(g, 0, s0) *
                         stats::pcauchy(value - g, 0, .25, lower.tail = FALSE),
                       pieces[i], pieces[i + 1L], rel.tol = 1e-12)$value
    }, numeric(1)))
  }
  for(value in c(-1, .1, .4)){
    reference <- .25 * stats::dnorm(value, 0, s0) + .25 * stats::dcauchy(value, 0, .25) +
      .25 * convolution(value)
    ordinate <- prior_density_ordinate(density, value)
    expect_true(ordinate$exact)
    expect_lt(abs(exp(ordinate$log_density) / reference - 1), 1e-8)
  }
  atom <- prior_density_ordinate(density, 0)
  expect_identical(atom$behavior, "point_mass")
  expect_equal(atom$point_mass, .25, tolerance = 1e-15)
  expect_identical(atom$provenance$continuous_behavior, "regular")
  side <- hypothesis_parse("theta > 0.4")$statements[[1L]]$left
  expect_lt(abs(.hypothesis_prior_density_prob(density, side, "theta") -
                  (.25 * stats::pnorm(.4, 0, s0, lower.tail = FALSE) +
                     .25 * stats::pcauchy(.4, 0, .25, lower.tail = FALSE) +
                     .25 * convolution_above(.4))), 1e-10)

  # two Cauchy terms sum to a Cauchy term, so a Cauchy intercept makes every
  # component closed-form: 0, Cauchy(0, .5), -0.5 Cauchy(0, .5) and their sum
  context$prior_list$mu_intercept <- prior_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("cauchy", list(0, .5), prior_weights = 1)
  ), is_null = c(TRUE, FALSE))
  density <- .prior_density_from_context(context, c(mu_intercept = 1))
  for(value in c(-1, .1)){
    ordinate <- prior_density_ordinate(density, value)
    expect_true(ordinate$exact)
    expect_equal(exp(ordinate$log_density),
                 .25 * (stats::dcauchy(value, 0, .5) + stats::dcauchy(value, 0, .25) +
                          stats::dcauchy(value, 0, .75)),
                 tolerance = 1e-14)
  }

  # a leaf without a structural route (three t terms) leaves the mixture
  # unknown, as for any combination
  slab <- function() prior_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("t", list(0, .5, 3), prior_weights = 1)
  ), is_null = c(TRUE, FALSE))
  context <- .prior_density_context(
    prior_list = list(mu_intercept = slab(), mu_x = slab(), mu_z = prior("t", list(0, .5, 3))),
    column_names = c("mu_intercept", "mu_x", "mu_z"),
    formula_scale = list(mu = formula_scale_for_test(~ x + z, list(x = list(mean = 1, sd = 2), z = list(mean = 1, sd = 2))))
  )
  ordinate <- prior_density_ordinate(.prior_density_from_context(context, c(mu_intercept = 1)), .1)
  expect_identical(ordinate$behavior, "unknown")
})

test_that("point-mass ordinates record the behavior of their continuous part", {

  # provenance$continuous_behavior is present exactly for point-mass results
  point <- prior_density_ordinate(prior("point", list(location = .3)), .3)
  expect_identical(point$behavior, "point_mass")
  expect_identical(point$provenance$continuous_behavior, "zero")
  expect_identical(point$log_density, -Inf)
  regular <- prior_density_ordinate(prior("normal", list(0, 1)), .3)
  expect_null(regular$provenance$continuous_behavior)

  # an atom next to a continuous part without a structural route (the
  # component with three t terms)
  spike_t <- prior_spike_and_slab(prior("t", list(0, .5, 3)),
                                  prior_inclusion = prior("spike", list(.5)))
  density <- .prior_linear_combination_density(
    list(a = spike_t, b = spike_t, c = spike_t), c(a = 1, b = -.5, c = .25)
  )
  atom <- prior_density_ordinate(density, 0)
  expect_identical(atom$behavior, "point_mass")
  expect_equal(atom$point_mass, .125, tolerance = 1e-15)
  expect_identical(atom$provenance$continuous_behavior, "unknown")
  away <- prior_density_ordinate(density, .1)
  expect_identical(away$behavior, "unknown")
  expect_null(away$provenance$continuous_behavior)
})

test_that("a regular ordinate without a value is never reported as exact", {

  # exact = TRUE always comes with an available regular log density. The
  # exp_lin image y = x^-2 of an inverse-gamma(2, 1) source has a positive
  # finite limit at y = 0 (the source density ~ x^-3 at infinity over the
  # Jacobian 2 x^-3), which is not available structurally for this family.
  image <- .prior_linear_combination_density(
    list(x = prior("invgamma", list(2, 1))), c(x = 1),
    output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = -2)
  )
  limit <- prior_density_ordinate(image, 0)
  expect_identical(limit$behavior, "regular")
  expect_false(limit$exact)
  expect_true(is.na(limit$log_density))
  expect_identical(limit$reason,
                   "The regular prior-density ordinate has no structural value at the requested value.")

  # a pure scale mixture N(0, 1e-3) * Cauchy(0, 1) whose quadrature at .3 is
  # rejected by its diagnostics, alone and inside a finite mixture
  priors <- list(beta = prior("normal", list(0, 1e-3)), sigma = prior("cauchy", list(0, 1)))
  attr(priors$beta, "multiply_by") <- "sigma"
  rejected <- prior_density_ordinate(.prior_linear_combination_density(priors, c(beta = 1)), .3)
  expect_identical(rejected$behavior, "regular")
  expect_false(rejected$exact)
  expect_true(is.na(rejected$log_density))
  expect_false(rejected$provenance$integration$converged)
  expect_match(rejected$reason, "The prior-density quadrature was rejected by its diagnostics",
               fixed = TRUE)
  beta_mixture <- prior_mixture(list(prior("normal", list(0, 1e-3), prior_weights = 1),
                                     prior("normal", list(0, 1), prior_weights = 1)),
                                is_null = c(FALSE, FALSE))
  attr(beta_mixture, "multiply_by") <- "sigma"
  mixture <- prior_density_ordinate(.prior_linear_combination_density(
    list(beta = beta_mixture, sigma = priors$sigma), c(beta = 1)
  ), .3)
  expect_identical(mixture$behavior, "regular")
  expect_false(mixture$exact)
  expect_true(is.na(mixture$log_density))
})

test_that("square-root share scale products match 50-digit references", {

  # Y = T sqrt(k S), S ~ Beta(alpha, beta): the component SDs and partial-set
  # totals of variance allocations. References: mpmath (50 digits) of
  # f(y) = int_0^1 Beta(s; alpha, beta) f_T(y / sqrt(k s)) / sqrt(k s) ds and
  # P(l < Y < u) = int_0^1 Beta(s; alpha, beta) [F_T(u / sqrt(k s)) - F_T(l / sqrt(k s))] ds
  # for T half-normal(0.5) and gamma(2, 2); K = 4 components with a common
  # concentration a in {0.5, 1, 2}, a set of m in {1, 2, 3} components has
  # the share Beta(m a, (4 - m) a), k in {1, 4}, and y in {1e-3, 0.06, the
  # median, the 0.999 quantile} (6 digits). Two quadrature splits of every
  # reference agree to 7e-28. The
  # quadratures accept a relative error of 1e-4 by their QUADPACK estimate;
  # every reference lies within that estimate and within 1e-8 relative.
  densities <- utils::read.csv(text = "
family,alpha,beta,kappa,y,density
halfnormal,0.5,1.5,1.0,0.001,12.1211464118084787
halfnormal,0.5,1.5,1.0,0.06,3.8305507780385271834
halfnormal,0.5,1.5,1.0,0.106452,2.7130794107251911593
halfnormal,0.5,1.5,1.0,1.18029,0.0070365258775834877193
halfnormal,0.5,1.5,4.0,0.001,6.7647342683777474047
halfnormal,0.5,1.5,4.0,0.06,2.6094674471018446359
halfnormal,0.5,1.5,4.0,0.212904,1.3565397053625955796
halfnormal,0.5,1.5,4.0,2.36058,0.0035182629387917438597
halfnormal,1.0,1.0,1.0,0.001,3.1835446262858201551
halfnormal,1.0,1.0,1.0,0.06,2.7344897833290089887
halfnormal,1.0,1.0,1.0,0.202617,1.8291228411874105379
halfnormal,1.0,1.0,1.0,1.37694,0.0071017626401189283671
halfnormal,1.0,1.0,4.0,0.001,1.5937699194902250243
halfnormal,1.0,1.0,4.0,0.06,1.4786406446194131923
halfnormal,1.0,1.0,4.0,0.405234,0.91456142059370526896
halfnormal,1.0,1.0,4.0,2.75388,0.0035508813200594641835
halfnormal,1.5,0.5,1.0,0.001,2.0317660122608680356
halfnormal,1.5,0.5,1.0,0.06,1.9825587475815563224
halfnormal,1.5,0.5,1.0,0.276089,1.4784110380248839041
halfnormal,1.5,0.5,1.0,1.52301,0.0071139412564858591465
halfnormal,1.5,0.5,4.0,0.001,1.0158940306586918455
halfnormal,1.5,0.5,4.0,0.06,1.0084688564968874622
halfnormal,1.5,0.5,4.0,0.552178,0.73920551901244195203
halfnormal,1.5,0.5,4.0,3.04602,0.0035569706282429295732
halfnormal,1.0,3.0,1.0,0.001,5.0825121898012686797
halfnormal,1.0,3.0,1.0,0.06,3.8371120003385229215
halfnormal,1.0,3.0,1.0,0.132132,2.6937536314980462943
halfnormal,1.0,3.0,1.0,1.08341,0.0077809242556820838564
halfnormal,1.0,3.0,4.0,0.001,2.5472369736472509309
halfnormal,1.0,3.0,4.0,0.06,2.215366042529688288
halfnormal,1.0,3.0,4.0,0.264264,1.3468768157490231471
halfnormal,1.0,3.0,4.0,2.16682,0.0038904621278410419282
halfnormal,2.0,2.0,1.0,0.001,2.5532051261866777123
halfnormal,2.0,2.0,1.0,0.06,2.4741846163583865645
halfnormal,2.0,2.0,1.0,0.219144,1.847957469493536721
halfnormal,2.0,2.0,1.0,1.3152,0.0075097290545527150945
halfnormal,2.0,2.0,4.0,0.001,1.2766121097439485041
halfnormal,2.0,2.0,4.0,0.06,1.2659593537335072535
halfnormal,2.0,2.0,4.0,0.438287,0.92397989855525983118
halfnormal,2.0,2.0,4.0,2.63039,0.0037549971960245362402
halfnormal,3.0,1.0,1.0,0.001,1.9149165628695140734
halfnormal,3.0,1.0,1.0,0.06,1.8921727332901174801
halfnormal,3.0,1.0,1.0,0.283861,1.484332459165976421
halfnormal,3.0,1.0,1.0,1.49361,0.0072920278190469234122
halfnormal,3.0,1.0,4.0,0.001,0.95746067507947563771
halfnormal,3.0,1.0,4.0,0.06,0.95459653759504403555
halfnormal,3.0,1.0,4.0,0.567722,0.74216622958298821052
halfnormal,3.0,1.0,4.0,2.98723,0.0036458899627828528681
halfnormal,2.0,6.0,1.0,0.001,3.8089167517310002891
halfnormal,2.0,6.0,1.0,0.06,3.5364384525362891962
halfnormal,2.0,6.0,1.0,0.148864,2.6675656704451543217
halfnormal,2.0,6.0,1.0,0.99669,0.0088525314931313670204
halfnormal,2.0,6.0,4.0,0.001,1.9044953182748367045
halfnormal,2.0,6.0,4.0,0.06,1.8655699924396668258
halfnormal,2.0,6.0,4.0,0.297728,1.3337828352225771609
halfnormal,2.0,6.0,4.0,1.99338,0.0044262657465656835102
halfnormal,4.0,4.0,1.0,0.001,2.3806222077564183538
halfnormal,4.0,4.0,1.0,0.06,2.3366464250896933922
halfnormal,4.0,4.0,1.0,0.228707,1.834871124939731248
halfnormal,4.0,4.0,1.0,1.26097,0.0080717576658517429955
halfnormal,4.0,4.0,4.0,0.001,1.1903157460943770602
halfnormal,4.0,4.0,4.0,0.06,1.1847648997235440524
halfnormal,4.0,4.0,4.0,0.457415,0.91743455574186607092
halfnormal,4.0,4.0,4.0,2.52193,0.0040360326218048408291
halfnormal,6.0,2.0,1.0,0.001,1.8747443213358935575
halfnormal,6.0,2.0,1.0,0.06,1.8553621650590443726
halfnormal,6.0,2.0,1.0,0.288032,1.4786771888723587968
halfnormal,6.0,2.0,1.0,1.46805,0.0075286982455601021565
halfnormal,6.0,2.0,4.0,0.001,0.93737419164283940461
halfnormal,6.0,2.0,4.0,0.06,0.93494113711967618907
halfnormal,6.0,2.0,4.0,0.576065,0.73933799084857796352
halfnormal,6.0,2.0,4.0,2.93611,0.0037642168147291586259
gamma,1.0,3.0,1.0,0.001,0.11742473312962439688
gamma,1.0,3.0,1.0,0.06,1.5386369814599611014
gamma,1.0,3.0,1.0,0.340895,1.2310692321110704506
gamma,1.0,3.0,1.0,2.88881,0.0023406901542596959513
gamma,2.0,2.0,1.0,0.001,0.011936637317410546226
gamma,2.0,2.0,1.0,0.06,0.54388844791819300311
gamma,2.0,2.0,1.0,0.552177,0.88208849156700690014
gamma,2.0,2.0,1.0,3.59397,0.002095261673286114544
gamma,3.0,1.0,1.0,0.001,0.0059840239681235380112
gamma,3.0,1.0,1.0,0.06,0.30721469363055722183
gamma,3.0,1.0,1.0,0.709247,0.72300674910655449583
gamma,3.0,1.0,1.0,4.14471,0.0019298905897387106128
")
  regions <- utils::read.csv(text = "
family,alpha,beta,kappa,lower,upper,probability
halfnormal,0.5,1.5,1.0,0.06,0.106452,0.14934213771794214658
halfnormal,0.5,1.5,4.0,0.06,0.212904,0.2826389574003450586
halfnormal,1.0,1.0,1.0,0.06,0.202617,0.32244863063755080775
halfnormal,1.0,1.0,4.0,0.06,0.405234,0.40779659153179269564
halfnormal,1.5,0.5,1.0,0.06,0.276089,0.37917479265006984894
halfnormal,1.5,0.5,4.0,0.06,0.552178,0.43920661247627180811
halfnormal,1.0,3.0,1.0,0.06,0.132132,0.23333470715250768631
halfnormal,1.0,3.0,4.0,0.06,0.264264,0.35715871351605066342
halfnormal,2.0,2.0,1.0,0.06,0.219144,0.34844956645837690886
halfnormal,2.0,2.0,4.0,0.06,0.438287,0.42362027010833819594
halfnormal,3.0,1.0,1.0,0.06,0.283861,0.3855616860044431745
halfnormal,3.0,1.0,4.0,0.06,0.567722,0.44260993469434792294
halfnormal,2.0,6.0,1.0,0.06,0.148864,0.27729422178891449124
halfnormal,2.0,6.0,4.0,0.06,0.297728,0.38653508263892599338
halfnormal,4.0,4.0,1.0,0.06,0.228707,0.35804536027558892689
halfnormal,4.0,4.0,4.0,0.06,0.457415,0.42869216320545956428
halfnormal,6.0,2.0,1.0,0.06,0.288032,0.38790290524704623139
halfnormal,6.0,2.0,4.0,0.06,0.576065,0.44380621044460805701
gamma,1.0,3.0,1.0,0.06,0.340895,0.4359609283926595674
gamma,2.0,2.0,1.0,0.06,0.552177,0.48212037473443138236
gamma,3.0,1.0,1.0,0.06,0.709247,0.49027885728844261084
")
  factor_of <- function(family){
    if(family == "halfnormal") prior("normal", list(0, .5), list(0, Inf)) else prior("gamma", list(2, 2))
  }
  leaf <- function(row){
    spec <- BayesTools:::.prior_scale_product_spec(
      offset = 0, scale = 1, factor = factor_of(row$family),
      multiplier = prior("beta", list(row$alpha, row$beta)), sources = list(),
      map = list(type = "sqrt", scale = row$kappa)
    )
    list(type = "scale_product", spec = spec, n_grid = 1024L)
  }

  # Every reference row is checked. Each bound holds for every row when its
  # largest violation does, and the rows that violate a bound are named, so
  # one expectation per bound replaces one per row.
  cases <- paste(densities$family, densities$alpha, densities$beta, densities$kappa, densities$y)
  ordinates <- lapply(seq_len(nrow(densities)), function(i){
    BayesTools:::.prior_density_route_ordinate(leaf(densities[i, ]), densities$y[[i]])
  })
  heights <- vapply(ordinates, function(ordinate) exp(ordinate$log_density), numeric(1))
  absolute_errors <- vapply(ordinates, function(ordinate){
    ordinate$provenance$integration$absolute_error
  }, numeric(1))
  expect_identical(
    cases[vapply(ordinates, function(ordinate) !identical(ordinate$behavior, "regular"), logical(1))],
    character()
  )
  expect_identical(
    cases[!vapply(ordinates, function(ordinate) isTRUE(ordinate$exact), logical(1))],
    character()
  )
  expect_identical(cases[!(abs(heights / densities$density - 1) <= 1e-8)], character())
  expect_identical(cases[!(abs(heights - densities$density) <= absolute_errors)], character())
  # the batched quadrature of plotted densities (relative 1e-8 acceptance)
  keys <- paste(densities$family, densities$alpha, densities$beta, densities$kappa)
  plotted_errors <- vapply(unique(keys), function(key){
    rows <- densities[keys == key, ]
    plotted <- BayesTools:::.prior_density_route_density(leaf(rows[1L, ]), rows$y)
    max(abs(plotted / rows$density - 1))
  }, numeric(1))
  expect_identical(names(plotted_errors)[!(plotted_errors <= 1e-8)], character())
  probabilities <- lapply(seq_len(nrow(regions)), function(i){
    row <- regions[i, ]
    region <- list(
      intervals = matrix(c(row$lower, row$upper), 1L),
      indicator = function(x) x > row$lower & x < row$upper
    )
    BayesTools:::.prior_density_route_region(leaf(row), region)
  })
  region_cases <- paste(regions$family, regions$alpha, regions$beta, regions$kappa, regions$lower, regions$upper)
  expect_identical(
    region_cases[!vapply(probabilities, function(probability) isTRUE(probability$converged), logical(1))],
    character()
  )
  region_errors <- vapply(seq_along(probabilities), function(i){
    abs(probabilities[[i]]$probability / regions$probability[[i]] - 1)
  }, numeric(1))
  expect_identical(region_cases[!(region_errors <= 1e-8)], character())

  # the offset 0: f_T(0) E[(k S)^(-1/2)] = f_T(0) B(alpha - 1/2, beta) /
  # (sqrt(k) B(alpha, beta)) for alpha > 1/2; infinite for alpha <= 1/2, where
  # the mapped share has a positive or infinite density at 0
  problems <- expectation_problems({
    for(i in which(densities$family == "halfnormal" & densities$y == 1e-3)){
      row <- densities[i, ]
      ordinate <- BayesTools:::.prior_density_route_ordinate(leaf(row), 0)
      if(row$alpha > 1 / 2){
        expect_identical(ordinate$behavior, "regular", info = paste(row$alpha, row$beta, row$kappa))
        expect_equal(
          ordinate$log_density,
          log(2 * stats::dnorm(0, sd = .5)) + lbeta(row$alpha - 1 / 2, row$beta) -
            lbeta(row$alpha, row$beta) - log(row$kappa) / 2,
          tolerance = 1e-12,
          info = paste(row$alpha, row$beta, row$kappa)
        )
      }else{
        expect_identical(ordinate$behavior, "infinite", info = paste(row$alpha, row$beta, row$kappa))
      }
    }
  })
  expect_identical(problems, character())
  # outside the support
  expect_identical(
    BayesTools:::.prior_density_route_ordinate(leaf(densities[1L, ]), -.1)$behavior,
    "zero"
  )

  # the square of the product (a variance) through an 'exp_lin' node: its
  # density at v is f(sqrt(v)) / (2 sqrt(v)), and infinite at 0 where f(0+)
  # is finite and positive
  row <- densities[densities$family == "halfnormal" & densities$alpha == 2 &
                     densities$beta == 2 & densities$kappa == 1 & densities$y == .06, ]
  squared <- BayesTools:::.prior_density_route_transform(
    leaf(row), "exp_lin", list(a = 0, b = 2),
    function() BayesTools:::.prior_scale_product_hull(leaf(row)$spec)
  )
  ordinate <- BayesTools:::.prior_density_route_ordinate(squared, .06^2)
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density), row$density / (2 * .06), tolerance = 1e-8)
  expect_identical(BayesTools:::.prior_density_route_ordinate(squared, 0)$behavior, "infinite")
})

test_that("quadrature ordinates are exact only within a relative error bound", {

  # half-t3(0.5) times the square root of a Beta(2, 0.7) share at y = 1e3: a
  # density of 5.5e-13, which the former absolute floor (1e-12 + 1e-4 x value)
  # accepted with a reported error of 0.2 %; the quadrature is now refined
  # against the value. Reference: mpmath at 50 digits
  # (.work/tmp/pr58-decisions/P4-G1-review/mp_refs.csv).
  spec <- BayesTools:::.prior_scale_product_spec(
    offset = 0, scale = 1, factor = prior("t", list(0, .5, 3), list(0, Inf)),
    multiplier = prior("beta", list(2, .7)), sources = list(),
    map = list(type = "sqrt", scale = 1)
  )
  route <- list(type = "scale_product", spec = spec, n_grid = 1024L)
  for(case in list(c(1e3, 5.47320150023372434679774351794e-13),
                   c(1e5, 5.47320834105333987681610471592e-21))){
    ordinate <- BayesTools:::.prior_density_route_ordinate(route, case[[1L]])
    integration <- ordinate$provenance$integration
    expect_true(ordinate$exact)
    expect_true(integration$refined)
    expect_lte(integration$absolute_error, 1e-4 * exp(ordinate$log_density))
    expect_lte(integration$error_bound, 1e-4 * exp(ordinate$log_density) * (1 + 1e-12))
    expect_lte(abs(exp(ordinate$log_density) / case[[2L]] - 1), 1e-6)
  }
  # ordinary values keep their first-pass quadrature (no refinement)
  ordinary <- BayesTools:::.prior_density_route_ordinate(route, .5)$provenance$integration
  expect_false(ordinary$refined)
  expect_lte(ordinary$absolute_error, 1e-4 * exp(BayesTools:::.prior_density_route_ordinate(route, .5)$log_density))

  # a total whose reported error misses the relative criterion is rejected
  # with that reason, and only a plotted curve draws its estimate
  quadrature <- BayesTools:::.prior_conditional_normal_quadrature(
    function(x) rep(1e-3, length(x)), c(0, 1), 64L, zero_message = "zero"
  )
  expect_true(quadrature$integration$converged)
  # a pure scale mixture b * s, b ~ N(1, .5), s ~ gamma(2, 1) on (0, 10), at
  # its offset 0: the ordinate is dnorm(2) / .5 * E[1 / s], with E[1 / s] by
  # quadrature (the truncated gamma has no closed form)
  slope <- prior("normal", list(1, .5))
  attr(slope, "multiply_by") <- "s"
  scale_mixture <- BayesTools:::.prior_linear_combination_density(
    list(b = slope, s = prior("gamma", list(2, 1), list(0, 10))), c(b = 1), n_grid = 1024
  )
  scale_route <- BayesTools:::.prior_density_route_from_adaptive(
    attr(scale_mixture, "adaptive_evaluation", exact = TRUE)
  )
  offset <- BayesTools:::.prior_density_route_ordinate(scale_route, 0)
  expect_true(offset$exact)
  expect_identical(offset$provenance$inverse_moment$method, "quadrature")
  expect_equal(exp(offset$log_density),
               stats::dnorm(2) / .5 * offset$provenance$inverse_moment$value)
  expect_equal(BayesTools:::.prior_density_route_quadrature_density(scale_route, 0),
               exp(offset$log_density))
  local_mocked_bindings(
    .prior_conditional_normal_piece = function(integrand, lower, upper, n_grid, relative, absolute){
      list(value = 1e-3, abs.error = 1e-6, message = "OK", evaluations = 21L)
    }
  )
  rejected <- BayesTools:::.prior_conditional_normal_quadrature(
    function(x) rep(1e-3, length(x)), c(0, 1), 64L, zero_message = "zero"
  )
  expect_true(is.na(rejected$value))
  expect_false(rejected$integration$converged)
  expect_true(rejected$integration$refined)
  # the message states the observed relative error (1e-6 / 1e-3), the
  # criterion stays in the error bound
  expect_identical(rejected$integration$message, "relative error estimate 0.001")
  expect_equal(rejected$integration$error_bound, 1e-4 * 1e-3)
  expect_identical(rejected$integration$estimate, 1e-3)
  ordinate <- BayesTools:::.prior_density_route_ordinate(route, 1e3)
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "integration reported 'relative error estimate 0.001'", fixed = TRUE)
  # the mocked pieces sum to the estimate the plotted curve draws
  estimate <- ordinate$provenance$integration$estimate
  expect_equal(estimate, 1e-3 * (length(ordinate$provenance$integration$breakpoints) - 1L))
  expect_equal(BayesTools:::.prior_density_route_quadrature_density(route, 1e3), estimate)
  # the estimate of a rejected offset quadrature is the inverse moment, not
  # the density, so the plotted curve omits the offset
  offset <- BayesTools:::.prior_density_route_ordinate(scale_route, 0)
  expect_false(offset$exact)
  expect_true(is.numeric(offset$provenance$integration$estimate))
  expect_true(is.na(BayesTools:::.prior_density_route_quadrature_density(scale_route, 0)))
})
test_that("R116 N06 nested named transformations retain declared support", {
  normal <- .prior_linear_combination_density(list(x = prior("normal", list(0, 1))), c(x = 1), n_grid = 128)
  exponential <- .prior_density_output_transform(normal, .bt_posterior_transformation("exp"))
  squared <- .prior_density_output_transform(exponential, .bt_posterior_transformation("exp_lin", list(b = 2)))
  ordinate <- prior_density_ordinate(squared, 1)
  expect_identical(ordinate$behavior, "regular")
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density), stats::dlnorm(1, 0, 2), tolerance = 1e-12)
  route <- .prior_density_route_from_adaptive(attr(squared, "adaptive_evaluation"))
  provenance <- .prior_density_route_provenance(route)
  expect_identical(provenance$arguments, list(a = 0, b = 2))
  expect_equal(.prior_density_ordinate_provenance_support(provenance), c(lower = 0, upper = Inf))
  expect_identical(.prior_density_ordinate_provenance_atoms(provenance), numeric())
  expect_false(any(vapply(provenance, is.environment, logical(1))))
  expect_false(any(vapply(provenance, is.function, logical(1))))
  positive <- .prior_linear_combination_density(list(x = prior("exp", list(1)), y = prior("exp", list(1))), c(x = 1, y = 1), n_grid = 128)
  square <- .prior_density_output_transform(positive, .bt_posterior_transformation("exp_lin", list(b = 2)))
  expect_equal(exp(prior_density_ordinate(square, 1)$log_density), exp(-1) / 2, tolerance = 1e-8)
  convolution <- .prior_density_route_provenance(.prior_density_route_from_adaptive(attr(positive, "adaptive_evaluation")))
  expect_identical(convolution$kind, "convolution")
  expect_length(convolution$terms, 2L)
  expect_equal(.prior_density_ordinate_provenance_support(convolution), c(lower = 0, upper = Inf))
  convolution$weights <- c(-2, 0)
  convolution$offset <- 3
  convolution$terms[[2L]] <- list(kind = "unsupported_provenance")
  expect_equal(.prior_density_ordinate_provenance_support(convolution), c(lower = -Inf, upper = 3))
  convolution$weights[2L] <- 1
  expect_null(.prior_density_ordinate_provenance_support(convolution))
  unknown <- .prior_density_ordinate_named_transform(function(x) NULL, list(kind = "unsupported_provenance"), "exp_lin", list(b = 2), 1)
  expect_identical(unknown$behavior, "unknown")
  expect_false(unknown$exact)
  expect_match(unknown$reason, "structural-domain")
  negative <- .prior_linear_combination_density(list(x = prior("normal", list(0, 1), list(0, Inf))), c(x = 1), n_grid = 128)
  negative <- .prior_density_output_transform(negative, .bt_posterior_transformation("lin", list(b = -1)))
  negative_provenance <- .prior_density_route_provenance(.prior_density_route_from_adaptive(attr(negative, "adaptive_evaluation")))
  expect_equal(.prior_density_ordinate_provenance_support(negative_provenance), c(lower = -Inf, upper = 0))
  refused <- .prior_density_ordinate_named_transform(function(x) NULL, negative_provenance, "exp_lin", list(b = 2), 1)
  expect_identical(refused$behavior, "undefined")
  expect_true(refused$exact)
})

test_that("R116 N06 constants and finite images have conservative provenance", {
  constant <- list(kind = "scalar_affine", offset = 2, scale = 0)
  expect_equal(.prior_density_ordinate_provenance_support(constant), c(lower = 2, upper = 2))
  expect_identical(.prior_density_ordinate_provenance_atoms(constant), 2)
  expect_identical(.prior_density_ordinate_exp_lin_boundary(constant, 2), "zero")
  expect_identical(.prior_density_ordinate_lower_log_coefficient(constant, 2), -Inf)
  for(location in c(0, 2)){
    density <- .prior_linear_combination_density(list(x = prior("point", list(location))), c(x = 1), n_grid = 128)
    source <- .prior_density_route_provenance(.prior_density_route_from_adaptive(attr(density, "adaptive_evaluation")))
    expect_equal(.prior_density_ordinate_provenance_support(source), c(lower = location, upper = location))
    expect_equal(.prior_density_ordinate_provenance_atoms(source), location)
    expect_identical(prior_density_ordinate(density, location)$provenance$continuous_behavior, "zero")
  }
  point <- .prior_linear_combination_density(list(x = prior("point", list(2))), c(x = 1), n_grid = 128, output_transformation = "exp_lin", output_transformation_arguments = list(b = 2))
  expect_identical(prior_density_ordinate(point, 4)$behavior, "point_mass")
  expect_identical(prior_density_ordinate(point, 3)$behavior, "zero")
  for(transformation in c("exp", "exp_lin", "tanh", "lin")){
    value <- if(transformation == "lin") 1e308 else 1000
    source <- list(kind = "scalar_affine", offset = value, scale = 0)
    mapped <- list(kind = "named_transform", transformation = transformation, arguments = list(a = 0, b = 2), source = source)
    expect_null(.prior_density_ordinate_provenance_support(mapped))
    expect_null(.prior_density_ordinate_provenance_atoms(mapped))
    source$offset <- -1000
    mapped$source <- source
    if(transformation == "exp"){
      expect_null(.prior_density_ordinate_provenance_support(mapped))
      expect_null(.prior_density_ordinate_provenance_atoms(mapped))
    }
  }
})

test_that("R116 N89 point-only densities are exact off their stored atom", {

  point <- .prior_linear_density_point(.5)
  at <- prior_density_ordinate(point, .5)
  off <- prior_density_ordinate(point, .3)
  expect_identical(at$behavior, "point_mass")
  expect_identical(at$provenance$continuous_behavior, "zero")
  expect_identical(at$point_mass, 1)
  expect_identical(at$log_density, -Inf)
  expect_true(at$exact)
  expect_identical(off$behavior, "zero")
  expect_true(off$exact)
  expect_identical(off$log_density, -Inf)
  expect_identical(off$point_mass, 0)
  grid <- point
  grid$density <- list(x = c(0, 1), y = c(1, 1), mass = 1)
  grid$points <- NULL
  positive <- prior_density_ordinate(grid, .3)
  expect_identical(positive$behavior, "unknown")
  expect_false(positive$exact)
  deferred <- point
  deferred$density <- list()
  deferred$points <- NULL
  unavailable <- prior_density_ordinate(deferred, .3)
  expect_identical(unavailable$behavior, "unknown")
  expect_false(unavailable$exact)
})
