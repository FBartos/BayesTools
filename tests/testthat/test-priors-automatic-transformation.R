skip_if_not_test_profile("unit")

automatic_log_map <- list(fun = log, inv = exp, jac = function(x) 1 / x)
automatic_neglog_map <- list(fun = function(x) -log(x), inv = function(x) exp(-x),
                             jac = function(x) -1 / x)
automatic_halfnormal <- function(sd = 1, weight = 1){
  prior("normal", list(0, sd), list(0, Inf), prior_weights = weight)
}
automatic_list_data <- function(priors, ...){
  args <- list(prior_list = priors, x_seq = NULL, x_range = NULL,
               x_range_quant = NULL, n_points = 80, n_samples = 500,
               force_samples = FALSE, individual = TRUE,
               transformation = automatic_log_map, transformation_arguments = NULL,
               transformation_settings = TRUE)
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(BayesTools:::.plot_data_prior_list.simple, args)
}
automatic_expect_halfnormal <- function(d, alpha = .005, decreasing = FALSE){
  source <- if(decreasing) exp(-d$x) else exp(d$x)
  expect_equal(d$y, 2 * dnorm(source) * source, tolerance = 1e-13)
  expect_equal(min(source), qnorm((1 + alpha) / 2), tolerance = 1e-13)
  expect_equal(2 * pnorm(min(source)) - 1, alpha, tolerance = 1e-13)
  expect_true(all(is.finite(attr(d, "x_range"))))
  expect_true(diff(attr(d, "x_range")) > 0)
}

for(settings in c(FALSE, TRUE)){
  test_that(paste("automatic density log range uses source coordinates", settings), {
    p <- automatic_halfnormal()
    original <- serialize(p, NULL)
    d <- density(p, transformation = automatic_log_map,
                 transformation_settings = settings, n_points = 80)
    automatic_expect_halfnormal(d)
    expect_equal(attr(d, "x_range"), range(d$x))
    expect_true(any(d$x < -5))
    expect_identical(serialize(p, NULL), original)
  })
}

test_that("automatic decreasing ranges preserve paired ordinates and requested tail mass", {
  p <- automatic_halfnormal()
  d <- density(p, transformation = automatic_neglog_map,
               transformation_settings = TRUE, x_range_quant = .02, n_points = 80)
  automatic_expect_halfnormal(d, alpha = .02, decreasing = TRUE)
  expect_true(all(diff(d$x) < 0))
  expect_equal(attr(d, "x_range"), range(d$x))
})

for(route in c("plot", "plot_list", "lines", "lines_list", "geom", "geom_list")){
  for(settings in c(FALSE, TRUE)){
    test_that(paste("automatic public prior range reaches negative log coordinates", route, settings), {
      p <- automatic_halfnormal()
      args <- list(transformation = automatic_log_map,
                   transformation_settings = settings, n_points = 80)
      captured <- NULL
      testthat::local_mocked_bindings(
        .plot.prior.simple = function(x, plot_type, plot_data, ...) captured <<- plot_data,
        .plot_prior_list.both = function(plot_data, ...) captured <<- plot_data$density,
        .lines.prior.simple = function(plot_data, ...) captured <<- plot_data,
        .package = "BayesTools"
      )
      if(route == "plot") do.call(plot, c(list(x = p, plot_type = "base"), args))
      if(route == "plot_list") do.call(plot_prior_list, c(list(prior_list = list(p), plot_type = "base"), args))
      if(route == "lines") do.call(lines, c(list(x = p), args))
      if(route == "lines_list") do.call(lines_prior_list, c(list(prior_list = list(p)), args))
      if(route == "geom"){
        g <- ggplot2::ggplot() + do.call(geom_prior, c(list(x = p), args))
        captured <- ggplot2::ggplot_build(g)$data[[1L]]
      }
      if(route == "geom_list"){
        g <- ggplot2::ggplot() + do.call(geom_prior_list, c(list(prior_list = list(p)), args))
        captured <- ggplot2::ggplot_build(g)$data[[1L]]
      }
      expect_true(min(captured$x) < -5)
      expect_equal(min(exp(captured$x)), mquant(p, .005), tolerance = 1e-13)
      expect_equal(captured$y, 2 * dnorm(exp(captured$x)) * exp(captured$x), tolerance = 1e-13)
    })
  }
}

test_that("actual direct and list ggplots have the automatic transformed curve", {
  p <- automatic_halfnormal()
  for(make in list(function(...) plot(p, plot_type = "ggplot", ...),
                  function(...) plot_prior_list(list(p), plot_type = "ggplot", ...))){
    g <- make(transformation = automatic_log_map, transformation_settings = TRUE,
              n_points = 80, x_range_quant = .02)
    layer <- ggplot2::ggplot_build(g)$data[[1L]]
    expect_equal(min(exp(layer$x)), mquant(p, .02), tolerance = 1e-13)
    expect_equal(layer$y, 2 * dnorm(exp(layer$x)) * exp(layer$x), tolerance = 1e-13)
  }
})

test_that("explicit transformed and source coordinates remain authoritative", {
  p <- automatic_halfnormal()
  for(settings in c(FALSE, TRUE)){
    source <- seq(.1, 2, length.out = 30)
    knots <- if(settings) log(source) else source
    bounds <- if(settings) c(-3, 1) else c(.04, 3)
    for(mode in c("seq", "range", "both")){
      args <- list(transformation = automatic_log_map,
                   transformation_settings = settings, truncate_end = FALSE, n_points = 30)
      if(mode != "range") args$x_seq <- knots
      if(mode != "seq") args$x_range <- bounds
      d <- do.call(density, c(list(x = p), args))
      expected_x <- if(mode == "range"){
        if(settings) seq(bounds[1], bounds[2], length.out = 30) else
          log(seq(bounds[1], bounds[2], length.out = 30))
      }else log(source)
      expect_equal(d$x, expected_x, tolerance = 1e-13)
      expect_equal(d$y, 2 * dnorm(exp(d$x)) * exp(d$x), tolerance = 1e-13)
      expect_equal(attr(d, "x_range"), if(settings){
        if(mode == "seq") range(knots) else bounds
      }else log(if(mode == "seq") range(source) else bounds))
    }
    for(route in c("direct", "list", "geom", "geom_list")){
      args <- list(transformation = automatic_log_map, transformation_settings = settings,
                   x_seq = knots, xlim = bounds, n_points = 30)
      g <- switch(route,
        direct = do.call(plot, c(list(x = p, plot_type = "ggplot"), args)),
        list = do.call(plot_prior_list, c(list(prior_list = list(p), plot_type = "ggplot"), args)),
        geom = ggplot2::ggplot() + do.call(geom_prior, c(list(x = p), args)),
        geom_list = ggplot2::ggplot() + do.call(geom_prior_list, c(list(prior_list = list(p)), args)))
      layer <- ggplot2::ggplot_build(g)$data[[1L]]
      expect_equal(layer$x, log(source), tolerance = 1e-13)
      expect_equal(layer$y, 2 * dnorm(source) * source, tolerance = 1e-13)
    }
  }
})

test_that("finite automatic source grids remain unchanged for named exp", {
  p <- automatic_halfnormal()
  automatic <- density(p, transformation = "exp", n_points = 80)
  source <- density(p, x_range = range(p), transformation = "exp", n_points = 80)
  expect_identical(automatic$x, source$x)
  expect_identical(automatic$y, source$y)
  expect_identical(attr(automatic, "x_range"), attr(source, "x_range"))
  flagged <- density(p, transformation = "exp", transformation_settings = TRUE, n_points = 80)
  expect_identical(flagged$x, automatic$x)
  expect_identical(flagged$y, automatic$y)
})

test_that("automatic mixtures use common marginal quantile knots with original weights and atoms", {
  priors <- list(automatic_halfnormal(weight = 1), automatic_halfnormal(2, 3),
                 prior("point", list(2), prior_weights = 2))
  original <- serialize(priors, NULL)
  d <- automatic_list_data(priors)
  x <- exp(d$density$x)
  expected <- (2 * dnorm(x) + 3 * 2 * dnorm(x, sd = 2)) * x / 6
  expect_equal(d$density$y, expected, tolerance = 1e-13)
  expect_equal(min(x), mquant(priors[[1]], .005), tolerance = 1e-13)
  expect_equal(d$points1$x, log(2))
  expect_equal(d$points1$y, 1 / 3)
  expect_equal(attr(d$density, "x_range"), range(d$density$x))
  expect_equal(attr(d$points1, "x_range"), rep(log(2), 2))
  decreasing <- automatic_list_data(priors, transformation = automatic_neglog_map)
  x_decreasing <- exp(-decreasing$density$x)
  expect_true(all(diff(decreasing$density$x) > 0))
  expect_equal(decreasing$density$y,
    (2 * dnorm(x_decreasing) + 3 * 2 * dnorm(x_decreasing, sd = 2)) * x_decreasing / 6,
    tolerance = 1e-13)
  expect_true(diff(attr(decreasing$density, "x_range")) > 0)
  expect_identical(serialize(priors, NULL), original)
})

test_that("automatic endpoint selection never moves an invalid atom", {
  point <- prior("point", list(0))
  p <- automatic_halfnormal()
  spike <- prior_spike_and_slab(p, prior("bernoulli", list(.5)))
  for(object in list(point, spike)){
    error <- tryCatch(density(object, transformation = automatic_log_map,
                              transformation_settings = TRUE, n_points = 30), error = identity)
    expect_s3_class(error, "BayesTools_transformation_image_unavailable")
    expect_equal(error$source_values, 0)
    expect_equal(error$images, -Inf)
  }
  expect_error(automatic_list_data(list(p, point)),
               class = "BayesTools_transformation_image_unavailable")
  expect_error(automatic_list_data(list(prior("bernoulli", list(.5)))),
               class = "BayesTools_transformation_image_unavailable")
})

test_that("automatic derivation validates supplied transformation settings", {
  p <- automatic_halfnormal()
  for(bad in list(NA, 1, c(TRUE, FALSE))){
    for(call in list(
      function() density(p, transformation = automatic_log_map, transformation_settings = bad),
      function() plot(p, transformation = automatic_log_map, transformation_settings = bad),
      function() lines(p, transformation = automatic_log_map, transformation_settings = bad),
      function() geom_prior(p, transformation = automatic_log_map, transformation_settings = bad),
      function() plot_prior_list(list(p), transformation = automatic_log_map, transformation_settings = bad),
      function() lines_prior_list(list(p), transformation = automatic_log_map, transformation_settings = bad),
      function() geom_prior_list(list(p), transformation = automatic_log_map, transformation_settings = bad))){
      expect_error(call(), "transformation_settings", fixed = TRUE)
    }
  }
})

test_that("seeded automatic transformed samples use the existing source KDE", {
  p <- automatic_halfnormal()
  bounds <- c(mquant(p, .005), range(p)[2])
  set.seed(728)
  automatic <- density(p, transformation = automatic_log_map, transformation_settings = TRUE,
                       force_samples = TRUE, n_samples = 500, n_points = 80)
  set.seed(728)
  source <- density(p, x_range = bounds, force_samples = TRUE,
                    n_samples = 500, n_points = 80)
  expect_identical(automatic$x, log(source$x))
  expect_equal(automatic$y, source$y * source$x, tolerance = 1e-13)
  expect_identical(automatic$samples, log(source$samples))
  expect_identical(attr(automatic, "boundary_reflection"), attr(source, "boundary_reflection"))
})

test_that("finite automatic public source grids and explicit decreasing limits stay authoritative", {
  p <- automatic_halfnormal()
  g <- plot(p, plot_type = "ggplot", transformation = "exp", n_points = 80)
  expected <- density(p, x_range = range(pretty(range(p))),
                      transformation = "exp", n_points = 80)
  layer <- ggplot2::ggplot_build(g)$data[[1L]]
  expect_identical(layer$x, expected$x)
  expect_identical(layer$y, expected$y)
  decreasing <- density(p, transformation = automatic_neglog_map,
                         transformation_settings = TRUE, x_range = c(-2, 4), n_points = 80)
  expect_identical(attr(decreasing, "x_range"), c(-2, 4))
  expect_equal(decreasing$y, 2 * dnorm(exp(-decreasing$x)) * exp(-decreasing$x), tolerance = 1e-13)

  captured <- NULL
  testthat::local_mocked_bindings(
    .lines.prior.simple = function(plot_data, ...) captured <<- plot_data,
    .package = "BayesTools"
  )
  for(settings in c(FALSE, TRUE)){
    source <- seq(.1, 2, length.out = 30)
    knots <- if(settings) log(source) else source
    bounds <- if(settings) c(-3, 1) else c(.04, 3)
    for(route in list(function(...) lines(p, ...), function(...) lines_prior_list(list(p), ...))){
      route(x_seq = knots, xlim = bounds, transformation = automatic_log_map,
            transformation_settings = settings)
      expect_equal(captured$x, log(source), tolerance = 1e-13)
      expect_equal(captured$y, 2 * dnorm(source) * source, tolerance = 1e-13)
    }
  }
})

test_that("common automatic source knots include valid atoms outside component tail windows", {
  priors <- list(automatic_halfnormal(), automatic_halfnormal(2, 3),
                 prior("point", list(.001), prior_weights = 2))
  d <- automatic_list_data(priors)
  expect_equal(min(exp(d$density$x)), .001, tolerance = 1e-13)
  expect_equal(d$points1$x, log(.001))
  bounds <- range(exp(d$density$x))
  components <- lapply(priors[1:2], density, x_range = bounds,
    transformation = automatic_log_map, n_points = 80, truncate_end = FALSE)
  expect_identical(components[[1]]$x, components[[2]]$x)
  expect_equal(d$density$x, components[[1]]$x, tolerance = 1e-13)
  expect_equal(d$density$y, (components[[1]]$y + 3 * components[[2]]$y) / 6,
               tolerance = 1e-13)
})

test_that("automatic ranges retain specialized axes and individual coefficient dispatch", {
  pet <- prior_PET("normal", list(0, 1), list(0, Inf))
  axis <- seq(0, 1, length.out = 30)
  specialized <- density(pet, transformation = "exp", n_points = 30)
  expect_identical(specialized$x, axis)
  individual <- density(pet, individual = TRUE, transformation = automatic_log_map,
                         transformation_settings = TRUE, n_points = 30)
  automatic_expect_halfnormal(individual)
  weight <- prior_weightfunction("one-sided", steps = c(.05, .5),
                                 weights = wf_fixed(c(1, .6, .2)))
  original <- density(weight, n_points = 30)
  mapped <- density(weight, transformation = "exp", n_points = 30)
  expect_identical(mapped$x, original$x)
  expect_identical(attr(mapped, "x_range"), c(0, 1))
  for(settings in c(FALSE, TRUE)){
    figure <- plot_prior_list(list(pet), prior_list_mu = list(prior("point", list(1))),
                              transformation = "exp", transformation_settings = settings,
                              plot_type = "ggplot", n_points = 8)
    layer <- ggplot2::ggplot_build(figure)$data[[1L]]
    expect_equal(range(layer$x), c(0, 1))
  }
})

test_that("automatic endpoint quantiles do not repair invalid or unrepresentable maps", {
  p <- automatic_halfnormal()
  invalid <- list(fun = function(x) log(x - 1), inv = function(x) exp(x) + 1,
                   jac = function(x) 1 / (x - 1))
  warnings <- character()
  error <- tryCatch(withCallingHandlers(density(p, transformation = invalid), warning = function(w){
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  expect_s3_class(error, "BayesTools_transformation_domain")
  expect_true(length(warnings) > 0L && all(warnings == "NaNs produced"))
  unavailable <- automatic_log_map
  unavailable$fun <- function(x) ifelse(x < .1, -Inf, log(x))
  expect_error(density(p, transformation = unavailable),
               class = "BayesTools_transformation_domain")
  expect_error(density(p, transformation = unavailable),
               class = "BayesTools_prior_range_unavailable")
})

test_that("an unavailable automatic transformed prior range is a classed refusal that names explicit limits", {
  wide <- prior("normal", list(0, 1000))
  # the 0.995 quantile, 2576, has no finite image under exp
  expect_identical(exp(mquant(wide, .995)), Inf)
  error <- tryCatch(density(wide, transformation = "exp"), error = identity)
  expect_s3_class(error, "BayesTools_prior_range_unavailable")
  expect_s3_class(error, "BayesTools_transformation_domain")
  expect_s3_class(error, "BayesTools_transformation")
  expect_null(error$call)
  expect_match(conditionMessage(error), "x_range", fixed = TRUE)
  expect_match(conditionMessage(error), "x_seq", fixed = TRUE)
  expect_match(conditionMessage(error), "xlim", fixed = TRUE)
  expect_false(grepl("posterior draws", conditionMessage(error), fixed = TRUE))
  # the prior plots refuse in the same way
  expect_error(plot_prior_list(list(wide), plot_type = "ggplot", transformation = "exp"),
               class = "BayesTools_prior_range_unavailable")
  expect_error(plot(wide, plot_type = "ggplot", transformation = "exp"),
               class = "BayesTools_prior_range_unavailable")
  # explicit limits avoid the refusal and keep their interpretation
  explicit <- density(wide, transformation = "exp", x_range = c(-5, 5), n_points = 50)
  expect_equal(explicit$x, exp(seq(-5, 5, length.out = 50)))
  # a prior that stays within the representable images keeps its automatic range
  moderate <- density(prior("normal", list(0, 100)), transformation = "exp", n_points = 50)
  expect_true(all(is.finite(moderate$x)) && all(is.finite(moderate$y)))
})

test_that("automatic list endpoints use only existing participating components", {
  p <- automatic_halfnormal()
  captured <- NULL
  testthat::local_mocked_bindings(
    .plot_prior_list.both = function(plot_data, ...) captured <<- plot_data,
    .lines.prior.simple = function(plot_data, ...) captured <<- list(density = plot_data),
    .package = "BayesTools"
  )
  for(force_samples in c(FALSE, TRUE)){
    atom <- prior("point", list(0), prior_weights = .001)
    if(!force_samples) atom$prior_weights <- 0
    priors <- list(p, atom)
    args <- list(transformation = automatic_log_map, transformation_settings = TRUE,
                 n_points = 30, n_samples = 100, force_samples = force_samples)
    for(route in c("plot", "lines", "geom")){
      set.seed(638)
      if(route == "plot"){
        expect_silent(do.call(plot_prior_list, c(list(prior_list = priors), args)))
        expect_named(captured, "density")
      }else if(route == "lines"){
        expect_silent(do.call(lines_prior_list, c(list(prior_list = priors), args)))
        expect_named(captured, "density")
      }else{
        expect_silent(g <- ggplot2::ggplot() + do.call(geom_prior_list,
          c(list(prior_list = priors), args)))
        expect_length(ggplot2::ggplot_build(g)$data, 1L)
      }
    }
  }
  # A direct point's external model weight cannot remove its own atom.
  atom$prior_weights <- 0
  expect_error(density(atom, transformation = automatic_log_map),
               class = "BayesTools_transformation_image_unavailable")
})

test_that("automatic decreasing spike bounds retain nested paired densities", {
  p <- prior_spike_and_slab(automatic_halfnormal(), prior("point", list(.5)))
  d <- density(p, transformation = "lin", transformation_arguments = list(a = 0, b = -1),
               transformation_settings = TRUE, n_points = 30)
  expect_true(diff(attr(d, "x_range")) > 0)
  expect_true(diff(attr(d$variable, "x_range")) > 0)
  expect_true(all(diff(d$variable$x) <= 0))
  expect_equal(d$inclusion$x, 0)
  expect_equal(d$inclusion$y, .5)
})

for(settings in c(FALSE, TRUE)){
  test_that(paste("individual weightfunction automatic ranges retain their source route", settings), {
    weights <- list(wf_fixed(c(1, .6, .2)), wf_fixed(c(1, .6, 0)), wf_cumulative(c(1, 1, 1)))
    for(weight in weights){
      p <- prior_weightfunction("one-sided", steps = c(.05, .5), weights = weight)
      original <- serialize(p, NULL)
      source <- density(p, individual = TRUE, n_points = 8, n_samples = 32)
      mapped <- tryCatch(density(p, individual = TRUE, transformation = automatic_log_map,
        transformation_settings = settings, n_points = 8, n_samples = 32), error = identity)
      expect_false(inherits(mapped, "error"))
      if(inherits(mapped, "error")) next
      expect_identical(lapply(mapped, function(d) list(x = d$x, y = d$y,
        x_range = attr(d, "x_range"), y_range = attr(d, "y_range"))),
        lapply(source, function(d) list(x = d$x, y = d$y,
          x_range = attr(d, "x_range"), y_range = attr(d, "y_range"))))
      expect_identical(attr(mapped[[2L]], "x_range"), c(0, 1))
      expect_identical(serialize(p, NULL), original)
    }
  })

  test_that(paste("individual weightfunction wrappers keep automatic source data", settings), {
    coordinates <- function(g) lapply(ggplot2::ggplot_build(g)$data, function(d){
      d[intersect(c("x", "y", "xend", "yend"), names(d))]
    })
    captured <- NULL
    testthat::local_mocked_bindings(
      .lines.prior.simple = function(plot_data, ...) captured <<- plot_data,
      .lines.prior.point = function(plot_data, ...) captured <<- plot_data,
      .package = "BayesTools"
    )
    for(weight in list(wf_fixed(c(1, .6, .2)), wf_fixed(c(1, .6, 0)), wf_cumulative(c(1, 1, 1)))){
      p <- prior_weightfunction("one-sided", steps = c(.05, .5), weights = weight)
      for(route in c("plot", "lines", "geom")){
        draw <- function(transformation){
          args <- list(individual = TRUE, transformation = transformation,
                       transformation_settings = settings, n_points = 8, n_samples = 32)
          if(route == "plot") return(coordinates(do.call(plot,
            c(list(x = p, plot_type = "ggplot", show_figures = 3L), args))))
          if(route == "geom") return(coordinates(ggplot2::ggplot() + do.call(geom_prior,
            c(list(x = p, show_parameter = 3L), args))))
          do.call(lines, c(list(x = p, show_parameter = 3L), args))
          list(x = captured$x, y = captured$y, x_range = attr(captured, "x_range"))
        }
        source <- draw(NULL)
        mapped <- tryCatch(draw(automatic_log_map), error = identity)
        expect_identical(mapped, source)
      }
    }
  })
}

test_that("individual weightfunction explicit grids retain their established interpretation", {
  p <- prior_weightfunction("one-sided", steps = c(.05, .5), weights = wf_cumulative(c(1, 1, 1)))
  for(settings in c(FALSE, TRUE)){
    bounds <- if(settings) c(-2, -.1) else c(.1, .9)
    supplied <- seq(bounds[1], bounds[2], length.out = 8)
    source <- if(settings) exp(supplied) else supplied
    for(mode in c("range", "seq", "both")){
      args <- list(individual = TRUE, transformation = automatic_log_map,
                   transformation_settings = settings, n_points = 8)
      if(mode != "range") args$x_seq <- supplied
      if(mode != "seq") args$x_range <- bounds
      d <- do.call(density, c(list(x = p), args))
      expect_equal(d[[2L]]$x, source, tolerance = 1e-14)
      expect_equal(d[[2L]]$y, 2 * source, tolerance = 1e-14)
      expect_equal(d[[3L]]$x, source, tolerance = 1e-14)
      expect_equal(d[[3L]]$y, 2 * (1 - source), tolerance = 1e-14)
      expect_identical(attr(d[[2L]], "x_range"), bounds)
      expect_identical(c(d[[1L]]$x, d[[1L]]$y), c(1, 1))
    }
  }
})

test_that("automatic ordered ranges preserve the typed reference atom refusal", {
  levels <- c("a", "b", "c")
  data <- data.frame(f = ordered(rep(levels, 3L), levels = levels))
  p <- JAGS_formula(~f, "mu", data, list(intercept = prior("normal", list(0, 1)),
    f = prior_ordered(prior("gamma", list(2, 1)),
      allocation = prior("dirichlet", list(c(1, 1))))))$prior_list$mu_f
  original <- serialize(p, NULL)
  expected_message <- paste0("The 'custom' transformation is numerically unavailable: finite source values ",
    "do not have representable images under this transformation. Use a transformation ",
    "whose images remain representable.")
  for(settings in c(FALSE, TRUE)){
    for(alpha in list(NULL, .005)){
      for(route in c("density", "plot", "lines", "geom")){
        args <- list(transformation = automatic_log_map, transformation_settings = settings,
                     x_range_quant = alpha, n_points = 8, n_samples = 32)
        set.seed(813)
        error <- tryCatch(switch(route,
          density = do.call(density, c(list(x = p), args)),
          plot = do.call(plot, c(list(x = p, show_figures = 2L), args)),
          lines = do.call(lines, c(list(x = p, show_parameter = 2L), args)),
          geom = do.call(geom_prior, c(list(x = p, show_parameter = 2L), args))), error = identity)
        expect_identical(class(error), c("BayesTools_transformation_image_unavailable",
          "BayesTools_transformation", "error", "condition"))
        expect_identical(list(message = conditionMessage(error), call = conditionCall(error),
          transformation = error$transformation, source_values = error$source_values, images = error$images),
          list(message = expected_message, call = NULL, transformation = "custom", source_values = 0, images = -Inf))
      }
    }
    # Explicit windows retain the same established atom refusal; positive model
    # weight cannot suppress an invalid scalar atom either.
    for(object in list(p, prior("point", list(0), prior_weights = 2))){
      set.seed(813)
      error <- tryCatch(density(object, x_range = if(settings) log(c(.1, 2)) else c(.1, 2),
        transformation = automatic_log_map, transformation_settings = settings,
        n_points = 8, n_samples = 32), error = identity)
      expect_identical(class(error), c("BayesTools_transformation_image_unavailable",
        "BayesTools_transformation", "error", "condition"))
      expect_identical(list(message = conditionMessage(error), call = conditionCall(error),
        transformation = error$transformation, source_values = error$source_values, images = error$images),
        list(message = expected_message, call = NULL, transformation = "custom", source_values = 0, images = -Inf))
    }
  }
  expect_identical(serialize(p, NULL), original)
})
