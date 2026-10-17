skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

test_that("transformed whole ordered factors plot attached current-scale laws", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))
  factor <- transform_factor_samples(fixture$samples)
  factor_exp <- factor
  factor_exp$mu_f <- posterior_transform(factor$mu_f, "exp")
  expect_identical(as.numeric(factor_exp$mu_f), exp(as.numeric(factor$mu_f)))
  expect_identical(posterior_metadata(factor_exp$mu_f, "ordered_source")$primitives,
    posterior_metadata(factor$mu_f, "ordered_source")$primitives)
  expect_identical(posterior_metadata(factor_exp$mu_f, "linear_weight_space"),
    posterior_metadata(factor$mu_f, "linear_weight_space"))
  laws <- posterior_metadata(factor_exp$mu_f, "prior_densities")
  expect_identical(names(laws), colnames(factor_exp$mu_f))
  reference <- prior_density_ordinate(laws[[1L]], 1)
  expect_identical(reference$behavior, "point_mass")
  expect_identical(reference$point_mass, 1)
  plotted <- plot_posterior(factor_exp, "mu_f", prior = TRUE,
    plot_type = "ggplot", n_points = 64L, n_samples = 128L)
  layers <- ggplot2::ggplot_build(plotted)$data
  expect_true(all(vapply(layers[1:2], function(layer) all(layer$x > 0), logical(1))))
  expect_equal_each(layers[[2L]]$y, stats::dnorm(log(layers[[2L]]$x), 0, .5) / layers[[2L]]$x,
    tolerance = 2e-12)
  mapped <- marginal_posterior(fixture$samples, "mu_f", formula = ~0 + f,
    prior_samples = TRUE, transformation = "exp", n_samples = 128L)
  marginal <- plot_marginal(list(mu_f = mapped), "mu_f", prior = TRUE,
    plot_type = "ggplot", n_points = 64L, n_samples = 128L)
  marginal_layers <- ggplot2::ggplot_build(marginal)$data
  expect_equal(range(layers[[2L]]$x), range(marginal_layers[[2L]]$x), tolerance = 1e-14)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot_posterior(factor_exp, "mu_f", prior = TRUE,
    n_points = 64L, n_samples = 128L))
  unavailable <- factor_exp
  unavailable$mu_f <- .bt_formula_measure_mark(unavailable$mu_f,
    colnames(unavailable$mu_f)[2L], "prior_density", "The current-scale law is unavailable.",
    cause = "structural_target_law_unavailable")
  expect_s3_class(plot_posterior(unavailable, "mu_f", prior = FALSE,
    plot_type = "ggplot", n_points = 64L), "ggplot")
  expect_error(plot_posterior(unavailable, "mu_f", prior = TRUE,
    plot_type = "ggplot", n_points = 64L), class = "BayesTools_formula_measure_unavailable")
})

test_that("ordered direct displays retain both weighted non-reference priors", {
  fixture <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))
  p <- fixture$prior
  x <- c(-.7, -.2, -.05, 0, .05, .2, .7)
  original <- density(p, x_seq = x, n_points = length(x))
  physical <- vapply(x[x != 0], function(value){
    integral <- stats::integrate(function(share){
      stats::dnorm(value / share, 0, .5) / share
    }, 0, 1, rel.tol = 1e-10, abs.tol = 1e-12)
    expect_lt(integral$abs.error, 1e-10)
    integral$value
  }, numeric(1))
  expect_equal_each(original[[2L]]$y[x != 0], physical, tolerance = 2e-10)
  expect_identical(original[[2L]]$y[x == 0], Inf)
  expect_equal_each(original[[3L]]$y, stats::dnorm(x, 0, .5), tolerance = 1e-14)
  expect_equal(stats::integrate(function(share) .25 * share^2, 0, 1)$value, .25 / 3)
  columns <- .JAGS_prior_factor_names("mu_f", p)
  context <- posterior_metadata(fixture$samples, "prior_context")
  partial <- .prior_density_from_context(context, stats::setNames(c(1, 0), columns))
  complete <- .prior_density_from_context(context, stats::setNames(c(1, 1), columns))
  expect_identical(prior_density_ordinate(partial, 0)$behavior, "infinite")
  expect_identical(prior_density_ordinate(partial, 0)$point_mass, 0)
  expect_true(prior_density_ordinate(partial, 0)$exact)
  expect_identical(prior_density_ordinate(complete, 0)$behavior, "regular")
  displayed <- .plot_data_ordered_prior_display(original)
  expect_identical(displayed[[2L]]$x, x[x != 0])
  expect_identical(displayed[[2L]]$y, original[[2L]]$y[x != 0])
  expect_identical(displayed[[3L]], original[[3L]])
  expect_identical(original, density(p, x_seq = x, n_points = length(x)))
  plots <- plot(p, plot_type = "ggplot", x_seq = x, n_points = length(x))
  expect_length(plots, 3L)
  expect_null(plots[[1L]])
  expect_true(all(vapply(plots[2:3], inherits, logical(1), "ggplot")))
  expect_equal(ggplot2::ggplot_build(plots[[2L]])$data[[1L]]$y, physical, tolerance = 2e-10)
  expect_equal_each(ggplot2::ggplot_build(plots[[3L]])$data[[1L]]$y,
    stats::dnorm(x, 0, .5), tolerance = 1e-14)
  expect_s3_class(plot(p, show_figures = 2L, plot_type = "ggplot", x_seq = x), "ggplot")
  expect_s3_class(plot(p, show_figures = 1L, plot_type = "ggplot", x_seq = x), "ggplot")
  selected <- plot(p, show_figures = -1L, plot_type = "ggplot", x_seq = x)
  expect_length(selected, 3L)
  expect_null(selected[[1L]])
  expect_true(all(vapply(selected[2:3], inherits, logical(1), "ggplot")))
  expect_length(geom_prior(p, x_seq = x), 2L)
  expect_length(geom_prior(p, show_parameter = 2L, x_seq = x), 1L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  original_lines <- graphics::lines.default
  original_arrows <- graphics::arrows
  curves <- arrows <- list()
  testthat::local_mocked_bindings(
    lines.default = function(x, y = NULL, ...){
      curves[[length(curves) + 1L]] <<- list(x = x, y = y)
      original_lines(x, y, ...)
    },
    arrows = function(x0, y0, x1 = x0, y1 = y0, ...){
      arrows[[length(arrows) + 1L]] <<- list(x = x1, y = y1)
      original_arrows(x0, y0, x1, y1, ...)
    }, .package = "graphics")
  expect_no_warning(plot(p, x_seq = x))
  expect_length(curves, 2L)
  expect_length(arrows, 0L)
  expect_identical(curves[[1L]]$x, x[x != 0])
  expect_equal_each(curves[[1L]]$y, physical, tolerance = 2e-10)
  expect_equal_each(curves[[2L]]$y, stats::dnorm(x, 0, .5), tolerance = 1e-14)
  curves <- list()
  expect_no_warning(lines(p, x_seq = x))
  expect_length(curves, 2L)
  expect_length(arrows, 0L)
  expect_equal_each(curves[[1L]]$y, physical, tolerance = 2e-10)
  expect_equal_each(curves[[2L]]$y, stats::dnorm(x, 0, .5), tolerance = 1e-14)
  curves <- list()
  expect_no_warning(plot(p, show_figures = 1L, x_seq = x))
  expect_length(curves, 0L)
  expect_length(arrows, 1L)
  expect_identical(arrows[[1L]], list(x = 0, y = 1))
})

test_that("visible ordered priors preserve genuine atoms and honest density failures", {
  for(total in list(prior("gamma", list(.5, 1)), prior("point", list(2)),
    prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))))){
    p <- ordered_plot_test_fixture(total, prior("dirichlet", list(c(2, 2, .5))))$prior
    original <- density(p, n_points = 64L)
    displayed <- .plot_data_ordered_prior_display(original)
    expect_true(any(displayed[[3L]]$y > 0))
    if(!is.null(original[[3L]]$atoms)){
      expect_identical(displayed[[3L]]$atoms, original[[3L]]$atoms)
      expect_identical(displayed[[3L]]$continuous, original[[3L]]$continuous)
    }
    expect_s3_class(plot(p, show_figures = 3L, plot_type = "ggplot", n_points = 64L), "ggplot")
  }
  p <- ordered_plot_test_fixture(prior("normal", list(0, 1)),
    prior("dirichlet", list(c(.5, .25, .25))))$prior
  analytic <- density(p, x_range = c(-3, 3), n_points = 64L)
  expect_true(all(vapply(analytic[2:4], function(x) all(is.finite(x$y)), logical(1))))
  expect_equal_each(analytic[[4L]]$y, stats::dnorm(analytic[[4L]]$x), tolerance = 0)
  expect_no_warning(plot(p, plot_type = "ggplot", xlim = c(-3, 3), n_points = 64L))
  set.seed(600)
  original <- density(p, force_samples = TRUE, n_points = 64L, n_samples = 128L)
  original_rng <- .Random.seed
  set.seed(600)
  displayed <- .plot_data_ordered_prior_display(
    density(p, force_samples = TRUE, n_points = 64L, n_samples = 128L))
  expect_identical(displayed, original)
  expect_identical(.Random.seed, original_rng)
})

test_that("ordered points-only totals retain exact scalar measures and sampled metadata", {
  cases <- list(
    positive = list(location = 2, probability = 1, alpha = c(2, 2, .5)),
    slab = list(location = 2, probability = .5, alpha = c(2, 2, .5)),
    finite = list(location = 2, probability = .5, alpha = c(2, 2, 2)),
    negative = list(location = -2, probability = 1, alpha = c(2, 2, .5)),
    negative_slab = list(location = -2, probability = .5, alpha = c(2, 2, .5)),
    fixed = list(location = 2, probability = 1, allocation = c(.2, .3, .5)),
    fixed_zero = list(location = 2, probability = .5, allocation = c(0, .3, .7)))
  transformations <- list(identity = list(name = NULL, arguments = NULL),
    exp = list(name = "exp", arguments = NULL),
    tanh = list(name = "tanh", arguments = NULL),
    lin = list(name = "lin", arguments = list(a = 3, b = -2)))
  for(case in cases){
    total <- prior("point", list(case$location))
    if(case$probability < 1) total <- prior_spike_and_slab(total, prior("point", list(case$probability)))
    allocation <- if(is.null(case$alpha)) case$allocation else prior("dirichlet", list(case$alpha))
    p <- ordered_plot_test_fixture(total, allocation)$prior
    for(transformation in transformations){
      forward <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = x, exp = exp(x), tanh = tanh(x), lin = 3 - 2 * x)
      }
      inverse <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = x, exp = log(x), tanh = atanh(x), lin = (x - 3) / -2)
      }
      jacobian <- function(x){
        switch(if(is.null(transformation$name)) "identity" else transformation$name,
          identity = rep(1, length(x)), exp = 1 / x, tanh = 1 / (1 - x^2), lin = rep(.5, length(x)))
      }
      args <- list(x = p, n_points = 64L, n_samples = 128L,
        transformation = transformation$name, transformation_arguments = transformation$arguments)
      set.seed(600)
      analytic <- do.call(density, args)
      set.seed(600)
      forced <- do.call(density, c(args, list(force_samples = TRUE)))
      set.seed(600)
      samples <- rng(p, 128L, transform_factor_samples = TRUE)
      expect_identical(attr(analytic, "method"), "analytic_mixed_measure")
      for(i in seq_len(4L)){
        m <- i - 1L
        fractional <- !is.null(case$alpha) && m > 0L && m < 3L
        share <- if(m == 0L) 0 else if(is.null(case$alpha)) sum(case$allocation[seq_len(m)]) else 1
        expected_atoms <- if(fractional){
          if(case$probability == 1) data.frame(location = numeric(), mass = numeric()) else
            data.frame(location = forward(0), mass = 1 - case$probability)
        }else if(share == 0){
          data.frame(location = forward(0), mass = 1)
        }else if(case$probability == 1){
          data.frame(location = forward(case$location * share), mass = 1)
        }else{
          data.frame(location = forward(c(0, case$location * share)), mass = c(1 - case$probability, case$probability))
        }
        atoms <- analytic[[i]]$atoms
        expect_equal(unname(as.matrix(atoms[order(atoms$location), , drop = FALSE])),
          unname(as.matrix(expected_atoms[order(expected_atoms$location), , drop = FALSE])), tolerance = 2e-12)
        expect_null(analytic[[i]]$samples)
        expect_identical(forced[[i]]$samples, forward(samples[, i]))
        forced[[i]]["samples"] <- list(NULL)
        expect_identical(forced[[i]], analytic[[i]])
        expect_equal(analytic[[i]]$diagnostics$continuous_mass, if(fractional) case$probability else 0)
        if(fractional){
          a <- sum(case$alpha[seq_len(m)])
          b <- sum(case$alpha[(m + 1L):3L])
          curve <- analytic[[i]]$continuous
          source <- inverse(curve$x)
          expected <- case$probability * stats::dbeta(source / case$location, a, b) /
            abs(case$location) * jacobian(curve$x)
          expect_equal_each(curve$density, expected, tolerance = 2e-12)
          expect_equal(case$probability * diff(stats::pbeta(c(0, 1), a, b)),
            analytic[[i]]$diagnostics$continuous_mass, tolerance = 0)
          expect_equal(analytic[[i]]$diagnostics$continuous_integral,
            sum(diff(curve$x) * (head(curve$density, -1L) + tail(curve$density, -1L)) / 2), tolerance = 0)
          expect_true(all(is.finite(curve$x) & is.finite(curve$density)))
          if(is.null(transformation$name) && b >= 1){
            expect_equal(curve$density[c(1L, nrow(curve))],
              case$probability * stats::dbeta(if(case$location > 0) c(0, 1) else c(1, 0), a, b) /
                abs(case$location), tolerance = 2e-12)
          }
        }else{
          expect_null(analytic[[i]]$continuous)
        }
      }
    }
  }
})

test_that("ordered points-only family guard certifies literal mixtures and preserves discrete sampling", {
  total <- prior_mixture(list(prior("point", list(0), prior_weights = .3),
    prior("point", list(2), prior_weights = .7)))
  for(allocation in list(c(.2, .3, .5), prior("dirichlet", list(c(2, 2, .5))))){
    p <- ordered_plot_test_fixture(total, allocation)$prior
    set.seed(600)
    initial_rng <- .Random.seed
    analytic <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L)
    expect_identical(.Random.seed, initial_rng)
    expect_identical(attr(analytic, "method"), "analytic_mixed_measure")
    set.seed(600)
    forced <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = TRUE)
    forced_rng <- .Random.seed
    set.seed(600)
    samples <- rng(p, 128L, transform_factor_samples = TRUE)
    expect_identical(.Random.seed, forced_rng)
    random <- is.prior(allocation)
    for(i in seq_len(4L)){
      partial <- random && i %in% 2:3
      share <- if(i == 1L) 0 else if(random) 1 else sum(allocation[seq_len(i - 1L)])
      expected_atoms <- if(partial) data.frame(location = 0, mass = .3) else if(share == 0)
        data.frame(location = 0, mass = 1) else data.frame(location = c(0, 2 * share), mass = c(.3, .7))
      expect_equal(analytic[[i]]$atoms, expected_atoms, tolerance = 2e-12)
      expect_equal(sum(analytic[[i]]$atoms$mass) + analytic[[i]]$diagnostics$continuous_mass, 1, tolerance = 2e-12)
      expect_null(analytic[[i]]$samples)
      expect_identical(forced[[i]]$samples, samples[, i])
      forced[[i]]["samples"] <- list(NULL)
      expect_identical(forced[[i]], analytic[[i]])
      if(partial){
        curve <- analytic[[i]]$continuous
        shapes <- if(i == 2L) c(2, 2.5) else c(4, .5)
        expect_equal_each(curve$density, .7 * stats::dbeta(curve$x / 2, shapes[1L], shapes[2L]) / 2, tolerance = 2e-12)
      }
    }
  }
  mixed_producer <- .density.prior.ordered_mixed
  unsupported_totals <- list(prior("bernoulli", list(.5)),
    prior("bernoulli", list(.5), truncation = list(lower = .5)),
    prior_mixture(list(prior("bernoulli", list(.5), prior_weights = .7),
      prior("point", list(2), prior_weights = .3))))
  for(total in unsupported_totals){
    for(allocation in list(c(.2, .3, .5), prior("dirichlet", list(c(2, 2, .5))))){
      p <- ordered_plot_test_fixture(total, allocation)$prior
      for(force in c(FALSE, TRUE)){
        set.seed(600)
        actual <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
        actual_rng <- .Random.seed
        expect_null(attr(actual, "method", exact = TRUE))
        expect_identical(unname(vapply(actual, function(component) length(component$samples), integer(1))), rep(128L, 4L))
        local_mocked_bindings(.density.prior.ordered_mixed = function(...) NULL)
        set.seed(600)
        fallback <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
        expect_identical(actual, fallback)
        expect_identical(.Random.seed, actual_rng)
        local_mocked_bindings(.density.prior.ordered_mixed = mixed_producer)
      }
    }
  }
})

test_that("ordered mixed-measure labels describe the actual displayed measure", {
  allocation <- prior("dirichlet", list(c(2, 2, 2)))
  point <- ordered_plot_test_fixture(prior("point", list(2)), allocation)$prior
  slab <- ordered_plot_test_fixture(prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))), allocation)$prior
  expect_identical(plot(point, show_figures = 4L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Probability")
  expect_identical(plot(point, show_figures = 2L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Density")
  expect_identical(plot(slab, show_figures = 2L, plot_type = "ggplot", n_points = 64L)$scales$get_scales("y")$name, "Density / probability mass")
  expect_identical(plot(point, show_figures = 4L, plot_type = "ggplot", n_points = 64L, ylab = "Custom")$scales$get_scales("y")$name, "Custom")
})

test_that("ordered complex points-only totals retain their existing sampled fallback", {
  for(two_ordered in c(FALSE, TRUE)){
    data <- expand.grid(f = ordered(c("early", "middle", "late", "last"),
      levels = c("early", "middle", "late", "last")),
      g = if(two_ordered) ordered(c("low", "mid", "high"), levels = c("low", "mid", "high")) else factor(c("a", "b", "c")))
    build <- function(total){
      JAGS_formula(~f*g, "mu", data, list(intercept = prior("point", list(0)),
        f = prior_ordered(prior("point", list(2))),
        g = if(two_ordered) prior_ordered(prior("point", list(2))) else
          prior_factor("normal", list(0, 1), contrast = "treatment"),
        "f:g" = prior_ordered(total)))$prior_list$mu_f__xXx__g
    }
    p <- build(prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))))
    for(force in c(FALSE, TRUE)){
      set.seed(600)
      sampled <- density(p, x_range = c(0, 2), n_points = 64L, n_samples = 128L, force_samples = force)
      expect_null(attr(sampled, "method", exact = TRUE))
      expect_identical(unname(vapply(sampled, function(component) length(component$samples), integer(1))),
        rep(128L, length(sampled)))
    }
    mixed <- build(prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5))))
    expect_error(density(mixed, x_range = c(-3, 3), n_points = 64L, n_samples = 128L),
      paste0("Mixed-measure ordered densities currently require one ordered term and one scalar total. ",
        "Split the interaction into explicitly named terms before requesting its density."), fixed = TRUE)
  }
})

test_that("centered Normal ordered products agree with independent physical references", {
  # Original-share and endpoint-u integrals at 80/120 decimal working precision,
  # with the binary64 x and sd supplied exactly (directed verification103).
  # Columns: a, b, sd, x, physical reference, maximum reference disagreement.
  references <- do.call(rbind, list(
    c(0.5, 0.5, 0x1.999999999999ap-4, -0x1.15fc9e0f116p-5, 3.82664948847685648811616724671593738978563012867555788981120015191638636345972557371953508314546277420636228737462055993, 5.14950967501997117072376690992153516737195350831454627742063622873746205599299427394583068037783165924111621866535782289e-44),
    c(0.5, 0.5, 0x1p-1, 0x1.5b7bc592d5b78p-3, 0.765329897695372059391787209198367646747723821414882275749266129483775200210612028044434868388801866023925621646590308466, 1.02990193500399443993555123637552477780444348683888018660239256216465903084660036800763463646240873892160296766652246011e-44),
    c(0.5, 0.5, 0x1p+0, -0x1.5b7bc592d5b8p-2, 0.382664948847685654122954439818181314962583216624471974199940069095727325089152950714646418194691677300218123007851995122, 5.14950967501997142366241040253664346071464641819469167730021812300785199512203533753017052006199134422343303391822699493e-45),
    c(0.5, 0.5, 0x1p+0, -0x1.8p+1, 0.000804763347585775120878231132982508083782070633154401095347239288941597495508371156887066077560779082122815929537469833784, 6.05962114626367121098912740822456849468870660775607790821228159295374698337840270046982605605675419500221128173398236068e-47),
    c(0.5, 0.5, 0x1p+0, -0x1.899c0f60188p-9, 4.9884551753737355577465243717037154796608897590451287382032763021649513772643276289561286310224175969038900721922615266, 5.45467252829873728990866204223036952895612863102241759690389007219226152659974166527545910420273851905318227545993641256e-45),
    c(0.5, 0.5, 0x1.8p+1, 0x1.049cd42e204ap+0, 0.127554982949228551374318146606060438320861072208157324733313356365242441696384316904882139398230559100072707669283998374, 1.71650322500665714122080346751221448690488213939823055910007270766928399837399564337707830119452703961350565150671078906e-45),
    c(0.75, 0.25, 0x1.999999999999ap-4, -0x1.25fa2848925fcp-4, 2.64913791375915945389854934289773344822513386471215108354787882261097032486641564289721860567627857087165474032506337973, 0.000000000000000000000413400781750339390547573523043728340494451987848949114046042897218605676278570871654740325063379730064351350652960976058),
    c(0.75, 0.25, 0x1p-1, -0x1.6f78b25ab6f78p-2, 0.52982758275183211324801622097263115186591539301802857857589691922101969546620321232852270228795921404878192350502598507, 0.0000000000000000000000826801563500679000994877996023430657650724037701117240945523285227022879592140487819235050259850699215944349939686836828),
    c(0.75, 0.25, 0x1p+0, -0x1.6f78b25ab6f78p-1, 0.264913791375916056624008110486315575932957696509014289287948459610509847733101606164261351143979607024390961752512992535, 0.0000000000000000000000413400781750339500497438998011715328825362018850558620472761642613511439796070243909617525129925349607972174969843418414),
    c(0.75, 0.25, 0x1p+0, -0x1.8p+1, 0.00205720662386305884592253044389522750488491311836516837085318197874629835317433137743352575562789887510430950070450671816, 0.000000000000000000000000594159074336491242562095975003080574788183977695529517880777433525755627898875104309500704506718159949182820754910149454),
    c(0.75, 0.25, 0x1p+0, -0x1.899c0f60188p-9, 1.57443399441549927674075184694688379309054025636855464969742186369460085297896066469536031479946230696878792198488849421, 0.0000000000000000000000534842542461753350594465961670203812764597843399403131423646953603147994623069687879219848884942099867135750824378819547),
    c(0.75, 0.25, 0x1.8p+1, 0x1.139a85c40939cp+1, 0.0883045971253053197208723369176916362245971104263617251445791515689554495834758272076494025675097392341590340060272663507, 0.0000000000000000000000137800260583446470885031192184324034205150426539805853049212076494025675097392341590340060272663506958180515271566003618),
    c(1, 1, 0x1.999999999999ap-4, -0x1.2c20988612c2p-6, 7.03569956540695001996521759448343176368067730432539912291621969796392187576405827622694485996292287947545985466632604344, 2.37730551400370771205245401453336739565603096370471388852106711200749690005087371367221332628445030128928227054236911137e-80),
    c(1, 1, 0x1p-1, 0x1.7728bea79772p-4, 1.40713991308139098961078317919234647397016750904456297508610174936043692640797696323861283093554161137174768396094657198, 3.67613871690644583886282523160390534280200226130511275954183938979032511717754917016640084619102127219646703783702084092e-80),
    c(1, 1, 0x1p+0, -0x1.7728bea79773p-3, 0.703569956540694543745042931265904320097821918817647992969189706925588230666630754106768275554324372275546684734587690971, 4.10676827555432437227554668473458769097103259174227720019175375243614748233790350241489318232591632047276233723895468387e-81),
    c(1, 1, 0x1p+0, -0x1.8p+1, 0.000413583612635949924292039355179629910338655050758448307225229430909085175493709932221090543671376742578569913702751560149, 2.22109054367137674257856991370275156014899754927901763986932593263997437695374042373621775052512879768829305013880909065e-84),
    c(1, 1, 0x1p+0, -0x1.899c0f60188p-9, 2.34023950087991183006503142744383603056289674883078964269434646595775640641334577075701848890509423260172032445450585338, 2.92429815110949057673982796755454941466202021977139202581242031605289210604705513850342760837698994059448882586251700432e-80),
    c(1, 1, 0x1.8p+1, -0x1.195e8efdb196p-1, 0.234523318846898286921719716903020837656889315571108244772443114528238370041455195312043765026944801019527299640205636062, 4.68795623497305519898047270035979436393800229230043878951774922350247509800819692713360198334766425333490590066134748165e-81),
    c(1.5, 1.5, 0x1.999999999999ap-4, 0x1.b4a9210e9b4a8p-4, 0.937310402969806386980469169213031390332904262904398970215686637217385379335384683942473565118500375143576920156080393624, 3.94247356511850037514357692015608039362402516025521109918279241462663396328008029011110738886773012776210351986543665946e-81),
    c(1.5, 1.5, 0x1p-1, 0x1.10e9b4a9210e8p-1, 0.187462080593961362331899879064136978022099287977439271683453302673746169034673008647999496034647814629715625793967904612, 1.3520005039653521853702843742060320953880090670117517825638504317071743845431615710420932356374882759096225826683388953e-81),
    c(1.5, 1.5, 0x1p+0, -0x1.8p+1, 0.000248719254429097984865515793499500575237053788927667589834399955662890521859038914550789150352468342839706148081164024322, 5.44921084964753165716029385191883597567795810370837888632145816563307351484075450275295986492520633641122955772711724165e-84),
    c(100, 100, 0x1p+0, -0x1.8p+1, 0.0000000844410411096973926758485025623859281330427528705225673275710127871995313771315334773545521909228388290936620909690160086, 1.23560062511308633213364988352605954355162772142802391079495640495821158399827059181259242014943745861377246933776231088e-72),
    c(0.5, 0.5, 0x1p+0, -0x1.7f3b31f84ff3bp+1, 0.000820924597906134914515244382924453233718075159622709465063248029347575049872951279484180264395089332560391391707488511683, 6.16968179079763948869190559072999001094841802643950893325603913917074885116829503209020966484515778028830172335271668568e-47),
    c(1, 1, 0x1p-1, -0x1.dc08bb712893ap+0, 0.0000506787494603685233392044161821033698430515442360507324241141053574773586958238996727629524881711959674572755370748168135, 3.27237047511828804032542724462925183186500318958186312565643605764122106924644452318472868612798488639803728962635542431e-85)
  ))
  for(i in seq_len(nrow(references))){
    ref <- references[i, ]
    total <- prior("normal", list(0, ref[3L]))
    value <- .density.prior.ordered_dirichlet_product(total,
      ref[1L], ref[2L], ref[4L], 1L)
    error <- attr(value, "quadrature")[[1L]]$abs.error
    # Eight relative ulps cover conversion of the reference to binary64 and
    # arithmetic in the comparison; there is no absolute tolerance floor.
    allowance <- ref[6L] + 8 * .Machine$double.eps * abs(ref[5L])
    expect_lte(abs(value$y - ref[5L]), error + allowance)
    expect_identical(attr(value, "quadrature_tolerance"),
      c(relative = 1e-7, absolute = 1e-10))
    zero <- .density.prior.ordered_dirichlet_product(total,
      ref[1L], ref[2L], 0, 1L)
    if(ref[1L] <= 1){
      expect_identical(zero$y, Inf)
    }else{
      expect_equal(zero$y, stats::dnorm(0, 0, ref[3L]) *
        (ref[1L] + ref[2L] - 1) / (ref[1L] - 1), tolerance = 0)
    }
  }
})

test_that("ordered product integration failures return no partial curve", {
  calls <- 0L
  testthat::local_mocked_bindings(integrate = function(...){
    calls <<- calls + 1L
    list(value = 1, abs.error = 0, message = if(calls == 1L) "OK" else "injected failure")
  }, .package = "stats")
  result <- NULL
  expect_error(result <- .density.prior.ordered_dirichlet_product(
    prior("normal", list(0, 1)), .75, .25, c(.1, .2, .3), 3L),
    "Ordered-prior product density quadrature failed at x = 0.20000000000000001: injected failure; absolute error 0.",
    fixed = TRUE)
  expect_identical(calls, 2L)
  expect_null(result)
})

test_that("ordered overlays map and clip only selected atom probabilities", {
  ordinary <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))$prior
  point <- ordered_plot_test_fixture(prior("point", list(2)))$prior
  mixed <- ordered_plot_test_fixture(prior_spike_and_slab(
    prior("point", list(-2)), prior("point", list(.5))),
    prior("dirichlet", list(c(2, 2, .5))))$prior
  priors <- list(prior("normal", list(0, .5), prior_weights = .5),
    prior("point", list(0), prior_weights = .5))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_prior_list(priors, ylim = c(0, 2), ylim2 = c(0, .6))
  expect_no_warning(lines(ordinary, x_seq = c(-.2, .2)))
  expect_warning(lines(ordinary, show_parameter = 1L, x_seq = c(-.2, .2)),
    "Point-mass probabilities outside the active secondary-axis limits will be clipped.", fixed = TRUE)
  expect_warning(lines(point, show_parameter = 4L),
    "Point-mass probabilities outside the active secondary-axis limits will be clipped.", fixed = TRUE)
  original_arrows <- graphics::arrows
  captured <- list()
  testthat::local_mocked_bindings(arrows = function(x0, y0, x1 = x0, y1 = y0, ...){
    captured[[length(captured) + 1L]] <<- list(x = x1, y = y1)
    original_arrows(x0, y0, x1, y1, ...)
  }, .package = "graphics")
  for(transformation in list(list(name = "exp", arguments = NULL, location = 1),
    list(name = "tanh", arguments = NULL, location = 0),
    list(name = "lin", arguments = list(a = 3, b = -2), location = 3))){
    captured <- list()
    state <- .plot_scale_y2_state_current()
    expect_no_warning(lines(mixed, show_parameter = 3L,
      transformation = transformation$name, transformation_arguments = transformation$arguments))
    expect_length(captured, 1L)
    expect_equal(captured[[1L]], list(x = transformation$location, y = .5 * state$scale_y2))
    base <- plot_prior_list(priors, plot_type = "ggplot", ylim = c(0, 2), ylim2 = c(0, .6))
    expect_no_warning(base + geom_prior(ordinary, x_seq = c(-.2, .2)))
    expect_warning(base + geom_prior(ordinary, show_parameter = 1L, x_seq = c(-.2, .2)),
      "Point-mass probabilities outside the active secondary-axis limits will be clipped.", fixed = TRUE)
    expect_warning(base + geom_prior(point, show_parameter = 4L),
      "Point-mass probabilities outside the active secondary-axis limits will be clipped.", fixed = TRUE)
    expect_no_warning(overlay <- base + geom_prior(mixed, show_parameter = 3L,
      transformation = transformation$name, transformation_arguments = transformation$arguments))
    data <- ggplot2::ggplot_build(overlay)$data
    expect_equal(tail(data, 1L)[[1L]]$x, transformation$location)
    expect_equal(tail(data, 1L)[[1L]]$yend, .5 * base$bt_scale_y2_state$scale_y2)
  }
})

test_that("ordered display with no finite ordinate gives a controlled remedy", {
  p <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))$prior
  unavailable <- density(p, x_seq = c(-1, 0, 1))
  unavailable[[2L]]$x <- 0
  unavailable[[2L]]$y <- Inf
  testthat::local_mocked_bindings(density.prior = function(...) unavailable)
  expect_error(plot(p, show_figures = 2L, x_seq = 0, plot_type = "ggplot"),
    "The ordered prior density curve is unavailable: the requested 'x_seq' contains no finite plotting ordinate. Supply 'x_seq' with finite density ordinates.",
    fixed = TRUE, class = "BayesTools_prior_curve_unavailable")
  expect_no_warning(plot(p, show_figures = 1L, x_seq = 0, plot_type = "ggplot"))
})

test_that("ordered overlay warnings retain source coordinates and available measures", {
  capture <- function(f){
    warnings <- list()
    value <- withCallingHandlers(f(), warning = function(w){
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
  }
  fixture <- ordered_plot_test_fixture(prior_spike_and_slab(
    prior("normal", list(0, .5)), prior("point", list(.5))),
    levels = c("systematic", "alternate", "random"))
  samples <- as_mixed_posteriors(fixture$fit, "mu_f")
  marginal <- marginal_posterior(samples, "mu_f", use_formula = FALSE, prior_samples = TRUE)
  original_density <- .prior_density_route_density
  original_ordinate <- .prior_density_route_ordinate
  original_plot_data <- .prior_linear_density_to_plot_data
  injected <- FALSE
  refusal_mode <- "single"
  refused_values <- c(.125, -.375)
  records <- list()
  zero_heights <- numeric()
  refused <- list()
  testthat::local_mocked_bindings(
    .prior_density_route_density = function(route, x, batch_singular = FALSE){
      y <- original_density(route, x, batch_singular)
      if(identical(route$type, "conditional_normal")){
        zero_heights <<- c(zero_heights, y[x == 0])
      }
      continuous_leaf <- !route$type %in% c("mixture", "atom")
      at_refusal <- is.finite(x) & is.finite(y) & y > 0 & if(refusal_mode == "plural"){
        continuous_leaf & x %in% refused_values
      }else if(refusal_mode == "whole"){
        identical(route$type, "conditional_normal")
      }else{
        identical(route$type, "conditional_normal") & x == .125
      }
      if(injected && any(at_refusal)){
        original_results <- lapply(x[at_refusal], function(value) original_ordinate(route, value))
        expect_true(all(vapply(original_results, function(result){
          identical(result$behavior, "regular") && result$exact && is.finite(result$log_density)
        }, logical(1))))
        unavailable <- vapply(original_results, function(result){
          result$log_density <- NA_real_
          result$exact <- FALSE
          result$reason <- "Injected numerical quadrature refusal."
          result$provenance$integration <- list(converged = FALSE)
          refused[[length(refused) + 1L]] <<- result
          .prior_density_ordinate_height_value(result)
        }, numeric(1))
        y[at_refusal] <- unavailable
      }
      y
    },
    .prior_linear_density_to_plot_data = function(...){
      args <- list(...)
      if(refusal_mode == "plural"){
        route <- .prior_density_route_from_adaptive(attr(args[[1L]], "adaptive_evaluation", exact = TRUE))
        # Quadrature merges one display value; 130/129 give binary 2^-6 steps
        # placing .125 and -.375 exactly in both real baseline/injected curves.
        args$n_points <- if(.prior_density_route_has_quadrature(route)) 130L else 129L
      }
      value <- do.call(original_plot_data, args)
      records[[length(records) + 1L]] <<- value
      value
    }
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for(transformation in list(list(name = NULL, arguments = NULL, refused = .125),
    list(name = "lin", arguments = list(a = 2, b = 2), refused = 2.25))){
    for(backend in c("base", "ggplot")){
      for(marginal_plot in c(FALSE, TRUE)){
        draw <- function(){
          if(marginal_plot) plot_marginal(list(mu_f = marginal), "mu_f", prior = TRUE,
            plot_type = backend, xlim = c(-.875, 1.125),
            transformation = transformation$name, transformation_arguments = transformation$arguments) else
            plot_posterior(samples, "mu_f", prior = TRUE, plot_type = backend,
              xlim = c(-.875, 1.125), transformation = transformation$name,
              transformation_arguments = transformation$arguments)
        }
        injected <- FALSE
        records <- list()
        baseline <- capture(draw)
        expect_length(baseline$warnings, 0L)
        expected <- records
        expect_length(expected, 2L)
        partial <- which(vapply(expected, function(record){
          any(record$density$x == transformation$refused)
        }, logical(1)))
        expect_identical(partial, 1L)
        expect_length(expected[[partial]]$points1$x, 1L)
        expect_identical(expected[[partial]]$points1$y, .5)
        baseline_curve <- expected[[partial]]$density
        available <- expected[[partial]]$density$x != transformation$refused
        expected[[partial]]$density$x <- expected[[partial]]$density$x[available]
        expected[[partial]]$density$y <- expected[[partial]]$density$y[available]
        injected <- TRUE
        records <- list()
        output <- capture(draw)
        expect_length(output$warnings, 1L)
        warning <- output$warnings[[1L]]
        expect_s3_class(warning, "BayesTools_prior_curve_unavailable")
        expect_s3_class(warning, "BayesTools_plot_condition")
        expect_identical(warning$unresolved_values, .125)
        expect_identical(conditionMessage(warning), paste0(
          "The prior density curve is partially unavailable: numerical evaluations were unresolved at 1 evaluation coordinate on the source scale. ",
          "Available curve points and declared atoms are retained; use 'prior = FALSE' to draw the posterior alone."))
        expect_identical(records, expected)
        if(!is.null(transformation$name)){
          expect_false(identical(warning$unresolved_values, transformation$refused))
        }
        if(backend == "ggplot"){
          expected_render <- ggplot2::ggplot_build(baseline$value)$data
          curve_layer <- which(vapply(expected_render, function(layer){
            identical(layer$x, baseline_curve$x) && identical(layer$y, baseline_curve$y)
          }, logical(1)))
          expect_length(curve_layer, 1L)
          expected_render[[curve_layer]] <- expected_render[[curve_layer]][available, , drop = FALSE]
          rownames(expected_render[[curve_layer]]) <- NULL
          expect_identical(ggplot2::ggplot_build(output$value)$data, expected_render)
        }
      }
    }
  }
  for(refusal_mode in c("whole", "plural")){
    for(backend in c("base", "ggplot")){
      for(marginal_plot in c(FALSE, TRUE)){
        draw <- function(){
          if(marginal_plot) plot_marginal(list(mu_f = marginal), "mu_f", prior = TRUE,
            plot_type = backend, xlim = c(-.875, 1.125)) else
            plot_posterior(samples, "mu_f", prior = TRUE, plot_type = backend,
              xlim = c(-.875, 1.125))
        }
        injected <- FALSE
        records <- list()
        baseline <- capture(draw)
        expect_length(baseline$warnings, 0L)
        expected <- records
        expect_length(expected, 2L)
        baseline_curves <- lapply(expected, `[[`, "density")
        expect_true(all(vapply(baseline_curves, function(curve){
          all(is.finite(curve$x)) && all(is.finite(curve$y)) && all(curve$y > 0)
        }, logical(1))))
        if(refusal_mode == "whole"){
          expected_values <- unique(baseline_curves[[1L]]$x)
          expect_false(0 %in% expected_values)
          expected[[1L]]$density <- NULL
        }else{
          expect_true(all(vapply(baseline_curves, function(curve){
            all(refused_values %in% curve$x)
          }, logical(1))))
          available <- lapply(baseline_curves, function(curve) !curve$x %in% refused_values)
          for(i in seq_along(expected)){
            expected[[i]]$density$x <- baseline_curves[[i]]$x[available[[i]]]
            expected[[i]]$density$y <- baseline_curves[[i]]$y[available[[i]]]
          }
          expected_values <- refused_values
        }
        injected <- TRUE
        records <- list()
        output <- capture(draw)
        expect_length(output$warnings, if(refusal_mode == "whole") 1L else 2L)
        for(warning in output$warnings){
          expect_s3_class(warning, "BayesTools_prior_curve_unavailable")
          expect_s3_class(warning, "BayesTools_plot_condition")
          expect_type(warning$unresolved_values, "double")
          expect_true(all(is.finite(warning$unresolved_values)))
          expect_identical(sort(warning$unresolved_values), sort(expected_values))
          expect_identical(conditionMessage(warning), paste0(
            "The prior density curve is ",
            if(refusal_mode == "whole") "unavailable" else "partially unavailable",
            ": numerical evaluations were unresolved at ", length(expected_values),
            " evaluation coordinates on the source scale. ",
            if(refusal_mode == "whole") "No continuous curve points are available; declared atoms are retained; " else
              "Available curve points and declared atoms are retained; ",
            "use 'prior = FALSE' to draw the posterior alone."))
        }
        expect_identical(records, expected)
        expect_identical(lapply(records, `[[`, "points1"), lapply(expected, `[[`, "points1"))
        if(refusal_mode == "whole"){
          expect_null(records[[1L]]$density)
          expect_identical(records[[2L]]$density, baseline_curves[[2L]])
          if(backend == "ggplot") expect_error(ggplot2::ggplot_build(output$value), NA)
        }else if(backend == "ggplot"){
          expected_render <- ggplot2::ggplot_build(baseline$value)$data
          for(i in seq_along(baseline_curves)){
            curve_layer <- which(vapply(expected_render, function(layer){
              identical(layer$x, baseline_curves[[i]]$x) && identical(layer$y, baseline_curves[[i]]$y)
            }, logical(1)))
            expect_length(curve_layer, 1L)
            expected_render[[curve_layer]] <- expected_render[[curve_layer]][available[[i]], , drop = FALSE]
            rownames(expected_render[[curve_layer]]) <- NULL
          }
          expect_identical(ggplot2::ggplot_build(output$value)$data, expected_render)
        }
      }
    }
  }
  expect_true(length(zero_heights) > 0L)
  expect_true(all(zero_heights == Inf))
  expect_true(length(refused) > 0L)
  expect_true(all(vapply(refused, function(result){
    identical(result$behavior, "regular") && !result$exact && is.na(result$log_density)
  }, logical(1))))
  ordinary <- ordered_plot_test_fixture(prior("normal", list(0, .5)),
    levels = c("systematic", "alternate", "random"))
  ordinary$samples <- as_mixed_posteriors(ordinary$fit, "mu_f")
  injected <- FALSE
  refusal_mode <- "single"
  expect_no_warning(plot_posterior(ordinary$samples, "mu_f", prior = TRUE, plot_type = "ggplot"))
})

test_that("ordered direct displays of a total away from zero follow the structural route", {
  # density of total * share, share ~ Beta(a1, a2): int f_T(x / c) dbeta(c) / c dc over the band of c
  # where the total is supported, in 200 pieces
  reference <- function(total_pdf, lower, upper, a1, a2, x){
    lo <- max(0, x / upper)
    hi <- min(1, x / lower)
    if(hi <= lo) return(0)
    knots <- seq(lo, hi, length.out = 201L)
    sum(vapply(seq_len(200L), function(i){
      stats::integrate(function(c) total_pdf(x / c) * stats::dbeta(c, a1, a2) / c,
        knots[i], knots[i + 1L], rel.tol = 1e-10, abs.tol = 0)$value
    }, numeric(1)))
  }
  truncated_normal <- prior("normal", list(2, .1), list(1.5, 2.5))
  cases <- list(
    list(total = prior("uniform", list(1, 2)), pdf = function(t) stats::dunif(t, 1, 2), lower = 1, upper = 2),
    list(total = prior("normal", list(1, .01)), pdf = function(t) stats::dnorm(t, 1, .01), lower = .88, upper = 1.12),
    list(total = truncated_normal,
         pdf = function(t) stats::dnorm(t, 2, .1) / diff(stats::pnorm(c(1.5, 2.5), 2, .1)), lower = 1.5, upper = 2.5)
  )
  x <- c(.005, .05, .15, .5, .9)
  for(case in cases){
    p <- ordered_plot_test_fixture(case$total, levels = c("a", "b", "c", "d"))$prior
    metadata <- .prior_ordered_metadata(.prior_ordered_default_bound(p))
    design <- .factor_term_design_from_metadata(p)$design
    colnames(design) <- .JAGS_prior_factor_names(metadata$parameter_name, p)
    displayed <- density(p, x_seq = x)
    expect_identical(attr(displayed, "method"), "direct")
    for(level in 2:3){
      # cumulative shares of three equal Dirichlet increments: Beta(level - 1, 4 - level)
      expected <- vapply(x, function(value){
        reference(case$pdf, case$lower, case$upper, level - 1, 4 - level, value)
      }, numeric(1))
      structural <- .prior_linear_combination_density(
        stats::setNames(list(p), metadata$parameter_name), design[level, design[level, ] != 0], n_grid = 512L)
      ordinates <- vapply(x, function(value) exp(prior_density_ordinate(structural, value)$log_density), numeric(1))
      expect_identical(displayed[[level]]$x, x)
      expect_equal(displayed[[level]]$y, expected, tolerance = 1e-6)
      expect_equal(displayed[[level]]$y, ordinates, tolerance = 1e-6)
    }
  }

  # closed forms for the uniform(1, 2) total: 2 (log 2 - x / 2) and x on (0, 1)
  p <- ordered_plot_test_fixture(cases[[1L]]$total, levels = c("a", "b", "c", "d"))$prior
  displayed <- density(p, x_seq = x)
  expect_equal(displayed[[2L]]$y, 2 * (log(2) - x / 2), tolerance = 1e-8)
  expect_equal(displayed[[3L]]$y, x, tolerance = 1e-8)
  # the default grid has no zero height inside the support and each level keeps its unit mass
  default <- density(p)
  for(level in 2:3){
    inside <- default[[level]]$x > 0 & default[[level]]$x < 2
    expect_true(all(default[[level]]$y[inside] > 0))
    grid <- default[[level]]$x
    height <- default[[level]]$y
    expect_equal(sum(diff(grid) * (head(height, -1L) + tail(height, -1L)) / 2), 1, tolerance = 1e-4)
  }

  # a Gamma(2, 1) total with a Dirichlet(1, 1) share: the intermediate level is Exponential(1);
  # at 0 the height is the one-sided limit 1
  gamma <- ordered_plot_test_fixture(prior("gamma", list(2, 1)), prior("dirichlet", list(c(1, 1))),
    levels = c("a", "b", "c"))$prior
  at <- c(0, .001, .1, 1)
  intermediate <- density(gamma, x_seq = at, truncate_end = FALSE)[[2L]]
  expect_identical(intermediate$x, at)
  expect_equal(intermediate$y, stats::dexp(at), tolerance = 1e-6)
})

test_that("ordered direct displays keep the Dirichlet-product quadrature and route only the values it leaves open", {
  route_density <- .prior_density_route_density
  levels <- c("a", "b", "c", "d")
  # independent reference: the density of total * share, share ~ Beta(a1, a2)
  reference <- function(total_pdf, a1, a2, x){
    stats::integrate(function(c) total_pdf(x / c) * stats::dbeta(c, a1, a2) / c,
      0, 1, rel.tol = 1e-10, abs.tol = 0, subdivisions = 1000L)$value
  }

  # a lognormal total: the quadrature resolves every positive display value, so the structural route
  # is only asked for the support end 0 (a height below the absolute tolerance), and no numerical
  # grid of a level is built
  p <- ordered_plot_test_fixture(prior("lognormal", list(0, .5)), prior("dirichlet", list(c(2, 1, .5))),
    levels = levels)$prior
  evaluated <- numeric()
  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...) stop("the numerical grid of a level was built", call. = FALSE),
    .prior_density_route_density = function(route, x, ...){
      evaluated <<- c(evaluated, x)
      route_density(route, x, ...)
    })
  x <- c(0, .3, 1, 2)
  displayed <- density(p, x_seq = x)
  expect_identical(attr(displayed, "method"), "direct")
  expect_true(length(evaluated) > 0L)
  expect_true(all(evaluated == 0))
  for(level in 2:3){
    # cumulative shares of the Dirichlet(2, 1, .5) increments: Beta(sum(first m), sum(rest))
    a <- c(2, 3)[level - 1L]
    expect_identical(displayed[[level]]$x, x)
    expect_equal(displayed[[level]]$y[-1L],
      vapply(x[-1L], function(value) reference(function(t) stats::dlnorm(t, 0, .5), a, 3.5 - a, value), numeric(1)),
      tolerance = 1e-6)
    expect_equal(displayed[[level]]$y[1L], 0)
  }
  # the default display is the same quadrature: no value of the level is left to the route but 0
  evaluated <- numeric()
  default <- density(p)
  expect_true(all(evaluated == 0))
  for(level in 2:3){
    expect_true(all(is.finite(default[[level]]$y)))
    expect_length(default[[level]]$x, 1000L)
  }
})

test_that("the structural route takes the display values that the quadrature misses", {
  route_density <- .prior_density_route_density
  # uniform(1, 2) total: the integrand over logit(c) is zero outside a band that the quadrature
  # may not visit; its height is then 0 (below the absolute tolerance), not a resolved value
  p <- ordered_plot_test_fixture(prior("uniform", list(1, 2)), levels = c("a", "b", "c", "d"))$prior
  x <- c(.005, .05, .15, .5, .9)
  evaluated <- numeric()
  testthat::local_mocked_bindings(
    .prior_density_route_density = function(route, x, ...){
      evaluated <<- c(evaluated, x)
      route_density(route, x, ...)
    })
  displayed <- density(p, x_seq = x)
  expect_true(length(evaluated) > 0L)
  expect_true(all(evaluated %in% c(.005, .05)))
  expect_equal(displayed[[2L]]$y, 2 * (log(2) - x / 2), tolerance = 1e-8)
  expect_equal(displayed[[3L]]$y, x, tolerance = 1e-8)

  # the quadrature on its own leaves a missed value open (NA) for a total other than the
  # certified centered Normal
  product <- .density.prior.ordered_dirichlet_product(prior("uniform", list(1, 2)), 1, 2, x, 5L)
  expect_true(is.na(product$y[1L]))
  expect_identical(attr(product, "quadrature")[[1L]]$message, "value not above the absolute tolerance")
  expect_equal(product$y[4L], 2 * (log(2) - x[4L] / 2), tolerance = 1e-6)
})

test_that("a quadrature that does not converge is not repeated along the level of an uncertified total", {
  calls <- 0L
  testthat::local_mocked_bindings(integrate = function(...){
    calls <<- calls + 1L
    list(value = 1, abs.error = 0, message = "maximum number of subdivisions reached")
  }, .package = "stats")
  level <- .density.prior.ordered_dirichlet_product(prior("lognormal", list(0, 1)), 1, 2, c(.1, .2, .3, .4), 4L)
  expect_identical(calls, 1L)
  expect_true(all(is.na(level$y)))
  messages <- vapply(attr(level, "quadrature"), `[[`, character(1), "message")
  expect_identical(messages[1L], "maximum number of subdivisions reached; absolute error 0")
  expect_true(all(startsWith(messages[-1L], "not attempted")))
  # the certified centered Normal total keeps its refusal of an unresolved value
  expect_error(.density.prior.ordered_dirichlet_product(prior("normal", list(0, 1)), 1, 2, c(.1, .2), 2L),
    "Ordered-prior product density quadrature failed at x = 0.10000000000000001: maximum number of subdivisions reached; absolute error 0.",
    fixed = TRUE)
  expect_identical(calls, 2L)
})

test_that("unresolved ordered display values are reported, not dropped silently", {
  levels <- c("a", "b", "c", "d")
  # a lognormal total with a Dirichlet(1, 1, .05) allocation: the quadrature of the third level does not
  # converge and its structural route resolves 1 of 200 display values (0.3.1.185 stopped with the
  # quadrature failure); the other levels are complete
  p <- ordered_plot_test_fixture(prior("lognormal", list(0, 1)), prior("dirichlet", list(c(1, 1, .05))),
    levels = levels)$prior
  condition <- NULL
  displayed <- withCallingHandlers(density(p, n_points = 200L),
    warning = function(w){ condition <<- w; invokeRestart("muffleWarning") })
  expect_s3_class(condition, "BayesTools_prior_curve_unavailable")
  expect_s3_class(condition, "BayesTools_plot_condition")
  unresolved <- condition$unresolved_values
  expect_gt(length(unresolved), 100L)
  expect_identical(condition$levels, names(displayed)[3L])
  expect_match(conditionMessage(condition), "partially unavailable", fixed = TRUE)
  expect_match(conditionMessage(condition), "'x_seq'", fixed = TRUE)
  expect_length(displayed[[2L]]$x, 200L)
  expect_false(any(unresolved %in% displayed[[3L]]$x))
  expect_identical(sort(c(displayed[[3L]]$x, unresolved)), displayed[[2L]]$x)
  expect_true(all(is.finite(displayed[[3L]]$y)))
  expect_warning(plot(p, n_points = 200L, plot_type = "ggplot"), class = "BayesTools_prior_curve_unavailable")
  # the default display of 1000 values resolves every value of the level and does not warn
  default <- NULL
  expect_no_warning(default <- density(p))
  expect_length(default[[3L]]$x, 1000L)

  # a display with every value resolved does not warn
  resolved <- ordered_plot_test_fixture(prior("lognormal", list(0, .5)), prior("dirichlet", list(c(2, 1, .5))),
    levels = levels)$prior
  expect_no_warning(density(resolved, n_points = 200L))
})

test_that("unresolved values of a mixed-measure ordered level are reported", {
  route_density <- .prior_density_route_density
  # a spike-and-slab total has atoms, so its levels take the mixed-measure display
  p <- ordered_plot_test_fixture(prior_spike_and_slab(prior("gamma", list(2, 1)), prior("point", list(.5))),
    prior("dirichlet", list(c(2, 2, .5))))$prior
  x <- seq(.2, 4, length.out = 20)
  testthat::local_mocked_bindings(
    .prior_density_route_density = function(route, x, ...){
      out <- route_density(route, x, ...)
      out[x > 3] <- NA_real_
      out
    })
  condition <- NULL
  displayed <- withCallingHandlers(density(p, x_seq = x, n_points = 20L),
    warning = function(w){ condition <<- w; invokeRestart("muffleWarning") })
  expect_s3_class(condition, "BayesTools_prior_curve_unavailable")
  # the total (last level) is the one with a structural route here
  expect_identical(condition$levels, names(displayed)[4L])
  expect_identical(condition$unresolved_values, x[x > 3])
  expect_identical(displayed[[4L]]$continuous$x, x[x <= 3])
})

test_that("ordered mixed-measure plots honour ylim2 and ylab2 and store the probability mapping", {

  allocation <- prior("dirichlet", list(c(2, 2)))
  mixed <- ordered_plot_test_fixture(prior_spike_and_slab(prior("point", list(2)), prior("point", list(.5))),
    allocation, levels = c("a", "b", "c"))$prior
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  scale_y2 <- .plot.prior_settings()[["scale_y2"]] * .1 / 1

  # without ylim2 the masses keep their unit scale and no mapping is stored
  default <- plot(mixed, show_figures = 2L, plot_type = "ggplot", n_points = 32L)
  default_layers <- ggplot2::ggplot_build(default)$data
  expect_null(default$bt_scale_y2_state)
  expect_identical(default$scales$get_scales("y")$name, "Density / probability mass")
  expect_equal(default_layers[[2L]]$yend, .5)

  # ggplot: the secondary axis spans the requested limits and carries the label
  figure <- plot(mixed, show_figures = 2L, plot_type = "ggplot", n_points = 32L,
    ylim = c(0, .1), ylim2 = c(0, 1), ylab2 = "P")
  layers <- ggplot2::ggplot_build(figure)$data
  y_scale <- figure$scales$get_scales("y")
  expect_identical(y_scale$secondary.axis$name, "P")
  expect_equal(range(y_scale$secondary.axis$labels), c(0, 1))
  expect_equal(y_scale$limits, c(0, scale_y2))
  expect_equal(figure$bt_scale_y2_state$scale_y2, scale_y2)
  expect_equal(figure$bt_scale_y2_state$ylim2, c(0, 1))
  expect_identical(y_scale$name, "Density")
  # the atom of mass .5 is not clipped; the density heights and masses are unchanged
  expect_equal(layers[[2L]]$yend, .5 * scale_y2)
  expect_lt(layers[[2L]]$yend, y_scale$limits[2L])
  expect_equal(layers[[1L]]$y, default_layers[[1L]]$y)
  expect_equal(layers[[1L]]$x, default_layers[[1L]]$x)
  # the default label of the secondary axis
  expect_identical(plot(mixed, show_figures = 2L, plot_type = "ggplot", n_points = 32L,
    ylim = c(0, .1), ylim2 = c(0, 1))$scales$get_scales("y")$secondary.axis$name, "Probability")

  # the stored mapping is used by an overlay (the atom is scaled like the plot)
  overlay <- figure + geom_prior(mixed, show_parameter = 2L, n_points = 32L)
  expect_equal(tail(ggplot2::ggplot_build(overlay)$data, 1L)[[1L]]$yend, .5 * scale_y2)

  # base graphics: the probability mapping is remembered for overlays
  expect_no_error(plot(mixed, show_figures = 2L, n_points = 32L, ylim = c(0, .1), ylim2 = c(0, 1), ylab2 = "P"))
  state <- .plot_scale_y2_state_current()
  expect_false(is.null(state))
  expect_equal(state$scale_y2, scale_y2)
  expect_equal(state$ylim2, c(0, 1))
  expect_gte(graphics::par("usr")[4L], .5 * scale_y2)
  captured <- list()
  original_arrows <- graphics::arrows
  testthat::local_mocked_bindings(arrows = function(x0, y0, x1 = x0, y1 = y0, ...){
    captured[[length(captured) + 1L]] <<- list(x = x1, y = y1)
    original_arrows(x0, y0, x1, y1, ...)
  }, .package = "graphics")
  expect_no_warning(lines(mixed, show_parameter = 2L, n_points = 32L))
  expect_equal(captured[[1L]]$y, .5 * scale_y2)
})

test_that("a constant transformation turns every ordered level into one atom", {

  data <- data.frame(f = ordered(c("a", "b", "c"), levels = c("a", "b", "c")))
  build <- function(total, allocation = prior("dirichlet", list(c(1, 1)))){
    JAGS_formula(~f, "mu", data, list(intercept = prior("point", list(0)),
      f = prior_ordered(total, allocation = allocation)))$prior_list$mu_f
  }
  # a 50% Normal slab and a 50% spike at zero: the non-reference levels carry
  # atoms and continuous mass on the structural route
  mixed <- build(prior_spike_and_slab(prior("normal", list(0, 1)), prior("point", list(.5))))

  cases <- list(
    list(transformation = "lin",     arguments = list(a = 3, b = 0), x_seq = c(-1, 0, 1),  location = 3),
    list(transformation = "exp_lin", arguments = list(a = 1, b = 0), x_seq = c(.5, 1, 2), location = exp(1))
  )
  for(case in cases){
    constant <- density(mixed, x_seq = case$x_seq, n_points = 16, transformation = case$transformation,
                        transformation_arguments = case$arguments)
    expect_length(constant, 3L)
    for(level in constant){
      expect_identical(level$atoms$location, case$location)
      expect_equal(level$atoms$mass, 1)
      expect_null(level$continuous)
      expect_equal(level$x, case$location)
      expect_equal(level$y, 1)
      expect_equal(level$diagnostics$atom_mass, 1)
      expect_equal(level$diagnostics$continuous_mass, 0)
    }
  }

  # a non-constant map keeps the atom and the continuous part of the levels
  shifted <- density(mixed, x_seq = c(-1, 0, 1), n_points = 16, transformation = "lin",
                     transformation_arguments = list(a = 3, b = 1))
  for(level in shifted[-1L]){
    expect_equal(level$atoms$location, 3)
    expect_equal(level$atoms$mass, .5)
    expect_false(is.null(level$continuous))
    expect_equal(level$diagnostics$continuous_mass, .5)
  }
})

test_that("the atoms of an ordered total follow from its structure", {

  point <- prior("point", list(1))
  normal <- prior("normal", list(0, 1))
  expect_true(.density.prior_ordered_total_has_atoms(point))
  expect_true(.density.prior_ordered_total_has_atoms(prior_none()))
  expect_true(.density.prior_ordered_total_has_atoms(prior("bernoulli", list(.5))))
  expect_true(.density.prior_ordered_total_has_atoms(prior_spike_and_slab(normal)))
  expect_true(.density.prior_ordered_total_has_atoms(prior_mixture(list(normal, point))))
  expect_false(.density.prior_ordered_total_has_atoms(normal))
  expect_false(.density.prior_ordered_total_has_atoms(prior("gamma", list(2, 1))))
  expect_false(.density.prior_ordered_total_has_atoms(prior_mixture(list(normal, prior("gamma", list(2, 1))))))
})

test_that("an ordered prior with a continuous total is displayed on a narrow range without building its total grid", {

  data <- data.frame(f = ordered(rep(c("a", "b", "c", "d"), 3L), levels = c("a", "b", "c", "d")))
  p <- JAGS_formula(~f, "mu", data, list(intercept = prior("normal", list(0, 1)),
    f = prior_ordered(prior("normal", list(0, 1)), allocation = prior("dirichlet", list(c(5, 5, 5))))))$prior_list$mu_f

  # the heights at points shared with a wide display range (the display
  # points -2e-6, 0, 2e-6 are points 500 to 502 of the wide grid)
  wide   <- density(p, x_range = c(-1e-3, 1e-3), n_points = 1001L)
  narrow <- density(p, x_range = c(-2e-6, 2e-6), n_points = 3L)
  expect_identical(names(narrow), names(wide))
  for(level in 2:4){
    expect_equal(narrow[[level]]$y, wide[[level]]$y[500:502], tolerance = 1e-10)
  }

  # a range far narrower than the display spacing of the total returns the
  # level heights, finite and positive
  extreme <- density(p, x_range = c(-1e-6, 1e-6))
  for(level in 2:4){
    expect_true(all(is.finite(extreme[[level]]$y)) && all(extreme[[level]]$y > 0))
  }
})
