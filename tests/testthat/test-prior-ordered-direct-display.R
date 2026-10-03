skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

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
  records <- list()
  zero_heights <- numeric()
  refused <- list()
  testthat::local_mocked_bindings(
    .prior_density_route_density = function(route, x, batch_singular = FALSE){
      y <- original_density(route, x, batch_singular)
      if(identical(route$type, "conditional_normal")){
        zero_heights <<- c(zero_heights, y[x == 0])
        at_refusal <- x == .125
        if(injected && any(at_refusal)){
          result <- original_ordinate(route, .125)
          expect_identical(result$behavior, "regular")
          expect_true(result$exact)
          result$log_density <- NA_real_
          result$exact <- FALSE
          result$reason <- "Injected numerical quadrature refusal."
          result$provenance$integration <- list(converged = FALSE)
          refused[[length(refused) + 1L]] <<- result
          y[at_refusal] <- .prior_density_ordinate_height_value(result)
        }
      }
      y
    },
    .prior_linear_density_to_plot_data = function(...){
      value <- original_plot_data(...)
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
  expect_no_warning(plot_posterior(ordinary$samples, "mu_f", prior = TRUE, plot_type = "ggplot"))
})
