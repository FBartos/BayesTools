# Mockup: deterministic original-scale prior densities after formula scaling.
#
# Run from the package root:
#   Rscript dev/prior-unscale-density-mockup.R
#
# The prototype targets the scalar marginal density of each original-scale
# coefficient. If beta_orig[j] = sum_i M[j, i] * beta_z[i] and the standardized
# priors are independent, the density is the convolution of the scaled marginal
# priors M[j, i] * beta_z[i]. This avoids Monte Carlo bumps in transformed prior
# overlays while still allowing mixtures, spikes, truncation, and non-normal
# priors with a usable mpdf()/range() implementation.

load_bayestools_here <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }

  for (file in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    sys.source(file, envir = .GlobalEnv)
  }
  invisible(TRUE)
}

bt_get <- function(name) {
  if ("BayesTools" %in% loadedNamespaces() &&
      exists(name, envir = asNamespace("BayesTools"), inherits = FALSE)) {
    get(name, envir = asNamespace("BayesTools"), inherits = FALSE)
  } else {
    get(name, envir = .GlobalEnv, inherits = FALSE)
  }
}

fft_convolve <- function(a, b) {
  n <- length(a) + length(b) - 1L
  n_fft <- 2^ceiling(log2(n))
  out <- fft(fft(c(a, rep(0, n_fft - length(a)))) *
               fft(c(b, rep(0, n_fft - length(b)))),
             inverse = TRUE) / n_fft
  pmax(Re(out[seq_len(n)]), 0)
}

aggregate_points <- function(points, dx) {
  if (is.null(points) || nrow(points) == 0) {
    return(data.frame(x = numeric(), p = numeric()))
  }

  key <- round(points$x / dx)
  data.frame(
    x = as.numeric(names(tapply(points$x, key, mean))),
    p = as.numeric(tapply(points$p, key, sum)),
    row.names = NULL
  ) |>
    transform(x = x * dx)
}

coalesce_distribution <- function(dist, dx) {
  points <- aggregate_points(dist$points, dx)
  densities <- dist$densities

  if (length(densities) == 0) {
    return(list(points = points, density = NULL))
  }

  x_min <- min(vapply(densities, function(d) min(d$x), numeric(1)))
  x_max <- max(vapply(densities, function(d) max(d$x), numeric(1)))
  x <- seq(x_min, x_max, by = dx)
  y <- numeric(length(x))
  mass <- 0

  for (d in densities) {
    y <- y + d$mass * stats::approx(d$x, d$y, xout = x, yleft = 0, yright = 0)$y
    mass <- mass + d$mass
  }

  area <- sum(y) * dx
  if (is.finite(area) && area > 0 && mass > 0) {
    y <- y / area
  }

  list(points = points, density = list(x = x, y = y, mass = mass))
}

raw_density_components <- function(prior, scale, mass = 1, tail_prob = 1e-4) {
  is_prior_mixture <- bt_get("is.prior.mixture")
  is_prior_spike_and_slab <- bt_get("is.prior.spike_and_slab")
  is_prior_point <- bt_get("is.prior.point")
  is_prior_simple <- bt_get("is.prior.simple")
  prior_mean <- bt_get("mean.prior")

  if (abs(scale) < sqrt(.Machine$double.eps)) {
    return(list(list(type = "point", x = 0, mass = mass)))
  }

  if (is_prior_spike_and_slab(prior)) {
    variable_prior <- bt_get(".get_spike_and_slab_variable")(prior)
    inclusion_prior <- bt_get(".get_spike_and_slab_inclusion")(prior)
    inclusion <- as.numeric(prior_mean(inclusion_prior))
    inclusion <- min(max(inclusion, 0), 1)

    return(c(
      raw_density_components(variable_prior, scale, mass * inclusion, tail_prob),
      list(list(type = "point", x = 0, mass = mass * (1 - inclusion)))
    ))
  }

  if (is_prior_mixture(prior)) {
    weights <- attr(prior, "prior_weights")
    weights <- weights / sum(weights)
    return(unlist(
      Map(function(component, weight) {
        raw_density_components(component, scale, mass * weight, tail_prob)
      }, prior, weights),
      recursive = FALSE
    ))
  }

  if (is_prior_point(prior)) {
    return(list(list(
      type = "point",
      x = scale * prior$parameters[["location"]],
      mass = mass
    )))
  }

  if (!is_prior_simple(prior)) {
    stop("This mockup handles simple, mixture, and spike-and-slab priors.")
  }

  prior_range <- range(prior, quantiles = tail_prob)
  transformed_range <- sort(scale * prior_range)
  list(list(
    type = "continuous",
    prior = prior,
    scale = scale,
    range = transformed_range,
    mass = mass
  ))
}

source_distribution <- function(components, dx) {
  points <- data.frame(x = numeric(), p = numeric())
  densities <- list()
  mpdf_prior <- bt_get("mpdf")

  for (component in components) {
    if (component$type == "point") {
      points <- rbind(points, data.frame(x = component$x, p = component$mass))
      next
    }

    x <- seq(component$range[1], component$range[2], by = dx)
    if (length(x) < 3) {
      x <- seq(component$range[1], component$range[2], length.out = 3)
    }

    y <- mpdf_prior(component$prior, x / component$scale) / abs(component$scale)
    y[!is.finite(y)] <- 0
    area <- sum(y) * (x[2] - x[1])
    if (area <= 0) {
      stop("A continuous prior component evaluated to zero mass on its grid.")
    }
    y <- y / area
    densities[[length(densities) + 1L]] <- list(
      x = x,
      y = y,
      mass = component$mass
    )
  }

  coalesce_distribution(list(points = points, densities = densities), dx)
}

shift_density <- function(density, shift, mass) {
  list(x = density$x + shift, y = density$y, mass = density$mass * mass)
}

convolve_distributions <- function(lhs, rhs, dx) {
  densities <- list()
  points <- data.frame(x = numeric(), p = numeric())

  if (nrow(lhs$points) > 0 && nrow(rhs$points) > 0) {
    grid <- merge(lhs$points, rhs$points, by = NULL)
    points <- rbind(points, data.frame(
      x = grid$x.x + grid$x.y,
      p = grid$p.x * grid$p.y
    ))
  }

  if (!is.null(lhs$density) && !is.null(rhs$density)) {
    y <- fft_convolve(lhs$density$y, rhs$density$y) * dx
    area <- sum(y) * dx
    if (area > 0) {
      y <- y / area
    }
    x <- lhs$density$x[1] + rhs$density$x[1] + dx * (seq_along(y) - 1)
    densities[[length(densities) + 1L]] <- list(
      x = x,
      y = y,
      mass = lhs$density$mass * rhs$density$mass
    )
  }

  if (!is.null(lhs$density) && nrow(rhs$points) > 0) {
    for (i in seq_len(nrow(rhs$points))) {
      densities[[length(densities) + 1L]] <- shift_density(
        lhs$density, rhs$points$x[i], rhs$points$p[i]
      )
    }
  }

  if (!is.null(rhs$density) && nrow(lhs$points) > 0) {
    for (i in seq_len(nrow(lhs$points))) {
      densities[[length(densities) + 1L]] <- shift_density(
        rhs$density, lhs$points$x[i], lhs$points$p[i]
      )
    }
  }

  coalesce_distribution(list(points = points, densities = densities), dx)
}

linear_combination_density <- function(priors, weights, n_grid = 2^12, tail_prob = 1e-4) {
  keep <- abs(weights) > sqrt(.Machine$double.eps)
  priors <- priors[keep]
  weights <- weights[keep]

  if (length(weights) == 0) {
    return(list(x = 0, y = NA_real_, points = data.frame(x = 0, p = 1)))
  }

  raw_components <- Map(raw_density_components, priors, weights,
                        MoreArgs = list(tail_prob = tail_prob))
  source_ranges <- vapply(raw_components, function(components) {
    x_min <- 0
    x_max <- 0
    for (component in components) {
      if (component$type == "point") {
        x_min <- min(x_min, component$x)
        x_max <- max(x_max, component$x)
      } else {
        x_min <- min(x_min, component$range[1])
        x_max <- max(x_max, component$range[2])
      }
    }
    c(min = x_min, max = x_max)
  }, numeric(2))

  target_width <- diff(range(c(sum(source_ranges["min", ]), sum(source_ranges["max", ]))))
  dx <- target_width / (n_grid - 1)
  if (!is.finite(dx) || dx <= 0) {
    dx <- 1
  }

  dist <- list(points = data.frame(x = 0, p = 1), density = NULL)
  for (components in raw_components) {
    source_dist <- source_distribution(components, dx)
    dist <- convolve_distributions(dist, source_dist, dx)
  }

  x_range <- range(
    if (!is.null(dist$density)) dist$density$x else numeric(),
    dist$points$x,
    finite = TRUE
  )
  x <- seq(x_range[1], x_range[2], length.out = n_grid)
  y <- rep(0, length(x))
  if (!is.null(dist$density) && dist$density$mass > 0) {
    y <- dist$density$mass *
      stats::approx(dist$density$x, dist$density$y, xout = x, yleft = 0, yright = 0)$y
  }

  list(x = x, y = y, points = dist$points, continuous_mass = sum(y) * (x[2] - x[1]))
}

load_bayestools_here()

set.seed(1)
df <- data.frame(
  x1 = rnorm(80, mean = 10, sd = 2.5),
  x2 = rgamma(80, shape = 2, rate = 0.4)
)

formula_result <- JAGS_formula(
  formula = ~ x1 * x2,
  parameter = "mu",
  data = df,
  formula_scale = list(x1 = TRUE, x2 = TRUE),
  prior_list = list(
    intercept = prior_mixture(list(
      prior("normal", list(0, 1), prior_weights = 3),
      prior("t", list(0, 1.5, 3), prior_weights = 1)
    )),
    x1 = prior("t", list(location = 0, scale = 0.7, df = 3)),
    x2 = prior_spike_and_slab(
      prior("normal", list(mean = 0.3, sd = 0.5)),
      prior_inclusion = prior("beta", list(alpha = 2, beta = 3))
    ),
    "x1:x2" = prior_mixture(list(
      prior("normal", list(0, 0.4), prior_weights = 2),
      prior("gamma", list(shape = 2, rate = 4), prior_weights = 1)
    ))
  )
)

prior_list <- formula_result$prior_list
formula_scale <- list(mu = formula_result$formula_scale)
term_names <- names(prior_list)
M <- bt_get(".build_unscale_matrix")(term_names, formula_result$formula_scale, "mu")

deterministic <- lapply(term_names, function(term) {
  linear_combination_density(
    priors = prior_list,
    weights = M[term, ],
    n_grid = 2^12,
    tail_prob = 1e-4
  )
})
names(deterministic) <- term_names

set.seed(11)
prior_samples <- bt_get(".generate_transformed_prior_samples")(
  prior_list = prior_list,
  column_names = term_names,
  n_samples = 200000,
  formula_scale = formula_scale
)

comparison <- do.call(rbind, lapply(term_names, function(term) {
  sample_density <- stats::density(prior_samples[, term], n = 2^12)
  deterministic_y <- stats::approx(
    deterministic[[term]]$x,
    deterministic[[term]]$y,
    xout = sample_density$x,
    yleft = 0,
    yright = 0
  )$y
  data.frame(
    term = term,
    sample_mean = mean(prior_samples[, term]),
    sample_sd = stats::sd(prior_samples[, term]),
    l1_density_gap = mean(abs(sample_density$y - deterministic_y)),
    continuous_mass = deterministic[[term]]$continuous_mass,
    point_mass = sum(deterministic[[term]]$points$p)
  )
}))

comparison_print <- comparison
numeric_columns <- vapply(comparison_print, is.numeric, logical(1))
comparison_print[numeric_columns] <- lapply(comparison_print[numeric_columns], round, digits = 5)
print(comparison_print)

if (interactive() || identical(Sys.getenv("BAYESTOOLS_MOCKUP_PLOT"), "true")) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  for (term in term_names) {
    plot(deterministic[[term]]$x, deterministic[[term]]$y, type = "l",
         main = term, xlab = "original scale", ylab = "density")
    lines(stats::density(prior_samples[, term], n = 2^12), col = "grey50")
    if (nrow(deterministic[[term]]$points) > 0) {
      rug(deterministic[[term]]$points$x, col = "red")
    }
  }
}

invisible(list(
  formula_scale = formula_scale,
  unscale_matrix = M,
  deterministic = deterministic,
  prior_samples = prior_samples,
  comparison = comparison
))
