# Run with an installed BayesTools and JAGS 5:
# Rscript tools/verify-allocation-relationships.R
library(BayesTools)

set.seed(789)
d <- expand.grid(x_unit = c(-1, 1), f = factor(c("a", "b", "c")), replicate = 1:8)
d$x <- 4 + 3 * d$x_unit
d$study <- factor(rep(1:8, each = 6))
d$paper <- factor(rep(1:6, 8))
d$lab <- factor(rep(1:4, 12))
d$sei <- rep(c(0.8, 1, 1.2), 16)
d$y <- 0.10 * (d$x - mean(d$x)) + rnorm(nrow(d), sd = d$sei)

# One external SD and one realized R² govern all moderator and random branches.
# The exact law on these orthogonal term blocks allocates empirical variances.
half <- prior("spike", list(0.5))
budget <- prior_R2(
  terms = c("x", "f"), allocation = R2_dirichlet(c(x = 3, f = 2)),
  total = prior_spike_and_slab(prior("beta", list(2, 3)), prior_inclusion = half),
  model = R2_model(law = "exact", sigma = "random_sd", contrast = "orthonormal")
)
random <- prior_random(
  random_variance_allocation("random_tree", terms = c(cluster = "intermediate", lab = "lab"),
    sd_source = random_sd_source("random_sd"), weights = prior("dirichlet", list(c(2, 3)))),
  random_variance_allocation("cluster_tree", terms = c(study = "study", paper = "paper"),
    parent = allocation_ref("random_tree", "cluster"), weights = prior("dirichlet", list(c(3, 2))))
)
fit <- JAGS_fit(
  model_syntax = "model{
    random_sd <- tau_total * sqrt(1 - fixed_Signal)
    for(i in 1:N) { y[i] ~ dnorm(fixed[i] + mu[i], pow(sei[i], -2)) }
  }",
  data = list(y = d$y, sei = d$sei, N = nrow(d)),
  prior_list = list(tau_total = prior_spike_and_slab(prior("gamma", list(2, 4)), prior_inclusion = half)),
  formula_list = list(fixed = ~ x + f, mu = ~ 1 + random(1 | study, name = "study") +
    random(1 | paper, name = "paper") + random(1 | lab, name = "lab")),
  formula_data_list = list(fixed = d, mu = d), formula_scale_list = list(fixed = TRUE),
  formula_prior_list = list(fixed = list(intercept = prior("spike", list(0)), Signal = budget),
                            mu = list(intercept = prior("normal", list(0, 1)))),
  formula_random_prior_list = list(mu = random), add_parameters = "random_sd",
  chains = 2, adapt = 150, burnin = 200, sample = 400, seed = 789, silent = TRUE
)
if (inherits(fit, "error")) stop(conditionMessage(fit))
catalog <- parameter_catalog(fit)
select <- function(alias, namespace = NULL) parameter_catalog_resolve(catalog, alias, namespace = namespace)
declarations <- data.frame(
  source_id = c(select("tau_total")$quantity_id, select("Signal{total}", "fixed")$quantity_id),
  consumer_id = select("random_sd")$quantity_id,
  relation = c("sd_source", "complement_variance_fraction")
)
native <- parameter_relationships(fit)
registered <- register_parameter_relationships(fit, declarations, "RoBMA")
restored <- unserialize(serialize(registered, NULL))
relationships <- parameter_relationships(restored)
stopifnot(identical(native, parameter_relationships(fit)),
          identical(relationships, parameter_relationships(registered)),
          sum(relationships$provider == "RoBMA") == 2L)
stopifnot(any(relationships$relation == "parent_component" & relationships$allocation == "cluster_tree"),
          any(relationships$relation == "budget_total" & relationships$source_id == declarations$source_id[2L]))
draw <- function(alias, namespace = NULL) as.numeric(as.matrix(parameter_draws(restored, select(alias, namespace))))
r <- draw("Signal{total}", "fixed")
share_x <- draw("Signal[share: x]", "fixed")
tau <- draw("tau_total")
residual_sd <- draw("random_sd")
on <- r > 0
stopifnot(all(r < 1), any(!on & tau > 0), any(on & tau == 0), any(on & tau > 0),
          all(is.na(share_x[!on])), all(is.finite(share_x[on])))

# Read only terminal variance quantities. Aggregate, gate and share rows
# describe other aspects of the same tree and would count branches twice.
leaves <- relationships[relationships$relation == "block_var" & relationships$formula_parameter == "mu", ]
q <- catalog$quantities
indices <- match(leaves$consumer_id, q$quantity_id)
stopifnot(nrow(leaves) == 3L, all(leaves$consumer_type == "quantity"), !anyNA(indices),
          !anyDuplicated(indices), all(q$quantity[indices] == "var"))
random_variances <- vapply(indices, function(i) draw(q$canonical_name[i], q$namespace[i]), numeric(length(tau)))

# Recover the outcome-scale slope and factor levels independently from predictions.
grid <- data.frame(x = c(3, 5), f = factor(c("a", "a"), levels = levels(d$f)))
predict_fixed <- function(grid) JAGS_predict_formula(restored, "fixed", data = grid, formula_target = "fixed")$value
p <- predict_fixed(grid)
slope <- (p[2, ] - p[1, ]) / diff(grid$x)
coefficient_map <- JAGS_formula_coefficient_transform(restored, "fixed")$matrix
fitted_slope <- draw("x", "fixed")
original_slope <- unname(coefficient_map["fixed_x", "fixed_x"]) * fitted_slope
levels_f <- levels(d$f)
grid <- data.frame(x = mean(d$x), f = factor(levels_f, levels = levels_f))
centred <- function(m) sweep(m, 2, colMeans(m))
effects <- centred(predict_fixed(grid))
cells <- t(vapply(levels_f, function(f) draw(paste0("f[", f, "]"), "fixed"), numeric(length(tau))))
v_x <- slope^2 * mean((d$x - mean(d$x))^2)
v_f <- colMeans(effects^2)
expected_x <- expected_f <- numeric(length(tau))
expected_x[on] <- tau[on]^2 * r[on] * share_x[on]
expected_f[on] <- tau[on]^2 * r[on] * (1 - share_x[on])
errors <- c(
  slope = max(abs(slope - residual_sd * original_slope)),
  cells = max(abs(effects - sweep(centred(cells), 2, residual_sd, `*`))),
  moderator_x = max(abs(v_x - expected_x)), moderator_f = max(abs(v_f - expected_f)),
  random = max(abs(rowSums(random_variances) - tau^2 * (1 - r))),
  total = max(abs(v_x + v_f + rowSums(random_variances) - tau^2))
)
stopifnot(all(errors < 1e-10))
# r = 0 may come from the budget spike. A point total at r = 1 cannot use
# this residual-SD parameterization because its r / (1 - r) is undefined.
boundary <- tryCatch(prior_R2(total = prior("spike", list(1))), error = identity)
stopifnot(inherits(boundary, "BayesTools_R2_invalid"))
cat(length(tau), "draws; R² off and tau positive:", sum(!on & tau > 0),
    "; R² on and tau zero:", sum(on & tau == 0), "\n")
print(signif(errors, 3))

# A row source with a declared scalar input uses the new rows' exposure.
d$exposure <- rep(c(0.5, 1, 1.5), length.out = nrow(d))
source <- parameter_source("row_sd", shape = "row", inputs = "tau_total",
  values = function(parameters, data, n_rows) parameters$tau_total * data$exposure[seq_len(n_rows)])
pool <- random_variance_allocation("row_pool", terms = c(a = "study", b = "lab"),
  sd_source = random_sd_source(source), weights = prior("dirichlet", list(c(2, 3))))
row_fit <- JAGS_fit(
  model_syntax = "model{ for(i in 1:N) {
    row_sd[i] <- tau_total * exposure[i]
    y[i] ~ dnorm(mu[i], 1)
  } }", data = list(y = d$y, N = nrow(d), exposure = d$exposure),
  prior_list = list(tau_total = prior("gamma", list(2, 4))),
  formula_list = list(mu = ~ 1 + random(1 | study, name = "study") + random(1 | lab, name = "lab")),
  formula_data_list = list(mu = d), formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
  formula_random_prior_list = list(mu = prior_random(pool)),
  chains = 1, adapt = 100, burnin = 100, sample = 200, seed = 787, silent = TRUE
)
if (inherits(row_fit, "error")) stop(conditionMessage(row_fit))
e <- parameter_relationships(row_fit)
input <- parameter_catalog_resolve(parameter_catalog(row_fit), "tau_total")$quantity_id
stopifnot(any(e$relation == "source_input" & e$source_type == "quantity" & e$source_id == input &
              e$consumer_type == "row_source" & e$consumer_id == "row_sd"),
          any(e$relation == "sd_source" & e$source_type == "row_source" & e$source_id == "row_sd" &
              e$consumer_type == "allocation" & e$allocation == "row_pool"))
original <- d[c(2, 7, 11, 19), , drop = FALSE]
changed <- transform(original, exposure = 2 * exposure)
predict_rows <- function(grid, target) JAGS_predict_formula(row_fit, "mu", data = grid, formula_target = target)$value
fixed <- predict_rows(original, "fixed")
before <- predict_rows(original, "conditional") - fixed
after <- predict_rows(changed, "conditional") - fixed
row_error <- max(abs(after - 2 * before))
stopifnot(max(abs(before)) > 0.01, row_error < 1e-10)
cat("New-row random-effect rescaling error:", format(row_error, digits = 3), "\n")
