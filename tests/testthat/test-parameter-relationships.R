skip_if_not_test_profile("unit")

test_that("relationships join declared public IDs and survive catalog and map round trips", {

  names <- c("tau_total", "random_sd", "moderator_sd", "rho")
  priors <- stats::setNames(lapply(names, function(x) prior("normal", list(0, 1))), names)
  map <- .bt_build_parameter_map(names, prior_list = priors)
  catalog <- .bt_parameter_map_catalog(map)
  expect_identical(parameter_relationships(map), parameter_relationships(catalog))
  expect_identical(names(parameter_relationships(map)), parameter_map_schema()$relationships$field)
  id <- function(x) parameter_catalog_resolve(catalog, x)$quantity_id
  declarations <- data.frame(
    source_id = c(id("tau_total"), id("rho"), id("tau_total"), id("rho")),
    consumer_id = c(rep(id("random_sd"), 2), rep(id("moderator_sd"), 2)),
    relation = c("sd_source", "complement_variance_fraction", "sd_source", "variance_fraction"),
    component = rep(c("random", "moderators"), each = 2)
  )
  updated <- register_parameter_relationships(map, declarations, "RoBMA")
  expect_identical(nrow(parameter_relationships(map)), 0L)
  expect_identical(parameter_relationships(updated)$provider, rep("RoBMA", 4))
  expect_identical(parameter_relationships(unserialize(serialize(updated, NULL))), parameter_relationships(updated))
  expect_identical(register_parameter_relationships(updated, declarations, "RoBMA"), updated)
  expect_identical(parameter_relationships(register_parameter_relationships(catalog, declarations, "RoBMA")),
                   parameter_relationships(updated))
  expect_identical(.bt_parameter_map_catalog(updated)$relationships, updated$relationships)
  expect_identical(updated$quantities, map$quantities)
  expect_error(register_parameter_relationships(updated, declarations, "Other"),
               class = "BayesTools_parameter_relationship_error")
  cycle <- data.frame(source_id = id("random_sd"), consumer_id = id("tau_total"), relation = "model_input")
  expect_error(register_parameter_relationships(updated, cycle, "RoBMA"), "cycle",
               class = "BayesTools_parameter_relationship_error")
  duplicate_root <- transform(cycle, source_id = id("rho"), consumer_id = id("random_sd"), relation = "sd_source")
  expect_error(register_parameter_relationships(updated, duplicate_root, "RoBMA"), "one declared SD source",
               class = "BayesTools_parameter_relationship_error")
  unknown <- transform(cycle, source_id = "absent")
  expect_error(register_parameter_relationships(updated, unknown, "RoBMA"),
               class = "BayesTools_parameter_relationship_error")
  expect_error(register_parameter_relationships(map, declarations, "BayesTools"),
               class = "BayesTools_parameter_relationship_error")
})

test_that("a nested allocation links the consumed component and retains the external source", {

  dat <- data.frame(study = factor(rep(1:4, each = 6)), paper = factor(rep(1:3, 8)),
                    lab = factor(rep(1:2, 12)), y = rep(c(-1, 1, 0), 8))
  random <- prior_random(
    random_variance_allocation("pool", terms = c(branch = "nested", other = "lab"),
      sd_source = random_sd_source("tau_total"), weights = prior("dirichlet", list(c(2, 3))),
      display_name = "Root"),
    random_variance_allocation("split", parent = allocation_ref("pool", "branch"),
      terms = c(first = "study", second = "paper"), weights = prior("dirichlet", list(c(3, 2))),
      display_name = "Studies", component_names = c("Study effect", "Paper effect"))
  )
  fit <- JAGS_fit(model_syntax = "model{for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) }}",
    data = list(y = dat$y, N = nrow(dat)), prior_list = list(tau_total = prior("gamma", list(2, 3))),
    formula_list = list(mu = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | lab, name = "lab", covariance = "diag")),
    formula_data_list = list(mu = dat), formula_random_prior_list = list(mu = random),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
    chains = 1, sample = 100, burnin = 100, adapt = 100, seed = 24, silent = TRUE)
  catalog <- parameter_catalog(fit)
  edges <- parameter_relationships(fit)
  root <- parameter_catalog_resolve(catalog, "tau_total")
  scale <- edges[edges$relation == "sd_source" & edges$allocation == "pool", ]
  expect_identical(scale$source_id, root$quantity_id)
  branches <- edges[edges$source_type == "allocation" & edges$source_id == scale$consumer_id &
                      edges$relation == "allocates", ]
  child <- subset(edges, relation == "parent_component" & allocation == "split")
  expect_identical(child$source_id, branches$consumer_id[branches$component == "branch"])
  expect_identical(child$component, "branch")
  leaves <- subset(edges, relation == "block_sd")
  expect_setequal(leaves$component, c("other", "first", "second"))
  ids <- catalog$quantities$quantity_id
  variances <- vapply(leaves$consumer_id, function(id){
    selection <- parameter_catalog_resolve(catalog, catalog$quantities$canonical_name[match(id, ids)])
    as.numeric(as.matrix(parameter_draws(fit, selection)))^2
  }, numeric(100))
  expect_equal(rowSums(variances), as.numeric(as.matrix(parameter_draws(fit, root)))^2,
               tolerance = 1e-12)
  expect_identical(parameter_relationships(parameter_catalog(fit)), edges)
  expect_identical(parameter_relationships(unserialize(serialize(fit, NULL))), edges)
})
