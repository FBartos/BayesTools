skip_if_not_test_profile("unit")

test_that("relationships join declared public IDs and survive catalog and map round trips", {

  priors <- list(tau_total = prior("gamma", list(2, 3)), rho = prior("beta", list(2, 2)))
  map <- .bt_build_parameter_map(c(names(priors), "random_sd", "moderator_sd"), prior_list = priors)
  catalog <- .bt_parameter_map_catalog(map)
  expect_identical(parameter_relationships(map), parameter_relationships(catalog))
  expect_identical(names(parameter_relationships(map)), parameter_map_schema()$relationships$field)
  id <- function(x) parameter_catalog_resolve(catalog, x)$quantity_id
  native <- parameter_relationships(map)
  expect_identical(native$source_type, rep("prior", 2))
  expect_identical(native$source_id, names(priors))
  expect_identical(native$relation, rep("prior_owns", 2))
  expect_identical(native$consumer_id, c(id("tau_total"), id("rho")))
  declarations <- data.frame(
    source_id = c(id("tau_total"), id("rho"), id("tau_total"), id("rho")),
    consumer_id = c(rep(id("random_sd"), 2), rep(id("moderator_sd"), 2)),
    relation = c("sd_source", "complement_variance_fraction", "sd_source", "variance_fraction"),
    component = rep(c("random", "moderators"), each = 2)
  )
  updated <- register_parameter_relationships(map, declarations, "RoBMA")
  expect_identical(parameter_relationships(map), native)
  expect_identical(parameter_relationships(updated)$provider, c(native$provider, rep("RoBMA", 4)))
  expect_identical(parameter_relationships(unserialize(serialize(updated, NULL))), parameter_relationships(updated))
  expect_identical(register_parameter_relationships(updated, declarations, "RoBMA"), updated)
  expect_identical(parameter_relationships(register_parameter_relationships(catalog, declarations, "RoBMA")),
                   parameter_relationships(updated))
  expect_identical(.bt_parameter_map_catalog(updated)$relationships, updated$relationships)
  expect_identical(updated$quantities, map$quantities)
  expect_error(register_parameter_relationships(updated, declarations, "Other"),
               class = "BayesTools_parameter_relationship_error")
  supplied <- data.frame(source_id = id("random_sd"), consumer_id = id("moderator_sd"), relation = "model_input")
  linked <- register_parameter_relationships(updated, supplied, "RoBMA")
  cycle <- data.frame(source_id = id("moderator_sd"), consumer_id = id("random_sd"), relation = "model_input")
  expect_error(register_parameter_relationships(linked, cycle, "RoBMA"), "cycle",
               class = "BayesTools_parameter_relationship_error")
  duplicate_root <- transform(cycle, source_id = id("rho"), consumer_id = id("random_sd"), relation = "sd_source")
  expect_error(register_parameter_relationships(updated, duplicate_root, "RoBMA"), "one declared SD source",
               class = "BayesTools_parameter_relationship_error")
  unknown <- transform(cycle, source_id = "absent")
  expect_error(register_parameter_relationships(updated, unknown, "RoBMA"),
               class = "BayesTools_parameter_relationship_error")
  expect_error(register_parameter_relationships(map, declarations, "BayesTools"),
               class = "BayesTools_parameter_relationship_error")
  owned <- data.frame(source_id = id("rho"), consumer_id = id("tau_total"), relation = "model_input")
  for(object in list(map, catalog, unserialize(serialize(updated, NULL)))){
    expect_error(register_parameter_relationships(object, owned, "RoBMA"), "fitted prior",
                 class = "BayesTools_parameter_relationship_error")
  }
  forged <- transform(owned, source_type = "prior", source_id = "tau_total",
                      consumer_id = id("random_sd"), relation = "prior_owns")
  expect_error(register_parameter_relationships(map, forged, "RoBMA"), "native prior",
               class = "BayesTools_parameter_relationship_error")
})

test_that("an indexed ordinary coordinate keeps the exact fitted prior-list owner", {

  priors <- list(p = prior("dirichlet", list(alpha = c(2, 3))))
  columns <- c("p[1]", "p[2]", "p_copy", "prior_par_eta_p[1]", "prior_par_eta_p[2]")
  map <- .bt_build_parameter_map(columns, prior_list = priors)
  catalog <- .bt_parameter_map_catalog(map)
  id <- function(x) parameter_catalog_resolve(catalog, x)$quantity_id
  native <- parameter_relationships(map)
  expect_identical(native$source_id, c("p", "p"))
  expect_identical(native$consumer_id, c(id("p[1]"), id("p[2]")))
  input <- data.frame(source_id = id("p[1]"), consumer_id = id("p_copy"), relation = "model_input")
  registered <- register_parameter_relationships(map, input, "RoBMA")
  expect_equal(sum(parameter_relationships(registered)$provider == "RoBMA"), 1L)
  expect_error(register_parameter_relationships(map, transform(input, consumer_id = id("p[2]")), "RoBMA"),
               "fitted prior", class = "BayesTools_parameter_relationship_error")
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
  expect_setequal(leaves$block, c("lab", "study", "paper"))
  forged <- data.frame(source_id = root$quantity_id, consumer_id = leaves$consumer_id[1],
                       relation = "sd_source")
  expect_error(register_parameter_relationships(fit, forged, "RoBMA"), "own quantities",
               class = "BayesTools_parameter_relationship_error")
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

test_that("direct random blocks retain their external SD source and row inputs", {

  dat <- data.frame(id = factor(rep(1:2, each = 3)), lab = factor(rep(1:3, 2)))
  bind <- function(parameter, source){
    JAGS_formula(~ 1 + random(1 | id, name = "same", covariance = "diag") +
                     random(1 | lab, name = "ordinary", covariance = "diag"),
      parameter, dat, list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(same = random_block(sd_source = random_sd_source(source)),
                                  ordinary = random_block(sd = prior("gamma", list(2, 3)))))
  }
  scalar <- bind("mu", parameter_source("tau_total"))
  row <- bind("eta", parameter_source("row_sd", shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters$tau_total * exp(parameters[["offset[2]"]]), n_rows)
    }, inputs = c("tau_total", "offset[2]")))
  priors <- c(list(tau_total = prior("gamma", list(2, 3))), scalar$prior_list, row$prior_list)
  map <- .bt_build_parameter_map(names(priors), prior_list = priors,
    formula_design = list(mu = scalar$formula_design, eta = row$formula_design))
  edges <- parameter_relationships(map)
  root <- parameter_catalog_resolve(.bt_parameter_map_catalog(map), "tau_total")$quantity_id

  scales <- edges[edges$relation == "sd_source", , drop = FALSE]
  expect_identical(scales$source_type, c("quantity", "row_source"))
  expect_identical(scales$source_id, c(root, "row_sd"))
  expect_identical(scales$consumer_type, rep("random_block", 2))
  expect_identical(scales$formula_parameter, c("mu", "eta"))
  expect_identical(scales$block, rep("same", 2))
  expect_length(unique(scales$consumer_id), 2L)
  expect_identical(scales$allocation, rep("", 2))
  expect_false(any(edges$block == "ordinary"))

  inputs <- edges[edges$relation == "source_input", , drop = FALSE]
  expect_identical(inputs$source_type, c("quantity", "coordinate"))
  expect_identical(inputs$source_id, c(root, "offset[2]"))
  expect_identical(inputs$consumer_type, rep("row_source", 2))
  expect_identical(inputs$consumer_id, rep("row_sd", 2))
  expect_identical(inputs$formula_parameter, rep("eta", 2))
  expect_identical(inputs$block, rep("same", 2))
})

test_that("a provider connects an unmonitored random-effect SD to its public variance inputs", {

  dat <- data.frame(study = factor(rep(1:4, each = 6)), lab = factor(rep(1:3, 8)),
                    y = rep(c(-1, 0, 1), 8))
  random <- prior_random(random_variance_allocation("pool", terms = c(study = "study", lab = "lab"),
    sd_source = random_sd_source("random_sd"), weights = prior("dirichlet", list(c(2, 3)))))
  fit <- JAGS_fit(
    model_syntax = "model{random_sd <- tau_total * sqrt(1 - r)
                          for(i in 1:N){y[i] ~ dnorm(mu[i], 1)}}",
    data = list(y = dat$y, N = nrow(dat)),
    prior_list = list(tau_total = prior("gamma", list(2, 3)), r = prior("beta", list(2, 2))),
    formula_list = list(mu = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
                        random(1 | lab, name = "lab", covariance = "diag")),
    formula_data_list = list(mu = dat), formula_random_prior_list = list(mu = random),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
    chains = 1, sample = 100, burnin = 100, adapt = 100, seed = 35, silent = TRUE
  )
  catalog <- parameter_catalog(fit)
  id <- function(name) parameter_catalog_resolve(catalog, name)$quantity_id
  declarations <- data.frame(
    source_id = c(id("tau_total"), id("r")), consumer_id = "random_sd", consumer_type = "coordinate",
    relation = c("sd_source", "complement_variance_fraction")
  )
  registered <- register_parameter_relationships(fit, declarations, "RoBMA")
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(registered, saved)
  restored <- readRDS(saved)
  graph <- parameter_relationships(restored)
  provider <- graph[graph$provider == "RoBMA", ]
  expect_identical(provider$consumer_type, rep("coordinate", 2))
  expect_identical(provider$source_type, rep("quantity", 2))
  expect_identical(provider$source_id, declarations$source_id)
  expect_identical(parameter_relationships(register_parameter_relationships(catalog, declarations, "RoBMA")), graph)
  expect_identical(parameter_relationships(register_parameter_relationships(registered, declarations, "RoBMA")), graph)
  expect_identical(parameter_catalog(restored)$quantities, catalog$quantities)
  expect_false("random_sd" %in% parameter_coordinates(restored)$coordinate_name)
  expect_error(parameter_catalog_resolve(parameter_catalog(restored), "random_sd"),
               class = "BayesTools_parameter_not_found")
  # The native edge starts at the very same typed endpoint that this model supplies.
  root <- graph[graph$provider == "BayesTools" & graph$relation == "sd_source", ]
  expect_identical(unique(provider$consumer_type), root$source_type)
  expect_identical(unique(provider$consumer_id), root$source_id)
  split <- graph[graph$source_type == root$consumer_type & graph$source_id == root$consumer_id &
                 graph$relation == "allocates", ]
  blocks <- graph[graph$source_type == "allocation_component" & graph$source_id %in% split$consumer_id &
                  graph$relation == "scales_block", ]
  expect_setequal(blocks$block, c("study", "lab"))
  expect_identical(parameter_relationships(fit), graph[graph$provider == "BayesTools", , drop = FALSE])

  private_weight <- graph$source_id[graph$relation == "allocation_weight"][[1L]]
  for(target in c(private_weight, "undeclared_sd")){
    expect_error(register_parameter_relationships(fit, transform(declarations, consumer_id = target), "RoBMA"),
                 class = "BayesTools_parameter_relationship_error")
  }
  expect_error(register_parameter_relationships(fit, transform(declarations, consumer_type = "row_source"), "RoBMA"),
               class = "BayesTools_parameter_relationship_error")
  expect_error(register_parameter_relationships(registered, declarations, "Other"), "one provider",
               class = "BayesTools_parameter_relationship_error")
  inverse <- data.frame(source_id = "random_sd", source_type = "coordinate", consumer_id = id("tau_total"),
                        relation = "model_input")
  expect_error(register_parameter_relationships(registered, inverse, "RoBMA"), "fitted prior",
               class = "BayesTools_parameter_relationship_error")
})

test_that("provider row endpoints use their fitted shape and respect native callback inputs", {

  dat <- data.frame(group = factor(rep(letters[1:3], each = 2)))
  bind <- function(name, source){
    JAGS_formula(~ 1 + random(1 | group, name = "group", covariance = "diag"),
      name, dat, list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(group = random_block(sd_source = random_sd_source(source))))
  }
  raw <- bind("mu", parameter_source("row_sd", shape = "row"))
  callback <- bind("eta", parameter_source("callback_sd", shape = "row",
    values = function(parameters, data, n_rows) rep(parameters[["offset[2]"]], n_rows), inputs = "offset[2]"))
  data_only <- bind("theta", parameter_source("known_sd", shape = "row",
    values = function(parameters, data, n_rows) rep(1, n_rows), inputs = character()))
  priors <- c(list(tau_total = prior("gamma", list(2, 2))), raw$prior_list, callback$prior_list, data_only$prior_list)
  map <- .bt_build_parameter_map(names(priors), prior_list = priors,
    formula_design = list(mu = raw$formula_design, eta = callback$formula_design, theta = data_only$formula_design))
  tau <- parameter_catalog_resolve(.bt_parameter_map_catalog(map), "tau_total")$quantity_id
  declarations <- data.frame(source_id = tau, consumer_type = c("row_source", "coordinate"),
                             consumer_id = c("row_sd", "offset[2]"), relation = "model_input")
  registered <- register_parameter_relationships(map, declarations, "RoBMA")
  added <- subset(parameter_relationships(registered), provider == "RoBMA")
  expect_identical(added$consumer_type, declarations$consumer_type)
  expect_identical(added$consumer_id, declarations$consumer_id)
  callback_input <- subset(parameter_relationships(registered), relation == "source_input")
  expect_identical(callback_input$source_id, added$consumer_id[2L])
  expect_identical(callback_input$source_type, added$consumer_type[2L])
  expect_error(register_parameter_relationships(map,
      transform(declarations[1L, ], consumer_id = "callback_sd"), "RoBMA"), "callback",
    class = "BayesTools_parameter_relationship_error")
  data_edge <- subset(parameter_relationships(map), relation == "callback_data")
  expect_identical(data_edge$source_type, "row_data")
  expect_identical(data_edge$source_id, "known_sd")
  expect_identical(data_edge$consumer_type, "row_source")
  expect_identical(data_edge$consumer_id, "known_sd")
  expect_error(register_parameter_relationships(map,
      transform(declarations[1L, ], consumer_id = "known_sd"), "RoBMA"), "callback",
    class = "BayesTools_parameter_relationship_error")
  expect_error(register_parameter_relationships(map,
      transform(declarations[1L, ], consumer_type = "row_data", consumer_id = "known_sd"), "RoBMA"),
    class = "BayesTools_parameter_relationship_error")
  # Including a row source as an input is distinct from redefining its callback.
  cycle <- data.frame(source_type = "row_source", source_id = "callback_sd",
                      consumer_type = "coordinate", consumer_id = "offset[2]", relation = "model_input")
  expect_error(register_parameter_relationships(registered, cycle, "RoBMA"), "cycle",
               class = "BayesTools_parameter_relationship_error")
})
