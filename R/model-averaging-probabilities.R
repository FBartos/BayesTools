# Original model declarations and their visible numeric probability boundary.
.model_probability_range_stop <- function(indices){

  .prior_numerical_signal("model-probability calculation", "model", "log", indices,
    "Finite model-score arithmetic is outside the supported logarithmic range",
    error = TRUE)
}

.model_probability_log_sum <- function(logs){

  if(!any(is.finite(logs))) return(-Inf)
  center <- max(logs)
  result <- center + log(sum(exp(logs - center)))
  if(!is.finite(result)) .model_probability_range_stop(which(is.finite(logs)))
  result
}

.model_probability_normalize_logs <- function(logs){

  if(!any(is.finite(logs))) stop("At least one model probability must be positive.", call. = FALSE)
  centered <- logs - max(logs)
  if(any(!is.finite(centered[is.finite(logs)]))) .model_probability_range_stop(which(is.finite(logs)))
  normalized <- centered - log(sum(exp(centered)))
  if(any(!is.finite(normalized[is.finite(logs)]))) .model_probability_range_stop(which(is.finite(logs)))
  unname(normalized)
}

.model_probability_eta <- function(raw_logs, logs){

  64 * .Machine$double.eps * max(1, abs(raw_logs[is.finite(raw_logs)]), abs(logs[is.finite(logs)]))
}

.model_probability_safe <- function(probabilities, logs, eta, intermediates = probabilities){

  positive <- is.finite(logs)
  all(is.finite(intermediates)) &&
    all(intermediates[intermediates != 0] >= .Machine$double.xmin) &&
    all(probabilities[positive] >= .Machine$double.xmin) &&
    all(probabilities[!positive] == 0) &&
    all(abs(log(probabilities[positive]) - logs[positive]) <= eta)
}

.model_probability_pair <- function(probabilities, logs, stage, route = "stabilized",
                                    eta = .model_probability_eta(logs, logs), model_indices = seq_along(logs)){

  declaration <- list(schema_version = 1L, model_indices = as.integer(model_indices),
    stage = stage, route = route, eta = eta)
  .model_probability_validate(probabilities, logs, declaration)
  list(probabilities = probabilities, logs = logs, declaration = declaration)
}

.model_probability_validate <- function(probabilities, logs, declaration, normalized = FALSE){

  valid <- is.list(declaration) && identical(names(declaration),
    c("schema_version", "model_indices", "stage", "route", "eta")) &&
    identical(declaration$schema_version, 1L) &&
    is.integer(declaration$model_indices) && length(declaration$model_indices) == length(logs) &&
    !anyNA(declaration$model_indices) && all(declaration$model_indices > 0L) &&
    !anyDuplicated(declaration$model_indices) &&
    is.character(declaration$stage) && length(declaration$stage) == 1L &&
    declaration$stage %in% c("prior", "posterior", "conditional_prior", "conditional_posterior", "event", "component") &&
    is.character(declaration$route) && length(declaration$route) == 1L &&
    declaration$route %in% c("ordinary", "stabilized", "raw") &&
    is.numeric(declaration$eta) && length(declaration$eta) == 1L &&
    is.finite(declaration$eta) && declaration$eta >= 0 &&
    is.numeric(probabilities) && is.numeric(logs) && length(probabilities) == length(logs) &&
    !anyNA(probabilities) && all(is.finite(probabilities)) && all(probabilities >= 0) &&
    !anyNA(logs) && all(is.finite(logs) | logs == -Inf)
  if(valid){
    positive <- is.finite(logs)
    visible <- probabilities > 0
    subnormal <- visible & probabilities < .Machine$double.xmin
    compare_logs <- if(declaration$route == "stabilized") visible & !subnormal else visible
    valid <- all(probabilities[!positive] == 0) &&
      all(exp(logs[positive & !visible]) == 0) &&
      (declaration$route == "stabilized" || !any(positive & !visible)) &&
      (declaration$route != "ordinary" || !any(subnormal)) &&
      (declaration$route != "stabilized" || all(probabilities[subnormal] == exp(logs[subnormal]))) &&
      all(abs(log(probabilities[compare_logs]) - logs[compare_logs]) <= declaration$eta)
  }
  if(!isTRUE(valid)) stop("Model probability ownership is missing or malformed. Recompute or refit with the current BayesTools version.", call. = FALSE)
  if(normalized){
    .inclusion_BF_check_probs(probabilities, "model_probabilities")
    if(abs(.model_probability_log_sum(logs)) > declaration$eta) stop("Model log probabilities are not normalized.", call. = FALSE)
  }
  invisible(TRUE)
}

.model_probability_prior <- function(weights, logs = log(weights), stage = "prior",
                                      ordinary = NULL, prior_eta = 0){

  if(any(!is.finite(weights))) stop("'prior_weights' must be finite.", call. = FALSE)
  if(!any(is.finite(logs))) stop("At least one prior model weight must be positive.", call. = FALSE)
  canonical <- .model_probability_normalize_logs(logs)
  eta <- max(prior_eta, .model_probability_eta(logs, canonical))
  if(is.null(ordinary)){
    scaled <- weights / max(weights)
    ordinary <- scaled / sum(scaled)
  }else scaled <- ordinary
  safe <- .model_probability_safe(ordinary, canonical, eta, scaled) &&
    all(scaled[is.finite(logs)] >= .Machine$double.xmin)
  .model_probability_pair(if(safe) ordinary else exp(canonical), canonical, stage,
    if(safe) "ordinary" else "stabilized", eta)
}

.model_probability_plot_priors <- function(priors, allow_weightfunction_null = FALSE, validated_pair = NULL){

  if(is.prior(priors)) return(priors)
  pair <- validated_pair
  if(is.null(pair)){
    pair <- .prior_model_probability_pair(priors)
    .model_probability_plot_check(priors, allow_weightfunction_null = allow_weightfunction_null)
  }else .model_probability_validate(pair$probabilities, pair$logs, pair$declaration, normalized = TRUE)
  weights <- vapply(priors, .prior_model_weight, numeric(1))
  if(.model_probability_safe(weights / sum(weights), pair$logs, pair$declaration$eta)) return(priors)
  for(i in seq_along(priors)){
    owner <- pair$declaration
    owner$model_indices <- as.integer(i)
    priors[[i]] <- .set_prior_model_probability(priors[[i]], pair$probabilities[[i]], pair$logs[[i]], owner)
  }
  priors
}

.model_probability_posterior <- function(evidence, prior){

  active <- is.finite(prior$logs) & is.finite(evidence)
  if(!any(active)) stop("No finite marginal likelihoods are available for models with positive prior probability.", call. = FALSE)
  centered <- evidence[active] - max(evidence[active])
  if(any(!is.finite(centered))) .model_probability_range_stop(which(active)[!is.finite(centered)])
  scores <- rep(-Inf, length(evidence))
  scores[active] <- prior$logs[active] + centered
  if(any(!is.finite(scores[active]))) .model_probability_range_stop(which(active)[!is.finite(scores[active])])
  canonical <- .model_probability_normalize_logs(scores)
  eta <- max(prior$declaration$eta, .model_probability_eta(prior$logs, canonical))
  old_scores <- rep(-Inf, length(evidence))
  old_scores[active] <- log(prior$probabilities[active]) + evidence[active]
  old_weights <- exp(old_scores - max(old_scores))
  ordinary <- unname(old_weights / sum(old_weights))
  safe <- .model_probability_safe(ordinary, canonical, eta, old_weights) &&
    all(old_weights[active] >= .Machine$double.xmin)
  .model_probability_pair(if(safe) ordinary else exp(canonical), canonical, "posterior",
    if(safe) "ordinary" else "stabilized", eta)
}

.model_probability_inference_set <- function(inference, prior, posterior){

  .model_probability_validate(prior$probabilities, prior$logs, prior$declaration)
  .model_probability_validate(posterior$probabilities, posterior$logs, posterior$declaration)
  inference$prior_probs <- prior$probabilities
  inference$post_probs <- posterior$probabilities
  attr(inference, "log_prior_probs") <- prior$logs
  attr(inference, "log_post_probs") <- posterior$logs
  attr(inference, "model_probability_declaration") <- list(prior = prior$declaration,
    posterior = posterior$declaration)
  inference
}

.model_probability_inference_get <- function(inference, stage){

  prior <- stage == "prior"
  probabilities <- inference[[if(prior) "prior_probs" else "post_probs"]]
  logs <- attr(inference, if(prior) "log_prior_probs" else "log_post_probs", exact = TRUE)
  declaration <- attr(inference, "model_probability_declaration", exact = TRUE)[[if(prior) "prior" else "posterior"]]
  .model_probability_validate(probabilities, logs, declaration, normalized = TRUE)
  stages <- if(prior) c("prior", "conditional_prior") else c("posterior", "conditional_posterior")
  if(!declaration$stage %in% stages ||
     !identical(declaration$model_indices, seq_along(logs))) stop("Inference model probability stage or identities are malformed.", call. = FALSE)
  list(probabilities = probabilities, logs = logs, declaration = declaration)
}

.model_probability_scalar_inference_validate <- function(inference){

  declarations <- attr(inference, "model_probability_declaration", exact = TRUE)
  prior <- attr(inference, "log_prior_prob", exact = TRUE)
  posterior <- attr(inference, "log_post_prob", exact = TRUE)
  if(is.null(declarations) && is.null(prior) && is.null(posterior)) return(invisible(TRUE))
  if(!is.list(declarations) || !identical(names(declarations), c("prior", "posterior"))) stop("Model inference probability ownership is malformed.", call. = FALSE)
  .model_probability_validate(inference$prior_prob, prior, declarations$prior)
  .model_probability_validate(inference$post_prob, posterior, declarations$posterior)
  if(!declarations$prior$stage %in% c("prior", "conditional_prior") ||
     !declarations$posterior$stage %in% c("posterior", "conditional_posterior")) stop("Model inference probability stages are malformed.", call. = FALSE)
  if(!identical(declarations$prior$model_indices, as.integer(inference$m_number)) ||
     !identical(declarations$posterior$model_indices, as.integer(inference$m_number))){
    stop("Model inference probability identities are malformed.", call. = FALSE)
  }
  invisible(TRUE)
}

.model_probability_condition <- function(pair, event_logs, stage, event_probabilities = exp(event_logs)){

  logs <- pair$logs + event_logs
  finite <- is.finite(pair$logs) & is.finite(event_logs)
  if(any(!is.finite(logs[finite]))) .model_probability_range_stop(which(finite)[!is.finite(logs[finite])])
  canonical <- .model_probability_normalize_logs(logs)
  eta <- max(pair$declaration$eta, .model_probability_eta(logs, canonical))
  product <- pair$probabilities * event_probabilities
  ordinary <- product / sum(product)
  if(stage != "prior") ordinary <- unname(ordinary)
  safe <- .model_probability_safe(ordinary, canonical, eta, product) &&
    all(product[finite] >= .Machine$double.xmin)
  .model_probability_pair(if(safe) ordinary else exp(canonical), canonical, stage,
    if(safe) "ordinary" else "stabilized", eta, pair$declaration$model_indices)
}

.model_probability_diagnostics <- function(prior = NULL, posterior = NULL, stage){

  pair <- if(!is.null(prior)) prior else posterior
  list(model_indices = pair$declaration$model_indices,
    log_prior_probabilities = if(!is.null(prior)) prior$logs else NULL,
    log_posterior_probabilities = if(!is.null(posterior)) posterior$logs else NULL,
    stage = stage)
}

.model_probability_measure_stop <- function(pair, measure = "prior_density", prior = TRUE){

  .model_probability_validate(pair$probabilities, pair$logs, pair$declaration)
  if(all(pair$probabilities[is.finite(pair$logs)] >= .Machine$double.xmin)) return(invisible(TRUE))
  diagnostics <- .model_probability_diagnostics(if(prior) pair else NULL, if(!prior) pair else NULL,
    pair$declaration$stage)
  message <- paste0("The declared ", measure,
    " is unavailable because active model probabilities are not representable at full precision.")
  if(measure == "prior_density") .bt_formula_density_stop(message,
    reason = "numerical_model_probability_unavailable", detail = message, diagnostics = diagnostics)
  stop(errorCondition(message, call = NULL,
    class = c(paste0("BayesTools_formula_", measure, "_unavailable"), "BayesTools_formula_measure_unavailable"),
    reason = "numerical_model_probability_unavailable", detail = message, diagnostics = diagnostics))
}

.model_probability_prior_atoms <- function(priors, pair, n_columns, column_names,
                                           source, null_location, exclusion_probabilities, point_locations = NULL){

  .model_probability_validate(pair$probabilities, pair$logs, pair$declaration)
  locations <- matrix(numeric(), 0L, n_columns)
  log_masses <- numeric()
  continuous <- FALSE
  for(i in which(is.finite(pair$logs))){
    location <- if(!is.null(point_locations) && !is.null(point_locations[[i]])) point_locations[[i]] else
      .posterior_atoms_point_location(priors[[i]], n_columns)
    fraction <- 1
    if(is.null(location) && !is.null(null_location) && .is_prior_weightfunction_null(priors[[i]])) location <- rep(null_location, n_columns)
    if(is.null(location)){
      fraction <- if(is.null(exclusion_probabilities)) 0 else exclusion_probabilities[[i]]
      continuous <- continuous || fraction < 1
      if(fraction == 0) next
      location <- rep(0, n_columns)
    }
    locations <- rbind(locations, location)
    log_masses <- c(log_masses, pair$logs[[i]] + log(fraction))
  }
  unique_locations <- unique(locations)
  logs <- vapply(seq_len(nrow(unique_locations)), function(row){
    selected <- rep(TRUE, nrow(locations))
    for(column in seq_len(n_columns)) selected <- selected & locations[, column] == unique_locations[row, column]
    .model_probability_log_sum(log_masses[selected])
  }, numeric(1))
  masses <- exp(logs)
  lost <- any(masses == 0 | (masses > 0 & masses < .Machine$double.xmin)) ||
    (continuous && length(masses) > 0L && sum(masses) == 1)
  if(lost){
    diagnostics <- .model_probability_diagnostics(posterior = pair, stage = pair$declaration$stage)
    stop(errorCondition("Posterior atoms are unavailable because positive model mass or a continuous remainder is not representable at full precision.",
      call = NULL, class = c("BayesTools_formula_atoms_unavailable", "BayesTools_formula_measure_unavailable"),
      reason = "numerical_model_probability_unavailable", diagnostics = diagnostics))
  }
  if(all(pair$probabilities[is.finite(pair$logs)] >= .Machine$double.xmin)){
    ordinary <- .posterior_atoms_from_priors(priors, pair$probabilities, n_columns, column_names,
      source, null_location, exclusion_probabilities, point_locations = point_locations)
    if(all(ordinary$mass >= .Machine$double.xmin) && !(continuous && sum(ordinary$mass) == 1)){
      ordinary$component_log_probabilities <- pair$logs
      ordinary$model_probability_declaration <- pair$declaration
      return(.posterior_atoms_from_attribute(ordinary))
    }
  }
  .posterior_atoms_new(unique_locations, masses, column_names, source,
    component_probabilities = pair$probabilities, component_log_probabilities = pair$logs,
    model_probability_declaration = pair$declaration)
}

.model_probability_mixed_measures <- function(samples, prior, posterior, priors, parameter){

  samples <- .bt_meta_set(samples, "model_probabilities", list(prior = prior, posterior = posterior))
  constant <- .model_probability_point_declarations(lapply(which(is.finite(prior$logs)), function(i){
    list(.bt_prior_without_multiply_by(priors[[i]]))
  }))
  constant <- constant || (inherits(samples, "mixed_posteriors.weightfunction") &&
    all(vapply(priors[is.finite(prior$logs)], .is_prior_weightfunction_null, logical(1))))
  if(inherits(samples, "mixed_posteriors.weightfunction")){
    locations <- .model_probability_weightfunction_points(priors, weightfunctions_mapping(priors), ncol(samples))
    locations <- locations[is.finite(prior$logs)]
    constant <- constant || (!any(vapply(locations, is.null, logical(1))) &&
      all(vapply(locations, identical, logical(1), locations[[1L]])))
  }
  unavailable <- if(constant) NULL else tryCatch(.model_probability_measure_stop(prior),
    BayesTools_formula_measure_unavailable = function(condition) condition)
  if(inherits(unavailable, "BayesTools_formula_measure_unavailable")){
    samples <- .bt_meta_assign(samples, list(prior_density = NULL, prior_densities = NULL, prior_context = NULL))
    columns <- if(is.matrix(samples)) colnames(samples) else parameter
    for(column in columns) samples <- .bt_formula_measure_mark(samples, column, "prior_density",
      unavailable$detail, cause = unavailable$reason, diagnostics = unavailable$diagnostics)
  }
  samples
}

.model_probability_weightfunction_points <- function(priors, mapping, n_columns){

  lapply(seq_along(priors), function(i){
    prior <- priors[[i]]
    if(.is_prior_weightfunction_null(prior)) return(rep(1, n_columns))
    .weightfunction_validate_steps(prior$steps)
    .weightfunction_validate_weights(prior$weights, .weightfunction_n_bins(prior), prior$reference)
    if(!identical(prior$weights$type, "fixed")) return(NULL)
    indices <- mapping[[i]]
    check_int(indices, "omega_mapping", lower = 1, upper = length(prior$weights$omega),
      check_length = n_columns, allow_NA = FALSE)
    unname(prior$weights$omega[indices])
  })
}

.model_probability_component_set <- function(samples, component, priors, probabilities, posterior = NULL){

  if(is.null(posterior)) posterior <- .model_probability_pair(probabilities, log(probabilities), "posterior", "raw")
  samples <- .bt_meta_set(samples, "model_probabilities", list(prior = .prior_model_probability_pair(priors), posterior = posterior))
  .bt_draws_set_component(samples, component, "model")
}

.model_probability_context_validate <- function(context){

  if(!inherits(context, "prior_density_model_mixture_context") &&
     !inherits(context, "prior_density_conditional_context")) return(invisible(TRUE))
  if(!.bt_formula_prior_density_context_valid(context)) stop("The model prior context is incomplete or malformed. Recompute or refit with the current BayesTools version.", call. = FALSE)
  if(!identical(context$schema_version, 2L)) stop("Model prior contexts require current probability ownership. Recompute or refit with the current BayesTools version.", call. = FALSE)
  .model_probability_validate(context$model_weights, context$model_log_weights,
    context$model_probability_declaration, normalized = TRUE)
  if(!context$model_probability_declaration$stage %in% c("prior", "conditional_prior", "event", "component")){
    stop("The model prior context probability stage is malformed.", call. = FALSE)
  }
  invisible(TRUE)
}

.model_probability_context_pair <- function(context){

  .model_probability_context_validate(context)
  list(probabilities = context$model_weights, logs = context$model_log_weights,
    declaration = context$model_probability_declaration)
}

.model_probability_context_measure_check <- function(context, weights = NULL){

  if(inherits(context, "prior_density_model_mixture_context") ||
     inherits(context, "prior_density_conditional_context")){
    .model_probability_context_validate(context)
    .bt_formula_context_check(context)
    if(!is.null(weights) && all(is.finite(weights)) && all(weights == 0)) return(invisible(TRUE))
    indices <- which(is.finite(context$model_log_weights))
    priors <- if(inherits(context, "prior_density_model_mixture_context")) lapply(indices, function(i){
      .prior_density_model_prior_list(context$prior_list, i)
    }) else context$prior_lists[indices]
    if(.model_probability_point_declarations(priors)) return(invisible(TRUE))
    .model_probability_measure_stop(.model_probability_context_pair(context))
  }
  invisible(TRUE)
}

.model_probability_plot_check <- function(priors, allow_weightfunction_null = FALSE){

  if(is.prior(priors)) priors <- list(priors)
  if(length(priors)){
    pair <- .prior_model_probability_pair(priors)
    if(allow_weightfunction_null && all(vapply(priors[is.finite(pair$logs)], .is_prior_weightfunction_null, logical(1)))) return(invisible(TRUE))
    if(.model_probability_point_declarations(lapply(priors[is.finite(pair$logs)], list))) return(invisible(TRUE))
    .model_probability_measure_stop(pair)
  }
  invisible(TRUE)
}

.model_probability_petpeese_plot_check <- function(priors, mu_priors){

  check_list(priors, "prior_list")
  check_list(mu_priors, "prior_list_mu", check_length = length(priors))
  for(prior in c(priors, mu_priors)) .check_prior(prior)
  pair <- .prior_model_probability_pair(priors)
  mu_pair <- .prior_model_probability_pair(mu_priors)
  mu_owned <- any(vapply(mu_priors, function(prior){
    !is.null(attr(prior, "model_probability_declaration", exact = TRUE))
  }, logical(1)))
  if(mu_owned && (!identical(mu_pair$probabilities, pair$probabilities) ||
     !identical(mu_pair$logs, pair$logs) ||
     !identical(mu_pair$declaration$model_indices, pair$declaration$model_indices))){
    stop("The PET-PEESE model prior distributions are not aligned across parameters.", call. = FALSE)
  }
  point <- function(prior){
    if(is.prior.none(prior)) return(0)
    if(!is.null(attr(prior, "multiply_by", exact = TRUE))) return(NULL)
    .posterior_atoms_point_location(prior, 1L)
  }
  declarations <- lapply(which(is.finite(pair$logs)), function(i){
    mu <- point(mu_priors[[i]])
    type <- if(is.prior.PET(priors[[i]])) "PET" else if(is.prior.PEESE(priors[[i]])) "PEESE" else "none"
    bias <- if(type == "none") 0 else point(priors[[i]])
    if(is.null(mu) || is.null(bias)) return(NULL)
    list(mu = mu, bias = bias, type = type)
  })
  complete <- !any(vapply(declarations, is.null, logical(1)))
  if(complete){
    mus <- vapply(declarations, `[[`, numeric(1), "mu")
    biases <- vapply(declarations, `[[`, numeric(1), "bias")
    types <- vapply(declarations, `[[`, character(1), "type")
    complete <- all(mus == mus[[1L]]) &&
      (all(biases == 0) || (all(types == types[[1L]]) && all(biases == biases[[1L]])))
  }
  if(!complete) .model_probability_measure_stop(pair)
  pair
}

.model_probability_point_declarations <- function(prior_lists){

  if(!length(prior_lists)) return(FALSE)
  points <- lapply(prior_lists, function(priors){
    lapply(priors, function(prior){
      if(!is.null(attr(prior, "multiply_by", exact = TRUE))) return(NULL)
      .posterior_atoms_point_location(prior, .prior_linear_prior_dimension(prior))
    })
  })
  if(any(vapply(points, function(model) any(vapply(model, is.null, logical(1))), logical(1)))) return(FALSE)
  all(vapply(points, identical, logical(1), points[[1L]]))
}

.model_probability_split_prior <- function(parent, child, fraction){

  weight <- .prior_model_weight(parent)
  if(is.null(weight)) weight <- 1
  owner <- attr(parent, "model_probability_declaration", exact = TRUE)
  if(is.null(owner)) return(.set_prior_model_weight(child, weight * fraction))
  logs <- .prior_model_log_weight(parent) + log(fraction)
  child_owner <- attr(child, "model_probability_declaration", exact = TRUE)
  if(!is.null(child_owner) && abs(.prior_model_log_weight(child) - logs) > max(owner$eta, child_owner$eta)){
    stop("An independent child model probability owner contradicts the requested parent split.", call. = FALSE)
  }
  owner$stage <- "component"
  owner$eta <- max(owner$eta, .model_probability_eta(logs, logs))
  value <- weight * fraction
  if(!.model_probability_safe(value, logs, owner$eta)){
    value <- exp(logs)
    owner$route <- "stabilized"
  }
  .set_prior_model_probability(child, value, logs, owner)
}

.model_probability_reindex_priors <- function(priors){

  owners <- lapply(priors, attr, which = "model_probability_declaration", exact = TRUE)
  if(!length(priors) || any(vapply(owners, is.null, logical(1)))) return(priors)
  indices <- vapply(owners, function(owner) owner$model_indices, integer(1))
  if(!anyDuplicated(indices)) return(priors)
  for(i in seq_along(priors)){
    owner <- owners[[i]]
    owner$model_indices <- as.integer(i)
    owner$stage <- "component"
    priors[[i]] <- .set_prior_model_probability(priors[[i]], .prior_model_weight(priors[[i]]),
      .prior_model_log_weight(priors[[i]]), owner)
  }
  priors
}

.model_probability_atom_table <- function(points, continuous, pair){

  if(is.null(points) || !nrow(points)) return(data.frame(x = numeric(), mass = numeric()))
  locations <- unique(points$x)
  logs <- vapply(locations, function(location) .model_probability_log_sum(points$log_mass[points$x == location]), numeric(1))
  mass <- exp(logs)
  if(any(mass < .Machine$double.xmin) || (continuous && sum(mass) == 1)){
    return(errorCondition("Posterior atoms are unavailable because positive model mass or a continuous remainder is not representable at full precision.",
      call = NULL, class = c("BayesTools_formula_atoms_unavailable", "BayesTools_formula_measure_unavailable"),
      reason = "numerical_model_probability_unavailable",
      diagnostics = .model_probability_diagnostics(posterior = pair, stage = pair$declaration$stage)))
  }
  ordinary <- .posterior_atoms_point_mass_table(points[c("x", "mass")])
  if(!is.null(ordinary) && all(points$mass >= .Machine$double.xmin) &&
     all(pair$probabilities[is.finite(pair$logs)] >= .Machine$double.xmin)) return(ordinary)
  .posterior_atoms_point_mass_table(data.frame(x = locations, mass = mass))
}
