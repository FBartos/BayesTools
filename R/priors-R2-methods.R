# Distribution methods of bound R2 term priors: random generation through the
# registered coefficient node, and marginal densities, distribution functions,
# quantiles, and moments of single coefficients or factor levels from the
# deterministic quadrature of R/priors-R2-density.R.

# Specification of the implied law of a linear functional L'beta of the term
# coefficients (in residual-SD units).
.bt_R2_term_functional <- function(prior, L){

  .bt_R2_require_bound(prior)
  L <- as.numeric(L)
  if(length(L) != prior$d){
    stop("Internal R2 functional dimension does not match the term.", call. = FALSE)
  }
  shapes <- prior$shapes
  if(identical(prior$law, "exact")){
    scale <- sqrt(sum((t(prior$map) %*% L)^2))
    unit <- prior$units[1L]
    d <- prior$d
  }else{
    coordinate <- which(L != 0)
    if(length(coordinate) != 1L){
      .bt_R2_unavailable(paste0(
        "The marginal prior distribution of a combination of several coded coefficients is unavailable for the '",
        prior$law, "' kernel; it is available for single coded coefficients, and for factor levels under the 'exact' law."),
        class = "BayesTools_R2_density_unavailable")
    }
    scale <- abs(L[coordinate] * prior$map[coordinate, coordinate])
    unit <- prior$units[coordinate]
    d <- 1L
  }
  .bt_R2_functional_spec(
    law     = prior$law,
    a       = shapes$a,
    b       = shapes$b,
    alpha   = shapes$alpha[unit],
    A       = shapes$A,
    n_units = prior$n_units,
    d       = d,
    scale   = scale
  )
}

# The rows whose marginals the marginal methods describe: the coefficient of a
# continuous term, the level effects of a mean-difference or orthonormal
# factor (rows of its declared design), and the coded coefficients of a
# treatment factor. Rows that vanish identically are atoms at zero, not
# continuous levels.
.bt_R2_marginal_rows <- function(prior){

  if(!is.prior.factor(prior)){
    return(matrix(1, 1, 1))
  }
  if(is.prior.orthonormal(prior) || is.prior.meandif(prior)){
    design <- attr(prior, "factor_design", exact = TRUE)
    if(is.matrix(design)){
      return(design)
    }
    contrast <- if(is.prior.orthonormal(prior)) contr.orthonormal else contr.meandif
    return(contrast(prior$d + 1L))
  }
  diag(prior$d)
}

# One common marginal of the term's rows, or a refusal when the rows differ.
.bt_R2_marginal_spec <- function(prior){

  rows <- .bt_R2_marginal_rows(prior)
  rows <- rows[rowSums(abs(rows)) > 0, , drop = FALSE]
  specs <- lapply(seq_len(nrow(rows)), function(i) .bt_R2_term_functional(prior, rows[i, ]))
  scales <- vapply(specs, `[[`, numeric(1), "scale")
  # the allocation shape of the units (absent with a single unit)
  shapes <- vapply(specs, function(spec) if(spec$n_units > 1L) spec$alpha else 1, numeric(1))
  if(length(specs) > 1L && (max(abs(scales - scales[1L])) > 1e-10 * scales[1L] ||
                            max(abs(shapes - shapes[1L])) > 1e-12 * shapes[1L])){
    .bt_R2_unavailable(paste0(
      "The levels of the R-squared budget term '", gsub("__xXx__", ":", prior$term, fixed = TRUE),
      "' have different marginal prior distributions (unequal cell counts or allocations), so there is no common marginal. ",
      "Use rng() for the joint distribution of the levels."),
      class = "BayesTools_R2_density_unavailable")
  }
  specs[[1L]]
}

# Draws of the term coefficients (rows) through the registered node, from
# draws of the budget total, the allocation, and the direction coordinates.
.bt_R2_term_rng <- function(prior, n, parameter_name = "R2_term"){

  .bt_R2_require_bound(prior)
  shapes <- prior$shapes
  total <- stats::rbeta(n, shapes$a, shapes$b)
  draws <- matrix(total, ncol = 1L, dimnames = list(NULL, prior$total_name))
  budget_priors <- .bt_R2_budget_prior_list(prior)
  if(prior$n_units > 1L){
    gammas <- vapply(shapes$alpha, function(alpha) stats::rgamma(n, shape = alpha), numeric(n))
    gammas <- matrix(gammas, nrow = n, dimnames = list(NULL, .bt_R2_allocation_eta_coordinates(prior)))
    draws <- cbind(draws, gammas)
  }
  directions <- matrix(.bt_R2_direction_rng(prior$law, n * prior$d), nrow = n)
  colnames(directions) <- .bt_R2_direction_coordinates(prior, parameter_name)
  draws <- cbind(draws, directions)
  node <- .bt_dnode_R2_coefficient(parameter_name, prior)
  values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(draws, budget_priors))
  unname(values)
}

.bt_R2_rng <- function(prior, n, transform_factor_samples = TRUE){

  values <- .bt_R2_term_rng(prior, n)
  if(!is.prior.factor(prior)){
    return(as.vector(values))
  }
  if(transform_factor_samples && (is.prior.orthonormal(prior) || is.prior.meandif(prior))){
    rows <- .bt_R2_marginal_rows(prior)
    return(values %*% t(rows))
  }
  values
}

.bt_R2_joint_refusal <- function(prior, method){

  if(is.prior.factor(prior)){
    stop("Joint ", method, " are unavailable for factor priors. Use the marginal ",
         "methods (e.g., 'mquant()') for the marginal distributions.", call. = FALSE)
  }
  invisible(NULL)
}

# The classed refusal of the linear-combination densities (formula prior
# densities, transformed and scaled prior densities, contexts of plots and
# hypotheses) for R2 budget parameters: a formula-measure condition, so that
# the consumers report the measure as unavailable.
.bt_R2_linear_density_unavailable <- function(parameter){

  .bt_formula_density_stop(
    paste0("The prior density of a linear combination involving '", parameter,
           "' is unavailable: the terms of an R-squared budget are dependent through the budget. ",
           "Use the marginal prior methods of the term (e.g., 'mpdf()') or prior draws."),
    reason = "unsupported_contribution_measure", prior_family = "R2"
  )
}

# The groups of a linear combination are convolved as independent terms:
# parameters of one R-squared budget (its terms, total, and allocation) are
# dependent, and so are the coded coefficients of an R2 factor (combinations
# of its levels have their marginal methods instead). A single R2 scalar term,
# combined with terms outside its budget, keeps its exact marginal law.
.bt_R2_check_linear_groups <- function(groups){

  budget <- vapply(groups, function(group){
    prior <- group$prior
    if(is.prior.R2(prior)){
      return(prior$total_name)
    }
    if(.bt_R2_is_budget_prior(prior)){
      return(sub("_allocation$", "", group$parameter))
    }
    NA_character_
  }, character(1))
  shared <- unique(budget[!is.na(budget) & duplicated(budget)])
  if(length(shared) > 0L){
    .bt_R2_linear_density_unavailable(names(groups)[budget %in% shared][1L])
  }
  for(group in groups[!is.na(budget)]){
    # several coded coefficients of a term (or of the allocation) are dependent
    if(is.prior.vector(group$prior) || length(group$weights) > 1L){
      .bt_R2_linear_density_unavailable(group$parameter)
    }
  }
  invisible(TRUE)
}

# Dispatcher of the distribution methods of an R2 term prior.
.bt_R2_prior_method <- function(prior, method, value = NULL, ...){

  .bt_R2_require_bound(prior)
  if(method %in% c("cdf", "ccdf", "pdf", "lpdf", "quant")){
    .bt_R2_joint_refusal(prior, switch(method, quant = "quantile functions", paste0(method, "s")))
  }
  switch(
    method,
    rng = .bt_R2_rng(prior, value, ...),
    cdf = , mcdf = .bt_R2_functional_cdf(.bt_R2_marginal_spec(prior), as.vector(value)),
    ccdf = , mccdf = .bt_R2_functional_cdf(.bt_R2_marginal_spec(prior), as.vector(value), lower.tail = FALSE),
    pdf = , mpdf = .bt_R2_functional_pdf(.bt_R2_marginal_spec(prior), as.vector(value)),
    lpdf = , mlpdf = .bt_R2_functional_log_pdf(.bt_R2_marginal_spec(prior), as.vector(value)),
    quant = , mquant = .bt_R2_functional_quantile(.bt_R2_marginal_spec(prior), as.vector(value)),
    mean = if(prior$shapes$b > 1 / 2) 0 else NaN,
    var = {
      variance <- .bt_R2_functional_variance(.bt_R2_marginal_spec(prior))
      if(is.finite(variance)) variance else NaN
    },
    stop("Internal R2 prior method '", method, "' is unknown.", call. = FALSE)
  )
}

# Joint prior draws of the parameters of a prior list that belong to an
# R-squared budget: the budget total and allocation are drawn once per budget
# (kept in 'registry'), and every term's coefficients are its registered node
# evaluated on those draws, so that the columns keep the joint prior law of
# the fitted model.
.bt_R2_budget_registry_draws <- function(registry, prior_list, total_name, n){

  if(exists(total_name, envir = registry, inherits = FALSE)){
    return(get(total_name, envir = registry, inherits = FALSE))
  }
  terms <- prior_list[vapply(prior_list, function(prior){
    is.prior.R2(prior) && isTRUE(prior$bound) && identical(prior$total_name, total_name)
  }, logical(1))]
  budget <- prior_list[vapply(prior_list, .bt_R2_is_budget_prior, logical(1)) &
    names(prior_list) == total_name]
  if(length(terms) > 0L){
    shapes <- terms[[1L]]$shapes
    n_units <- terms[[1L]]$n_units
    allocation_name <- terms[[1L]]$allocation_name
  }else if(length(budget) == 1L){
    shapes <- list(a = budget[[1L]]$parameters[["alpha"]], b = budget[[1L]]$parameters[["beta"]])
    allocation_prior <- prior_list[[paste0(total_name, "_allocation")]]
    shapes$alpha <- if(is.null(allocation_prior)) 1 else allocation_prior$parameters[["alpha"]]
    n_units <- length(shapes$alpha)
    allocation_name <- paste0(total_name, "_allocation")
  }else{
    return(NULL)
  }
  draws <- matrix(stats::rbeta(n, shapes$a, shapes$b), ncol = 1L, dimnames = list(NULL, total_name))
  if(n_units > 1L){
    gammas <- matrix(vapply(shapes$alpha, function(alpha) stats::rgamma(n, shape = alpha), numeric(n)), nrow = n)
    allocation <- gammas / rowSums(gammas)
    colnames(allocation) <- paste0(allocation_name, "[", seq_len(n_units), "]")
    colnames(gammas) <- paste0(.JAGS_prior_dirichlet_eta_name(allocation_name), "[", seq_len(n_units), "]")
    draws <- cbind(draws, allocation, gammas)
  }
  for(term_name in names(terms)){
    directions <- matrix(.bt_R2_direction_rng(terms[[term_name]]$law, n * terms[[term_name]]$d), nrow = n)
    colnames(directions) <- .bt_R2_direction_coordinates(terms[[term_name]], term_name)
    draws <- cbind(draws, directions)
  }
  assign(total_name, draws, envir = registry)
  draws
}

.bt_R2_prior_sample_columns <- function(prior, param_name, prior_list, n, registry){

  if(.bt_R2_is_budget_prior(prior)){
    role <- attr(prior, "R2_budget", exact = TRUE)$role
    total_name <- if(identical(role, "allocation")) sub("_allocation$", "", param_name) else param_name
    draws <- .bt_R2_budget_registry_draws(registry, prior_list, total_name, n)
    columns <- if(identical(role, "allocation")){
      K <- prior$parameters[["K"]]
      c(paste0(param_name, "[", seq_len(K), "]"),
        paste0(.JAGS_prior_dirichlet_eta_name(param_name), "[", seq_len(K), "]"))
    }else{
      param_name
    }
    return(draws[, columns, drop = FALSE])
  }
  .bt_R2_require_bound(prior, param_name)
  draws <- .bt_R2_budget_registry_draws(registry, prior_list, prior$total_name, n)
  if(!all(.bt_R2_direction_coordinates(prior, param_name) %in% colnames(draws))){
    directions <- matrix(.bt_R2_direction_rng(prior$law, n * prior$d), nrow = n)
    colnames(directions) <- .bt_R2_direction_coordinates(prior, param_name)
    draws <- cbind(draws, directions)
    assign(prior$total_name, draws, envir = registry)
  }
  node <- .bt_dnode_R2_coefficient(param_name, prior)
  values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(draws, .bt_R2_budget_prior_list(prior)))
  cbind(values, draws[, .bt_R2_direction_coordinates(prior, param_name), drop = FALSE])
}
