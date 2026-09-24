#' @title Generated deterministic nodes
#'
#' @description \code{JAGS_deterministic_nodes()} lists the deterministic
#' nodes that BayesTools generates in the JAGS model of a fit, with the family
#' that defines them, the coordinates they produce, and the coordinates they
#' are computed from. \code{JAGS_evaluate_deterministic()} recomputes their
#' values from draws of those coordinates.
#'
#' @param fit a model fitted with [JAGS_fit()].
#' @param draws draws to evaluate the nodes on: a numeric matrix with one
#' column per coordinate (for example posterior draws with some columns
#' replaced, or prior draws from [transform_prior_samples()] with
#' \code{formula_scale = list()}), an \code{mcmc} or \code{mcmc.list} object,
#' or a named numeric vector holding a single draw. Coordinates that have a
#' point prior do not need a column. Defaults to the posterior draws of
#' \code{fit}.
#' @param nodes optional character vector of node names (the \code{node}
#' column of \code{JAGS_deterministic_nodes()}) to evaluate. By default, every
#' node whose dependencies are available in \code{draws} is evaluated.
#'
#' @details Each generated node belongs to one family, which defines both the
#' JAGS syntax of the node and the R evaluator used everywhere BayesTools needs
#' its value (prior draws, parameter-catalog quantities, bridge sampling,
#' marginal-likelihood parameters, prediction, marginal posteriors of formula
#' parameters, random-effect unscaling, and convergence roles):
#' \describe{
#'   \item{\code{"random_sd"}}{random-effect SDs derived from a variance
#'   allocation ([random_variance_allocation()]): the source SD (the
#'   allocation's own SD or an external scalar SD source) times the chain of
#'   allocation factors from the root allocation, \code{sqrt(w[i])} (total
#'   variance) or \code{sqrt(K * w[i])} (mean variance) of the Dirichlet
#'   weights and the inclusion gates. The family also contains the parent
#'   allocation components that a child allocation splits further.}
#'   \item{\code{"random_rho"}}{scalar correlations of structured random-effect
#'   blocks sampled on the Fisher-z (\code{rho = tanh(z)}) or logit scale
#'   (\code{rho = lower + (upper - lower) * plogis(z)}).}
#'   \item{\code{"lkj"}}{the Cholesky factor \code{L}, the correlation matrix
#'   \eqn{R = L L^T} (with an exact unit diagonal), and the monitored partial
#'   correlations \code{cpc = 2 u - 1} of an LKJ-Cholesky block, computed from
#'   its primitives \code{u} with the kernel of the BayesTools JAGS module.}
#'   \item{\code{"omega"}}{the publication weights of a weight-function prior
#'   (alone, as the selection of [prior_bias()], or as branches of a
#'   publication-bias mixture) on the one-sided global p-value bins, with the
#'   bins of two-sided weight functions mirrored. The dependencies are the free
#'   coordinates of the weights (\code{eta}, the monitored \code{omega[2]} of
#'   a binary cumulative weight function, or the independent \code{omega} or
#'   \code{log_omega}); mixture branches are read from their own component
#'   nodes and the \code{bias_indicator}.}
#'   \item{\code{"prior_mixture"}}{parameters with spike-and-slab
#'   (\code{p = p_variable * p_indicator}) or mixture priors
#'   (\code{p = sum_k p_component_k * (p_indicator == k)}), and the PET and
#'   PEESE terms of publication-bias mixtures
#'   (\code{PET <- PET_1 * equals(bias_indicator, k)}), per coefficient of
#'   factor priors. Only the components that are active in some draw need to
#'   be available; point components are constants.}
#'   \item{\code{"linear_predictor"}}{the linear predictor of a formula
#'   parameter on the fitted rows: the intercept (or its logarithm), the
#'   continuous and factor terms with their \code{multiply_by} multipliers,
#'   the \code{expression()} terms, and the contributions of the sampled
#'   random-effect blocks, evaluated with [JAGS_evaluate_formula()]. It is
#'   evaluated only by \code{JAGS_evaluate_deterministic()}, and only when the
#'   standardized latent effects of the sampled blocks are among the draws.
#'   A block whose SDs come from a row-shaped external source
#'   (\code{random_sd_source(..., shape = "row")}) depends on the source rows
#'   \code{source[1]}, ..., \code{source[N]}, or, when the source's
#'   \code{values} function reconstructs them, on the inputs that function
#'   declares (\code{parameter_source(inputs = )}; a function without declared
#'   inputs adds no dependencies and is evaluated whenever the other
#'   dependencies are available), and on the Dirichlet weights and inclusion
#'   gates of its allocation.}
#' }
#' Nodes are evaluated with the arithmetic of the R evaluator, which reproduces
#' the JAGS monitors exactly or to the last bits of floating-point rounding.
#'
#' @return \code{JAGS_deterministic_nodes()} returns a data frame with one row
#' per node and the columns \code{node} (node name), \code{family},
#' \code{parameter} (the formula parameter or prior name the node belongs
#' to), \code{block} (the random-effect block, \code{NA} otherwise),
#' \code{coordinates} and \code{dependencies} (list columns of coordinate
#' names), and \code{monitored} (whether all coordinates of the node are
#' monitored in \code{fit}).
#'
#' \code{JAGS_evaluate_deterministic()} returns a numeric matrix with one row
#' per draw and one column per coordinate of the evaluated nodes.
#'
#' @seealso [JAGS_fit()] [parameter_coordinates()] [transform_prior_samples()]
#' @export
JAGS_deterministic_nodes <- function(fit){

  .bt_require_fit_contract(fit, "fit")
  nodes <- .bt_deterministic_nodes_fit(fit)
  monitored_columns <- colnames(as.matrix(.fit_to_posterior(fit)))

  .bt_deterministic_node_table(nodes, monitored_columns)
}

#' @rdname JAGS_deterministic_nodes
#' @export
JAGS_evaluate_deterministic <- function(fit, draws = NULL, nodes = NULL){

  .bt_require_fit_contract(fit, "fit")
  if(is.null(draws)){
    draws <- .fit_to_posterior(fit)
  }
  draws <- .bt_deterministic_draws_matrix(draws)
  check_char(nodes, "nodes", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  if(is.null(prior_list)){
    prior_list <- list()
  }
  all_nodes <- .bt_deterministic_nodes_fit(fit)
  requested <- !is.null(nodes)
  if(requested){
    unknown <- setdiff(nodes, names(all_nodes))
    if(length(unknown) > 0L){
      stop(
        "'nodes' contains names that are not generated deterministic nodes of 'fit': ",
        paste0("'", unknown, "'", collapse = ", "),
        ". See JAGS_deterministic_nodes(fit) for the available nodes.",
        call. = FALSE
      )
    }
    all_nodes <- all_nodes[unique(nodes)]
  }

  lookup <- .bt_deterministic_lookup(draws, prior_list)
  lookup$fit <- fit
  values <- list()
  for(node in all_nodes){
    node_values <- .bt_deterministic_node_evaluate(node, lookup)
    if(is.null(node_values)){
      if(requested){
        stop(
          "Deterministic node '", node$node, "' is unavailable from 'draws': ",
          "its dependencies ",
          paste0("'", node$dependencies, "'", collapse = ", "),
          " are not all available as columns of 'draws' or as point priors.",
          call. = FALSE
        )
      }
      next
    }
    values[[length(values) + 1L]] <- node_values
  }

  if(length(values) == 0L){
    return(matrix(numeric(), nrow = nrow(draws), ncol = 0L))
  }
  out <- do.call(cbind, values)
  rownames(out) <- NULL
  out
}


# Node specifications ----------------------------------------------------------
#
# Every deterministic node that BayesTools generates in a JAGS model is
# described by one specification of a registered family. The family defines,
# from that specification alone,
#   - the JAGS syntax of the node (.bt_deterministic_node_emit()),
#   - its R evaluator over draws (.bt_deterministic_node_evaluate()),
#   - its declared dependencies and the coordinates it produces.
# The model syntax is emitted from the same specification that the evaluators
# of prior draws, catalog quantities, bridge sampling, marginal-likelihood
# parameters, prediction, and convergence roles use, so the definitions cannot
# drift apart.

.bt_deterministic_node_families <- c(
  "random_sd",
  "random_rho",
  "lkj",
  "omega",
  "prior_mixture",
  "linear_predictor"
)

.bt_deterministic_node <- function(family, node, coordinates,
                                   dependencies = character(),
                                   parameter = NA_character_,
                                   block = NA_character_,
                                   spec = list()){

  check_char(family, "family", allow_values = .bt_deterministic_node_families,
             allow_NA = FALSE)
  check_char(node, "node", allow_NA = FALSE)
  check_char(coordinates, "coordinates", check_length = 0, allow_NA = FALSE)
  check_char(dependencies, "dependencies", check_length = 0, allow_NULL = TRUE,
             allow_NA = FALSE)
  if(anyDuplicated(coordinates)){
    stop("Deterministic node '", node, "' must define unique coordinates.",
         call. = FALSE)
  }

  out <- list(
    family       = family,
    node         = node,
    coordinates  = coordinates,
    dependencies = unique(dependencies),
    parameter    = as.character(parameter),
    block        = as.character(block),
    spec         = spec
  )
  class(out) <- c("BayesTools_deterministic_node", "list")
  out
}

# JAGS syntax lines of a node.
.bt_deterministic_node_emit <- function(node){

  switch(
    node$family,
    random_sd = .bt_dnode_random_sd_emit(node),
    random_rho = .bt_dnode_rho_emit(node),
    lkj = .bt_dnode_lkj_emit(node),
    omega = .bt_dnode_omega_emit(node),
    prior_mixture = .bt_dnode_prior_mixture_emit(node),
    linear_predictor = .bt_dnode_linear_predictor_emit(node),
    stop("Unsupported deterministic node family '", node$family, "'.", call. = FALSE)
  )
}

# Values of a node on the draws of a lookup: a matrix with one row per draw
# and one column per coordinate, or NULL when a dependency is unavailable.
.bt_deterministic_node_evaluate <- function(node, lookup){

  values <- switch(
    node$family,
    random_sd = .bt_dnode_random_sd_evaluate(node, lookup),
    random_rho = .bt_dnode_rho_evaluate(node, lookup),
    lkj = .bt_dnode_lkj_evaluate(node, lookup),
    omega = .bt_dnode_omega_evaluate(node, lookup),
    prior_mixture = .bt_dnode_prior_mixture_evaluate(node, lookup),
    linear_predictor = .bt_dnode_linear_predictor_evaluate(node, lookup),
    stop("Unsupported deterministic node family '", node$family, "'.", call. = FALSE)
  )
  if(is.null(values)){
    return(NULL)
  }
  values <- matrix(values, nrow = lookup$n, ncol = length(node$coordinates))
  colnames(values) <- node$coordinates
  values
}

# The generated deterministic nodes of a fit, keyed by node name.
.bt_deterministic_nodes_fit <- function(fit){

  .bt_deterministic_nodes(
    prior_list = attr(fit, "prior_list", exact = TRUE),
    formula_design = attr(fit, "formula_design", exact = TRUE)
  )
}

.bt_deterministic_nodes <- function(prior_list = NULL, formula_design = NULL){

  nodes <- .bt_deterministic_nodes_prior_list(prior_list)
  if(inherits(formula_design, "BayesTools_formula_design")){
    formula_design <- list(formula_design)
  }
  if(is.list(formula_design)){
    for(design in formula_design){
      for(random_term in .bt_formula_design_random_effects(design)){
        nodes <- c(nodes, .bt_deterministic_nodes_random_term(
          random_term,
          parameter = design$parameter
        ))
      }
      nodes <- c(nodes, .bt_dnode_random_sd_components(design))
      if(.bt_dnode_linear_predictor_is_defined(design)){
        nodes[[length(nodes) + 1L]] <- .bt_dnode_linear_predictor(design)
      }
    }
  }

  names(nodes) <- vapply(nodes, function(node) node$node, character(1))
  if(anyDuplicated(names(nodes))){
    stop(
      "Generated deterministic nodes must have unique names; duplicated: ",
      paste0("'", unique(names(nodes)[duplicated(names(nodes))]), "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  nodes
}

# Nodes that the priors of a prior list define.
.bt_deterministic_nodes_prior_list <- function(prior_list){

  nodes <- list()
  prior_names <- names(prior_list)
  for(i in seq_along(prior_list)){
    prior <- prior_list[[i]]
    prior_nodes <- if(is.prior.weightfunction(prior) || is_prior_bias(prior)){
      list(.bt_dnode_omega(prior_names[[i]], prior))
    }else if(inherits(prior, "prior.bias_mixture")){
      list(
        .bt_dnode_omega(prior_names[[i]], prior),
        .bt_dnode_prior_mixture_bias_term(prior_names[[i]], prior, "PET"),
        .bt_dnode_prior_mixture_bias_term(prior_names[[i]], prior, "PEESE")
      )
    }else if(is.prior.mixture(prior)){
      list(.bt_dnode_prior_mixture(prior_names[[i]], prior))
    }else{
      list()
    }
    for(node in prior_nodes){
      if(!is.null(node)){
        nodes[[length(nodes) + 1L]] <- node
      }
    }
  }

  nodes
}

.bt_deterministic_nodes_random_term <- function(random_term,
                                                parameter = NA_character_){

  nodes <- .bt_dnode_random_sd_from_random_term(random_term, parameter = parameter)
  for(node in list(
    .bt_dnode_rho_from_random_term(random_term, parameter = parameter),
    .bt_dnode_lkj_from_random_term(random_term, parameter = parameter)
  )){
    if(!is.null(node)){
      nodes[[length(nodes) + 1L]] <- node
    }
  }

  nodes
}

# The parents of the generated nodes by JAGS node name (coordinate base): the
# bases of their declared dependencies. The convergence roles read them.
.bt_deterministic_node_parent_bases <- function(nodes){

  parents <- list()
  for(node in nodes){
    dependency_bases <- unique(.bt_parameter_coordinates_base(node$dependencies))
    for(base in unique(.bt_parameter_coordinates_base(node$coordinates))){
      parents[[base]] <- setdiff(
        unique(c(parents[[base]], dependency_bases)),
        base
      )
    }
  }

  parents
}

.bt_deterministic_node_table <- function(nodes, monitored_columns = character()){

  nodes <- unname(nodes)
  out <- data.frame(
    node = vapply(nodes, function(node) node$node, character(1)),
    family = vapply(nodes, function(node) node$family, character(1)),
    parameter = vapply(nodes, function(node) node$parameter, character(1)),
    block = vapply(nodes, function(node) node$block, character(1)),
    stringsAsFactors = FALSE
  )
  out$coordinates <- lapply(nodes, function(node) node$coordinates)
  out$dependencies <- lapply(nodes, function(node) node$dependencies)
  out$monitored <- vapply(nodes, function(node){
    all(node$coordinates %in% monitored_columns)
  }, logical(1))
  rownames(out) <- NULL
  out
}


# Draws and their lookup ---------------------------------------------------------

.bt_deterministic_draws_matrix <- function(draws){

  if(inherits(draws, "mcmc.list") || inherits(draws, "mcmc")){
    draws <- as.matrix(draws)
  }else if(is.data.frame(draws)){
    draws <- as.matrix(draws)
  }else if(is.numeric(draws) && is.null(dim(draws))){
    if(is.null(names(draws))){
      stop("A single draw in 'draws' must be a named numeric vector.", call. = FALSE)
    }
    draws <- matrix(draws, nrow = 1L, dimnames = list(NULL, names(draws)))
  }
  if(!is.matrix(draws) || !is.numeric(draws) || is.null(colnames(draws)) ||
     anyNA(colnames(draws)) || anyDuplicated(colnames(draws))){
    stop(
      "'draws' must be a numeric matrix, 'mcmc', or 'mcmc.list' object with ",
      "unique column names, or a named numeric vector.",
      call. = FALSE
    )
  }

  draws
}

# The values a node is computed from: draws with coordinate columns, and the
# prior list whose point priors supply constant coordinates without a column.
.bt_deterministic_lookup <- function(draws, prior_list = list()){

  if(is.null(prior_list)){
    prior_list <- list()
  }
  list(
    draws = draws,
    prior_list = prior_list,
    n = nrow(draws)
  )
}

# A lookup of one draw of bridge sampling or the marginal-likelihood
# parameters: a named numeric vector or list of coordinates.
.bt_deterministic_row_lookup <- function(samples, prior_list = list()){

  if(is.list(samples)){
    samples <- unlist(samples)
  }
  .bt_deterministic_lookup(
    .bt_JAGS_marglik_random_effect_posterior_row(samples),
    prior_list
  )
}

# Draws of one scalar coordinate: its column, or the location of its point
# prior; NULL when neither is available.
.bt_deterministic_lookup_value <- function(lookup, name){

  .bt_random_effect_parameter_draws(
    parameter_name = name,
    posterior = lookup$draws,
    prior_list = lookup$prior_list
  )
}

# Draws of the Dirichlet weights of an allocation factor (a draws x K matrix,
# from the weight columns or the normalized auxiliary gamma draws); NULL when
# unavailable.
.bt_deterministic_lookup_simplex <- function(lookup, factor){

  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = factor$weight_name,
    posterior = lookup$draws,
    prior_list = lookup$prior_list
  )
  if(is.null(weights)){
    return(NULL)
  }
  if(ncol(weights) != factor$n_targets){
    stop(
      "Random-effect allocation factor metadata for '",
      factor$weight_name,
      "' expected ", factor$n_targets,
      " Dirichlet coordinate(s), but found ", ncol(weights), ".",
      call. = FALSE
    )
  }
  if(factor$index > ncol(weights)){
    stop(
      "Random-effect allocation factor metadata for '",
      factor$weight_name,
      "' reference coordinate ", factor$index,
      ", but only ", ncol(weights), " coordinate(s) are available.",
      call. = FALSE
    )
  }

  weights
}

# Draws of several coordinates as a matrix with one column per name; NULL when
# any of them is unavailable.
.bt_deterministic_lookup_values <- function(lookup, names){

  indices <- match(names, colnames(lookup$draws))
  if(!anyNA(indices)){
    return(lookup$draws[, indices, drop = FALSE])
  }
  values <- matrix(NA_real_, nrow = lookup$n, ncol = length(names))
  for(i in seq_along(names)){
    value <- .bt_deterministic_lookup_value(lookup, names[[i]])
    if(is.null(value)){
      return(NULL)
    }
    values[, i] <- value
  }
  colnames(values) <- names
  values
}
