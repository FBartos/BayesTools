# Marginal-likelihood and bridge-sampling support of R-squared budget terms.
#
# The bridge coordinates of a budget term are its direction coordinates;
# the budget total and allocation are ordinary Beta and Dirichlet priors of
# the formula and contribute their own densities. The coefficients are the
# registered 'R2_coefficient' node evaluated on the draw.

.bt_R2_direction_values <- function(samples, prior, parameter_name, signal = TRUE){

  names_z <- .bt_R2_direction_coordinates(prior, parameter_name)
  if(!all(names_z %in% names(samples))){
    if(signal){
      .bt_JAGS_marglik_missing_columns("'samples' does not contain all monitored R-squared budget direction coordinates.")
    }
    return(NULL)
  }
  unname(unlist(samples[names_z], use.names = FALSE))
}

.JAGS_marglik_priors.R2 <- function(samples, prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  z <- .bt_R2_direction_values(samples, prior, parameter_name)
  .bt_R2_direction_lpdf(prior$law, z)
}

.bt_JAGS_marglik_compile_R2_prior_rows <- function(prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  names_z <- .bt_R2_direction_coordinates(prior, parameter_name)
  law <- prior$law
  function(samples){
    if(!all(names_z %in% colnames(samples)))
      .bt_JAGS_marglik_missing_columns("'samples' does not contain all monitored R-squared budget direction coordinates.")
    z <- as.matrix(samples[, names_z, drop = FALSE])
    if(identical(law, "laplace")){
      rowSums(log(sqrt(2) / 2) - sqrt(2) * abs(z))
    }else{
      rowSums(stats::dnorm(z, log = TRUE))
    }
  }
}

# Coefficient values of a term on one draw (a named vector or list of
# coordinates), from the budget total, the allocation (or its gamma
# coordinates), and the direction coordinates.
.bt_JAGS_marglik_compile_R2_parameter_values <- function(prior, parameter_name){

  node <- .bt_dnode_R2_coefficient(parameter_name, prior)
  budget_priors <- .bt_R2_budget_prior_list(prior)
  function(samples){
    lookup <- .bt_deterministic_row_lookup(samples, budget_priors)
    values <- .bt_deterministic_node_evaluate(node, lookup)
    if(is.null(values)){
      .bt_JAGS_marglik_missing_columns("'samples' does not contain all monitored R-squared budget parameters.")
    }
    as.vector(values)
  }
}

.JAGS_marglik_parameters.R2 <- function(samples, prior, parameter_name){

  parameter <- list()
  parameter[[parameter_name]] <- .bt_JAGS_marglik_compile_R2_parameter_values(prior, parameter_name)(samples)
  parameter
}

.bt_JAGS_bridge_compile_R2_evaluator <- function(prior, parameter_name){

  force(prior)
  force(parameter_name)
  .bt_R2_require_bound(prior, parameter_name)
  parameter_values <- .bt_JAGS_marglik_compile_R2_parameter_values(prior, parameter_name)
  list(
    log_prior = function(samples){
      .JAGS_marglik_priors.R2(samples, prior, parameter_name)
    },
    parameters = function(samples){
      parameter <- list()
      parameter[[parameter_name]] <- parameter_values(samples)
      parameter
    }
  )
}

.bt_JAGS_bridge_owned_R2_names <- function(prior, parameter_name){

  c(parameter_name,
    .bt_R2_coefficient_names(prior, parameter_name),
    .bt_R2_direction_name(parameter_name),
    .bt_R2_direction_coordinates(prior, parameter_name))
}
