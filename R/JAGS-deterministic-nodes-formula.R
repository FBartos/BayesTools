# Formula linear-predictor deterministic node family.


# Formula linear predictors ("linear_predictor") ------------------------------
#
# JAGS_formula() defines each formula parameter 'p' as
#   for(i in 1:N_p){ p[i] = <term> + <term> + ... }
# with, in this order,
#   the intercept    p_intercept, or log(p_intercept) for log(intercept)
#                    formulas (a mean-centered random intercept replaces it by
#                    its group location),
#   continuous terms [m *] p_<term> * p_data_<term>[i],
#   factor terms     [m *] inprod(p_<term>, p_data_<term>[i,]),
#   expression() terms,
#   sampled random-effect blocks  p__xREx__<block>[i] (or the group location
#                    'p__xREx__<block>_xRE_MEANx[...]' of a mean-centered block),
# where m is the term prior's 'multiply_by' (a number or a parameter). The
# fixed terms are evaluated from one term specification per intercept or model
# term, with the arithmetic of the JAGS terms: coefficient times multiplier
# times data column, inner products of factor coefficients, and the terms
# summed in model order.

.bt_dnode_linear_predictor_term <- function(parameter, model_term, type,
                                            columns, prior,
                                            log = FALSE){

  list(
    parameter = parameter,
    model_term = model_term,
    type = type,
    name = paste0(parameter, "_", model_term),
    data_name = paste0(parameter, "_data_", model_term),
    columns = columns,
    prior = prior,
    multiply_by = attr(prior, "multiply_by", exact = TRUE),
    log = isTRUE(log)
  )
}

# Coordinate names of a fixed term's coefficients.
.bt_dnode_linear_predictor_coefficient_names <- function(term){

  if(identical(term$type, "factor") && .get_prior_factor_levels(term$prior) > 1L){
    return(paste0(term$name, "[", seq_len(.get_prior_factor_levels(term$prior)), "]"))
  }

  term$name
}

# JAGS syntax of one fixed term.
.bt_dnode_linear_predictor_term_syntax <- function(term){

  if(identical(term$type, "intercept")){
    # The intercept term does not emit its prior's 'multiply_by'.
    return(if(isTRUE(term$log)) paste0("log(", term$name, ")") else term$name)
  }
  multiplier <- if(!is.null(term$multiply_by)) paste0(term$multiply_by, " * ")
  if(identical(term$type, "continuous")){
    return(paste0(multiplier, term$name, " * ", term$data_name, "[i]"))
  }

  paste0(multiplier, "inprod(", term$name, ", ", term$data_name, "[i,])")
}

# JAGS syntax of the linear predictor of formula parameter 'parameter' from its
# ordered term syntax.
.bt_dnode_linear_predictor_syntax <- function(parameter, terms_syntax){

  paste0(
    "for(i in 1:N_", parameter, "){\n",
    "  ", parameter, "[i] = ", paste0(terms_syntax, collapse = " + "), "\n",
    "}\n"
  )
}

# The fixed terms of a formula parameter in model order from a fitted design
# and the formula prior list.
.bt_dnode_linear_predictor_design_terms <- function(design, prior_list,
                                                    log_intercept = design$log_intercept,
                                                    require_priors = TRUE){

  parameter <- design$parameter
  terms <- list()
  intercept_name <- paste0(parameter, "_intercept")
  if(0L %in% design$assign && (require_priors || intercept_name %in% names(prior_list))){
    prior <- prior_list[[intercept_name]]
    if(is.null(prior)){
      stop(
        "Stored formula design for parameter '", parameter,
        "' is missing the prior of its intercept.",
        call. = FALSE
      )
    }
    terms[[1L]] <- .bt_dnode_linear_predictor_term(
      parameter = parameter,
      model_term = "intercept",
      type = "intercept",
      columns = which(design$assign == 0L),
      prior = prior,
      log = isTRUE(log_intercept)
    )
  }
  for(i in seq_along(design$model_terms)){
    model_term <- design$model_terms[[i]]
    if(identical(model_term, "intercept")){
      next
    }
    prior <- prior_list[[paste0(parameter, "_", model_term)]]
    if(is.null(prior) && !require_priors){
      next
    }
    if(is.null(prior)){
      stop(
        "Stored formula design for parameter '", parameter,
        "' is missing the prior of model term '", model_term, "'.",
        call. = FALSE
      )
    }
    terms[[length(terms) + 1L]] <- .bt_dnode_linear_predictor_term(
      parameter = parameter,
      model_term = model_term,
      type = if(identical(design$model_terms_type[[i]], "factor")) "factor" else "continuous",
      columns = .bt_JAGS_formula_design_term_columns(design, model_term),
      prior = prior
    )
  }

  terms
}

# Whether a formula design defines its linear predictor node: a fitted design
# that carries the priors of its intercept and all its model terms.
.bt_dnode_linear_predictor_is_defined <- function(design){

  if(!.bt_JAGS_formula_design_can_reconstruct(design) || !is.list(design$prior_list)){
    return(FALSE)
  }
  model_terms <- setdiff(design$model_terms, "intercept")
  term_names <- character()
  if(length(model_terms) > 0L){
    term_names <- paste0(design$parameter, "_", model_terms)
  }
  if(0L %in% design$assign){
    term_names <- c(term_names, paste0(design$parameter, "_intercept"))
  }

  all(term_names %in% names(design$prior_list))
}

# The linear predictor node of a fitted formula design. Its coordinates are the
# rows of the fitted data; its dependencies are the fixed coefficients, the
# coefficient multipliers, the parameters of expression terms, and for every
# sampled random-effect block the standardized latent effects and the SD and
# correlation monitors the block contribution is reconstructed from.
.bt_dnode_linear_predictor <- function(design){

  parameter <- design$parameter
  n_rows <- nrow(design$model_matrix)
  terms <- .bt_dnode_linear_predictor_design_terms(design, design$prior_list)
  random_terms <- .bt_formula_design_sampled_random_effects(design)
  mean_translated <- vapply(random_terms, function(random_term){
    !is.null(random_term$mean_translation)
  }, logical(1))

  terms_syntax <- vapply(terms, .bt_dnode_linear_predictor_term_syntax, character(1))
  if(any(mean_translated)){
    # A mean-centered random intercept replaces the fixed intercept.
    terms_syntax <- terms_syntax[vapply(terms, `[[`, character(1), "type") != "intercept"]
  }
  expressions <- design$transformed_terms
  terms_syntax <- c(
    terms_syntax,
    vapply(expressions, .clean_from_expression, character(1), USE.NAMES = FALSE),
    vapply(random_terms, function(random_term){
      translation <- random_term$mean_translation
      if(!is.null(translation)){
        return(paste0(translation$location_name, "[", translation$group_map_name, "[i],1]"))
      }
      paste0(random_term$parameter_stem, "[i]")
    }, character(1))
  )

  coefficient_dependencies <- unlist(lapply(terms, function(term){
    if(is.prior.point(term$prior)){
      return(character())
    }
    .bt_dnode_linear_predictor_coefficient_names(term)
  }), use.names = FALSE)
  random_dependencies <- unlist(lapply(random_terms, .bt_dnode_linear_predictor_random_dependencies), use.names = FALSE)

  .bt_deterministic_node(
    family = "linear_predictor",
    node = parameter,
    coordinates = paste0(parameter, "[", seq_len(n_rows), "]"),
    dependencies = c(
      coefficient_dependencies,
      .bt_formula_predictor_multiplier_dependencies(design),
      unlist(lapply(design$expression_specs, `[[`, "parameter_dependencies"), use.names = FALSE),
      random_dependencies
    ),
    parameter = parameter,
    spec = list(
      terms_syntax = terms_syntax,
      sampled_blocks = vapply(random_terms, function(random_term) random_term$block_name, character(1)),
      has_random = .bt_formula_design_has_any_random_effects(design)
    )
  )
}

.bt_dnode_linear_predictor_random_dependencies <- function(random_term){

  translation <- random_term$mean_translation
  if(!is.null(translation)){
    return(paste0(translation$location_name, "[", seq_len(random_term$n_groups), ",1]"))
  }
  latent <- as.vector(.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  sd_names <- unique(random_term$sd_parameter_names)
  correlation <- random_term$correlation
  correlation_names <- if(!is.list(correlation) || random_term$n_columns < 2L){
    character()
  }else if(identical(correlation$type, "rho")){
    correlation$rho_name
  }else if(identical(correlation$type, "lkj")){
    as.vector(.bt_random_effect_cholesky_names(random_term, random_term$n_columns))
  }else{
    character()
  }

  c(latent, sd_names[!is.na(sd_names)], correlation_names)
}

.bt_dnode_linear_predictor_emit <- function(node){

  .bt_dnode_linear_predictor_syntax(node$node, node$spec$terms_syntax)
}

# The linear predictor on the fitted rows: the formula evaluator of
# JAGS_evaluate_formula() on the draws, with the random-effect contributions of
# the sampled blocks (marginalized blocks are not part of the node). It needs
# the fit, which only JAGS_evaluate_deterministic() supplies.
.bt_dnode_linear_predictor_evaluate <- function(node, lookup){

  if(is.null(lookup$fit)){
    return(NULL)
  }
  for(name in node$dependencies){
    if(is.null(.bt_deterministic_lookup_value(lookup, name))){
      return(NULL)
    }
  }
  spec <- node$spec
  sampled <- length(spec$sampled_blocks) > 0L
  values <- .bt_JAGS_evaluate_formula(
    fit = lookup$fit,
    parameter = node$node,
    formula_target = if(sampled) "conditional" else if(isTRUE(spec$has_random)) "fixed",
    blocks = if(sampled) spec$sampled_blocks,
    posterior = lookup$draws
  )

  t(unname(values))
}

# The fixed part of a linear predictor: a rows x draws matrix. 'values_of(term)'
# returns the term's coefficient draws (draws x coefficients) and
# 'multiplier_of(term)' its multiplier draws (a vector, or NULL without one).
.bt_dnode_linear_predictor_fixed <- function(terms, model_matrix, n_draws,
                                             values_of, multiplier_of){

  n_rows <- nrow(model_matrix)
  output <- matrix(0, nrow = n_rows, ncol = n_draws)
  for(term in terms){
    values <- values_of(term)
    multiplier <- multiplier_of(term)
    if(identical(term$type, "intercept")){
      value <- values[, 1L]
      if(isTRUE(term$log)){
        value <- log(value)
      }
      if(!is.null(multiplier)){
        value <- multiplier * value
      }
      contribution <- matrix(value, nrow = n_rows, ncol = n_draws, byrow = TRUE)
    }else if(identical(term$type, "continuous")){
      # (multiplier * coefficient) * data, one exact product per cell.
      coefficient <- values[, 1L]
      if(!is.null(multiplier)){
        coefficient <- multiplier * coefficient
      }
      contribution <- model_matrix[, term$columns, drop = FALSE] %*%
        matrix(coefficient, nrow = 1L)
    }else{
      contribution <- model_matrix[, term$columns, drop = FALSE] %*% t(values)
      if(!is.null(multiplier)){
        contribution <- contribution * rep(multiplier, each = n_rows)
      }
    }
    output <- output + contribution
  }

  output
}

# The fixed part and expression terms of a fitted linear predictor for one
# draw of bridge sampling or the marginal-likelihood parameters, compiled once
# per design.
.bt_dnode_linear_predictor_fixed_plan <- function(design, formula_prior_list,
                                                  log_intercept = FALSE,
                                                  context = "Bridge reconstruction"){

  force(design)
  parameter <- design$parameter
  context <- paste0(context, " for parameter '", parameter, "'")

  # Every prior of the formula prior list must be a term of the design; the
  # fixed part consists of the terms that have a prior.
  intercept_name <- paste0(parameter, "_intercept")
  if(intercept_name %in% names(formula_prior_list)){
    .bt_validate_formula_reconstruction_prior(
      formula_prior_list[[intercept_name]],
      intercept_name
    )
  }
  for(name in setdiff(names(formula_prior_list), intercept_name)){
    .bt_JAGS_formula_design_term_columns(
      design,
      sub(paste0("^", JAGS_regex_escape(parameter), "_"), "", name)
    )
    .bt_validate_formula_reconstruction_prior(formula_prior_list[[name]], name)
  }
  terms <- .bt_dnode_linear_predictor_design_terms(
    design = design,
    prior_list = formula_prior_list,
    log_intercept = isTRUE(log_intercept) || isTRUE(design$log_intercept),
    require_priors = FALSE
  )
  for(i in seq_along(terms)){
    terms[[i]]$index <- i
  }
  value_evaluators <- lapply(terms, function(term){
    names <- .bt_dnode_linear_predictor_coefficient_names(term)
    if(is.prior.point(term$prior)){
      location <- term$prior[["parameters"]][["location"]]
      n_values <- length(term$columns)
      return(function(samples) rep(location, n_values))
    }
    .bt_JAGS_bridge_compile_parameter_values(term$prior, names)
  })
  multiplier_evaluators <- lapply(terms, function(term){
    .bt_JAGS_bridge_compile_prior_multiply_by(term$prior)
  })
  has_multiplier <- vapply(terms, function(term) !is.null(term$multiply_by), logical(1))

  expressions <- design$expression_specs
  if(is.null(expressions)){
    expressions <- design$transformed_terms
  }
  expression_data <- NULL
  if(length(expressions) > 0L){
    expression_data <- .bt_formula_expression_merge_data(
      design$source_data,
      design$expression_data,
      context = context
    )
  }
  model_matrix <- design$model_matrix

  list(
    value = function(samples, prior_list_parameters){
      output <- .bt_dnode_linear_predictor_fixed(
        terms = terms,
        model_matrix = model_matrix,
        n_draws = 1L,
        values_of = function(term){
          matrix(value_evaluators[[term$index]](samples), nrow = 1L)
        },
        multiplier_of = function(term){
          if(has_multiplier[[term$index]]){
            multiplier_evaluators[[term$index]](prior_list_parameters)
          }else{
            NULL
          }
        }
      )
      output <- as.vector(output)
      if(length(expressions) > 0L){
        output <- output + .bt_formula_expression_row_values(
          expressions = expressions,
          data = expression_data,
          n_rows = nrow(model_matrix),
          context = context,
          samples = samples,
          parameters = prior_list_parameters
        )
      }
      output
    },
    terms = terms
  )
}
