.bt_hypothesis_linear_target_name <- ".BayesTools_linear_target"

#' Compile a hypothesis linear in the levels of one parameter
#'
#' @description
#' When structured quantity labels are available, hypothesis tables display
#' the compiled combination's label rather than its internal scalar name.
#' Compiles point and simple region statements that all refer to one linear
#' combination \eqn{t = \sum_l c_l L_l} of levels \eqn{L_l} of the same
#' parameter of one marginal posterior, e.g. a difference of two levels, a
#' scaled level, or a weighted average of several levels. The returned scalar
#' marginal posterior of \eqn{t} retains the exact joint-prior density and the
#' structural conditioning metadata needed by downstream inference methods
#' (e.g. precomputed posterior-ordinate estimators that need the combined
#' linear weights).
#'
#' Every referenced level must have one fixed named linear-weight row over the
#' fitted coordinates and a declared posterior-atom status, and the levels
#' must share a joint prior context and conditioning event. The prior density
#' of \eqn{t} is the joint-context density at the combined weights
#' \eqn{\sum_l c_l w_l} shifted by the combined level offsets: on the fitted
#' scale it equals [parameter_prior_density()] at those coordinate weights,
#' and on the original scale of scaled formula coefficients
#' [JAGS_formula_prior_density()] with the corresponding `weights`. It must be
#' atom-free; by absolute continuity the posterior of \eqn{t} is then
#' atom-free as well. Affine transformations already applied to the marginal
#' levels retain their scale and offset; nonlinear transformed marginal levels,
#' nonlinear expressions, row-varying level weights, and statements about
#' different combinations are rejected.
#'
#' @param posterior a factor-like marginal posterior containing the referenced
#'   levels.
#' @param hypothesis hypothesis text or a `BayesTools_hypothesis_ast`. Each
#'   side is a point statement (`=`, `!=`) or a simple comparison (`<`, `<=`,
#'   `>`, `>=`) whose expressions are linear in the level references: sums and
#'   differences of level references and numbers, level references multiplied
#'   or divided by numbers, and parentheses.
#' @param parameter exact parameter root used by the level references.
#'
#' @return A list with `posterior`, `hypothesis`, `parameter`, and `weights`.
#'   `posterior` is a scalar `marginal_posterior` of the linear combination
#'   whose draw metadata ([posterior_metadata()]) hold the combined
#'   `linear_weights` and `linear_offset`, its `prior_density` and
#'   `prior_context`, declared (empty) `atoms`, the `condition` of the
#'   levels, and, when the levels carry label parts, `quantities` labelling
#'   the target by the combination of the levels' labels (e.g.
#'   `(mu) 2*f[A] - f[C]`, see [parameter_labels()]); `hypothesis` is an equivalent AST written against that scalar
#'   target, `parameter` its name, and `weights` the combined linear weights.
#'   A target that cannot be certified stops with an error of class
#'   `BayesTools_linear_target_unavailable` (also
#'   `BayesTools_hypothesis_target`) whose field `reason` is
#'   `"posterior_atoms"` (the target's prior is not structurally atom-free),
#'   `"atom_declarations"` (a level lacks its posterior-atom declaration), or
#'   `"prior_context"` (no valid joint prior context, including level weights
#'   on columns the joint prior context does not contain). A reference to a
#'   level that `posterior` does not contain stops with an error of class
#'   `BayesTools_parameter_not_found` (also
#'   `BayesTools_parameter_resolution_error`), as in [hypothesis_BF()].
#'
#' @export
hypothesis_linear_target <- function(posterior, hypothesis, parameter){

  if(!.hypothesis_inherits_marginal_posterior(posterior) ||
     !is.list(posterior)){
    stop("'posterior' must be a factor-like marginal posterior.", call. = FALSE)
  }
  check_char(parameter, "parameter", check_length = 1L, allow_NA = FALSE)
  ast <- if(inherits(hypothesis, "BayesTools_hypothesis_ast")){
    .bt_validate_hypothesis_ast(hypothesis)
    hypothesis
  }else{
    hypothesis_parse(hypothesis)
  }

  target_name <- .bt_hypothesis_linear_target_name
  sides <- unlist(lapply(ast$statements, function(statement){
    list(statement$left, statement$right)
  }), recursive = FALSE)
  forms <- lapply(sides, .hypothesis_linear_target_side, parameter = parameter)
  coefficients <- forms[[1L]]$coefficients
  same_target <- vapply(forms, function(form){
    identical(form$coefficients, coefficients)
  }, logical(1))
  if(!all(same_target)){
    stop(
      "Hypothesis statements must all use the same linear combination of levels.",
      call. = FALSE
    )
  }
  if(length(coefficients) == 0L){
    .hypothesis_linear_target_stop("The linear target prior is a structural point measure.", "posterior_atoms")
  }

  levels <- .hypothesis_match_level_names(
    levels    = names(coefficients),
    available = names(posterior),
    parameter = parameter
  )
  if(anyDuplicated(levels)){
    stop("A linear target must reference each level of '", parameter,
         "' once.", call. = FALSE)
  }
  names(coefficients) <- levels
  numerator_coefficients <- forms[[1L]]$numerator_coefficients
  names(numerator_coefficients) <- levels

  .hypothesis_validate_level_conditionals(posterior, parameter, levels)
  level_draws <- lapply(levels, function(level){
    values <- as.numeric(posterior[[level]])
    if(length(values) == 0L || any(!is.finite(values))){
      stop("Finite posterior draws are required for level '", level, "'.",
           call. = FALSE)
    }
    values
  })
  names(level_draws) <- levels
  draw_lengths <- vapply(level_draws, length, integer(1))
  if(length(unique(draw_lengths)) != 1L){
    stop("Level comparisons require equal-length posterior draws.",
         call. = FALSE)
  }

  level_weights <- lapply(levels, function(level){
    weights <- .hypothesis_level_linear_weights(posterior[[level]])
    if(is.null(weights)){
      .hypothesis_linear_target_stop(
        paste0("Linear prior weights are missing for level '", level, "'."),
        "prior_context"
      )
    }
    weights <- .hypothesis_prepare_level_weights(weights)
    if(nrow(weights) != 1L){
      stop("Linear targets require one fixed linear-weight row per level.",
           call. = FALSE)
    }
    weights[1L, ]
  })
  names(level_weights) <- levels
  level_spaces <- vapply(posterior[levels], .bt_linear_weight_space, character(1))
  if(length(unique(level_spaces)) != 1L){
    .hypothesis_linear_target_stop("Linear target levels use incompatible weight spaces.", "weight_space")
  }

  context <- tryCatch(
    .hypothesis_child_prior_context(posterior, levels),
    error = function(e){
      .hypothesis_linear_target_stop(conditionMessage(e), "prior_context")
    }
  )
  if(is.null(context)){
    context <- .bt_meta_get(posterior, "prior_context")
  }
  if(is.null(context) || !.hypothesis_is_prior_density_context(context)){
    .hypothesis_linear_target_stop(
      "A valid joint prior context is required for a linear target.",
      "prior_context"
    )
  }
  # level weights the joint context cannot evaluate (e.g. columns missing
  # from it) leave no valid joint prior context for the target
  tryCatch(
    .hypothesis_validate_level_weights_context(level_weights, context),
    error = function(e){
      .hypothesis_linear_target_stop(conditionMessage(e), "prior_context")
    }
  )

  weights <- .hypothesis_linear_target_combine_weights(
    level_weights,
    coefficients
  )
  if(length(weights) == 0L || all(weights == 0)){
    .hypothesis_linear_target_stop("The linear target prior is a structural point measure.", "posterior_atoms")
  }
  numerator_offset <- 0
  for(level in levels){
    numerator_offset <- .hypothesis_affine_sum(numerator_offset,
      .hypothesis_affine_product(numerator_coefficients[[level]],
        .hypothesis_level_linear_offset(posterior[[level]])))
  }
  offset <- .hypothesis_linear_target_divide(numerator_offset, forms[[1L]]$divisor)
  prior_density <- .prior_density_from_context(
    context, weights,
    output_transformation = if(offset != 0) "lin" else NULL,
    output_transformation_arguments = if(offset != 0) list(a = offset, b = 1) else NULL
  )
  if(!.hypothesis_linear_target_prior_atom_free(prior_density)){
    .hypothesis_linear_target_stop(
      "The linear target prior is not structurally atom-free.",
      "posterior_atoms"
    )
  }
  declared_atoms <- vapply(
    posterior[levels],
    .hypothesis_linear_target_posterior_atoms_declared,
    logical(1)
  )
  if(!all(declared_atoms)){
    .hypothesis_linear_target_stop(
      paste0("Complete structural posterior-atom declarations are required for ",
             "a linear target."),
      "atom_declarations"
    )
  }

  values <- rep(0, draw_lengths[[1L]])
  for(level in levels){
    values <- .hypothesis_affine_sum(values,
      .hypothesis_affine_product(numerator_coefficients[[level]], level_draws[[level]]))
  }
  numerator_values <- values
  numerator_weights <- .hypothesis_linear_target_combine_weights(level_weights, numerator_coefficients)
  values <- values / forms[[1L]]$divisor
  if(any(!is.finite(values)) || any(numerator_values != 0 & values == 0)){
    .hypothesis_numerical_stop("coefficient_range", "public linear target values")
  }
  class(values) <- c(
    "numeric", "marginal_posterior.linear_target", "marginal_posterior"
  )
  attr(values, "parameter")             <- target_name
  values <- .bt_meta_assign(values, list(linear_weights = weights,
    linear_offset = offset, prior_density = prior_density, prior_context = context,
    linear_weight_space = if(identical(context$linear_weight_space, "formula_contribution"))
      "formula_contribution" else "coefficient"))
  values <- .bt_meta_set(values, "hypothesis_evaluation", list(
    numerator = numerator_values, divisor = forms[[1L]]$divisor,
    weights = numerator_weights, offset = numerator_offset))
  values <- .bt_meta_set(values, "atoms", .posterior_atoms_new(
    column_names = target_name,
    source       = "linear_target"
  ))
  values <- .posterior_support_set(values, .posterior_support_from_prior_context_weights(
    context, weights, "lin", list(a = offset, b = 1)))
  values <- .posterior_components_set(values, .hypothesis_linear_components(
    posterior[levels], context, weights, list(a = offset, b = 1), draw_lengths[[1L]]))
  condition <- .bt_meta_get(posterior[[levels[[1L]]]], "condition")
  condition <- condition[intersect(names(condition), c(
    "conditional", "conditional_rule", "condition_key", "condition_event",
    "resolved_condition_event", "averaged"
  ))]
  values <- .bt_meta_set(values, "condition", if(length(condition) > 0L) condition)
  values <- .bt_meta_set(values, "quantities", .hypothesis_linear_target_quantities(
    posterior    = posterior,
    levels       = levels,
    coefficients = coefficients,
    target_name  = target_name
  ))

  rewritten <- .hypothesis_linear_target_rewrite(
    ast         = ast,
    forms       = forms,
    target_name = target_name
  )

  return(list(
    posterior  = values,
    hypothesis = rewritten,
    parameter  = target_name,
    weights    = weights
  ))
}


# The column table of a linear target: no catalog quantity, and no fitted
# coordinates declared (its levels may be estimated marginal means, which are
# predictions), labelled by the combination of the rendered labels of its
# levels (`2*f[A] - f[C]`) under their common formula parameter. NULL when a
# level carries no label parts.
.hypothesis_linear_target_quantities <- function(posterior, levels,
                                                 coefficients, target_name){

  parts <- lapply(levels, function(level){
    quantities <- .bt_draws_quantities(posterior[[level]])
    if(is.null(quantities) || nrow(quantities) != 1L){
      return(NULL)
    }
    quantities$label_parts[[1L]]
  })
  if(length(parts) == 0L || any(vapply(parts, is.null, logical(1)))){
    return(NULL)
  }
  formula_parameters <- unique(vapply(parts, `[[`, character(1), "formula_parameter"))
  formula_parameter <- if(length(formula_parameters) == 1L) formula_parameters else ""
  labels <- .bt_label(parts, style = "table", formula_prefix = !nzchar(formula_parameter))
  terms <- vapply(seq_along(levels), function(i){
    coefficient <- coefficients[[levels[[i]]]]
    magnitude <- abs(coefficient)
    term <- if(magnitude == 1){
      labels[[i]]
    }else{
      paste0(format(magnitude, digits = 15), "*", labels[[i]])
    }
    if(i == 1L){
      paste0(if(coefficient < 0) "-", term)
    }else{
      paste0(if(coefficient < 0) " - " else " + ", term)
    }
  }, character(1))

  .bt_draws_quantity_table(
    columns      = target_name,
    quantity_ids = "",
    dependencies = list(character()),
    weights      = list(numeric()),
    label_parts  = list(.bt_label_parts(
      components        = paste(terms, collapse = ""),
      formula_parameter = formula_parameter,
      selector          = target_name
    ))
  )
}


# A linear target that cannot be certified: class
# BayesTools_linear_target_unavailable (parent BayesTools_hypothesis_target)
# with the condition field 'reason': "posterior_atoms" (the target's prior is
# not structurally atom-free, so its posterior may have atoms),
# "atom_declarations" (a referenced level lacks its posterior-atom
# declaration), or "prior_context" (no valid joint prior context). Callers
# match the class and reason, never the message.
.hypothesis_linear_target_stop <- function(message, reason){

  stop(structure(
    class = c("BayesTools_linear_target_unavailable",
              "BayesTools_hypothesis_target", "error", "condition"),
    list(message = message, call = NULL, reason = reason)
  ))
}


# One side of a statement as 'coefficients %*% levels  <operator>  value':
# point sides move the expression's constant to the value, comparison sides
# subtract the right expression from the left.
.hypothesis_linear_target_side <- function(side, parameter){

  if(side$type %in% c("point", "not_point")){
    form <- .hypothesis_linear_target_form(side$expression, parameter)
    return(list(
      operator     = if(side$type == "point") "=" else "!=",
      coefficients = .hypothesis_linear_target_scaled_coefficients(form),
      numerator_coefficients = form$coefficients,
      value        = .hypothesis_linear_target_divide(.hypothesis_affine_null(form, side$value), form$divisor),
      divisor      = form$divisor
    ))
  }
  if(!identical(side$type, "region") ||
     !identical(side$expression$type, "comparison") ||
     !side$expression$operator %in% c("<", "<=", ">", ">=")){
    stop("Linear targets support only point and simple region statements.",
         call. = FALSE)
  }

  left  <- .hypothesis_linear_target_form(side$expression$left, parameter)
  right <- .hypothesis_linear_target_form(side$expression$right, parameter)
  form  <- .hypothesis_linear_target_add(left, right, -1)
  list(
    operator     = side$expression$operator,
    coefficients = .hypothesis_linear_target_scaled_coefficients(form),
    numerator_coefficients = form$coefficients,
    value        = .hypothesis_linear_target_divide(-form$constant, form$divisor),
    divisor      = form$divisor
  )
}


# An expression linear in the levels of 'parameter' as a constant 'offset'
# and level 'coefficients' (named by level, sorted, without zeros).
.hypothesis_linear_target_form <- function(node, parameter){

  expression <- .bt_hypothesis_node_language(node)
  symbols <- .hypothesis_expression_symbols(.hypothesis_parse_expression(expression))
  if(length(symbols) == 0L) return(.hypothesis_affine_read(expression, symbols))
  references <- hypothesis_parse_level_reference(vapply(symbols, function(symbol){
    .hypothesis_expression_text(as.name(symbol))
  }, character(1)))
  if(all(references$direct) && any(references$parameter != parameter)){
    stop("A linear target may reference levels of only '", parameter, "'.", call. = FALSE)
  }
  if(any(!references$direct) || any(references$parameter != parameter)){
    stop("A linear target must be a linear combination of levels of '", parameter,
         "' and numbers.", call. = FALSE)
  }
  form <- .hypothesis_affine_read(expression, symbols)
  if(is.null(form)){
    stop("A linear target must be a linear combination of levels of '", parameter,
         "' and numbers.", call. = FALSE)
  }
  names(form$coefficients) <- references$level[match(names(form$coefficients), references$symbol)]
  form
}

.hypothesis_linear_target_add <- function(left, right, sign){

  .hypothesis_affine_add(left, right, sign)
}

.hypothesis_linear_target_scale <- function(form, factor){

  .hypothesis_affine_scale(form, .hypothesis_affine_new(factor))
}


.hypothesis_linear_target_combine_weights <- function(level_weights,
                                                      coefficients){

  columns <- sort(unique(unlist(lapply(level_weights, names), use.names = FALSE)))
  out <- stats::setNames(numeric(length(columns)), columns)
  for(level in names(level_weights)){
    columns_i <- names(level_weights[[level]])
    out[columns_i] <-
      .hypothesis_affine_sum(out[columns_i],
        .hypothesis_affine_product(coefficients[[level]], level_weights[[level]]))
  }
  out
}


.hypothesis_linear_target_prior_atom_free <- function(prior_density){

  points <- prior_density$points
  inherits(prior_density, "prior_density") &&
    is.data.frame(points) &&
    all(c("x", "p") %in% names(points)) &&
    nrow(points) == 0L
}


.hypothesis_linear_target_posterior_atoms_declared <- function(posterior){

  atoms <- .posterior_atoms_get(posterior)
  !is.null(atoms) && isTRUE(atoms$declared)
}


.hypothesis_linear_target_rewrite <- function(ast, forms, target_name){

  form_i <- 0L
  statements <- vapply(ast$statements, function(statement){
    side_text <- vapply(c("left", "right"), function(side_name){
      form_i <<- form_i + 1L
      form <- forms[[form_i]]
      paste(
        target_name,
        form$operator,
        .hypothesis_number_label(form$value)
      )
    }, character(1))
    if(isTRUE(statement$explicit)){
      paste(side_text[[1L]], "vs", side_text[[2L]])
    }else{
      side_text[[1L]]
    }
  }, character(1))
  out <- hypothesis_parse(statements)
  for(i in seq_along(out$statements)){
    out$statements[[i]]$left$label  <- ast$statements[[i]]$left$label
    out$statements[[i]]$right$label <- ast$statements[[i]]$right$label
  }
  .bt_validate_hypothesis_ast(out)
  out
}

.hypothesis_linear_target_scaled_coefficients <- function(form){

  coefficients <- form$coefficients / form$divisor
  if(any(!is.finite(coefficients)) || any(form$coefficients != 0 & coefficients == 0)){
    .hypothesis_numerical_stop("coefficient_range", "public linear target units")
  }
  coefficients
}

.hypothesis_linear_target_divide <- function(value, divisor){

  out <- value / divisor
  if(!is.finite(out) || (value != 0 && out == 0) || out * divisor != value){
    .hypothesis_numerical_stop("boundary_range", "compiled target threshold")
  }
  out
}
