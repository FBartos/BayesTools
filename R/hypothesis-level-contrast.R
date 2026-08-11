#' Compile a two-level hypothesis contrast
#'
#' @description
#' Compiles point and simple region statements that all refer to one ordered
#' difference between two levels of the same parameter. The returned scalar
#' marginal posterior retains the exact joint-prior density and structural
#' conditioning metadata needed by downstream inference methods.
#'
#' This is a deliberately narrow cross-package interface. Both levels must
#' have one fixed named linear-weight row, share a joint prior context and
#' conditioning event, and carry complete posterior-atom declarations. The
#' exact induced contrast prior must be atom-free; by absolute continuity, its
#' posterior is then atom-free as well. Scaled, nonlinear, row-varying, and
#' multiple-target hypotheses are rejected.
#'
#' @param posterior a factor-like marginal posterior containing the referenced
#'   levels.
#' @param hypothesis hypothesis text or a `BayesTools_hypothesis_ast`.
#' @param parameter exact parameter root used by the level references.
#'
#' @return A list with `posterior`, `hypothesis`, `parameter`, and `weights`.
#'   `posterior` is a scalar `marginal_posterior`; `hypothesis` is an equivalent
#'   AST written against that scalar contrast.
#'
#' @export
hypothesis_level_contrast <- function(posterior, hypothesis, parameter){

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

  target_name <- ".BayesTools_level_contrast"
  sides <- unlist(lapply(ast$statements, function(statement){
    list(statement$left, statement$right)
  }), recursive = FALSE)
  forms <- lapply(sides, .hypothesis_level_contrast_side)
  coefficients <- forms[[1L]]$coefficients
  same_target <- vapply(forms, function(form){
    identical(form$coefficients, coefficients)
  }, logical(1))
  if(!all(same_target)){
    stop(
      "Hypothesis statements must all use the same ordered level contrast.",
      call. = FALSE
    )
  }
  if(length(coefficients) != 2L ||
     !identical(sort(unname(coefficients)), c(-1, 1))){
    stop("A level contrast must be one unscaled difference between two levels.",
         call. = FALSE)
  }

  expected_prefix <- paste0(parameter, "[")
  valid_symbols <- startsWith(names(coefficients), expected_prefix) &
    endsWith(names(coefficients), "]")
  if(!all(valid_symbols)){
    stop("A level contrast may reference levels of only '", parameter, "'.",
         call. = FALSE)
  }
  levels <- substring(
    names(coefficients),
    nchar(expected_prefix) + 1L,
    nchar(names(coefficients)) - 1L
  )
  names(coefficients) <- levels
  missing <- setdiff(levels, names(posterior))
  if(length(missing) > 0L){
    stop(
      "Hypothesis references unknown level '", paste(missing, collapse = "', '"),
      "' for parameter '", parameter, "'.",
      call. = FALSE
    )
  }

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
    weights <- .hypothesis_prepare_level_weights(
      attr(posterior[[level]], "linear_weights", exact = TRUE)
    )
    if(nrow(weights) != 1L){
      stop("Level contrasts require one fixed linear-weight row per level.",
           call. = FALSE)
    }
    weights[1L, ]
  })
  names(level_weights) <- levels

  context <- .hypothesis_child_prior_context(posterior, levels)
  if(is.null(context)){
    context <- attr(posterior, "prior_density_context", exact = TRUE)
  }
  if(is.null(context) || !.hypothesis_is_prior_density_context(context)){
    stop("A valid joint prior context is required for a level contrast.",
         call. = FALSE)
  }
  .hypothesis_validate_level_weights_context(level_weights, context)

  weights <- .hypothesis_level_contrast_combine_weights(
    level_weights,
    coefficients
  )
  if(length(weights) == 0L || all(weights == 0)){
    stop("The level contrast has zero combined linear weight.", call. = FALSE)
  }
  prior_density <- .prior_density_from_context(context, weights)
  if(!.hypothesis_level_contrast_prior_atom_free(prior_density)){
    stop("The level contrast prior is not structurally atom-free.",
         call. = FALSE)
  }
  declared_atoms <- vapply(
    posterior[levels],
    .hypothesis_level_contrast_posterior_atoms_declared,
    logical(1)
  )
  if(!all(declared_atoms)){
    stop("Complete structural posterior-atom declarations are required for ",
         "a level contrast.", call. = FALSE)
  }

  values <- rep(0, draw_lengths[[1L]])
  for(level in levels){
    values <- values + coefficients[[level]] * level_draws[[level]]
  }
  class(values) <- c(
    "numeric", "marginal_posterior.level_contrast", "marginal_posterior"
  )
  attr(values, "parameter")             <- target_name
  attr(values, "linear_weights")        <- weights
  attr(values, "prior_density")         <- prior_density
  attr(values, "prior_density_context") <- context
  attr(values, "posterior_atoms") <- .posterior_atoms_new(
    column_names = target_name,
    source       = "level_contrast"
  )
  for(attribute in c(
    "conditional", "conditional_rule", "condition_key", "condition_event",
    "resolved_condition_event"
  )){
    attr(values, attribute) <- attr(
      posterior[[levels[[1L]]]],
      attribute,
      exact = TRUE
    )
  }

  rewritten <- .hypothesis_level_contrast_rewrite(
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


.hypothesis_level_contrast_side <- function(side){

  if(side$type %in% c("point", "not_point")){
    form <- .hypothesis_level_contrast_form(side$expression)
    return(list(
      operator     = if(side$type == "point") "=" else "!=",
      coefficients = form$coefficients,
      value        = side$value - form$offset
    ))
  }
  if(!identical(side$type, "region") ||
     !identical(side$expression$type, "comparison") ||
     !side$expression$operator %in% c("<", "<=", ">", ">=")){
    stop("Level contrasts support only point and simple region statements.",
         call. = FALSE)
  }

  left  <- .hypothesis_level_contrast_form(side$expression$left)
  right <- .hypothesis_level_contrast_form(side$expression$right)
  form  <- .hypothesis_level_contrast_subtract(left, right)
  list(
    operator     = side$expression$operator,
    coefficients = form$coefficients,
    value        = -form$offset
  )
}


.hypothesis_level_contrast_form <- function(node){

  if(identical(node$type, "literal")){
    return(list(offset = node$value, coefficients = numeric()))
  }
  if(identical(node$type, "level_reference")){
    name <- paste0(node$parameter, "[", node$level, "]")
    return(list(
      offset      = 0,
      coefficients = stats::setNames(1, name)
    ))
  }
  if(identical(node$type, "parentheses")){
    return(.hypothesis_level_contrast_form(node$expression))
  }
  if(identical(node$type, "arithmetic") &&
     identical(node$operator, "-") &&
     length(node$arguments) == 2L){
    left  <- .hypothesis_level_contrast_form(node$arguments[[1L]])
    right <- .hypothesis_level_contrast_form(node$arguments[[2L]])
    return(.hypothesis_level_contrast_subtract(left, right))
  }

  stop("A level contrast must be one unscaled difference between two levels.",
       call. = FALSE)
}


.hypothesis_level_contrast_subtract <- function(left, right){

  names_all <- sort(unique(c(
    names(left$coefficients),
    names(right$coefficients)
  )))
  coefficients <- stats::setNames(numeric(length(names_all)), names_all)
  if(length(left$coefficients) > 0L){
    coefficients[names(left$coefficients)] <- left$coefficients
  }
  if(length(right$coefficients) > 0L){
    coefficients[names(right$coefficients)] <-
      coefficients[names(right$coefficients)] - right$coefficients
  }
  coefficients <- coefficients[coefficients != 0]
  list(
    offset       = left$offset - right$offset,
    coefficients = coefficients
  )
}


.hypothesis_level_contrast_combine_weights <- function(level_weights,
                                                       coefficients){

  columns <- sort(unique(unlist(lapply(level_weights, names), use.names = FALSE)))
  out <- stats::setNames(numeric(length(columns)), columns)
  for(level in names(level_weights)){
    columns_i <- names(level_weights[[level]])
    out[columns_i] <-
      out[columns_i] + coefficients[[level]] * level_weights[[level]]
  }
  out
}


.hypothesis_level_contrast_prior_atom_free <- function(prior_density){

  points <- prior_density$points
  inherits(prior_density, "prior_density") &&
    is.data.frame(points) &&
    all(c("x", "p") %in% names(points)) &&
    nrow(points) == 0L
}


.hypothesis_level_contrast_posterior_atoms_declared <- function(posterior){

  atoms <- .posterior_atoms_get(posterior)
  !is.null(atoms) && isTRUE(atoms$declared)
}


.hypothesis_level_contrast_rewrite <- function(ast, forms, target_name){

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
