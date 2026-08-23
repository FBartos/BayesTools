# Internal random-effect formula adapter.

.bt_random_effect_specials <- function(){
  c("random", "re", "us", "un", "diag", "id", "cs", "hcs", "ar1", "ar", "car", "har")
}

.bt_require_reformulas <- function(){
  if(!requireNamespace("reformulas", quietly = TRUE)){
    stop("The 'reformulas' package is required for random-effect formulas.", call. = FALSE)
  }
}

.bt_deparse_expr <- function(x){
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

.bt_rhs_formula <- function(rhs, env){
  stats::as.formula(call("~", rhs), env = env)
}

.bt_random_effect_normalize_special_term_parens <- function(x){

  if(inherits(x, "formula")){
    rhs_index <- if(length(x) == 3L) 3L else 2L
    x[[rhs_index]] <- .bt_random_effect_normalize_special_term_parens(x[[rhs_index]])
    return(x)
  }
  if(!is.call(x)){
    return(x)
  }

  if(length(x) >= 2L){
    for(i in seq.int(2L, length(x))){
      x[[i]] <- .bt_random_effect_normalize_special_term_parens(x[[i]])
    }
  }
  if(!is.symbol(x[[1L]]) ||
     !as.character(x[[1L]]) %in% setdiff(.bt_random_effect_specials(), c("random", "re"))){
    return(x)
  }

  args <- as.list(x)
  term_position <- .bt_random_effect_unnamed_term_positions(args)
  if(length(term_position) != 1L){
    return(x)
  }
  term_arg <- .bt_random_effect_strip_group_parens(args[[term_position]])
  if(!.bt_random_effect_is_bar_call(term_arg)){
    return(x)
  }

  args[[term_position]] <- term_arg
  as.call(args)
}

.bt_parse_random_effects <- function(formula, double_bar = "diag"){

  .bt_require_reformulas()

  if(!inherits(formula, "formula")){
    formula <- stats::as.formula(formula)
  }
  parsed_formula <- .bt_random_effect_normalize_special_term_parens(formula)

  expand_method <- switch(
    double_bar,
    diag = "diag_special",
    split = "split",
    stop("'double_bar' must be either 'diag' or 'split'.", call. = FALSE)
  )

  if(length(.bt_find_random_wrapper_calls(parsed_formula)) > 0L){
    bars <- .bt_find_random_effect_calls_ordered(parsed_formula, expand_method)
  }else{
    bars <- reformulas::findbars_x(
      parsed_formula,
      specials = .bt_random_effect_specials(),
      default.special = NULL,
      expand_doublevert_method = expand_method
    )
  }

  terms <- lapply(seq_along(bars), function(i){
    .bt_random_effect_term_from_call(bars[[i]], index = i, env = environment(formula))
  })

  out <- list(
    terms = terms,
    fixed_formula = .bt_fixed_formula(parsed_formula),
    full_formula = formula,
    policy = list(double_bar = double_bar)
  )
  class(out) <- c("BayesTools_random_effects", "list")

  return(out)
}

# TRUE for a `expr | group` or `expr || group` bar call.
.bt_random_effect_is_bar_call <- function(x){
  is.call(x) && length(x) == 3L &&
    (identical(x[[1L]], as.name("|")) || identical(x[[1L]], as.name("||")))
}

# TRUE for a random()/re() wrapper or a covariance special (diag/cs/us/...) that
# encloses a bar. These are the calls whose grouping the ordered walker has to
# expand itself (reformulas::findbars_x cannot, see .bt_random_effect_expand_call_nesting).
.bt_random_effect_is_wrapper_call <- function(x){
  if(!is.call(x)){
    return(FALSE)
  }
  call_name <- as.character(x[[1L]])
  any(call_name %in% c("random", "re")) ||
    (any(call_name %in% setdiff(.bt_random_effect_specials(), c("random", "re"))) &&
       .bt_random_effect_has_bar_arg(x))
}

.bt_find_random_effect_calls_ordered <- function(x, expand_method,
                                                env = parent.frame()){

  if(inherits(x, "formula")){
    rhs_index <- if(length(x) == 3L) 3L else 2L
    return(.bt_find_random_effect_calls_ordered(
      x[[rhs_index]],
      expand_method,
      env = environment(x)
    ))
  }

  if(!is.call(x)){
    return(list())
  }

  if(.bt_random_effect_is_wrapper_call(x)){
    return(.bt_random_effect_expand_call_nesting(x, expand_method, env))
  }
  if(.bt_random_effect_is_bar_call(x)){
    term_formula <- stats::as.formula(call("~", x), env = env)
    return(reformulas::findbars_x(
      term_formula,
      specials = .bt_random_effect_specials(),
      default.special = NULL,
      expand_doublevert_method = expand_method
    ))
  }

  out <- list()
  if(length(x) >= 2L){
    for(i in seq.int(2L, length(x))){
      out <- c(out, .bt_find_random_effect_calls_ordered(x[[i]], expand_method, env = env))
    }
  }

  out
}

# Expand a covariance-special or random()/re() wrapper call whose grouping is a
# nested expression ('g1/g2') into one call per nesting level, preserving the
# wrapper structure and all named arguments. Returns 'list(x)' unchanged when the
# grouping is not nested, so non-nested wrappers behave exactly as before.
#
# reformulas::findbars_x already expands nested grouping for *un-named* specials
# and wrappers, but it drops or fails to expand wrappers carrying named arguments
# such as name =/covariance =/hom =. The ordered walker therefore has to perform
# the nesting expansion itself while keeping those named arguments intact.
.bt_random_effect_expand_call_nesting <- function(x, expand_method, env){

  args <- as.list(x)
  term_position <- .bt_random_effect_unnamed_term_positions(args)
  # A malformed wrapper (not exactly one unnamed term) is left for
  # .bt_random_effect_term_from_call to reject with a context-rich message.
  if(length(term_position) != 1L){
    return(list(x))
  }

  inner <- args[[term_position]]
  if(!.bt_random_effect_is_bar_call(inner)){
    return(list(x))
  }

  group_levels <- .bt_random_effect_nested_group_levels(inner[[3L]], expand_method, env)
  if(is.null(group_levels) || length(group_levels) <= 1L){
    return(list(x))
  }

  # An explicit block name becomes a per-level prefix so the user's naming intent
  # is preserved while keeping the expanded blocks distinct (e.g.
  # random(1 | site/plot, name = "spatial") -> "spatial_plot_site", "spatial_site").
  base_name <- if("name" %in% names(args)){
    .bt_random_effect_eval_name(args[["name"]], env)
  }else{
    NULL
  }

  lapply(group_levels, function(level_group){
    new_args <- args
    new_args[[term_position]] <- as.call(list(inner[[1L]], inner[[2L]], level_group))
    if(!is.null(base_name)){
      new_args[["name"]] <- paste0(base_name, "_", .bt_deparse_expr(level_group))
    }
    as.call(new_args)
  })
}

# Return the per-level grouping expressions for a nested grouping expression
# ('g1/g2/...'), in the exact order and interaction convention that
# reformulas::findbars_x uses for plain bars (so special/wrapper and plain-bar
# nesting expand identically). Returns NULL when the expression is not nested.
.bt_random_effect_nested_group_levels <- function(group_expr, expand_method, env){

  group_expr <- .bt_random_effect_strip_group_parens(group_expr)
  if(!(is.call(group_expr) && identical(group_expr[[1L]], as.name("/")))){
    return(NULL)
  }

  dummy_formula <- stats::as.formula(call("~", call("|", 1, group_expr)), env = env)
  expanded <- reformulas::findbars_x(
    dummy_formula,
    specials = .bt_random_effect_specials(),
    default.special = NULL,
    expand_doublevert_method = expand_method
  )

  lapply(expanded, function(bar) bar[[3L]])
}

.bt_random_effect_strip_group_parens <- function(expr){

  while(is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L){
    expr <- expr[[2L]]
  }

  expr
}

.bt_random_effect_validate_group_expr <- function(expr){

  if(.bt_random_effect_group_expr_allowed(expr)){
    return(invisible(TRUE))
  }

  stop(
    "Unsupported random-effect grouping expression '",
    .bt_deparse_expr(expr),
    "'. Random-effect grouping expressions must be variables, ':' interactions, or '/' ",
    "nested grouping. Create an explicit data column for other groupings.",
    call. = FALSE
  )
}

.bt_random_effect_group_expr_allowed <- function(expr){

  expr <- .bt_random_effect_strip_group_parens(expr)
  if(is.symbol(expr)){
    return(TRUE)
  }
  operator <- if(is.call(expr) && is.symbol(expr[[1L]])) as.character(expr[[1L]]) else ""
  if(is.call(expr) && length(expr) == 3L && operator %in% c(":", "/")){
    return(
      .bt_random_effect_group_expr_allowed(expr[[2L]]) &&
        .bt_random_effect_group_expr_allowed(expr[[3L]])
    )
  }

  FALSE
}

.bt_random_effect_validate_predictor_expr <- function(expr, block_name){

  validate_expression <- function(expression){
    if(is.symbol(expression)){
      return(invisible(TRUE))
    }
    if(is.numeric(expression) && length(expression) == 1L &&
       is.finite(expression) && expression %in% c(0, 1)){
      return(invisible(TRUE))
    }
    if(is.call(expression)){
      call_name <- if(is.symbol(expression[[1L]])){
        as.character(expression[[1L]])
      }else{
        ""
      }
      if(call_name %in% c("+", "-", "*", ":", "/", "^", "(")){
        for(argument in as.list(expression)[-1L]){
          validate_expression(argument)
        }
        return(invisible(TRUE))
      }
    }

    stop(
      "Unsupported random-effect predictor call '",
      .bt_deparse_expr(expression),
      "' in block '", block_name,
      "'. Create the transformed predictor as an explicit data column and ",
      "reference that column by name.",
      call. = FALSE
    )
  }

  validate_expression(expr)
  invisible(TRUE)
}

.bt_find_random_wrapper_calls <- function(x){

  if(inherits(x, "formula")){
    rhs_index <- if(length(x) == 3L) 3L else 2L
    return(.bt_find_random_wrapper_calls(x[[rhs_index]]))
  }

  if(!is.call(x)){
    return(list())
  }

  if(.bt_random_effect_is_wrapper_call(x)){
    return(list(x))
  }

  out <- list()
  if(length(x) >= 2L){
    for(i in seq.int(2L, length(x))){
      out <- c(out, .bt_find_random_wrapper_calls(x[[i]]))
    }
  }

  out
}

.bt_random_effect_has_bar_arg <- function(x){

  args <- as.list(x)
  if(length(args) < 2L){
    return(FALSE)
  }
  any(vapply(args[-1L], .bt_random_effect_is_bar_call, logical(1)))
}

.bt_random_effect_term_from_call <- function(x, index, env){

  specials <- .bt_random_effect_specials()
  structure <- "us"
  explicit_special <- FALSE
  special_call <- x
  hom <- NULL
  block_name <- NULL
  has_explicit_name <- FALSE
  extra_args <- character()

  if(is.call(x) && as.character(x[[1L]]) %in% c("random", "re")){
    wrapper_name <- as.character(x[[1L]])
    wrapper_args <- as.list(x)
    term_arg <- .bt_random_effect_call_term_arg(
      wrapper_args,
      "Random-effect wrappers"
    )
    .bt_random_effect_reject_wrapped_special(wrapper_name, term_arg)
    if("name" %in% names(wrapper_args)){
      block_name <- .bt_random_effect_eval_name(wrapper_args[["name"]], env)
      has_explicit_name <- TRUE
    }
    if("covariance" %in% names(wrapper_args)){
      structure <- .bt_random_covariance_from_arg(wrapper_args[["covariance"]], env)
      explicit_special <- TRUE
    }
    if("hom" %in% names(wrapper_args)){
      hom <- .bt_random_effect_eval_hom(wrapper_args[["hom"]], env)
      resolved_hom <- .bt_random_effect_resolve_wrapper_hom(structure, hom)
      structure <- resolved_hom$structure
      hom <- resolved_hom$hom
    }
    extra_args <- .bt_random_effect_extra_named_args(
      wrapper_args,
      c("", "name", "covariance", "hom")
    )
    x <- term_arg
    special_call <- x
  }

  if(is.call(x) && as.character(x[[1L]]) %in% setdiff(specials, c("random", "re"))){
    structure <- .bt_random_covariance_from_special(as.character(x[[1L]]))
    explicit_special <- TRUE
    special_args <- as.list(x)
    term_arg <- .bt_random_effect_call_term_arg(
      special_args,
      paste0("The '", structure, "' random-effect covariance structure")
    )
    if("hom" %in% names(special_args)){
      hom <- .bt_random_effect_eval_hom(special_args[["hom"]], env)
    }
    if("name" %in% names(special_args)){
      block_name <- .bt_random_effect_eval_name(special_args[["name"]], env)
      has_explicit_name <- TRUE
    }
    extra_args <- c(
      extra_args,
      .bt_random_effect_extra_named_args(special_args, c("", "hom", "name"))
    )
    x <- term_arg
  }

  if(is.call(x) && identical(x[[1L]], as.name("||")) && length(x) == 3L){
    if(isTRUE(explicit_special) && !identical(structure, "diag")){
      stop(
        "Random-effect covariance structure '", structure,
        "' cannot be combined with '||' syntax. Use covariance = 'diag' or replace '||' with '|'.",
        call. = FALSE
      )
    }
    structure <- "diag"
    independent <- TRUE
    x[[1L]] <- as.name("|")
  }

  if(!is.call(x) || !identical(x[[1L]], as.name("|")) || length(x) != 3L){
    stop("Random-effect terms must have the form 'expr | group'.", call. = FALSE)
  }

  expr <- x[[2L]]
  group_expr <- .bt_random_effect_strip_group_parens(x[[3L]])
  x[[3L]] <- group_expr
  group_label <- .bt_deparse_expr(group_expr)
  .bt_random_effect_validate_group_expr(group_expr)
  independent <- structure %in% c("diag", "id")
  if(is.null(block_name)){
    block_name <- group_label
  }
  block_name <- .bt_random_effect_sanitize_name(block_name)
  .bt_random_effect_validate_predictor_expr(expr, block_name)

  term <- list(
    id = paste0("random_", index),
    block_name = block_name,
    has_explicit_name = has_explicit_name,
    original_call = special_call,
    bar_call = x,
    expr = expr,
    group_expr = group_expr,
    group_label = group_label,
    group_is_symbol = is.symbol(group_expr),
    structure = structure,
    explicit_special = explicit_special,
    independent = independent,
    hom = hom,
    extra_args = extra_args,
    term_formula = .bt_rhs_formula(expr, env = env),
    index = index
  )

  attr(term, "grouping_factor") <- group_label
  attr(term, "random_block") <- block_name
  attr(term, "independent") <- independent
  attr(term, "structure") <- structure
  class(term) <- c("BayesTools_random_effect_term", "list")

  return(term)
}

.bt_random_effect_reject_wrapped_special <- function(wrapper_name, term_arg){

  if(!.bt_random_effect_is_covariance_special_call(term_arg)){
    return(invisible(TRUE))
  }

  special <- as.character(term_arg[[1L]])
  structure <- .bt_random_covariance_from_special(special)
  stop(
    "Do not wrap covariance-special random-effect calls in '", wrapper_name,
    "()'. Use direct '", special, "(expr | group)' syntax or '",
    wrapper_name, "(expr | group, covariance = \"", structure, "\")'.",
    call. = FALSE
  )
}

.bt_random_effect_is_covariance_special_call <- function(x){

  is.call(x) &&
    as.character(x[[1L]]) %in% setdiff(.bt_random_effect_specials(), c("random", "re"))
}

.bt_random_effect_resolve_wrapper_hom <- function(structure, hom){

  if(!identical(hom, FALSE)){
    return(list(structure = structure, hom = hom))
  }
  if(identical(structure, "cs")){
    return(list(structure = "hcs", hom = NULL))
  }
  if(identical(structure, "ar1")){
    return(list(structure = "har", hom = NULL))
  }

  list(structure = structure, hom = hom)
}

.bt_random_effect_metadata_block_detail <- function(random_term){

  block <- random_term$block_name
  if(!is.null(block) && length(block) == 1L && nzchar(block)){
    return(paste0(" for block '", block, "'"))
  }

  ""
}

.bt_random_effect_structure <- function(random_term,
                                        context = "Random-effect metadata"){

  structure <- random_term$structure
  if(!is.null(structure) && length(structure) == 1L){
    structure <- as.character(structure)
    if(!is.na(structure) && nzchar(structure)){
      return(tolower(.bt_random_covariance_normalize(structure)))
    }
  }

  stop(
    context,
    .bt_random_effect_metadata_block_detail(random_term),
    " is missing canonical 'random_term$structure'.",
    call. = FALSE
  )
}

.bt_random_effect_homogeneous_sd_metadata <- function(
    random_term,
    context = "Random-effect metadata"){

  homogeneous_sd <- random_term$homogeneous_sd
  if(is.logical(homogeneous_sd) && length(homogeneous_sd) == 1L &&
     !is.na(homogeneous_sd)){
    return(homogeneous_sd)
  }

  stop(
    context,
    .bt_random_effect_metadata_block_detail(random_term),
    " is missing canonical 'random_term$homogeneous_sd'.",
    call. = FALSE
  )
}

.bt_random_effect_correlation_metadata <- function(
    random_term,
    structure = NULL,
    context = "Random-effect metadata"){

  if(is.null(structure)){
    structure <- .bt_random_effect_structure(random_term, context = context)
  }
  n_columns <- random_term$n_columns
  requires_correlation <- is.numeric(n_columns) &&
    length(n_columns) == 1L &&
    !is.na(n_columns) &&
    n_columns > 1L &&
    structure %in% c("us", "cs", "hcs", "ar1", "car", "har")

  correlation <- random_term$correlation
  if(!isTRUE(requires_correlation)){
    return(correlation)
  }
  if(is.list(correlation) && !is.null(correlation$type) &&
     length(correlation$type) == 1L && !is.na(correlation$type) &&
     nzchar(correlation$type)){
    return(correlation)
  }

  stop(
    context,
    .bt_random_effect_metadata_block_detail(random_term),
    " is missing canonical 'random_term$correlation'.",
    call. = FALSE
  )
}

# Positions of the unnamed arguments of a wrapper/special call (the term slots),
# excluding the call head. A well-formed random-effect call has exactly one.
.bt_random_effect_unnamed_term_positions <- function(args){

  arg_names <- names(args)
  if(is.null(arg_names)){
    arg_names <- rep("", length(args))
  }
  which(!nzchar(arg_names) & seq_along(args) > 1L)
}

.bt_random_effect_call_term_arg <- function(args, label){

  term_position <- .bt_random_effect_unnamed_term_positions(args)
  if(length(term_position) != 1L){
    stop(
      label,
      " must contain exactly one unnamed term of the form 'expr | group'.",
      call. = FALSE
    )
  }

  args[[term_position]]
}

.bt_random_effect_extra_named_args <- function(args, allowed){

  arg_names <- names(args)
  if(is.null(arg_names)){
    return(character())
  }
  setdiff(arg_names[nzchar(arg_names)], allowed)
}

.bt_random_covariance_from_special <- function(x){

  switch(
    x,
    us   = "us",
    un   = "us",
    diag = "diag",
    id   = "id",
    cs   = "cs",
    hcs  = "hcs",
    ar1  = "ar1",
    ar   = "ar1",
    car  = "car",
    har  = "har",
    stop("Unknown random-effect covariance special '", x, "'.", call. = FALSE)
  )
}

.bt_random_covariance_from_arg <- function(x, env){

  value <- try(eval(x, envir = env), silent = TRUE)
  if(inherits(value, "try-error")){
    value <- .bt_deparse_expr(x)
  }
  check_char(value, "covariance", allow_NA = FALSE)
  tolower(.bt_random_covariance_normalize(value))
}

.bt_random_effect_eval_name <- function(x, env){

  value <- try(eval(x, envir = env), silent = TRUE)
  if(inherits(value, "try-error")){
    value <- .bt_deparse_expr(x)
  }
  check_char(value, "name", allow_NA = FALSE)
  if(!nzchar(value)){
    stop("Random-effect block names cannot be empty.", call. = FALSE)
  }
  value
}

.bt_random_effect_eval_hom <- function(x, env){

  value <- try(eval(x, envir = env), silent = TRUE)
  if(inherits(value, "try-error")){
    value <- x
  }
  if(!is.logical(value) || length(value) != 1L || is.na(value)){
    stop("'hom' must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.bt_random_effect_sanitize_name <- function(x){

  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if(!nzchar(x)){
    stop("Random-effect block names must contain at least one letter, digit, or underscore.", call. = FALSE)
  }
  if(!grepl("^[A-Za-z]", x)){
    x <- paste0("RE_", x)
  }
  x
}

.bt_random_effect_reserved_terms <- function(){

  c(
    "__xXx__", "__xREx__", "xRE_ALLOCx", "xRE_PRECx", "xRE_CORx",
    "xRE_Zx", "xRE_STDx", "xRE_COEFx", "xRE_MAPx", "xRE_DATAx",
    "__xRE_SUMMARY__", "__default_factor", "__default_continuous"
  )
}

.bt_validate_random_effect_reserved_name <- function(x,
                                                    context = "naming random-effect blocks"){

  for(reserved_term in .bt_random_effect_reserved_terms()){
    if(any(grepl(reserved_term, x, fixed = TRUE))){
      stop(
        "'", reserved_term,
        "' string is internally used by the BayesTools package and can't be used for ",
        context,
        ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_validate_categorical_level_names <- function(level_names, variable_name,
                                                 context = "Categorical variable"){

  level_names <- unique(as.character(level_names))
  level_names <- level_names[!is.na(level_names)]
  for(reserved_term in .bt_random_effect_reserved_terms()){
    offending <- level_names[
      grepl(reserved_term, level_names, fixed = TRUE)
    ]
    if(length(offending) > 0L){
      stop(
        context, " '", variable_name,
        "' contains the internally reserved token '", reserved_term,
        "' in level '", offending[[1L]],
        "'. Rename the level before fitting or prediction.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_validate_categorical_values <- function(value, variable_name,
                                            context = "Categorical variable"){

  level_names <- if(is.factor(value)){
    levels(value)
  }else{
    unique(as.character(value))
  }
  .bt_validate_categorical_level_names(
    level_names,
    variable_name,
    context = context
  )
}

.bt_fixed_formula <- function(formula){

  .bt_require_reformulas()

  fixed_formula <- reformulas::nobars(formula)
  environment(fixed_formula) <- environment(formula)

  return(fixed_formula)
}

.bt_validate_random_effect_block_names <- function(terms, prior_random = NULL){

  if(!is.null(prior_random)){
    .bt_check_prior_random(prior_random)
  }

  block_names <- character()
  if(length(terms) > 0L){
    block_names <- vapply(terms, function(term) term$block_name, character(1))
    .bt_validate_random_effect_reserved_name(block_names)
    if(anyDuplicated(block_names)){
      stop("Random-effect block names must be unique.", call. = FALSE)
    }
  }

  if(!is.null(prior_random)){
    unknown_blocks <- setdiff(.bt_random_prior_block_names(prior_random), block_names)
    if(length(unknown_blocks) > 0L){
      stop(
        "The following 'prior_random' block override names were not found in the formula: ",
        paste(unknown_blocks, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_as_random_effect_term <- function(x){

  if(inherits(x, "BayesTools_random_effect_term")){
    return(x)
  }

  if(is.character(x) && length(x) == 1L){
    formula <- stats::as.formula(paste0("~ (", x, ")"))
    terms <- .bt_parse_random_effects(formula)$terms
    if(length(terms) != 1L){
      stop("Expected exactly one random-effect term.", call. = FALSE)
    }
    return(terms[[1L]])
  }

  stop("Unrecognized random-effect term representation.", call. = FALSE)
}

.bt_validate_random_effect_term_supported <- function(term){

  structure <- .bt_random_effect_structure(term)

  if(!is.null(term$hom) && !identical(term$hom, FALSE) && !identical(term$hom, TRUE)){
    stop("'hom' must be TRUE or FALSE.", call. = FALSE)
  }

  if(!is.null(term$hom) && structure %in% c("hcs", "har")){
    stop("The 'hom' argument is not supported for already heteroscedastic random-effect structures.", call. = FALSE)
  }
  if(!is.null(term$hom) && identical(term$hom, TRUE) &&
     !structure %in% c("diag", "id", "cs", "ar1", "car")){
    stop("Homogeneous random-effect standard deviations are not supported yet.", call. = FALSE)
  }
  if(!is.null(term$hom) && identical(term$hom, FALSE) &&
     structure == "car"){
    stop(
      "Structure 'car' has homogeneous random-effect standard deviations; heteroscedastic CAR is not supported yet.",
      call. = FALSE
    )
  }
  if(!is.null(term$hom) && identical(term$hom, FALSE) &&
     structure %in% c("id", "cs", "ar1")){
    replacement <- switch(
      structure,
      id = "diag",
      cs = "hcs",
      ar1 = "har"
    )
    stop(
      "Structure '", structure, "' has homogeneous random-effect standard deviations. ",
      "Use '", replacement, "' for heteroscedastic standard deviations.",
      call. = FALSE
    )
  }

  if(length(term$extra_args) > 0L){
    stop(
      "The '", structure,
      "' random-effect covariance structure does not support extra arguments yet.",
      call. = FALSE
    )
  }

  if(structure %in% c("us", "diag", "id", "cs", "hcs", "ar1", "car", "har")){
    return(invisible(TRUE))
  }

  if(!isTRUE(term$independent)){
    stop(
      "The '", structure,
      "' random-effect covariance structure is not supported yet.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_group_is_colon_expr <- function(expr){

  expr <- .bt_random_effect_strip_group_parens(expr)
  is.call(expr) && identical(expr[[1L]], as.name(":")) && length(expr) == 3L
}

.bt_random_group_colon_terms <- function(expr){

  expr <- .bt_random_effect_strip_group_parens(expr)
  if(.bt_random_group_is_colon_expr(expr)){
    return(c(
      .bt_random_group_colon_terms(expr[[2L]]),
      .bt_random_group_colon_terms(expr[[3L]])
    ))
  }

  list(expr)
}

.bt_random_group_component_names <- function(term){

  components <- .bt_random_group_colon_terms(term$group_expr)
  component_names <- vapply(components, function(component){
    if(!is.symbol(component)){
      stop(
        "Random-effect grouping metadata for block '", term$block_name,
        "' contains a non-variable interaction component.",
        call. = FALSE
      )
    }
    as.character(component)
  }, character(1))
  if(anyDuplicated(component_names)){
    stop(
      "Random-effect grouping expression '", term$group_label,
      "' repeats a grouping variable.",
      call. = FALSE
    )
  }

  component_names
}

.bt_random_group_tuple_key <- function(values){

  values <- enc2utf8(as.character(values))
  encoded <- paste0(
    nchar(values, type = "bytes"),
    ":",
    values
  )
  paste(encoded, collapse = "|")
}

.bt_random_group_observations <- function(term, data){

  component_names <- .bt_random_group_component_names(term)
  missing_components <- component_names[!component_names %in% colnames(data)]
  if(length(missing_components) > 0L){
    stop(
      "The ",
      paste0("'", missing_components, "'", collapse = ", "),
      " random-effect grouping variable",
      if(length(missing_components) > 1L) "s are" else " is",
      " missing in the data set.",
      call. = FALSE
    )
  }

  component_values <- lapply(component_names, function(component_name){
    value <- .bt_validate_random_group_values(
      data[[component_name]],
      term,
      data
    )
    .bt_validate_categorical_values(
      value,
      component_name,
      context = "Random-effect grouping variable"
    )
    value
  })
  names(component_values) <- component_names
  tuple_values <- do.call(cbind, lapply(component_values, as.character))
  if(length(component_names) == 1L){
    tuple_values <- matrix(
      tuple_values,
      ncol = 1L,
      dimnames = list(NULL, component_names)
    )
  }else{
    colnames(tuple_values) <- component_names
  }
  tuple_keys <- apply(tuple_values, 1L, .bt_random_group_tuple_key)
  display_labels <- apply(tuple_values, 1L, paste, collapse = ":")

  list(
    component_names = component_names,
    component_values = component_values,
    component_levels = lapply(component_values, function(value){
      levels(as.factor(value))
    }),
    tuple_values = tuple_values,
    tuple_keys = unname(tuple_keys),
    display_labels = unname(display_labels)
  )
}

.bt_random_group_metadata <- function(term, data){

  observations <- .bt_random_group_observations(term, data)
  if(length(observations$component_names) == 1L){
    group_tuples <- matrix(
      observations$component_levels[[1L]],
      ncol = 1L,
      dimnames = list(NULL, observations$component_names)
    )
    group_tuple_keys <- apply(group_tuples, 1L, .bt_random_group_tuple_key)
  }else{
    first_rows <- which(!duplicated(observations$tuple_keys))
    component_codes <- vapply(
      seq_along(observations$component_names),
      function(i){
        match(
          observations$tuple_values[, i],
          observations$component_levels[[i]]
        )
      },
      integer(nrow(observations$tuple_values))
    )
    tuple_order <- do.call(
      order,
      as.data.frame(component_codes[first_rows, , drop = FALSE])
    )
    first_rows <- first_rows[tuple_order]
    group_tuples <- observations$tuple_values[first_rows, , drop = FALSE]
    group_tuple_keys <- observations$tuple_keys[first_rows]
  }
  group_tuple_index <- stats::setNames(
    seq_along(group_tuple_keys),
    group_tuple_keys
  )
  group_map <- unname(group_tuple_index[observations$tuple_keys])
  group_labels <- apply(group_tuples, 1L, paste, collapse = ":")
  group_levels <- .bt_random_group_unique_labels(
    group_labels,
    group_tuple_keys
  )

  list(
    values = group_levels[group_map],
    levels = group_levels,
    map = group_map,
    components = observations$component_names,
    component_levels = observations$component_levels,
    tuples = group_tuples,
    labels = unname(group_labels),
    tuple_keys = unname(group_tuple_keys),
    tuple_index = group_tuple_index
  )
}

.bt_random_group_unique_labels <- function(labels, tuple_keys,
                                           existing = character()){

  combined <- c(existing, labels)
  duplicate_labels <- duplicated(combined) |
    duplicated(combined, fromLast = TRUE)
  duplicate_labels <- tail(duplicate_labels, length(labels))
  labels[duplicate_labels] <- paste0(
    labels[duplicate_labels],
    " [",
    tuple_keys[duplicate_labels],
    "]"
  )
  labels
}

.bt_validate_random_group_values <- function(value, term, data){

  if(is.data.frame(value) || is.matrix(value) || is.list(value)){
    stop(
      "Random-effect grouping expression '",
      term$group_label,
      "' must evaluate to one atomic value per row of data.",
      call. = FALSE
    )
  }
  if(length(value) != nrow(data)){
    stop(
      "Random-effect grouping expression '",
      term$group_label,
      "' must evaluate to one value per row of data.",
      call. = FALSE
    )
  }
  if(anyNA(value) || (is.factor(value) && anyNA(levels(value)))){
    stop(
      "Random-effect grouping expression '",
      term$group_label,
      "' must not contain missing values.",
      call. = FALSE
    )
  }

  value
}
