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

.bt_parse_random_effects <- function(formula, double_bar = "diag"){

  .bt_require_reformulas()

  if(!inherits(formula, "formula")){
    formula <- stats::as.formula(formula)
  }

  expand_method <- switch(
    double_bar,
    diag = "diag_special",
    split = "split",
    stop("'double_bar' must be either 'diag' or 'split'.", call. = FALSE)
  )

  bars <- reformulas::findbars_x(
    formula,
    specials = .bt_random_effect_specials(),
    default.special = NULL,
    expand_doublevert_method = expand_method
  )
  bars <- .bt_merge_random_wrapper_bars(formula, bars, expand_method)

  terms <- lapply(seq_along(bars), function(i){
    .bt_random_effect_term_from_call(bars[[i]], index = i, env = environment(formula))
  })

  out <- list(
    terms = terms,
    fixed_formula = .bt_fixed_formula(formula),
    full_formula = formula,
    policy = list(double_bar = double_bar)
  )
  class(out) <- c("BayesTools_random_effects", "list")

  return(out)
}

.bt_merge_random_wrapper_bars <- function(formula, bars, expand_method){

  wrapper_bars <- .bt_find_random_wrapper_calls(formula)
  if(length(wrapper_bars) == 0L){
    return(bars)
  }

  .bt_find_random_effect_calls_ordered(formula, expand_method)
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

  call_name <- as.character(x[[1L]])
  if(call_name %in% c("random", "re")){
    return(list(x))
  }
  if(call_name %in% setdiff(.bt_random_effect_specials(), c("random", "re")) &&
     .bt_random_effect_has_bar_arg(x)){
    return(list(x))
  }
  if(identical(x[[1L]], as.name("|")) || identical(x[[1L]], as.name("||"))){
    term_formula <- stats::as.formula(call("~", x), env = env)
    return(reformulas::findbars_x(
      term_formula,
      specials = .bt_random_effect_specials(),
      default.special = NULL,
      expand_doublevert_method = expand_method
    ))
  }

  out <- list()
  for(i in seq.int(2L, length(x))){
    out <- c(out, .bt_find_random_effect_calls_ordered(x[[i]], expand_method, env = env))
  }

  out
}

.bt_find_random_wrapper_calls <- function(x){

  if(inherits(x, "formula")){
    rhs_index <- if(length(x) == 3L) 3L else 2L
    return(.bt_find_random_wrapper_calls(x[[rhs_index]]))
  }

  if(!is.call(x)){
    return(list())
  }

  call_name <- as.character(x[[1L]])
  if(call_name %in% c("random", "re")){
    return(list(x))
  }
  if(call_name %in% setdiff(.bt_random_effect_specials(), c("random", "re")) &&
     .bt_random_effect_has_bar_arg(x)){
    return(list(x))
  }

  out <- list()
  for(i in seq.int(2L, length(x))){
    out <- c(out, .bt_find_random_wrapper_calls(x[[i]]))
  }

  out
}

.bt_random_effect_has_bar_arg <- function(x){

  args <- as.list(x)
  if(length(args) < 2L){
    return(FALSE)
  }
  any(vapply(args[-1L], function(arg){
    is.call(arg) &&
      (identical(arg[[1L]], as.name("|")) ||
         identical(arg[[1L]], as.name("||")))
  }, logical(1)))
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
    wrapper_args <- as.list(x)
    term_arg <- .bt_random_effect_call_term_arg(
      wrapper_args,
      "Random-effect wrappers"
    )
    if("name" %in% names(wrapper_args)){
      block_name <- .bt_random_effect_eval_name(wrapper_args[["name"]], env)
      has_explicit_name <- TRUE
    }
    if("covariance" %in% names(wrapper_args)){
      structure <- .bt_random_covariance_from_arg(wrapper_args[["covariance"]], env)
      explicit_special <- TRUE
    }
    extra_args <- .bt_random_effect_extra_named_args(
      wrapper_args,
      c("", "name", "covariance")
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
  group_expr <- x[[3L]]
  group_label <- .bt_deparse_expr(group_expr)
  independent <- structure %in% c("diag", "id")
  if(is.null(block_name)){
    block_name <- group_label
  }
  block_name <- .bt_random_effect_sanitize_name(block_name)

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

.bt_random_effect_call_term_arg <- function(args, label){

  arg_names <- names(args)
  if(is.null(arg_names)){
    arg_names <- rep("", length(args))
  }
  unnamed_positions <- which(!nzchar(arg_names) & seq_along(args) > 1L)
  if(length(unnamed_positions) != 1L){
    stop(
      label,
      " must contain exactly one unnamed term of the form 'expr | group'.",
      call. = FALSE
    )
  }

  args[[unnamed_positions]]
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
    "__default_factor", "__default_continuous"
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

.bt_fixed_formula <- function(formula){

  .bt_require_reformulas()

  fixed_formula <- reformulas::nobars(formula)
  environment(fixed_formula) <- environment(formula)

  return(fixed_formula)
}

.bt_validate_random_effect_block_names <- function(terms, prior_random = NULL){

  if(length(terms) == 0L){
    return(invisible(TRUE))
  }

  block_names <- vapply(terms, function(term) term$block_name, character(1))
  .bt_validate_random_effect_reserved_name(block_names)
  if(anyDuplicated(block_names)){
    stop("Random-effect block names must be unique.", call. = FALSE)
  }

  if(!is.null(prior_random)){
    .bt_check_prior_random(prior_random)
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

.bt_random_group_values <- function(term, data){

  if(isTRUE(term$group_is_symbol)){
    if(!term$group_label %in% colnames(data)){
      stop("The '", term$group_label, "' grouping factor is missing in the data set.", call. = FALSE)
    }
    return(.bt_validate_random_group_values(data[[term$group_label]], term, data))
  }

  term_env <- environment(term$term_formula)
  if(is.null(term_env)){
    term_env <- parent.frame()
  }
  if(.bt_random_group_is_colon_expr(term$group_expr)){
    value <- try(
      .bt_random_group_interaction_values(term$group_expr, data, term_env),
      silent = TRUE
    )
  }else{
    value <- try(eval(term$group_expr, envir = data, enclos = term_env), silent = TRUE)
  }
  if(inherits(value, "try-error")){
    stop(
      "Could not evaluate random-effect grouping expression '",
      term$group_label,
      "' in the data set.",
      call. = FALSE
    )
  }

  .bt_validate_random_group_values(value, term, data)
}

.bt_random_group_is_colon_expr <- function(expr){

  is.call(expr) && identical(expr[[1L]], as.name(":")) && length(expr) == 3L
}

.bt_random_group_colon_terms <- function(expr){

  if(.bt_random_group_is_colon_expr(expr)){
    return(c(
      .bt_random_group_colon_terms(expr[[2L]]),
      .bt_random_group_colon_terms(expr[[3L]])
    ))
  }

  list(expr)
}

.bt_random_group_interaction_values <- function(expr, data, env){

  components <- .bt_random_group_colon_terms(expr)
  values <- lapply(components, function(component){
    eval(component, envir = data, enclos = env)
  })

  do.call(interaction, c(values, list(drop = TRUE, sep = ":")))
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
  if(anyNA(value)){
    stop(
      "Random-effect grouping expression '",
      term$group_label,
      "' must not contain missing values.",
      call. = FALSE
    )
  }

  value
}
