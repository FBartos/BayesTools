# Replay a selected set of existing terms with its original factor coding.
.bt_formula_selected_terms <- function(formula, design){

  requested <- stats::terms(formula)
  fitted <- design$terms
  requested_factors <- attr(requested, "factors")
  fitted_factors <- attr(fitted, "factors")
  fitted_labels <- attr(fitted, "term.labels")
  requested_labels <- attr(requested, "term.labels")
  if(length(fitted_labels) == 0L && length(fitted_factors) == 0L){
    fitted_factors <- matrix(integer(), 0L, 0L,
      dimnames = list(character(), character()))
  }
  if(!inherits(fitted, "terms") || !is.matrix(fitted_factors) ||
     length(fitted_labels) != ncol(fitted_factors)){
    .bt_stop_refit_required("Fitted formula terms are malformed. Refit the model with this version of BayesTools.")
  }
  identity <- function(factors, column){
    sort(rownames(factors)[factors[, column] != 0L])
  }
  matched <- vapply(seq_along(requested_labels), function(column){
    candidates <- which(vapply(seq_along(fitted_labels), function(fitted_column){
      identical(identity(requested_factors, column), identity(fitted_factors, fitted_column))
    }, logical(1)))
    if(length(candidates) != 1L){
      stop("The selected formula term '", requested_labels[[column]],
           "' is unavailable in the fitted formula.", call. = FALSE)
    }
    candidates[[1L]]
  }, integer(1))
  matched <- sort(unique(matched))
  needed <- if(length(matched)){
    which(rowSums(fitted_factors[, matched, drop = FALSE] != 0L) > 0L)
  }else{
    integer()
  }
  selected_formula <- if(length(matched)){
    stats::reformulate(fitted_labels[matched], intercept = TRUE, env = baseenv())
  }else{
    stats::as.formula("~ 1", env = baseenv())
  }
  selected <- stats::terms(selected_formula, keep.order = TRUE)
  variables <- attr(fitted, "variables")
  if(length(variables) - 1L != nrow(fitted_factors) ||
     (nrow(fitted_factors) > 0L && !identical(vapply(as.list(variables)[-1L], function(x) paste(deparse(x), collapse = " "), character(1)), rownames(fitted_factors)))){
    .bt_stop_refit_required("Fitted formula variables disagree with their term coding. Refit the model with this version of BayesTools.")
  }
  attr(selected, "variables") <- as.call(c(list(as.name("list")), as.list(variables)[needed + 1L]))
  predvars <- attr(fitted, "predvars")
  if(!is.null(predvars)){
    attr(selected, "predvars") <- as.call(c(list(as.name("list")), as.list(predvars)[needed + 1L]))
  }
  attr(selected, "factors") <- fitted_factors[needed, matched, drop = FALSE]
  attr(selected, "term.labels") <- fitted_labels[matched]
  attr(selected, "order") <- attr(fitted, "order", exact = TRUE)[matched]
  attr(selected, "intercept") <- 1L
  attr(selected, "response") <- 0L
  attr(selected, ".Environment") <- baseenv()
  data_classes <- attr(fitted, "dataClasses")
  if(!is.null(data_classes)) attr(selected, "dataClasses") <- data_classes[needed]
  columns <- which(design$assign %in% c(0L, matched))
  list(terms = selected, raw_column_names = design$raw_column_names[columns])
}
