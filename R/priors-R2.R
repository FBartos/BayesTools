# R-squared variance-budget priors for formula terms.
#
# A budget declared by prior_R2() is placed in a formula prior list under the
# reserved name "__R2". JAGS_formula() binds it to the terms of the formula:
# the total R-squared r ~ Beta(a_M, b_M) and the allocation phi ~
# Dirichlet(alpha_M) over the retained allocation units become ordinary priors
# of the formula (named after the budget), and every retained term becomes an
# R2 term container whose coefficients are deterministic functions of
# W = r / (1 - r), its allocation, and its own direction coordinates
# prior_par_z_<coefficient>. The residual scale enters through the
# 'multiply_by' attribute of the term, so coefficients are in residual-SD
# units. See ?prior_R2 for the coefficient laws and cross-model conventions.

.bt_R2_reserved_name    <- "__R2"
.bt_R2_laws             <- c("exact", "normal", "laplace")
.bt_R2_conventions      <- c("common", "linked", "conditioning")
.bt_R2_concentrations   <- c("inherited", "fixed_total")
.bt_R2_contrasts        <- c("orthonormal", "meandif", "treatment")

#' @title R-squared variance-budget prior for formula terms
#'
#' @description Declares a joint prior distribution for the coefficients of
#' several formula terms through a budget for the explained variance. The
#' total R-squared of the terms is \eqn{r \sim \textrm{Beta}(a, b)}, its
#' signal-to-noise ratio \eqn{W = r / (1 - r)} is allocated to the terms by
#' Dirichlet weights, and each term turns its share into coefficients by the
#' selected coefficient law. Place the object in a formula prior list under
#' the reserved name \code{"__R2"}; see Details.
#'
#' @param alpha named numeric vector of Dirichlet concentrations, one per term
#' of the full candidate model. The names are formula term labels (e.g.,
#' \code{c(x = 1, f = 1, "x:f" = 1)}). Terms of the budget that a model's
#' formula does not contain are deleted, and the budget is redistributed in
#' full to the retained terms.
#' @param a first shape of the Beta prior on the total R-squared. Required
#' unless \code{convention = "linked"}, which determines it from the retained
#' concentrations (and then it must be \code{NULL}). Under
#' \code{convention = "conditioning"}, \code{a} and \code{b} are the shapes of
#' the full candidate model.
#' @param b second shape of the Beta prior on the total R-squared.
#' @param law coefficient law: \code{"exact"} allocates squared signal
#' magnitudes exactly, with a fair sign for one-column terms and a uniform
#' direction in the term's coefficient space (in the empirical design metric)
#' for multi-column terms; \code{"normal"} is the normal kernel of the R2D2M2
#' prior and \code{"laplace"} the double-exponential kernel of the R2D2
#' prior, both allocating the budget over the individual coded coefficients.
#' @param convention cross-model convention for the Beta shapes in a model
#' with retained concentration total \eqn{A_M}: \code{"common"} keeps
#' \eqn{(a, b)} in every model; \code{"linked"} uses \eqn{(A_M, b)};
#' \code{"conditioning"} keeps \eqn{a_M - A_M} and \eqn{a_M + b_M} equal to
#' their values in the full candidate model.
#' @param concentration concentrations of the retained units:
#' \code{"inherited"} keeps the declared values; \code{"fixed_total"}
#' rescales them to sum to \code{total} in every model.
#' @param sigma residual standard deviation that multiplies the coefficients:
#' the name of a parameter of the model (e.g., \code{"sigma"}) or a positive
#' number.
#' @param total total concentration for \code{concentration = "fixed_total"}.
#' @param contrast contrast of the factor terms of the budget:
#' \code{"orthonormal"}, \code{"meandif"}, or \code{"treatment"}. Required
#' when the budget contains a factor term.
#' @param name name of the budget. The total R-squared is the formula
#' parameter \code{<parameter>_<name>} and the allocation the formula
#' parameter \code{<parameter>_<name>_allocation}, whose elements follow the
#' retained terms in the order of the formula (for the kernels, the coded
#' coefficients of each term in turn). Several budgets of one formula need
#' distinct names.
#'
#' @details
#' A formula prior list declares budgets with the entry
#' \code{"__R2" = prior_R2(...)}, or \code{"__R2" = list(prior_R2(...), ...)}
#' for several budgets (for example, a budget for compulsory covariates and a
#' separate budget for the tested terms, whose R-squared is then the
#' incremental partial R-squared of the tested terms). Every term named in a
#' budget's \code{alpha} must not have its own prior in the list, except a
#' point prior at zero (\code{prior("spike", list(0))} or
#' \code{prior_none()}), which excludes it from the model and from the
#' budget. Terms that no budget names keep their own priors.
#'
#' In a model \eqn{M} with retained allocation units \eqn{u} (the terms for
#' \code{law = "exact"}, the coded coefficients for the kernels, each with
#' its term's concentration) the prior is
#' \deqn{r \sim \textrm{Beta}(a_M, b_M), \quad W = r / (1 - r), \quad
#' \phi \sim \textrm{Dirichlet}(\alpha_M).}
#' With the centered design block \eqn{X_t} of a term and its empirical metric
#' \eqn{G_t = X_t^\top C X_t / n}, the coefficients in residual-SD units are
#' \itemize{
#'   \item \code{"exact"}: \eqn{\beta_t = \sqrt{W \phi_t}\, R_t^{-1} u_t},
#'   with \eqn{R_t^\top R_t = G_t} and \eqn{u_t} uniform on the unit sphere,
#'   so the term's signal variance \eqn{\beta_t^\top G_t \beta_t} equals
#'   \eqn{W \phi_t} for every draw and the law of its fitted contribution does
#'   not depend on the contrast;
#'   \item \code{"normal"}: \eqn{\beta_j \sim N(0, W \phi_j / G_{t,jj})};
#'   \item \code{"laplace"}: \eqn{\beta_j \sim
#'   \textrm{Laplace}(0, \sqrt{W \phi_j / (2 G_{t,jj})})}.
#' }
#' When the term blocks are mutually orthogonal, the realized R-squared of the
#' budget equals \eqn{r} for every draw under \code{"exact"}. With correlated
#' terms, the signal variance of the whole budget differs from \eqn{W} by the
#' cross products of the terms' contributions. A model with a single
#' retained unit has \eqn{\phi = 1}; a model whose budget retains no term has
#' no signal from it.
#'
#' The prior targets the empirical variance of the fitted contributions,
#' computed with denominator \eqn{n} in the fitted data; unequal cell counts
#' therefore weight the levels of a factor by their frequencies. With
#' treatment contrasts, the intercept prior applies to the reference level
#' rather than to the grand mean.
#'
#' @return \code{prior_R2} returns an object of class
#' \code{"BayesTools_R2_budget"}.
#'
#' @examples
#' budget <- prior_R2(
#'   alpha = c(x = 1, f = 1), a = 1, b = 1, law = "exact",
#'   convention = "common", concentration = "inherited", sigma = "sigma",
#'   contrast = "orthonormal"
#' )
#' budget
#'
#' # formula prior list of a model 'mu ~ x + f'
#' formula_prior_list <- list(
#'   intercept = prior("normal", list(0, 10)),
#'   "__R2"    = budget
#' )
#'
#' \dontrun{
#' # ANCOVA: compare the models with and without the factor f
#' set.seed(1)
#' data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 20)), x = rnorm(60))
#' data$y <- 0.5 * (data$f == "b") + 0.3 * data$x + rnorm(60)
#'
#' fit_model <- function(formula){
#'   arguments <- list(
#'     model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
#'     data = list(y = data$y, N = nrow(data)),
#'     prior_list = list(sigma = prior("invgamma", list(1, 1))),
#'     formula_list = list(mu = formula),
#'     formula_data_list = list(mu = data),
#'     formula_prior_list = list(mu = formula_prior_list)
#'   )
#'   fit <- do.call(JAGS_fit, arguments)
#'   marglik <- do.call(JAGS_bridgesampling, c(list(fit = fit,
#'     log_posterior = function(parameters, data){
#'       sum(dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
#'     }), arguments[-1]))
#'   list(fit = fit, marglik = marglik, prior_weights = 1)
#' }
#'
#' # the same budget in both models: without f, x receives the whole budget
#' models <- models_inference(list(fit_model(~ x), fit_model(~ x + f)))
#' inference <- ensemble_inference(models, parameters = "mu_f",
#'   is_null_list = list(mu_f = c(TRUE, FALSE)), conditional = FALSE)
#' inference$mu_f$BF
#'
#' # posterior summaries in residual-SD units, with factor levels
#' runjags_estimates_table(models[[2]]$fit, transform_factors = TRUE)
#'
#' # implied prior distributions: the total R-squared and a coefficient
#' priors <- attr(models[[2]]$fit, "prior_list")
#' plot(priors$mu_R2)
#' plot(priors$mu_x)
#' }
#'
#' @seealso [JAGS_formula()], [prior_factor()]
#' @export
prior_R2 <- function(alpha, a = NULL, b, law, convention, concentration, sigma,
                     total = NULL, contrast = NULL, name = "R2"){

  required <- c(alpha = missing(alpha), b = missing(b), law = missing(law),
                convention = missing(convention),
                concentration = missing(concentration), sigma = missing(sigma))
  if(any(required)){
    stop("The ", paste0("'", names(required)[required], "'", collapse = ", "),
         " argument", if(sum(required) > 1) "s are" else " is",
         " required: the R-squared budget prior has no default for ",
         if(sum(required) > 1) "them" else "it", ".", call. = FALSE)
  }

  check_char(law, "law", allow_values = .bt_R2_laws, allow_NA = FALSE)
  check_char(convention, "convention", allow_values = .bt_R2_conventions, allow_NA = FALSE)
  check_char(concentration, "concentration", allow_values = .bt_R2_concentrations, allow_NA = FALSE)
  check_char(contrast, "contrast", allow_values = .bt_R2_contrasts, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(name, "name", allow_NA = FALSE)
  .bt_check_jags_node_name(name, "name")
  if(grepl("intercept", name, fixed = TRUE)){
    stop("The budget 'name' cannot contain 'intercept'.", call. = FALSE)
  }

  if(!is.numeric(alpha) || length(alpha) == 0L || is.null(names(alpha)) ||
     anyNA(names(alpha)) || any(!nzchar(trimws(names(alpha))))){
    stop("'alpha' must be a named numeric vector with one concentration per formula term.", call. = FALSE)
  }
  check_real(alpha, "alpha", lower = 0.01, check_length = 0, allow_NA = FALSE)
  term_keys <- .bt_R2_term_key(names(alpha))
  if(anyDuplicated(term_keys)){
    stop("'alpha' names each formula term once; duplicated: ",
         paste0("'", names(alpha)[duplicated(term_keys)], "'", collapse = ", "), ".", call. = FALSE)
  }
  if(any(grepl("intercept", names(alpha), fixed = TRUE))){
    stop("The intercept cannot be part of an R-squared budget.", call. = FALSE)
  }
  alpha <- stats::setNames(as.numeric(alpha), trimws(names(alpha)))

  check_real(b, "b", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if(identical(convention, "linked")){
    if(!is.null(a)){
      stop("'a' must be NULL with convention = \"linked\", which sets it to the total retained concentration of each model.", call. = FALSE)
    }
  }else{
    if(is.null(a)){
      stop("The 'a' argument is required with convention = \"", convention, "\".", call. = FALSE)
    }
    check_real(a, "a", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  }

  if(identical(concentration, "fixed_total")){
    if(is.null(total)){
      stop("The 'total' argument is required with concentration = \"fixed_total\".", call. = FALSE)
    }
    check_real(total, "total", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  }else if(!is.null(total)){
    stop("'total' applies only to concentration = \"fixed_total\".", call. = FALSE)
  }

  if(is.character(sigma)){
    check_char(sigma, "sigma", allow_NA = FALSE)
    .bt_check_jags_node_name(sigma, "sigma")
  }else{
    check_real(sigma, "sigma", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  }

  if(identical(convention, "conditioning") && identical(concentration, "inherited") &&
     identical(law, "exact") && a - sum(alpha) + min(alpha) <= 0){
    stop("The 'conditioning' convention needs a - sum(alpha) + min(alpha) > 0, so that every model has a positive first shape.", call. = FALSE)
  }

  output <- list(
    alpha         = alpha,
    a             = a,
    b             = b,
    law           = law,
    convention    = convention,
    concentration = concentration,
    total         = total,
    sigma         = sigma,
    contrast      = contrast,
    name          = name
  )
  class(output) <- c("BayesTools_R2_budget", "list")
  output
}

#' @export
print.BayesTools_R2_budget <- function(x, ...){

  law_label <- switch(
    x$law,
    exact   = "exact allocation",
    normal  = "normal kernel (R2D2M2)",
    laplace = "Laplace kernel (R2D2)"
  )
  shapes <- switch(
    x$convention,
    common       = paste0("Beta(", format(x$a), ", ", format(x$b), ") in every model"),
    linked       = paste0("Beta(A_M, ", format(x$b), "), A_M the retained concentration"),
    conditioning = paste0("Beta(", format(x$a), ", ", format(x$b),
                          ") in the full model; a - A and a + b fixed across models")
  )
  concentration <- if(identical(x$concentration, "fixed_total")){
    paste0("rescaled to total ", format(x$total))
  }else{
    "inherited"
  }
  cat("R-squared budget '", x$name, "' (", law_label, ")\n", sep = "")
  cat("  R2 ~ ", shapes, "\n", sep = "")
  cat("  allocation: Dirichlet(", paste0(names(x$alpha), " = ", format(x$alpha), collapse = ", "),
      "), ", concentration, " concentrations\n", sep = "")
  cat("  coefficients multiplied by ", if(is.character(x$sigma)) x$sigma else format(x$sigma), "\n", sep = "")
  if(!is.null(x$contrast)){
    cat("  factor contrast: ", x$contrast, "\n", sep = "")
  }
  invisible(x)
}

#' @title Check whether an object is an R-squared budget or an R2 term prior
#'
#' @param x an object.
#' @return \code{is.prior.R2} returns \code{TRUE} for the term priors that
#' \code{JAGS_formula()} binds from a [prior_R2()] budget, and
#' \code{is.R2_budget} for budgets created by [prior_R2()].
#' @export
is.prior.R2 <- function(x){
  inherits(x, "prior.R2")
}
#' @rdname is.prior.R2
#' @export
is.R2_budget <- function(x){
  inherits(x, "BayesTools_R2_budget")
}

.bt_R2_stop <- function(message, class = "BayesTools_R2_invalid"){
  stop(structure(list(message = message, call = NULL),
    class = unique(c(class, "BayesTools_R2_condition", "error", "condition"))))
}

# A coefficient state outside the support of the budget (a total outside
# [0, 1) or invalid allocation coordinates): bridge sampling treats it as a
# state of zero density.
.bt_R2_stop_state <- function(message){
  .bt_R2_stop(message, class = c("BayesTools_R2_invalid_state", "BayesTools_marglik_out_of_support"))
}

.bt_R2_unavailable <- function(message, class = NULL){
  .bt_R2_stop(message, class = c(class, "BayesTools_R2_unavailable"))
}

# Canonical key of a term label: components sorted, so that 'f:x' and 'x:f'
# name the same term.
.bt_R2_term_key <- function(labels){

  labels <- gsub("__xXx__", ":", labels, fixed = TRUE)
  vapply(labels, function(label){
    paste(sort(trimws(strsplit(label, ":", fixed = TRUE)[[1]])), collapse = ":")
  }, character(1), USE.NAMES = FALSE)
}

.bt_R2_direction_name <- function(coefficient_name){
  paste0("prior_par_z_", coefficient_name)
}

# ---------------------------------------------------------------------------
# Formula binding

# The budgets of a formula prior list (NULL when there are none).
.bt_R2_formula_budgets <- function(entry){

  if(is.null(entry)){
    return(NULL)
  }
  budgets <- if(is.R2_budget(entry)) list(entry) else entry
  if(!is.list(budgets) || length(budgets) == 0L || !all(vapply(budgets, is.R2_budget, logical(1)))){
    .bt_R2_stop("The 'prior_list[[\"__R2\"]]' entry must be a prior_R2() budget or a list of prior_R2() budgets.")
  }
  budget_names <- vapply(budgets, `[[`, character(1), "name")
  if(anyDuplicated(budget_names)){
    .bt_R2_stop(paste0("The R-squared budgets of a formula need distinct names; duplicated: ",
      paste0("'", unique(budget_names[duplicated(budget_names)]), "'", collapse = ", "), "."))
  }
  keys <- unlist(lapply(budgets, function(budget) .bt_R2_term_key(names(budget$alpha))))
  if(anyDuplicated(keys)){
    .bt_R2_stop(paste0("A formula term can belong to only one R-squared budget; repeated: ",
      paste0("'", unique(keys[duplicated(keys)]), "'", collapse = ", "), "."))
  }
  names(budgets) <- budget_names
  budgets
}

.bt_R2_is_exclusion_prior <- function(prior){

  if(is.prior.none(prior)){
    return(TRUE)
  }
  location <- prior$parameters[["location"]]
  is.prior.point(prior) && !.is_prior_expression(prior) && is.numeric(location) &&
    length(location) >= 1L && all(location == 0)
}

# An excluded factor term keeps a zero spike coded by the budget's contrast
# (a scalar spike or prior_none() is not a factor prior).
.bt_R2_exclusion_prior <- function(prior, budget, term, type){

  if(!identical(type, "factor") || is.prior.factor(prior)){
    return(prior)
  }
  if(is.null(budget$contrast)){
    .bt_R2_stop(paste0("The factor term '", term, "' of the R-squared budget '", budget$name,
      "' is excluded by a scalar point prior; give it a factor spike (prior_factor(\"spike\", list(0), contrast = ...)) ",
      "or set the budget's 'contrast'."))
  }
  prior_factor("spike", list(0), contrast = budget$contrast, prior_weights = .prior_model_weight(prior))
}

# Assigns the formula terms to the budgets. Retained budget terms get an
# unbound R2 term prior in 'prior_list' (so that defaults do not fill them);
# terms with a point prior at zero stay excluded.
.bt_R2_assign_terms <- function(budgets, prior_list, model_terms, model_terms_type){

  terms <- setdiff(model_terms, "intercept")
  term_keys <- .bt_R2_term_key(terms)
  for(budget in budgets){
    budget_keys <- .bt_R2_term_key(names(budget$alpha))
    for(i in which(term_keys %in% budget_keys)){
      term <- terms[[i]]
      declared <- names(budget$alpha)[match(term_keys[[i]], budget_keys)]
      if(term %in% names(prior_list)){
        if(.bt_R2_is_exclusion_prior(prior_list[[term]])){
          prior_list[[term]] <- .bt_R2_exclusion_prior(prior_list[[term]], budget, term, model_terms_type[[term]])
          next
        }
        .bt_R2_stop(paste0("The formula term '", term, "' belongs to the R-squared budget '", budget$name,
          "' and cannot have its own prior. Remove its prior, or give it a point prior at zero ",
          "(prior(\"spike\", list(0)) or prior_none()) to exclude it from the model and the budget."))
      }
      type <- model_terms_type[[term]]
      if(identical(type, "factor") && is.null(budget$contrast)){
        .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' contains the factor term '", term,
          "' and therefore needs a 'contrast' (\"orthonormal\", \"meandif\", or \"treatment\")."))
      }
      prior_list[[term]] <- .bt_R2_term_stub(budget, term, declared, type)
    }
  }
  prior_list
}

.bt_R2_term_stub <- function(budget, term, declared, type){

  output <- list(
    distribution = "R2",
    parameters   = list(),
    truncation   = list(lower = -Inf, upper = Inf),
    budget       = budget,
    term         = term,
    declared     = declared,
    alpha        = unname(budget$alpha[[declared]]),
    law          = budget$law,
    bound        = FALSE,
    prior_weights = 1
  )
  # Like ordered priors, R2 term priors carry the classes of the ordinary
  # priors of their term (scalar; treatment factors are scalar factor priors,
  # orthonormal and mean-difference factors vector factor priors) for the
  # consumers of term metadata (labels, coordinates, factor transforms,
  # mixing); every distribution-specific consumer dispatches on 'prior.R2'
  # first.
  class(output) <- if(identical(type, "factor")){
    c("prior", if(identical(budget$contrast, "treatment")) "prior.simple" else "prior.vector",
      "prior.factor", paste0("prior.", budget$contrast), "prior.R2", "prior.R2_factor")
  }else{
    c("prior", "prior.simple", "prior.R2")
  }
  attr(output, "multiply_by") <- budget$sigma
  output
}

# Number of design columns of every declared term of a budget in the full
# candidate model (all formula terms and all declared terms).
.bt_R2_declared_dimensions <- function(budget, formula, data){

  declared <- names(budget$alpha)
  formula_labels <- attr(stats::terms(formula), "term.labels")
  missing <- setdiff(.bt_R2_term_key(declared), .bt_R2_term_key(formula_labels))
  full_labels <- c(formula_labels, declared[.bt_R2_term_key(declared) %in% missing])
  variables <- unique(unlist(strsplit(full_labels, ":", fixed = TRUE)))
  variables <- trimws(variables)
  if(any(!variables %in% colnames(data))){
    .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' uses the 'conditioning' convention with a kernel law, ",
      "which needs the dimensions of all its declared terms; the data lack ",
      paste0("'", setdiff(variables, colnames(data)), "'", collapse = ", "), "."))
  }
  full_formula <- stats::reformulate(full_labels)
  model_frame <- stats::model.frame(full_formula, data = data, na.action = stats::na.pass)
  model_matrix <- stats::model.matrix(full_formula, model_frame)
  assign <- attr(model_matrix, "assign")
  full_keys <- .bt_R2_term_key(attr(stats::terms(full_formula), "term.labels"))
  vapply(.bt_R2_term_key(declared), function(key){
    sum(assign == match(key, full_keys))
  }, numeric(1))
}

# Beta shapes and unit concentrations of a budget in a model with the
# retained units 'unit_alpha' (declared concentrations of the units).
.bt_R2_model_shapes <- function(budget, unit_alpha, reference_total){

  alpha_M <- if(identical(budget$concentration, "fixed_total")){
    unit_alpha * budget$total / sum(unit_alpha)
  }else{
    unit_alpha
  }
  A_M <- sum(alpha_M)
  shapes <- switch(
    budget$convention,
    common = c(budget$a, budget$b),
    linked = c(A_M, budget$b),
    conditioning = {
      A_ref <- if(identical(budget$concentration, "fixed_total")) budget$total else reference_total
      c(budget$a - (A_ref - A_M), budget$b + (A_ref - A_M))
    }
  )
  if(!all(is.finite(shapes)) || any(shapes <= 0)){
    .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' gives the Beta shapes (",
      format(shapes[1]), ", ", format(shapes[2]), ") in this model; both must be positive. ",
      "Increase 'a' or change the 'convention'."))
  }
  list(a = shapes[[1]], b = shapes[[2]], alpha = alpha_M, A = A_M)
}

# Centered empirical metric of a design block.
.bt_R2_metric <- function(X){

  X <- as.matrix(X)
  centered <- sweep(X, 2, colMeans(X), "-")
  crossprod(centered) / nrow(X)
}

# Binds the R2 term stubs of a formula once its design matrix is known:
# metric, coefficient maps, budget shapes, and the injected total and
# allocation priors. Called before the prior names get the parameter prefix.
.bt_R2_bind_formula <- function(prior_list, budgets, parameter, model_matrix,
                                terms_indexes, model_terms, formula, data){

  r2_terms <- names(prior_list)[vapply(prior_list, is.prior.R2, logical(1))]
  for(budget in budgets){

    budget_terms <- r2_terms[vapply(prior_list[r2_terms], function(prior) identical(prior$budget$name, budget$name), logical(1))]
    if(length(budget_terms) == 0L){
      next
    }
    if(any(c(budget$name, paste0(budget$name, "_allocation")) %in% names(prior_list))){
      .bt_R2_stop(paste0("The R-squared budget name '", budget$name, "' clashes with a formula term; choose another 'name'."))
    }
    if(identical(budget$law, "exact")){
      unit_reference <- sum(budget$alpha)
    }else if(identical(budget$convention, "conditioning") && identical(budget$concentration, "inherited")){
      dims <- .bt_R2_declared_dimensions(budget, formula, data)
      unit_reference <- sum(budget$alpha * dims)
    }else{
      unit_reference <- NA_real_
    }

    maps <- list()
    unit_alpha <- numeric()
    unit_term <- character()
    for(term in budget_terms){
      columns <- which(terms_indexes == match(term, model_terms))
      X <- model_matrix[, columns, drop = FALSE]
      metric <- .bt_R2_metric(X)
      d <- ncol(X)
      label <- gsub("__xXx__", ":", term, fixed = TRUE)
      scales <- sqrt(diag(metric))
      if(any(!is.finite(scales)) || any(scales <= 0)){
        .bt_R2_stop(paste0("A design column of the R-squared budget term '", label,
          "' is constant after centering; its signal variance cannot carry a share of the budget."))
      }
      if(identical(budget$law, "exact")){
        # rank is judged on the correlation matrix, so that it does not depend
        # on the units of the predictors
        correlation <- metric / (scales %o% scales)
        R_correlation <- tryCatch(chol(correlation), error = function(e) NULL)
        if(is.null(R_correlation) || min(diag(R_correlation)) <= sqrt(.Machine$double.eps)){
          .bt_R2_stop(paste0("The design of the R-squared budget term '", label,
            "' is rank deficient after centering; the exact law needs linearly independent centered columns."))
        }
        map <- backsolve(sweep(R_correlation, 2, scales, "*"), diag(d))
        units <- 1L
      }else{
        map <- diag(1 / scales, nrow = d)
        units <- d
      }
      maps[[term]] <- list(map = map, metric = metric, d = d, units = units)
      unit_alpha <- c(unit_alpha, rep(prior_list[[term]]$alpha, units))
      unit_term <- c(unit_term, rep(term, units))
    }

    shapes <- .bt_R2_model_shapes(budget, unit_alpha, unit_reference)
    n_units <- length(unit_alpha)
    total_name <- paste0(parameter, "_", budget$name)
    allocation_name <- if(n_units > 1L) paste0(parameter, "_", budget$name, "_allocation") else NA_character_

    total_prior <- prior("beta", list(alpha = shapes$a, beta = shapes$b))
    attr(total_prior, "R2_budget") <- list(name = budget$name, role = "total", law = budget$law)
    prior_list[[budget$name]] <- total_prior
    if(n_units > 1L){
      allocation_prior <- prior("dirichlet", list(alpha = shapes$alpha))
      attr(allocation_prior, "R2_budget") <- list(name = budget$name, role = "allocation", law = budget$law,
        units = .bt_R2_unit_labels(unit_term, maps))
      prior_list[[paste0(budget$name, "_allocation")]] <- allocation_prior
    }

    for(term in budget_terms){
      this_prior <- prior_list[[term]]
      this_prior$bound <- TRUE
      if(is.prior.vector(this_prior)){
        this_prior$parameters <- list(K = maps[[term]]$d)
      }
      this_prior$parameter <- parameter
      this_prior$d <- maps[[term]]$d
      this_prior$map <- maps[[term]]$map
      this_prior$metric <- maps[[term]]$metric
      this_prior$units <- which(unit_term == term)
      this_prior$n_units <- n_units
      this_prior$shapes <- shapes
      this_prior$total_name <- total_name
      this_prior$allocation_name <- allocation_name
      this_prior$budget_terms <- budget_terms
      prior_list[[term]] <- this_prior
    }
  }
  prior_list
}

# The total and allocation priors that the coefficient node of a term reads.
.bt_R2_budget_prior_list <- function(prior){

  .bt_R2_require_bound(prior)
  out <- list()
  out[[prior$total_name]] <- prior("beta", list(alpha = prior$shapes$a, beta = prior$shapes$b))
  if(prior$n_units > 1L){
    out[[prior$allocation_name]] <- prior("dirichlet", list(alpha = prior$shapes$alpha))
  }
  out
}

.bt_R2_unit_labels <- function(unit_term, maps){

  labels <- gsub("__xXx__", ":", unit_term, fixed = TRUE)
  for(term in unique(unit_term)){
    if(maps[[term]]$units > 1L){
      labels[unit_term == term] <- paste0(labels[unit_term == term], "[", seq_len(maps[[term]]$units), "]")
    }
  }
  labels
}

# The Beta total and Dirichlet allocation priors that a budget adds to a
# formula prior list (they are parameters of the formula, not terms).
.bt_R2_is_budget_prior <- function(prior){
  !is.null(attr(prior, "R2_budget", exact = TRUE))
}

.bt_R2_require_bound <- function(prior, parameter_name = NULL){

  if(!is.prior.R2(prior) || !isTRUE(prior$bound)){
    .bt_R2_unavailable(paste0(
      "The R-squared budget prior", if(!is.null(parameter_name)) paste0(" of '", parameter_name, "'"),
      " is unavailable before JAGS_formula() binds it to a formula design."),
      class = "BayesTools_R2_unbound")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Coordinates and JAGS syntax of a bound term

.bt_R2_coefficient_names <- function(prior, parameter_name){

  if(prior$d > 1L){
    paste0(parameter_name, "[", seq_len(prior$d), "]")
  }else{
    parameter_name
  }
}

is.prior.R2_factor <- function(x){
  inherits(x, "prior.R2_factor")
}

.bt_R2_direction_coordinates <- function(prior, parameter_name){

  name <- .bt_R2_direction_name(parameter_name)
  if(prior$d > 1L){
    paste0(name, "[", seq_len(prior$d), "]")
  }else{
    name
  }
}

.bt_R2_allocation_coordinates <- function(prior){

  if(prior$n_units > 1L){
    paste0(prior$allocation_name, "[", prior$units, "]")
  }else{
    character()
  }
}

.bt_R2_allocation_eta_coordinates <- function(prior){

  paste0(.JAGS_prior_dirichlet_eta_name(prior$allocation_name), "[", seq_len(prior$n_units), "]")
}

.bt_dnode_R2_coefficient <- function(parameter_name, prior){

  .bt_R2_require_bound(prior, parameter_name)
  .bt_deterministic_node(
    "R2_coefficient",
    parameter_name,
    .bt_R2_coefficient_names(prior, parameter_name),
    dependencies = c(
      prior$total_name,
      if(prior$n_units > 1L) .bt_R2_allocation_eta_coordinates(prior),
      .bt_R2_direction_coordinates(prior, parameter_name)
    ),
    parameter = parameter_name,
    spec = list(
      law             = prior$law,
      d               = prior$d,
      map             = prior$map,
      units           = prior$units,
      n_units         = prior$n_units,
      total_name      = prior$total_name,
      allocation_name = prior$allocation_name,
      directions      = .bt_R2_direction_coordinates(prior, parameter_name),
      eta             = if(prior$n_units > 1L) .bt_R2_allocation_eta_coordinates(prior) else character(),
      coefficients    = .bt_R2_coefficient_names(prior, parameter_name)
    )
  )
}

.bt_dnode_R2_coefficient_emit <- function(node){

  spec <- node$spec
  share <- function(unit){
    if(spec$n_units > 1L){
      paste0(spec$total_name, " / (1 - ", spec$total_name, ") * ", spec$allocation_name, "[", unit, "]")
    }else{
      paste0(spec$total_name, " / (1 - ", spec$total_name, ")")
    }
  }
  z <- spec$directions
  if(identical(spec$law, "exact")){
    radius <- paste0("sqrt(", share(spec$units[1]), ")")
    if(spec$d == 1L){
      return(paste0(spec$coefficients, " <- ", radius, " * ", .JAGS_prior_number(spec$map[1, 1]),
        " * ", z, " / abs(", z, ")"))
    }
    norm <- paste0("sqrt(", paste0(z, "^2", collapse = " + "), ")")
    vapply(seq_len(spec$d), function(i){
      columns <- which(spec$map[i, ] != 0)
      combination <- paste0(vapply(columns, function(j){
        paste0(.JAGS_prior_number(spec$map[i, j]), " * ", z[j])
      }, character(1)), collapse = " + ")
      paste0(spec$coefficients[i], " <- ", radius, " * (", combination, ") / ", norm)
    }, character(1))
  }else{
    vapply(seq_len(spec$d), function(i){
      paste0(spec$coefficients[i], " <- sqrt(", share(spec$units[i]), ") * ",
        .JAGS_prior_number(spec$map[i, i]), " * ", z[i])
    }, character(1))
  }
}

.bt_dnode_R2_coefficient_evaluate <- function(node, lookup){

  spec <- node$spec
  total <- .bt_deterministic_lookup_value(lookup, spec$total_name)
  if(is.null(total)){
    return(NULL)
  }
  if(spec$n_units > 1L){
    # the allocation is normalized from its gamma coordinates (the declared
    # dependencies), so that a replay with changed coordinates never reads
    # stale normalized weights
    if(!all(spec$eta %in% colnames(lookup$draws))){
      return(NULL)
    }
    eta <- lookup$draws[, spec$eta, drop = FALSE]
    eta_sum <- rowSums(eta)
    if(any(!is.finite(eta)) || any(eta < 0) || any(!is.finite(eta_sum)) || any(eta_sum <= 0)){
      .bt_R2_stop_state("R-squared budget allocation gamma coordinates must be finite and nonnegative with positive row sums.")
    }
    allocation <- (eta / eta_sum)[, spec$units, drop = FALSE]
  }else{
    allocation <- matrix(1, nrow = lookup$n, ncol = 1L)
  }
  if(!all(spec$directions %in% colnames(lookup$draws))){
    return(NULL)
  }
  z <- lookup$draws[, spec$directions, drop = FALSE]
  if(any(!is.finite(total)) || any(total < 0 | total >= 1)){
    .bt_R2_stop_state("R-squared budget totals must lie in [0, 1).")
  }
  W <- total / (1 - total)
  if(identical(spec$law, "exact")){
    norm <- sqrt(rowSums(z^2))
    direction <- z / norm
    values <- sqrt(W * allocation[, 1]) * (direction %*% t(spec$map))
  }else{
    values <- sqrt(W * allocation) * sweep(z, 2, diag(spec$map), "*")
  }
  values <- matrix(values, nrow = lookup$n, ncol = spec$d)
  values
}

.JAGS_prior.R2 <- function(prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  directions <- .bt_R2_direction_coordinates(prior, parameter_name)
  kernel <- if(identical(prior$law, "laplace")){
    paste0("ddexp(0, ", .JAGS_prior_number(sqrt(2)), ")")
  }else{
    "dnorm(0, 1)"
  }
  paste0(
    paste0(directions, " ~ ", kernel, "\n", collapse = ""),
    paste0(.bt_deterministic_node_emit(.bt_dnode_R2_coefficient(parameter_name, prior)), "\n", collapse = "")
  )
}

.JAGS_init.R2 <- function(prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  values <- .bt_R2_direction_rng(prior$law, prior$d)
  init <- list()
  init[[.bt_R2_direction_name(parameter_name)]] <- values
  init
}

.JAGS_monitor.R2 <- function(prior, parameter_name){

  c(parameter_name, .bt_R2_direction_name(parameter_name))
}

# Direction coordinates drawn from their prior: standard normal for the exact
# law and the normal kernel, unit-variance Laplace for the Laplace kernel.
.bt_R2_direction_rng <- function(law, n){

  if(identical(law, "laplace")){
    stats::rexp(n, rate = sqrt(2)) * sample(c(-1, 1), n, replace = TRUE)
  }else{
    stats::rnorm(n)
  }
}

# Log density of the direction coordinates of one draw.
.bt_R2_direction_lpdf <- function(law, z){

  if(identical(law, "laplace")){
    sum(log(sqrt(2) / 2) - sqrt(2) * abs(z))
  }else{
    sum(stats::dnorm(z, log = TRUE))
  }
}
