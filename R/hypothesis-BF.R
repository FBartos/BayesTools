# ============================================================================ #
# hypothesis-BF.R
# ============================================================================ #
#
# Generic Bayes factors for scalar prior/posterior quantities.
#
# ============================================================================ #


#' @title Hypothesis Bayes Factors
#'
#' @description Computes Bayes factors for scalar hypotheses written as
#' expressions, for example \code{"theta = 0"}, \code{"theta > 0"}, or
#' \code{"theta = 0 vs theta > 0"}. Factor/posterior-list levels can be
#' referenced as \code{"mu_alloc[alternate] > mu_alloc[random]"} when the
#' marginal posterior carries joint prior information. Level names that
#' contain brackets, such as the \code{cut()} level \code{(0,1]}, are
#' referenced by quoting the whole reference in backticks,
#' \code{"`mu[(0,1]]` > `mu[(1,2]]`"}, or with the parameter catalog's escaped
#' level form \code{mu["(0,1\%5D"]}. Point-null hypotheses use
#' \code{\link{Savage_Dickey_BF}} when a \code{marginal_posterior} object is
#' supplied, including precomputed qCMDE/IWMDE posterior ordinates when
#' \code{density_method = "precomputed"}. Level-specific
#' \code{marginal_posterior.*} subclass vectors with attached density
#' metadata are accepted as marginal posterior inputs; parent-level
#' \code{posterior_density}, \code{posterior_densities},
#' \code{posterior_ordinate}, or \code{posterior_ordinates} metadata on a
#' marginal-posterior list is reused for matching level hypotheses.
#'
#' @param posterior posterior draws, a \code{marginal_posterior}, a
#' \code{marginal_inference} object, or a data frame/matrix of posterior draws.
#' Posterior draws must be finite and quantity names must be unique.
#' @param prior prior draws for numeric/data-frame inputs, or for numeric
#' draws a scalar BayesTools prior object (including mixture and
#' spike-and-slab priors). Quantity names in draw tables must be unique.
#' Ignored when \code{posterior} already contains deterministic prior density
#' information.
#' @param hypothesis character vector with scalar hypothesis statements written
#' in the restricted grammar described in Details, or a validated object from
#' [hypothesis_parse()].
#' @param parameter optional scalar quantity name for numeric vectors or
#' \code{marginal_posterior} objects.
#' @param logBF whether to display the Bayes factor on the log scale.
#' @param BF01 whether to display the inverse Bayes factor.
#' @param seed optional finite seed used by downstream helpers that sample.
#' The caller's random-number generator state is restored after the call.
#' @param density_method posterior density source for point-null tests.
#' \code{"KDE"} uses kernel density estimates. \code{"normal"} uses a normal
#' approximation to the posterior density at the null. \code{"precomputed"}
#' requires valid \code{posterior_ordinate} or \code{posterior_density}
#' metadata ([posterior_metadata()]) and errors if no valid precomputed
#' source is available. KDE
#' point-null tests use boundary reflection when exact support metadata is
#' attached to the marginal posterior or to matched \code{posterior_density}
#' metadata; numeric and data-frame expression tests use standard Gaussian
#' sample KDE. A point hypothesis on a linear expression of marginal-posterior
#' parameters or levels (sums, differences, and constant multiples, e.g.,
#' \code{mu[B] - mu[A] = 0}; functions and powers only of constants) is tested
#' as a Savage-Dickey ratio of that linear combination: the prior density and
#' exact support come from the joint prior context and the posterior ordinate
#' is the normal approximation (\code{"normal"}) or the KDE with boundary
#' reflection at bounded supports (per mixture component, see
#' [Savage_Dickey_BF]). An affine expression of a numeric quantity with a
#' prior object uses the exact density of that transformed prior. Only
#' user-supplied prior draws use the sample KDE of the prior expression draws,
#' with a classed warning (see Value).
#' Finite sample and KDE evaluation ranges are not treated as exact support,
#' so finite point hypotheses outside those ranges use kernel-tail density
#' estimates.
#' @param columns output columns. \code{"default"} returns \code{Alternative},
#' \code{Null}, \code{BF}, and \code{BF_error}. \code{"all"} also returns
#' \code{prior}, \code{posterior}, and \code{method} columns. The
#' \code{prior} and \code{posterior} columns are diagnostics, not always
#' probabilities: point tests report density heights, region tests report odds,
#' and transitive point-vs-region tests report \code{NA}. The default
#' \code{BF_error} column is printed as \code{error\%(BF)}. Individual columns
#' are selected by these column names; other spellings are rejected.
#' @param ... unused.
#'
#' @details The hypothesis language deliberately accepts only a small,
#' non-programmable subset of R expressions:
#'
#' \preformatted{
#' hypothesis := statement [ "vs" statement ]
#' statement  := point | region
#' point      := arithmetic ( "=" | "==" | "!=" ) finite_number
#' region     := relation
#'             | "(" region ")"
#'             | "!" region
#'             | region ( "&" | "|" ) region
#' relation   := arithmetic ( "<" | "<=" | ">" | ">=" ) arithmetic
#' }
#'
#' Arithmetic expressions may contain parameter identifiers, finite numeric
#' literals, parentheses, \code{+}, \code{-}, \code{*}, \code{/}, \code{^},
#' and the one-argument functions \code{abs()}, \code{exp()}, \code{log()},
#' \code{sqrt()}, \code{plogis()}, and \code{qlogis()}. Every constant
#' arithmetic subexpression must evaluate to one finite number. A point value
#' may be one finite numeric literal, optionally preceded by one unary sign, or
#' another parameter expression. Symbolic right-hand sides are normalized to a
#' difference from zero. Relations with a constant on the left and a parameter
#' expression on the right are normalized to the equivalent parameter-left
#' form.
#'
#' Equality defines a point hypothesis and may appear only as the top-level
#' relation of a statement. It cannot be parenthesized inside, negated, or
#' combined as a region. Region relations may be parenthesized, negated, and
#' combined with the element-wise operators \code{&} and \code{|};
#' \code{&&} and \code{||} are not supported. Backticks permit
#' non-syntactic parameter names. In particular, escaped identifiers such as
#' \code{`Inf`} refer to parameters, while unescaped \code{Inf}, \code{NaN},
#' \code{NA}, \code{TRUE}, and \code{FALSE} are reserved literals and are
#' rejected.
#'
#' \tabular{lll}{
#' \strong{Form} \tab \strong{Status} \tab \strong{Reason} \cr
#' \code{theta = -0.5} \tab accepted \tab point with a finite literal \cr
#' \code{theta = phi} \tab accepted \tab normalized to \code{theta - phi = 0} \cr
#' \code{0 > theta} \tab accepted \tab normalized to \code{theta < 0} \cr
#' \code{(theta > 0)} \tab accepted \tab parenthesized region \cr
#' \code{!(theta > 0)} \tab accepted \tab negated region \cr
#' \code{theta > 0 & abs(phi) < 2} \tab accepted \tab combined regions \cr
#' \code{`Inf` > 0} \tab accepted \tab escaped parameter identifier \cr
#' \code{!(theta == 0)} \tab rejected \tab point equalities cannot be negated \cr
#' \code{theta = 1 + 1} \tab rejected \tab constant point value is not one literal \cr
#' \code{theta > 0 && phi < 1} \tab rejected \tab scalar boolean operator \cr
#' \code{sin(theta) > 0} \tab rejected \tab function outside the whitelist
#' }
#'
#' Prior region masses of deterministic prior densities use the structure that
#' [prior_density_ordinate()] classifies, for regions whose relations are
#' linear in the quantity (e.g., \code{theta > 0}, \code{2 * theta < 1}, and
#' their combinations with \code{&}, \code{|}, and \code{!}): point masses
#' contribute their exact probability, scalar and normal-sum priors their
#' exact distribution functions, and Gaussian convolutions and
#' conditional-normal scale mixtures the conditional-normal quadrature of the
#' ordinate over the other term, split at the same breakpoints, with the
#' Gaussian peak (where the ordinate uses one) at every finite region bound;
#' mixture, model, and design-row components are summed with their
#' probabilities. A quadrature that fails its diagnostics, or evaluates to
#' exactly zero, stops with an error. Other prior densities and regions use the
#' prior-density grid, which resolves the region boundaries on that grid: a
#' boundary is located by bisection within each grid cell where the condition
#' changes, but region features narrower than the grid spacing may be missed.
#'
#' @return A BayesTools table of class \code{BayesTools_hypothesis_BF}. The
#' \code{BF_error} column reports approximate relative Monte Carlo error
#' percentage when available. Region odds errors are computed on
#' \code{log(BF)} from prior/posterior region indicators using an iid
#' delta-method approximation on the flattened draws; they do not adjust for
#' MCMC autocorrelation. Prior region masses computed exactly from a prior
#' object's distribution function contribute no Monte Carlo error.
#' Point-vs-region errors combine the available
#' point-density and region-mass errors on the \code{log(BF)} scale.
#' Point-null tests require a regular (positive and finite) prior density at
#' the null value that [prior_density_ordinate()] classifies exactly from the
#' prior's structure, on every route. Otherwise they stop with a classed error
#' condition, which callers should match by class rather than by message:
#' \code{BayesTools_point_mass_at_null} (a prior point mass at the null),
#' \code{BayesTools_infinite_ordinate}, \code{BayesTools_zero_ordinate},
#' \code{BayesTools_undefined_ordinate}, or \code{BayesTools_inexact_ordinate}
#' (no exact structural ordinate, e.g., a nonlinear expression of a
#' deterministic prior, a prior-density combination evaluated only on a
#' numerical grid, a structural quadrature rejected by its diagnostics,
#' or a density grid without recorded provenance). Each of
#' these conditions also has class \code{BayesTools_hypothesis_ordinate}.
#' For a whole-factor statement expanded into several quantities, known
#' ordinate errors produce an \code{NA} row with method \code{"unavailable"}
#' and the original reason; valid siblings are computed. Direct scalar and
#' explicitly indexed tests remain strict. Known numerical prior-region errors
#' behave likewise, with classes \code{BayesTools_prior_region_mass_unavailable},
#' \code{BayesTools_prior_region_route_unavailable},
#' \code{BayesTools_prior_region_grid_unavailable}, or
#' \code{BayesTools_prior_region_probability_rejected}, all also
#' \code{BayesTools_hypothesis_region}. Missing deterministic provenance for a
#' region remains a metadata error, including in whole-factor tables.
#' The specific \code{BayesTools_transformation_image_unavailable} leaf also
#' produces an unavailable whole-factor row with its reason; other
#' transformation and input errors remain strict. Scalar tests retain the
#' original refusal.
#' User-supplied prior draws (numeric or data-frame inputs without a prior
#' object) have no structural prior density: their prior ordinate is the
#' kernel (or normal) estimate of the prior expression draws, returned with a
#' warning of classes \code{BayesTools_inexact_ordinate} and
#' \code{BayesTools_hypothesis_ordinate}.
#' Region tests require positive prior mass for every compared region. An
#' implicit region statement is compared with its complement, which therefore
#' also needs positive prior mass; explicit comparisons such as
#' \code{"theta > 0.5 vs theta > 0"} accept an encompassing region with prior
#' mass one. Rows are labelled by quantity; when several statements refer to
#' the same quantity, the statement number is appended, e.g.
#' \code{theta (2)}.
#' A symbolic comparison on both sides uses its centered scalar expression:
#' \code{"theta = phi vs theta > phi"} and
#' \code{"theta - phi = 0 vs theta - phi > 0"} are supported. Reversed spelling
#' \code{"theta = phi vs phi < theta"} remains conservatively unavailable;
#' use the supported orientation or centered spelling.
#' A reference to a quantity or level that \code{posterior} does not contain,
#' and a \code{parameter} that it does not contain, stop with an error of
#' class \code{BayesTools_parameter_not_found} (also
#' \code{BayesTools_parameter_resolution_error}), as unresolved selectors of
#' [parameter_catalog_resolve()] do; its fields \code{alias} and
#' \code{available} name the unresolved and the available references.
#'
#' @export
hypothesis_BF <- function(posterior, prior = NULL, hypothesis, parameter = NULL,
                          logBF = FALSE, BF01 = FALSE, seed = NULL,
                          density_method = c("KDE", "normal", "precomputed"),
                          columns = "default", ...) {

  hypothesis_ast <- if(inherits(hypothesis, "BayesTools_hypothesis_ast")){
    .bt_validate_hypothesis_ast(hypothesis)
    hypothesis
  }else{
    hypothesis_parse(hypothesis)
  }
  check_char(parameter, "parameter", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01, "BF01", allow_NA = FALSE)
  check_real(seed, "seed", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  if(!is.null(seed) && !is.finite(seed)){
    stop("'seed' must be finite.", call. = FALSE)
  }
  check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  density_method <- .hypothesis_density_method(density_method)
  columns        <- .hypothesis_BF_output_columns(columns)
  .hypothesis_validate_input_integrity(posterior, prior)

  seed_state <- .bt_formula_prediction_seed(seed)
  on.exit(.bt_formula_prediction_restore_seed(seed_state), add = TRUE)

  statements <- hypothesis_ast$statements
  quantities <- .as_hypothesis_quantities(
    posterior  = posterior,
    prior      = prior,
    statements = statements,
    parameter  = parameter
  )

  rows <- list()
  log_results <- numeric()
  numerical_diagnostics <- prior_numerical_diagnostics <- list()
  row_labels <- character()
  row_statements <- integer()
  inexact_priors <- character()
  row_i <- 1L
  for(hyp_i in seq_along(statements)){
    for(quantity_i in seq_along(quantities)){
      result <- .hypothesis_BF_compute(
        quantity       = quantities[[quantity_i]],
        statement      = statements[[hyp_i]],
        density_method = density_method,
        allow_unavailable = length(quantities) > 1L
      )
      log_results[[row_i]] <- if(is.null(result$log_BF)) log(result$BF) else result$log_BF
      numerical_diagnostics[row_i] <- list(result$numerical_diagnostics)
      prior_numerical_diagnostics[row_i] <- list(result$prior_numerical_diagnostics)
      rows[[row_i]] <- .hypothesis_BF_row(
        quantity = quantities[[quantity_i]],
        result   = result
      )
      inexact_priors <- c(inexact_priors, result[["inexact_prior"]])
      row_labels[[row_i]]     <- quantities[[quantity_i]][["label"]]
      row_statements[[row_i]] <- hyp_i
      row_i <- row_i + 1L
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- .hypothesis_BF_row_names(row_labels, row_statements)
  raw_BF <- out[["BF"]]
  out[["BF"]] <- .format_BF_from_log(log_results, logBF = logBF, BF01 = BF01,
    BF = raw_BF, diagnostics = numerical_diagnostics)
  attr(out[["BF_error"]], "name") <- "error%(BF)"

  warnings <- .hypothesis_BF_table_warnings(out)
  out      <- out[, columns, drop = FALSE]

  attr(out, "raw_BF")   <- raw_BF
  attr(out, "raw_log_BF") <- log_results
  attr(out, "numerical_diagnostics") <- numerical_diagnostics
  attr(out, "prior_numerical_diagnostics") <- prior_numerical_diagnostics
  attr(out, "hypothesis_ast") <- hypothesis_ast
  attr(out, "logBF")    <- logBF
  attr(out, "BF01")     <- BF01
  attr(out, "type")      <- .hypothesis_BF_table_types(colnames(out))
  attr(out, "footnotes") <- .hypothesis_BF_table_footnotes(colnames(out))
  if(any(!vapply(prior_numerical_diagnostics, is.null, logical(1)))){
    attr(out, "footnotes") <- c(attr(out, "footnotes"),
      "Note: Deterministic prior-region diagnostics are stored in 'prior_numerical_diagnostics'; 'BF_error' reports Monte Carlo error only.")
  }
  attr(out, "warnings")  <- warnings
  attr(out, "rownames")  <- TRUE
  class(out) <- c("BayesTools_table", "BayesTools_hypothesis_BF", "data.frame")

  if(length(inexact_priors) > 0L){
    .hypothesis_warn_inexact_ordinate(unique(inexact_priors))
  }

  return(out)
}


.hypothesis_density_method <- function(density_method){

  return(posterior_density_method_match(
    density_method,
    allowed = c("KDE", "normal", "precomputed"),
    name    = "density_method"
  ))
}


#' @export
print.BayesTools_hypothesis_BF <- function(x, ...) {

  if(!inherits(x, "BayesTools_table")){
    class(x) <- c("BayesTools_table", class(x))
  }
  print(x, ...)

  return(invisible(x))
}

#' @rdname hypothesis_BF
#' @param x a hypothesis Bayes factor table.
#' @param row.names optional row names for the complete exported frame.
#' @param optional logical; exported column names are always syntactic and unique.
#' @details Data-frame coercion returns a plain long frame with leading
#' \code{component}, \code{parameter} (the displayed row identity), and
#' \code{warning} fields. Statistical values retain full precision. Bayes-factor
#' names are \code{BF10}, \code{BF01}, \code{logBF10}, or \code{logBF01}
#' according to the displayed column attributes; \code{BF_error} retains its
#' existing percent relative-error units. Declared nonmissing bounds are in
#' \code{BF_bound_operator}. Unmatched/global warnings have component
#' \code{"hypothesis/warnings"}. Compact printing retains its existing layout.
#' @exportS3Method
as.data.frame.BayesTools_hypothesis_BF <- function(x, row.names = NULL,
                                                  optional = FALSE, ...){

  parameters <- rownames(x)
  warnings <- attr(x, "warnings", exact = TRUE)
  warning_names <- names(warnings)
  if(is.null(warning_names)) warning_names <- rep(NA_character_, length(warnings))
  row_warnings <- vapply(parameters, function(parameter){
    matches <- !is.na(warning_names) & warning_names == parameter
    .hypothesis_collapse_warning(warnings[matches])
  }, character(1))
  out <- data.frame(component = rep("hypothesis", nrow(x)),
                    parameter = parameters, warning = unname(row_warnings),
                    stringsAsFactors = FALSE)
  visible <- lapply(x, function(column){
    if(inherits(column, "BayesTools_BF")) column <- as.numeric(column)
    for(attribute in c("name", "logBF", "BF01", "bound_operator")){
      attr(column, attribute) <- NULL
    }
    column
  })
  visible_names <- names(x)
  bf <- which(visible_names == "BF")
  bf_name <- character()
  bound_index <- integer()
  if(length(bf) == 1L){
    bf_column <- x[[bf]]
    visible_names[bf] <- paste0(if(isTRUE(attr(bf_column, "logBF"))) "log" else "",
                                if(isTRUE(attr(bf_column, "BF01"))) "BF01" else "BF10")
    bf_name <- visible_names[bf]
    bound <- attr(bf_column, "bound_operator", exact = TRUE)
    if(!is.null(bound) && any(!is.na(bound) & nzchar(bound))){
      visible <- c(visible, list(bound))
      visible_names <- c(visible_names, "BF_bound_operator")
      bound_index <- length(visible)
    }
  }
  other_columns <- setdiff(seq_along(visible), c(bf, bound_index))
  normalized <- make.names(c(names(out), bf_name, visible_names[bound_index],
                             visible_names[other_columns]), unique = TRUE)
  visible_names[other_columns] <- tail(normalized, length(other_columns))
  names(visible) <- visible_names
  for(name in names(visible)) out[[name]] <- visible[[name]]
  unmatched <- is.na(warning_names) | !nzchar(warning_names) |
    !warning_names %in% parameters
  if(any(unmatched)){
    extra <- out[rep(NA_integer_, sum(unmatched)), , drop = FALSE]
    extra$component <- "hypothesis/warnings"
    extra$parameter <- warning_names[unmatched]
    extra$parameter[is.na(extra$parameter) | !nzchar(extra$parameter)] <- NA_character_
    extra$warning <- as.character(warnings[unmatched])
    out <- rbind(out, extra)
  }
  rownames(out) <- row.names
  out
}
