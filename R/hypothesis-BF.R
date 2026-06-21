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
#' marginal posterior carries joint prior information. Point-null hypotheses use
#' \code{\link{Savage_Dickey_BF}} when a \code{marginal_posterior} object is
#' supplied, including precomputed qCMDE/IWMDE posterior ordinates when
#' \code{density_method = "precomputed"}. Level-specific
#' \code{marginal_posterior.*} subclass vectors with attached density
#' attributes are accepted as marginal posterior inputs; parent-level
#' \code{posterior_density}, \code{posterior_densities},
#' \code{posterior_ordinate}, or \code{posterior_ordinates} metadata on a
#' marginal-posterior list is reused for matching level hypotheses.
#'
#' @param posterior posterior draws, a \code{marginal_posterior}, a
#' \code{marginal_inference} object, or a data frame/matrix of posterior draws.
#' @param prior prior draws for numeric/data-frame inputs. Ignored when
#' \code{posterior} already contains deterministic prior density information.
#' @param hypothesis character vector with scalar hypothesis statements.
#' @param parameter optional scalar quantity name for numeric vectors or
#' \code{marginal_posterior} objects.
#' @param logBF whether to display the Bayes factor on the log scale.
#' @param BF01 whether to display the inverse Bayes factor.
#' @param seed optional seed used only by downstream helpers that sample.
#' @param density_method posterior density source for point-null tests.
#' \code{"KDE"} uses kernel density estimates. \code{"normal"} uses a normal
#' approximation to the posterior density at the null. \code{"precomputed"}
#' requires valid \code{posterior_ordinate} or \code{posterior_density}
#' attributes and errors if no valid precomputed source is available. KDE
#' point-null tests use boundary reflection when exact support metadata is
#' attached to the marginal posterior or to a matched \code{posterior_density}
#' attribute; numeric and data-frame expression tests use standard sample KDE.
#' @param columns output columns. \code{"default"} returns \code{Alternative},
#' \code{Null}, \code{BF}, and \code{BF_error}. \code{"all"} also returns
#' \code{prior}, \code{posterior}, and \code{method} columns. The
#' \code{prior} and \code{posterior} columns are diagnostics, not always
#' probabilities: point tests report density heights, region tests report odds,
#' and transitive point-vs-region tests report \code{NA}. The default
#' \code{BF_error} column is printed as \code{error\%(BF)}.
#' @param ... unused.
#'
#' @return A BayesTools table of class \code{BayesTools_hypothesis_BF}. The
#' \code{BF_error} column reports approximate relative Monte Carlo error
#' percentage when available. Region odds errors are computed on
#' \code{log(BF)} from prior/posterior region indicators using an iid
#' delta-method approximation on the flattened draws; they do not adjust for
#' MCMC autocorrelation. Point-vs-region errors combine the available
#' point-density and region-mass errors on the \code{log(BF)} scale.
#' Point-null tests require a positive finite prior density at the null value.
#'
#' @export
hypothesis_BF <- function(posterior, prior = NULL, hypothesis, parameter = NULL,
                          logBF = FALSE, BF01 = FALSE, seed = NULL,
                          density_method = c("KDE", "normal", "precomputed"),
                          columns = "default", ...) {

  check_char(hypothesis, "hypothesis", check_length = 0, allow_NA = FALSE)
  check_char(parameter, "parameter", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01, "BF01", allow_NA = FALSE)
  check_real(seed, "seed", check_length = 1, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(columns, "columns", check_length = 0, allow_NA = FALSE)
  density_method <- .hypothesis_density_method(density_method)
  columns        <- .hypothesis_BF_output_columns(columns)

  if(!is.null(seed)){
    set.seed(seed)
  }

  parsed <- lapply(hypothesis, .parse_hypothesis_BF)
  quantities <- .as_hypothesis_quantities(
    posterior  = posterior,
    prior      = prior,
    parsed     = parsed,
    parameter  = parameter
  )

  rows <- list()
  row_i <- 1L
  for(hyp_i in seq_along(parsed)){
    for(quantity_i in seq_along(quantities)){
      result <- .hypothesis_BF_compute(
        quantity       = quantities[[quantity_i]],
        parsed         = parsed[[hyp_i]],
        density_method = density_method
      )
      rows[[row_i]] <- .hypothesis_BF_row(
        quantity = quantities[[quantity_i]],
        parsed   = parsed[[hyp_i]],
        result   = result
      )
      row_i <- row_i + 1L
    }
  }

  out <- do.call(rbind, rows)
  raw_BF <- out[["BF"]]
  out[["BF"]] <- format_BF(raw_BF, logBF = logBF, BF01 = BF01)
  attr(out[["BF_error"]], "name") <- "error%(BF)"

  warnings <- .hypothesis_BF_table_warnings(out)
  out      <- out[, columns, drop = FALSE]

  attr(out, "raw_BF")   <- raw_BF
  attr(out, "parsed")   <- parsed
  attr(out, "logBF")    <- logBF
  attr(out, "BF01")     <- BF01
  attr(out, "type")      <- .hypothesis_BF_table_types(colnames(out))
  attr(out, "footnotes") <- .hypothesis_BF_table_footnotes(colnames(out))
  attr(out, "warnings")  <- warnings
  attr(out, "rownames")  <- TRUE
  class(out) <- c("BayesTools_table", "BayesTools_hypothesis_BF", "data.frame")

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


