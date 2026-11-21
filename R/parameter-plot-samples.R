#' Plot-ready samples of one selected catalog quantity
#'
#' Prepare a scalar quantity selected by [parameter_catalog_resolve()] for
#' [plot_posterior()]. The list key is the canonical name. The selected
#' posterior carries its prior law, support, structural atoms and defined-draw
#' counts when available. R2 shares are conditional on the term being open.
#'
#' @param fit Fit produced by [JAGS_fit()].
#' @param selection A selection of one BayesTools quantity.
#' @param conditional Whether to condition on its declared inclusion event.
#' @param transform_scaled Whether to request the original predictor scale.
#'   A pure rescaling of an R2 term transforms the draws and prior together.
#'   When an original-scale R2 coefficient mixes dependent fitted terms with
#'   no supported prior law, the method raises `BayesTools_R2_density_unavailable`.
#'   An ordinary scalar fixed coefficient is transformed using the formula
#'   posterior interface. An affected fixed factor cell without a scalar
#'   display route raises `BayesTools_parameter_plot_scale_unavailable`.
#'   Use `FALSE` for the fitted scale. Budget shares, totals and rates are unitless.
#' @return A named `as_mixed_posteriors` list with one scalar posterior.
#' @seealso [parameter_mixed_posterior()], [parameter_relationships()]
#' @export
parameter_plot_samples <- function(fit, selection, conditional = FALSE, transform_scaled = TRUE){

  values <- parameter_mixed_posterior(fit, selection, conditional = conditional,
                                      transform_scaled = transform_scaled)
  q <- selection$quantities
  key <- q$extraction_key[[1L]]
  name <- q$canonical_name
  if(transform_scaled && identical(q$owner_type, "r2_budget")){
    values <- .bt_R2_original_scale_mixed_posterior(fit, selection, values, conditional)
  }else if(transform_scaled && identical(q$role, "fixed_coefficient") &&
           length(attr(fit, "formula_scale", exact = TRUE)[[q$formula_parameter]]) &&
           length(key$dependencies)){
    prior_list <- attr(fit, "prior_list", exact = TRUE)
    owner <- if(identical(key$type, "coordinate") && length(key$dependencies) == 1L) key$dependencies
    if(length(owner) && owner %in% names(prior_list) && !is.prior.factor(prior_list[[owner]]) &&
       !is.prior.vector(prior_list[[owner]])){
      return(as_mixed_posteriors(fit, parameters = owner,
        conditional = if(conditional) owner, transform_scaled = TRUE))
    }
    factor <- .bt_R2_original_scale_factor(fit, q$formula_parameter, key$dependencies)
    if(!identical(factor, 1)){
      stop(errorCondition(paste0("The original-scale plot of '", name,
        "' is unavailable for this fixed coefficient selection. Use 'transform_scaled = FALSE'."),
        class = "BayesTools_parameter_plot_scale_unavailable", call = NULL))
    }
  }
  out <- stats::setNames(list(values), name)
  attr(out, "prior_list") <- stats::setNames(list(prior_none()), name)
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}
