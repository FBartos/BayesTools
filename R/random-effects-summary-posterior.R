#' Extract random-effect summary posterior distributions
#'
#' @description Extracts posterior draws for derived random-effect summary
#' quantities with their canonical prior densities. For allocation summaries,
#' the helper keeps raw Dirichlet allocation weights internal and exposes
#' interpretable scalar summaries such as mean-variance SD-component variance
#' multipliers. Gated total-variance proportions are normalized over active
#' components and omit draws on which every component is excluded because a
#' variance share is undefined there.
#'
#' Each summary is the [parameter_mixed_posterior()] of its catalog quantity:
#' its prior density is [parameter_prior_density()] (for gated allocations the
#' mixed measure with its inclusion-gate atoms), and its posterior atoms are
#' declared from its structure: none when the prior density has no point
#' mass; the inclusion-gate atoms of gated total-variance proportions (at 0
#' and 1) and allocation totals (at 0), and the point components of a mixture
#' scale prior whose indicator is monitored, with masses from the posterior
#' inclusion- and mixture-indicator draws. Point masses without such a source
#' leave the atom status undeclared, and posterior plots of such summaries stop
#' instead of inferring point masses from the draws.
#'
#' @param fit model fit created by [JAGS_fit] or a derived-draw view from
#' [JAGS_with_draws()]. Summaries use its supplied draws and local analysis metadata.
#' @param summary semantic quantity to extract: `"var_mult"`, `"var_prop"`,
#'   `"sd_mult"`, `"sd_total"`, `"var_total"`, `"sd_common"`, or
#'   `"var_common"`. The default `"var_mult"` exists only for mean-variance
#'   SD-component allocations: for a fit with total-variance allocations it
#'   stops (request `"var_prop"` or a total instead).
#' @param allocation optional allocation filter. It matches the stable
#'   allocation `name` given to `random_variance_allocation()` (the
#'   identifier the backend uses), not the public `display_name`.
#' @param component optional allocation component label filter, matched
#'   against the component labels shown in the quantity names.
#' @param formula_parameter optional formula parameter filter.
#' @param simplify_names whether returned quantities should use centrally
#'   generated simplified display labels. Defaults to `FALSE`.
#' @param n_prior_points grid size (`n_grid` of [parameter_prior_density()]) of
#'   the attached prior densities.
#'
#' @return A named list of posterior vectors with class `mixed_posteriors` and
#'   `marginal_posterior`. The list can be passed to [plot_posterior], and the
#'   individual vectors carry their `prior_density`, `support`, and declared
#'   `atoms` as draw metadata (see [posterior_metadata()]).
#'
#' @export
random_effects_summary_posterior <- function(
    fit,
    summary = c(
      "var_mult", "var_prop", "sd_mult",
      "sd_total", "var_total", "sd_common", "var_common"
    ),
    allocation = NULL,
    component = NULL,
    formula_parameter = NULL,
    simplify_names = FALSE,
    n_prior_points = 4096){

  if(!.bt_is_jags_analysis_fit(fit) || !inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a BayesTools JAGS fit.", call. = FALSE)
  }
  summary <- match.arg(summary)
  check_char(allocation, "allocation", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(component, "component", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
  check_int(n_prior_points, "n_prior_points", lower = 16)

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  check_list(prior_list, "prior_list")
  if(!all(vapply(prior_list, is.prior, logical(1)))){
    stop("'prior_list' must be a list of priors.", call. = FALSE)
  }

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  keep <- startsWith(quantities$role, "random_") &
    !quantities$internal & quantities$status != "unavailable" &
    quantities$quantity == summary
  keys <- quantities$extraction_key
  if(!is.null(allocation)){
    keep <- keep & vapply(keys, function(key){
      !is.null(key$allocation_label) && key$allocation_label %in% allocation
    }, logical(1))
  }
  if(!is.null(component)){
    keep <- keep & quantities$component %in% component
  }
  if(!is.null(formula_parameter)){
    keep <- keep & quantities$formula_parameter %in% formula_parameter
  }
  quantities <- quantities[keep, , drop = FALSE]
  if(nrow(quantities) == 0L && !any(startsWith(
    catalog$quantities$role,
    "random_"
  ))){
    stop("No random-effect summaries are available in 'fit'.", call. = FALSE)
  }
  if(nrow(quantities) == 0L){
    .bt_random_effect_summary_posterior_no_match_stop(
      summary = summary,
      allocation = allocation,
      component = component,
      formula_parameter = formula_parameter
    )
  }

  display_names <- if(simplify_names){
    quantities$display_label
  }else{
    quantities$canonical_name
  }

  out <- vector("list", nrow(quantities))
  names(out) <- display_names
  out_prior_list <- vector("list", nrow(quantities))
  names(out_prior_list) <- display_names

  for(i in seq_len(nrow(quantities))){
    quantity  <- quantities[i, , drop = FALSE]
    key       <- quantity$extraction_key[[1L]]
    selection <- list(
      schema_version        = .bt_parameter_selection_version,
      parameter_map_version = catalog$schema_version,
      quantity_id           = quantity$quantity_id,
      quantities            = quantity
    )
    class(selection) <- c("BayesTools_parameter_selection", "list")
    # the catalog quantity's mixed posterior: defined draws, catalog support,
    # canonical prior density, and atoms declared from the gate states
    values <- .bt_parameter_mixed_posterior(
      fit            = fit,
      selection      = selection,
      n_grid         = n_prior_points,
      simplify_label = simplify_names
    )
    attr(values, "parameter") <- display_names[i]
    attr(values, "summary_name") <- key$summary_name
    attr(values, "random_summary") <- summary
    attr(values, "random_summary_label") <- display_names[i]
    attr(values, "random_allocation") <- key$allocation_label
    attr(values, "random_component") <- quantity$component

    out[[i]] <- values
    out_prior_list[[i]] <- attr(values, "prior_list", exact = TRUE)
  }

  attr(out, "prior_list") <- out_prior_list
  attr(out, "summary_names") <- vapply(
    quantities$extraction_key,
    `[[`,
    character(1),
    "summary_name"
  )
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}

.bt_random_effect_summary_posterior_no_match_stop <- function(summary,
                                                              allocation,
                                                              component,
                                                              formula_parameter){

  quoted <- function(x) paste0("'", x, "'", collapse = ", ")
  filters <- c(
    if(!is.null(allocation)) paste0("allocation = ", quoted(allocation)),
    if(!is.null(component)) paste0("component = ", quoted(component)),
    if(!is.null(formula_parameter)) paste0("formula_parameter = ", quoted(formula_parameter))
  )
  filter_text <- if(length(filters) > 0L){
    paste0(" matched the requested filters (", paste(filters, collapse = "; "), ")")
  }else{
    " are available"
  }

  detail <- switch(
    summary,
    "var_mult" = paste0(
      "Variance-multiplier summaries are created only for ",
      "random_variance_allocation(..., target = \"sd_component\", ",
      "scale = \"mean_variance\"). Total-variance allocations are returned ",
      "by summary = \"var_prop\"."
    ),
    "var_prop" = paste0(
      "Variance-proportion summaries are created for true total-variance ",
      "allocations. Mean-variance SD-component allocations are returned by ",
      "summary = \"var_mult\"."
    ),
    "sd_mult" = paste0(
      "SD-multiplier summaries are available only for SD-component ",
      "variance allocations."
    ),
    "Requested random-effect summaries are not available."
  )

  stop(
    "No random-effect ",
    summary,
    " summaries",
    filter_text,
    ". ",
    detail,
    call. = FALSE
  )
}
