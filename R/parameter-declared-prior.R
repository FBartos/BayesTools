#' @title Declared scalar prior of a random-effect quantity
#'
#' @description Returns the scalar prior that the fitted model assigns to a
#' selected random-effect SD or correlation on the selection's predictor scale.
#' For a treatment-coded or independently coded factor SD this is the marginal
#' of the selected coefficient. It retains the mixture components, their weights,
#' and any spike-and-slab inclusion prior.
#'
#' @param object a `BayesTools_fit`.
#' @param selection one quantity resolved by [parameter_catalog_resolve()].
#' @param source logical. With `FALSE`, the selected quantity must be the
#' identity of a scalar prior. With `TRUE`, a one-to-one transformation such as
#' a variance or a rescaled slope SD returns its source prior in fitted units.
#' For a fitted two-column LKJ correlation this is the Beta prior on its primitive
#' in \code{[0, 1]}. A composite quantity has no single source prior.
#'
#' @return A scalar `prior` object, including a scalar mixture or spike-and-slab
#' prior. Returns `NULL` for a quantity without a supported declaration of this
#' kind. [parameter_prior_density()] gives the law on the selected quantity's
#' own scale when available, including for derived quantities.
#' @export
parameter_declared_prior <- function(object, selection, source = FALSE){

  if(!inherits(object, "BayesTools_fit")){
    stop("'object' must be a BayesTools fit.", call. = FALSE)
  }
  check_bool(source, "source", allow_NA = FALSE)
  object <- .bt_parameter_selection_fit(object, selection)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.", call. = FALSE)
  }
  quantity <- selection$quantities
  key <- quantity$extraction_key[[1L]]
  if(!identical(quantity$provider, "BayesTools") || !identical(key$type, "random_summary")){
    return(NULL)
  }
  if(!source && !(identical(key$source_type, "identity") &&
                  identical(key$source_transform, "identity") &&
                  (identical(key$source_scale, NA_real_) || isTRUE(key$source_scale == 1)))){
    return(NULL)
  }

  .bt_random_scalar_source_prior(object, key)
}

.bt_random_scalar_source_prior <- function(object, key){

  if(!key$source_type %in% c("identity", "one_to_one_transform")){
    return(NULL)
  }
  if(identical(key$source_transform, "lkj2")){
    term <- .bt_parameter_catalog_find_random_term(object, key)
    correlation <- .bt_random_effect_correlation_metadata(term)
    eta <- correlation$eta
    if(!identical(correlation$type, "lkj") || !isTRUE(term$n_columns == 2L) ||
       !identical(key$source_parameter, correlation$primitive_names) ||
       !is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta <= 0){
      return(NULL)
    }
    return(prior("beta", list(alpha = eta, beta = eta)))
  }
  if(!key$evaluator %in% c("sd", "sd_variance", "allocation_sd", "allocation_var", "rho")){
    return(NULL)
  }
  owner <- attr(object, "prior_list", exact = TRUE)[[key$source_prior]]
  if(!is.prior(owner)){
    return(NULL)
  }
  columns <- .prior_linear_prior_columns(key$source_prior, owner)
  index <- match(key$source_parameter, columns)
  if(is.na(index)){
    return(NULL)
  }
  if(!.bt_prior_is_factor_family(owner) && !is.prior.vector(owner)){
    if(is.prior.simple(owner) || inherits(owner, "prior.simple_mixture")) return(owner)
    return(NULL)
  }
  if(is.prior.vector(owner) && is.prior.point(owner)){
    locations <- .posterior_atoms_point_location(owner, length(columns))
    if(is.null(locations)) return(NULL)
    return(prior("point", list(location = locations[[index]]), prior_weights = .prior_model_weight(owner)))
  }
  if(is.prior.ordered(owner) && !is.prior.mixture(owner) && length(columns) == 1L){
    # A two-level cumulative SD with a single increment owns its total.
    spec <- .bt_ordered_spec(key$source_prior, owner)
    increment <- .bt_ordered_tensor(spec, spec$slice_index[[index]], 1)
    if(!length(increment$records) && identical(as.numeric(increment$coefficients), 1)){
      return(owner$total)
    }
  }
  .bt_random_factor_sd_scalar_prior(owner)
}

# The compiler binds a scalar SD prior to each treatment or independent
# coefficient, with one shared component indicator for the whole term. Its
# coordinate marginal keeps that indicator's law and the component metadata.
.bt_random_factor_sd_scalar_prior <- function(x){

  if(!is.prior.treatment(x) && !is.prior.independent(x)) return(NULL)
  if(is.prior.mixture(x)){
    members <- lapply(x, .bt_random_factor_sd_scalar_prior)
    if(any(vapply(members, is.null, logical(1)))) return(NULL)
    x[] <- members
    # Restore the constructors' method dispatch order. The formula compiler
    # can add its factor wrapper after the existing scalar mixture classes.
    class(x) <- c("prior",
      if(is.prior.spike_and_slab(x)) c("prior.spike_and_slab", "prior.simple_spike_and_slab"),
      "prior.simple_mixture", "prior.mixture")
  }else if(!is.prior.simple(x) || is.prior.vector(x)){
    return(NULL)
  }else{
    class(x) <- setdiff(class(x), c("prior.factor", "prior.treatment", "prior.independent"))
  }
  for(field in c(.bt_factor_metadata_names, "K")){
    attr(x, field) <- NULL
  }
  x
}
