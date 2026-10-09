# One bound container component: product scalar coordinates or its joint law,
# followed by the exact declared contrast image. Standalone rng() is unchanged.
.bt_rng_bound_factor_component <- function(component, n, container,
                                            transform_factor_samples){

  design <- attr(container, "factor_design", exact = TRUE)
  component_design <- attr(component, "factor_design", exact = TRUE)
  if(!is.matrix(design) || !is.numeric(design) || any(!is.finite(design)) ||
     ncol(design) < 1L || !identical(design, component_design) ||
     !identical(attr(container, "factor_contrasts", exact = TRUE),
                attr(component, "factor_contrasts", exact = TRUE))){
    stop("Bound factor components must share their declared coding and dimensions.", call. = FALSE)
  }
  coefficient_dim <- ncol(design)
  declared_K <- component$parameters[["K"]]
  container_K <- attr(container, "K", exact = TRUE)
  for(value in list(declared_K, container_K)){
    if(!is.null(value) && (length(value) != 1L ||
       !(is.numeric(value) || (is.logical(value) && is.na(value))) ||
       (!is.na(value) && (!is.finite(value) || value != coefficient_dim)))){
      stop("Bound factor dimensions disagree with an explicit 'K' declaration.", call. = FALSE)
    }
  }
  if(is.prior.ordered(component)){
    stop("This bound component sampler is unavailable for ordered factor priors.", call. = FALSE)
  }
  if(is.prior.vector(component)){
    component$parameters[["K"]] <- coefficient_dim
    raw <- rng(component, n, transform_factor_samples = FALSE)
  }else if(is.prior.simple(component)){
    # Column-major stream: a complete independent scalar stream per coordinate.
    raw <- matrix(.prior_simple_rng(component, n * coefficient_dim),
                  nrow = n, ncol = coefficient_dim)
  }else{
    stop("The bound factor component has no supported coefficient sampler.", call. = FALSE)
  }
  raw <- matrix(raw, nrow = n, ncol = coefficient_dim)
  if(transform_factor_samples) raw %*% t(design) else raw
}
