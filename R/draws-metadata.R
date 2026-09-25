# ============================================================================ #
# draws-metadata.R
# ============================================================================ #
#
# One metadata container for posterior and prior draws: mixed and marginal
# posteriors, their elements, and extracted parameter draws carry their
# supports, atoms, mixture components, undefined draws, prior densities and
# their context, precomputed posterior densities and ordinates, formula
# flags, linear weights, and conditioning in the single validated attribute
# 'bayestools_meta'. Every read and write goes through .bt_meta_get() and
# .bt_meta_set() (exact field access, validated values); posterior_metadata()
# exposes the fields downstream packages attach or read. Free attributes of
# earlier development versions are not read.
#
# ============================================================================ #

.bt_meta_attribute <- "bayestools_meta"

# Free attributes that held draw metadata before the container; they must not
# be used on draws (see the unit lint test).
.bt_meta_legacy_names <- c(
  "posterior_support", "posterior_atoms", "posterior_components",
  "models_ind", "sample_ind", "ordered_total_indicator",
  "ordered_total_component", "undefined_draws",
  "prior_density", "prior_density_context", "prior_densities",
  "posterior_density", "posterior_densities", "posterior_ordinate",
  "posterior_ordinates", "formula_parameter", "formula_log_intercept",
  "formula_scale", "transform_scaled", "linear_weights", "linear_offset",
  "joint_prior_transformation", "conditional", "conditional_rule",
  "condition_key", "condition_event", "resolved_condition_event",
  "effective_conditional", "effective_conditional_rule"
)

# Elements of the 'condition' field (the conditioning of the draws).
.bt_meta_condition_names <- c(
  "conditional", "conditional_rule", "condition_key", "condition_event",
  "resolved_condition_event", "effective_conditional",
  "effective_conditional_rule"
)

# Validators of the draw-metadata fields: each returns NULL for a valid value
# and the reason otherwise.
.bt_meta_validators <- function(){

  list(
    support = function(value){
      # a single support or a list of them keyed by column; malformed support
      # stops with the constructor message
      supports <- if(.posterior_metadata_is_container(value)) value else list(value)
      for(support in supports){
        if(!is.null(support) && is.null(.posterior_support_from_attribute(support))){
          return("it must be created with 'posterior_support_attribute()'")
        }
      }
      NULL
    },
    atoms = function(value){
      .posterior_atoms_from_attribute(value)
      NULL
    },
    components = function(value){
      if(inherits(value, "BayesTools_posterior_components")) NULL else
        "it must be a posterior component table"
    },
    component = function(value){
      if(.bt_meta_is_index(value)) NULL else
        "it must be a vector of positive integer component indices"
    },
    component_source = function(value){
      if(is.character(value) && length(value) == 1L &&
         value %in% .bt_meta_component_sources){
        return(NULL)
      }
      "it must be 'model', 'mixture', or 'spike_and_slab'"
    },
    draw_index = function(value){
      if(.bt_meta_is_index(value)) NULL else
        "it must be a vector of positive integer draw indices"
    },
    ordered_total_component = function(value){
      if(.bt_meta_is_index(value, allow_NA = TRUE)) NULL else
        "it must hold positive integer component indices (NA for models without a total spike)"
    },
    undefined_draws = function(value){
      if(is.character(value) && !anyNA(value)) NULL else
        "it must be a character vector naming the undefined quantities"
    },
    prior_density = function(value){
      if(inherits(value, "prior_density") || is.prior(value)) NULL else
        "it must be a prior density or a prior distribution"
    },
    prior_context = function(value){
      if(.bt_meta_is_prior_context(value)) NULL else
        "it must be a prior-density context"
    },
    prior_densities = function(value){
      if(is.list(value) && (!is.object(value) || inherits(value, "prior_density_list")) &&
         all(vapply(value, function(density){
           is.null(density) || inherits(density, "prior_density") || is.prior(density)
         }, logical(1)))){
        return(NULL)
      }
      "it must be a list of prior densities"
    },
    posterior_density = function(value){
      .posterior_density_kind(value)
      NULL
    },
    posterior_densities = function(value){
      if(!is.list(value) || is.object(value)){
        return("it must be a list of posterior density metadata")
      }
      lapply(value, .posterior_density_kind)
      NULL
    },
    posterior_ordinate = function(value){
      .posterior_ordinate_kind(value)
      NULL
    },
    posterior_ordinates = function(value){
      if(!is.list(value) || is.object(value)){
        return("it must be a list of posterior ordinate metadata")
      }
      lapply(value, .posterior_ordinate_kind)
      NULL
    },
    formula_parameter = function(value){
      if(is.character(value) && length(value) > 0L && !anyNA(value)) NULL else
        "it must name the formula parameter"
    },
    log_intercept = function(value){
      if(is.logical(value) && length(value) == 1L) NULL else
        "it must be TRUE, FALSE, or NA"
    },
    formula_scale = function(value){
      if(is.list(value)) NULL else "it must be a list of formula scaling information"
    },
    transform_scaled = function(value){
      if(isTRUE(value) || isFALSE(value)) NULL else "it must be TRUE or FALSE"
    },
    condition = function(value){
      if(!is.list(value) || is.null(names(value)) ||
         !all(names(value) %in% .bt_meta_condition_names)){
        return("it must be a named list of condition fields")
      }
      for(name in c("conditional", "effective_conditional")){
        if(!is.null(value[[name]]) && !is.character(value[[name]])){
          return(paste0("'", name, "' must be a character vector"))
        }
      }
      for(name in c("conditional_rule", "effective_conditional_rule")){
        if(!is.null(value[[name]]) &&
           !(is.character(value[[name]]) && length(value[[name]]) == 1L &&
             value[[name]] %in% c("AND", "OR"))){
          return(paste0("'", name, "' must be 'AND' or 'OR'"))
        }
      }
      if(!is.null(value[["condition_key"]]) &&
         !(is.character(value[["condition_key"]]) &&
           length(value[["condition_key"]]) == 1L)){
        return("'condition_key' must be a character string")
      }
      for(name in c("condition_event", "resolved_condition_event")){
        if(!is.null(value[[name]]) && !is.list(value[[name]])){
          return(paste0("'", name, "' must be a condition event"))
        }
      }
      NULL
    },
    linear_weights = function(value){
      if(is.numeric(value)) NULL else "it must be a numeric vector or matrix"
    },
    linear_offset = function(value){
      if(is.numeric(value) && length(value) == 1L) NULL else
        "it must be a number"
    },
    joint_prior_transformation = function(value){
      if(is.character(value) && length(value) == 1L) NULL else
        "it must be a character string"
    }
  )
}

# Sources of the per-draw component index: the model of a model-averaged
# ensemble (mix_posteriors()), or the component of a mixture or spike-and-slab
# prior of a single fit (as_mixed_posteriors()).
.bt_meta_component_sources <- c("model", "mixture", "spike_and_slab")

.bt_meta_is_index <- function(x, allow_NA = FALSE){

  if(is.logical(x) && allow_NA && all(is.na(x))){
    return(TRUE)
  }
  if(!is.numeric(x) || (!allow_NA && anyNA(x))){
    return(FALSE)
  }
  x <- x[!is.na(x)]
  all(x >= 1 & x == round(x))
}

.bt_meta_is_prior_context <- function(x){

  inherits(x, "prior_density_context") ||
    inherits(x, "prior_density_model_mixture_context") ||
    inherits(x, "prior_density_conditional_context")
}

.bt_meta_field_names <- c(
  "support", "atoms", "components", "component", "component_source",
  "draw_index", "ordered_total_component", "undefined_draws", "prior_density",
  "prior_context", "prior_densities", "posterior_density",
  "posterior_densities", "posterior_ordinate", "posterior_ordinates",
  "formula_parameter", "log_intercept", "formula_scale", "transform_scaled",
  "condition", "linear_weights", "linear_offset", "joint_prior_transformation"
)

.bt_meta_fields <- function(){
  .bt_meta_field_names
}

.bt_meta_check_field <- function(field){

  if(!is.character(field) || length(field) != 1L || is.na(field) ||
     !field %in% .bt_meta_fields()){
    stop(
      "'", paste(field, collapse = ", "), "' is not a draw-metadata field.",
      call. = FALSE
    )
  }
  invisible(field)
}

.bt_meta_validate <- function(field, value){

  reason <- .bt_meta_validators()[[field]](value)
  if(!is.null(reason)){
    stop("Draw metadata '", field, "' is invalid: ", reason, ".", call. = FALSE)
  }
  invisible(value)
}

# The draw-metadata container of 'x' (NULL when 'x' carries none).
.bt_meta_container <- function(x){

  meta <- attr(x, .bt_meta_attribute, exact = TRUE)
  if(is.null(meta)){
    return(NULL)
  }
  if(!inherits(meta, "BayesTools_draw_metadata") || !is.list(meta) ||
     is.null(names(meta)) ||
     !all(names(meta) %in% c(.bt_meta_fields(), .bt_meta_fingerprint_field))){
    stop("The draw-metadata container is invalid.", call. = FALSE)
  }
  meta
}

# The value of a draw-metadata field, or NULL when it is absent. The metadata
# of draws must describe their current values (.bt_meta_check_current()).
.bt_meta_get <- function(x, field){

  .bt_meta_check_field(field)
  .bt_meta_current_container(x)[[field]]
}

# The container of 'x' (NULL when 'x' carries none), checked to describe the
# current values of draws.
.bt_meta_current_container <- function(x){

  meta <- .bt_meta_container(x)
  if(!is.null(meta) && .bt_meta_is_draws(x)){
    .bt_meta_check_current(meta, .bt_meta_fingerprint(x))
  }
  meta
}

# 'x' with the draw-metadata field set to 'value'; NULL removes the field.
# Setting a field of draws whose values changed after their metadata were
# attached stops: producers that change the values of draws rebuild the
# container first (.bt_meta_refresh()).
.bt_meta_set <- function(x, field, value){

  .bt_meta_assign(x, stats::setNames(list(value), field))
}

# 'x' with the fields of the named list 'fields' set (NULL removes a field),
# with one check and one fingerprint of the draws.
.bt_meta_assign <- function(x, fields){

  for(field in names(fields)){
    .bt_meta_check_field(field)
    if(!is.null(fields[[field]])){
      .bt_meta_validate(field, fields[[field]])
    }
  }
  meta <- .bt_meta_container(x)
  if(is.null(meta) && all(vapply(fields, is.null, logical(1)))){
    return(x)
  }
  fingerprint <- if(.bt_meta_is_draws(x)) .bt_meta_fingerprint(x)
  if(is.null(meta)){
    meta <- structure(list(), names = character(), class = c("BayesTools_draw_metadata", "list"))
  }else if(!is.null(fingerprint)){
    .bt_meta_check_current(meta, fingerprint)
  }
  for(field in names(fields)){
    meta[[field]] <- fields[[field]]
  }
  .bt_meta_write(x, meta, fingerprint)
}

# Draw metadata go stale when the values of the draws change after the
# metadata were attached (e.g., by 'x[] <- ', 'x[i] <- ', or pmin(), which
# keep the attributes of 'x'). The container of draws therefore stores a
# fingerprint of the values it describes: their number, the number of missing
# values, and the sums of the observed values and of the observed values
# weighted by their positions. The fingerprint is computed in one native pass
# (src/r-draw-fingerprint.c).
.bt_meta_fingerprint_field <- "fingerprint"

.bt_meta_is_draws <- function(x){

  is.atomic(x) && typeof(x) %in% c("double", "integer", "logical")
}

# The fingerprint of the values of draws 'x' ('value', stored in the
# container) with the scales of its rounding-error bound ('scale').
.bt_meta_fingerprint <- function(x){

  out <- .Call("BayesTools_draw_fingerprint", x, PACKAGE = "BayesTools")
  list(
    value = c(length = out[[1L]], missing = out[[2L]], sum = out[[3L]], weighted_sum = out[[4L]]),
    scale = c(sum = out[[5L]], weighted_sum = out[[6L]])
  )
}

# Whether a stored fingerprint describes the values of the current one. The
# number of values and of missing values must be equal. The native pass sums
# in double precision in a fixed order, four interleaved partial sums of at
# most m = ceiling(n / 4) summands combined pairwise, so a computed sum is
# within gamma_k * sum(|t|) of the exact one, with k = m + 2 (the weighted
# sum's products are rounded too), gamma_k = k u / (1 - k u), u = 2^-53, and
# t the summands. Two computations of the same values - on the same build
# bit-identical, on other platforms possibly with fused multiply-adds - thus
# agree within 2 gamma_k (1 + gamma_k) times the computed sum(|t|), the last
# factor bounding the rounding of that scale. Differences within this bound
# (for example changes of single draws far below their magnitude) are not
# detected.
.bt_meta_fingerprint_matches <- function(stored, current){

  if(!is.numeric(stored) || !identical(names(stored), names(current$value)) ||
     !identical(stored[c("length", "missing")], current$value[c("length", "missing")])){
    return(FALSE)
  }
  if(identical(stored, current$value)){
    return(TRUE)
  }
  k <- ceiling(current$value[["length"]] / 4) + 2
  gamma <- k * .Machine$double.eps / 2
  gamma <- gamma / (1 - gamma)
  all(vapply(c("sum", "weighted_sum"), function(name){
    stored_sum  <- stored[[name]]
    current_sum <- current$value[[name]]
    if(!is.finite(stored_sum) || !is.finite(current_sum)){
      return(identical(stored_sum, current_sum))
    }
    abs(stored_sum - current_sum) <= 2 * gamma * (1 + gamma) * current$scale[[name]]
  }, logical(1)))
}

.bt_meta_check_current <- function(meta, current){

  if(!.bt_meta_fingerprint_matches(meta[[.bt_meta_fingerprint_field]], current)){
    stop(errorCondition(
      paste0(
        "The metadata of these posterior draws are unavailable: their values ",
        "changed after the metadata were attached (for example by 'x[] <- ', ",
        "'x[i] <- ', or 'pmin()'), so the supports, atoms, and prior densities ",
        "no longer describe them. Use 'marginal_posterior(transformation = )' ",
        "for transformed posterior distributions."
      ),
      class = c("BayesTools_stale_metadata", "BayesTools_metadata"),
      call = NULL
    ))
  }
  invisible(TRUE)
}

# 'x' with the container 'meta' (and, for draws, the fingerprint of their
# values); a container without fields is removed.
.bt_meta_write <- function(x, meta, fingerprint = NULL){

  meta[[.bt_meta_fingerprint_field]] <- NULL
  if(length(meta) == 0L){
    attr(x, .bt_meta_attribute) <- NULL
    return(x)
  }
  if(!is.null(fingerprint)){
    meta[[.bt_meta_fingerprint_field]] <- fingerprint$value
  }
  attr(x, .bt_meta_attribute) <- meta
  x
}

# Rebuilds the container of draws whose values a producer changed together
# with the metadata that describe them (the producer transforms, replaces, or
# removes the affected fields).
.bt_meta_refresh <- function(x){

  meta <- .bt_meta_container(x)
  if(is.null(meta)){
    return(x)
  }
  .bt_meta_write(x, meta, if(.bt_meta_is_draws(x)) .bt_meta_fingerprint(x))
}

# 'x' with several draw-metadata fields set, given as named arguments (one
# check and one fingerprint of the draws).
.bt_meta_update <- function(x, ...){

  fields <- list(...)
  if(length(fields) == 0L){
    return(x)
  }
  if(is.null(names(fields)) || any(!nzchar(names(fields)))){
    stop("Draw-metadata fields must be named.", call. = FALSE)
  }
  .bt_meta_assign(x, fields)
}

# One element of the 'condition' field.
.bt_meta_condition <- function(x, name){

  if(!name %in% .bt_meta_condition_names){
    stop("'", name, "' is not a draw-conditioning field.", call. = FALSE)
  }
  .bt_meta_get(x, "condition")[[name]]
}

# Fields of the public accessor posterior_metadata().
.bt_meta_public_fields <- c(
  "support", "atoms", "undefined_draws", "prior_density", "prior_densities",
  "prior_context", "posterior_density", "posterior_densities",
  "posterior_ordinate", "posterior_ordinates", "condition", "linear_weights"
)

#' @title Metadata of BayesTools posterior draws
#'
#' @description Reads and sets the metadata that BayesTools attaches to
#' posterior draws: the elements and lists returned by [mix_posteriors()],
#' [as_mixed_posteriors()], [marginal_posterior()], and
#' [random_effects_summary_posterior()], and the draws of [parameter_draws()].
#' All metadata are stored in one attribute and validated when they are set;
#' use this accessor instead of setting attributes.
#'
#' @param x posterior draws (an element or a list of them).
#' @param field metadata field:
#' \describe{
#'   \item{\code{"support"}}{exact support, created with
#'   [posterior_support_attribute()] (or a list of such objects keyed by
#'   column).}
#'   \item{\code{"atoms"}}{declared posterior point masses, created with
#'   [posterior_atom_attribute()]; point masses are read only from this
#'   field.}
#'   \item{\code{"undefined_draws"}}{a character vector naming quantities
#'   whose draws may be \code{NA} because the quantity is undefined in those
#'   draws (e.g. \code{"correlation"}), as returned by [parameter_draws()].}
#'   \item{\code{"prior_density"}}{the prior density of the quantity.}
#'   \item{\code{"prior_densities"}}{a list of prior densities keyed by
#'   parameter, attached to the list returned by [as_mixed_posteriors()]
#'   with \code{transform_scaled = TRUE}.}
#'   \item{\code{"prior_context"}}{the joint prior-density context of the
#'   draws.}
#'   \item{\code{"posterior_density"}, \code{"posterior_densities"}}{
#'   precomputed posterior densities created with
#'   [posterior_density_attribute()].}
#'   \item{\code{"posterior_ordinate"}, \code{"posterior_ordinates"}}{
#'   precomputed posterior ordinates created with
#'   [posterior_ordinate_attribute()] or [posterior_ordinate_append()].}
#'   \item{\code{"condition"}}{the conditioning of the draws: a list with
#'   \code{conditional}, \code{conditional_rule}, \code{condition_key},
#'   \code{condition_event}, \code{resolved_condition_event}, and, for
#'   levels whose conditioning was resolved per level,
#'   \code{effective_conditional} and \code{effective_conditional_rule}.}
#'   \item{\code{"linear_weights"}}{the weights of the fitted coordinates
#'   that form a level of a [marginal_posterior()] (a named numeric vector,
#'   or a matrix with one row per draw).}
#' }
#' @param value the new value of the field; \code{NULL} removes it.
#'
#' @details Arithmetic and mathematical functions of posterior draws
#' (\code{Ops} and \code{Math} group generics), \code{c()},
#' \code{as.numeric()}, and subsetting return plain numeric draws without
#' metadata, because supports, atoms, and prior densities do not follow the
#' transformed values. Functions that need the metadata (e.g.
#' [Savage_Dickey_BF()], [plot_posterior()], [marginal_posterior()]) stop on
#' such draws; [marginal_posterior()] with \code{transformation} transforms
#' the draws together with their metadata.
#'
#' Operations that replace values but keep the attributes of the draws
#' (\code{x[] <- }, \code{x[i] <- }, \code{pmin()}, \code{pmax()}) leave
#' metadata that no longer describe the draws. The metadata record a
#' fingerprint of the values they describe (their number, the number of
#' missing values, and two sums of the observed values), and reading or
#' setting the metadata of such draws, here or in any function that uses
#' them, stops with an error of class \code{BayesTools_stale_metadata}
#' (parent class \code{BayesTools_metadata}). Subsetting a list of mixed
#' posteriors with \code{[} keeps the list's metadata, with the prior
#' densities of the omitted parameters removed.
#'
#' The fingerprint cannot detect every change: the two sums are compared
#' within their rounding-error bound, so changes that keep both sums (e.g.,
#' adding d, -2d, and d to three consecutive draws), changes of single draws
#' smaller than about \eqn{5.6 \cdot 10^{-17} n^2} times the draws' mean
#' absolute value for \eqn{n} draws (about \eqn{6 \cdot 10^{-5}} for
#' \eqn{10^6} draws), and reorderings that move draws by only a few positions
#' keep the metadata.
#'
#' @return \code{posterior_metadata()} returns the value of the field or
#' \code{NULL}; the replacement form returns \code{x} with the field set.
#'
#' @examples
#' draws <- stats::rnorm(100)
#' posterior_metadata(draws, "atoms") <- posterior_atom_attribute()
#' posterior_metadata(draws, "atoms")
#'
#' @export
posterior_metadata <- function(x, field){

  .bt_meta_check_public_field(field)
  .bt_meta_get(x, field)
}

#' @rdname posterior_metadata
#' @export
`posterior_metadata<-` <- function(x, field, value){

  .bt_meta_check_public_field(field)
  .bt_meta_set(x, field, value)
}

.bt_meta_check_public_field <- function(field){

  check_char(field, "field", allow_values = .bt_meta_public_fields)
  invisible(field)
}

# Plain numeric draws: 'x' without its class and metadata (dimensions and
# names are kept).
.bt_draws_plain <- function(x){

  keep <- intersect(names(attributes(x)), c("dim", "dimnames", "names"))
  attributes(x) <- attributes(x)[keep]
  x
}

# Arithmetic and mathematical functions of posterior draws return plain
# numeric draws: supports, atoms, prior densities and the other draw metadata
# do not follow the transformed values. marginal_posterior(transformation = )
# transforms the draws together with their metadata.
.bt_draws_ops <- function(e1, e2){

  if(missing(e2)){
    return(get(.Generic)(.bt_draws_plain(e1)))
  }
  get(.Generic)(.bt_draws_plain(e1), .bt_draws_plain(e2))
}

.bt_draws_math <- function(x, ...){

  get(.Generic)(.bt_draws_plain(x), ...)
}

# One method object per group generic for every draw class, so that
# operations combining draws of different classes dispatch to it (R uses the
# internal method, which keeps attributes, when the two methods differ).
#' @exportS3Method Ops mixed_posteriors
Ops.mixed_posteriors <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior
Ops.marginal_posterior <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior.simple
Ops.marginal_posterior.simple <- .bt_draws_ops
#' @exportS3Method Ops marginal_posterior.factor
Ops.marginal_posterior.factor <- .bt_draws_ops
#' @exportS3Method Math mixed_posteriors
Math.mixed_posteriors <- .bt_draws_math
#' @exportS3Method Math marginal_posterior
Math.marginal_posterior <- .bt_draws_math
#' @exportS3Method Math marginal_posterior.simple
Math.marginal_posterior.simple <- .bt_draws_math
#' @exportS3Method Math marginal_posterior.factor
Math.marginal_posterior.factor <- .bt_draws_math

# Subsetting a list of mixed posteriors keeps the list and its metadata: the
# per-element prior densities are subset to the kept elements, and the joint
# fields (prior context, scaling, conditioning, posterior density sources)
# are kept. Subsetting draws returns plain numeric draws.
#' @exportS3Method "[" mixed_posteriors
`[.mixed_posteriors` <- function(x, ...){

  out <- NextMethod()
  if(!is.list(x)){
    return(out)
  }
  kept <- names(out)
  x_attributes <- attributes(x)
  x_attributes <- x_attributes[setdiff(names(x_attributes), c("names", .bt_meta_attribute))]
  x_attributes[["names"]] <- kept
  attributes(out) <- x_attributes

  meta <- .bt_meta_container(x)
  if(is.null(meta)){
    return(out)
  }
  prior_densities <- meta[["prior_densities"]]
  if(!is.null(prior_densities)){
    keep <- vapply(names(prior_densities), function(key){
      any(key == kept | startsWith(key, paste0(kept, "[")))
    }, logical(1))
    subset <- unclass(prior_densities)[keep]
    density_attributes <- attributes(prior_densities)
    density_attributes[["names"]] <- names(subset)
    attributes(subset) <- density_attributes
    meta[["prior_densities"]] <- if(any(keep)) subset
  }
  .bt_meta_write(out, meta)
}

# 'fun' applied to the values of draws 'x', keeping every attribute of 'x'
# (for producers that transform the metadata explicitly).
.bt_draws_transform_values <- function(x, fun){

  meta <- .bt_meta_container(x)
  if(!is.null(meta) && .bt_meta_is_draws(x)){
    .bt_meta_check_current(meta, .bt_meta_fingerprint(x))
  }
  out <- fun(.bt_draws_plain(x))
  attributes(out) <- attributes(x)
  .bt_meta_refresh(out)
}

# Error for plain numeric draws passed where BayesTools posterior draws with
# their metadata are required.
.bt_draws_stop_plain <- function(what){

  stop(
    what, ": arithmetic, mathematical functions, and subsetting of posterior ",
    "draws return plain numeric draws without their supports, atoms, and ",
    "prior densities. Use marginal_posterior(transformation = ) for ",
    "transformed posterior distributions.",
    call. = FALSE
  )
}

# Component indices (into the declared component list of a mixture or
# spike-and-slab prior) of fitted indicator draws: the 'dcat' index of a
# mixture, and for a spike-and-slab prior the 'dbern' inclusion indicator
# (1 = slab) mapped to the positions of its slab and spike components.
.bt_component_from_indicator <- function(prior, indicator){

  indicator <- as.numeric(indicator)
  if(anyNA(indicator)){
    stop("Mixture component indicator draws must not be missing.", call. = FALSE)
  }
  if(is.prior.spike_and_slab(prior)){
    if(any(!indicator %in% c(0, 1))){
      stop("Spike-and-slab indicator draws must be 0 or 1.", call. = FALSE)
    }
    components <- attr(prior, "components", exact = TRUE)
    slab  <- which(components == "alternative")
    spike <- which(components == "null")
    return(ifelse(indicator == 1, slab, spike))
  }
  if(is.prior.mixture(prior)){
    if(any(!indicator %in% seq_along(prior))){
      stop("Mixture indicator draws must index the mixture components.", call. = FALSE)
    }
    return(as.integer(indicator))
  }

  stop("Component indicators require a mixture or spike-and-slab prior.", call. = FALSE)
}

# Whether 'component' indexes the spike (the 'null' component) of a
# spike-and-slab prior.
.bt_component_is_spike <- function(prior, component){

  identical(attr(prior, "components", exact = TRUE)[component], "null")
}

# 'x' with its per-draw component index and the source of that index.
.bt_draws_set_component <- function(x, component, source){

  .bt_meta_update(x, component = as.integer(component), component_source = source)
}

# The per-draw component index of draws. By explicit rule, draws without a
# mixture (parameters of a single fit without mixture priors) form one
# component.
.bt_draws_component <- function(x, n = NROW(x)){

  component <- .bt_meta_get(x, "component")
  if(is.null(component)){
    return(rep(1L, n))
  }
  component
}

# The model of each draw of a model-averaged ensemble (mix_posteriors()), or
# NULL for draws of a single fit.
.bt_draws_model_component <- function(x){

  if(!identical(.bt_meta_get(x, "component_source"), "model")){
    return(NULL)
  }
  .bt_meta_get(x, "component")
}
