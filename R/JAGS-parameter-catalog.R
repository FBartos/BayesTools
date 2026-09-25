# Metadata-only semantic parameter catalog and deferred draw extraction.

.bt_parameter_selection_version <- 1L

.bt_parameter_catalog_quantity_columns <- c(
  "quantity_id", "canonical_name", "provider", "namespace", "role",
  "formula_parameter", "owner_type", "owner_name", "quantity", "scale_role",
  "parent_quantity_id", "arguments",
  "term", "component", "display_label",
  "fitted_scale", "display_scale", "status", "fixed_value", "internal",
  "source_type", "support", "definedness", "extraction_key"
)
.bt_parameter_catalog_alias_columns <- c(
  "alias", "quantity_id", "namespace", "component", "simplified"
)

#' Semantic parameter catalogs and deferred draw extraction
#'
#' @description
#' `parameter_catalog()` returns the semantic quantities and aliases from the
#' fitted object's versioned [parameter_map()]. The catalog contains only
#' metadata: selectable quantities, exact aliases, and serializable extraction
#' keys. Square brackets after a factor term always hold a level label, never a
#' coordinate position. Fixed factor terms, including ordinary factor priors
#' whose levels are labelled `1, ..., K`, expose every level or interaction cell
#' as a quantity named `<parameter>[<level>]`, so `term[level]`,
#' `<parameter>[level]`, and `term` plus `component = "level"` all select the
#' level with that label. Treatment and independent cells and the first ordered
#' coordinate are structurally the fitted coordinate; reference cells are
#' structural zeroes; and the remaining cells are reconstructed from the
#' persisted term-only design matrix. A coordinate that is not a level cell (a
#' mean-difference or orthonormal coefficient, or a later ordered increment) is
#' coefficient `j` of the contrast coding, named `<parameter>{j}` and displayed
#' as `term{j}`. The level names of transformed summaries,
#' `<parameter>[dif: level]`, and the row labels those summaries display,
#' `(formula) term[dif: level]` and `term[dif: level]`, are aliases of the
#' level quantity. JAGS
#' coordinate names such as `mu_g[1]` remain backend columns, used by
#' [JAGS_materialize_draws()], and do not select factor quantities. Ordinary
#' level labels remain unchanged; syntax-sensitive characters (including curly
#' braces) are percent-escaped and ambiguous interaction tokens are quoted so
#' that every component remains hypothesis-safe and injective. These are
#' coefficient-level quantities, distinct from estimated marginal means based
#' on full predictions.
#' [parameter_coordinates()] is the linked concrete posterior-coordinate view;
#' the catalog is the semantic view of the same fitted map. Random-effect
#' canonical names follow `(formula) owner: quantity(arguments)`, with `owner: `
#' omitted for a bare or unnamed one-entry random formula, retained for an
#' explicitly named one-entry list, and required when multiple blocks need
#' disambiguation. Parentheses contain parameter or
#' coefficient names, while square brackets inside an argument contain its
#' factor or index level. Examples include `(mu) study: sd(intercept)`,
#' `(mu) study: cor(group[sensitivity],group[specificity])`, and
#' `(mu) var_prop(study)`. Formula-prefix omission is accepted as
#' an alias. Additional aliases are limited to genuine semantic equivalences,
#' such as a CS/HCS pairwise correlation referring to its shared `cor`.
#' Random-effect quantities also carry centrally generated simplified aliases.
#' These remove a sole intercept argument and may omit a redundant owner, but
#' are considered by resolvers only when `simplify_names = TRUE`.
#'
#' Identity random-effect summaries reuse their sampled or structural
#' coordinate. Summaries requiring a scale or covariance transformation expose
#' only the transformed quantity; their fitted-scale inputs and all other
#' private implementation coordinates remain dependencies rather than public
#' catalog rows. Each random extraction key records whether its source mapping
#' is an identity, one-to-one transform, or composite, together with the source
#' coordinate, prior, and transform when those are defined.
#' Constructing or resolving the catalog never accesses posterior draws.
#'
#' `parameter_catalog_extend()` adds plain-data quantities and aliases owned by
#' another provider; the `"BayesTools"` provider name is reserved for native
#' quantities. `parameter_catalog_resolve()` applies optional namespace
#' and component filters and returns a versioned selection only when the match
#' is unique.
#'
#' `parameter_draws()` is the deferred extraction boundary. The BayesTools fit
#' method reads only the coordinates declared by the selected
#' extraction key. Downstream packages can provide methods for package-owned
#' derived quantities. For gated total-variance allocations, realized
#' `sd_total` and `var_total` draws include the all-off zero branch.
#' `var_prop(...)` draws are normalized over active components and are `NA` on
#' draws where the realized allocation total is zero. Original-scale
#' random-effect correlations (`cor(...)` of LKJ blocks) are `NA` on draws
#' where the correlation is undefined, i.e. where one of its SDs is zero.
#' Each catalog quantity declares its `definedness` (`"always"`, or the reason
#' of possibly undefined draws: `"correlation"`, or `"allocation_active"` for
#' the variance shares of gated total-variance allocations) and its exact
#' `support` from the prior provenance of its source coordinates (`NULL` when
#' it is not derivable). The `mcmc.list` returned by `parameter_draws()`
#' carries these as draw metadata ([posterior_metadata()]): `support`, keyed by
#' canonical name, and `undefined_draws`, a character vector named by the
#' canonical names of the possibly undefined quantities with their reasons.
#' Summaries such as [ensemble_estimates_table()] accept missing draws only for
#' columns carrying this declaration: callers that extract a numeric vector
#' keep it by copying the element with
#' `posterior_metadata(x, "undefined_draws") <- `.
#'
#' `parameter_prior_density()` constructs a deterministic
#' `prior_linear_density` for supported map-defined quantities: fitted
#' coordinates and factor levels (the weighted sum of their fitted
#' coordinates under the fitted prior, whose `multiply_by` scales enter only
#' the linear predictor; structural levels are point masses at their fixed
#' value), one-to-one transformations, Dirichlet allocation marginals, and
#' allocation-derived component SDs. Variance proportions are conditional on
#' positive allocation variance. A shared parent inclusion gate cancels from
#' that conditional law, leaving the Dirichlet marginal
#' `Beta(α_i, α_•−α_i)`. Independently gated components keep the mixed measure
#' over the realized active set: atoms at 0 and 1, and a Beta mixture over
#' nonempty sets of other active components. It returns
#' `NULL` when the fitted map does
#' not declare a supported deterministic prior composition; for fitted
#' coordinates and factor levels, only when no fitted prior owns their
#' coordinates (an owning prior-list entry that is not a BayesTools prior
#' stops). The density records
#' the source prior and transform it was built from, so its heights and region
#' probabilities are evaluated on that prior's structural route (bounded-logit
#' correlations on the refined numerical grid). Allocation-derived component
#' SDs are products of density grids without such provenance: they are
#' returned for plotting, and their heights and region probabilities are
#' unavailable. Allocation totals (`sd_total`, `var_total`) and component SDs
#' and variances of gated allocations carry the exact gate atom at zero (no
#' active component, or the component's gate off); the continuous part of a
#' total is the scale prior itself when every component is active and a
#' product of density grids otherwise. Totals of nested allocations and
#' inclusion indicators return `NULL`.
#'
#' `parameter_transform()` returns the one-to-one map from the selected source
#' coordinate to its public semantic quantity when that map exists. Composite
#' quantities return `NULL`. The forward, inverse, and Jacobian helpers are the
#' authoritative evaluators for these serializable descriptors, so downstream
#' consumers do not need to reproduce backend correlation, allocation, or
#' formula-scale transformations.
#'
#' @param object fitted object.
#' @param catalog a `BayesTools_parameter_catalog` object.
#' @param quantities quantity rows matching the quantity schema.
#' @param aliases alias rows matching the alias schema.
#' @param provider scalar provider name owning every added quantity.
#' @param alias scalar exact canonical name or alias.
#' @param namespace optional exact namespace filter.
#' @param component optional exact component filter.
#' @param simplify_names whether to accept centrally generated simplified
#'   random-effect aliases. Defaults to `FALSE`.
#' @param selection a `BayesTools_parameter_selection` returned by
#'   `parameter_catalog_resolve()`.
#' @param model_samples optional numeric matrix containing the declared source
#'   coordinates. This lets downstream summaries evaluate a selected semantic
#'   quantity on an already materialized posterior sample.
#' @param transform a serializable transform descriptor returned by
#'   `parameter_transform()`.
#' @param n_grid number of grid points used for deterministic induced prior
#'   densities.
#' @param tail_prob probability omitted from each continuous source-prior tail
#'   when constructing a finite numerical grid.
#' @param values numeric values on the source scale for
#'   `parameter_transform_forward()` and `parameter_transform_jacobian()`, or
#'   on the public quantity scale for `parameter_transform_inverse()`.
#' @param ... arguments for methods.
#'
#' @return `parameter_catalog()` and `parameter_catalog_extend()` return a
#' `BayesTools_parameter_catalog`. `parameter_catalog_schema()` returns schema
#' descriptions. `parameter_catalog_resolve()` returns a
#' `BayesTools_parameter_selection`. `parameter_draws()` returns a
#' `coda::mcmc.list` for BayesTools-owned quantities.
#' `parameter_prior_density()` returns a `prior_linear_density` or `NULL`.
#' `parameter_transform()` returns a serializable transform descriptor or
#' `NULL`; the transform helpers return numeric values.
#'
#' @export parameter_catalog
#' @export parameter_catalog_schema
#' @export parameter_catalog_extend
#' @export parameter_catalog_resolve
#' @export parameter_draws
#' @export parameter_prior_density
#' @export parameter_transform
#' @export parameter_transform_forward
#' @export parameter_transform_inverse
#' @export parameter_transform_jacobian
#' @name parameter_catalog
NULL

#' @rdname parameter_catalog
parameter_catalog <- function(object, ...){

  UseMethod("parameter_catalog")
}

#' @rdname parameter_catalog
#' @exportS3Method parameter_catalog BayesTools_fit
parameter_catalog.BayesTools_fit <- function(object, ...){

  JAGS_validate_fit_contract(object, requires = "parameter_map")
  .bt_parameter_map_catalog(parameter_map(object))
}

#' @rdname parameter_catalog
parameter_catalog_schema <- function(){

  quantities <- data.frame(
    field = .bt_parameter_catalog_quantity_columns,
    type = c(
      rep("character", 11L), "list", rep("character", 6L),
      "numeric", "logical", "character", "list", "character", "list"
    ),
    description = c(
      "Stable provider-namespaced quantity identifier.",
      "Exact canonical semantic name.",
      "Package or subsystem owning extraction.",
      "Exact resolver namespace.",
      "Semantic role.",
      "Owning formula output parameter, or an empty string.",
      "Semantic owner type, such as random_block or variance_allocation, or an empty string.",
      "Stable semantic owner name, or an empty string.",
      "Semantic quantity such as sd, cor, var_prop, var_mult, or sd_mult, or an empty string.",
      "Allocation scale role: total, common, or an empty string.",
      "Owning aggregate quantity identifier for a nested quantity, or an empty string.",
      "Ordered semantic parameter arguments; factor/index levels use square brackets inside each argument.",
      "Formula term or semantic subterm, or an empty string.",
      "Quantity component, or an empty string.",
      "Unambiguous default display label.",
      "Scale of stored dependency coordinates.",
      "Scale of the selected quantity.",
      "One of sampled, structural, derived, or unavailable.",
      "Exact structural value; otherwise NA.",
      "Whether the quantity is private implementation metadata.",
      "Relationship to concrete source coordinates: identity, one_to_one_transform, composite, or none.",
      "Exact support of the quantity from the prior provenance of its source coordinates (a posterior_support_attribute() object), or NULL when it is not derivable.",
      "Draws where the quantity is defined: 'always', or the reason of possibly undefined (NA) draws ('correlation': an SD is zero; 'allocation_active': no allocation component is active).",
      "Serializable plain-data extraction recipe."
    ),
    stringsAsFactors = FALSE
  )
  aliases <- data.frame(
    field = .bt_parameter_catalog_alias_columns,
    type = c(rep("character", 4L), "logical"),
    description = c(
      "Exact accepted alias.",
      "Quantity identifier targeted by the alias.",
      "Exact resolver namespace.",
      "Optional component filter, or an empty string.",
      "Whether the alias requires 'simplify_names = TRUE'."
    ),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = .bt_parameter_map_version,
    quantities = quantities,
    aliases = aliases
  )
}

#' @rdname parameter_catalog
parameter_catalog_extend <- function(catalog, quantities, aliases,
                                     provider){

  .bt_validate_parameter_catalog(catalog)
  check_char(provider, "provider", check_length = 1L, allow_NA = FALSE)
  if(!grepl("^[A-Za-z][A-Za-z0-9.]*$", provider)){
    stop("'provider' must start with a letter and contain only letters, numbers, and dots.",
         call. = FALSE)
  }
  if(identical(provider, "BayesTools")){
    stop("'BayesTools' is reserved for quantities created by BayesTools.",
         call. = FALSE)
  }
  .bt_validate_parameter_catalog_tables(
    quantities,
    aliases,
    known_quantity_ids = c(
      catalog$quantities$quantity_id,
      quantities$quantity_id
    )
  )
  if(nrow(quantities) == 0L && nrow(aliases) == 0L){
    stop("At least one extension quantity or alias must be supplied.",
         call. = FALSE)
  }
  if(any(quantities$provider != provider)){
    stop("Every extension quantity must be owned by 'provider'.", call. = FALSE)
  }
  provider_prefix <- paste0(provider, "::")
  if(any(!startsWith(quantities$quantity_id, provider_prefix))){
    stop("Every extension 'quantity_id' must start with the provider namespace '",
         provider_prefix, "'.", call. = FALSE)
  }
  collisions <- intersect(
    catalog$quantities$quantity_id,
    quantities$quantity_id
  )
  if(length(collisions) > 0L){
    stop(
      "Extension quantity IDs already exist: ",
      paste0("'", collisions, "'", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out <- catalog
  out$quantities <- rbind(out$quantities, quantities)
  out$aliases <- rbind(out$aliases, aliases)
  rownames(out$quantities) <- NULL
  rownames(out$aliases) <- NULL
  .bt_validate_parameter_catalog(out)
  out
}

#' @rdname parameter_catalog
parameter_catalog_resolve <- function(catalog, alias, namespace = NULL,
                                      component = NULL,
                                      simplify_names = FALSE){

  .bt_validate_parameter_catalog(catalog)
  check_char(alias, "alias", check_length = 1L, allow_NA = FALSE)
  check_char(namespace, "namespace", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(component, "component", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)

  quantities <- catalog$quantities
  public <- !quantities$internal
  canonical_rows <- public & quantities$canonical_name == alias
  if(!is.null(namespace)){
    canonical_rows <- canonical_rows & quantities$namespace == namespace
  }
  if(!is.null(component)){
    canonical_rows <- canonical_rows & quantities$component == component
  }
  canonical_ids <- quantities$quantity_id[canonical_rows]
  alias_rows <- catalog$aliases$alias == alias &
    (!catalog$aliases$simplified | simplify_names)
  if(!is.null(namespace)){
    alias_rows <- alias_rows & catalog$aliases$namespace == namespace
  }
  if(!is.null(component)){
    alias_rows <- alias_rows & catalog$aliases$component == component
  }
  candidate_ids <- unique(c(
    canonical_ids,
    catalog$aliases$quantity_id[alias_rows]
  ))
  candidates <- quantities[
    quantities$quantity_id %in% candidate_ids & public,
    ,
    drop = FALSE
  ]

  if(nrow(candidates) == 0L){
    available <- sort(unique(c(
      quantities$canonical_name[public],
      catalog$aliases$alias[
        catalog$aliases$quantity_id %in% quantities$quantity_id[public] &
          (!catalog$aliases$simplified | simplify_names)
      ]
    )))
    .bt_parameter_catalog_stop(
      class = "BayesTools_parameter_not_found",
      message = paste0("No public parameter quantity matches '", alias, "'."),
      alias = alias,
      available = available
    )
  }
  if(nrow(candidates) > 1L){
    candidates <- candidates[order(candidates$quantity_id), , drop = FALSE]
    .bt_parameter_catalog_stop(
      class = "BayesTools_parameter_ambiguous",
      message = paste0(
        "Parameter alias '", alias, "' is ambiguous; use 'namespace' or ",
        "'component' to select one quantity."
      ),
      alias = alias,
      candidates = candidates
    )
  }

  out <- list(
    schema_version = .bt_parameter_selection_version,
    parameter_map_version = catalog$schema_version,
    quantity_id = candidates$quantity_id,
    quantities = candidates
  )
  class(out) <- c("BayesTools_parameter_selection", "list")
  .bt_validate_parameter_selection(out)
  out
}

#' @rdname parameter_catalog
parameter_draws <- function(object, selection, ...){

  UseMethod("parameter_draws")
}

#' @rdname parameter_catalog
parameter_prior_density <- function(object, selection, ...){

  UseMethod("parameter_prior_density")
}

#' @rdname parameter_catalog
#' @exportS3Method parameter_draws BayesTools_fit
parameter_draws.BayesTools_fit <- function(object, selection,
                                            model_samples = NULL, ...){

  catalog <- parameter_catalog(object)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  quantities <- selection$quantities

  .bt_parameter_draws_from_quantities(
    object        = object,
    quantities    = quantities,
    model_samples = model_samples
  )
}

# Draws of catalog quantities with the catalog's draw metadata: the exact
# supports of the quantities, keyed by canonical name, and the quantities
# whose draws may be undefined (NA) with the reason of their definedness.
.bt_parameter_draws_from_quantities <- function(
    object, quantities, model_samples = NULL){

  draws <- .bt_parameter_draws_values(object, quantities, model_samples)
  supports <- stats::setNames(quantities$support, quantities$canonical_name)
  supports <- supports[!vapply(supports, is.null, logical(1))]
  if(length(supports) > 0L){
    draws <- .bt_meta_set(draws, "support", unclass(supports))
  }
  undefined <- quantities$definedness != "always"
  if(any(undefined)){
    draws <- .bt_meta_set(draws, "undefined_draws", stats::setNames(
      quantities$definedness[undefined],
      quantities$canonical_name[undefined]
    ))
  }
  draws
}

.bt_parameter_draws_values <- function(object, quantities, model_samples = NULL){

  if(any(quantities$provider != "BayesTools")){
    stop(
      "The selection contains quantities owned by another provider; use that provider's 'parameter_draws()' method.",
      call. = FALSE
    )
  }
  unavailable <- quantities$status == "unavailable"
  if(any(unavailable)){
    stop(
      "The parameter quantity '", quantities$canonical_name[unavailable][[1L]],
      "' is unavailable in this fit: its source coordinates are not part of ",
      "the posterior draws. Refit the model with those coordinates monitored.",
      call. = FALSE
    )
  }

  ordinary <- vapply(
    quantities$extraction_key,
    function(key) identical(key$type, "coordinate"),
    logical(1)
  )
  if(all(ordinary)){
    if(!is.null(model_samples)){
      dependencies <- vapply(
        quantities$extraction_key,
        function(key) key$dependencies[[1L]],
        character(1)
      )
      out <- .bt_parameter_draw_supplied_dependencies(
        model_samples,
        dependencies
      )
      out <- lapply(out, function(chain){
        values <- as.matrix(chain)
        colnames(values) <- quantities$canonical_name
        coda::mcmc(values)
      })
      return(coda::mcmc.list(out))
    }
    return(JAGS_materialize_draws(
      object,
      parameters = quantities$canonical_name,
      include_internal = FALSE
    ))
  }
  if(nrow(quantities) != 1L){
    stop("Mixed or multiple derived selections are not supported in one extraction call.",
         call. = FALSE)
  }

  key <- quantities$extraction_key[[1L]]
  dependencies <- if(is.null(model_samples)){
    .bt_parameter_draw_dependencies(object, key$dependencies)
  }else{
    .bt_parameter_draw_supplied_dependencies(
      model_samples,
      key$dependencies
    )
  }
  out <- vector("list", length(dependencies))
  for(chain_i in seq_along(dependencies)){
    chain <- dependencies[[chain_i]]
    model_samples <- if(identical(key$type, "factor_level") &&
                        length(key$dependencies) == 0L){
      matrix(
        numeric(),
        nrow = nrow(chain),
        ncol = 0L,
        dimnames = list(NULL, character())
      )
    }else{
      as.matrix(chain)
    }
    values <- if(identical(key$type, "factor_level")){
      .bt_parameter_draw_factor_level(key, model_samples)
    }else{
      .bt_parameter_draw_random_summary(
        fit = object,
        key = key,
        model_samples = model_samples
      )
    }
    values <- matrix(
      as.numeric(values),
      ncol = 1L,
      dimnames = list(NULL, quantities$canonical_name)
    )
    mcpar <- attr(chain, "mcpar", exact = TRUE)
    out[[chain_i]] <- coda::mcmc(
      values,
      start = mcpar[1L],
      end = mcpar[2L],
      thin = mcpar[3L]
    )
  }
  coda::mcmc.list(out)
}

#' @rdname parameter_catalog
#' @exportS3Method parameter_prior_density BayesTools_fit
parameter_prior_density.BayesTools_fit <- function(
    object, selection, n_grid = .prior_linear_density_default_grid(),
    tail_prob = .prior_linear_density_tail_prob(), ...){

  catalog <- parameter_catalog(object)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.",
         call. = FALSE)
  }
  check_int(n_grid, "n_grid", lower = 16)
  check_real(tail_prob, "tail_prob", lower = 0, upper = 0.5,
             allow_bound = FALSE)

  .bt_parameter_prior_density_quantity(
    object    = object,
    selection = selection,
    n_grid    = n_grid,
    tail_prob = tail_prob
  )
}

# The prior density of the one quantity of a validated selection. With
# 'conditional', the density of a gated random-effect allocation quantity is
# conditional on its inclusion event (parameter_mixed_posterior()): the
# component's gates on for a component SD or variance, its own gate on for a
# variance proportion, and any component active for an allocation total.
.bt_parameter_prior_density_quantity <- function(object, selection, n_grid,
                                                 tail_prob,
                                                 conditional = FALSE){

  quantity <- selection$quantities[1L, , drop = FALSE]
  if(!identical(quantity$provider, "BayesTools")){
    stop("The selected quantity is owned by another provider.",
         call. = FALSE)
  }
  key <- quantity$extraction_key[[1L]]
  if(key$type %in% c("coordinate", "factor_level")){
    out <- .bt_parameter_prior_density_coordinates(
      object = object,
      quantity = quantity,
      key = key,
      n_grid = n_grid,
      tail_prob = tail_prob
    )
    if(!is.null(out)){
      attr(out, "parameter_prior_density") <- list(
        quantity_id = quantity$quantity_id,
        source = "fitted_parameter_map"
      )
    }
    return(out)
  }
  if(!identical(key$type, "random_summary")){
    return(NULL)
  }

  out <- if(identical(key$evaluator, "sd") &&
             isTRUE(key$allocation_derived)){
    .bt_parameter_prior_density_random_component_sd(
      object = object,
      key = key,
      n_grid = n_grid,
      tail_prob = tail_prob,
      conditional = conditional
    )
  }else if(identical(key$evaluator, "sd_variance") &&
           isTRUE(key$allocation_derived)){
    .bt_parameter_prior_density_square(
      .bt_parameter_prior_density_random_component_sd(
        object = object,
        key = key,
        n_grid = n_grid,
        tail_prob = tail_prob,
        conditional = conditional
      )
    )
  }else if(identical(key$evaluator, "allocation")){
    .bt_parameter_prior_density_allocation_quantity(
      object = object,
      selection = selection,
      key = key,
      n_grid = n_grid,
      tail_prob = tail_prob,
      conditional = conditional
    )
  }else if(key$evaluator %in% c("allocation_sd", "allocation_var") &&
           identical(key$source_type, "composite")){
    .bt_parameter_prior_density_allocation_total(
      object = object,
      key = key,
      n_grid = n_grid,
      tail_prob = tail_prob,
      conditional = conditional
    )
  }else{
    .bt_parameter_prior_density_direct_quantity(
      object = object,
      selection = selection,
      key = key,
      n_grid = n_grid,
      tail_prob = tail_prob
    )
  }
  if(!is.null(out)){
    attr(out, "parameter_prior_density") <- c(
      list(
        quantity_id = quantity$quantity_id,
        source = "fitted_parameter_map"
      ),
      if(isTRUE(conditional)) list(conditional = "inclusion")
    )
  }
  out
}

# Prior density of a fitted coordinate or a factor level, the weighted sum
# sum_j w_j beta_j of its fitted coordinates: the prior-density context of
# the priors owning those coordinates (a coefficient's 'multiply_by' scales
# only its linear-predictor contribution, not the coefficient) evaluated at
# the level's weights. Structural quantities are their fixed value. NULL only
# when no fitted prior owns the coordinates; an owner that is not a BayesTools
# prior (or a list of them) stops, and failures to build the prior-density
# context propagate.
.bt_parameter_prior_density_coordinates <- function(object, quantity, key,
                                                    n_grid, tail_prob){

  dependencies <- key$dependencies
  if(length(dependencies) == 0L || identical(quantity$status, "structural")){
    fixed <- quantity$fixed_value
    if(!is.numeric(fixed) || length(fixed) != 1L || !is.finite(fixed)){
      return(NULL)
    }
    return(.bt_parameter_prior_density_scalar(
      prior("point", list(location = fixed)),
      n_grid = n_grid,
      tail_prob = tail_prob
    ))
  }
  weights <- if(identical(key$type, "coordinate")){
    rep(1, length(dependencies))
  }else{
    key$weights
  }

  prior_list <- .marginal_posterior_strip_multiply_by(
    attr(object, "prior_list", exact = TRUE)
  )
  if(!is.list(prior_list) || length(prior_list) == 0L){
    return(NULL)
  }
  owner_columns <- lapply(names(prior_list), function(parameter){
    columns <- .prior_linear_prior_columns(parameter, prior_list[[parameter]])
    if(any(dependencies %in% columns)) columns else NULL
  })
  names(owner_columns) <- names(prior_list)
  owners <- names(owner_columns)[!vapply(owner_columns, is.null, logical(1))]
  columns <- unlist(owner_columns[owners], use.names = FALSE)
  if(length(owners) == 0L || !all(dependencies %in% columns)){
    return(NULL)
  }
  unsupported <- owners[!vapply(
    prior_list[owners], .bt_parameter_prior_density_is_prior, logical(1)
  )]
  if(length(unsupported) > 0L){
    stop(
      "The prior density of '", quantity$canonical_name, "' is unavailable: ",
      "the prior distribution of '", unsupported[1L], "' is not a BayesTools prior.",
      call. = FALSE
    )
  }

  context <- .prior_density_build_context(
    prior_list   = prior_list[owners],
    column_names = columns,
    n_grid       = n_grid,
    tail_prob    = tail_prob
  )
  full_weights <- stats::setNames(numeric(length(columns)), columns)
  full_weights[dependencies] <- weights
  .prior_density_from_context(context, full_weights)
}

# Whether a prior-list entry is a BayesTools prior or a model list of them
# (the prior-list entries a prior-density context accepts).
.bt_parameter_prior_density_is_prior <- function(prior){

  if(is.prior(prior)){
    return(TRUE)
  }
  is.list(prior) && !is.object(prior) && length(prior) > 0L &&
    all(vapply(prior, is.prior, logical(1)))
}

.bt_parameter_prior_density_direct_quantity <- function(
    object, selection, key, n_grid, tail_prob){

  if(!key$source_type %in% c("identity", "one_to_one_transform")){
    return(NULL)
  }
  prior_list <- attr(object, "prior_list", exact = TRUE)
  source_prior <- if(nzchar(key$source_prior) &&
                     key$source_prior %in% names(prior_list)){
    prior_list[[key$source_prior]]
  }else{
    NULL
  }
  if(is.null(source_prior) && identical(key$source_transform, "lkj2")){
    random_term <- .bt_parameter_catalog_find_random_term(object, key)
    eta <- random_term$correlation$eta
    if(is.numeric(eta) && length(eta) == 1L && is.finite(eta) && eta > 0){
      source_prior <- prior("beta", list(alpha = eta, beta = eta))
    }
  }
  if(is.null(source_prior) || !is.prior(source_prior) ||
     .prior_linear_prior_dimension(source_prior) != 1L){
    return(NULL)
  }
  transform <- .bt_parameter_transform_from_quantity(
    object,
    selection$quantities[1L, , drop = FALSE]
  )
  .bt_parameter_prior_density_transformed(
    source_prior,
    transform,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
}

.bt_parameter_prior_density_allocation_quantity <- function(
    object, selection, key, n_grid, tail_prob, conditional = FALSE){

  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(object, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(
    object,
    key,
    random_term
  )
  prior_list <- attr(object, "prior_list", exact = TRUE)
  source_prior <- prior_list[[allocation$weight_name]]
  transform <- .bt_parameter_transform_from_quantity(
    object,
    selection$quantities[1L, , drop = FALSE]
  )
  if(is.null(transform)){
    if(!identical(selection$quantities$quantity, "var_prop") ||
       !identical(allocation$scale, "total_variance")){
      return(NULL)
    }
    return(.bt_parameter_prior_density_gated_var_prop(
      allocation = allocation,
      source_prior = source_prior,
      prior_list = prior_list,
      index = key$index,
      n_grid = n_grid,
      tail_prob = tail_prob,
      conditional = conditional
    ))
  }
  index <- key$index
  if(identical(selection$quantities$quantity, "sd_mult") &&
     index > allocation$n_targets){
    index <- index - allocation$n_targets
  }
  beta_prior <- .bt_parameter_prior_density_simplex_marginal(
    source_prior,
    index
  )
  if(is.null(beta_prior)){
    return(NULL)
  }
  .bt_parameter_prior_density_transformed(
    beta_prior,
    transform,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
}

.bt_parameter_prior_density_gated_var_prop_free_gate_limit <- function(){
  20L
}

# Prior probability that a scalar scale prior, or a finite mixture of such
# priors (e.g. a model-averaged scale), is positive.
.bt_parameter_prior_density_positive_probability <- function(prior_object){

  if(is.null(prior_object) || !is.prior(prior_object)){
    return(NA_real_)
  }
  if(is.prior.mixture(prior_object) || is.prior.spike_and_slab(prior_object)){
    weights <- .prior_density_ordinate_mixture_weights(prior_object)
    if(is.null(weights)){
      return(NA_real_)
    }
    components <- vapply(seq_along(weights), function(i){
      .bt_parameter_prior_density_positive_probability(prior_object[[i]])
    }, numeric(1))
    return(sum(weights * components))
  }
  if(!is.prior.simple(prior_object)){
    return(NA_real_)
  }
  probability <- tryCatch(ccdf(prior_object, 0), error = function(e) NA_real_)
  if(!is.numeric(probability) || length(probability) != 1L ||
     !is.finite(probability)){
    return(NA_real_)
  }
  as.numeric(probability)
}

# Prior inclusion probability of the allocation gate 'indicator_name'.
.bt_parameter_prior_density_gate_probability <- function(prior_list,
                                                         indicator_name){

  gate_priors <- Filter(function(prior){
    identical(attr(prior, "random_allocation_indicator", exact = TRUE),
              indicator_name)
  }, prior_list)
  if(length(gate_priors) != 1L){
    return(NA_real_)
  }
  .bt_parameter_prior_density_inclusion_probability(gate_priors[[1L]])
}

.bt_parameter_prior_density_inclusion_probability <- function(prior_object){

  if(is.null(prior_object) || !is.prior(prior_object)){
    return(NA_real_)
  }
  probability <- mean(prior_object)
  if(!is.numeric(probability) || length(probability) != 1L ||
     !is.finite(probability) || probability < 0 || probability > 1){
    return(NA_real_)
  }
  probability
}

.bt_parameter_prior_density_allocation_parent_gate_probability <- function(
    allocation, prior_list){

  parent_factors <- allocation$parent_factors
  if(!is.list(parent_factors) || length(parent_factors) == 0L){
    return(1)
  }

  probability <- 1
  for(parent in parent_factors){
    if(is.null(parent$inclusion_name) || !nzchar(parent$inclusion_name)){
      next
    }
    gate_priors <- Filter(function(prior){
      identical(attr(prior, "random_allocation_indicator", exact = TRUE),
                parent$inclusion_name)
    }, prior_list)
    if(length(gate_priors) != 1L){
      return(NA_real_)
    }
    parent_probability <- .bt_parameter_prior_density_inclusion_probability(
      gate_priors[[1L]]
    )
    if(!is.finite(parent_probability) || parent_probability <= 0){
      return(parent_probability)
    }
    probability <- probability * parent_probability
  }
  probability
}

.bt_parameter_prior_density_allocation_component_probabilities <- function(
    allocation, K){

  probability <- rep(1, K)
  inclusion <- allocation$inclusion
  if(!is.list(inclusion) || length(inclusion) == 0L){
    return(probability)
  }

  for(record in inclusion){
    index <- record$index
    if(!is.numeric(index) || length(index) != 1L || is.na(index) ||
       index != as.integer(index) || index < 1L || index > K){
      return(NULL)
    }
    component_probability <- .bt_parameter_prior_density_inclusion_probability(
      record$prior
    )
    if(!is.finite(component_probability)){
      return(NULL)
    }
    probability[[as.integer(index)]] <- component_probability
  }
  probability
}

.bt_parameter_prior_density_gated_var_prop <- function(
    allocation, source_prior, prior_list, index, n_grid, tail_prob,
    conditional = FALSE){

  positive_scale <- .bt_parameter_prior_density_positive_probability(
    allocation$source$prior
  )
  if(!is.finite(positive_scale) || positive_scale <= 0){
    return(NULL)
  }

  parent_probability <- .bt_parameter_prior_density_allocation_parent_gate_probability(
    allocation,
    prior_list
  )
  if(!is.numeric(parent_probability) || length(parent_probability) != 1L ||
     !is.finite(parent_probability) || parent_probability <= 0 ||
     parent_probability > 1){
    return(NULL)
  }

  if(!is.prior.simplex(source_prior) ||
     !identical(source_prior$distribution, "dirichlet")){
    return(NULL)
  }
  alpha <- source_prior$parameters$alpha
  if(!is.numeric(alpha) || length(alpha) < 2L ||
     any(!is.finite(alpha)) || any(alpha <= 0) ||
     !is.numeric(index) || length(index) != 1L || is.na(index) ||
     index != as.integer(index) || index < 1L || index > length(alpha)){
    return(NULL)
  }
  index <- as.integer(index)
  K <- length(alpha)
  probability <- .bt_parameter_prior_density_allocation_component_probabilities(
    allocation,
    K
  )
  if(is.null(probability)){
    return(NULL)
  }
  if(isTRUE(conditional)){
    # conditional on the component's own gate: the proportion given that it
    # is active (the atom at 0 of an inactive component drops out)
    probability[[index]] <- 1
  }

  .bt_parameter_prior_density_gated_var_prop_mixture(
    alpha = alpha,
    index = index,
    probability = probability,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
}

.bt_parameter_prior_density_gated_var_prop_mixture <- function(
    alpha, index, probability, n_grid, tail_prob){

  K <- length(alpha)
  p_i <- probability[[index]]
  others <- seq_len(K)[-index]
  p_others <- probability[others]
  alpha_i <- alpha[[index]]
  alpha_others <- alpha[others]

  p_all_off <- prod(1 - probability)
  p_positive <- 1 - p_all_off
  if(!is.finite(p_positive) || p_positive <= 0){
    return(NULL)
  }

  always_off <- p_others <= 0
  always_on <- p_others >= 1
  free <- !always_off & !always_on
  fixed_beta <- sum(alpha_others[always_on])
  free_alpha <- alpha_others[free]
  free_p <- p_others[free]
  n_free <- length(free_p)
  if(n_free > .bt_parameter_prior_density_gated_var_prop_free_gate_limit()){
    stop(
      "Independently gated variance-proportion prior density is unavailable for more than 20 free inclusion gates.",
      call. = FALSE
    )
  }

  p_others_all_off <- prod(1 - p_others)
  p_atom0 <- (1 - p_i) * (1 - p_others_all_off) / p_positive
  p_atom1 <- p_i * p_others_all_off / p_positive

  mix_weight <- 1
  mix_beta <- fixed_beta
  if(n_free > 0L){
    for(j in seq_len(n_free)){
      mix_weight <- c(mix_weight * (1 - free_p[[j]]), mix_weight * free_p[[j]])
      mix_beta <- c(mix_beta, mix_beta + free_alpha[[j]])
    }
  }

  cont_keep <- mix_beta > 0
  cont_weight <- (p_i / p_positive) * mix_weight[cont_keep]
  cont_beta <- mix_beta[cont_keep]
  if(length(cont_beta) > 0L){
    key <- sprintf("%a", cont_beta)
    unique_key <- unique(key)
    grouped_beta <- cont_beta[match(unique_key, key)]
    grouped_weight <- vapply(unique_key, function(k){
      sum(cont_weight[key == k])
    }, numeric(1), USE.NAMES = FALSE)
  }else{
    grouped_beta <- numeric()
    grouped_weight <- numeric()
  }

  # The measure is one mixture prior (atoms at 0 and 1 and Beta components),
  # so its density carries the deterministic provenance of that prior.
  components <- list()
  if(p_atom0 > 0){
    components[[length(components) + 1L]] <- prior(
      "point", list(location = 0), prior_weights = p_atom0
    )
  }
  if(p_atom1 > 0){
    components[[length(components) + 1L]] <- prior(
      "point", list(location = 1), prior_weights = p_atom1
    )
  }
  for(j in seq_along(grouped_beta)){
    if(grouped_weight[[j]] > 0){
      components[[length(components) + 1L]] <- prior(
        "beta",
        list(alpha = alpha_i, beta = grouped_beta[[j]]),
        prior_weights = grouped_weight[[j]]
      )
    }
  }
  if(length(components) == 0L){
    return(NULL)
  }
  measure <- if(length(components) == 1L){
    components[[1L]]
  }else{
    prior_mixture(
      components,
      is_null = vapply(components, is.prior.point, logical(1))
    )
  }
  .bt_parameter_prior_density_scalar(
    measure,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
}

# The allocation chain of an allocation-derived component SD: its
# allocation, the scalar source prior, and the factors (Dirichlet weight,
# index, scale and optional inclusion gate) from the root allocation to the
# component. NULL when the component SD is not a direct allocation product.
.bt_parameter_prior_density_component_chain <- function(object, key){

  random_term <- .bt_parameter_catalog_find_random_term(object, key)
  if(!.bt_parameter_catalog_random_sd_is_direct(
    random_term = random_term,
    parameter = key$formula_parameter,
    formula_scale = attr(object, "formula_scale", exact = TRUE)
  )){
    return(NULL)
  }
  binding <- random_term$sd_binding
  if(is.null(binding) || !isTRUE(binding$true_allocation) ||
     length(binding$allocations) != 1L ||
     !is.numeric(key$index) || length(key$index) != 1L ||
     is.na(key$index)){
    return(NULL)
  }
  allocation <- binding$allocations[[1L]]
  target <- .bt_random_effect_allocation_target_metadata(allocation)
  source <- allocation$source
  source_prior <- source$prior
  if(is.null(source_prior) || !is.prior(source_prior) ||
     .prior_linear_prior_dimension(source_prior) != 1L){
    return(NULL)
  }
  if(identical(target, "sd_component")){
    leaf_index <- allocation$leaf_index_by_column
    if(length(leaf_index) < key$index || is.na(leaf_index[[key$index]])){
      return(NULL)
    }
    factors <- allocation$parent_factors
    factors[[length(factors) + 1L]] <- list(
      weight_name = allocation$weight_name,
      index = leaf_index[[key$index]],
      scale = allocation$scale,
      n_targets = allocation$n_targets
    )
  }else if(identical(target, "block")){
    factors <- allocation$factors
  }else{
    return(NULL)
  }

  list(allocation = allocation, source_prior = source_prior, factors = factors)
}

.bt_parameter_prior_density_random_component_sd <- function(
    object, key, n_grid, tail_prob, conditional = FALSE){

  chain <- .bt_parameter_prior_density_component_chain(object, key)
  if(is.null(chain)){
    return(NULL)
  }
  source_prior <- chain$source_prior
  factors <- chain$factors
  prior_list <- attr(object, "prior_list", exact = TRUE)
  dist <- .bt_parameter_prior_density_scalar(
    source_prior,
    n_grid = n_grid,
    tail_prob = tail_prob
  )
  # the component SD is the gated product sd * prod_f sqrt(scale_f w_f) *
  # gate_f: the gates, independent of the rest, add an atom at zero
  active_probability <- 1
  for(factor in factors){
    if(!is.null(factor$inclusion_name) &&
       is.character(factor$inclusion_name) &&
       length(factor$inclusion_name) == 1L &&
       !is.na(factor$inclusion_name) && nzchar(factor$inclusion_name)){
      gate_probability <- .bt_parameter_prior_density_gate_probability(
        prior_list,
        factor$inclusion_name
      )
      if(!is.finite(gate_probability)){
        return(NULL)
      }
      active_probability <- active_probability * gate_probability
    }
    factor_prior <- prior_list[[factor$weight_name]]
    beta_prior <- .bt_parameter_prior_density_simplex_marginal(
      factor_prior,
      factor$index
    )
    if(is.null(beta_prior)){
      return(NULL)
    }
    scale <- if(identical(factor$scale, "mean_variance")){
      factor$n_targets
    }else if(identical(factor$scale, "total_variance")){
      1
    }else{
      return(NULL)
    }
    factor_dist <- .bt_parameter_prior_density_transformed(
      beta_prior,
      list(type = "sqrt_scale", scale = as.numeric(scale)),
      n_grid = n_grid,
      tail_prob = tail_prob
    )
    dist <- .prior_linear_density_product(
      dist,
      factor_dist,
      n_grid = n_grid
    )
  }
  # conditional on the inclusion event (every gate of the chain on), the
  # gate atom drops out
  .bt_parameter_prior_density_with_zero_atom(
    dist,
    if(isTRUE(conditional)) 0 else 1 - active_probability,
    n_grid,
    tail_prob
  )
}

# 'dist' with an additional atom at zero of probability 'zero' (the rest of
# the mass is scaled by 1 - zero): one mixture prior with provenance when
# 'dist' is the density of a single (mixture) prior, otherwise the mixed
# grid. Returns NULL for an unusable probability.
.bt_parameter_prior_density_with_zero_atom <- function(dist, zero, n_grid,
                                                       tail_prob){

  if(is.null(dist) || !is.numeric(zero) || length(zero) != 1L ||
     !is.finite(zero) || zero < 0 || zero > 1){
    return(NULL)
  }
  if(zero == 0){
    return(dist)
  }
  if(zero == 1){
    return(.bt_parameter_prior_density_scalar(
      prior("point", list(location = 0)), n_grid = n_grid, tail_prob = tail_prob
    ))
  }
  adaptive <- attr(dist, "adaptive_evaluation", exact = TRUE)
  arguments <- adaptive$arguments
  if(identical(adaptive$kind, "linear_combination") &&
     is.null(arguments$output_transformation) &&
     length(arguments$weights) == 1L && isTRUE(unname(arguments$weights) == 1) &&
     is.prior(arguments$prior_list[[names(arguments$weights)]])){
    source <- arguments$prior_list[[names(arguments$weights)]]
    if(is.prior.mixture(source) || is.prior.spike_and_slab(source)){
      source_weights <- .prior_density_ordinate_mixture_weights(source)
      components <- if(is.null(source_weights)) NULL else lapply(seq_along(source), function(i) source[[i]])
    }else if(is.prior.simple(source)){
      source_weights <- 1
      components <- list(source)
    }else{
      components <- NULL
    }
    if(!is.null(components) &&
       all(vapply(components, is.prior.simple, logical(1)))){
      components <- lapply(seq_along(components), function(i){
        component <- components[[i]]
        attr(component, "model_prior_weights") <- NULL
        component[["prior_weights"]] <- (1 - zero) * source_weights[[i]]
        component
      })
      components <- c(
        list(prior("point", list(location = 0), prior_weights = zero)),
        components
      )
      return(.bt_parameter_prior_density_scalar(
        prior_mixture(components, is_null = vapply(components, is.prior.point, logical(1))),
        n_grid = n_grid,
        tail_prob = tail_prob
      ))
    }
  }
  dx <- .prior_linear_density_dx(dist)
  if(!is.finite(dx) || dx <= 0){
    dx <- 1 / max(n_grid - 1L, 1L)
  }
  .prior_linear_density_normalize(
    .prior_linear_density_mix(
      list(.prior_linear_density_point(0), dist),
      c(zero, 1 - zero),
      dx = dx,
      n_grid = n_grid
    ),
    warn = TRUE
  )
}

# The square of a nonnegative density (a variance from an SD): exact through
# 'exp_lin' when the density records a transformation-free linear
# combination without point masses, otherwise its transformed grid (atoms at
# zero stay at zero).
.bt_parameter_prior_density_square <- function(dist){

  if(is.null(dist)){
    return(NULL)
  }
  adaptive <- attr(dist, "adaptive_evaluation", exact = TRUE)
  if(identical(adaptive$kind, "linear_combination") &&
     is.null(adaptive$arguments$output_transformation) &&
     (is.null(dist$points) || nrow(dist$points) == 0L)){
    hull <- .prior_linear_combination_support_hull(
      adaptive$arguments$prior_list,
      adaptive$arguments$weights,
      adaptive$arguments$source_transforms
    )
    if(!is.null(hull) && hull[1L] >= 0){
      arguments <- adaptive$arguments
      return(.prior_linear_combination_density(
        prior_list        = arguments$prior_list,
        weights           = arguments$weights,
        n_grid            = arguments$n_grid,
        tail_prob         = arguments$tail_prob,
        source_transforms = arguments$source_transforms,
        output_transformation = "exp_lin",
        output_transformation_arguments = list(a = 0, b = 2)
      ))
    }
  }
  .prior_linear_density_transform(dist, "exp_lin", list(a = 0, b = 2))
}

# Total SD (or variance) of an allocation with component inclusion gates and
# no nested parent shares: sd * sqrt(T), where T sums the Dirichlet weights of
# the active components (T = 1 without component gates or for other
# allocation scales). With no active component the total is exactly zero.
# The active sets S mix: sd itself when every component is active, and
# sd * sqrt(W_S), W_S ~ Beta(alpha_S, alpha_{-S}), otherwise. Only a total
# whose continuous part is the scale prior itself is exact; the products
# are plotting densities (grids without provenance), and the atom at zero is
# exact in both cases. A total variance is the square of the total SD.
.bt_parameter_prior_density_allocation_total <- function(object, key,
                                                         n_grid, tail_prob,
                                                         conditional = FALSE){

  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(object, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(object, key, random_term)
  source <- allocation$source
  if(!is.list(source) || !identical(source$shape, "scalar") ||
     !is.prior(source$prior) ||
     .prior_linear_prior_dimension(source$prior) != 1L ||
     length(allocation$parent_factors) > 0L){
    return(NULL)
  }
  K <- allocation$n_targets
  probability <- if(identical(allocation$scale, "total_variance")){
    .bt_parameter_prior_density_allocation_component_probabilities(allocation, K)
  }else{
    rep(1, K)
  }
  prior_list <- attr(object, "prior_list", exact = TRUE)
  weight_prior <- prior_list[[allocation$weight_name]]
  if(is.null(probability) || !is.prior.simplex(weight_prior) ||
     !identical(weight_prior$distribution, "dirichlet") ||
     length(weight_prior$parameters$alpha) != K){
    return(NULL)
  }
  alpha <- weight_prior$parameters$alpha
  # the total SD; a total variance is its square (squaring the SD products
  # keeps them on uniform source grids)
  finish <- if(identical(key$evaluator, "allocation_var")){
    .bt_parameter_prior_density_square
  }else{
    identity
  }
  scale_density <- .bt_parameter_prior_density_scalar(
    source$prior,
    n_grid = n_grid,
    tail_prob = tail_prob
  )

  free <- which(probability > 0 & probability < 1)
  if(length(free) > .bt_parameter_prior_density_gated_var_prop_free_gate_limit()){
    stop(
      "Independently gated allocation-total prior density is unavailable for more than 20 free inclusion gates.",
      call. = FALSE
    )
  }
  # active sets: always-on components plus every subset of the free gates
  active <- matrix(probability >= 1, nrow = 1L)
  set_probability <- 1
  for(k in free){
    off <- active
    on <- active
    on[, k] <- TRUE
    active <- rbind(off, on)
    set_probability <- c(set_probability * (1 - probability[[k]]),
                         set_probability * probability[[k]])
  }
  zero <- sum(set_probability[rowSums(active) == 0L])
  if(isTRUE(conditional)){
    # conditional on an active component, the empty-set atom drops out (the
    # set mixture below is normalized over the nonempty sets)
    if(zero >= 1){
      return(NULL)
    }
    zero <- 0
  }
  full <- rowSums(active) == K
  partial <- !full & rowSums(active) > 0L & set_probability > 0

  if(!any(partial)){
    return(finish(.bt_parameter_prior_density_with_zero_atom(
      scale_density, zero, n_grid, tail_prob
    )))
  }

  dists <- list()
  weights <- numeric()
  if(any(full & set_probability > 0)){
    dists[[1L]] <- scale_density
    weights <- sum(set_probability[full])
  }
  for(i in which(partial)){
    share <- prior("beta", list(alpha = sum(alpha[active[i, ]]),
                                beta  = sum(alpha[!active[i, ]])))
    share_density <- .bt_parameter_prior_density_transformed(
      share,
      list(type = "sqrt_scale", scale = 1),
      n_grid = n_grid,
      tail_prob = tail_prob
    )
    dists[[length(dists) + 1L]] <- .prior_linear_density_product(
      scale_density,
      share_density,
      n_grid = n_grid
    )
    weights <- c(weights, set_probability[[i]])
  }
  dx_values <- vapply(dists, .prior_linear_density_dx, numeric(1))
  dx_values <- dx_values[is.finite(dx_values) & dx_values > 0]
  dx <- if(length(dx_values) > 0L) min(dx_values) else 1 / max(n_grid - 1L, 1L)
  continuous <- .prior_linear_density_normalize(
    .prior_linear_density_mix(dists, weights, dx = dx, n_grid = n_grid),
    warn = TRUE
  )
  finish(.bt_parameter_prior_density_with_zero_atom(continuous, zero, n_grid, tail_prob))
}

.bt_parameter_prior_density_simplex_marginal <- function(prior_object,
                                                         index){

  if(!is.prior.simplex(prior_object) ||
     !identical(prior_object$distribution, "dirichlet")){
    return(NULL)
  }
  alpha <- prior_object$parameters$alpha
  if(!is.numeric(alpha) || length(alpha) < 2L ||
     any(!is.finite(alpha)) || any(alpha <= 0) ||
     !is.numeric(index) || length(index) != 1L || is.na(index) ||
     index != as.integer(index) || index < 1L || index > length(alpha)){
    return(NULL)
  }
  index <- as.integer(index)
  prior(
    "beta",
    list(alpha = alpha[[index]], beta = sum(alpha[-index]))
  )
}

.bt_parameter_prior_density_scalar <- function(
    prior_object, n_grid, tail_prob, output_transformation = NULL,
    output_transformation_arguments = NULL){

  .prior_linear_combination_density(
    prior_list = list(source = prior_object),
    weights = c(source = 1),
    n_grid = n_grid,
    tail_prob = tail_prob,
    output_transformation = output_transformation,
    output_transformation_arguments = output_transformation_arguments
  )
}

# Density of the semantic transform of a scalar source prior, built with the
# transform as the output transformation so that the density records its
# deterministic provenance: affine maps are 'lin', 'tanh' is named, and the
# square-root and square maps of a nonnegative source are 'exp_lin'
# (sqrt(s x) = exp(log(s) / 2) x^(1/2), (s x)^2 = exp(2 log(s)) x^2).
# Bounded-logit maps have no named equivalent: their provenance records the
# map itself, and heights use the refined numerical grid.
.bt_parameter_prior_density_transformed <- function(prior_object, transform,
                                                    n_grid, tail_prob){

  if(is.null(transform)){
    return(NULL)
  }
  scalar <- function(output_transformation = NULL, arguments = NULL){
    .bt_parameter_prior_density_scalar(
      prior_object,
      n_grid = n_grid,
      tail_prob = tail_prob,
      output_transformation = output_transformation,
      output_transformation_arguments = arguments
    )
  }
  nonnegative <- function(){
    lower <- prior_object$truncation$lower
    is.numeric(lower) && length(lower) == 1L && !is.na(lower) && lower >= 0
  }
  switch(
    transform$type,
    "identity" = scalar(),
    "affine" = scalar("lin", list(a = transform$offset, b = transform$scale)),
    "tanh" = scalar("tanh"),
    "sqrt_scale" = if(nonnegative()){
      scalar("exp_lin", list(a = log(transform$scale) / 2, b = 1 / 2))
    },
    "square" = if(nonnegative()){
      scalar("exp_lin", list(
        a = 2 * log(.bt_parameter_transform_square_scale(transform)),
        b = 2
      ))
    },
    "bounded_logit" = scalar(.bt_parameter_prior_density_bounded_logit(
      transform$lower, transform$upper
    )),
    NULL
  )
}

# lower + (upper - lower) * plogis(x) as a transformation list whose closures
# capture only the two bounds.
.bt_parameter_prior_density_bounded_logit <- function(lower, upper){

  force(lower)
  width <- upper - lower
  list(
    fun = function(x) lower + width * stats::plogis(x),
    inv = function(x) stats::qlogis((x - lower) / width),
    jac = function(x) width * stats::plogis(x) * (1 - stats::plogis(x))
  )
}

#' @rdname parameter_catalog
parameter_transform <- function(object, selection){

  if(!inherits(object, "BayesTools_fit")){
    stop("'object' must be a BayesTools fit.", call. = FALSE)
  }
  catalog <- parameter_catalog(object)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.",
         call. = FALSE)
  }
  quantity <- selection$quantities[1L, , drop = FALSE]

  .bt_parameter_transform_from_quantity(object, quantity)
}

.bt_parameter_transform_from_quantity <- function(object, quantity){

  if(!identical(quantity$provider, "BayesTools")){
    stop("The selected quantity is owned by another provider.",
         call. = FALSE)
  }
  key <- quantity$extraction_key[[1L]]
  if(!identical(key$type, "random_summary") ||
     identical(key$source_type, "composite")){
    return(NULL)
  }

  source_transform <- key$source_transform
  transform <- if(identical(source_transform, "identity") ||
                  identical(source_transform, "var_prop")){
    list(type = "identity")
  }else if(identical(source_transform, "lkj2")){
    list(type = "affine", offset = -1, scale = 2)
  }else if(identical(source_transform, "fisher_z")){
    list(type = "tanh")
  }else if(identical(source_transform, "logit")){
    random_term <- .bt_parameter_catalog_find_random_term(object, key)
    bounds <- random_term$correlation$bounds
    if(!is.numeric(bounds) || length(bounds) != 2L ||
       any(!is.finite(bounds)) || bounds[1L] >= bounds[2L]){
      return(NULL)
    }
    list(
      type = "bounded_logit",
      lower = unname(bounds[1L]),
      upper = unname(bounds[2L])
    )
  }else if(identical(source_transform, "square")){
    list(type = "square")
  }else if(identical(source_transform, "random_var")){
    formula_scale <- attr(object, "formula_scale", exact = TRUE)
    if(is.null(formula_scale) || length(formula_scale) == 0L){
      list(type = "square")
    }else if(is.numeric(key$source_scale) &&
             length(key$source_scale) == 1L &&
             is.finite(key$source_scale) && key$source_scale > 0){
      list(type = "square", scale = key$source_scale)
    }else{
      NULL
    }
  }else if(identical(source_transform, "random_sd")){
    formula_scale <- attr(object, "formula_scale", exact = TRUE)
    if(is.null(formula_scale) || length(formula_scale) == 0L){
      list(type = "identity")
    }else if(is.numeric(key$source_scale) &&
             length(key$source_scale) == 1L &&
             is.finite(key$source_scale) && key$source_scale > 0){
      list(type = "affine", offset = 0, scale = key$source_scale)
    }else{
      NULL
    }
  }else if(source_transform %in% c("var_mult", "sd_mult")){
    random_term <- if(nzchar(key$random_block)){
      .bt_parameter_catalog_find_random_term(object, key)
    }else{
      NULL
    }
    allocation <- .bt_parameter_catalog_find_allocation(
      object,
      key,
      random_term
    )
    allocation_scale <- .bt_random_effect_allocation_scale_metadata(
      allocation,
      context = "Parameter transform"
    )
    n_targets <- .bt_random_effect_summary_allocation_n_targets(
      allocation,
      K = allocation$n_targets
    )
    variance_scale <- if(identical(allocation_scale, "mean_variance")){
      as.numeric(n_targets)
    }else{
      1
    }
    if(identical(source_transform, "var_mult")){
      list(type = "affine", offset = 0, scale = variance_scale)
    }else{
      list(type = "sqrt_scale", scale = variance_scale)
    }
  }else{
    NULL
  }

  if(!is.null(transform)){
    .bt_validate_parameter_transform(transform)
  }
  transform
}

#' @rdname parameter_catalog
parameter_transform_forward <- function(values, transform){

  .bt_validate_parameter_transform(transform)
  if(identical(transform$type, "identity")){
    return(values)
  }
  if(identical(transform$type, "affine")){
    return(transform$offset + transform$scale * values)
  }
  if(identical(transform$type, "tanh")){
    return(tanh(values))
  }
  if(identical(transform$type, "bounded_logit")){
    return(transform$lower +
      (transform$upper - transform$lower) * stats::plogis(values))
  }
  if(identical(transform$type, "sqrt_scale")){
    scaled <- transform$scale * values
    out <- rep(NaN, length(scaled))
    valid <- !is.na(scaled) & scaled >= 0
    out[valid] <- sqrt(scaled[valid])
    return(out)
  }
  if(identical(transform$type, "square")){
    return((.bt_parameter_transform_square_scale(transform) * values)^2)
  }

  stop("Unsupported semantic parameter transform.", call. = FALSE)
}

#' @rdname parameter_catalog
parameter_transform_inverse <- function(values, transform){

  .bt_validate_parameter_transform(transform)
  if(identical(transform$type, "identity")){
    return(values)
  }
  if(identical(transform$type, "affine")){
    return((values - transform$offset) / transform$scale)
  }
  if(identical(transform$type, "tanh")){
    return(atanh(values))
  }
  if(identical(transform$type, "bounded_logit")){
    probability <- (values - transform$lower) /
      (transform$upper - transform$lower)
    return(stats::qlogis(probability))
  }
  if(identical(transform$type, "sqrt_scale")){
    return(values^2 / transform$scale)
  }
  if(identical(transform$type, "square")){
    return(sqrt(values) / .bt_parameter_transform_square_scale(transform))
  }

  stop("Unsupported semantic parameter transform.", call. = FALSE)
}

#' @rdname parameter_catalog
parameter_transform_jacobian <- function(values, transform){

  .bt_validate_parameter_transform(transform)
  if(identical(transform$type, "identity")){
    return(rep(1, length(values)))
  }
  if(identical(transform$type, "affine")){
    return(rep(abs(transform$scale), length(values)))
  }
  if(identical(transform$type, "tanh")){
    return(1 - tanh(values)^2)
  }
  if(identical(transform$type, "bounded_logit")){
    probability <- stats::plogis(values)
    return((transform$upper - transform$lower) *
      probability * (1 - probability))
  }
  if(identical(transform$type, "sqrt_scale")){
    scaled <- transform$scale * values
    out <- rep(NaN, length(scaled))
    valid <- !is.na(scaled) & scaled >= 0
    out[valid] <- transform$scale / (2 * sqrt(scaled[valid]))
    return(out)
  }
  if(identical(transform$type, "square")){
    return(2 * .bt_parameter_transform_square_scale(transform)^2 * abs(values))
  }

  stop("Unsupported semantic parameter transform.", call. = FALSE)
}

.bt_validate_parameter_transform <- function(transform){

  scalar_number <- function(value, finite = TRUE){
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      (!finite || is.finite(value))
  }
  valid <- is.list(transform) &&
    is.character(transform$type) && length(transform$type) == 1L &&
    !is.na(transform$type) && transform$type %in% c(
      "identity", "affine", "tanh", "bounded_logit", "sqrt_scale",
      "square"
    )
  if(valid && identical(transform$type, "affine")){
    valid <- scalar_number(transform$offset) &&
      scalar_number(transform$scale) && transform$scale != 0
  }
  if(valid && identical(transform$type, "bounded_logit")){
    valid <- scalar_number(transform$lower) &&
      scalar_number(transform$upper) && transform$lower < transform$upper
  }
  if(valid && identical(transform$type, "sqrt_scale")){
    valid <- scalar_number(transform$scale) && transform$scale > 0
  }
  if(valid && identical(transform$type, "square") && !is.null(transform$scale)){
    valid <- scalar_number(transform$scale) && transform$scale > 0
  }
  if(!valid){
    stop("Unsupported semantic parameter transform.", call. = FALSE)
  }

  invisible(TRUE)
}

# The optional 'scale' of a square transform: (scale * x)^2.
.bt_parameter_transform_square_scale <- function(transform){

  if(is.null(transform$scale)) 1 else transform$scale
}

.bt_parameter_catalog_empty_quantities <- function(){

  out <- data.frame(
    quantity_id = character(),
    canonical_name = character(),
    provider = character(),
    namespace = character(),
    role = character(),
    formula_parameter = character(),
    owner_type = character(),
    owner_name = character(),
    quantity = character(),
    scale_role = character(),
    parent_quantity_id = character(),
    term = character(),
    component = character(),
    display_label = character(),
    fitted_scale = character(),
    display_scale = character(),
    status = character(),
    fixed_value = numeric(),
    internal = logical(),
    source_type = character(),
    definedness = character(),
    stringsAsFactors = FALSE
  )
  out$arguments <- I(list())
  out$support <- I(list())
  out <- out[, setdiff(.bt_parameter_catalog_quantity_columns, "extraction_key"),
             drop = FALSE]
  out$extraction_key <- I(list())
  out
}

.bt_parameter_catalog_empty_aliases <- function(){

  data.frame(
    alias = character(),
    quantity_id = character(),
    namespace = character(),
    component = character(),
    simplified = logical(),
    stringsAsFactors = FALSE
  )
}

.bt_parameter_catalog_new <- function(quantities, aliases){

  out <- list(
    schema_version = .bt_parameter_map_version,
    quantities = quantities,
    aliases = aliases
  )
  class(out) <- c("BayesTools_parameter_catalog", "list")
  .bt_validate_parameter_catalog(out)
  out
}

.bt_parameter_catalog_quantity_id <- function(canonical_name, namespace,
                                               role){

  encoded <- JAGS_parameter_encode(list(
    kind = "catalog",
    formula_parameter = namespace,
    term = canonical_name,
    role = role
  ))
  paste0("BayesTools::", encoded)
}

.bt_parameter_catalog_quantity <- function(
    canonical_name, namespace, role, formula_parameter = "", term = "",
    component = "", display_label = canonical_name,
    fitted_scale = "fitted_original", display_scale = fitted_scale,
    status = "derived", fixed_value = NA_real_, internal = FALSE,
    owner_type = "", owner_name = "", quantity = "",
    scale_role = "", parent_quantity_id = "",
    arguments = character(), source_type = "none", extraction_key){

  out <- .bt_parameter_catalog_empty_quantities()
  scalar_columns <- setdiff(names(out), c("arguments", "support", "extraction_key"))
  out[1L, scalar_columns] <- list(
    .bt_parameter_catalog_quantity_id(canonical_name, namespace, role),
    canonical_name,
    "BayesTools",
    namespace,
    role,
    formula_parameter,
    owner_type,
    owner_name,
    quantity,
    scale_role,
    parent_quantity_id,
    term,
    component,
    display_label,
    fitted_scale,
    display_scale,
    status,
    fixed_value,
    internal,
    source_type,
    "always"
  )
  out$arguments <- I(list(as.character(arguments)))
  # declared by .bt_parameter_catalog_add_support() when the catalog is built
  out$support <- I(list(NULL))
  out$extraction_key <- I(list(extraction_key))
  out
}

.bt_parameter_catalog_factor_component_token <- function(x){

  replacements <- c(
    "%" = "%25",
    "[" = "%5B",
    "]" = "%5D",
    "`" = "%60",
    "\\" = "%5C",
    "\"" = "%22",
    "{" = "%7B",
    "}" = "%7D",
    "\r" = "%0D",
    "\n" = "%0A",
    "\t" = "%09",
    "\f" = "%0C",
    "\b" = "%08",
    "\a" = "%07",
    "\v" = "%0B"
  )
  for(token in names(replacements)){
    x <- gsub(token, replacements[[token]], x, fixed = TRUE)
  }
  reserved <- c(",", "=")
  needs_quotes <- !nzchar(x) ||
    grepl("^[[:space:]]|[[:space:]]$", x) ||
    any(vapply(reserved, function(token){
    grepl(token, x, fixed = TRUE)
  }, logical(1)))
  if(needs_quotes){
    encodeString(x, quote = "\"")
  }else{
    x
  }
}

.bt_parameter_catalog_factor_cell_names <- function(design_info){

  level_names <- design_info$level_names
  if(is.null(level_names)){
    return(design_info$cell_names)
  }
  if(length(level_names) == 1L){
    cell_names <- vapply(
      design_info$cell_names,
      .bt_parameter_catalog_factor_component_token,
      character(1)
    )
    if(anyDuplicated(cell_names)){
      stop(
        "Parameter catalog factor metadata do not identify factor levels uniquely.",
        call. = FALSE
      )
    }
    return(cell_names)
  }
  level_grid <- .factor_cell_grid(level_names)
  if(nrow(level_grid) != length(design_info$cell_names) ||
     is.null(names(level_grid)) || any(!nzchar(names(level_grid)))){
    stop(
      "Parameter catalog factor metadata do not identify every interaction cell.",
      call. = FALSE
    )
  }
  cell_names <- vapply(seq_len(nrow(level_grid)), function(cell){
    terms <- vapply(seq_along(level_grid), function(term){
      paste0(
        .bt_parameter_catalog_factor_component_token(names(level_grid)[term]),
        "=",
        .bt_parameter_catalog_factor_component_token(
          as.character(level_grid[[term]][cell])
        )
      )
    }, character(1))
    paste0(terms, collapse = ", ")
  }, character(1))
  if(anyDuplicated(cell_names)){
    stop(
      "Parameter catalog factor metadata do not identify interaction cells uniquely.",
      call. = FALSE
    )
  }
  cell_names
}

# Square brackets after a factor term always hold a level (or cell) label. A
# fitted coordinate that is not itself a level cell - a mean-difference or
# orthonormal coefficient, or an ordered increment - is coefficient `j` of the
# contrast coding and is written with curly braces, `term{j}`.
.bt_parameter_catalog_factor_label <- function(formula_parameter, term,
                                               component){

  paste0(
    .bt_random_effect_summary_formula_prefix(formula_parameter, TRUE),
    term,
    if(startsWith(component, "{")) component else paste0("[", component, "]")
  )
}

.bt_parameter_catalog_factor_coefficient_component <- function(index){

  paste0("{", index, "}")
}

# Contrasts whose term design rows are exact 0/1 structure. For these, a design
# row equal to a unit vector states that the level cell is that coordinate.
# Mean-difference and orthonormal rows come from an eigendecomposition, where a
# unit row is a floating-point coincidence, never a structural identity.
.bt_parameter_catalog_structural_contrasts <- function(){

  c(
    "contr.treatment",
    "contr.independent",
    "contr.ordered_cumulative",
    "contr.ordered_cumulative_levels"
  )
}

.bt_parameter_catalog_factor_prior <- function(parameter, prior){

  if(!.bt_prior_is_factor_family(prior) ||
     .bt_is_random_effect_prior(prior)){
    return(NULL)
  }
  .complete_factor_metadata(prior, parameter)
}

# For each fitted coordinate (design column) of a factor term, the design row
# (level cell) that the coordinate structurally is, or NA.
.bt_factor_direct_cells <- function(prior, design){

  .bt_factor_direct_cells_design(
    design,
    attr(prior, "factor_contrasts", exact = TRUE)
  )
}

.bt_factor_direct_cells_design <- function(design, contrasts){

  direct <- rep.int(NA_integer_, ncol(design))
  if(length(contrasts) == 0L ||
     !all(as.character(contrasts) %in%
            .bt_parameter_catalog_structural_contrasts())){
    return(direct)
  }
  # A unit row: exactly one nonzero entry, equal to one.
  nonzero <- design != 0
  unit_rows <- which(rowSums(nonzero) == 1L)
  if(length(unit_rows) == 0L){
    return(direct)
  }
  unit_columns <- max.col(nonzero[unit_rows, , drop = FALSE], ties.method = "first")
  unit <- design[cbind(unit_rows, unit_columns)] == 1
  unit_rows <- unit_rows[unit]
  unit_columns <- unit_columns[unit]
  for(coordinate in seq_len(ncol(design))){
    matches <- unit_rows[unit_columns == coordinate]
    if(length(matches) == 1L){
      direct[coordinate] <- matches
    }
  }
  direct
}

# Display names of a fixed factor prior's fitted coordinates, in coordinate
# order: a coordinate that is structurally one level cell is named by that
# cell's level labels, as summary tables name level cells (`mu_g[10]`,
# `mu_f[b]__xXx__g[v]`); any other coordinate is coefficient `j` of the
# contrast coding, `<parameter>{j}`. NULL when `prior` is not a fixed factor
# prior with design metadata.
.bt_factor_coordinate_display_names <- function(parameter, prior){

  prior <- .bt_parameter_catalog_factor_prior(parameter, prior)
  if(is.null(prior)){
    return(NULL)
  }
  design_info <- tryCatch(
    .factor_term_design_from_metadata(prior),
    error = function(error) NULL
  )
  if(is.null(design_info) || is.null(design_info$level_names)){
    return(NULL)
  }
  design <- as.matrix(design_info$design)
  out <- paste0(
    parameter,
    .bt_parameter_catalog_factor_coefficient_component(seq_len(ncol(design)))
  )
  direct <- .bt_factor_direct_cells(prior, design)
  if(any(!is.na(direct))){
    cell_names <- tryCatch(
      .format_factor_level_parameter_names(
        parameter,
        design_info$level_names,
        nrow(design)
      ),
      error = function(error){
        paste0(parameter, "[", design_info$cell_names, "]")
      }
    )
    out[!is.na(direct)] <- cell_names[direct[!is.na(direct)]]
  }
  out
}

# Maps every fixed factor term (formula terms of all contrasts and ordinary
# factor priors, whose level labels are 1..K by construction) to label-keyed
# quantities: one `<parameter>[<level token>]` quantity per level or cell, and
# one `<parameter>{j}` quantity per coordinate that is not a structural level
# cell. JAGS coordinate names are backend columns, never factor selectors.
.bt_parameter_catalog_factor_map <- function(coordinates, prior_list,
                                             derived = character()){

  out <- list(
    quantities = .bt_parameter_catalog_empty_quantities(),
    coordinates = character(),
    # Displayed estimates-table rows of these coordinates are level cells.
    representatives = stats::setNames(character(), character()),
    aliases = data.frame(
      alias = character(),
      quantity_id = character(),
      stringsAsFactors = FALSE
    )
  )
  if(length(prior_list) == 0L || is.null(names(prior_list))){
    return(out)
  }
  quantity_rows <- list()
  alias_rows <- list()
  for(parameter in names(prior_list)){
    prior <- .bt_parameter_catalog_factor_prior(
      parameter,
      prior_list[[parameter]]
    )
    if(is.null(prior)){
      next
    }
    expected_coordinate_names <- .JAGS_prior_factor_names(parameter, prior)
    coordinate_names <- .bt_parameter_catalog_coordinates(
      coordinates,
      expected_coordinate_names
    )
    coordinate_index <- suppressWarnings(as.integer(
      .bt_parameter_coordinates_index(coordinate_names)
    ))
    if(length(coordinate_names) > 1L){
      coordinate_names <- coordinate_names[order(coordinate_index)]
    }
    coordinate_rows <- match(coordinate_names, coordinates$coordinate_name)
    if(length(coordinate_names) != length(expected_coordinate_names) ||
       anyNA(coordinate_rows) ||
       (length(coordinate_names) > 1L && anyNA(coordinate_index)) ||
       length(unique(coordinates$role[coordinate_rows])) != 1L ||
       !coordinates$role[coordinate_rows][1L] %in%
         c("fixed_coefficient", "parameter")){
      stop(
        "Parameter catalog factor coordinates are missing or malformed for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    coordinate_metadata <- coordinates[coordinate_rows, , drop = FALSE]
    design_info <- .factor_term_design_from_metadata(prior)
    design <- design_info$design
    if(is.null(design_info$level_names)){
      stop(
        "Parameter catalog factor levels are missing for '", parameter,
        "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    cell_names <- .bt_parameter_catalog_factor_cell_names(design_info)
    if(ncol(design) != length(coordinate_names) ||
       nrow(design) != length(design_info$cell_names) ||
       length(cell_names) != nrow(design) ||
       any(!is.finite(design))){
      stop(
        "Parameter catalog factor metadata disagree with parameter coordinates for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    owner_fields <- c("formula_parameter", "term", "fitted_scale")
    if(any(vapply(owner_fields, function(field){
      length(unique(coordinate_metadata[[field]])) != 1L
    }, logical(1)))){
      stop(
        "Parameter catalog factor coordinates have inconsistent ownership for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    dif_names <- .factor_contrast_parameter_names(
      parameter,
      design_info$level_names,
      design_info$cell_names
    )
    if(length(dif_names) != nrow(design) ||
       anyNA(dif_names) || any(!nzchar(dif_names)) ||
       anyDuplicated(dif_names)){
      stop(
        "Parameter catalog factor semantic names are malformed for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }

    formula_parameter <- coordinate_metadata$formula_parameter[1L]
    namespace <- if(nzchar(formula_parameter)) formula_parameter else "model"
    role <- coordinate_metadata$role[1L]
    # Ordinary factor priors own no formula term; the prior name is the term.
    term <- if(nzchar(coordinate_metadata$term[1L])){
      coordinate_metadata$term[1L]
    }else{
      parameter
    }
    fitted_scale <- coordinate_metadata$fitted_scale[1L]
    coordinate_status <- ifelse(
      coordinate_names %in% derived,
      "derived",
      coordinate_metadata$monitor_status
    )
    factor_quantity <- function(canonical_name, component, status,
                                fixed_value, source_type, dependencies,
                                weights){
      .bt_parameter_catalog_quantity(
        canonical_name = canonical_name,
        namespace = namespace,
        role = role,
        formula_parameter = formula_parameter,
        term = term,
        component = component,
        display_label = .bt_parameter_catalog_factor_label(
          formula_parameter,
          term,
          component
        ),
        fitted_scale = fitted_scale,
        display_scale = fitted_scale,
        status = status,
        fixed_value = fixed_value,
        internal = FALSE,
        source_type = source_type,
        extraction_key = list(
          type = "factor_level",
          dependencies = dependencies,
          weights = weights
        )
      )
    }

    direct_cells <- .bt_factor_direct_cells(prior, design)

    for(cell in seq_len(nrow(design))){
      component <- cell_names[cell]
      direct <- match(cell, direct_cells)
      if(!is.na(direct)){
        status <- coordinate_status[direct]
        quantity <- factor_quantity(
          canonical_name = paste0(parameter, "[", component, "]"),
          component = component,
          status = status,
          fixed_value = if(identical(status, "structural")){
            coordinate_metadata$fixed_value[direct]
          }else{
            NA_real_
          },
          source_type = "identity",
          dependencies = coordinate_names[direct],
          weights = 1
        )
        out$representatives[[coordinate_names[direct]]] <- quantity$quantity_id
      }else{
        nonzero <- which(design[cell, ] != 0)
        dependencies <- coordinate_names[nonzero]
        weights <- unname(as.numeric(design[cell, nonzero]))
        dependency_status <- coordinate_status[nonzero]
        structural <- length(nonzero) == 0L ||
          all(dependency_status == "structural")
        unavailable <- length(nonzero) > 0L &&
          any(dependency_status == "unavailable")
        quantity <- factor_quantity(
          canonical_name = paste0(parameter, "[", component, "]"),
          component = component,
          status = if(structural){
            "structural"
          }else if(unavailable){
            "unavailable"
          }else{
            "derived"
          },
          fixed_value = if(!structural){
            NA_real_
          }else if(length(nonzero) == 0L){
            0
          }else{
            sum(weights * coordinate_metadata$fixed_value[nonzero])
          },
          source_type = if(length(dependencies) == 0L){
            "structural_zero"
          }else if(length(dependencies) == 1L){
            "identity"
          }else{
            "composite"
          },
          dependencies = dependencies,
          weights = weights
        )
      }
      quantity_rows[[length(quantity_rows) + 1L]] <- quantity
      # The level name of transformed summaries, `<parameter>[dif: <level>]`,
      # and the row labels those summaries display with and without the
      # formula prefix, `(mu) term[dif: <level>]` and `term[dif: <level>]`.
      dif_labels <- dif_names[cell]
      if(nzchar(formula_parameter)){
        dif_labels <- c(dif_labels, vapply(c(TRUE, FALSE), function(prefix){
          format_parameter_names(
            dif_names[cell],
            formula_parameters = formula_parameter,
            formula_prefix = prefix
          )
        }, character(1)))
      }
      alias_rows[[length(alias_rows) + 1L]] <- data.frame(
        alias = unique(dif_labels),
        quantity_id = quantity$quantity_id,
        stringsAsFactors = FALSE
      )
    }

    for(coordinate in which(is.na(direct_cells))){
      component <- .bt_parameter_catalog_factor_coefficient_component(
        coordinate
      )
      status <- coordinate_status[coordinate]
      quantity <- factor_quantity(
        canonical_name = paste0(parameter, component),
        component = component,
        status = status,
        fixed_value = if(identical(status, "structural")){
          coordinate_metadata$fixed_value[coordinate]
        }else{
          NA_real_
        },
        source_type = "identity",
        dependencies = coordinate_names[coordinate],
        weights = 1
      )
      quantity_rows[[length(quantity_rows) + 1L]] <- quantity
    }
    out$coordinates <- c(out$coordinates, coordinate_names)
  }
  if(length(quantity_rows) > 0L){
    out$quantities <- do.call(rbind, quantity_rows)
    rownames(out$quantities) <- NULL
  }
  if(length(alias_rows) > 0L){
    out$aliases <- do.call(rbind, alias_rows)
    rownames(out$aliases) <- NULL
  }
  out
}

# Coordinates owned by a point prior with an expression location are
# deterministic functions of other nodes: derived, never structural.
.bt_parameter_catalog_derived_coordinates <- function(coordinates, prior_list){

  if(length(prior_list) == 0L || is.null(names(prior_list))){
    return(character())
  }
  expression_points <- names(prior_list)[vapply(prior_list, function(prior){
    is.prior.point(prior) && .is_prior_expression(prior)
  }, logical(1))]
  bases <- .bt_parameter_coordinates_base(coordinates$coordinate_name)
  coordinates$coordinate_name[
    bases %in% expression_points & coordinates$monitor_status == "sampled"
  ]
}

.bt_parameter_catalog_coordinate_quantities <- function(
    coordinates, suppress = character(), derived = character()){

  out <- .bt_parameter_catalog_empty_quantities()
  keep <- coordinates$role != "backend_anchor" & !coordinates$internal &
    !coordinates$coordinate_name %in% suppress
  coordinates <- coordinates[keep, , drop = FALSE]
  if(nrow(coordinates) == 0L){
    return(out)
  }
  rows <- vector("list", nrow(coordinates))
  for(i in seq_len(nrow(coordinates))){
    row <- coordinates[i, , drop = FALSE]
    namespace <- if(nzchar(row$formula_parameter)){
      row$formula_parameter
    }else{
      "model"
    }
    rows[[i]] <- .bt_parameter_catalog_quantity(
      canonical_name = row$coordinate_name,
      namespace = namespace,
      role = row$role,
      formula_parameter = row$formula_parameter,
      term = row$term,
      component = row$column,
      display_label = row$display_label,
      fitted_scale = row$fitted_scale,
      display_scale = row$fitted_scale,
      status = if(row$coordinate_name %in% derived){
        "derived"
      }else{
        row$monitor_status
      },
      fixed_value = row$fixed_value,
      internal = FALSE,
      source_type = "identity",
      extraction_key = list(
        type = "coordinate",
        dependencies = row$coordinate_name
      )
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_level_alias <- function(term, component){

  if(startsWith(component, "{")){
    paste0(term, component)
  }else{
    paste0(term, "[", component, "]")
  }
}

.bt_parameter_catalog_aliases <- function(quantities, formula_design = NULL,
                                          table_labels = NULL,
                                          representatives = NULL,
                                          secondary = NULL){

  out <- .bt_parameter_catalog_empty_aliases()
  public <- quantities[!quantities$internal, , drop = FALSE]
  if(nrow(public) == 0L){
    return(out)
  }
  rows <- list()
  add_aliases <- function(quantity, values, simplified){
    values <- unique(values[!is.na(values) & nzchar(values)])
    if(length(values) == 0L){
      return(invisible(NULL))
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      alias = values,
      quantity_id = rep(quantity$quantity_id, length(values)),
      namespace = rep(quantity$namespace, length(values)),
      component = rep(quantity$component, length(values)),
      simplified = rep(simplified, length(values)),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  for(i in seq_len(nrow(public))){
    quantity <- public[i, , drop = FALSE]
    semantic_label <- character()
    if(startsWith(quantity$role, "random_") &&
       !is.na(quantity$formula_parameter)){
      prefix <- .bt_random_effect_summary_formula_prefix(
        quantity$formula_parameter,
        TRUE
      )
      if(nzchar(prefix) && startsWith(quantity$canonical_name, prefix)){
        semantic_label <- substring(
          quantity$canonical_name,
          nchar(prefix) + 1L
        )
      }
    }
    values <- if(startsWith(quantity$role, "random_")){
      unique(c(
        semantic_label,
        .bt_parameter_catalog_random_correlation_aliases(
          quantity,
          formula_design
        )
      ))
    }else{
      unique(c(
        quantity$canonical_name,
        quantity$display_label,
        quantity$term,
        if(nzchar(quantity$term) && nzchar(quantity$component) &&
           (identical(quantity$role, "fixed_coefficient") ||
              .bt_parameter_catalog_is_factor_quantity(quantity))){
          .bt_parameter_catalog_level_alias(
            quantity$term,
            quantity$component
          )
        }else{
          character()
        }
      ))
    }
    add_aliases(quantity, values, simplified = FALSE)
    if(startsWith(quantity$role, "random_")){
      add_aliases(
        quantity,
        .bt_parameter_catalog_random_simplified_aliases(quantity),
        simplified = TRUE
      )
    }
  }
  if(length(rows) == 0L){
    return(out)
  }
  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  .bt_parameter_catalog_table_label_aliases(
    out,
    public,
    table_labels,
    representatives = representatives,
    secondary = secondary
  )
}

.bt_parameter_catalog_is_factor_quantity <- function(quantity){

  key <- quantity$extraction_key[[1L]]
  is.list(key) && identical(key$type, "factor_level")
}

# Adds secondary labels as aliases unless the label already names a different
# quantity in the same namespace: displayed estimates-table labels of public
# coordinates (through the quantity representing each coordinate) and other
# supplied label-to-quantity rows, such as transformed factor-level names.
.bt_parameter_catalog_table_label_aliases <- function(aliases, public,
                                                      table_labels,
                                                      representatives = NULL,
                                                      secondary = NULL){

  added <- data.frame(
    alias = character(),
    quantity_id = character(),
    stringsAsFactors = FALSE
  )
  if(!is.null(table_labels) && nrow(table_labels) > 0L){
    coordinate_rows <- vapply(public$extraction_key, function(key){
      if(identical(key$type, "coordinate") && length(key$dependencies) == 1L){
        key$dependencies
      }else{
        NA_character_
      }
    }, character(1))
    owner <- public$quantity_id[
      match(table_labels$coordinate_name, coordinate_rows)
    ]
    if(length(representatives) > 0L){
      # A factor coordinate's table row names its level cell only when the
      # table relabelled it; an unrenamed indexed row shows a position.
      positional <- table_labels$coordinate_name == table_labels$renamed &
        grepl("\\[[^]]*\\]$", table_labels$coordinate_name)
      represented <- is.na(owner) & !positional &
        table_labels$coordinate_name %in% names(representatives)
      owner[represented] <- unname(
        representatives[table_labels$coordinate_name[represented]]
      )
    }
    keep <- !is.na(owner)
    added <- rbind(added, data.frame(
      alias = table_labels$alias[keep],
      quantity_id = owner[keep],
      stringsAsFactors = FALSE
    ))
  }
  if(!is.null(secondary) && nrow(secondary) > 0L){
    added <- rbind(added, secondary[c("alias", "quantity_id")])
  }
  quantity_rows <- match(added$quantity_id, public$quantity_id)
  added <- added[!is.na(quantity_rows), , drop = FALSE]
  quantity_rows <- quantity_rows[!is.na(quantity_rows)]
  if(nrow(added) == 0L){
    return(aliases)
  }
  added <- unique(data.frame(
    alias = added$alias,
    quantity_id = added$quantity_id,
    namespace = public$namespace[quantity_rows],
    component = public$component[quantity_rows],
    simplified = FALSE,
    stringsAsFactors = FALSE
  ))
  added <- added[!is.na(added$alias) & nzchar(added$alias), , drop = FALSE]
  if(nrow(added) == 0L){
    return(aliases)
  }

  label_key <- function(alias, namespace){
    paste(alias, namespace, sep = "\r")
  }
  # A label shown for several coordinates identifies none of them.
  added_key <- label_key(added$alias, added$namespace)
  added_owners <- tapply(added$quantity_id, added_key, function(x){
    length(unique(x))
  })
  added <- added[added_owners[added_key] == 1L, , drop = FALSE]

  # Existing selectors of another quantity in the namespace take precedence.
  existing <- unique(data.frame(
    key = c(
      label_key(aliases$alias, aliases$namespace),
      label_key(public$canonical_name, public$namespace)
    ),
    quantity_id = c(aliases$quantity_id, public$quantity_id),
    stringsAsFactors = FALSE
  ))
  added_key <- label_key(added$alias, added$namespace)
  existing_owners <- table(existing$key)
  first_owner <- existing$quantity_id[match(added_key, existing$key)]
  free <- is.na(first_owner) |
    (as.integer(existing_owners[added_key]) == 1L &
       first_owner == added$quantity_id)
  added <- added[free, , drop = FALSE]
  if(nrow(added) == 0L){
    return(aliases)
  }

  out <- unique(rbind(aliases, added))
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_random_simplified_aliases <- function(quantity){

  label <- quantity$display_label
  prefix <- .bt_random_effect_summary_formula_prefix(
    quantity$formula_parameter,
    TRUE
  )
  without_prefix <- if(nzchar(prefix) && startsWith(label, prefix)){
    substring(label, nchar(prefix) + 1L)
  }else{
    label
  }
  without_owner <- sub("^.*: ", "", without_prefix)

  unique(c(label, without_prefix, without_owner))
}

.bt_parameter_catalog_random_correlation_aliases <- function(
    quantity, formula_design){

  if(!identical(quantity$role, "random_correlation")){
    return(character())
  }
  key <- quantity$extraction_key[[1L]]
  random_term <- .bt_parameter_catalog_alias_random_term(
    formula_design,
    key
  )
  if(!is.list(random_term)){
    return(character())
  }
  owners <- .bt_parameter_catalog_random_public_block_owner(
    formula_design,
    quantity$formula_parameter,
    random_term
  )
  if(!is.character(owners) || anyNA(owners)){
    return(character())
  }

  if(!identical(key$evaluator, "rho")){
    return(character())
  }
  structure <- .bt_random_effect_summary_term_structure(random_term)
  if(!structure %in% c("cs", "hcs")){
    return(character())
  }
  components <- .bt_random_effect_summary_column_components(random_term)
  if(length(components) < 2L){
    return(character())
  }
  pairs <- utils::combn(components, 2L)
  correlations <- unlist(lapply(owners, function(owner){
    apply(pairs, 2L, function(pair){
      .bt_random_effect_semantic_name(
        parameter = quantity$formula_parameter,
        owner = owner,
        quantity = "cor",
        arguments = pair,
        formula_prefix = TRUE
      )
    })
  }), use.names = FALSE)
  prefix <- .bt_random_effect_summary_formula_prefix(
    quantity$formula_parameter,
    TRUE
  )
  without_prefix <- if(nzchar(prefix)){
    substring(correlations, nchar(prefix) + 1L)
  }else{
    correlations
  }
  unique(c(correlations, without_prefix))
}

.bt_parameter_catalog_random_public_block_owner <- function(
    formula_design, parameter, random_term){

  designs <- .bt_random_effect_summary_designs(formula_design)
  designs <- Filter(function(design){
    identical(design$parameter, parameter)
  }, designs)
  terms <- unlist(lapply(designs, `[[`, "random_effects"), recursive = FALSE)
  blocks <- unique(vapply(terms, `[[`, character(1), "block_name"))
  if(length(blocks) <= 1L && !isTRUE(random_term$component_visible)){
    return("")
  }

  .bt_random_effect_public_name(random_term)
}

.bt_parameter_catalog_alias_random_term <- function(formula_design, key){

  designs <- .bt_random_effect_summary_designs(formula_design)
  designs <- Filter(function(design){
    identical(design$parameter, key$formula_parameter)
  }, designs)
  terms <- unlist(lapply(designs, `[[`, "random_effects"), recursive = FALSE)
  terms <- Filter(function(term){
    identical(term$block_name, key$random_block)
  }, terms)

  if(length(terms) == 1L) terms[[1L]] else NULL
}

.bt_parameter_catalog_coordinates <- function(coordinates, names){

  names <- unique(names[!is.na(names) & nzchar(names)])
  if(length(names) == 0L){
    return(character())
  }
  bases <- .bt_parameter_coordinates_base(coordinates$coordinate_name)
  unique(coordinates$coordinate_name[
    coordinates$coordinate_name %in% names |
      coordinates$monitor_name %in% names |
      bases %in% names
  ])
}

.bt_parameter_catalog_random_block_dependencies <- function(
    coordinates, formula_parameter, random_block,
    roles = c("random_sd", "random_correlation")){

  coordinates$coordinate_name[
    coordinates$formula_parameter == formula_parameter &
      coordinates$random_block == random_block &
      coordinates$role %in% roles
  ]
}

.bt_parameter_catalog_random_correlation_sources <- function(coordinates,
                                                             random_term){

  structure <- .bt_random_effect_summary_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Parameter catalog"
  )
  if(is.null(correlation)){
    return(character())
  }
  if(identical(correlation$type, "lkj")){
    primitive_dependencies <- .bt_parameter_catalog_coordinates(
      coordinates,
      correlation$primitive_names
    )
    n_pairs <- random_term$n_columns * (random_term$n_columns - 1L) / 2L
    if(length(primitive_dependencies) == n_pairs){
      return(primitive_dependencies)
    }
    return(.bt_parameter_catalog_coordinates(
      coordinates,
      as.vector(.bt_random_effect_cholesky_names(
        random_term,
        random_term$n_columns
      ))
    ))
  }
  if(identical(correlation$type, "rho")){
    return(.bt_parameter_catalog_coordinates(
      coordinates,
      c(correlation$rho_name, correlation$sample_name)
    ))
  }
  character()
}

.bt_parameter_catalog_allocation_names <- function(allocation){

  if(is.null(allocation)){
    return(character())
  }
  source_name <- if(is.null(allocation$source)){
    character()
  }else{
    .bt_random_sd_binding_source_name(allocation$source)
  }
  factors <- c(allocation$factors, allocation$parent_factors)
  factor_names <- unlist(lapply(factors, function(factor){
    c(factor$weight_name, factor$inclusion_name)
  }), use.names = FALSE)
  inclusion_names <- unlist(lapply(allocation$inclusion, function(record){
    record$indicator_name
  }), use.names = FALSE)
  unique(c(
    allocation$weight_name,
    allocation$scale_name,
    allocation$source_node,
    source_name,
    factor_names,
    inclusion_names
  ))
}

.bt_parameter_catalog_allocation_scale_names <- function(allocation){

  if(is.null(allocation) || !is.list(allocation$source) ||
     !identical(allocation$source$shape, "scalar")){
    return(character())
  }
  factor_names <- unlist(lapply(allocation$parent_factors, function(factor){
    c(factor$weight_name, factor$inclusion_name)
  }), use.names = FALSE)
  unique(c(
    .bt_random_sd_binding_source_name(allocation$source),
    factor_names
  ))
}

# A scalar allocation source that is absent from the parameter map (for
# example, a source node that was not monitored) leaves every quantity scaled
# by it unevaluable.
.bt_parameter_catalog_allocation_source_missing <- function(allocation,
                                                             coordinates){

  source <- allocation$source
  is.list(source) && identical(source$shape, "scalar") &&
    length(.bt_parameter_catalog_coordinates(
      coordinates,
      .bt_random_sd_binding_source_name(source)
    )) == 0L
}

.bt_parameter_catalog_random_status <- function(
    key, coordinates, prior_list, formula_design, formula_scale){

  dependency_rows <- match(key$dependencies, coordinates$coordinate_name)
  if(anyNA(dependency_rows)){
    stop(
      "Parameter catalog random-summary dependencies are missing from the parameter map. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  dependency_status <- coordinates$monitor_status[dependency_rows]
  if(any(dependency_status == "unavailable")){
    return(list(status = "unavailable", fixed_value = NA_real_))
  }
  if(any(dependency_status == "sampled")){
    return(list(status = "sampled", fixed_value = NA_real_))
  }

  values <- matrix(
    coordinates$fixed_value[dependency_rows],
    nrow = 1L,
    dimnames = list(NULL, key$dependencies)
  )
  fit <- structure(
    list(),
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  fixed_value <- tryCatch(
    .bt_parameter_draw_random_summary(fit, key, values),
    error = function(error) error
  )
  if(inherits(fixed_value, "error") || length(fixed_value) != 1L ||
     !is.finite(fixed_value)){
    stop(
      "Parameter catalog could not evaluate a structural random summary from its declared dependencies. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  list(status = "structural", fixed_value = as.numeric(fixed_value))
}

# A one-to-one random summary (identity or square of one source coordinate) is
# structural exactly when that coordinate is structural, whatever the status
# of the other coordinates its block evaluator reads.
.bt_parameter_catalog_random_source_status <- function(source, coordinates,
                                                       source_transform){

  row <- match(source, coordinates$coordinate_name)
  if(length(source) != 1L || is.na(row) ||
     !source_transform %in% c("identity", "square")){
    stop(
      "Parameter catalog one-to-one random-summary source metadata are malformed. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  status <- coordinates$monitor_status[row]
  if(!identical(status, "structural")){
    return(list(status = status, fixed_value = NA_real_))
  }
  value <- coordinates$fixed_value[row]
  if(identical(source_transform, "square")){
    value <- value^2
  }

  list(status = "structural", fixed_value = as.numeric(value))
}

.bt_parameter_catalog_random_sd_is_direct <- function(
    random_term, parameter, formula_scale){

  if(is.null(formula_scale) || length(formula_scale) == 0L ||
     is.null(formula_scale[[parameter]]) ||
     length(formula_scale[[parameter]]) == 0L){
    return(TRUE)
  }
  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  parameter_scale <- formula_scale[[parameter]]
  column_groups <- .random_sd_column_unscale_groups(
    random_sd_cols = sd_names,
    formula_scale = parameter_scale,
    prefix = parameter
  )
  term_map <- if(is.null(column_groups)){
    .random_sd_term_map(sd_names, parameter_scale, parameter)
  }else{
    character()
  }
  (is.null(column_groups) || length(column_groups) == 0L) &&
    length(term_map) == 0L
}

.bt_parameter_catalog_random_sd_source_scale <- function(
    random_term, index, source_parameter, prior_list, parameter,
    formula_scale){

  # The caller has already proved a one-coordinate homogeneous SD transform
  # from the parameter map. Evaluate its unit response through the same
  # unscaling engine used for posterior draws.
  model_samples <- matrix(
    1,
    ncol = 1L,
    dimnames = list(NULL, source_parameter)
  )
  summary <- .bt_random_effect_summary_sd_samples(
    random_term = random_term,
    model_samples = model_samples,
    prior_list = prior_list,
    parameter = parameter,
    formula_scale = formula_scale
  )
  if(is.null(summary) || ncol(summary$values) < index){
    stop(
      "Parameter catalog could not resolve a one-coordinate random-SD transform. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  scale <- as.numeric(summary$values[, index])
  if(length(scale) != 1L || !is.finite(scale) || scale <= 0){
    stop(
      "Parameter catalog random-SD transform must have one finite positive source scale. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  unname(scale)
}

.bt_parameter_catalog_random_definitions <- function(coordinates, prior_list,
                                                     formula_design,
                                                     formula_scale = NULL){

  out <- list(
    derived = .bt_parameter_catalog_empty_quantities(),
    suppress = character()
  )
  if(is.null(prior_list)){
    prior_list <- list()
  }
  out$suppress <- unique(c(
    out$suppress,
    names(.bt_parameter_coordinates_random_prior_auxiliary_owners(prior_list)),
    .bt_random_variance_allocation_inclusion_indicator_names(formula_design)
  ))
  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(out)
  }

  rows <- list()
  used_names <- character()
  add_definition <- function(raw_name, role, parameter, label,
                             evaluator, dependencies, metadata = NULL,
                             block = "",
                             term = "", component = "",
                             fitted_scale = "fitted_covariance",
                             display_scale = "original",
                             owner_type, owner_name, quantity,
                             public_owner = owner_name,
                              scale_role = "", parent_quantity_id = "",
                              arguments = character(),
                              display_arguments = arguments,
                              source_type, source_parameter = "",
                              source_prior = "",
                              source_transform = "identity",
                              source_scale = NA_real_,
                              allocation_derived = FALSE,
                              status_source = NULL,
                              unavailable = FALSE){
    canonical_name <- .bt_random_effect_semantic_name(
      parameter = parameter,
      owner = public_owner,
      quantity = quantity,
      arguments = arguments,
      formula_prefix = TRUE
    )
    if(canonical_name %in% used_names){
      stop(
        "Random-effect semantic parameter names are not unique: '",
        canonical_name,
        "'. Assign unique random-block or variance-allocation names.",
        call. = FALSE
      )
    }
    used_names <<- c(used_names, canonical_name)
    namespace <- if(nzchar(parameter)) parameter else "model"
    dependencies <- unique(dependencies)
    if(identical(source_transform, "identity") && is.na(source_scale)){
      source_scale <- 1
    }
    key <- c(
      list(
        type = "random_summary",
        evaluator = evaluator,
        formula_parameter = parameter,
        random_block = block,
        summary_name = raw_name,
        dependencies = dependencies,
        source_type = source_type,
        source_parameter = source_parameter,
        source_prior = source_prior,
        source_transform = source_transform,
        source_scale = source_scale,
        allocation_derived = allocation_derived
      ),
      metadata[names(metadata) %in% c(
        "prior_name", "allocation_label", "parent_allocation", "index"
      )]
    )
    display_label <- .bt_random_effect_semantic_name(
      parameter = parameter,
      owner = public_owner,
      quantity = quantity,
      arguments = display_arguments,
      formula_prefix = TRUE
    )
    state <- if(isTRUE(unavailable)){
      list(status = "unavailable", fixed_value = NA_real_)
    }else{
      .bt_parameter_catalog_random_status(
        key = key,
        coordinates = coordinates,
        prior_list = prior_list,
        formula_design = formula_design,
        formula_scale = formula_scale
      )
    }
    if(!is.null(status_source) && !identical(state$status, "unavailable")){
      state <- .bt_parameter_catalog_random_source_status(
        source = status_source,
        coordinates = coordinates,
        source_transform = source_transform
      )
    }
    rows[[length(rows) + 1L]] <<- .bt_parameter_catalog_quantity(
      canonical_name = canonical_name,
      namespace = namespace,
      role = role,
      formula_parameter = parameter,
      term = term,
      component = component,
      display_label = display_label,
      fitted_scale = fitted_scale,
      display_scale = display_scale,
      status = state$status,
      fixed_value = state$fixed_value,
      owner_type = owner_type,
      owner_name = owner_name,
      quantity = quantity,
      scale_role = scale_role,
      parent_quantity_id = parent_quantity_id,
      arguments = arguments,
      source_type = source_type,
      extraction_key = key
    )
    invisible(NULL)
  }

  add_allocation <- function(allocation, parameter, random_term = NULL){
    if(is.null(allocation)){
      return(invisible(NULL))
    }
    K <- allocation$n_targets
    gate_only <- isTRUE(allocation$gate_only)
    if(!is.numeric(K) || length(K) != 1L || is.na(K) ||
       K < if(gate_only) 1L else 2L){
      stop("Random-effect allocation metadata have no valid 'n_targets'. Refit the model with this version of BayesTools.",
           call. = FALSE)
    }
    K <- as.integer(K)
    allocation_type <- .bt_random_effect_summary_allocation_type(allocation)
    allocation_owner <- .bt_random_effect_allocation_public_name(allocation)
    allocation_public_owner <- allocation_owner
    if(!is.null(random_term) && !nzchar(allocation_owner) && identical(
      .bt_random_effect_summary_allocation_target(allocation),
      "sd_component"
    )){
      allocation_owner <- .bt_random_effect_public_name(random_term)
      allocation_public_owner <-
        .bt_parameter_catalog_random_public_block_owner(
          formula_design,
          parameter,
          random_term
        )
    }
    components <- .bt_random_effect_summary_allocation_components(
      allocation,
      K = K,
      random_term = random_term
    )
    block <- if(is.null(random_term)) "" else random_term$block_name
    metadata <- list(
      allocation_label = allocation$label,
      parent_allocation = if(is.null(allocation$parent)) "" else
        allocation$parent$allocation,
      allocation = allocation,
      random_term = random_term
    )
    scale_role <- .bt_random_effect_allocation_scale_role(allocation)
    scale_names <- .bt_parameter_catalog_allocation_scale_names(allocation)
    # The allocation scale is evaluable only from its scalar source
    # coordinate; without it (a row-shaped or unmonitored source) the
    # sd/var totals are unavailable, and gate or weight coordinates must not
    # stand in as their dependencies.
    source_available <- length(scale_names) > 0L &&
      !.bt_parameter_catalog_allocation_source_missing(allocation, coordinates)
    allocation_gate_names <-
      .bt_random_effect_summary_allocation_gate_names(allocation)
    component_gate_names <-
      .bt_random_effect_summary_allocation_gate_names(
        allocation,
        include_parents = FALSE
      )
    realized_total <- identical(scale_role, "total") &&
      length(component_gate_names) > 0L
    if(realized_total && source_available){
      scale_names <- unique(c(
        scale_names,
        allocation$weight_name,
        component_gate_names
      ))
    }
    scale_dependencies <- .bt_parameter_catalog_coordinates(
      coordinates,
      scale_names
    )
    if(source_available && !gate_only){
      source_name <- .bt_random_sd_binding_source_name(allocation$source)
      parent_factors <- allocation$parent_factors
      if(is.null(parent_factors)){
        parent_factors <- list()
      }
      direct_scale <- length(parent_factors) == 0L && !realized_total
      source_prior <- if(direct_scale && source_name %in% names(prior_list)){
        source_name
      }else{
        ""
      }
      sd_quantity <- .bt_random_effect_allocation_sd_quantity(allocation)
      add_definition(
        raw_name = .bt_random_effect_summary_name(
          parameter = parameter,
          type = sd_quantity,
          parts = allocation$label
        ),
        role = paste0("random_", sd_quantity),
        parameter = parameter,
        label = allocation$label,
        evaluator = "allocation_sd",
        dependencies = scale_dependencies,
        metadata = metadata,
        block = block,
        term = allocation$label,
        component = allocation$label,
        owner_type = "variance_allocation",
        owner_name = allocation_owner,
        public_owner = allocation_public_owner,
        quantity = sd_quantity,
        scale_role = scale_role,
        source_type = if(direct_scale) "identity" else "composite",
        source_parameter = if(direct_scale) source_name else "",
        source_prior = source_prior,
        source_transform = "identity"
      )
      var_quantity <- .bt_random_effect_allocation_var_quantity(allocation)
      add_definition(
        raw_name = .bt_random_effect_summary_name(
          parameter = parameter,
          type = var_quantity,
          parts = allocation$label
        ),
        role = paste0("random_", var_quantity),
        parameter = parameter,
        label = allocation$label,
        evaluator = "allocation_var",
        dependencies = scale_dependencies,
        metadata = metadata,
        block = block,
        term = allocation$label,
        component = allocation$label,
        owner_type = "variance_allocation",
        owner_name = allocation_owner,
        public_owner = allocation_public_owner,
        quantity = var_quantity,
        scale_role = scale_role,
        source_type = if(direct_scale){
          "one_to_one_transform"
        }else{
          "composite"
        },
        source_parameter = if(direct_scale) source_name else "",
        source_prior = source_prior,
        source_transform = "square"
      )
    }
    out$suppress <<- unique(c(out$suppress, scale_dependencies))
    dependencies <- .bt_parameter_catalog_coordinates(
      coordinates,
      c(allocation$weight_name, allocation_gate_names)
    )
    out$suppress <<- unique(c(out$suppress, dependencies))
    if(!gate_only) for(i in seq_len(K)){
      raw_name <- .bt_random_effect_summary_name(
        parameter = parameter,
        type = allocation_type$name,
        parts = c(allocation$label, components[i])
      )
      add_definition(
        raw_name = raw_name,
        role = paste0("random_", allocation_type$summary),
        parameter = parameter,
        label = paste0(allocation_type$label, "(", allocation$label,
                       ": ", components[i], ")"),
        evaluator = "allocation",
        dependencies = dependencies,
        metadata = c(metadata, list(index = i)),
        block = block,
        term = allocation$label,
        component = components[i],
        fitted_scale = "unitless",
        display_scale = "unitless",
        owner_type = "variance_allocation",
        owner_name = allocation_owner,
        public_owner = allocation_public_owner,
        quantity = allocation_type$summary,
        scale_role = scale_role,
        arguments = components[i],
        source_type = if(identical(allocation_type$summary, "var_prop") &&
                         length(allocation_gate_names) == 0L){
          "identity"
        }else if(identical(allocation_type$summary, "var_prop")){
          "composite"
        }else{
          "one_to_one_transform"
        },
        source_parameter = if(identical(allocation_type$summary, "var_prop") &&
                              length(allocation_gate_names) > 0L){
          ""
        }else{
          allocation$weight_name
        },
        source_prior = if(identical(allocation_type$summary, "var_prop") &&
                          length(allocation_gate_names) > 0L){
          ""
        }else{
          allocation$weight_name
        },
        source_transform = allocation_type$summary
      )
      if(identical(.bt_random_effect_summary_allocation_target(allocation),
                   "sd_component")){
        multiplier_name <- .bt_random_effect_summary_name(
          parameter = parameter,
          type = "sd_mult",
          parts = c(allocation$label, components[i])
        )
        add_definition(
          raw_name = multiplier_name,
          role = "random_sd_mult",
          parameter = parameter,
          label = allocation$label,
          evaluator = "allocation",
          dependencies = dependencies,
          metadata = c(metadata, list(index = K + i)),
          block = block,
          term = allocation$label,
          component = components[i],
          fitted_scale = "unitless",
          display_scale = "unitless",
          owner_type = "variance_allocation",
          owner_name = allocation_owner,
          public_owner = allocation_public_owner,
          quantity = "sd_mult",
          scale_role = scale_role,
          arguments = components[i],
          source_type = "one_to_one_transform",
          source_parameter = allocation$weight_name,
          source_prior = allocation$weight_name,
          source_transform = "sd_mult"
        )
      }
    }
    inclusion <- allocation$inclusion
    if(!is.null(inclusion) && length(inclusion) > 0L){
      inclusion_i <- 0L
      for(component_label in names(inclusion)){
        inclusion_i <- inclusion_i + 1L
        component_name <- .bt_random_effect_allocation_component_name(
          allocation,
          component_label
        )
        raw_name <- .bt_random_effect_summary_name(
          parameter = parameter,
          type = "inclusion",
          parts = c(allocation$label, component_label)
        )
        inclusion_dependencies <- .bt_parameter_catalog_coordinates(
          coordinates,
          inclusion[[component_label]]$indicator_name
        )
        add_definition(
          raw_name = raw_name,
          role = "random_inclusion",
          parameter = parameter,
          label = paste0("inclusion(", allocation$label, ": ",
                         component_name, ")"),
          evaluator = "allocation_inclusion",
          dependencies = inclusion_dependencies,
          metadata = c(metadata, list(index = inclusion_i)),
          block = block,
          term = allocation$label,
          component = component_name,
          fitted_scale = "unitless",
          display_scale = "unitless",
          owner_type = "variance_allocation",
          owner_name = allocation_owner,
          public_owner = allocation_public_owner,
          quantity = "inclusion",
          arguments = component_name,
          source_type = "identity",
          source_parameter = inclusion[[component_label]]$indicator_name,
          source_prior = inclusion[[component_label]]$indicator_name
        )
      }
    }
    invisible(NULL)
  }

  seen_allocations <- character()
  for(design in random_design){
    parameter <- design$parameter
    for(random_term in design$random_effects){
      block <- random_term$block_name
      group <- .bt_random_effect_summary_group_label(random_term)
      owner <- .bt_random_effect_public_name(random_term)
      public_owner <- .bt_parameter_catalog_random_public_block_owner(
        formula_design,
        parameter,
        random_term
      )
      structure <- .bt_random_effect_summary_term_structure(random_term)
      term_metadata <- list(random_term = random_term)
      external_sd <- !is.null(random_term$sd_binding) &&
        .bt_random_sd_binding_has_external_source(random_term$sd_binding) &&
        !isTRUE(random_term$sd_binding$true_allocation)
      sd_names <- unique(random_term$sd_parameter_names)
      sd_names <- sd_names[!is.na(sd_names)]
      if(!external_sd && length(sd_names) > 0L){
        components <- .bt_random_effect_summary_sd_components(
          random_term,
          sd_names
        )
        direct_sd <- .bt_parameter_catalog_random_sd_is_direct(
          random_term = random_term,
          parameter = parameter,
          formula_scale = formula_scale
        )
        allocation_sd <- !is.null(random_term$sd_binding) &&
          length(random_term$sd_binding$allocations) > 0L
        allocation_derived <- allocation_sd && !all(vapply(
          random_term$sd_binding$allocations,
          function(allocation) isTRUE(allocation$gate_only),
          logical(1)
        ))
        direct_source <- direct_sd && !allocation_sd
        # SD rows scaled by an allocation source missing from the map cannot
        # be evaluated; they are unavailable, as are the allocation totals.
        allocation_source_missing <- allocation_sd && any(vapply(
          random_term$sd_binding$allocations,
          .bt_parameter_catalog_allocation_source_missing,
          logical(1),
          coordinates = coordinates
        ))
        allocation_dependencies <- .bt_parameter_catalog_coordinates(
          coordinates,
          unlist(lapply(
            random_term$sd_binding$allocations,
            .bt_parameter_catalog_allocation_names
          ), use.names = FALSE)
        )
        correlation_source_dependencies <-
          .bt_parameter_catalog_random_correlation_sources(
            coordinates,
            random_term
          )
        sd_dependencies <- if(allocation_sd){
          unique(c(
            allocation_dependencies,
            if(!direct_sd) correlation_source_dependencies else character()
          ))
        }else{
          unique(c(
            .bt_parameter_catalog_coordinates(coordinates, sd_names),
            if(!direct_sd) correlation_source_dependencies else character()
          ))
        }
        for(i in seq_along(sd_names)){
          sd_quantity <- "sd"
          sd_arguments <- .bt_random_effect_semantic_sd_arguments(
            components[i]
          )
          sd_display_arguments <-
            .bt_random_effect_semantic_sd_display_arguments(
              random_term,
              components[i]
            )
          raw_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "sd",
            parts = c(block, components[i])
          )
          label <- .bt_random_effect_sd_summary_label(
            component = components[i],
            random_term = random_term
          )
          out$suppress <- unique(c(
            out$suppress,
            intersect(sd_names[i], coordinates$coordinate_name)
          ))
          dependencies <- if(direct_source){
            .bt_parameter_catalog_coordinates(coordinates, sd_names)
          }else{
            sd_dependencies
          }
          source_coordinate <- .bt_parameter_catalog_coordinates(
            coordinates,
            sd_names[i]
          )
          source_type <- if(direct_source){
            "identity"
          }else if(length(dependencies) == 1L){
            if(allocation_sd) "composite" else "one_to_one_transform"
          }else{
            "composite"
          }
          source_prior <- .bt_parameter_coordinates_base(sd_names[i])
          if(!source_prior %in% names(prior_list)){
            source_prior <- ""
          }
          source_parameter <- if(direct_source &&
                                 length(source_coordinate) == 1L){
            source_coordinate
          }else if(!allocation_sd && length(dependencies) == 1L){
            dependencies
          }else{
            ""
          }
          source_transform <- if(direct_source) "identity" else "random_sd"
          status_source <- if(direct_source &&
                              length(source_coordinate) == 1L){
            source_coordinate
          }else{
            NULL
          }
          source_scale <- if(identical(source_type, "one_to_one_transform") &&
                             nzchar(source_parameter)){
            .bt_parameter_catalog_random_sd_source_scale(
              random_term = random_term,
              index = i,
              source_parameter = source_parameter,
              prior_list = prior_list,
              parameter = parameter,
              formula_scale = formula_scale
            )
          }else{
            NA_real_
          }
          add_definition(
            raw_name = raw_name,
            role = "random_sd",
            parameter = parameter,
            label = label,
            evaluator = "sd",
            dependencies = dependencies,
            metadata = c(term_metadata, list(index = i)),
            block = block,
            term = block,
            component = components[i],
            owner_type = "random_block",
            owner_name = owner,
            public_owner = public_owner,
            quantity = sd_quantity,
            arguments = sd_arguments,
            display_arguments = sd_display_arguments,
            source_type = source_type,
            source_parameter = source_parameter,
            source_prior = source_prior,
            source_transform = source_transform,
            source_scale = source_scale,
            allocation_derived = allocation_derived,
            status_source = status_source,
            unavailable = allocation_source_missing
          )
          var_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "var",
            parts = c(block, components[i])
          )
          add_definition(
            raw_name = var_name,
            role = "random_var",
            parameter = parameter,
            label = label,
            evaluator = "sd_variance",
            dependencies = dependencies,
            metadata = c(term_metadata, list(index = i)),
            block = block,
            term = block,
            component = components[i],
            owner_type = "random_block",
            owner_name = owner,
            public_owner = public_owner,
            quantity = "var",
            arguments = sd_arguments,
            display_arguments = sd_display_arguments,
            # var = sd^2 is one-to-one exactly when sd is; a scaled SD
            # (sd = s * source) gives var = (s * source)^2.
            source_type = if(direct_source ||
                             identical(source_type, "one_to_one_transform")){
              "one_to_one_transform"
            }else{
              "composite"
            },
            source_parameter = if(direct_source ||
                                  identical(source_type, "one_to_one_transform")){
              source_parameter
            }else{
              ""
            },
            source_prior = source_prior,
            source_transform = if(!direct_source &&
                                  identical(source_type, "one_to_one_transform")){
              "random_var"
            }else{
              "square"
            },
            source_scale = if(!direct_source &&
                              identical(source_type, "one_to_one_transform")){
              source_scale
            }else{
              NA_real_
            },
            allocation_derived = allocation_derived,
            status_source = status_source,
            unavailable = allocation_source_missing
          )
        }
      }

      inclusion_names <- .bt_random_effect_summary_inclusion_prior_names(
        random_term,
        prior_list
      )
      inclusion_i <- 0L
      for(prior_name in inclusion_names){
        prior <- prior_list[[prior_name]]
        prior_components <- attr(prior, "components", exact = TRUE)
        summary_components <- if(is.prior.spike_and_slab(prior)){
          "alternative"
        }else{
          unique(prior_components)
        }
        for(summary_component in summary_components){
          inclusion_i <- inclusion_i + 1L
          parts <- c(block, .bt_random_effect_summary_safe_label(prior_name))
          if(!is.prior.spike_and_slab(prior)){
            parts <- c(parts, summary_component)
          }
          raw_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "inclusion",
            parts = parts
          )
          effect <- .bt_random_effect_summary_inclusion_effect_label(
            random_term,
            prior_name
          )
          effect_parameter <- if(identical(effect, owner)){
            "sd"
          }else{
            paste0("sd(", effect, ")")
          }
          argument <- if(is.prior.spike_and_slab(prior)){
            effect_parameter
          }else{
            paste0(effect_parameter, "[", summary_component, "]")
          }
          label <- if(is.prior.spike_and_slab(prior)){
            paste0(effect, " | ", group, " (inclusion)")
          }else{
            paste0(effect, " | ", group, " (inclusion: ",
                   summary_component, ")")
          }
          add_definition(
            raw_name = raw_name,
            role = "random_inclusion",
            parameter = parameter,
            label = label,
            evaluator = "inclusion",
            dependencies = .bt_parameter_catalog_coordinates(
              coordinates,
              paste0(prior_name, "_indicator")
            ),
            metadata = c(term_metadata, list(index = inclusion_i)),
            block = block,
            term = block,
            component = summary_component,
            fitted_scale = "unitless",
            display_scale = "unitless",
            owner_type = "random_block",
            owner_name = owner,
            public_owner = public_owner,
            quantity = "inclusion",
            arguments = argument,
            source_type = "identity",
            source_parameter = paste0(prior_name, "_indicator"),
            source_prior = prior_name
          )
        }
      }

      correlation <- .bt_random_effect_correlation_metadata(
        random_term,
        structure = structure,
        context = "Parameter catalog"
      )
      if(!is.null(correlation) && identical(correlation$type, "rho")){
        raw_name <- .bt_random_effect_summary_name(
          parameter = parameter,
          type = "cor",
          parts = block
        )
        rho_scale <- .bt_random_effect_rho_scale_metadata(
          correlation,
          random_term,
          context = "Parameter catalog"
        )
        source_parameter <- if(identical(rho_scale, "rho")){
          correlation$rho_name
        }else{
          correlation$sample_name
        }
        dependencies <- .bt_parameter_catalog_coordinates(
          coordinates,
          c(correlation$rho_name, correlation$sample_name)
        )
        out$suppress <- unique(c(out$suppress, dependencies))
        source_prior <- correlation$prior_name
        if(is.null(source_prior) || length(source_prior) != 1L ||
           is.na(source_prior) || !source_prior %in% names(prior_list)){
          source_prior <- .bt_parameter_coordinates_base(source_parameter)
        }
        if(!source_prior %in% names(prior_list)){
          source_prior <- ""
        }
        add_definition(
          raw_name = raw_name,
          role = "random_correlation",
          parameter = parameter,
          label = owner,
          evaluator = "rho",
          dependencies = dependencies,
          metadata = term_metadata,
          block = block,
          term = block,
          fitted_scale = "fitted_covariance",
          display_scale = "unitless",
          owner_type = "random_block",
          owner_name = owner,
          public_owner = public_owner,
          quantity = "cor",
          source_type = if(identical(rho_scale, "rho")){
            "identity"
          }else{
            "one_to_one_transform"
          },
          source_parameter = source_parameter,
          source_prior = source_prior,
          source_transform = if(identical(rho_scale, "rho")){
            "identity"
          }else{
            rho_scale
          }
        )
      }
      if(!is.null(correlation) && identical(correlation$type, "lkj") &&
         random_term$n_columns > 1L){
        pairs <- utils::combn(seq_len(random_term$n_columns), 2L)
        components <- .bt_random_effect_summary_column_components(random_term)
        scaled_correlation <- !.bt_parameter_catalog_random_sd_is_direct(
          random_term = random_term,
          parameter = parameter,
          formula_scale = formula_scale
        )
        correlation_suppress <- .bt_parameter_catalog_random_block_dependencies(
          coordinates,
          formula_parameter = parameter,
          random_block = block,
          roles = c("random_sd", "random_correlation")
        )
        correlation_dependencies <-
          .bt_parameter_catalog_random_correlation_sources(
            coordinates,
            random_term
          )
        if(scaled_correlation){
          correlation_dependencies <- unique(c(
            correlation_dependencies,
            .bt_parameter_catalog_coordinates(coordinates, sd_names),
            .bt_parameter_catalog_coordinates(
              coordinates,
              unlist(lapply(
                random_term$sd_binding$allocations,
                .bt_parameter_catalog_allocation_names
              ), use.names = FALSE)
            )
          ))
        }
        out$suppress <- unique(c(out$suppress, correlation_suppress))
        one_to_one <- random_term$n_columns == 2L && !scaled_correlation &&
          length(correlation$primitive_names) == 1L
        for(i in seq_len(ncol(pairs))){
          pair <- components[pairs[, i]]
          raw_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "cor",
            parts = c(block, pair)
          )
          add_definition(
            raw_name = raw_name,
            role = "random_correlation",
            parameter = parameter,
            label = paste0("cor(", pair[1L], ",", pair[2L],
                           " | ", group, ")"),
            evaluator = "correlation",
            dependencies = correlation_dependencies,
            metadata = c(term_metadata, list(index = i)),
            block = block,
            term = block,
            component = paste0(pair[1L], ",", pair[2L]),
            fitted_scale = "fitted_covariance",
            display_scale = "unitless",
            owner_type = "random_block",
            owner_name = owner,
            public_owner = public_owner,
            quantity = "cor",
            arguments = pair,
            source_type = if(one_to_one){
              "one_to_one_transform"
            }else{
              "composite"
            },
            source_parameter = if(one_to_one){
              correlation$primitive_names
            }else{
              ""
            },
            source_transform = if(one_to_one) "lkj2" else "lkj"
          )
        }
      }

      allocations <- if(is.null(random_term$sd_binding)){
        list()
      }else{
        random_term$sd_binding$allocations
      }
      for(allocation in allocations){
        allocation_key <- paste(parameter, allocation$label, sep = "::")
        if(!allocation_key %in% seen_allocations){
          add_allocation(allocation, parameter, random_term)
          seen_allocations <- c(seen_allocations, allocation_key)
        }
      }
    }
    for(allocation in design$random_allocations){
      allocation_key <- paste(parameter, allocation$label, sep = "::")
      if(!allocation_key %in% seen_allocations){
        add_allocation(allocation, parameter)
        seen_allocations <- c(seen_allocations, allocation_key)
      }
    }
  }

  if(length(rows) > 0L){
    out$derived <- do.call(rbind, rows)
    rownames(out$derived) <- NULL
    allocation_rows <- which(
      out$derived$owner_type == "variance_allocation"
    )
    allocation_labels <- vapply(allocation_rows, function(i){
      key <- out$derived$extraction_key[[i]]
      if(!is.character(key$allocation_label) ||
         length(key$allocation_label) != 1L ||
         is.na(key$allocation_label) || !nzchar(key$allocation_label)){
        stop(
          "Parameter catalog variance-allocation metadata require one non-empty allocation label. Refit the model with this version of BayesTools.",
          call. = FALSE
        )
      }
      key$allocation_label
    }, character(1))
    # Allocation labels are unique only within one formula parameter.
    allocation_parameters <- out$derived$formula_parameter[allocation_rows]
    allocation_owner_key <- function(parameter, label){
      .bt_random_group_tuple_key(c(parameter, label))
    }
    allocation_keys <- vapply(seq_along(allocation_rows), function(row_i){
      allocation_owner_key(
        allocation_parameters[row_i],
        allocation_labels[row_i]
      )
    }, character(1))
    allocation_owners <- unique(allocation_keys)
    allocation_sd_ids <- stats::setNames(rep("", length(allocation_owners)),
                                          allocation_owners)
    for(owner in allocation_owners){
      candidates <- allocation_rows[
        allocation_keys == owner &
          out$derived$quantity[allocation_rows] %in% c("sd_total", "sd_common")
      ]
      if(length(candidates) == 1L){
        allocation_sd_ids[[owner]] <- out$derived$quantity_id[candidates]
      }
    }
    for(row_i in seq_along(allocation_rows)){
      i <- allocation_rows[row_i]
      owner <- allocation_keys[row_i]
      if(!out$derived$quantity[i] %in% c("sd_total", "sd_common")){
        out$derived$parent_quantity_id[i] <- allocation_sd_ids[[owner]]
        next
      }
      key <- out$derived$extraction_key[[i]]
      parent <- key$parent_allocation
      if(is.character(parent) && length(parent) == 1L && nzchar(parent)){
        parent_owner <- allocation_owner_key(allocation_parameters[row_i], parent)
        if(parent_owner %in% names(allocation_sd_ids)){
          out$derived$parent_quantity_id[i] <- allocation_sd_ids[[parent_owner]]
        }
      }
    }
  }
  out
}

.bt_build_parameter_catalog <- function(coordinates, prior_list = NULL,
                                        formula_design = NULL,
                                        formula_scale = NULL){

  .bt_validate_parameter_coordinates(coordinates)
  derived <- .bt_parameter_catalog_derived_coordinates(
    coordinates = coordinates,
    prior_list = prior_list
  )
  factor_map <- .bt_parameter_catalog_factor_map(
    coordinates = coordinates,
    prior_list = prior_list,
    derived = derived
  )
  random_map <- .bt_parameter_catalog_random_definitions(
    coordinates = coordinates,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  base <- .bt_parameter_catalog_coordinate_quantities(
    coordinates = coordinates,
    suppress = c(random_map$suppress, factor_map$coordinates),
    derived = derived
  )
  quantities <- rbind(base, factor_map$quantities, random_map$derived)
  rownames(quantities) <- NULL
  quantities <- .bt_parameter_catalog_add_support(
    quantities     = quantities,
    prior_list     = prior_list,
    formula_design = formula_design,
    formula_scale  = formula_scale
  )
  aliases <- .bt_parameter_catalog_aliases(
    quantities,
    formula_design,
    table_labels = .bt_parameter_catalog_table_labels(
      coordinates = coordinates,
      prior_list = prior_list,
      formula_scale = formula_scale
    ),
    representatives = factor_map$representatives,
    secondary = factor_map$aliases
  )
  .bt_parameter_catalog_new(quantities, aliases)
}

# Exact support and definedness of the catalog quantities, declared when the
# catalog is built from the prior provenance of the quantities' source
# coordinates (never from draws). A support that is not derivable is NULL.
.bt_parameter_catalog_add_support <- function(quantities, prior_list,
                                              formula_design = NULL,
                                              formula_scale = NULL){

  if(nrow(quantities) == 0L){
    return(quantities)
  }
  object <- structure(
    list(),
    prior_list     = prior_list,
    formula_design = formula_design,
    formula_scale  = formula_scale
  )
  quantities$support <- I(lapply(seq_len(nrow(quantities)), function(i){
    .bt_parameter_catalog_quantity_support(
      object   = object,
      quantity = quantities[i, , drop = FALSE]
    )
  }))
  quantities$definedness <- vapply(seq_len(nrow(quantities)), function(i){
    .bt_parameter_catalog_quantity_definedness(
      object   = object,
      quantity = quantities[i, , drop = FALSE]
    )
  }, character(1))
  quantities
}

.bt_parameter_catalog_quantity_support <- function(object, quantity){

  if(identical(quantity$status, "unavailable")){
    return(NULL)
  }
  if(identical(quantity$status, "structural")){
    return(.posterior_support_point(quantity$fixed_value, source = "catalog"))
  }
  key <- quantity$extraction_key[[1L]]
  if(key$type %in% c("coordinate", "factor_level")){
    weights <- if(identical(key$type, "coordinate")){
      rep(1, length(key$dependencies))
    }else{
      key$weights
    }
    return(.bt_parameter_catalog_linear_support(
      prior_list = attr(object, "prior_list", exact = TRUE),
      weights    = stats::setNames(as.numeric(weights), key$dependencies)
    ))
  }
  if(!identical(key$type, "random_summary")){
    return(NULL)
  }

  # one-to-one quantities: the source prior's support mapped by the transform
  if(key$source_type %in% c("identity", "one_to_one_transform")){
    source_support <- .bt_parameter_catalog_source_support(object, key)
    transform <- .bt_parameter_transform_from_quantity(object, quantity)
    if(!is.null(source_support) && !is.null(transform)){
      return(.bt_parameter_catalog_transform_support(source_support, transform))
    }
  }

  .bt_parameter_catalog_random_support(object, quantity, key)
}

# Support of a linear combination of fitted coordinates under their priors
# (a coefficient's 'multiply_by' scales only its linear-predictor
# contribution); NULL when a coordinate has no owning prior.
.bt_parameter_catalog_linear_support <- function(prior_list, weights){

  if(!is.list(prior_list) || length(prior_list) == 0L || length(weights) == 0L){
    return(NULL)
  }
  prior_list <- .marginal_posterior_strip_multiply_by(prior_list)
  owned <- unlist(lapply(names(prior_list), function(parameter){
    if(!is.prior(prior_list[[parameter]])){
      return(NULL)
    }
    .prior_linear_prior_columns(parameter, prior_list[[parameter]])
  }), use.names = FALSE)
  if(!all(names(weights) %in% owned)){
    return(NULL)
  }
  .posterior_support_from_prior_list_weights(prior_list, weights, source = "catalog")
}

# Support of the source prior of a one-to-one random-effect quantity: the
# named prior, or the Beta(eta, eta) of an LKJ pairwise correlation.
.bt_parameter_catalog_source_support <- function(object, key){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  if(nzchar(key$source_prior) && key$source_prior %in% names(prior_list)){
    source_prior <- prior_list[[key$source_prior]]
    if(is.prior.simplex(source_prior)){
      # a Dirichlet weight lies in (0, 1)
      return(.posterior_support_new(c(0, 1), source = "catalog", type = "interval"))
    }
    if(!is.prior(source_prior) || .prior_linear_prior_dimension(source_prior) != 1L){
      return(NULL)
    }
    return(.posterior_support_from_prior(source_prior, source = "catalog"))
  }
  if(identical(key$source_transform, "lkj2")){
    return(.posterior_support_new(c(0, 1), source = "catalog", type = "interval"))
  }

  NULL
}

# A support mapped by a monotone transform descriptor (parameter_transform()).
.bt_parameter_catalog_transform_support <- function(support, transform){

  map <- function(values){
    if(length(values) == 0L){
      return(values)
    }
    parameter_transform_forward(values, transform)
  }
  if(identical(transform$type, "square") && support$bounds[1L] < 0){
    # a square is monotone only on a nonnegative source
    return(NULL)
  }
  bounds <- map(support$bounds)
  points <- map(support$points)
  if(anyNA(bounds) || anyNA(points)){
    return(NULL)
  }

  .posterior_support_new(
    bounds = range(bounds),
    points = points,
    exact  = support$exact,
    source = "catalog",
    type   = support$type
  )
}

# Supports of composite random-effect quantities: correlations lie in
# [-1, 1], variance shares in [0, 1], multipliers between zero and the
# allocation scale, and inclusion indicators at 0 and 1. Composite SDs and
# variances are nonnegative; their hull [0, Inf) is exact only when every
# scale prior is supported on [0, Inf).
.bt_parameter_catalog_random_support <- function(object, quantity, key){

  interval <- function(lower, upper, exact = TRUE){
    .posterior_support_new(c(lower, upper), exact = exact, source = "catalog",
                           type = "interval")
  }
  switch(
    quantity$quantity,
    "cor"       = interval(-1, 1),
    "var_prop"  = interval(0, 1),
    "inclusion" = .posterior_support_new(c(0, 1), points = c(0, 1), source = "catalog",
                                         type = "points"),
    "var_mult"  = ,
    "sd_mult"   = {
      scale <- .bt_parameter_catalog_allocation_variance_scale(object, key)
      if(is.null(scale)){
        NULL
      }else if(identical(quantity$quantity, "var_mult")){
        interval(0, scale)
      }else{
        interval(0, sqrt(scale))
      }
    },
    "sd"        = ,
    "var"       = ,
    "sd_total"  = ,
    "var_total" = ,
    "sd_common" = ,
    "var_common" = interval(
      0, Inf,
      exact = .bt_parameter_catalog_unbounded_scale(object, key)
    ),
    NULL
  )
}

# The variance scale of an allocation's multipliers: the number of targets
# for mean-variance allocations, 1 for total-variance allocations.
.bt_parameter_catalog_allocation_variance_scale <- function(object, key){

  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(object, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(object, key, random_term)
  if(is.null(allocation)){
    return(NULL)
  }
  scale <- .bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Parameter catalog support"
  )
  if(identical(scale, "mean_variance")){
    return(as.numeric(.bt_random_effect_summary_allocation_n_targets(
      allocation,
      K = allocation$n_targets
    )))
  }
  1
}

# Whether the scale priors of a composite random-effect SD or variance are all
# supported on [0, Inf): a composite SD is positively homogeneous in the
# scales, so it then reaches every value in (0, Inf). A scale bounded away
# from zero can bound the composite away from zero as well (e.g. the
# original-scale intercept SD sqrt(sd_0^2 + c^2 sd_1^2) of a scaled
# predictor), and a scale bounded above bounds it above, so the hull is then
# not exact.
.bt_parameter_catalog_unbounded_scale <- function(object, key){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  scale_priors <- prior_list[intersect(key$dependencies, names(prior_list))]
  scale_priors <- Filter(function(prior){
    is.prior(prior) && !is.prior.simplex(prior) && !is.prior.discrete(prior) &&
      .prior_linear_prior_dimension(prior) == 1L
  }, scale_priors)
  if(length(scale_priors) == 0L){
    return(FALSE)
  }
  all(vapply(scale_priors, function(prior){
    support <- .posterior_support_from_prior(prior, source = "catalog")
    !is.null(support) && isTRUE(support$exact) && support$bounds[1L] == 0 &&
      is.infinite(support$bounds[2L])
  }, logical(1)))
}

# Draws where a catalog quantity is defined: original-scale correlations are
# undefined where an SD is zero, and variance shares of gated total-variance
# allocations where no component is active.
.bt_parameter_catalog_quantity_definedness <- function(object, quantity){

  key <- quantity$extraction_key[[1L]]
  if(!identical(key$type, "random_summary")){
    return("always")
  }
  if(identical(key$evaluator, "correlation")){
    return("correlation")
  }
  if(identical(key$evaluator, "allocation") &&
     identical(quantity$quantity, "var_prop")){
    random_term <- if(nzchar(key$random_block)){
      .bt_parameter_catalog_find_random_term(object, key)
    }else{
      NULL
    }
    allocation <- .bt_parameter_catalog_find_allocation(object, key, random_term)
    if(!is.null(allocation) &&
       identical(.bt_random_effect_allocation_scale_metadata(
         allocation,
         context = "Parameter catalog definedness"
       ), "total_variance") &&
       length(.bt_random_effect_summary_allocation_gate_names(allocation)) > 0L){
      return("allocation_active")
    }
  }

  "always"
}

# Row labels that JAGS_estimates_table() and runjags_estimates_table() display
# for public coordinates: factor levels renamed as those tables do, with and
# without the formula prefix, and with the log-intercept exp(intercept) label
# of transformed tables. A label that cannot be computed is simply not added.
.bt_parameter_catalog_table_labels <- function(coordinates, prior_list,
                                               formula_scale = NULL){

  empty <- data.frame(
    coordinate_name = character(),
    renamed = character(),
    alias = character(),
    stringsAsFactors = FALSE
  )
  public <- !coordinates$internal & coordinates$role != "backend_anchor"
  coordinate_names <- coordinates$coordinate_name[public]
  if(length(coordinate_names) == 0L || length(prior_list) == 0L){
    return(empty)
  }
  renamed <- tryCatch(
    .bt_random_effect_summary_renamed_parameter_names(
      coordinate_names,
      prior_list
    ),
    error = function(error) NULL
  )
  if(!is.character(renamed) || length(renamed) != length(coordinate_names)){
    return(empty)
  }
  formula_parameters <- unique(unlist(
    lapply(prior_list, attr, which = "parameter", exact = TRUE),
    use.names = FALSE
  ))
  formula_random <- unique(unlist(
    lapply(prior_list, attr, which = "random_factor"),
    use.names = FALSE
  ))
  settings <- expand.grid(
    prefix = c(TRUE, FALSE),
    scaled = c(FALSE, TRUE)
  )
  labels <- lapply(seq_len(nrow(settings)), function(i){
    format_parameter_names(
      renamed,
      formula_parameters = formula_parameters,
      formula_random = formula_random,
      formula_prefix = settings$prefix[i],
      formula_scale = if(settings$scaled[i]) formula_scale else NULL
    )
  })
  out <- data.frame(
    coordinate_name = rep(coordinate_names, length(labels)),
    renamed = rep(renamed, length(labels)),
    alias = unlist(labels, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  out <- unique(out[!is.na(out$alias) & nzchar(out$alias), , drop = FALSE])
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_valid_native_key <- function(key, quantity){

  if(identical(key$type, "coordinate")){
    return(
      length(key$dependencies) == 1L &&
        identical(key$dependencies, quantity$canonical_name)
    )
  }
  if(identical(key$type, "factor_level")){
    return(
      is.numeric(key$weights) && !anyNA(key$weights) &&
        all(is.finite(key$weights)) &&
        length(key$weights) == length(key$dependencies)
    )
  }
  if(!identical(key$type, "random_summary")){
    return(FALSE)
  }
  scalar_character <- function(value, allow_empty = FALSE){
    is.character(value) && length(value) == 1L && !is.na(value) &&
      (allow_empty || nzchar(value))
  }
  evaluators <- c(
    "allocation_sd", "allocation_var", "sd", "sd_variance", "inclusion",
    "rho", "correlation",
    "allocation", "allocation_inclusion"
  )
  if(!scalar_character(key$evaluator) ||
     !key$evaluator %in% evaluators ||
     !scalar_character(key$formula_parameter, allow_empty = TRUE) ||
     !scalar_character(key$random_block, allow_empty = TRUE) ||
     !scalar_character(key$summary_name) ||
     !scalar_character(key$source_type) ||
     !key$source_type %in% c(
       "identity", "one_to_one_transform", "composite"
     ) ||
     !scalar_character(key$source_parameter, allow_empty = TRUE) ||
     !scalar_character(key$source_prior, allow_empty = TRUE) ||
     !scalar_character(key$source_transform) ||
     !is.numeric(key$source_scale) || length(key$source_scale) != 1L){
    return(FALSE)
  }
  if(key$evaluator %in% c("sd", "sd_variance", "inclusion", "correlation")){
    return(
      is.numeric(key$index) && length(key$index) == 1L &&
        !is.na(key$index) && key$index == as.integer(key$index) &&
        key$index >= 1L
    )
  }
  if(key$evaluator %in% c(
    "allocation_sd", "allocation_var", "allocation", "allocation_inclusion"
  )){
    return(
      scalar_character(key$allocation_label) &&
        (key$evaluator %in% c("allocation_sd", "allocation_var") ||
           (is.numeric(key$index) && length(key$index) == 1L &&
              !is.na(key$index) && key$index == as.integer(key$index) &&
              key$index >= 1L))
    )
  }
  identical(key$evaluator, "rho")
}

.bt_validate_parameter_catalog_tables <- function(
    quantities, aliases,
    known_quantity_ids = quantities$quantity_id){

  valid <- is.data.frame(quantities) && is.data.frame(aliases) &&
    identical(names(quantities), .bt_parameter_catalog_quantity_columns) &&
    identical(names(aliases), .bt_parameter_catalog_alias_columns)
  if(!valid){
    stop("Parameter catalog tables have a missing or unsupported schema. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  character_columns <- setdiff(
    .bt_parameter_catalog_quantity_columns,
    c("arguments", "fixed_value", "internal", "support", "extraction_key")
  )
  if(!all(vapply(quantities[character_columns], is.character, logical(1))) ||
     !is.list(quantities$arguments) ||
     !all(vapply(quantities$arguments, function(arguments){
       is.character(arguments) && !anyNA(arguments)
     }, logical(1))) ||
     !is.numeric(quantities$fixed_value) ||
     !is.logical(quantities$internal) ||
     !is.list(quantities$extraction_key) ||
     !all(vapply(
       aliases[setdiff(names(aliases), "simplified")],
       is.character,
       logical(1)
     )) ||
     !is.logical(aliases$simplified) ||
     !is.list(quantities$support) ||
     !all(vapply(quantities$support, function(support){
       is.null(support) || inherits(support, "BayesTools_posterior_support")
     }, logical(1))) ||
     anyNA(quantities[setdiff(names(quantities),
                             c("fixed_value", "support", "extraction_key"))]) ||
     anyNA(aliases)){
    stop("Parameter catalog tables contain malformed field types or missing metadata. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  required_nonempty <- c(
    "quantity_id", "canonical_name", "provider", "namespace", "role",
    "display_label", "fitted_scale", "display_scale", "status",
    "source_type"
  )
  if(any(!nzchar(as.matrix(quantities[required_nonempty]))) ||
     any(!nzchar(aliases$alias)) || any(!nzchar(aliases$quantity_id)) ||
     any(!nzchar(aliases$namespace)) ||
     anyDuplicated(quantities$quantity_id) ||
     any(!quantities$source_type %in%
           c("identity", "one_to_one_transform", "composite",
             "structural_zero", "none")) ||
     any(!quantities$status %in%
           c("sampled", "structural", "derived", "unavailable")) ||
     any(!quantities$definedness %in%
           c("always", names(.bt_undefined_draws_reasons))) ||
     any(!is.na(quantities$fixed_value[quantities$status != "structural"])) ||
     any(!is.finite(quantities$fixed_value[quantities$status == "structural"]))){
    stop("Parameter catalog tables contain invalid names, statuses, or structural values. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  # `quantity_id` is the key; `canonical_name` is a selector, and
  # `parameter_catalog_resolve()` narrows it by namespace and component before
  # reporting a typed ambiguity. Two providers describing one term under the
  # same public name is therefore legitimate - that is what
  # `parameter_catalog_extend()` exists to produce - and only rows that the
  # resolver could not tell apart are rejected here.
  selector <- duplicated(quantities[c("canonical_name", "namespace", "component")])
  if(any(selector)){
    stop(
      "Parameter catalog quantities repeat the selector '",
      quantities$canonical_name[selector][[1L]],
      "' within one namespace and component, so it cannot be resolved. ",
      "Give the extending provider a distinct canonical name or component.",
      call. = FALSE
    )
  }
  provider_prefix <- paste0(quantities$provider, "::")
  if(any(!startsWith(quantities$quantity_id, provider_prefix))){
    stop("Parameter catalog quantity IDs do not match their providers. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  if(any(!aliases$quantity_id %in% known_quantity_ids)){
    stop("Parameter catalog aliases reference unknown quantity IDs. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  if(any(!quantities$scale_role %in% c("", "total", "common")) ||
     any(nzchar(quantities$parent_quantity_id) &
           !quantities$parent_quantity_id %in% known_quantity_ids) ||
     any(quantities$parent_quantity_id == quantities$quantity_id)){
    stop("Parameter catalog quantities contain invalid scale hierarchy metadata. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  random <- startsWith(quantities$role, "random_")
  if(any(random & (
    !nzchar(quantities$owner_type) |
      (!nzchar(quantities$owner_name) &
         quantities$owner_type != "variance_allocation") |
      !nzchar(quantities$quantity) |
      quantities$source_type == "none"
  ))){
    stop("Parameter catalog random quantities contain incomplete semantic ownership or source metadata. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  valid_keys <- vapply(seq_len(nrow(quantities)), function(i){
    key <- quantities$extraction_key[[i]]
    valid <- is.list(key) && is.character(key$type) && length(key$type) == 1L &&
      !is.na(key$type) && nzchar(key$type) &&
      is.character(key$dependencies) && !anyNA(key$dependencies) &&
      !anyDuplicated(key$dependencies)
    if(!isTRUE(valid)){
      return(FALSE)
    }
    if(!identical(quantities$provider[i], "BayesTools")){
      return(TRUE)
    }
    .bt_parameter_catalog_valid_native_key(
      key,
      quantities[i, , drop = FALSE]
    )
  }, logical(1))
  if(!all(valid_keys)){
    stop("Parameter catalog extraction keys are malformed. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.bt_validate_parameter_catalog <- function(catalog){

  valid <- inherits(catalog, "BayesTools_parameter_catalog") &&
    is.list(catalog) &&
    identical(names(catalog), c("schema_version", "quantities", "aliases")) &&
    identical(catalog$schema_version, .bt_parameter_map_version)
  if(!valid){
    stop("Parameter catalog metadata are missing or unsupported. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  .bt_validate_parameter_catalog_tables(catalog$quantities, catalog$aliases)
  invisible(TRUE)
}

.bt_parameter_catalog_stop <- function(class, message, ...){

  condition <- structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(class, "BayesTools_parameter_resolution_error", "error",
              "condition")
  )
  stop(condition)
}

.bt_validate_parameter_selection <- function(selection, catalog = NULL){

  valid <- inherits(selection, "BayesTools_parameter_selection") &&
    is.list(selection) &&
    identical(
      names(selection),
      c("schema_version", "parameter_map_version", "quantity_id",
        "quantities")
    ) &&
    identical(selection$schema_version, .bt_parameter_selection_version) &&
    identical(selection$parameter_map_version,
              .bt_parameter_map_version) &&
    is.character(selection$quantity_id) &&
    length(selection$quantity_id) == nrow(selection$quantities)
  if(!valid){
    stop("Parameter selection metadata are missing or unsupported. Resolve the parameter again with this version of BayesTools.",
         call. = FALSE)
  }
  .bt_validate_parameter_catalog_tables(
    selection$quantities,
    .bt_parameter_catalog_empty_aliases(),
    known_quantity_ids = unique(c(
      selection$quantities$quantity_id,
      selection$quantities$parent_quantity_id[
        nzchar(selection$quantities$parent_quantity_id)
      ]
    ))
  )
  if(!identical(selection$quantity_id, selection$quantities$quantity_id)){
    stop("Parameter selection IDs and quantity metadata disagree. Resolve the parameter again.",
         call. = FALSE)
  }
  if(!is.null(catalog)){
    expected <- match(selection$quantity_id, catalog$quantities$quantity_id)
    if(anyNA(expected)){
      stop("The parameter selection does not belong to this fitted catalog.",
           call. = FALSE)
    }
    catalog_rows <- catalog$quantities[expected, , drop = FALSE]
    if(!identical(selection$quantities, catalog_rows)){
      stop("The parameter selection is stale or does not belong to this fitted catalog.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

.bt_parameter_draw_dependencies <- function(fit, dependencies){

  JAGS_materialize_draws(
    fit,
    parameters = dependencies,
    include_internal = TRUE
  )
}

.bt_parameter_draw_supplied_dependencies <- function(model_samples,
                                                      dependencies){

  if(!is.matrix(model_samples) || !is.numeric(model_samples) ||
     is.null(colnames(model_samples)) || anyDuplicated(colnames(model_samples))){
    stop("'model_samples' must be a numeric matrix with unique column names.",
         call. = FALSE)
  }
  missing <- setdiff(dependencies, colnames(model_samples))
  if(length(missing) > 0L){
    stop(
      "The supplied 'model_samples' are missing declared source coordinates: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  values <- model_samples[, dependencies, drop = FALSE]

  coda::mcmc.list(coda::mcmc(values))
}

.bt_parameter_draw_factor_level <- function(key, model_samples){

  dependency_names <- colnames(model_samples)
  valid_dependencies <- if(length(key$dependencies) == 0L){
    ncol(model_samples) == 0L
  }else{
    identical(dependency_names, key$dependencies)
  }
  if(!isTRUE(valid_dependencies)){
    stop(
      "Selected factor level is unavailable from its declared dependencies.",
      call. = FALSE
    )
  }
  as.vector(model_samples %*% key$weights)
}

.bt_parameter_catalog_find_random_term <- function(fit, key){

  designs <- attr(fit, "formula_design", exact = TRUE)
  design <- designs[[key$formula_parameter]]
  if(is.null(design)){
    stop("Selected random summary has no matching formula design. Refit the model with this version of BayesTools.",
         call. = FALSE)
  }
  matches <- vapply(design$random_effects, function(term){
    identical(term$block_name, key$random_block)
  }, logical(1))
  if(sum(matches) != 1L){
    stop("Selected random summary has no unique random-effect block. Refit the model with this version of BayesTools.",
         call. = FALSE)
  }
  design$random_effects[[which(matches)]]
}

.bt_parameter_catalog_find_allocation <- function(fit, key, random_term){

  designs <- attr(fit, "formula_design", exact = TRUE)
  design <- designs[[key$formula_parameter]]
  allocations <- list()
  if(!is.null(random_term) && !is.null(random_term$sd_binding)){
    allocations <- random_term$sd_binding$allocations
  }
  allocations <- c(allocations, design$random_allocations)
  allocation_keys <- vapply(allocations, function(allocation){
    allocation$label
  }, character(1))
  allocations <- allocations[!duplicated(allocation_keys)]
  matches <- vapply(allocations, function(allocation){
    identical(allocation$label, key$allocation_label)
  }, logical(1))
  if(sum(matches) != 1L){
    stop("Selected random summary has no unique variance-allocation definition. Refit the model with this version of BayesTools.",
         call. = FALSE)
  }
  allocations[[which(matches)]]
}

.bt_parameter_catalog_complete_cholesky_triangle <- function(random_term,
                                                             model_samples){

  structure <- .bt_random_effect_summary_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Parameter catalog"
  )
  if(is.null(correlation) || !identical(correlation$type, "lkj") ||
     random_term$n_columns < 2L){
    return(model_samples)
  }
  names <- .bt_random_effect_cholesky_names(
    random_term,
    random_term$n_columns
  )
  upper_names <- names[upper.tri(names)]
  missing <- upper_names[!upper_names %in% colnames(model_samples)]
  if(length(missing) == 0L){
    return(model_samples)
  }
  zeroes <- matrix(
    0,
    nrow = nrow(model_samples),
    ncol = length(missing),
    dimnames = list(NULL, missing)
  )
  cbind(model_samples, zeroes)
}

.bt_parameter_draw_random_summary <- function(fit, key, model_samples){

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  evaluator <- key$evaluator
  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  if(!is.null(random_term)){
    model_samples <- .bt_parameter_catalog_complete_cholesky_triangle(
      random_term,
      model_samples
    )
  }
  if(evaluator %in% c("allocation_sd", "allocation_var")){
    allocation <- .bt_parameter_catalog_find_allocation(
      fit,
      key,
      random_term
    )
    values <- .bt_random_effect_summary_allocation_scale_samples(
      allocation = allocation,
      model_samples = model_samples,
      prior_list = prior_list
    )
    if(is.null(values)){
      stop("Selected random allocation scale is unavailable from its declared dependencies.",
           call. = FALSE)
    }
    return(if(identical(evaluator, "allocation_var")) values^2 else values)
  }else if(evaluator %in% c("sd", "sd_variance")){
    summary <- .bt_random_effect_summary_sd_samples(
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      parameter = key$formula_parameter,
      formula_scale = attr(fit, "formula_scale", exact = TRUE)
    )
    values <- summary$values[, key$index]
    return(if(identical(evaluator, "sd_variance")) values^2 else values)
  }
  if(identical(evaluator, "inclusion")){
    summary <- .bt_random_effect_summary_inclusion_samples(
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      parameter = key$formula_parameter
    )
    match <- match(key$summary_name, summary$names)
    if(is.na(match)){
      stop("Selected random inclusion summary is unavailable from its declared dependencies.",
           call. = FALSE)
    }
    return(summary$values[, match])
  }
  if(identical(evaluator, "rho")){
    return(.bt_random_effect_summary_rho_samples(random_term, model_samples))
  }
  if(identical(evaluator, "correlation")){
    complete <- .bt_random_effect_summary_complete_scaled_samples(
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      parameter = key$formula_parameter,
      formula_scale = attr(fit, "formula_scale", exact = TRUE)
    )
    summary <- .bt_random_effect_summary_correlation_samples(
      random_term,
      complete
    )
    return(summary$values[, key$index])
  }
  if(identical(evaluator, "allocation_inclusion")){
    values <- .bt_random_effect_allocation_gate_draws(
      parameter_name = key$source_parameter,
      posterior = model_samples
    )
    if(is.null(values)){
      stop(
        "Selected random allocation inclusion is unavailable from its declared dependency.",
        call. = FALSE
      )
    }
    return(values)
  }
  if(identical(evaluator, "allocation")){
    allocation <- .bt_parameter_catalog_find_allocation(
      fit,
      key,
      random_term
    )
    summary <- .bt_random_effect_summary_allocation_samples(
      allocation = allocation,
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      include_multipliers = TRUE
    )
    match <- match(key$summary_name, summary$names)
    if(is.na(match)){
      stop("Selected random allocation summary is unavailable from its declared dependencies.",
           call. = FALSE)
    }
    return(summary$values[, match])
  }

  stop("Unsupported BayesTools parameter-catalog extraction key.", call. = FALSE)
}
