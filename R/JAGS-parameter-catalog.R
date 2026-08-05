# Metadata-only semantic parameter catalog and deferred draw extraction.

.bt_parameter_catalog_version <- 1L
.bt_parameter_selection_version <- 1L

.bt_parameter_catalog_quantity_columns <- c(
  "quantity_id", "canonical_name", "provider", "namespace", "role",
  "formula_parameter", "term", "component", "display_label",
  "fitted_scale", "display_scale", "status", "fixed_value", "internal",
  "extraction_key"
)
.bt_parameter_catalog_alias_columns <- c(
  "alias", "quantity_id", "namespace", "component"
)

#' Semantic parameter catalogs and deferred draw extraction
#'
#' @description
#' `parameter_catalog()` returns the versioned semantic catalog cached on a
#' fitted object. The catalog contains only metadata: selectable quantities,
#' exact aliases, and serializable extraction keys. Constructing or resolving
#' it never accesses posterior draws.
#'
#' `parameter_catalog_extend()` adds plain-data quantities and aliases owned by
#' another provider. `parameter_catalog_resolve()` applies optional namespace
#' and component filters and returns a versioned selection only when the match
#' is unique.
#'
#' `parameter_draws()` is the deferred extraction boundary. The BayesTools fit
#' method reads only the registry coordinates declared by the selected
#' extraction key. Downstream packages can provide methods for package-owned
#' derived quantities.
#'
#' @param object fitted object.
#' @param catalog a `BayesTools_parameter_catalog` object.
#' @param quantities quantity rows matching the quantity schema.
#' @param aliases alias rows matching the alias schema.
#' @param provider scalar provider name owning every added quantity.
#' @param alias scalar exact canonical name or alias.
#' @param namespace optional exact namespace filter.
#' @param component optional exact component filter.
#' @param selection a `BayesTools_parameter_selection` returned by
#'   `parameter_catalog_resolve()`.
#' @param ... arguments for methods.
#'
#' @return `parameter_catalog()` and `parameter_catalog_extend()` return a
#' `BayesTools_parameter_catalog`. `parameter_catalog_schema()` returns schema
#' descriptions. `parameter_catalog_resolve()` returns a
#' `BayesTools_parameter_selection`. `parameter_draws()` returns a
#' `coda::mcmc.list` for BayesTools-owned quantities.
#'
#' @export parameter_catalog
#' @export parameter_catalog_schema
#' @export parameter_catalog_extend
#' @export parameter_catalog_resolve
#' @export parameter_draws
#' @name parameter_catalog
NULL

#' @rdname parameter_catalog
parameter_catalog <- function(object, ...){

  UseMethod("parameter_catalog")
}

#' @rdname parameter_catalog
#' @exportS3Method parameter_catalog BayesTools_fit
parameter_catalog.BayesTools_fit <- function(object, ...){

  JAGS_validate_fit_contract(object, requires = "parameter_catalog")
  catalog <- attr(object, "parameter_catalog", exact = TRUE)
  .bt_validate_parameter_catalog(catalog)
  catalog
}

#' @rdname parameter_catalog
parameter_catalog_schema <- function(){

  quantities <- data.frame(
    field = .bt_parameter_catalog_quantity_columns,
    type = c(rep("character", 12L), "numeric", "logical", "list"),
    description = c(
      "Stable provider-namespaced quantity identifier.",
      "Exact canonical semantic name.",
      "Package or subsystem owning extraction.",
      "Exact resolver namespace.",
      "Semantic role.",
      "Owning formula output parameter, or an empty string.",
      "Formula term or semantic subterm, or an empty string.",
      "Quantity component, or an empty string.",
      "Unambiguous default display label.",
      "Scale of stored dependency coordinates.",
      "Scale of the selected quantity.",
      "One of sampled, structural, derived, or unavailable.",
      "Exact structural value; otherwise NA.",
      "Whether the quantity is private implementation metadata.",
      "Serializable plain-data extraction recipe."
    ),
    stringsAsFactors = FALSE
  )
  aliases <- data.frame(
    field = .bt_parameter_catalog_alias_columns,
    type = rep("character", 4L),
    description = c(
      "Exact accepted alias.",
      "Quantity identifier targeted by the alias.",
      "Exact resolver namespace.",
      "Optional component filter, or an empty string."
    ),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = .bt_parameter_catalog_version,
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
                                      component = NULL){

  .bt_validate_parameter_catalog(catalog)
  check_char(alias, "alias", check_length = 1L, allow_NA = FALSE)
  check_char(namespace, "namespace", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_char(component, "component", check_length = 1L, allow_NULL = TRUE,
             allow_NA = FALSE)

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
  alias_rows <- catalog$aliases$alias == alias
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
        catalog$aliases$quantity_id %in% quantities$quantity_id[public]
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
    catalog_schema_version = catalog$schema_version,
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
#' @exportS3Method parameter_draws BayesTools_fit
parameter_draws.BayesTools_fit <- function(object, selection, ...){

  catalog <- parameter_catalog(object)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  quantities <- selection$quantities
  if(any(quantities$provider != "BayesTools")){
    stop(
      "The selection contains quantities owned by another provider; use that provider's 'parameter_draws()' method.",
      call. = FALSE
    )
  }

  ordinary <- vapply(
    quantities$extraction_key,
    function(key) identical(key$type, "registry"),
    logical(1)
  )
  if(all(ordinary)){
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
  dependencies <- .bt_parameter_draw_dependencies(object, key$dependencies)
  out <- vector("list", length(dependencies))
  for(chain_i in seq_along(dependencies)){
    chain <- dependencies[[chain_i]]
    values <- .bt_parameter_draw_random_summary(
      fit = object,
      key = key,
      model_samples = as.matrix(chain)
    )
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

.bt_parameter_catalog_empty_quantities <- function(){

  out <- data.frame(
    quantity_id = character(),
    canonical_name = character(),
    provider = character(),
    namespace = character(),
    role = character(),
    formula_parameter = character(),
    term = character(),
    component = character(),
    display_label = character(),
    fitted_scale = character(),
    display_scale = character(),
    status = character(),
    fixed_value = numeric(),
    internal = logical(),
    stringsAsFactors = FALSE
  )
  out$extraction_key <- I(list())
  out
}

.bt_parameter_catalog_empty_aliases <- function(){

  data.frame(
    alias = character(),
    quantity_id = character(),
    namespace = character(),
    component = character(),
    stringsAsFactors = FALSE
  )
}

.bt_parameter_catalog_new <- function(quantities, aliases){

  out <- list(
    schema_version = .bt_parameter_catalog_version,
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
    extraction_key){

  out <- .bt_parameter_catalog_empty_quantities()
  out[1L, setdiff(names(out), "extraction_key")] <- list(
    .bt_parameter_catalog_quantity_id(canonical_name, namespace, role),
    canonical_name,
    "BayesTools",
    namespace,
    role,
    formula_parameter,
    term,
    component,
    display_label,
    fitted_scale,
    display_scale,
    status,
    fixed_value,
    internal
  )
  out$extraction_key <- I(list(extraction_key))
  out
}

.bt_parameter_catalog_registry_quantities <- function(registry){

  out <- .bt_parameter_catalog_empty_quantities()
  keep <- registry$role != "backend_anchor"
  registry <- registry[keep, , drop = FALSE]
  if(nrow(registry) == 0L){
    return(out)
  }
  rows <- vector("list", nrow(registry))
  for(i in seq_len(nrow(registry))){
    row <- registry[i, , drop = FALSE]
    namespace <- if(nzchar(row$formula_parameter)){
      row$formula_parameter
    }else{
      "model"
    }
    rows[[i]] <- .bt_parameter_catalog_quantity(
      canonical_name = row$canonical_name,
      namespace = namespace,
      role = row$role,
      formula_parameter = row$formula_parameter,
      term = row$term,
      component = row$column,
      display_label = row$display_label,
      fitted_scale = row$fitted_scale,
      display_scale = row$fitted_scale,
      status = row$monitor_status,
      fixed_value = row$fixed_value,
      internal = row$internal,
      extraction_key = list(
        type = "registry",
        dependencies = row$canonical_name
      )
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_aliases <- function(quantities){

  out <- .bt_parameter_catalog_empty_aliases()
  public <- quantities[!quantities$internal, , drop = FALSE]
  if(nrow(public) == 0L){
    return(out)
  }
  rows <- vector("list", nrow(public))
  for(i in seq_len(nrow(public))){
    quantity <- public[i, , drop = FALSE]
    semantic_label <- character()
    if(startsWith(quantity$role, "random_") &&
       !is.na(quantity$formula_parameter)){
      prefix <- .bt_random_effect_summary_formula_prefix(
        quantity$formula_parameter,
        TRUE
      )
      if(nzchar(prefix) && startsWith(quantity$display_label, prefix)){
        semantic_label <- substring(
          quantity$display_label,
          nchar(prefix) + 1L
        )
      }
    }
    values <- unique(c(
      quantity$canonical_name,
      quantity$display_label,
      semantic_label,
      quantity$term,
      quantity$component
    ))
    values <- values[!is.na(values) & nzchar(values)]
    rows[[i]] <- data.frame(
      alias = values,
      quantity_id = rep(quantity$quantity_id, length(values)),
      namespace = rep(quantity$namespace, length(values)),
      component = rep(quantity$component, length(values)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_character_values <- function(x){

  if(is.character(x)){
    return(x[!is.na(x) & nzchar(x)])
  }
  if(!is.list(x) || length(x) == 0L){
    return(character())
  }
  unique(unlist(lapply(x, .bt_parameter_catalog_character_values),
                use.names = FALSE))
}

.bt_parameter_catalog_dependencies <- function(registry, metadata = NULL,
                                               formula_parameter = "",
                                               random_block = ""){

  available <- registry$monitor_status %in% c("sampled", "structural")
  tokens <- .bt_parameter_catalog_character_values(metadata)
  bases <- .bt_parameter_registry_base(registry$canonical_name)
  token_match <- registry$canonical_name %in% tokens |
    registry$monitor_name %in% tokens |
    bases %in% tokens
  block_match <- nzchar(random_block) &
    registry$random_block == random_block &
    (!nzchar(formula_parameter) |
       registry$formula_parameter == formula_parameter)
  unique(registry$canonical_name[available & (token_match | block_match)])
}

.bt_parameter_catalog_random_definitions <- function(registry, prior_list,
                                                     formula_design,
                                                     formula_scale = NULL){

  out <- .bt_parameter_catalog_empty_quantities()
  if(is.null(prior_list)){
    prior_list <- list()
  }
  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(out)
  }

  rows <- list()
  used_names <- registry$canonical_name
  add_definition <- function(raw_name, role, parameter, label,
                             evaluator, metadata = NULL, block = "",
                             term = "", component = "",
                             fitted_scale = "fitted_covariance",
                             display_scale = "original"){
    canonical_name <- .bt_random_effect_summary_unique_name(
      raw_name,
      used_names
    )
    used_names <<- c(used_names, canonical_name)
    namespace <- if(nzchar(parameter)) parameter else "model"
    dependencies <- .bt_parameter_catalog_dependencies(
      registry = registry,
      metadata = metadata,
      formula_parameter = parameter,
      random_block = block
    )
    key <- c(
      list(
        type = "random_summary",
        evaluator = evaluator,
        formula_parameter = parameter,
        random_block = block,
        summary_name = raw_name,
        dependencies = dependencies
      ),
      metadata[names(metadata) %in% c("prior_name", "allocation_label", "index")]
    )
    display_label <- paste0(
      .bt_random_effect_summary_formula_prefix(parameter, TRUE),
      label
    )
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
      extraction_key = key
    )
    invisible(NULL)
  }

  for(prior_name in names(prior_list)){
    prior <- prior_list[[prior_name]]
    if(!isTRUE(attr(prior, "random_sd_total", exact = TRUE))){
      next
    }
    parameter <- attr(prior, "parameter", exact = TRUE)
    allocation <- attr(prior, "random_allocation", exact = TRUE)
    raw_name <- .bt_random_effect_summary_name(
      parameter = parameter,
      type = "sd_total",
      parts = allocation
    )
    add_definition(
      raw_name = raw_name,
      role = "random_sd_total",
      parameter = parameter,
      label = paste0("sd_total(", allocation, ")"),
      evaluator = "sd_total",
      metadata = list(prior_name = prior_name),
      term = allocation,
      component = allocation
    )
  }

  add_allocation <- function(allocation, parameter, random_term = NULL){
    if(is.null(allocation)){
      return(invisible(NULL))
    }
    K <- allocation$n_targets
    if(!is.numeric(K) || length(K) != 1L || is.na(K) || K < 2L){
      stop("Random-effect allocation metadata have no valid 'n_targets'. Refit the model with this version of BayesTools.",
           call. = FALSE)
    }
    K <- as.integer(K)
    allocation_type <- .bt_random_effect_summary_allocation_type(allocation)
    components <- .bt_random_effect_summary_allocation_components(
      allocation,
      K = K,
      random_term = random_term
    )
    block <- if(is.null(random_term)) "" else random_term$block_name
    metadata <- list(
      allocation_label = allocation$label,
      allocation = allocation,
      random_term = random_term
    )
    for(i in seq_len(K)){
      raw_name <- .bt_random_effect_summary_name(
        parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
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
        metadata = c(metadata, list(index = i)),
        block = block,
        term = allocation$label,
        component = components[i],
        fitted_scale = "unitless",
        display_scale = "unitless"
      )
      if(identical(.bt_random_effect_summary_allocation_target(allocation),
                   "sd_component")){
        multiplier_name <- .bt_random_effect_summary_name(
          parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
          type = "sd_mult",
          parts = c(allocation$label, components[i])
        )
        add_definition(
          raw_name = multiplier_name,
          role = "random_sd_multiplier",
          parameter = parameter,
          label = paste0("sd_mult(", allocation$label, ": ",
                         components[i], ")"),
          evaluator = "allocation",
          metadata = c(metadata, list(index = K + i)),
          block = block,
          term = allocation$label,
          component = components[i],
          fitted_scale = "unitless",
          display_scale = "unitless"
        )
      }
    }
    inclusion <- allocation$inclusion
    if(!is.null(inclusion) && length(inclusion) > 0L){
      inclusion_i <- 0L
      for(component_label in names(inclusion)){
        inclusion_i <- inclusion_i + 1L
        raw_name <- .bt_random_effect_summary_name(
          parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
          type = "inclusion",
          parts = c(allocation$label, component_label)
        )
        add_definition(
          raw_name = raw_name,
          role = "random_inclusion",
          parameter = parameter,
          label = paste0("inclusion(", allocation$label, ": ",
                         component_label, ")"),
          evaluator = "allocation_inclusion",
          metadata = c(metadata, list(index = inclusion_i)),
          block = block,
          term = allocation$label,
          component = component_label,
          fitted_scale = "unitless",
          display_scale = "unitless"
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
        for(i in seq_along(sd_names)){
          raw_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "sd",
            parts = c(block, components[i])
          )
          add_definition(
            raw_name = raw_name,
            role = "random_sd",
            parameter = parameter,
            label = .bt_random_effect_sd_summary_label(
              component = components[i],
              group = group,
              random_term = random_term
            ),
            evaluator = "sd",
            metadata = c(term_metadata, list(index = i)),
            block = block,
            term = block,
            component = components[i]
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
            metadata = c(term_metadata, list(index = inclusion_i)),
            block = block,
            term = block,
            component = summary_component,
            fitted_scale = "unitless",
            display_scale = "unitless"
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
          type = "rho",
          parts = block
        )
        add_definition(
          raw_name = raw_name,
          role = "random_correlation",
          parameter = parameter,
          label = paste0("rho(", group, ")"),
          evaluator = "rho",
          metadata = term_metadata,
          block = block,
          term = block,
          fitted_scale = "fitted_covariance",
          display_scale = "unitless"
        )
      }
      if(!is.null(correlation) && identical(correlation$type, "lkj") &&
         random_term$n_columns > 1L){
        pairs <- utils::combn(seq_len(random_term$n_columns), 2L)
        components <- .bt_random_effect_summary_column_components(random_term)
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
            metadata = c(term_metadata, list(index = i)),
            block = block,
            term = block,
            component = paste0(pair[1L], ",", pair[2L]),
            fitted_scale = "fitted_covariance",
            display_scale = "unitless"
          )
        }
      }

      allocations <- if(is.null(random_term$sd_binding)){
        list()
      }else{
        random_term$sd_binding$allocations
      }
      for(allocation in allocations){
        if(!allocation$weight_name %in% seen_allocations){
          add_allocation(allocation, parameter, random_term)
          seen_allocations <- c(seen_allocations, allocation$weight_name)
        }
      }
    }
    for(allocation in design$random_allocations){
      if(!allocation$weight_name %in% seen_allocations){
        add_allocation(allocation, parameter)
        seen_allocations <- c(seen_allocations, allocation$weight_name)
      }
    }
  }

  if(length(rows) == 0L){
    return(out)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_build_parameter_catalog <- function(registry, prior_list = NULL,
                                        formula_design = NULL,
                                        formula_scale = NULL){

  .bt_validate_parameter_registry(registry)
  base <- .bt_parameter_catalog_registry_quantities(registry)
  derived <- .bt_parameter_catalog_random_definitions(
    registry = registry,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  quantities <- rbind(base, derived)
  rownames(quantities) <- NULL
  aliases <- .bt_parameter_catalog_aliases(quantities)
  .bt_parameter_catalog_new(quantities, aliases)
}

.bt_attach_parameter_catalog <- function(fit){

  if(inherits(fit, "error")){
    return(fit)
  }
  registry <- JAGS_parameter_registry(fit)
  attr(fit, "parameter_catalog") <- .bt_build_parameter_catalog(
    registry = registry,
    prior_list = attr(fit, "prior_list", exact = TRUE),
    formula_design = attr(fit, "formula_design", exact = TRUE),
    formula_scale = attr(fit, "formula_scale", exact = TRUE)
  )
  fit
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
    c("fixed_value", "internal", "extraction_key")
  )
  if(!all(vapply(quantities[character_columns], is.character, logical(1))) ||
     !is.numeric(quantities$fixed_value) ||
     !is.logical(quantities$internal) ||
     !is.list(quantities$extraction_key) ||
     !all(vapply(aliases, is.character, logical(1))) ||
     anyNA(quantities[setdiff(names(quantities),
                             c("fixed_value", "extraction_key"))]) ||
     anyNA(aliases)){
    stop("Parameter catalog tables contain malformed field types or missing metadata. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
  }
  required_nonempty <- c(
    "quantity_id", "canonical_name", "provider", "namespace", "role",
    "display_label", "fitted_scale", "display_scale", "status"
  )
  if(any(!nzchar(as.matrix(quantities[required_nonempty]))) ||
     any(!nzchar(aliases$alias)) || any(!nzchar(aliases$quantity_id)) ||
     any(!nzchar(aliases$namespace)) ||
     anyDuplicated(quantities$quantity_id) ||
     any(!quantities$status %in%
           c("sampled", "structural", "derived", "unavailable")) ||
     any(!is.na(quantities$fixed_value[quantities$status != "structural"])) ||
     any(!is.finite(quantities$fixed_value[quantities$status == "structural"]))){
    stop("Parameter catalog tables contain invalid names, statuses, or structural values. Refit or rebuild the catalog with this version of BayesTools.",
         call. = FALSE)
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
  valid_keys <- vapply(quantities$extraction_key, function(key){
    is.list(key) && is.character(key$type) && length(key$type) == 1L &&
      !is.na(key$type) && nzchar(key$type) &&
      is.character(key$dependencies) && !anyNA(key$dependencies) &&
      !anyDuplicated(key$dependencies)
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
    identical(catalog$schema_version, .bt_parameter_catalog_version)
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
      c("schema_version", "catalog_schema_version", "quantity_id",
        "quantities")
    ) &&
    identical(selection$schema_version, .bt_parameter_selection_version) &&
    identical(selection$catalog_schema_version,
              .bt_parameter_catalog_version) &&
    is.character(selection$quantity_id) &&
    length(selection$quantity_id) == nrow(selection$quantities)
  if(!valid){
    stop("Parameter selection metadata are missing or unsupported. Resolve the parameter again with this version of BayesTools.",
         call. = FALSE)
  }
  .bt_validate_parameter_catalog_tables(
    selection$quantities,
    .bt_parameter_catalog_empty_aliases()
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
  allocations <- design$random_allocations
  if(!is.null(random_term) && !is.null(random_term$sd_binding)){
    allocations <- c(allocations, random_term$sd_binding$allocations)
  }
  allocation_keys <- vapply(allocations, function(allocation){
    allocation$weight_name
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

.bt_parameter_draw_random_summary <- function(fit, key, model_samples){

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  evaluator <- key$evaluator
  if(identical(evaluator, "sd_total")){
    return(.bt_random_effect_parameter_draws(
      key$prior_name,
      model_samples,
      prior_list
    ))
  }

  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  if(identical(evaluator, "sd")){
    summary <- .bt_random_effect_summary_sd_samples(
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      parameter = key$formula_parameter,
      formula_scale = attr(fit, "formula_scale", exact = TRUE)
    )
    return(summary$values[, key$index])
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
  if(evaluator %in% c("allocation", "allocation_inclusion")){
    allocation <- .bt_parameter_catalog_find_allocation(
      fit,
      key,
      random_term
    )
    summary <- if(identical(evaluator, "allocation")){
      .bt_random_effect_summary_allocation_samples(
        allocation = allocation,
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        include_multipliers = TRUE
      )
    }else{
      .bt_random_effect_summary_allocation_inclusion_samples(
        allocation,
        model_samples
      )
    }
    match <- match(key$summary_name, summary$names)
    if(is.na(match)){
      stop("Selected random allocation summary is unavailable from its declared dependencies.",
           call. = FALSE)
    }
    return(summary$values[, match])
  }

  stop("Unsupported BayesTools parameter-catalog extraction key.", call. = FALSE)
}
