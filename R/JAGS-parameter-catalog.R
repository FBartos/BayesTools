# Metadata-only semantic parameter catalog and deferred draw extraction.

.bt_parameter_catalog_version <- 5L
.bt_parameter_selection_version <- 1L

.bt_parameter_catalog_quantity_columns <- c(
  "quantity_id", "canonical_name", "provider", "namespace", "role",
  "formula_parameter", "owner_type", "owner_name", "quantity", "scale_role",
  "parent_quantity_id", "arguments",
  "term", "component", "display_label",
  "fitted_scale", "display_scale", "status", "fixed_value", "internal",
  "source_type", "extraction_key"
)
.bt_parameter_catalog_alias_columns <- c(
  "alias", "quantity_id", "namespace", "component"
)

#' Semantic parameter catalogs and deferred draw extraction
#'
#' @description
#' `parameter_catalog()` returns the versioned semantic catalog cached on a
#' fitted object. The catalog contains only metadata: selectable quantities,
#' exact aliases, and serializable extraction keys. Fixed factor terms expose
#' every fitted level or interaction cell as a catalog component, so both
#' `term[level]` and `term` plus `component = "level"` resolve without parsing
#' backend coordinate names. Direct coordinates are reused; reference cells are
#' structural zeroes; and contrast-coded cells are reconstructed from the
#' persisted term-only design matrix. Ordinary level labels remain unchanged;
#' syntax-sensitive characters are percent-escaped and ambiguous interaction
#' tokens are quoted so that every component remains hypothesis-safe and
#' injective. These are coefficient-level quantities, distinct from estimated
#' marginal means based on full predictions.
#' [JAGS_parameter_registry()] remains the concrete posterior-coordinate map;
#' the catalog is the separate semantic layer. Random-effect canonical names
#' follow `(formula) owner: quantity(arguments)`: parentheses contain parameter
#' or coefficient names, while square brackets inside an argument contain its
#' factor or index level. Examples include `(mu) study: sd(intercept)`,
#' `(mu) study: cor(group[sensitivity],group[specificity])`, and
#' `(mu) heterogeneity: var_prop(study)`. Formula-prefix omission is accepted as
#' an alias. Additional aliases are limited to genuine semantic equivalences,
#' such as a CS/HCS pairwise correlation referring to its shared `cor`.
#'
#' Identity random-effect summaries reuse their sampled or structural registry
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
    type = c(
      rep("character", 11L), "list", rep("character", 6L),
      "numeric", "logical", "character", "list"
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
      "Semantic quantity such as sd, cor, var_prop, var_ratio, or sd_ratio, or an empty string.",
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
    stringsAsFactors = FALSE
  )
  out$arguments <- I(list())
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
    owner_type = "", owner_name = "", quantity = "",
    scale_role = "", parent_quantity_id = "",
    arguments = character(), source_type = "none", extraction_key){

  out <- .bt_parameter_catalog_empty_quantities()
  scalar_columns <- setdiff(names(out), c("arguments", "extraction_key"))
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
    source_type
  )
  out$arguments <- I(list(as.character(arguments)))
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

.bt_parameter_catalog_empty_overrides <- function(){

  data.frame(
    canonical_name = character(),
    role = character(),
    term = character(),
    component = character(),
    display_label = character(),
    display_scale = character(),
    stringsAsFactors = FALSE
  )
}

.bt_parameter_catalog_factor_label <- function(formula_parameter, term,
                                               component){

  paste0(
    .bt_random_effect_summary_formula_prefix(formula_parameter, TRUE),
    term,
    "[",
    component,
    "]"
  )
}

.bt_parameter_catalog_factor_map <- function(registry, prior_list){

  out <- list(
    direct = .bt_parameter_catalog_empty_overrides(),
    derived = .bt_parameter_catalog_empty_quantities()
  )
  if(length(prior_list) == 0L || is.null(names(prior_list))){
    return(out)
  }
  direct_rows <- list()
  derived_rows <- list()
  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    formula_factor <- .bt_formula_prior_is_factor(prior) &&
      !is.null(attr(prior, "term_components", exact = TRUE)) &&
      !isTRUE(attr(prior, "random_sd", exact = TRUE))
    if(!formula_factor){
      next
    }
    coordinates <- .JAGS_prior_factor_names(parameter, prior)
    registry_rows <- match(coordinates, registry$coordinate_name)
    if(anyNA(registry_rows) ||
       any(registry$role[registry_rows] != "fixed_coefficient")){
      stop(
        "Parameter catalog factor coordinates are missing or malformed for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    coordinate_metadata <- registry[registry_rows, , drop = FALSE]
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
    if(ncol(design) != length(coordinates) ||
       nrow(design) != length(design_info$cell_names) ||
       length(cell_names) != nrow(design) ||
       any(!is.finite(design))){
      stop(
        "Parameter catalog factor metadata disagree with registry coordinates for '",
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
    semantic_names <- .factor_contrast_parameter_names(
      parameter,
      design_info$level_names,
      design_info$cell_names
    )
    if(length(semantic_names) != nrow(design) ||
       anyNA(semantic_names) || any(!nzchar(semantic_names)) ||
       anyDuplicated(semantic_names)){
      stop(
        "Parameter catalog factor semantic names are malformed for '",
        parameter, "'. Refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    direct_cells <- rep.int(NA_integer_, length(coordinates))
    for(coordinate in seq_along(coordinates)){
      identity_row <- rep.int(0, ncol(design))
      identity_row[coordinate] <- 1
      matches <- which(vapply(seq_len(nrow(design)), function(cell){
        identical(unname(design[cell, ]), identity_row)
      }, logical(1)))
      if(length(matches) == 1L){
        direct_cells[coordinate] <- matches
        metadata <- coordinate_metadata[coordinate, , drop = FALSE]
        component <- cell_names[matches]
        direct_rows[[length(direct_rows) + 1L]] <- data.frame(
          canonical_name = coordinates[coordinate],
          role = metadata$role,
          term = metadata$term,
          component = component,
          display_label = .bt_parameter_catalog_factor_label(
            metadata$formula_parameter,
            metadata$term,
            component
          ),
          display_scale = metadata$fitted_scale,
          stringsAsFactors = FALSE
        )
      }
    }
    for(cell in seq_len(nrow(design))){
      component <- cell_names[cell]
      if(cell %in% direct_cells){
        next
      }
      nonzero <- which(design[cell, ] != 0)
      dependencies <- coordinates[nonzero]
      weights <- unname(design[cell, nonzero])
      dependency_metadata <- coordinate_metadata[nonzero, , drop = FALSE]
      structural <- length(nonzero) == 0L ||
        all(dependency_metadata$monitor_status == "structural")
      unavailable <- length(nonzero) > 0L &&
        any(dependency_metadata$monitor_status == "unavailable")
      status <- if(structural){
        "structural"
      }else if(unavailable){
        "unavailable"
      }else{
        "derived"
      }
      fixed_value <- if(structural){
        if(length(nonzero) == 0L){
          0
        }else{
          sum(weights * dependency_metadata$fixed_value)
        }
      }else{
        NA_real_
      }
      derived_rows[[length(derived_rows) + 1L]] <- .bt_parameter_catalog_quantity(
        canonical_name = semantic_names[cell],
        namespace = if(nzchar(coordinate_metadata$formula_parameter[1L])){
          coordinate_metadata$formula_parameter[1L]
        }else{
          "model"
        },
        role = "fixed_coefficient",
        formula_parameter = coordinate_metadata$formula_parameter[1L],
        term = coordinate_metadata$term[1L],
        component = component,
        display_label = .bt_parameter_catalog_factor_label(
          coordinate_metadata$formula_parameter[1L],
          coordinate_metadata$term[1L],
          component
        ),
        fitted_scale = coordinate_metadata$fitted_scale[1L],
        display_scale = coordinate_metadata$fitted_scale[1L],
        status = status,
        fixed_value = fixed_value,
        internal = FALSE,
        source_type = if(length(dependencies) <= 1L){
          "identity"
        }else{
          "composite"
        },
        extraction_key = list(
          type = "factor_level",
          dependencies = dependencies,
          weights = weights
        )
      )
    }
  }
  if(length(direct_rows) > 0L){
    out$direct <- do.call(rbind, direct_rows)
    rownames(out$direct) <- NULL
  }
  if(length(derived_rows) > 0L){
    out$derived <- do.call(rbind, derived_rows)
    rownames(out$derived) <- NULL
  }
  out
}

.bt_parameter_catalog_registry_quantities <- function(
    registry, overrides = .bt_parameter_catalog_empty_overrides(),
    suppress = character()){

  out <- .bt_parameter_catalog_empty_quantities()
  if(anyDuplicated(overrides$canonical_name)){
    stop(
      "Parameter catalog direct semantic mappings are not unique. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  keep <- registry$role != "backend_anchor" & !registry$internal &
    !registry$coordinate_name %in% suppress
  registry <- registry[keep, , drop = FALSE]
  if(nrow(registry) == 0L){
    return(out)
  }
  rows <- vector("list", nrow(registry))
  for(i in seq_len(nrow(registry))){
    row <- registry[i, , drop = FALSE]
    override <- match(row$coordinate_name, overrides$canonical_name)
    role <- row$role
    term <- row$term
    component <- row$column
    display_label <- row$display_label
    display_scale <- row$fitted_scale
    if(!is.na(override)){
      role <- overrides$role[override]
      term <- overrides$term[override]
      component <- overrides$component[override]
      display_label <- overrides$display_label[override]
      display_scale <- overrides$display_scale[override]
    }
    namespace <- if(nzchar(row$formula_parameter)){
      row$formula_parameter
    }else{
      "model"
    }
    rows[[i]] <- .bt_parameter_catalog_quantity(
      canonical_name = row$coordinate_name,
      namespace = namespace,
      role = role,
      formula_parameter = row$formula_parameter,
      term = term,
      component = component,
      display_label = display_label,
      fitted_scale = row$fitted_scale,
      display_scale = display_scale,
      status = row$monitor_status,
      fixed_value = row$fixed_value,
      internal = FALSE,
      source_type = "identity",
      extraction_key = list(
        type = "registry",
        dependencies = row$coordinate_name
      )
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.bt_parameter_catalog_aliases <- function(quantities, formula_design = NULL){

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
        if(identical(quantity$role, "fixed_coefficient") &&
           nzchar(quantity$term) && nzchar(quantity$component)){
          paste0(quantity$term, "[", quantity$component, "]")
        }else{
          character()
        }
      ))
    }
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
  owner <- .bt_random_effect_public_name(random_term)
  if(!is.character(owner) || length(owner) != 1L || is.na(owner) ||
     !nzchar(owner)){
    return(character())
  }

  if(identical(key$evaluator, "correlation")){
    components <- .bt_random_effect_summary_column_components(random_term)
    index <- key$index
    if(length(components) < 2L || !is.numeric(index) || length(index) != 1L ||
       is.na(index)){
      return(character())
    }
    pairs <- utils::combn(seq_along(components), 2L)
    if(index < 1L || index > ncol(pairs)){
      return(character())
    }
    pair <- components[pairs[, as.integer(index)]]
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
  correlations <- apply(pairs, 2L, function(pair){
    .bt_random_effect_semantic_name(
      parameter = quantity$formula_parameter,
      owner = owner,
      quantity = "cor",
      arguments = pair,
      formula_prefix = TRUE
    )
  })
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

.bt_parameter_catalog_coordinates <- function(registry, names){

  names <- unique(names[!is.na(names) & nzchar(names)])
  if(length(names) == 0L){
    return(character())
  }
  bases <- .bt_parameter_registry_base(registry$coordinate_name)
  unique(registry$coordinate_name[
    registry$coordinate_name %in% names |
      registry$monitor_name %in% names |
      bases %in% names
  ])
}

.bt_parameter_catalog_random_block_dependencies <- function(
    registry, formula_parameter, random_block,
    roles = c("random_sd", "random_correlation")){

  registry$coordinate_name[
    registry$formula_parameter == formula_parameter &
      registry$random_block == random_block &
      registry$role %in% roles
  ]
}

.bt_parameter_catalog_random_correlation_sources <- function(registry,
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
      registry,
      correlation$primitive_names
    )
    n_pairs <- random_term$n_columns * (random_term$n_columns - 1L) / 2L
    if(length(primitive_dependencies) == n_pairs){
      return(primitive_dependencies)
    }
    return(.bt_parameter_catalog_coordinates(
      registry,
      as.vector(.bt_random_effect_cholesky_names(
        random_term,
        random_term$n_columns
      ))
    ))
  }
  if(identical(correlation$type, "rho")){
    return(.bt_parameter_catalog_coordinates(
      registry,
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

.bt_parameter_catalog_random_status <- function(
    key, registry, prior_list, formula_design, formula_scale){

  dependency_rows <- match(key$dependencies, registry$coordinate_name)
  if(anyNA(dependency_rows)){
    stop(
      "Parameter catalog random-summary dependencies are missing from the registry. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }
  dependency_status <- registry$monitor_status[dependency_rows]
  if(any(dependency_status == "unavailable")){
    return(list(status = "unavailable", fixed_value = NA_real_))
  }
  if(any(dependency_status == "sampled")){
    return(list(status = "sampled", fixed_value = NA_real_))
  }

  values <- matrix(
    registry$fixed_value[dependency_rows],
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

  source_values <- c(1, 2)
  model_samples <- matrix(
    source_values,
    ncol = 1L,
    dimnames = list(NULL, source_parameter)
  )
  summary <- tryCatch(
    .bt_random_effect_summary_sd_samples(
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      parameter = parameter,
      formula_scale = formula_scale
    ),
    error = function(error) NULL
  )
  if(is.null(summary) || ncol(summary$values) < index){
    return(NA_real_)
  }
  scale <- as.numeric(summary$values[, index]) / source_values
  if(length(scale) != 2L || any(!is.finite(scale) | scale <= 0) ||
     abs(scale[1L] - scale[2L]) >
       sqrt(.Machine$double.eps) * max(1, abs(scale))){
    return(NA_real_)
  }
  unname(scale[1L])
}

.bt_parameter_catalog_random_definitions <- function(registry, prior_list,
                                                     formula_design,
                                                     formula_scale = NULL){

  out <- list(
    direct = .bt_parameter_catalog_empty_overrides(),
    derived = .bt_parameter_catalog_empty_quantities(),
    suppress = character()
  )
  if(is.null(prior_list)){
    prior_list <- list()
  }
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
                             scale_role = "", parent_quantity_id = "",
                             arguments = character(),
                              source_type, source_parameter = "",
                              source_prior = "",
                              source_transform = "identity",
                              source_scale = NA_real_){
    canonical_name <- .bt_random_effect_semantic_name(
      parameter = parameter,
      owner = owner_name,
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
        source_scale = source_scale
      ),
      metadata[names(metadata) %in% c(
        "prior_name", "allocation_label", "parent_allocation", "index"
      )]
    )
    display_label <- canonical_name
    state <- .bt_parameter_catalog_random_status(
      key = key,
      registry = registry,
      prior_list = prior_list,
      formula_design = formula_design,
      formula_scale = formula_scale
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
      parent_allocation = if(is.null(allocation$parent)) "" else
        allocation$parent$allocation,
      allocation = allocation,
      random_term = random_term
    )
    scale_role <- .bt_random_effect_allocation_scale_role(allocation)
    scale_names <- .bt_parameter_catalog_allocation_scale_names(allocation)
    scale_dependencies <- .bt_parameter_catalog_coordinates(
      registry,
      scale_names
    )
    if(length(scale_names) > 0L && length(scale_dependencies) > 0L){
      source_name <- .bt_random_sd_binding_source_name(allocation$source)
      parent_factors <- allocation$parent_factors
      if(is.null(parent_factors)){
        parent_factors <- list()
      }
      direct_scale <- length(parent_factors) == 0L
      source_prior <- if(source_name %in% names(prior_list)) source_name else ""
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
        owner_name = allocation$label,
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
        owner_name = allocation$label,
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
      out$suppress <<- unique(c(out$suppress, scale_dependencies))
    }
    dependencies <- .bt_parameter_catalog_coordinates(
      registry,
      allocation$weight_name
    )
    out$suppress <<- unique(c(out$suppress, dependencies))
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
        dependencies = dependencies,
        metadata = c(metadata, list(index = i)),
        block = block,
        term = allocation$label,
        component = components[i],
        fitted_scale = "unitless",
        display_scale = "unitless",
        owner_type = "variance_allocation",
        owner_name = allocation$label,
        quantity = allocation_type$summary,
        scale_role = scale_role,
        arguments = components[i],
        source_type = if(identical(allocation_type$summary, "var_prop")){
          "identity"
        }else{
          "one_to_one_transform"
        },
        source_parameter = allocation$weight_name,
        source_prior = allocation$weight_name,
        source_transform = allocation_type$summary
      )
      if(identical(.bt_random_effect_summary_allocation_target(allocation),
                   "sd_component")){
        multiplier_name <- .bt_random_effect_summary_name(
          parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
          type = "sd_ratio",
          parts = c(allocation$label, components[i])
        )
        add_definition(
          raw_name = multiplier_name,
          role = "random_sd_ratio",
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
          owner_name = allocation$label,
          quantity = "sd_ratio",
          scale_role = scale_role,
          arguments = components[i],
          source_type = "one_to_one_transform",
          source_parameter = allocation$weight_name,
          source_prior = allocation$weight_name,
          source_transform = "sd_ratio"
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
        inclusion_dependencies <- .bt_parameter_catalog_coordinates(
          registry,
          inclusion[[component_label]]$indicator_name
        )
        add_definition(
          raw_name = raw_name,
          role = "random_inclusion",
          parameter = parameter,
          label = paste0("inclusion(", allocation$label, ": ",
                         component_label, ")"),
          evaluator = "allocation_inclusion",
          dependencies = inclusion_dependencies,
          metadata = c(metadata, list(index = inclusion_i)),
          block = block,
          term = allocation$label,
          component = component_label,
          fitted_scale = "unitless",
          display_scale = "unitless",
          owner_type = "variance_allocation",
          owner_name = allocation$label,
          quantity = "inclusion",
          arguments = component_label,
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
        direct_source <- direct_sd && !allocation_sd
        allocation_dependencies <- .bt_parameter_catalog_coordinates(
          registry,
          unlist(lapply(
            random_term$sd_binding$allocations,
            .bt_parameter_catalog_allocation_names
          ), use.names = FALSE)
        )
        correlation_source_dependencies <-
          .bt_parameter_catalog_random_correlation_sources(
            registry,
            random_term
          )
        sd_dependencies <- if(allocation_sd){
          unique(c(
            allocation_dependencies,
            if(!direct_sd) correlation_source_dependencies else character()
          ))
        }else{
          unique(c(
            .bt_parameter_catalog_coordinates(registry, sd_names),
            if(!direct_sd) correlation_source_dependencies else character()
          ))
        }
        for(i in seq_along(sd_names)){
          sd_quantity <- .bt_random_effect_semantic_sd_quantity(random_term)
          sd_arguments <- .bt_random_effect_semantic_sd_arguments(
            components[i]
          )
          raw_name <- .bt_random_effect_summary_name(
            parameter = parameter,
            type = "sd",
            parts = c(block, components[i])
          )
          label <- .bt_random_effect_sd_summary_label(
            component = components[i],
            group = group,
            random_term = random_term
          )
          out$suppress <- unique(c(
            out$suppress,
            intersect(sd_names[i], registry$coordinate_name)
          ))
          dependencies <- if(direct_source){
            .bt_parameter_catalog_coordinates(registry, sd_names)
          }else{
            sd_dependencies
          }
          source_coordinate <- .bt_parameter_catalog_coordinates(
            registry,
            sd_names[i]
          )
          source_type <- if(direct_source){
            "identity"
          }else if(length(dependencies) == 1L){
            if(allocation_sd) "composite" else "one_to_one_transform"
          }else{
            "composite"
          }
          source_prior <- .bt_parameter_registry_base(sd_names[i])
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
            role = if(identical(sd_quantity, "sd")){
              "random_sd"
            }else{
              "random_sd_ratio"
            },
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
            quantity = sd_quantity,
            arguments = sd_arguments,
            source_type = source_type,
            source_parameter = source_parameter,
            source_prior = source_prior,
            source_transform = source_transform,
            source_scale = source_scale
          )
          if(identical(sd_quantity, "sd_ratio")){
            var_name <- .bt_random_effect_summary_name(
              parameter = parameter,
              type = "var_ratio",
              parts = c(block, components[i])
            )
            add_definition(
              raw_name = var_name,
              role = "random_var_ratio",
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
              quantity = "var_ratio",
              arguments = sd_arguments,
              source_type = if(direct_source){
                "one_to_one_transform"
              }else{
                "composite"
              },
              source_parameter = if(direct_source) source_parameter else "",
              source_prior = source_prior,
              source_transform = "square"
            )
          }
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
              registry,
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
          registry,
          c(correlation$rho_name, correlation$sample_name)
        )
        out$suppress <- unique(c(out$suppress, dependencies))
        source_prior <- correlation$prior_name
        if(is.null(source_prior) || length(source_prior) != 1L ||
           is.na(source_prior) || !source_prior %in% names(prior_list)){
          source_prior <- .bt_parameter_registry_base(source_parameter)
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
          registry,
          formula_parameter = parameter,
          random_block = block,
          roles = c("random_sd", "random_correlation")
        )
        correlation_dependencies <-
          .bt_parameter_catalog_random_correlation_sources(
            registry,
            random_term
          )
        if(scaled_correlation){
          correlation_dependencies <- unique(c(
            correlation_dependencies,
            .bt_parameter_catalog_coordinates(registry, sd_names),
            .bt_parameter_catalog_coordinates(
              registry,
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

  if(length(rows) > 0L){
    out$derived <- do.call(rbind, rows)
    rownames(out$derived) <- NULL
    allocation_rows <- which(
      out$derived$owner_type == "variance_allocation"
    )
    allocation_owners <- unique(out$derived$owner_name[allocation_rows])
    allocation_sd_ids <- stats::setNames(rep("", length(allocation_owners)),
                                          allocation_owners)
    for(owner in allocation_owners){
      candidates <- allocation_rows[
        out$derived$owner_name[allocation_rows] == owner &
          out$derived$quantity[allocation_rows] %in% c("sd_total", "sd_common")
      ]
      if(length(candidates) == 1L){
        allocation_sd_ids[[owner]] <- out$derived$quantity_id[candidates]
      }
    }
    for(i in allocation_rows){
      owner <- out$derived$owner_name[i]
      if(!out$derived$quantity[i] %in% c("sd_total", "sd_common")){
        out$derived$parent_quantity_id[i] <- allocation_sd_ids[[owner]]
        next
      }
      key <- out$derived$extraction_key[[i]]
      parent <- key$parent_allocation
      if(is.character(parent) && length(parent) == 1L && nzchar(parent) &&
         parent %in% names(allocation_sd_ids)){
        out$derived$parent_quantity_id[i] <- allocation_sd_ids[[parent]]
      }
    }
  }
  out
}

.bt_build_parameter_catalog <- function(registry, prior_list = NULL,
                                        formula_design = NULL,
                                        formula_scale = NULL){

  .bt_validate_parameter_registry(registry)
  factor_map <- .bt_parameter_catalog_factor_map(
    registry = registry,
    prior_list = prior_list
  )
  random_map <- .bt_parameter_catalog_random_definitions(
    registry = registry,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  overrides <- rbind(factor_map$direct, random_map$direct)
  base <- .bt_parameter_catalog_registry_quantities(
    registry = registry,
    overrides = overrides,
    suppress = random_map$suppress
  )
  quantities <- rbind(base, factor_map$derived, random_map$derived)
  rownames(quantities) <- NULL
  aliases <- .bt_parameter_catalog_aliases(quantities, formula_design)
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

.bt_parameter_catalog_valid_native_key <- function(key, quantity){

  if(identical(key$type, "registry")){
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
    c("arguments", "fixed_value", "internal", "extraction_key")
  )
  if(!all(vapply(quantities[character_columns], is.character, logical(1))) ||
     !is.list(quantities$arguments) ||
     !all(vapply(quantities$arguments, function(arguments){
       is.character(arguments) && !anyNA(arguments)
     }, logical(1))) ||
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
    "display_label", "fitted_scale", "display_scale", "status",
    "source_type"
  )
  if(any(!nzchar(as.matrix(quantities[required_nonempty]))) ||
     any(!nzchar(aliases$alias)) || any(!nzchar(aliases$quantity_id)) ||
     any(!nzchar(aliases$namespace)) ||
     anyDuplicated(quantities$quantity_id) ||
     any(!quantities$source_type %in%
           c("identity", "one_to_one_transform", "composite", "none")) ||
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
      !nzchar(quantities$owner_name) |
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
