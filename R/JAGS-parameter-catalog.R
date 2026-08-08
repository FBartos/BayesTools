# Metadata-only semantic parameter catalog and deferred draw extraction.

.bt_parameter_catalog_version <- 3L
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
#' Identity random-effect summaries reuse their sampled or structural registry
#' coordinate. Summaries requiring a scale or covariance transformation expose
#' only the transformed quantity; their fitted-scale inputs and all other
#' private implementation coordinates remain dependencies rather than public
#' catalog rows.
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
    registry_rows <- match(coordinates, registry$canonical_name)
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
    !registry$canonical_name %in% suppress
  registry <- registry[keep, , drop = FALSE]
  if(nrow(registry) == 0L){
    return(out)
  }
  rows <- vector("list", nrow(registry))
  for(i in seq_len(nrow(registry))){
    row <- registry[i, , drop = FALSE]
    override <- match(row$canonical_name, overrides$canonical_name)
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
      canonical_name = row$canonical_name,
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
      if(identical(quantity$role, "fixed_coefficient") &&
         nzchar(quantity$term) && nzchar(quantity$component)){
        paste0(quantity$term, "[", quantity$component, "]")
      }else{
        character()
      }
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

.bt_parameter_catalog_coordinates <- function(registry, names){

  names <- unique(names[!is.na(names) & nzchar(names)])
  if(length(names) == 0L){
    return(character())
  }
  bases <- .bt_parameter_registry_base(registry$canonical_name)
  unique(registry$canonical_name[
    registry$canonical_name %in% names |
      registry$monitor_name %in% names |
      bases %in% names
  ])
}

.bt_parameter_catalog_random_block_dependencies <- function(
    registry, formula_parameter, random_block,
    roles = c("random_sd", "random_correlation")){

  registry$canonical_name[
    registry$formula_parameter == formula_parameter &
      registry$random_block == random_block &
      registry$role %in% roles
  ]
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
    allocation$total_name,
    allocation$source_node,
    source_name,
    factor_names,
    inclusion_names
  ))
}

.bt_parameter_catalog_random_status <- function(
    key, registry, prior_list, formula_design, formula_scale){

  dependency_rows <- match(key$dependencies, registry$canonical_name)
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
    return(list(status = "derived", fixed_value = NA_real_))
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
  direct_rows <- list()
  used_names <- registry$canonical_name
  add_direct <- function(source_name, role, parameter, label,
                         term = "", component = "",
                         display_scale = "original"){
    source <- match(source_name, registry$canonical_name)
    if(is.na(source) || isTRUE(registry$internal[source])){
      return(FALSE)
    }
    direct_rows[[length(direct_rows) + 1L]] <<- data.frame(
      canonical_name = source_name,
      role = role,
      term = term,
      component = component,
      display_label = paste0(
        .bt_random_effect_summary_formula_prefix(parameter, TRUE),
        label
      ),
      display_scale = display_scale,
      stringsAsFactors = FALSE
    )
    TRUE
  }
  add_definition <- function(raw_name, role, parameter, label,
                             evaluator, dependencies, metadata = NULL,
                             block = "",
                             term = "", component = "",
                             fitted_scale = "fitted_covariance",
                             display_scale = "original"){
    canonical_name <- .bt_random_effect_summary_unique_name(
      raw_name,
      used_names
    )
    used_names <<- c(used_names, canonical_name)
    namespace <- if(nzchar(parameter)) parameter else "model"
    dependencies <- unique(dependencies)
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
    label <- paste0("sd_total(", allocation, ")")
    if(!add_direct(
      source_name = prior_name,
      role = "random_sd_total",
      parameter = parameter,
      label = label,
      term = allocation,
      component = allocation
    )){
      add_definition(
        raw_name = raw_name,
        role = "random_sd_total",
        parameter = parameter,
        label = label,
        evaluator = "sd_total",
        dependencies = .bt_parameter_catalog_coordinates(
          registry,
          prior_name
        ),
        metadata = list(prior_name = prior_name),
        term = allocation,
        component = allocation
      )
    }
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
    dependencies <- .bt_parameter_catalog_coordinates(
      registry,
      allocation$weight_name
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
        dependencies = dependencies,
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
          dependencies = dependencies,
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
        direct_sd <- .bt_parameter_catalog_random_sd_is_direct(
          random_term = random_term,
          parameter = parameter,
          formula_scale = formula_scale
        )
        sd_dependencies <- unique(c(
          .bt_parameter_catalog_coordinates(registry, sd_names),
          .bt_parameter_catalog_random_block_dependencies(
            registry,
            formula_parameter = parameter,
            random_block = block
          ),
          .bt_parameter_catalog_coordinates(
            registry,
            unlist(lapply(
              random_term$sd_binding$allocations,
              .bt_parameter_catalog_allocation_names
            ), use.names = FALSE)
          )
        ))
        for(i in seq_along(sd_names)){
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
          if(!direct_sd || !add_direct(
            source_name = sd_names[i],
            role = "random_sd",
            parameter = parameter,
            label = label,
            term = block,
            component = components[i]
          )){
            out$suppress <- unique(c(
              out$suppress,
              intersect(sd_names[i], registry$canonical_name)
            ))
            add_definition(
              raw_name = raw_name,
              role = "random_sd",
              parameter = parameter,
              label = label,
              evaluator = "sd",
              dependencies = sd_dependencies,
              metadata = c(term_metadata, list(index = i)),
              block = block,
              term = block,
              component = components[i]
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
        label <- paste0("rho(", group, ")")
        rho_scale <- .bt_random_effect_rho_scale_metadata(
          correlation,
          random_term,
          context = "Parameter catalog"
        )
        direct_rho <- identical(rho_scale, "rho") && add_direct(
          source_name = correlation$rho_name,
          role = "random_correlation",
          parameter = parameter,
          label = label,
          term = block,
          display_scale = "unitless"
        )
        if(!direct_rho){
          out$suppress <- unique(c(
            out$suppress,
            intersect(
              c(correlation$rho_name, correlation$sample_name),
              registry$canonical_name[!registry$internal]
            )
          ))
          add_definition(
            raw_name = raw_name,
            role = "random_correlation",
            parameter = parameter,
            label = label,
            evaluator = "rho",
            dependencies = .bt_parameter_catalog_coordinates(
              registry,
              c(correlation$rho_name, correlation$sample_name)
            ),
            metadata = term_metadata,
            block = block,
            term = block,
            fitted_scale = "fitted_covariance",
            display_scale = "unitless"
          )
        }
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
        correlation_roles <- if(scaled_correlation){
          c("random_sd", "random_correlation")
        }else{
          "random_correlation"
        }
        correlation_dependencies <- .bt_parameter_catalog_random_block_dependencies(
          registry,
          formula_parameter = parameter,
          random_block = block,
          roles = correlation_roles
        )
        if(scaled_correlation){
          correlation_dependencies <- unique(c(
            correlation_dependencies,
            .bt_parameter_catalog_coordinates(
              registry,
              unlist(lapply(
                random_term$sd_binding$allocations,
                .bt_parameter_catalog_allocation_names
              ), use.names = FALSE)
            )
          ))
        }
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

  if(length(direct_rows) > 0L){
    out$direct <- do.call(rbind, direct_rows)
    rownames(out$direct) <- NULL
  }
  if(length(rows) > 0L){
    out$derived <- do.call(rbind, rows)
    rownames(out$derived) <- NULL
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
    "sd_total", "sd", "inclusion", "rho", "correlation",
    "allocation", "allocation_inclusion"
  )
  if(!scalar_character(key$evaluator) ||
     !key$evaluator %in% evaluators ||
     !scalar_character(key$formula_parameter, allow_empty = TRUE) ||
     !scalar_character(key$random_block, allow_empty = TRUE) ||
     !scalar_character(key$summary_name)){
    return(FALSE)
  }
  if(identical(key$evaluator, "sd_total")){
    return(scalar_character(key$prior_name))
  }
  if(key$evaluator %in% c("sd", "inclusion", "correlation")){
    return(
      is.numeric(key$index) && length(key$index) == 1L &&
        !is.na(key$index) && key$index == as.integer(key$index) &&
        key$index >= 1L
    )
  }
  if(key$evaluator %in% c("allocation", "allocation_inclusion")){
    return(
      scalar_character(key$allocation_label) &&
        is.numeric(key$index) && length(key$index) == 1L &&
        !is.na(key$index) && key$index == as.integer(key$index) &&
        key$index >= 1L
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
  if(!is.null(random_term)){
    model_samples <- .bt_parameter_catalog_complete_cholesky_triangle(
      random_term,
      model_samples
    )
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
