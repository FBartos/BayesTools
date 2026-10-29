# ============================================================================ #
# hypothesis-model-provider.R
# ============================================================================ #
#
# The provider of the hypothesis engine (hypothesis_BF_model()): one object of
# callbacks that tells the engine what the tested model offers. The method for
# a BayesTools_fit is the default provider; a downstream model class starts
# from it and replaces the callbacks that differ. See
# .agents/instructions/hypothesis.md.
#
# ============================================================================ #


#' @title Hypothesis provider of a fitted model
#'
#' @description \code{hypothesis_provider()} is the single hook through which
#' [hypothesis_BF_model()] reaches a model. Its method returns a provider, a
#' list of the callbacks and values that the engine needs for the statements
#' it plans: the parameter catalog and its entries, the draws of the tested
#' quantities, the random-effect quantities, and the density backends that can
#' evaluate them. The method for a \code{BayesTools_fit} returns the default
#' provider, which implements every callback with the functions of BayesTools
#' ([parameter_catalog()], [as_mixed_posteriors()], [marginal_posterior()],
#' [parameter_mixed_posterior()], [JAGS_formula_coefficient_transform()], and
#' [JAGS_formula_prior_density()]). A downstream model class that wraps a
#' fit starts from \code{hypothesis_provider(<its fit>)} and replaces the
#' callbacks that differ with [utils::modifyList()], so that it overrides only
#' what is specific to its model.
#'
#' @details A provider is a list of class
#' \code{BayesTools_hypothesis_provider}; the engine checks it before use. Its
#' elements are the following. The callbacks close over the model object.
#' \describe{
#'   \item{\code{kind}}{\code{"model"} (default) or \code{"marginal_means"}.}
#'   \item{\code{object}}{the model object. It is hashed for the cache of
#'   plans and passed to the density backends.}
#'   \item{\code{fit}}{its \code{BayesTools_fit}, from which the transforms of
#'   formula coefficients and their prior densities are read.}
#'   \item{\code{capability_object}}{the object the density backends receive
#'   in \code{available()} and in their requests (default \code{object}).}
#'   \item{\code{build}, \code{refit_class}}{the name of the build in the
#'   requests to refit a model with outdated metadata ("Refit the model with
#'   the current <build> build."), and classes that precede
#'   \code{BayesTools_refit_required} in those conditions.}
#'   \item{\code{caller}, \code{random_point_mass_advice},
#'   \code{transform_scaled_default}}{the name of the entry point in the
#'   reasons of refusals (default \code{"hypothesis_BF()"}); the sentence that
#'   ends the reason of a refused point null at a prior point mass of a
#'   random-effect quantity; and the scale of [hypothesis_quantities()] when
#'   \code{transform_scaled} is \code{NULL}.}
#'   \item{\code{ordinate_label}}{how the messages of refusals name the
#'   ordinate backends (default \code{"ordinate density methods"}).}
#'   \item{\code{component_argument}}{the name of the argument through which
#'   the caller selects a value of \code{components} (default
#'   \code{"namespace"}, the argument of [hypothesis_BF_model()] and of the
#'   method of [hypothesis_BF()] for a \code{BayesTools_fit}). The messages of
#'   the engine that name that argument ("The 'hypothesis' argument selects
#'   namespace = 'mu' but 'namespace' was set to 'random'.") use it; a model
#'   class that passes the selection under another name sets it, for example
#'   \code{"component"}.}
#'   \item{\code{density_methods}}{the names that \code{density_method} can
#'   match, in order (default: the registered backends that are not of type
#'   \code{"stored"}).}
#'   \item{\code{components}}{the selection vocabulary of the \code{namespace}
#'   argument as data: \code{values} (the accepted values, in the order of the
#'   error message), \code{synonyms} (a named character vector, for example
#'   \code{c(location = "mods")}), \code{default} (\code{"auto"}),
#'   \code{refused} (a named character vector of the reasons for which a
#'   component cannot be tested), \code{random} (the value that selects the
#'   random-effect quantities), and \code{advice} (the sentence that follows
#'   the message of an ambiguous statement).}
#'   \item{\code{metadata()}}{a list with the \code{catalog} (a
#'   [parameter_catalog()], possibly extended with
#'   [parameter_catalog_extend()]; the engine parses and resolves with it) and
#'   the \code{entries}: a data frame with one row per tested parameter and
#'   the columns \code{quantity_id}, \code{parameter} (the name of the tested
#'   parameter), \code{component} (a value of \code{components$values}),
#'   \code{term}, \code{source}, \code{formula_parameter}, \code{role}
#'   (\code{"formula_coefficient_group"} for a whole factor term, whose levels
#'   are its \code{member_quantity_ids}), \code{quantity}, \code{status},
#'   \code{fixed_value}, \code{aliases} (a list column), and
#'   \code{member_quantity_ids} (a list column). The alias rows of the catalog
#'   of the entries carry the component as their \code{namespace} and
#'   \code{component}.}
#'   \item{\code{uncovered(resolution, quantity_ids, symbol, metadata)}}{the
#'   refusal for resolved quantities that no entry covers: a list with
#'   \code{reason} and either \code{type} (\code{"target"} or
#'   \code{"statement"}) or the complete condition \code{class}.}
#'   \item{\code{sample_parameters(parameter)}, \code{mixed_posteriors(parameters,
#'   conditional, conditional_rule, transform_scaled, n_prior_samples)},
#'   \code{density_parameter(samples, parameter, sample_parameters)},
#'   \code{formula_name_map(formula_parameter)}}{the draws of a fixed
#'   coefficient: the prior-list names whose draws are mixed (the parameter and
#'   the owners of the fitted coordinates that its prior density combines), the
#'   mixed posteriors ([as_mixed_posteriors()]), the element of the samples
#'   that is the tested parameter, and the map of the formula coefficients to
#'   the fitted names.}
#'   \item{\code{random}}{the callbacks for the random-effect quantities:
#'   \code{select(parameter, transform_scaled)} (a list with the element
#'   \code{spec}, a list with the \code{label} and \code{status} of the
#'   quantity, and any elements the other callbacks need),
#'   \code{mixed_posterior(parameter, transform_scaled, conditional, selected)}
#'   (the mixed posterior as the one-element list returned by
#'   [as_mixed_posteriors()]), \code{support(selected)} (the exact bounds on
#'   the scale requested from \code{select()}, which are those of the draws),
#'   \code{density_target(parameter, transform_scaled, operation)} (a list with the
#'   quantity that a density backend evaluates, its \code{parameter_spec}, and
#'   \code{display_transform}, or the refusal of a quantity without one: its
#'   \code{reason}, and optionally its \code{type} (\code{"target"}, the
#'   default, \code{"method"}, or \code{"fixed"}: the engine's refusal
#'   \code{BayesTools_hypothesis_<type>} and its parent
#'   \code{BayesTools_hypothesis_unavailable}) and its \code{cause}
#'   (\code{"random_target"}, the default, or \code{"original_scale"}), for
#'   which the density backend appends its classes
#'   (\code{refusal_classes(cause)})), \code{sd_pair(selected)} (the standard deviation
#'   of which a variance is the square, or \code{NULL}),
#'   \code{zero_boundary(selected)} (\code{NULL}, or a list with the
#'   statement that tests the omission of a component),
#'   \code{defined_share(values, selected, conditional)},
#'   \code{defined_footnote(label, samples, posterior_defined, prior_defined)},
#'   \code{defined_draws(values, samples, label)} (the treatment of the draws
#'   at which a quantity is undefined), and \code{prior_draws(parameter,
#'   selected, transform_scaled, n_samples, seed)} (prior draws of a quantity
#'   without a prior density with provenance).}
#'   \item{\code{point_mass_reason()}}{\code{NULL}, or the text that replaces
#'   the reason of a point hypothesis at a prior point mass.}
#'   \item{\code{full_ensemble(parameter)}}{whether omitting
#'   \code{conditional} tests the parameter within the full ensemble,
#'   including the models in which it is fixed by its null component (the
#'   engine then warns); \code{combination_hint(plan)}: \code{NULL}, or text
#'   appended to the refusal of a linear combination of the levels of such a
#'   term.}
#'   \item{\code{check_plans(plans)}}{\code{NULL}, or a function that
#'   validates the plans of all statements of a request and stops with a
#'   refusal.}
#'   \item{\code{check_density_control(density_method, density_control)}}{
#'   \code{NULL}, or a function that checks the \code{density_control} of a
#'   call before the \code{control()} of the density backend sees it. It
#'   receives the name of the selected backend and the settings as supplied (a
#'   list, possibly empty, or \code{NULL}) wherever \code{control()} is
#'   called: on every call with an ordinate backend, and on every call that
#'   supplies settings to another backend (for stored marginal means also
#'   where the stored settings are normalized). It stops with the condition
#'   that refuses the settings; its value is ignored. Without it the
#'   built-in backends accept \code{NULL} and an empty list and refuse other
#'   settings with the message \code{'density_control' is not used by
#'   density_method = 'KDE'.}, and the \code{control()} of a backend
#'   validates its own settings.}
#' }
#' A provider of \code{kind = "marginal_means"} describes the stored marginal
#' means of a model instead: \code{inference} (the marginal inference with the
#' elements \code{averaged} and \code{conditional}), \code{model_averaged}
#' (logical; model-averaged point statements are tested on the conditional
#' draws and cannot mix with region statements), \code{capability_object} (the
#' fit that the backends compute ordinates from, or \code{NULL}),
#' \code{select(hypothesis, parameter)} (the marginal-means parameter and its
#' aliases that a statement names), \code{density_method_default()} and
#' \code{density_control_default(density_method)} (used when the arguments are
#' \code{NULL}), and \code{metadata()} for [hypothesis_quantities()].
#'
#' @param object a model object.
#' @param ... further arguments for methods.
#'
#' @return A provider, a list of class \code{BayesTools_hypothesis_provider}.
#'
#' @seealso [hypothesis_BF_model()], [hypothesis_density_backend()]
#' @rdname hypothesis_provider
#' @export
hypothesis_provider <- function(object, ...){

  UseMethod("hypothesis_provider")
}


#' @rdname hypothesis_provider
#' @exportS3Method hypothesis_provider default
hypothesis_provider.default <- function(object, ...){

  stop(
    "No hypothesis provider is available for objects of class '",
    class(object)[[1L]], "'. Define a 'hypothesis_provider()' method for it.",
    call. = FALSE
  )
}


#' @rdname hypothesis_provider
#' @exportS3Method hypothesis_provider BayesTools_fit
hypothesis_provider.BayesTools_fit <- function(object, ...){

  if(...length() > 0L){
    stop("Unused argument(s) in hypothesis_provider(): ",
         paste0("'", names(list(...)), "'", collapse = ", "), ".", call. = FALSE)
  }
  fit      <- object
  metadata <- function() .hypothesis_model_default_metadata(fit)
  entries  <- metadata()[["entries"]]
  values   <- unique(entries[["component"]])
  values   <- c("auto", setdiff(values, "random"), if("random" %in% values) "random")

  provider <- list(
    kind              = "model",
    object            = fit,
    fit               = fit,
    capability_object = fit,
    build             = "BayesTools",
    refit_class       = NULL,
    caller            = "hypothesis_BF()",
    component_argument = "namespace",
    transform_scaled_default = FALSE,
    ordinate_label    = "ordinate density methods",
    random_point_mass_advice = paste0(
      "Use a region or directional hypothesis, or 'conditional = TRUE' to ",
      "test the quantity within its inclusion event."
    ),
    density_methods   = NULL,
    components        = list(
      values   = values,
      synonyms = character(),
      default  = "auto",
      refused  = character(),
      random   = "random",
      advice   = paste0(
        "Set 'namespace' to ",
        paste0("'", setdiff(values, "auto"), "'", collapse = ", "), "."
      )
    ),
    metadata          = metadata,
    uncovered         = function(resolution, quantity_ids, symbol, metadata){
      .hypothesis_model_default_uncovered(
        resolution, quantity_ids, symbol, metadata
      )
    },
    sample_parameters = function(parameter) parameter,
    mixed_posteriors  = function(parameters, conditional, conditional_rule,
                                 transform_scaled, n_prior_samples){
      .hypothesis_model_default_mixed_posteriors(
        fit, parameters, conditional, conditional_rule, transform_scaled,
        n_prior_samples
      )
    },
    density_parameter = function(samples, parameter, sample_parameters){
      if(!is.null(samples[[parameter]])) parameter else sample_parameters[[1L]]
    },
    formula_name_map  = function(formula_parameter){
      JAGS_formula_name_map(fit, formula_parameter)
    },
    random            = .hypothesis_model_default_random(fit, metadata),
    point_mass_reason = NULL,
    full_ensemble     = function(parameter){
      .hypothesis_model_default_full_ensemble(fit, parameter)
    },
    combination_hint  = function(plan){
      if(isTRUE(plan[["conditional"]]) ||
         !.hypothesis_model_default_full_ensemble(fit, plan[["parameter"]])){
        return(NULL)
      }
      paste0(
        " The levels of '", plan[["label"]], "' share the null component of ",
        "the prior; test the combination within the models that include the ",
        "term with conditional = TRUE."
      )
    },
    check_plans       = NULL,
    check_density_control = NULL
  )
  class(provider) <- c("BayesTools_hypothesis_provider", "list")

  provider
}


# Checks the provider that a method returns before the engine uses it.
.hypothesis_model_validate_provider <- function(provider){

  fail <- function(...){
    stop("Invalid hypothesis provider: ", ..., call. = FALSE)
  }
  if(!inherits(provider, "BayesTools_hypothesis_provider") || !is.list(provider)){
    fail("'hypothesis_provider()' must return an object of class ",
         "'BayesTools_hypothesis_provider'.")
  }
  need_functions <- function(x, names, where = "provider"){
    missing <- names[!vapply(names, function(name) is.function(x[[name]]), logical(1))]
    if(length(missing) > 0L){
      fail("the ", where, " lacks the function(s) ",
           paste0("'", missing, "'", collapse = ", "), ".")
    }
  }
  kind <- provider[["kind"]]
  if(!is.character(kind) || length(kind) != 1L ||
     !kind %in% c("model", "marginal_means")){
    fail("'kind' must be 'model' or 'marginal_means'.")
  }
  if(is.null(provider[["object"]])){
    fail("'object' is missing.")
  }
  for(name in c("build", "caller", "ordinate_label")){
    value <- provider[[name]]
    if(!is.character(value) || length(value) != 1L || is.na(value)){
      fail("'", name, "' must be a string.")
    }
  }
  if(!is.null(provider[["component_argument"]]) &&
     (!is.character(provider[["component_argument"]]) ||
      length(provider[["component_argument"]]) != 1L ||
      is.na(provider[["component_argument"]]) ||
      !nzchar(provider[["component_argument"]]))){
    fail("'component_argument' must be a string or NULL.")
  }
  if(!is.null(provider[["density_methods"]]) &&
     (!is.character(provider[["density_methods"]]) ||
      length(provider[["density_methods"]]) == 0L)){
    fail("'density_methods' must be NULL or a character vector.")
  }
  for(name in c("point_mass_reason", "full_ensemble", "combination_hint",
                "check_plans", "check_density_control")){
    if(!is.null(provider[[name]]) && !is.function(provider[[name]])){
      fail("'", name, "' must be a function or NULL.")
    }
  }
  if(identical(kind, "marginal_means")){
    need_functions(provider, c("select", "density_method_default",
                               "density_control_default", "metadata"))
    if(is.null(provider[["capability_object"]]) &&
       (!is.character(provider[["source_missing_reason"]]) ||
        length(provider[["source_missing_reason"]]) != 1L)){
      fail("'source_missing_reason' is needed without a 'capability_object'.")
    }
    if(!is.logical(provider[["model_averaged"]]) ||
       length(provider[["model_averaged"]]) != 1L ||
       is.na(provider[["model_averaged"]])){
      fail("'model_averaged' must be TRUE or FALSE.")
    }
    if(!is.list(provider[["inference"]])){
      fail("'inference' must be the stored marginal inference.")
    }
    return(invisible(provider))
  }

  components <- provider[["components"]]
  if(!is.list(components) || !is.character(components[["values"]]) ||
     !is.character(components[["default"]]) ||
     length(components[["default"]]) != 1L ||
     !components[["default"]] %in% components[["values"]] ||
     !is.character(components[["random"]]) ||
     length(components[["random"]]) != 1L ||
     !is.character(components[["advice"]]) ||
     length(components[["advice"]]) != 1L){
    fail("'components' needs the elements 'values', 'default', 'random', ",
         "and 'advice' (and optionally 'synonyms' and 'refused').")
  }
  if(!inherits(provider[["fit"]], "BayesTools_fit")){
    fail("'fit' must be a 'BayesTools_fit'.")
  }
  if(!is.character(provider[["random_point_mass_advice"]]) ||
     length(provider[["random_point_mass_advice"]]) != 1L){
    fail("'random_point_mass_advice' must be a string.")
  }
  need_functions(provider, c(
    "metadata", "uncovered", "sample_parameters", "mixed_posteriors",
    "density_parameter", "formula_name_map"
  ))
  need_functions(provider[["random"]], c(
    "select", "mixed_posterior", "support", "density_target", "sd_pair",
    "zero_boundary", "defined_share", "defined_footnote", "defined_draws",
    "prior_draws"
  ), where = "'random' element")

  invisible(provider)
}


# ---------------------------------------------------------------------------- #
# The default provider: metadata
# ---------------------------------------------------------------------------- #

# The catalog and the entries of the tested parameters of a fit: the fixed
# coefficients (a whole factor term is one entry whose members are its level
# cells), the scalar parameters with a prior of their own, and the supported
# random-effect quantities. The alias rows of the entries carry their
# component as namespace and component, so that the component of a request
# filters the resolution.
.hypothesis_model_default_metadata <- function(fit){

  parameter_map_cache(
    map      = parameter_map(fit),
    provider = "BayesTools.hypothesis",
    key      = list(version = 1L),
    compute  = function() .hypothesis_model_default_metadata_compute(fit)
  )
}


.hypothesis_model_random_quantities <- function(){

  c("sd", "var", "sd_total", "var_total", "sd_common", "var_common", "cor",
    "var_prop", "sd_mult", "var_mult")
}


.hypothesis_model_default_metadata_compute <- function(fit){

  catalog    <- parameter_catalog(fit)
  quantities <- catalog[["quantities"]]
  public     <- !quantities[["internal"]]
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  entries    <- list()
  aliases    <- list()
  extension_quantities <- list()

  native_aliases <- function(quantity_ids){
    rows <- catalog[["aliases"]][
      catalog[["aliases"]][["quantity_id"]] %in% quantity_ids, , drop = FALSE
    ]
    unique(rows[, c("alias", "simplified"), drop = FALSE])
  }
  add_entry <- function(quantity, parameter, component, term, source,
                        formula_parameter, alias_rows,
                        member_quantity_ids = character(),
                        entry_aliases = NULL){
    keep <- !is.na(alias_rows[["alias"]]) & nzchar(alias_rows[["alias"]])
    alias_rows <- unique(alias_rows[keep, c("alias", "simplified"), drop = FALSE])
    entry <- data.frame(
      quantity_id       = quantity[["quantity_id"]],
      parameter         = parameter,
      component         = component,
      term              = term,
      source            = source,
      formula_parameter = formula_parameter,
      role              = quantity[["role"]],
      quantity          = quantity[["quantity"]],
      status            = quantity[["status"]],
      fixed_value       = quantity[["fixed_value"]],
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )
    entry[["aliases"]] <- I(list(
      if(is.null(entry_aliases)) unique(alias_rows[["alias"]]) else entry_aliases
    ))
    entry[["member_quantity_ids"]] <- I(list(unique(member_quantity_ids)))
    entries[[length(entries) + 1L]] <<- entry
    aliases[[length(aliases) + 1L]] <<- data.frame(
      alias       = alias_rows[["alias"]],
      quantity_id = rep(quantity[["quantity_id"]], nrow(alias_rows)),
      namespace   = rep(component, nrow(alias_rows)),
      component   = rep(component, nrow(alias_rows)),
      simplified  = alias_rows[["simplified"]],
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  # fixed coefficients of the formulas
  designs <- JAGS_formula_design(fit)
  for(formula_parameter in names(designs)){
    name_map   <- JAGS_formula_name_map(fit, formula_parameter)
    fixed_rows <- name_map[
      name_map[["kind"]] == "fixed" & name_map[["role"]] == "coefficient", ,
      drop = FALSE
    ]
    for(row in seq_len(nrow(fixed_rows))){
      map_row <- fixed_rows[row, , drop = FALSE]
      term    <- map_row[["term"]]
      coordinate_rows <- which(
        public &
          quantities[["formula_parameter"]] == formula_parameter &
          quantities[["role"]] == "fixed_coefficient" &
          quantities[["term"]] == term &
          !.hypothesis_model_contrast_coefficient(quantities)
      )
      if(length(coordinate_rows) == 0L){
        next
      }
      if(length(coordinate_rows) == 1L){
        quantity <- as.list(quantities[coordinate_rows, , drop = FALSE])
        member_quantity_ids <- character()
        alias_rows <- native_aliases(quantity[["quantity_id"]])
      }else{
        quantity <- .hypothesis_model_default_group_quantity(
          catalog     = catalog,
          map_row     = map_row,
          coordinates = quantities[coordinate_rows, , drop = FALSE],
          component   = formula_parameter
        )
        extension_quantities[[length(extension_quantities) + 1L]] <- quantity
        quantity <- as.list(quantity)
        member_quantity_ids <- quantities[["quantity_id"]][coordinate_rows]
        alias_rows <- data.frame(
          alias = c(
            map_row[["jags_name"]], term,
            paste0("(", formula_parameter, ") ", term)
          ),
          simplified = FALSE,
          stringsAsFactors = FALSE
        )
      }
      add_entry(
        quantity            = quantity,
        parameter           = map_row[["jags_name"]],
        component           = formula_parameter,
        term                = term,
        source              = formula_parameter,
        formula_parameter   = formula_parameter,
        alias_rows          = alias_rows,
        member_quantity_ids = member_quantity_ids
      )
    }
  }

  # scalar parameters with a prior of their own (not the indicators and
  # inclusions of mixture priors, and not the components of vector priors)
  scalar <- which(
    public & quantities[["role"]] == "parameter" &
      quantities[["status"]] %in% c("sampled", "structural") &
      quantities[["canonical_name"]] %in% names(prior_list)
  )
  for(row in scalar){
    quantity <- as.list(quantities[row, , drop = FALSE])
    add_entry(
      quantity          = quantity,
      parameter         = quantity[["canonical_name"]],
      component         = quantity[["namespace"]],
      term              = quantity[["canonical_name"]],
      source            = quantity[["namespace"]],
      formula_parameter = NA_character_,
      alias_rows        = native_aliases(quantity[["quantity_id"]])
    )
  }

  # random-effect quantities
  random <- which(
    public & startsWith(quantities[["role"]], "random_") &
      quantities[["status"]] != "unavailable" &
      quantities[["quantity"]] %in% .hypothesis_model_random_quantities()
  )
  for(row in random){
    quantity <- as.list(quantities[row, , drop = FALSE])
    add_entry(
      quantity          = quantity,
      parameter         = quantity[["canonical_name"]],
      component         = "random",
      term              = quantity[["canonical_name"]],
      source            = "random",
      formula_parameter = quantity[["formula_parameter"]],
      alias_rows        = rbind(
        data.frame(alias = quantity[["canonical_name"]], simplified = FALSE,
                   stringsAsFactors = FALSE),
        native_aliases(quantity[["quantity_id"]])
      ),
      # the rows of the table are named as the estimates tables name the
      # quantity; its other names still resolve through the catalog
      entry_aliases     = quantity[["canonical_name"]]
    )
  }

  if(length(entries) == 0L){
    entries <- data.frame(
      quantity_id = character(), parameter = character(),
      component = character(), term = character(), source = character(),
      formula_parameter = character(), role = character(),
      quantity = character(), status = character(),
      fixed_value = numeric(), stringsAsFactors = FALSE
    )
    entries[["aliases"]] <- I(list())
    entries[["member_quantity_ids"]] <- I(list())
    return(list(catalog = catalog, entries = entries))
  }
  entries <- do.call(rbind, entries)
  rownames(entries) <- NULL
  extension_aliases <- do.call(rbind, aliases)
  extension_quantities <- if(length(extension_quantities) == 0L){
    quantities[FALSE, , drop = FALSE]
  }else{
    do.call(rbind, extension_quantities)
  }
  catalog <- parameter_catalog_extend(
    catalog    = catalog,
    quantities = extension_quantities,
    aliases    = extension_aliases,
    provider   = "BayesToolsHypothesis"
  )

  list(catalog = catalog, entries = entries)
}


# Catalog quantities that are factor contrast coefficients ('g{1}'): their
# label parts name a contrast coefficient rather than a level cell.
.hypothesis_model_contrast_coefficient <- function(quantities){

  vapply(quantities[["label_parts"]], function(parts){
    !is.null(parts) && !is.na(parts[["coefficient"]])
  }, logical(1))
}


# The catalog quantity of a whole factor term, composed of the level cells.
.hypothesis_model_default_group_quantity <- function(catalog, map_row,
                                                     coordinates, component){

  fitted_scale  <- unique(coordinates[["fitted_scale"]])
  display_scale <- unique(coordinates[["display_scale"]])
  if(length(fitted_scale) != 1L || length(display_scale) != 1L){
    .bt_stop_refit_required(
      "Grouped formula coefficient coordinates have inconsistent scale ",
      "metadata. Refit the model with the current BayesTools build."
    )
  }
  dependencies <- unique(unlist(lapply(
    coordinates[["extraction_key"]],
    function(key) key[["dependencies"]]
  ), use.names = FALSE))
  if(!is.character(dependencies) || length(dependencies) == 0L ||
     anyNA(dependencies) || any(!nzchar(dependencies))){
    .bt_stop_refit_required(
      "Grouped formula coefficient cells have invalid extraction metadata. ",
      "Refit the model with the current BayesTools build."
    )
  }

  out <- data.frame(
    quantity_id        = paste0("BayesToolsHypothesis::formula_", map_row[["encoded_name"]]),
    canonical_name     = map_row[["jags_name"]],
    provider           = "BayesToolsHypothesis",
    namespace          = map_row[["formula_parameter"]],
    role               = "formula_coefficient_group",
    formula_parameter  = map_row[["formula_parameter"]],
    owner_type         = "",
    owner_name         = "",
    quantity           = "",
    scale_role         = "",
    parent_quantity_id = "",
    term               = map_row[["term"]],
    component          = component,
    display_label      = map_row[["term"]],
    fitted_scale       = fitted_scale,
    display_scale      = display_scale,
    status             = "derived",
    fixed_value        = NA_real_,
    internal           = FALSE,
    source_type        = "composite",
    stringsAsFactors   = FALSE,
    check.names        = FALSE
  )
  out[["arguments"]] <- I(list(character()))
  out[["extraction_key"]] <- I(list(list(
    type              = "hypothesis_formula_group",
    dependencies      = dependencies,
    formula_parameter = map_row[["formula_parameter"]],
    term              = map_row[["term"]]
  )))
  # A coefficient group has no scalar support of its own.
  out[["support"]]     <- I(list(NULL))
  out[["definedness"]] <- "always"
  out[["label_parts"]] <- I(list(NULL))

  out[, names(catalog[["quantities"]]), drop = FALSE]
}


# The refusal of resolved quantities that no entry covers: the contrast
# coefficients of a factor term are to be restated on its levels (a statement
# error); other quantities (inclusion indicators, the components of vector
# priors, latent effects) are no tested parameters.
.hypothesis_model_default_uncovered <- function(resolution, quantity_ids,
                                                symbol, metadata){

  quantities  <- metadata[["catalog"]][["quantities"]]
  rows        <- quantities[quantities[["quantity_id"]] %in% quantity_ids, , drop = FALSE]
  entries     <- metadata[["entries"]]
  covered     <- rows[["quantity_id"]] %in% c(
    entries[["quantity_id"]],
    unlist(entries[["member_quantity_ids"]], use.names = FALSE)
  )
  factor_level <- vapply(rows[["extraction_key"]], function(key){
    is.list(key) && identical(key[["type"]], "factor_level")
  }, logical(1))
  coefficient <- .hypothesis_model_contrast_coefficient(rows)
  formula     <- nzchar(rows[["formula_parameter"]])
  coefficients <- rows[(factor_level | coefficient) & formula & !covered, , drop = FALSE]
  if(nrow(coefficients) > 0L){
    coefficient <- coefficients[1L, , drop = FALSE]
    groups <- entries[
      entries[["role"]] == "formula_coefficient_group" &
        entries[["term"]] == coefficient[["term"]] &
        entries[["formula_parameter"]] == coefficient[["formula_parameter"]],
      , drop = FALSE
    ]
    example <- "their labels"
    if(nrow(groups) == 1L){
      aliases <- as.list(rep(groups[["parameter"]], length(groups[["aliases"]][[1L]])))
      names(aliases) <- groups[["aliases"]][[1L]]
      label  <- .hypothesis_model_alias_label(aliases, groups[["parameter"]])
      levels <- quantities[
        quantities[["quantity_id"]] %in% unlist(groups[["member_quantity_ids"]]) &
          quantities[["status"]] != "structural", "component"
      ]
      if(length(levels) > 0L){
        example <- paste0("their labels, such as '", label, "[", levels[[1L]], "]'")
      }
    }
    return(list(
      reason = paste0(
        "Hypotheses on factor contrast coefficients such as '",
        coefficient[["canonical_name"]], "' are not supported. State them on ",
        "factor levels by ", example, "."
      ),
      type = "statement"
    ))
  }

  list(
    reason = paste0(
      "Hypothesis tests are unavailable for '", symbol, "'. Use ",
      "hypothesis_quantities() to list the quantities that hypothesis_BF() tests."
    ),
    type = "target"
  )
}


# ---------------------------------------------------------------------------- #
# The default provider: draws
# ---------------------------------------------------------------------------- #

# The mixed posteriors of 'parameters' of a fit. The vector priors of random
# effects of the fit are left out unless they are requested, as they do not
# take part in the posterior of a fixed coefficient.
.hypothesis_model_default_mixed_posteriors <- function(
    fit, parameters, conditional, conditional_rule, transform_scaled,
    n_prior_samples){

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  if(!is.null(prior_list)){
    protected <- unique(c(parameters, conditional))
    drop <- vapply(prior_list, function(prior){
      prior_attributes <- names(attributes(prior))
      is.prior.vector(prior) && any(startsWith(prior_attributes, "random_"))
    }, logical(1))
    drop[names(prior_list) %in% protected] <- FALSE
    if(any(drop)){
      attr(fit, "prior_list") <- prior_list[!drop]
    }
  }

  as_mixed_posteriors(
    model            = fit,
    parameters       = parameters,
    conditional      = conditional,
    conditional_rule = conditional_rule,
    transform_scaled = transform_scaled,
    n_prior_samples  = n_prior_samples
  )
}


# Whether the prior of 'parameter' has a null and an alternative component.
.hypothesis_model_default_full_ensemble <- function(fit, parameter){

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  if(!parameter %in% names(prior_list) ||
     !is.prior.mixture(prior_list[[parameter]])){
    return(FALSE)
  }
  prior_components <- attr(prior_list[[parameter]], "components", exact = TRUE)

  any(prior_components == "null") && any(prior_components == "alternative")
}


# ---------------------------------------------------------------------------- #
# The default provider: random-effect quantities
# ---------------------------------------------------------------------------- #

.hypothesis_model_default_random <- function(fit, metadata){

  entry_of <- function(parameter){
    entries <- metadata()[["entries"]]
    entry   <- entries[entries[["component"]] == "random" &
                         entries[["parameter"]] == parameter, , drop = FALSE]
    if(nrow(entry) != 1L){
      stop("Random-effect quantity '", parameter, "' is not available.",
           call. = FALSE)
    }

    entry
  }
  selection_of <- function(entry){
    catalog   <- parameter_catalog(fit)
    quantity  <- catalog[["quantities"]][
      match(entry[["quantity_id"]], catalog[["quantities"]][["quantity_id"]]), ,
      drop = FALSE
    ]
    parameter_catalog_resolve(
      catalog,
      alias     = quantity[["canonical_name"]],
      namespace = quantity[["namespace"]]
    )
  }

  list(
    select = function(parameter, transform_scaled){
      entry     <- entry_of(parameter)
      selection <- selection_of(entry)
      quantity  <- selection[["quantities"]]
      # the quantity on the scale of the request (its support is the one of the
      # draws): the catalog quantity of the fitted scale for transform_scaled =
      # FALSE
      on_scale  <- if(transform_scaled) selection else
        .bt_parameter_fitted_scale(fit, selection)[["selection"]]
      list(
        entry = list(
          parameter   = entry[["parameter"]],
          quantity_id = entry[["quantity_id"]],
          selection   = selection,
          on_scale    = on_scale
        ),
        spec  = list(
          label   = quantity[["canonical_name"]],
          status  = quantity[["status"]],
          quantity = quantity[["quantity"]]
        )
      )
    },
    mixed_posterior = function(parameter, transform_scaled, conditional, selected){
      if(is.null(selected)){
        selected <- list(entry = list(selection = selection_of(entry_of(parameter))))
      }
      values <- parameter_mixed_posterior(
        fit              = fit,
        selection        = selected[["entry"]][["selection"]],
        conditional      = conditional,
        transform_scaled = transform_scaled
      )
      attr(values, "parameter") <- parameter
      out <- list(values)
      names(out) <- parameter
      attr(out, "prior_list") <- stats::setNames(list(prior_none()), parameter)
      attr(out, "random_parameter_label") <- parameter
      class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")

      out
    },
    support = function(selected){
      support <- selected[["entry"]][["on_scale"]][["quantities"]][["support"]][[1L]]
      if(is.null(support) || !isTRUE(support[["exact"]])){
        return(c(-Inf, Inf))
      }

      as.numeric(support[["bounds"]])
    },
    density_target = function(parameter, transform_scaled, operation){
      list(
        parameter         = parameter,
        parameter_spec    = list(type = "primitive"),
        display_transform = NULL
      )
    },
    sd_pair = function(selected){
      .hypothesis_model_default_sd_pair(fit, selected)
    },
    zero_boundary = function(selected){
      quantity <- selected[["entry"]][["selection"]][["quantities"]]
      key      <- quantity[["extraction_key"]][[1L]]
      if(!identical(quantity[["quantity"]], "sd") ||
         !isTRUE(key[["allocation_derived"]])){
        return(NULL)
      }

      list(alternative = NULL)
    },
    defined_share = function(values, selected, conditional){
      if(is.null(posterior_metadata(values, "undefined_draws"))){
        return(NULL)
      }
      n_draws <- nrow(as.matrix(.fit_to_posterior(fit)))
      if(isTRUE(conditional)){
        states  <- parameter_gate_states(fit, selected[["entry"]][["selection"]])
        n_draws <- sum(states[["event"]])
      }
      if(length(values) >= n_draws) return(NULL)

      list(defined = length(values), eligible = n_draws)
    },
    defined_footnote = .hypothesis_model_defined_footnote,
    defined_draws = .hypothesis_model_defined_draws,
    prior_draws = function(parameter, selected, transform_scaled, n_samples, seed){
      raw_samples <- transform_prior_samples(
        fit           = fit,
        n_samples     = n_samples,
        seed          = seed,
        formula_scale = list()
      )
      extraction_fit <- JAGS_with_draws(
        fit, coda::mcmc.list(coda::mcmc(raw_samples))
      )
      if(!transform_scaled){
        attr(extraction_fit, "formula_scale") <- list()
      }
      draws   <- parameter_draws(extraction_fit, selected[["entry"]][["selection"]])
      samples <- as.matrix(draws)
      colnames(samples) <- parameter
      declared <- posterior_metadata(draws, "undefined_draws")
      if(!is.null(declared)){
        posterior_metadata(samples, "undefined_draws") <-
          stats::setNames(unname(declared[[1L]]), parameter)
      }

      list(samples = samples)
    }
  )
}


# The standard deviation whose square is the selected variance quantity (the
# catalog quantity with the same owner and component, and the standard
# deviation evaluator of the variance's extraction key), or NULL.
.hypothesis_model_default_sd_pair <- function(fit, selected){

  quantity <- selected[["entry"]][["selection"]][["quantities"]]
  key      <- quantity[["extraction_key"]][[1L]]
  sd_evaluator <- switch(
    as.character(key[["evaluator"]]),
    sd_variance    = "sd",
    allocation_var = "allocation_sd",
    NULL
  )
  if(is.null(sd_evaluator)){
    return(NULL)
  }
  quantities <- parameter_catalog(fit)[["quantities"]]
  same <- function(field){
    vapply(quantities[[field]], identical, logical(1), y = quantity[[field]][[1L]])
  }
  candidates <- which(
    !quantities[["internal"]] &
      same("formula_parameter") & same("owner_type") & same("owner_name") &
      same("component") &
      vapply(quantities[["extraction_key"]], function(other){
        is.list(other) && identical(other[["type"]], key[["type"]]) &&
          identical(other[["evaluator"]], sd_evaluator) &&
          identical(other[["random_block"]], key[["random_block"]]) &&
          identical(other[["dependencies"]], key[["dependencies"]])
      }, logical(1))
  )
  if(length(candidates) != 1L){
    return(NULL)
  }

  quantities[["canonical_name"]][candidates]
}


# Draws of a quantity declared as possibly undefined (an original-scale
# correlation with a zero SD; 'undefined_draws' metadata of the extracted
# samples) are left out where the quantity is undefined. Any other missing
# draw is an error. Returns the mask of defined draws.
.hypothesis_model_defined_draws <- function(values, samples, label){

  defined <- !is.na(values)
  if(all(defined)){
    return(defined)
  }
  if(is.null(posterior_metadata(samples, "undefined_draws"))){
    stop(
      "The draws of random-effect quantity '", label, "' contain missing ",
      "values. Missing draws are accepted only for quantities declared as ",
      "possibly undefined (original-scale random-effect correlations).",
      call. = FALSE
    )
  }
  if(!any(defined)){
    stop(
      "Random-effect quantity '", label, "' is undefined in every draw.",
      call. = FALSE
    )
  }

  defined
}


# Footnote for results computed over the defined draws only.
.hypothesis_model_defined_footnote <- function(label, samples, posterior_defined,
                                               prior_defined = NULL){

  counts <- function(defined, what){
    if(is.null(defined)) return(NULL)
    if(is.list(defined)){
      if(defined[["defined"]] == defined[["eligible"]]) return(NULL)
      return(paste(defined[["defined"]], "of", defined[["eligible"]], what, "draws"))
    }
    if(all(defined)) return(NULL)
    paste(sum(defined), "of", length(defined), what, "draws")
  }
  parts <- c(
    counts(posterior_defined, "posterior"),
    counts(prior_defined, "prior")
  )
  if(length(parts) == 0L){
    return(NULL)
  }
  reason <- posterior_metadata(samples, "undefined_draws")
  condition <- if(identical(unname(reason[[1L]]), "correlation")){
    "where the correlation is defined, i.e. both SDs are positive."
  }else{
    "where the quantity is defined."
  }

  paste0(
    label, ": computed from ", paste(parts, collapse = " and "), " ",
    condition
  )
}
