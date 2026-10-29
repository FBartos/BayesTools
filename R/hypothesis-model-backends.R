# ============================================================================ #
# hypothesis-model-backends.R
# ============================================================================ #
#
# Registry of the density backends of the hypothesis engine
# (hypothesis_BF_model()): the methods a point hypothesis can be evaluated
# with. The built-ins evaluate the posterior draws ("KDE", "normal") or draws
# that carry their ordinates ("precomputed"); a downstream package registers an
# ordinate backend that computes its own posterior ordinates. See
# .agents/instructions/hypothesis.md.
#
# ============================================================================ #


#' @title Density backends of fit-level hypothesis Bayes factors
#'
#' @description The density backend of a point hypothesis (the
#' \code{density_method} of [hypothesis_BF()] and [hypothesis_BF_model()])
#' decides how the posterior ordinate at the null value is obtained. The
#' built-in backends \code{"KDE"} and \code{"normal"} estimate it from the
#' posterior draws, and \code{"precomputed"} uses ordinates that are stored on
#' the draws (it is for draw inputs only, not for fitted models). A package
#' that computes posterior ordinates itself, for example by a likelihood-aware
#' estimator that cannot be written as a function of the draws, registers an
#' \emph{ordinate} backend. [hypothesis_density_backend()] creates a backend
#' object, [register_hypothesis_density_backend()] adds it to the registry (a
#' package registers its backends in \code{.onLoad()} and removes them with
#' [unregister_hypothesis_density_backend()] in \code{.onUnload()}), and
#' [hypothesis_density_backends()] lists the registered backends.
#'
#' @details A backend has a \code{type}: \code{"draws"} (evaluated from the
#' posterior draws by [hypothesis_BF()] with the density method
#' \code{method}), \code{"ordinate"} (the backend attaches the posterior
#' ordinates of the tested values to the draws; they are then evaluated as
#' \code{"precomputed"} ordinates), or \code{"stored"} (the ordinates are
#' already attached to the draws). Matching of \code{density_method} to the
#' registered names ignores case and accepts unambiguous partial names.
#'
#' The engine calls the callbacks of a backend as follows; each is optional
#' except \code{attach} of an ordinate backend.
#' \describe{
#'   \item{\code{available(object, plan, control)}}{whether the backend can
#'   evaluate the statement \code{plan} of the model \code{object} (the
#'   \code{capability_object} of the provider) with the normalized
#'   \code{control}. It is evaluated on every call and never cached. It
#'   returns \code{NULL} when the backend can, and otherwise a refusal: a list
#'   with the \code{reason} (a string) and optionally \code{class} (additional
#'   condition classes appended to \code{BayesTools_hypothesis_method} and
#'   \code{BayesTools_hypothesis_unavailable}) and \code{fields} (a named list
#'   of further fields of the condition). The engine stops with the refusal
#'   before any statement is evaluated.}
#'   \item{\code{control(control)}}{the normalized \code{density_control} of
#'   the backend (defaults filled in), stopping on invalid settings. The
#'   default accepts no settings: \code{NULL} and an empty list give
#'   \code{NULL}, anything else is an error.}
#'   \item{\code{attach(posterior, request, control, cache)}}{the whole
#'   ordinate step of an ordinate backend. \code{posterior} is the marginal
#'   posterior of the tested quantity (a vector, or for the request of a
#'   factor level the list of the posteriors of the levels of the factor
#'   term); it returns it with the ordinates of the requested values stored in
#'   its \code{posterior_ordinate} metadata ([posterior_metadata()],
#'   [posterior_ordinate_attribute()], [posterior_ordinate_append()]), one
#'   entry per value, and the \code{parameter}, \code{conditional},
#'   \code{conditional_rule}, and \code{condition_key} metadata of the draws,
#'   or stops with its own typed condition when an ordinate cannot be
#'   constructed. \code{cache} is an environment that lives for the
#'   evaluation of one group of statements, for objects that the requests of
#'   the group share. The \code{request} is a list with the elements
#'   \code{kind} (\code{"scalar"}: the values of the vector \code{posterior};
#'   \code{"level"}: the values of the level \code{level}), \code{route} (the
#'   kind of the statement plan: \code{"scalar"}, \code{"linear"},
#'   \code{"random"}, \code{"exp_affine"}, or \code{"marginal_mean"}),
#'   \code{object} (the model object, or the object that the provider names as
#'   its \code{capability_object}), \code{model} (the object of the provider),
#'   \code{parameter}, \code{label} (the tested quantity, as the table names
#'   it), \code{level}, \code{values} (the point values to compute),
#'   \code{conditional} (the conditioning of the draws: the parameter whose
#'   inclusion event the draws are conditioned on, or \code{NULL}),
#'   \code{spec} (the structure of the tested quantity that the provider
#'   declares: for example \code{type = "primitive"} or \code{"linear"} with
#'   the \code{weights} on the fitted coordinates and the \code{prior_density}),
#'   \code{display_transform} (the transformation from the evaluated to the
#'   tested scale, or \code{NULL}), \code{density_method} (the name of the
#'   backend), and \code{inference_type} (marginal means: the marginal
#'   inference the draws belong to).}
#'   \item{\code{reuse(entry, request, control)}}{whether the ordinate
#'   \code{entry} that is stored on the draws at a requested value is
#'   current for the request (the method, control, source of the fit, target,
#'   and conditioning agree), so that it is not computed again. The default
#'   reuses nothing: stored ordinates at the requested values are replaced.}
#'   \item{\code{warnings(posterior, request)}}{the warnings of the ordinates
#'   stored on \code{posterior} (recursively over the levels of a list), as a
#'   data frame with the columns \code{warning}, \code{parameter},
#'   \code{level}, \code{label}, and \code{value}. The engine puts each
#'   warning on the rows of the statements it belongs to. The default reads
#'   the \code{warning} element of the diagnostics of each stored ordinate.}
#'   \item{\code{diagnostics(posterior, request)}}{the public diagnostics of
#'   the stored ordinates (a data frame), stored as the
#'   \code{density_diagnostics} attribute of the table, and
#'   \code{bind_diagnostics(list)}, which combines the diagnostics of several
#'   groups of statements. Without \code{diagnostics} the table has no such
#'   attribute; the default of \code{bind_diagnostics} row-binds the data
#'   frames.}
#'   \item{\code{rescale_failure(condition, requested, evaluated, label,
#'   density_method)}}{restates a condition that \code{attach} raised for the
#'   evaluated values \code{evaluated} of a request for the values
#'   \code{requested} of a quantity that is evaluated through a
#'   transformation of itself (a variance through its standard deviation); it
#'   returns the condition to signal, or \code{NULL} to signal the original
#'   one.}
#'   \item{\code{refusal_classes(cause)}}{additional condition classes for the
#'   refusals of the ordinate backends for a generic cause
#'   (\code{"original_scale"}, \code{"conditional_random"}, or
#'   \code{"random_target"}).}
#'   \item{\code{clear_cache()}}{drops the caches of the backend; the engine
#'   calls it when it drops its own caches.}
#' }
#' \code{policy} is a list of the settings that the backend declares for its
#' acceptance of ordinates (thresholds on its error estimates, for example). It
#' is reported by [hypothesis_density_backends()] and not used by the engine:
#' the backend applies its own acceptance rules in \code{attach}, and the
#' engine only checks that the attached ordinates are valid
#' ([posterior_ordinate_supports_bf()]).
#'
#' @param name the name of the density method (a single string that starts
#'   with a letter).
#' @param type \code{"ordinate"}, \code{"draws"}, or \code{"stored"}.
#' @param available,control,attach,reuse,warnings,diagnostics,bind_diagnostics,rescale_failure,refusal_classes,clear_cache
#'   the callbacks of the backend, functions or \code{NULL} for the default;
#'   see Details.
#' @param method for a \code{"draws"} backend, the evaluation of
#'   [hypothesis_BF()] that it stands for: \code{"KDE"} or \code{"normal"}.
#' @param policy a list of declared settings of the backend.
#' @param advertised whether [hypothesis_quantities()] lists the backend among
#'   the methods that evaluate a quantity.
#' @param package the name of the package that provides the backend.
#' @param backend a backend created with [hypothesis_density_backend()].
#' @param overwrite whether registering replaces a registered backend of the
#'   same name. Replacing a built-in backend is not possible.
#'
#' @return \code{hypothesis_density_backend()} returns an object of class
#'   \code{BayesTools_hypothesis_density_backend}.
#'   \code{register_hypothesis_density_backend()} returns the backend
#'   invisibly, and \code{unregister_hypothesis_density_backend()} returns
#'   whether a backend was removed (invisibly).
#'   \code{hypothesis_density_backends()} returns a data frame with one row per
#'   registered backend in the order of registration: \code{name},
#'   \code{type}, \code{package}, \code{advertised}, and the list column
#'   \code{policy}.
#'
#' @examples
#' hypothesis_density_backends()
#'
#' # a toy backend that evaluates the ordinate of a normal posterior exactly
#' normal_ordinate <- hypothesis_density_backend(
#'   name   = "exactNormal",
#'   type   = "ordinate",
#'   attach = function(posterior, request, control, cache) {
#'     draws    <- as.numeric(posterior)
#'     ordinate <- stats::dnorm(request$values, mean(draws), stats::sd(draws))
#'     posterior_metadata(posterior, "posterior_ordinate") <-
#'       posterior_ordinate_attribute(
#'         value = request$values, ordinate = ordinate,
#'         method = "exact normal", density_method = "exactNormal",
#'         parameter = attr(posterior, "parameter")
#'       )
#'     posterior
#'   }
#' )
#' register_hypothesis_density_backend(normal_ordinate)
#' hypothesis_density_backends()$name
#' unregister_hypothesis_density_backend("exactNormal")
#'
#' @seealso [hypothesis_BF()], [hypothesis_BF_model()],
#'   [hypothesis_provider()]
#' @rdname hypothesis_density_backend
#' @export
hypothesis_density_backend <- function(name, type = c("ordinate", "draws", "stored"),
                                       available = NULL, control = NULL,
                                       attach = NULL, reuse = NULL,
                                       warnings = NULL, diagnostics = NULL,
                                       bind_diagnostics = NULL,
                                       rescale_failure = NULL,
                                       refusal_classes = NULL,
                                       clear_cache = NULL, method = NULL,
                                       policy = list(), advertised = TRUE,
                                       package = NULL){

  check_char(name, "name", check_length = 1, allow_NA = FALSE)
  if(!grepl("^[A-Za-z][A-Za-z0-9._]*$", name)){
    stop("'name' must start with a letter and contain only letters, digits, '.' and '_'.",
         call. = FALSE)
  }
  type <- match.arg(type)
  callbacks <- list(
    available = available, control = control, attach = attach, reuse = reuse,
    warnings = warnings, diagnostics = diagnostics,
    bind_diagnostics = bind_diagnostics, rescale_failure = rescale_failure,
    refusal_classes = refusal_classes, clear_cache = clear_cache
  )
  for(callback in names(callbacks)){
    if(!is.null(callbacks[[callback]]) && !is.function(callbacks[[callback]])){
      stop("'", callback, "' must be a function or NULL.", call. = FALSE)
    }
  }
  if(identical(type, "ordinate") && is.null(attach)){
    stop("An ordinate backend needs an 'attach' function.", call. = FALSE)
  }
  if(!identical(type, "ordinate") && !is.null(attach)){
    stop("Only an ordinate backend has an 'attach' function.", call. = FALSE)
  }
  if(identical(type, "draws")){
    if(is.null(method)){
      method <- name
    }
    method <- match.arg(method, c("KDE", "normal"))
  }else if(!is.null(method)){
    stop("'method' is used by the backends of type 'draws' only.", call. = FALSE)
  }
  if(!is.list(policy)){
    stop("'policy' must be a list.", call. = FALSE)
  }
  check_bool(advertised, "advertised", allow_NA = FALSE)
  check_char(package, "package", check_length = 1, allow_NA = FALSE,
             allow_NULL = TRUE)

  out <- c(
    list(name = name, type = type, method = method),
    callbacks,
    list(policy = policy, advertised = advertised, package = package)
  )
  class(out) <- c("BayesTools_hypothesis_density_backend", "list")

  out
}


#' @rdname hypothesis_density_backend
#' @export
register_hypothesis_density_backend <- function(backend, overwrite = FALSE){

  if(!inherits(backend, "BayesTools_hypothesis_density_backend")){
    stop("'backend' must be created with 'hypothesis_density_backend()'.",
         call. = FALSE)
  }
  check_bool(overwrite, "overwrite", allow_NA = FALSE)
  registry <- .hypothesis_model_backend_registry()
  name     <- backend[["name"]]
  if(!is.null(registry[[name]])){
    if(isTRUE(registry[[name]][["builtin"]])){
      stop("The built-in density backend '", name, "' cannot be replaced.",
           call. = FALSE)
    }
    if(!overwrite){
      stop("A density backend named '", name, "' is registered already; use ",
           "'overwrite = TRUE' to replace it.", call. = FALSE)
    }
  }
  # Names that differ only in case or by a prefix would make the matching of
  # 'density_method' ambiguous in exactly the cases it is meant to resolve.
  others <- setdiff(names(registry), name)
  clash  <- others[tolower(others) == tolower(name)]
  if(length(clash) > 0L){
    stop("The name '", name, "' differs from the registered density backend '",
         clash[[1L]], "' by case only.", call. = FALSE)
  }
  if(!is.null(registry[[name]])){
    .hypothesis_model_backend_clear_cache(registry[[name]][["backend"]])
  }
  entry <- list(backend = backend, builtin = FALSE)
  registry[[name]] <- entry
  .hypothesis_model_backend_registry_order(add = name)

  invisible(backend)
}


#' @rdname hypothesis_density_backend
#' @export
unregister_hypothesis_density_backend <- function(name){

  check_char(name, "name", check_length = 1, allow_NA = FALSE)
  registry <- .hypothesis_model_backend_registry()
  entry    <- registry[[name]]
  if(is.null(entry)){
    return(invisible(FALSE))
  }
  if(isTRUE(entry[["builtin"]])){
    stop("The built-in density backend '", name, "' cannot be removed.",
         call. = FALSE)
  }
  .hypothesis_model_backend_clear_cache(entry[["backend"]])
  rm(list = name, envir = registry)
  .hypothesis_model_backend_registry_order(remove = name)

  invisible(TRUE)
}


#' @rdname hypothesis_density_backend
#' @export
hypothesis_density_backends <- function(){

  backends <- .hypothesis_model_backends()
  out <- data.frame(
    name       = vapply(backends, `[[`, character(1), "name"),
    type       = vapply(backends, `[[`, character(1), "type"),
    package    = vapply(backends, function(backend) {
      if(is.null(backend[["package"]])) NA_character_ else backend[["package"]]
    }, character(1)),
    advertised = vapply(backends, `[[`, logical(1), "advertised"),
    stringsAsFactors = FALSE
  )
  out[["policy"]] <- I(lapply(backends, `[[`, "policy"))
  rownames(out) <- NULL

  out
}


# ---------------------------------------------------------------------------- #
# Registry
# ---------------------------------------------------------------------------- #

# The registry is an environment of the package: the registered backends by
# name, and the order of their registration (the order of the methods that the
# engine offers).
.hypothesis_model_backend_registry <- function(){

  registry <- .BayesTools_private$hypothesis_backends
  if(is.null(registry)){
    registry <- .hypothesis_model_backend_registry_reset()
  }

  registry
}


.hypothesis_model_backend_registry_reset <- function(){

  registry <- new.env(parent = emptyenv())
  .BayesTools_private$hypothesis_backends      <- registry
  .BayesTools_private$hypothesis_backend_order <- character()
  for(backend in .hypothesis_model_builtin_backends()){
    registry[[backend[["name"]]]] <- list(backend = backend, builtin = TRUE)
    .hypothesis_model_backend_registry_order(add = backend[["name"]])
  }

  registry
}


.hypothesis_model_backend_registry_order <- function(add = NULL, remove = NULL){

  order <- .BayesTools_private$hypothesis_backend_order
  if(is.null(order)){
    order <- character()
  }
  order <- c(setdiff(order, c(add, remove)), add)
  .BayesTools_private$hypothesis_backend_order <- order

  invisible(order)
}


# The registered backends in the order of registration.
.hypothesis_model_backends <- function(){

  registry <- .hypothesis_model_backend_registry()
  order    <- .BayesTools_private$hypothesis_backend_order
  order    <- order[order %in% names(registry)]

  stats::setNames(lapply(order, function(name) registry[[name]][["backend"]]), order)
}


.hypothesis_model_backend_clear_cache <- function(backend){

  if(is.function(backend[["clear_cache"]])){
    backend[["clear_cache"]]()
  }

  invisible(NULL)
}


# The names a request for 'density_method' can match: the registered backends
# that a fitted model can use (not the stored-ordinate method of draw inputs),
# unless the provider names its own.
.hypothesis_model_density_methods <- function(provider = NULL){

  if(!is.null(provider[["density_methods"]])){
    return(provider[["density_methods"]])
  }
  backends <- .hypothesis_model_backends()
  names(backends)[vapply(backends, function(backend) {
    !identical(backend[["type"]], "stored")
  }, logical(1))]
}


# The registered name that 'density_method' stands for. Matching ignores case
# and accepts an unambiguous partial name; NULL selects the first allowed
# name.
.hypothesis_model_density_method <- function(density_method, allowed){

  if(is.null(density_method)){
    density_method <- allowed[[1L]]
  }
  if(is.character(density_method) && length(density_method) > 0L &&
     !anyNA(density_method)){
    density_method_lower <- tolower(density_method)
    allowed_lower        <- tolower(allowed)
    if(length(density_method) > 1L){
      if(identical(density_method_lower, allowed_lower)){
        density_method <- allowed
      }
    }else{
      index <- pmatch(density_method_lower, allowed_lower, duplicates.ok = TRUE)
      if(!is.na(index)){
        density_method <- allowed[[index]]
      }
    }
  }

  posterior_density_method_match(
    method  = density_method,
    allowed = allowed,
    name    = "density_method"
  )
}


# The registered backend named 'name'.
.hypothesis_model_backend <- function(name){

  backend <- .hypothesis_model_backend_registry()[[name]][["backend"]]
  if(is.null(backend)){
    stop("The density backend '", name, "' is not registered.", call. = FALSE)
  }

  backend
}


# The normalized control of a backend for a request, with the checks that
# apply to every backend: settings are validated by the backend, and a backend
# that has none refuses them.
.hypothesis_model_backend_control <- function(backend, density_control){

  if(is.function(backend[["control"]])){
    return(backend[["control"]](density_control))
  }
  if(is.null(density_control)){
    return(NULL)
  }
  if(!is.list(density_control)){
    stop("'density_control' must be a named list.", call. = FALSE)
  }
  if(length(density_control) == 0L){
    return(NULL)
  }

  stop("'density_control' is not used by density_method = '",
       backend[["name"]], "'.", call. = FALSE)
}


# The additional condition classes with which the backend refuses a generic
# cause (see hypothesis_density_backend()).
.hypothesis_model_backend_cause_classes <- function(backend, cause){

  if(is.null(cause) || !is.function(backend[["refusal_classes"]])){
    return(NULL)
  }

  backend[["refusal_classes"]](cause)
}


# ---------------------------------------------------------------------------- #
# Built-in backends
# ---------------------------------------------------------------------------- #

.hypothesis_model_builtin_backends <- function(){

  list(
    hypothesis_density_backend(
      name       = "KDE",
      type       = "draws",
      method     = "KDE",
      package    = "BayesTools"
    ),
    hypothesis_density_backend(
      name       = "normal",
      type       = "draws",
      method     = "normal",
      advertised = FALSE,
      package    = "BayesTools"
    ),
    hypothesis_density_backend(
      name       = "precomputed",
      type       = "stored",
      available  = function(object, plan, control){
        list(
          reason = paste0(
            "density_method = 'precomputed' evaluates ordinates that are ",
            "stored on posterior draws; it is not available for fitted models."
          ),
          type = "method"
        )
      },
      advertised = FALSE,
      package    = "BayesTools"
    )
  )
}
