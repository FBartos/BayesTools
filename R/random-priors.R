#' Random-effect prior specification
#'
#' @description
#' `prior_random()` defines priors and policies for formula random effects.
#' Random-effect formulas own the grouping variables, slopes, intercepts, and
#' covariance structure. `prior_random()` supplies priors and monitoring policy
#' for the already parsed random-effect blocks.
#'
#' @details
#' Formula random effects in [JAGS_formula()] and [JAGS_fit()] require an
#' explicit `prior_random()` object. The top-level arguments are defaults used
#' by every random-effect block unless a named `random_block()` override is
#' supplied in `...`.
#'
#' Each override name must match a stable random-effect block label generated
#' from the formula parser or supplied by the user through named random terms.
#' This keeps prior specification separate from formula syntax. For example,
#' if the formula contains a block named `"study"`, then
#' `prior_random(study = random_block(...))` overrides only that block.
#'
#' `random_term()` is an alias of `random_block()` retained for readability in
#' user-facing code.
#'
#' @section Standard-deviation priors:
#' `sd` must be a BayesTools prior object with nonnegative support, such as a
#' truncated normal, gamma, point prior at zero, or mixture prior. When a
#' structure needs multiple marginal standard deviations, a single `sd` prior is
#' replicated over those SD components unless a more specific override is
#' supplied by the random-effect resolver.
#'
#' For ordinary estimation, use continuous positive SD priors. Use point or
#' mixture priors at zero only when a random-effect component is intentionally
#' fixed or model-averaged.
#'
#' @section Covariance structures:
#' `random_covariance()` stores structure-specific covariance settings.
#' Supported structure labels are:
#'
#' * `"ID"`: independent columns with one shared SD.
#' * `"DIAG"`: independent columns with one SD per generated column.
#' * `"CS"` / `"HCS"`: compound symmetry over one or more index variables,
#'   with one scalar correlation `rho`; `"HCS"` has one SD per index level.
#'   Multiple index variables in `cs(index1 + index2 | group)` are combined
#'   with `interaction(..., drop = TRUE)`.
#' * `"AR"` / `"AR1"` / `"HAR"`: discrete autoregressive structure over a
#'   single ordered index variable, with one scalar correlation `rho`; `"HAR"`
#'   has one SD per index level.
#' * `"CAR"`: continuous-time AR(1), currently for `car(time | group)`,
#'   using `rho ^ abs(t_i - t_j)` with `0 < rho < 1`.
#' * `"UN"` / `"US"`: unstructured covariance with SD priors and `prior_lkj()`.
#'
#' The formula parser owns the actual structure used by a block. Supplying a
#' conflicting `structure` in `random_covariance()` is rejected when the formula
#' is resolved.
#'
#' @section Correlation priors:
#' Use `prior_lkj(eta = ...)` for full unstructured correlation matrices. The
#' compiled JAGS backend is used by default, with a syntax fallback available
#' through `backend = "syntax"`.
#'
#' Use `rho` for scalar-correlation structures (`"CS"`, `"HCS"`, `"AR"`,
#' `"HAR"`, and `"CAR"`). `rho_scale` controls the scale on which the prior is
#' placed:
#'
#' * `"fisher_z"` places the prior on `atanh(rho)` and is the default. A point
#'   prior at zero corresponds to independence. For one-sided structures such
#'   as `"CAR"`, the Fisher-z prior is truncated to the valid raw-rho interval.
#' * `"logit"` places the prior on a logit transform of the valid raw-rho
#'   interval.
#' * `"rho"` places the prior directly on raw `rho`; the prior is truncated to
#'   the structure-specific valid interval.
#'
#' @section Variance allocation:
#' `random_variance_allocation()` expresses one total SD prior plus a Dirichlet
#' prior over variance fractions. Multiple named allocations can be supplied in
#' `prior_random()`; child allocations can inherit one component of an earlier
#' allocation through `parent = allocation_ref(...)`.
#'
#' \deqn{\sigma_j = \sigma_{\mathrm{total}}\sqrt{w_j}, \quad
#'       w \sim \mathrm{Dirichlet}(\alpha).}
#'
#' This is useful for nested or crossed random intercepts when the prior should
#' control total heterogeneity separately from how that heterogeneity is
#' allocated across components. The allocation targets named random-effect
#' blocks through `terms`. If `terms = NULL`, the allocation targets all
#' remaining unallocated random-effect blocks at resolution time. Generated JAGS
#' syntax supports allocated blocks with one SD component, homogeneous
#' structures where one SD controls the block, and explicit `components = "sd"`
#' allocations that split the resolved SD leaves of one heterogeneous block.
#' Use `scale = "total_variance"` to preserve summed variance and
#' `scale = "mean_variance"` to preserve average variance inside a
#' heterogeneous block.
#'
#' @section Monitoring and prediction:
#' `random_monitor()` controls which random-effect quantities are saved from
#' JAGS. The default monitors standardized latent random effects because they
#' are required for bridge sampling and for reconstructing group-level
#' coefficients in prediction. Turning `latent = FALSE` can reduce memory use,
#' but bridge sampling and conditional prediction may require refitting.
#'
#' `random_new_levels()` is reserved for future prediction support for new
#' grouping levels. It currently requires `allow = FALSE`.
#'
#' @param ... named `random_block()` specifications and, optionally,
#'   top-level `random_variance_allocation()` specifications.
#' @param sd prior distribution for random-effect standard deviations.
#' @param covariance optional `random_covariance()` specification.
#' @param cor optional correlation prior. Use `prior_lkj()` for unstructured
#'   correlation matrices.
#' @param rho optional scalar correlation prior for structured correlation
#'   models. The default scale is Fisher's z.
#' @param monitor `random_monitor()` object.
#' @param new_levels `random_new_levels()` object. New-level prediction is
#'   reserved for a future release and currently must leave `allow = FALSE`.
#' @param allocation optional `random_variance_allocation()` specification, or
#'   a list of such specifications, defining total-SD plus Dirichlet variance
#'   allocation across named random-effect blocks. Block-local allocation is
#'   reserved for a future release.
#'
#' @return A list-like S3 object describing random-effect priors, covariance
#'   settings, allocation settings, and monitoring policy.
#'
#' @examples
#' sd_prior <- prior(
#'   "normal", list(mean = 0, sd = 0.5),
#'   truncation = list(lower = 0, upper = Inf)
#' )
#' rho_prior <- prior("normal", list(mean = 0, sd = 0.5))
#'
#' # One SD prior used for every random-effect block.
#' prior_random(sd = sd_prior)
#'
#' # Override a named block while keeping a global default.
#' prior_random(
#'   sd = sd_prior,
#'   study = random_block(sd = prior("gamma", list(shape = 2, rate = 2)))
#' )
#'
#' # Unstructured covariance: SD prior plus LKJ correlation prior.
#' prior_random(
#'   study = random_term(
#'     covariance = random_covariance(
#'       structure = "UN",
#'       sd = sd_prior,
#'       cor = prior_lkj(eta = 2)
#'     )
#'   )
#' )
#'
#' # Scalar-rho structure with a Fisher-z prior.
#' prior_random(
#'   study_time = random_term(
#'     covariance = random_covariance(
#'       structure = "AR",
#'       sd = sd_prior,
#'       rho = rho_prior,
#'       rho_scale = "fisher_z"
#'     )
#'   )
#' )
#'
#' # Total SD plus Dirichlet variance allocation.
#' prior_random(
#'   random_variance_allocation(
#'     terms = c("study", "drug"),
#'     sd = sd_prior,
#'     allocation = prior("dirichlet", list(alpha = c(1, 1)))
#'   )
#' )
#'
#' # A child allocation can split one inherited variance budget.
#' prior_random(
#'   random_variance_allocation(
#'     name = "total_re",
#'     terms = c(nested = "nested", simple = "study"),
#'     sd = sd_prior,
#'     allocation = prior("dirichlet", list(alpha = c(1, 1)))
#'   ),
#'   random_variance_allocation(
#'     name = "nested_split",
#'     parent = allocation_ref("total_re", "nested"),
#'     terms = c("paper", "country"),
#'     allocation = prior("dirichlet", list(alpha = c(3, 1)))
#'   )
#' )
#' @export
prior_random <- function(..., sd = NULL, covariance = NULL, cor = NULL,
                         rho = NULL, monitor = random_monitor(),
                         new_levels = random_new_levels(),
                         allocation = NULL){

  blocks <- list(...)
  dot_allocations <- vapply(blocks, inherits, logical(1), what = "random_variance_allocation")
  if(any(dot_allocations)){
    if(!is.null(allocation)){
      stop("Variance allocation priors must be supplied either in '...' or in 'allocation', not both.", call. = FALSE)
    }
    allocation <- blocks[dot_allocations]
    blocks <- blocks[!dot_allocations]
    if(length(allocation) == 1L &&
       (is.null(names(allocation)) || !nzchar(names(allocation)))){
      allocation <- allocation[[1L]]
    }
  }
  .bt_check_random_prior_blocks(blocks)
  .bt_check_random_sd_prior(sd, allow_NULL = TRUE)

  if(is.null(covariance)){
    covariance <- random_covariance(cor = cor, rho = rho)
  }else{
    .bt_check_random_covariance(covariance)
    if(!is.null(cor) || !is.null(rho)){
      stop("'cor' and 'rho' cannot be supplied together with 'covariance'.", call. = FALSE)
    }
  }

  .bt_check_random_monitor(monitor)
  .bt_check_random_new_levels(new_levels)
  .bt_check_random_allocation(allocation)

  out <- list(
    sd         = sd,
    covariance = covariance,
    monitor    = monitor,
    new_levels = new_levels,
    allocation = allocation,
    blocks     = blocks
  )
  class(out) <- c("prior_random", "list")

  out
}

#' @rdname prior_random
#' @param terms for `random_block()`, an optional named list of term-specific
#'   SD overrides. Each entry must be either a prior object or
#'   `random_block(sd = ...)`, and names must match resolved model terms such
#'   as `"intercept"` or a slope term label. For
#'   `random_variance_allocation()`, optional random-effect block names targeted
#'   by the shared allocation. If `NULL`, the allocation targets all remaining
#'   unallocated random-effect blocks at resolution time.
#' @export
random_block <- function(sd = NULL, covariance = NULL, cor = NULL, rho = NULL,
                         monitor = NULL, new_levels = NULL, allocation = NULL,
                         terms = NULL){

  .bt_check_random_sd_prior(sd, allow_NULL = TRUE)

  if(is.null(covariance)){
    covariance <- if(is.null(cor) && is.null(rho)){
      NULL
    }else{
      random_covariance(cor = cor, rho = rho)
    }
  }else{
    .bt_check_random_covariance(covariance)
    if(!is.null(cor) || !is.null(rho)){
      stop("'cor' and 'rho' cannot be supplied together with 'covariance'.", call. = FALSE)
    }
  }
  if(!is.null(monitor)){
    .bt_check_random_monitor(monitor)
  }
  if(!is.null(new_levels)){
    .bt_check_random_new_levels(new_levels)
  }
  if(!is.null(allocation)){
    stop(
      "Block-local variance allocation priors are not implemented yet; use top-level prior_random(allocation = random_variance_allocation(...)).",
      call. = FALSE
    )
  }

  out <- list(
    sd         = sd,
    covariance = covariance,
    monitor    = monitor,
    new_levels = new_levels,
    allocation = allocation,
    terms      = terms
  )
  class(out) <- c("random_block", "list")

  out
}

#' @rdname prior_random
#' @export
random_term <- random_block

#' @rdname prior_random
#' @param name optional allocation label. Required when another allocation uses
#'   this allocation as `parent`.
#' @param parent optional `allocation_ref()` object selecting a component of an
#'   earlier named allocation. Child allocations inherit that component's SD
#'   budget and must not specify `sd`.
#' @param components allocation target type. `"block"` allocates across named
#'   random-effect blocks or symbolic intermediate components. `"sd"` expands
#'   one selected heterogeneous block to its resolved SD leaves.
#' @param scale variance normalization. `"total_variance"` uses
#'   `sd_child = sd_parent * sqrt(w)`. `"mean_variance"` uses
#'   `sd_child = sd_parent * sqrt(K * w)` and is intended for
#'   `components = "sd"` heterogeneity tests.
#' @export
random_variance_allocation <- function(terms = NULL, sd = NULL,
                                       allocation = NULL,
                                       name = NULL,
                                       parent = NULL,
                                       components = c("block", "sd"),
                                       scale = c("total_variance", "mean_variance")){

  check_char(terms, "terms", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  if(!is.null(terms) && anyDuplicated(terms)){
    stop("'terms' in random_variance_allocation() must be unique.", call. = FALSE)
  }
  .bt_check_random_allocation_term_labels(terms)
  components <- match.arg(components)
  scale <- match.arg(scale)
  if(identical(scale, "mean_variance") && !identical(components, "sd")){
    stop("'scale = \"mean_variance\"' is supported only with components = \"sd\".", call. = FALSE)
  }
  if(identical(components, "sd") && (is.null(terms) || length(terms) != 1L)){
    stop("'components = \"sd\"' requires exactly one random-effect block in 'terms'.", call. = FALSE)
  }
  .bt_check_random_allocation_ref(parent, allow_NULL = TRUE)
  if(is.null(parent)){
    .bt_check_random_sd_prior(sd, allow_NULL = FALSE)
  }else{
    if(!is.null(sd)){
      stop("Child variance allocations inherit their SD budget from 'parent' and must not specify 'sd'.", call. = FALSE)
    }
  }
  if(is.null(allocation)){
    if(identical(components, "block") && !is.null(terms) && length(terms) < 2L){
      stop(
        "A default Dirichlet allocation prior can be used only when at least two 'terms' are supplied.",
        call. = FALSE
      )
    }
    if(identical(components, "block") && !is.null(terms)){
      allocation <- prior("dirichlet", list(alpha = rep(1, length(terms))))
    }
  }else{
    .bt_check_random_allocation_prior(allocation)
  }
  check_char(name, "name", allow_NULL = TRUE, allow_NA = FALSE)
  if(!is.null(name) && !grepl("^[A-Za-z][A-Za-z0-9_]*$", name)){
    stop("'name' must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
  }

  if(identical(components, "block") && !is.null(terms) &&
     !is.null(allocation) && length(terms) != allocation$parameters[["K"]]){
    stop(
      "The Dirichlet allocation dimension must match the number of targeted random-effect terms.",
      call. = FALSE
    )
  }

  out <- list(
    terms      = terms,
    sd         = sd,
    allocation = allocation,
    name       = name,
    parent     = parent,
    components = components,
    scale      = scale
  )
  class(out) <- c("random_variance_allocation", "list")

  out
}

#' @rdname prior_random
#' @param allocation_name allocation name referenced by `allocation_ref()`.
#' @param component component label inside the referenced allocation. Named
#'   `terms` supply component labels; unnamed `terms` use sanitized term names.
#' @export
allocation_ref <- function(allocation_name, component){

  check_char(allocation_name, "allocation_name", allow_NA = FALSE)
  check_char(component, "component", allow_NA = FALSE)
  if(!grepl("^[A-Za-z][A-Za-z0-9_]*$", allocation_name)){
    stop("'allocation_name' must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
  }
  if(!grepl("^[A-Za-z][A-Za-z0-9_]*$", component)){
    stop("'component' must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
  }

  out <- list(
    allocation = allocation_name,
    component  = component
  )
  class(out) <- c("random_allocation_ref", "list")

  out
}

#' @rdname prior_random
#' @param structure covariance structure label. Accepted labels are `"ID"`,
#'   `"DIAG"`, `"CS"`, `"HCS"`, `"AR"` / `"AR1"`, `"HAR"`, `"CAR"`, and
#'   `"UN"` / `"US"`; matching is case-insensitive.
#' @param eta LKJ concentration parameter used when `cor` is omitted for
#'   unstructured covariance. Values above one regularize correlations toward
#'   zero; `eta = 1` is uniform over correlation matrices.
#' @param rho_scale scale for scalar correlation priors. `"fisher_z"` is the
#'   default, `"logit"` uses a logit transform of the valid raw-rho interval,
#'   and `"rho"` uses the raw correlation parameter.
#' @param backend LKJ backend, either the compiled JAGS module (`"module"`) or
#'   the pure-JAGS syntax fallback (`"syntax"`).
#' @export
random_covariance <- function(structure = NULL, sd = NULL, cor = NULL,
                              rho = NULL, eta = NULL,
                              rho_scale = c("fisher_z", "rho", "logit"),
                              backend = c("module", "syntax")){

  explicit_fields <- character()
  if(!missing(structure)) explicit_fields <- c(explicit_fields, "structure")
  if(!missing(sd)) explicit_fields <- c(explicit_fields, "sd")
  if(!missing(cor)) explicit_fields <- c(explicit_fields, "cor")
  if(!missing(rho)) explicit_fields <- c(explicit_fields, "rho")
  if(!missing(rho_scale)) explicit_fields <- c(explicit_fields, "rho_scale")
  if(!missing(backend)) explicit_fields <- c(explicit_fields, "backend")

  rho_scale <- match.arg(rho_scale)
  backend <- match.arg(backend)
  .bt_check_random_sd_prior(sd, allow_NULL = TRUE)

  if(!is.null(cor) && !is.null(rho)){
    stop("'cor' and 'rho' cannot both be supplied.", call. = FALSE)
  }
  if(!is.null(eta)){
    if(!is.null(cor) || !is.null(rho)){
      stop("'eta' cannot be supplied together with explicit 'cor' or 'rho'.", call. = FALSE)
    }
    check_real(eta, "eta", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
    if(!is.finite(eta)){
      stop("'eta' must be finite.", call. = FALSE)
    }
    if(is.null(cor)){
      cor <- prior_lkj(eta = eta, backend = backend)
    }
  }
  if(!is.null(cor)){
    .bt_check_random_cor_prior(cor)
  }
  if(!is.null(rho)){
    .bt_check_random_rho_prior(rho)
  }
  if(!is.null(structure)){
    check_char(structure, "structure", allow_NA = FALSE)
    structure <- .bt_random_covariance_normalize(structure)
  }

  out <- list(
    structure = structure,
    sd        = sd,
    cor       = cor,
    rho       = rho,
    rho_scale = rho_scale,
    backend   = backend
  )
  if(!is.null(out$cor) && "backend" %in% explicit_fields){
    out$cor$backend <- backend
  }
  class(out) <- c("random_covariance", "list")
  attr(out, "explicit_fields") <- explicit_fields

  if(!is.null(structure)){
    .bt_validate_random_covariance_for_structure(
      out,
      structure = tolower(structure),
      label = "random_covariance()"
    )
  }

  out
}

#' @rdname prior_random
#' @param include_correlation whether to monitor the correlation matrix for
#'   `prior_lkj()`.
#' @param include_primitives whether to monitor LKJ primitive beta coordinates
#'   for `prior_lkj()`.
#' @export
prior_lkj <- function(eta = 1, backend = c("module", "syntax"),
                      include_correlation = TRUE,
                      include_primitives = FALSE){

  backend_explicit <- !missing(backend)
  backend <- match.arg(backend)
  check_real(eta, "eta", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if(!is.finite(eta)){
    stop("'eta' must be finite.", call. = FALSE)
  }
  check_bool(include_correlation, "include_correlation")
  check_bool(include_primitives, "include_primitives")

  out <- list(
    eta = eta,
    backend = backend,
    include_correlation = include_correlation,
    include_primitives = include_primitives
  )
  class(out) <- c("prior_lkj", "list")
  attr(out, "explicit_fields") <- if(backend_explicit) "backend" else character()

  out
}

#' @rdname prior_random
#' @param latent whether to monitor standardized latent random effects.
#'   Required for bridge sampling and coefficient reconstruction.
#' @param coefficients whether to monitor realized group-level coefficients.
#'   This is convenient for inspection but can be memory intensive.
#' @param correlation whether to monitor correlation matrices.
#' @param lkj_primitives whether to monitor LKJ primitive beta coordinates.
#' @export
random_monitor <- function(latent = TRUE, coefficients = FALSE,
                           correlation = TRUE, lkj_primitives = FALSE){

  check_bool(latent, "latent")
  check_bool(coefficients, "coefficients")
  check_bool(correlation, "correlation")
  check_bool(lkj_primitives, "lkj_primitives")

  out <- list(
    latent = latent,
    coefficients = coefficients,
    correlation = correlation,
    lkj_primitives = lkj_primitives
  )
  class(out) <- c("random_monitor", "list")

  out
}

#' @rdname prior_random
#' @param allow whether new levels are allowed during prediction.
#' @param method new-level prediction method.
#' @export
random_new_levels <- function(allow = FALSE, method = c("zero", "sample")){

  check_bool(allow, "allow")
  method <- match.arg(method)

  out <- list(allow = allow, method = method)
  class(out) <- c("random_new_levels", "list")

  out
}

is.prior_random <- function(x){
  inherits(x, "prior_random")
}

.bt_check_prior_random <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!is.prior_random(x)){
    stop("'prior_random' must be created with prior_random().", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_prior_blocks <- function(blocks){

  if(length(blocks) == 0L){
    return(invisible(TRUE))
  }
  if(is.null(names(blocks)) || any(!nzchar(names(blocks)))){
    stop("Random-effect block overrides in '...' must be named.", call. = FALSE)
  }
  if(anyDuplicated(names(blocks))){
    stop("Random-effect block override names must be unique.", call. = FALSE)
  }
  for(i in seq_along(blocks)){
    if(!inherits(blocks[[i]], "random_block")){
      stop("Random-effect block override '", names(blocks)[i], "' must be created with random_block().", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.bt_check_random_sd_prior <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!is.prior(x)){
    stop("Random-effect SD specifications must be prior objects.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_covariance <- function(x){

  if(!inherits(x, "random_covariance")){
    stop("'covariance' must be created with random_covariance().", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_cor_prior <- function(x){

  if(!inherits(x, "prior_lkj")){
    stop("Correlation priors currently must be created with prior_lkj().", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_rho_prior <- function(x){

  if(!is.prior(x)){
    stop("Scalar correlation priors must be prior objects.", call. = FALSE)
  }
  if(is.prior.none(x)){
    stop("Scalar correlation priors cannot use prior_none().", call. = FALSE)
  }
  if(.bt_random_prior_has_non_scalar_family(x)){
    stop("Scalar correlation priors must be ordinary scalar prior objects.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_random_prior_is_factor_family <- function(x){

  is.prior.factor(x) ||
    inherits(x, "prior.factor_mixture") ||
    inherits(x, "prior.factor_spike_and_slab")
}

.bt_random_prior_has_non_scalar_family <- function(x){

  if(.bt_random_prior_is_factor_family(x) ||
     is.prior.simplex(x) ||
     is.prior.vector(x) ||
     (is.prior.discrete(x) && !is.prior.point(x)) ||
     is.prior.PET(x) ||
     is.prior.PEESE(x) ||
     is.prior.weightfunction(x) ||
     inherits(x, "prior.bias_mixture")){
    return(TRUE)
  }
  if(is.prior.spike_and_slab(x) || is.prior.mixture(x)){
    return(any(vapply(x, .bt_random_prior_has_non_scalar_family, logical(1))))
  }

  FALSE
}

.bt_check_random_monitor <- function(x){

  if(!inherits(x, "random_monitor")){
    stop("'monitor' must be created with random_monitor().", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_new_levels <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!inherits(x, "random_new_levels")){
    stop("'new_levels' must be created with random_new_levels().", call. = FALSE)
  }
  if(isTRUE(x$allow)){
    stop(
      "New-level random-effect prediction is not implemented yet; leave 'new_levels' at random_new_levels(allow = FALSE).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_check_random_allocation <- function(x){

  if(is.null(x)){
    return(invisible(TRUE))
  }

  allocations <- .bt_random_allocation_list(x)
  allocation_names <- names(allocations)
  if(!is.null(allocation_names) && any(nzchar(allocation_names))){
    if(any(!nzchar(allocation_names))){
      stop("Variance allocation list names must be either all named or all unnamed.", call. = FALSE)
    }
    if(anyDuplicated(allocation_names)){
      stop("Variance allocation list names must be unique.", call. = FALSE)
    }
    bad <- !grepl("^[A-Za-z][A-Za-z0-9_]*$", allocation_names)
    if(any(bad)){
      stop("Variance allocation list names must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
    }
  }
  if(length(allocations) == 0L){
    stop("'allocation' must contain at least one random_variance_allocation() object.", call. = FALSE)
  }
  for(i in seq_along(allocations)){
    allocation <- allocations[[i]]
    if(!inherits(allocation, "random_variance_allocation")){
      stop("'allocation' entries must be created with random_variance_allocation().", call. = FALSE)
    }
    components <- allocation$components
    if(is.null(components)){
      components <- "block"
    }
    scale <- allocation$scale
    if(is.null(scale)){
      scale <- "total_variance"
    }
    check_char(components, "components", allow_values = c("block", "sd"), allow_NA = FALSE)
    check_char(scale, "scale", allow_values = c("total_variance", "mean_variance"), allow_NA = FALSE)
    if(identical(scale, "mean_variance") && !identical(components, "sd")){
      stop("'scale = \"mean_variance\"' is supported only with components = \"sd\".", call. = FALSE)
    }
    .bt_check_random_allocation_term_labels(allocation$terms)
    .bt_check_random_allocation_ref(allocation$parent, allow_NULL = TRUE)
    if(is.null(allocation$parent)){
      .bt_check_random_sd_prior(allocation$sd, allow_NULL = FALSE)
    }else if(!is.null(allocation$sd)){
      stop("Child variance allocations inherit their SD budget from 'parent' and must not specify 'sd'.", call. = FALSE)
    }
    if(!is.null(allocation$allocation)){
      .bt_check_random_allocation_prior(allocation$allocation)
    }
    if(identical(components, "sd") &&
       (is.null(allocation$terms) || length(allocation$terms) != 1L)){
      stop("'components = \"sd\"' requires exactly one random-effect block in 'terms'.", call. = FALSE)
    }
    if(identical(components, "block") &&
       !is.null(allocation$terms) &&
       !is.null(allocation$allocation) &&
       length(allocation$terms) != allocation$allocation$parameters[["K"]]){
      stop(
        "The Dirichlet allocation dimension must match the number of targeted random-effect terms.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_check_random_allocation_ref <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!inherits(x, "random_allocation_ref")){
    stop("'parent' must be created with allocation_ref().", call. = FALSE)
  }
  check_char(x$allocation, "allocation", allow_NA = FALSE)
  check_char(x$component, "component", allow_NA = FALSE)

  invisible(TRUE)
}

.bt_check_random_allocation_term_labels <- function(terms){

  if(is.null(terms)){
    return(invisible(TRUE))
  }
  labels <- names(terms)
  if(is.null(labels) || all(!nzchar(labels))){
    return(invisible(TRUE))
  }
  if(any(!nzchar(labels))){
    stop("Variance allocation 'terms' must be either all named or all unnamed.", call. = FALSE)
  }
  if(anyDuplicated(labels)){
    stop("Variance allocation term labels must be unique.", call. = FALSE)
  }
  bad <- !grepl("^[A-Za-z][A-Za-z0-9_]*$", labels)
  if(any(bad)){
    stop("Variance allocation term labels must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_allocation_prior <- function(x){

  if(!is.prior.simplex(x) || !identical(x$distribution, "dirichlet")){
    stop(
      "Variance allocation priors must be Dirichlet simplex priors created with prior('dirichlet', list(alpha = ...)).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_allocation_list <- function(x){

  if(is.null(x)){
    return(list())
  }
  if(inherits(x, "random_variance_allocation")){
    return(list(x))
  }
  if(!is.list(x)){
    stop("'allocation' must be a random_variance_allocation() object or a list of them.", call. = FALSE)
  }

  x
}

.bt_random_covariance_normalize <- function(x){

  x <- toupper(x)
  aliases <- c(
    US = "US",
    UN = "US",
    DIAG = "DIAG",
    ID = "ID",
    CS = "CS",
    HCS = "HCS",
    AR = "AR1",
    AR1 = "AR1",
    CAR = "CAR",
    HAR = "HAR"
  )
  if(!x %in% names(aliases)){
    stop("Unknown random-effect covariance structure '", x, "'.", call. = FALSE)
  }

  unname(aliases[[x]])
}

.bt_random_prior_block_names <- function(prior_random){

  .bt_check_prior_random(prior_random)
  names(prior_random$blocks)
}

.bt_random_prior_for_block <- function(prior_random, block_name){

  .bt_check_prior_random(prior_random)
  check_char(block_name, "block_name", allow_NA = FALSE)

  block <- list(
    sd         = prior_random$sd,
    covariance = prior_random$covariance,
    monitor    = prior_random$monitor,
    new_levels = prior_random$new_levels,
    allocation = NULL,
    terms      = NULL
  )

  if(block_name %in% names(prior_random$blocks)){
    override <- prior_random$blocks[[block_name]]
    for(field in names(override)){
      if(!is.null(override[[field]])){
        if(identical(field, "covariance")){
          block[[field]] <- .bt_random_merge_covariance(block[[field]], override[[field]])
        }else{
          block[[field]] <- override[[field]]
        }
      }
    }
  }

  class(block) <- c("random_block", "list")
  block
}

.bt_random_block_lkj_prior <- function(block){

  if(is.null(block$covariance)){
    block$covariance <- random_covariance()
  }
  .bt_check_random_covariance(block$covariance)
  cor_prior <- block$covariance$cor
  if(is.null(cor_prior)){
    cor_prior <- prior_lkj(eta = 1, backend = block$covariance$backend)
  }

  cor_prior
}

.bt_random_merge_covariance <- function(base, override){

  .bt_check_random_covariance(base)
  .bt_check_random_covariance(override)

  out <- base
  override_fields <- .bt_random_covariance_explicit_fields(override)
  for(field in c("structure", "sd", "cor", "rho")){
    if(!is.null(override[[field]])){
      out[[field]] <- override[[field]]
    }
  }
  if("rho_scale" %in% override_fields){
    out$rho_scale <- override$rho_scale
  }
  if("backend" %in% override_fields){
    out$backend <- override$backend
    if(!is.null(out$cor)){
      out$cor$backend <- override$backend
    }
  }
  if(!is.null(override$rho)){
    out["cor"] <- list(NULL)
  }
  if(!is.null(override$cor)){
    if("backend" %in% override_fields){
      out$backend <- override$backend
      out$cor$backend <- override$backend
    }else if("backend" %in% .bt_random_covariance_explicit_fields(base) &&
             !.bt_prior_lkj_backend_explicit(override$cor)){
      out$backend <- base$backend
      out$cor$backend <- base$backend
    }else{
      out$backend <- out$cor$backend
    }
    out["rho"] <- list(NULL)
  }
  if(!is.null(out$structure)){
    structure <- tolower(.bt_random_covariance_normalize(out$structure))
    .bt_random_validate_explicit_covariance_override(
      override = override,
      override_fields = override_fields,
      structure = structure
    )
    if(.bt_random_structure_uses_no_correlation(structure)){
      out["cor"] <- list(NULL)
      out["rho"] <- list(NULL)
    }else if(.bt_random_structure_uses_lkj(structure)){
      out["rho"] <- list(NULL)
    }else if(.bt_random_structure_uses_scalar_rho(structure)){
      out["cor"] <- list(NULL)
    }
  }

  class(out) <- c("random_covariance", "list")
  attr(out, "explicit_fields") <- unique(c(
    .bt_random_covariance_explicit_fields(base),
    override_fields
  ))
  out
}

.bt_random_validate_explicit_covariance_override <- function(override,
                                                            override_fields,
                                                            structure){

  if(.bt_random_structure_uses_no_correlation(structure)){
    if("cor" %in% override_fields && !is.null(override$cor)){
      stop(
        "Block covariance override supplies an LKJ correlation prior, but structure '",
        structure, "' has no correlation parameter.",
        call. = FALSE
      )
    }
    if("rho" %in% override_fields && !is.null(override$rho)){
      stop(
        "Block covariance override supplies a scalar correlation prior, but structure '",
        structure, "' has no correlation parameter.",
        call. = FALSE
      )
    }
  }else if(.bt_random_structure_uses_lkj(structure)){
    if("rho" %in% override_fields && !is.null(override$rho)){
      stop(
        "Block covariance override supplies a scalar correlation prior, but structure '",
        structure, "' uses an LKJ correlation prior.",
        call. = FALSE
      )
    }
  }else if(.bt_random_structure_uses_scalar_rho(structure)){
    if("cor" %in% override_fields && !is.null(override$cor)){
      stop(
        "Block covariance override supplies an LKJ correlation prior, but structure '",
        structure, "' uses a scalar correlation prior.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_random_covariance_explicit_fields <- function(x){

  fields <- attr(x, "explicit_fields")
  if(is.null(fields)){
    fields <- setdiff(names(x)[!vapply(x, is.null, logical(1))], c("rho_scale", "backend"))
    if(!is.null(x$rho)){
      fields <- c(fields, "rho_scale")
    }
    if(!is.null(x$cor)){
      fields <- c(fields, "backend")
    }
  }

  fields
}

.bt_prior_lkj_backend_explicit <- function(x){

  fields <- attr(x, "explicit_fields")
  if(!is.null(fields)){
    return("backend" %in% fields)
  }

  !is.null(x$backend) && !identical(x$backend, "module")
}

.bt_random_structure_uses_lkj <- function(structure){

  structure %in% c("us")
}

.bt_random_structure_uses_scalar_rho <- function(structure){

  structure %in% c("cs", "hcs", "ar1", "car", "har")
}

.bt_random_structure_uses_no_correlation <- function(structure){

  structure %in% c("diag", "id")
}

.bt_validate_random_covariance_for_structure <- function(covariance, structure,
                                                         label = "random-effect covariance"){

  .bt_check_random_covariance(covariance)
  check_char(structure, "structure", allow_NA = FALSE)
  structure <- tolower(.bt_random_covariance_normalize(structure))

  if(.bt_random_structure_uses_no_correlation(structure)){
    if(!is.null(covariance$cor)){
      stop(
        label, " supplies an LKJ correlation prior, but structure '",
        structure, "' has no correlation parameter.",
        call. = FALSE
      )
    }
    if(!is.null(covariance$rho)){
      stop(
        label, " supplies a scalar correlation prior, but structure '",
        structure, "' has no correlation parameter.",
        call. = FALSE
      )
    }
  }else if(.bt_random_structure_uses_lkj(structure)){
    if(!is.null(covariance$rho)){
      stop(
        label, " supplies a scalar correlation prior, but structure '",
        structure, "' uses an LKJ correlation prior.",
        call. = FALSE
      )
    }
  }else if(.bt_random_structure_uses_scalar_rho(structure)){
    if(!is.null(covariance$cor)){
      stop(
        label, " supplies an LKJ correlation prior, but structure '",
        structure, "' uses a scalar correlation prior.",
        call. = FALSE
      )
    }
  }else{
    stop("Unsupported random-effect covariance structure '", structure, "'.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_validate_random_block_for_structure <- function(block, structure,
                                                    block_name = NULL){

  if(!inherits(block, "random_block")){
    stop("'block' must be a random_block object.", call. = FALSE)
  }
  label <- "random-effect block"
  if(!is.null(block_name) && nzchar(block_name)){
    label <- paste0("random-effect block '", block_name, "'")
  }

  .bt_validate_random_covariance_for_structure(
    block$covariance,
    structure = structure,
    label = label
  )
  .bt_check_random_new_levels(block$new_levels, allow_NULL = TRUE)
  .bt_check_random_allocation(block$allocation)

  invisible(TRUE)
}
