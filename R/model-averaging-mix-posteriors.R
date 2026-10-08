#' @title Model-average posterior distributions
#'
#' @description Model-averages posterior distributions based on
#' a list of models, vector of parameters, and a list of
#' indicators of the null or alternative hypothesis models
#' for each parameter.
#'
#' @details For simplex parameters such as Dirichlet priors, absent model
#' parameters are not filled with an implicit spike at zero because the all-zero
#' vector is off the simplex. Every model with positive prior or posterior
#' probability must provide an explicit compatible simplex prior or point prior
#' on the simplex.
#'
#' Mixture draws are aligned across parameters by restarting from one shared
#' sampling seed and requiring identical posterior model probabilities for every
#' requested parameter. That holds for unconditional averaging even when
#' \code{is_null_list} differs by parameter. Conditional averaging
#' (\code{conditional = TRUE}) with different null indicators per parameter
#' produces different \code{post_probs} and therefore different mixture
#' allocations; \code{mix_posteriors()} hard-errors in that case. Mix each
#' conditional parameter in a separate call instead.
#'
#' For an ordered total with a declared inclusion event, conditional averaging
#' also conditions within each model. Posterior model probabilities are
#' multiplied by the event's fraction of fitted draws; prior probabilities are
#' multiplied by its declared prior probability. Resampling uses eligible
#' original draw indices, and conditional prior contexts preserve the original
#' fitted source specifications. A parameterized total without an inclusion
#' event, including a Bernoulli total, retains every state. A compound event
#' must be applied once to full aligned sources; separately conditioned ordered
#' parameters cannot be combined as joint draws by this function.
#'
#' Positive model declarations remain active when their visible probability
#' rounds to zero. Numeric mixed draws remain usable. A prior law or atom law
#' requiring model probabilities unavailable at full precision is explicitly
#' refused with \code{BayesTools_formula_measure_unavailable} and the cause
#' \code{numerical_model_probability_unavailable}; structured model/log
#' diagnostics are retained in [posterior_metadata()]. Descriptive
#' \code{marginal_posterior(use_formula = FALSE, prior_samples = FALSE)} remains
#' available. An ordered projection with no allocated rows for a positive
#' model requires a complete declared atom certificate; otherwise that measure
#' has cause \code{structural_target_law_unavailable}. Independently known
#' support and valid sibling measures remain available.
#'
#' @param seed integer specifying seed for sampling posteriors for
#' model averaging. The caller's random-number state (\code{.Random.seed} and
#' \code{RNGkind()}) is restored afterwards. Defaults to \code{NULL}, which
#' draws the shared sampling seed from the caller's random-number stream (one
#' draw).
#' @param n_samples number of samples to be drawn for the model-averaged
#' posterior distribution
#' @inheritParams ensemble_inference
#'
#' @return \code{mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [ensemble_inference] [BayesTools_ensemble_tables] [as_mixed_posteriors]
#'
#' @name mix_posteriors
#' @export
mix_posteriors <- function(model_list, parameters, is_null_list,
                           conditional = FALSE, seed = NULL,
                           n_samples = 10000,
                           on_failure = c("error", "drop", "zero")){

  # A seeded call leaves the caller's random-number state as it found it.
  if(!is.null(seed)){
    rng_state <- .bt_rng_state()
    on.exit(.bt_rng_restore(rng_state), add = TRUE)
  }

  # check input
  check_list(model_list, "model_list")
  check_char(parameters, "parameters", check_length = FALSE)
  check_list(is_null_list, "is_null_list", check_length = length(parameters))
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  on_failure <- match.arg(on_failure)
  sapply(model_list, function(m)check_list(m, "model_list:model", check_names = c("fit", "marglik", "prior_weights"), all_objects = TRUE, allow_other = TRUE))
  if(!all(sapply(model_list, function(m) inherits(m[["fit"]], what = "runjags")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "stanfit")) | sapply(model_list, function(m)inherits(m[["fit"]], what = "null_model"))))
    stop("model_list:fit must contain 'runjags' or 'rstan' models")
  for(m in model_list){
    if(inherits(m[["fit"]], "runjags")){
      .bt_require_fit_contract(m[["fit"]], "model_list:fit")
    }
  }
  if(!all(unlist(sapply(model_list, function(m) sapply(attr(m[["fit"]], "prior_list"), function(p) is.prior(p))))))
    stop("model_list:priors must contain 'BayesTools' priors")
  sapply(model_list, function(m) check_real(m[["prior_weights"]], "model_list:prior_weights", lower = 0))


  # extract the object
  fits           <- lapply(model_list, function(m) m[["fit"]])
  priors         <- lapply(model_list, function(m) attr(m[["fit"]], "prior_list"))

  inference <- ensemble_inference(
    model_list,
    parameters,
    is_null_list,
    conditional,
    on_failure = on_failure
  )

  # Use one shared sampling seed so formula terms are drawn from aligned
  # posterior rows across parameters.
  if(!is.null(seed)){
    set.seed(seed)
    common_sample_seed <- seed
  }else{
    # One draw from the caller's stream; the sampling below is scoped.
    common_sample_seed <- sample(.Machine$integer.max, 1)
    rng_state <- .bt_rng_state()
    on.exit(.bt_rng_restore(rng_state), add = TRUE)
  }

  .mix_posteriors_assert_aligned_post_probs(inference, parameters)
  ordered_conditions <- list()
  if(conditional){
    for(parameter in parameters){
      parameter_priors <- lapply(priors,function(prior) prior[[parameter]])
      if(!any(vapply(parameter_priors,is.prior.ordered,logical(1)))) next
      condition <- .mix_posteriors_ordered_condition(fits,parameter_priors,parameter,inference[[parameter]])
      inference[[parameter]] <- .model_probability_inference_set(inference[[parameter]],
        condition$prior, condition$posterior)
      ordered_conditions[[parameter]] <- condition
    }
    .mix_posteriors_assert_aligned_post_probs(inference,parameters)
    if(length(parameters)>1L && any(vapply(ordered_conditions,function(condition){
      isTRUE(condition$has_nested_event)
    },logical(1)))){
      stop("Joint conditional mixed draws are unavailable for independently requested ordered events. Mix one conditional parameter per call or use 'conditional = FALSE'; a compound 'AND' event must be applied once to full aligned sources.",call.=FALSE)
    }
  }

  out <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter    <- parameters[p]
    temp_inference    <- inference[[temp_parameter]]
    prior_pair <- .model_probability_inference_get(temp_inference, "prior")
    posterior_pair <- .model_probability_inference_get(temp_inference, "posterior")
    temp_priors       <- lapply(priors, function(p) p[[temp_parameter]])

    if(all(sapply(temp_priors, is.null))){
      stop(
        "The parameter '", temp_parameter,
        "' is not available in any model prior list.",
        call. = FALSE
      )
    }

    if(any(sapply(temp_priors, is.prior.weightfunction)) && all(sapply(temp_priors, is.prior.weightfunction) | sapply(temp_priors, .is_prior_weightfunction_null) | sapply(temp_priors, is.null))){
      # weightfunctions:

      # replace missing priors with default prior: none
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior_none()
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        declaration <- prior_pair$declaration
        declaration$model_indices <- as.integer(i)
        temp_priors[[i]] <- .set_prior_model_probability(temp_priors[[i]], prior_pair$probabilities[i], prior_pair$logs[i], declaration)
      }

      out[[temp_parameter]] <- .mix_posteriors.weightfunction(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples, posterior_pair = posterior_pair)

    }else if(any(sapply(temp_priors, is.prior.factor)) && all(sapply(temp_priors, is.prior.factor) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # factor priors

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        declaration <- prior_pair$declaration
        declaration$model_indices <- as.integer(i)
        temp_priors[[i]] <- .set_prior_model_probability(temp_priors[[i]], prior_pair$probabilities[i], prior_pair$logs[i], declaration)
      }

      out[[temp_parameter]] <- .mix_posteriors.factor(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples,
        ordered_condition=ordered_conditions[[temp_parameter]], posterior_pair = posterior_pair)

    }else if(any(sapply(temp_priors, is.prior.vector)) && all(sapply(temp_priors, is.prior.vector) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # vector priors:

      temp_priors <- .mix_posteriors_validate_simplex_priors(
        temp_priors,
        temp_parameter,
        temp_inference$prior_probs,
        temp_inference$post_probs,
        log_prior_probs = prior_pair$logs, log_post_probs = posterior_pair$logs
      )

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        declaration <- prior_pair$declaration
        declaration$model_indices <- as.integer(i)
        temp_priors[[i]] <- .set_prior_model_probability(temp_priors[[i]], prior_pair$probabilities[i], prior_pair$logs[i], declaration)
      }

      out[[temp_parameter]] <- .mix_posteriors.vector(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples, posterior_pair = posterior_pair)

    }else if(all(sapply(temp_priors, is.prior.simple) | sapply(temp_priors, is.prior.point) | sapply(temp_priors, is.null))){
      # simple priors:

      # replace missing priors with default prior: spike(0)
      for(i in 1:length(fits)){
        if(is.null(temp_priors[[i]])){
          temp_priors[[i]] <- prior("spike", parameters = list("location" = 0))
        }
      }

      # replace prior odds with the corresponding prior model odds
      for(i in seq_along(temp_priors)){
        declaration <- prior_pair$declaration
        declaration$model_indices <- as.integer(i)
        temp_priors[[i]] <- .set_prior_model_probability(temp_priors[[i]], prior_pair$probabilities[i], prior_pair$logs[i], declaration)
      }

      out[[temp_parameter]] <- .mix_posteriors.simple(fits, temp_priors, temp_parameter, temp_inference$post_probs, common_sample_seed, n_samples, posterior_pair = posterior_pair)

    }else{
      stop("The posterior samples cannot be mixed: unsupported mixture of prior distributions.")
    }
    out[[temp_parameter]] <- .model_probability_mixed_measures(out[[temp_parameter]],
      prior_pair, posterior_pair, temp_priors, temp_parameter)

    # the fitted coordinate and label parts of every column
    if(is.null(.bt_meta_get(out[[temp_parameter]], "quantities"))){
      out[[temp_parameter]] <- .mix_posteriors_set_quantities(
        samples   = out[[temp_parameter]],
        parameter = temp_parameter,
        priors    = temp_priors,
        fits      = fits
      )
    }

    # add formula relevant information
    if(!is.null(unique(unlist(lapply(temp_priors, attr, which = "parameter", exact = TRUE))))){
      class(out[[temp_parameter]]) <- c(class(out[[temp_parameter]]), "mixed_posteriors.formula")
      out[[temp_parameter]] <- .bt_meta_set(out[[temp_parameter]], "formula_parameter", unique(unlist(lapply(temp_priors, attr, which = "parameter", exact = TRUE))))
      out[[temp_parameter]] <- .bt_meta_set(out[[temp_parameter]], "log_intercept", .mixed_posteriors_formula_log_intercept(
        fits,
        .bt_meta_get(out[[temp_parameter]], "formula_parameter")
      ))
    }

  }

  class(out) <- c(class(out), "mixed_posteriors")
  formula_owners <- parameters[vapply(out, function(x) !is.null(.bt_meta_get(x, "formula_parameter")), logical(1))]
  declared_owners <- unique(unlist(lapply(fits, function(fit){
    unlist(lapply(attr(fit, "formula_design", exact = TRUE), function(design){
      intersect(names(design$prior_list), paste0(design$parameter, "_", design$model_terms))
    }), use.names = FALSE)
  }), use.names = FALSE))
  formula_owners <- intersect(formula_owners, declared_owners)
  if(length(formula_owners)){
    reference <- out[[formula_owners[[1L]]]]
    model_index <- .bt_draws_model_component(reference)
    row_index <- .bt_meta_get(reference, "draw_index")
    records <- vector("list", length(fits))
    blocks <- vector("list", length(fits))
    for(model in seq_along(fits)){
      declared_priors <- lapply(formula_owners, function(owner) attr(out[[owner]], "prior_list", exact = TRUE)[[model]])
      names(declared_priors) <- formula_owners
      all_zero <- all(vapply(declared_priors, .posterior_atoms_is_zero_point, logical(1)))
      if(inherits(fits[[model]], "null_model") && all_zero){
        zero_columns <- unique(unlist(lapply(formula_owners, function(owner){
          .posterior_atoms_coefficient_columns(out[[owner]], owner)
        }), use.names = FALSE))
        records[[model]] <- list(prior_list = declared_priors, formula_scale = list(), required = character(),
          prior_probability = inference[[formula_owners[[1L]]]]$prior_probs[[model]],
          prior_log_probability = attr(inference[[formula_owners[[1L]]]], "log_prior_probs", exact = TRUE)[[model]],
          condition_event = NULL, absent = TRUE, zero_coordinates = zero_columns,
          gate_plan = list(components = matrix(rep(1L, length(formula_owners)), nrow = 1L,
            dimnames = list(NULL, formula_owners)), probabilities = 1, index = 1L,
            draw_index = 1L, model_mixture = FALSE), eligible_n = 1L)
        blocks[[model]] <- matrix(numeric(), sum(model_index == model), 0L,
          dimnames = list(NULL, character()))
        next
      }
      original <- .extract_posterior_samples(fits[[model]], as_list = FALSE)
      original_rows <- seq_len(nrow(original))
      condition <- ordered_conditions[[formula_owners[[1L]]]]
      if(!is.null(condition)) original_rows <- condition$eligible[[model]]
      fitted_state <- .bt_formula_state_new(fits[[model]], original[original_rows, , drop = FALSE],
        formula_owners, original_rows, condition_event = if(!is.null(condition))
          .condition_event(attr(fits[[model]], "prior_list", exact = TRUE), formula_owners[[1L]], "AND"))
      if(is.null(fitted_state)) .bt_formula_transform_stop(
        "A model in the formula mixture has no retained declaration owner.", reason = "missing_multiplier_state", model = model)
      records[[model]] <- fitted_state$models[[1L]]
      missing_owners <- setdiff(formula_owners, names(attr(fits[[model]], "prior_list", exact = TRUE)))
      if(length(missing_owners)){
        if(!all(vapply(declared_priors[missing_owners], .posterior_atoms_is_zero_point, logical(1)))){
          .bt_formula_transform_stop("Missing formula owners have no declared zero contribution.",
            reason = "missing_multiplier_state", model = model, missing = missing_owners)
        }
        records[[model]]$prior_list <- c(records[[model]]$prior_list, declared_priors[missing_owners])
        records[[model]]$zero_owners <- missing_owners
        records[[model]]$zero_coordinates <- unique(unlist(lapply(missing_owners, function(owner){
          .posterior_atoms_coefficient_columns(out[[owner]], owner)
        }), use.names = FALSE))
      }
      records[[model]]$prior_probability <- inference[[formula_owners[[1L]]]]$prior_probs[[model]]
      records[[model]]$prior_log_probability <- attr(inference[[formula_owners[[1L]]]], "log_prior_probs", exact = TRUE)[[model]]
      requested <- which(model_index == model)
      block_rows <- match(row_index[requested], original_rows)
      if(anyNA(block_rows)) stop("Formula mixture draw rows do not belong to their eligible model population.", call. = FALSE)
      blocks[[model]] <- fitted_state$values[block_rows, , drop = FALSE]
    }
    columns <- unique(unlist(lapply(blocks, colnames), use.names = FALSE))
    values <- matrix(NA_real_, length(model_index), length(columns), dimnames = list(NULL, columns))
    for(model in seq_along(blocks)) values[model_index == model, colnames(blocks[[model]])] <- blocks[[model]]
    state <- list(schema_version = 2L, models = records, model = as.integer(model_index),
      draw_index = as.integer(row_index), values = values,
      posterior_model_probabilities = as.numeric(inference[[formula_owners[[1L]]]]$post_probs),
      posterior_log_model_probabilities = attr(inference[[formula_owners[[1L]]]], "log_post_probs", exact = TRUE))
    attr(state, "model_probability_declaration") <- attr(inference[[formula_owners[[1L]]]], "model_probability_declaration", exact = TRUE)
    out <- .bt_formula_state_attach(out, state)
  }
  if(length(parameters)==1L && isTRUE(ordered_conditions[[parameters[[1L]]]]$has_nested_event)){
    out <- .bt_meta_set(out,"prior_context",.bt_meta_get(out[[parameters[[1L]]]],"prior_context"))
    out <- .bt_meta_set(out,"condition",.bt_meta_get(out[[parameters[[1L]]]],"condition"))
  }
  return(out)
}

# Mixed draws of one parameter with the quantities of their columns: label
# parts and fitted coordinates from a model prior that defines it, and catalog
# quantities from the first fit whose parameter map contains it.
.mix_posteriors_set_quantities <- function(samples, parameter, priors, fits){

  defined <- vapply(priors, function(prior){
    !is.null(prior) && !is.prior.point(prior) && !is.prior.none(prior)
  }, logical(1))
  if(!any(defined)){
    defined <- !vapply(priors, is.null, logical(1))
  }
  prior <- priors[[which(defined)[1L]]]
  catalog <- NULL
  for(i in which(defined)){
    if(inherits(fits[[i]], "BayesTools_fit")){
      catalog <- parameter_catalog(fits[[i]])
      break
    }
  }

  .bt_mixed_set_quantities(
    samples,
    parameter = parameter,
    prior     = prior,
    columns   = if(is.null(dim(samples))) parameter else colnames(samples),
    catalog   = catalog
  )
}

# Persisted log(intercept) flag of the fitted formula; NULL when no fit stores
# a formula design and NA when the fitted designs disagree.
.mixed_posteriors_formula_log_intercept <- function(fits, formula_parameter){

  if(length(formula_parameter) != 1L){
    return(NULL)
  }

  flags <- unlist(lapply(fits, function(fit){
    design <- attr(fit, "formula_design", exact = TRUE)
    if(!is.list(design) || !is.list(design[[formula_parameter]])){
      return(NULL)
    }
    isTRUE(design[[formula_parameter]][["log_intercept"]])
  }), use.names = FALSE)

  if(length(flags) == 0L){
    return(NULL)
  }
  if(length(unique(flags)) != 1L){
    return(NA)
  }

  flags[[1L]]
}

.mix_posteriors_assert_aligned_post_probs <- function(inference, parameters){

  if(length(parameters) <= 1L){
    return(invisible(TRUE))
  }

  reference <- as.numeric(inference[[parameters[[1L]]]][["post_probs"]])
  reference_pair <- .model_probability_inference_get(inference[[parameters[[1L]]]], "posterior")
  for(parameter in parameters[-1L]){
    other <- as.numeric(inference[[parameter]][["post_probs"]])
    other_pair <- .model_probability_inference_get(inference[[parameter]], "posterior")
    if(!identical(reference_pair$logs, other_pair$logs) ||
       !identical(reference_pair$declaration$model_indices, other_pair$declaration$model_indices) ||
       !isTRUE(all.equal(reference, other, tolerance = 0, check.attributes = FALSE))){
      stop(
        "mix_posteriors() requires identical posterior model probabilities ",
        "across parameters so mixture draws stay aligned. Differing ",
        "post_probs usually arise from conditional = TRUE with different ",
        "null indicators per parameter. Mix each conditional parameter ",
        "separately, or use unconditional averaging.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.mix_posteriors_is_dirichlet_simplex <- function(prior){
  is.prior.simplex(prior) && identical(prior[["distribution"]], "dirichlet")
}

.mix_posteriors_canonicalize_simplex_point <- function(prior, K){

  if(!is.prior.point(prior)){
    return(NULL)
  }

  location <- prior$parameters[["location"]]
  if(is.numeric(location) && length(location) == 1L){
    location <- rep(location, K)
  }
  if(!is.numeric(location) || length(location) != K){
    return(NULL)
  }

  canonical <- tryCatch(
    .canonicalize_simplex(
      location,
      name = "point-prior location",
      diagnostics = TRUE
    ),
    error = function(e) NULL
  )
  if(is.null(canonical)){
    return(NULL)
  }

  prior$parameters[["location"]] <- canonical$values
  attr(prior, "simplex_canonicalization") <- canonical$diagnostics
  prior
}

.mix_posteriors_validate_simplex_priors <- function(priors, parameter, prior_probs, post_probs,
                                                  log_prior_probs = log(prior_probs), log_post_probs = log(post_probs)){

  simplex <- vapply(priors, .mix_posteriors_is_dirichlet_simplex, logical(1))
  if(!any(simplex)){
    return(priors)
  }

  K <- unique(vapply(priors[simplex], function(prior){
    prior$parameters[["K"]]
  }, numeric(1)))
  if(length(K) != 1L){
    stop(
      "The simplex parameter '", parameter,
      "' cannot be mixed across Dirichlet priors with different dimensions.",
      call. = FALSE
    )
  }

  contributing <- is.finite(log_prior_probs) | is.finite(log_post_probs)
  missing <- vapply(priors, is.null, logical(1))
  if(any(missing & contributing)){
    stop(
      "The simplex parameter '", parameter,
      "' has no implicit spike-at-zero null. Every model with positive prior ",
      "or posterior probability must provide an explicit compatible Dirichlet ",
      "simplex prior or point prior on the simplex.",
      call. = FALSE
    )
  }

  canonical_points <- lapply(priors, function(prior){
    if(is.null(prior)){
      return(NULL)
    }
    .mix_posteriors_canonicalize_simplex_point(prior, K)
  })
  valid_point <- !vapply(canonical_points, is.null, logical(1))
  compatible <- missing | simplex | valid_point

  if(any(!compatible & contributing)){
    stop(
      "The simplex parameter '", parameter,
      "' can only mix Dirichlet simplex priors with the same dimension or ",
      "explicit point priors on the simplex.",
      call. = FALSE
    )
  }

  for(i in which(valid_point)){
    priors[[i]] <- canonical_points[[i]]
  }

  priors
}

.mix_posteriors.simple         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000, posterior_pair = NULL){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  if(is.null(posterior_pair)) posterior_pair <- .model_probability_pair(post_probs, log(post_probs), "posterior", "raw")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.simple) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of simple priors")

  # gather and check compatibility of prior distributions
  priors_info <- lapply(priors, function(p){
    if(is.prior.point(p) | is.prior.none(p)){
      return(FALSE)
    }else{
      list(
        "interaction"       = .is_prior_interaction(p),
        "interaction_terms" = attr(p, "interaction_terms")
      )
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!vapply(
    priors_info,
    function(i) isTRUE(all.equal(i, priors_info[[1]])),
    logical(1)
  ))){
    stop("non-matching prior factor type specifications")
  }else if(length(priors_info) != 0){
    priors_info <- priors_info[[1]]
  }



  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  samples <- NULL
  draw_index <- NULL
  model_component <- NULL

  # mix samples
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

    if(is.prior.point(priors[[i]])){
      # not sampling the priors as the samples would be already transformed
      samples <- c(samples, rep(priors[[i]]$parameters[["location"]], length(temp_ind)))
    }else{
      samples <- c(samples, model_samples[temp_ind, parameter])
    }

    draw_index <- c(draw_index, temp_ind)
    model_component <- c(model_component, rep(i, length(temp_ind)))
  }

  samples <- unname(samples)
  samples <- .bt_meta_set(samples, "draw_index", draw_index)
  samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  attr(samples, "interaction")       <- if(length(priors_info) == 0) FALSE else priors_info[["interaction"]]
  attr(samples, "interaction_terms") <- priors_info[["interaction_terms"]]
  samples <- .posterior_support_set_from_prior_list(samples, priors)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors,
      post_probs,
      posterior_pair = posterior_pair,
      n_columns = 1L,
      column_names = parameter
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.mix_posteriors.vector         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000,
                                           column_names = NULL, posterior_pair = NULL){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  if(is.null(posterior_pair)) posterior_pair <- .model_probability_pair(post_probs, log(post_probs), "posterior", "raw")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.vector) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of vector priors")

  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
  }

  # prepare output objects
  K <- unique(sapply(priors[sapply(priors, is.prior.vector)], function(p) p$parameters[["K"]]))
  if(length(K) != 1)
    stop("all vector priors must be of the same length")

  samples    <- matrix(nrow = 0, ncol = K)
  draw_index <- NULL
  model_component <- NULL

  # mix samples
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

    if(is.prior.point(priors[[i]])){
      # not sampling the priors as the samples would be already transformed
      samples <- rbind(
        samples,
        matrix(
          rep(priors[[i]]$parameters[["location"]], times = length(temp_ind)),
          nrow = length(temp_ind),
          ncol = K,
          byrow = TRUE
        )
      )
    }else if(K == 1){
      samples <- rbind(samples, matrix(model_samples[temp_ind, parameter], nrow = length(temp_ind), ncol = K))
    }else{
      samples <- rbind(samples, model_samples[temp_ind, paste0(parameter,"[",1:K,"]")])
    }

    draw_index <- c(draw_index, temp_ind)
    model_component <- c(model_component, rep(i, length(temp_ind)))
  }

  rownames(samples) <- NULL
  colnames(samples) <- if(is.null(column_names)) paste0(parameter,"[",1:K,"]") else column_names
  samples <- .bt_meta_set(samples, "draw_index", draw_index)
  samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  samples <- .posterior_support_set_columns_from_prior_list(samples, priors)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors,
      post_probs,
      posterior_pair = posterior_pair,
      n_columns = K,
      column_names = colnames(samples)
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.mix_posteriors.factor         <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000,
                                           ordered_condition = NULL, posterior_pair = NULL){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  if(is.null(posterior_pair)) posterior_pair <- .model_probability_pair(post_probs, log(post_probs), "posterior", "raw")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.factor) | sapply(priors, is.prior.point)))
    stop("'priors' must be a list of factor priors")
  priors <- .complete_factor_metadata_prior_list(priors)

  # check the prior levels
  levels <- unique(sapply(priors[sapply(priors, is.prior.factor)], .get_prior_factor_levels))
  if(length(levels) != 1)
    stop("all factor priors must be of the same number of levels")

  # gather and check compatibility of prior distributions
  priors_info <- lapply(priors, function(p){
    if(is.prior.factor(p)){
      return(list(
        "levels"            = .get_prior_factor_levels(p),
        "level_names"       = .get_prior_factor_level_names(p),
        "interaction"       = .is_prior_interaction(p),
        "interaction_terms" = attr(p, "interaction_terms"),
        "term_components"   = attr(p, "term_components"),
        "factor_terms"      = attr(p, "factor_terms"),
        "factor_contrasts"  = attr(p, "factor_contrasts"),
        "factor_design"     = attr(p, "factor_design"),
        "factor_cell_names" = attr(p, "factor_cell_names"),
        "treatment"         = is.prior.treatment(p),
        "independent"       = is.prior.independent(p),
        "orthonormal"       = is.prior.orthonormal(p),
        "meandif"           = is.prior.meandif(p),
        "ordered"           = is.prior.ordered(p)
      ))
    }else if(is.prior.point(p) | is.prior.none(p)){
      return(FALSE)
    }else{
      stop("unsupported prior type")
    }
  })
  priors_info <- priors_info[!sapply(priors_info, isFALSE)]
  if(length(priors_info) >= 2 && any(!vapply(
    priors_info,
    function(i) isTRUE(all.equal(i, priors_info[[1]])),
    logical(1)
  ))){
    stop("non-matching prior factor type specifications")
  }
  if(length(priors_info) != 0){
    priors_info <- priors_info[[1]]
  }


  if(priors_info[["ordered"]]){

    if(!is.null(seed)){
      set.seed(seed)
    }

    ordered_prior <- priors[[which(vapply(priors, is.prior.ordered, logical(1)))[1]]]
    coefficient_names <- .JAGS_prior_factor_names(parameter, ordered_prior)
    indicator_name    <- paste0(.prior_ordered_total_name(parameter), "_indicator")
    samples    <- matrix(nrow = 0, ncol = levels)
    draw_index <- NULL
    model_component <- NULL
    # per-draw component of the ordered total of models whose total has a
    # spike at zero (NA for the other models); formula-level atoms split by it
    total_indicator <- NULL
    source_models <- lapply(seq_along(priors),function(i){
      .bt_ordered_source_model(priors[[i]],parameter,if(inherits(fits[[i]],"BayesTools_fit")) fits[[i]])
    })
    source_samples <- vector("list", length(priors))

    sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
    for(i in seq_along(fits)[sample_counts > 0]){

      if(inherits(fits[[i]], "null_model")){
        model_samples <- matrix()
      }else if(inherits(fits[[i]], "runjags")){
        model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
        if(!is.matrix(model_samples)){
          model_samples <- matrix(model_samples, ncol = 1)
          colnames(model_samples) <- fits[[i]]$monitor
        }
      }else if(inherits(fits[[i]], "stanfit")){
        .check_rstan()
        model_samples <- .extract_stan(fits[[i]])
      }

      eligible <- if(is.null(ordered_condition)) seq_len(nrow(model_samples)) else ordered_condition$eligible[[i]]
      temp_ind <- eligible[sample(length(eligible),sample_counts[i],replace=TRUE)]
      if(.mix_posteriors_ordered_total_has_indicator(priors[[i]]) && !indicator_name %in% colnames(model_samples)){
        .mix_posteriors_stop_missing_total_indicator(parameter,indicator_name)
      }
      source_samples[[i]] <- .bt_ordered_source_rows(source_models[[i]], model_samples, temp_ind)

      if(is.prior.point(priors[[i]])){
        samples <- rbind(
          samples,
          matrix(
            rep(priors[[i]]$parameters[["location"]], times = length(temp_ind)),
            nrow = length(temp_ind),
            ncol = levels,
            byrow = TRUE
          )
        )
      }else{
        temp_names <- .JAGS_prior_factor_names(parameter, priors[[i]])
        samples <- rbind(samples, model_samples[temp_ind, temp_names, drop = FALSE])
      }

      temp_total_indicator <- rep(NA_integer_, length(temp_ind))
      if(.mix_posteriors_ordered_total_has_indicator(priors[[i]])){
        if(!indicator_name %in% colnames(model_samples)){
          .mix_posteriors_stop_missing_total_indicator(parameter, indicator_name)
        }
        temp_total_indicator <- .bt_component_from_indicator(
          priors[[i]]$total,
          model_samples[temp_ind, indicator_name]
        )
      }

      draw_index <- c(draw_index, temp_ind)
      model_component <- c(model_component, rep(i, length(temp_ind)))
      total_indicator <- c(total_indicator, temp_total_indicator)
    }

    rownames(samples) <- NULL
    # the first ordered coordinate is a level cell; later ones are contrast
    # coefficients `{j}`, never bracketed positions
    colnames(samples) <- .bt_label_prior_column_names(parameter, ordered_prior)
    samples <- .bt_meta_set(samples, "draw_index", draw_index)
    samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
    samples <- .bt_meta_set(samples, "ordered_source",
      .bt_ordered_source_new(parameter, source_models, source_samples, model_component, draw_index,
        posterior_pair = posterior_pair))
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    if(any(!is.na(total_indicator))){
      samples <- .bt_meta_set(samples, "ordered_total_component", total_indicator)
    }
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["treatment"]]){

    if(levels == 1){

      samples <- .mix_posteriors.simple(fits, priors, parameter, post_probs, seed, n_samples, posterior_pair = posterior_pair)

      draw_index <- .bt_meta_get(samples, "draw_index")
      model_component <- .bt_meta_get(samples, "component")

      samples <- matrix(samples, ncol = 1)

    }else{

      # keep the same seed across levels
      if(is.null(seed)){
        seed <- sample(666666, 1)
      }

      samples <- lapply(1:levels, function(i) .mix_posteriors.simple(fits, priors, paste0(parameter, "[", i, "]"), post_probs, seed, n_samples, posterior_pair = posterior_pair))

      draw_index <- .bt_meta_get(samples[[1]], "draw_index")
      model_component <- .bt_meta_get(samples[[1]], "component")

      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    # level cells from the term's design (a full-rank interaction such as
    # `~ g + g:x` includes the first level); cumulative increments of an
    # interaction with an ordered factor are contrast coefficients `{j}`
    factor_prior <- priors[vapply(priors, is.prior.factor, logical(1))][[1]]
    colnames(samples) <- .bt_label_prior_column_names(parameter, factor_prior)
    samples <- .bt_meta_set(samples, "draw_index", draw_index)
    samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["independent"]]){

    if(levels == 1){

      samples <- .mix_posteriors.simple(fits, priors, parameter, post_probs, seed, n_samples)

      draw_index <- .bt_meta_get(samples, "draw_index")
      model_component <- .bt_meta_get(samples, "component")

      samples <- matrix(samples, ncol = 1)

    }else{

      # keep the same seed across levels
      if(is.null(seed)){
        seed <- sample(666666, 1)
      }

      samples <- lapply(1:levels, function(i) .mix_posteriors.simple(fits, priors, paste0(parameter, "[", i, "]"), post_probs, seed, n_samples))

      draw_index <- .bt_meta_get(samples[[1]], "draw_index")
      model_component <- .bt_meta_get(samples[[1]], "component")

      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    factor_prior <- priors[vapply(priors, is.prior.factor, logical(1))][[1]]
    colnames(samples) <- .bt_label_prior_column_names(parameter, factor_prior)
    samples <- .bt_meta_set(samples, "draw_index", draw_index)
    samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- priors
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(priors_info[["orthonormal"]] | priors_info[["meandif"]]){

    for(i in seq_along(priors)){
      if(is.prior.factor(priors[[i]])){
        priors[[i]]$parameters[["K"]] <- levels
      }
    }

    factor_prior <- priors[vapply(priors, is.prior.factor, logical(1))][[1]]
    samples <- .mix_posteriors.vector(
      fits, priors, parameter, post_probs, seed, n_samples,
      column_names = .bt_label_prior_column_names(parameter, factor_prior), posterior_pair = posterior_pair
    )
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")            <- priors_info[["levels"]]
  attr(samples, "level_names")       <- priors_info[["level_names"]]
  attr(samples, "interaction")       <- if(length(priors_info) == 0) FALSE else priors_info[["interaction"]]
  attr(samples, "interaction_terms") <- priors_info[["interaction_terms"]]
  attr(samples, "term_components")   <- priors_info[["term_components"]]
  attr(samples, "factor_terms")      <- priors_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- priors_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- priors_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- priors_info[["factor_cell_names"]]
  attr(samples, "treatment")         <- priors_info[["treatment"]]
  attr(samples, "independent")       <- priors_info[["independent"]]
  attr(samples, "orthonormal")       <- priors_info[["orthonormal"]]
  attr(samples, "meandif")           <- priors_info[["meandif"]]
  attr(samples, "ordered")           <- priors_info[["ordered"]]
  if(isTRUE(priors_info[["ordered"]])){
    ordered_prior <- priors[[which(vapply(priors, is.prior.ordered, logical(1)))[1]]]
    attr(samples, "ordered_metadata") <- attr(ordered_prior, "ordered_metadata")
  }

  if(isTRUE(priors_info[["treatment"]]) || isTRUE(priors_info[["independent"]])){
    factor_support <- .posterior_support_from_prior_list(priors)
    if(!is.null(factor_support) && !is.null(colnames(samples))){
      samples <- .bt_meta_set(samples, "support", stats::setNames(
        rep(list(factor_support), ncol(samples)),
        colnames(samples)
      ))
    }
  }

  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors,
      post_probs,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      posterior_pair = posterior_pair,
      exclusion_probabilities = if(isTRUE(priors_info[["ordered"]])){
        .mix_posteriors_ordered_exclusion_probabilities(fits, priors, parameter, post_probs,
          log_post_probs=if(!is.null(posterior_pair)) posterior_pair$logs else log(post_probs),
          eligible=if(!is.null(ordered_condition)) ordered_condition$eligible)
      }
    )
  )
  if(isTRUE(priors_info[["ordered"]])){
    samples <- .bt_ordered_source_semantics(samples,diag(ncol(samples)),colnames(samples))
    if(!is.null(ordered_condition)){
      source <- .bt_meta_get(samples,"ordered_source")
      source$conditioning <- ordered_condition[c("prior_probs","post_probs","log_prior_probs","log_post_probs","prior_fractions","prior_log_fractions","posterior_fractions")]
      source$conditioning$model_probability_declaration <- list(prior=ordered_condition$prior$declaration,posterior=ordered_condition$posterior$declaration)
      samples <- .bt_meta_set(samples,"ordered_source",source)
      if(isTRUE(ordered_condition$has_nested_event)){
        samples <- .bt_meta_set(samples,"prior_context",.mix_posteriors_ordered_prior_context(
          priors,parameter,ordered_condition))
        ordered_model <- which(vapply(source$models,function(spec) is.null(spec$parameterization),logical(1)))[1L]
        event <- .condition_event(setNames(list(source$models[[ordered_model]]$prior),parameter),parameter)
        samples <- .condition_event_set_attributes(samples,event)
      }
    }
  }

  return(samples)
}

.mix_posteriors_ordered_condition <- function(fits, priors, parameter, inference){

  eligible <- options <- vector("list",length(priors))
  posterior <- prior_probability <- numeric(length(priors))
  prior_log_probability <- rep(-Inf, length(priors))
  nested <- FALSE
  for(i in seq_along(priors)){
    prior <- priors[[i]]
    if(is.null(prior) || is.prior.none(prior) || .posterior_atoms_is_zero_point(prior) || inherits(fits[[i]],"null_model")) next
    draws <- if(inherits(fits[[i]],"stanfit")) .extract_stan(fits[[i]]) else as.matrix(.fit_to_posterior(fits[[i]]))
    if(is.prior.ordered(prior) && is.prior.mixture(prior$total)){
      nested <- TRUE
      own <- setNames(list(prior),parameter)
      event <- .condition_event(own,parameter)
      mask <- .condition_event_posterior_mask(event,own,draws)
      options[[i]] <- .condition_event_model_options(own,event)
      prior_probability[[i]] <- options[[i]]$event_probability
      prior_log_probability[[i]] <- options[[i]]$log_event_probability
    }else{
      mask <- rep(TRUE,nrow(draws))
      options[[i]] <- list(prior_lists=list(setNames(list(prior),parameter)),weights=1,log_weights=0,event_probability=1,log_event_probability=0)
      prior_probability[[i]] <- 1
      prior_log_probability[[i]] <- 0
    }
    eligible[[i]] <- which(mask)
    posterior[[i]] <- mean(mask)
  }
  prior <- .model_probability_inference_get(inference, "prior")
  post <- .model_probability_inference_get(inference, "posterior")
  if(!any(is.finite(prior$logs) & is.finite(prior_log_probability))) stop("Conditional inference requires at least one non-null model.",call.=FALSE)
  if(!any(is.finite(post$logs) & posterior > 0)){
    .bt_ordered_stop(paste0("The conditional posterior of '",parameter,
      "' is unavailable: no fitted draw lies in its declared inclusion event. Obtain more upstream posterior draws or use 'conditional = FALSE'."))
  }
  prior <- .model_probability_condition(prior, prior_log_probability, "conditional_prior", prior_probability)
  post <- .model_probability_condition(post, log(posterior), "conditional_posterior", posterior)
  list(eligible=eligible,prior_options=options,prior_fractions=prior_probability,
    prior_log_fractions=prior_log_probability, posterior_fractions=posterior,has_nested_event=nested,
    prior_probs=prior$probabilities, post_probs=post$probabilities,
    log_prior_probs=prior$logs,log_post_probs=post$logs,prior=prior,posterior=post)
}

.mix_posteriors_ordered_prior_context <- function(priors, parameter, condition){

  prior_lists <- list()
  weights <- numeric()
  log_weights <- numeric()
  for(i in seq_along(priors)){
    options <- condition$prior_options[[i]]
    if(is.null(options) || !is.finite(condition$log_prior_probs[[i]])) next
    prior_lists <- c(prior_lists,options$prior_lists)
    weights <- c(weights,condition$prior_probs[[i]] * options$weights)
    log_weights <- c(log_weights,condition$log_prior_probs[[i]] + options$log_weights)
  }
  ordered <- priors[vapply(priors,is.prior.ordered,logical(1))][[1L]]
  event <- .condition_event(setNames(list(ordered),parameter),parameter)
  pair <- .model_probability_prior(weights, log_weights, "conditional_prior")
  structure(list(schema_version=2L,linear_weight_space="coefficient",prior_list=setNames(list(priors),parameter),
    column_names=.JAGS_prior_factor_names(parameter,ordered),formula_scale=NULL,transforms=list(),
    conditional=parameter,conditional_rule="AND",condition_event=event,condition_key=event$condition_key,
    prior_lists=prior_lists,model_weights=pair$probabilities,model_log_weights=pair$logs,
    model_probability_declaration=pair$declaration,n_grid=.prior_linear_density_default_grid(),
    tail_prob=.prior_linear_density_tail_prob()),class="prior_density_conditional_context")
}
# Posterior probability that an ordered total (spike-and-slab, or a mixture
# with point(0) components) is excluded within each model, from the fitted
# total-prior indicator.
.mix_posteriors_ordered_exclusion_probabilities <- function(fits, priors, parameter, post_probs,eligible=NULL,
                                                           log_post_probs=log(post_probs)){

  indicator_name <- paste0(.prior_ordered_total_name(parameter), "_indicator")
  vapply(seq_along(priors), function(i){
    if(!is.finite(log_post_probs[i]) || !.mix_posteriors_ordered_total_has_indicator(priors[[i]])){
      return(0)
    }
    model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
    if(!is.matrix(model_samples)){
      model_samples <- matrix(model_samples, ncol = 1)
      colnames(model_samples) <- fits[[i]]$monitor
    }
    if(!indicator_name %in% colnames(model_samples)){
      .mix_posteriors_stop_missing_total_indicator(parameter, indicator_name)
    }
    component <- .bt_component_from_indicator(priors[[i]]$total,model_samples[,indicator_name])
    if(!is.null(eligible)) component <- component[eligible[[i]]]
    .posterior_atoms_ordered_exclusion(priors[[i]]$total,component)
  }, numeric(1))
}
# Whether a model's ordered prior has a total with a within-model spike at zero
# whose posterior share is read from the fitted total-prior indicator.
.mix_posteriors_ordered_total_has_indicator <- function(prior){

  is.prior.ordered(prior) &&
    !.posterior_atoms_is_ordered_zero_total(prior) &&
    .posterior_atoms_ordered_total_has_spike(prior$total)
}
.mix_posteriors_stop_missing_total_indicator <- function(parameter, indicator_name){

  .bt_stop_refit_required(
    "The fitted samples for ordered factor '", parameter,
    "' do not contain the required total-prior indicator '",
    indicator_name, "'. Refit the model with this package version."
  )
}
.mix_posteriors.weightfunction <- function(fits, priors, parameter, post_probs, seed = NULL, n_samples = 10000, posterior_pair = NULL){

  # check input
  check_list(fits, "fits")
  check_list(priors, "priors", check_length = length(fits))
  check_char(parameter, "parameter")
  check_real(post_probs, "post_probs", lower = 0, upper = 1, check_length = length(fits))
  if(is.null(posterior_pair)) posterior_pair <- .model_probability_pair(post_probs, log(post_probs), "posterior", "raw")
  check_real(seed, "seed", allow_NULL = TRUE)
  check_int(n_samples, "n_samples")
  if(!all(sapply(fits, inherits, what = "runjags") | sapply(fits, inherits, what = "stanfit") | sapply(fits, inherits, what = "null_model")))
    stop("'fits' must be a list of 'runjags' or 'rstan' models")
  if(!all(sapply(priors, is.prior.weightfunction) | sapply(priors, .is_prior_weightfunction_null)))
    stop("'priors' must be a list of weightfunction priors or point(1)/none null priors")


  # When seed is NULL, continue using the current RNG state so mix_posteriors()
  # can coordinate draws across parameters.
  if(!is.null(seed)){
    set.seed(seed)
  }

  # obtain mapping for the weight coefficients
  omega_mapping <- weightfunctions_mapping(priors)
  omega_cuts    <- weightfunctions_mapping(priors, cuts_only = TRUE)
  omega_names   <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))

  # prepare output objects
  samples    <- matrix(nrow = 0, ncol = length(omega_cuts) - 1)
  draw_index <- NULL
  model_component <- NULL

  # mix samples
  sample_counts <- .posterior_mixture_sample_counts(post_probs, n_samples)
  for(i in seq_along(fits)[sample_counts > 0]){

    # obtain posterior samples
    if(inherits(fits[[i]], "null_model")){
      # deal with a possibility of completely null model
      model_samples <- matrix()
    }else if(inherits(fits[[i]], "runjags")){
      model_samples <- .extract_posterior_samples(fits[[i]], as_list = FALSE)
      if(!is.matrix(model_samples)){
        # deal with automatic coercion into a vector in case of a single predictor
        model_samples <- matrix(model_samples, ncol = 1)
        colnames(model_samples) <- fits[[i]]$monitor
      }
    }else if(inherits(fits[[i]], "stanfit")){
      .check_rstan()
      model_samples <- .extract_stan(fits[[i]])
    }

    # sample indexes
    temp_ind <- sample(nrow(model_samples), sample_counts[i], replace = TRUE)

    if(.is_prior_weightfunction_null(priors[[i]])){
      samples <- rbind(samples, matrix(1, ncol = length(omega_cuts) - 1, nrow = length(temp_ind)))
    }else{
      samples <- rbind(samples, model_samples[temp_ind, paste0("omega[",omega_mapping[[i]],"]")])
    }

    draw_index <- c(draw_index, temp_ind)
    model_component <- c(model_component, rep(i, length(temp_ind)))
  }

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  # each weight mixes different fitted coordinates across the models
  samples <- .bt_meta_set(samples, "quantities", .bt_verbatim_quantities(omega_names))
  samples <- .bt_meta_set(samples, "draw_index", draw_index)
  samples <- .model_probability_component_set(samples, model_component, priors, post_probs, posterior_pair)
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- priors
  samples <- .posterior_support_set_weightfunction_columns(
    samples,
    priors,
    .weightfunction_mapping_info(priors)
  )
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      priors,
      post_probs,
      posterior_pair = posterior_pair,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      null_location = 1,
      point_locations = .model_probability_weightfunction_points(priors, omega_mapping, ncol(samples))
    )
  )
  samples <- .posterior_weightfunction_declarations(samples, priors, post_probs)
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}
