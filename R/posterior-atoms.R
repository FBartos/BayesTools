#' @title Declare posterior point-mass metadata
#'
#' @description Creates exact posterior measure metadata used to distinguish a
#' continuous posterior from one containing atoms. Supply \code{NULL} to
#' explicitly declare that the posterior has no point masses.
#'
#' @param point_masses optional data frame or list with \code{x} and
#' \code{mass} entries.
#' @param source short label describing the source of the declaration.
#'
#' @return A posterior atom metadata object for the \code{atoms} metadata
#' of posterior draws (\code{posterior_metadata(x, "atoms") <- }). Posterior
#' point masses are read only from these metadata: precomputed posterior
#' densities ([posterior_density_attribute()]) describe the continuous part
#' and carry no point masses.
#' BayesTools producers can additionally declare positional scalar
#' \code{marginals}, named in the full draw-column order. Scalar extraction
#' prefers these declarations; a continuous joint measure can have a point
#' marginal, for example the last level of an ordered point-total prior.
#' When a producer can certify only scalar marginals, the object records
#' \code{joint_declared = FALSE} and retains its typed joint refusal in
#' \code{joint_unavailable}. Its empty joint table is dimensional storage.
#' It does not declare an atom-free joint law. Scalar selection uses only
#' certified marginals; whole-law inference retains the joint refusal.
#'
#' @export
posterior_atom_attribute <- function(point_masses = NULL, source = "user"){

  check_char(source, "source", check_length = 1L, allow_NA = FALSE)
  point_masses <- .posterior_atoms_point_mass_table(point_masses)
  if(is.null(point_masses)){
    stop("Posterior 'point_masses' metadata is invalid.", call. = FALSE)
  }

  locations <- matrix(point_masses$x, ncol = 1L)
  colnames(locations) <- "value"
  .posterior_atoms_new(
    locations = locations,
    mass = point_masses$mass,
    source = source,
    declared = TRUE
  )
}

#' @title Whether posterior draws are declared atom-free
#'
#' @description Reads the declared posterior atoms of BayesTools posterior
#' draws (the \code{atoms} metadata, see [posterior_metadata()]) and reports
#' whether the draws are structurally atom-free.
#'
#' @param x BayesTools posterior draws, e.g. an element returned by
#' [as_mixed_posteriors()], [marginal_posterior()], or
#' [parameter_mixed_posterior()].
#'
#' @return \code{TRUE} when atoms are declared and all declared scalar marginals
#' and the joint measure have no positive atomic mass;
#' \code{FALSE} when the draws declare a point mass or do not declare their
#' atom status. Plain numeric draws, which carry no metadata, are an error.
#' An explicitly unavailable joint law signals its retained typed condition,
#' even when one or more scalar marginals have valid certificates.
#'
#' @examples
#' draws <- structure(stats::rnorm(10), class = c("marginal_posterior.simple", "marginal_posterior"))
#' posterior_atoms_free(draws)
#' posterior_metadata(draws, "atoms") <- posterior_atom_attribute()
#' posterior_atoms_free(draws)
#'
#' @export
posterior_atoms_free <- function(x){

  .bt_formula_measure_check(x, "atoms")
  if(is.null(.bt_meta_container(x)) &&
     !inherits(x, c("mixed_posteriors", "marginal_posterior"))){
    .bt_draws_stop_plain(
      "'posterior_atoms_free' requires BayesTools posterior draws, not plain numeric draws"
    )
  }
  atoms <- .posterior_atoms_get(x)

  !is.null(atoms) && nrow(atoms$locations) == 0L &&
    (is.null(atoms$marginals) || all(vapply(atoms$marginals, function(marginal){
      !is.null(marginal) && nrow(marginal$locations) == 0L
    }, logical(1))))
}

# A validated point-mass table (data frame with 'x' and 'mass', masses at
# exactly equal locations merged), an empty table for NULL, or NULL when
# 'point_masses' is malformed.
.posterior_atoms_point_mass_table <- function(point_masses){

  empty <- data.frame(x = numeric(), mass = numeric())
  if(is.null(point_masses)){
    return(empty)
  }

  if(is.data.frame(point_masses)){
    if(!all(c("x", "mass") %in% colnames(point_masses))){
      return(NULL)
    }
  }else if(!is.list(point_masses) ||
           is.null(point_masses[["x"]]) || is.null(point_masses[["mass"]])){
    return(NULL)
  }
  x    <- point_masses[["x"]]
  mass <- point_masses[["mass"]]

  if(length(x) != length(mass)){
    return(NULL)
  }
  if(length(x) == 0L){
    return(empty)
  }

  out <- data.frame(
    x    = suppressWarnings(as.numeric(x)),
    mass = suppressWarnings(as.numeric(mass))
  )
  if(any(!is.finite(out[["x"]])) ||
     any(!is.finite(out[["mass"]])) ||
     any(out[["mass"]] <= 0)){
    return(NULL)
  }

  if(nrow(out) > 0L && anyDuplicated(out[["x"]])){
    # merge atoms at exactly equal locations (character keys would round
    # distinct locations to 15 digits)
    unique_x <- unique(out[["x"]])
    index    <- match(out[["x"]], unique_x)
    out <- data.frame(
      x    = unique_x,
      mass = as.numeric(rowsum(out[["mass"]], index, reorder = TRUE))
    )
    out <- out[order(out[["x"]]), , drop = FALSE]
  }
  mass_bound <- .Machine$double.eps * max(8, nrow(out))
  if(nrow(out) > 0L && sum(out[["mass"]]) > 1 + mass_bound){
    return(NULL)
  }
  rownames(out) <- NULL

  return(out)
}

.posterior_atoms_new <- function(locations = NULL, mass = numeric(),
                                 column_names = NULL,
                                 source = "structural",
                                 declared = TRUE,
                                 component_probabilities = NULL,
                                 marginals = NULL,
                                 component_log_probabilities = NULL,
                                 model_probability_declaration = NULL,
                                 joint_declared = TRUE, joint_unavailable = NULL){

  check_bool(declared, "declared", allow_NA = FALSE)
  check_bool(joint_declared, "joint_declared", allow_NA = FALSE)

  if(is.null(locations)){
    n_columns <- if(is.null(column_names)) 1L else length(column_names)
    locations <- matrix(numeric(), nrow = 0L, ncol = n_columns)
  }else if(!is.matrix(locations)){
    locations <- matrix(locations, ncol = 1L)
  }
  if(!is.numeric(locations) || any(!is.finite(locations))){
    stop("Posterior atom locations must be a finite numeric matrix.", call. = FALSE)
  }
  mass <- as.numeric(mass)
  if(length(mass) != nrow(locations) ||
     any(!is.finite(mass)) || any(mass <= 0)){
    stop("Posterior atom masses must be finite, positive, and match the locations.",
         call. = FALSE)
  }
  mass_bound <- .Machine$double.eps * max(8, length(mass))
  if(sum(mass) > 1 + mass_bound){
    stop("Posterior atom masses cannot sum to more than one.", call. = FALSE)
  }
  if(!is.null(component_probabilities)){
    component_probabilities <- as.numeric(component_probabilities)
    if(length(component_probabilities) == 0L ||
       any(!is.finite(component_probabilities)) ||
       any(component_probabilities < 0)){
      stop("Posterior component probabilities must be finite and nonnegative.",
           call. = FALSE)
    }
    probability_sum <- sum(component_probabilities)
    if(probability_sum <= 0){
      stop("Posterior component probabilities must have positive total mass.",
           call. = FALSE)
    }
    if(is.null(model_probability_declaration) && is.null(component_log_probabilities)){
      component_probabilities <- component_probabilities / probability_sum
    }else .model_probability_validate(component_probabilities, component_log_probabilities, model_probability_declaration)
  }
  if(!is.null(column_names)){
    if(length(column_names) != ncol(locations)){
      stop("Posterior atom column names do not match the location matrix.",
           call. = FALSE)
    }
    colnames(locations) <- column_names
  }
  if(!is.null(marginals)){
    if(!is.list(marginals) || length(marginals) != ncol(locations) ||
       is.null(names(marginals)) || !identical(names(marginals), colnames(locations))){
      stop("Posterior scalar marginals must name every location-matrix column in order.", call. = FALSE)
    }
    marginals <- lapply(marginals, function(marginal){
      if(is.null(marginal)) return(NULL)
      marginal <- .posterior_atoms_from_attribute(marginal)
      if(!isTRUE(marginal$joint_declared) || ncol(marginal$locations) != 1L || !is.null(marginal$marginals)){
        stop("Each posterior scalar marginal must declare a single column without nested marginals.", call. = FALSE)
      }
      marginal
    })
  }

  if(!joint_declared){
    if(nrow(locations) != 0L || length(mass) != 0L || is.null(marginals) ||
       is.null(colnames(locations)) || anyNA(colnames(locations)) || any(!nzchar(colnames(locations))) || anyDuplicated(colnames(locations)) ||
       all(vapply(marginals, is.null, logical(1))) ||
       !inherits(joint_unavailable, "condition") ||
       !inherits(joint_unavailable, c("BayesTools_formula_atoms_unavailable", "BayesTools_formula_measure_unavailable"))){
      stop("Partial posterior atoms require certified scalar marginals and a typed joint refusal, with no declared joint mass.", call. = FALSE)
    }
  }else if(!is.null(joint_unavailable)){
    stop("Complete posterior atoms cannot carry a joint-unavailable condition.", call. = FALSE)
  }

  out <- list(
    declared = isTRUE(declared),
    joint_declared = joint_declared,
    locations = locations,
    mass = mass,
    source = source,
    component_probabilities = component_probabilities
  )
  if(!joint_declared) out$joint_unavailable <- joint_unavailable
  if(!is.null(marginals)) out$marginals <- marginals
  if(!is.null(model_probability_declaration) || !is.null(component_log_probabilities)){
    out$component_log_probabilities <- component_log_probabilities
    out$model_probability_declaration <- model_probability_declaration
  }
  class(out) <- c("BayesTools_posterior_atoms", "list")
  out
}

# Validated posterior atoms: NULL when absent (undeclared atom status); a
# malformed declaration stops.
.posterior_atoms_from_attribute <- function(atoms){

  if(is.null(atoms)){
    return(NULL)
  }
  if(!inherits(atoms, "BayesTools_posterior_atoms") || !isTRUE(atoms$declared)){
    stop("Posterior atom metadata must be created with 'posterior_atom_attribute()'.",
         call. = FALSE)
  }

  .posterior_atoms_new(
    locations = atoms$locations,
    mass = atoms$mass,
    column_names = colnames(atoms$locations),
    source = if(is.null(atoms$source)) "unknown" else atoms$source,
    declared = TRUE,
    component_probabilities = atoms$component_probabilities,
    component_log_probabilities = atoms$component_log_probabilities,
    model_probability_declaration = atoms$model_probability_declaration,
    marginals = atoms$marginals,
    joint_declared = if(is.null(atoms$joint_declared)) TRUE else atoms$joint_declared,
    joint_unavailable = atoms$joint_unavailable
  )
}

.posterior_atoms_point_location <- function(prior, n_columns){

  if(is.prior.ordered(prior) && is.prior.point(prior$total) && !.is_prior_expression(prior$total)){
    spec <- .bt_ordered_spec(.prior_ordered_metadata(prior)$parameter_name,prior)
    if(n_columns!=length(spec$coefficient_names)) return(NULL)
    values <- vapply(seq_along(spec$coefficient_names),function(i){
      weights <- rep(0,length(spec$coefficient_names))
      weights[[i]] <- 1
      tensor <- .bt_ordered_tensor(spec,spec$slice_index[[i]],weights)
      if(length(tensor$records) && prior$total$parameters$location!=0) return(NA_real_)
      prior$total$parameters$location * if(length(tensor$records)) 0 else tensor$coefficients[[1L]]
    },numeric(1))
    if(!anyNA(values)) return(values)
    return(NULL)
  }
  if(.posterior_atoms_is_ordered_zero_total(prior)){
    # every ordered coefficient is total x allocation share = 0
    return(rep(0, n_columns))
  }
  if(!is.prior.point(prior)){
    return(NULL)
  }
  location <- prior$parameters[["location"]]
  if(!is.numeric(location) || any(!is.finite(location))){
    return(NULL)
  }
  if(length(location) == 1L){
    return(rep(location, n_columns))
  }
  if(length(location) == n_columns){
    return(as.numeric(location))
  }

  NULL
}

.posterior_atoms_is_zero_point <- function(prior){

  if(!is.prior(prior) || !is.prior.point(prior)){
    return(FALSE)
  }
  location <- prior$parameters[["location"]]
  is.numeric(location) && length(location) == 1L && isTRUE(location == 0)
}

.posterior_atoms_is_ordered_zero_total <- function(prior){

  if(!is.prior.ordered(prior)){
    return(FALSE)
  }
  total <- prior$total
  if(.posterior_atoms_is_zero_point(total)){
    return(TRUE)
  }
  zero_components <- .posterior_atoms_ordered_total_zero_components(total)
  length(zero_components) > 0L && length(zero_components) == length(total)
}

# Mixture components of an ordered total that are a point at zero (the
# total's component index selects them).
.posterior_atoms_ordered_total_zero_components <- function(total){

  if(!is.prior.mixture(total) || is.prior.spike_and_slab(total)){
    return(integer())
  }

  which(vapply(seq_along(total), function(i){
    .posterior_atoms_is_zero_point(total[[i]])
  }, logical(1)))
}

# Whether an ordered total has a spike at zero whose posterior mass is read
# from the fitted total-prior indicator (as the total's component index).
.posterior_atoms_ordered_total_has_spike <- function(total){

  is.prior.spike_and_slab(total) ||
    length(.posterior_atoms_ordered_total_zero_components(total)) > 0L
}

# Probability that an ordered total is exactly zero: the spike of a
# spike-and-slab total, or the point(0) components of a mixture total, from
# the posterior component indices of the total (indices into its component
# list). Without components, the prior probability of the spike.
.posterior_atoms_ordered_exclusion <- function(total, component = NULL){

  if(is.prior.spike_and_slab(total)){
    if(is.null(component)){
      return(1 - mean(.get_spike_and_slab_inclusion(total)))
    }
    component <- .posterior_atoms_check_components(total, component)
    return(mean(vapply(component, .bt_component_is_spike, logical(1), prior = total)))
  }

  zero_components <- .posterior_atoms_ordered_total_zero_components(total)
  if(length(zero_components) == 0L){
    return(0)
  }
  if(is.null(component)){
    prior_weights <- attr(total, "prior_weights", exact = TRUE)
    return(sum(prior_weights[zero_components]) / sum(prior_weights))
  }
  component <- .posterior_atoms_check_components(total, component)
  mean(component %in% zero_components)
}

.posterior_atoms_check_components <- function(prior, component){

  component <- as.integer(component)
  if(length(component) == 0L || anyNA(component) ||
     any(!component %in% seq_along(prior))){
    stop(
      "Posterior component indices must index the components of the ",
      "mixture or spike-and-slab prior.",
      call. = FALSE
    )
  }
  component
}

.posterior_atoms_from_priors <- function(priors, probabilities, n_columns = 1L,
                                         column_names = NULL,
                                         source = "model_probabilities",
                                         null_location = NULL,
                                         exclusion_probabilities = NULL, posterior_pair = NULL, point_locations = NULL){

  if(is.prior(priors) && length(probabilities) == 1L){
    priors <- list(priors)
  }else if(is.prior(priors)){
    priors <- unclass(priors)
  }
  probabilities <- as.numeric(probabilities)
  if(length(priors) != length(probabilities)){
    stop("Prior components and atom probabilities must have the same length.",
         call. = FALSE)
  }
  if(!is.null(exclusion_probabilities) &&
     length(exclusion_probabilities) != length(priors)){
    stop("Within-model exclusion probabilities must match the prior components.",
         call. = FALSE)
  }

  if(!is.null(point_locations) && (length(point_locations) != length(priors) ||
     !all(vapply(point_locations, function(location){
       is.null(location) || (is.numeric(location) && length(location) == n_columns && all(is.finite(location)))
     }, logical(1))))){
    stop("Declared model point locations must align with every model and target column.", call. = FALSE)
  }
  if(!is.null(posterior_pair)) return(tryCatch(
    .model_probability_prior_atoms(priors, posterior_pair, n_columns, column_names,
      source, null_location, exclusion_probabilities, point_locations),
    BayesTools_formula_measure_unavailable = function(condition) condition))

  locations <- matrix(numeric(), nrow = 0L, ncol = n_columns)
  masses <- numeric()
  for(i in seq_along(priors)){
    if(probabilities[i] <= 0){
      next
    }
    location <- if(!is.null(point_locations) && !is.null(point_locations[[i]])) point_locations[[i]] else
      .posterior_atoms_point_location(priors[[i]], n_columns)
    mass <- probabilities[i]
    if(is.null(location) && !is.null(null_location) &&
       .is_prior_weightfunction_null(priors[[i]])){
      location <- rep(null_location, n_columns)
    }
    if(is.null(location) && !is.null(exclusion_probabilities) &&
       exclusion_probabilities[i] > 0){
      # within-model spike at zero (e.g., a spike-and-slab ordered total)
      location <- rep(0, n_columns)
      mass <- probabilities[i] * exclusion_probabilities[i]
    }
    if(!is.null(location)){
      locations <- rbind(locations, location)
      masses <- c(masses, mass)
    }
  }

  .posterior_atoms_new(
    locations = locations,
    mass = masses,
    column_names = column_names,
    source = source,
    declared = TRUE,
    component_probabilities = probabilities,
    marginals = if(!is.null(column_names) &&
      any(vapply(priors, is.prior.weightfunction, logical(1))) &&
      all(vapply(priors, function(prior) is.prior.weightfunction(prior) || .is_prior_weightfunction_null(prior), logical(1)))){
      .posterior_weightfunction_scalar_laws(priors, probabilities, column_names, source)$marginals
    }else NULL
  )
}

.posterior_atoms_set <- function(samples, atoms){

  if(inherits(atoms, c("BayesTools_formula_atoms_unavailable", "BayesTools_formula_measure_unavailable"))){
    columns <- if(is.matrix(samples)) colnames(samples) else attr(samples, "parameter", exact = TRUE)
    ordered <- .bt_meta_get(samples, "ordered_source")
    if(!is.null(ordered)){
      for(column in columns) samples <- .bt_formula_measure_mark(samples, column, "atoms",
        if(is.null(atoms$detail)) atoms$message else atoms$detail, cause = atoms$reason, diagnostics = atoms$diagnostics)
      marginals <- lapply(seq_along(columns), function(i){
        weights <- diag(length(columns))[i, ]
        projection <- .bt_ordered_source_project(ordered, weights)
        .bt_ordered_source_projection_atoms(ordered, projection, columns[[i]], weights)
      })
      for(i in seq_along(marginals)){
        if(inherits(marginals[[i]], "BayesTools_formula_measure_unavailable")){
          condition <- marginals[[i]]
          samples <- .bt_formula_measure_mark(samples, columns[[i]], "atoms", conditionMessage(condition),
            cause = condition$reason, diagnostics = condition$diagnostics)
          marginals[i] <- list(NULL)
        }else if(is.null(marginals[[i]])){
          samples <- .bt_formula_measure_mark(samples, columns[[i]], "atoms", conditionMessage(atoms),
            cause = atoms$reason, diagnostics = atoms$diagnostics)
        }
      }
      names(marginals) <- columns
      if(all(vapply(marginals, is.null, logical(1)))) return(.bt_meta_set(samples, "atoms", NULL))
      return(.bt_meta_set(samples, "atoms", .posterior_atoms_new(column_names = columns, marginals = marginals,
        component_probabilities = ordered$model_probabilities, component_log_probabilities = ordered$model_log_probabilities,
        model_probability_declaration = ordered$model_probability_declaration,
        joint_declared = FALSE, joint_unavailable = atoms)))
    }
    samples <- .bt_meta_set(samples, "atoms", NULL)
    for(column in columns) samples <- .bt_formula_measure_mark(samples, column, "atoms",
      if(is.null(atoms$detail)) atoms$message else atoms$detail, cause = atoms$reason, diagnostics = atoms$diagnostics)
    return(samples)
  }
  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    stop("Cannot attach invalid posterior atom metadata.", call. = FALSE)
  }
  if(is.matrix(samples) && ncol(samples)==ncol(atoms$locations) && !is.null(colnames(samples))){
    atoms <- .posterior_atoms_rename_columns(atoms,colnames(samples))
  }
  samples <- .bt_meta_set(samples, "atoms", atoms)
  samples
}

# Declare each retained scalar omega from its model/branch law. Joint atoms
# remain separate; an unclassified non-omega column has a NULL marginal.
.posterior_weightfunction_scalar_laws <- function(priors, probabilities, columns,
                                                  source = "model_probabilities", posterior_pair = NULL){

  expanded <- .weightfunction_expand_bias_mixture_priors(priors)
  if(length(expanded) != length(probabilities)){
    stop("Weightfunction branches and posterior probabilities do not align.", call. = FALSE)
  }
  if(is.null(posterior_pair)) posterior_pair <- .model_probability_pair(
    probabilities, log(probabilities), "component", "raw")
  .model_probability_validate(probabilities, posterior_pair$logs,
    posterior_pair$declaration, normalized = TRUE)
  # All declarations participate in the global bin map, including a branch
  # with zero posterior probability; mapping never uses sampled occupancy.
  mapped_priors <- lapply(expanded, function(prior) .set_prior_model_weight(prior, 1))
  attr(mapped_priors, "omega_context") <- attr(priors, "omega_context", exact = TRUE)
  context <- .weightfunction_prior_list_context(mapped_priors, merge = FALSE)
  marginals <- supports <- stats::setNames(rep(list(NULL), length(columns)), columns)
  for(column in intersect(columns, context$omega_names)){
    components <- .weightfunction_prior_entry_components(context, match(column, context$omega_names))
    locations <- rep(list(NULL), length(components))
    component_priors <- vector("list", length(components))
    component_supports <- list()
    for(i in seq_along(components)){
      component <- components[[i]]
      location <- if(component$type == "point") component$location else NULL
      if(component$type == "prior" && is.prior.point(component$prior)){
        location <- component$prior$parameters$location
        if(component$scale != "omega") location <- .density.prior_transformation_checked_x(location, "exp")
      }
      if(!is.null(location)) location <- unname(location)
      locations[i] <- list(location)
      component_priors[[i]] <- if(!is.null(location)) prior("point", list(location)) else
        if(component$type == "beta") prior("beta", list(component$alpha, component$beta)) else component$prior
      if(is.finite(posterior_pair$logs[[i]])) component_supports[[length(component_supports) + 1L]] <-
        .posterior_support_from_weightfunction_component(component, source)
    }
    marginal <- .posterior_atoms_from_priors(component_priors, probabilities,
      n_columns = 1L, column_names = column, source = source, posterior_pair = posterior_pair,
      point_locations = locations)
    if(!inherits(marginal, "BayesTools_formula_measure_unavailable")){
      points <- .posterior_atoms_point_mass_table(list(x = marginal$locations[, 1L], mass = marginal$mass))
      point_values <- vapply(locations, function(location){
        if(is.null(location)) NA_real_ else location[[1L]]
      }, numeric(1))
      continuous <- any(is.finite(posterior_pair$logs) & is.na(point_values))
      grouped_logs <- vapply(points$x, function(location){
        .model_probability_log_sum(posterior_pair$logs[!is.na(point_values) & point_values == location])
      }, numeric(1))
      grouped_masses <- exp(grouped_logs)
      lost <- any(grouped_masses == 0 | grouped_masses < .Machine$double.xmin) ||
        (continuous && length(grouped_masses) > 0L && sum(grouped_masses) == 1)
      if(lost){
        marginal <- errorCondition(
          "Posterior atoms are unavailable because positive model mass or a continuous remainder is not representable at full precision.",
          call = NULL, class = c("BayesTools_formula_atoms_unavailable", "BayesTools_formula_measure_unavailable"),
          reason = "numerical_model_probability_unavailable",
          diagnostics = .model_probability_diagnostics(posterior = posterior_pair, stage = posterior_pair$declaration$stage))
      }else{
        if(continuous && sum(points$mass) == 1) points$mass <- grouped_masses
        marginal <- .posterior_atoms_new(matrix(points$x, ncol = 1L), points$mass,
          column_names = column, source = marginal$source,
          component_probabilities = marginal$component_probabilities,
          component_log_probabilities = marginal$component_log_probabilities,
          model_probability_declaration = marginal$model_probability_declaration)
      }
    }
    marginals[column] <- list(marginal)
    supports[column] <- list(.posterior_support_union(component_supports, source))
  }
  list(marginals = marginals, supports = supports)
}

.posterior_weightfunction_declarations <- function(samples, priors, probabilities,
                                                   source = "model_probabilities", posterior_pair = NULL){

  if(!is.matrix(samples) || ncol(samples) == 0L) return(samples)
  unavailable <- .bt_meta_get(samples, "measure_unavailable")
  joint_unavailable <- if(!is.null(unavailable) && any(unavailable$measure == "atoms")){
    tryCatch(.bt_formula_measure_check(samples, "atoms", unavailable$column[unavailable$measure == "atoms"]),
      BayesTools_formula_atoms_unavailable = function(condition) condition)
  }else NULL
  atoms <- .posterior_atoms_get(samples, allow_partial = TRUE)
  if(!is.null(atoms) && !atoms$joint_declared) joint_unavailable <- atoms$joint_unavailable
  if(is.null(posterior_pair)){
    owners <- .bt_meta_get(samples, "model_probabilities")
    if(!is.null(owners)) posterior_pair <- owners$posterior else
      if(!is.null(atoms$model_probability_declaration)) posterior_pair <- list(
        probabilities = atoms$component_probabilities, logs = atoms$component_log_probabilities,
        declaration = atoms$model_probability_declaration)
  }
  declarations <- .posterior_weightfunction_scalar_laws(priors, probabilities, colnames(samples), source, posterior_pair)
  atoms <- .posterior_atoms_get(samples, allow_partial = TRUE)
  marginals <- if(is.null(atoms$marginals)) declarations$marginals else atoms$marginals
  supports <- .bt_meta_get(samples, "support")
  if(!.posterior_metadata_is_container(supports)) supports <- declarations$supports
  for(column in names(declarations$marginals)){
    condition <- declarations$marginals[[column]]
    if(inherits(condition, "BayesTools_formula_measure_unavailable")){
      samples <- .bt_formula_measure_mark(samples, column, "atoms", conditionMessage(condition),
        cause = condition$reason, diagnostics = condition$diagnostics)
      declarations$marginals[column] <- list(NULL)
      marginals[column] <- list(NULL)
    }
  }
  declared <- !vapply(declarations$marginals, is.null, logical(1))
  marginals[declared] <- declarations$marginals[declared]
  supported <- !vapply(declarations$supports, is.null, logical(1))
  supports[supported] <- declarations$supports[supported]
  if(is.null(joint_unavailable) || any(!vapply(marginals, is.null, logical(1)))){
    atoms <- .posterior_atoms_new(if(is.null(joint_unavailable)) atoms$locations else NULL,
      if(is.null(joint_unavailable) && !is.null(atoms)) atoms$mass else numeric(),
      column_names = colnames(samples), source = if(is.null(atoms)) source else atoms$source,
      component_probabilities = if(is.null(posterior_pair)) probabilities else posterior_pair$probabilities,
      component_log_probabilities = posterior_pair$logs,
      model_probability_declaration = posterior_pair$declaration, marginals = marginals,
      joint_declared = is.null(joint_unavailable), joint_unavailable = joint_unavailable)
    samples <- .posterior_atoms_set(samples, atoms)
  }else samples <- .bt_meta_set(samples, "atoms", NULL)
  .posterior_support_set(samples, supports)
}

# The declared posterior atoms of 'samples' (their only source), or NULL when
# the atom status is undeclared.
.posterior_atoms_get <- function(samples, allow_partial = FALSE){

  check_bool(allow_partial, "allow_partial", allow_NA = FALSE)
  atoms <- .posterior_atoms_from_attribute(.bt_meta_get(samples, "atoms"))
  if(!is.null(atoms) && !atoms$joint_declared && !allow_partial) stop(atoms$joint_unavailable)
  if(!is.null(atoms) && is.matrix(samples) && ncol(samples)==ncol(atoms$locations) && !is.null(colnames(samples))){
    atoms <- .posterior_atoms_rename_columns(atoms,colnames(samples))
  }
  atoms
}

.posterior_atoms_rename_columns <- function(atoms, column_names){

  colnames(atoms$locations) <- column_names
  if(!is.null(atoms$marginals)){
    names(atoms$marginals) <- column_names
    for(i in seq_along(atoms$marginals)){
      if(!is.null(atoms$marginals[[i]])) colnames(atoms$marginals[[i]]$locations) <- column_names[i]
    }
  }
  atoms
}

.posterior_atoms_for_column <- function(atoms, column){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    return(NULL)
  }
  if(is.character(column)){
    column <- match(column, colnames(atoms$locations))
  }
  if(length(column) != 1L || is.na(column) ||
     column < 1L || column > ncol(atoms$locations)){
    return(NULL)
  }
  if(!is.null(atoms$marginals)) return(atoms$marginals[[column]])
  if(!atoms$joint_declared) return(NULL)

  x <- atoms$locations[, column]
  point_masses <- data.frame(x = x, mass = atoms$mass)
  point_masses <- .posterior_atoms_point_mass_table(point_masses)
  if(is.null(point_masses)){
    return(NULL)
  }

  .posterior_atoms_new(
    locations = matrix(point_masses$x, ncol = 1L),
    mass = point_masses$mass,
    column_names = colnames(atoms$locations)[column],
    source = atoms$source,
    declared = TRUE,
    component_probabilities = atoms$component_probabilities,
    component_log_probabilities = atoms$component_log_probabilities,
    model_probability_declaration = atoms$model_probability_declaration
  )
}

.posterior_atoms_transform <- function(atoms, transformation,
                                       transformation_arguments = NULL){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    return(NULL)
  }
  # the scalar transformation maps every coordinate of the atom locations
  locations <- atoms$locations
  locations[] <- .density.prior_transformation_checked_x(
    as.numeric(atoms$locations),
    transformation,
    transformation_arguments
  )
  if(any(!is.finite(locations))){
    stop("The posterior transformation maps an atom to a non-finite location.",
         call. = FALSE)
  }

  .posterior_atoms_new(
    locations = locations,
    mass = atoms$mass,
    column_names = colnames(atoms$locations),
    source = paste0(atoms$source, ":transformed"),
    declared = TRUE,
    component_probabilities = atoms$component_probabilities,
    component_log_probabilities = atoms$component_log_probabilities,
    model_probability_declaration = atoms$model_probability_declaration,
    joint_declared = atoms$joint_declared,
    joint_unavailable = atoms$joint_unavailable,
    marginals = if(!is.null(atoms$marginals)) lapply(atoms$marginals, function(marginal){
      if(is.null(marginal)) NULL else .posterior_atoms_transform(marginal, transformation, transformation_arguments)
    })
  )
}

.posterior_atoms_linear_transform <- function(atoms, design,
                                              column_names = NULL){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    return(NULL)
  }
  if(!is.matrix(design) || !is.numeric(design) || any(!is.finite(design)) || ncol(design) != ncol(atoms$locations)){
    stop("The atom transformation design does not match the joint atom locations.",
         call. = FALSE)
  }
  if(is.null(column_names)) column_names <- rownames(design)
  if(is.null(column_names) && !is.null(atoms$marginals)) column_names <- paste0("value",seq_len(nrow(design)))

  marginals <- if(!is.null(atoms$marginals)){
    stats::setNames(lapply(seq_len(nrow(design)), function(i){
      active <- which(design[i, ] != 0)
      if(length(active) == 0L) return(.posterior_atoms_new(matrix(0, 1L, 1L), 1, column_names = column_names[[i]]))
      if(length(active) != 1L) return(NULL)
      marginal <- atoms$marginals[[active]]
      if(is.null(marginal)) return(NULL)
      .posterior_atoms_new(marginal$locations * design[i, active], marginal$mass,
        column_names = column_names[[i]], source = paste0(marginal$source, ":linear_transform"),
        component_probabilities = marginal$component_probabilities,
        component_log_probabilities = marginal$component_log_probabilities,
        model_probability_declaration = marginal$model_probability_declaration)
    }), column_names)
  }else NULL
  if(!atoms$joint_declared){
    if(nrow(design) == 1L && !is.null(marginals[[1L]])) return(marginals[[1L]])
    if(is.null(marginals) || all(vapply(marginals, is.null, logical(1)))) return(atoms$joint_unavailable)
  }

  .posterior_atoms_new(
    locations = atoms$locations %*% t(design),
    mass = atoms$mass,
    column_names = column_names,
    source = paste0(atoms$source, ":linear_transform"),
    declared = TRUE,
    component_probabilities = atoms$component_probabilities,
    component_log_probabilities = atoms$component_log_probabilities,
    model_probability_declaration = atoms$model_probability_declaration,
    marginals = marginals,
    joint_declared = atoms$joint_declared, joint_unavailable = atoms$joint_unavailable
  )
}

# Posterior atoms of a mixture or spike-and-slab parameter of a single fit
# from the posterior component index of its draws (an index into the prior's
# component list).
.posterior_atoms_from_components <- function(prior, component, n_columns,
                                             column_names = NULL){

  component <- .posterior_atoms_check_components(prior, component)

  components <- as.list(prior)
  probabilities <- vapply(seq_along(components), function(i){
    mean(component == i)
  }, numeric(1))
  if(is.prior.spike_and_slab(prior) &&
     .posterior_atoms_is_zero_point(.get_spike_and_slab_variable(prior))){
    return(.posterior_atoms_new(
      locations = matrix(0, nrow = 1L, ncol = n_columns), mass = 1,
      column_names = column_names, source = "posterior_indicator",
      component_probabilities = probabilities
    ))
  }
  .posterior_atoms_from_priors(
    priors = components,
    probabilities = probabilities,
    n_columns = n_columns,
    column_names = column_names,
    source = "posterior_indicator"
  )
}

.posterior_atoms_from_ordered_total <- function(
    prior, n_columns, column_names = NULL, component = NULL,
    source = "ordered_total_structure"){

  if(!is.prior.ordered(prior) ||
     !.posterior_atoms_ordered_total_has_spike(prior$total)){
    return(NULL)
  }
  prior <- .prior_ordered_default_bound(prior)
  metadata <- .prior_ordered_metadata(prior)

  exclusion_mass <- .posterior_atoms_ordered_exclusion(prior$total, component)
  inclusion_mass <- 1 - exclusion_mass
  locations <- if(exclusion_mass > 0){
    matrix(0, nrow = 1L, ncol = n_columns)
  }else{
    matrix(numeric(), nrow = 0L, ncol = n_columns)
  }
  masses <- if(exclusion_mass > 0) exclusion_mass else numeric()

  .posterior_atoms_new(
    locations = locations,
    mass = masses,
    column_names = column_names,
    source = source,
    declared = TRUE,
    component_probabilities = c(
      excluded = exclusion_mass,
      included = inclusion_mass
    )
  )
}

# The prior of one mixture component of a coefficient: 'component' indexes
# the declared component list (the model of a model-averaged ensemble, the
# component of a mixture or spike-and-slab prior, or for an ordered prior the
# component of its total); 'total_component' is the component of the total of
# an ordered prior selected within a model.
.posterior_atoms_component_prior <- function(prior_entry, component,
                                             model_mixture,
                                             total_component = NA_integer_){

  if(is.prior.ordered(prior_entry)){
    total <- prior_entry$total
    if(.posterior_atoms_is_ordered_zero_total(prior_entry) ||
       (is.prior.spike_and_slab(total) && .bt_component_is_spike(total, component)) ||
       (!is.prior.spike_and_slab(total) &&
        component %in% .posterior_atoms_ordered_total_zero_components(total))){
      # the total's component is its spike (or a point(0) mixture component)
      return(prior("point", list(location = 0)))
    }
    if(is.prior.mixture(total) && !is.na(component)){
      return(.bt_ordered_localize_total(prior_entry,total[[component]]))
    }
    return(prior_entry)
  }

  if(model_mixture){
    if(is.prior(prior_entry)){
      return(prior_entry)
    }
    if(component < 1L || component > length(prior_entry)){
      return(NULL)
    }
    prior <- prior_entry[[component]]
    if(is.null(prior)){
      return(prior("point", list(location = 0)))
    }
    if(!is.na(total_component) && is.prior.ordered(prior)){
      # the model's total component selects its spike at zero or its slab
      return(.posterior_atoms_component_prior(
        prior, total_component, model_mixture = FALSE
      ))
    }
    return(prior)
  }

  if(is.prior.spike_and_slab(prior_entry)){
    if(.bt_component_is_spike(prior_entry, component)){
      return(prior_entry[[component]])
    }
    return(.get_spike_and_slab_variable(prior_entry))
  }
  if(is.prior.mixture(prior_entry)){
    if(component < 1L || component > length(prior_entry)){
      return(NULL)
    }
    return(prior_entry[[component]])
  }

  prior_entry
}

.posterior_atoms_formula_plan <- function(samples, prior_list){

  parameter_names <- intersect(names(prior_list), names(samples))
  if(length(parameter_names) == 0L){
    return(NULL)
  }

  if(inherits(samples, "as_mixed_posteriors")){
    # the component of each draw: of an ordered prior's total, or of a
    # mixture or spike-and-slab prior (one component otherwise)
    indicators <- lapply(parameter_names, function(parameter){
      total_component <- .bt_meta_get(samples[[parameter]], "ordered_total_component")
      if(!is.null(total_component)){
        return(as.integer(total_component))
      }
      as.integer(.bt_draws_component(samples[[parameter]]))
    })
    lengths <- vapply(indicators, length, integer(1))
    if(length(unique(lengths)) != 1L || lengths[1L] == 0L){
      return(NULL)
    }
    indicators <- do.call(cbind, indicators)
    colnames(indicators) <- parameter_names
    keys <- apply(indicators, 1L, paste0, collapse = "\r")
    unique_keys <- unique(keys)
    rows <- match(unique_keys, keys)
    probabilities <- tabulate(
      match(keys, unique_keys),
      nbins = length(unique_keys)
    ) / length(keys)

    return(list(
      components = indicators[rows, , drop = FALSE],
      probabilities = probabilities,
      index = match(keys, unique_keys),
      model_mixture = FALSE
    ))
  }

  atom_metadata <- lapply(parameter_names, function(parameter){
    .posterior_atoms_get(samples[[parameter]])
  })
  probabilities <- lapply(atom_metadata, function(atoms){
    if(is.null(atoms)) NULL else atoms$component_probabilities
  })
  probabilities <- probabilities[!vapply(probabilities, is.null, logical(1))]
  if(length(probabilities) == 0L){
    return(NULL)
  }
  reference <- probabilities[[1L]]
  aligned <- vapply(probabilities, function(probability){
    isTRUE(all.equal(probability, reference, tolerance = 0))
  }, logical(1))
  if(!all(aligned)){
    return(NULL)
  }

  components <- matrix(
    rep(seq_along(reference), length(parameter_names)),
    nrow = length(reference),
    ncol = length(parameter_names)
  )
  colnames(components) <- parameter_names
  .posterior_atoms_split_model_ordered_totals(
    list(
      components = components,
      probabilities = reference,
      model_mixture = TRUE
    ),
    samples,
    parameter_names
  )
}

# Model-averaged ordered factors whose total has a within-model spike at zero
# carry the per-draw total component (mix_posteriors()). Each model component is
# split by the joint total components of its draws, with the posterior model
# probability times the within-model draw frequency (the single-model plan of
# as_mixed_posteriors() uses the draw frequencies of the components).
.posterior_atoms_split_model_ordered_totals <- function(plan, samples, parameter_names){

  indicators <- lapply(parameter_names, function(parameter){
    .bt_meta_get(samples[[parameter]], "ordered_total_component")
  })
  names(indicators) <- parameter_names
  carrying <- parameter_names[!vapply(indicators, is.null, logical(1))]
  if(length(carrying) == 0L){
    return(plan)
  }

  # the indicators must describe the same mixture draws
  model_component <- as.integer(.bt_meta_get(samples[[carrying[1L]]], "component"))
  aligned <- vapply(carrying, function(parameter){
    parameter_models <- as.integer(.bt_meta_get(samples[[parameter]], "component"))
    identical(parameter_models, model_component) &&
      length(indicators[[parameter]]) == length(model_component)
  }, logical(1))
  if(length(model_component) == 0L || !all(aligned)){
    return(NULL)
  }
  indicator_matrix <- do.call(cbind, lapply(indicators[carrying], as.integer))
  colnames(indicator_matrix) <- carrying

  n_components <- nrow(plan$components)
  components <- vector("list", n_components)
  total_indicators <- vector("list", n_components)
  probabilities <- vector("list", n_components)
  for(i in seq_len(n_components)){
    draws <- which(model_component == plan$components[i, 1L])
    if(length(draws) == 0L){
      # a model without mixture draws keeps its unsplit component
      rows <- integer()
      frequencies <- 1
      row_indicators <- matrix(NA_integer_, 1L, length(carrying))
    }else{
      keys <- apply(indicator_matrix[draws, , drop = FALSE], 1L, paste0, collapse = "\r")
      unique_keys <- unique(keys)
      rows <- draws[match(unique_keys, keys)]
      frequencies <- tabulate(
        match(keys, unique_keys),
        nbins = length(unique_keys)
      ) / length(keys)
      row_indicators <- indicator_matrix[rows, , drop = FALSE]
    }
    component_indicators <- matrix(
      NA_integer_, length(frequencies), length(parameter_names),
      dimnames = list(NULL, parameter_names)
    )
    component_indicators[, carrying] <- row_indicators
    components[[i]] <- plan$components[rep(i, length(frequencies)), , drop = FALSE]
    total_indicators[[i]] <- component_indicators
    probabilities[[i]] <- plan$probabilities[i] * frequencies
  }

  plan$components <- do.call(rbind, components)
  plan$total_indicators <- do.call(rbind, total_indicators)
  plan$probabilities <- unlist(probabilities, use.names = FALSE)
  plan
}

.posterior_atoms_plan_total_indicator <- function(plan, row, parameter){

  if(is.null(plan$total_indicators)){
    return(NA_integer_)
  }
  plan$total_indicators[row, parameter]
}

.posterior_atoms_formula <- function(samples, prior_list, weights,
                                     column_name = "value",
                                     source_transforms = NULL){

  weights <- as.matrix(weights)
  log_columns <- names(source_transforms)[source_transforms %in% "log"]
  if(nrow(weights) == 0L || is.null(colnames(weights))){
    return(NULL)
  }
  if(isTRUE(.bt_meta_get(samples, "transform_scaled"))){
    context <- .prior_density_context(
      prior_list, colnames(weights),
      formula_scale = .bt_meta_get(samples, "formula_scale")
    )
    # the log of the unscaled intercept of a log-intercept formula scaling is
    # linear in the fitted coefficients with the log of the fitted intercept
    # (its 'log' source transformation, applied to its atom locations below)
    standardized <- matrix(0, nrow(weights), length(context$column_names),
      dimnames = list(rownames(weights), context$column_names))
    for(i in seq_len(nrow(weights))){
      row <- .prior_density_context_standardized_weights(context, weights[i, ], source_transforms)
      standardized[i, names(row)] <- row
    }
    weights <- standardized
  }
  ordered <- .bt_ordered_formula_projections(samples,weights,source_transforms)
  if(!is.null(ordered)){
    for(projection in ordered) .bt_ordered_require_measure(projection)
    states <- list(atom=as.vector(do.call(rbind,lapply(ordered,`[[`,"atom"))),
      state=as.vector(do.call(rbind,lapply(ordered,`[[`,"state"))))
    model <- rep(attr(ordered,"model",exact=TRUE),each=length(ordered))
    return(.bt_ordered_projection_atoms(states,column_name,model,attr(ordered,"probabilities",exact=TRUE)))
  }
  plan <- .posterior_atoms_formula_plan(samples, prior_list)
  if(is.null(plan)){
    return(NULL)
  }

  active_columns <- colnames(weights)[colSums(abs(weights)) != 0]
  active_owners <- names(prior_list)[vapply(names(prior_list), function(parameter){
    any(active_columns == parameter |
      startsWith(active_columns, paste0(parameter, "[")))
  }, logical(1))]
  owned_columns <- vapply(active_columns, function(column){
    any(column == active_owners | startsWith(column, paste0(active_owners, "[")))
  }, logical(1))
  if(!all(owned_columns)){
    return(NULL)
  }
  missing_owners <- setdiff(active_owners, colnames(plan$components))
  ordinary <- vapply(prior_list[missing_owners], function(prior){
    is.prior(prior) && !is.prior.mixture(prior) &&
      !is.prior.spike_and_slab(prior) && !is.prior.ordered(prior)
  }, logical(1))
  if(!all(ordinary)){
    return(NULL)
  }

  atom_locations <- numeric()
  atom_masses <- numeric()
  parameter_names <- union(colnames(plan$components), missing_owners)
  for(component_i in seq_len(nrow(plan$components))){
    coefficient_locations <- rep(NA_real_, ncol(weights))
    names(coefficient_locations) <- colnames(weights)

    for(parameter in parameter_names){
      parameter_columns <- colnames(weights) == parameter |
        startsWith(colnames(weights), paste0(parameter, "["))
      if(!any(parameter_columns)){
        next
      }
      component_prior <- if(parameter %in% missing_owners){
        # An ordinary declared prior has no unobserved component selection.
        prior_list[[parameter]]
      }else{
        .posterior_atoms_component_prior(
          prior_list[[parameter]],
          plan$components[component_i, parameter],
          model_mixture = plan$model_mixture,
          total_component = .posterior_atoms_plan_total_indicator(plan, component_i, parameter)
        )
      }
      if(is.null(component_prior)){
        next
      }
      location <- .posterior_atoms_point_location(
        component_prior,
        sum(parameter_columns)
      )
      if(parameter %in% missing_owners && is.prior.point(component_prior) && is.null(location)){
        return(NULL)
      }
      if(!is.null(location)){
        coefficient_locations[parameter_columns] <- location
      }
    }
    # log(intercept) formulas enter the linear predictor through log(location)
    logged <- intersect(log_columns, names(coefficient_locations))
    logged <- logged[!is.na(coefficient_locations[logged])]
    if(length(logged) > 0L){
      if(any(coefficient_locations[logged] <= 0)){
        stop("A log(intercept) formula has a non-positive point-prior intercept.",
             call. = FALSE)
      }
      coefficient_locations[logged] <- log(coefficient_locations[logged])
    }

    for(weight_i in seq_len(nrow(weights))){
      active <- weights[weight_i, ] != 0
      if(any(active) && anyNA(coefficient_locations[active])){
        next
      }
      location <- if(any(active)){
        sum(weights[weight_i, active] * coefficient_locations[active])
      }else{
        0
      }
      atom_locations <- c(atom_locations, location)
      atom_masses <- c(
        atom_masses,
        plan$probabilities[component_i] / nrow(weights)
      )
    }
  }

  point_masses <- .posterior_atoms_point_mass_table(
    data.frame(x = atom_locations, mass = atom_masses)
  )
  if(is.null(point_masses)){
    return(NULL)
  }
  .posterior_atoms_new(
    locations = matrix(point_masses$x, ncol = 1L),
    mass = point_masses$mass,
    column_names = column_name,
    source = "formula_structure",
    declared = TRUE
  )
}

.posterior_atoms_joint_linear <- function(prior_list, plan, design,
                                           source_transforms = NULL,
                                           output_transforms = NULL,
                                           samples = NULL){

  active <- colSums(abs(design)) != 0
  design <- design[, active, drop = FALSE]
  source_transforms <- if(is.null(source_transforms)){
    rep("identity", ncol(design))
  }else{
    source_transforms[colnames(design)]
  }
  if(anyNA(source_transforms) || any(!source_transforms %in% c("identity", "log"))){
    stop("Coefficient source-transform metadata are incomplete or unsupported.", call. = FALSE)
  }
  locations <- matrix(numeric(), 0L, ncol(design),
                      dimnames = list(NULL, colnames(design)))
  masses <- numeric()
  for(i in seq_len(nrow(plan$components))){
    component_locations <- stats::setNames(rep(NA_real_, ncol(design)), colnames(design))
    for(parameter in colnames(plan$components)){
      columns <- colnames(design) == parameter |
        startsWith(colnames(design), paste0(parameter, "["))
      if(!any(columns)) next
      component_prior <- .posterior_atoms_component_prior(
        prior_list[[parameter]], plan$components[i, parameter], plan$model_mixture,
        total_component = .posterior_atoms_plan_total_indicator(plan, i, parameter)
      )
      point <- .posterior_atoms_point_location(component_prior, sum(columns))
      if(!is.null(point)) component_locations[columns] <- point
    }
    if(anyNA(component_locations)) next
    logged <- source_transforms == "log"
    component_locations[logged] <- log(component_locations[logged])
    locations <- rbind(locations, component_locations)
    masses <- c(masses, plan$probabilities[i])
  }
  atoms <- .posterior_atoms_new(
    locations, masses, source = "joint_coefficient_structure"
  )
  atoms <- .posterior_atoms_linear_transform(atoms, design, rownames(design))
  if(!is.null(output_transforms)){
    output_transforms <- output_transforms[colnames(atoms$locations)]
    if(anyNA(output_transforms) || any(!output_transforms %in% c("identity", "exp"))){
      stop("Coefficient output-transform metadata are incomplete or unsupported.", call. = FALSE)
    }
    exponentiated <- output_transforms == "exp"
    atoms$locations[, exponentiated] <- .density.prior_transformation_checked_x(
      atoms$locations[, exponentiated, drop = FALSE], "exp"
    )
  }
  if(!is.null(samples)){
    ordered <- .bt_ordered_formula_projections(samples,design,
      source_transforms[source_transforms!="identity"])
    if(!is.null(ordered)){
      atoms$marginals <- stats::setNames(lapply(seq_along(ordered),function(i){
        margin <- .bt_ordered_projection_atoms(ordered[[i]],rownames(design)[[i]])
        if(!is.null(margin) && !is.null(output_transforms) && identical(output_transforms[[i]],"exp")){
          margin <- .posterior_atoms_transform(margin,"exp")
        }
        margin
      }),rownames(design))
    }
  }
  atoms
}

.posterior_atoms_unscale_joint <- function(atoms, scalar_marginals, columns, unavailable = NULL){

  atoms <- .posterior_atoms_rename_columns(.posterior_atoms_from_attribute(atoms), columns)
  joint_marginals <- atoms$marginals
  certified <- character()
  marginals <- scalar_marginals
  if(!is.null(joint_marginals)){
    for(column in columns){
      marginal <- joint_marginals[[column]]
      if(is.null(marginal)) next
      marginal <- .posterior_atoms_from_attribute(marginal)
      marginals[column] <- list(marginal)
      if(isTRUE(marginal$declared)) certified <- c(certified, column)
    }
  }
  atoms$marginals <- marginals
  atoms <- .posterior_atoms_from_attribute(atoms)
  if(!is.null(unavailable)){
    .bt_meta_validate("measure_unavailable", unavailable)
    unavailable <- unavailable[!(unavailable$measure == "atoms" &
      unavailable$column %in% certified), , drop = FALSE]
    if(nrow(unavailable) == 0L) unavailable <- NULL
  }
  list(atoms = atoms, certified_columns = certified, unavailable = unavailable)
}

.posterior_atoms_unscale_mixed <- function(
    samples, model, model_samples, prior_list, formula_scale,
    conditional, conditional_rule, n_grid = .prior_linear_density_default_grid()){

  state <- .bt_formula_state_get(samples)
  if(is.null(state)) return(samples)
  attach_ordered <- function(x, transform, columns){
    if(any(!is.finite(transform$matrix[columns, , drop = FALSE]))) return(x)
    fixed <- intersect(names(prior_list), names(JAGS_formula_design(model, transform$parameter)$prior_list))
    if(!any(vapply(prior_list[fixed], is.prior.ordered, logical(1)))) return(x)
    weights <- .bt_formula_static_rows(transform, columns)
    active <- colnames(weights)[colSums(abs(weights)) != 0]
    contributors <- fixed[vapply(fixed, function(parameter){
      any(active %in% .prior_linear_prior_columns(parameter, prior_list[[parameter]]))
    }, logical(1))]
    requested <- fixed[fixed %in% union(contributors, intersect(names(samples), fixed))]
    raw <- as_mixed_posteriors(model, unique(c(requested, intersect(conditional, names(prior_list)))),
      conditional, conditional_rule, transform_scaled = FALSE, n_prior_samples = n_grid)
    projections <- .bt_ordered_formula_projections(raw, weights,
      transform$source_transforms[transform$source_transforms != "identity"], weight_space = "coefficient")
    .bt_ordered_attach_linear_view(x, projections, weights, raw)
  }
  for(owner in names(samples)){
    prefix <- .bt_meta_get(samples[[owner]], "formula_parameter")
    if(is.null(prefix) || length(prefix) != 1L || !prefix %in% names(formula_scale)) next
    columns <- .posterior_atoms_coefficient_columns(samples[[owner]], owner)
    scale <- state$models[[1L]]$formula_scale[[prefix]]
    spec <- attr(scale, "unscale_design", exact = TRUE)
    transform <- .bt_formula_coefficient_transform(names(spec$multipliers), scale, prefix)
    if(!all(columns %in% transform$target_names)) next
    samples[[owner]] <- .bt_formula_dynamic_quantities(samples[[owner]], transform, columns)
    source <- .bt_meta_get(samples[[owner]], "ordered_source")
    if(!is.null(source) && length(scale) == 0L){
      samples[[owner]] <- .bt_ordered_source_semantics(samples[[owner]],
        diag(length(columns)), colnames(samples[[owner]]))
      next
    }
    cached <- .bt_meta_get(samples, "prior_densities")
    if(!is.matrix(samples[[owner]])){
      samples[[owner]] <- .bt_formula_contribution_metadata(samples[[owner]], state,
        NULL, NULL, prefix, owner, prior_samples = TRUE, n_grid = n_grid, coefficient = columns[[1L]],
        prior_density = cached[[columns[[1L]]]])
      context <- .bt_meta_get(samples, "prior_context")
      weights <- stats::setNames(rep(0, length(context$column_names)), context$column_names)
      weights[[columns[[1L]]]] <- 1
      samples[[owner]] <- .bt_meta_assign(samples[[owner]], list(prior_context = context,
        linear_weights = weights, linear_weight_space = "coefficient"))
      samples[[owner]] <- attach_ordered(samples[[owner]], transform, columns)
      next
    }
    supports <- atom_marginals <- densities <- vector("list", length(columns))
    names(supports) <- names(atom_marginals) <- names(densities) <- colnames(samples[[owner]])
    unavailable <- NULL
    for(i in seq_along(columns)){
      column <- colnames(samples[[owner]])[[i]]
      scalar <- as.numeric(samples[[owner]][, i])
      attr(scalar, "parameter") <- column
      scalar <- .bt_formula_contribution_metadata(scalar, state, NULL, NULL, prefix,
        column, prior_samples = TRUE, n_grid = n_grid, coefficient = columns[[i]],
        prior_density = cached[[columns[[i]]]])
      supports[i] <- list(.bt_meta_get(scalar, "support"))
      atom_marginals[i] <- list(.bt_meta_get(scalar, "atoms"))
      densities[i] <- list(.bt_meta_get(scalar, "prior_density"))
      unavailable <- rbind(unavailable, .bt_meta_get(scalar, "measure_unavailable"))
    }
    samples[[owner]] <- .bt_meta_assign(samples[[owner]], list(support = supports,
      atoms = .posterior_atoms_new(column_names = colnames(samples[[owner]]),
        source = "coefficient_recipe", marginals = atom_marginals),
      components = NULL, prior_densities = densities, measure_unavailable = unavailable))
    if(all(is.finite(transform$matrix[columns, , drop = FALSE]))){
      weights <- .bt_formula_static_rows(transform, columns)
      active <- colnames(weights)[colSums(abs(weights)) != 0]
      contributors <- names(prior_list)[vapply(names(prior_list), function(parameter){
        any(active %in% .prior_linear_prior_columns(parameter, prior_list[[parameter]]))
      }, logical(1))]
      raw <- as_mixed_posteriors(model, unique(c(contributors, intersect(conditional, names(prior_list)))),
        conditional, conditional_rule, transform_scaled = FALSE, n_prior_samples = n_grid)
      plan <- .posterior_atoms_formula_plan(raw, prior_list)
      if(!is.null(plan)){
        atoms <- .posterior_atoms_joint_linear(prior_list, plan, weights,
          source_transforms = transform$source_transforms, output_transforms = transform$output_transforms,
          samples = raw)
        joint <- .posterior_atoms_unscale_joint(atoms, atom_marginals, colnames(samples[[owner]]),
          .bt_meta_get(samples[[owner]], "measure_unavailable"))
        samples[[owner]] <- .posterior_atoms_set(samples[[owner]], joint$atoms)
        samples[[owner]] <- .bt_meta_set(samples[[owner]], "measure_unavailable", joint$unavailable)
      }
    }
    if(!is.null(source)){
      if(any(!is.finite(transform$matrix[columns, , drop = FALSE]))){
        source$view_transformations <- c(source$view_transformations,
          list(list(transformation = "unavailable", arguments = list())))
        samples[[owner]] <- .bt_meta_set(samples[[owner]], "ordered_source", source)
      }else{
        samples[[owner]] <- attach_ordered(samples[[owner]], transform, columns)
      }
    }
  }
  samples
}

# Coefficient (index) names of mixed posterior columns. Factor samples are
# labelled by factor level or contrast coefficient `{j}`, while formula-scale
# transformations name coefficients by index.
.posterior_atoms_coefficient_columns <- function(parameter_samples, parameter){

  if(!is.matrix(parameter_samples)){
    return(parameter)
  }
  if(inherits(parameter_samples, "mixed_posteriors.factor") ||
     isTRUE(attr(parameter_samples, "treatment", exact = TRUE)) ||
     isTRUE(attr(parameter_samples, "independent", exact = TRUE))){
    n_columns <- ncol(parameter_samples)
    if(n_columns == 1L){
      return(parameter)
    }
    return(paste0(parameter, "[", seq_len(n_columns), "]"))
  }

  colnames(parameter_samples)
}
