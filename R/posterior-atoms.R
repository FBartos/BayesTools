#' @title Declare posterior point-mass metadata
#'
#' @description Creates exact posterior measure metadata used to distinguish a
#' continuous posterior from one containing atoms. Supply \code{NULL} to
#' explicitly declare that the posterior has no point masses.
#'
#' @param point_masses optional data frame or list with \code{x} (or
#' \code{location}) and \code{mass} (or \code{p}) entries.
#' @param source short label describing the source of the declaration.
#'
#' @return A posterior atom metadata object suitable for the
#' \code{posterior_atoms} attribute.
#'
#' @export
posterior_atom_attribute <- function(point_masses = NULL, source = "user"){

  check_char(source, "source", check_length = 1L, allow_NA = FALSE)
  point_masses <- .posterior_density_point_masses(point_masses)
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

.posterior_atoms_new <- function(locations = NULL, mass = numeric(),
                                 column_names = NULL,
                                 source = "structural",
                                 declared = TRUE,
                                 component_probabilities = NULL){

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
    component_probabilities <- component_probabilities / probability_sum
  }
  if(!is.null(column_names)){
    if(length(column_names) != ncol(locations)){
      stop("Posterior atom column names do not match the location matrix.",
           call. = FALSE)
    }
    colnames(locations) <- column_names
  }

  out <- list(
    declared = isTRUE(declared),
    locations = locations,
    mass = mass,
    source = source,
    component_probabilities = component_probabilities
  )
  class(out) <- c("BayesTools_posterior_atoms", "list")
  out
}

.posterior_atoms_from_attribute <- function(atoms){

  if(is.null(atoms) || !is.list(atoms) || !isTRUE(atoms$declared)){
    return(NULL)
  }

  tryCatch(
    .posterior_atoms_new(
      locations = atoms$locations,
      mass = atoms$mass,
      column_names = colnames(atoms$locations),
      source = if(is.null(atoms$source)) "unknown" else atoms$source,
      declared = TRUE,
      component_probabilities = atoms$component_probabilities
    ),
    error = function(e) NULL
  )
}

.posterior_atoms_point_location <- function(prior, n_columns){

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

.posterior_atoms_is_ordered_zero_total <- function(prior){

  if(!is.prior.ordered(prior) || !is.prior.point(prior$total)){
    return(FALSE)
  }
  location <- prior$total$parameters[["location"]]
  is.numeric(location) && length(location) == 1L && isTRUE(location == 0)
}

.posterior_atoms_from_priors <- function(priors, probabilities, n_columns = 1L,
                                         column_names = NULL,
                                         source = "model_probabilities",
                                         null_location = NULL,
                                         exclusion_probabilities = NULL){

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

  locations <- matrix(numeric(), nrow = 0L, ncol = n_columns)
  masses <- numeric()
  for(i in seq_along(priors)){
    if(probabilities[i] <= 0){
      next
    }
    location <- .posterior_atoms_point_location(priors[[i]], n_columns)
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
    component_probabilities = probabilities
  )
}

.posterior_atoms_set <- function(samples, atoms){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    stop("Cannot attach invalid posterior atom metadata.", call. = FALSE)
  }
  attr(samples, "posterior_atoms") <- atoms
  samples
}

.posterior_atoms_get <- function(samples){

  atoms <- .posterior_atoms_from_attribute(
    attr(samples, "posterior_atoms", exact = TRUE)
  )
  if(!is.null(atoms)){
    return(atoms)
  }

  density <- .posterior_density_from_attribute(
    attr(samples, "posterior_density", exact = TRUE)
  )
  if(!is.null(density) &&
     .posterior_density_point_masses_declared(density)){
    point_masses <- density$point_masses
    return(.posterior_atoms_new(
      locations = matrix(point_masses$x, ncol = 1L),
      mass = point_masses$mass,
      source = "posterior_density",
      declared = TRUE
    ))
  }

  NULL
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

  x <- atoms$locations[, column]
  point_masses <- data.frame(x = x, mass = atoms$mass)
  point_masses <- .posterior_density_point_masses(point_masses)
  if(is.null(point_masses)){
    return(NULL)
  }

  .posterior_atoms_new(
    locations = matrix(point_masses$x, ncol = 1L),
    mass = point_masses$mass,
    column_names = colnames(atoms$locations)[column],
    source = atoms$source,
    declared = TRUE,
    component_probabilities = atoms$component_probabilities
  )
}

.posterior_atoms_transform <- function(atoms, transformation,
                                       transformation_arguments = NULL){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    return(NULL)
  }
  if(ncol(atoms$locations) != 1L){
    stop("A scalar transformation requires scalar posterior atom metadata.",
         call. = FALSE)
  }
  locations <- .density.prior_transformation_x(
    atoms$locations[, 1L],
    transformation,
    transformation_arguments
  )
  if(any(!is.finite(locations))){
    stop("The posterior transformation maps an atom to a non-finite location.",
         call. = FALSE)
  }

  .posterior_atoms_new(
    locations = matrix(locations, ncol = 1L),
    mass = atoms$mass,
    column_names = colnames(atoms$locations),
    source = paste0(atoms$source, ":transformed"),
    declared = TRUE,
    component_probabilities = atoms$component_probabilities
  )
}

.posterior_atoms_linear_transform <- function(atoms, design,
                                              column_names = NULL){

  atoms <- .posterior_atoms_from_attribute(atoms)
  if(is.null(atoms)){
    return(NULL)
  }
  if(!is.matrix(design) || ncol(design) != ncol(atoms$locations)){
    stop("The atom transformation design does not match the joint atom locations.",
         call. = FALSE)
  }

  .posterior_atoms_new(
    locations = atoms$locations %*% t(design),
    mass = atoms$mass,
    column_names = column_names,
    source = paste0(atoms$source, ":linear_transform"),
    declared = TRUE,
    component_probabilities = atoms$component_probabilities
  )
}

.posterior_atoms_from_indicator <- function(prior, indicator, n_columns,
                                            column_names = NULL,
                                            spike_and_slab = FALSE){

  indicator <- as.integer(indicator)
  if(anyNA(indicator)){
    stop("Posterior component indicators must not be missing.", call. = FALSE)
  }

  if(spike_and_slab){
    exclusion_mass <- mean(indicator == 0L)
    locations <- if(exclusion_mass > 0){
      matrix(0, nrow = 1L, ncol = n_columns)
    }else{
      matrix(numeric(), nrow = 0L, ncol = n_columns)
    }
    masses <- if(exclusion_mass > 0) exclusion_mass else numeric()
    return(.posterior_atoms_new(
      locations = locations,
      mass = masses,
      column_names = column_names,
      source = "posterior_indicator",
      declared = TRUE,
      component_probabilities = c(excluded = exclusion_mass,
                                  included = 1 - exclusion_mass)
    ))
  }

  components <- if(is.prior(prior)) as.list(prior) else prior
  probabilities <- vapply(seq_along(components), function(i){
    mean(indicator == i)
  }, numeric(1))
  .posterior_atoms_from_priors(
    priors = components,
    probabilities = probabilities,
    n_columns = n_columns,
    column_names = column_names,
    source = "posterior_indicator"
  )
}

.posterior_atoms_from_ordered_total <- function(
    prior, n_columns, column_names = NULL, indicator = NULL,
    source = "ordered_total_structure"){

  if(!is.prior.ordered(prior) || !is.prior.spike_and_slab(prior$total)){
    return(NULL)
  }
  prior <- .prior_ordered_default_bound(prior)
  metadata <- .prior_ordered_metadata(prior)
  if(metadata$theta_dim != 1L){
    return(NULL)
  }

  if(is.null(indicator)){
    inclusion_mass <- mean(.get_spike_and_slab_inclusion(prior$total))
  }else{
    indicator <- as.integer(indicator)
    if(length(indicator) == 0L || anyNA(indicator) ||
       any(!indicator %in% c(0L, 1L))){
      stop(
        "Ordered-total posterior indicators must contain only zero and one.",
        call. = FALSE
      )
    }
    inclusion_mass <- mean(indicator == 1L)
  }
  exclusion_mass <- 1 - inclusion_mass
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

.posterior_atoms_component_prior <- function(prior_entry, component,
                                             model_mixture){

  if(is.prior.ordered(prior_entry)){
    if(.posterior_atoms_is_ordered_zero_total(prior_entry)){
      return(prior("point", list(location = 0)))
    }
    if(is.prior.spike_and_slab(prior_entry$total)){
      excluded_component <- if(model_mixture) 1L else 0L
      if(component == excluded_component){
        return(prior("point", list(location = 0)))
      }
    }
    return(NULL)
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
    return(prior)
  }

  if(is.prior.spike_and_slab(prior_entry)){
    if(component == 0L){
      return(prior("point", list(location = 0)))
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
    indicators <- lapply(parameter_names, function(parameter){
      ordered_indicator <- attr(
        samples[[parameter]],
        "ordered_total_indicator",
        exact = TRUE
      )
      if(!is.null(ordered_indicator)){
        return(as.integer(ordered_indicator))
      }
      indicator <- attr(samples[[parameter]], "models_ind", exact = TRUE)
      if(is.null(indicator)){
        return(rep(1L, NROW(samples[[parameter]])))
      }
      as.integer(indicator)
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
  list(
    components = components,
    probabilities = reference,
    model_mixture = TRUE
  )
}

.posterior_atoms_formula <- function(samples, prior_list, weights,
                                     transformation = NULL,
                                     transformation_arguments = NULL,
                                     column_name = "value",
                                     source_transforms = NULL){

  weights <- as.matrix(weights)
  log_columns <- names(source_transforms)[source_transforms %in% "log"]
  if(nrow(weights) == 0L || is.null(colnames(weights))){
    return(NULL)
  }
  if(isTRUE(attr(samples, "transform_scaled", exact = TRUE))){
    context <- .prior_density_context(
      prior_list, colnames(weights),
      formula_scale = attr(samples, "formula_scale", exact = TRUE)
    )
    standardized <- matrix(0, nrow(weights), ncol(weights), dimnames = dimnames(weights))
    for(i in seq_len(nrow(weights))){
      row <- .prior_density_context_standardized_weights(context, weights[i, ])
      standardized[i, names(row)] <- row
    }
    weights <- standardized
  }
  plan <- .posterior_atoms_formula_plan(samples, prior_list)
  if(is.null(plan)){
    return(NULL)
  }

  atom_locations <- numeric()
  atom_masses <- numeric()
  parameter_names <- colnames(plan$components)
  for(component_i in seq_len(nrow(plan$components))){
    coefficient_locations <- rep(NA_real_, ncol(weights))
    names(coefficient_locations) <- colnames(weights)

    for(parameter in parameter_names){
      parameter_columns <- colnames(weights) == parameter |
        startsWith(colnames(weights), paste0(parameter, "["))
      if(!any(parameter_columns)){
        next
      }
      component_prior <- .posterior_atoms_component_prior(
        prior_list[[parameter]],
        plan$components[component_i, parameter],
        model_mixture = plan$model_mixture
      )
      if(is.null(component_prior)){
        next
      }
      location <- .posterior_atoms_point_location(
        component_prior,
        sum(parameter_columns)
      )
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

  point_masses <- .posterior_density_point_masses(
    data.frame(x = atom_locations, mass = atom_masses)
  )
  if(is.null(point_masses)){
    return(NULL)
  }
  atoms <- .posterior_atoms_new(
    locations = matrix(point_masses$x, ncol = 1L),
    mass = point_masses$mass,
    column_names = column_name,
    source = "formula_structure",
    declared = TRUE
  )
  if(!is.null(transformation)){
    atoms <- .posterior_atoms_transform(
      atoms,
      transformation,
      transformation_arguments
    )
  }
  atoms
}

.posterior_atoms_joint_linear <- function(prior_list, plan, design,
                                           source_transforms = NULL,
                                           output_transforms = NULL){

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
        prior_list[[parameter]], plan$components[i, parameter], plan$model_mixture
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
    atoms$locations[, exponentiated] <- exp(atoms$locations[, exponentiated, drop = FALSE])
  }
  atoms
}

.posterior_atoms_unscale_mixed <- function(
    samples, model, model_samples, prior_list, formula_scale,
    conditional, conditional_rule){

  for(prefix in names(formula_scale)){
    columns <- colnames(model_samples)[
      .formula_scale_matches_prefix(colnames(model_samples), prefix)
    ]
    columns <- columns[
      !.formula_scale_matches_prefix(columns, prefix, "__xREx__") &
      !.formula_scale_matches_prefix(columns, prefix, "__xRE_ALLOCx") &
      !.formula_scale_matches_prefix(columns, prefix, "__xRE_SUMMARY__")
    ]
    if(length(columns) == 0L) next
    transform <- .bt_formula_coefficient_transform(
      source_names = columns, formula_scale = formula_scale[[prefix]],
      parameter = prefix
    )
    requested_columns <- unique(unlist(lapply(names(samples), function(parameter){
      if(is.matrix(samples[[parameter]])) colnames(samples[[parameter]]) else parameter
    }), use.names = FALSE))
    requested_columns <- intersect(requested_columns, transform$target_names)
    if(length(requested_columns) == 0L) next
    requested_design <- transform$matrix[requested_columns, , drop = FALSE]
    active_columns <- colnames(requested_design)[colSums(abs(requested_design)) != 0]
    contributors <- names(prior_list)[vapply(names(prior_list), function(name){
      any(active_columns == name | startsWith(active_columns, paste0(name, "[")))
    }, logical(1))]
    continuous <- vapply(prior_list[contributors], function(prior){
      is.prior.simple(prior) && !is.prior.point(prior) &&
        !is.prior.discrete(prior) && !is.prior.spike_and_slab(prior) &&
        !is.prior.mixture(prior) && is.null(attr(prior, "multiply_by", exact = TRUE))
    }, logical(1))
    if(length(contributors) > 0L && all(continuous) &&
       all(rowSums(abs(requested_design)) > 0)){
      # Nonconstant combinations of independent continuous coefficients remain
      # atom-free, so their existing empty declarations require no draw replay.
      next
    }
    raw_samples <- as_mixed_posteriors(
      model, parameters = unique(c(contributors, intersect(conditional, names(prior_list)))),
      conditional = conditional, conditional_rule = conditional_rule,
      transform_scaled = FALSE
    )
    plan <- .posterior_atoms_formula_plan(raw_samples, prior_list)
    if(is.null(plan)){
      stop("Joint posterior atom metadata are unavailable for unscaled coefficients of '",
           prefix, "'.", call. = FALSE)
    }
    for(parameter in names(samples)){
      target_columns <- if(is.matrix(samples[[parameter]])){
        colnames(samples[[parameter]])
      }else{
        parameter
      }
      if(!all(target_columns %in% transform$target_names)) next
      design <- transform$matrix[target_columns, , drop = FALSE]
      samples[[parameter]] <- .posterior_atoms_set(
        samples[[parameter]],
        .posterior_atoms_joint_linear(
          prior_list, plan, design,
          source_transforms = transform$source_transforms,
          output_transforms = transform$output_transforms
        )
      )
    }
  }
  samples
}
