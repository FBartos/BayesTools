# A named regularization group supplies the coefficient priors of its terms
# before formula defaults and contrasts are resolved. The ordinary formula
# compiler determines their fitted basis and dimension.
.bt_regularized_assign_terms <- function(prior_list, model_terms, model_terms_type,
                                         R2_budgets = list()){

  labels <- names(prior_list)
  declarations <- vapply(prior_list, is.prior.regularized, logical(1))
  # Fixed-scale ridge/lasso are ordinary scalar priors when assigned to a
  # parameter. Under a group name the same constructor can supply all slopes.
  if(!is.null(labels)){
    for(term in intersect(labels, model_terms[model_terms_type == "factor"])){
      p <- prior_list[[term]]
      declaration <- attr(p, "regularized_spec", exact = TRUE)
      if(is.null(declaration) || is.prior.factor(p)) next
      prior_list[[term]] <- .bt_regularized_term_prior(declaration, term, term, TRUE,
        attr(p, "regularized_contrast", exact = TRUE), multiplier = 1, prior_weights = .prior_model_weight(p))
    }
    scalar_groups <- !labels %in% c(model_terms, .JAGS_formula_default_prior_names()) &
      vapply(prior_list, function(p) !is.null(attr(p, "regularized_spec", exact = TRUE)), logical(1))
    declarations <- declarations | scalar_groups
  }
  if(!any(declarations)) return(prior_list)
  if(is.null(labels) || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)){
    stop("Regularization groups must be uniquely named entries in the formula prior list.", call. = FALSE)
  }
  groups <- prior_list[declarations]
  prior_list <- prior_list[!declarations]
  terms <- setdiff(model_terms, "intercept")
  keys <- function(x) vapply(strsplit(x, ":", fixed = TRUE), function(parts){
    paste(sort(parts), collapse = ":")
  }, character(1), USE.NAMES = FALSE)
  if(any(names(groups) %in% c(model_terms, .JAGS_formula_default_prior_names()))){
    stop("A regularization group needs a name distinct from the formula terms and defaults.", call. = FALSE)
  }
  specs <- lapply(groups, function(group){
    if(is.prior.regularized(group)) group else attr(group, "regularized_spec", exact = TRUE)
  })
  members <- lapply(specs, function(spec){
    if(is.null(spec$terms)) return(NULL)
    indices <- match(keys(spec$terms), keys(terms))
    if(anyNA(indices)){
      stop("Unknown regularization term '", spec$terms[which(is.na(indices))[1L]],
           "'. Supply a non-intercept term of this formula.", call. = FALSE)
    }
    terms[indices]
  })
  explicit <- unlist(members, use.names = FALSE)
  if(anyDuplicated(explicit) || any(explicit %in% names(prior_list))){
    stop("A regularized term cannot have its own prior or belong to multiple groups.", call. = FALSE)
  }
  implicit <- which(vapply(members, is.null, logical(1)))
  if(length(implicit) > 1L){
    stop("At most one regularization group per formula may omit 'terms'.", call. = FALSE)
  }
  # Reserve explicitly budgeted R2 terms before expanding a regularization
  # default. The R2 binder runs next and sees the resolved regularized terms.
  reserved <- character()
  if(length(R2_budgets)){
    if(any(names(groups) %in% names(R2_budgets))){
      .bt_R2_stop("An R-squared budget and a regularization group must have distinct names.")
    }
    R2_terms <- lapply(R2_budgets, `[[`, "terms")
    if(length(implicit) && any(vapply(R2_terms, is.null, logical(1)))){
      .bt_R2_stop("An R-squared budget and a regularization group cannot both omit 'terms'.")
    }
    reserved <- terms[keys(terms) %in% keys(as.character(unlist(R2_terms, use.names = FALSE)))]
    overlap <- intersect(explicit, reserved)
    if(length(overlap)){
      .bt_R2_stop(paste0("Term '", overlap[1L],
        "' belongs to both an R-squared budget and a regularization group. Name each term in only one group."))
    }
  }
  if(length(implicit)) members[[implicit]] <- setdiff(terms, c(names(prior_list), explicit, reserved))

  for(group in names(groups)){
    .bt_check_jags_node_name(group, "regularization group name")
    if(!length(members[[group]])){
      stop("Regularization group '", group, "' has no remaining slopes.", call. = FALSE)
    }
    spec <- specs[[group]]
    multiplier <- attr(groups[[group]], "regularized_scale_multiplier", exact = TRUE)
    if(is.null(multiplier)) multiplier <- 1
    .regularized_check_positive(multiplier, "regularized_scale_multiplier")
    fallback_contrast <- attr(groups[[group]], "regularized_contrast", exact = TRUE)
    for(term in members[[group]]){
      factor_term <- identical(model_terms_type[[term]], "factor")
      prior_list[[term]] <- .bt_regularized_term_prior(spec, group, term, factor_term, fallback_contrast, multiplier,
        owned = intersect(strsplit(term, ":", fixed = TRUE)[[1L]], members[[group]]))
    }
  }
  prior_list
}

.bt_regularized_term_prior <- function(spec, group, term, factor_term, fallback_contrast, multiplier, prior_weights = 1,
                                       owned = character()){

  contrast <- spec$contrast
  if(is.null(contrast)) contrast <- fallback_contrast
  if(is.null(contrast)) contrast <- "treatment"
  centered <- factor_term && contrast %in% c("meandif", "orthonormal")
  binding <- if(factor_term && grepl(":", term, fixed = TRUE)){
    list(group = group, requested = spec$contrast, owned = owned)
  }
  # In an interaction the main-effect priors can supply a different basis
  # for each factor; only the fitted design can identify its lasso units.
  needs_lasso_units <- identical(spec$family, "lasso") && factor_term &&
    (centered || grepl(":", term, fixed = TRUE))
  if(is.numeric(spec$scale) && spec$family %in% c("ridge", "lasso") && !needs_lasso_units){
    law <- if(spec$family == "lasso") "laplace" else if(centered) "mnormal" else "normal"
    parameters <- list(0, spec$scale * multiplier)
    if(!factor_term) return(prior(law, parameters, prior_weights = prior_weights))
    p <- prior_factor(law, parameters, contrast = contrast, prior_weights = prior_weights)
    attr(p, "regularized_factor_binding") <- binding
    return(p)
  }
  p <- list(
    distribution = "regularized", parameters = list(),
    truncation = list(lower = -Inf, upper = Inf), prior_weights = prior_weights,
    regularization = list(
      group = group, term = term, family = spec$family, scale = spec$scale,
      slab_scale = spec$slab_scale, slab_df = spec$slab_df, multiplier = multiplier
    )
  )
  class(p) <- c("prior", "prior.regularized_term", "prior.simple",
                if(factor_term) c("prior.factor", paste0("prior.", contrast)))
  attr(p, "regularized_factor_binding") <- binding
  p
}

# A term has one factor class for reporting. Its stored design and the
# separate coding of each constituent define the transformation to cells.
# Choose a centered class whenever a mixed design requires that transformation.
.bt_regularized_factor_parts <- function(prior){
  if(is.prior.mixture(prior)) prior else list(prior)
}

.bt_regularized_interaction <- function(prior){
  any(vapply(.bt_regularized_factor_parts(prior), function(p){
    !is.null(attr(p, "regularized_factor_binding", exact = TRUE))
  }, logical(1)))
}

.bt_regularized_resolve_factor <- function(prior, term, contrasts, indicators){

  if(!.bt_regularized_interaction(prior)) return(prior)
  actual <- sub("^contr\\.", "", contrasts)
  labels <- paste0("'", names(actual), "' = '", actual, "'", collapse = ", ")
  conflict <- function(requested, reason, factors){
    .bt_formula_transform_stop(paste0("Regularized interaction '", term, "' ", reason,
      ". Its fitted factor codings are ", labels, ". Set the main-effect contrasts to match the request, ",
      "or omit 'contrast' to inherit their fitted codings."),
      reason = "regularized_contrast_conflict", terms = term, predictors = factors,
      requested = requested, factor_contrasts = contrasts)
  }
  parts <- .bt_regularized_factor_parts(prior)
  for(p in parts){
    binding <- attr(p, "regularized_factor_binding", exact = TRUE)
    requested <- binding$requested
    if(is.null(requested)) next
    if(length(indicators) == length(actual) && requested %in% c("meandif", "orthonormal")){
      conflict(requested, paste0("codes ", paste0("'", indicators, "'", collapse = " and "),
        " by level indicators, so it cannot use '", requested, "' contrasts"), indicators)
    }
    # A group can supply a main effect and its interactions. A partner
    # with a separately specified main effect keeps its own coding.
    owned <- intersect(names(actual), binding$owned)
    constrained <- setdiff(if(length(owned)) owned else names(actual), indicators)
    different <- constrained[actual[constrained] != requested]
    if(length(different)) conflict(requested, paste0("requested '", requested, "' contrasts"), different)
  }
  representative <- intersect(c("meandif", "orthonormal", "treatment", "independent"), actual)[1L]
  if(is.na(representative)) conflict(NULL, "has no supported factor coding", names(actual))
  reclass <- function(p){
    old <- .factor_object_contrast_name(p)
    if(identical(old, paste0("contr.", representative))) return(p)
    if(is.prior.mixture(p)){
      for(i in seq_along(p)) p[[i]] <- reclass(p[[i]])
    }else if(.posterior_atoms_is_zero_point(p) ||
             (!is.null(attr(p, "regularized_factor_binding", exact = TRUE)) &&
              p$distribution %in% c("normal", "mnormal"))){
      # JAGS uses different point/Normal prior objects for centered and
      # directly coded factors. Recreate the law on the same raw coordinates.
      point <- is.prior.point(p)
      family <- if(point) "point" else if(representative %in% c("meandif", "orthonormal")) "mnormal" else "normal"
      parameters <- if(point) list(0) else list(0, p$parameters$sd)
      replacement <- prior_factor(family, parameters, contrast = representative,
                                  prior_weights = .prior_model_weight(p))
      metadata <- setdiff(names(attributes(p)), c("class", "names"))
      for(name in metadata) attr(replacement, name) <- attr(p, name, exact = TRUE)
      return(replacement)
    }else if(!is.prior.regularized_term(p)){
      conflict(sub("^contr\\.", "", old), "has an incompatible mixture component", names(actual))
    }
    class(p)[class(p) == sub("^contr\\.", "prior.", old)] <- paste0("prior.", representative)
    p
  }
  reclass(prior)
}

# Rebind from the fitted design, including groups inside ordinary mixtures.
# This also elects a single owner for the shared nodes when only a subset of
# a previously compiled group's members is present in a submodel.
.bt_regularized_bind_terms <- function(prior_list, parameter){

  groups <- list()
  # 'node' names the emitted coefficient: a mixture writes each component as
  # its own node (.JAGS_prior.mixture()), so the local scales of components
  # from different groups stay distinct.
  bound <- function(prior, term, node = term){
    if(is.prior.mixture(prior)){
      for(i in seq_along(prior)){
        prior[[i]] <- bound(prior[[i]], term,
                            if(is.prior.spike_and_slab(prior)) node else paste0(node, "_component_", i))
      }
      return(prior)
    }
    if(!is.prior.regularized_term(prior)) return(prior)
    spec <- prior$regularization
    if(!identical(gsub(":", "__xXx__", spec$term, fixed = TRUE), term)){
      stop("The regularization declaration for '", spec$term, "' cannot be bound to '", term, "'.", call. = FALSE)
    }
    d <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else 1L
    spec$dimension <- d
    units <- NULL
    if(is.prior.factor(prior)){
      design <- .factor_term_design_from_metadata(prior)
      design$factor_contrasts <- attr(prior, "factor_contrasts", exact = TRUE)
      units <- .bt_factor_exchangeable_units(term, seq_len(d), design)
    }
    if(is.null(units) && is.numeric(spec$scale) && identical(spec$family, "lasso") &&
       (is.prior.treatment(prior) || is.prior.independent(prior))){
      # A fixed lasso interaction without any centered constituent is the
      # ordinary product Laplace factor, with no extra scale coordinates.
      prior$distribution <- "laplace"
      prior$parameters <- list(location = 0, scale = spec$scale * spec$multiplier)
      prior$regularization <- NULL
      class(prior) <- setdiff(class(prior), "prior.regularized_term")
      return(prior)
    }
    group <- spec$group
    record <- groups[[group]]
    definition <- spec[c("family", "scale", "slab_scale", "slab_df", "multiplier")]
    if(!is.null(record) && (!identical(definition, record$definition) || term %in% record$terms)){
      stop("Incompatible or duplicate regularized members for group '", group, "'.", call. = FALSE)
    }
    groups[[group]] <<- list(definition = definition, terms = c(record$terms, term))
    stem <- paste(parameter, group, sep = "_")
    spec$owner <- is.null(record)
    spec$global_name <- if(!is.numeric(spec$scale)) paste0(stem, "_global")
    spec$slab_name <- if(!is.null(spec$slab_scale) && is.finite(spec$slab_scale)) paste0(stem, "_slab2")
    spec$local_names <- if(identical(spec$family, "horseshoe")){
      paste0(parameter, "_", node, "_local_", seq_len(if(is.null(units)) d else length(units)))
    }else character()
    # Ordinary coordinates have their own local scales (the original fitted
    # declaration). Store the unit map only for actual centered blocks.
    spec[c("unit_index", "unit_dimension", "unit_centered", "variance_names")] <- NULL
    if(!is.null(units)){
      spec$unit_dimension <- vapply(units, function(unit) length(unit$columns), integer(1))
      spec$unit_index <- integer(d)
      for(i in seq_along(units)) spec$unit_index[units[[i]]$columns] <- i
      spec$unit_centered <- rep(TRUE, length(units))
      if(identical(spec$family, "lasso")) spec$variance_names <- paste0(parameter, "_", node, "_variance_", seq_along(units))
    }
    prior$regularization <- spec
    prior
  }
  for(term in names(prior_list)) prior_list[[term]] <- bound(prior_list[[term]], term)
  prior_list
}

.bt_regularized_primitives <- function(prior, include_shared = prior$regularization$owner){

  spec <- prior$regularization
  primitives <- list()
  if(include_shared){
    if(!is.null(spec$global_name)) primitives[[spec$global_name]] <- spec$scale
    if(!is.null(spec$slab_name)){
      shape <- spec$slab_df / 2
      scale <- shape * spec$slab_scale^2
      if(!is.finite(scale) || scale == 0){
        scale <- exp(log(shape) + 2 * log(spec$slab_scale))
      }
      if(!is.finite(scale) || scale == 0 || shape == 0){
        .prior_numerical_signal("horseshoe slab prior", "invgamma", "finite", 1L,
          "The inverse-gamma parameters of the slab variance are not representable in double precision.",
          error = TRUE)
      }
      primitives[[spec$slab_name]] <- prior("invgamma", list(shape, scale))
    }
  }
  if(length(spec$local_names)){
    local <- prior("cauchy", list(0, 1), truncation = list(lower = 0))
    primitives[spec$local_names] <- rep(list(local), length(spec$local_names))
  }
  if(length(spec$variance_names)){
    primitives[spec$variance_names] <- lapply(spec$unit_dimension, function(d){
      prior("gamma", list(shape = (d + 1) / 2, rate = 0.5))
    })
  }
  primitives
}

# At a public compilation boundary the whole set of prior owners is available.
# Compare actual emitted roots, including ordinary mixture and factor roots,
# so that a user's parameter cannot silently become a group's hyperparameter.
.bt_regularized_validate_roots <- function(prior_list){

  members <- function(p){
    if(is.prior.regularized_term(p)) return(list(p))
    if(is.prior.mixture(p)) return(unlist(lapply(p, members), recursive = FALSE, use.names = FALSE))
    list()
  }
  all_members <- unlist(lapply(prior_list, members), recursive = FALSE, use.names = FALSE)
  if(!length(all_members)) return(invisible(TRUE))
  primitives <- unlist(lapply(all_members, function(p){
    names(.bt_regularized_primitives(p))
  }), use.names = FALSE)
  references <- unlist(lapply(all_members, function(p){
    c(p$regularization$global_name, p$regularization$slab_name)
  }), use.names = FALSE)
  missing <- setdiff(references, primitives)
  if(length(missing)){
    stop("Regularization scale '", missing[[1L]], "' has no owner in the supplied prior list.", call. = FALSE)
  }
  all_roots <- unlist(lapply(names(prior_list), function(name){
    names(.bt_prior_emitted_roots(prior_list[[name]], name))
  }), use.names = FALSE)
  collisions <- intersect(primitives, c(primitives[duplicated(primitives)], names(prior_list),
                                        all_roots[duplicated(all_roots)]))
  if(length(collisions)){
    stop("Regularization node '", collisions[[1L]], "' conflicts with another JAGS parameter. Rename the group or term.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# The finite slab adds its precision to the ordinary horseshoe precision.
# Work with log scales so a very wide local scale can still yield a finite
# regularized scale, including when global * local overflows.
.bt_regularized_scale <- function(prior, global, local = 1, slab = NULL){

  log_sd <- log(global) + log(local)
  if(!is.null(slab)){
    log_c <- 0.5 * log(slab)
    log_sd <- pmin(log_sd, log_c) - 0.5 * log1p(exp(-2 * abs(log_sd - log_c)))
  }
  exp(log(prior$regularization$multiplier) + log_sd)
}

# The extension continues the original JAGS graph. Bridge evaluation uses
# the compiled formula's copy of the priors. Both must keep the same joint
# law as the fit-level copy, including the ownership of shared scales.
.bt_regularized_check_fitted_priors <- function(design, fit_priors){

  for(name in names(design$prior_list)){
    compiled <- design$prior_list[[name]]
    fitted <- fit_priors[[name]]
    if((.bt_has_regularized_term(compiled) || .bt_has_regularized_term(fitted)) &&
       !identical(.bt_JAGS_bridge_prior_semantics(compiled), .bt_JAGS_bridge_prior_semantics(fitted))){
      .bt_stop_refit_required("The fitted regularization priors disagree for '", name,
                             "'. Refit the model with the intended group and scale priors.")
    }
  }
  invisible(TRUE)
}

.bt_regularized_jags <- function(prior, name){

  spec <- prior$regularization
  scalar <- .prior_ordered_format_number
  global <- if(is.null(spec$global_name)) scalar(spec$scale) else spec$global_name
  coefficient <- if(is.prior.factor(prior)) paste0(name, "[", seq_len(spec$dimension), "]") else name
  log_scale <- paste0("log(", global, ")")
  if(spec$multiplier != 1) log_scale <- paste0("log(", scalar(spec$multiplier), ") + ", log_scale)
  unit_index <- spec$unit_index
  if(is.null(unit_index)) unit_index <- seq_len(spec$dimension)
  if(length(spec$local_names)) log_scale <- paste0(log_scale, " + log(", spec$local_names[unit_index], ")")
  if(length(spec$variance_names)) log_scale <- paste0(log_scale, " + 0.5 * log(", spec$variance_names[unit_index], ")")
  if(identical(spec$family, "lasso") && !length(spec$variance_names)){
    law <- paste0("ddexp(0, exp(-(", log_scale, ")))")
  }else{
    precision <- paste0("exp(-2 * (", log_scale, "))")
    if(!is.null(spec$slab_name)){
      precision <- paste0(precision, " + exp(-log(", spec$slab_name, ") - 2 * log(", scalar(spec$multiplier), "))")
    }
    law <- paste0("dnorm(0, ", precision, ")")
  }
  primitives <- .bt_regularized_primitives(prior)
  paste0(paste(c(.JAGS_add_priors.fun(primitives, numeric_literals = TRUE),
                 paste0(coefficient, " ~ ", law)), collapse = "\n"), "\n")
}
