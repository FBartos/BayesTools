# A named regularization group supplies the coefficient priors of its terms
# before formula defaults and contrasts are resolved. The ordinary formula
# compiler determines their fitted basis and dimension.
.bt_regularized_assign_terms <- function(prior_list, model_terms, model_terms_type){

  labels <- names(prior_list)
  declarations <- vapply(prior_list, inherits, logical(1), "prior.regularized")
  # Fixed-scale ridge/lasso are ordinary scalar priors when assigned to a
  # parameter. Under a group name the same constructor can supply all slopes.
  if(!is.null(labels)){
    for(term in intersect(labels, model_terms[model_terms_type == "factor"])){
      p <- prior_list[[term]]
      declaration <- attr(p, "regularized_spec", exact = TRUE)
      if(!is.null(declaration$contrast) && !is.prior.factor(p)){
        prior_list[[term]] <- prior_factor(p$distribution, p$parameters, truncation = p$truncation,
          prior_weights = .prior_model_weight(p), contrast = declaration$contrast)
      }
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
    if(inherits(group, "prior.regularized")) group else attr(group, "regularized_spec", exact = TRUE)
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
  if(length(implicit)) members[[implicit]] <- setdiff(terms, c(names(prior_list), explicit))

  for(group in names(groups)){
    .bt_check_jags_node_name(group, "regularization group name")
    if(!length(members[[group]])){
      stop("Regularization group '", group, "' has no remaining slopes.", call. = FALSE)
    }
    spec <- specs[[group]]
    multiplier <- attr(groups[[group]], "regularized_scale_multiplier", exact = TRUE)
    if(is.null(multiplier)) multiplier <- 1
    .regularized_check_positive(multiplier, "regularized_scale_multiplier")
    contrast <- spec$contrast
    if(is.null(contrast)) contrast <- attr(groups[[group]], "regularized_contrast", exact = TRUE)
    if(is.null(contrast)) contrast <- "treatment"
    fixed <- is.numeric(spec$scale) && spec$family %in% c("ridge", "lasso")
    for(term in members[[group]]){
      factor_term <- identical(model_terms_type[[term]], "factor")
      if(factor_term && !contrast %in% c("treatment", "independent") &&
         !(fixed && identical(spec$family, "ridge") && contrast %in% c("meandif", "orthonormal"))){
        stop("Regularized factor term '", term, "' requires treatment or independent contrasts.",
             call. = FALSE)
      }
      if(fixed){
        law <- if(identical(spec$family, "ridge")) "normal" else "laplace"
        if(factor_term && contrast %in% c("meandif", "orthonormal")) law <- "mnormal"
        par <- list(0, spec$scale * multiplier)
        p <- if(factor_term) prior_factor(law, par, contrast = contrast) else prior(law, par)
      }else{
        p <- list(
          distribution = "regularized", parameters = list(),
          truncation = list(lower = -Inf, upper = Inf), prior_weights = 1,
          regularization = list(
            group = group, term = term, family = spec$family, scale = spec$scale,
            slab_scale = spec$slab_scale, slab_df = spec$slab_df, multiplier = multiplier
          )
        )
        class(p) <- c("prior", "prior.regularized_term", "prior.simple",
                      if(factor_term) c("prior.factor", paste0("prior.", contrast)))
      }
      prior_list[[term]] <- p
    }
  }
  prior_list
}

# Rebind from the fitted design, including groups inside ordinary mixtures.
# This also elects a single owner for the shared nodes when only a subset of
# a previously compiled group's members is present in a submodel.
.bt_regularized_bind_terms <- function(prior_list, parameter){

  groups <- list()
  bound <- function(prior, term){
    if(is.prior.mixture(prior)){
      for(i in seq_along(prior)) prior[[i]] <- bound(prior[[i]], term)
      return(prior)
    }
    if(!inherits(prior, "prior.regularized_term")) return(prior)
    spec <- prior$regularization
    if(is.prior.factor(prior) && any(!attr(prior, "factor_contrasts", exact = TRUE) %in%
                                     c("contr.treatment", "contr.independent"))){
      stop("Regularized term '", spec$term,
           "' uses a fitted factor basis other than treatment or independent contrasts.", call. = FALSE)
    }
    if(!identical(gsub(":", "__xXx__", spec$term, fixed = TRUE), term)){
      stop("The regularization declaration for '", spec$term, "' cannot be bound to '", term, "'.", call. = FALSE)
    }
    group <- spec$group
    record <- groups[[group]]
    definition <- spec[c("family", "scale", "slab_scale", "slab_df", "multiplier")]
    if(!is.null(record) && (!identical(definition, record$definition) || term %in% record$terms)){
      stop("Incompatible or duplicate regularized members for group '", group, "'.", call. = FALSE)
    }
    groups[[group]] <<- list(definition = definition, terms = c(record$terms, term))
    stem <- paste(parameter, group, sep = "_")
    d <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else 1L
    spec$dimension <- d
    spec$owner <- is.null(record)
    spec$global_name <- if(!is.numeric(spec$scale)) paste0(stem, "_global")
    spec$slab_name <- if(!is.null(spec$slab_scale) && is.finite(spec$slab_scale)) paste0(stem, "_slab2")
    spec$local_names <- if(identical(spec$family, "horseshoe")){
      paste0(parameter, "_", term, "_local_", seq_len(d))
    }else character()
    prior$regularization <- spec
    prior
  }
  for(term in names(prior_list)) prior_list[[term]] <- bound(prior_list[[term]], term)
  prior_list
}

.bt_regularized_primitives <- function(prior, shared = TRUE){

  spec <- prior$regularization
  primitives <- list()
  if(shared && spec$owner){
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
  primitives
}

# At a public compilation boundary the whole set of prior owners is available.
# Compare actual emitted roots, including ordinary mixture and factor roots,
# so that a user's parameter cannot silently become a group's hyperparameter.
.bt_regularized_validate_roots <- function(prior_list){

  members <- function(p){
    if(inherits(p, "prior.regularized_term")) return(list(p))
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

  has_member <- function(p){
    inherits(p, "prior.regularized_term") ||
      (is.prior.mixture(p) && any(vapply(p, has_member, logical(1))))
  }
  for(name in names(design$prior_list)){
    compiled <- design$prior_list[[name]]
    fitted <- fit_priors[[name]]
    if((has_member(compiled) || has_member(fitted)) &&
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
  if(length(spec$local_names)) log_scale <- paste0(log_scale, " + log(", spec$local_names, ")")
  if(identical(spec$family, "lasso")){
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
