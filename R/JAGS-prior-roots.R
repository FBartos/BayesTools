# Declared root shapes mirror the prior emitters; coordinates of one array
# occupy one root. No generated syntax or posterior draws are inspected.
.bt_prior_emitted_roots <- function(prior, parameter){

  roots <- list()
  add <- function(name, shape = integer()){
    roots[[name]] <<- as.integer(shape)
  }
  append <- function(p, name){
    nested <- .bt_prior_emitted_roots(p, name)
    roots[names(nested)] <<- nested
  }
  if(is.null(prior) || is.prior.none(prior)) return(roots)
  if(is.prior.weightfunction(prior) || is_prior_phacking(prior) || is_prior_bias(prior) ||
     inherits(prior, "prior.bias_mixture")){
    branches <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
    has_selection <- vapply(branches, function(x) !is.null(x$selection), logical(1))
    has_phacking <- vapply(branches, function(x) !is.null(x$phacking), logical(1))
    if(any(has_selection) || any(has_phacking)){
      spec <- .bt_dnode_omega_spec(prior)
      multiple <- spec$uses_indicator
      add("sel_vector_rule")
      if(multiple) add("bias_indicator")
      for(k in seq_along(branches)){
        component <- if(multiple) k else NULL
        selection <- branches[[k]]$selection
        if(!is.null(selection)){
          nodes <- .weightfunction_component_node_names(selection, component,
            spec$global_cuts, force_one_sided = TRUE)
          add(nodes$omega_local, nodes$n_bins)
          if(identical(selection$weights$type, "cumulative")){
            if(nodes$n_bins == 2L) add(nodes$omega_ratio) else{
              add(nodes$eta, nodes$n_bins)
              add(nodes$std_eta, nodes$n_bins)
            }
          }else if(identical(selection$weights$type, "independent") &&
                   identical(selection$weights$scale, "log_omega") && nodes$n_bins > 1L){
            add(nodes$log_omega, nodes$n_bins)
          }
          if(multiple || nodes$needs_mapping) add(nodes$omega_target, spec$n_bins)
        }else if(multiple || !is.null(branches[[k]]$phacking)){
          add(if(multiple) paste0("omega_component_", k) else "omega", spec$n_bins)
        }
        if(multiple || !is.null(branches[[k]]$phacking)){
          for(field in c("alpha", "phack_kind", "pi_null", "beta_null")){
            add(.selection_phacking_node_name(field, component))
          }
          for(field in c("phack_z_source", "phack_z_dest")){
            add(.selection_phacking_node_name(field, component), 2L)
          }
        }
      }
      if(multiple){
        add("omega", spec$n_bins)
        for(field in c("alpha", "phack_kind", "pi_null", "beta_null")) add(field)
        if(any(has_phacking)){
          add("phack_z_source", 2L)
          add("phack_z_dest", 2L)
        }
      }
    }else if(inherits(prior, "prior.bias_mixture")) add("bias_indicator")
    if(inherits(prior, "prior.bias_mixture")){
      for(field in c("PET", "PEESE")){
        predicate <- if(field == "PET") is.prior.PET else is.prior.PEESE
        if(any(vapply(prior, predicate, logical(1)))){
          add(paste0(field, "_1"))
          add(field)
        }
      }
    }
    return(roots)
  }
  if(is.prior.PET(prior) || is.prior.PEESE(prior)){
    add(if(is.prior.PET(prior)) "PET" else "PEESE")
  }else if(is.prior.spike_and_slab(prior)){
    append(.get_spike_and_slab_variable(prior), paste0(parameter, "_variable"))
    append(.get_spike_and_slab_inclusion(prior), paste0(parameter, "_inclusion"))
    add(paste0(parameter, "_indicator"))
    shape <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else if(is.prior.vector(prior)) prior$parameters$K else integer()
    add(parameter, shape)
  }else if(is.prior.mixture(prior)){
    add(paste0(parameter, "_indicator"))
    for(k in seq_along(prior)) append(prior[[k]], paste0(parameter, "_component_", k))
    shape <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else if(is.prior.vector(prior)) prior[[1L]]$parameters$K else integer()
    add(parameter, shape)
  }else if(is.prior.factor(prior) || is.prior.vector(prior)){
    K <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else prior$parameters$K
    add(parameter, K)
    if(is.prior.factor(prior) && (is.prior.treatment(prior) || is.prior.independent(prior))) return(roots)
    if(identical(prior$distribution, "dirichlet")){
      add(.JAGS_prior_dirichlet_eta_name(parameter), K)
    }else if(prior$distribution %in% c("mnormal", "mt")){
      add(paste0("prior_par1_", parameter), K)
      add(paste0("prior_par2_", parameter), c(K, K))
      if(identical(prior$distribution, "mt")){
        add(paste0("prior_par_s_", parameter))
        add(paste0("prior_par_z_", parameter), K)
      }
    }
  }else if(is.prior.simple(prior)) add(parameter)
  roots
}

.bt_ordered_total_emitted_roots <- function(prior, parameter){

  metadata <- .prior_ordered_metadata(prior)
  total_name <- .prior_ordered_total_name(parameter)
  if(metadata$theta_dim == 1L) return(.bt_prior_emitted_roots(prior$total, total_name))
  roots <- list()
  roots[[total_name]] <- as.integer(metadata$theta_dim)
  if(is.prior.spike_and_slab(prior$total)){
    roots[[paste0(total_name, "_variable")]] <- as.integer(metadata$theta_dim)
    inclusion <- .bt_prior_emitted_roots(.get_spike_and_slab_inclusion(prior$total), paste0(total_name, "_inclusion"))
    roots[names(inclusion)] <- inclusion
    roots[[paste0(total_name, "_indicator")]] <- integer()
  }
  roots
}
