# Fit-local coordinates for one nested pair of cumulative weight functions.
# An active branch keeps its two weights below the common cutoffs. The finer
# branch also has an independent spacing within the split at the leading end.
.bt_nested_weights_plan <- function(prior){

  attr(prior, "BayesTools_nested_weights", exact = TRUE)
}

.bt_nested_weights_component_names <- function(prior, plan, fields){

  unlist(lapply(c(plan$coarse, plan$fine), function(k){
    selection <- .selection_branch_info(prior[[k]])$selection
    .weightfunction_component_node_names(selection, component_id = k)[fields]
  }), use.names = FALSE)
}

.bt_nested_weights_retired <- function(prior, plan = .bt_nested_weights_plan(prior)){

  if(is.null(plan)){
    return(character())
  }
  .bt_nested_weights_component_names(prior, plan, c("eta", "std_eta"))
}

.bt_nested_weights_bind <- function(prior_list, model_syntax, add_parameters,
                                    sampler_nested_weights, data = NULL){

  check_bool(sampler_nested_weights, "sampler_nested_weights", allow_NA = FALSE)
  bias <- prior_list[["bias"]]
  if(!is.null(.bt_nested_weights_plan(bias))){
    attr(bias, "BayesTools_nested_weights") <- NULL
    prior_list[["bias"]] <- bias
  }
  if(!sampler_nested_weights){
    return(list(prior_list = prior_list, status = NULL))
  }

  decline <- function(reason, branches = integer()){
    list(prior_list = prior_list,
         status = list(realized = FALSE, branches = branches, reason = reason))
  }
  if(!inherits(bias, "prior.bias_mixture")){
    return(decline("The fit has no publication-bias mixture."))
  }

  branch_info <- lapply(bias, .selection_branch_info)
  steps <- lapply(branch_info, function(branch){
    prior <- branch$selection
    if(is.null(prior) || !is.null(branch$phacking) ||
       !identical(prior$side, "one-sided") ||
       !identical(prior$weights$type, "cumulative") ||
       length(prior$weights$alpha) != length(prior$steps) + 1L ||
       !all(prior$weights$alpha == 1)){
      return(NULL)
    }
    prior$steps
  })
  coarse <- which(lengths(steps) == 2L)
  fine <- which(lengths(steps) == 3L)
  pairs <- list()
  for(i in coarse){
    for(j in fine){
      if(all(steps[[i]] == steps[[j]][-1L])){
        pairs[[length(pairs) + 1L]] <- list(coarse = i, fine = j)
      }
    }
  }
  if(length(pairs) == 0L){
    return(decline("The bias mixture has no eligible nested one-sided cumulative pair."))
  }
  if(length(pairs) > 1L){
    return(decline("The bias mixture has multiple eligible nested pairs."))
  }

  plan <- pairs[[1L]]
  branches <- c(coarse = plan$coarse, fine = plan$fine)
  plan$nodes <- c(a = "BayesTools_nested_a", v = "BayesTools_nested_v",
                  u = "BayesTools_nested_u")
  # Raw eta, standardized eta, and the branch-local weight aliases may be used
  # in caller syntax. Only the composed weights use these leaves in our chart.
  affected <- .bt_nested_weights_component_names(
    bias, plan, c("eta", "std_eta", "omega_local", "omega_target")
  )
  # Expression-valued parameters in other priors enter the graph only when
  # their priors are emitted, so include that source in the dependency guard.
  other_priors <- prior_list[names(prior_list) != "bias"]
  other_branches <- .bt_sampler_bias_branches(bias, "bias")$branches[-unname(branches)]
  other_branch_syntax <- unlist(lapply(other_branches, `[[`, "own"), use.names = FALSE)
  sources <- unique(c(names(data), names(prior_list),
                      .bt_sampler_tokens(c(model_syntax, add_parameters,
                                           .JAGS_add_priors.fun(other_priors), other_branch_syntax))))
  used <- intersect(sources, c(affected, plan$nodes))
  if(length(used) > 0L){
    return(decline(paste0("The fit uses or monitors the private weight node '", used[1L], "'."),
                   branches))
  }

  attr(bias, "BayesTools_nested_weights") <- plan
  prior_list[["bias"]] <- bias
  list(prior_list = prior_list,
       status = list(realized = TRUE, branches = branches, reason = NULL))
}

.bt_nested_weights_activity <- function(plan){

  c(coarse = paste0("equals(bias_indicator, ", plan$coarse, ")"),
    fine = paste0("equals(bias_indicator, ", plan$fine, ")"))
}

.bt_nested_weights_syntax <- function(plan){

  nodes <- plan$nodes
  active <- .bt_nested_weights_activity(plan)
  c(
    paste0(nodes[["a"]], " ~ dbeta(1 + ", active[["coarse"]], " + ", active[["fine"]],
           ", 1 + ", active[["fine"]], ")"),
    paste0(nodes[["v"]], " ~ dbeta(1, 1)"),
    paste0(nodes[["u"]], " ~ dbeta(1, 1)")
  )
}

.bt_nested_weights_local <- function(plan, branch, a, v, u = NULL){

  if(branch == plan$fine){
    return(cbind(1, a + (1 - a) * u, a, a * v))
  }
  cbind(1, a, a * v)
}

.bt_nested_weights_local_syntax <- function(plan, branch){

  a <- plan$nodes[["a"]]
  v <- plan$nodes[["v"]]
  u <- plan$nodes[["u"]]
  if(branch == plan$fine){
    return(c("1", paste0(a, " + (1 - ", a, ") * ", u), a, paste0(a, " * ", v)))
  }
  c("1", a, paste0(a, " * ", v))
}

.bt_nested_weights_init <- function(plan, init, branch_info){

  eta_name <- function(k){
    .weightfunction_component_node_names(branch_info[[k]]$selection, component_id = k)$eta
  }
  coarse_name <- eta_name(plan$coarse)
  fine_name <- eta_name(plan$fine)
  coarse <- init[[coarse_name]]
  fine <- init[[fine_name]]
  tail <- if(init[["bias_indicator"]] == plan$fine) fine[3:4] else coarse[2:3]
  values <- if(init[["bias_indicator"]] == plan$fine) fine else coarse
  a <- sum(tail) / sum(values)
  if(!init[["bias_indicator"]] %in% c(plan$coarse, plan$fine)){
    a <- a^2
  }
  v <- tail[2L] / sum(tail)
  u <- fine[2L] / sum(fine[1:2])
  init[c(coarse_name, fine_name)] <- NULL
  init[unname(plan$nodes)] <- list(a, v, u)
  init
}
