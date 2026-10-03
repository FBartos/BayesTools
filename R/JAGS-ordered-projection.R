# Finite tensor contractions of persisted ordered recipes. Fixed simplex axes
# contract with declared weights. A random simplex axis disappears exactly when
# its coordinate coefficient tensors agree; no draw variability decides this.
.bt_ordered_tensor_reduce <- function(tensor){

  for(key in names(tensor$records)){
    axis <- match(key, names(tensor$records))
    dims <- vapply(tensor$records, `[[`, integer(1), "dim")
    permutation <- c(axis, setdiff(seq_along(dims), axis))
    coefficients <- matrix(aperm(array(tensor$coefficients, dim=dims), permutation), nrow=dims[[axis]])
    constant <- all(coefficients == matrix(rep(coefficients[1L,], each=dims[[axis]]), nrow=dims[[axis]]))
    record <- tensor$records[[key]]
    if(!constant && !identical(record$spec$type, "fixed")) next
    reduced <- if(constant) coefficients[1L,] else as.vector(crossprod(record$spec$weights, coefficients))
    tensor$records[[key]] <- NULL
    tensor$coefficients <- reduced
    # Re-start after removing an axis; the remaining axes keep their order.
    return(.bt_ordered_tensor_reduce(tensor))
  }
  tensor
}

.bt_ordered_tensor <- function(spec, slice, weights){

  indices <- which(spec$slice_index == slice)
  records <- lapply(spec$metadata$ordered_terms, function(factor){
    record <- .prior_ordered_allocation_for_coefficient(spec$metadata, factor, slice)
    spec$allocations[[record$key]]
  })
  names(records) <- vapply(records, `[[`, character(1), "key")
  .bt_ordered_tensor_reduce(list(records=records, coefficients=as.numeric(weights[indices])))
}

.bt_ordered_tensor_values <- function(tensor, allocations, n){

  if(!length(tensor$records)) return(rep(tensor$coefficients[[1L]], n))
  grid <- expand.grid(lapply(tensor$records, function(record) seq_len(record$dim)), KEEP.OUT.ATTRS=FALSE)
  out <- numeric(n)
  for(i in seq_along(tensor$coefficients)){
    if(tensor$coefficients[[i]] == 0) next
    value <- rep(tensor$coefficients[[i]], n)
    for(j in seq_along(tensor$records)){
      value <- value * allocations[[names(tensor$records)[[j]]]][, grid[i,j]]
    }
    out <- out + value
  }
  out
}

.bt_ordered_total_structure <- function(spec, draws, totals){

  states <- .bt_ordered_total_components(spec, draws)
  n <- nrow(draws)
  point <- matrix(FALSE, n, length(spec$total_names))
  unavailable <- rep(FALSE, n)
  location <- matrix(NA_real_, n, length(spec$total_names))
  for(k in unique(states$component)){
    prior <- states$priors[[k]]
    rows <- which(states$component == k)
    if(is.prior.point(prior) && .is_prior_expression(prior)){
      unavailable[rows] <- TRUE
    }else if(is.prior.point(prior)){
      point[rows,] <- TRUE
      location[rows,] <- prior$parameters$location
    }else if(is.prior.discrete(prior)){
      if(!identical(prior$distribution, "bernoulli")){
        .bt_ordered_stop("Ordered discrete total structure is unavailable for this distribution. Use a supported scalar total prior.")
      }
      values <- totals[rows,,drop=FALSE]
      if(!is.numeric(prior$truncation$lower) || !is.numeric(prior$truncation$upper)){
        unavailable[rows] <- TRUE
        next
      }
      support <- c(0, 1)
      support <- support[support >= prior$truncation$lower & support <= prior$truncation$upper]
      invalid <- any(!is.finite(values)) || any(!values %in% support)
      if(!.is_prior_expression(prior)) invalid <- invalid || any(lpdf(prior, as.vector(values)) == -Inf)
      if(invalid){
        .bt_ordered_stop("Ordered discrete total draws do not belong to the declared support.", "BayesTools_ordered_invalid_state")
      }
      point[rows,] <- TRUE
      location[rows,] <- values
    }
  }
  list(point=point, location=location, unavailable=unavailable, component=states$component)
}

.bt_ordered_projection <- function(specs, weights, draws = NULL, prior_list = list()){

  check_real(weights, "weights", check_length=0, allow_NA=FALSE)
  if(is.null(names(weights)) || anyDuplicated(names(weights)) || any(!is.finite(weights))){
    stop("'weights' must be finite named fitted-coordinate weights.", call. = FALSE)
  }
  requested_weights <- weights
  weights <- .bt_ordered_expand_total_weights(specs,weights)
  tensors <- list()
  for(spec in specs){
    local_weights <- stats::setNames(rep(0, length(spec$coefficient_names)), spec$coefficient_names)
    present <- intersect(names(weights), spec$coefficient_names)
    local_weights[present] <- weights[present]
    for(slice in seq_along(spec$total_names)){
      tensors[[paste0(spec$parameter,"|",slice)]] <- list(parameter=spec$parameter, slice=slice,
        identity=length(unique(local_weights[spec$slice_index==slice]))==1L,
        tensor=.bt_ordered_tensor(spec,slice,local_weights))
    }
  }
  ordinary <- weights[!names(weights) %in% unlist(lapply(specs, `[[`, "coefficient_names"), use.names=FALSE) & weights != 0]
  out <- list(specs=specs, weights=weights, requested_weights=requested_weights,
    contractions=tensors, ordinary_weights=ordinary)
  class(out) <- c("BayesTools_ordered_projection", "list")
  if(is.null(draws)) return(out)
  draws <- .bt_deterministic_draws_matrix(draws)
  n <- nrow(draws)
  values <- numeric(n)
  structures <- list()
  allocations <- list()
  totals <- list()
  total_contractions <- list()
  allocation_contractions <- list()
  for(spec in specs){
    total <- .bt_ordered_total_values(spec, .bt_deterministic_lookup(draws,prior_list))
    if(is.null(total)) .bt_ordered_stop(paste0("Ordered total sources for '",spec$parameter,"' are unavailable in 'draws'. Include its declared source coordinates."))
    totals[[spec$parameter]] <- total
    total_contractions[[spec$parameter]] <- matrix(0,n,length(spec$total_names),
      dimnames=list(NULL,spec$total_names))
    structures[[spec$parameter]] <- .bt_ordered_total_structure(spec,draws,total)
    for(record in spec$allocations){
      needed <- any(vapply(tensors,function(term) record$key %in% names(term$tensor$records),logical(1)))
      if(!needed) next
      allocation <- .bt_ordered_allocation_values(record,draws)
      if(is.null(allocation)) .bt_ordered_stop(paste0("Ordered allocation sources for '",spec$parameter,"' are unavailable in 'draws'. Include its declared source coordinates."))
      allocations[[record$key]] <- allocation
    }
  }
  for(term in tensors){
    contraction <- .bt_ordered_tensor_values(term$tensor, allocations, n)
    total_contractions[[term$parameter]][,term$slice] <- contraction
    values <- values + totals[[term$parameter]][,term$slice] * contraction
    for(key in names(term$tensor$records)){
      records <- term$tensor$records
      axis <- match(key,names(records))
      dimensions <- vapply(records,`[[`,integer(1),"dim")
      coefficients <- matrix(aperm(array(term$tensor$coefficients,dimensions),
        c(axis,setdiff(seq_along(dimensions),axis))),nrow=dimensions[[axis]])
      record <- records[[key]]
      if(is.null(allocation_contractions[[key]])) allocation_contractions[[key]] <- matrix(0,n,record$dim,
        dimnames=list(NULL,record$coordinates))
      for(j in seq_len(record$dim)){
        remaining <- records
        remaining[[key]] <- NULL
        tensor <- .bt_ordered_tensor_reduce(list(records=remaining,coefficients=coefficients[j,]))
        allocation_contractions[[key]][,j] <- allocation_contractions[[key]][,j] +
          totals[[term$parameter]][,term$slice] * .bt_ordered_tensor_values(tensor,allocations,n)
      }
    }
  }
  ordinary_point <- rep(TRUE,n)
  ordinary_location <- numeric(n)
  ordinary_unavailable <- rep(FALSE,n)
  for(name in names(ordinary)){
    owner <- .bt_parameter_coordinate_owner(prior_list,name)
    prior <- if(!is.null(owner)) prior_list[[owner]] else NULL
    if(is.null(prior)) .bt_ordered_stop(paste0("Projection source '",name,
      "' has unavailable structural provenance. Supply weights on fitted coefficients with declared priors."))
    source <- .bt_deterministic_lookup_value(.bt_deterministic_lookup(draws,prior_list),name)
    if(is.null(source)) .bt_ordered_stop(paste0("Projection source '",name,"' is unavailable in 'draws'. Include its fitted coordinate."))
    values <- values + ordinary[[name]] * source
    components <- rep(1L,n)
    if(!is.null(prior) && is.prior.mixture(prior)){
      indicator <- paste0(owner,"_indicator")
      if(!indicator %in% colnames(draws)) .bt_ordered_stop(paste0("Projection component indicator '",indicator,"' is unavailable. Include its monitored coordinate."))
      components <- .bt_component_from_indicator(prior,draws[,indicator])
    }
    for(component in unique(components)){
      rows <- which(components==component)
      active_prior <- if(!is.null(prior)) .posterior_atoms_component_prior(prior,component,FALSE) else NULL
      point <- if(!is.null(active_prior)) .posterior_atoms_point_location(active_prior,1L) else NULL
      if(is.null(point)) ordinary_point[rows] <- FALSE else{
        ordinary_location[rows] <- ordinary_location[rows] + ordinary[[name]] * point
      }
      if(!is.null(active_prior) && is.prior.point(active_prior) && .is_prior_expression(active_prior)){
        ordinary_unavailable[rows] <- TRUE
      }
    }
  }
  # State groups depend only on declared component/discrete sources. Tensor
  # cancellation across shared allocations is evaluated once per such group.
  signatures <- paste(ordinary_point,match(ordinary_location,unique(ordinary_location)),ordinary_unavailable,sep="|")
  for(structure in structures){
    signatures <- paste(signatures,structure$component,structure$unavailable,sep="|")
    for(i in seq_len(ncol(structure$point))){
      location <- structure$location[,i]
      signatures <- paste(signatures,structure$point[,i],match(location,unique(location)),sep="|")
    }
  }
  atom <- rep(NA_real_,n)
  unavailable <- rep(FALSE,n)
  for(signature in unique(signatures)){
    rows <- which(signatures==signature)
    first <- rows[[1L]]
    fixed <- list()
    continuous <- !ordinary_point[[first]]
    unclassified <- ordinary_unavailable[[first]]
    for(term in tensors){
      tensor <- term$tensor
      if(all(tensor$coefficients==0)) next
      structure <- structures[[term$parameter]]
      if(structure$unavailable[[first]]){
        unclassified <- TRUE
        next
      }
      if(!structure$point[first,term$slice]){
        continuous <- TRUE
        next
      }
      scalar <- structure$location[first,term$slice]
      if(scalar==0) next
      key <- paste(names(tensor$records),collapse="\r")
      if(!nzchar(key)) key <- "<constant>"
      tensor$coefficients <- scalar * tensor$coefficients
      if(is.null(fixed[[key]])) fixed[[key]] <- tensor else{
        fixed[[key]]$coefficients <- fixed[[key]]$coefficients + tensor$coefficients
      }
    }
    location <- ordinary_location[[first]]
    for(tensor in fixed){
      tensor <- .bt_ordered_tensor_reduce(tensor)
      if(length(tensor$records)) continuous <- TRUE else location <- location + tensor$coefficients[[1L]]
    }
    unavailable[rows] <- unclassified
    if(!continuous && !unclassified) atom[rows] <- location
  }
  on_atom <- !is.na(atom)
  values[on_atom] <- atom[on_atom]
  out$values <- values
  out$total_contraction_weights <- total_contractions
  out$allocation_contractions <- allocation_contractions
  out$exact <- ordinary_point & all(vapply(tensors,`[[`,logical(1),"identity"))
  out$atom <- atom
  out$state <- ifelse(unavailable,"unavailable",ifelse(on_atom,"point","continuous"))
  out$reason <- if(any(unavailable)) .bt_ordered_expression_reason() else NULL
  out
}

.bt_ordered_expression_reason <- function(){

  structure(list(message=paste0("Ordered scalar measure is unavailable for expression totals with unclassified stochastic ancestry. ",
    "Use a supported scalar total prior with numeric point locations, or inspect fitted snapshot values with 'parameter_draws()'."),call=NULL),
    class=c("BayesTools_ordered_expression_unavailable","BayesTools_ordered_unavailable","error","condition"))
}

.bt_ordered_require_measure <- function(projection){

  if(!is.null(projection$reason)) stop(projection$reason)
  invisible(NULL)
}

.bt_ordered_source_require_measure <- function(x){

  source <- .bt_meta_get(x,"ordered_source")
  if(is.null(source)) return(invisible(NULL))
  templates <- which(vapply(source$models,function(spec) is.null(spec$parameterization),logical(1)))
  if(!length(templates)) return(invisible(NULL))
  columns <- if(is.null(source$projection_design)){
    source$models[[templates[[1L]]]]$coefficient_names
  }else rownames(source$projection_design)
  for(i in seq_along(columns)){
    .bt_ordered_require_measure(.bt_ordered_source_project(source,diag(length(columns))[i,]))
  }
  invisible(NULL)
}

.bt_ordered_expand_total_weights <- function(specs, weights){

  for(spec in specs){
    selected <- intersect(names(weights),spec$total_names)
    for(total in selected){
      slice <- match(total,spec$total_names)
      coefficients <- spec$coefficient_names[spec$slice_index==slice]
      existing <- intersect(names(weights),coefficients)
      added <- setdiff(coefficients,names(weights))
      weights[existing] <- weights[existing] + weights[[total]]
      weights[added] <- weights[[total]]
      weights <- weights[names(weights)!=total]
    }
  }
  weights
}

.bt_ordered_projection_atoms <- function(projection, column, model = NULL, probabilities = NULL){

  if(any(projection$state=="unavailable")) return(NULL)
  locations <- sort(unique(projection$atom[!is.na(projection$atom)]))
  .posterior_atoms_new(matrix(locations,ncol=1L),
    vapply(locations,function(value){
      rows <- !is.na(projection$atom) & projection$atom==value
      .bt_ordered_state_mass(rows,model,probabilities)
    },numeric(1)),
    column_names=column,source="ordered_projection_structure")
}

.bt_ordered_state_mass <- function(rows, model=NULL, probabilities=NULL){

  if(is.null(model)) return(mean(rows))
  sum(vapply(seq_along(probabilities),function(i){
    selected <- model==i
    if(!any(selected)) 0 else probabilities[[i]] * mean(rows[selected])
  },numeric(1)))
}

.bt_ordered_quantity_plan <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  if(!key$type %in% c("coordinate","factor_level")) return(NULL)
  priors <- attr(fit,"prior_list",exact=TRUE)
  weights <- stats::setNames(if(identical(key$type,"coordinate")) 1 else key$weights,key$dependencies)
  ordered <- priors[vapply(priors,is.prior.ordered,logical(1))]
  relevant <- vapply(names(ordered),function(parameter){
    any(names(weights)[weights!=0] %in% c(.JAGS_prior_factor_names(parameter,ordered[[parameter]]),
      .prior_ordered_total_monitor_names(ordered[[parameter]],parameter)))
  },logical(1))
  if(!any(relevant)) return(NULL)
  specs <- JAGS_ordered_parameter_spec(fit,names(ordered)[relevant])
  weights <- .bt_ordered_expand_total_weights(specs,weights)
  events <- names(ordered)[relevant & vapply(ordered,function(prior){
    is.prior.mixture(prior$total) && all(attr(prior$total,"components",exact=TRUE) %in% c("null","alternative"))
  },logical(1))]
  list(kind="ordered_projection",specs=specs,weights=weights,prior_list=priors,
    chain_gates=character(),event_gates=events,event_rule="AND")
}

.bt_ordered_quantity_projection <- function(fit, quantity, draws=NULL){

  plan <- .bt_ordered_quantity_plan(fit,quantity)
  if(is.null(plan)) return(NULL)
  if(is.null(draws)) draws <- as.matrix(.fit_to_posterior(fit))
  .bt_ordered_projection(plan$specs,plan$weights,draws,plan$prior_list)
}

.bt_ordered_source_project <- function(source, weights){

  if(!is.null(source$projection_design)){
    weights <- as.vector(as.numeric(weights) %*% source$projection_design)
    names(weights) <- colnames(source$projection_design)
  }
  values <- atom <- rep(NA_real_,length(source$model))
  exact <- rep(FALSE,length(source$model))
  state <- rep("unavailable",length(source$model))
  reason <- NULL
  for(i in seq_along(source$models)){
    rows <- which(source$model==i)
    if(!length(rows)) next
    spec <- source$models[[i]]
    if(!is.null(source$projection_context)){
      context <- source$projection_context
      model <- context$models[[i]]
      projection <- .bt_ordered_projection(model$specs,weights,context$primitives[rows,,drop=FALSE],model$priors)
      values[rows] <- projection$values
      exact[rows] <- projection$exact
      atom[rows] <- projection$atom
      state[rows] <- projection$state
      if(!is.null(projection$reason)) reason <- projection$reason
      next
    }
    if(identical(spec$parameterization,"absent")){
      values[rows] <- atom[rows] <- 0
      state[rows] <- "point"
      next
    }
    if(!is.null(spec$parameterization)){
      point <- .posterior_atoms_point_location(spec$prior,length(weights))
      if(!is.null(point)){
        values[rows] <- atom[rows] <- sum(as.numeric(weights)*point)
        state[rows] <- "point"
      }
      next
    }
    own <- if(is.null(names(weights))) stats::setNames(as.numeric(weights),spec$coefficient_names) else weights
    projection <- .bt_ordered_projection(setNames(list(spec),spec$parameter),own,
      source$primitives[rows,,drop=FALSE])
    values[rows] <- projection$values
    exact[rows] <- projection$exact
    atom[rows] <- projection$atom
    state[rows] <- projection$state
    if(!is.null(projection$reason)) reason <- projection$reason
  }
  for(map in source$view_transformations){
    if(identical(map$transformation,"unavailable")){
      state[] <- "unavailable"
      atom[] <- NA_real_
      exact[] <- FALSE
      next
    }
    values <- .density.prior_transformation_x(values,map$transformation,map$arguments)
    atom <- .density.prior_transformation_x(atom,map$transformation,map$arguments)
  }
  list(values=values,atom=atom,state=state,exact=exact,reason=reason)
}

.bt_ordered_projection_defined_rows <- function(projection, values){

  projection$state %in% c("point", "continuous") &
    !is.na(projection$values) & !is.na(values)
}

.bt_ordered_source_semantics <- function(x, design, columns = colnames(x)){

  source <- .bt_meta_get(x,"ordered_source")
  if(is.null(source)) return(x)
  existing <- .posterior_atoms_get(x)
  if(is.null(source$model_probabilities)){
    probabilities <- existing$component_probabilities
    source$model_probabilities <- if(length(probabilities)==length(source$models)) probabilities else{
      tabulate(source$model,nbins=length(source$models))/length(source$model)
    }
  }
  projections <- lapply(seq_len(nrow(design)),function(i) .bt_ordered_source_project(source,design[i,]))
  marginals <- lapply(seq_along(projections),function(i) .bt_ordered_projection_atoms(projections[[i]],columns[[i]],source$model,source$model_probabilities))
  names(marginals) <- columns
  values <- .bt_draws_plain(x)
  for(i in seq_along(projections)){
    # Every available semantic value follows the primitive projection in its
    # declared view, including continuous contractions. Raw monitors stay intact.
    replace <- .bt_ordered_projection_defined_rows(projections[[i]],
      if(is.null(dim(values))) values else values[,i])
    if(is.null(dim(values))) values[replace] <- projections[[i]]$values[replace] else values[replace,i] <- projections[[i]]$values[replace]
  }
  x <- .bt_draws_transform_values(x,function(old) values)
  old_design <- source$projection_design
  if(is.null(old_design)){
    template <- source$models[[which(vapply(source$models,function(spec) is.null(spec$parameterization),logical(1)))[1L]]]
    old_design <- diag(length(template$coefficient_names))
    colnames(old_design) <- template$coefficient_names
  }
  source$projection_design <- design %*% old_design
  rownames(source$projection_design) <- columns
  x <- .bt_meta_set(x,"ordered_source",source)
  atoms <- .posterior_atoms_get(x)
  if(is.null(atoms)) atoms <- .posterior_atoms_new(column_names=columns)
  joint_locations <- do.call(cbind,lapply(projections,`[[`,"atom"))
  joint_rows <- rowSums(is.na(joint_locations))==0
  if(all(!vapply(marginals,is.null,logical(1))) && nrow(atoms$locations)==0L){
    locations <- unique(joint_locations[joint_rows,,drop=FALSE])
    mass <- vapply(seq_len(nrow(locations)),function(i){
      matching <- joint_rows
      for(j in seq_len(ncol(joint_locations))){
        matching <- matching & !is.na(joint_locations[,j]) & joint_locations[,j]==locations[i,j]
      }
      .bt_ordered_state_mass(matching,
        source$model,source$model_probabilities)
    },numeric(1))
    atoms <- .posterior_atoms_new(locations,mass,column_names=columns,source="ordered_projection_joint_structure",
      component_probabilities=atoms$component_probabilities)
  }
  if(nrow(atoms$locations)>0L){
    current_locations <- unique(joint_locations[joint_rows,,drop=FALSE])
    if(nrow(current_locations)==1L && nrow(unique(atoms$locations))==1L){
      atoms$locations[,] <- rep(current_locations[1,],each=nrow(atoms$locations))
    }
    joint <- atoms
    joint$marginals <- NULL
    for(i in seq_along(marginals)){
      if(!is.null(marginals[[i]]) && identical(!is.na(projections[[i]]$atom),joint_rows)){
        marginals[i] <- list(.posterior_atoms_for_column(joint,i))
      }
    }
  }
  atoms$marginals <- marginals
  .posterior_atoms_set(x,atoms)
}

.bt_ordered_formula_projections <- function(samples, weights, source_transforms=NULL){

  sources <- lapply(samples,function(x) .bt_meta_get(x,"ordered_source"))
  ordered_parameters <- names(sources)[!vapply(sources,is.null,logical(1))]
  if(!length(ordered_parameters)) return(NULL)
  if(any(vapply(sources[ordered_parameters],function(source) !is.null(source$view_transformations),logical(1)))) return(NULL)
  if(length(source_transforms)){
    # Existing log-source atom recipes own this unsupported combined view.
    return(NULL)
  }
  weights <- as.matrix(weights)
  source <- sources[[ordered_parameters[[1L]]]]
  n <- length(source$model)
  projections <- lapply(seq_len(nrow(weights)),function(i){
    list(values=rep(NA_real_,n),atom=rep(NA_real_,n),state=rep("unavailable",n),exact=rep(FALSE,n))
  })
  contexts <- vector("list",length(source$models))
  source_inputs <- vector("list",length(source$models))
  for(model in unique(source$model)){
    rows <- which(source$model==model)
    specs <- priors <- list()
    inputs <- list()
    absent <- character()
    for(parameter in names(samples)){
      x <- samples[[parameter]]
      prior <- attr(x,"prior_list",exact=TRUE)
      priors[[parameter]] <- if(is.prior(prior)) prior else prior[[model]]
      if(!is.null(sources[[parameter]])){
        own <- sources[[parameter]]
        if(!identical(own$model,source$model)) .bt_ordered_stop("Ordered formula sources are unavailable because their model rows do not align. Recreate mixed posteriors with aligned sampling.")
        spec <- own$models[[model]]
        if(!is.null(spec$parameterization)){
          if(identical(spec$parameterization,"absent")) absent <- c(absent,.posterior_atoms_coefficient_columns(x,parameter)) else return(NULL)
        }else{
          specs[[parameter]] <- spec
          coordinates <- unique(c(spec$total_names,
            unlist(lapply(spec$allocations,`[[`,"coordinates"),use.names=FALSE),spec$total_node$spec$indicator))
          inputs[[length(inputs)+1L]] <- own$primitives[rows,coordinates,drop=FALSE]
        }
      }else{
        names <- .posterior_atoms_coefficient_columns(x,parameter)
        values <- matrix(as.numeric(x),nrow=NROW(x))
        if(length(names)!=ncol(values) || nrow(values)!=n) return(NULL)
        colnames(values) <- names
        inputs[[length(inputs)+1L]] <- values[rows,,drop=FALSE]
        if(is.prior.mixture(priors[[parameter]]) && .bt_meta_get(x,"component_source") %in% c("mixture","spike_and_slab")){
          component <- .bt_meta_get(x,"component")[rows]
          indicator <- if(is.prior.spike_and_slab(priors[[parameter]])){
            as.numeric(attr(priors[[parameter]],"components",exact=TRUE)[component]=="alternative")
          }else component
          inputs[[length(inputs)+1L]] <- matrix(indicator,ncol=1L,
            dimnames=list(NULL,paste0(parameter,"_indicator")))
        }
      }
    }
    if(!length(inputs)) next
    draws <- do.call(cbind,inputs)
    draws <- draws[,!duplicated(colnames(draws)),drop=FALSE]
    contexts[[model]] <- list(specs=specs,priors=priors)
    source_inputs[[model]] <- draws
    for(i in seq_len(nrow(weights))){
      target <- stats::setNames(as.numeric(weights[i,]),colnames(weights))
      target[names(target) %in% absent] <- 0
      projection <- .bt_ordered_projection(specs,target,draws,priors)
      for(field in c("values","atom","state","exact")) projections[[i]][[field]][rows] <- projection[[field]]
      if(!is.null(projection$reason)) projections[[i]]$reason <- projection$reason
    }
  }
  columns <- unique(unlist(lapply(source_inputs,colnames),use.names=FALSE))
  primitives <- matrix(NA_real_,n,length(columns),dimnames=list(NULL,columns))
  for(model in unique(source$model)){
    rows <- which(source$model==model)
    if(!is.null(source_inputs[[model]])) primitives[rows,colnames(source_inputs[[model]])] <- source_inputs[[model]]
  }
  attr(projections,"context") <- list(models=contexts,primitives=primitives)
  counts <- tabulate(source$model,nbins=length(source$models))
  probabilities <- if(is.null(source$model_probabilities)) counts/sum(counts) else source$model_probabilities
  attr(projections,"model") <- source$model
  attr(projections,"probabilities") <- probabilities
  projections
}

.bt_ordered_attach_linear_view <- function(x, projections, weights, samples){

  if(is.null(projections)) return(x)
  sources <- lapply(samples,function(draws) .bt_meta_get(draws,"ordered_source"))
  source <- sources[[which(!vapply(sources,is.null,logical(1)))[[1L]]]]
  columns <- if(is.null(dim(x))) attr(x,"parameter",exact=TRUE) else colnames(x)
  if(is.null(columns) || !length(columns)) columns <- "value"
  source$projection_context <- attr(projections,"context",exact=TRUE)
  source$projection_design <- weights
  rownames(source$projection_design) <- columns
  x <- .bt_meta_set(x,"ordered_source",source)
  .bt_ordered_source_semantics(x,diag(length(columns)),columns)
}
