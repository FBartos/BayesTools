.bt_random_relationships <- function(add, quantities, design, formula_scale = NULL){

  parameter <- design$parameter
  vertex <- function(name, component = NULL){
    .bt_relationship_vertex("variance_allocation", parameter, name, component)
  }
  q <- quantities[quantities$formula_parameter == parameter & !quantities$internal, , drop = FALSE]
  blocks <- design$random_effects
  names(blocks) <- vapply(blocks, `[[`, character(1), "block_name")
  for(declaration in design$random_allocations){
    name <- declaration$label
    owner <- vertex(name)
    normalization <- declaration$scale
    link <- function(from_type, from, to_type, to, relation, component = ""){
      add(from_type, from, to_type, to, relation, parameter, name, component, normalization)
    }
    if(is.null(declaration$parent)){
      source <- declaration$source
      if(identical(source$shape, "row")){
        root <- list(type = "row_source", id = source$name)
        for(input in .bt_parameter_source_inputs(source)){
          ref <- .bt_relationship_source(quantities, input)
          link(ref$type, ref$id, root$type, root$id, "source_input")
        }
      }else{
        root <- .bt_relationship_source(quantities, source$name)
      }
      link(root$type, root$id, "allocation", owner,
        if(isTRUE(source$owned)) "owned_sd_source" else "sd_source")
    }else{
      parent <- declaration$parent
      link("allocation_component", vertex(parent$allocation, parent$component),
        "allocation", owner, "parent_component", parent$component)
    }

    local <- declaration
    per_sd <- identical(local$target, "sd_component")
    if(per_sd){
      # Only the fitted block knows which SD units a column allocation owns.
      block <- blocks[[unname(local$terms[[1L]])]]
      bound <- block$sd_binding$allocations
      match_owner <- vapply(bound, function(x) identical(x$label, name), logical(1))
      if(sum(match_owner) != 1L){
        .bt_relationship_stop("An SD-component allocation has no unique bound block record.")
      }
      local <- bound[[which(match_owner)]]
      components <- unname(local$leaf_terms)
    }else{
      components <- local$component_labels
    }
    owner_rows <- which(q$owner_type == "variance_allocation" &
      vapply(q$extraction_key, function(key) identical(key$allocation_label, name), logical(1)))
    for(row in owner_rows){
      key <- q$extraction_key[[row]]
      kind <- q$quantity[[row]]
      index <- switch(kind,
        var_prop = key$index, var_mult = key$index,
        sd_mult = key$index - length(components),
        inclusion = match(names(local$inclusion)[key$index], components))
      if(is.null(index)){
        link("allocation", owner, "quantity", q$quantity_id[row], kind)
      }else{
        component <- components[[index]]
        link("allocation_component", vertex(name, component), "quantity", q$quantity_id[row], kind, component)
      }
    }
    for(i in seq_along(components)){
      component <- components[[i]]
      port <- vertex(name, component)
      link("allocation", owner, "allocation_component", port, "allocates", component)
      if(length(local$weight_name)){
        weight <- .bt_dnode_prior_vector_names(local$weight_name, length(components))[[i]]
        link("coordinate", weight, "allocation_component", port, "allocation_weight", component)
      }
      gate <- local$inclusion[[component]]
      if(!is.null(gate)){
        link("coordinate", gate$indicator_name, "allocation_component", port, "inclusion_gate", component)
        link("coordinate", gate$prob_name, "allocation_component", port, "inclusion_probability", component)
      }
      target <- if(per_sd) unname(local$terms[[1L]]) else unname(local$components[[component]]$term)
      block <- blocks[[target]]
      if(is.null(block)) next
      # A parent can name the very same block that a nested SD allocation owns.
      # Only the child's declaration directly governs those SD units.
      owns <- vapply(block$sd_binding$allocations, function(x) identical(x$label, name), logical(1))
      if(!any(owns)) next
      block_id <- .bt_relationship_vertex("random_block", parameter, target)
      link("allocation_component", port, "random_block", block_id, "scales_block", component)
      block_rows <- which(q$owner_type == "random_block" & q$quantity %in% c("sd", "var") &
        vapply(q$extraction_key, function(key) identical(key$random_block, target), logical(1)))
      direct <- .bt_parameter_catalog_random_sd_is_direct(block, parameter, formula_scale)
      if(per_sd && direct){
        sd_names <- unique(block$sd_parameter_names[!is.na(block$sd_parameter_names)])
        block_rows <- block_rows[vapply(q$extraction_key[block_rows], function(key){
          identical(unname(sd_names[key$index]), unname(local$leaf_names[i]))
        }, logical(1))]
      }
      if(direct){
        for(kind in c("sd", "var")){
          leaves <- block_rows[q$quantity[block_rows] == kind]
          link("allocation_component", port, "quantity", q$quantity_id[leaves],
            paste0("block_", kind), component)
        }
      }else{
        # Back-transformed summaries belong to the covariance block. Its fitted
        # component shares cannot in general be assigned to individual marginals.
        link("random_block", block_id, "quantity", q$quantity_id[block_rows], "block_covariance")
      }
    }
  }
}
