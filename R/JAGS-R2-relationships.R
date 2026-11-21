# Fitted provenance for R-squared formula budgets. The normalized weights of
# a kernel term may have several units, while an exact-law term owns one.
.bt_R2_relationships <- function(add, quantities, prior_list){

  for(record in .bt_R2_budget_records(prior_list)){
    owned <- quantities[quantities$provider == "BayesTools" & quantities$owner_type == "r2_budget" &
                        quantities$formula_parameter == record$parameter & quantities$owner_name == record$name, , drop = FALSE]
    allocation <- .bt_relationship_vertex("r2_budget", record$parameter, record$name)
    budget <- prior_list[[record$terms[[1L]]$parameter]]$budget
    edge <- function(source_type, source, consumer_type, consumer, relation, term = ""){
      add(source_type, source, consumer_type, consumer, relation,
          formula_parameter = record$parameter, allocation = record$name,
          component = term, normalization = budget$law)
    }
    id <- function(kind) owned$quantity_id[owned$quantity == kind]
    total <- id("total")
    edge("quantity", total, "allocation", allocation, "budget_total")
    scale <- budget$sigma
    if(is.numeric(scale)){
      edge("constant", sprintf("%.17g", scale), "allocation", allocation, "residual_sd")
    }else{
      source <- .bt_relationship_source(quantities, scale)
      edge(source$type, source$id, "allocation", allocation, "residual_sd")
    }
    # All Gamma draws and realized term gates participate in normalization.
    edge("coordinate", record$gamma_coordinates, "allocation", allocation, "simplex_source")
    nodes <- if(!is.null(record$core)) .bt_R2_core_nodes(record$core)
    edge("coordinate", .bt_R2_record_gating(record), "allocation", allocation, "active_set")
    edge("coordinate", nodes$variable, "quantity", total, "slab_total")
    edge("coordinate", nodes$indicator, "quantity", total, "budget_gate")
    budget_event <- id("budget_inclusion")
    if(length(budget_event)){
      # Without a spike, the budget indicator is the OR of its primitive gates.
      event_sources <- if(is.null(record$core$gates$spike)) nodes$gates else nodes$indicator
      edge("coordinate", event_sources, "quantity", budget_event, "budget_inclusion")
      edge("coordinate", nodes$inclusion, "quantity", budget_event, "inclusion_probability")
    }
    rate <- id("inclusion_rate")
    # The collapsed gate law draws this rate conditional on the raw gates.
    edge("coordinate", nodes$rate_u, "quantity", rate, "rate_draw")
    edge("coordinate", nodes$gates, "quantity", rate, "conditions_rate")

    for(term in record$terms){
      component <- .bt_relationship_vertex("r2_budget", record$parameter, record$name, term$key)
      member <- owned[owned$term == term$label, , drop = FALSE]
      coefficient <- prior_list[[term$parameter]]
      edge("allocation", allocation, "allocation_component", component, "allocates", term$label)
      weights <- if(!is.null(record$allocation_name)){
        paste0(record$allocation_name, "[", term$units, "]")
      }
      edge("coordinate", weights, "allocation_component", component, "allocation_weight", term$label)
      edge("coordinate", .bt_R2_direction_coordinates(coefficient, term$parameter),
           "allocation_component", component, "coefficient_direction", term$label)
      gate <- if(!is.null(record$core)) .bt_R2_term_open_column(coefficient)
      edge("coordinate", gate, "allocation_component", component, "inclusion_gate", term$label)
      edge("allocation_component", component, "quantity", member$quantity_id[member$quantity == "coefficient"],
           "budget_coefficient", term$label)
      edge("allocation_component", component, "quantity", member$quantity_id[member$quantity == "share"],
           "share", term$label)
      event <- member$quantity_id[member$quantity == "term_inclusion"]
      edge("allocation_component", component, "quantity", event, "term_inclusion", term$label)
      if(length(event)){
        index <- record$core$term_gate[[term$key]]
        edge("coordinate", nodes$gates[[index]], "quantity", event, "raw_gate", term$label)
        if(!is.null(record$core$gates$spike)){
          edge("coordinate", nodes$indicator, "quantity", event, "budget_gate", term$label)
        }
      }
    }
  }
  invisible(NULL)
}
