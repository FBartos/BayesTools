# Internal random-effect summary helpers.

.bt_random_effect_summary_designs <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    formula_design <- list(formula_design)
  }
  if(!is.list(formula_design)){
    return(list())
  }

  formula_design[vapply(
    formula_design,
    .bt_formula_design_has_any_random_effects,
    logical(1)
  )]
}
