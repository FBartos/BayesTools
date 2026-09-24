# Random-effect deterministic node families.


# Scalar correlations ("random_rho") ------------------------------------------
#
# A structured random-effect block (cs, hcs, ar1, har, car) with a correlation
# prior on the Fisher-z or logit scale samples 'sample_name' and defines
#   rho <- tanh(sample_name)                                   (Fisher-z)
#   rho <- lower + (upper - lower) * ilogit(sample_name)       (logit)
# with the structure's correlation bounds. A raw ('rho' scale) correlation is
# the sampled node itself and has no generated node.

.bt_dnode_rho <- function(rho_name, sample_name, rho_scale, bounds,
                          sample_fixed = NULL, parameter = NA_character_,
                          block = NA_character_){

  if(!rho_scale %in% c("fisher_z", "logit")){
    stop("Scalar correlation nodes are generated only on the Fisher-z or logit scale.",
         call. = FALSE)
  }

  .bt_deterministic_node(
    family = "random_rho",
    node = rho_name,
    coordinates = rho_name,
    dependencies = sample_name,
    parameter = parameter,
    block = block,
    spec = list(
      sample_name = sample_name,
      rho_scale = rho_scale,
      bounds = c(lower = bounds[["lower"]], upper = bounds[["upper"]]),
      sample_fixed = sample_fixed
    )
  )
}

.bt_dnode_rho_from_random_term <- function(random_term, parameter = NA_character_){

  correlation <- random_term$correlation
  if(!is.list(correlation) || !identical(correlation$type, "rho") ||
     !correlation$rho_scale %in% c("fisher_z", "logit")){
    return(NULL)
  }
  plan <- .bt_random_effect_compile_rho_draw_plan(
    random_term = random_term,
    context = "Random-effect deterministic node metadata"
  )

  .bt_dnode_rho(
    rho_name = correlation$rho_name,
    sample_name = correlation$sample_name,
    rho_scale = plan$rho_scale,
    bounds = plan$bounds,
    sample_fixed = plan$sample_fixed,
    parameter = parameter,
    block = random_term$block_name
  )
}

.bt_dnode_rho_emit <- function(node){

  spec <- node$spec
  if(identical(spec$rho_scale, "fisher_z")){
    return(paste0(node$node, " <- tanh(", spec$sample_name, ")"))
  }

  paste0(
    node$node, " <- ",
    .bt_JAGS_numeric_literal(spec$bounds[["lower"]]), " + ",
    .bt_JAGS_numeric_literal(spec$bounds[["upper"]] - spec$bounds[["lower"]]),
    " * ilogit(", spec$sample_name, ")"
  )
}

# The map from the sampled coordinate to the correlation.
.bt_dnode_rho_transform <- function(value, rho_scale, bounds){

  if(identical(rho_scale, "fisher_z")){
    return(tanh(value))
  }
  if(identical(rho_scale, "logit")){
    return(
      bounds[["lower"]] +
        (bounds[["upper"]] - bounds[["lower"]]) * stats::plogis(value)
    )
  }

  value
}

.bt_dnode_rho_evaluate <- function(node, lookup){

  spec <- node$spec
  value <- .bt_deterministic_lookup_value(lookup, spec$sample_name)
  if(is.null(value) && !is.null(spec$sample_fixed)){
    value <- rep(spec$sample_fixed, lookup$n)
  }
  if(is.null(value)){
    return(NULL)
  }

  .bt_dnode_rho_transform(value, spec$rho_scale, spec$bounds)
}

