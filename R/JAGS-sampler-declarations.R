# Sampler declarations of JAGS_fit() ('sampler_gates', 'sampler_blocks').
#
# A declaration is an observed marker node with a degenerate distribution at
# zero, added to the final model text with data value 0:
#   BayesTools_gate_<k>  ~ dbt_gate(<node>, <active>)
#   BayesTools_block_<k> ~ dbt_block(<node>, <id>)
# The density of a marker is constant in its parents, so the posterior is the
# same with and without markers; the samplers of the BayesTools JAGS module
# (src/samplers/) read them. A gate lets 'BayesTools::GatedSlice' draw a
# switched-off leaf of a publication-bias mixture branch directly from its
# prior; a block lets 'BayesTools::FactorSlice' update the scalar leaves of the
# listed priors jointly. Branches, leaves, node names, and activity tests are
# read from the emitters' declared structure, never from node-name patterns.

.bt_sampler_gate_prefix  <- "BayesTools_gate_"
.bt_sampler_block_prefix <- "BayesTools_block_"

.bt_sampler_declaration_stop <- function(...){

  stop(structure(
    class = c("BayesTools_sampler_declaration_error", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

.bt_check_sampler_declaration_arguments <- function(sampler_gates, sampler_blocks){

  if(!is.logical(sampler_gates) || length(sampler_gates) != 1L || is.na(sampler_gates)){
    .bt_sampler_declaration_stop("The 'sampler_gates' argument must be TRUE or FALSE.")
  }
  if(is.null(sampler_blocks) || is.function(sampler_blocks)){
    return(invisible(TRUE))
  }
  if(!.bt_sampler_valid_blocks(sampler_blocks)){
    .bt_sampler_declaration_stop(
      "The 'sampler_blocks' argument must be NULL, a list of character vectors of prior names, ",
      "or a function of the combined prior list that returns such a list or NULL."
    )
  }

  invisible(TRUE)
}

.bt_sampler_valid_blocks <- function(sampler_blocks){

  valid_block <- function(block){
    is.character(block) && length(block) > 0L && !anyNA(block) && all(nzchar(block))
  }
  is.list(sampler_blocks) && !is.prior(sampler_blocks) && length(sampler_blocks) > 0L &&
    all(vapply(sampler_blocks, valid_block, logical(1)))
}

# The blocks of 'sampler_blocks': the list itself, or the value of the function
# called once with the combined prior list (formula priors under their final
# names, with their attributes), validated like a list given directly.
.bt_sampler_resolve_blocks <- function(sampler_blocks, prior_list){

  if(!is.function(sampler_blocks)){
    return(sampler_blocks)
  }
  blocks <- sampler_blocks(prior_list)
  if(!is.null(blocks) && !.bt_sampler_valid_blocks(blocks)){
    .bt_sampler_declaration_stop(
      "The 'sampler_blocks' function returned an invalid value; it must return NULL or a list of ",
      "character vectors of prior names."
    )
  }

  blocks
}

# The sampler declarations of a fit: the model text and data with the marker
# nodes, the record stored as attr(fit, "sampler_declarations"), and whether a
# marker was emitted (the fit then requires the BayesTools JAGS module). With
# both arguments at their defaults the text and data are returned unchanged.
.bt_sampler_declarations <- function(model_text, data, prior_list, add_parameters,
                                     sampler_gates = FALSE, sampler_blocks = NULL){

  .bt_check_sampler_declaration_arguments(sampler_gates, sampler_blocks)
  if(!isTRUE(sampler_gates) && is.null(sampler_blocks)){
    return(list(model = model_text, data = data, declarations = NULL, uses_module = FALSE))
  }
  .bt_sampler_check_reserved_names(model_text, data, prior_list, add_parameters)

  gates <- if(isTRUE(sampler_gates)){
    .bt_sampler_gate_declarations(model_text, prior_list)
  }else{
    list(declared = .bt_sampler_gate_table(), declined = .bt_sampler_declined_table())
  }
  blocks <- .bt_sampler_block_declarations(.bt_sampler_resolve_blocks(sampler_blocks, prior_list), prior_list)

  gate_markers <- if(nrow(gates$declared) > 0L) paste0(.bt_sampler_gate_prefix, seq_len(nrow(gates$declared))) else character()
  block_markers <- if(nrow(blocks$declared) > 0L) paste0(.bt_sampler_block_prefix, seq_len(nrow(blocks$declared))) else character()
  gates$declared$marker <- gate_markers
  blocks$declared$marker <- block_markers

  marker_lines <- c(
    paste0(gate_markers, " ~ dbt_gate(", gates$declared$node, ", ", gates$declared$active, ")", recycle0 = TRUE),
    paste0(block_markers, " ~ dbt_block(", blocks$declared$node, ", ", blocks$declared$block, ")", recycle0 = TRUE)
  )
  markers <- c(gate_markers, block_markers)
  if(length(markers) > 0L){
    model_text <- .bt_sampler_insert_markers(model_text, marker_lines)
    if(is.null(data)){
      data <- list()
    }
    data <- c(data, stats::setNames(as.list(rep(0, length(markers))), markers))
  }

  list(
    model = model_text,
    data = data,
    declarations = list(
      gates = gates$declared[, c("marker", "node", "active", "prior", "branch")],
      blocks = blocks$declared[, c("marker", "node", "block", "prior")],
      declined = rbind(gates$declined, blocks$declined)
    ),
    uses_module = length(markers) > 0L
  )
}

.bt_sampler_gate_table <- function(){

  data.frame(node = character(), active = character(), prior = character(),
             branch = integer(), stringsAsFactors = FALSE)
}

# Declined declarations: a gate row names its prior and branch, a block row its
# block (position in 'sampler_blocks') and its priors; 'nodes' lists the leaves.
.bt_sampler_declined_table <- function(){

  data.frame(kind = character(), prior = character(), branch = integer(), block = integer(),
             nodes = character(), reason = character(), stringsAsFactors = FALSE)
}

# Marker node names are reserved: no data, prior, monitored parameter, or node
# of the model may use them.
.bt_sampler_check_reserved_names <- function(model_text, data, prior_list, add_parameters){

  prefixes <- c(.bt_sampler_gate_prefix, .bt_sampler_block_prefix)
  names <- unique(c(
    names(data), names(prior_list), add_parameters,
    .bt_sampler_tokens(model_text)
  ))
  reserved <- names[vapply(names, function(name) any(startsWith(name, prefixes)), logical(1))]
  if(length(reserved) > 0L){
    .bt_sampler_declaration_stop(
      "The model uses reserved BayesTools sampler node name '", reserved[1L], "'."
    )
  }

  invisible(TRUE)
}

# The JAGS identifiers of a syntax, comments excluded.
.bt_sampler_tokens <- function(syntax){

  syntax <- .JAGS_mask_comments(paste(syntax, collapse = "\n"))
  unlist(regmatches(
    syntax,
    gregexpr("(?<![A-Za-z0-9_.])[A-Za-z][A-Za-z0-9_.]*", syntax, perl = TRUE)
  ), use.names = FALSE)
}

.bt_sampler_lines <- function(syntax){

  lines <- trimws(unlist(strsplit(paste(syntax, collapse = "\n"), "\n", fixed = TRUE), use.names = FALSE))
  lines[nzchar(lines)]
}

# Inserts the marker statements before the closing brace of the model block,
# the last block of the syntax.
.bt_sampler_insert_markers <- function(model_text, marker_lines){

  closing <- max(gregexpr("}", .JAGS_mask_comments(model_text), fixed = TRUE)[[1L]])
  if(closing < 1L){
    stop("The JAGS model syntax has no closing brace.", call. = FALSE)
  }

  paste0(
    substr(model_text, 1L, closing - 1L),
    paste0("  ", marker_lines, "\n", collapse = ""),
    substr(model_text, closing, nchar(model_text))
  )
}


# Gates ---------------------------------------------------------------------------

# Whether a prior gives a sampled continuous node (the samplers take scalar
# continuous nodes only; a point prior gives a deterministic node).
.bt_sampler_continuous_prior <- function(prior){

  !is.prior.point(prior) && !is.prior.discrete(prior)
}

# The stochastic leaves of one weight-function component as its emitter writes
# them (.bt_dnode_omega_component_syntax()): the binary cumulative ratio, the
# cumulative Gamma auxiliaries, or the free independent weights when their prior
# is continuous.
.bt_sampler_weightfunction_leaves <- function(prior, component_id, global_cuts){

  node_names <- .weightfunction_component_node_names(
    prior = prior, component_id = component_id,
    global_cuts = global_cuts, force_one_sided = TRUE
  )
  J <- node_names$n_bins
  type <- prior$weights$type
  if(identical(type, "cumulative")){
    if(J == 2L){
      return(node_names$omega_ratio)
    }
    return(paste0(node_names$eta, "[", seq_len(J), "]"))
  }
  if(identical(type, "independent") && J > 1L && .bt_sampler_continuous_prior(prior$weights$prior)){
    name <- if(identical(prior$weights$scale, "log_omega")) node_names$log_omega else node_names$omega_local
    return(paste0(name, "[", 2:J, "]"))
  }

  character()
}

# The branches of a publication-bias mixture with gate candidates. Each branch
# lists its stochastic leaves, every node name it emits before the composition
# (leaves and deterministic intermediates), its own emitted statements, and the
# activity test the composition multiplies it by; 'composition' holds the
# composition statements of the mixture.
.bt_sampler_bias_branches <- function(prior, parameter){

  branch_info <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
  branches <- lapply(seq_along(branch_info), function(k){
    list(branch = k, leaves = character(), names = character(), own = character(), active = NULL)
  })
  composition <- character()

  if(any(has_selection) || any(has_phacking)){
    omega_spec <- .bt_dnode_omega_spec(prior)
    phacking_spec <- .bt_dnode_phacking_spec(omega_spec$branches)
    composition <- c(
      .bt_dnode_phacking_emit_vector_rule(phacking_spec),
      .bt_dnode_omega_emit_composition(omega_spec),
      .bt_dnode_phacking_emit_composition(phacking_spec)
    )
    for(k in seq_along(branch_info)){
      branch <- branches[[k]]
      branch$active <- phacking_spec$indicator_terms[k]
      branch$own <- c(
        .bt_dnode_omega_emit_branch(omega_spec, k),
        .bt_dnode_phacking_emit_branch(phacking_spec, k)
      )
      selection <- branch_info[[k]]$selection
      if(!is.null(selection)){
        node_names <- .weightfunction_component_node_names(
          prior = selection, component_id = k,
          global_cuts = omega_spec$global_cuts, force_one_sided = TRUE
        )
        branch$names <- c(branch$names, unlist(node_names[c(
          "omega_local", "omega_target", "log_omega", "omega_ratio", "eta", "std_eta"
        )], use.names = FALSE))
        branch$leaves <- c(branch$leaves, .bt_sampler_weightfunction_leaves(selection, k, omega_spec$global_cuts))
      }
      phacking <- branch_info[[k]]$phacking
      if(!is.null(phacking)){
        branch$names <- c(branch$names, vapply(
          c("alpha", "phack_kind", "pi_null", "beta_null", "phack_z_source", "phack_z_dest"),
          .selection_phacking_node_name, character(1),
          component_id = k, names = phacking_spec$names
        ))
        if(.bt_sampler_continuous_prior(phacking$alpha)){
          branch$leaves <- c(branch$leaves, .selection_phacking_node_name("alpha", k, phacking_spec$names))
        }
      }
      branches[[k]] <- branch
    }
  }

  for(term in c("PET", "PEESE")){
    node <- .bt_dnode_prior_mixture_bias_term(parameter, prior, term)
    if(is.null(node)){
      next
    }
    k <- node$spec$branch
    branch <- branches[[k]]
    branch$active <- paste0("equals(", node$spec$indicator, ", ", k, ")")
    branch$own <- c(branch$own, .JAGS_prior.bias_term_source(prior, term))
    branch$names <- c(branch$names, node$spec$source)
    if(.bt_sampler_continuous_prior(prior[[k]])){
      branch$leaves <- c(branch$leaves, node$spec$source)
    }
    composition <- c(composition, .bt_deterministic_node_emit(node))
    branches[[k]] <- branch
  }

  list(branches = branches, composition = .bt_sampler_lines(composition))
}

# Whether every composition statement that names a node of the branch uses it
# only in a term '<node> * <active>' or '<node>[i] * <active>'.
.bt_sampler_composition_gated <- function(line, names, active){

  parts <- strsplit(line, " <- ", fixed = TRUE)[[1L]]
  if(length(parts) != 2L || any(.bt_sampler_tokens(parts[1L]) %in% names)){
    return(FALSE)
  }
  for(term in strsplit(parts[2L], " + ", fixed = TRUE)[[1L]]){
    if(!any(.bt_sampler_tokens(term) %in% names)){
      next
    }
    factors <- strsplit(term, " * ", fixed = TRUE)[[1L]]
    operand <- sub("\\[[^]]*\\]$", "", factors[1L])
    if(length(factors) != 2L || !operand %in% names || !identical(factors[2L], active)){
      return(FALSE)
    }
  }

  TRUE
}

# Gate declarations for the stochastic leaves of every branch of each
# publication-bias mixture (with more than one branch) in the prior list. A
# branch is gated only when its nodes occur in the final model text only in its
# own statements and, multiplied by its activity test, in the composition
# statements (condition G); otherwise its gates are declined with the reason.
.bt_sampler_gate_declarations <- function(model_text, prior_list){

  declared <- .bt_sampler_gate_table()
  declined <- .bt_sampler_declined_table()
  decline <- function(name, branch, reason){
    declined[nrow(declined) + 1L, ] <<- list(
      "gate", name, as.integer(branch$branch), NA_integer_, paste(branch$leaves, collapse = ", "), reason
    )
  }

  for(name in names(prior_list)){
    prior <- prior_list[[name]]
    if(!inherits(prior, "prior.bias_mixture") || length(prior) < 2L){
      next
    }
    structure <- .bt_sampler_bias_branches(prior, name)
    candidates <- Filter(function(branch) length(branch$leaves) > 0L, structure$branches)
    if(length(candidates) == 0L){
      next
    }

    # the mixture's statements as JAGS_add_priors() emits them, once in the text
    piece <- .JAGS_add_priors.fun(prior_list[name])
    occurrences <- gregexpr(piece, model_text, fixed = TRUE)[[1L]]
    if(length(occurrences) != 1L || occurrences[1L] < 1L){
      for(branch in candidates){
        decline(name, branch, "the emitted statements of the mixture were not found once in the model syntax")
      }
      next
    }
    remainder <- paste0(
      substr(model_text, 1L, occurrences[1L] - 1L), " ",
      substr(model_text, occurrences[1L] + nchar(piece), nchar(model_text))
    )
    remainder_tokens <- unique(.bt_sampler_tokens(remainder))
    piece_lines <- .bt_sampler_lines(piece)

    for(branch in candidates){
      outside <- intersect(branch$names, remainder_tokens)
      if(length(outside) > 0L){
        decline(name, branch, paste0("node '", outside[1L], "' occurs outside the statements of the mixture"))
        next
      }
      own <- .bt_sampler_lines(branch$own)
      gated <- TRUE
      for(line in piece_lines){
        if(!any(.bt_sampler_tokens(line) %in% branch$names) || line %in% own){
          next
        }
        if(!(line %in% structure$composition &&
             .bt_sampler_composition_gated(line, branch$names, branch$active))){
          gated <- FALSE
          decline(name, branch, paste0("statement '", line, "' uses the branch outside its gated composition term"))
          break
        }
      }
      if(!gated){
        next
      }
      for(leaf in branch$leaves){
        declared[nrow(declared) + 1L, ] <- list(leaf, branch$active, name, as.integer(branch$branch))
      }
    }
  }

  list(declared = declared, declined = declined)
}


# Blocks --------------------------------------------------------------------------

# The scalar stochastic leaves of a prior that may join a factor slice block,
# or the reason it may not. The class guards come first, because a prior keeps
# its class when it also carries random-effect metadata.
.bt_sampler_block_leaves <- function(prior, name){

  reason <- function(text) list(leaves = NULL, reason = text)
  if(is.prior.spike_and_slab(prior)) return(reason("it is a spike-and-slab prior"))
  if(is.prior.mixture(prior))        return(reason("it is a mixture prior"))
  if(is.prior.ordered(prior))        return(reason("it is an ordered prior"))
  if(is.prior.R2(prior) || is.prior.R2_budget(prior)){
    return(reason("its R-squared budget value is deterministic"))
  }
  if(is.prior.point(prior))          return(reason("it is a point prior"))
  if(is.prior.discrete(prior))       return(reason("it is a discrete prior"))
  if(!is.null(attr(prior, "random_allocation_inclusion", exact = TRUE))){
    return(reason("it is an inclusion prior"))
  }

  # plain allocation primitives of random-effect variance allocation
  if(isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) &&
     is.prior.simple(prior) && !is.prior.factor(prior) && !is.prior.vector(prior)){
    return(list(leaves = name, reason = NULL))
  }
  if(!is.null(attr(prior, "random_allocation", exact = TRUE)) &&
     !isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) &&
     is.prior.vector(prior) && identical(prior[["distribution"]], "dirichlet")){
    return(list(
      leaves = paste0(.JAGS_prior_dirichlet_eta_name(name), "[", seq_len(prior$parameters[["K"]]), "]"),
      reason = NULL
    ))
  }
  # a direct random-effect SD or a scalar structured correlation with a simple
  # continuous prior is a leaf under its own name (without a summary and
  # without allocation metadata nothing else carries random-effect metadata)
  if(.bt_is_random_effect_prior(prior)){
    metadata <- .bt_random_effect_metadata(prior)
    if(is.prior.simple(prior) && !is.prior.factor(prior) && !is.prior.vector(prior) &&
       identical(metadata$summary, "") && !metadata$allocation){
      return(list(leaves = name, reason = NULL))
    }
    return(reason(paste0(
      "it is a random-effect prior other than a simple continuous scalar prior, ",
      "a plain allocation SD or Dirichlet allocation weights"
    )))
  }
  if(is.prior.vector(prior)) return(reason("it is a vector prior"))

  if(is.prior.factor(prior) && (is.prior.treatment(prior) || is.prior.independent(prior))){
    return(list(leaves = paste0(name, "[", seq_len(.get_prior_factor_levels(prior)), "]"), reason = NULL))
  }
  if(is.prior.weightfunction(prior)){
    spec <- .bt_dnode_omega_spec(prior)
    return(list(leaves = .bt_sampler_weightfunction_leaves(prior, NULL, spec$global_cuts), reason = NULL))
  }
  if(is.prior.PET(prior))   return(list(leaves = "PET", reason = NULL))
  if(is.prior.PEESE(prior)) return(list(leaves = "PEESE", reason = NULL))
  if(is.prior.simple(prior) && !is.prior.factor(prior)){
    return(list(leaves = name, reason = NULL))
  }

  reason("its prior type is not supported in a block")
}

# Block declarations: block k lists prior names of the combined prior list,
# each expanded to its scalar stochastic leaves. A block of fewer than two
# leaves is not declared but recorded as declined (a caller cannot count the
# leaves of its priors); guarded priors, unknown names and a prior in two
# blocks are errors.
.bt_sampler_block_declarations <- function(sampler_blocks, prior_list){

  blocks <- data.frame(node = character(), block = integer(), prior = character(), stringsAsFactors = FALSE)
  declined <- .bt_sampler_declined_table()
  if(is.null(sampler_blocks)){
    return(list(declared = blocks, declined = declined))
  }

  listed <- unlist(sampler_blocks, use.names = FALSE)
  if(anyDuplicated(listed)){
    .bt_sampler_declaration_stop(
      "Prior '", listed[duplicated(listed)][1L], "' appears more than once in 'sampler_blocks'; a prior may be in one block only."
    )
  }
  for(k in seq_along(sampler_blocks)){
    leaves <- character()
    block <- blocks[0L, ]
    for(name in sampler_blocks[[k]]){
      if(!name %in% names(prior_list)){
        .bt_sampler_declaration_stop(
          "Prior '", name, "' in 'sampler_blocks' is not in the combined prior list of the model."
        )
      }
      expansion <- .bt_sampler_block_leaves(prior_list[[name]], name)
      if(!is.null(expansion$reason)){
        .bt_sampler_declaration_stop(
          "Prior '", name, "' cannot join a block of 'sampler_blocks': ", expansion$reason, "."
        )
      }
      leaves <- c(leaves, expansion$leaves)
      block <- rbind(block, data.frame(
        node = expansion$leaves, block = rep(k, length(expansion$leaves)),
        prior = rep(name, length(expansion$leaves)), stringsAsFactors = FALSE
      ))
    }
    if(length(leaves) < 2L){
      declined[nrow(declined) + 1L, ] <- list(
        "block", paste(sampler_blocks[[k]], collapse = ", "), NA_integer_, as.integer(k),
        paste(leaves, collapse = ", "), "the block expands to fewer than two sampled nodes"
      )
      next
    }
    blocks <- rbind(blocks, block)
  }
  blocks$block <- as.integer(blocks$block)

  list(declared = blocks, declined = declined)
}
