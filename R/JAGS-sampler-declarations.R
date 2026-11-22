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

#' Resolve prior sources for JAGS sampler declarations
#'
#' @param prior_list Named list of priors with formula priors already bound.
#'   A `sampler_blocks` callback in [JAGS_fit()] receives this list.
#' @param priors Character vector of names in `prior_list`, in the desired
#'   order. Defaults to all priors.
#' @param kind Type of sampler declaration. Currently `"block"`.
#' @return A data frame with columns `prior`, `node`, `status` and `reason`.
#'   A status of `"eligible"` identifies an emitted scalar stochastic node
#'   of a supported prior. A status of `"empty"` identifies a supported prior
#'   with no such nodes, and `"declined"` identifies a prior that cannot
#'   join a block. In those two cases `node` is missing and `reason` explains
#'   the status.
#' @details Sources are resolved before model compilation using the same
#'   rules as `sampler_blocks`. A source name is the emitted node name. The
#'   compiled graph determines the actual sampler and its reported node name.
#' @export
JAGS_sampler_sources <- function(prior_list, priors = names(prior_list), kind = "block"){

  if(!identical(kind, "block")){
    .bt_sampler_declaration_stop("The 'kind' argument must be 'block'.")
  }
  if(!is.list(prior_list) || is.prior(prior_list) ||
     (length(prior_list) > 0L && (is.null(names(prior_list)) ||
      anyNA(names(prior_list)) || any(!nzchar(names(prior_list))) ||
      anyDuplicated(names(prior_list)) || !all(vapply(prior_list, is.prior, logical(1)))))){
    .bt_sampler_declaration_stop("The 'prior_list' argument must be a named list of priors with unique names.")
  }
  if(is.null(priors) && length(prior_list) == 0L){
    priors <- character()
  }
  if(!is.character(priors) || anyNA(priors) || anyDuplicated(priors) ||
     any(!priors %in% names(prior_list))){
    .bt_sampler_declaration_stop("The 'priors' argument must contain distinct names in 'prior_list'.")
  }

  if(length(priors) == 0L){
    return(data.frame(prior = character(), node = character(), status = character(), reason = character(),
                      stringsAsFactors = FALSE))
  }
  rows <- lapply(unname(priors), function(name){
    expansion <- .bt_sampler_block_leaves(prior_list[[name]], name)
    if(!is.null(expansion$reason)){
      return(data.frame(prior = name, node = NA_character_, status = "declined", reason = expansion$reason,
                        stringsAsFactors = FALSE))
    }
    if(length(expansion$leaves) == 0L){
      return(data.frame(prior = name, node = NA_character_, status = "empty",
                        reason = "it has no supported scalar stochastic nodes", stringsAsFactors = FALSE))
    }
    data.frame(prior = name, node = expansion$leaves, status = "eligible", reason = NA_character_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Resolve the emitted sources of a JAGS mixture
#'
#' @param prior_list Named list of priors with formula priors already bound.
#' @param name Name of an ordinary or publication-bias mixture in `prior_list`.
#' @return A list with the emitted `gate` node name, the branch
#'   `probabilities`, and `branches` in the order indexed by the gate.
#'   Each branch contains its `prior`, `family`, `source`, `coordinates`,
#'   `stochastic`, `location`, and `omega`.
#' @details For ordinary mixtures the family identifies a `"point"`,
#'   `"scalar"`, `"vector"`, `"scalar_array"` (separately sampled factor
#'   coefficients), or `"unsupported"` source. For publication-bias mixtures
#'   it identifies `"none"`, `"selection"`, `"PET"`, `"PEESE"`, `"phack"`,
#'   or `"combined"`. The `source` and ordered `coordinates` refer to the
#'   emitted component before applying the mixture gate. The source name is
#'   missing for bias branches without a raw slope. The `stochastic` field
#'   is true when the source is emitted as a single stochastic node, false
#'   for constant and aggregate sources, and missing for unsupported sources.
#'   A literal point value appears in `location`; a prior expression has
#'   no literal location. For selection branches, `omega` gives the exact
#'   physical weight references in the order of the global bins from
#'   [selection_backend_spec()]. Other branches have no weight references.
#' @export
JAGS_mixture_sources <- function(prior_list, name){

  if(!is.list(prior_list) || !is.character(name) || length(name) != 1L ||
     is.na(name) || !name %in% names(prior_list)){
    .bt_sampler_declaration_stop("'name' must identify a prior in 'prior_list'.")
  }
  prior <- prior_list[[name]]
  if(!is.prior.mixture(prior) || is.prior.spike_and_slab(prior) || length(prior) < 2L){
    .bt_sampler_declaration_stop("Prior '", name, "' must be an ordinary or publication-bias mixture of at least two branches.")
  }

  weights <- unname(attr(prior, "prior_weights", exact = TRUE))
  weights <- weights / max(weights)
  probabilities <- weights / sum(weights)
  if(inherits(prior, "prior.bias_mixture")){
    if(name != "bias"){
      .bt_sampler_declaration_stop("A publication-bias mixture must be named 'bias' in the prior list.")
    }
    omega_spec <- .bt_dnode_omega_spec(prior)
    branches <- lapply(seq_along(prior), function(k){
      component <- prior[[k]]
      info <- omega_spec$branches[[k]]
      family <- if(info$type == "weightfunction") "selection" else info$type
      source <- NA_character_
      coordinates <- character()
      stochastic <- FALSE
      location <- NULL
      omega <- character()
      if(!is.null(info$selection)){
        names <- .weightfunction_component_node_names(
          info$selection, component_id = k, global_cuts = omega_spec$global_cuts, force_one_sided = TRUE
        )
        omega <- paste0(names$omega_target, "[", seq_len(omega_spec$n_bins), "]")
      }
      if(family %in% c("PET", "PEESE")){
        bias_term <- .bt_dnode_prior_mixture_bias_term(name, prior, family)
        if(!identical(unname(bias_term$spec$branch), k)){
          .bt_sampler_declaration_stop("Only one ", family, " style publication bias adjustment is allowed.")
        }
        source <- bias_term$spec$source
        coordinates <- source
        stochastic <- !is.prior.point(component)
        if(!stochastic && !.is_prior_expression(component)){
          location <- component$parameters[["location"]]
        }
      }
      list(prior = component, family = family, source = source, coordinates = coordinates,
           stochastic = stochastic, location = location, omega = omega)
    })
    # The bias emitter uses one fixed categorical node for all bias families.
    return(list(gate = "bias_indicator", probabilities = probabilities, branches = branches))
  }

  spec <- .bt_dnode_prior_mixture(name, prior)$spec
  branches <- lapply(seq_along(prior), function(k){
    component <- spec$components[[k]]
    branch_prior <- prior[[k]]
    family <- .bt_sampler_mixture_family(branch_prior)
    stochastic <- if(family == "unsupported") NA else
      family %in% c("scalar", "vector") ||
      (family == "scalar_array" && length(component$coordinates) == 1L)
    list(prior = branch_prior, family = family, source = component$name,
         coordinates = component$coordinates, stochastic = stochastic,
         location = component$location, omega = character())
  })
  list(gate = spec$indicator, probabilities = probabilities, branches = branches)
}

.bt_sampler_mixture_family <- function(prior){

  if(is.prior.point(prior)) return("point")
  if(is.prior.factor(prior) && (is.prior.treatment(prior) || is.prior.independent(prior))){
    return("scalar_array")
  }
  if(is.prior.vector(prior) && identical(prior$distribution, "mnormal")) return("vector")
  if(is.prior.simple(prior) && !is.prior.factor(prior) && !is.prior.vector(prior)) return("scalar")
  "unsupported"
}


# The scalar stochastic leaves of a prior that may join a factor slice block,
# or the reason it may not. The class guards come first, because a prior keeps
# its class when it also carries random-effect metadata.
.bt_sampler_block_leaves <- function(prior, name){

  reason <- function(text) list(leaves = NULL, reason = text)
  if(is.prior.spike_and_slab(prior)) return(reason("it is a spike-and-slab prior"))
  if(is.prior.mixture(prior))        return(reason("it is a mixture prior"))
  if(is.prior.ordered(prior))        return(reason("it is an ordered prior"))
  if(is.prior.R2(prior)){
    if(!isTRUE(prior$bound)){
      return(reason("its R-squared term has not been bound to the formula"))
    }
    return(reason("its bound R-squared coefficient is deterministic"))
  }
  if(is.prior.R2_budget(prior)){
    return(reason("its realized R-squared budget total is deterministic"))
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
  missing <- setdiff(listed, names(prior_list))
  if(length(missing) > 0L){
    .bt_sampler_declaration_stop(
      "Prior '", missing[1L], "' in 'sampler_blocks' is not in the combined prior list of the model."
    )
  }
  sources <- JAGS_sampler_sources(prior_list, listed, kind = "block")
  for(k in seq_along(sampler_blocks)){
    members <- sources[sources$prior %in% sampler_blocks[[k]], , drop = FALSE]
    first_declined <- match("declined", members$status)
    if(!is.na(first_declined)){
      .bt_sampler_declaration_stop(
        "Prior '", members$prior[first_declined], "' cannot join a block of 'sampler_blocks': ",
        members$reason[first_declined], "."
      )
    }
    members <- members[members$status == "eligible", , drop = FALSE]
    leaves <- members$node
    if(length(leaves) < 2L){
      declined[nrow(declined) + 1L, ] <- list(
        "block", paste(sampler_blocks[[k]], collapse = ", "), NA_integer_, as.integer(k),
        paste(leaves, collapse = ", "), "the block expands to fewer than two sampled nodes"
      )
      next
    }
    block <- data.frame(node = leaves, block = rep(k, length(leaves)), prior = members$prior, stringsAsFactors = FALSE)
    blocks <- rbind(blocks, block)
  }
  blocks$block <- as.integer(blocks$block)

  list(declared = blocks, declined = declined)
}
