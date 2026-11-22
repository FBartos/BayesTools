# Package-supplied declarations add observed zero markers or terminal
# deterministic tags, with optional numeric data for their arguments.
.bt_sampler_finite_prefix <- "BayesTools_finite_"
.bt_sampler_finite_readout_prefix <- "BayesTools_finite_readout_"
.bt_sampler_node_pattern <- "[A-Za-z][A-Za-z0-9_.]*"

.bt_sampler_check_plan_arguments <- function(sampler_joint, sampler_probabilities){

  if(!is.null(sampler_joint) && !is.function(sampler_joint)){
    .bt_sampler_declaration_stop("The 'sampler_joint' argument must be NULL or a function.")
  }
  if(is.null(sampler_probabilities)){
    return(character())
  }
  nodes <- .bt_sampler_node_references(sampler_probabilities, ranges = FALSE)
  if(is.null(nodes) || anyDuplicated(nodes)){
    .bt_sampler_declaration_stop(
      "The 'sampler_probabilities' argument must contain distinct JAGS node names or scalar indexed references."
    )
  }
  nodes
}

.bt_sampler_node_references <- function(value, ranges = TRUE){

  if(!is.character(value) || !is.null(dim(value)) || anyNA(value)){
    return(NULL)
  }
  # Spaces around the brackets and separators are valid. Spaces within a
  # symbol or index are invalid rather than silently joined together.
  index <- if(ranges) "[1-9][0-9]*(?:[ \\t]*:[ \\t]*[1-9][0-9]*)?" else "[1-9][0-9]*"
  pattern <- paste0("\\A[ \\t]*", .bt_sampler_node_pattern,
    "(?:[ \\t]*\\[[ \\t]*", index, "(?:[ \\t]*,[ \\t]*", index,
    ")*[ \\t]*\\])?[ \\t]*\\z")
  if(!all(grepl(pattern, value, perl = TRUE))){
    return(NULL)
  }
  unname(gsub("[ \\t]", "", value, perl = TRUE))
}

.bt_sampler_joint_fields <- function(object, required, optional = character()){

  is.list(object) && !is.object(object) &&
    length(names(object)) == length(object) && !anyDuplicated(names(object)) &&
    all(required %in% names(object)) && all(names(object) %in% c(required, optional))
}

.bt_sampler_joint_numeric <- function(value){

  typeof(value) %in% c("double", "integer") && !is.object(value) &&
    length(value) > 0L && all(is.finite(value))
}

.bt_sampler_joint_names <- function(value){

  (is.character(value) || is.null(value)) && !anyNA(value) &&
    all(grepl(paste0("\\A", .bt_sampler_node_pattern, "\\z"), value, perl = TRUE))
}

.bt_sampler_joint_markers <- function(joint, existing, jags_modules, retired = character(),
                                      model_nodes = existing, ordinary_blocks = FALSE){

  fail <- function(...){
    .bt_sampler_declaration_stop("The 'sampler_joint' function returned an invalid declaration: ", ..., ".")
  }
  if(!.bt_sampler_joint_fields(joint, "markers", "data") ||
     !is.list(joint[["markers"]]) || is.object(joint[["markers"]]) || !length(joint[["markers"]])){
    fail("expected a nonempty 'markers' list and optional 'data'")
  }
  payload <- joint[["data"]]
  if(!is.null(payload) && (!is.list(payload) || is.object(payload) ||
     length(names(payload)) != length(payload) || !.bt_sampler_joint_names(names(payload)) ||
     anyDuplicated(names(payload)) || !all(vapply(payload, .bt_sampler_joint_numeric, logical(1))))){
    fail("'data' must be a uniquely named list of finite numeric constants")
  }
  payload <- lapply(payload, function(value){
    if(is.null(dim(value))) as.vector(value) else array(as.vector(value), dim(value))
  })

  reserved <- c(.bt_sampler_gate_prefix, .bt_sampler_block_prefix, .bt_sampler_finite_prefix)
  tag_prefix <- paste0(tolower(jags_modules), "_", recycle0 = TRUE)
  distribution_prefix <- paste0("d", tag_prefix, recycle0 = TRUE)
  marked <- names(payload)
  check_nested <- function(nodes){
    conflict <- intersect(nodes, retired)
    if(length(conflict)){
      fail("node '", conflict[[1L]],
        "' was replaced by nested weight sharing; use 'sampler_nested_weights' = FALSE to access this node")
    }
  }
  check_fresh <- function(names){
    check_nested(names)
    collision <- names %in% existing |
      vapply(names, function(name) any(startsWith(name, reserved)), logical(1))
    if(any(collision)){
      fail("name '", names[which(collision)[1L]], "' is already used or reserved")
    }
  }
  check_fresh(marked)
  known <- c(model_nodes, marked)
  markers <- lapply(unname(joint[["markers"]]), function(marker){
    type <- intersect(names(marker), c("distribution", "tag"))
    if(!.bt_sampler_joint_fields(marker, c("name", "arguments"), c("distribution", "tag")) ||
       length(type) != 1L ||
       length(marker$name) != 1L || !.bt_sampler_joint_names(marker$name) ||
       length(marker[[type]]) != 1L || !.bt_sampler_joint_names(marker[[type]]) ||
       !is.list(marker$arguments) || is.object(marker$arguments)){
      fail("each marker must have a JAGS 'name', exactly one of 'distribution' and 'tag', and an ordered 'arguments' list")
    }
    symbol <- as.vector(marker[[type]])
    if(type == "tag"){
      if(!any(startsWith(symbol, tag_prefix))){
        fail("tag '", symbol, "' must use the lowercase prefix of a provider listed in 'jags_modules'")
      }
    }else if(symbol == "dbt_gate"){
      fail("distribution 'dbt_gate' is reserved for gates verified from the model. Use 'sampler_gates' to enable those gates")
    }else if(symbol == "dbt_block" && ordinary_blocks){
      fail("callback 'dbt_block' markers cannot be combined with blocks from 'sampler_blocks'. Declare BayesTools blocks through one argument")
    }else if(symbol != "dbt_block" && !any(startsWith(symbol, distribution_prefix))){
      fail("distribution '", symbol,
        "' must be a BayesTools marker or use the prefix of a provider listed in 'jags_modules'")
    }
    check_fresh(marker$name)
    if(marker$name %in% marked){
      fail("name '", marker$name, "' is used more than once")
    }
    marked <<- c(marked, marker$name)
    arguments <- lapply(unname(marker$arguments), function(value){
      references <- .bt_sampler_node_references(value)
      if(length(references)){
        bases <- .bt_jags_parameter_base(references)
        check_nested(bases)
        unknown <- setdiff(bases, known)
        if(length(unknown)){
          fail("unknown node '", unknown[[1L]], "' in marker '", marker$name, "'")
        }
        return(references)
      }
      if(.bt_sampler_joint_numeric(value) && is.null(dim(value))){
        return(as.vector(value))
      }
      fail("each argument must contain JAGS node references or finite numeric literals")
    })
    stats::setNames(list(as.vector(marker$name), symbol, arguments), c("name", type, "arguments"))
  })
  arguments <- unlist(lapply(markers, "[[", "arguments"), recursive = FALSE, use.names = FALSE)
  references <- unlist(Filter(is.character, arguments), use.names = FALSE)
  unused <- setdiff(names(payload), .bt_jags_parameter_base(references))
  if(length(unused)){
    fail("data '", unused[[1L]], "' is not referenced by any marker")
  }
  lines <- vapply(markers, function(marker){
    args <- vapply(marker$arguments, function(value){
      if(is.numeric(value)) value <- trimws(formatC(value, digits = 17L, format = "g", decimal.mark = "."))
      if(length(value) == 1L) value else paste0("c(", paste(value, collapse = ", "), ")")
    }, character(1))
    if(!is.null(marker$tag)){
      paste0(marker$name, " <- ", marker$tag, "(", paste(args, collapse = ", "), ")")
    }else{
      paste0(marker$name, " ~ ", marker$distribution, "(", paste(args, collapse = ", "), ")")
    }
  }, character(1))
  observed <- Filter(function(marker) !is.null(marker$distribution), markers)
  observed <- vapply(observed, "[[", character(1), "name")
  list(lines = lines, data = c(payload, stats::setNames(rep(list(0), length(observed)), observed)),
    declaration = list(markers = markers, data = names(payload)))
}

.bt_sampler_fit_gate_claims <- function(gates, prior_list, fit){

  workers <- fit[["samplers_by_worker"]]
  if(is.null(workers)){
    return(.bt_sampler_gate_claims(gates, prior_list, fit$samplers))
  }
  claims <- lapply(workers, function(worker){
    .bt_sampler_gate_claims(gates, prior_list, worker$samplers)
  })
  result <- claims[[1L]]
  mixed <- rep(FALSE, nrow(result))
  for(field in c("reported_node", "sampler", "assignment", "status")){
    values <- do.call(cbind, lapply(claims, "[[", field))
    differs <- apply(values, 1L, function(value) length(unique(value)) > 1L)
    if(field %in% c("sampler", "status")){
      mixed <- mixed | differs
    }
    if(field != "status"){
      result[[field]][differs] <- NA
    }
  }
  result$status[mixed] <- "mixed"
  result
}

.bt_sampler_apply_plan <- function(legacy, model_names, prior_list, formula_design, data,
                                   sampler_joint, sampler_probabilities, jags_modules = NULL){

  joint <- if(!is.null(sampler_joint)){
    sampler_joint(prior_list = prior_list, formula_design = formula_design, data = data)
  }
  if(is.null(joint) && !length(sampler_probabilities)){
    return(legacy)
  }
  existing <- unique(c(.bt_sampler_tokens(legacy$model), names(legacy$data),
    names(prior_list), trimws(.bt_jags_parameter_base(trimws(model_names)))))
  graph <- .bt_jags_syntax_graph(legacy$model)
  model_nodes <- unique(c(graph$stochastic, names(graph$deterministic), graph$data, names(legacy$data)))
  if(!is.null(joint)){
    bias <- prior_list[["bias"]]
    nested <- .bt_nested_weights_plan(bias)
    retired <- if(!is.null(nested)) .bt_nested_weights_retired(bias, nested)
    package_markers <- .bt_sampler_joint_markers(joint, existing, jags_modules, retired,
      model_nodes = model_nodes, ordinary_blocks = NROW(legacy$declarations$blocks) > 0L)
    legacy$model <- .bt_sampler_insert_markers(legacy$model, package_markers$lines)
    if(length(package_markers$data)){
      legacy$data <- c(legacy$data, package_markers$data)
    }
    legacy$declarations$joint <- package_markers$declaration
    legacy$uses_module <- TRUE
  }

  if(length(sampler_probabilities)){
    reserved <- existing[startsWith(existing, .bt_sampler_finite_prefix)]
    if(length(reserved)){
      .bt_sampler_declaration_stop("The model uses reserved BayesTools sampler node name '", reserved[[1L]], "'.")
    }
    unknown <- setdiff(.bt_jags_parameter_base(sampler_probabilities), model_nodes)
    if(length(unknown)){
      .bt_sampler_declaration_stop("Node '", unknown[[1L]], "' in 'sampler_probabilities' is unavailable in the model.")
    }
    id <- seq_along(sampler_probabilities)
    plan <- data.frame(
      node = sampler_probabilities,
      marker = paste0(.bt_sampler_finite_prefix, id),
      readout = paste0(.bt_sampler_finite_readout_prefix, id),
      stringsAsFactors = FALSE
    )
    legacy$model <- .bt_sampler_insert_markers(legacy$model,
      paste0(plan$marker, " <- bt_finite_tag(", plan$node, ")"))
    legacy$declarations$probabilities <- plan
    legacy$uses_module <- TRUE
  }

  legacy
}
