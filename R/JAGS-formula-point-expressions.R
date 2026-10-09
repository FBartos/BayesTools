# Scalar expression-point coefficients retain declarations, never caller state.
.bt_formula_numeric_point <- function(prior){

  if(!is.prior.point(prior)) return(NULL)
  location <- prior$parameters[["location"]]
  if(is.numeric(location) && is.null(dim(location)) && length(location) == 1L && is.finite(location)){
    as.numeric(location)
  }else NULL
}

.bt_formula_point_stop <- function(name, reason, detail){

  condition <- structure(list(message = paste0("Expression point '", name,
    "' is unavailable: ", detail, "."), call = NULL, parameter = name,
    reason = reason), class = c("BayesTools_formula_point_unavailable",
      "BayesTools_formula_transform_unavailable", "error", "condition"))
  stop(condition)
}

.bt_formula_point_dimension <- function(prior){

  dimension <- if(.bt_prior_is_factor_family(prior)) .get_prior_factor_levels(prior) else if(is.prior.vector(prior)) prior$parameters[["K"]] else 1L
  if(!is.numeric(dimension) || length(dimension) != 1L || !is.finite(dimension) ||
     dimension < 1L || dimension != floor(dimension)){
    .bt_formula_point_stop("dimension", "malformed_point_state", "the declared coordinate dimension is unavailable")
  }
  as.integer(dimension)
}

.bt_formula_point_finalize <- function(design, prior_list, formula_data,
                                       model_data = NULL,
                                       parameter_names = character()){

  fixed <- intersect(names(design$prior_list), paste0(design$parameter, "_", design$model_terms))
  targets <- fixed[vapply(design$prior_list[fixed], function(prior){
    is.prior.point(prior) && is.expression(prior$parameters[["location"]])
  }, logical(1))]
  if(!length(targets)) return(design)
  if(!is.list(prior_list) || is.null(names(prior_list)) || anyDuplicated(names(prior_list)) ||
     any(!vapply(prior_list, is.prior, logical(1)))){
    .bt_formula_point_stop(design$parameter, "malformed_point_owner", "prior declarations are malformed")
  }
  formula_data <- .bt_formula_expression_data_list(formula_data, "Expression point data")
  model_data <- .bt_formula_expression_data_list(model_data, "Expression point data")
  for(source in list(formula_data, model_data)){
    if(length(source) && (is.null(names(source)) || anyNA(names(source)) ||
       any(!nzchar(names(source))) || anyDuplicated(names(source)))){
      .bt_formula_point_stop(design$parameter, "malformed_point_data", "owned data must have unique nonempty names")
    }
  }
  data <- formula_data
  data[setdiff(names(model_data), names(data))] <- model_data[setdiff(names(model_data), names(data))]
  parameters <- unique(c(names(prior_list), .bt_parameter_coordinates_base(parameter_names)))
  points <- list()
  retained <- character()
  needed_data <- character()
  pending <- targets
  while(length(pending)){
    name <- pending[[1L]]
    pending <- pending[-1L]
    if(name %in% retained) next
    retained <- c(retained, name)
    prior <- prior_list[[name]]
    if(is.null(prior)) next
    if(!is.prior.point(prior) || !is.expression(prior$parameters[["location"]])) next
    location <- prior$parameters[["location"]]
    if(length(location) != 1L){
      .bt_formula_point_stop(name, "malformed_point_declaration", "the location must contain one expression")
    }
    original <- location[[1L]]
    label <- .bt_formula_expression_label(original)
    dependencies <- setdiff(.bt_formula_expression_symbols(original), "i")
    overlap <- intersect(dependencies, intersect(names(data), parameters))
    if(length(overlap)) .bt_formula_point_stop(name, "contradictory_point_owner",
      "a parent is declared as both data and a parameter")
    spec <- tryCatch(.bt_formula_expression_specs(list(label), names(data), parameters,
      allow_unresolved = TRUE)[[1L]], error = identity)
    reason <- detail <- NULL
    if(inherits(spec, "error")){
      reason <- "unsupported_point_expression"
      detail <- conditionMessage(spec)
      spec <- NULL
    }else if("i" %in% .bt_formula_expression_symbols(spec$parsed)){
      reason <- "unsupported_point_row_index"
      detail <- "loop index 'i' has no scalar coefficient replay"
    }else if(length(spec$unresolved_dependencies)){
      reason <- "unresolved_point_parent"
      detail <- paste0("undeclared parent ", paste(spec$unresolved_dependencies, collapse = ", "))
    }
    points[[name]] <- list(location = original, spec = spec, reason = reason, detail = detail,
      parameter_dependencies = intersect(dependencies, parameters))
    needed_data <- union(needed_data, intersect(dependencies, names(data)))
    pending <- union(pending, intersect(dependencies, parameters))
  }
  data <- data[needed_data]
  for(name in intersect(needed_data, intersect(names(formula_data), names(model_data)))){
    if(!identical(formula_data[[name]], model_data[[name]])){
      .bt_formula_point_stop(name, "contradictory_point_data", "formula and model data disagree")
    }
  }
  for(name in names(data)){
    value <- data[[name]]
    if((!is.numeric(value) && !is.logical(value)) || anyNA(value) || any(!is.finite(value))){
      .bt_formula_point_stop(name, "malformed_point_data", "owned data must be finite numeric or logical values")
    }
  }
  design$point_expression_owner <- list(version = 1L, targets = targets,
    points = points, prior_list = prior_list[intersect(retained, names(prior_list))],
    data = data, extra_parents = setdiff(retained, names(prior_list)))
  design
}

.bt_formula_point_owner <- function(design, name){

  owner <- design$point_expression_owner
  if(is.null(owner)){
    .bt_stop_refit_required("Expression-point formula declarations are missing. Refit the model with this version of BayesTools.")
  }
  if(!is.list(owner) || !identical(owner$version, 1L) ||
     !identical(names(owner), c("version", "targets", "points", "prior_list", "data", "extra_parents")) ||
     !is.character(owner$targets) || !is.list(owner$points) ||
     !is.list(owner$prior_list) || !is.list(owner$data) ||
     !is.character(owner$extra_parents) ||
     !name %in% names(owner$points)){
    .bt_formula_point_stop(name, "malformed_point_owner", "the retained expression owner is malformed")
  }
  for(point_name in names(owner$points)){
    prior <- owner$prior_list[[point_name]]
    point <- owner$points[[point_name]]
    if(!is.prior.point(prior) || !is.expression(prior$parameters[["location"]]) ||
       length(prior$parameters[["location"]]) != 1L ||
       !identical(prior$parameters[["location"]][[1L]], point$location)){
      .bt_formula_point_stop(point_name, "contradictory_point_owner", "the retained location disagrees with its declaration")
    }
  }
  .bt_validate_once("formula_point_owner", owner, function(){
    for(point_name in names(owner$points)){
      point <- owner$points[[point_name]]
      if(!is.list(point) || !is.character(point$parameter_dependencies)){
        .bt_formula_point_stop(point_name, "malformed_point_owner", "the replay declaration is malformed")
      }
      if(!is.null(point$spec)){
        validated <- tryCatch(.bt_parse_formula_expression(point$location), error = identity)
        if(inherits(validated, "error") || !identical(validated, point$spec$parsed) ||
           !identical(.bt_formula_expression_symbols(validated),
             .bt_formula_expression_symbols(point$location))){
          .bt_formula_point_stop(point_name, "malformed_point_owner", "the replay AST disagrees with its location declaration")
        }
        dependencies <- setdiff(.bt_formula_expression_symbols(validated), "i")
        if(!identical(dependencies, point$spec$dependencies) ||
           !setequal(point$parameter_dependencies, point$spec$parameter_dependencies) ||
           any(!point$spec$data_dependencies %in% names(owner$data))){
          .bt_formula_point_stop(point_name, "malformed_point_owner", "the replay dependency declaration is malformed")
        }
      }else if(!identical(point$reason, "unsupported_point_expression")){
        .bt_formula_point_stop(point_name, "malformed_point_owner", "an unavailable replay declaration lacks its syntax reason")
      }
    }
    for(data_name in names(owner$data)){
      value <- owner$data[[data_name]]
      if((!is.numeric(value) && !is.logical(value)) || anyNA(value) || any(!is.finite(value))){
        .bt_formula_point_stop(data_name, "malformed_point_data", "owned replay data must remain finite")
      }
    }
  })
  owner
}

.bt_formula_point_values <- function(name, prior, design = NULL, samples = NULL,
                                     parameters = NULL, n_draws = 1L,
                                     n_values = 1L, recompute = FALSE){

  literal <- .bt_formula_numeric_point(prior)
  if(!is.null(literal)) return(matrix(literal, n_draws, n_values))
  if(!is.expression(prior$parameters[["location"]])){
    .bt_formula_point_stop(name, "malformed_point_declaration", "the location is neither a finite scalar nor a single expression")
  }
  owner <- .bt_formula_point_owner(design, name)
  if(length(prior$parameters[["location"]]) != 1L ||
     !identical(prior$parameters[["location"]][[1L]], owner$points[[name]]$location)){
    .bt_formula_point_stop(name, "contradictory_point_owner", "the coefficient declaration disagrees with its retained owner")
  }
  context <- paste0("Expression point '", name, "'")
  monitor <- function(parameter){
    sample_names <- .bt_formula_expression_sample_names(samples)
    if(!parameter %in% sample_names &&
       !any(startsWith(sample_names, paste0(parameter, "["))) &&
       is.null(parameters[[parameter]])) return(NULL)
    draws <- tryCatch(.bt_formula_expression_parameter_draws(samples, parameters, parameter, n_draws, context), error = identity)
    if(inherits(draws, "error")) .bt_formula_point_stop(parameter, "malformed_point_state", conditionMessage(draws))
    indexed <- inherits(draws, "BayesTools_formula_expression_indexed_draws")
    if(indexed){
      declaration <- owner$prior_list[[parameter]]
      if(.bt_prior_is_factor_family(declaration) || is.prior.vector(declaration)){
        draws$size <- .bt_formula_point_dimension(declaration)
        if(any(draws$indices > draws$size)) .bt_formula_point_stop(parameter,
          "malformed_point_state", "supplied indices exceed the declared coordinate dimension")
      }
      if(any(!is.finite(draws$values))) .bt_formula_point_stop(parameter,
        "nonfinite_point_state", "supplied coordinates must be finite")
    }
    values <- lapply(seq_len(n_draws), function(draw){
      .bt_formula_expression_parameter_value(draws, draw)
    })
    declaration <- owner$prior_list[[parameter]]
    if(!indexed && (.bt_prior_is_factor_family(declaration) || is.prior.vector(declaration)) &&
       any(lengths(values) != .bt_formula_point_dimension(declaration))){
      .bt_formula_point_stop(parameter, "malformed_point_state", "supplied values disagree with the declared coordinate dimension")
    }
    if(!indexed && any(vapply(values, function(value) !is.numeric(value) || any(!is.finite(value)), logical(1)))){
      .bt_formula_point_stop(parameter, "nonfinite_point_state", "supplied values must be finite numeric values")
    }
    values
  }
  resolve <- function(parameter, path = character(), prefer_replay = FALSE){
    if(parameter %in% path){
      .bt_formula_point_stop(parameter, "cyclic_point_parent", "declared parents contain a cycle")
    }
    declaration <- owner$prior_list[[parameter]]
    fixed <- .bt_formula_numeric_point(declaration)
    if(!is.null(fixed)) return(rep(list(rep(fixed, .bt_formula_point_dimension(declaration))), n_draws))
    point <- owner$points[[parameter]]
    if(!prefer_replay || is.null(point)){
      supplied <- monitor(parameter)
      if(!is.null(supplied)) return(supplied)
    }
    if(is.null(point)){
      .bt_formula_point_stop(parameter, "missing_point_parent", "its declared source has no supplied coordinate")
    }
    replay <- function(){
      if(!is.null(point$reason)) .bt_formula_point_stop(parameter, point$reason, point$detail)
      spec <- point$spec
      parents <- lapply(spec$parameter_dependencies, function(parent){
        resolve(parent, c(path, parameter), prefer_replay)
      })
      names(parents) <- spec$parameter_dependencies
      lapply(seq_len(n_draws), function(draw){
        values <- lapply(parents, `[[`, draw)
        value <- tryCatch(.bt_formula_expression_eval(spec, owner$data, 1L,
          values, context, mode = "point"), error = identity)
        if(inherits(value, "error")){
          if(inherits(value, "BayesTools_formula_point_unavailable")) stop(value)
          .bt_formula_point_stop(parameter, "point_expression_evaluation", conditionMessage(value))
        }
        rep(value, .bt_formula_point_dimension(declaration))
      })
    }
    result <- tryCatch(replay(), BayesTools_formula_point_unavailable = identity)
    if(inherits(result, "BayesTools_formula_point_unavailable")){
      if(prefer_replay && length(path) > 0L && identical(result$reason, "missing_point_parent")){
        supplied <- monitor(parameter)
        if(!is.null(supplied)) return(supplied)
      }
      stop(result)
    }
    result
  }
  values <- resolve(name, prefer_replay = recompute)
  if(any(lengths(values) != n_values) && any(lengths(values) != 1L)){
    .bt_formula_point_stop(name, "nonscalar_point_expression", "the coefficient result must be scalar")
  }
  output <- do.call(rbind, lapply(values, function(value){
    if(length(value) == 1L) rep(value, n_values) else value
  }))
  matrix(output, n_draws, n_values)
}

.bt_formula_point_index <- function(x, ...){

  indices <- list(...)
  if(!is.null(dim(x)) && length(dim(x)) != length(indices)){
    .bt_formula_point_stop("index", "point_data_shape", "one index is required for every owned array dimension")
  }
  if(length(indices) != 1L){
    return(tryCatch(.bt_formula_expression_index(x, ...), error = function(e){
      .bt_formula_point_stop("index", "invalid_point_index", conditionMessage(e))
    }))
  }
  index <- indices[[1L]]
  if(!is.numeric(index) || !is.null(dim(index)) || !length(index) ||
     any(!is.finite(index)) || any(index <= 0 | index != floor(index) | index > length(x))){
    .bt_formula_point_stop("index", "invalid_point_index", "indices must be finite positive integers within their declared dimensions")
  }
  value <- x[index]
  if(anyNA(value)) .bt_formula_point_stop("index", "missing_point_coordinate", "an indexed parameter coordinate is unavailable")
  value
}

.bt_dnode_point_expression <- function(name, prior, design){

  owner <- .bt_formula_point_owner(design, name)
  point <- owner$points[[name]]
  n_values <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else 1L
  coordinates <- if(n_values > 1L) paste0(name, "[", seq_len(n_values), "]") else name
  .bt_deterministic_node("point_expression", name, coordinates,
    dependencies = point$parameter_dependencies,
    parameter = design$parameter,
    spec = list(prior = prior, design = list(point_expression_owner = owner), n_values = n_values,
      syntax = if(is.prior.factor(prior)) .JAGS_prior.factor(prior, name) else .JAGS_prior.simple(prior, name)))
}

.bt_dnode_point_expression_evaluate <- function(node, lookup){

  .bt_formula_point_values(node$node, node$spec$prior,
    design = node$spec$design, samples = lookup$draws, n_draws = lookup$n,
    n_values = node$spec$n_values, recompute = TRUE)
}
