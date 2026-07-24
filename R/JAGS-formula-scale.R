# Helper: Parse a term name into its component variable names
# e.g., "mu_x1__xXx__x2" with prefix "mu" -> c("x1", "x2")
# e.g., "mu_intercept" -> character(0) (intercept has no components)
# e.g., "mu_x1" -> c("x1")
.parse_term_components <- function(term_name, prefix) {
  # Remove prefix
  term_part <- sub(paste0("^", prefix, "_"), "", term_name)

  # Check if it's the intercept
  if (term_part == "intercept") {
    return(character(0))
  }

  # Indexed factor terms attach the index to the full term name. Strip the
  # trailing index so scaled continuous components inside interactions are
  # still detected correctly (e.g., alloc__xXx__year[1] -> alloc, year).
  term_part <- sub("\\[[^]]+\\]$", "", term_part)

  # Split by interaction separator
  components <- strsplit(term_part, "__xXx__", fixed = TRUE)[[1]]
  return(components)
}


# Helper: Parse a term into scaled components and an unscaled identity string
# used to determine which coefficients can contribute to each other during
# de-standardization.
.parse_unscale_term_structure <- function(term_name, prefix, scaled_vars) {

  term_part <- sub(paste0("^", prefix, "_"), "", term_name)

  if (term_part == "intercept") {
    return(list(
      components  = character(0),
      scaled      = character(0),
      unscaled_id = ""
    ))
  }

  term_index <- ""
  if (grepl("\\[[^]]+\\]$", term_part)) {
    term_index <- sub("^.*(\\[[^]]+\\])$", "\\1", term_part)
    term_core  <- sub("\\[[^]]+\\]$", "", term_part)
  } else {
    term_core <- term_part
  }

  components <- strsplit(term_core, "__xXx__", fixed = TRUE)[[1]]
  is_scaled  <- components %in% scaled_vars

  unscaled_components <- components[!is_scaled]
  unscaled_id <- paste(unscaled_components, collapse = "__xXx__")

  if (nzchar(term_index) && length(unscaled_components) > 0) {
    unscaled_id <- paste0(unscaled_id, term_index)
  }

  return(list(
    components  = components,
    scaled      = components[is_scaled],
    unscaled_id = unscaled_id
  ))
}


# Helper: Check if set A is a subset of set B (including equality)
.is_subset <- function(A, B) {

  length(A) == 0 || all(A %in% B)
}


# Helper: Compare lower-order unscaled term identities. Two-level factor main
# effects are stored without [1], while their interactions are indexed.
.unscale_ids_match <- function(target_id, source_id) {

  identical(target_id, source_id) ||
    (nzchar(target_id) && !grepl("\\[[^]]+\\]$", target_id) && identical(paste0(target_id, "[1]"), source_id))
}


# Helper: Refuse incomplete coefficient transformations. Centered interactions
# induce lower-order terms, which cannot be represented by a square transform
# when the corresponding posterior columns are absent.
.check_unscale_term_closure <- function(term_names, term_scaled, term_unscaled,
                                        formula_scale, prefix) {

  for(source_name in term_names){
    source_scaled <- unique(term_scaled[[source_name]])
    if(length(source_scaled) == 0L){
      next
    }

    source_unscaled <- term_unscaled[[source_name]]
    has_returned_lower_order <- any(vapply(
      setdiff(term_names, source_name),
      function(target_name){
        target_scaled <- term_scaled[[target_name]]
        .unscale_ids_match(term_unscaled[[target_name]], source_unscaled) &&
          length(target_scaled) < length(source_scaled) &&
          .is_subset(target_scaled, source_scaled)
      },
      logical(1)
    ))
    if(!has_returned_lower_order){
      # Transforming an interaction coefficient by itself is complete: the
      # omitted lower-order coefficients are outside the requested output.
      next
    }

    missing_targets <- character()
    for(subset_size in 0:(length(source_scaled) - 1L)){
      target_sets <- if(subset_size == 0L){
        list(character())
      }else{
        utils::combn(source_scaled, subset_size, simplify = FALSE)
      }

      for(target_scaled in target_sets){
        extra_scaled <- setdiff(source_scaled, target_scaled)
        extra_params <- paste0(prefix, "_", extra_scaled)
        mean_product <- prod(vapply(
          extra_params,
          function(parameter) formula_scale[[parameter]][["mean"]],
          numeric(1)
        ))
        if(mean_product == 0){
          next
        }

        target_exists <- any(vapply(
          term_names,
          function(target_name){
            target_components <- term_scaled[[target_name]]
            .unscale_ids_match(term_unscaled[[target_name]], source_unscaled) &&
              length(target_components) == length(target_scaled) &&
              all(target_components %in% target_scaled)
          },
          logical(1)
        ))
        if(!target_exists){
          missing_targets <- c(
            missing_targets,
            if(length(target_scaled) == 0L){
              "<no scaled components>"
            }else{
              paste(target_scaled, collapse = ":")
            }
          )
        }
      }
    }

    if(length(missing_targets) > 0L){
      stop(
        "Cannot unscale posterior term '", source_name,
        "' because centering induces missing lower-order coefficient(s) for ",
        "scaled component set(s): ",
        paste0(unique(missing_targets), collapse = ", "),
        ". Include the corresponding lower-order posterior columns.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


# Helper: Validate the nested formula_scale structure used for unscaling.
.check_formula_scale_info <- function(formula_scale, name = "formula_scale") {

  check_list(formula_scale, name)

  if(is.null(names(formula_scale)) || anyNA(names(formula_scale)) || any(names(formula_scale) == ""))
    stop(paste0("The '", name, "' argument must be a named nested list keyed by parameter name."), call. = FALSE)
  if(anyDuplicated(names(formula_scale)))
    stop(paste0("The '", name, "' argument must not contain duplicate parameter names."), call. = FALSE)

  for(param_name in names(formula_scale)){
    param_scale <- formula_scale[[param_name]]
    check_list(param_scale, paste0(name, "[['", param_name, "]]"))

    if(length(param_scale) == 0)
      next

    if(is.null(names(param_scale)) || anyNA(names(param_scale)) || any(names(param_scale) == ""))
      stop(paste0("The '", name, "[['", param_name, "]]" ,"' entry must be a named list keyed by parameter term."), call. = FALSE)
    if(anyDuplicated(names(param_scale)))
      stop(paste0("The '", name, "[['", param_name, "']]' entry must not contain duplicate term names."), call. = FALSE)

    for(term_name in names(param_scale)){
      term_scale <- param_scale[[term_name]]

      check_list(
        term_scale,
        paste0(name, "[['", param_name, "]][['", term_name, "]]"),
        check_names = c("mean", "sd"),
        all_objects = TRUE,
        allow_other = TRUE
      )
      check_real(
        term_scale[["mean"]],
        paste0(name, "[['", param_name, "]][['", term_name, "]][['mean'] ]"),
        allow_NA = FALSE
      )
      if(!is.finite(term_scale[["mean"]])){
        stop(
          "The '", name, "[['", param_name, "]][['", term_name,
          "]][['mean']]' entry must be finite.",
          call. = FALSE
        )
      }
      check_real(
        term_scale[["sd"]],
        paste0(name, "[['", param_name, "]][['", term_name, "]][['sd'] ]"),
        lower = 0,
        allow_bound = FALSE,
        allow_NA = FALSE
      )
      if(!is.finite(term_scale[["sd"]])){
        stop(
          "The '", name, "[['", param_name, "]][['", term_name,
          "]][['sd']]' entry must be finite.",
          call. = FALSE
        )
      }
    }
  }

  invisible(NULL)
}


# Helper: warn when formula_scale entries do not map to posterior terms.
.warn_unused_formula_scale_terms <- function(term_names, formula_scale, prefix) {

  scaled_terms <- names(formula_scale)
  if(length(scaled_terms) == 0)
    return(invisible(NULL))

  scaled_vars <- sub(paste0("^", prefix, "_"), "", scaled_terms)
  term_components <- unique(unlist(lapply(term_names, .parse_term_components, prefix = prefix), use.names = FALSE))
  random_scaled_vars <- .formula_scale_random_scaled_vars(
    formula_scale,
    scaled_vars = scaled_vars
  )
  unused_terms <- scaled_terms[
    !scaled_vars %in% term_components &
      !scaled_vars %in% random_scaled_vars
  ]

  if(length(unused_terms) == 0)
    return(invisible(NULL))

  if(length(unused_terms) == length(scaled_terms)){
    warning(
      paste0(
        "Ignoring all formula_scale[['", prefix, "]] entries because none match posterior terms: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Samples for this parameter prefix are returned unchanged."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }else{
    warning(
      paste0(
        "Ignoring unused formula_scale[['", prefix, "]] entries: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Matched entries are still applied."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  invisible(NULL)
}

.formula_scale_random_scaled_vars <- function(formula_scale, scaled_vars = NULL){

  if(is.null(scaled_vars)){
    prefix <- attr(formula_scale, "parameter")
    if(!is.null(prefix) && length(prefix) == 1L && !is.na(prefix) && nzchar(prefix)){
      scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))
    }else{
      scaled_vars <- names(formula_scale)
    }
  }

  random_terms <- character()

  metadata <- attr(formula_scale, "random_effect_terms")
  if(!is.null(metadata)){
    metadata <- unname(as.character(metadata))
    metadata <- metadata[is.na(metadata) | metadata != "sd"]
    random_terms <- c(random_terms, metadata)
  }

  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves")
  if(!is.null(sd_leaves) && length(sd_leaves) > 0L){
    for(leaves in sd_leaves){
      random_terms <- c(
        random_terms,
        .formula_scale_random_leaf_terms(
          leaves = leaves,
          scaled_vars = scaled_vars
        )
      )
    }
  }

  random_terms <- random_terms[!is.na(random_terms) & nzchar(random_terms)]
  random_terms <- setdiff(random_terms, "intercept")
  if(length(random_terms) == 0L){
    return(character())
  }

  random_components <- unique(unlist(lapply(random_terms, .bt_random_effect_term_components), use.names = FALSE))
  random_components <- setdiff(random_components, "intercept")
  if(length(scaled_vars) > 0L){
    random_components <- random_components[random_components %in% scaled_vars]
  }
  unique(random_components)
}

.formula_scale_random_leaf_terms <- function(leaves, scaled_vars){

  if(.random_sd_structured_leaves(leaves)){
    return(character())
  }

  leaf_terms <- character()
  if(!is.null(leaves$leaf_terms)){
    leaf_terms <- unname(as.character(leaves$leaf_terms))
  }

  homogeneous_sd <- length(leaf_terms) > 0L &&
    all(!is.na(leaf_terms) & leaf_terms == "sd")
  if(!homogeneous_sd){
    if(!is.null(leaves$leaf_terms_by_column)){
      leaf_terms <- c(leaf_terms, unname(as.character(leaves$leaf_terms_by_column)))
    }
    return(leaf_terms)
  }

  if(is.null(leaves$column_names)){
    return(character())
  }

  column_terms <- vapply(
    leaves$column_names,
    .random_sd_term_from_column_name,
    character(1)
  )
  scaled_columns <- vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  )

  column_terms[scaled_columns]
}


# Helper: Build the transformation matrix for unscaling coefficients
#
# For each target term T and source term S, computes the coefficient M[T,S] such that:
#   coef_orig[T] = sum over S of M[T,S] * coef_z[S]
#
# The formula is based on expanding products of (x_i - mu_i)/sigma_i terms.
# For S to contribute to T:
#   1. T_unscaled == S_unscaled (unscaled components must match exactly)
#   2. T_scaled is a subset of S_scaled
#
# The contribution is: (-1)^|extra| * prod(mu_extra) / prod(sigma_S_scaled)
# where extra = S_scaled \ T_scaled
#
# @param term_names Character vector of all term names in the posterior
# @param formula_scale Named list with scaling info (mean, sd) for scaled predictors
# @param prefix The parameter prefix (e.g., "mu")
# @param require_closure Whether missing induced coefficient terms should fail.
#   Random-effect SD transforms operate in covariance space and set this to FALSE.
# @return A square transformation matrix
.build_unscale_matrix <- function(term_names, formula_scale, prefix,
                                  require_closure = TRUE) {

  n_terms <- length(term_names)
  M <- diag(n_terms)  # Start with identity matrix
  rownames(M) <- colnames(M) <- term_names

  # Extract the variable names that are scaled (without prefix)
  scaled_vars <- sub(paste0("^", prefix, "_"), "", names(formula_scale))

  # Parse all terms into their scaled components and unscaled identity.
  term_structure <- lapply(
    term_names,
    .parse_unscale_term_structure,
    prefix = prefix,
    scaled_vars = scaled_vars
  )
  names(term_structure) <- term_names

  term_components <- lapply(term_structure, `[[`, "components")
  term_scaled <- lapply(term_structure, `[[`, "scaled")
  term_unscaled <- vapply(term_structure, `[[`, character(1), "unscaled_id")

  if(isTRUE(require_closure)){
    .check_unscale_term_closure(
      term_names = term_names,
      term_scaled = term_scaled,
      term_unscaled = term_unscaled,
      formula_scale = formula_scale,
      prefix = prefix
    )
  }

  # Warn about high-order interactions
  max_order <- max(sapply(term_components, length))
  if (max_order >= 5) {
    warning("Model contains ", max_order, "-way or higher interactions. ",
            "Unscaling transformation may be computationally intensive.",
            immediate. = TRUE)
  }

  # Build the transformation matrix
  for (t_idx in seq_along(term_names)) {
    T_name <- term_names[t_idx]
    T_scaled <- term_scaled[[T_name]]
    T_unscaled <- term_unscaled[[T_name]]

    for (s_idx in seq_along(term_names)) {
      S_name <- term_names[s_idx]
      S_scaled <- term_scaled[[S_name]]
      S_unscaled <- term_unscaled[[S_name]]

      # Check contribution conditions
      # 1. Unscaled parts must match exactly
      if (!.unscale_ids_match(T_unscaled, S_unscaled)) next

      # 2. T_scaled must be a subset of S_scaled
      if (!.is_subset(T_scaled, S_scaled)) next

      # 3. S must have at least one scaled component (otherwise no transformation needed)
      if (length(S_scaled) == 0) {
        # No scaling for this source term - keep identity (already set)
        next
      }

      # Compute the coefficient
      extra_scaled <- setdiff(S_scaled, T_scaled)

      # Sign: (-1)^|extra|
      sign <- (-1)^length(extra_scaled)

      # Product of means for extra scaled components
      if (length(extra_scaled) > 0) {
        extra_params <- paste0(prefix, "_", extra_scaled)
        mean_product <- prod(sapply(extra_params, function(p) formula_scale[[p]]$mean))
      } else {
        mean_product <- 1
      }

      # Product of SDs for all scaled components in S
      S_scaled_params <- paste0(prefix, "_", S_scaled)
      sd_product <- prod(sapply(S_scaled_params, function(p) formula_scale[[p]]$sd))

      # Contribution coefficient
      M[t_idx, s_idx] <- sign * mean_product / sd_product
    }
  }

  return(M)
}


# Helper: Apply unscaling transformation to a matrix of posterior samples
#
# @param posterior Matrix with samples in rows, parameters in columns
# Apply the unscaling transformation to posterior samples
#
# @param posterior Matrix of posterior samples with parameter names as column names
# @param formula_scale Nested list with scaling info keyed by parameter name:
#   list(mu = list(mu_x1 = list(mean, sd)), log_sigma = list(log_sigma_x = list(mean, sd)))
# @return Transformed posterior matrix
.apply_unscale_transform <- function(posterior, formula_scale) {

  if (is.null(formula_scale) || length(formula_scale) == 0) {
    return(posterior)
  }

  # Handle nested structure: iterate over each parameter
  matched_prefix <- FALSE
  for (param_name in names(formula_scale)) {
    param_scale <- formula_scale[[param_name]]
    affected_cols <- grep(paste0("^", param_name, "_"), colnames(posterior), value = TRUE)
    if(length(affected_cols) == 0)
      next

    matched_prefix <- TRUE
    posterior <- .apply_unscale_transform_single(posterior, param_scale, prefix = param_name)
  }

  if(!matched_prefix){
    warning(
      "Ignoring formula_scale because none of its parameter prefixes match the posterior columns.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  return(posterior)
}

# Helper: Apply unscaling for a single parameter's predictors
# @param posterior Matrix of posterior samples
# @param formula_scale Flat list of scaling info: list(mu_x1 = list(mean, sd), mu_x2 = list(mean, sd))
# @param prefix Parameter prefix (e.g., "mu")
# @return Transformed posterior matrix
.apply_unscale_transform_single <- function(posterior, formula_scale, prefix) {

  if (is.null(formula_scale) || length(formula_scale) == 0) {
    return(posterior)
  }

  # Check if this parameter uses log(intercept)
  log_intercept <- isTRUE(attr(formula_scale, "log_intercept"))
  intercept_col <- paste0(prefix, "_intercept")

  # Identify which columns are affected by the transformation
  affected_cols <- grep(paste0("^", prefix, "_"), colnames(posterior), value = TRUE)
  if (length(affected_cols) > 0) {
    posterior <- .materialize_formula_scale_point_terms(
      posterior = posterior,
      formula_scale = formula_scale,
      prefix = prefix
    )
    affected_cols <- grep(paste0("^", prefix, "_"), colnames(posterior), value = TRUE)
  }
  random_sd_cols  <- grep(paste0("^", prefix, "__xREx__"), affected_cols, value = TRUE)
  random_aux_cols <- grep(paste0("^", prefix, "__(xRE_ALLOCx|xRE_SUMMARY__)"), affected_cols, value = TRUE)
  fixed_cols      <- setdiff(affected_cols, c(random_sd_cols, random_aux_cols))

  if (length(affected_cols) == 0) {
    return(posterior)
  }

  if(length(fixed_cols) > 0){
    .warn_unused_formula_scale_terms(fixed_cols, formula_scale, prefix)
  }

  # For log(intercept): transform to log scale before unscaling, then exp() back
  # This works because: log_sigma = log(intercept) + beta * x_z
  # is equivalent to: log_sigma = log_int + beta * x_z (standard additive form)
  # where log_int = log(intercept)
  if (length(fixed_cols) > 0 && log_intercept && intercept_col %in% colnames(posterior)) {
    posterior[, intercept_col] <- log(posterior[, intercept_col])
  }

  # Build and apply standard transformation matrix
  if(length(fixed_cols) > 0){
    M <- .build_unscale_matrix(fixed_cols, formula_scale, prefix)
    posterior[, fixed_cols] <- posterior[, fixed_cols, drop = FALSE] %*% t(M)
  }

  # Transform intercept back from log scale
  if (length(fixed_cols) > 0 && log_intercept && intercept_col %in% colnames(posterior)) {
    posterior[, intercept_col] <- exp(posterior[, intercept_col])
  }

  posterior <- .apply_random_sd_unscale(posterior, random_sd_cols, formula_scale, prefix)

  return(posterior)
}
