.formula_scale_matches_prefix <- function(x, prefix, separator = "_"){

  stem <- paste0(prefix, separator)
  !is.na(x) & startsWith(x, stem)
}


.formula_scale_strip_prefix <- function(x, prefix, separator = "_"){

  matched <- .formula_scale_matches_prefix(x, prefix, separator)
  if(any(matched)){
    stem_length <- nchar(paste0(prefix, separator))
    x[matched] <- substring(x[matched], stem_length + 1L)
  }

  x
}


# Helper: Parse a term name into its component variable names
# e.g., "mu_x1__xXx__x2" with prefix "mu" -> c("x1", "x2")
# e.g., "mu_intercept" -> character(0) (intercept has no components)
# e.g., "mu_x1" -> c("x1")
.parse_term_components <- function(term_name, prefix) {
  # Remove prefix
  term_part <- .formula_scale_strip_prefix(term_name, prefix)

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

  term_part <- .formula_scale_strip_prefix(term_name, prefix)

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

  scaled_vars <- .formula_scale_strip_prefix(scaled_terms, prefix)
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
      scaled_vars <- .formula_scale_strip_prefix(names(formula_scale), prefix)
    }else{
      scaled_vars <- names(formula_scale)
    }
  }

  random_terms <- character()

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


# Helper: Build the transformation matrix for unscaling fixed coefficients
#
# For each target term T and source term S, computes M[T,S] such that:
#   coef_orig[T] = sum over S of M[T,S] * coef_z[S]
#
# When the formula-scale metadata carry the fitted fixed-effect design
# (attribute "unscale_design"), the matrix is derived from that design and
# verified exactly (see .bt_formula_unscale_design_transform()). Otherwise, the
# matrix is paired by coefficient names (.build_unscale_matrix_by_names()).
#
# @param term_names Character vector of all term names in the posterior
# @param formula_scale Named list with scaling info (mean, sd) for scaled predictors
# @param prefix The parameter prefix (e.g., "mu")
# @param require_closure Whether missing induced coefficient terms should fail.
# @return A square transformation matrix
.build_unscale_matrix <- function(term_names, formula_scale, prefix,
                                  require_closure = TRUE) {

  design_spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  if(is.null(design_spec) || !isTRUE(require_closure)){
    return(.build_unscale_matrix_by_names(
      term_names = term_names,
      formula_scale = formula_scale,
      prefix = prefix,
      require_closure = require_closure
    ))
  }

  design_transform <- .bt_formula_unscale_design_transform(
    spec = design_spec,
    formula_scale = formula_scale,
    prefix = prefix
  )
  .bt_formula_unscale_design_submatrix(
    design_transform = design_transform,
    term_names = term_names,
    prefix = prefix
  )
}


# Helper: Build the unscaling matrix by pairing coefficient names
#
# The formula is based on expanding products of (x_i - mu_i)/sigma_i terms.
# For S to contribute to T:
#   1. T_unscaled == S_unscaled (unscaled components must match exactly)
#   2. T_scaled is a subset of S_scaled
#
# The contribution is: (-1)^|extra| * prod(mu_extra) / prod(sigma_S_scaled)
# where extra = S_scaled \ T_scaled
#
# The pairing assumes that column k of an interaction with a factor is coded
# like column k of the lower-order factor term. Fixed-effect transforms of
# fitted formulas are verified against the fitted design instead; this
# name-based map is used directly only for covariance-space random-effect SD
# transforms and for formula-scale metadata without a stored design.
#
# @param require_closure Whether missing induced coefficient terms should fail.
#   Random-effect SD transforms operate in covariance space and set this to FALSE.
.build_unscale_matrix_by_names <- function(term_names, formula_scale, prefix,
                                           require_closure = TRUE) {

  n_terms <- length(term_names)
  M <- diag(n_terms)  # Start with identity matrix
  rownames(M) <- colnames(M) <- term_names

  # Extract the variable names that are scaled (without prefix)
  scaled_vars <- .formula_scale_strip_prefix(names(formula_scale), prefix)

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


.bt_formula_unscale_design_spec_version <- 1L

# Helper: Minimal fitted-design metadata for exact fixed-effect unscaling: the
# fixed-effect formula, persisted factor levels and concrete contrasts, and the
# fitted column layout that defines the coefficient names.
.bt_formula_unscale_design_spec <- function(design){

  if(!inherits(design, "BayesTools_formula_design")){
    return(NULL)
  }

  predictors <- as.character(design$predictors)
  predictor_types <- design$predictor_types[predictors]
  factors <- predictors[predictor_types == "factor"]
  factor_levels <- lapply(factors, function(factor_name){
    as.character(design$xlevels[[factor_name]])
  })
  names(factor_levels) <- factors
  factor_ordered <- vapply(factors, function(factor_name){
    is.ordered(design$model_frame[[factor_name]])
  }, logical(1))
  names(factor_ordered) <- factors
  contrast_matrices <- design$contrast_matrices[factors]
  names(contrast_matrices) <- factors
  formula <- design$formula
  environment(formula) <- emptyenv()

  list(
    schema_version    = .bt_formula_unscale_design_spec_version,
    formula           = formula,
    continuous        = predictors[predictor_types == "continuous"],
    factor_levels     = factor_levels,
    factor_ordered    = factor_ordered,
    contrast_matrices = contrast_matrices,
    model_terms       = as.character(design$model_terms),
    assign            = as.integer(design$assign),
    raw_column_names  = as.character(design$raw_column_names)
  )
}

.bt_formula_unscale_design_spec_check <- function(spec, prefix){

  factor_names <- names(spec$factor_levels)
  valid <- is.list(spec) &&
    identical(spec$schema_version, .bt_formula_unscale_design_spec_version) &&
    inherits(spec$formula, "formula") &&
    is.character(spec$continuous) && !anyNA(spec$continuous) &&
    is.list(spec$factor_levels) &&
    (length(spec$factor_levels) == 0L || (
      !is.null(factor_names) && !anyNA(factor_names) &&
        all(nzchar(factor_names)) && !anyDuplicated(factor_names))) &&
    is.logical(spec$factor_ordered) &&
    identical(names(spec$factor_ordered), factor_names) &&
    is.list(spec$contrast_matrices) &&
    identical(names(spec$contrast_matrices), factor_names) &&
    all(vapply(factor_names, function(factor_name){
      levels <- spec$factor_levels[[factor_name]]
      contrast <- spec$contrast_matrices[[factor_name]]
      is.character(levels) && length(levels) > 0L && !anyNA(levels) &&
        is.matrix(contrast) && is.numeric(contrast) &&
        nrow(contrast) == length(levels) && all(is.finite(contrast))
    }, logical(1))) &&
    is.character(spec$model_terms) && !anyNA(spec$model_terms) &&
    is.integer(spec$assign) && !anyNA(spec$assign) &&
    is.character(spec$raw_column_names) &&
    length(spec$raw_column_names) == length(spec$assign) &&
    length(spec$assign) > 0L
  if(!isTRUE(valid)){
    stop(
      "Formula-scale design metadata for parameter '", prefix,
      "' are malformed. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Helper: Attach the fitted fixed-effect design to formula-scale metadata so
# that fixed-coefficient unscaling can be derived from the design rather than
# from coefficient names.
.bt_formula_scale_with_unscale_design <- function(formula_scale, design){

  if(is.null(formula_scale) || length(formula_scale) == 0L ||
     !is.null(attr(formula_scale, "unscale_design", exact = TRUE))){
    return(formula_scale)
  }

  spec <- .bt_formula_unscale_design_spec(design)
  if(!is.null(spec)){
    attr(formula_scale, "unscale_design") <- spec
  }

  formula_scale
}

.bt_formula_scale_list_with_unscale_designs <- function(formula_scale, designs){

  if(is.null(formula_scale) || length(formula_scale) == 0L ||
     !is.list(designs) || length(designs) == 0L){
    return(formula_scale)
  }

  for(parameter in intersect(names(formula_scale), names(designs))){
    formula_scale[[parameter]] <- .bt_formula_scale_with_unscale_design(
      formula_scale[[parameter]],
      designs[[parameter]]
    )
  }

  formula_scale
}

# Helper: Synthetic standardized and original-scale data for the fitted
# fixed-effect design.
#
# Rows cross every combination of the persisted factor levels with one point
# per product set of continuous predictors: every subset of the continuous
# predictors of each formula term. Within one factor cell, all columns of both
# designs are multilinear polynomials spanned by these products (centering a
# product only adds its subsets). Evaluating that span at two levels per
# predictor, raising exactly the predictors of one product set per row, is
# unisolvent (the evaluation matrix is triangular under set inclusion). A
# linear identity between the designs on these rows therefore holds for all
# predictor values, and the original-scale design has full column rank on these
# rows exactly when its columns are linearly independent functions.
#
# Scaled predictors take the levels m and m + s (standardized values 0 and 1);
# unscaled continuous predictors take the levels 0 and 1.
.bt_formula_unscale_design_data <- function(spec, formula_scale, prefix){

  continuous <- spec$continuous
  scaled_names <- paste0(prefix, "_", continuous)
  is_scaled <- scaled_names %in% names(formula_scale)

  factor_table <- attr(stats::terms(spec$formula), "factors")
  product_sets <- list(character())
  if(length(factor_table) > 0L && length(continuous) > 0L){
    for(term_i in seq_len(ncol(factor_table))){
      term_variables <- rownames(factor_table)[factor_table[, term_i] > 0]
      term_continuous <- intersect(term_variables, continuous)
      for(size in seq_along(term_continuous)){
        product_sets <- c(
          product_sets,
          utils::combn(term_continuous, size, simplify = FALSE)
        )
      }
    }
  }
  set_keys <- vapply(
    product_sets,
    function(set) paste(sort(set), collapse = "\r"),
    character(1)
  )
  product_sets <- product_sets[!duplicated(set_keys)]

  if(length(spec$factor_levels) > 0L){
    cells <- expand.grid(
      spec$factor_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }else{
    cells <- data.frame(row.names = 1L)
  }
  n_cells <- nrow(cells)
  n_rows <- n_cells * length(product_sets)
  n_columns <- length(spec$assign)
  if(n_rows * n_columns > 5e7){
    stop(
      "Cannot verify the original-scale coefficient transformation for formula parameter '",
      prefix, "' because its design has too many factor-level combinations.",
      call. = FALSE
    )
  }
  cell_index <- rep(seq_len(n_cells), times = length(product_sets))
  set_index <- rep(seq_along(product_sets), each = n_cells)

  standardized <- data.frame(row.names = seq_len(n_rows))
  original <- data.frame(row.names = seq_len(n_rows))
  for(factor_name in names(spec$factor_levels)){
    levels <- spec$factor_levels[[factor_name]]
    values <- cells[[factor_name]][cell_index]
    column <- if(isTRUE(spec$factor_ordered[[factor_name]])){
      ordered(values, levels = levels)
    }else{
      factor(values, levels = levels)
    }
    attr(column, "contrasts") <- spec$contrast_matrices[[factor_name]]
    standardized[[factor_name]] <- column
    original[[factor_name]] <- column
  }
  for(j in seq_along(continuous)){
    raised <- vapply(
      product_sets[set_index],
      function(set) continuous[j] %in% set,
      logical(1)
    )
    level <- as.numeric(raised)
    standardized[[continuous[j]]] <- level
    if(is_scaled[j]){
      scale_info <- formula_scale[[scaled_names[j]]]
      original[[continuous[j]]] <- scale_info[["mean"]] + scale_info[["sd"]] * level
    }else{
      original[[continuous[j]]] <- level
    }
  }

  list(standardized = standardized, original = original)
}

.bt_formula_unscale_model_matrix <- function(spec, data, prefix){

  # Stored formulas drop their environment; fixed formulas contain only
  # data-column names, so base functions suffice for model-frame evaluation.
  formula <- spec$formula
  environment(formula) <- baseenv()
  model_frame <- stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.fail
  )
  model_matrix <- .bt_model_matrix(model_frame, formula = formula, data = data)
  if(!identical(colnames(model_matrix), spec$raw_column_names) ||
     !identical(as.integer(attr(model_matrix, "assign")), spec$assign)){
    stop(
      "Formula-scale design metadata for parameter '", prefix,
      "' do not reproduce the fitted design. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }

  matrix(
    as.numeric(model_matrix),
    nrow = nrow(model_matrix),
    ncol = ncol(model_matrix)
  )
}

# Helper: Fitted coefficient names of the fixed-effect design columns. Every
# formula term is fitted as one coefficient vector (inprod over its design
# columns), indexed unless the term has a single column.
.bt_formula_unscale_coefficient_names <- function(spec, prefix){

  assign <- spec$assign
  has_intercept <- any(assign == 0L)
  term_index <- assign + if(has_intercept) 1L else 0L
  if(any(term_index < 1L) || any(term_index > length(spec$model_terms)) ||
     (has_intercept && !identical(spec$model_terms[1L], "intercept"))){
    stop(
      "Formula-scale design metadata for parameter '", prefix,
      "' do not match the fitted formula terms. Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }

  out <- character(length(assign))
  for(term in unique(assign)){
    columns <- which(assign == term)
    base <- paste0(prefix, "_", spec$model_terms[term_index[columns[1L]]])
    out[columns] <- if(length(columns) == 1L){
      base
    }else{
      paste0(base, "[", seq_along(columns), "]")
    }
  }

  out
}

.bt_formula_unscale_term_labels <- function(spec, columns){

  assign <- spec$assign[columns]
  has_intercept <- any(spec$assign == 0L)
  terms <- spec$model_terms[assign + if(has_intercept) 1L else 0L]
  unique(gsub("__xXx__", ":", terms, fixed = TRUE))
}

# Helper: Exact fixed-effect unscaling matrix derived from the fitted design.
#
# The standardized design X_s (as fitted) and the original-scale design X_o
# (same terms and contrasts, scaling disabled) are evaluated on the synthetic
# rows of .bt_formula_unscale_design_data(). The linear predictor is preserved,
# X_s b_s = X_o b_o, exactly when X_s = X_o A and b_o = A b_s. The name-paired
# matrix is used when it satisfies this identity (it is then the unique
# solution for full-rank designs and keeps its exact structural zeros);
# otherwise the least-squares solution is used when it is exact. A formula whose
# centered terms induce effects outside the fitted design has no exact
# solution and is rejected.
.bt_formula_unscale_design_transform <- function(spec, formula_scale, prefix){

  .bt_formula_unscale_design_spec_check(spec, prefix)
  data <- .bt_formula_unscale_design_data(spec, formula_scale, prefix)
  X_s <- .bt_formula_unscale_model_matrix(spec, data$standardized, prefix)
  X_o <- .bt_formula_unscale_model_matrix(spec, data$original, prefix)
  coefficient_names <- .bt_formula_unscale_coefficient_names(spec, prefix)

  tolerance <- 1e-8 * max(1, max(abs(X_s)))
  residual <- function(A){
    abs(X_o %*% A - X_s)
  }

  name_paired <- .build_unscale_matrix_by_names(
    term_names = coefficient_names,
    formula_scale = formula_scale,
    prefix = prefix,
    require_closure = FALSE
  )
  if(max(residual(name_paired)) <= tolerance){
    return(list(
      matrix = name_paired,
      method = "name_paired",
      coefficient_names = coefficient_names
    ))
  }

  qr_original <- qr(X_o)
  if(qr_original$rank < ncol(X_o)){
    .bt_formula_transform_stop(
      paste0(
        "Cannot transform the coefficients of formula parameter '", prefix,
        "' to the original predictor scale: the fitted formula terms are linearly ",
        "dependent, so the original-scale coefficients are not uniquely determined."
      ),
      parameter = prefix,
      reason = "original_scale_not_identified"
    )
  }
  least_squares <- qr.coef(qr_original, X_s)
  column_scale <- pmax(1, apply(abs(least_squares), 2L, max))
  zapped <- least_squares
  zapped[abs(zapped) <= 1e-10 * rep(column_scale, each = nrow(zapped))] <- 0
  if(max(residual(zapped)) <= tolerance){
    least_squares <- zapped
  }
  column_residual <- apply(residual(least_squares), 2L, max)
  if(any(column_residual > tolerance)){
    terms <- .bt_formula_unscale_term_labels(
      spec,
      which(column_residual > tolerance)
    )
    .bt_formula_transform_stop(
      paste0(
        "Cannot transform the coefficients of formula parameter '", prefix,
        "' to the original predictor scale: centering the standardized predictors in term",
        if(length(terms) > 1L) "s " else " ",
        paste0("'", terms, "'", collapse = ", "),
        " induces lower-order effects that the formula does not contain. ",
        "Include the corresponding lower-order terms or keep the coefficients on ",
        "the fitted, standardized scale."
      ),
      parameter = prefix,
      reason = "original_scale_not_representable",
      terms = terms
    )
  }
  dimnames(least_squares) <- list(coefficient_names, coefficient_names)

  list(
    matrix = least_squares,
    method = "least_squares",
    coefficient_names = coefficient_names
  )
}

# Helper: Restrict an exact design-level unscaling matrix to the requested
# posterior columns. Columns that are not fixed-design coefficients keep the
# identity transformation.
.bt_formula_unscale_design_submatrix <- function(design_transform, term_names,
                                                 prefix){

  A <- design_transform$matrix
  design_names <- rownames(A)
  out <- diag(length(term_names))
  dimnames(out) <- list(term_names, term_names)
  present <- term_names[term_names %in% design_names]
  if(length(present) == 0L){
    return(out)
  }

  missing <- setdiff(design_names, term_names)
  if(length(missing) > 0L){
    for(target in present){
      absent <- intersect(design_names[A[target, ] != 0], missing)
      if(length(absent) > 0L){
        stop(
          "Cannot unscale posterior term '", target,
          "' because its original-scale value depends on the missing fitted ",
          "coefficient(s) ", paste0("'", absent, "'", collapse = ", "),
          ". Include the corresponding posterior columns.",
          call. = FALSE
        )
      }
    }
    for(source in present){
      induced <- setdiff(design_names[A[, source] != 0], source)
      if(any(induced %in% present) && any(induced %in% missing)){
        stop(
          "Cannot unscale posterior term '", source,
          "' because centering induces missing lower-order coefficient(s) ",
          paste0("'", intersect(induced, missing), "'", collapse = ", "),
          ". Include the corresponding lower-order posterior columns.",
          call. = FALSE
        )
      }
    }
  }

  out[present, present] <- A[present, present]
  out
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
    affected_cols <- colnames(posterior)[
      .formula_scale_matches_prefix(colnames(posterior), param_name)
    ]
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

  # Identify which columns are affected by the transformation
  affected_cols <- colnames(posterior)[
    .formula_scale_matches_prefix(colnames(posterior), prefix)
  ]
  if (length(affected_cols) > 0) {
    posterior <- .materialize_formula_scale_point_terms(
      posterior = posterior,
      formula_scale = formula_scale,
      prefix = prefix
    )
    affected_cols <- colnames(posterior)[
      .formula_scale_matches_prefix(colnames(posterior), prefix)
    ]
  }
  random_sd_cols <- affected_cols[
    .formula_scale_matches_prefix(affected_cols, prefix, "__xREx__")
  ]
  random_aux_cols <- affected_cols[
    .formula_scale_matches_prefix(affected_cols, prefix, "__xRE_ALLOCx") |
      .formula_scale_matches_prefix(affected_cols, prefix, "__xRE_SUMMARY__")
  ]
  fixed_cols      <- setdiff(affected_cols, c(random_sd_cols, random_aux_cols))

  if (length(affected_cols) == 0) {
    return(posterior)
  }

  if(length(fixed_cols) > 0){
    .warn_unused_formula_scale_terms(fixed_cols, formula_scale, prefix)
  }

  if(length(fixed_cols) > 0){
    transform <- .bt_formula_coefficient_transform(
      source_names = fixed_cols,
      formula_scale = formula_scale,
      parameter = prefix
    )
    posterior <- .bt_apply_formula_coefficient_transform(
      posterior,
      transform
    )
  }

  posterior <- .apply_random_sd_unscale(posterior, random_sd_cols, formula_scale, prefix)

  return(posterior)
}
