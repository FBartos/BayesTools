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
  # trailing index, a level `[..]` or a contrast coefficient `{j}`, so scaled
  # continuous components inside interactions are still detected correctly
  # (e.g., alloc__xXx__year[1] -> alloc, year).
  term_part <- sub("(\\[[^]]+\\]|\\{[^}]+\\})$", "", term_part)

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

  # A trailing level `[..]` or contrast coefficient `{j}` index.
  term_index <- ""
  if (grepl("(\\[[^]]+\\]|\\{[^}]+\\})$", term_part)) {
    term_index <- sub("^.*(\\[[^]]+\\]|\\{[^}]+\\})$", "\\1", term_part)
    term_core  <- sub("(\\[[^]]+\\]|\\{[^}]+\\})$", "", term_part)
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
    (nzchar(target_id) && !grepl("(\\[[^]]+\\]|\\{[^}]+\\})$", target_id) &&
       source_id %in% paste0(target_id, c("[1]", "{1}")))
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
    if(is.list(param_scale) && length(param_scale) == 0L){
      spec <- attr(param_scale, "unscale_design", exact = TRUE)
      if(!is.null(spec)) .bt_formula_unscale_design_spec_check(spec, param_name)
      next
    }
    entry <- paste0(name, "[['", param_name, "']]")
    check_list(param_scale, entry)

    if(length(param_scale) == 0)
      next

    if(is.null(names(param_scale)) || anyNA(names(param_scale)) || any(names(param_scale) == ""))
      stop(paste0("The '", entry, "' entry must be a named list keyed by parameter term."), call. = FALSE)
    if(anyDuplicated(names(param_scale)))
      stop(paste0("The '", entry, "' entry must not contain duplicate term names."), call. = FALSE)

    for(term_name in names(param_scale)){
      term_scale <- param_scale[[term_name]]
      term_entry <- paste0(entry, "[['", term_name, "']]")

      check_list(
        term_scale,
        term_entry,
        check_names = c("mean", "sd"),
        all_objects = TRUE,
        allow_other = TRUE
      )
      check_real(
        term_scale[["mean"]],
        paste0(term_entry, "[['mean']]"),
        allow_NA = FALSE
      )
      if(!is.finite(term_scale[["mean"]])){
        stop(
          "The '", term_entry, "[['mean']]' entry must be finite.",
          call. = FALSE
        )
      }
      check_real(
        term_scale[["sd"]],
        paste0(term_entry, "[['sd']]"),
        lower = 0,
        allow_bound = FALSE,
        allow_NA = FALSE
      )
      if(!is.finite(term_scale[["sd"]])){
        stop(
          "The '", term_entry, "[['sd']]' entry must be finite.",
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
        "Ignoring all formula_scale[['", prefix, "']] entries because none match posterior terms: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Samples for this parameter prefix are returned unchanged."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }else{
    warning(
      paste0(
        "Ignoring unused formula_scale[['", prefix, "']] entries: '",
        paste0(unused_terms, collapse = "', '"),
        "'. Matched entries are still applied."
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  invisible(NULL)
}

.formula_scale_random_scaled_vars <- function(formula_scale, scaled_vars){

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
# The matrix is derived from the fitted fixed-effect design that the
# formula-scale metadata carry (attribute "unscale_design") and verified
# exactly (see .bt_formula_unscale_design_transform()); the columns must be
# fitted coefficient coordinates (callers map labelled columns to their
# coordinates first). Formula-scale metadata without the fitted design (built
# by hand) cannot be transformed: pairing coefficient names misses the
# contributions of terms whose coding differs from their names (for example,
# the level slopes of `~ g + g:x` to the intercept and the levels).
#
# @param term_names Character vector of all term names in the posterior
# @param formula_scale Named list with scaling info (mean, sd) for scaled predictors
# @param prefix The parameter prefix (e.g., "mu")
# @return A square transformation matrix
.build_unscale_matrix <- function(term_names, formula_scale, prefix) {

  design_spec <- .bt_formula_scale_require_design(formula_scale, prefix)
  if(!.bt_formula_unscale_design_columns(design_spec, term_names, prefix)){
    stop(
      "Cannot transform the coefficients of formula parameter '", prefix,
      "' to the original predictor scale: the posterior columns are not the ",
      "fitted coefficient coordinates of its design.",
      call. = FALSE
    )
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

# Helper: The fitted design that the formula-scale metadata of formula
# parameter 'prefix' carry (attribute "unscale_design"). Every original-scale
# transform of its coefficients and random-effect SDs requires it: the
# formula_scale of a fitted model carries it together with the rest of the
# fitted structure (log intercept, point terms, random-effect SD structure),
# while standardization information built by hand carries none of it and
# stops.
.bt_formula_scale_require_design <- function(formula_scale, prefix){

  design_spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  if(is.null(design_spec)){
    .bt_formula_transform_stop(
      paste0(
        "Cannot transform the coefficients of formula parameter '", prefix,
        "' to the original predictor scale: 'formula_scale' does not carry ",
        "the fitted design of the formula. Use the 'formula_scale' attribute ",
        "of the fitted model (attr(fit, \"formula_scale\")), which carries it."
      ),
      parameter = prefix,
      reason    = "missing_fitted_design"
    )
  }

  design_spec
}

# Helper: Standardization information passed for the fitted model 'fit',
# completed with the fitted formula-scale metadata of the model that it does
# not carry itself: the fitted design, the log-intercept flag, the point
# terms, and the random-effect SD structure. A formula_scale list built by
# hand carries only means and SDs; without the rest of the fitted structure,
# random-effect SDs would stay standardized and a log intercept would be
# transformed as a linear one.
.bt_formula_scale_list_complete <- function(formula_scale, fit){

  if(is.null(formula_scale) || length(formula_scale) == 0L){
    return(formula_scale)
  }
  .check_formula_scale_info(formula_scale)

  fitted_scale <- attr(fit, "formula_scale", exact = TRUE)
  # A formula parameter that the fitted model does not have is an argument
  # error, not a missing declaration of the fit: also for a fit without any
  # formula parameter, for which a refit cannot help.
  fitted_parameters <- unique(c(
    names(attr(fit, "formula_design", exact = TRUE)),
    names(fitted_scale)
  ))
  unknown <- setdiff(names(formula_scale), fitted_parameters)
  if(length(unknown) > 0L){
    stop(
      "The 'formula_scale' argument names unknown formula parameter",
      if(length(unknown) > 1L) "s", " ", paste0("'", unknown, "'", collapse = ", "),
      if(length(fitted_parameters) > 0L){
        paste0(
          ". The fitted model has the formula parameter",
          if(length(fitted_parameters) > 1L) "s", " ",
          paste0("'", fitted_parameters, "'", collapse = ", "), "."
        )
      }else{
        ". The fitted model has no formula parameters."
      },
      call. = FALSE
    )
  }
  if(is.list(fitted_scale)){
    for(parameter in intersect(names(formula_scale), names(fitted_scale))){
      # The override replaces the fitted scale: a predictor it omits would be
      # treated as unscaled, giving a result of mixed scales.
      omitted <- setdiff(names(fitted_scale[[parameter]]), names(formula_scale[[parameter]]))
      if(length(omitted) > 0L){
        .bt_formula_transform_stop(
          paste0(
            "The 'formula_scale' entry for parameter '", parameter, "' omits the standardized predictor",
            if(length(omitted) > 1L) "s", " ", paste0("'", omitted, "'", collapse = ", "),
            " of the fitted model. Supply the mean and SD of every standardized predictor, ",
            "or leave 'formula_scale' unspecified to use the fitted values."
          ),
          parameter = parameter, reason = "incomplete_formula_scale",
          predictors = omitted
        )
      }
      fitted_attributes <- attributes(fitted_scale[[parameter]])
      .bt_formula_scale_finalized_check(fitted_scale[[parameter]], parameter)
      attributes(formula_scale[[parameter]]) <- c(
        list(names = names(formula_scale[[parameter]])),
        fitted_attributes[setdiff(names(fitted_attributes), "names")]
      )
    }
  }

  for(parameter in names(formula_scale)){
    .bt_formula_scale_finalized_check(formula_scale[[parameter]], parameter)
    design <- attr(fit, "formula_design", exact = TRUE)[[parameter]]
    if(is.null(design) || !identical(attr(design$formula_scale, "unscale_design", exact = TRUE),
                                    attr(fitted_scale[[parameter]], "unscale_design", exact = TRUE))){
      .bt_stop_refit_required("Formula declaration carriers disagree. Refit the model with this version of BayesTools.")
    }
    # A fit without automatic scaling stores only the fitted design: the
    # structure follows from it, as if the predictors had been scaled.
    if(is.null(attr(formula_scale[[parameter]], "parameter", exact = TRUE))){
      formula_scale[[parameter]] <- .bt_formula_scale_design_structure(formula_scale[[parameter]], design)
    }
  }
  formula_scale
}

# The formula-scale structure beyond the means and SDs: the formula parameter,
# the log-intercept flag, the numeric point terms, and the random-effect SD
# leaves of each block (keyed '__xREx__<block>') with the blocks whose
# correlation the SD transform needs. JAGS_formula() attaches it to the scale
# it builds; 'random_sd_leaves = NULL' attaches no random-effect structure.
.bt_formula_scale_attach_structure <- function(formula_scale, parameter, log_intercept,
                                               point_terms = NULL, random_sd_leaves = NULL,
                                               random_correlation_required = character()){

  attr(formula_scale, "parameter") <- parameter
  attr(formula_scale, "log_intercept") <- log_intercept
  if(length(point_terms) > 0L){
    attr(formula_scale, "point_terms") <- point_terms
  }
  if(!is.null(random_sd_leaves)){
    attr(formula_scale, "random_effect_sd_leaves") <- random_sd_leaves
    if(length(random_correlation_required) > 0L){
      attr(formula_scale, "random_effect_correlation_required") <- unique(random_correlation_required)
    }
  }
  formula_scale
}

# The structure JAGS_formula() would have attached had the predictors of the
# fitted design been scaled automatically, taken from the stored design.
.bt_formula_scale_design_structure <- function(formula_scale, design){

  random_sd_leaves <- NULL
  random_correlation_required <- character()
  random_terms <- .bt_formula_design_random_effects(design)
  if(length(random_terms) > 0L){
    random_sd_leaves <- list()
    for(random_term in random_terms){
      random_key <- paste0("__xREx__", random_term$block_name)
      random_sd_leaves[[random_key]] <- random_term$sd_leaves
      n_columns <- random_term$n_columns
      if(.bt_random_effect_structure(random_term) %in% c("us", "cs", "hcs", "ar1", "car", "har") &&
         is.numeric(n_columns) && length(n_columns) == 1L && !is.na(n_columns) && n_columns > 1L){
        random_correlation_required <- c(random_correlation_required, random_key)
      }
    }
  }
  .bt_formula_scale_attach_structure(
    formula_scale,
    parameter = design$parameter,
    log_intercept = isTRUE(design$log_intercept),
    point_terms = .formula_scale_point_terms(
      prior_list = design$prior_list,
      parameter = design$parameter,
      model_terms = design$model_terms
    ),
    random_sd_leaves = random_sd_leaves,
    random_correlation_required = random_correlation_required
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
# like column k of the lower-order factor term. Fixed-effect transforms are
# derived from the fitted design (.build_unscale_matrix()); this name-based map
# is used only for covariance-space random-effect SD transforms and as the
# first candidate that the design verification accepts or rejects.
.build_unscale_matrix_by_names <- function(term_names, formula_scale, prefix) {

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

  # Build the transformation matrix
  for (t_idx in seq_along(term_names)) {
    T_name <- term_names[t_idx]
    T_scaled <- term_scaled[[T_name]]
    T_unscaled <- term_unscaled[[T_name]]

    for (s_idx in seq_along(term_names)) {
      S_name <- term_names[s_idx]
      S_scaled <- term_scaled[[S_name]]
      S_unscaled <- term_unscaled[[S_name]]

      # Columns map to coordinates by position: another column of the same term
      # name (the columns of one shared SD leaf) is its own coordinate and puts
      # no weight on this one
      if (s_idx != t_idx && identical(S_name, T_name)) next

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

      M[t_idx, s_idx] <- .bt_formula_unscale_affine_coefficient(
        S_scaled, extra_scaled, formula_scale, prefix, S_name, T_name
      )
    }
  }

  return(M)
}


# Expand one numeric monomial without turning range loss into a structural zero.
# Keep the ordinary grouping when its intermediate products are normal, then
# try a whole product of ratios. Base prod() can preserve a representable whole
# product even when a sequential double product would lose an intermediate.
.bt_formula_unscale_affine_coefficient <- function(scaled, dropped,
                                                   formula_scale, prefix,
                                                   source, target){

  means <- vapply(dropped, function(v) formula_scale[[paste0(prefix, "_", v)]]$mean, numeric(1))
  if(any(means == 0)) return(0)
  sds <- vapply(scaled, function(v) formula_scale[[paste0(prefix, "_", v)]]$sd, numeric(1))
  normal <- function(x) is.finite(x) & abs(x) >= .Machine$double.xmin
  final <- function(x) is.finite(x) & x != 0
  numerator <- prod(-means)
  denominator <- prod(sds)
  direct <- numerator / denominator
  if(normal(numerator) && normal(denominator) && final(direct)) return(direct)

  ratios <- vapply(scaled, function(v){
    info <- formula_scale[[paste0(prefix, "_", v)]]
    if(v %in% dropped) -info$mean / info$sd else 1 / info$sd
  }, numeric(1))
  alternate <- prod(ratios)
  ratios_valid <- if(length(ratios) == 1L) final(ratios) else all(normal(ratios))
  if(ratios_valid && final(alternate)) return(alternate)
  .bt_formula_transform_stop(
    "Original coefficient affine arithmetic is numerically unavailable.",
    parameter = prefix, source = source, target = target,
    reason = "nonfinite_transform",
    observed = list(numerator = numerator, denominator = denominator,
      direct = direct, ratios = ratios, alternate = alternate)
  )
}


.bt_formula_unscale_design_spec_version <- 2L

# The compiler's top-level term multiplier is the sole semantic owner. Child
# mixture/slab attributes are deliberately not consulted.
.bt_formula_compiled_multipliers <- function(design){

  out <- list()
  for(owner in intersect(names(design$prior_list), paste0(design$parameter, "_", design$model_terms))){
    prior <- design$prior_list[[owner]]
    columns <- .prior_linear_prior_columns(owner, prior)
    multiplier <- attr(prior, "multiply_by", exact = TRUE)
    if(identical(owner, paste0(design$parameter, "_intercept")) || is.null(multiplier)){
      multiplier <- 1
    }
    declaration <- if(is.numeric(multiplier) && length(multiplier) == 1L && is.finite(multiplier)){
      list(type = "constant", value = as.numeric(multiplier))
    }else if(is.character(multiplier) && length(multiplier) == 1L && !is.na(multiplier) && nzchar(multiplier)){
      list(type = "state", name = multiplier)
    }else{
      .bt_formula_transform_stop("Compiled formula multiplier is unavailable.",
        parameter = design$parameter, reason = "unsupported_multiplier_declaration")
    }
    for(column in columns) out[[column]] <- declaration
  }
  out
}

.bt_formula_multiplier_valid <- function(x){

  is.list(x) && (
    (identical(names(x), c("type", "value")) && identical(x$type, "constant") &&
       is.numeric(x$value) && length(x$value) == 1L && is.finite(x$value)) ||
    (identical(names(x), c("type", "name")) && identical(x$type, "state") &&
       is.character(x$name) && length(x$name) == 1L && !is.na(x$name) && nzchar(x$name))
  )
}

# Complete a fresh current-format owner copy; this is never a fit-reader
# migration path. Formula literal points keep their existing point_terms owner.
.bt_formula_scale_finalize <- function(formula_scale, design = NULL,
                                       prior_list = list(), model_data = NULL,
                                       owner_scope = "fit"){

  if(is.null(formula_scale)) return(NULL)
  spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  if(is.null(spec) && !is.null(design)) spec <- .bt_formula_unscale_design_spec(design)
  if(is.null(spec)){
    if(length(formula_scale) == 0L) return(formula_scale)
    .bt_stop_refit_required("Formula declarations are missing. Refit the model with this version of BayesTools.")
  }
  # The formula parameter named in refusals: the design's, else the one the
  # scale was built for.
  prefix <- if(!is.null(design)){
    design$parameter
  }else{
    scale_parameter <- attr(formula_scale, "parameter", exact = TRUE)
    if(is.character(scale_parameter) && length(scale_parameter) == 1L &&
       !is.na(scale_parameter) && nzchar(scale_parameter)){
      scale_parameter
    }else{
      "(unnamed)"
    }
  }
  .bt_formula_unscale_design_spec_check(spec, prefix)
  if(!owner_scope %in% c("fit", "prior_only")) stop("Unsupported formula declaration owner.", call. = FALSE)
  constants <- spec$state_constants
  formula_points <- attr(formula_scale, "point_terms", exact = TRUE)
  needed <- unique(vapply(Filter(function(x) identical(x$type, "state"), spec$multipliers), `[[`, character(1), "name"))
  if(!is.null(design)){
    scalar_sd <- .bt_formula_design_random_effects(design)
    scalar_sd <- Filter(function(term) !.bt_random_effect_has_row_indexed_external_sd(term), scalar_sd)
    needed <- union(needed, unlist(lapply(scalar_sd, `[[`, "sd_parameter_names"), use.names = FALSE))
    needed <- needed[!is.na(needed)]
  }
  for(name in needed){
    candidates <- numeric()
    if(name %in% names(formula_points)) candidates <- c(candidates, formula_points[[name]])
    prior <- prior_list[[name]]
    if(!is.null(prior) && is.prior.point(prior)){
      point <- prior$parameters[["location"]]
      if(is.numeric(point) && length(point) == 1L && is.finite(point)) candidates <- c(candidates, point)
    }
    data_value <- model_data[[name]]
    if(is.numeric(data_value) && length(data_value) == 1L && is.finite(data_value)) candidates <- c(candidates, data_value)
    if(name %in% names(constants)) candidates <- c(candidates, constants[[name]])
    if(length(candidates)){
      if(any(candidates != candidates[[1L]])) .bt_formula_transform_stop(
        "Formula scalar declarations have contradictory owners.", parameter = prefix,
        reason = "contradictory_state_constants", state = name)
      # A formula coefficient point remains owned by point_terms.
      if(!name %in% names(formula_points)) constants[name] <- candidates[[1L]]
    }
  }
  spec$state_constants <- constants
  spec$owner_scope <- owner_scope
  spec$complete <- TRUE
  attr(formula_scale, "unscale_design") <- spec
  formula_scale
}

#' @rdname JAGS_formula_point_finalize
#' @param formula_scale the \code{formula_scale} element returned by
#' [JAGS_formula()] (\code{NULL} is returned unchanged).
#' @param owner_scope the owner of the completed scale: \code{"fit"} for the
#' scale of a fitted model, as [JAGS_fit()] stores it, or \code{"prior_only"}
#' for a design used for prior densities without a fit.
#'
#' @return \code{JAGS_formula_scale_finalize()} returns \code{formula_scale}
#' with the completed \code{unscale_design} attribute (the resolved
#' \code{state_constants}, \code{owner_scope}, and \code{complete = TRUE}), or
#' \code{NULL} for a \code{NULL} scale.
#'
#' @export
JAGS_formula_scale_finalize <- function(formula_scale, design = NULL,
                                        prior_list = list(), model_data = NULL,
                                        owner_scope = "fit"){

  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)
  .bt_formula_finalize_design_check(design, allow_NULL = TRUE)
  check_list(prior_list, "prior_list", allow_NULL = TRUE)
  check_list(model_data, "model_data", allow_NULL = TRUE)
  check_char(owner_scope, "owner_scope", check_length = 1, allow_NA = FALSE,
             allow_values = c("fit", "prior_only"))

  .bt_formula_scale_finalize(formula_scale, design, prior_list, model_data,
                             owner_scope)
}

.bt_formula_scale_finalized_check <- function(formula_scale, prefix, require_owner = FALSE){

  if(!require_owner && (is.null(formula_scale) || length(formula_scale) == 0L)) return(invisible(TRUE))
  spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  .bt_formula_unscale_design_spec_check(spec, prefix)
  if(!isTRUE(spec$complete) || !identical(spec$owner_scope, "fit")){
    .bt_stop_refit_required("Fitted formula declarations are incomplete. Refit the model with this version of BayesTools.")
  }
  invisible(TRUE)
}

.bt_formula_state_constants <- function(formula_scale){

  spec <- attr(formula_scale, "unscale_design", exact = TRUE)
  points <- attr(formula_scale, "point_terms", exact = TRUE)
  out <- c(spec$state_constants, points)
  if(is.null(out)) out <- stats::setNames(numeric(), character())
  if(anyDuplicated(names(out))) .bt_formula_transform_stop(
    "Formula scalar declarations have contradictory owners.", reason = "contradictory_state_constants")
  out
}

.bt_formula_bridge_state_constants <- function(designs){

  out <- list()
  for(design in designs){
    constants <- .bt_formula_state_constants(design$formula_scale)
    for(name in names(constants)){
      value <- constants[[name]]
      if(!is.numeric(value) || length(value) != 1L || !is.finite(value)){
        .bt_formula_transform_stop("Retained formula scalar state is malformed.", reason = "malformed_state_constants", state = name)
      }
      if(!is.null(out[[name]]) && out[[name]] != value){
        .bt_formula_transform_stop("Formula scalar declarations have contradictory owners.", reason = "contradictory_state_constants", state = name)
      }
      out[[name]] <- value
    }
  }
  out
}

.bt_formula_materialize_state <- function(samples, constants, validate = character()){

  for(name in names(constants)){
    if(name %in% colnames(samples)){
      if(name %in% validate && (any(!is.finite(samples[, name])) || any(samples[, name] != constants[[name]]))){
        .bt_formula_transform_stop("Supplied fitted state disagrees with its declaration.",
          reason = "contradictory_state_constants", state = name)
      }
    }else{
      samples <- cbind(samples, rep(constants[[name]], nrow(samples)))
      colnames(samples)[ncol(samples)] <- name
    }
  }
  samples
}

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
  # Logical predictors are coded as indicators (e.g., column 'dTRUE') unless
  # they were standardized before the design was built.
  continuous <- predictors[predictor_types == "continuous"]
  continuous_logical <- vapply(continuous, function(predictor){
    is.logical(design$model_frame[[predictor]])
  }, logical(1))
  names(continuous_logical) <- continuous
  formula <- design$formula
  environment(formula) <- emptyenv()

  list(
    schema_version    = .bt_formula_unscale_design_spec_version,
    formula           = formula,
    continuous        = continuous,
    continuous_logical = continuous_logical,
    factor_levels     = factor_levels,
    factor_ordered    = factor_ordered,
    contrast_matrices = contrast_matrices,
    model_terms       = as.character(design$model_terms),
    assign            = as.integer(design$assign),
    raw_column_names  = as.character(design$raw_column_names),
    multipliers       = .bt_formula_compiled_multipliers(design),
    state_constants   = stats::setNames(numeric(), character()),
    owner_scope       = "compiler",
    complete          = FALSE
  )
}

.bt_formula_unscale_design_spec_check <- function(spec, prefix){

  factor_names <- names(spec$factor_levels)
  valid <- is.list(spec) &&
    identical(spec$schema_version, .bt_formula_unscale_design_spec_version) &&
    inherits(spec$formula, "formula") &&
    is.character(spec$continuous) && !anyNA(spec$continuous) &&
    is.logical(spec$continuous_logical) && !anyNA(spec$continuous_logical) &&
    length(spec$continuous_logical) == length(spec$continuous) &&
    (length(spec$continuous) == 0L ||
       identical(names(spec$continuous_logical), spec$continuous)) &&
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
    length(spec$assign) > 0L &&
    is.list(spec$multipliers) && !is.null(names(spec$multipliers)) &&
    !anyDuplicated(names(spec$multipliers)) &&
    all(vapply(spec$multipliers, .bt_formula_multiplier_valid, logical(1))) &&
    is.numeric(spec$state_constants) && all(is.finite(spec$state_constants)) &&
    (length(spec$state_constants) == 0L ||
       (!is.null(names(spec$state_constants)) && !anyNA(names(spec$state_constants)) &&
          all(nzchar(names(spec$state_constants))) && !anyDuplicated(names(spec$state_constants)))) &&
    spec$owner_scope %in% c("compiler", "fit", "prior_only") &&
    (isTRUE(spec$complete) || isFALSE(spec$complete))
  if(!isTRUE(valid)){
    .bt_stop_refit_required(
      "Formula-scale design metadata for parameter '", prefix,
      "' are malformed. Refit the model with this version of BayesTools."
    )
  }

  invisible(TRUE)
}

# Helper: Attach the fitted fixed-effect design to formula-scale metadata so
# that fixed-coefficient unscaling can be derived from the design rather than
# from coefficient names.
.bt_formula_scale_with_unscale_design <- function(formula_scale, design){

  if(!is.null(attr(formula_scale, "unscale_design", exact = TRUE))){
    return(formula_scale)
  }
  if(is.null(formula_scale)) formula_scale <- list()

  spec <- .bt_formula_unscale_design_spec(design)
  if(!is.null(spec)){
    attr(formula_scale, "unscale_design") <- spec
  }

  formula_scale
}

# Whether any formula parameter has a scaled predictor.
.bt_formula_scale_has_scaled <- function(formula_scale){

  is.list(formula_scale) && any(lengths(formula_scale) > 0L)
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
# unscaled continuous predictors take the levels 0 and 1 (FALSE and TRUE for
# logical predictors).
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

  # The guard comes before expand.grid(), which would allocate every
  # combination of the factor levels.
  n_cells <- prod(as.numeric(lengths(spec$factor_levels)))
  n_rows <- n_cells * length(product_sets)
  n_columns <- length(spec$assign)
  if(n_rows * n_columns > 5e7){
    stop(
      "Cannot verify the original-scale coefficient transformation for formula parameter '",
      prefix, "' because its design has too many factor-level combinations.",
      call. = FALSE
    )
  }
  if(length(spec$factor_levels) > 0L){
    cells <- expand.grid(
      spec$factor_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }else{
    cells <- data.frame(row.names = 1L)
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
    }else if(isTRUE(spec$continuous_logical[[continuous[j]]])){
      # unstandardized logical predictors keep their indicator coding
      standardized[[continuous[j]]] <- raised
      original[[continuous[j]]] <- raised
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
  if(attr(stats::terms(formula), "intercept") == 0L){
    # A random-effect block without an intercept keeps the contrast coding of
    # its factors: it is built as the block's own design matrix is.
    model_matrix <- .bt_random_effect_design_matrix(
      formula,
      data,
      preserve_no_intercept_contrasts = TRUE
    )$model_matrix
  }else{
    model_frame <- stats::model.frame(
      formula,
      data = data,
      na.action = stats::na.fail
    )
    model_matrix <- .bt_model_matrix(model_frame, formula = formula, data = data)
  }
  if(!identical(colnames(model_matrix), spec$raw_column_names) ||
     !identical(as.integer(attr(model_matrix, "assign")), spec$assign)){
    .bt_stop_refit_required(
      "Formula-scale design metadata for parameter '", prefix,
      "' do not reproduce the fitted design. Refit the model with this version of BayesTools."
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
    .bt_stop_refit_required(
      "Formula-scale design metadata for parameter '", prefix,
      "' do not match the fitted formula terms. Refit the model with this version of BayesTools."
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

# Helper: Whether the requested columns are in the fitted coefficient
# coordinates of the design. Columns that belong to a design term but carry
# other names (factor-level labels, contrast coefficients `{j}` of mixed
# posteriors, or more levels than fitted coefficients) are not the fitted
# coefficient vector and cannot be transformed.
.bt_formula_unscale_design_columns <- function(spec, term_names, prefix){

  .bt_formula_unscale_design_spec_check(spec, prefix)
  design_names <- .bt_formula_unscale_coefficient_names(spec, prefix)
  strip_index <- function(x) sub("(\\[[^]]+\\]|\\{[^}]+\\})$", "", x)
  foreign <- strip_index(term_names) %in% strip_index(design_names) &
    !term_names %in% design_names

  !any(foreign)
}

.bt_formula_unscale_term_labels <- function(spec, columns){

  assign <- spec$assign[columns]
  has_intercept <- any(spec$assign == 0L)
  terms <- spec$model_terms[assign + if(has_intercept) 1L else 0L]
  unique(gsub("__xXx__", ":", terms, fixed = TRUE))
}

# Number of declared categorical cells, including both values of a logical.
.bt_formula_unscale_cell_count <- function(data, variables){

  prod(vapply(variables, function(v){
    if(is.factor(data[[v]])) nlevels(data[[v]]) else 2L
  }, integer(1)))
}

# Exact observed-cell identities use persisted factor levels/logical values.
.bt_formula_unscale_cell_keys <- function(data, variables){

  if(!length(variables)) return(rep("", nrow(data)))
  codes <- lapply(data[variables], as.integer)
  do.call(paste, c(codes, sep = "\r"))
}

# A certified 0/1 lift into functions of the source categorical variables.
# Source cells must be complete. Outside cells need only be the observed cells,
# with every source/outside pair present exactly once after deduplication.
# Groups from one formula term may partition those outside cells; copying their
# common nonzero column avoids cancellation and certifies all excluded zeros.
.bt_formula_unscale_partition_basis <- function(data, C, columns, source_variables,
                                                categorical_variables, assign){

  source_count <- .bt_formula_unscale_cell_count(data, source_variables)
  if(!is.finite(source_count) || source_count > length(columns)) return(NULL)
  outside <- setdiff(categorical_variables, source_variables)
  cell_keys <- .bt_formula_unscale_cell_keys(data, categorical_variables)
  first <- which(!duplicated(cell_keys))
  source_keys <- .bt_formula_unscale_cell_keys(data, source_variables)
  outside_keys <- .bt_formula_unscale_cell_keys(data, outside)
  source_levels <- unique(source_keys[first])
  outside_levels <- unique(outside_keys[first])
  ns <- length(source_levels)
  no <- length(outside_levels)
  if(ns != source_count || length(first) != ns * no) return(NULL)
  sid <- match(source_keys[first], source_levels)
  oid <- match(outside_keys[first], outside_levels)
  if(anyDuplicated((oid - 1L) * ns + sid)) return(NULL)

  values <- C[first, columns, drop = FALSE]
  functions <- lapply(seq_along(columns), function(k){
    out <- matrix(0, ns, no)
    out[cbind(sid, oid)] <- values[, k]
    out
  })
  lifts <- list()
  basis <- list()
  separable <- list()
  for(k in seq_along(columns)){
    F <- functions[[k]]
    invariant <- all(vapply(seq_len(no), function(o) identical(F[, o], F[, 1L]), logical(1)))
    if(invariant){
      lift <- numeric(length(columns))
      lift[k] <- 1
      lifts[[length(lifts) + 1L]] <- lift
      basis[[length(basis) + 1L]] <- F[, 1L]
    }else{
      first_nonzero <- which(colSums(F != 0) > 0L)[1L]
      pattern <- F[, first_nonzero]
      mask <- vapply(seq_len(no), function(o) identical(F[, o], pattern), logical(1))
      zeros <- vapply(seq_len(no), function(o) all(F[, o] == 0), logical(1))
      if(all(mask | zeros)){
        separable[[length(separable) + 1L]] <- list(
          column = k, pattern = pattern, mask = mask, term = assign[columns[k]]
        )
      }
    }
  }
  pending <- seq_along(separable)
  while(length(pending)){
    candidate <- separable[[pending[1L]]]
    group <- pending[vapply(pending, function(k){
      identical(separable[[k]]$term, candidate$term) &&
        identical(separable[[k]]$pattern, candidate$pattern)
    }, logical(1))]
    masks <- do.call(rbind, lapply(separable[group], `[[`, "mask"))
    if(all(colSums(masks) == 1L)){
      lift <- numeric(length(columns))
      lift[vapply(separable[group], `[[`, integer(1), "column")] <- 1
      lifts[[length(lifts) + 1L]] <- lift
      basis[[length(basis) + 1L]] <- candidate$pattern
    }
    pending <- setdiff(pending, group)
  }
  if(length(lifts) != ns) return(NULL)
  L <- do.call(cbind, lifts)
  B <- do.call(cbind, basis)
  magnitude <- apply(abs(B), 2L, max)
  if(any(magnitude == 0) || qr(sweep(B, 2L, magnitude, `/`))$rank < ns) return(NULL)
  list(matrix = B[match(source_keys, source_levels), , drop = FALSE], lift = L)
}

.bt_formula_unscale_missing_terms <- function(spec, prefix, columns){

  terms <- .bt_formula_unscale_term_labels(spec, columns)
  .bt_formula_transform_stop(
    paste0(
      "Cannot transform the coefficients of formula parameter '", prefix,
      "' to the original predictor scale: centering the standardized predictors in term",
      if(length(terms) > 1L) "s " else " ", paste0("'", terms, "'", collapse = ", "),
      " induces lower-order effects that the formula does not contain. ",
      "Include the corresponding lower-order terms or keep the coefficients on ",
      "the fitted, standardized scale."
    ),
    parameter = prefix, reason = "original_scale_not_representable", terms = terms
  )
}

# Solve only the categorical representation. Numeric affine coefficients are
# applied afterwards, so their size cannot hide a missing lower-order effect.
.bt_formula_unscale_categorical_coefficients <- function(C, columns, source,
                                                         categorical_sets, data,
                                                         categorical_variables,
                                                         spec, prefix, names){

  out <- numeric(length(columns))
  identical_columns <- which(vapply(columns, function(k) identical(C[, k], C[, source]), logical(1)))
  if(source %in% columns){
    out[match(source, columns)] <- 1
    return(out)
  }
  if(length(identical_columns) == 1L){
    out[identical_columns] <- 1
    return(out)
  }

  source_variables <- categorical_sets[[source]]
  keep <- which(vapply(categorical_sets[columns], function(x) all(x %in% source_variables), logical(1)))
  source_count <- .bt_formula_unscale_cell_count(data, source_variables)
  target <- C[, columns, drop = FALSE]
  lift <- NULL
  if(length(keep) == source_count && source_count > 0L &&
     length(unique(.bt_formula_unscale_cell_keys(data, source_variables))) == source_count){
    restricted <- target[, keep, drop = FALSE]
    magnitude <- apply(abs(restricted), 2L, max)
    if(all(magnitude > 0) && qr(sweep(restricted, 2L, magnitude, `/`))$rank == length(keep)){
      target <- restricted
      lift <- diag(length(columns))[, keep, drop = FALSE]
    }
  }
  if(is.null(lift)){
    partition <- .bt_formula_unscale_partition_basis(
      data, C, columns, source_variables, categorical_variables, spec$assign
    )
    if(!is.null(partition)){
      target <- partition$matrix
      lift <- partition$lift
    }
  }

  source_max <- max(abs(C[, source]))
  target_max <- apply(abs(target), 2L, max)
  target_normalized <- sweep(target, 2L, target_max, `/`)
  source_normalized <- C[, source] / source_max
  rows <- qr(t(target_normalized))$pivot[seq_len(ncol(target))]
  solved <- tryCatch(
    as.numeric(solve(target_normalized[rows, , drop = FALSE], source_normalized[rows])),
    error = function(e) .bt_formula_transform_stop(
      "Original coefficient categorical arithmetic is numerically unavailable.",
      parameter = prefix, source = names[source], reason = "nonfinite_transform"
    )
  )
  residual <- as.numeric(target_normalized %*% solved - source_normalized)
  bound <- as.numeric(1e3 * .Machine$double.eps *
    (abs(target_normalized) %*% abs(solved) + max(abs(source_normalized))))
  if(any(!is.finite(solved)) || any(!is.finite(residual))){
    .bt_formula_transform_stop(
      "Original coefficient categorical arithmetic is numerically unavailable.",
      parameter = prefix, source = names[source], reason = "nonfinite_transform"
    )
  }
  if(any(abs(residual) > bound)) .bt_formula_unscale_missing_terms(spec, prefix, source)
  b <- numeric(length(solved))
  nonzero <- which(solved != 0)
  if(length(nonzero)) b[nonzero] <- .bt_formula_coefficient_ratio(
    solved[nonzero], numerator = source_max, denominator = target_max[nonzero],
    source = names[source], target = names[columns]
  )
  if(is.null(lift)) return(b)
  # Each original column belongs to at most one virtual column. Copying avoids
  # both cancellation and rounded reconstruction of its certified zeros.
  for(k in seq_along(b)) out[lift[, k] == 1] <- b[k]
  out
}

# X_s = X_o A and b_o = A b_s. Expand numeric monomials structurally and solve
# their categorical functions with the actual persisted coding. No magnitude
# threshold removes a coefficient or accepts a missing induced component.
.bt_formula_unscale_design_transform <- function(spec, formula_scale, prefix){

  .bt_formula_unscale_design_spec_check(spec, prefix)
  data <- .bt_formula_unscale_design_data(spec, formula_scale, prefix)$standardized
  # Rebuilding before grouping retains the raw-column/assignment guard,
  # including the existing boundary for matrix-valued predictors.
  .bt_formula_unscale_model_matrix(spec, data, prefix)
  coefficient_names <- .bt_formula_unscale_coefficient_names(spec, prefix)
  is_scaled <- stats::setNames(paste0(prefix, "_", spec$continuous) %in% names(formula_scale), spec$continuous)
  is_logical <- vapply(spec$continuous, function(v) isTRUE(spec$continuous_logical[[v]]), logical(1))
  numeric_names <- spec$continuous[!is_logical | is_scaled]
  categorical_names <- setdiff(names(data), numeric_names)
  factor_table <- attr(stats::terms(spec$formula), "factors")
  variables <- lapply(spec$assign, function(a){
    if(a == 0L) character() else rownames(factor_table)[factor_table[, a] > 0L]
  })
  subsets <- lapply(variables, function(v) sort(intersect(v, numeric_names), method = "radix"))
  categorical_sets <- lapply(variables, function(v) setdiff(v, numeric_names))
  key <- function(v) paste0("{", paste(v, collapse = "\r"), "}")
  keys <- vapply(subsets, key, character(1))
  groups <- split(seq_along(keys), keys)
  for(v in numeric_names) data[[v]] <- rep(1, nrow(data))
  C <- .bt_formula_unscale_model_matrix(spec, data, prefix)
  if(any(!is.finite(C))) .bt_formula_transform_stop(
    "Original coefficient categorical arithmetic is numerically unavailable.",
    parameter = prefix, reason = "nonfinite_transform")

  # Required support is derived before evaluating any numeric coefficient.
  components <- lapply(seq_along(subsets), function(j){
    scaled <- subsets[[j]][is_scaled[subsets[[j]]]]
    droppable <- scaled[vapply(scaled, function(v) formula_scale[[paste0(prefix, "_", v)]]$mean != 0, logical(1))]
    drops <- list(character())
    for(size in seq_along(droppable)) drops <- c(drops, utils::combn(droppable, size, simplify = FALSE))
    lapply(drops, function(d) list(scaled = scaled, dropped = d, target = key(setdiff(subsets[[j]], d))))
  })
  parsed <- lapply(coefficient_names, .parse_unscale_term_structure,
    prefix = prefix, scaled_vars = .formula_scale_strip_prefix(names(formula_scale), prefix))
  partners <- lapply(seq_along(components), function(j){
    lapply(components[[j]], function(component){
      if(!length(component$dropped)) return(j)
      possible <- groups[[component$target]]
      possible[vapply(possible, function(k){
        .unscale_ids_match(parsed[[k]]$unscaled_id, parsed[[j]]$unscaled_id) &&
          identical(C[, k], C[, j])
      }, logical(1))]
    })
  })
  early <- all(vapply(partners, function(x) all(lengths(x) == 1L), logical(1)))
  if(!early){
    # Rank precedes every missing-subset/residual refusal, so a random block
    # cannot route an unidentified design through its dropping convention.
    for(columns in groups){
      magnitude <- apply(abs(C[, columns, drop = FALSE]), 2L, max)
      if(any(!is.finite(magnitude))) .bt_formula_transform_stop(
        "Original coefficient categorical arithmetic is numerically unavailable.",
        parameter = prefix, reason = "nonfinite_transform")
      if(any(magnitude == 0) || qr(sweep(C[, columns, drop = FALSE], 2L, magnitude, `/`))$rank < length(columns)){
        .bt_formula_transform_stop(
          paste0("Cannot transform the coefficients of formula parameter '", prefix,
            "' to the original predictor scale: the fitted formula terms are linearly ",
            "dependent, so the original-scale coefficients are not uniquely determined."),
          parameter = prefix, reason = "original_scale_not_identified"
        )
      }
    }
  }
  A <- matrix(0, length(keys), length(keys), dimnames = list(coefficient_names, coefficient_names))
  for(j in seq_along(components)){
    for(i in seq_along(components[[j]])){
      component <- components[[j]][[i]]
      columns <- groups[[component$target]]
      if(!length(columns)) .bt_formula_unscale_missing_terms(spec, prefix, j)
      g <- .bt_formula_unscale_affine_coefficient(component$scaled, component$dropped,
        formula_scale, prefix, coefficient_names[j], coefficient_names[columns])
      if(early){
        A[partners[[j]][[i]], j] <- g
      }else{
        b <- .bt_formula_unscale_categorical_coefficients(C, columns, j,
          categorical_sets, data, categorical_names, spec, prefix, coefficient_names)
        A[columns, j] <- .bt_formula_checked_product(g, b, parameter = prefix,
          source = coefficient_names[j], target = coefficient_names[columns], transform = TRUE)
      }
    }
  }
  list(matrix = A, method = if(early) "name_paired" else "categorical_basis",
    coefficient_names = coefficient_names)
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
# @param formula_design Optional fitted formula designs keyed by parameter name;
#   the random-effect SDs of a parameter are unscaled with the design of its
#   blocks (see .apply_random_sd_unscale()).
# @return Transformed posterior matrix
.apply_unscale_transform <- function(posterior, formula_scale, targets = NULL,
                                     formula_design = NULL) {

  if (is.null(formula_scale) || length(formula_scale) == 0) {
    return(posterior)
  }

  # Every prefix reads the same original fitted snapshot, including states
  # whose displayed coefficient is changed by another prefix.
  fitted_state <- posterior
  for(parameter in names(formula_scale)){
    fitted_state <- .bt_formula_materialize_state(fitted_state,
      .bt_formula_state_constants(formula_scale[[parameter]]))
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
    posterior <- .apply_unscale_transform_single(posterior, param_scale, prefix = param_name,
      fitted_state = fitted_state, targets = targets,
      formula_design = formula_design[[param_name]])
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
# @param formula_scale Scaling info of one formula parameter, keyed by its
#   prefixed predictors: list(mu_x1 = list(mean, sd), mu_x2 = list(mean, sd)),
#   with the structure attributes of the fitted formula scale
# @param prefix Parameter prefix (e.g., "mu")
# @return Transformed posterior matrix
.apply_unscale_transform_single <- function(posterior, formula_scale, prefix,
                                             fitted_state = posterior, targets = NULL,
                                             formula_design = NULL) {

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
    spec <- attr(formula_scale, "unscale_design", exact = TRUE)
    source_names <- union(names(spec$multipliers), fixed_cols)
    transform <- .bt_formula_coefficient_transform(
      source_names = source_names,
      formula_scale = formula_scale,
      parameter = prefix
    )
    transformed <- .bt_apply_formula_coefficient_transform(
      fitted_state, transform, targets = if(is.null(targets)) transform$target_names else intersect(targets, transform$target_names)
    )
    selected <- if(is.null(targets)) transform$target_names else intersect(targets, transform$target_names)
    absent <- setdiff(selected, colnames(posterior))
    if(length(absent)) posterior <- cbind(posterior, transformed[, absent, drop = FALSE])
    posterior[, selected] <- transformed[, selected, drop = FALSE]
  }

  # random-effect SDs are transformed with the fitted random-effect structure,
  # which standardization information without the fitted design lacks
  if(length(random_sd_cols) > 0L){
    .bt_formula_scale_require_design(formula_scale, prefix)
  }
  posterior <- .apply_random_sd_unscale(posterior, random_sd_cols, formula_scale, prefix,
    formula_design = formula_design)

  return(posterior)
}
