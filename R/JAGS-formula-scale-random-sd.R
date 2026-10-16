.materialize_formula_scale_point_terms <- function(posterior, formula_scale,
                                                   prefix){

  point_terms <- attr(formula_scale, "point_terms", exact = TRUE)
  if(is.null(point_terms) || length(point_terms) == 0L){
    return(posterior)
  }

  point_terms <- point_terms[startsWith(names(point_terms), paste0(prefix, "_"))]
  point_terms <- point_terms[!names(point_terms) %in% colnames(posterior)]
  # Fits created before factor point terms were stored per coefficient carry
  # one unindexed name for a factor term; its indexed coefficients are present.
  indexed_present <- vapply(names(point_terms), function(name){
    !grepl("\\]$", name) &&
      any(startsWith(colnames(posterior), paste0(name, "[")))
  }, logical(1))
  point_terms <- point_terms[!indexed_present]
  if(length(point_terms) == 0L){
    return(posterior)
  }

  point_matrix <- matrix(
    rep(unname(point_terms), each = nrow(posterior)),
    nrow = nrow(posterior),
    dimnames = list(NULL, names(point_terms))
  )

  cbind(posterior, point_matrix)
}

# Original-scale random-effect SDs (and correlations) of the scaled blocks.
#
# @param formula_design the fitted formula design of the parameter (or NULL):
#   when it holds the block that owns an SD group, the block's own design
#   defines the map from its fitted-scale to its original-scale coefficients
#   (see .random_sd_unscale_matrix()); without it, hand-built SD-leaf
#   metadata are transformed by pairing coefficient names.
.apply_random_sd_unscale <- function(posterior, random_sd_cols, formula_scale,
                                     prefix,
                                     correlation_required_groups = NULL,
                                     formula_design = NULL){

  if(length(random_sd_cols) == 0){
    return(posterior)
  }

  column_groups <- .random_sd_column_unscale_groups(
    random_sd_cols = random_sd_cols,
    formula_scale = formula_scale,
    prefix = prefix
  )
  if(!is.null(column_groups)){
    return(.apply_random_sd_column_unscale(
      posterior = posterior,
      column_groups = column_groups,
      formula_scale = formula_scale,
      prefix = prefix,
      correlation_required_groups = correlation_required_groups,
      formula_design = formula_design
    ))
  }

  # Without the SD-leaf structure that the fitted formula scale carries there is
  # no random-effect SD to unscale: the name-based term map needs the same
  # structure, so it has no further route.
  posterior
}

.apply_random_sd_column_unscale <- function(posterior, column_groups,
                                            formula_scale, prefix,
                                            correlation_required_groups = NULL,
                                            formula_design = NULL){

  if(length(column_groups) == 0L){
    return(posterior)
  }

  required_groups <- .random_sd_correlation_required_groups(
    formula_scale = formula_scale,
    prefix = prefix,
    correlation_required_groups = correlation_required_groups
  )

  for(group in column_groups){
    group_key <- group$group_key
    M <- .random_sd_unscale_matrix(
      column_terms = group$column_terms,
      formula_scale = formula_scale,
      prefix = prefix,
      random_term = .random_sd_design_block(formula_design, group_key, group$leaf_names),
      source_data = formula_design$source_data
    )

    source_sd <- posterior[, group$leaf_names_by_column, drop = FALSE]
    source_cor <- .random_sd_correlation_draws(
      posterior = posterior,
      prefix = prefix,
      group_key = group_key,
      n_terms = length(group$column_terms),
      required = group_key %in% required_groups
    )
    transformed_sd_by_column <- matrix(
      NA_real_,
      nrow = nrow(source_sd),
      ncol = ncol(source_sd),
      dimnames = list(NULL, group$leaf_names_by_column)
    )

    if(is.null(source_cor)){
      for(target_i in seq_along(group$column_terms)){
        transformed_var <- rowSums(t(t(source_sd^2) * (M[target_i, ]^2)))
        transformed_sd_by_column[, target_i] <- sqrt(transformed_var)
      }
      .random_sd_check_range(source_sd, M, transformed_sd_by_column, group_key)
    }else{
      transformed_cor <- array(NA_real_, dim = dim(source_cor))
      for(draw_i in seq_len(nrow(source_sd))){
        source_cov <- diag(source_sd[draw_i, ], nrow = ncol(source_sd)) %*%
          source_cor[draw_i, , ] %*%
          diag(source_sd[draw_i, ], nrow = ncol(source_sd))
        transformed_cov <- M %*% source_cov %*% t(M)
        transformed_sd_by_column[draw_i, ] <- sqrt(diag(transformed_cov))
        transformed_cor[draw_i, , ] <- .random_sd_transformed_correlation(
          covariance = transformed_cov,
          sd = transformed_sd_by_column[draw_i, ],
          group_key = group_key
        )
      }
      .random_sd_check_range(source_sd, M, transformed_sd_by_column, group_key,
                             check_zero = FALSE)
      posterior <- .random_sd_assign_transformed_correlation(
        posterior = posterior,
        prefix = prefix,
        group_key = group_key,
        correlation = transformed_cor,
        nonzero_source_sd = .random_sd_nonzero_source_sd(source_sd)
      )
    }

    for(sd_col in group$leaf_names){
      column_index <- which(group$leaf_names_by_column == sd_col)
      if(length(column_index) == 1L){
        posterior[, sd_col] <- transformed_sd_by_column[, column_index]
        next
      }
      if(!.random_sd_shared_leaf_is_invariant(
        M = M,
        leaf_names_by_column = group$leaf_names_by_column,
        target_columns = column_index,
        correlated = !is.null(source_cor)
      )){
        .bt_formula_transform_stop(
          paste0(
            "Random-effect SD column '", sd_col,
            "' cannot be unscaled because its shared design columns transform to different SDs."
          ),
          parameter = prefix,
          reason = "original_scale_not_representable",
          coordinates = sd_col
        )
      }
      posterior[, sd_col] <- transformed_sd_by_column[, column_index[1L]]
    }
  }

  posterior
}

# The fitted random-effect block of 'formula_design' that owns the SD group:
# the block named by the group's key (the block name the SD coordinates are
# keyed by) whose SD coordinates include the group's leaves. NULL when there is
# no fitted design or it holds no such block (hand-built SD-leaf metadata).
.random_sd_design_block <- function(formula_design, group_key, leaf_names){

  blocks <- .bt_formula_design_random_effects(formula_design)
  if(length(blocks) == 0L || length(group_key) != 1L || is.na(group_key)){
    return(NULL)
  }

  owner <- vapply(blocks, function(block){
    identical(block$block_name, group_key) &&
      all(leaf_names %in% block$sd_parameter_names)
  }, logical(1))
  if(sum(owner) != 1L){
    return(NULL)
  }

  blocks[[which(owner)]]
}

# The part of a fitted formula design that the original-scale random-effect SDs
# need (.random_sd_design_block(), .bt_random_effect_unscale_design_spec()):
# per block its name, SD coordinates and the design to rebuild (formula without
# its environment, factor levels, contrasts, column layout), and the zero-row
# typed fitted data that identify logical predictors. Formula states keep it
# with the samples instead of the full design (design matrices, group maps).
# NULL when the design has no random-effect block.
.bt_formula_design_random_unscale_view <- function(design){

  blocks <- .bt_formula_design_random_effects(design)
  if(length(blocks) == 0L){
    return(NULL)
  }

  fields <- c("block_name", "sd_parameter_names", "term_formula", "xlevels",
              "contrast_matrices", "model_terms", "assign", "raw_column_names")
  blocks <- lapply(blocks, function(block){
    view <- block[intersect(fields, names(block))]
    if(inherits(view$term_formula, "formula")){
      environment(view$term_formula) <- emptyenv()
    }
    view
  })
  source_data <- design$source_data
  if(is.data.frame(source_data)){
    variables <- unique(unlist(lapply(blocks, function(block){
      if(inherits(block$term_formula, "formula")) all.vars(block$term_formula)
    }), use.names = FALSE))
    source_data <- source_data[0L, intersect(variables, colnames(source_data)), drop = FALSE]
  }else{
    source_data <- NULL
  }

  structure(list(random_effects = blocks, source_data = source_data), class = class(design))
}

# The views of the formula designs of a fit, by formula parameter, for those
# with a random-effect block; NULL when none has one.
.bt_formula_designs_random_unscale_view <- function(designs){

  views <- lapply(designs, .bt_formula_design_random_unscale_view)
  views <- views[!vapply(views, is.null, logical(1))]
  if(length(views) == 0L){
    return(NULL)
  }

  views
}

# Whether pairing coefficient names gives the exact original-scale map of a
# block that is known by its column terms only (no fitted design): every column
# with scaled components needs, for each subset of them whose means do not
# multiply to zero, a column with the same unscaled part and that subset (the
# intercept for the empty subset). Without such a partner, as for an
# interaction `f:x` of a factor with a scaled slope whose margin `x` is not in
# the block, centering induces terms the pairing drops (a fixed-effect map is
# derived from the fitted design instead). Two scaled columns that share one term
# name (the indicator columns of a logical predictor's interaction with a
# scaled slope) have different partners that one name cannot tell apart, so
# such a block is not exact either.
.random_sd_name_map_is_exact <- function(column_terms, formula_scale, prefix){

  scaled_vars <- .formula_scale_strip_prefix(names(formula_scale), prefix)
  structure <- lapply(
    paste0(prefix, "_", column_terms),
    .parse_unscale_term_structure,
    prefix = prefix,
    scaled_vars = scaled_vars
  )
  term_scaled <- lapply(structure, `[[`, "scaled")
  term_unscaled <- vapply(structure, `[[`, character(1), "unscaled_id")

  if(anyDuplicated(column_terms[lengths(term_scaled) > 0L]) > 0L){
    return(FALSE)
  }

  for(source_i in seq_along(structure)){
    source_scaled <- unique(term_scaled[[source_i]])
    for(subset_size in seq_len(length(source_scaled)) - 1L){
      target_sets <- if(subset_size == 0L){
        list(character())
      }else{
        utils::combn(source_scaled, subset_size, simplify = FALSE)
      }
      for(target_scaled in target_sets){
        extra_scaled <- setdiff(source_scaled, target_scaled)
        mean_product <- prod(vapply(paste0(prefix, "_", extra_scaled), function(parameter){
          formula_scale[[parameter]][["mean"]]
        }, numeric(1)))
        if(mean_product == 0){
          next
        }
        target_exists <- any(vapply(seq_along(structure), function(target_i){
          .unscale_ids_match(term_unscaled[[target_i]], term_unscaled[[source_i]]) &&
            length(term_scaled[[target_i]]) == length(target_scaled) &&
            all(term_scaled[[target_i]] %in% target_scaled)
        }, logical(1)))
        if(!target_exists){
          return(FALSE)
        }
      }
    }
  }

  TRUE
}

# The map A from the fitted-scale coefficients of a scaled random-effect block
# to its original-scale coefficients (rows and columns in block column order);
# the covariance of the block transforms as A Sigma A'. With the fitted block
# it is derived from the block's design exactly as the fixed-effect map is
# (.bt_formula_unscale_design_transform()), so a factor that a scaled slope
# interacts with is carried by the design and not by coefficient names. A
# block whose scaled design has no original-scale representation in its own
# columns (a term that would have to be induced, e.g. the implied intercept of
# `(0 + x | g)`) keeps the documented name-paired map, which drops the induced
# terms; any other failure propagates. A block whose design is linearly
# dependent has no unique original-scale map (`original_scale_not_identified`):
# a request for original-scale values stops with that condition
# (unidentified = "stop"), while a caller that classifies the block without
# transforming draws (the parameter map) asks for NULL (unidentified =
# "unavailable"). Without the fitted block the map pairs coefficient names. The
# returned matrix is the one map of the block: the draws, the prior law of an SD
# and the catalog read the same retained entries.
.random_sd_unscale_matrix <- function(column_terms, formula_scale, prefix,
                                      random_term = NULL, source_data = NULL,
                                      unidentified = c("stop", "unavailable")){

  unidentified <- match.arg(unidentified)

  if(!is.null(random_term)){
    spec <- .bt_random_effect_unscale_design_spec(random_term, source_data)
    if(!is.null(spec)){
      transform <- tryCatch(
        .bt_formula_unscale_design_transform(spec, formula_scale, prefix),
        BayesTools_formula_transform_unavailable = function(condition){
          if(identical(condition$reason, "original_scale_not_identified") &&
             identical(unidentified, "unavailable")){
            return("unidentified")
          }
          if(!identical(condition$reason, "original_scale_not_representable")){
            stop(condition)
          }
          NULL
        }
      )
      if(identical(transform, "unidentified")){
        return(NULL)
      }
      if(!is.null(transform)){
        if(nrow(transform$matrix) != length(column_terms)){
          .bt_stop_refit_required(
            "Random-effect design metadata for block '", random_term$block_name,
            "' disagree with its SD columns. Refit the model with this version of BayesTools."
          )
        }
        return(unname(transform$matrix))
      }
    }
  }

  .build_unscale_matrix_by_names(
    paste0(prefix, "_", column_terms),
    formula_scale,
    prefix
  )
}

# Whether a scaled block keeps the documented name-paired map that drops the
# lower-order terms of the scaled slope (see .random_sd_unscale_matrix()): its
# design is refused as `original_scale_not_representable`, e.g. the level
# intercepts that `x` induces in `diag(1 + f:x | g)`, which no column of the
# block carries.
.random_sd_block_drops_induced_terms <- function(formula_scale, prefix,
                                                 random_term, source_data = NULL){

  spec <- .bt_random_effect_unscale_design_spec(random_term, source_data)
  if(is.null(spec)){
    return(FALSE)
  }

  tryCatch(
    {
      .bt_formula_unscale_design_transform(spec, formula_scale, prefix)
      FALSE
    },
    BayesTools_formula_transform_unavailable = function(condition){
      identical(condition$reason, "original_scale_not_representable")
    }
  )
}

# The unscale-design specification of a fitted random-effect block, in the
# format of .bt_formula_unscale_design_spec(): the block's formula (with its
# own intercept control), factor levels and concrete contrasts, and its fitted
# column layout. The variables of the block that are not factors are the
# continuous ones; those named in 'formula_scale' are the scaled ones, and
# logical predictors are those that are logical in the fitted source data.
# NULL for a block record without a design to rebuild.
.bt_random_effect_unscale_design_spec <- function(random_term, source_data = NULL){

  if(!inherits(random_term$term_formula, "formula") ||
     is.null(random_term$model_terms) || is.null(random_term$assign) ||
     is.null(random_term$raw_column_names)){
    return(NULL)
  }

  variables <- all.vars(random_term$term_formula)
  factors <- intersect(names(random_term$xlevels), variables)
  factor_levels <- lapply(factors, function(factor_name){
    as.character(random_term$xlevels[[factor_name]])
  })
  names(factor_levels) <- factors
  factor_ordered <- vapply(factors, function(factor_name) FALSE, logical(1))
  names(factor_ordered) <- factors
  contrast_matrices <- random_term$contrast_matrices[factors]
  names(contrast_matrices) <- factors
  continuous <- setdiff(variables, factors)
  continuous_logical <- vapply(continuous, function(predictor){
    is.data.frame(source_data) && is.logical(source_data[[predictor]])
  }, logical(1))
  names(continuous_logical) <- continuous
  formula <- random_term$term_formula
  environment(formula) <- emptyenv()

  list(
    schema_version     = .bt_formula_unscale_design_spec_version,
    formula            = formula,
    continuous         = continuous,
    continuous_logical = continuous_logical,
    factor_levels      = factor_levels,
    factor_ordered     = factor_ordered,
    contrast_matrices  = contrast_matrices,
    model_terms        = as.character(random_term$model_terms),
    assign             = as.integer(random_term$assign),
    raw_column_names   = as.character(random_term$raw_column_names),
    multipliers        = stats::setNames(list(), character()),
    state_constants    = stats::setNames(numeric(), character()),
    owner_scope        = "compiler",
    complete           = FALSE
  )
}

# The one zero rule of the original-scale map of a random-effect block (the
# unscale matrix): an entry is zero exactly when the matrix that the draws use
# (.random_sd_unscale_matrix()) holds an exact zero there, so the draws and
# every classification of the map read the same zeros. Every retained nonzero
# entry, however small beside the largest of its row, enters the draws and so
# counts here. The rows of the matrix decide
# which fitted SDs an original-scale SD reads (the catalog support and its
# exactness, the prior law of a correlation) and whether it is a pure rescaling
# (.random_sd_original_scale_laws()). Returns a logical vector, all FALSE for an
# all-zero row.
.random_sd_unscale_row_nonzero <- function(row){

  row != 0
}

# The relation of the original-scale random-effect SDs of the scaled blocks of
# one formula parameter to their fitted SDs, for the prior law of an SD (see
# .prior_density_context_random_sd_scale()). The original-scale SD of a design
# column is a pure rescaling of its fitted SD when its row of the block's
# original-scale map (.random_sd_unscale_matrix(), the matrix the draws use)
# has exactly one nonzero entry, on its own coordinate; that entry c is the
# scale, and the original-scale law is the fitted SD's law rescaled by |c|.
# An entry counts as zero by the one rule of .random_sd_unscale_row_nonzero():
# exactly the zeros of the matrix the draws use. A retained nonzero entry of
# any size makes a composite, including one induced by a numerically nonzero
# mean after centering by subtraction. Small predictor units alone can make
# such an entry consequential, so no additional magnitude threshold is safe. A
# shared leaf is a pure rescaling only when every column of the leaf is, with
# the same |c|. Any other SD (e.g. an intercept SD combined with the SD of a
# non-centered scaled slope) is a composite of several fitted SDs.
#
# Returns the SD coordinates of the blocks with a scaled column: 'scale', the
# |c| of the pure ones (a named numeric), 'refused', the cause of the others (a
# named character): "random_effect_design_unavailable" for a block without a
# fitted design in 'design' (the part of the formula design of the parameter
# that .bt_formula_design_random_unscale_view() keeps), and
# "random_sd_prior_law_unavailable" for a composite, and
# "original_scale_not_identified" for a block whose linearly dependent design
# has no unique original-scale map, and 'blocks', the block
# name of every one of these coordinates. The SDs of blocks without a scaled
# column, and of structured blocks, are in none of them: their fitted law is
# the original-scale law.
.random_sd_original_scale_laws <- function(parameter_scale, prefix, design = NULL){

  # every SD of a block asks for the laws of its block: the session memo keeps
  # the laws of the inputs it recognises as identical()
  .bt_content_memo(
    "random_sd_original_scale_laws",
    list(parameter_scale, prefix, design),
    function() .random_sd_original_scale_laws_build(parameter_scale, prefix, design)
  )
}

.random_sd_original_scale_laws_build <- function(parameter_scale, prefix, design){

  scale <- stats::setNames(numeric(), character())
  refused <- blocks <- stats::setNames(character(), character())
  leaves <- attr(parameter_scale, "random_effect_sd_leaves", exact = TRUE)
  if(is.null(leaves) || length(leaves) == 0L){
    return(list(scale = scale, refused = refused, blocks = blocks))
  }

  leaf_names <- unique(unlist(lapply(leaves, function(block){
    as.character(block$leaf_names)
  }), use.names = FALSE))
  groups <- .random_sd_column_unscale_groups(leaf_names, parameter_scale, prefix)
  for(group in groups){
    blocks[group$leaf_names] <- group$group_key
    block <- .random_sd_design_block(design, group$group_key, group$leaf_names)
    if(is.null(block) ||
       is.null(.bt_random_effect_unscale_design_spec(block, design$source_data))){
      refused[group$leaf_names] <- "random_effect_design_unavailable"
      next
    }

    A <- .random_sd_unscale_matrix(
      column_terms = group$column_terms,
      formula_scale = parameter_scale,
      prefix = prefix,
      random_term = block,
      source_data = design$source_data,
      unidentified = "unavailable"
    )
    if(is.null(A)){
      # a linearly dependent block design: the SDs are composites without an
      # original-scale law, refused when original-scale values are requested
      refused[group$leaf_names] <- "original_scale_not_identified"
      next
    }
    column_scale <- vapply(seq_len(nrow(A)), function(row){
      if(identical(which(.random_sd_unscale_row_nonzero(A[row, ])), row)){
        abs(A[row, row])
      }else{
        NA_real_
      }
    }, numeric(1))
    for(leaf in unique(group$leaf_names_by_column)){
      leaf_scale <- column_scale[group$leaf_names_by_column == leaf]
      if(anyNA(leaf_scale) || !isTRUE(all.equal(
        leaf_scale, rep(leaf_scale[1L], length(leaf_scale)),
        tolerance = 1e-12, scale = max(leaf_scale)
      ))){
        refused[leaf] <- "random_sd_prior_law_unavailable"
      }else{
        scale[leaf] <- leaf_scale[1L]
      }
    }
  }

  list(scale = scale, refused = refused, blocks = blocks)
}

.random_sd_shared_leaf_is_invariant <- function(
    M, leaf_names_by_column, target_columns, correlated){

  if(isTRUE(correlated)){
    return(FALSE)
  }
  leaf_names <- unique(leaf_names_by_column)
  variance_weights <- vapply(leaf_names, function(leaf_name){
    rowSums(M[, leaf_names_by_column == leaf_name, drop = FALSE]^2)
  }, numeric(nrow(M)))
  variance_weights <- matrix(
    variance_weights,
    nrow = nrow(M),
    dimnames = list(NULL, leaf_names)
  )
  # The weights come from a least-squares map, so equal weights agree only up to
  # rounding of a few ulps relative to their magnitude. The comparison is
  # relative to the larger weight, whatever its size: all.equal() compares
  # absolutely when the differing elements have a mean magnitude below the
  # tolerance, which merges weights of order 1 / sd^2 of a predictor with a very
  # large SD (all-zero vectors are equal before the scale is used).
  target_weights <- variance_weights[target_columns, , drop = FALSE]
  first_weights <- unname(target_weights[1L, , drop = TRUE])
  all(vapply(seq_len(nrow(target_weights)), function(i){
    weights <- unname(target_weights[i, , drop = TRUE])
    isTRUE(all.equal(
      weights, first_weights,
      tolerance = 1e-12,
      scale = max(abs(c(weights, first_weights)))
    ))
  }, logical(1)))
}

# Original-scale correlation matrix of one draw from its transformed
# covariance. A correlation is defined whenever both SDs are positive,
# including a singular covariance (perfect correlation); entries of a
# coefficient with a zero SD are missing. Rounding can push a perfect
# correlation just beyond +-1, which is set to +-1; a larger excess signals an
# inconsistent covariance.
.random_sd_transformed_correlation <- function(covariance, sd, group_key,
                                               tolerance = 1e-8){

  n_terms <- length(sd)
  correlation <- matrix(NA_real_, nrow = n_terms, ncol = n_terms)
  positive <- !is.na(sd) & is.finite(sd) & sd > 0
  if(!any(positive)){
    return(correlation)
  }

  correlation[positive, positive] <- covariance[positive, positive, drop = FALSE] /
    tcrossprod(sd[positive])
  diag(correlation)[positive] <- 1
  excess <- abs(correlation) - 1
  if(any(excess > tolerance, na.rm = TRUE)){
    stop(
      "Internal error: an original-scale random-effect correlation of block '",
      group_key, "' exceeds 1 in absolute value by ",
      signif(max(excess, na.rm = TRUE), 3), ".",
      call. = FALSE
    )
  }
  rounded <- !is.na(excess) & excess > 0
  correlation[rounded] <- sign(correlation[rounded])

  correlation
}

# Double-precision range edge. The original-scale SDs are formed from the
# squares of the fitted-scale SDs, so finite positive SDs beyond about 1e154
# overflow to Inf and SDs below about 1e-154 underflow to an exact zero. Such a
# draw is refused instead of returning an infinite or zero SD; no rescaling is
# attempted. A zero cannot be told from an exact cancellation under a
# correlation, so only the uncorrelated transform checks for it
# (`check_zero`).
.random_sd_check_range <- function(source_sd, M, transformed_sd, group_key,
                                   check_zero = TRUE){

  finite_source <- rowSums(!is.finite(source_sd)) == 0L
  lost <- rowSums(!is.finite(transformed_sd)) > 0L
  if(isTRUE(check_zero)){
    declared <- (source_sd > 0) %*% t(M != 0) > 0
    declared[is.na(declared)] <- FALSE
    lost <- lost | rowSums(declared & !is.na(transformed_sd) &
                             transformed_sd == 0) > 0L
  }
  lost <- finite_source & lost
  if(any(lost)){
    .bt_formula_transform_stop(
      paste0(
        "Original-scale random-effect SDs of block '", group_key,
        "' are unavailable: the transformed SD of ", sum(lost), " draw",
        if(sum(lost) > 1L) "s" else "",
        " is not representable in double precision although the fitted-scale ",
        "SDs are finite. Rescale the response or predictors so that the ",
        "random-effect SDs are well within 1e-150 to 1e150."
      ),
      reason = "nonfinite_transform", block = group_key,
      draws = which(lost)
    )
  }

  invisible(TRUE)
}

# A zero source SD makes the original-scale covariance M S M' singular. The
# transformed correlation cannot decide this reliably: a perfect correlation
# can round to just inside +-1, where chol() succeeds.
.random_sd_nonzero_source_sd <- function(source_sd){

  rowSums(!is.finite(source_sd) | source_sd == 0) == 0L
}

.random_sd_assign_transformed_correlation <- function(posterior, prefix,
                                                      group_key,
                                                      correlation,
                                                      nonzero_source_sd){

  n_terms <- dim(correlation)[2L]

  R_names <- .random_sd_correlation_matrix_names(
    prefix = prefix,
    group_key = group_key,
    suffix = "_xRE_CORx_R",
    n_terms = n_terms
  )
  L_names <- .random_sd_correlation_matrix_names(
    prefix = prefix,
    group_key = group_key,
    suffix = "_xRE_CORx_L",
    n_terms = n_terms
  )

  u_names <- .random_sd_lkj_primitive_names(
    prefix = prefix,
    group_key = group_key,
    n_terms = n_terms
  )
  cpc_names <- sub("_lkj_u[", "_lkj_cpc[", u_names, fixed = TRUE)

  has_R <- all(as.vector(R_names) %in% colnames(posterior))
  has_L <- all(as.vector(L_names) %in% colnames(posterior))
  # Monitored LKJ primitives define the Cholesky factor preferred by the
  # reconstruction evaluators; they (and their monitored canonical partial
  # correlations, cpc = 2u - 1) are rewritten from the transformed factor.
  has_u <- length(u_names) > 0L && all(u_names %in% colnames(posterior))
  has_cpc <- length(cpc_names) > 0L && all(cpc_names %in% colnames(posterior))
  if(!has_R && !has_L && !has_u && !has_cpc){
    return(posterior)
  }

  # Correlation entries are defined whenever both SDs are positive, singular
  # draws included.
  if(has_R){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[, R_names[row, column]] <- correlation[, row, column]
      }
    }
  }
  if(!has_L && !has_u && !has_cpc){
    return(posterior)
  }

  # The Cholesky factor and the LKJ primitives parameterize positive-definite
  # correlation matrices. Draws with a zero SD, or with a singular (perfect)
  # correlation from a zero source SD, leave them missing.
  if(has_L){
    posterior[, as.vector(L_names)] <- NA_real_
  }
  if(has_u){
    posterior[, u_names] <- NA_real_
  }
  if(has_cpc){
    posterior[, cpc_names] <- NA_real_
  }
  valid_draw <- rep(FALSE, dim(correlation)[1L])
  L <- array(NA_real_, dim = dim(correlation))
  for(draw_i in seq_len(dim(correlation)[1L])){
    if(!nonzero_source_sd[draw_i] ||
       !all(is.finite(correlation[draw_i, , ]))){
      next
    }
    this_L <- try(t(chol(correlation[draw_i, , ])), silent = TRUE)
    if(inherits(this_L, "try-error")){
      next
    }
    L[draw_i, , ] <- this_L
    valid_draw[draw_i] <- TRUE
  }
  if(!any(valid_draw)){
    return(posterior)
  }

  if(has_L){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[valid_draw, L_names[row, column]] <-
          L[valid_draw, row, column]
      }
    }
  }
  if(has_u || has_cpc){
    u <- .bt_lkj_cholesky_L_to_cpc_u(
      L[valid_draw, , , drop = FALSE],
      K = n_terms
    )
    if(has_u){
      posterior[valid_draw, u_names] <- u
    }
    if(has_cpc){
      posterior[valid_draw, cpc_names] <- .bt_dnode_lkj_cpc(u)
    }
  }

  posterior
}

.random_sd_lkj_primitive_names <- function(prefix, group_key, n_terms){

  n_pairs <- n_terms * (n_terms - 1L) / 2L
  if(n_pairs < 1L){
    return(character())
  }

  paste0(prefix, "__xREx__", group_key, "_xRE_CORx_lkj_u[", seq_len(n_pairs), "]")
}

.random_sd_correlation_matrix_names <- function(prefix, group_key, suffix,
                                                n_terms){

  stem <- paste0(prefix, "__xREx__", group_key, suffix)
  outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "[", row, ",", column, "]")
    })
  )
}

.random_sd_column_unscale_groups <- function(random_sd_cols, formula_scale,
                                             prefix){

  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves", exact = TRUE)
  if(is.null(sd_leaves)){
    return(NULL)
  }

  scaled_vars <- .formula_scale_strip_prefix(names(formula_scale), prefix)
  groups <- list()
  leaf_keys <- names(sd_leaves)
  if(is.null(leaf_keys)){
    leaf_keys <- rep("", length(sd_leaves))
  }
  for(i in seq_along(sd_leaves)){
    group <- .random_sd_column_unscale_group(
      leaves = sd_leaves[[i]],
      leaf_key = leaf_keys[[i]],
      random_sd_cols = random_sd_cols,
      scaled_vars = scaled_vars,
      prefix = prefix
    )
    if(is.null(group)){
      next
    }
    groups[[length(groups) + 1L]] <- group
  }

  groups
}

.random_sd_column_unscale_group <- function(leaves, leaf_key, random_sd_cols,
                                            scaled_vars, prefix){

  if(.random_sd_structured_leaves(leaves)){
    return(NULL)
  }

  if(is.null(leaves$leaf_names_by_column) ||
     is.null(leaves$leaf_terms_by_column) ||
     is.null(leaves$leaf_names)){
    return(NULL)
  }

  leaf_names_by_column <- as.character(leaves$leaf_names_by_column)
  column_terms <- as.character(leaves$leaf_terms_by_column)
  leaf_names <- as.character(leaves$leaf_names)

  if(length(leaf_names_by_column) == 0L ||
     length(column_terms) != length(leaf_names_by_column)){
    return(NULL)
  }

  available_leaf_names <- leaf_names[leaf_names %in% random_sd_cols]
  if(length(available_leaf_names) == 0L){
    return(NULL)
  }

  homogeneous <- !is.null(leaves$leaf_terms) &&
    identical(unique(unname(as.character(leaves$leaf_terms))), "sd")
  if(.random_sd_homogeneous_scaled_block(leaves, scaled_vars)){
    .random_sd_homogeneous_scaled_stop(available_leaf_names, parameter = prefix)
  }
  if(isTRUE(homogeneous) && !is.null(leaves$column_names)){
    column_terms <- vapply(
      leaves$column_names,
      .random_sd_term_from_column_name,
      character(1)
    )
  }
  scaled_columns <- vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  )

  if(!any(scaled_columns)){
    return(NULL)
  }

  missing_leaf_names <- setdiff(leaf_names, random_sd_cols)
  group_key <- .random_sd_group_key_from_leaf_key(
    leaf_key = leaf_key,
    leaf_names = leaf_names,
    column_terms = column_terms,
    prefix = prefix
  )
  if(length(missing_leaf_names) > 0L){
    .bt_formula_transform_stop(
      paste0(
        "Random-effect SD unscaling for block '", group_key,
        "' requires all SD columns from the block. Missing: ",
        paste0("'", missing_leaf_names, "'", collapse = ", "),
        "."
      ),
      parameter   = prefix,
      reason      = "random_effect_block_incomplete",
      coordinates = missing_leaf_names
    )
  }

  list(
    group_key = group_key,
    leaf_names = leaf_names,
    leaf_names_by_column = leaf_names_by_column,
    column_terms = column_terms
  )
}

.random_sd_structured_leaves <- function(leaves){

  random_structure <- leaves$random_structure
  !is.null(random_structure) &&
    length(random_structure) == 1L &&
    random_structure %in% c("cs", "hcs", "ar1", "car", "har")
}

.random_sd_group_key_from_leaf_key <- function(leaf_key, leaf_names,
                                               column_terms, prefix){

  if(!is.null(leaf_key) && length(leaf_key) == 1L &&
     !is.na(leaf_key) && nzchar(leaf_key)){
    return(sub("^__xREx__", "", leaf_key))
  }

  term_map <- stats::setNames(column_terms[seq_along(leaf_names)], leaf_names)
  .random_sd_group_key(
    col = leaf_names[[1L]],
    term_map = term_map,
    prefix = prefix
  )
}

# A homogeneous block shares one SD leaf among all its columns. With a column
# of a scaled predictor the columns transform to different original-scale SDs,
# which the shared leaf cannot carry, so the block has no original-scale SD.
# Structured blocks (which use unscaled data) and blocks whose SD is an
# external or row source (their leaves are not the homogeneous 'sd' leaf) are
# not covered; neither is a single-column block.
.random_sd_homogeneous_scaled_block <- function(leaves, scaled_vars){

  if(.random_sd_structured_leaves(leaves) || is.null(leaves$leaf_terms) ||
     !identical(unique(unname(as.character(leaves$leaf_terms))), "sd")){
    return(FALSE)
  }
  column_terms <- if(!is.null(leaves$column_names)){
    vapply(leaves$column_names, .random_sd_term_from_column_name, character(1))
  }else{
    as.character(leaves$leaf_terms_by_column)
  }
  length(column_terms) > 1L && any(vapply(
    column_terms,
    .random_sd_term_uses_scaled_var,
    logical(1),
    scaled_vars = scaled_vars
  ))
}

# The refusal of JAGS_formula() and of the unscaling of hand-built metadata:
# BayesTools_formula_transform_unavailable, reason "homogeneous_scaled_block",
# coordinates the shared SD leaves. The remedy holds in both: a structure with
# one SD per column, or a slope that is not standardized.
.random_sd_homogeneous_scaled_stop <- function(leaf_names, parameter){

  leaf_names <- unique(leaf_names)
  .bt_formula_transform_stop(
    paste0(
      "Cannot unscale homogeneous random-effect SD '",
      paste(leaf_names, collapse = "', '"),
      "' because its block contains scaled random-slope columns, which one ",
      "shared SD cannot represent. Use a heterogeneous random-effect SD ",
      "structure (for example diag() or us()), or do not standardize the slope ",
      "with 'formula_scale'."
    ),
    parameter   = parameter,
    reason      = "homogeneous_scaled_block",
    coordinates = leaf_names
  )
}

.random_sd_term_from_column_name <- function(column_name){

  if(identical(column_name, "(Intercept)")){
    return("intercept")
  }

  column_name
}

.random_sd_term_uses_scaled_var <- function(term, scaled_vars){

  if(length(scaled_vars) == 0L || identical(term, "intercept")){
    return(FALSE)
  }

  term <- sub("(\\[[^]]+\\]|\\{[^}]+\\})$", "", term)
  components <- .bt_random_effect_term_components(term)
  any(components %in% scaled_vars)
}

.random_sd_group_key <- function(col, term_map, prefix){

  base_col <- sub("\\[[^]]+\\]$", "", col)
  rest <- .formula_scale_strip_prefix(base_col, prefix, "__xREx__")
  term <- unname(term_map[[col]])
  term_core <- gsub("\\[[^]]+\\]|\\{[^}]+\\}", "", term)
  .random_sd_strip_term_suffix(rest, term_core)
}

.random_sd_strip_term_suffix <- function(x, term){

  suffix <- paste0("_", term)
  if(endsWith(x, suffix)){
    return(substr(x, 1L, nchar(x) - nchar(suffix)))
  }

  x
}

.random_sd_correlation_required_groups <- function(formula_scale, prefix,
                                                   correlation_required_groups = NULL){

  required <- character()
  metadata <- attr(formula_scale, "random_effect_correlation_required")
  if(!is.null(metadata)){
    if(is.logical(metadata)){
      if(!is.null(names(metadata))){
        metadata <- names(metadata)[metadata]
      }else{
        metadata <- character()
      }
    }
    required <- c(required, as.character(metadata))
  }
  if(!is.null(correlation_required_groups)){
    required <- c(required, as.character(correlation_required_groups))
  }
  required <- required[!is.na(required) & nzchar(required)]
  if(length(required) == 0L){
    return(character())
  }

  required <- .formula_scale_strip_prefix(required, prefix, "__xREx__")
  required <- sub("^__xREx__", "", required)
  unique(required)
}

.random_sd_correlation_draws <- function(posterior, prefix, group_key, n_terms,
                                         required = FALSE){

  if(n_terms < 2L){
    return(NULL)
  }

  stem <- paste0(prefix, "__xREx__", group_key)
  R_names <- outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  R_present <- as.vector(R_names) %in% colnames(posterior)

  L_names <- outer(
    seq_len(n_terms),
    seq_len(n_terms),
    Vectorize(function(row, column){
      paste0(stem, "_xRE_CORx_L[", row, ",", column, "]")
    })
  )
  L_present <- as.vector(L_names) %in% colnames(posterior)

  if(any(R_present) && !all(R_present)){
    stop(
      "Random-effect correlation samples are incomplete for block '",
      group_key,
      "'. Expected a complete random-effect correlation matrix.",
      call. = FALSE
    )
  }
  if(any(L_present) && !all(L_present)){
    stop(
      "Random-effect Cholesky samples are incomplete for block '",
      group_key,
      "'. Expected a complete random-effect Cholesky matrix.",
      call. = FALSE
    )
  }

  if(all(R_present)){
    out <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        out[, row, column] <- posterior[, R_names[row, column]]
      }
    }
    return(out)
  }

  u_names <- .random_sd_lkj_primitive_names(
    prefix = prefix,
    group_key = group_key,
    n_terms = n_terms
  )
  u_present <- u_names %in% colnames(posterior)
  if(!any(L_present) && all(u_present)){
    # Only the LKJ primitives are available (e.g., prior draws or fits without
    # monitored correlation matrices): the correlation of the 'lkj' node.
    return(.bt_dnode_lkj_correlation(
      posterior[, u_names, drop = FALSE],
      K = n_terms
    ))
  }

  if(!all(L_present)){
    if(isTRUE(required)){
      stop(
        "Random-effect SD unscaling requires random-effect correlation samples for block '",
        group_key,
        "'. Expected monitored correlation or Cholesky coordinates.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  L <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
  for(row in seq_len(n_terms)){
    for(column in seq_len(n_terms)){
      L[, row, column] <- posterior[, L_names[row, column]]
    }
  }

  .bt_lkj_cholesky_L_to_R(L)
}
