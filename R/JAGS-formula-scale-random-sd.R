.materialize_formula_scale_point_terms <- function(posterior, formula_scale,
                                                   prefix){

  point_terms <- attr(formula_scale, "point_terms", exact = TRUE)
  if(is.null(point_terms) || length(point_terms) == 0L){
    return(posterior)
  }

  point_terms <- point_terms[startsWith(names(point_terms), paste0(prefix, "_"))]
  point_terms <- point_terms[!names(point_terms) %in% colnames(posterior)]
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

.apply_random_sd_unscale <- function(posterior, random_sd_cols, formula_scale,
                                     prefix,
                                     correlation_required_groups = NULL){

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
      correlation_required_groups = correlation_required_groups
    ))
  }

  term_map <- .random_sd_term_map(random_sd_cols, formula_scale, prefix)
  if(length(term_map) == 0){
    return(posterior)
  }

  random_sd_cols <- names(term_map)
  group_keys <- vapply(random_sd_cols, .random_sd_group_key,
                       character(1), term_map = term_map, prefix = prefix)
  required_groups <- .random_sd_correlation_required_groups(
    formula_scale = formula_scale,
    prefix = prefix,
    correlation_required_groups = correlation_required_groups
  )

  for(group_key in unique(group_keys)){
    group_cols <- random_sd_cols[group_keys == group_key]
    group_terms <- unname(term_map[group_cols])

    if(any(duplicated(group_terms))){
      duplicated_terms <- unique(group_terms[duplicated(group_terms)])
      stop(
        "Random-effect SD columns for group '", group_key,
        "' cannot be unscaled because multiple columns map to the same term: ",
        paste0("'", duplicated_terms, "'", collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    pseudo_terms <- paste0(prefix, "_", group_terms)
    names(pseudo_terms) <- group_cols
    M <- .build_unscale_matrix(
      unname(pseudo_terms),
      formula_scale,
      prefix,
      require_closure = FALSE
    )

    source_sd <- posterior[, group_cols, drop = FALSE]
    source_cor <- .random_sd_correlation_draws(
      posterior = posterior,
      prefix = prefix,
      group_key = group_key,
      n_terms = length(group_cols),
      required = group_key %in% required_groups
    )
    transformed_sd <- matrix(NA_real_, nrow = nrow(source_sd), ncol = ncol(source_sd))
    colnames(transformed_sd) <- group_cols

    if(is.null(source_cor)){
      for(target_i in seq_along(group_cols)){
        transformed_var <- rowSums(t(t(source_sd^2) * (M[target_i, ]^2)))
        transformed_sd[, target_i] <- sqrt(transformed_var)
      }
    }else{
      transformed_cor <- array(NA_real_, dim = dim(source_cor))
      valid_cor_draw <- rep(FALSE, nrow(source_sd))
      for(draw_i in seq_len(nrow(source_sd))){
        source_cov <- diag(source_sd[draw_i, ], nrow = ncol(source_sd)) %*%
          source_cor[draw_i, , ] %*%
          diag(source_sd[draw_i, ], nrow = ncol(source_sd))
        transformed_cov <- M %*% source_cov %*% t(M)
        transformed_sd[draw_i, ] <- sqrt(diag(transformed_cov))
        if(any(!is.finite(transformed_sd[draw_i, ]) | transformed_sd[draw_i, ] <= 0)){
          next
        }
        transformed_cor[draw_i, , ] <- transformed_cov /
          tcrossprod(transformed_sd[draw_i, ])
        diag(transformed_cor[draw_i, , ]) <- 1
        valid_cor_draw[draw_i] <- all(is.finite(transformed_cor[draw_i, , ]))
      }
      posterior <- .random_sd_assign_transformed_correlation(
        posterior = posterior,
        prefix = prefix,
        group_key = group_key,
        correlation = transformed_cor,
        valid_draw = valid_cor_draw
      )
    }

    posterior[, group_cols] <- transformed_sd
  }

  posterior
}

.apply_random_sd_column_unscale <- function(posterior, column_groups,
                                            formula_scale, prefix,
                                            correlation_required_groups = NULL){

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
    pseudo_terms <- paste0(prefix, "_", group$column_terms)
    M <- .build_unscale_matrix(
      pseudo_terms,
      formula_scale,
      prefix,
      require_closure = FALSE
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
    }else{
      transformed_cor <- array(NA_real_, dim = dim(source_cor))
      valid_cor_draw <- rep(FALSE, nrow(source_sd))
      for(draw_i in seq_len(nrow(source_sd))){
        source_cov <- diag(source_sd[draw_i, ], nrow = ncol(source_sd)) %*%
          source_cor[draw_i, , ] %*%
          diag(source_sd[draw_i, ], nrow = ncol(source_sd))
        transformed_cov <- M %*% source_cov %*% t(M)
        transformed_sd_by_column[draw_i, ] <- sqrt(diag(transformed_cov))
        if(any(!is.finite(transformed_sd_by_column[draw_i, ]) |
               transformed_sd_by_column[draw_i, ] <= 0)){
          next
        }
        transformed_cor[draw_i, , ] <- transformed_cov /
          tcrossprod(transformed_sd_by_column[draw_i, ])
        diag(transformed_cor[draw_i, , ]) <- 1
        valid_cor_draw[draw_i] <- all(is.finite(transformed_cor[draw_i, , ]))
      }
      posterior <- .random_sd_assign_transformed_correlation(
        posterior = posterior,
        prefix = prefix,
        group_key = group_key,
        correlation = transformed_cor,
        valid_draw = valid_cor_draw
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
        stop(
          "Random-effect SD column '", sd_col,
          "' cannot be unscaled because its shared design columns transform to different SDs.",
          call. = FALSE
        )
      }
      posterior[, sd_col] <- transformed_sd_by_column[, column_index[1L]]
    }
  }

  posterior
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
  target_weights <- variance_weights[target_columns, , drop = FALSE]
  all(vapply(seq_len(nrow(target_weights)), function(i){
    identical(
      unname(target_weights[i, , drop = TRUE]),
      unname(target_weights[1L, , drop = TRUE])
    )
  }, logical(1)))
}

.random_sd_assign_transformed_correlation <- function(posterior, prefix,
                                                      group_key,
                                                      correlation,
                                                      valid_draw = NULL){

  n_terms <- dim(correlation)[2L]
  if(is.null(valid_draw)){
    valid_draw <- rep(TRUE, dim(correlation)[1L])
  }
  valid_draw <- valid_draw & is.finite(vapply(
    seq_len(dim(correlation)[1L]),
    function(draw_i) sum(correlation[draw_i, , ]),
    numeric(1)
  ))

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

  has_R <- all(as.vector(R_names) %in% colnames(posterior))
  has_L <- all(as.vector(L_names) %in% colnames(posterior))
  if(!has_R && !has_L){
    return(posterior)
  }

  if(has_R){
    posterior[, as.vector(R_names)] <- NA_real_
  }
  if(has_L){
    posterior[, as.vector(L_names)] <- NA_real_
  }
  if(!any(valid_draw)){
    return(posterior)
  }

  L <- NULL
  if(has_L){
    L <- array(NA_real_, dim = dim(correlation))
    for(draw_i in which(valid_draw)){
      this_L <- try(t(chol(correlation[draw_i, , ])), silent = TRUE)
      if(inherits(this_L, "try-error")){
        valid_draw[draw_i] <- FALSE
        next
      }
      L[draw_i, , ] <- this_L
    }
    if(!any(valid_draw)){
      return(posterior)
    }
  }

  if(has_R){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[valid_draw, R_names[row, column]] <-
          correlation[valid_draw, row, column]
      }
    }
  }
  if(has_L){
    for(row in seq_len(n_terms)){
      for(column in seq_len(n_terms)){
        posterior[valid_draw, L_names[row, column]] <-
          L[valid_draw, row, column]
      }
    }
  }

  posterior
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

.random_sd_term_map <- function(random_sd_cols, formula_scale, prefix){

  scaled_vars <- .formula_scale_strip_prefix(names(formula_scale), prefix)
  sd_leaves <- attr(formula_scale, "random_effect_sd_leaves", exact = TRUE)
  if(!is.null(sd_leaves) && length(sd_leaves) > 0){
    descriptor_terms <- do.call(
      c,
      unname(lapply(sd_leaves, .random_sd_leaf_term_map, scaled_vars = scaled_vars))
    )
    term_map <- descriptor_terms[random_sd_cols]
    names(term_map) <- random_sd_cols
    term_map <- term_map[!is.na(term_map)]
    if(length(term_map) > 0){
      return(term_map)
    }
  }

  stats::setNames(character(), character())
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

  if(isTRUE(homogeneous)){
    if(length(column_terms) == 1L){
      column_terms <- .random_sd_term_from_column_name(leaves$column_names[[1L]])
    }else{
      stop(
        "Cannot unscale homogeneous random-effect SD '",
        paste(available_leaf_names, collapse = "', '"),
        "' because its block contains scaled random-slope columns. ",
        "Use a heterogeneous random-effect SD structure or leave samples on the fitted scale.",
        call. = FALSE
      )
    }
  }

  missing_leaf_names <- setdiff(leaf_names, random_sd_cols)
  group_key <- .random_sd_group_key_from_leaf_key(
    leaf_key = leaf_key,
    leaf_names = leaf_names,
    column_terms = column_terms,
    prefix = prefix
  )
  if(length(missing_leaf_names) > 0L){
    stop(
      "Random-effect SD unscaling for block '", group_key,
      "' requires all SD columns from the block. Missing: ",
      paste0("'", missing_leaf_names, "'", collapse = ", "),
      ".",
      call. = FALSE
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

.random_sd_leaf_term_map <- function(leaves, scaled_vars){

  if(is.null(leaves$leaf_terms)){
    return(stats::setNames(character(), character()))
  }

  if(.random_sd_structured_leaves(leaves)){
    return(stats::setNames(character(), character()))
  }

  out <- leaves$leaf_terms
  homogeneous <- identical(unique(unname(leaves$leaf_terms)), "sd")
  if(!isTRUE(homogeneous)){
    return(out)
  }

  if(is.null(leaves$leaf_names_by_column) || is.null(leaves$column_names)){
    return(out)
  }

  leaf_names <- unique(leaves$leaf_names_by_column)
  if(length(leaf_names) != 1L){
    return(out)
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

  if(length(column_terms) == 1L){
    if(isTRUE(scaled_columns)){
      out[[leaf_names]] <- column_terms
    }
    return(out)
  }

  if(any(scaled_columns)){
    stop(
      "Cannot unscale homogeneous random-effect SD '", leaf_names,
      "' because its block contains scaled random-slope columns. ",
      "Use a heterogeneous random-effect SD structure or leave samples on the fitted scale.",
      call. = FALSE
    )
  }

  out
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

  term <- sub("\\[[^]]+\\]$", "", term)
  components <- .bt_random_effect_term_components(term)
  any(components %in% scaled_vars)
}

.random_sd_group_key <- function(col, term_map, prefix){

  base_col <- sub("\\[[^]]+\\]$", "", col)
  rest <- .formula_scale_strip_prefix(base_col, prefix, "__xREx__")
  term <- unname(term_map[[col]])
  term_core <- gsub("\\[[^]]+\\]", "", term)
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

  out <- array(NA_real_, dim = c(nrow(posterior), n_terms, n_terms))
  for(draw_i in seq_len(nrow(posterior))){
    out[draw_i, , ] <- L[draw_i, , ] %*% t(L[draw_i, , ])
  }

  out
}
