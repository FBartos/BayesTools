# Hierarchical centering of mean-centered random-effect blocks.
#
# A block with parameterization = "mean_centered" samples its group effects
# around a location instead of zero,
#   eta_g = b_g + mu_g,   mu_g = A c_a + A_g c_ga + P eta_parent(g),
# where c are the effective coefficients of fixed terms that the block absorbs
# (the values their design columns are multiplied by in the linear predictor),
# A a constant map with Z A = X_a on the fitted rows up to rounding (the term's
# columns lie in the span of the block design Z), A_g = e_k0 x_g' a group map
# for terms whose columns are constant within the block's groups (they enter
# the location of the block's all-ones column k0), and P eta_parent(g) the
# location of a linked parent block whose groups contain the block's groups
# and whose design is exactly Z P. The fitted predictor drops the absorbed
# terms and the linked parent's term and reads Z_i eta_g(i), plus the rounding
# residual (X_a - Z A) c_a of an inexact constant map, so it equals the
# semantic predictor X c + Z b whatever the map. Blocks are translated parents
# first, so (b, c, ...) -> (eta, c, ...) is triangular with a unit diagonal:
# the posterior and every prior are unchanged (conditional on discrete prior
# states). The semantic coefficients b = eta - mu and the standardized latents
# z = (D L)^-1 b keep the evaluator, prediction, and bridge conventions.

# The fixed terms of a formula parameter that a mean-centered block can
# absorb: the JAGS design columns of each term and the effective coefficient
# expression of each column.
.bt_random_effect_mean_fixed_design <- function(parameter, model_matrix, terms_indexes,
                                                model_terms, model_terms_type,
                                                prior_list, has_intercept,
                                                log_intercept){

  terms <- list()
  if(has_intercept){
    terms[[1L]] <- .bt_dnode_linear_predictor_term(
      parameter = parameter, model_term = "intercept", type = "intercept",
      columns = 1L, prior = prior_list[["intercept"]], log = log_intercept
    )
  }
  for(i in unique(terms_indexes[terms_indexes > 0])){
    terms[[length(terms) + 1L]] <- .bt_dnode_linear_predictor_term(
      parameter = parameter, model_term = model_terms[i],
      type = if(identical(unname(model_terms_type[i]), "factor")) "factor" else "continuous",
      columns = which(terms_indexes == i), prior = prior_list[[model_terms[i]]]
    )
  }
  terms <- lapply(terms, function(term){
    term$expressions <- .bt_dnode_linear_predictor_term_column_expressions(term)
    term
  })

  list(matrix = unname(as.matrix(model_matrix)), terms = terms)
}

# The row group of every group of 'child' within 'parent' (both row group
# maps), or NULL when a child group has no rows or spans several parent groups.
.bt_random_effect_group_nesting <- function(child_map, child_groups, parent_map){

  parent_of <- lapply(split(parent_map, factor(child_map, levels = seq_len(child_groups))), unique)
  if(any(lengths(parent_of) != 1L)){
    return(NULL)
  }
  as.integer(unlist(parent_of, use.names = FALSE))
}

# Translation order: a topological order of proper nesting (a block whose
# groups nest in another mean-centered block's groups, and do not equal them,
# comes after it) that takes the earliest block in formula order among those
# ready.
.bt_random_effect_mean_order <- function(blocks){

  n_blocks <- length(blocks)
  nested_in <- matrix(FALSE, n_blocks, n_blocks)
  for(child in seq_len(n_blocks)){
    for(parent in seq_len(n_blocks)[-child]){
      nested_in[child, parent] <- !is.null(.bt_random_effect_group_nesting(
        blocks[[child]]$group_map, blocks[[child]]$n_groups, blocks[[parent]]$group_map
      ))
    }
  }
  strict <- nested_in & !t(nested_in)
  order <- integer()
  remaining <- seq_len(n_blocks)
  while(length(remaining) > 0L){
    ready <- remaining[vapply(remaining, function(block){
      !any(strict[block, setdiff(remaining, block)])
    }, logical(1))]
    order <- c(order, ready[[1L]])
    remaining <- setdiff(remaining, ready[[1L]])
  }
  order
}

# The constant map A (ncol(Z) x ncol(X)) of least squares, with entries
# within 1e-12 of an integer snapped to it, its residual X - Z A, and whether
# each column of X lies in the span of Z (a residual of at most 1e-9 times the
# column's largest absolute value). The map decides only which terms a block
# absorbs; the fitted predictor keeps the residual.
.bt_random_effect_mean_constant_map <- function(Z, X){

  A <- qr.coef(qr(Z), X)
  A[is.na(A)] <- 0
  A <- matrix(A, nrow = ncol(Z), ncol = ncol(X))
  near <- abs(A - round(A)) <= 1e-12 * pmax(1, abs(A))
  A[near] <- round(A[near])
  residual <- X - Z %*% A
  scale <- apply(abs(X), 2L, max)

  list(map = A, residual = residual, in_span = apply(abs(residual), 2L, max) <= 1e-9 * scale)
}

# The group values (groups x columns) of X when every column is constant
# within each group, or NULL (also when a declared group has no rows).
.bt_random_effect_mean_group_values <- function(X, group_map, n_groups){

  rows <- split(seq_len(nrow(X)), factor(group_map, levels = seq_len(n_groups)))
  if(any(lengths(rows) == 0L)){
    return(NULL)
  }
  constant <- all(vapply(rows, function(r){
    all(X[r, , drop = FALSE] == matrix(X[r[[1L]], ], nrow = length(r), ncol = ncol(X), byrow = TRUE))
  }, logical(1)))
  if(!constant){
    return(NULL)
  }
  X[vapply(rows, `[[`, integer(1), 1L), , drop = FALSE]
}

# A weighted sum of JAGS expressions with literal weights (zero weights are
# skipped; "0" when every weight is zero).
.bt_JAGS_weighted_sum <- function(weights, expressions){

  parts <- character()
  for(j in seq_along(weights)){
    w <- weights[[j]]
    if(w == 0){
      next
    }
    term <- if(w == 1){
      expressions[[j]]
    }else if(w == -1){
      paste0("-(", expressions[[j]], ")")
    }else{
      paste0(.bt_JAGS_numeric_literal(w), " * (", expressions[[j]], ")")
    }
    parts <- c(parts, term)
  }
  if(length(parts) == 0L){
    return("0")
  }
  out <- parts[[1L]]
  for(part in parts[-1L]){
    out <- if(startsWith(part, "-")) paste0(out, " - ", substring(part, 2L)) else paste0(out, " + ", part)
  }
  out
}

# Resolve the translations of the mean-centered blocks of one formula
# parameter. 'blocks' are list(term = <compiled random term>, centering =
# <.JAGS_random_effect_formula() centering metadata>); 'fixed' is
# .bt_random_effect_mean_fixed_design(). Returns, per block in input order,
# the translation specification, its JAGS syntax, and its JAGS data.
.bt_random_effect_mean_translations <- function(blocks, fixed){

  n_blocks <- length(blocks)
  terms <- lapply(blocks, `[[`, "term")
  block_names <- vapply(terms, `[[`, character(1), "block_name")
  specs <- vector("list", n_blocks)
  owned <- character()
  linked <- logical(n_blocks)
  done <- integer()

  for(b in .bt_random_effect_mean_order(terms)){
    term <- terms[[b]]
    Z <- unname(as.matrix(term$model_matrix))
    stem <- term$parameter_stem
    if(qr(Z)$rank < ncol(Z)){
      stop("Mean-centered parameterization is unavailable for random-effect block '",
           term$block_name, "': its random-effect design is rank deficient.", call. = FALSE)
    }
    ones <- which(colSums(Z != 1) == 0L)

    # whole terms only: a term partly in the span stays in the predictor
    constant <- list(terms = character(), columns = integer(), expressions = character(),
                     map = matrix(0, ncol(Z), 0L), residual = matrix(0, nrow(Z), 0L))
    group <- list(terms = character(), columns = integer(), expressions = character(),
                  values = matrix(0, term$n_groups, 0L))
    for(fixed_term in fixed$terms){
      if(fixed_term$model_term %in% owned){
        next
      }
      X <- fixed$matrix[, fixed_term$columns, drop = FALSE]
      solved <- .bt_random_effect_mean_constant_map(Z, X)
      if(all(solved$in_span)){
        constant$terms <- c(constant$terms, fixed_term$model_term)
        constant$columns <- c(constant$columns, fixed_term$columns)
        constant$expressions <- c(constant$expressions, fixed_term$expressions)
        constant$map <- cbind(constant$map, solved$map)
        constant$residual <- cbind(constant$residual, solved$residual)
        owned <- c(owned, fixed_term$model_term)
        next
      }
      values <- if(length(ones) > 0L){
        .bt_random_effect_mean_group_values(X, term$group_map, term$n_groups)
      }
      if(!is.null(values)){
        group$terms <- c(group$terms, fixed_term$model_term)
        group$columns <- c(group$columns, fixed_term$columns)
        group$expressions <- c(group$expressions, fixed_term$expressions)
        group$values <- cbind(group$values, values)
        owned <- c(owned, fixed_term$model_term)
      }
    }

    # the finest earlier block whose groups contain this block's groups (every
    # group observed) and whose design is exactly this block's design times a
    # constant map
    parent <- NULL
    candidates <- done[!linked[done]]
    candidates <- candidates[order(-vapply(terms[candidates], `[[`, numeric(1), "n_groups"))]
    for(p in candidates){
      index <- .bt_random_effect_group_nesting(term$group_map, term$n_groups, terms[[p]]$group_map)
      if(is.null(index)){
        next
      }
      solved <- .bt_random_effect_mean_constant_map(Z, unname(as.matrix(terms[[p]]$model_matrix)))
      if(any(solved$residual != 0)){
        next
      }
      parent <- list(block = terms[[p]]$block_name, map = solved$map, index = index,
                     location_name = specs[[p]]$location_name,
                     data_name = paste0(stem, "_xRE_PARENTx"))
      linked[p] <- TRUE
      break
    }
    if(length(constant$terms) == 0L && length(group$terms) == 0L && is.null(parent)){
      stop("Mean-centered parameterization is unavailable for random-effect block '",
           term$block_name, "': its random-effect design spans no fixed term that another ",
           "mean-centered block does not already absorb, and no mean-centered parent block. ",
           "Use parameterization = \"centered\" or \"noncentered\" for this block.",
           call. = FALSE)
    }

    residual_term <- NULL
    if(length(constant$terms) > 0L){
      if(any(constant$residual != 0)){
        # the rounding residual of the map stays in the predictor (also when
        # a child absorbs the block's locations)
        m <- length(constant$expressions)
        constant$data_name <- paste0(stem, "_xRE_RESIDx")
        constant$coefficient_name <- paste0(stem, "_xRE_ACOEFx")
        residual_term <- paste0("inprod(", constant$coefficient_name, "[1:", m, "], ",
                                constant$data_name, "[i,1:", m, "])")
      }else{
        constant$residual <- NULL
      }
    }else{
      constant <- NULL
    }
    if(length(group$terms) > 0L){
      group$column <- ones[[1L]]
      group$data_name <- paste0(stem, "_xRE_GLEVx")
      group$coefficient_name <- paste0(stem, "_xRE_GCOEFx")
      dimnames(group$values) <- NULL
    }else{
      group <- NULL
    }

    specs[[b]] <- list(
      schema_version = 2L,
      block_name = term$block_name,
      location_name = paste0(stem, "_xRE_MEANx"),
      coefficient_name = paste0(stem, "_xRE_COEFx"),
      latent_name = paste0(stem, "_xRE_Zx"),
      group_map_name = paste0(stem, "_xRE_MAPx"),
      constant = constant,
      group = group,
      parent = parent,
      absorbed_terms = c(constant$terms, group$terms),
      predictor_term = blocks[[b]]$centering$predictor_term,
      residual_term = residual_term
    )
    done <- c(done, b)
  }

  # a linked parent's contribution is part of its child's location
  for(b in seq_len(n_blocks)){
    if(!is.null(specs[[b]]$parent)){
      p <- match(specs[[b]]$parent$block, block_names)
      specs[[p]]$child <- block_names[[b]]
    }
  }

  lapply(seq_len(n_blocks), function(b){
    syntax <- .bt_JAGS_mean_centered_random(specs[[b]], terms[[b]], blocks[[b]]$centering)
    list(spec = specs[[b]], syntax = syntax$syntax, data = syntax$data)
  })
}

# The block's term of the fitted predictor: its group locations read through
# the block design (the column of an exact indicator layout, the single
# column of an all-ones design, or the inner product with the design row).
.bt_random_effect_mean_predictor_term <- function(stem, model_matrix, indicator_layout){

  location <- paste0(stem, "_xRE_MEANx")
  map <- paste0(stem, "_xRE_MAPx[i]")
  K <- ncol(model_matrix)
  if(isTRUE(indicator_layout)){
    return(paste0(location, "[", map, ",", stem, "_xRE_COLx[i]]"))
  }
  if(K == 1L && all(model_matrix == 1)){
    return(paste0(location, "[", map, ",1]"))
  }
  paste0("inprod(", location, "[", map, ",1:", K, "], ", stem, "_xRE_DATAx[i,1:", K, "])")
}

# JAGS syntax and data of a translated block: the locations, the group
# effects as scalar nodes, the semantic coefficients, and the standardized
# latents. Correlated blocks sample sequential conditionals in the arithmetic
# of their noncentered transforms (the lower LKJ factor, the compound-symmetry
# prefix recurrence, or the AR(1) recurrence): column i given the earlier
# latents has the standardized mean `offset(i)` and SD `scale(i)`. At the
# representability limit (an LKJ primitive or a correlation that rounds to a
# boundary) a conditional precision is not finite and JAGS stops; the
# arithmetic is not altered there.
.bt_JAGS_mean_centered_random <- function(spec, term, centering){

  K <- ncol(term$model_matrix)
  G <- term$n_groups
  stem <- term$parameter_stem
  sd_name <- centering$sd_name
  syntax <- character()
  data <- list()

  per_group <- !is.null(spec$group) || !is.null(spec$parent)
  constant_part <- vapply(seq_len(K), function(k){
    if(is.null(spec$constant)) "0" else .bt_JAGS_weighted_sum(spec$constant$map[k, ], spec$constant$expressions)
  }, character(1))
  if(!is.null(spec$constant$residual)){
    syntax <- c(syntax, paste0(spec$constant$coefficient_name, "[", seq_along(spec$constant$expressions),
                               "] <- ", spec$constant$expressions))
    data[[spec$constant$data_name]] <- spec$constant$residual
  }
  if(!is.null(spec$group)){
    syntax <- c(syntax, paste0(spec$group$coefficient_name, "[", seq_along(spec$group$expressions),
                               "] <- ", spec$group$expressions))
    data[[spec$group$data_name]] <- spec$group$values
  }
  if(!is.null(spec$parent)){
    data[[spec$parent$data_name]] <- spec$parent$index
  }

  # a scalar location that is one coefficient expression is written inline
  inline <- !per_group && K == 1L && !is.null(spec$constant) &&
    sum(spec$constant$map != 0) == 1L && all(spec$constant$map[spec$constant$map != 0] == 1)
  location_name <- paste0(stem, "_xRE_LOCx")
  location <- if(inline){
    function(k) constant_part[[k]]
  }else if(!per_group){
    syntax <- c(syntax, paste0(location_name, "[", seq_len(K), "] <- ", constant_part))
    function(k) paste0(location_name, "[", k, "]")
  }else{
    lines <- vapply(seq_len(K), function(k){
      parts <- if(constant_part[[k]] != "0") constant_part[[k]]
      if(!is.null(spec$group) && k == spec$group$column){
        m <- length(spec$group$expressions)
        parts <- c(parts, paste0("inprod(", spec$group$data_name, "[g,1:", m, "], ",
                                 spec$group$coefficient_name, "[1:", m, "])"))
      }
      if(!is.null(spec$parent)){
        parent_part <- .bt_JAGS_weighted_sum(spec$parent$map[k, ], paste0(
          spec$parent$location_name, "[", spec$parent$data_name, "[g],", seq_len(ncol(spec$parent$map)), "]"
        ))
        if(parent_part != "0") parts <- c(parts, parent_part)
      }
      paste0("   ", location_name, "[g,", k, "] <- ", if(length(parts) == 0L) "0" else paste(parts, collapse = " + "))
    }, character(1))
    syntax <- c(syntax, paste0(" for(g in 1:", G, "){\n", paste(lines, collapse = "\n"), "\n }"))
    function(k) paste0(location_name, "[g,", k, "]")
  }

  mean <- spec$location_name
  coef <- spec$coefficient_name
  z <- spec$latent_name
  if(inline){
    # the scalar random intercept (one coefficient expression)
    syntax <- c(syntax, paste0(
      " for(g in 1:", G, "){\n",
      "   ", mean, "[g,1] ~ dnorm(", location(1L), ", pow(", sd_name, "[1], -2))\n",
      "   ", coef, "[g,1] <- ", mean, "[g,1] - ", location(1L), "\n",
      "   ", z, "[g,1] <- ", coef, "[g,1] / ", sd_name, "[1]\n",
      " }"))
    return(list(syntax = paste0(paste(syntax, collapse = "\n"), "\n"), data = data))
  }

  # the standardized conditional mean and SD of column i given the earlier
  # latents of the group, and the recurrence nodes they read
  structured <- centering$structured
  offset <- function(i) NULL
  scale <- function(i) NULL
  recurrence <- function(i) character()
  if(!is.null(centering$cholesky_name) && K > 1L){
    L <- centering$cholesky_name
    offset <- function(i) if(i > 1L) paste0(L, "[", i, ",", seq_len(i - 1L), "] * ", z, "[g,", seq_len(i - 1L), "]", collapse = " + ")
    scale <- function(i) if(i > 1L) paste0(L, "[", i, ",", i, "]")
  }else if(!is.null(structured) && K > 1L && structured$structure %in% c("cs", "hcs")){
    rho <- structured$rho_name
    diagonal <- paste0(stem, "_xRE_CS_DIAGx")
    update <- paste0(stem, "_xRE_CS_UPDATEx")
    prefix <- paste0(stem, "_xRE_CS_PREFIXx")
    syntax <- c(syntax, paste0(
      " for(i in 2:", K, "){\n",
      "   ", diagonal, "[i] <- sqrt((1 - ", rho, ") * (1 + (i - 1) * ", rho, ") / (1 + (i - 2) * ", rho, "))\n",
      "   ", update, "[i] <- ", rho, " * sqrt((1 - ", rho, ") / ((1 + (i - 2) * ", rho, ") * (1 + (i - 1) * ", rho, ")))\n",
      " }"))
    offset <- function(i) if(i > 1L) paste0(prefix, "[g,", i - 1L, "]")
    scale <- function(i) if(i > 1L) paste0(diagonal, "[", i, "]")
    recurrence <- function(i){
      if(i == K) return(character())
      paste0("   ", prefix, "[g,", i, "] <- ",
             if(i == 1L) paste0(rho, " * ", z, "[g,1]") else paste0(prefix, "[g,", i - 1L, "] + ", update, "[", i, "] * ", z, "[g,", i, "]"))
    }
  }else if(!is.null(structured) && K > 1L){
    rho <- structured$rho_name
    innovation <- paste0(stem, "_xRE_AR_INNOVx")
    unit <- paste0(stem, "_xRE_UNIT_COEFx")
    syntax <- c(syntax, paste0(innovation, " <- sqrt(1 - pow(", rho, ", 2))"))
    offset <- function(i) if(i > 1L) paste0(rho, " * ", unit, "[g,", i - 1L, "]")
    scale <- function(i) if(i > 1L) innovation
    recurrence <- function(i) if(i < K) paste0("   ", unit, "[g,", i, "] <- ", coef, "[g,", i, "] / ", sd_name, "[", i, "]")
  }

  lines <- character()
  for(i in seq_len(K)){
    sd_i <- paste0(sd_name, "[", i, "]")
    conditional_mean <- if(is.null(offset(i))) location(i) else paste0(location(i), " + ", sd_i, " * (", offset(i), ")")
    conditional_sd <- if(is.null(scale(i))) sd_i else paste0(sd_i, " * ", scale(i))
    standardized <- paste0(coef, "[g,", i, "] / ", sd_i)
    if(!is.null(offset(i))) standardized <- paste0("(", standardized, " - (", offset(i), "))")
    if(!is.null(scale(i))) standardized <- paste0(standardized, " / ", scale(i))
    lines <- c(lines,
      paste0("   ", mean, "[g,", i, "] ~ dnorm(", conditional_mean, ", pow(", conditional_sd, ", -2))"),
      paste0("   ", coef, "[g,", i, "] <- ", mean, "[g,", i, "] - ", location(i)),
      paste0("   ", z, "[g,", i, "] <- ", standardized),
      recurrence(i))
  }
  syntax <- c(syntax, paste0(" for(g in 1:", G, "){\n", paste(lines, collapse = "\n"), "\n }"))

  list(syntax = paste0(paste(syntax, collapse = "\n"), "\n"), data = data)
}
