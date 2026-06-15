## Compact print methods for random-effect specification objects.

.bt_print_spec_lines <- function(lines, silent){

  check_bool(silent, "silent")
  if(!silent){
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }

  invisible(lines)
}

.bt_format_random_print_number <- function(x, digits_estimates){

  if(is.numeric(x)){
    return(paste(round(x, digits_estimates), collapse = ", "))
  }

  paste(as.character(x), collapse = ", ")
}

.bt_format_random_print_prior <- function(x, digits_estimates){

  if(is.null(x)){
    return("none")
  }
  if(is.prior(x)){
    return(print(
      x,
      digits_estimates = digits_estimates,
      silent = TRUE,
      inline = TRUE
    ))
  }
  if(inherits(x, "prior_lkj")){
    return(.bt_format_prior_lkj_inline(x, digits_estimates))
  }
  if(inherits(x, "random_sd_source")){
    return(.bt_format_random_sd_source_inline(x))
  }

  paste0("<", class(x)[[1L]], ">")
}

.bt_format_random_print_terms <- function(terms){

  if(is.null(terms) || length(terms) == 0L){
    return("none")
  }

  term_names <- names(terms)
  if(is.null(term_names) || all(!nzchar(term_names))){
    return(paste(as.character(terms), collapse = ", "))
  }
  if(length(term_names) == length(terms) && all(nzchar(term_names))){
    if(is.character(terms)){
      return(paste(paste0(term_names, " = ", terms), collapse = ", "))
    }
    return(paste(term_names, collapse = ", "))
  }

  paste(as.character(terms), collapse = ", ")
}

.bt_format_prior_lkj_inline <- function(x, digits_estimates){

  paste0(
    "prior_lkj(eta = ", .bt_format_random_print_number(x$eta, digits_estimates),
    ", include_correlation = ", as.character(x$include_correlation),
    ", include_primitives = ", as.character(x$include_primitives),
    ")"
  )
}

.bt_format_random_monitor_inline <- function(x){

  paste0(
    "random_monitor(latent = ", as.character(x$latent),
    ", coefficients = ", as.character(x$coefficients),
    ", correlation = ", as.character(x$correlation),
    ", lkj_primitives = ", as.character(x$lkj_primitives),
    ")"
  )
}

.bt_format_random_new_levels_inline <- function(x){

  paste0(
    "random_new_levels(allow = ", as.character(x$allow),
    ", method = \"", x$method, "\")"
  )
}

.bt_format_random_covariance_inline <- function(x, digits_estimates){

  if(is.null(x)){
    return("inherit")
  }

  parts <- paste0(
    "structure = ",
    if(is.null(x$structure)) "formula-owned" else x$structure
  )
  if(!is.null(x$sd)){
    parts <- c(parts, paste0(
      "sd = ", .bt_format_random_print_prior(x$sd, digits_estimates)
    ))
  }
  if(!is.null(x$cor)){
    parts <- c(parts, paste0(
      "cor = ", .bt_format_random_print_prior(x$cor, digits_estimates)
    ))
  }
  if(!is.null(x$rho)){
    parts <- c(parts, paste0(
      "rho = ", .bt_format_random_print_prior(x$rho, digits_estimates)
    ))
  }
  explicit_fields <- attr(x, "explicit_fields", exact = TRUE)
  if(!is.null(x$rho) || "rho_scale" %in% explicit_fields){
    parts <- c(parts, paste0("rho_scale = ", x$rho_scale))
  }

  paste0("random_covariance(", paste(parts, collapse = ", "), ")")
}

.bt_format_parameter_source_inline <- function(x){

  paste0(
    "parameter_source(name = \"", x$name,
    "\", shape = \"", x$shape,
    "\", values = ", if(is.null(x$values)) "none" else "function",
    ")"
  )
}

.bt_format_parameter_source_label <- function(x){

  if(is.null(x)){
    return("none")
  }
  if(identical(x$shape, "row")){
    return(paste0(x$name, "[row]"))
  }

  x$name
}

.bt_format_random_sd_source_inline <- function(x){

  paste0("random_sd_source(source = ", .bt_format_parameter_source_label(x$source), ")")
}

.bt_format_random_allocation_ref_inline <- function(x){

  paste0(
    "allocation_ref(allocation = \"", x$allocation,
    "\", component = \"", x$component,
    "\")"
  )
}

.bt_format_random_allocation_label <- function(allocations, i){

  allocation_names <- names(allocations)
  if(!is.null(allocation_names) && length(allocation_names) >= i &&
     !is.na(allocation_names[[i]]) && nzchar(allocation_names[[i]])){
    return(allocation_names[[i]])
  }

  allocation <- allocations[[i]]
  if(inherits(allocation, "random_variance_allocation") &&
     !is.null(allocation$name) && length(allocation$name) == 1L &&
     !is.na(allocation$name) && nzchar(allocation$name)){
    return(allocation$name)
  }

  paste0("#", i)
}

.bt_format_random_allocation_summary <- function(allocation){

  allocations <- .bt_random_allocation_list(allocation)
  if(length(allocations) == 0L){
    return("none")
  }

  labels <- vapply(
    seq_along(allocations),
    function(i) .bt_format_random_allocation_label(allocations, i),
    character(1)
  )

  paste0(
    length(allocations),
    if(length(allocations) == 1L) " allocation" else " allocations",
    " (", paste(labels, collapse = ", "), ")"
  )
}

.bt_format_random_block_inline <- function(x, digits_estimates){

  parts <- character()
  if(!is.null(x$sd)){
    parts <- c(parts, paste0("sd = ", .bt_format_random_print_prior(x$sd, digits_estimates)))
  }
  if(!is.null(x$sd_source)){
    parts <- c(parts, paste0("sd_source = ", .bt_format_random_sd_source_inline(x$sd_source)))
  }
  if(!is.null(x$covariance)){
    parts <- c(parts, paste0("covariance = ", .bt_format_random_covariance_inline(x$covariance, digits_estimates)))
  }
  if(!is.null(x$monitor)){
    parts <- c(parts, paste0("monitor = ", .bt_format_random_monitor_inline(x$monitor)))
  }
  if(!is.null(x$new_levels)){
    parts <- c(parts, paste0("new_levels = ", .bt_format_random_new_levels_inline(x$new_levels)))
  }
  if(!is.null(x$terms)){
    parts <- c(parts, paste0("terms = ", .bt_format_random_print_terms(x$terms)))
  }

  if(length(parts) == 0L){
    return("random_block(inherits top-level settings)")
  }

  paste0("random_block(", paste(parts, collapse = ", "), ")")
}

.bt_format_random_blocks_summary <- function(blocks, digits_estimates){

  if(length(blocks) == 0L){
    return("none")
  }

  block_names <- names(blocks)
  parts <- vapply(seq_along(blocks), function(i){
    paste0(
      block_names[[i]], " = ",
      .bt_format_random_block_inline(blocks[[i]], digits_estimates)
    )
  }, character(1))

  paste(parts, collapse = "; ")
}

#' @export
print.prior_lkj <- function(x, digits_estimates = 2, silent = FALSE, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  if(!inherits(x, "prior_lkj")){
    stop("'x' must be created with prior_lkj().", call. = FALSE)
  }

  .bt_print_spec_lines(c(
    "prior_lkj()",
    paste0("  eta: ", .bt_format_random_print_number(x$eta, digits_estimates)),
    paste0("  include_correlation: ", as.character(x$include_correlation)),
    paste0("  include_primitives: ", as.character(x$include_primitives))
  ), silent = silent)
}

#' @export
print.random_monitor <- function(x, silent = FALSE, ...){

  check_bool(silent, "silent")
  if(!inherits(x, "random_monitor")){
    stop("'x' must be created with random_monitor().", call. = FALSE)
  }

  .bt_print_spec_lines(c(
    "random_monitor()",
    paste0("  latent: ", as.character(x$latent)),
    paste0("  coefficients: ", as.character(x$coefficients)),
    paste0("  correlation: ", as.character(x$correlation)),
    paste0("  lkj_primitives: ", as.character(x$lkj_primitives))
  ), silent = silent)
}

#' @export
print.random_new_levels <- function(x, silent = FALSE, ...){

  check_bool(silent, "silent")
  if(!inherits(x, "random_new_levels")){
    stop("'x' must be created with random_new_levels().", call. = FALSE)
  }

  .bt_print_spec_lines(c(
    "random_new_levels()",
    paste0("  allow: ", as.character(x$allow)),
    paste0("  method: ", x$method)
  ), silent = silent)
}

#' @export
print.random_covariance <- function(x, digits_estimates = 2, silent = FALSE,
                                    ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  .bt_check_random_covariance(x)

  .bt_print_spec_lines(c(
    "random_covariance()",
    paste0("  structure: ", if(is.null(x$structure)) "formula-owned" else x$structure),
    paste0("  sd: ", .bt_format_random_print_prior(x$sd, digits_estimates)),
    paste0("  cor: ", .bt_format_random_print_prior(x$cor, digits_estimates)),
    paste0("  rho: ", .bt_format_random_print_prior(x$rho, digits_estimates)),
    paste0("  rho_scale: ", x$rho_scale)
  ), silent = silent)
}

#' @export
print.random_block <- function(x, digits_estimates = 2, silent = FALSE, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  if(!inherits(x, "random_block")){
    stop("'x' must be created with random_block().", call. = FALSE)
  }

  .bt_print_spec_lines(c(
    "random_block()",
    paste0(
      "  sd: ",
      if(is.null(x$sd)) "inherit" else .bt_format_random_print_prior(x$sd, digits_estimates)
    ),
    paste0(
      "  sd_source: ",
      if(is.null(x$sd_source)) "none" else .bt_format_random_sd_source_inline(x$sd_source)
    ),
    paste0(
      "  covariance: ",
      if(is.null(x$covariance)) "inherit" else .bt_format_random_covariance_inline(x$covariance, digits_estimates)
    ),
    paste0(
      "  monitor: ",
      if(is.null(x$monitor)) "inherit" else .bt_format_random_monitor_inline(x$monitor)
    ),
    paste0(
      "  new_levels: ",
      if(is.null(x$new_levels)) "inherit" else .bt_format_random_new_levels_inline(x$new_levels)
    ),
    paste0("  terms: ", .bt_format_random_print_terms(x$terms)),
    paste0("  allocation: ", if(is.null(x$allocation)) "none" else "<reserved>")
  ), silent = silent)
}

#' @export
print.random_variance_allocation <- function(x, digits_estimates = 2,
                                             silent = FALSE, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  if(!inherits(x, "random_variance_allocation")){
    stop("'x' must be created with random_variance_allocation().", call. = FALSE)
  }

  .bt_print_spec_lines(c(
    "random_variance_allocation()",
    paste0("  name: ", if(is.null(x$name)) "none" else x$name),
    paste0("  terms: ", .bt_format_random_print_terms(x$terms)),
    paste0("  sd: ", .bt_format_random_print_prior(x$sd, digits_estimates)),
    paste0(
      "  sd_source: ",
      if(is.null(x$sd_source)) "none" else .bt_format_random_sd_source_inline(x$sd_source)
    ),
    paste0("  weights: ", .bt_format_random_print_prior(x$weights, digits_estimates)),
    paste0(
      "  parent: ",
      if(is.null(x$parent)) "none" else .bt_format_random_allocation_ref_inline(x$parent)
    ),
    paste0("  target: ", x$target),
    paste0("  scale: ", x$scale)
  ), silent = silent)
}

#' @export
print.random_allocation_ref <- function(x, silent = FALSE, ...){

  check_bool(silent, "silent")
  .bt_check_random_allocation_ref(x)

  .bt_print_spec_lines(c(
    "allocation_ref()",
    paste0("  allocation: ", x$allocation),
    paste0("  component: ", x$component)
  ), silent = silent)
}

#' @export
print.prior_random <- function(x, digits_estimates = 2, silent = FALSE, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  .bt_check_prior_random(x)

  .bt_print_spec_lines(c(
    "prior_random()",
    paste0("  sd: ", .bt_format_random_print_prior(x$sd, digits_estimates)),
    paste0(
      "  covariance: ",
      .bt_format_random_covariance_inline(x$covariance, digits_estimates)
    ),
    paste0("  monitor: ", .bt_format_random_monitor_inline(x$monitor)),
    paste0("  new_levels: ", .bt_format_random_new_levels_inline(x$new_levels)),
    paste0("  allocation: ", .bt_format_random_allocation_summary(x$allocation)),
    paste0("  blocks: ", .bt_format_random_blocks_summary(x$blocks, digits_estimates))
  ), silent = silent)
}

#' @export
print.parameter_source <- function(x, silent = FALSE, ...){

  check_bool(silent, "silent")
  .bt_check_parameter_source(x)

  .bt_print_spec_lines(c(
    "parameter_source()",
    paste0("  name: ", x$name),
    paste0("  shape: ", x$shape),
    paste0("  values: ", if(is.null(x$values)) "none" else "function")
  ), silent = silent)
}

#' @export
print.random_sd_source <- function(x, silent = FALSE, ...){

  check_bool(silent, "silent")
  .bt_check_random_sd_source(x)

  .bt_print_spec_lines(c(
    "random_sd_source()",
    paste0("  source: ", .bt_format_parameter_source_label(x$source))
  ), silent = silent)
}
