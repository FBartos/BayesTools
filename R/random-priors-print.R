## Prior-style print methods for random-effect specification objects.

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

  paste0("<", class(x)[[1L]], ">")
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

.bt_format_prior_lkj_distribution <- function(x, digits_estimates){

  paste0("LKJ(eta = ", .bt_format_random_print_number(x$eta, digits_estimates), ")")
}

.bt_format_random_covariance_structure <- function(x){

  if(is.null(x) || is.null(x$structure)){
    return("formula-owned")
  }

  x$structure
}

.bt_format_random_sigma_name <- function(component = NULL, index = NULL){

  if(is.null(component)){
    out <- "sigma"
  }else if(grepl("^[A-Za-z][A-Za-z0-9_]*$", component)){
    out <- paste0("sigma_", component)
  }else{
    out <- paste0("sigma[", component, "]")
  }

  if(!is.null(index)){
    out <- paste0(out, "[", index, "]")
  }

  out
}

.bt_format_random_prior_equation <- function(lhs, rhs, operator = "~"){

  paste0(lhs, " ", operator, " ", rhs)
}

.bt_format_random_covariance_prior_lines <- function(x, digits_estimates,
                                                     include_structure = FALSE,
                                                     include_sd = TRUE){

  if(is.null(x)){
    return(if(include_structure) "covariance: inherit" else character())
  }
  .bt_check_random_covariance(x)

  header <- character()
  lines <- character()
  if(isTRUE(include_structure)){
    header <- paste0("covariance: ", .bt_format_random_covariance_structure(x))
  }
  if(isTRUE(include_sd) && !is.null(x$sd)){
    lines <- c(lines, .bt_format_random_prior_equation(
      "sigma",
      .bt_format_random_print_prior(x$sd, digits_estimates)
    ))
  }
  if(!is.null(x$cor) && inherits(x$cor, "prior_lkj")){
    lines <- c(lines, .bt_format_random_prior_equation(
      "R",
      .bt_format_prior_lkj_distribution(x$cor, digits_estimates)
    ))
  }
  if(!is.null(x$cor) && !inherits(x$cor, "prior_lkj")){
    lines <- c(lines, .bt_format_random_prior_equation(
      "cor",
      .bt_format_random_print_prior(x$cor, digits_estimates)
    ))
  }
  explicit_fields <- attr(x, "explicit_fields", exact = TRUE)
  if((!is.null(x$cor) && !inherits(x$cor, "prior_lkj")) ||
     "cor_scale" %in% explicit_fields){
    lines <- c(lines, paste0("cor_scale: ", x$cor_scale))
  }

  if(length(header) > 0L){
    if(length(lines) == 0L){
      return(header)
    }
    return(c(header, paste0("  ", lines)))
  }

  lines
}

.bt_format_random_block_header <- function(name = NULL, covariance = NULL){

  header <- if(is.null(name)){
    "block"
  }else{
    paste0("block: ", name)
  }
  if(!is.null(covariance) && !is.null(covariance$structure)){
    header <- paste0(header, " (", covariance$structure, ")")
  }

  header
}

.bt_format_random_monitor_settings <- function(x){

  paste0(
    "monitor: latent = ", as.character(x$latent),
    ", coefficients = ", as.character(x$coefficients),
    ", correlation = ", as.character(x$correlation),
    ", lkj_primitives = ", as.character(x$lkj_primitives)
  )
}

.bt_random_monitor_is_default <- function(x){

  isTRUE(x$latent) &&
    identical(x$coefficients, FALSE) &&
    isTRUE(x$correlation) &&
    identical(x$lkj_primitives, FALSE)
}

.bt_random_new_levels_is_default <- function(x){

  identical(x$method, "error")
}

.bt_format_random_block_term_prior <- function(x){

  if(is.prior(x)){
    return(x)
  }
  if(inherits(x, "random_block") && !is.null(x$sd)){
    return(x$sd)
  }

  NULL
}

.bt_format_random_block_prior_lines <- function(x, digits_estimates,
                                                name = NULL,
                                                include_empty = TRUE){

  header <- .bt_format_random_block_header(name, x$covariance)
  lines <- character()
  if(!is.null(x$sd)){
    lines <- c(lines, .bt_format_random_prior_equation(
      "sigma",
      .bt_format_random_print_prior(x$sd, digits_estimates)
    ))
  }
  if(!is.null(x$sd_source)){
    lines <- c(lines, .bt_format_random_prior_equation(
      "sigma",
      .bt_random_sd_source_label(x$sd_source),
      operator = "="
    ))
  }

  lines <- c(lines, .bt_format_random_covariance_prior_lines(
    x$covariance,
    digits_estimates = digits_estimates,
    include_structure = FALSE,
    include_sd = is.null(x$sd) && is.null(x$sd_source)
  ))

  if(!is.null(x$terms)){
    for(term in names(x$terms)){
      term_prior <- .bt_format_random_block_term_prior(x$terms[[term]])
      if(!is.null(term_prior)){
        lines <- c(lines, .bt_format_random_prior_equation(
          .bt_format_random_sigma_name(term),
          .bt_format_random_print_prior(term_prior, digits_estimates)
        ))
      }
    }
  }
  if(!is.null(x$contrasts)){
    lines <- c(
      lines,
      paste0(
        "contrasts: ",
        paste0(
          names(x$contrasts),
          "=",
          sub("^contr\\.", "", x$contrasts),
          collapse = ", "
        )
      )
    )
  }
  if(!is.null(x$monitor)){
    lines <- c(lines, .bt_format_random_monitor_settings(x$monitor))
  }
  if(!is.null(x$new_levels)){
    lines <- c(lines, paste0("new_levels: ", x$new_levels$method))
  }
  if(!is.null(x$parameterization)){
    lines <- c(lines, paste0("parameterization: ", x$parameterization))
  }

  if(length(lines) == 0L && isTRUE(include_empty)){
    lines <- "inherits defaults"
  }

  if(length(lines) == 0L){
    return(header)
  }

  c(header, paste0("  ", lines))
}

.bt_format_random_allocation_weight_prior <- function(x, digits_estimates){

  if(is.null(x$weights)){
    return(character())
  }

  .bt_format_random_print_prior(x$weights, digits_estimates)
}

.bt_format_random_allocation_inclusion_lines <- function(x, digits_estimates){

  if(!is.list(x$inclusion) || length(x$inclusion) == 0L){
    return(character())
  }

  lines <- character()
  for(component in names(x$inclusion)){
    prob_name <- paste0("p_", component)
    indicator_name <- paste0("I_", component)
    lines <- c(lines, .bt_format_random_prior_equation(
      prob_name,
      .bt_format_random_print_prior(x$inclusion[[component]], digits_estimates)
    ))
    lines <- c(lines, .bt_format_random_prior_equation(
      indicator_name,
      paste0("Bernoulli(", prob_name, ")")
    ))
  }

  lines
}

.bt_format_random_allocation_source_name <- function(x){

  if(is.null(x$parent)){
    role <- if(identical(
      .bt_random_variance_allocation_scale(x),
      "mean_variance"
    )) "common" else "total"
    return(paste0("sigma_", role))
  }

  .bt_format_random_sigma_name(x$parent$component)
}

.bt_format_random_allocation_component_labels <- function(x){

  if(is.null(x$terms)){
    return(character())
  }

  .bt_random_variance_allocation_component_labels(x$terms)
}

.bt_format_random_allocation_block_component_lines <- function(x){

  labels <- .bt_format_random_allocation_component_labels(x)
  source_name <- .bt_format_random_allocation_source_name(x)
  if(length(labels) == 0L){
    return(c(
      .bt_format_random_prior_equation(
        "sigma_block[k]",
        paste0(source_name, " * sqrt(w[k])"),
        operator = "="
      ),
      "terms: resolved from formula"
    ))
  }
  if(is.null(x$weights)){
    return(vapply(labels, function(label){
      gate <- if(is.list(x$inclusion) && label %in% names(x$inclusion)){
        paste0("I_", label)
      }else{
        ""
      }
      .bt_format_random_prior_equation(
        .bt_format_random_sigma_name(label),
        if(nzchar(gate)) paste0(source_name, " * ", gate) else source_name,
        operator = "="
      )
    }, character(1)))
  }

  vapply(seq_along(labels), function(i){
    gate <- if(is.list(x$inclusion) && labels[i] %in% names(x$inclusion)){
      paste0("I_", labels[i], " * ")
    }else{
      ""
    }
    .bt_format_random_prior_equation(
      .bt_format_random_sigma_name(labels[i]),
      paste0(source_name, " * ", gate, "sqrt(w[", i, "])"),
      operator = "="
    )
  }, character(1))
}

.bt_format_random_allocation_sd_component_lines <- function(x){

  source_name <- .bt_format_random_allocation_source_name(x)
  block <- if(is.null(x$terms)){
    "block"
  }else{
    unname(x$terms[[1L]])
  }

  K <- if(!is.null(x$weights)) x$weights$parameters[["K"]] else NA_integer_
  scale <- .bt_random_variance_allocation_scale(x)
  if(is.na(K)){
    multiplier <- if(identical(scale, "mean_variance")) "K * w[k]" else "w[k]"
    return(c(
      .bt_format_random_prior_equation(
        .bt_format_random_sigma_name(block, "k"),
        paste0(source_name, " * sqrt(", multiplier, ")"),
        operator = "="
      ),
      "components: resolved from formula"
    ))
  }

  vapply(seq_len(K), function(i){
    multiplier <- if(identical(scale, "mean_variance")){
      paste0(K, " * w[", i, "]")
    }else{
      paste0("w[", i, "]")
    }
    .bt_format_random_prior_equation(
      .bt_format_random_sigma_name(block, i),
      paste0(source_name, " * sqrt(", multiplier, ")"),
      operator = "="
    )
  }, character(1))
}

.bt_format_random_allocation_lines <- function(x, digits_estimates,
                                               label = NULL){

  if(is.null(label)){
    label <- if(is.null(x$name)) "#1" else x$name
  }
  target <- .bt_random_variance_allocation_target(x)
  scale <- .bt_random_variance_allocation_scale(x)

  lines <- character()
  if(is.null(x$parent)){
    if(!is.null(x$sd)){
      lines <- c(lines, .bt_format_random_prior_equation(
        .bt_format_random_allocation_source_name(x),
        .bt_format_random_print_prior(x$sd, digits_estimates)
      ))
    }else if(!is.null(x$sd_source)){
      lines <- c(lines, .bt_format_random_prior_equation(
        .bt_format_random_allocation_source_name(x),
        .bt_random_sd_source_label(x$sd_source),
        operator = "="
      ))
    }
  }
  if(!is.null(x$weights)){
    lines <- c(lines, .bt_format_random_prior_equation(
      "w",
      .bt_format_random_allocation_weight_prior(x, digits_estimates)
    ))
  }
  lines <- c(lines, .bt_format_random_allocation_inclusion_lines(
    x,
    digits_estimates
  ))

  if(identical(target, "block")){
    lines <- c(lines, .bt_format_random_allocation_block_component_lines(x))
  }else{
    lines <- c(lines, .bt_format_random_allocation_sd_component_lines(x))
  }
  if(!identical(scale, "total_variance")){
    lines <- c(lines, paste0("scale: ", scale))
  }

  c(paste0("allocation: ", label), paste0("  ", lines))
}

.bt_format_prior_random_defaults <- function(x, digits_estimates){

  lines <- character()
  if(!is.null(x$sd)){
    lines <- c(lines, .bt_format_random_prior_equation(
      "sigma",
      .bt_format_random_print_prior(x$sd, digits_estimates)
    ))
  }
  lines <- c(lines, .bt_format_random_covariance_prior_lines(
    x$covariance,
    digits_estimates = digits_estimates,
    include_structure = FALSE,
    include_sd = TRUE
  ))

  has_structure <- !is.null(x$covariance) && !is.null(x$covariance$structure)
  if(length(lines) == 0L && !isTRUE(has_structure)){
    return(character())
  }

  header <- "defaults"
  if(isTRUE(has_structure)){
    header <- paste0(header, " (", x$covariance$structure, ")")
  }

  c(header, paste0("  ", lines))
}

.bt_format_prior_random_settings <- function(x){

  lines <- character()
  if(!.bt_random_monitor_is_default(x$monitor)){
    lines <- c(lines, .bt_format_random_monitor_settings(x$monitor))
  }
  if(!.bt_random_new_levels_is_default(x$new_levels)){
    lines <- c(lines, paste0("new_levels: ", x$new_levels$method))
  }
  if(!identical(x$parameterization, "noncentered")){
    lines <- c(lines, paste0("parameterization: ", x$parameterization))
  }
  if(length(lines) == 0L){
    return(character())
  }

  c("settings", paste0("  ", lines))
}

#' @export
print.prior_lkj <- function(x, digits_estimates = 2, silent = FALSE, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  if(!inherits(x, "prior_lkj")){
    stop("'x' must be created with prior_lkj().", call. = FALSE)
  }

  lines <- .bt_format_random_prior_equation(
    "R",
    .bt_format_prior_lkj_distribution(x, digits_estimates)
  )
  if(!isTRUE(x$include_correlation)){
    lines <- c(lines, "  include_correlation: FALSE")
  }
  if(isTRUE(x$include_primitives)){
    lines <- c(lines, "  include_primitives: TRUE")
  }

  .bt_print_spec_lines(lines, silent = silent)
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
    paste0("  method: ", x$method)
  ), silent = silent)
}

#' @export
print.random_covariance <- function(x, digits_estimates = 2, silent = FALSE,
                                    ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  .bt_check_random_covariance(x)

  .bt_print_spec_lines(.bt_format_random_covariance_prior_lines(
    x,
    digits_estimates = digits_estimates,
    include_structure = TRUE,
    include_sd = TRUE
  ), silent = silent)
}

#' @export
print.random_block <- function(x, digits_estimates = 2, silent = FALSE,
                               name = NULL, ...){

  check_int(digits_estimates, "digits_estimates", lower = 0)
  check_bool(silent, "silent")
  check_char(name, "name", allow_NULL = TRUE, allow_NA = FALSE)
  if(!inherits(x, "random_block")){
    stop("'x' must be created with random_block().", call. = FALSE)
  }

  .bt_print_spec_lines(.bt_format_random_block_prior_lines(
    x,
    digits_estimates = digits_estimates,
    name = name,
    include_empty = TRUE
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

  .bt_print_spec_lines(.bt_format_random_allocation_lines(
    x,
    digits_estimates = digits_estimates
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

  lines <- character()
  lines <- c(lines, .bt_format_prior_random_defaults(x, digits_estimates))

  allocations <- .bt_random_allocation_list(x$allocation)
  if(length(allocations) > 0L){
    for(i in seq_along(allocations)){
      lines <- c(lines, .bt_format_random_allocation_lines(
        allocations[[i]],
        digits_estimates = digits_estimates,
        label = .bt_format_random_allocation_label(allocations, i)
      ))
    }
  }

  if(length(x$blocks) > 0L){
    block_names <- names(x$blocks)
    for(i in seq_along(x$blocks)){
      lines <- c(lines, .bt_format_random_block_prior_lines(
        x$blocks[[i]],
        digits_estimates = digits_estimates,
        name = block_names[[i]],
        include_empty = TRUE
      ))
    }
  }

  lines <- c(lines, .bt_format_prior_random_settings(x))
  if(length(lines) == 0L){
    lines <- "no random-effect priors specified"
  }

  .bt_print_spec_lines(lines, silent = silent)
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
    paste0("  source: ", .bt_random_sd_source_label(x))
  ), silent = silent)
}
