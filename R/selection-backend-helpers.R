.selection_null_default <- function(x, default){

  if(is.null(x)){
    return(default)
  }

  return(x)
}

.phack_validate_alpha_prior <- function(alpha){

  .check_prior(alpha)

  if(!is.prior.simple(alpha)){
    stop("'alpha' must be a simple prior distribution.", call. = FALSE)
  }
  if(is.prior.discrete(alpha) && !is.prior.point(alpha)){
    stop("'alpha' must be a continuous prior distribution or a point prior.", call. = FALSE)
  }

  if(is.prior.point(alpha)){
    if(alpha$parameters[["location"]] < 0 || alpha$parameters[["location"]] >= 1){
      stop("Point p-hacking alpha priors must be in the interval [0, 1).", call. = FALSE)
    }
  }else{
    if(alpha$truncation[["lower"]] < 0 || alpha$truncation[["upper"]] > 1){
      stop("P-hacking alpha priors must have support within [0, 1].", call. = FALSE)
    }
  }

  return()
}

.phack_kind <- function(form){
  switch(
    form,
    "linear"    = 1L,
    "quadratic" = 2L
  )
}

# Standard-normal probability of (lower, upper). pnorm(upper) - pnorm(lower)
# cancels to 0 (or a negative value) once both cut points lie far in the upper
# tail, which is where small p-value cuts put them; an interval in the upper
# half is therefore a difference of upper-tail probabilities, and one in the
# lower half a difference of lower-tail probabilities.
.phack_normal_interval_mass <- function(lower, upper){

  if(lower >= 0){
    return(stats::pnorm(lower, lower.tail = FALSE) - stats::pnorm(upper, lower.tail = FALSE))
  }
  if(upper <= 0){
    return(stats::pnorm(upper) - stats::pnorm(lower))
  }

  1 - stats::pnorm(lower) - stats::pnorm(upper, lower.tail = FALSE)
}

.phack_power_null_moment <- function(lower, upper, q, anchor, reverse){

  width <- upper - lower
  m0 <- .phack_normal_interval_mass(lower, upper)
  m1 <- stats::dnorm(lower) - stats::dnorm(upper)

  if(q == 1L){
    if(reverse){
      return((anchor * m0 - m1) / width)
    }else{
      return((m1 - anchor * m0) / width)
    }
  }

  m2 <- m0 + lower * stats::dnorm(lower) - upper * stats::dnorm(upper)
  if(reverse){
    return((m2 - 2 * anchor * m1 + anchor^2 * m0) / width^2)
  }

  return((m2 - 2 * anchor * m1 + anchor^2 * m0) / width^2)
}

.selection_backend_names <- function(names){

  backend_names <- names
  backend_name_names <- base::names(backend_names)
  if(is.null(backend_name_names) ||
     any(!nzchar(backend_name_names))){
    stop("All entries in the 'names' argument must be named.", call. = FALSE)
  }

  defaults <- list(
    omega               = "omega",
    alpha               = "alpha",
    pi_null             = "pi_null",
    beta_null           = "beta_null",
    phack_kind          = "phack_kind",
    phack_z_source      = "phack_z_source",
    phack_z_dest        = "phack_z_dest"
  )
  defaults[base::names(backend_names)] <- backend_names

  for(i in seq_along(defaults)){
    check_char(defaults[[i]], paste0("names$", names(defaults)[i]))
  }

  # The names become JAGS nodes written by the generated code, so they must be
  # distinct JAGS identifiers that no generated internal node uses.
  node_names <- unlist(defaults, use.names = FALSE)
  if(!all(grepl("^[A-Za-z][A-Za-z0-9._]*$", node_names))){
    stop("All entries in the 'names' argument must be valid JAGS node names.", call. = FALSE)
  }
  if(anyDuplicated(node_names)){
    stop("All entries in the 'names' argument must be distinct.", call. = FALSE)
  }
  reserved <- c("bias_indicator", "sel_vector_rule", "omega_local", "omega_ratio", "log_omega", "eta", "std_eta")
  if(any(node_names %in% reserved | grepl("_component_[0-9]+$", node_names))){
    stop(
      "The 'names' argument must not use the internal selection node names: ",
      paste0("'", reserved, "'", collapse = ", "), " or '*_component_<k>'.",
      call. = FALSE
    )
  }

  return(defaults)
}

# Renames the indexed JAGS node `from` to `to` in generated syntax. A single
# selection branch writes its public nodes directly, and the weight-function
# syntax names that node by its default name.
.selection_rename_jags_node <- function(code, from, to){

  if(identical(from, to) || length(code) == 0L){
    return(code)
  }

  gsub(paste0("(?<![A-Za-z0-9._])", from, "(?=\\[)"), to, code, perl = TRUE)
}

.selection_normalize_priors <- function(priors){

  if(is.null(priors)){
    priors <- list(prior_none())
  }else if(is.prior(priors) && is.prior.mixture(priors)){
    priors <- as.list(priors)
  }else if(is.prior(priors)){
    priors <- list(priors)
  }else if(is.list(priors)){
    if(length(priors) == 0L){
      priors <- list(prior_none())
    }
  }else{
    stop("'priors' must be a prior, prior mixture, list of priors, or NULL.", call. = FALSE)
  }

  if(!all(vapply(priors, is.prior, logical(1)))){
    stop("'priors' must contain only prior objects.", call. = FALSE)
  }

  allowed <- vapply(priors, function(x){
    is.prior.none(x) || is.prior.PET(x) || is.prior.PEESE(x) ||
      is.prior.weightfunction(x) || is_prior_phacking(x) || is_prior_bias(x)
  }, logical(1))

  if(!all(allowed)){
    stop("'priors' contains unsupported selection prior objects.", call. = FALSE)
  }

  return(priors)
}

.selection_branch_info <- function(prior){

  if(is.prior.none(prior)){
    return(list(type = "none", selection = NULL, phacking = NULL))
  }
  if(is.prior.PET(prior)){
    return(list(type = "PET", selection = NULL, phacking = NULL))
  }
  if(is.prior.PEESE(prior)){
    return(list(type = "PEESE", selection = NULL, phacking = NULL))
  }
  if(is.prior.weightfunction(prior)){
    return(list(type = "weightfunction", selection = prior, phacking = NULL))
  }
  if(is_prior_phacking(prior)){
    return(list(type = "phack", selection = NULL, phacking = prior))
  }
  if(is_prior_bias(prior)){
    has_selection <- !is.null(prior$selection)
    has_phacking  <- !is.null(prior$phacking)
    type <- if(has_selection && has_phacking){
      "combined"
    }else if(has_selection){
      "weightfunction"
    }else{
      "phack"
    }
    return(list(type = type, selection = prior$selection, phacking = prior$phacking))
  }

  stop("Unsupported selection prior object.", call. = FALSE)
}

.selection_backend_mode <- function(has_selection, has_phacking){

  if(has_selection && has_phacking){
    return("step_phack_power")
  }
  if(has_selection){
    return("step")
  }
  if(has_phacking){
    return("phack_power")
  }

  return("none")
}

.selection_mode_code <- function(mode){
  switch(
    mode,
    "none"             = 0L,
    "step"             = 1L,
    "phack_power"      = 2L,
    "step_phack_power" = 3L
  )
}

.selection_validate_global_breaks <- function(global_breaks){

  check_real(global_breaks, "global_breaks", check_length = 0, allow_NA = FALSE)
  if(length(global_breaks) < 2L){
    stop("'global_breaks' must contain at least 0 and 1.", call. = FALSE)
  }
  if(any(global_breaks < 0 | global_breaks > 1)){
    stop("'global_breaks' must be within [0, 1].", call. = FALSE)
  }
  if(anyDuplicated(global_breaks)){
    stop("'global_breaks' must not contain duplicate values.", call. = FALSE)
  }
  if(!all(global_breaks == cummax(global_breaks))){
    stop("'global_breaks' must be monotonically increasing.", call. = FALSE)
  }
  if(global_breaks[1] != 0 || global_breaks[length(global_breaks)] != 1){
    stop("'global_breaks' must start at 0 and end at 1.", call. = FALSE)
  }

  return(global_breaks)
}

.selection_jags_step_component_code <- function(selection, component_id, n_bins, global_cuts){

  if(is.null(selection)){
    return(.JAGS_weightfunction_none_component_syntax(component_id = component_id, n_bins = n_bins))
  }

  return(.JAGS_weightfunction_component_syntax(
    prior           = selection,
    component_id    = component_id,
    global_cuts     = global_cuts,
    force_one_sided = TRUE
  ))
}

.selection_jags_active_scalar <- function(prefix, target, indicator_terms, component_ids){

  paste0(
    target, " <- ",
    paste0(prefix, component_ids, " * ", indicator_terms, collapse = " + ")
  )
}

.selection_jags_active_vector <- function(prefix, target, index, indicator_terms, component_ids){

  paste0(
    target, "[", index, "] <- ",
    paste0(prefix, component_ids, "[", index, "] * ", indicator_terms, collapse = " + ")
  )
}

.selection_backend_phacking_info <- function(phacking_priors, names){

  if(length(phacking_priors) == 0L){
    return(list(
      form                 = "none",
      q                    = 0L,
      z_source             = c(0, 0),
      z_destination        = c(0, 0),
      coefficient          = names$alpha,
      coefficient_ids      = names$alpha,
      branch_form          = character(),
      branch_q             = integer(),
      branch_phack_kind    = integer(),
      branch_beta_null_per_alpha = numeric(),
      branch_z_source      = matrix(nrow = 0L, ncol = 2L),
      branch_z_destination = matrix(nrow = 0L, ncol = 2L)
    ))
  }

  constants <- lapply(phacking_priors, function(x){
    phack_backend_constants(x$form, x$source, x$destination, target = x$target)
  })

  forms <- vapply(constants, function(x) x$form, character(1))
  q <- vapply(constants, function(x) x$q, integer(1))
  z_source <- do.call(rbind, lapply(constants, function(x) x$z_source))
  z_destination <- do.call(rbind, lapply(constants, function(x) x$z_destination))

  return(list(
    form                 = if(length(unique(forms)) == 1L) forms[1] else unique(forms),
    q                    = if(length(unique(q)) == 1L) q[1] else unique(q),
    z_source             = if(nrow(z_source) == 1L) as.numeric(z_source[1,]) else z_source,
    z_destination        = if(nrow(z_destination) == 1L) as.numeric(z_destination[1,]) else z_destination,
    coefficient          = names$alpha,
    coefficient_ids      = names$alpha,
    branch_form          = forms,
    branch_q             = q,
    branch_phack_kind    = vapply(constants, function(x) x$phack_kind, integer(1)),
    branch_beta_null_per_alpha = vapply(constants, function(x) x$beta_null_per_alpha, numeric(1)),
    branch_z_source      = z_source,
    branch_z_destination = z_destination
  ))
}

.selection_backend_init <- function(branch_info, prior_weights, uses_indicator, global_cuts = NULL, names = NULL){

  active_branch <- which.max(prior_weights)

  init <- list()
  for(i in seq_along(branch_info)){
    component_id <- if(uses_indicator) i else NULL
    if(!is.null(branch_info[[i]]$selection)){
      selection_init <- .selection_JAGS_init_weightfunction_component(branch_info[[i]]$selection, component_id = component_id, global_cuts = global_cuts)
      if(is.null(component_id) && !is.null(names) && "omega" %in% base::names(selection_init)){
        # A single unmapped independent weight function samples the public
        # omega node itself, which the syntax writes under names$omega.
        base::names(selection_init)[base::names(selection_init) == "omega"] <- names$omega
      }
      init <- c(init, selection_init)
    }
    if(!is.null(branch_info[[i]]$phacking)){
      init <- c(init, .selection_JAGS_init_phacking_component(branch_info[[i]]$phacking, component_id = component_id, names = names))
    }
  }
  if(uses_indicator){
    init[["bias_indicator"]] <- active_branch
  }

  return(init)
}

.selection_JAGS_init_weightfunction_component <- function(prior, component_id = NULL, global_cuts = NULL){

  init <- list()
  # Same resolver as the model syntax: initial values must name the stochastic
  # node, never the deterministic `omega_target` expanded onto the global cuts.
  node_names <- .weightfunction_component_node_names(
    prior           = prior,
    component_id    = component_id,
    global_cuts     = global_cuts,
    force_one_sided = TRUE
  )
  if(prior$weights$type == "fixed"){
    return()
  }else if(prior$weights$type == "cumulative"){
    label <- if(is.null(component_id)){
      "cumulative weight function"
    }else{
      paste0("cumulative weight-function component '", component_id, "'")
    }
    if(node_names$n_bins == 2L){
      init[[node_names$omega_ratio]] <- .JAGS_binary_cumulative_initialization(
        alpha = prior$weights[["alpha"]],
        label = label
      )
    }else{
      init[[node_names$eta]] <- .JAGS_positive_gamma_initialization(
        shape = prior$weights[["alpha"]],
        label = label
      )
    }
  }else if(prior$weights$type == "independent"){
    n_bins <- node_names$n_bins
    if(n_bins > 1L){
      n_free <- n_bins - 1L
      draws <- as.numeric(rng(prior$weights$prior, n_free))
      if(length(draws) != n_free || any(!is.finite(draws))){
        stop(
          "Independent weight-function initialization produced a non-finite draw.",
          call. = FALSE
        )
      }
      values <- rep(NA_real_, n_bins)
      values[seq.int(2L, n_bins)] <- draws
      node_name <- if(identical(prior$weights$scale, "log_omega")){
        node_names$log_omega
      }else{
        node_names$omega_local
      }
      init[[node_name]] <- values
    }
  }

  return(init)
}

.selection_JAGS_init_phacking_component <- function(prior, component_id = NULL, names = NULL){

  .JAGS_init.simple(prior$alpha, .selection_phacking_node_name("alpha", component_id, names))
}

# JAGS node of a p-hacking field: `<field>_component_<k>` inside a mixture,
# otherwise the public name from `names` (default: the field name).
.selection_phacking_node_name <- function(field, component_id = NULL, names = NULL){

  if(!is.null(component_id)){
    return(paste0(field, "_component_", component_id))
  }
  if(!is.null(names[[field]])){
    return(names[[field]])
  }

  field
}

.selection_prior_branch_info <- function(prior){

  if(is.prior.mixture(prior)){
    return(lapply(prior, .selection_branch_info))
  }

  list(.selection_branch_info(prior))
}

.selection_prior_has_selection <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  any(vapply(branch_info, function(x) !is.null(x$selection), logical(1)))
}

.selection_prior_has_phacking <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  any(vapply(branch_info, function(x) !is.null(x$phacking), logical(1)))
}

.selection_prior_selection_priors <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  selection_priors <- lapply(branch_info, function(x) x$selection)
  selection_priors[!vapply(selection_priors, is.null, logical(1))]
}

.selection_prior_phacking_priors <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  phacking_priors <- lapply(branch_info, function(x) x$phacking)
  phacking_priors[!vapply(phacking_priors, is.null, logical(1))]
}

.selection_prior_stop_unsupported_generic <- function(generic, prior){

  if(is_prior_bias(prior)){
    selection_backend_spec(prior, include_init = FALSE)
    stop(
      sprintf(
        "No %s is implemented for composed bias priors; use the selection or p-hacking component explicitly.",
        generic
      ),
      call. = FALSE
    )
  }

  if(is_prior_phacking(prior)){
    selection_backend_spec(prior, include_init = FALSE)
    stop(
      sprintf(
        "No %s is implemented for p-hacking priors; use the alpha prior component explicitly.",
        generic
      ),
      call. = FALSE
    )
  }

  stop(sprintf("No %s is implemented for this prior.", generic), call. = FALSE)
}

.selection_prior_has_PET <- function(prior){

  if(is.prior.mixture(prior)){
    return(any(sapply(prior, is.prior.PET)))
  }

  is.prior.PET(prior)
}

.selection_prior_has_PEESE <- function(prior){

  if(is.prior.mixture(prior)){
    return(any(sapply(prior, is.prior.PEESE)))
  }

  is.prior.PEESE(prior)
}

.selection_bias_parameter_names <- function(prior, include_kind = TRUE){

  c(
    if(.selection_prior_has_PET(prior)) "PET",
    if(.selection_prior_has_PEESE(prior)) "PEESE",
    if(.selection_prior_has_selection(prior)) "omega",
    if(.selection_prior_has_phacking(prior)) c(
      .selection_prior_phacking_report_parameters(prior),
      if(include_kind) c("beta_null", "phack_kind")
    )
  )
}

.selection_format_number <- function(x){
  format(x, scientific = FALSE, digits = 16, trim = TRUE)
}

.JAGS_phacking_component_syntax <- function(prior, component_id = NULL, names = NULL){

  alpha_name         <- .selection_phacking_node_name("alpha", component_id, names)
  kind_name          <- .selection_phacking_node_name("phack_kind", component_id, names)
  pi_null_name       <- .selection_phacking_node_name("pi_null", component_id, names)
  beta_null_name     <- .selection_phacking_node_name("beta_null", component_id, names)
  z_source_name      <- .selection_phacking_node_name("phack_z_source", component_id, names)
  z_destination_name <- .selection_phacking_node_name("phack_z_dest", component_id, names)

  if(is.null(prior)){
    return(paste0(
      alpha_name, " <- 0\n",
      kind_name, " <- 0\n",
      pi_null_name, " <- 0\n",
      beta_null_name, " <- 0\n",
      z_source_name, "[1] <- 0\n",
      z_source_name, "[2] <- 0\n",
      z_destination_name, "[1] <- 0\n",
      z_destination_name, "[2] <- 0\n"
    ))
  }

  if(!is_prior_phacking(prior)){
    stop("'prior' must be a p-hacking prior.", call. = FALSE)
  }

  constants <- phack_backend_constants(prior$form, prior$source, prior$destination, target = prior$target)

  paste0(
    .JAGS_prior.simple(prior$alpha, alpha_name),
    kind_name, " <- ", constants$phack_kind, "\n",
    pi_null_name, " <- ", alpha_name, " * ", .selection_format_number(constants$pi_null_per_alpha), "\n",
    beta_null_name, " <- ", alpha_name, " * ", .selection_format_number(constants$beta_null_per_alpha), "\n",
    z_source_name, "[1] <- ", .selection_format_number(constants$z_source[1]), "\n",
    z_source_name, "[2] <- ", .selection_format_number(constants$z_source[2]), "\n",
    z_destination_name, "[1] <- ", .selection_format_number(constants$z_destination[1]), "\n",
    z_destination_name, "[2] <- ", .selection_format_number(constants$z_destination[2]), "\n"
  )
}
