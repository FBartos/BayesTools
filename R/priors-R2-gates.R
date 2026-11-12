# Inclusion gates and spike-and-slab totals of R-squared budgets.
#
# A gated or spiked budget is one model-averaging fit over the models of its
# terms. A configuration S of open terms (the open gated terms and the
# ungated terms) is a model, and given S the coefficients follow the prior of
# that model as a separate fit: the allocation is redistributed over the open
# terms (a kernel gate closes every unit of its term), and the total has the
# Beta prior of the convention for S (common: Beta(a, b); linked: Beta(A_S, b),
# sampled as term signals g_t eta_t / H; conditioning: Beta(a - (A - A_S),
# b + (A - A_S))), with the concentrations rescaled per S under fixed_total.
# The empty configuration is the budget null: no signal from the budget.
#
# The gate prior is exact: the gates are Bernoulli nodes sampled in sequence,
# each with its conditional probability given the gates before it, so that
# the joint law of the gates is the declared model prior (independent fixed
# probabilities, or the beta-binomial law of a shared Beta rate), restricted
# to nonempty configurations when a spike-and-slab total makes the null its
# own event (P(S empty) = the spike weight; under the slab P(S) is the gate
# prior conditioned on S nonempty). No configuration is enumerated, so the
# number of gated terms is not limited.

# The gate structure of 'budget' in a model that retains 'budget_terms' (the
# prior-list names of its retained terms, with their declared labels
# 'declared'); NULL for a budget without gates and spike.
.bt_R2_gate_structure <- function(budget, budget_terms, declared){

  spike <- budget$total_inclusion
  gated <- if(is.null(budget$inclusion)){
    character()
  }else if(identical(budget$inclusion$form, "common")){
    budget_terms
  }else{
    budget_terms[.bt_R2_term_key(declared) %in% .bt_R2_term_key(names(budget$inclusion$laws))]
  }
  if(length(gated) == 0L && is.null(spike)){
    return(NULL)
  }
  K <- length(gated)
  law <- if(K == 0L){
    NULL
  }else if(identical(budget$inclusion$form, "common") && identical(budget$inclusion$law$type, "beta")){
    list(type = "shared", a = budget$inclusion$law$a, b = budget$inclusion$law$b,
         pi = rep(budget$inclusion$law$p, K))
  }else if(identical(budget$inclusion$form, "common")){
    list(type = "fixed", pi = rep(budget$inclusion$law$p, K))
  }else{
    gated_declared <- declared[match(gated, budget_terms)]
    laws <- budget$inclusion$laws[match(.bt_R2_term_key(gated_declared), .bt_R2_term_key(names(budget$inclusion$laws)))]
    list(type = "fixed", pi = unname(vapply(laws, `[[`, numeric(1), "p")))
  }
  all_gated <- K == length(budget_terms)
  truncated <- !is.null(spike) && all_gated && K > 0L

  list(
    gated     = gated,
    K         = K,
    law       = law,
    spike     = spike,
    all_gated = all_gated,
    truncated = truncated,
    p_trunc   = if(truncated) .bt_R2_gate_truncated_probabilities(law, K)
  )
}

# The log probability that the gates t..K are all closed given that the gates
# before t are closed, under the gate law before the spike (t = 1: the log
# probability of the empty configuration). Under a shared Beta(a, b) rate it is
# the product over the closed gates j = t - 1, ..., K - 1 of the beta-binomial
# conditionals (b + j) / (a + b + j) = 1 - a / (a + b + j), summed in logs
# without differences of large lbeta() values (which lose huge shapes).
.bt_R2_gate_log_closed_from <- function(law, K, t = 1L){

  if(t > K){
    return(0)
  }
  if(identical(law$type, "shared")){
    j <- (t - 1):(K - 1)
    return(sum(.bt_R2_log_fraction(law$b + j, law$a)))
  }
  sum(log1p(-law$pi[t:K]))
}

# log(x / (x + y)) for x, y > 0 without cancellation and without forming
# x + y: log(x / y) - log1p(x / y) when the fraction is small, -log1p(y / x)
# when it is near one (both stay accurate for very unequal x and y).
.bt_R2_log_fraction <- function(x, y){
  ifelse(x < y, log(x) - log(y) - log1p(x / y), -log1p(y / x))
}

# The probability of a nonempty configuration of the gates before the spike,
# 1 - P(empty), computed without cancellation (tiny gate probabilities keep
# their value).
.bt_R2_gate_nonempty_probability <- function(gates){
  -expm1(.bt_R2_gate_log_closed_from(gates$law, gates$K))
}

# P(gate t open | gates 1..t-1 closed, at least one gate open), t = 1..K.
.bt_R2_gate_truncated_probabilities <- function(law, K){

  p <- vapply(seq_len(K), function(t){
    first <- if(identical(law$type, "shared")) law$a / (law$a + law$b + (t - 1)) else law$pi[[t]]
    # P(some of t..K open | none before), without cancellation
    first / -expm1(.bt_R2_gate_log_closed_from(law, K, t))
  }, numeric(1))
  if(any(!is.finite(p)) || any(p > 1 + 1e-10) || any(p <= 0)){
    .bt_R2_stop(paste0("The inclusion probabilities of the gates are numerically degenerate: the conditional ",
      "probabilities of the gates given a nonempty configuration cannot be computed."))
  }
  # with every gate before it closed, the last gate is open (exactly; the
  # ratio above rounds)
  p[[K]] <- 1

  pmin(1, p)
}

# The JAGS expression of P(gate t open | s open gates before it), with 's' the
# JAGS expression of that count ("0" for t = 1).
.bt_R2_gate_probability_syntax <- function(gates, t, s){

  regular <- if(identical(gates$law$type, "shared")){
    paste0("(", .JAGS_prior_number(gates$law$a), " + ", s, ") / ",
           .JAGS_prior_number(gates$law$a + gates$law$b + (t - 1)))
  }else{
    .JAGS_prior_number(gates$law$pi[[t]])
  }
  if(!gates$truncated){
    return(regular)
  }
  if(t == 1L){
    return(.JAGS_prior_number(gates$p_trunc[[1L]]))
  }

  paste0("ifelse(", s, " > 0, ", regular, ", ", .JAGS_prior_number(gates$p_trunc[[t]]), ")")
}

# The R-side conditional gate probability (the same law as the syntax).
.bt_R2_gate_probability <- function(gates, t, s){

  regular <- if(identical(gates$law$type, "shared")){
    (gates$law$a + s) / (gates$law$a + gates$law$b + (t - 1))
  }else{
    rep(gates$law$pi[[t]], length(s))
  }
  if(gates$truncated){
    regular[s == 0] <- gates$p_trunc[[t]]
  }

  regular
}

# The probability that gate t is closed given 's' open gates before it, before
# the spike and the truncation: under a shared rate (b + t - 1 - s) / (a + b +
# t - 1), not 1 minus the open probability (which cancels when it is near 1).
.bt_R2_gate_closed_probability <- function(gates, t, s){

  if(identical(gates$law$type, "shared")){
    return((gates$law$b + (t - 1 - s)) / (gates$law$a + gates$law$b + (t - 1)))
  }
  rep(1 - gates$law$pi[[t]], length(s))
}

# The log prior probability of the gate configurations 'open' (a logical
# matrix, one row per configuration of the K gates) under the gate law,
# before the spike and the truncation: the probability that exactly these
# gates are open. Under a shared rate it depends only on the number s of open
# gates (exchangeable), the product of the sequential beta-binomial
# conditionals with the open gates first, summed in logs without differences
# of large lbeta() values.
.bt_R2_gate_log_configuration_probability <- function(gates, open){

  open <- matrix(as.logical(open), ncol = gates$K)
  if(!identical(gates$law$type, "shared")){
    log_open <- matrix(log(gates$law$pi), nrow(open), gates$K, byrow = TRUE)
    log_closed <- matrix(log1p(-gates$law$pi), nrow(open), gates$K, byrow = TRUE)
    return(rowSums(ifelse(open, log_open, log_closed)))
  }
  a <- gates$law$a
  b <- gates$law$b
  K <- gates$K
  s <- rowSums(open)
  by_s <- vapply(0:K, function(count){
    i <- seq_len(count) - 1
    j <- seq_len(K - count) - 1
    sum(.bt_R2_log_fraction(a + i, b)) + sum(.bt_R2_log_fraction(b + j, a + count))
  }, numeric(1))

  by_s[s + 1L]
}

# Prior probability of the gate configurations 'open' (see above), with the
# truncation of a spiked budget whose terms are all gated (the empty
# configuration then has probability zero under the slab).
.bt_R2_gate_configuration_probability <- function(gates, open){

  open <- matrix(as.logical(open), ncol = gates$K)
  base <- exp(.bt_R2_gate_log_configuration_probability(gates, open))
  if(gates$truncated){
    base <- ifelse(rowSums(open) == 0, 0, base / .bt_R2_gate_nonempty_probability(gates))
  }

  base
}

# Marginal prior probabilities of the budget being on (some signal) and of
# each retained term being open (named by the prior-list term names).
.bt_R2_gate_marginals <- function(gates, budget_terms){

  on_spike <- if(is.null(gates$spike)) 1 else gates$spike$p
  open <- stats::setNames(rep(on_spike, length(budget_terms)), budget_terms)
  if(gates$K > 0L){
    pi <- gates$law$pi
    if(gates$truncated){
      pi <- pi / .bt_R2_gate_nonempty_probability(gates)
    }
    open[gates$gated] <- on_spike * pi
  }
  budget <- if(gates$K > 0L && gates$all_gated && is.null(gates$spike)){
    .bt_R2_gate_nonempty_probability(gates)
  }else{
    on_spike
  }

  list(budget = budget, open = open)
}

# The core prior of a gated or spiked budget: the prior-list entry of its
# total, which emits the gates, the budget indicator, the allocation
# coordinates, the total and the gated allocation (see the header).
.bt_R2_core_prior <- function(budget, gates, parameter, budget_terms, unit_term, unit_alpha,
                              unit_reference, shapes){

  total_name <- paste0(parameter, "_", budget$name)
  n_units <- length(unit_alpha)
  core <- list(
    distribution    = "R2_budget",
    parameters      = list(),
    truncation      = list(lower = 0, upper = 1),
    prior_weights   = 1,
    name            = budget$name,
    parameter       = parameter,
    total_name      = total_name,
    allocation_name = if(n_units > 1L) paste0(total_name, "_allocation") else NA_character_,
    law             = budget$law,
    convention      = budget$convention,
    concentration   = budget$concentration,
    a               = budget$a,
    b               = budget$b,
    A_reference     = unit_reference,
    fixed_total     = sum(budget$alpha),
    unit_alpha      = unit_alpha,
    unit_term       = unit_term,
    n_units         = n_units,
    terms           = budget_terms,
    gates           = gates,
    term_gate       = stats::setNames(match(budget_terms, gates$gated, nomatch = 0L), budget_terms),
    shapes          = shapes
  )
  class(core) <- c("prior", "prior.R2_budget")
  attr(core, "R2_budget") <- list(name = budget$name, role = "total", law = budget$law, gated = TRUE)

  core
}

is.prior.R2_budget <- function(x){
  inherits(x, "prior.R2_budget")
}

# Node names of a core.
.bt_R2_core_nodes <- function(core){

  T <- core$total_name
  list(
    total      = T,
    variable   = paste0(T, "_variable"),
    indicator  = if(!is.null(core$gates$spike) || (core$gates$K > 0L && core$gates$all_gated)) paste0(T, "_indicator"),
    inclusion  = if(identical(core$gates$spike$type, "beta")) paste0(T, "_inclusion"),
    # (JAGS names a vector of length one by its scalar name: one gate or one
    # unit is a scalar node)
    gates      = if(core$gates$K > 0L) .bt_dnode_prior_vector_names(paste0(T, "_gate"), core$gates$K),
    open       = if(core$gates$K > 0L) .bt_dnode_prior_vector_names(paste0(T, "_open"), core$gates$K),
    allocation = if(core$n_units > 1L) paste0(core$allocation_name, "[", seq_len(core$n_units), "]"),
    eta        = .bt_dnode_prior_vector_names(.JAGS_prior_dirichlet_eta_name(
      if(core$n_units > 1L) core$allocation_name else paste0(T, "_allocation")), core$n_units),
    H          = if(identical(core$convention, "linked")) paste0("prior_par_H_", T),
    rate       = if(.bt_R2_core_has_rate(core)) paste0(T, "_rate"),
    rate_u     = if(.bt_R2_core_has_rate(core)) paste0(T, "_rate_u")
  )
}

# Whether a core has a shared Beta inclusion rate theta, a node drawn given the
# gates: theta | S ~ Beta(a + s, b + K - s), with s the number of open gates.
# It is the quantile of an independent uniform at that Beta law, so it does not
# feed back into the collapsed sampler of the gates (whose law integrates
# theta out); under a spike-and-slab total the gates keep their slab law when
# the spike is selected, and so does theta.
.bt_R2_core_has_rate <- function(core){
  core$gates$K > 0L && identical(core$gates$law$type, "shared")
}

# Open flag of every unit as a JAGS expression ("1" for an always open unit).
.bt_R2_core_unit_open_syntax <- function(core, nodes){

  vapply(seq_len(core$n_units), function(u){
    gate <- core$term_gate[[core$unit_term[[u]]]]
    factors <- c(
      if(!is.null(core$gates$spike)) nodes$indicator,
      if(gate > 0L) nodes$gates[[gate]]
    )
    if(length(factors) == 0L) "1" else paste(factors, collapse = " * ")
  }, character(1))
}

.bt_R2_core_syntax <- function(core){

  nodes <- .bt_R2_core_nodes(core)
  gates <- core$gates
  num <- .JAGS_prior_number
  lines <- character()
  # the number of open gates among the first k (a single gate is a scalar)
  gate_sum <- function(k) if(gates$K == 1L) nodes$gates[[1L]] else paste0("sum(", core$total_name, "_gate[1:", k, "])")

  # the gates, each given the number of open gates before it
  for(t in seq_len(gates$K)){
    s <- if(t == 1L) "0" else gate_sum(t - 1L)
    lines <- c(lines, paste0(nodes$gates[[t]], " ~ dbern(", .bt_R2_gate_probability_syntax(gates, t, s), ")"))
  }
  # the budget indicator: the spike-and-slab total, or any gate open
  if(!is.null(gates$spike)){
    if(identical(gates$spike$type, "beta")){
      lines <- c(lines, paste0(nodes$inclusion, " ~ dbeta(", num(gates$spike$a), ", ", num(gates$spike$b), ")"),
                 paste0(nodes$indicator, " ~ dbern(", nodes$inclusion, ")"))
    }else{
      lines <- c(lines, paste0(nodes$indicator, " ~ dbern(", num(gates$spike$p), ")"))
    }
  }else if(!is.null(nodes$indicator)){
    lines <- c(lines, paste0(nodes$indicator, " <- step(", gate_sum(gates$K), " - 0.5)"))
  }

  # the shared inclusion rate given the gates
  if(!is.null(nodes$rate)){
    s_all <- gate_sum(gates$K)
    lines <- c(lines, paste0(nodes$rate_u, " ~ dunif(0, 1)"),
               paste0(nodes$rate, " <- qbeta(", nodes$rate_u, ", ", num(gates$law$a), " + ", s_all, ", ",
                      num(gates$law$b), " + (", gates$K, " - ", s_all, "))"))
  }
  # the open flag of every gated term (its gate times the budget indicator of
  # a spike-and-slab total): the inclusion of the term
  for(t in seq_len(gates$K)){
    lines <- c(lines, paste0(nodes$open[[t]], " <- ", if(!is.null(gates$spike)) paste0(nodes$indicator, " * "), nodes$gates[[t]]))
  }
  open <- .bt_R2_core_unit_open_syntax(core, nodes)
  open_alpha <- paste0("(", paste0(open, " * ", vapply(core$unit_alpha, num, character(1)), collapse = " + "), ")")
  any_open <- paste0("step(", paste0(open, collapse = " + "), " - 0.5)")

  # allocation gamma coordinates (shapes rescaled per configuration under
  # fixed_total; a closed unit keeps its declared shape)
  for(u in seq_len(core$n_units)){
    shape <- if(identical(core$concentration, "fixed_total")){
      paste0("ifelse(", open[[u]], " > 0, ", num(core$unit_alpha[[u]] * core$fixed_total), " / (", open_alpha,
             " + equals(", open_alpha, ", 0)), ", num(core$unit_alpha[[u]]), ")")
    }else{
      num(core$unit_alpha[[u]])
    }
    lines <- c(lines, paste0(nodes$eta[[u]], " ~ dgamma(", shape, ", 1)"))
  }
  open_eta <- paste0("(", paste0(open, " * ", nodes$eta, collapse = " + "), ")")
  if(core$n_units > 1L){
    for(u in seq_len(core$n_units)){
      lines <- c(lines, paste0(nodes$allocation[[u]], " <- ", open[[u]], " * ", nodes$eta[[u]], " / (",
                               open_eta, " + equals(", open_eta, ", 0))"))
    }
  }

  # the slab total of the open configuration
  A_open <- if(identical(core$concentration, "fixed_total")) num(core$fixed_total) else open_alpha
  if(identical(core$convention, "linked")){
    lines <- c(lines, paste0(nodes$H, " ~ dgamma(", num(core$b), ", 1)"),
               paste0(nodes$variable, " <- ", open_eta, " / (", open_eta, " + ", nodes$H, ")"))
  }else if(identical(core$convention, "conditioning")){
    A <- if(identical(core$concentration, "fixed_total")) num(core$fixed_total) else num(core$A_reference)
    lines <- c(lines, paste0(nodes$variable, " ~ dbeta(ifelse(", any_open, " > 0, ", num(core$a), " - (", A, " - ", A_open,
      "), ", num(core$a), "), ifelse(", any_open, " > 0, ", num(core$b), " + (", A, " - ", A_open, "), ", num(core$b), "))"))
  }else{
    lines <- c(lines, paste0(nodes$variable, " ~ dbeta(", num(core$a), ", ", num(core$b), ")"))
  }
  # the realized total: no signal from the budget when no term is open
  lines <- c(lines, paste0(nodes$total, " <- ", nodes$variable, " * ", any_open))

  paste0(paste0(lines, "\n", collapse = ""))
}

.bt_R2_core_monitor <- function(core){

  nodes <- .bt_R2_core_nodes(core)
  c(nodes$total, nodes$indicator, nodes$inclusion, nodes$rate,
    if(core$gates$K > 0L) paste0(core$total_name, c("_gate", "_open")),
    if(core$n_units > 1L) core$allocation_name,
    sub("\\[.*$", "", nodes$eta[[1L]]), nodes$H,
    if(!identical(core$convention, "linked")) nodes$variable)
}

# Initial values: every gate and the budget open, slab draws of the full
# configuration.
.bt_R2_core_inits <- function(core){

  nodes <- .bt_R2_core_nodes(core)
  init <- list()
  if(core$gates$K > 0L){
    init[[paste0(core$total_name, "_gate")]] <- rep(1, core$gates$K)
  }
  if(!is.null(nodes$rate_u)){
    init[[nodes$rate_u]] <- stats::runif(1)
  }
  if(!is.null(core$gates$spike)){
    init[[nodes$indicator]] <- 1
    if(!is.null(nodes$inclusion)){
      init[[nodes$inclusion]] <- stats::rbeta(1, core$gates$spike$a, core$gates$spike$b)
    }
  }
  init[[sub("\\[.*$", "", nodes$eta[[1L]])]] <- stats::rgamma(core$n_units, shape = core$shapes$alpha)
  if(identical(core$convention, "linked")){
    init[[nodes$H]] <- stats::rgamma(1, shape = core$b)
  }else{
    init[[nodes$variable]] <- stats::rbeta(1, core$shapes$a, core$shapes$b)
  }

  init
}

# Joint prior draws of a core: gates in sequence, the budget indicator, the
# allocation coordinates with the shapes of the configuration, the slab total
# of the configuration, the realized total and the gated allocation.
.bt_R2_core_draws <- function(core, n){

  if(!is.null(core$pattern)){
    return(.bt_R2_core_pattern_draws(core, n))
  }
  nodes <- .bt_R2_core_nodes(core)
  gates <- core$gates
  out <- list()
  g <- matrix(0, n, gates$K)
  for(t in seq_len(gates$K)){
    s <- if(t == 1L) rep(0, n) else rowSums(g[, seq_len(t - 1L), drop = FALSE])
    g[, t] <- stats::rbinom(n, 1, .bt_R2_gate_probability(gates, t, s))
  }
  indicator <- rep(1, n)
  if(!is.null(gates$spike)){
    inclusion <- if(identical(gates$spike$type, "beta")) stats::rbeta(n, gates$spike$a, gates$spike$b) else rep(gates$spike$p, n)
    indicator <- stats::rbinom(n, 1, inclusion)
    if(!is.null(nodes$inclusion)) out[[nodes$inclusion]] <- inclusion
  }
  open_term <- vapply(core$terms, function(term){
    gate <- core$term_gate[[term]]
    indicator * (if(gate > 0L) g[, gate] else 1)
  }, numeric(n))
  open_term <- matrix(open_term, nrow = n, dimnames = list(NULL, core$terms))
  open <- open_term[, core$unit_term, drop = FALSE]
  open_alpha <- as.vector(open %*% core$unit_alpha)
  any_open <- as.numeric(rowSums(open) > 0)

  shape <- matrix(core$unit_alpha, n, core$n_units, byrow = TRUE)
  if(identical(core$concentration, "fixed_total")){
    rescaled <- shape * core$fixed_total / ifelse(open_alpha > 0, open_alpha, 1)
    shape <- ifelse(open > 0, rescaled, shape)
  }
  eta <- matrix(stats::rgamma(n * core$n_units, shape = as.vector(shape)), nrow = n)
  open_eta <- rowSums(open * eta)
  A_open <- if(identical(core$concentration, "fixed_total")) rep(core$fixed_total, n) else open_alpha
  if(identical(core$convention, "linked")){
    H <- stats::rgamma(n, shape = core$b)
    out[[nodes$H]] <- H
    variable <- open_eta / (open_eta + H)
  }else if(identical(core$convention, "conditioning")){
    A <- if(identical(core$concentration, "fixed_total")) core$fixed_total else core$A_reference
    variable <- stats::rbeta(n, ifelse(any_open > 0, core$a - (A - A_open), core$a),
                                ifelse(any_open > 0, core$b + (A - A_open), core$b))
  }else{
    variable <- stats::rbeta(n, core$a, core$b)
  }

  columns <- list()
  columns[[nodes$total]] <- variable * any_open
  if(!is.null(nodes$indicator)){
    columns[[nodes$indicator]] <- if(is.null(gates$spike)) any_open else indicator
  }
  columns <- c(columns, out)
  for(t in seq_len(gates$K)){
    columns[[nodes$gates[[t]]]] <- g[, t]
    columns[[nodes$open[[t]]]] <- g[, t] * indicator
  }
  if(!is.null(nodes$rate)){
    s <- rowSums(g)
    columns[[nodes$rate]] <- stats::rbeta(n, gates$law$a + s, gates$law$b + (gates$K - s))
  }
  if(core$n_units > 1L){
    allocation <- open * eta / (open_eta + (open_eta == 0))
    for(u in seq_len(core$n_units)){
      columns[[nodes$allocation[[u]]]] <- allocation[, u]
    }
  }
  for(u in seq_len(core$n_units)){
    columns[[nodes$eta[[u]]]] <- eta[, u]
  }
  if(!identical(core$convention, "linked")){
    columns[[nodes$variable]] <- variable
  }

  do.call(cbind, columns)
}

# Joint prior draws of a core conditioned on a pattern of open flags (a
# condition event, see .bt_R2_condition_core()): the draws of the pattern.
.bt_R2_core_pattern_draws <- function(core, n){

  pattern <- core$pattern
  free <- core
  free$pattern <- NULL
  nodes <- .bt_R2_core_nodes(free)
  out <- NULL
  for(attempt in seq_len(1000L)){
    draws <- .bt_R2_core_draws(free, max(n, 1000L))
    indicator <- if(!is.null(free$gates$spike)) draws[, nodes$indicator] else rep(1, nrow(draws))
    open <- vapply(free$terms, function(term){
      gate <- free$term_gate[[term]]
      indicator * (if(gate > 0L) draws[, nodes$gates[[gate]]] else 1)
    }, numeric(nrow(draws)))
    open <- matrix(open, nrow = nrow(draws), dimnames = list(NULL, free$terms))
    keep <- rep(TRUE, nrow(draws))
    for(term in names(pattern$flags)){
      keep <- keep & open[, term] == as.numeric(pattern$flags[[term]])
    }
    if(!is.na(pattern$total)){
      keep <- keep & (rowSums(open) > 0) == pattern$total
    }
    out <- rbind(out, draws[keep, , drop = FALSE])
    if(NROW(out) >= n){
      return(out[seq_len(n), , drop = FALSE])
    }
  }
  .bt_R2_unavailable("The condition event of the R-squared budget has too small a prior probability to sample its draws.",
    class = "BayesTools_R2_density_unavailable")
}

# The open flag of every unit of a term's budget from the gate and indicator
# draws (for the coefficient evaluator).
.bt_R2_unit_open_values <- function(gating, lookup, n){

  indicator <- if(!is.null(gating$indicator)) .bt_deterministic_lookup_value(lookup, gating$indicator) else rep(1, n)
  vapply(seq_along(gating$unit_gate), function(u){
    gate <- gating$unit_gate[[u]]
    indicator * (if(gate > 0L) .bt_deterministic_lookup_value(lookup, gating$gate_names[[gate]]) else 1)
  }, numeric(n))
}

.bt_R2_gated_unavailable <- function(what){
  .bt_R2_unavailable(paste0(
    what, " is unavailable for an R-squared budget with inclusion gates or a spike-and-slab total: ",
    "the single fit already averages over the models of the budget; use its inclusion summaries."),
    class = "BayesTools_R2_gated_unavailable")
}

# ---------------------------------------------------------------------------
# Prior laws of a gated or spiked budget. A configuration S of open terms
# enters every law of the budget only through its summed open concentration
# A_S (the Beta shapes of the total and, under fixed_total, the rescaling of
# the concentrations), whether a second term and a second unit are open, and,
# for a shared Beta rate, the number of open gates (its probability). The
# configurations are therefore aggregated into such states by an exact forward
# pass over the gates with the sequential gate law of the sampler, merging
# equal states; the number of gated terms is not limited, only the number of
# distinct states (at most K + 1 per open-concentration value with a shared
# rate, and K + 1 in all when the gated terms have equal concentrations).

.bt_R2_law_max_states <- 65536L

# The raw concentration of the term 'term' (a core term name): the sum of the
# concentrations of its units.
.bt_R2_core_term_alpha <- function(core, term){
  sum(core$unit_alpha[core$unit_term == term])
}

# The states of the open configurations of a core. 'fixed' (a named logical
# vector over gated core terms) restricts the configurations to the gates
# open (TRUE) or closed (FALSE); the result is then the joint law of the
# configuration and that pattern. 'target' (a core term open in every
# configuration of the result) is left out of 'A', so that 'A' is the open
# concentration of the other terms; without 'concentrations', 'A' is not
# tracked (probabilities of events only). Returns the states with a nonempty
# open set: their raw open concentration 'A', the number of open units and
# terms (both capped at 2), the number of open gates 's' (shared rate only)
# and the absolute probability 'prob' (with the spike weight and the
# truncation of a spiked budget whose terms are all gated); 'empty', the
# probability that the gates leave no term open while the spike is on (0
# under truncation), and 'off', the spike probability. The budget null has
# the probability off + empty. A core conditioned on a pattern of open flags
# ('pattern', see .bt_R2_condition_core()) gives the states of the pattern,
# normalized by its probability.
.bt_R2_core_states <- function(core, fixed = NULL, target = NULL, concentrations = TRUE){

  gates <- core$gates
  shared <- identical(gates$law$type, "shared")
  on <- if(is.null(gates$spike)) 1 else gates$spike$p
  ungated <- core$terms[core$term_gate == 0L]
  units_of <- function(term) sum(core$unit_term == term)
  alpha_of <- function(term) if(concentrations && !identical(term, target)) .bt_R2_core_term_alpha(core, term) else 0
  # the pattern: its gated flags fix the gates (the open flags given the
  # spike), its ungated flags need the spike open (TRUE) or closed (FALSE)
  pattern <- core$pattern
  allow <- c(on = TRUE, empty = TRUE, off = TRUE)
  if(!is.null(pattern)){
    flags <- pattern$flags
    gated_flags <- flags[core$term_gate[names(flags)] > 0L]
    both <- intersect(names(fixed), names(gated_flags))
    if(any(fixed[both] != gated_flags[both])){
      allow[] <- FALSE
    }
    fixed <- c(fixed, gated_flags[setdiff(names(gated_flags), names(fixed))])
    spike_open <- all(flags[core$term_gate[names(flags)] == 0L])
    allow["on"] <- allow[["on"]] && spike_open && !isFALSE(pattern$total)
    allow["empty"] <- allow[["empty"]] && spike_open && !isTRUE(pattern$total)
    allow["off"] <- allow[["off"]] && !any(flags) && !isTRUE(pattern$total)
  }
  states <- list(
    s     = 0,
    A     = sum(vapply(ungated, alpha_of, numeric(1))),
    units = min(2, sum(vapply(ungated, units_of, numeric(1)))),
    terms = min(2, length(ungated)),
    prob  = 1
  )
  for(t in seq_len(gates$K)){
    term <- gates$gated[[t]]
    p_open <- if(shared) (gates$law$a + states$s) / (gates$law$a + gates$law$b + (t - 1)) else gates$law$pi[[t]]
    force <- if(term %in% names(fixed)) fixed[[term]] else NA
    alpha <- alpha_of(term)
    open <- if(!isFALSE(force)) list(
      s     = states$s + 1,
      A     = states$A + alpha,
      units = pmin(2, states$units + units_of(term)),
      terms = pmin(2, states$terms + 1),
      prob  = states$prob * p_open
    )
    closed <- if(!isTRUE(force)) utils::modifyList(states, list(prob = states$prob *
      .bt_R2_gate_closed_probability(gates, t, states$s)))
    branches <- Filter(Negate(is.null), list(open, closed))
    states <- .bt_R2_states_merge(if(length(branches) == 2L) Map(c, open, closed) else branches[[1L]], shared, core)
  }
  empty <- states$terms == 0
  empty_mass <- sum(states$prob[empty])
  keep <- !empty & states$prob > 0 & allow[["on"]]
  prob <- states$prob[keep]
  if(gates$truncated){
    # the gate prior under the slab is conditioned on a nonempty configuration
    prob <- prob / .bt_R2_gate_nonempty_probability(gates)
    empty_mass <- 0
  }
  mass <- if(is.null(pattern)) 1 else pattern$mass

  list(
    A     = states$A[keep],
    units = states$units[keep],
    terms = states$terms[keep],
    s     = states$s[keep],
    prob  = on * prob / mass,
    empty = if(allow[["empty"]]) on * empty_mass / mass else 0,
    off   = if(allow[["off"]]) (1 - on) / mass else 0
  )
}

# A key that identifies a double exactly (17 significant digits round-trip):
# states merge only when their concentrations are equal, never when they are
# merely close (close raw concentrations can give very different shapes).
.bt_R2_exact_key <- function(x){
  formatC(x, digits = 17L, format = "g")
}

# Merges equal states (the number of open gates only under a shared rate,
# where it sets the probability of the later gates).
.bt_R2_states_merge <- function(states, shared, core){

  if(!shared){
    states$s <- rep(0, length(states$prob))
  }
  key <- paste(states$s, .bt_R2_exact_key(states$A), states$units, states$terms)
  first <- !duplicated(key)
  if(sum(first) > .bt_R2_law_max_states){
    .bt_R2_unavailable(paste0(
      "The prior laws of the R-squared budget '", core$name, "' need more than ", .bt_R2_law_max_states,
      " distinct states of its gate configurations (", sum(first), " at least; the gated terms have many ",
      "distinct concentrations). Use prior draws instead."),
      class = "BayesTools_R2_density_unavailable")
  }
  group <- match(key, key[first])
  list(
    s     = states$s[first],
    A     = states$A[first],
    units = states$units[first],
    terms = states$terms[first],
    prob  = as.vector(rowsum(states$prob, group, reorder = FALSE))
  )
}

# The Beta shapes of the total and the concentration rescaling of the states
# with the raw open concentrations 'A_raw': list(a, b, A, rescale), A the
# total concentration of the open units after the rescaling of fixed_total.
.bt_R2_core_state_shapes <- function(core, A_raw){

  fixed_total <- identical(core$concentration, "fixed_total")
  A <- if(fixed_total) rep(core$fixed_total, length(A_raw)) else A_raw
  A_ref <- if(fixed_total) core$fixed_total else core$A_reference
  a <- switch(core$convention, common = rep(core$a, length(A)), linked = A, conditioning = core$a - (A_ref - A))
  b <- switch(core$convention, common = , linked = rep(core$b, length(A)), conditioning = core$b + (A_ref - A))

  list(a = a, b = b, A = A, rescale = if(fixed_total) core$fixed_total / A_raw else rep(1, length(A_raw)))
}

# The probability that the gated budget term 'term' is closed (the budget
# null or its gate closed) and the open states of the term (with 'target').
.bt_R2_core_term_states <- function(core, term){

  gated <- core$term_gate[[term]] > 0L
  open <- .bt_R2_core_states(core, fixed = if(gated) stats::setNames(TRUE, term), target = term)
  closed <- if(gated){
    states <- .bt_R2_core_states(core, fixed = stats::setNames(FALSE, term))
    states$off + states$empty + sum(states$prob)
  }else{
    # an ungated term is closed only by the spike (it keeps a configuration
    # nonempty)
    open$off
  }

  list(open = open, closed = closed)
}

# A finite mixture of priors from components with weights (NULL components
# are dropped); a single component is returned as it is.
.bt_R2_mixture <- function(components, weights, is_null){

  keep <- weights > 0
  components <- components[keep]
  weights <- weights[keep] / sum(weights[keep])
  is_null <- is_null[keep]
  if(length(components) == 1L){
    return(components[[1L]])
  }
  components <- lapply(seq_along(components), function(i){
    component <- components[[i]]
    component$prior_weights <- weights[[i]]
    component
  })

  prior_mixture(components, is_null = is_null)
}

# Group Beta components with equal shapes.
.bt_R2_beta_components <- function(shapes, weights){

  key <- paste(.bt_R2_exact_key(shapes[, 1L]), .bt_R2_exact_key(shapes[, 2L]))
  groups <- split(seq_len(nrow(shapes)), factor(key, levels = unique(key)))
  list(
    components = lapply(groups, function(rows) prior("beta", list(alpha = shapes[rows[[1L]], 1L], beta = shapes[rows[[1L]], 2L]))),
    weights    = vapply(groups, function(rows) sum(weights[rows]), numeric(1))
  )
}

# The prior law of the realized total of a core: an atom at zero with the
# probability of the budget null and the mixture of the Beta totals of the
# open configurations; with 'conditional', the law given the budget is on.
.bt_R2_core_total_law <- function(core, conditional = FALSE){

  components <- .bt_R2_core_total_components(core)
  betas <- .bt_R2_beta_components(components$shapes, components$weights)
  if(conditional || components$atom <= 0){
    return(.bt_R2_mixture(betas$components, betas$weights, rep(FALSE, length(betas$weights))))
  }

  .bt_R2_mixture(c(list(prior("spike", list(0))), betas$components), c(components$atom, betas$weights),
                 c(TRUE, rep(FALSE, length(betas$weights))))
}

# The prior law of the share of the term 'term' (a core term name) given that
# it is open: the mixture over the open configurations that contain it of the
# Beta marginal of its concentration, an atom at one when it is the only open
# term.
.bt_R2_core_share_law <- function(core, term){

  states <- .bt_R2_core_term_states(core, term)$open
  alpha <- .bt_R2_core_term_alpha(core, term)
  # 'A' is the raw concentration of the other open terms (no cancellation)
  rescale <- .bt_R2_core_state_shapes(core, states$A + alpha)$rescale
  alone <- states$terms == 1
  only <- sum(states$prob[alone])
  shapes <- cbind(alpha * rescale, states$A * rescale)[!alone, , drop = FALSE]
  weights <- states$prob[!alone]
  betas <- if(nrow(shapes)) .bt_R2_beta_components(shapes, weights) else list(components = list(), weights = numeric())

  .bt_R2_mixture(c(if(only > 0) list(prior("spike", list(1))), betas$components),
                 c(if(only > 0) only, betas$weights),
                 c(if(only > 0) FALSE, rep(FALSE, length(betas$weights))))
}

# The Bernoulli law of an inclusion indicator with probability 'p'.
.bt_R2_bernoulli_law <- function(p){

  if(p >= 1) return(prior("spike", list(1)))
  if(p <= 0) return(prior("spike", list(0)))
  .bt_R2_mixture(list(prior("spike", list(0)), prior("spike", list(1))), c(1 - p, p), c(TRUE, FALSE))
}

# The gated core whose realized total is the coordinate 'coordinate' (NULL
# otherwise).
.bt_R2_gated_total_core <- function(prior_list, coordinate){

  if(length(coordinate) != 1L){
    return(NULL)
  }
  core <- prior_list[[coordinate]]
  if(is.prior.R2_budget(core)) core
}

# The prior law of the shared inclusion rate theta of a core: its Beta(a, b)
# hyperprior, and under the truncation of a spiked budget whose terms are all
# gated the law of theta given a nonempty configuration,
# Beta(theta; a, b) (1 - (1 - theta)^K) / (1 - c) with c = B(a, b + K) / B(a, b),
# which is the mixture over the first open gate j + 1 (j = 0, ..., K - 1) of
# Beta(a + 1, b + j) with the weights B(a + 1, b + j) / (B(a, b) (1 - c)).
.bt_R2_core_rate_law <- function(core){

  law <- core$gates$law
  if(!core$gates$truncated){
    return(prior("beta", list(alpha = law$a, beta = law$b)))
  }
  j <- seq_len(core$gates$K) - 1
  # B(a + 1, b + j) / B(a, b): the first j gates closed and gate j + 1 open,
  # as products of the sequential conditionals (no lbeta() differences)
  log_w <- vapply(j, function(first) .bt_R2_gate_log_configuration_probability(
    list(law = law, K = first + 1L), matrix(c(rep(FALSE, first), TRUE), nrow = 1L)), numeric(1))
  .bt_R2_mixture(lapply(j, function(first) prior("beta", list(alpha = law$a + 1, beta = law$b + first))),
                 exp(log_w - max(log_w)), rep(FALSE, length(j)))
}

# The gated core whose shared inclusion rate is the coordinate 'coordinate'.
.bt_R2_rate_core <- function(prior_list, coordinate){

  if(length(coordinate) != 1L){
    return(NULL)
  }
  for(core in prior_list[vapply(prior_list, is.prior.R2_budget, logical(1))]){
    if(identical(.bt_R2_core_nodes(core)$rate, coordinate)) return(core)
  }
  NULL
}

# The gated core whose budget indicator is the coordinate 'coordinate'.
.bt_R2_indicator_core <- function(prior_list, coordinate){

  if(length(coordinate) != 1L){
    return(NULL)
  }
  for(core in prior_list[vapply(prior_list, is.prior.R2_budget, logical(1))]){
    if(identical(.bt_R2_core_nodes(core)$indicator, coordinate)) return(core)
  }
  NULL
}

# ---------------------------------------------------------------------------
# Marginal prior methods of the terms of a gated or spiked budget: an atom at
# zero (the term closed) and the mixture over the open configurations that
# contain the term of the functional of model S.

.bt_R2_term_is_gated <- function(prior){
  !is.null(prior$core)
}

# The gated marginal of the row 'L' of a term: list(atom, specs, weights) with
# the absolute probabilities (atom + sum(weights) = 1).
.bt_R2_term_gated_specs <- function(prior, L){

  core <- prior$core
  found <- .bt_R2_term_states(prior)
  specs <- lapply(seq_along(found$weights), function(k){
    model <- prior
    model$core <- NULL
    model$shapes <- list(a = found$a[[k]], b = found$b[[k]], alpha = core$unit_alpha * found$rescale[[k]],
                         A = found$A[[k]])
    # two or more open units: the functional needs only whether the term
    # shares the budget
    model$n_units <- found$units[[k]]
    .bt_R2_term_functional(model, L)
  })

  structure(list(atom = found$atom, specs = specs, weights = found$weights), class = "BayesTools_R2_functional_mixture")
}

# The open states of a gated term grouped into the distinct models of its
# coefficients (Beta shapes, concentrations and whether it shares the budget),
# with their absolute probabilities, and the atom at zero (the term closed).
.bt_R2_term_states <- function(prior){

  core <- prior$core
  states <- .bt_R2_core_term_states(core, prior$term)
  open <- states$open
  shapes <- .bt_R2_core_state_shapes(core, open$A + .bt_R2_core_term_alpha(core, prior$term))
  key <- paste(.bt_R2_exact_key(shapes$a), .bt_R2_exact_key(shapes$b), .bt_R2_exact_key(shapes$A),
               .bt_R2_exact_key(shapes$rescale), open$units)
  first <- !duplicated(key)
  group <- match(key, key[first])

  list(
    atom    = states$closed,
    a       = shapes$a[first],
    b       = shapes$b[first],
    A       = shapes$A[first],
    rescale = shapes$rescale[first],
    units   = open$units[first],
    weights = as.vector(rowsum(open$prob, group, reorder = FALSE))
  )
}

.bt_R2_spec_cdf <- function(spec, q, lower.tail = TRUE){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    return(.bt_R2_functional_cdf(spec, q, lower.tail = lower.tail))
  }
  q <- as.numeric(q)
  continuous <- Reduce(`+`, lapply(seq_along(spec$specs), function(k){
    spec$weights[[k]] * .bt_R2_functional_cdf(spec$specs[[k]], q, lower.tail = lower.tail)
  }), rep(0, length(q)))
  atom <- if(lower.tail) spec$atom * (q >= 0) else spec$atom * (q < 0)

  continuous + atom
}

# The density of the continuous part (the atom at zero is a point mass that
# the cumulative distribution and the quantiles carry).
.bt_R2_spec_pdf <- function(spec, x){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    return(.bt_R2_functional_pdf(spec, x))
  }
  x <- as.numeric(x)
  Reduce(`+`, lapply(seq_along(spec$specs), function(k){
    spec$weights[[k]] * .bt_R2_functional_pdf(spec$specs[[k]], x)
  }), rep(0, length(x)))
}

.bt_R2_spec_log_pdf <- function(spec, x){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    return(.bt_R2_functional_log_pdf(spec, x))
  }
  # log-sum-exp of the weighted component log densities
  x <- as.numeric(x)
  terms <- vapply(seq_along(spec$specs), function(k){
    log(spec$weights[[k]]) + .bt_R2_functional_log_pdf(spec$specs[[k]], x)
  }, numeric(length(x)))
  terms <- matrix(terms, nrow = length(x))
  top <- apply(terms, 1, max)
  out <- top + log(rowSums(exp(terms - ifelse(is.finite(top), top, 0))))
  out[!is.finite(top)] <- top[!is.finite(top)]
  out
}

.bt_R2_spec_quantile <- function(spec, p){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    return(.bt_R2_functional_quantile(spec, p))
  }
  p <- as.numeric(p)
  # the continuous part is symmetric about zero (no evaluation near zero,
  # which a large scale would underflow)
  W <- sum(spec$weights)
  below <- W / 2
  at <- below + spec$atom
  vapply(p, function(prob){
    if(is.na(prob)) return(NA_real_)
    if(prob < 0 || prob > 1) .bt_R2_stop("The R-squared functional quantile needs probabilities in [0, 1].")
    if(prob == 0) return(-Inf)
    if(prob == 1) return(Inf)
    if(prob >= below && prob <= at) return(0)
    # the continuous part is a mixture of symmetric laws without atoms, so
    # the quantile lies between the component quantiles at the level of the
    # normalized continuous part; it is solved on log|x|, so that quantiles
    # near zero keep their relative precision
    side <- if(prob > at) 1 else -1
    level <- if(side > 0) (prob - spec$atom) / W else prob / W
    # (a level at the median of the continuous part is the atom itself)
    if((side > 0 && level <= 1 / 2) || (side < 0 && level >= 1 / 2)) return(0)
    if(level >= 1) return(Inf)
    if(level <= 0) return(-Inf)
    ends <- abs(vapply(spec$specs, .bt_R2_functional_quantile, numeric(1), p = level))
    finite <- ends[is.finite(ends) & ends > 0]
    unavailable <- function(){
      .bt_R2_unavailable("The quantile of the gated R-squared budget term is outside the representable range.",
        class = "BayesTools_R2_density_unavailable")
    }
    if(length(finite) == 0L) unavailable()
    f <- function(u) side * (.bt_R2_spec_cdf(spec, side * exp(u)) - prob)
    lo <- log(min(finite))
    hi <- log(max(finite))
    # a component quantile that underflows or overflows widens the bracket
    while(any(ends == 0) && f(lo) > 0 && lo > log(.Machine$double.xmin)) lo <- max(lo - 10, log(.Machine$double.xmin))
    while(any(!is.finite(ends)) && f(hi) < 0 && hi < log(.Machine$double.xmax)) hi <- min(hi + 10, log(.Machine$double.xmax))
    f_lo <- f(lo)
    f_hi <- f(hi)
    # the bracket is exact, so outside the range only where it was widened
    if((f_lo > 0 && any(ends == 0)) || (f_hi < 0 && any(!is.finite(ends)))) unavailable()
    if(hi - lo <= 0 || f_lo >= 0) return(side * exp(lo))
    if(f_hi <= 0) return(side * exp(hi))
    side * exp(stats::uniroot(f, c(lo, hi), f.lower = f_lo, f.upper = f_hi, tol = 1e-13)$root)
  }, numeric(1))
}

.bt_R2_spec_mean <- function(spec){

  specs <- if(inherits(spec, "BayesTools_R2_functional_mixture")) spec$specs else list(spec)
  if(all(vapply(specs, function(component) component$b > 1 / 2, logical(1)))) 0 else NaN
}

.bt_R2_spec_variance <- function(spec){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    variance <- .bt_R2_functional_variance(spec)
    return(if(is.finite(variance)) variance else NaN)
  }
  # zero mean: the second moment of the mixture
  variance <- sum(spec$weights * vapply(spec$specs, .bt_R2_functional_variance, numeric(1)))
  if(is.finite(variance)) variance else NaN
}

# Draws of the coefficients of a gated term from the joint prior of its
# budget.
.bt_R2_gated_term_rng <- function(prior, n){

  draw <- function(m){
    draws <- .bt_R2_core_draws(prior$core, m)
    directions <- matrix(.bt_R2_direction_rng(prior$law, m * prior$d), nrow = m)
    colnames(directions) <- .bt_R2_direction_coordinates(prior, "R2_term")
    draws <- cbind(draws, directions)
    node <- .bt_dnode_R2_coefficient("R2_term", prior)
    values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(draws, .bt_R2_budget_prior_list(prior)))
    list(values = matrix(values, nrow = m), open = .bt_R2_term_open_draws(prior, draws))
  }
  if(!isTRUE(prior$gated_slab)){
    return(draw(n)$values)
  }
  # the continuous part of the term, given that it is open
  out <- matrix(numeric(), 0L, prior$d)
  for(attempt in seq_len(1000L)){
    batch <- draw(max(n, 1000L))
    out <- rbind(out, batch$values[batch$open == 1, , drop = FALSE])
    if(nrow(out) >= n){
      return(out[seq_len(n), , drop = FALSE])
    }
  }
  .bt_R2_unavailable("The R-squared budget term is open with too small a prior probability to sample it given inclusion.",
    class = "BayesTools_R2_density_unavailable")
}

# The plot of a gated term: the continuous part and the atom at zero.
.bt_R2_gated_plot_components <- function(prior){

  spec <- .bt_R2_marginal_spec(prior)
  slab <- prior
  slab$gated_slab <- TRUE
  slab$prior_weights <- 1 - spec$atom
  list(slab, prior("spike", list(0), prior_weights = spec$atom))
}

# ---------------------------------------------------------------------------
# Mixed posteriors of gated budgets.

# The open flag of a gated term per draw from its gate and budget indicator.
.bt_R2_term_open_draws <- function(prior, samples){

  gate <- prior$core$term_gate[[prior$term]]
  indicator <- if(!is.null(prior$gating$indicator)) samples[, prior$gating$indicator] else rep(1, nrow(samples))
  as.numeric(indicator * (if(gate > 0L) samples[, prior$gating$gate_names[[gate]]] else 1))
}

# The atom at zero of the coefficients of a gated term (the term closed): the
# posterior mass of the closed draws and the mixed support with the point 0.
# 'law' is the prior law of the coefficient given the condition event (NULL:
# its unconditional law).
.bt_R2_set_gated_atoms <- function(samples, model_samples, parameter, prior, law = NULL){

  open <- .bt_R2_term_open_draws(prior, model_samples)
  columns <- if(is.null(dim(samples))) parameter else colnames(samples)
  closed <- mean(open == 0)
  # the prior density of a scalar coefficient given the conditioning event
  # (factor terms have no catalog density here; a law with too many states
  # leaves the draws without one)
  if(prior$d == 1L && !is.prior.factor(prior)){
    law <- tryCatch({
      if(is.null(law)){
        .bt_R2_gated_coefficient_law(stats::setNames(list(prior), parameter),
                                     list(type = "coordinate", dependencies = parameter))
      }else{
        # evaluate the law now so that an unavailable law leaves no density
        .bt_R2_check_law(law)
      }
    }, BayesTools_R2_density_unavailable = function(condition) NULL)
    samples <- .bt_R2_set_prior_density(samples, law)
  }
  # conditioned on its inclusion the term has no atom
  if(closed == 0){
    return(samples)
  }
  atoms <- .posterior_atoms_new(
    locations = matrix(0, nrow = 1L, ncol = length(columns), dimnames = list("location", columns)),
    mass = closed,
    column_names = columns,
    source = "posterior_indicator",
    component_probabilities = c(1 - closed, closed)
  )
  samples <- .posterior_atoms_set(samples, atoms)
  support <- .posterior_support_new(c(-Inf, Inf), points = 0, source = "prior_list", type = "mixed")
  samples <- if(is.null(dim(samples))){
    .posterior_support_set(samples, support)
  }else{
    .bt_meta_set(samples, "support", stats::setNames(rep(list(support), length(columns)), columns))
  }

  samples
}

# ---------------------------------------------------------------------------
# Structural supports, gate plans and prior laws of the quantities of gated
# budgets (no configuration enumeration for supports and plans).

# The gated R2 term (prior-list name) that owns the coordinates 'dependencies'
# (NULL when they are not all coordinates of one gated term).
.bt_R2_gated_term_owner <- function(prior_list, dependencies){

  if(length(dependencies) == 0L){
    return(NULL)
  }
  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.prior.R2(prior) && .bt_R2_term_is_gated(prior) &&
       all(dependencies %in% .bt_R2_coefficient_names(prior, parameter))){
      return(parameter)
    }
  }
  NULL
}

# The column whose draws are 1 where the gated term 'prior' is open (its open
# node, or the budget indicator of an ungated term of a spiked budget); NULL
# for a term that is always open.
.bt_R2_term_open_column <- function(prior){

  nodes <- .bt_R2_core_nodes(prior$core)
  gate <- prior$core$term_gate[[prior$term]]
  if(gate > 0L) nodes$open[[gate]] else if(!is.null(prior$core$gates$spike)) nodes$indicator
}

# Whether the term 'term' (a core term name) can be the only open term: every
# other term is gated.
.bt_R2_core_term_can_be_alone <- function(core, term){
  all(core$term_gate[setdiff(core$terms, term)] > 0L)
}

# The structural support of a quantity of a gated or spiked budget (NULL for
# other quantities).
.bt_R2_quantity_support <- function(object, key){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  continuous <- prior("beta", list(1, 1))
  if(identical(key$type, "R2_open")){
    return(.posterior_support_from_prior(.bt_R2_bernoulli_law(0.5), source = "catalog"))
  }
  if(identical(key$type, "R2_share")){
    record <- .bt_R2_budget_records(prior_list)[[key$total]]
    if(is.null(record$core)){
      return(NULL)
    }
    alone <- .bt_R2_core_term_can_be_alone(record$core, record$terms[[key$term]]$key)
    return(.posterior_support_from_prior(if(alone){
      .bt_R2_mixture(list(prior("spike", list(1)), continuous), c(0.5, 0.5), c(FALSE, FALSE))
    }else{
      continuous
    }, source = "catalog"))
  }
  if(!key$type %in% c("coordinate", "factor_level")){
    return(NULL)
  }
  if(identical(key$type, "coordinate") && length(key$dependencies) == 1L){
    core <- .bt_R2_gated_total_core(prior_list, key$dependencies)
    if(!is.null(core)){
      null <- !is.null(.bt_R2_core_nodes(core)$indicator) && .bt_R2_core_can_be_null(core)
      return(.posterior_support_from_prior(if(null){
        .bt_R2_mixture(list(prior("spike", list(0)), continuous), c(0.5, 0.5), c(TRUE, FALSE))
      }else{
        continuous
      }, source = "catalog"))
    }
    if(!is.null(.bt_R2_indicator_core(prior_list, key$dependencies))){
      return(.posterior_support_from_prior(.bt_R2_bernoulli_law(0.5), source = "catalog"))
    }
    if(!is.null(.bt_R2_rate_core(prior_list, key$dependencies))){
      return(.posterior_support_from_prior(continuous, source = "catalog"))
    }
  }
  owner <- .bt_R2_gated_term_owner(prior_list, key$dependencies)
  if(!is.null(owner) && !is.null(.bt_R2_term_open_column(prior_list[[owner]])) &&
     .bt_R2_term_can_close(prior_list[[owner]])){
    return(.posterior_support_new(c(-Inf, Inf), points = 0, source = "catalog", type = "mixed"))
  }
  if(!is.null(owner) && !is.null(.bt_R2_term_open_column(prior_list[[owner]]))){
    # a term that the condition pattern keeps open has no atom
    return(.posterior_support_new(c(-Inf, Inf), source = "catalog", type = "interval"))
  }

  NULL
}

# Whether the budget null (no term open) has a positive probability under the
# condition pattern of a core (always without a pattern); gate law only.
.bt_R2_core_can_be_null <- function(core){

  if(is.null(core$pattern)){
    return(TRUE)
  }
  states <- .bt_R2_core_states(core, concentrations = FALSE)
  states$off + states$empty > 0
}

# Whether a gated term can be closed under the condition pattern of its core
# (always without a pattern); gate law only.
.bt_R2_term_can_close <- function(prior){

  core <- prior$core
  if(is.null(core$pattern)){
    return(TRUE)
  }
  if(core$term_gate[[prior$term]] == 0L){
    return(.bt_R2_core_states(core, concentrations = FALSE)$off > 0)
  }
  states <- .bt_R2_core_states(core, fixed = stats::setNames(FALSE, prior$term), concentrations = FALSE)
  states$off + states$empty + sum(states$prob) > 0
}

# The prior law of a coefficient of a gated term: the atom at zero of the
# closed term and its continuous part given that it is open. A factor level
# (or a coded coefficient) is the linear combination of the term's
# coefficients of its catalog key.
.bt_R2_gated_coefficient_law <- function(prior_list, key, conditional = FALSE){

  owner <- .bt_R2_gated_term_owner(prior_list, key$dependencies)
  if(is.null(owner)){
    return(NULL)
  }
  slab <- .bt_R2_functional_law(prior_list[[owner]], owner, key)
  if(is.prior.point(slab)){
    return(slab)
  }
  slab$gated_slab <- TRUE
  atom <- .bt_R2_term_states(slab)$atom
  if(atom <= 0 || isTRUE(conditional)){
    return(slab)
  }

  .bt_R2_mixture(list(prior("spike", list(0)), slab), c(atom, 1 - atom), c(TRUE, FALSE))
}

# The bound R2 term (prior-list name) whose coefficients include all the
# coordinates 'dependencies' (NULL otherwise).
.bt_R2_term_owner <- function(prior_list, dependencies){

  if(length(dependencies) == 0L){
    return(NULL)
  }
  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.prior.R2(prior) && isTRUE(prior$bound) &&
       all(dependencies %in% .bt_R2_coefficient_names(prior, parameter))){
      return(parameter)
    }
  }
  NULL
}

# The scalar prior law, in residual-SD units (the fitted coordinates), of the
# linear combination of the coefficients of the R2 term 'parameter' that a
# catalog key declares (a coordinate, or a factor level with its weights):
# the term prior restricted to that combination ('functional_L'), whose
# methods give the combination's own law (an exact law for any combination; a
# kernel law for single coded coefficients). A combination that vanishes is a
# point at zero.
.bt_R2_functional_law <- function(prior, parameter, key){

  coefficients <- .bt_R2_coefficient_names(prior, parameter)
  weights <- if(is.null(key$weights)) rep(1, length(key$dependencies)) else as.numeric(key$weights)
  L <- numeric(prior$d)
  L[match(key$dependencies, coefficients)] <- weights
  if(all(L == 0)){
    return(prior("spike", list(0)))
  }
  law <- .bt_prior_without_multiply_by(prior)
  law$functional_L <- L
  attr(law, "factor_design") <- NULL
  class(law) <- c("prior", "prior.simple", "prior.R2")

  law
}

# The gate plan of a quantity of a gated or spiked budget: the point plan of
# a gated coefficient (atom at zero when its term is closed), the inclusion
# plan of an inclusion quantity, and the plan of a share (atom at one when its
# term is the only open term, undefined when it is closed). NULL otherwise.
.bt_R2_gate_plan <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  point_plan <- function(owner, indicator, on){
    plan_prior <- prior_spike_and_slab(prior("beta", list(1, 1)), prior_inclusion = prior("spike", list(on)))
    list(
      kind        = "point",
      origin      = "prior",
      owner       = owner,
      prior       = plan_prior,
      locations   = vapply(attr(plan_prior, "components"), function(component) if(component == "null") 0 else NA_real_, numeric(1)),
      prior_probabilities = .bt_prior_component_probabilities(plan_prior),
      indicator   = indicator,
      unknown     = FALSE,
      chain_gates = character(),
      event_components = attr(plan_prior, "components") == "alternative",
      event_reason     = NULL,
      event_gates = indicator,
      event_rule  = "AND"
    )
  }
  if(identical(key$type, "R2_open") ||
     (identical(key$type, "coordinate") && !is.null(.bt_R2_indicator_core(prior_list, key$dependencies)))){
    # the inclusion of a term, or of the budget (its indicator): the 0/1
    # draws are their own gate state
    return(list(kind = "inclusion", chain_gates = key$dependencies[[1L]], component_gates = character(),
                parent_gates = character(), event_gates = character(), event_rule = "AND"))
  }
  if(identical(key$type, "R2_share")){
    record <- .bt_R2_budget_records(prior_list)[[key$total]]
    if(is.null(record$core)){
      return(NULL)
    }
    return(list(kind = "R2_share", record = record, term = key$term, chain_gates = character(),
                event_gates = character(), event_rule = "AND"))
  }
  if(!key$type %in% c("coordinate", "factor_level")){
    return(NULL)
  }
  owner <- .bt_R2_gated_term_owner(prior_list, key$dependencies)
  if(!is.null(owner)){
    prior <- prior_list[[owner]]
    column <- .bt_R2_term_open_column(prior)
    if(is.null(column)){
      return(NULL)
    }
    marginals <- .bt_R2_gate_marginals(prior$core$gates, prior$core$terms)
    return(point_plan(owner, column, marginals$open[[prior$term]]))
  }

  NULL
}

# The per-draw states of an R2 share plan: an atom at one where its term is
# the only open term; undefined where it is closed.
.bt_R2_share_states <- function(plan, n, model_samples){

  open <- .bt_R2_record_open(plan$record, model_samples)
  if(is.null(open)){
    return(list(atom = rep(NA_real_, n), defined = rep(TRUE, n), event = NULL, known = FALSE))
  }
  alone <- open[, plan$term] == 1 & rowSums(open) == 1
  list(atom = ifelse(alone, 1, NA_real_), defined = open[, plan$term] == 1, event = NULL, known = TRUE)
}

# ---------------------------------------------------------------------------
# Prior methods of the realized total of a gated or spiked budget (the core
# entry): an atom at zero and a mixture of Beta totals.

.bt_R2_core_total_components <- function(core){

  states <- .bt_R2_core_states(core)
  shapes <- .bt_R2_core_state_shapes(core, states$A)
  betas <- .bt_R2_beta_components(cbind(shapes$a, shapes$b), states$prob)
  list(atom = states$off + states$empty,
       shapes = do.call(rbind, lapply(betas$components, function(component){
         c(component$parameters[["alpha"]], component$parameters[["beta"]])
       })),
       weights = unname(betas$weights))
}

# The root x of sum(w * pbeta(x, a, b)) = mass, at or below the median of the
# normalized mixture: bracketed by the component quantiles at the level of the
# mass and solved on the logit scale (relative precision near zero), with
# the bracket widened through the subnormal range when a component quantile
# underflows; unavailable below the double range.
.bt_R2_beta_mixture_lower_root <- function(a, b, w, mass){

  # (the inverse logit as exp(u) / (1 + exp(u)) below zero keeps subnormal
  # values that 1 / (1 + exp(-u)) underflows)
  inverse <- function(u) if(u < 0) exp(u) / (1 + exp(u)) else 1 - exp(-u) / (1 + exp(-u))
  f <- function(u) sum(w * stats::pbeta(inverse(u), a, b)) - mass
  ends <- stats::qlogis(range(stats::qbeta(min(max(mass / sum(w), 0), 1), a, b)))
  lowest <- log(.Machine$double.xmin) - 52 * log(2)
  if(!is.finite(ends[[1L]])){
    ends[[1L]] <- stats::qlogis(.Machine$double.xmin)
    while(f(ends[[1L]]) > 0 && ends[[1L]] > lowest) ends[[1L]] <- max(ends[[1L]] - 4, lowest)
    if(f(ends[[1L]]) > 0){
      .bt_R2_unavailable("The quantile of the gated R-squared budget total is beyond the double range.",
        class = "BayesTools_R2_density_unavailable")
    }
  }
  if(!is.finite(ends[[2L]])){
    ends[[2L]] <- 0
  }
  f_ends <- c(f(ends[[1L]]), f(ends[[2L]]))
  if(f_ends[[1L]] >= 0) return(inverse(ends[[1L]]))
  if(f_ends[[2L]] <= 0) return(inverse(ends[[2L]]))
  inverse(stats::uniroot(f, ends, f.lower = f_ends[[1L]], f.upper = f_ends[[2L]], tol = 1e-13)$root)
}

.bt_R2_core_method <- function(prior, method, value = NULL, ...){

  if(identical(method, "rng")){
    nodes <- .bt_R2_core_nodes(prior)
    return(unname(.bt_R2_core_draws(prior, value)[, nodes$total]))
  }
  if(method %in% c("cdf", "ccdf", "pdf", "lpdf", "quant")){
    method <- paste0("m", method)
  }
  components <- .bt_R2_core_total_components(prior)
  a <- components$shapes[, 1L]
  b <- components$shapes[, 2L]
  w <- components$weights
  cdf <- function(q, lower.tail = TRUE){
    q <- as.numeric(q)
    continuous <- vapply(q, function(v) sum(w * stats::pbeta(v, a, b, lower.tail = lower.tail)), numeric(1))
    continuous + if(lower.tail) components$atom * (q >= 0) else components$atom * (q < 0)
  }
  switch(
    method,
    mcdf  = cdf(value),
    mccdf = cdf(value, lower.tail = FALSE),
    mpdf  = vapply(as.numeric(value), function(v) sum(w * stats::dbeta(v, a, b)), numeric(1)),
    mlpdf = vapply(as.numeric(value), function(v){
      terms <- log(w) + stats::dbeta(v, a, b, log = TRUE)
      top <- max(terms)
      if(is.finite(top)) top + log(sum(exp(terms - top))) else top
    }, numeric(1)),
    mquant = vapply(as.numeric(value), function(p){
      if(is.na(p)) return(NA_real_)
      if(p < 0 || p > 1) .bt_R2_stop("Quantiles need probabilities in [0, 1].")
      # (the upper end first: a tiny continuous part can round the atom to 1)
      if(p >= 1) return(if(length(w) > 0L && any(w > 0)) 1 else 0)
      if(p <= components$atom) return(0)
      # the quantile of the Beta mixture lies between the Beta quantiles at the
      # level of its normalized continuous part (exact for one component);
      # it is solved on the logit scale, which keeps the relative precision
      # of quantiles near 0 and 1
      level <- min(max((p - components$atom) / sum(w), 0), 1)
      ends <- stats::qbeta(level, a, b)
      if(length(w) == 1L || max(ends) <= min(ends)){
        return(ends[[1L]])
      }
      # above the median of the continuous part the complement 1 - r is
      # solved (its shapes swapped, its mass 1 - p), so that a quantile next
      # to one rounds to its nearest double
      if(level > 1 / 2){
        return(1 - .bt_R2_beta_mixture_lower_root(b, a, w, 1 - p))
      }
      .bt_R2_beta_mixture_lower_root(a, b, w, p - components$atom)
    }, numeric(1)),
    mean = sum(w * a / (a + b)),
    var = {
      # the atom at zero and the Beta components around the mixture mean (no
      # difference of raw second moments, which cancels for concentrated
      # totals)
      # (weights normalized with the atom; the spread between the components
      # and the atom at zero as pairwise differences of their means, so that
      # one component has none whatever its rounded weight; the variance of a
      # component as ratios, as the square of the shape sum can overflow)
      total <- components$atom + sum(w)
      p <- c(components$atom, w) / total
      # (the shape sum is a finite double, .bt_R2_check_beta_sum())
      mu <- a / (a + b)
      within <- mu * (b / (a + b)) / (a + b + 1)
      means <- c(0, mu)
      spread <- outer(means, means, "-")^2
      sum(p[-1L] * within) + sum(outer(p, p) * spread) / 2
    },
    .bt_R2_unavailable(paste0("The '", method, "' method is unavailable for the total of a gated R-squared budget."))
  )
}

# The mixed posterior of the realized total of a gated or spiked budget from
# the (conditioned) draws: its prior law, the mixed support with the atom at
# zero, and the posterior mass of the atom.
.bt_R2_as_mixed_total <- function(model_samples, core, parameter, law = NULL){

  law <- tryCatch(if(is.null(law)) .bt_R2_core_total_law(core) else law,
                  BayesTools_R2_density_unavailable = function(condition) NULL)
  samples <- .as_mixed_posteriors.simple(model_samples, if(is.null(law)) prior("beta", list(1, 1)) else law, parameter)
  attr(samples, "prior_list") <- if(is.null(law)) core else law
  if(!is.null(.bt_R2_core_nodes(core)$indicator) && any(samples == 0)){
    closed <- mean(samples == 0)
    samples <- .posterior_atoms_set(samples, .posterior_atoms_new(
      locations = matrix(0, nrow = 1L, ncol = 1L, dimnames = list("location", parameter)),
      mass = closed, column_names = parameter, source = "posterior_indicator",
      component_probabilities = c(1 - closed, closed)
    ))
    samples <- .posterior_support_set(samples,
      .posterior_support_new(c(0, 1), points = 0, source = "prior_list", type = "mixed"))
  }

  .bt_R2_set_prior_density(samples, law)
}

# The prior mean of a term's coefficients: zero when the second Beta shape of
# the total exceeds 1/2 in every model where the term is open (each
# coefficient is symmetric about zero), NaN otherwise. No marginal density is
# needed, nor the configurations: the second shape is the same in every
# model (common, linked, fixed_total) or decreases with the open
# concentration (conditioning), so its smallest value is that of the model
# with every retained term open that the condition pattern of the core allows
# (it has a positive probability).
.bt_R2_term_mean <- function(prior){

  if(!.bt_R2_term_is_gated(prior)){
    return(if(prior$shapes$b > 1 / 2) 0 else NaN)
  }
  core <- prior$core
  flags <- core$pattern$flags
  open <- if(is.null(flags)) rep(TRUE, core$n_units) else !core$unit_term %in% names(flags)[!flags]
  if(!any(open)){
    return(0)
  }
  b <- .bt_R2_core_state_shapes(core, sum(core$unit_alpha[open]))$b
  if(b > 1 / 2) 0 else NaN
}

# ---------------------------------------------------------------------------
# Condition events of gated budgets on the prior side.

# The gated budget (its total name) of a prior-list entry: a gated term or the
# total entry; NULL otherwise.
.bt_R2_condition_budget <- function(prior){

  if(is.prior.R2(prior) && .bt_R2_term_is_gated(prior)){
    return(prior$core$total_name)
  }
  if(is.prior.R2_budget(prior)){
    return(prior$total_name)
  }
  NULL
}

# The labels of a condition event that belong to one gated budget are
# dependent (shared rate, truncation, spike): they form one joint family whose
# options are the truth patterns of their open flags (a term open; the budget
# on for its total), each with its exact joint prior probability, and whose
# prior lists replace every entry of the budget (labelled or not) by its prior
# with the core conditioned on the pattern (a term closed by the pattern is a
# point at zero): the laws, the draws and the dependence of the entries all
# follow the pattern. Families of different budgets and other priors stay
# independent.

# The prior-list names of the entries of the gated budget 'total_name'.
.bt_R2_budget_members <- function(prior_list, total_name){

  names(prior_list)[vapply(prior_list, function(prior){
    identical(.bt_R2_condition_budget(prior), total_name)
  }, logical(1))]
}

# The probability of a pattern of a gated budget: 'flags', the open flags of
# labelled terms (named by core term names), and 'total', the budget-on flag
# (NA when it is not labelled). With the spike, every term is closed and the
# budget off; otherwise the ungated terms are open and the gated terms follow
# their gates. It needs the gate law only (no concentrations).
.bt_R2_pattern_mass <- function(core, flags, total = NA){

  core$pattern <- NULL
  on <- if(is.null(core$gates$spike)) 1 else core$gates$spike$p
  off <- if(!any(flags) && !isTRUE(total)) 1 - on else 0
  gated <- core$term_gate[names(flags)] > 0L
  if(!all(flags[!gated])){
    return(off)
  }
  states <- .bt_R2_core_states(core, fixed = flags[gated], concentrations = FALSE)

  off + (if(!isFALSE(total)) sum(states$prob) else 0) + (if(!isTRUE(total)) states$empty else 0)
}

# The core 'core' conditioned on the pattern ('flags', 'total').
.bt_R2_condition_core <- function(core, flags, total = NA){

  core$pattern <- NULL
  core$pattern <- list(flags = flags, total = total, mass = .bt_R2_pattern_mass(core, flags, total))
  core
}

# The priors of the entries 'members' of a gated budget given the pattern
# ('flags', 'total').
.bt_R2_pattern_laws <- function(prior_list, members, core, flags, total){

  conditioned <- .bt_R2_condition_core(core, flags, total)
  laws <- lapply(members, function(member){
    prior <- prior_list[[member]]
    if(is.prior.R2_budget(prior)){
      on <- if(isFALSE(total)) 0 else .bt_R2_pattern_mass(core, flags, TRUE)
      return(if(on > 0) .bt_R2_condition_core(prior, flags, total) else prior("spike", list(0)))
    }
    term <- prior$term
    open <- if(term %in% names(flags)) as.numeric(flags[[term]]) else
      .bt_R2_pattern_mass(core, c(flags, stats::setNames(TRUE, term)), total)
    if(open <= 0){
      # a closed factor keeps its coded coefficients (a factor point)
      return(if(is.prior.factor(prior)) prior_factor("spike", list(0), contrast = prior$budget$contrast) else
        prior("spike", list(0)))
    }
    prior$core <- conditioned
    prior
  })
  names(laws) <- members

  laws
}

# The joint condition families of the gated budgets of 'families' (the
# families of a condition event): list(members, options), with 'members' the
# family names they replace and 'options' one list of options per budget in
# the format of .condition_event_family_options(), with 'replacements' (the
# priors of every entry of the budget given the pattern).
.bt_R2_condition_joint_options <- function(prior_list, families){

  budgets <- vapply(names(families), function(name){
    budget <- .bt_R2_condition_budget(prior_list[[name]])
    if(is.null(budget)) NA_character_ else budget
  }, character(1))
  out <- list(members = character(), options = list())
  for(budget in unique(budgets[!is.na(budgets)])){
    labelled <- names(families)[budgets %in% budget]
    core <- prior_list[[budget]]
    if(!is.prior.R2_budget(core)){
      core <- prior_list[[labelled[[1L]]]]$core
    }
    core$pattern <- NULL
    members <- .bt_R2_budget_members(prior_list, budget)
    term_of <- vapply(labelled, function(name){
      prior <- prior_list[[name]]
      if(is.prior.R2_budget(prior)) NA_character_ else prior$term
    }, character(1))
    grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(labelled))))
    options <- list()
    for(i in seq_len(nrow(grid))){
      pattern <- grid[i, ]
      flags <- stats::setNames(pattern[!is.na(term_of)], term_of[!is.na(term_of)])
      total <- if(any(is.na(term_of))) pattern[is.na(term_of)][[1L]] else NA
      mass <- .bt_R2_pattern_mass(core, flags, total)
      if(!(mass > 0)){
        next
      }
      values <- unlist(lapply(seq_along(labelled), function(j){
        labels <- unique(families[[labelled[[j]]]][["labels"]])
        stats::setNames(rep(pattern[[j]], length(labels)), labels)
      }))
      options[[length(options) + 1L]] <- list(
        family       = budget,
        prior        = NULL,
        probability  = mass,
        log_probability = log(mass),
        values       = values,
        replace      = TRUE,
        replacements = .bt_R2_pattern_laws(prior_list, members, core, flags, total)
      )
    }
    out$members <- c(out$members, labelled)
    out$options[[budget]] <- options
  }

  out
}

# The prior density of a mixed posterior of a gated quantity, on the
# conditioning event when its own inclusion is conditioned.
.bt_R2_set_prior_density <- function(samples, law){

  if(is.null(law)){
    return(samples)
  }
  density <- .bt_parameter_prior_density_scalar(.bt_R2_grid_law(law), n_grid = .prior_linear_density_default_grid(),
                                                tail_prob = .prior_linear_density_tail_prob())
  attr(samples, "prior_list") <- law
  if(!is.null(density)){
    samples <- .bt_meta_set(samples, "prior_density", density)
  }
  samples
}

# The prior law of the gated quantity 'parameter' (a gated term or a total)
# given the condition event of as_mixed_posteriors(): NULL when the event does
# not involve its budget (it is independent of the other conditioned priors),
# and otherwise the mixture over the options of the event of its law in each
# option, decomposed into the atom at zero (merged) and the continuous parts
# (a slab of the term, the Beta components of the total), on the residual-SD
# scale of the catalog laws.
.bt_R2_conditioned_law <- function(prior_list, condition_event, parameter){

  budget <- .bt_R2_condition_budget(prior_list[[parameter]])
  families <- condition_event[["families"]]
  involved <- vapply(names(families), function(name){
    identical(.bt_R2_condition_budget(prior_list[[name]]), budget)
  }, logical(1))
  if(is.null(budget) || !any(involved)){
    return(NULL)
  }
  options <- .condition_event_model_options(prior_list, condition_event)
  components <- list()
  weights <- numeric()
  atom <- 0
  for(i in seq_along(options$prior_lists)){
    law <- options$prior_lists[[i]][[parameter]]
    weight <- options$weights[[i]]
    if(is.prior.point(law)){
      atom <- atom + weight
    }else if(is.prior.R2_budget(law)){
      total <- .bt_R2_core_total_components(law)
      atom <- atom + weight * total$atom
      for(k in seq_along(total$weights)){
        components[[length(components) + 1L]] <- prior("beta", list(alpha = total$shapes[k, 1L], beta = total$shapes[k, 2L]))
        weights <- c(weights, weight * total$weights[[k]])
      }
    }else{
      closed <- .bt_R2_term_states(law)$atom
      atom <- atom + weight * closed
      if(closed < 1){
        slab <- .bt_prior_without_multiply_by(law)
        slab$gated_slab <- TRUE
        slab$prior_weights <- 1
        components[[length(components) + 1L]] <- slab
        weights <- c(weights, weight * (1 - closed))
      }
    }
  }

  .bt_R2_mixture(c(list(prior("spike", list(0))), components), c(atom, weights),
                 c(TRUE, rep(FALSE, length(components))))
}

# A conditioned law whose slab components need the states of their budget:
# the states are computed once so that a law beyond the state limit is
# refused here (BayesTools_R2_density_unavailable).
.bt_R2_check_law <- function(law){

  parts <- if(is.prior.mixture(law)) unclass(law) else list(law)
  for(part in parts){
    if(is.prior.R2(part) && .bt_R2_term_is_gated(part)){
      .bt_R2_term_states(part)
    }
  }
  law
}

# A law whose R2 coefficient components evaluate the densities of prior-density
# grids by interpolation (.bt_R2_spec_log_pdf_grid()) and cache their
# quantiles: every exact point is a deterministic quadrature over the budget
# (nested for factor levels), while prior-density grids, their refinements
# for region probabilities and their ranges need thousands of points. The
# cache is an environment shared by the copies of the law.
.bt_R2_grid_law <- function(law){

  if(is.prior.R2(law)){
    law$grid_interpolation <- new.env(parent = emptyenv())
    return(law)
  }
  if(is.prior.mixture(law)){
    for(k in seq_along(law)){
      if(is.prior.R2(law[[k]])) law[[k]]$grid_interpolation <- new.env(parent = emptyenv())
    }
  }
  law
}

# The log density of the continuous part of a (symmetric) R2 coefficient law
# 'spec' for prior-density grids: a cubic spline in log |x| through the exact
# log density at 256 nodes uniform in log |x| (the log densities of these laws
# are smooth in log |x|, power laws near zero and in the tails), cached in
# 'cache' and rebuilt when the values leave its domain. The domain extends
# the values by a factor of 10 on both sides, within 1e-8 to 1e8 times the
# scale of the law; zero, values outside the domain and the first short
# vectors are exact.
.bt_R2_spec_log_pdf_grid <- function(spec, x, cache = new.env(parent = emptyenv())){

  x <- as.numeric(x)
  y <- abs(x)
  out <- rep(NA_real_, length(x))
  scale <- if(inherits(spec, "BayesTools_R2_functional_mixture")){
    max(vapply(spec$specs, `[[`, numeric(1), "scale"))
  }else{
    spec$scale
  }
  inner <- !is.na(y) & is.finite(y) & y > scale * 1e-8 & y < scale * 1e8
  if(is.null(cache$spline) && length(unique(y[inner])) < 32L){
    inner[] <- FALSE
  }
  if(any(inner)){
    limits <- range(y[inner])
    if(is.null(cache$spline) || limits[[1L]] < cache$lower || limits[[2L]] > cache$upper){
      lower <- max(min(limits[[1L]], cache$lower) / 10, scale * 1e-8)
      upper <- min(max(limits[[2L]], cache$upper) * 10, scale * 1e8)
      cache$spline <- .bt_R2_log_pdf_spline(spec, lower, upper)
      cache$lower  <- lower
      cache$upper  <- upper
    }
    if(is.null(cache$spline)){
      inner[] <- FALSE
    }
  }
  if(any(inner)){
    out[inner] <- cache$spline(log(y[inner]))
  }
  exact <- !inner & !is.na(x)
  if(any(exact)){
    out[exact] <- .bt_R2_spec_log_pdf(spec, x[exact])
  }
  out
}

# The spline of the log density of 'spec' in log |x| on [lower, upper]: 256
# nodes uniform in log |x|, validated at every midpoint between nodes (the
# exact log density there within 1e-6) and refined by adding the midpoints as
# nodes (up to three times); NULL when it does not validate, and the
# densities are then evaluated exactly.
.bt_R2_log_pdf_spline <- function(spec, lower, upper){

  log_nodes <- seq(log(lower), log(upper), length.out = 256L)
  log_f <- .bt_R2_spec_log_pdf(spec, exp(log_nodes))
  for(round in seq_len(4L)){
    finite <- is.finite(log_f)
    if(!all(finite)){
      return(NULL)
    }
    spline <- stats::splinefun(log_nodes, log_f, method = "fmm")
    log_mid <- (log_nodes[-1L] + log_nodes[-length(log_nodes)]) / 2
    log_f_mid <- .bt_R2_spec_log_pdf(spec, exp(log_mid))
    if(all(is.finite(log_f_mid)) && max(abs(spline(log_mid) - log_f_mid)) <= 1e-6){
      return(spline)
    }
    if(round == 4L){
      return(NULL)
    }
    order <- order(c(log_nodes, log_mid))
    log_nodes <- c(log_nodes, log_mid)[order]
    log_f <- c(log_f, log_f_mid)[order]
  }
  NULL
}

# The quantiles of an R2 coefficient law, cached in 'cache' by probability.
.bt_R2_spec_quantile_cached <- function(spec, p, cache){

  if(is.null(cache$quantiles)){
    cache$quantiles <- numeric()
  }
  keys <- formatC(p, digits = 17L, format = "g")
  missing <- !keys %in% names(cache$quantiles)
  if(any(missing)){
    cache$quantiles[keys[missing]] <- .bt_R2_spec_quantile(spec, p[missing])
  }
  unname(cache$quantiles[keys])
}
