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

# P(gate t open | gates 1..t-1 closed, at least one gate open), t = 1..K.
.bt_R2_gate_truncated_probabilities <- function(law, K){

  p <- vapply(seq_len(K), function(t){
    if(identical(law$type, "shared")){
      # P(g_t = 1 | none before) / P(some of t..K open | none before) under BB
      first <- law$a / (law$a + law$b + t - 1)
      none <- exp(lbeta(law$a, law$b + K) - lbeta(law$a, law$b + t - 1))
      first / (1 - none)
    }else{
      law$pi[[t]] / (1 - prod(1 - law$pi[t:K]))
    }
  }, numeric(1))

  pmin(1, p)
}

# The JAGS expression of P(gate t open | s open gates before it), with 's' the
# JAGS expression of that count ("0" for t = 1).
.bt_R2_gate_probability_syntax <- function(gates, t, s){

  regular <- if(identical(gates$law$type, "shared")){
    paste0("(", .JAGS_prior_number(gates$law$a), " + ", s, ") / ",
           .JAGS_prior_number(gates$law$a + gates$law$b + t - 1))
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
    (gates$law$a + s) / (gates$law$a + gates$law$b + t - 1)
  }else{
    rep(gates$law$pi[[t]], length(s))
  }
  if(gates$truncated){
    regular[s == 0] <- gates$p_trunc[[t]]
  }

  regular
}

# Prior probability of the gate configurations 'open' (a logical matrix, one
# row per configuration of the K gates) under the gate law, before the spike:
# the probability that exactly these gates are open.
.bt_R2_gate_configuration_probability <- function(gates, open){

  open <- matrix(as.logical(open), ncol = gates$K)
  s <- rowSums(open)
  base <- if(identical(gates$law$type, "shared")){
    exp(lbeta(gates$law$a + s, gates$law$b + gates$K - s) - lbeta(gates$law$a, gates$law$b))
  }else{
    apply(open, 1, function(row) prod(ifelse(row, gates$law$pi, 1 - gates$law$pi)))
  }
  if(gates$truncated){
    empty <- .bt_R2_gate_configuration_probability(
      utils::modifyList(gates, list(truncated = FALSE)), matrix(FALSE, 1, gates$K))
    base <- ifelse(s == 0, 0, base / (1 - empty))
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
      empty <- .bt_R2_gate_configuration_probability(
        utils::modifyList(gates, list(truncated = FALSE)), matrix(FALSE, 1, gates$K))
      pi <- pi / (1 - empty)
    }
    open[gates$gated] <- on_spike * pi
  }
  budget <- if(gates$K > 0L && gates$all_gated && is.null(gates$spike)){
    1 - .bt_R2_gate_configuration_probability(gates, matrix(FALSE, 1, gates$K))
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
    gates      = if(core$gates$K > 0L) paste0(T, "_gate[", seq_len(core$gates$K), "]"),
    allocation = if(core$n_units > 1L) paste0(core$allocation_name, "[", seq_len(core$n_units), "]"),
    eta        = paste0(.JAGS_prior_dirichlet_eta_name(if(core$n_units > 1L) core$allocation_name else paste0(T, "_allocation")),
                        "[", seq_len(core$n_units), "]"),
    H          = if(identical(core$convention, "linked")) paste0("prior_par_H_", T)
  )
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

  # the gates, each given the number of open gates before it
  for(t in seq_len(gates$K)){
    s <- if(t == 1L) "0" else paste0("sum(", core$total_name, "_gate[1:", t - 1L, "])")
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
    lines <- c(lines, paste0(nodes$indicator, " <- step(sum(", core$total_name, "_gate[1:", gates$K, "]) - 0.5)"))
  }

  open <- .bt_R2_core_unit_open_syntax(core, nodes)
  open_alpha <- paste0("(", paste0(open, " * ", vapply(core$unit_alpha, num, character(1)), collapse = " + "), ")")
  any_open <- paste0("step(", paste0(open, collapse = " + "), " - 0.5)")

  # allocation gamma coordinates (shapes rescaled per configuration under
  # fixed_total; a closed unit keeps its declared shape)
  for(u in seq_len(core$n_units)){
    shape <- if(identical(core$concentration, "fixed_total") && core$n_units > 1L){
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
  c(nodes$total, nodes$indicator, nodes$inclusion,
    if(core$gates$K > 0L) paste0(core$total_name, "_gate"),
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
  if(identical(core$concentration, "fixed_total") && core$n_units > 1L){
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
# Prior laws of a gated or spiked budget, by enumeration of the gate
# configurations (the sampler is not limited; these analytic laws are).

.bt_R2_law_max_gates <- 16L

# The open configurations of a core with their prior probabilities, and the
# probability of the budget null.
.bt_R2_core_configurations <- function(core){

  K <- core$gates$K
  if(K > .bt_R2_law_max_gates){
    .bt_R2_unavailable(paste0(
      "The prior laws of the R-squared budget '", core$name, "' enumerate its gate configurations, ",
      "which is unavailable for more than ", .bt_R2_law_max_gates, " gated terms."),
      class = "BayesTools_R2_density_unavailable")
  }
  grid <- if(K > 0L) as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K))) else matrix(FALSE, 1L, 0L)
  p_gate <- if(K > 0L) .bt_R2_gate_configuration_probability(core$gates, grid) else 1
  on <- if(is.null(core$gates$spike)) 1 else core$gates$spike$p
  ungated <- core$terms[core$term_gate == 0L]
  configurations <- list()
  null <- 1 - on
  for(i in seq_len(nrow(grid))){
    open <- c(ungated, core$gates$gated[grid[i, ]])
    if(length(open) == 0L){
      null <- null + on * p_gate[[i]]
      next
    }
    if(on * p_gate[[i]] > 0){
      configurations[[length(configurations) + 1L]] <- list(open = open, prob = on * p_gate[[i]])
    }
  }

  list(configurations = configurations, null = null)
}

# The Beta shapes of the total and the concentrations of the open terms in the
# configuration with the open terms 'open'.
.bt_R2_core_configuration_shapes <- function(core, open){

  units <- core$unit_term %in% open
  unit_alpha <- core$unit_alpha
  if(identical(core$concentration, "fixed_total")){
    unit_alpha[units] <- unit_alpha[units] * core$fixed_total / sum(unit_alpha[units])
  }
  A <- sum(unit_alpha[units])
  shapes <- switch(
    core$convention,
    common = c(core$a, core$b),
    linked = c(A, core$b),
    conditioning = {
      A_ref <- if(identical(core$concentration, "fixed_total")) core$fixed_total else core$A_reference
      c(core$a - (A_ref - A), core$b + (A_ref - A))
    }
  )

  list(a = shapes[[1L]], b = shapes[[2L]], A = A,
       term_alpha = vapply(open, function(term) sum(unit_alpha[core$unit_term == term]), numeric(1)))
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

  key <- paste(signif(shapes[, 1L], 14), signif(shapes[, 2L], 14))
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

  configurations <- .bt_R2_core_configurations(core)
  shapes <- do.call(rbind, lapply(configurations$configurations, function(configuration){
    s <- .bt_R2_core_configuration_shapes(core, configuration$open)
    c(s$a, s$b)
  }))
  weights <- vapply(configurations$configurations, `[[`, numeric(1), "prob")
  betas <- .bt_R2_beta_components(shapes, weights)
  if(conditional || configurations$null <= 0){
    return(.bt_R2_mixture(betas$components, betas$weights, rep(FALSE, length(betas$weights))))
  }

  .bt_R2_mixture(c(list(prior("spike", list(0))), betas$components), c(configurations$null, betas$weights),
                 c(TRUE, rep(FALSE, length(betas$weights))))
}

# The prior law of the share of the term 'term' (a core term name) given that
# it is open: the mixture over the open configurations that contain it of the
# Beta marginal of its concentration, an atom at one when it is the only open
# term.
.bt_R2_core_share_law <- function(core, term){

  configurations <- Filter(function(configuration) term %in% configuration$open,
                           .bt_R2_core_configurations(core)$configurations)
  components <- list()
  weights <- numeric()
  only <- 0
  shapes <- matrix(numeric(), 0L, 2L)
  for(configuration in configurations){
    s <- .bt_R2_core_configuration_shapes(core, configuration$open)
    alpha <- s$term_alpha[[term]]
    rest <- sum(s$term_alpha[names(s$term_alpha) != term])
    if(rest <= 0){
      only <- only + configuration$prob
    }else{
      shapes <- rbind(shapes, c(alpha, rest))
      weights <- c(weights, configuration$prob)
    }
  }
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
  configurations <- .bt_R2_core_configurations(core)
  atom <- configurations$null
  specs <- list()
  weights <- numeric()
  for(configuration in configurations$configurations){
    if(!prior$term %in% configuration$open){
      atom <- atom + configuration$prob
      next
    }
    shapes <- .bt_R2_core_configuration_shapes(core, configuration$open)
    open_units <- core$unit_term %in% configuration$open
    unit_alpha <- core$unit_alpha
    if(identical(core$concentration, "fixed_total")){
      unit_alpha[open_units] <- unit_alpha[open_units] * core$fixed_total / sum(unit_alpha[open_units])
    }
    model <- prior
    model$core <- NULL
    model$shapes <- list(a = shapes$a, b = shapes$b, alpha = unit_alpha, A = shapes$A)
    model$n_units <- sum(open_units)
    specs[[length(specs) + 1L]] <- .bt_R2_term_functional(model, L)
    weights <- c(weights, configuration$prob)
  }

  structure(list(atom = atom, specs = specs, weights = weights), class = "BayesTools_R2_functional_mixture")
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
  log(.bt_R2_spec_pdf(spec, x))
}

.bt_R2_spec_quantile <- function(spec, p){

  if(!inherits(spec, "BayesTools_R2_functional_mixture")){
    return(.bt_R2_functional_quantile(spec, p))
  }
  p <- as.numeric(p)
  below <- .bt_R2_spec_cdf(spec, -.Machine$double.xmin)
  at <- below + spec$atom
  vapply(p, function(prob){
    if(is.na(prob)) return(NA_real_)
    if(prob < 0 || prob > 1) .bt_R2_stop("The R-squared functional quantile needs probabilities in [0, 1].")
    if(prob == 0) return(-Inf)
    if(prob == 1) return(Inf)
    if(prob > below && prob <= at) return(0)
    # bracket from the component quantiles of the continuous part
    level <- min(max(prob, 1e-12), 1 - 1e-12)
    ends <- vapply(spec$specs, function(component){
      .bt_R2_functional_quantile(component, c(min(level, 1 - level) / 2, 1 - min(level, 1 - level) / 2))
    }, numeric(2))
    lower <- min(ends[1L, ], -1e-12)
    upper <- max(ends[2L, ], 1e-12)
    while(.bt_R2_spec_cdf(spec, lower) > prob) lower <- 2 * lower - 1
    while(.bt_R2_spec_cdf(spec, upper) < prob) upper <- 2 * upper + 1
    stats::uniroot(function(x) .bt_R2_spec_cdf(spec, x) - prob, c(lower, upper), tol = 1e-12)$root
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

  draws <- .bt_R2_core_draws(prior$core, n)
  directions <- matrix(.bt_R2_direction_rng(prior$law, n * prior$d), nrow = n)
  colnames(directions) <- .bt_R2_direction_coordinates(prior, "R2_term")
  draws <- cbind(draws, directions)
  node <- .bt_dnode_R2_coefficient("R2_term", prior)

  .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(draws, .bt_R2_budget_prior_list(prior)))
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

# The mixed posterior of the realized total of a gated or spiked budget: the
# mixed posterior of its catalog quantity (atom at zero, prior density).
.bt_R2_as_mixed_total <- function(fit, parameter, conditional = FALSE){

  selection <- parameter_catalog_resolve(parameter_catalog(fit), alias = parameter)
  out <- parameter_mixed_posterior(fit, selection, conditional = isTRUE(conditional))
  attr(out, "parameter") <- parameter

  out
}

# The open flag of a gated term per draw from its gate and budget indicator.
.bt_R2_term_open_draws <- function(prior, samples){

  gate <- prior$core$term_gate[[prior$term]]
  indicator <- if(!is.null(prior$gating$indicator)) samples[, prior$gating$indicator] else rep(1, nrow(samples))
  indicator * (if(gate > 0L) samples[, prior$gating$gate_names[[gate]]] else 1)
}

# The atom at zero of the coefficients of a gated term (the term closed): the
# posterior mass of the closed draws and the mixed support with the point 0.
.bt_R2_set_gated_atoms <- function(samples, model_samples, parameter, prior){

  open <- .bt_R2_term_open_draws(prior, model_samples)
  columns <- if(is.null(dim(samples))) parameter else colnames(samples)
  closed <- mean(open == 0)
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
