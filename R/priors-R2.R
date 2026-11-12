# R-squared variance-budget priors for formula terms.
#
# A budget declared by prior_R2() is an entry of a formula prior list whose
# name is not a formula term; the name of the entry names the budget.
# JAGS_formula() binds it to the terms of the formula: the total R-squared r ~
# Beta(a_M, b_M) and the allocation phi ~ Dirichlet(alpha_M) over the
# retained allocation units become ordinary priors of the formula (named after
# the budget), and every retained term becomes an R2 term container whose
# coefficients are deterministic functions of W = r / (1 - r), its
# allocation, and its own direction coordinates prior_par_z_<coefficient>. The
# residual scale enters through the 'multiply_by' attribute of the term, so
# coefficients are in residual-SD units. See ?prior_R2 for the coefficient
# laws and cross-model conventions.

.bt_R2_laws             <- c("exact", "normal", "laplace")
.bt_R2_conventions      <- c("common", "linked", "conditioning")
.bt_R2_concentrations   <- c("inherited", "fixed_total")
.bt_R2_contrasts        <- c("orthonormal", "meandif", "treatment")

#' @title R-squared variance-budget prior for formula terms
#'
#' @description Declares a joint prior distribution for the coefficients of
#' several formula terms through a budget for the explained variance. The
#' total R-squared of the terms has the Beta prior \code{total}, its
#' signal-to-noise ratio \eqn{W = r / (1 - r)} is allocated to the terms by
#' the Dirichlet prior \code{allocation}, and each term turns its share into
#' coefficients by the selected coefficient law. Place the object in a formula
#' prior list under a name that is not a formula term, which names the budget
#' (e.g., \code{R2 = prior_R2(...)}); see Details.
#'
#' @param terms character vector of the formula term labels of the budget in
#' the full candidate model (e.g., \code{c("x", "f", "x:f")}; \code{"x:f"} and
#' \code{"f:x"} name the same term), or \code{NULL} (default) for every
#' fixed-effect term of the formula that has no prior of its own and that no
#' other budget names (at most one budget of a formula; a term with a point
#' prior at zero stays in it, excluded from the model). Terms of the budget
#' that a model's formula does not contain are deleted, and the budget is
#' redistributed in full to the retained terms.
#' @param total Beta prior on the total R-squared of the full candidate model,
#' \code{prior("beta", list(alpha = a, beta = b))} (default Beta(1, 1)), or a
#' spike-and-slab prior with that Beta slab,
#' \code{prior_spike_and_slab(prior("beta", ...), prior_inclusion = ...)},
#' whose spike is zero signal in the budget (the budget null) and whose
#' inclusion prior is a probability prior (\code{prior("spike", list(p))} or a
#' Beta prior). Other mixtures are not supported.
#' @param allocation Dirichlet allocation of the budget of the full candidate
#' model to its terms, from [R2_dirichlet()] (default flat: concentration 1
#' for every term). A budget of one term has the single concentration of its
#' term, which the linked convention and the kernel laws use.
#' @param inclusion inclusion gates of the terms, which make one fit average
#' over the models of the budget: \code{NULL} (default) for no gates; one
#' probability prior for every term, \code{prior("spike", list(p))} for a
#' common fixed inclusion probability or \code{prior("beta", list(a, b))} for
#' a shared inclusion rate (a beta-binomial prior over the models); or a list
#' of probability priors named by term labels for independent gates of the
#' listed terms (the other terms are always included). A Beta prior in the
#' list is marginally its mean probability for its term; a shared rate needs
#' the single-prior form.
#' @param model the fixed model choices of the budget, from [R2_model()]:
#' the coefficient law, the cross-model convention and concentration, the
#' residual scale \code{sigma} and the factor contrast (defaults: exact law,
#' common convention, inherited concentrations, \code{"sigma"}, orthonormal
#' contrasts).
#'
#' @details
#' A formula prior list declares a budget with an entry such as
#' \code{R2 = prior_R2(...)}. The name of the entry names the budget: the
#' total R-squared is the formula parameter \code{<parameter>_<name>}
#' (displayed as \code{<name>{total}}), and the share of the budget allocated
#' to a retained term is displayed as \code{<name>[share: <term>]}; both are
#' quantities of [parameter_catalog()] that the estimates tables show and
#' \code{hypothesis_BF()} tests (e.g., \code{"R2{total} > 0.1"} or
#' \code{"R2[share: f] > 0.5"}), with the Beta prior of the total in the model
#' and the Beta marginal of the share under the allocation. The name
#' must be a valid JAGS node name that contains no \code{"intercept"}, does
#' not end in a suffix of the nodes that BayesTools generates
#' (\code{"_ordered"}, \code{"_indicator"}, \code{"_inclusion"},
#' \code{"_variable"}, \code{"_allocation"}, \code{"_gate"}, \code{"_open"},
#' \code{"_rate"}, \code{"_rate_u"}), and neither it nor the names of its
#' nodes is a term label.
#' Several budgets are several entries (for example, a
#' budget for compulsory covariates and a separate budget for the tested
#' terms, whose R-squared is then the incremental partial R-squared of the
#' tested terms when their design blocks are orthogonal to the covariates;
#' with correlated blocks it is the share of the tested terms' own signal
#' variance, not their partial R-squared). Every term of a budget must not
#' have its own prior in the list, except a point prior at zero
#' (\code{prior("spike", list(0))} or \code{prior_none()}), which excludes it
#' from the model and from the budget. Terms that no budget names keep their
#' own priors.
#'
#' \code{total} and \code{allocation} describe the full candidate model. In a
#' model \eqn{M} with retained allocation units \eqn{u} (the terms for
#' \code{law = "exact"}, the coded coefficients for the kernels, each with
#' its term's concentration) the prior is
#' \deqn{r \sim \textrm{Beta}(a_M, b_M), \quad W = r / (1 - r), \quad
#' \phi \sim \textrm{Dirichlet}(\alpha_M),}
#' with \eqn{(a_M, b_M)} given by the \code{convention} and \eqn{\alpha_M} by
#' the \code{concentration} of \code{model}. With \code{convention = "linked"},
#' the first shape of \code{total} must equal the total concentration
#' \eqn{A} of the full model (the default total is not changed: the error
#' names the required total):
#' the sum of the concentrations of \code{allocation} for the exact law or
#' with \code{concentration = "fixed_total"}, and the sum of the
#' concentrations times the numbers of coded coefficients of the terms for
#' the kernels with inherited concentrations (checked when
#' \code{JAGS_formula()} binds the budget, which then needs the data of all
#' terms of the budget).
#'
#' With the centered design block \eqn{X_t} of a term and its empirical metric
#' \eqn{G_t = X_t^\top C X_t / n}, the coefficients in residual-SD units are
#' \itemize{
#'   \item \code{"exact"}: \eqn{\beta_t = \sqrt{W \phi_t}\, R_t^{-1} u_t},
#'   with \eqn{R_t^\top R_t = G_t} and \eqn{u_t} uniform on the unit sphere,
#'   so the term's signal variance \eqn{\beta_t^\top G_t \beta_t} equals
#'   \eqn{W \phi_t} for every draw and the law of its fitted contribution does
#'   not depend on the contrast;
#'   \item \code{"normal"}: \eqn{\beta_j \sim N(0, W \phi_j / G_{t,jj})};
#'   \item \code{"laplace"}: \eqn{\beta_j \sim
#'   \textrm{Laplace}(0, \sqrt{W \phi_j / (2 G_{t,jj})})}.
#' }
#' When the term blocks are mutually orthogonal, the realized R-squared of the
#' budget equals \eqn{r} for every draw under \code{"exact"}. With correlated
#' terms, the signal variance of the whole budget differs from \eqn{W} by the
#' cross products of the terms' contributions. A model with a single
#' retained unit has \eqn{\phi = 1}; a model whose budget retains no term has
#' no signal from it.
#'
#' The prior targets the empirical variance of the fitted contributions,
#' computed with denominator \eqn{n} in the fitted data; unequal cell counts
#' therefore weight the levels of a factor by their frequencies. With
#' treatment contrasts, the intercept prior applies to the reference level
#' rather than to the grand mean.
#'
#' A budget without \code{terms} takes the terms of each formula it binds to,
#' so that formula is its full candidate model. This is what a single fit
#' with inclusion gates needs; models of a model space fitted separately
#' under the \code{"linked"} or \code{"conditioning"} convention or with
#' \code{concentration = "fixed_total"} need explicit \code{terms} naming
#' the full model, because without them each fit is its own full model.
#' Random-effect terms are never part of a budget without \code{terms}.
#'
#' The R2 term priors that \code{JAGS_formula()} binds support the prior
#' methods (\code{rng()}, \code{mpdf()}, \code{mcdf()}, \code{mquant()},
#' \code{mean()}, \code{var()}, \code{plot()}) of single coefficients, and of
#' factor levels under the exact law, by deterministic quadrature over the
#' budget. The levels of a factor term are also catalog quantities
#' ([parameter_catalog()]) with the prior law of each level (any level under
#' the exact law, also with unequal cell counts; single coded coefficients
#' under the kernels), which \code{parameter_prior_density()} and
#' \code{hypothesis_BF()} use. On the original scale of declared predictor
#' scales (\code{transform_scaled = TRUE}), a coefficient that the scale
#' transformation multiplies by a constant has its rescaled law. The prior
#' densities of combinations of several parameters of one budget (an
#' original-scale coefficient that combines several fitted coefficients
#' among them), ordinates at a point (Savage-Dickey ratios and point
#' hypotheses of coefficients) and the common level marginal of the term
#' methods for the kernels with three or more levels or for unequal marginals
#' are unavailable; point hypotheses at a bound of the support of a share, a
#' total or an inclusion rate are refused (the posterior density of the draws
#' has no exact ordinate there).
#'
#' With \code{inclusion} gates or a spike-and-slab \code{total}, one fit
#' averages over the models of the budget. A configuration \eqn{S} of open
#' terms (the open gated terms and the ungated terms) is a model, and given
#' \eqn{S} the coefficients follow exactly the prior of model \eqn{S} as a
#' separate fit under the same \code{convention}, \code{concentration} and
#' \code{law}: the allocation is redistributed over the open terms (a gate of
#' a kernel term closes all its coded coefficients), the total is
#' Beta(a, b) (\code{"common"}), Beta(\eqn{A_S}, b) (\code{"linked"}, sampled
#' as term signals \eqn{g_t G_t / H} with \eqn{G_t} Gamma(\eqn{\alpha_t}) and
#' \eqn{H} Gamma(b)) or Beta(\eqn{a - (A - A_S)}, \eqn{b + (A - A_S)})
#' (\code{"conditioning"}), and \code{"fixed_total"} rescales the
#' concentrations of the open units. With gates only, the empty configuration
#' is the budget null, with prior probability \eqn{\prod_t (1 - \pi_t)} (or
#' its beta-binomial value). A spike-and-slab total makes the null its own
#' event with the spike weight \eqn{\pi_0}: when every term is gated, the gate
#' prior under the slab is conditioned on a nonempty configuration,
#' \eqn{P(S) = (1 - \pi_0) \pi(S) / (1 - \pi(\emptyset))}; with an ungated term
#' no conditioning is needed. The gates are sampled exactly as a sequence of
#' Bernoulli nodes, each with its probability given the gates before it, so
#' the number of gated terms is not limited. The realized total is zero under
#' the null. The prior methods of the term priors (\code{rng()},
#' \code{mcdf()}, \code{mquant()}, \code{mean()}, \code{var()},
#' \code{plot()}) include the atom at zero of a closed term (\code{mpdf()}
#' is the density of the continuous part), and those of the total entry the
#' atom at zero of the budget null. These prior laws sum over the gate
#' configurations through their open concentrations, exactly and without a
#' limit on the number of gated terms (the configurations are aggregated into
#' at most 65536 distinct states, which few distinct concentrations keep
#' small). A shared Beta(a, b) inclusion rate \eqn{\theta} is the node
#' \code{<parameter>_<name>_rate} and the catalog quantity
#' \code{<name>\{inclusion rate\}}, drawn given the gates as
#' \eqn{\theta | S \sim} Beta(\eqn{a + s}, \eqn{b + K - s}) with \eqn{s} of
#' the \eqn{K} gates open (without changing the sampler of the gates); its
#' prior is Beta(a, b), or, when the gate prior is conditioned on a nonempty
#' configuration, the law of \eqn{\theta} given that condition. Condition
#' events over several terms (and the total) of one gated budget, as in
#' \code{as_mixed_posteriors(conditional = ...)}, use their exact joint prior
#' law. The point mass at zero of a gated coefficient or a total is the
#' exclusion of its term or the null of the budget, tested by the inclusion
#' Bayes factors of \code{runjags_inference_table()} (or, for a total, the
#' region hypothesis \code{"<name>\{total\} > 0"}); the atom at one of a share
#' (its term is the only open term) is tested by region hypotheses such as
#' \code{"<name>[share: x] >= 1"}.
#' Bridge sampling of a gated or spiked budget is unavailable (the
#' fit is itself the model average) and stops with
#' \code{BayesTools_R2_gated_unavailable}.
#'
#' @section Conditions:
#' Invalid budgets and bindings stop with \code{BayesTools_R2_invalid};
#' unavailable results with \code{BayesTools_R2_unavailable} (with
#' \code{BayesTools_R2_unbound} for a term prior used before
#' \code{JAGS_formula()} binds it, and \code{BayesTools_R2_density_unavailable}
#' for an implied distribution); a coefficient state outside the support of
#' the budget with \code{BayesTools_R2_invalid_state} (also
#' \code{BayesTools_marglik_out_of_support}). All inherit
#' \code{BayesTools_R2_condition}. Linear combinations of several parameters of
#' one budget stop with \code{BayesTools_formula_measure_unavailable}.
#'
#' @return \code{prior_R2} returns an object of class
#' \code{"BayesTools_R2_budget"}.
#'
#' @examples
#' # the defaults: every term of the formula without its own prior, a
#' # uniform total R-squared, a flat allocation, exact allocation
#' prior_R2()
#'
#' budget <- prior_R2(
#'   terms      = c("x", "f"),
#'   total      = prior("beta", list(alpha = 1, beta = 1)),
#'   allocation = R2_dirichlet(c(x = 1, f = 2)),
#'   model      = R2_model(law = "exact", convention = "common")
#' )
#' budget
#'
#' # formula prior list of a model 'mu ~ x + f': the entry name 'R2' names
#' # the budget
#' formula_prior_list <- list(
#'   intercept = prior("normal", list(0, 10)),
#'   R2        = budget
#' )
#'
#' \dontrun{
#' # ANCOVA: compare the models with and without the factor f
#' set.seed(1)
#' data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 20)), x = rnorm(60))
#' data$y <- 0.5 * (data$f == "b") + 0.3 * data$x + rnorm(60)
#'
#' fit_model <- function(formula){
#'   arguments <- list(
#'     model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
#'     data = list(y = data$y, N = nrow(data)),
#'     prior_list = list(sigma = prior("invgamma", list(1, 1))),
#'     formula_list = list(mu = formula),
#'     formula_data_list = list(mu = data),
#'     formula_prior_list = list(mu = formula_prior_list)
#'   )
#'   fit <- do.call(JAGS_fit, arguments)
#'   marglik <- do.call(JAGS_bridgesampling, c(list(fit = fit,
#'     log_posterior = function(parameters, data){
#'       sum(dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
#'     }), arguments[-1]))
#'   list(fit = fit, marglik = marglik, prior_weights = 1)
#' }
#'
#' # the same budget in both models: without f, x receives the whole budget
#' models <- models_inference(list(fit_model(~ x), fit_model(~ x + f)))
#' inference <- ensemble_inference(models, parameters = "mu_f",
#'   is_null_list = list(mu_f = c(TRUE, FALSE)), conditional = FALSE)
#' inference$mu_f$BF
#'
#' # posterior summaries in residual-SD units, with factor levels
#' runjags_estimates_table(models[[2]]$fit, transform_factors = TRUE)
#'
#' # implied prior distributions: the total R-squared and a coefficient
#' priors <- attr(models[[2]]$fit, "prior_list")
#' plot(priors$mu_R2)
#' plot(priors$mu_x)
#'
#' # one fit averaging over the models of the budget: a shared inclusion rate
#' # (beta-binomial model prior) and a spike-and-slab total for the budget null
#' gated <- prior_R2(
#'   total     = prior_spike_and_slab(prior("beta", list(1, 1)),
#'                                    prior_inclusion = prior("spike", list(0.5))),
#'   inclusion = prior("beta", list(1, 1)),
#'   model     = R2_model(contrast = "meandif")
#' )
#' fit <- JAGS_fit(
#'   model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
#'   data = list(y = data$y, N = nrow(data)),
#'   prior_list = list(sigma = prior("invgamma", list(1, 1))),
#'   formula_list = list(mu = ~ x + f), formula_data_list = list(mu = data),
#'   formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 10)), R2 = gated))
#' )
#' runjags_inference_table(fit)   # inclusion of the budget and of each term
#' runjags_estimates_table(fit, conditional = TRUE)
#' }
#'
#' @seealso [R2_model()], [R2_dirichlet()], [JAGS_formula()], [prior_factor()]
#' @export
prior_R2 <- function(terms = NULL, total = prior("beta", list(1, 1)), allocation = R2_dirichlet(),
                     inclusion = NULL, model = R2_model()){

  check_R2_model(model)
  if(!is.null(terms)){
    terms <- .bt_R2_check_terms(terms)
  }
  shapes <- .bt_R2_check_total(total)
  .bt_R2_check_allocation_spec(allocation)
  .bt_R2_check_inclusion(inclusion, terms)

  output <- list(
    terms         = terms,
    alpha         = NULL,
    a             = shapes[["a"]],
    b             = shapes[["b"]],
    total_inclusion = shapes[["inclusion"]],
    inclusion     = NULL,
    allocation    = allocation,
    inclusion_spec = inclusion,
    implicit      = is.null(terms),
    law           = model$law,
    convention    = model$convention,
    concentration = model$concentration,
    sigma         = model$sigma,
    contrast      = model$contrast,
    model         = model,
    name          = NULL
  )
  class(output) <- c("BayesTools_R2_budget", "list")
  if(is.null(terms)){
    return(output)
  }

  .bt_R2_resolve_budget(output, terms)
}

# The budget 'budget' over the terms 'terms' (given to prior_R2() or resolved
# when JAGS_formula() binds a budget without terms): the concentrations of the
# allocation and the inclusion gates of the terms, and the checks of the
# shapes that the concentrations decide.
.bt_R2_resolve_budget <- function(budget, terms){

  budget$terms <- terms
  budget$alpha <- stats::setNames(.bt_R2_allocation_alpha(budget$allocation, terms, positional = !budget$implicit), terms)
  budget$inclusion <- .bt_R2_check_inclusion(budget$inclusion_spec, terms)

  # the full-model total concentration is known here unless a kernel law
  # inherits the concentrations (then it needs the term dimensions)
  if(identical(budget$convention, "linked") &&
     (identical(budget$law, "exact") || identical(budget$concentration, "fixed_total"))){
    .bt_R2_check_linked_shape(budget, sum(budget$alpha))
  }
  alpha <- budget$alpha
  if(identical(budget$convention, "conditioning") && identical(budget$concentration, "inherited") &&
     identical(budget$law, "exact") && budget$a - sum(alpha) + min(alpha) <= 0){
    .bt_R2_stop(paste0("The 'conditioning' convention needs a - sum(alpha) + min(alpha) > 0, with a the first shape of ",
         "'total' and alpha the concentrations of 'allocation', so that every model has a positive first shape",
         if(!is.null(budget$name)) paste0(" (budget '", budget$name, "')"), "."))
  }

  budget
}

#' @title Fixed model choices of an R-squared budget prior
#'
#' @description \code{R2_model()} declares the fixed choices of a
#' [prior_R2()] budget: the coefficient law, how a model of a model space
#' derives its prior from the full candidate model, the residual scale, and
#' the factor contrast. These are fixed model choices, not parameters with
#' prior distributions; one specification can serve several budgets.
#'
#' @param law coefficient law: \code{"exact"} (default) allocates squared
#' signal magnitudes exactly, with a fair sign for one-column terms and a
#' uniform direction in the term's coefficient space (in the empirical design
#' metric) for multi-column terms; \code{"normal"} is the normal kernel of the
#' R2D2M2 prior and \code{"laplace"} the double-exponential kernel of the R2D2
#' prior, both allocating the budget over the individual coded coefficients,
#' each with the concentration of its term.
#' @param convention how a model derives the Beta shapes of its total from
#' the \code{total} = Beta(a, b) of the full candidate model, with \eqn{A} the
#' total concentration of the full model and \eqn{A_M} the total concentration
#' retained by model \eqn{M}: \code{"common"} (default) uses Beta(a, b) in
#' every model; \code{"linked"} uses Beta(\eqn{A_M}, b), so the first shape of
#' \code{total} must equal \eqn{A}; \code{"conditioning"} keeps
#' \eqn{a_M - A_M} and \eqn{a_M + b_M} equal to \eqn{a - A} and \eqn{a + b}.
#' @param concentration concentrations of the retained units:
#' \code{"inherited"} (default) keeps the concentrations of the allocation;
#' \code{"fixed_total"} rescales them in every model to the sum of the
#' concentrations of the allocation.
#' @param sigma residual standard deviation that multiplies the coefficients,
#' as the \code{"multiply_by"} attribute of a formula prior: the name of a
#' parameter of the model (default \code{"sigma"}) or a positive number.
#' @param contrast contrast of the factor terms of the budget:
#' \code{"orthonormal"} (default), \code{"meandif"}, or \code{"treatment"};
#' used only when the budget contains a factor term.
#' @param model an R-squared budget specification to validate.
#' @param name argument name used in validation messages.
#' @param x an R-squared budget specification.
#' @param ... additional arguments, currently unused.
#'
#' @return \code{R2_model()} returns an object of class \code{"R2_model"}
#' with the five choices; \code{check_R2_model()} and \code{print()} return
#' their validated input invisibly.
#'
#' @examples
#' model <- R2_model(law = "normal", convention = "linked")
#' model
#' # one specification for two budgets of a formula prior list
#' list(
#'   covariates = prior_R2(terms = c("age", "sex"), model = model,
#'                         total = prior("beta", list(2, 2))),
#'   tested     = prior_R2(terms = "treatment", model = model)
#' )
#'
#' @seealso [prior_R2()], [R2_dirichlet()]
#' @export
R2_model <- function(law = "exact", convention = "common", concentration = "inherited",
                     sigma = "sigma", contrast = "orthonormal"){

  output <- structure(list(
    law           = law,
    convention    = convention,
    concentration = concentration,
    sigma         = sigma,
    contrast      = contrast
  ), class = "R2_model")
  check_R2_model(output)
  output
}

#' @rdname R2_model
#' @export
check_R2_model <- function(model, name = "model"){

  fields <- c("law", "convention", "concentration", "sigma", "contrast")
  if(!is.list(model) || !inherits(model, "R2_model") || !identical(names(model), fields)){
    stop(paste0("'", name, "' must be a specification from 'R2_model()'."), call. = FALSE)
  }
  check_char(model$law, "law", allow_values = .bt_R2_laws, allow_NA = FALSE)
  check_char(model$convention, "convention", allow_values = .bt_R2_conventions, allow_NA = FALSE)
  check_char(model$concentration, "concentration", allow_values = .bt_R2_concentrations, allow_NA = FALSE)
  check_char(model$contrast, "contrast", allow_values = .bt_R2_contrasts, allow_NA = FALSE)
  sigma <- model$sigma
  if(is.character(sigma)){
    check_char(sigma, "sigma", allow_NA = FALSE)
    .bt_check_jags_node_name(sigma, "sigma")
  }else{
    check_real(sigma, "sigma", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  }
  invisible(model)
}

#' @rdname R2_model
#' @export
print.R2_model <- function(x, ...){

  check_R2_model(x, "x")
  cat("R-squared budget model\n")
  cat("  law:           ", x$law, "\n", sep = "")
  cat("  convention:    ", x$convention, "\n", sep = "")
  cat("  concentration: ", x$concentration, "\n", sep = "")
  cat("  sigma:         ", if(is.character(x$sigma)) paste0("'", x$sigma, "'") else format(x$sigma), "\n", sep = "")
  cat("  contrast:      ", x$contrast, " (factor terms only)\n", sep = "")
  invisible(x)
}

#' @title Dirichlet allocation of an R-squared budget
#'
#' @description \code{R2_dirichlet()} specifies the Dirichlet allocation of a
#' [prior_R2()] budget over its terms, with one concentration per term.
#'
#' @param alpha \code{NULL} (default) for the flat allocation (concentration 1
#' for every term, sized when the budget binds its terms); a single value for
#' a symmetric allocation; a vector named by the terms of the budget, which it
#' must cover exactly; or an unnamed vector in the order of the explicit
#' \code{terms} of the budget. A budget of one term takes a single value, its
#' term's concentration (which the linked convention and the kernel laws
#' use). Each value must be finite and at least \code{0.01}.
#' @param x an allocation specification.
#' @param ... additional arguments, currently unused.
#'
#' @return \code{R2_dirichlet()} returns an object of class
#' \code{"R2_dirichlet"}.
#'
#' @examples
#' R2_dirichlet()                      # flat
#' R2_dirichlet(0.5)                   # symmetric, sparser allocations
#' R2_dirichlet(c(x = 1, f = 2))       # per term
#'
#' @seealso [prior_R2()], [R2_model()]
#' @export
R2_dirichlet <- function(alpha = NULL){

  if(!is.null(alpha)){
    check_real(alpha, "alpha", lower = 0.01, allow_bound = TRUE, check_length = 0, allow_NA = FALSE)
    if(any(!is.finite(alpha))){
      stop("The 'alpha' argument must be finite.", call. = FALSE)
    }
    labels <- names(alpha)
    if(!is.null(labels)){
      if(anyNA(labels) || any(!nzchar(trimws(labels)))){
        stop("The concentrations in 'alpha' must be all named by terms or all unnamed.", call. = FALSE)
      }
      if(anyDuplicated(.bt_R2_term_key(trimws(labels)))){
        stop("'alpha' names each term once.", call. = FALSE)
      }
    }
  }

  structure(list(alpha = alpha), class = "R2_dirichlet")
}

#' @rdname R2_dirichlet
#' @export
print.R2_dirichlet <- function(x, ...){

  alpha <- x$alpha
  cat("R-squared budget allocation: ",
    if(is.null(alpha)){
      "flat Dirichlet (concentration 1 for every term)"
    }else if(length(alpha) == 1L && is.null(names(alpha))){
      paste0("symmetric Dirichlet (concentration ", format(alpha), " for every term)")
    }else if(is.null(names(alpha))){
      paste0("Dirichlet(", paste(format(alpha), collapse = ", "), ") in the order of the terms")
    }else{
      paste0("Dirichlet(", paste0(names(alpha), " = ", format(alpha), collapse = ", "), ")")
    }, "\n", sep = "")
  invisible(x)
}

.bt_R2_check_terms <- function(terms){

  if(!is.character(terms) || length(terms) == 0L || anyNA(terms) || any(!nzchar(trimws(terms)))){
    stop("'terms' must be a character vector of the formula term labels of the budget.", call. = FALSE)
  }
  terms <- trimws(terms)
  keys <- .bt_R2_term_key(terms)
  if(anyDuplicated(keys)){
    stop("'terms' names each formula term once; duplicated: ",
         paste0("'", terms[duplicated(keys)], "'", collapse = ", "), ".", call. = FALSE)
  }
  if(any(grepl("intercept", terms, fixed = TRUE))){
    stop("The intercept cannot be part of an R-squared budget.", call. = FALSE)
  }

  terms
}

# The Beta shapes of the total prior and, for a spike-and-slab total, the law
# of its inclusion probability.
.bt_R2_check_total <- function(total){

  if(is.prior(total) && is.prior.spike_and_slab(total)){
    shapes <- .bt_R2_check_beta_total(.get_spike_and_slab_variable(total))
    inclusion <- .bt_R2_probability_law(.get_spike_and_slab_inclusion(total),
      "The inclusion prior of the spike-and-slab 'total'")
    return(list(a = shapes[["a"]], b = shapes[["b"]], inclusion = inclusion))
  }
  if(is.prior(total) && is.prior.mixture(total)){
    stop("'total' must be a Beta prior or a spike-and-slab prior with a Beta slab, ",
         "prior_spike_and_slab(prior(\"beta\", ...)); other mixtures are not supported.", call. = FALSE)
  }
  shapes <- .bt_R2_check_beta_total(total)

  list(a = shapes[["a"]], b = shapes[["b"]], inclusion = NULL)
}

.bt_R2_check_beta_total <- function(total){

  if(!is.prior(total) || !identical(total[["distribution"]], "beta") || !is.prior.simple(total)){
    stop("'total' must be a Beta prior on the total R-squared, prior(\"beta\", list(alpha = a, beta = b)).", call. = FALSE)
  }
  if(.is_prior_expression(total)){
    stop("The Beta prior 'total' must have numeric shapes; expressions are not supported.", call. = FALSE)
  }
  truncation <- unname(unlist(total[["truncation"]]))
  if(length(truncation) != 2L || !identical(as.numeric(truncation), c(0, 1))){
    stop("The Beta prior 'total' cannot be truncated.", call. = FALSE)
  }

  c(a = unname(total$parameters[["alpha"]]), b = unname(total$parameters[["beta"]]))
}

# The law of an inclusion probability: a fixed probability in (0, 1) (a spike
# prior) or a Beta prior, whose mean is its marginal inclusion probability.
.bt_R2_probability_law <- function(prior, what){

  if(!is.prior(prior) || is.prior.mixture(prior) || !is.prior.simple(prior) || .is_prior_expression(prior)){
    stop(what, " must be a probability prior: prior(\"spike\", list(p)) or prior(\"beta\", list(a, b)).", call. = FALSE)
  }
  if(is.prior.point(prior)){
    p <- unname(prior$parameters[["location"]])
    if(!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1){
      stop(what, " must be a probability strictly between 0 and 1 (a term that is always included takes no gate).",
           call. = FALSE)
    }
    return(list(type = "spike", p = p))
  }
  if(!identical(prior[["distribution"]], "beta") ||
     !identical(as.numeric(unname(unlist(prior[["truncation"]]))), c(0, 1))){
    stop(what, " must be a probability prior: prior(\"spike\", list(p)) or an untruncated prior(\"beta\", list(a, b)).",
         call. = FALSE)
  }
  a <- unname(prior$parameters[["alpha"]])
  b <- unname(prior$parameters[["beta"]])

  list(type = "beta", a = a, b = b, p = a / (a + b))
}

# The inclusion gates of the budget terms: NULL (no gates), one probability
# law for every term ("common": a fixed probability, or a shared Beta rate
# that makes the model prior beta-binomial), or per-term laws ("terms":
# independent gates with their marginal probabilities; terms not listed are
# always included).
.bt_R2_check_inclusion <- function(inclusion, terms){

  if(is.null(inclusion)){
    return(NULL)
  }
  if(is.prior(inclusion)){
    return(list(form = "common", law = .bt_R2_probability_law(inclusion, "'inclusion'")))
  }
  if(!is.list(inclusion) || length(inclusion) == 0L || is.null(names(inclusion)) ||
     anyNA(names(inclusion)) || any(!nzchar(trimws(names(inclusion))))){
    stop("'inclusion' must be NULL, one probability prior for every term, or a list of probability priors ",
         "named by the term labels.", call. = FALSE)
  }
  keys <- .bt_R2_term_key(trimws(names(inclusion)))
  if(anyDuplicated(keys)){
    stop("'inclusion' names each term once; duplicated: ",
         paste0("'", names(inclusion)[duplicated(keys)], "'", collapse = ", "), ".", call. = FALSE)
  }
  # without explicit terms the names are checked when the budget binds
  if(is.null(terms)){
    for(i in seq_along(inclusion)){
      .bt_R2_probability_law(inclusion[[i]], paste0("The inclusion prior of '", names(inclusion)[[i]], "'"))
    }
    return(invisible(NULL))
  }
  position <- match(keys, .bt_R2_term_key(terms))
  if(anyNA(position)){
    .bt_R2_stop(paste0("'inclusion' names terms that are not terms of the budget: ",
         paste0("'", names(inclusion)[is.na(position)], "'", collapse = ", "), "."))
  }
  laws <- lapply(seq_along(inclusion), function(i){
    .bt_R2_probability_law(inclusion[[i]], paste0("The inclusion prior of '", terms[[position[[i]]]], "'"))
  })
  names(laws) <- terms[position]

  list(form = "terms", laws = laws[terms[sort(position)]])
}

# The allocation of a budget is an R2_dirichlet() specification.
.bt_R2_check_allocation_spec <- function(allocation){

  if(!inherits(allocation, "R2_dirichlet")){
    stop("'allocation' must be an allocation from 'R2_dirichlet()' (e.g., R2_dirichlet(c(x = 1, f = 2))); ",
         "prior(\"dirichlet\", ...) is not accepted for a budget.", call. = FALSE)
  }
  invisible(allocation)
}

# The concentrations of the allocation of a budget, one per term: flat for
# NULL, a single value for every term, per-term values matched by name (they
# must cover the terms exactly), or values in the order of explicit terms
# ('positional').
.bt_R2_allocation_alpha <- function(allocation, terms, positional = TRUE){

  alpha <- allocation$alpha
  if(is.null(alpha)){
    return(rep(1, length(terms)))
  }
  labels <- names(alpha)
  if(is.null(labels) && length(alpha) == 1L){
    return(rep(as.numeric(alpha), length(terms)))
  }
  if(!is.null(labels)){
    keys <- .bt_R2_term_key(trimws(labels))
    position <- match(.bt_R2_term_key(terms), keys)
    if(length(keys) != length(terms) || anyNA(position)){
      .bt_R2_stop(paste0("The concentrations of 'allocation' must be named by the terms of the budget, each once; ",
        "terms: ", paste0("'", terms, "'", collapse = ", "), "; names: ",
        paste0("'", labels, "'", collapse = ", "), "."))
    }
    return(as.numeric(alpha[position]))
  }
  if(!positional){
    .bt_R2_stop(paste0("Unnamed concentrations of 'allocation' need the explicit 'terms' of the budget; ",
      "name them by the terms (", paste0("'", terms, "'", collapse = ", "), ") or give one value for every term."))
  }
  if(length(alpha) != length(terms)){
    .bt_R2_stop(paste0("'allocation' has ", length(alpha), " concentrations, but 'terms' has ",
      length(terms), " term(s)."))
  }

  as.numeric(alpha)
}

# The first shape of a linked total equals the total concentration of the
# full candidate model.
.bt_R2_check_linked_shape <- function(budget, A){

  if(abs(budget$a - A) > sqrt(.Machine$double.eps) * max(1, A)){
    .bt_R2_stop(paste0(
      "With convention = \"linked\", the first shape of 'total' must equal the total concentration of the full model, ",
      as.character(A), " for this budget",
      if(!is.null(budget$name)) paste0(" '", budget$name, "'"),
      "; use total = prior(\"beta\", list(alpha = ", as.character(A), ", beta = ", as.character(budget$b), "))."))
  }
  invisible(TRUE)
}

#' @export
print.BayesTools_R2_budget <- function(x, ...){

  law_label <- switch(
    x$law,
    exact   = "exact allocation",
    normal  = "normal kernel (R2D2M2)",
    laplace = "Laplace kernel (R2D2)"
  )
  total <- print(prior("beta", list(alpha = x$a, beta = x$b)), silent = TRUE)
  law_text <- function(law){
    if(identical(law$type, "spike")) format(law$p) else paste0("Beta(", format(law$a), ", ", format(law$b), ")")
  }
  total_rule <- switch(
    x$convention,
    common       = "in every model",
    linked       = paste0("in the full model; Beta(A_M, ", format(x$b), ") in a model with retained concentration A_M"),
    conditioning = "in the full model; a - A and a + b kept in every model"
  )
  allocation <- if(!is.null(x$alpha)){
    paste0("Dirichlet(", paste(format(x$alpha), collapse = ", "), ")")
  }else if(is.null(x$allocation$alpha)){
    "flat Dirichlet"
  }else if(length(x$allocation$alpha) == 1L && is.null(names(x$allocation$alpha))){
    paste0("symmetric Dirichlet(", format(x$allocation$alpha), ")")
  }else{
    paste0("Dirichlet(", paste0(names(x$allocation$alpha), " = ", format(x$allocation$alpha), collapse = ", "), ")")
  }
  allocation_rule <- if(identical(x$concentration, "fixed_total")){
    paste0("retained concentrations rescaled to sum to ",
           if(is.null(x$alpha)) "their full-model total" else format(sum(x$alpha)), " in every model")
  }else{
    "retained terms keep their concentrations"
  }
  if(!identical(x$law, "exact")){
    allocation_rule <- paste0(allocation_rule, "; each coded coefficient has its term's concentration")
  }
  scale <- if(is.character(x$sigma)){
    paste0("coefficients in units of '", x$sigma, "'")
  }else{
    paste0("coefficients in units of a fixed residual SD of ", format(x$sigma))
  }

  cat("R-squared budget prior", if(!is.null(x$name)) paste0(" '", x$name, "'"), " (", law_label, ")\n", sep = "")
  cat("  terms:      ", if(is.null(x$terms)){
    "every formula term without its own prior (resolved when the formula binds the budget)"
  }else{
    paste0(paste(x$terms, collapse = ", "), if(isTRUE(x$implicit)) " (resolved from the formula)")
  }, "\n", sep = "")
  cat("  total:      ", total, " ", total_rule, "\n", sep = "")
  if(!is.null(x$total_inclusion)){
    cat("  spike:      zero signal in the budget with probability ",
      if(identical(x$total_inclusion$type, "spike")) format(1 - x$total_inclusion$p) else
        paste0("1 - ", law_text(x$total_inclusion), " (marginally ", format(1 - x$total_inclusion$p), ")"),
      "\n", sep = "")
  }
  if(!is.null(x$inclusion)){
    gates <- if(identical(x$inclusion$form, "common")){
      if(identical(x$inclusion$law$type, "spike")){
        paste0("every term with probability ", format(x$inclusion$law$p))
      }else{
        paste0("every term with a shared rate ~ ", law_text(x$inclusion$law), " (beta-binomial model prior)")
      }
    }else{
      paste0(paste0(names(x$inclusion$laws), ": ", vapply(x$inclusion$laws, law_text, character(1)), collapse = ", "),
        if(length(x$inclusion$laws) < length(x$terms)) "; other terms always included")
    }
    cat("  inclusion:  ", gates, "\n", sep = "")
  }else if(is.list(x$inclusion_spec) && !is.prior(x$inclusion_spec)){
    laws <- lapply(seq_along(x$inclusion_spec), function(i){
      .bt_R2_probability_law(x$inclusion_spec[[i]], "'inclusion'")
    })
    cat("  inclusion:  ", paste0(names(x$inclusion_spec), ": ", vapply(laws, law_text, character(1)), collapse = ", "),
        " (terms matched when the formula binds the budget)\n", sep = "")
  }else if(is.prior(x$inclusion_spec)){
    law <- .bt_R2_probability_law(x$inclusion_spec, "'inclusion'")
    cat("  inclusion:  every term with ", if(identical(law$type, "spike")) paste0("probability ", format(law$p)) else
      paste0("a shared rate ~ ", law_text(law), " (beta-binomial model prior)"), "\n", sep = "")
  }
  if(is.null(x$terms) || length(x$terms) > 1L){
    cat("  allocation: ", allocation, ", ", allocation_rule, "\n", sep = "")
  }else{
    cat("  allocation: concentration ", format(x$alpha), " of the single term",
        if(!identical(x$law, "exact")){
          if(identical(x$concentration, "fixed_total")) " (shared by its coded coefficients under the kernel laws)" else
            " (of each coded coefficient under the kernel laws)"
        }, "\n", sep = "")
  }
  cat("  scale:      ", scale, "\n", sep = "")
  cat("  contrast:   ", x$contrast, " (factor terms only)\n", sep = "")
  invisible(x)
}

#' @title Check whether an object is an R-squared budget or an R2 term prior
#'
#' @param x an object.
#' @return \code{is.prior.R2} returns \code{TRUE} for the term priors that
#' \code{JAGS_formula()} binds from a [prior_R2()] budget, and
#' \code{is.R2_budget} for budgets created by [prior_R2()].
#' @export
is.prior.R2 <- function(x){
  inherits(x, "prior.R2")
}
#' @rdname is.prior.R2
#' @export
is.R2_budget <- function(x){
  inherits(x, "BayesTools_R2_budget")
}

.bt_R2_stop <- function(message, class = "BayesTools_R2_invalid"){
  stop(structure(list(message = message, call = NULL),
    class = unique(c(class, "BayesTools_R2_condition", "error", "condition"))))
}

# A coefficient state outside the support of the budget (a total outside
# [0, 1) or invalid allocation coordinates): bridge sampling treats it as a
# state of zero density.
.bt_R2_stop_state <- function(message){
  .bt_R2_stop(message, class = c("BayesTools_R2_invalid_state", "BayesTools_marglik_out_of_support"))
}

.bt_R2_unavailable <- function(message, class = NULL){
  .bt_R2_stop(message, class = c(class, "BayesTools_R2_unavailable"))
}

# Canonical key of a term label: components sorted, so that 'f:x' and 'x:f'
# name the same term.
.bt_R2_term_key <- function(labels){

  labels <- gsub("__xXx__", ":", labels, fixed = TRUE)
  vapply(labels, function(label){
    paste(sort(trimws(strsplit(label, ":", fixed = TRUE)[[1]])), collapse = ":")
  }, character(1), USE.NAMES = FALSE)
}

.bt_R2_direction_name <- function(coefficient_name){
  paste0("prior_par_z_", coefficient_name)
}

# ---------------------------------------------------------------------------
# Formula binding

# Suffixes of the nodes generated for budgets and other priors, which a budget
# name cannot end in.
.bt_R2_generated_suffixes <- c("_ordered", "_indicator", "_inclusion", "_variable", "_allocation",
                               "_gate", "_open", "_rate", "_rate_u")

# The entries of a formula prior list that are R-squared budgets.
.bt_R2_budget_entries <- function(prior_list){

  vapply(prior_list, is.R2_budget, logical(1))
}

# The budgets of a formula prior list ('entries', the entries that are
# budgets), each named by its entry (NULL when there are none). An entry name
# names the budget's nodes, so it must be a JAGS node name that is neither a
# formula term nor a term of a budget. A budget without terms (at most one)
# takes every fixed-effect term of the formula that has no prior of its own in
# 'prior_list' (a point prior at zero keeps the term in the budget, excluded
# from the model) and that no other budget names.
.bt_R2_formula_budgets <- function(entries, model_terms, prior_list = list()){

  if(length(entries) == 0L){
    return(NULL)
  }
  budget_names <- names(entries)
  if(is.null(budget_names) || anyNA(budget_names) || any(!nzchar(budget_names))){
    .bt_R2_stop("An R-squared budget must be a named entry of the formula prior list; the name names the budget.")
  }
  if(anyDuplicated(budget_names)){
    .bt_R2_stop(paste0("The R-squared budgets of a formula need distinct names; duplicated: ",
      paste0("'", unique(budget_names[duplicated(budget_names)]), "'", collapse = ", "), "."))
  }
  implicit <- budget_names[vapply(entries, function(budget) is.null(budget$terms), logical(1))]
  if(length(implicit) > 1L){
    .bt_R2_stop(paste0("At most one R-squared budget of a formula can take the formula's terms without 'terms'; ",
      "name the terms of the budgets ", paste0("'", implicit, "'", collapse = ", "), "."))
  }
  if(length(implicit) == 1L){
    claimed <- unlist(lapply(entries[setdiff(budget_names, implicit)], function(budget) .bt_R2_term_key(budget$terms)))
    candidates <- setdiff(model_terms, "intercept")
    own <- vapply(candidates, function(term){
      term %in% names(prior_list) && !.bt_R2_is_exclusion_prior(prior_list[[term]])
    }, logical(1))
    terms <- candidates[!own & !.bt_R2_term_key(candidates) %in% claimed]
    if(length(terms) == 0L){
      .bt_R2_stop(paste0("The R-squared budget '", implicit, "' has no terms: every term of the formula has its own ",
        "prior or belongs to another budget."))
    }
    budget <- entries[[implicit]]
    budget$name <- implicit
    entries[[implicit]] <- .bt_R2_resolve_budget(budget, gsub("__xXx__", ":", terms, fixed = TRUE))
  }
  term_keys <- unlist(lapply(entries, function(budget) .bt_R2_term_key(budget$terms)))
  if(anyDuplicated(term_keys)){
    .bt_R2_stop(paste0("A formula term can belong to only one R-squared budget; repeated: ",
      paste0("'", unique(term_keys[duplicated(term_keys)]), "'", collapse = ", "), "."))
  }
  term_names <- c(.bt_R2_term_key(model_terms), term_keys)
  for(name in budget_names){
    tryCatch(.bt_check_jags_node_name(name, "name"), error = function(e){
      .bt_R2_stop(paste0("The name of the R-squared budget '", name, "' must be a valid JAGS node name: ",
        conditionMessage(e)))
    })
    if(grepl("intercept", name, fixed = TRUE)){
      .bt_R2_stop(paste0("The name of the R-squared budget '", name, "' cannot contain 'intercept'."))
    }
    # the suffixes of the nodes that BayesTools generates for other priors
    # (ordered totals and allocations, mixture indicators, allocations) and
    # for gated budgets (gates, open flags, the shared inclusion rate)
    reserved <- .bt_R2_generated_suffixes
    if(any(endsWith(name, reserved))){
      .bt_R2_stop(paste0("The name of the R-squared budget '", name, "' cannot end in ",
        paste0("'", reserved, "'", collapse = ", "), ", which name generated nodes."))
    }
    if(.bt_R2_term_key(name) %in% term_names){
      .bt_R2_stop(paste0("The R-squared budget '", name, "' is named like a formula term; ",
        "a budget is an entry of the formula prior list under a name that is not a term."))
    }
    generated <- paste0(name, .bt_R2_generated_suffixes)
    if(any(.bt_R2_term_key(generated) %in% term_names)){
      .bt_R2_stop(paste0("The nodes of the R-squared budget '", name, "' would be named like the formula term '",
        generated[.bt_R2_term_key(generated) %in% term_names][1L], "'; rename the budget."))
    }
  }
  budgets <- entries
  for(name in budget_names){
    budgets[[name]]$name <- name
  }

  budgets
}

.bt_R2_is_exclusion_prior <- function(prior){

  if(is.prior.none(prior)){
    return(TRUE)
  }
  location <- prior$parameters[["location"]]
  is.prior.point(prior) && !.is_prior_expression(prior) && is.numeric(location) &&
    length(location) >= 1L && all(location == 0)
}

# An excluded factor term keeps a zero spike coded by the budget's contrast
# (a scalar spike or prior_none() is not a factor prior).
.bt_R2_exclusion_prior <- function(prior, budget, term, type){

  if(!identical(type, "factor") || is.prior.factor(prior)){
    return(prior)
  }
  prior_factor("spike", list(0), contrast = budget$contrast, prior_weights = .prior_model_weight(prior))
}

# Assigns the formula terms to the budgets. Retained budget terms get an
# unbound R2 term prior in 'prior_list' (so that defaults do not fill them);
# terms with a point prior at zero stay excluded.
.bt_R2_assign_terms <- function(budgets, prior_list, model_terms, model_terms_type){

  terms <- setdiff(model_terms, "intercept")
  term_keys <- .bt_R2_term_key(terms)
  for(budget in budgets){
    budget_keys <- .bt_R2_term_key(budget$terms)
    for(i in which(term_keys %in% budget_keys)){
      term <- terms[[i]]
      declared <- budget$terms[match(term_keys[[i]], budget_keys)]
      if(term %in% names(prior_list)){
        if(.bt_R2_is_exclusion_prior(prior_list[[term]])){
          prior_list[[term]] <- .bt_R2_exclusion_prior(prior_list[[term]], budget, term, model_terms_type[[term]])
          next
        }
        .bt_R2_stop(paste0("The formula term '", term, "' belongs to the R-squared budget '", budget$name,
          "' and cannot have its own prior. Remove its prior, or give it a point prior at zero ",
          "(prior(\"spike\", list(0)) or prior_none()) to exclude it from the model and the budget."))
      }
      prior_list[[term]] <- .bt_R2_term_stub(budget, term, declared, model_terms_type[[term]])
    }
  }
  prior_list
}

.bt_R2_term_stub <- function(budget, term, declared, type){

  output <- list(
    distribution = "R2",
    parameters   = list(),
    truncation   = list(lower = -Inf, upper = Inf),
    budget       = budget,
    term         = term,
    declared     = declared,
    alpha        = unname(budget$alpha[[declared]]),
    law          = budget$law,
    bound        = FALSE,
    prior_weights = 1
  )
  # Like ordered priors, R2 term priors carry the classes of the ordinary
  # priors of their term (scalar; treatment factors are scalar factor priors,
  # orthonormal and mean-difference factors vector factor priors) for the
  # consumers of term metadata (labels, coordinates, factor transforms,
  # mixing); every distribution-specific consumer dispatches on 'prior.R2'
  # first.
  class(output) <- if(identical(type, "factor")){
    c("prior", if(identical(budget$contrast, "treatment")) "prior.simple" else "prior.vector",
      "prior.factor", paste0("prior.", budget$contrast), "prior.R2")
  }else{
    c("prior", "prior.simple", "prior.R2")
  }
  attr(output, "multiply_by") <- budget$sigma
  output
}

# Number of design columns of every declared term of a budget in the full
# candidate model (all formula terms and all declared terms).
.bt_R2_declared_dimensions <- function(budget, formula, data){

  declared <- budget$terms
  formula_labels <- attr(stats::terms(formula), "term.labels")
  missing <- setdiff(.bt_R2_term_key(declared), .bt_R2_term_key(formula_labels))
  full_labels <- c(formula_labels, declared[.bt_R2_term_key(declared) %in% missing])
  variables <- unique(unlist(strsplit(full_labels, ":", fixed = TRUE)))
  variables <- trimws(variables)
  if(any(!variables %in% colnames(data))){
    .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' uses the '", budget$convention,
      "' convention with a kernel law, which needs the dimensions of all its terms; the data lack ",
      paste0("'", setdiff(variables, colnames(data)), "'", collapse = ", "), "."))
  }
  full_formula <- stats::reformulate(full_labels)
  model_frame <- stats::model.frame(full_formula, data = data, na.action = stats::na.pass)
  # the package resolves its own contrast functions (namespace-only callers)
  model_matrix <- .bt_model_matrix(model_frame, full_formula)
  assign <- attr(model_matrix, "assign")
  full_keys <- .bt_R2_term_key(attr(stats::terms(full_formula), "term.labels"))
  vapply(.bt_R2_term_key(declared), function(key){
    sum(assign == match(key, full_keys))
  }, numeric(1))
}

# Beta shapes and unit concentrations of a budget in a model with the
# retained units 'unit_alpha' (declared concentrations of the units).
.bt_R2_model_shapes <- function(budget, unit_alpha, reference_total){

  alpha_M <- if(identical(budget$concentration, "fixed_total")){
    unit_alpha * sum(budget$alpha) / sum(unit_alpha)
  }else{
    unit_alpha
  }
  A_M <- sum(alpha_M)
  shapes <- switch(
    budget$convention,
    common = c(budget$a, budget$b),
    linked = c(A_M, budget$b),
    conditioning = {
      A_ref <- if(identical(budget$concentration, "fixed_total")) sum(budget$alpha) else reference_total
      c(budget$a - (A_ref - A_M), budget$b + (A_ref - A_M))
    }
  )
  if(!all(is.finite(shapes)) || any(shapes <= 0)){
    .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' gives the Beta shapes (",
      format(shapes[1]), ", ", format(shapes[2]), ") in this model; both must be positive. ",
      "Increase the first shape of 'total' or change the 'convention'."))
  }
  list(a = shapes[[1]], b = shapes[[2]], alpha = alpha_M, A = A_M)
}

# Under the conditioning convention every open configuration of a gated
# budget needs a positive first shape: the smallest retained concentration
# is that of the ungated terms (always open) or of one gated term.
.bt_R2_check_gated_shapes <- function(budget, gates, budget_terms, unit_alpha, unit_term, unit_reference){

  if(!identical(budget$convention, "conditioning") || identical(budget$concentration, "fixed_total")){
    return(invisible(TRUE))
  }
  term_alpha <- vapply(budget_terms, function(term) sum(unit_alpha[unit_term == term]), numeric(1))
  ungated <- setdiff(budget_terms, gates$gated)
  A_min <- if(length(ungated) > 0L) sum(term_alpha[ungated]) else min(term_alpha[gates$gated])
  if(budget$a - (unit_reference - A_min) <= 0){
    .bt_R2_stop(paste0("The R-squared budget '", budget$name, "' gives a nonpositive first Beta shape (",
      format(budget$a - (unit_reference - A_min)), ") in its smallest open configuration; ",
      "increase the first shape of 'total' or change the 'convention'."))
  }

  invisible(TRUE)
}

# Centered empirical metric of a design block.
.bt_R2_metric <- function(X){

  X <- as.matrix(X)
  centered <- sweep(X, 2, colMeans(X), "-")
  crossprod(centered) / nrow(X)
}

# Binds the R2 term stubs of a formula once its design matrix is known:
# metric, coefficient maps, budget shapes, and the injected total and
# allocation priors. Called before the prior names get the parameter prefix.
.bt_R2_bind_formula <- function(prior_list, budgets, parameter, model_matrix,
                                terms_indexes, model_terms, formula, data){

  r2_terms <- names(prior_list)[vapply(prior_list, is.prior.R2, logical(1))]
  for(budget in budgets){

    budget_terms <- r2_terms[vapply(prior_list[r2_terms], function(prior) identical(prior$budget$name, budget$name), logical(1))]
    if(length(budget_terms) == 0L){
      next
    }
    if(any(c(budget$name, paste0(budget$name, "_allocation")) %in% names(prior_list))){
      .bt_R2_stop(paste0("The R-squared budget name '", budget$name, "' clashes with a formula term; ",
        "name the budget entry differently."))
    }
    # total concentration of the full candidate model
    if(identical(budget$law, "exact") || identical(budget$concentration, "fixed_total")){
      unit_reference <- sum(budget$alpha)
    }else if(budget$convention %in% c("conditioning", "linked")){
      dims <- .bt_R2_declared_dimensions(budget, formula, data)
      unit_reference <- sum(budget$alpha * dims)
    }else{
      unit_reference <- NA_real_
    }
    if(identical(budget$convention, "linked")){
      .bt_R2_check_linked_shape(budget, unit_reference)
    }

    maps <- list()
    unit_alpha <- numeric()
    unit_term <- character()
    for(term in budget_terms){
      columns <- which(terms_indexes == match(term, model_terms))
      X <- model_matrix[, columns, drop = FALSE]
      metric <- .bt_R2_metric(X)
      d <- ncol(X)
      label <- gsub("__xXx__", ":", term, fixed = TRUE)
      scales <- sqrt(diag(metric))
      if(any(!is.finite(scales)) || any(scales <= 0)){
        .bt_R2_stop(paste0("A design column of the R-squared budget term '", label,
          "' is constant after centering; its signal variance cannot carry a share of the budget."))
      }
      if(identical(budget$law, "exact")){
        # rank is judged on the correlation matrix, so that it does not depend
        # on the units of the predictors
        correlation <- metric / (scales %o% scales)
        R_correlation <- tryCatch(chol(correlation), error = function(e) NULL)
        if(is.null(R_correlation) || min(diag(R_correlation)) <= sqrt(.Machine$double.eps)){
          .bt_R2_stop(paste0("The design of the R-squared budget term '", label,
            "' is rank deficient after centering; the exact law needs linearly independent centered columns."))
        }
        map <- backsolve(sweep(R_correlation, 2, scales, "*"), diag(d))
        units <- 1L
      }else{
        map <- diag(1 / scales, nrow = d)
        units <- d
      }
      maps[[term]] <- list(map = map, metric = metric, d = d, units = units)
      unit_alpha <- c(unit_alpha, rep(prior_list[[term]]$alpha, units))
      unit_term <- c(unit_term, rep(term, units))
    }

    shapes <- .bt_R2_model_shapes(budget, unit_alpha, unit_reference)
    n_units <- length(unit_alpha)
    total_name <- paste0(parameter, "_", budget$name)
    allocation_name <- if(n_units > 1L) paste0(parameter, "_", budget$name, "_allocation") else NA_character_
    declared <- vapply(budget_terms, function(term) prior_list[[term]]$declared, character(1))
    gates <- .bt_R2_gate_structure(budget, budget_terms, declared)

    if(is.null(gates)){
      total_prior <- prior("beta", list(alpha = shapes$a, beta = shapes$b))
      attr(total_prior, "R2_budget") <- list(name = budget$name, role = "total", law = budget$law)
      prior_list[[budget$name]] <- total_prior
      if(n_units > 1L){
        allocation_prior <- .bt_R2_dirichlet_prior(shapes$alpha)
        attr(allocation_prior, "R2_budget") <- list(name = budget$name, role = "allocation", law = budget$law,
          units = .bt_R2_unit_labels(unit_term, maps))
        prior_list[[paste0(budget$name, "_allocation")]] <- allocation_prior
      }
      gating <- NULL
    }else{
      # a gated or spiked budget: one core entry emits the gates, the budget
      # indicator, the allocation and the total of the open configuration
      .bt_R2_check_gated_shapes(budget, gates, budget_terms, unit_alpha, unit_term, unit_reference)
      core <- .bt_R2_core_prior(budget, gates, parameter, budget_terms, unit_term, unit_alpha,
                                unit_reference, shapes)
      attr(core, "R2_budget")$units <- .bt_R2_unit_labels(unit_term, maps)
      prior_list[[budget$name]] <- core
      nodes <- .bt_R2_core_nodes(core)
      gating <- list(
        gate_names = nodes$gates,
        unit_gate  = unname(core$term_gate[unit_term]),
        indicator  = if(!is.null(gates$spike)) nodes$indicator
      )
    }

    for(term in budget_terms){
      this_prior <- prior_list[[term]]
      this_prior$bound <- TRUE
      if(is.prior.vector(this_prior)){
        this_prior$parameters <- list(K = maps[[term]]$d)
      }
      this_prior$parameter <- parameter
      this_prior$d <- maps[[term]]$d
      this_prior$map <- maps[[term]]$map
      this_prior$metric <- maps[[term]]$metric
      this_prior$units <- which(unit_term == term)
      this_prior$n_units <- n_units
      this_prior$shapes <- shapes
      this_prior$total_name <- total_name
      this_prior$allocation_name <- allocation_name
      this_prior$budget_terms <- budget_terms
      this_prior$gating <- gating
      # the core of a gated budget, for the marginal methods of the term
      this_prior$core <- if(!is.null(gating)) prior_list[[budget$name]]
      prior_list[[term]] <- this_prior
    }
  }
  prior_list
}

# The total and allocation priors that the coefficient node of a term reads.
.bt_R2_budget_prior_list <- function(prior){

  .bt_R2_require_bound(prior)
  out <- list()
  out[[prior$total_name]] <- prior("beta", list(alpha = prior$shapes$a, beta = prior$shapes$b))
  if(prior$n_units > 1L){
    out[[prior$allocation_name]] <- .bt_R2_dirichlet_prior(prior$shapes$alpha)
  }
  out
}

# The Dirichlet prior of the concentrations of the units of a model: the
# declared concentrations of R2_dirichlet() are at least 0.01, but the
# rescaling of fixed_total and the kernel units can give smaller (positive)
# concentrations, which the bound allocation keeps.
.bt_R2_dirichlet_prior <- function(alpha){

  out <- prior("dirichlet", list(alpha = pmax(alpha, 0.01)))
  out$parameters$alpha <- alpha
  out
}

.bt_R2_unit_labels <- function(unit_term, maps){

  labels <- gsub("__xXx__", ":", unit_term, fixed = TRUE)
  for(term in unique(unit_term)){
    if(maps[[term]]$units > 1L){
      labels[unit_term == term] <- paste0(labels[unit_term == term], "[", seq_len(maps[[term]]$units), "]")
    }
  }
  labels
}

# The Beta total and Dirichlet allocation priors that a budget adds to a
# formula prior list (they are parameters of the formula, not terms).
.bt_R2_is_budget_prior <- function(prior){
  !is.null(attr(prior, "R2_budget", exact = TRUE)) && !is.prior.point(prior)
}

.bt_R2_require_bound <- function(prior, parameter_name = NULL){

  if(!is.prior.R2(prior) || !isTRUE(prior$bound)){
    .bt_R2_unavailable(paste0(
      "The R-squared budget prior", if(!is.null(parameter_name)) paste0(" of '", parameter_name, "'"),
      " is unavailable before JAGS_formula() binds it to a formula design."),
      class = "BayesTools_R2_unbound")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Coordinates and JAGS syntax of a bound term

.bt_R2_coefficient_names <- function(prior, parameter_name){

  if(prior$d > 1L){
    paste0(parameter_name, "[", seq_len(prior$d), "]")
  }else{
    parameter_name
  }
}

.bt_R2_direction_coordinates <- function(prior, parameter_name){

  name <- .bt_R2_direction_name(parameter_name)
  if(prior$d > 1L){
    paste0(name, "[", seq_len(prior$d), "]")
  }else{
    name
  }
}

.bt_R2_allocation_eta_coordinates <- function(prior){

  paste0(.JAGS_prior_dirichlet_eta_name(prior$allocation_name), "[", seq_len(prior$n_units), "]")
}

.bt_dnode_R2_coefficient <- function(parameter_name, prior){

  .bt_R2_require_bound(prior, parameter_name)
  .bt_deterministic_node(
    "R2_coefficient",
    parameter_name,
    .bt_R2_coefficient_names(prior, parameter_name),
    dependencies = c(
      prior$total_name,
      if(prior$n_units > 1L) .bt_R2_allocation_eta_coordinates(prior),
      if(prior$n_units > 1L) c(prior$gating$gate_names, prior$gating$indicator),
      .bt_R2_direction_coordinates(prior, parameter_name)
    ),
    parameter = parameter_name,
    spec = list(
      law             = prior$law,
      d               = prior$d,
      map             = prior$map,
      units           = prior$units,
      n_units         = prior$n_units,
      total_name      = prior$total_name,
      allocation_name = prior$allocation_name,
      directions      = .bt_R2_direction_coordinates(prior, parameter_name),
      eta             = if(prior$n_units > 1L) .bt_R2_allocation_eta_coordinates(prior) else character(),
      gating          = prior$gating,
      coefficients    = .bt_R2_coefficient_names(prior, parameter_name)
    )
  )
}

.bt_dnode_R2_coefficient_emit <- function(node){

  spec <- node$spec
  share <- function(unit){
    if(spec$n_units > 1L){
      paste0(spec$total_name, " / (1 - ", spec$total_name, ") * ", spec$allocation_name, "[", unit, "]")
    }else{
      paste0(spec$total_name, " / (1 - ", spec$total_name, ")")
    }
  }
  z <- spec$directions
  if(identical(spec$law, "exact")){
    radius <- paste0("sqrt(", share(spec$units[1]), ")")
    if(spec$d == 1L){
      return(paste0(spec$coefficients, " <- ", radius, " * ", .JAGS_prior_number(spec$map[1, 1]),
        " * ", z, " / abs(", z, ")"))
    }
    norm <- paste0("sqrt(", paste0(z, "^2", collapse = " + "), ")")
    vapply(seq_len(spec$d), function(i){
      columns <- which(spec$map[i, ] != 0)
      combination <- paste0(vapply(columns, function(j){
        paste0(.JAGS_prior_number(spec$map[i, j]), " * ", z[j])
      }, character(1)), collapse = " + ")
      paste0(spec$coefficients[i], " <- ", radius, " * (", combination, ") / ", norm)
    }, character(1))
  }else{
    vapply(seq_len(spec$d), function(i){
      paste0(spec$coefficients[i], " <- sqrt(", share(spec$units[i]), ") * ",
        .JAGS_prior_number(spec$map[i, i]), " * ", z[i])
    }, character(1))
  }
}

.bt_dnode_R2_coefficient_evaluate <- function(node, lookup){

  spec <- node$spec
  total <- .bt_deterministic_lookup_value(lookup, spec$total_name)
  if(is.null(total)){
    return(NULL)
  }
  if(spec$n_units > 1L){
    # the allocation is normalized from its gamma coordinates (the declared
    # dependencies), so that a replay with changed coordinates never reads
    # stale normalized weights
    if(!all(spec$eta %in% colnames(lookup$draws))){
      return(NULL)
    }
    eta <- lookup$draws[, spec$eta, drop = FALSE]
    eta_sum <- rowSums(eta)
    if(any(!is.finite(eta)) || any(eta < 0) || any(!is.finite(eta_sum)) || any(eta_sum <= 0)){
      .bt_R2_stop_state("R-squared budget allocation gamma coordinates must be finite and nonnegative with positive row sums.")
    }
    if(!is.null(spec$gating)){
      # the allocation of a gated budget is redistributed over the open units
      open <- matrix(.bt_R2_unit_open_values(spec$gating, lookup, lookup$n), nrow = lookup$n)
      eta <- open * eta
      eta_sum <- rowSums(eta)
      eta_sum[eta_sum == 0] <- 1
    }
    allocation <- (eta / eta_sum)[, spec$units, drop = FALSE]
  }else{
    allocation <- matrix(1, nrow = lookup$n, ncol = 1L)
  }
  if(!all(spec$directions %in% colnames(lookup$draws))){
    return(NULL)
  }
  z <- lookup$draws[, spec$directions, drop = FALSE]
  if(any(!is.finite(total)) || any(total < 0 | total >= 1)){
    .bt_R2_stop_state("R-squared budget totals must lie in [0, 1).")
  }
  W <- total / (1 - total)
  if(identical(spec$law, "exact")){
    norm <- sqrt(rowSums(z^2))
    direction <- z / norm
    values <- sqrt(W * allocation[, 1]) * (direction %*% t(spec$map))
  }else{
    values <- sqrt(W * allocation) * sweep(z, 2, diag(spec$map), "*")
  }
  values <- matrix(values, nrow = lookup$n, ncol = spec$d)
  values
}

.JAGS_prior.R2 <- function(prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  directions <- .bt_R2_direction_coordinates(prior, parameter_name)
  kernel <- if(identical(prior$law, "laplace")){
    paste0("ddexp(0, ", .JAGS_prior_number(sqrt(2)), ")")
  }else{
    "dnorm(0, 1)"
  }
  paste0(
    paste0(directions, " ~ ", kernel, "\n", collapse = ""),
    paste0(.bt_deterministic_node_emit(.bt_dnode_R2_coefficient(parameter_name, prior)), "\n", collapse = "")
  )
}

.JAGS_init.R2 <- function(prior, parameter_name){

  .bt_R2_require_bound(prior, parameter_name)
  values <- .bt_R2_direction_rng(prior$law, prior$d)
  init <- list()
  init[[.bt_R2_direction_name(parameter_name)]] <- values
  init
}

.JAGS_monitor.R2 <- function(prior, parameter_name){

  c(parameter_name, .bt_R2_direction_name(parameter_name))
}

# Direction coordinates drawn from their prior: standard normal for the exact
# law and the normal kernel, unit-variance Laplace for the Laplace kernel.
.bt_R2_direction_rng <- function(law, n){

  if(identical(law, "laplace")){
    stats::rexp(n, rate = sqrt(2)) * sample(c(-1, 1), n, replace = TRUE)
  }else{
    stats::rnorm(n)
  }
}

# Log density of the direction coordinates of one draw.
.bt_R2_direction_lpdf <- function(law, z){

  if(identical(law, "laplace")){
    sum(log(sqrt(2) / 2) - sqrt(2) * abs(z))
  }else{
    sum(stats::dnorm(z, log = TRUE))
  }
}
