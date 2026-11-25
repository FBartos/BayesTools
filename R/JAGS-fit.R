#' @title Fits a 'JAGS' model
#'
#' @description Compiles JAGS model syntax, data, formulas and priors into
#' a model and samples it using \pkg{rjags}.
#' @param model_syntax jags syntax for the model part
#' @param data list containing data to fit the model (not including data for the formulas)
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names) of parameters not specified within the
#' \code{formula_list}. Factor priors must carry their levels, set with
#' [prior_factor_levels()].
#' @param formula_list named list of formulas to be added to the model
#' (names correspond to the parameter name created by each of the formula)
#' @param formula_data_list named list of data frames containing data for each formula
#' (names of the lists correspond to the parameter name created by each of the formula)
#' @param formula_prior_list named list of named lists of prior distributions
#' (names of the lists correspond to the parameter name created by each of the formula and
#' the names of the prior distribution correspond to the parameter names) of parameters specified
#' within the \code{formula}
#' @param formula_scale_list named list of named lists for standardizing continuous predictors
#' (names of the lists correspond to the parameter name created by each of the formula).
#' Each entry should be a named list where continuous predictors with \code{TRUE} values will
#' be standardized. Defaults to \code{NULL} (no standardization).
#' @param formula_random_prior_list optional named list of `prior_random()`
#' objects for random effects in `formula_list`. Required for any formula that
#' contains random effects.
#' @param formula_random_effects_compile_list optional named list of
#' `random_effects_compile()` objects controlling which formula random-effect
#' blocks are sampled and which are compiled as marginalized structural blocks.
#' Defaults to \code{NULL}, which preserves the current all-sampled behavior.
#' @param chains number of chains to be run, defaults to \code{4}
#' @param adapt number of samples used for adapting the MCMC chains, defaults to
#' \code{500}. In \code{JAGS_extend()}, the number of adaptation iterations run
#' when the backend builds the sampler from the stored chain states; the default
#' \code{NULL} uses the setting of the initial fit. A value supplied to
#' \code{JAGS_extend()} applies to that call only: the returned fit keeps the
#' initial setting, which later calls without \code{adapt} use again.
#' @param burnin number of burnin iterations of the MCMC chains, defaults to
#' \code{1000}. Its last iterations are checked for agreement of the chains
#' before sampling starts, and with \code{autofit} burn-in can be extended (see
#' \code{burnin_max_Rhat} in \code{autofit_control}); \code{fit$burnin} holds
#' the burn-in that was run.
#' @param sample number of draws to retain per chain, defaults to \code{4000}.
#' Each chain runs \code{sample * thin} sampling iterations.
#' @param thin thinning interval for the MCMC samples, defaults to \code{1}
#' @param autofit whether the models should be refitted until convergence criteria
#' specified in \code{autofit_control}. Defaults to \code{FALSE}.
#' @param autofit_control a list of arguments controlling the autofit function.
#' Convergence thresholds and \code{max_time$time} must be finite and non-missing.
#' Set a convergence threshold or \code{max_time} to \code{NULL} to disable that check.
#' Possible options are:
#' \describe{
#'   \item{max_Rhat}{maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}. R-hat is the larger of the point estimate and
#'   the upper confidence limit of the potential scale reduction factor,
#'   whereas the summary tables show the point estimate (see
#'   [JAGS_check_convergence()]).}
#'   \item{min_ESS}{minimum effective sample size. Defaults to \code{500}.}
#'   \item{max_error}{maximum MCMC error. Defaults to \code{0.01}.}
#'   \item{max_SD_error}{maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.}
#'   \item{max_time}{list with fields \code{time} and \code{unit}, specifying
#'   when the automatic fitting function is stopped. The \code{unit} value
#'   corresponds to the \code{units} argument of \link[base]{difftime}.}
#'   \item{max_extend}{number of times after which the automatic fitting function is stopped.}
#'   \item{sample_extend}{number of samples between each convergence check. Defaults to
#'   \code{1000}. With \code{"auto"}, each extension is sized from the
#'   convergence check that failed, which then assesses every selected
#'   parameter: with \code{n} retained draws per chain and \code{f} the largest
#'   required growth over the parameters (\code{min_ESS / ESS},
#'   \code{(error / max_error)^2}, and \code{(SD_error / max_SD_error)^2} for
#'   the criteria that are set; at least \code{2} for an R-hat above
#'   \code{max_Rhat}, and for a parameter whose diagnostics are not assessable
#'   unless \code{allow_not_assessable}), the
#'   extension is \code{n * (f - 1) * 1.1} draws per chain rounded up to a
#'   multiple of 100, then at most \code{2 * n} rounded down to a multiple of
#'   100, then at least 1000 (the minimum wins when \code{2 * n} is smaller),
#'   and at most the largest multiple of 100 whose sampling iterations
#'   (draws times \code{thin}) stay within the backend's integer limit.
#'   With \code{max_time}, it is then at most the draws per chain that the
#'   time left allows at the rate observed so far (the time on the
#'   \code{max_time} clock divided by the retained draws per chain run since
#'   that clock started; in \code{JAGS_fit()} the clock starts before the
#'   initial run, so the rate includes its adaptation and burn-in), rounded
#'   down to a multiple of 100; when that is below 100, the fitting stops with
#'   the \code{max_time} warning. In \code{JAGS_extend()}, the first extension
#'   is sized from the check of the supplied fit and, with no rate observed
#'   yet, is not shortened. Without \code{max_time}, each extension of a fit
#'   that cannot converge multiplies the draws per chain by about 2.1 (a
#'   failed R-hat or a parameter whose diagnostics are not assessable) and up
#'   to 3 (large ESS or error deficits), up to \code{max_extend}. An
#'   extension, of either kind, whose last draw would have an iteration
#'   number beyond \code{.Machine$integer.max} is not started: the fitting
#'   stops with the last fit and a warning of its own.}
#'   \item{restarts}{maximum number of fitting attempts when the model fails
#'   to initialize; every failed attempt but the last is retried with new
#'   initial values (\code{restarts = 1} fits once, without a retry). A
#'   failed attempt that is retried is reported as a warning of the fit.
#'   Defaults to \code{10}.}
#'   \item{check_indicators}{whether model indicator variables should be included
#'   in convergence checks, including when \code{monitor} selects a narrower
#'   parameter set. Auxiliary inclusion-probability coordinates remain excluded
#'   unless requested explicitly. Defaults to \code{FALSE}.}
#'   \item{monitor}{optional character vector selecting parameters for
#'   convergence checks. Base names select all indexed elements. Requests are
#'   resolved against the monitored nodes that the summary tables show,
#'   including the derived and auxiliary ones that the default selection
#'   excludes. Base names that are not
#'   monitored are rejected before sampling (error class
#'   \code{BayesTools_parameter_not_found}); an indexed element such as
#'   \code{"mu[2]"} is resolved against the fitted columns after the first
#'   run and fails then when it does not exist. Defaults to \code{NULL}, which
#'   checks the default selection of [JAGS_check_convergence()]: the sampled
#'   parameters of the priors, \code{add_parameters} that depend on a
#'   stochastic node, and the stochastic nodes that formulas monitor (such as
#'   latent random effects and LKJ primitives). Deterministic nodes that
#'   formulas monitor (such as random-effect correlation matrices and their
#'   Cholesky factors) are derived and not checked by default. The role of
#'   every fitted coordinate is stored in the \code{convergence_role} column
#'   of [parameter_coordinates()].}
#'   \item{allow_not_assessable}{whether undefined diagnostics for requested
#'   sampled parameters may be ignored. Defaults to \code{FALSE}.}
#'   \item{burnin_max_Rhat}{the threshold of the readiness check between
#'   burn-in and sampling, which runs in every fit, also without
#'   \code{autofit}. Defaults to \code{1.5}; an omitted element takes its
#'   default and \code{NULL} turns the check off. The last
#'   \code{burnin_extend} iterations of burn-in (all of it when shorter) run
#'   with the monitors on; their draws are checked and then dropped. Each
#'   chain's block is split into halves, and the chains are ready when the
#'   point estimate of the potential scale reduction factor over the
#'   half-chains (\code{coda::gelman.diag()}) is at most \code{burnin_max_Rhat}
#'   for every quantity that the convergence check assesses (the selection of
#'   \code{monitor} and \code{check_indicators}; a requested state of an
#'   indicator, such as \code{"mu_indicator (state 3)"}, that the block does
#'   not visit is not assessable there). A quantity that is constant over the
#'   half-chains, or not finite in each of them, is not assessable; one that
#'   is not finite in some half-chains only fails the check. The check needs at least two
#'   chains and a block of at least 100 iterations. It is skipped in a model
#'   to which JAGS assigns no sampler (every node is drawn directly given its
#'   parents, as in a model without data): its chains cannot disagree, and its
#'   whole burn-in runs unmonitored, as with the check off. When the chains are not
#'   ready, a fit with \code{autofit} extends burn-in by blocks of
#'   \code{burnin_extend} iterations, each checked alone, until they are
#'   ready, \code{burnin_max_extend} blocks were added, or \code{max_time} is
#'   exceeded (checked before each block); a fit without \code{autofit} keeps
#'   its burn-in. Chains that are still not ready start sampling; a warning
#'   of class \code{BayesTools_burnin_not_converged} (also
#'   \code{BayesTools_fit_diagnostic}) is issued once that initial sampling
#'   of the attempt the fit keeps has finished (a restart attempt that is
#'   discarded warns nothing), before any automatic extension, whatever
#'   \code{silent} is, and kept in the fit's \code{warnings} attribute. When no monitored node lacks
#'   observed descendants, a fit whose chains are ready at once has the same
#'   draws, final states and samplers as with the check off: monitoring changes
#'   no sampler, but JAGS samples a monitored node without observed
#'   descendants (such as a predictive quantity) from the chain's random-number
#'   stream. The check cannot detect chains that are all in the same wrong
#'   state.}
#'   \item{burnin_extend}{the number of final burn-in iterations that the
#'   readiness check assesses and of each extension of burn-in; at least
#'   \code{100}. Defaults to \code{500}.}
#'   \item{burnin_max_extend}{the largest number of extensions of burn-in, with
#'   \code{autofit} only. Defaults to \code{10}. An extension after which the
#'   last retained draw would have an iteration number beyond
#'   \code{.Machine$integer.max} is not run; the warning then names that
#'   limit.}
#' }
#' @param parallel whether to distribute chains over PSOCK worker processes.
#' Defaults to \code{FALSE}, which runs all chains in the calling process.
#' Each sampling process uses one JAGS thread.
#' @param cores maximum number of workers when \code{parallel = TRUE},
#' defaults to \code{chains}.
#' @param silent whether the function should proceed silently, defaults to \code{TRUE}
#' @param seed seed for random number generation by \code{JAGS_fit()}: an
#'   integer between \code{-.Machine$integer.max} and \code{.Machine$integer.max}.
#'   Draws are reproducible for the same \code{seed}, \code{chains}, and R
#'   random-number generator settings (\code{RNGkind()}), and, with
#'   \code{parallel = TRUE}, the same \code{cores} and backend versions. Initial values are
#'   drawn from the priors after \code{set.seed(seed)}, and each chain's JAGS
#'   random-number seed is drawn from the same seed with
#'   \code{sample.int(.Machine$integer.max, chains)} (see
#'   \code{\link{JAGS_get_inits}}). A chain's JAGS seed therefore does not
#'   depend on the number of chains, but the JAGS generator does: with more
#'   than four chains it changes from \code{"base::Super-Duper"} to
#'   \code{"lecuyer::RngStream"}, so the draws of a chain differ between fits
#'   with up to four and with more chains. The chains of adjacent seeds (e.g.,
#'   \code{seed} and \code{seed + 1}) do not share random-number streams.
#'   Automatic restarts after a failed initialization (see \code{restarts} in
#'   \code{autofit_control}) draw their seeds from a separate stream of the
#'   same seed (the first L'Ecuyer-CMRG substream), so the i-th restart does
#'   not repeat the fit with \code{seed + i}. A seeded fit restores the
#'   caller's random-number state (\code{.Random.seed} and \code{RNGkind()})
#'   when it returns. With \code{seed = NULL}, the initial values are drawn
#'   from one seed taken from the caller's random-number stream (see
#'   \code{\link{JAGS_get_inits}}), which the fit then continues to use.
#' @param worker_output optional file path for parallel worker stdout and stderr.
#'   The parent directory must exist. Workers append to the same file, so messages
#'   can interleave. \code{NULL} retains the backend default of discarding worker
#'   output. This setting is call-specific and is not retained with the fit;
#'   supply it again to \code{JAGS_extend()} when needed. Connection failures stop
#'   automatic fitting retries and preserve the original backend error.
#' @param sampler_gates whether the branch parameters of a publication-bias
#'   mixture prior (a [prior_mixture()] of [prior_none()], [prior_PET()],
#'   [prior_PEESE()], [prior_weightfunction()], [prior_phacking()], and
#'   [prior_bias()] branches) are declared to the gated sampler of the
#'   BayesTools JAGS module: while its branch is not selected, a parameter is
#'   drawn directly from its prior instead of by slice sampling. A branch whose
#'   nodes the model uses outside the mixture composition (in
#'   \code{model_syntax}, another prior, or a formula) is not gated. Defaults to
#'   \code{FALSE}.
#' @param sampler_blocks optional list of character vectors of names of the
#'   combined prior list (\code{prior_list} and the formula priors). The scalar
#'   parameters of the priors of one vector are updated jointly by the factor
#'   slice sampler of the BayesTools JAGS module, which slice-samples along the
#'   principal directions of their posterior covariance learned during
#'   adaptation. A block may contain simple continuous scalar priors, treatment
#'   and independent factor priors, a single weight function, PET and PEESE
#'   priors, direct random-effect SDs and scalar correlations with simple
#'   continuous scalar priors, and the allocation SD and Dirichlet allocation
#'   weights of a random-effect variance allocation; a prior may be in one block
#'   only. A
#'   block that expands to fewer than two parameters is not declared and is
#'   recorded as declined. \code{sampler_blocks} may also be a function of one
#'   argument, the combined prior list (formula priors under their final names
#'   and with their attributes), called once by \code{JAGS_fit()}, that returns
#'   \code{NULL} or such a list. Defaults to \code{NULL}.
#'   Neither \code{sampler_gates} nor \code{sampler_blocks} changes the model or
#'   its posterior distribution, only the draws. Invalid declarations stop with
#'   an error of class \code{BayesTools_sampler_declaration_error}. The
#'   declarations are stored in \code{attr(fit, "sampler_declarations")}, which
#'   \code{JAGS_extend()} keeps: a list of the data frames \code{gates},
#'   \code{blocks} and \code{declined} (one row per declined gate branch or
#'   block, with its \code{kind}, \code{prior}, \code{branch} or \code{block},
#'   \code{nodes} and \code{reason}).
#'   For successful fits with gate declarations, \code{gate_claims} records the
#'   sampler assignment for each gate when the reported JAGS node can be matched.
#'   A gate whose assignment status or sampler differs across workers has
#'   status \code{"mixed"}.
#'   A density that cannot be calculated (NaN)
#'   at a point the factor slice sampler evaluates stops the fit with a JAGS
#'   error, as with the stock samplers.
#' @param sampler_joint optional package callback accepting `prior_list`,
#'   `formula_design` and `data`. It receives the combined bound prior list,
#'   the finalized formula designs and the combined model and formula data,
#'   once before compilation. Return `NULL` to decline, or a list with a
#'   nonempty `markers` list and optional `data` list of named finite numeric
#'   constants that are referenced by the markers. Each marker specifies a fresh
#'   `name`, an ordered `arguments` list and one of `tag` or `distribution`.
#'   Each argument contains one or more node references
#'   (bare names or positive integer subscripts or ranges) or finite numeric
#'   literals. Multiple entries in one argument are combined with `c()`.
#'   References can come from the existing model and data or the callback's
#'   numeric constants. Marker names are reserved to the declarations.
#'   A `tag` is a registered deterministic JAGS function with a name prefixed
#'   with `<lowercased-module>_` for a provider listed in `jags_modules`. It must
#'   return a fixed finite value for admissible arguments. Use a dedicated
#'   function without aliases. The model must reserve that function for tags
#'   to prevent JAGS from reusing a scientific node. The tag creates an
#'   unmonitored terminal node with no observed likelihood contribution.
#'   The BayesTools finite readout functions are reserved for `sampler_probabilities`.
#'   A `distribution` must be a normalized point mass at zero for every
#'   admissible parent value. Its marker is observed at zero.
#'   Use a distribution name prefixed with
#'   `d<lowercased-module>_` for a provider listed in `jags_modules`, or use
#'   the BayesTools block marker `dbt_block`. Declare BayesTools blocks through
#'   one of `sampler_joint` or `sampler_blocks`. Use `sampler_gates` for model-verified
#'   gates, which take precedence over `dbt_block` when both mark a node.
#'   Supply package and module requirements through `required_packages`
#'   and `jags_modules`. The
#'   validated declaration is stored in the `joint` entry of
#'   `attr(fit, "sampler_declarations")`. Extensions reuse the compiled model
#'   text and data.
#' @param sampler_probabilities optional character vector of sampled finite
#'   indicator node names or scalar indexed references in the compiled model.
#'   `fit$sampler_probabilities` gives one entry per requested indicator.
#'   The recorder can replace an active JAGS 5 `base::Finite` update for a
#'   full-rank scalar `dbern` or `dcat` with fixed integer support starting at
#'   0 or 1, respectively. It preserves competing sampler assignments. Its
#'   entry is an `mcmc.list` with one row per retained draw in the same chain
#'   and iteration order as the model draws. The columns `p[k]` and
#'   `log_p[k]` are conditional probabilities and log probabilities at the
#'   indicator update in that iteration for support label `k`. The columns
#'   `selected` and `lower` give the sampled label and the lowest support
#'   label. Compiled aliases of one node return the same trace. An entry is
#'   `NULL` when no recorder was assigned to that node. A probability trace
#'   requires the same availability and support labels
#'   on every worker.
#'   Extensions append the corresponding rows. By default no probabilities
#'   are recorded. An incomplete or incompatible recording raises a
#'   `BayesTools_sampler_probabilities_error` backend condition.
#' @param sampler_nested_weights logical. Share physical weight coordinates for eligible
#'   nested cumulative one-sided weights in a publication-bias mixture.
#'   The scientific priors and public weights are preserved. Unsupported pairs
#'   and mixtures whose private weight nodes are explicitly used or monitored
#'   keep their ordinary parameterization. The `nested_weights` entry in
#'   `attr(fit, "sampler_declarations")` records the two branch indices when
#'   sharing was used or the reason it was declined. Defaults to `FALSE`.
#' @param add_parameters vector of additional parameter names that should be
#' monitored but were not specified in the \code{prior_list}. Automatic fitting
#' checks their convergence like that of the other parameters. Fully observed
#' data and a node that the model syntax defines deterministically (\code{<-} or
#' \code{=}) from data and constants only are structural constants and do not
#' block convergence. A
#' node that is stochastic (\code{~}) or depends on a stochastic node is
#' checked, and its constant draws are not assessable (see
#' \code{autofit_control$allow_not_assessable}). A name that repeats a parameter
#' of \code{prior_list} is refused. A parameter that a formula needs monitored,
#' such as the scalar SD source of a random-effect block
#' (\code{random_sd_source("tau")}), is added automatically and, when
#' \code{prior_list} already holds its prior, is monitored through that prior.
#' @param required_packages character vector specifying list of packages containing
#' JAGS models required for sampling (in case that the function is run in parallel or in
#' detached R session). Defaults to \code{NULL}. Parallel workers must load
#' the same versions and native library builds as the parent session. When
#' developing packages with \code{load_all()}, install matching builds into
#' a library and set \code{R_LIBS_USER} before starting R and its workers.
#' @param jags_modules character vector specifying JAGS modules required by the
#' generated model syntax. Defaults to \code{NULL}.
#' @param runtime_setup optional function accepting one \code{context} list. After
#' required packages and JAGS modules load, it is called with \code{phase = "start"}
#' in the calling process, then each parallel worker. The context contains
#' \code{role} (\code{"local"}, \code{"coordinator"}, or \code{"worker"}),
#' \code{chains} (total chains), \code{processes} (sampling processes),
#' \code{process_id} (one-based worker index, or zero for the coordinator),
#' \code{process_chains} (chains assigned to this process; zero for the
#' coordinator), and \code{parallel}. Chains are divided round-robin among at
#' most \code{min(cores, chains)} workers. On exit, including fitting failures,
#' workers are stopped before the calling process receives \code{phase = "finish"}
#' with its original topology context. Worker processes do not receive a finish
#' callback. If stopping workers fails, a warning is issued and the finish
#' callback is not run; cleanup still attempts every remaining worker and
#' closes failed connections. This permits releasing coordinator resources before
#' parallel work and restoring them afterward. Return values are ignored. The callback must be
#' idempotent, must not change random-number generator state, and should capture
#' only small immutable settings. Packages it uses must be included in
#' \code{required_packages}. The callback is stored with the fit;
#' \code{JAGS_extend()} reuses it with the current chain and worker topology.
#' Supplying \code{NULL} explicitly disables replay.
#' @param runtime_cache optional function accepting \code{context} and
#' \code{state = NULL} for optional computational caches. It must not change
#' inference targets or random-number state and must capture only small settings.
#' After runtime setup, a fit or extension calls it with \code{context$phase = "restore"}
#' and a list of prior process shards assigned to this process (empty for new
#' fits or workers without an assigned shard). Each prior shard
#' is sent to at most one current process, in round-robin order; shards are not
#' duplicated when worker count increases. Before workers stop, a completed fit
#' calls it with \code{context$phase = "capture"} and \code{state = NULL}.
#' Capture returns one serializable shard, or \code{NULL} to retain no state.
#' The fitted object's \code{runtime_state} attribute stores process shards
#' separately from its small \code{runtime_cache} callback. Backend extension
#' calls receive neither attribute. \code{JAGS_extend()} restores old state
#' independently of whether its callback elects to capture fresh state.
#' Callback warnings and errors are reported without discarding valid draws;
#' a failed capture callback leaves a \code{NULL} shard for its process. If a
#' worker connection fails during capture, no shard is retained from any
#' process for that fit. Capture warnings, including such failures, are also
#' recorded in the fit's \code{warnings} attribute. Error and interrupt cleanup does
#' not attempt capture; ordinary runtime cleanup still runs. Retained R snapshots,
#' serialization buffers, and the caller's old fitted object consume memory in
#' addition to active cache storage. Callback packages must be listed in
#' \code{required_packages}. Passing \code{NULL} disables capture and restore.
#' @param runtime_report optional function accepting one \code{context} list,
#' as for \code{runtime_setup} but with \code{context$phase = "report"}, that
#' reads small values of a process, such as sampler counters. It is called in
#' every process that ran chains (the calling process of a serial fit, each
#' worker of a parallel fit) after every backend run: the initial run, each
#' restart attempt, each automatic extension, and each extension of
#' \code{JAGS_extend()}. JASP sampling batches share one report for the initial run.
#' An error in a JASP progress callback propagates to the caller. If warmup
#' has completed, the reporting callback runs before the error propagates.
#' The callback must not change
#' inference targets or random-number state, must capture only small settings,
#' and returns \code{NULL} or a small serializable value. The values of runs
#' that do not give the returned fit (failed attempts and failed extensions)
#' are dropped; the callback is still called after them, so that
#' read-and-clear state starts the next run clean, but not after a worker
#' connection failure. The \code{runtime_report_values} attribute of a fit is
#' the list of the kept runs of the call, in order, each a list with one
#' element per process (\code{NULL} elements kept in place); it is \code{NULL}
#' when no process reported a value, and a fit that failed has none. Callback
#' warnings and errors are signalled with the process they came from and
#' recorded in the fit's \code{warnings} attribute when the run is kept. A
#' callback error leaves \code{NULL} for its process; after a warning, the
#' value is kept. A failed transport of the report leaves \code{NULL} for every
#' process, with one warning. The callback is stored as the
#' \code{runtime_report} attribute; \code{JAGS_extend()} reuses it and reports
#' only its own runs (the values of the extended fit are not carried over).
#' Callback packages must be listed in \code{required_packages}. Passing
#' \code{NULL} disables the report.
#' @param fit a 'BayesTools_fit' object (created by \code{JAGS_fit()} function) to be
#' extended. A derived-draw view signals
#' \code{BayesTools_draws_view_sampling_unavailable}, also
#' \code{BayesTools_sampling_unavailable} and
#' \code{BayesTools_draws_view_unavailable}, before runtime arguments are forced;
#' extend its \code{original_fit} and regenerate it with [JAGS_with_draws()].
#' Fits without the parameter map and fit contract of this version
#' of BayesTools (such as fits created by BayesTools 0.3.0) must be refitted:
#' they stop with an error of class \code{BayesTools_refit_required} (see
#' [JAGS_validate_fit_contract()]).
#' @param ... additional hidden arguments
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#' }
#'
#' @return A successful fit is a \code{BayesTools_fit}. It contains the sampled
#' chains, model and data used for compilation, monitored nodes, and final
#' per-chain parameter and RNG states. Extract the chains with
#' \code{coda::as.mcmc.list(fit)}. The sampler table \code{fit$samplers} is from
#' the most recent successful backend run. For parallel fits it is the first
#' worker's table. If worker assignments differ, \code{fit$samplers_by_worker}
#' gives each worker's \code{chains} and \code{samplers}. A failed initial fit returns its error
#' condition.
#'
#' \code{JAGS_extend()} rebuilds the sampler from the stored states, adapts
#' it for \code{adapt} iterations when needed, then appends draws without
#' further burn-in. A saved and reloaded fit extends from the same state.
#' Within one call, automatic extensions continue the serial sampler. JASP fits
#' also keep parallel worker samplers through progress updates and automatic
#' extensions. Other parallel fits rebuild at each extension. The returned fit stores
#' ordinary R data and can be serialized with \code{saveRDS()}.
#'
#' The \code{adaptation} attribute of a fit records the adaptation settings:
#' \code{initial}, the \code{adapt} of \code{JAGS_fit()}, and
#' \code{extensions}, a data frame with one row
#' per backend extension (automatic extensions of \code{JAGS_fit()} and the
#' extensions of \code{JAGS_extend()}) giving the \code{adapt} passed to the
#' backend and whether the backend \code{rebuilt} the sampler, the only case in
#' which it adapts. Warnings of the backend that adaptation was incomplete are
#' issued and also kept in the fit's \code{warnings} attribute.
#'
#' The \code{burnin_check} attribute of a fit records the readiness check
#' between burn-in and sampling (see \code{burnin_max_Rhat} in
#' \code{autofit_control}); it is \code{NULL} when the check did not run. It
#' runs in the initial fit of \code{JAGS_fit()} only, afresh in each restart
#' attempt, and not in \code{JAGS_extend()}, the automatic extensions after
#' sampling, or a model without samplers. The record holds \code{burnin_max_Rhat}, \code{block} (the
#' iterations of the checked end of the requested burn-in),
#' \code{burnin_extend}, \code{blocks} (one row per checked block, with the
#' iterations \code{sweeps} it ran and its largest value \code{Rhat} and
#' \code{parameter}), \code{extension} (the iterations added to burn-in),
#' \code{status} (\code{"ready"}, \code{"not_ready"}, or
#' \code{"not_assessable"} when no quantity could be assessed, which makes no
#' claim of readiness and extends nothing), \code{stop_reason}
#' (\code{"ready"}, \code{"max_extend"}, \code{"max_time"},
#' \code{"iteration_limit"}, \code{"no_autofit"}, or \code{"not_assessable"}),
#' and \code{failing} (the
#' quantities of the last block above the threshold, with their values).
#' \code{fit$burnin} is the burn-in that was run, and the iteration numbers
#' of the draws and of \code{sampler_probabilities} count it.
#'
#' @seealso [JAGS_check_convergence()]
#'
#' @export JAGS_fit
#' @export JAGS_extend
#' @name JAGS_fit
NULL

# The adaptation record of a fit: the initial setting and, per backend
# extension, the requested adaptation length and whether it rebuilt the sampler.
.bt_adaptation_record <- function(initial){

  list(
    initial    = as.integer(initial),
    extensions = data.frame(adapt = integer(), rebuilt = logical())
  )
}

.bt_adaptation_add_extension <- function(record, adapt, rebuilt){

  record$extensions <- rbind(
    record$extensions,
    data.frame(adapt = as.integer(adapt), rebuilt = rebuilt)
  )
  record
}

# The adaptation policy belongs to the stored fit.
.bt_fit_adaptation_record <- function(fit){

  record <- attr(fit, "adaptation", exact = TRUE)
  valid <- is.list(record) && identical(names(record), c("initial", "extensions")) &&
    is.integer(record$initial) && length(record$initial) == 1L &&
    !is.na(record$initial) && record$initial >= 0L &&
    is.data.frame(record$extensions) &&
    identical(names(record$extensions), c("adapt", "rebuilt"))
  if(!valid){
    .bt_stop_refit_required(
      "The fitted object's adaptation record is malformed. Refit the model with this version of BayesTools."
    )
  }
  record
}

# Record warnings issued by the sampling process alongside the draws.
.bt_collect_adaptation_warnings <- function(expr){

  messages <- character()
  value <- withCallingHandlers(expr, warning = function(condition){
    if(inherits(condition, "BayesTools_adaptation_incomplete")){
      messages <<- c(messages, conditionMessage(condition))
    }
  })
  list(value = value, warnings = messages)
}

# Whether 'sample' more draws per chain keep the public iteration number of
# the fit's last draw within the integer range of the draw geometry (the
# chains' coda timing). An extension that would leave it is not started.
# Chains without timing have no iteration numbers to keep in range.
.bt_extension_fits_iterations <- function(fit, sample){

  timing <- lapply(fit[["mcmc"]], attr, which = "mcpar", exact = TRUE)
  if(length(timing) == 0L || any(lengths(timing) != 3L)){
    return(TRUE)
  }
  last <- max(vapply(timing, function(mcpar) as.double(mcpar[[2L]]), numeric(1)))
  last + as.double(sample) * timing[[1L]][[3L]] <= .Machine$integer.max
}

# The warning of automatic fitting stopped by the integer range of the
# iteration numbers.
.bt_iteration_limit_message <- function(sample){

  paste0("The automatic model fitting was terminated due to the iteration limit: an extension of ",
    sample, " draws per chain would number the draws beyond the largest iteration number (",
    .Machine$integer.max, "). Refit with a smaller 'adapt', 'burnin', or 'sample' times 'thin' ",
    "to leave room for extensions.")
}

# A JASP fit keeps its compiled samplers in its private workers through autofit.
.bt_run_extension <- function(fit, sample, adapt, cl, silent){

  rebuilt <- if(is.null(cl)) is.null(fit[["live_model"]]) else !isTRUE(fit[["live_model"]])
  run <- .bt_collect_adaptation_warnings(.JAGS_run_backend(
    .bt_rjags_extend,
    list(fit = fit, sample = sample, adapt = as.integer(adapt), cl = cl, silent = silent)
  ))
  list(
    extension = run$value,
    rebuilt   = rebuilt,
    warnings  = run$warnings
  )
}

# The runtime report after one backend run ('result': the run's fit or error).
# The callback runs after every attempt, so that a consumer with read-and-clear
# state starts the next attempt clean, but not when the workers are gone; the
# values and warnings of a run are retained only when the run succeeded.
.bt_runtime_report_after_run <- function(reports, runtime_report, result,
                                         chains, cl = NULL){

  if(is.null(runtime_report) ||
     inherits(result, "BayesTools_JAGS_worker_connection_error")){
    return(reports)
  }
  report <- .JAGS_capture_runtime_report(runtime_report, chains, cl)
  if(inherits(result, "error")){
    return(reports)
  }
  list(
    runs     = c(reports$runs, list(report$values)),
    warnings = c(reports$warnings, report$warnings)
  )
}

.bt_runtime_reports_empty <- function(){

  list(runs = list(), warnings = character())
}

.bt_append_fit_warnings <- function(fit, messages){

  if(length(messages) == 0L){
    return(fit)
  }
  attr(fit, "warnings") <- unique(c(
    attr(fit, "warnings", exact = TRUE),
    messages
  ))
  fit
}

#' @rdname JAGS_fit
JAGS_fit <- function(model_syntax, data = NULL, prior_list = NULL, formula_list = NULL, formula_data_list = NULL, formula_prior_list = NULL, formula_scale_list = NULL,
                     chains = 4, adapt = 500, burnin = 1000, sample = 4000, thin = 1,
                     autofit = FALSE, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10, check_indicators = FALSE, monitor = NULL, allow_not_assessable = FALSE, burnin_max_Rhat = 1.5, burnin_extend = 500, burnin_max_extend = 10),
                     parallel = FALSE, cores = chains, silent = TRUE, seed = NULL,
                     add_parameters = NULL, required_packages = NULL,
                     formula_random_prior_list = NULL, formula_random_effects_compile_list = NULL,
                     jags_modules = NULL, runtime_setup = NULL, runtime_cache = NULL, worker_output = NULL,
                     sampler_gates = FALSE, sampler_blocks = NULL, runtime_report = NULL,
                     sampler_joint = NULL, sampler_probabilities = NULL, sampler_nested_weights = FALSE, ...){

  # A seeded fit leaves the caller's random-number state as it found it.
  if(!is.null(seed)){
    rng_state <- .bt_rng_state()
    on.exit(.bt_rng_restore(rng_state), add = TRUE)
  }
  dots <- list(...)
  worker_output <- .JAGS_validate_worker_output(worker_output)

  ### check input
  model_syntax <- .check_JAGS_syntax(model_syntax)
  JAGS_check_and_list_fit_settings(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed)
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)
  # The readiness check between burn-in and sampling looks at the last
  # 'burnin_extend' sweeps of burn-in (all of it when shorter), if they are at
  # least 100 and there are at least two chains.
  burnin_settings <- .bt_burnin_check_settings(autofit_control)
  readiness_block <- min(burnin_settings[["burnin_extend"]], burnin)
  check_readiness <- !is.null(burnin_settings[["burnin_max_Rhat"]]) && chains >= 2 && readiness_block >= 100
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(required_packages, "required_packages", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(jags_modules, "jags_modules", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  .JAGS_validate_runtime_setup(runtime_setup)
  .JAGS_validate_runtime_cache(runtime_cache)
  .JAGS_validate_runtime_report(runtime_report)
  .bt_check_sampler_declaration_arguments(sampler_gates, sampler_blocks)
  sampler_probabilities <- .bt_sampler_check_plan_arguments(sampler_joint, sampler_probabilities)
  check_bool(sampler_nested_weights, "sampler_nested_weights", allow_NA = FALSE)
  check_list(formula_list, "formula_list", allow_NULL = TRUE)
  check_list(formula_data_list, "formula_data_list", allow_NULL = is.null(formula_list))
  check_list(formula_prior_list, "formula_prior_list", allow_NULL = is.null(formula_list))
  check_list(formula_random_prior_list, "formula_random_prior_list", allow_NULL = TRUE)
  check_list(formula_random_effects_compile_list, "formula_random_effects_compile_list", allow_NULL = TRUE)
  check_list(formula_scale_list, "formula_scale_list", allow_NULL = TRUE)
  .bt_validate_jags_formula_lists(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list,
    formula_scale_list = formula_scale_list
  )
  if(!is.null(formula_random_prior_list)){
    for(parameter in names(formula_random_prior_list)){
      .bt_check_prior_random(formula_random_prior_list[[parameter]])
    }
  }
  if(!is.null(formula_random_effects_compile_list)){
    for(parameter in names(formula_random_effects_compile_list)){
      .bt_check_random_effects_compile(formula_random_effects_compile_list[[parameter]])
    }
  }
  if(!is.null(formula_list)){
    for(parameter in names(formula_list)){
      if((is.language(formula_list[[parameter]]) ||
          inherits(formula_list[[parameter]], "BayesTools_random_effects")) &&
         .has_random_effects(formula_list[[parameter]]) &&
         (is.null(formula_random_prior_list) || is.null(formula_random_prior_list[[parameter]]))){
        stop(
          "JAGS_fit() requires 'formula_random_prior_list' with a prior_random() object for formula random effects in parameter '",
          parameter,
          "'.",
          call. = FALSE
        )
      }
    }
  }

  ### add formulas
  if(!is.null(formula_list)){

    # obtain settings for each formula
    formula_output <- list()
    for(parameter in names(formula_list)){
      formula_output[[parameter]] <- JAGS_formula(
        formula        = formula_list[[parameter]],
        parameter      = parameter,
        data           = formula_data_list[[parameter]],
        prior_list     = formula_prior_list[[parameter]],
        formula_scale  = if(!is.null(formula_scale_list)) formula_scale_list[[parameter]] else NULL,
        prior_random   = if(!is.null(formula_random_prior_list)) formula_random_prior_list[[parameter]] else NULL,
        random_effects_compile = if(!is.null(formula_random_effects_compile_list)) formula_random_effects_compile_list[[parameter]] else NULL)
    }

    .bt_check_formula_node_clashes(formula_output, prior_list)
    formula_prior_output <- do.call(c, unname(lapply(
      formula_output,
      function(output) output[["prior_list"]]
    )))
    formula_add_parameters <- unique(unlist(lapply(
      formula_output,
      function(output) output[["add_parameters"]]
    )), use.names = FALSE)
    combined_prior_list <- .complete_factor_metadata_prior_list(c(
      formula_prior_output,
      prior_list
    ))
    expression_parameter_names <- c(
      JAGS_to_monitor(combined_prior_list),
      add_parameters,
      formula_add_parameters
    )
    for(parameter in names(formula_output)){
      formula_output[[parameter]]$formula_design <- .bt_formula_point_finalize(
        formula_output[[parameter]]$formula_design, combined_prior_list,
        formula_data_list[[parameter]], data, expression_parameter_names
      )
      formula_output[[parameter]]$formula_design <-
        .bt_formula_expression_finalize_design(
          design = formula_output[[parameter]]$formula_design,
          formula_data = formula_data_list[[parameter]],
          model_data = data,
          parameter_names = expression_parameter_names,
          forbidden_parameters = names(formula_output),
          context = paste0(
            "JAGS_fit() expression for parameter '", parameter, "'"
          )
        )
    }

    # merge with the rest of the input
    prior_list <- combined_prior_list
    formula_generated_data <- list()
    for(parameter in names(formula_output)){
      formula_generated_data <- .bt_formula_expression_merge_jags_data(
        formula_generated_data,
        formula_output[[parameter]][["data"]],
        context = "JAGS_fit() formula data"
      )
    }
    data <- .bt_formula_expression_merge_jags_data(
      formula_generated_data,
      data,
      context = "JAGS_fit() model data"
    )
    formula_syntax <- paste0(lapply(formula_output, function(output) output[["formula_syntax"]]), collapse = "")

    for(parameter in names(formula_output)){
      completed_scale <- .bt_formula_scale_finalize(
        formula_output[[parameter]]$formula_scale,
        design = formula_output[[parameter]]$formula_design,
        prior_list = prior_list, model_data = data, owner_scope = "fit"
      )
      formula_output[[parameter]]$formula_scale <- completed_scale
      formula_output[[parameter]]$formula_design$formula_scale <- completed_scale
    }
    # collect formula_scale information
    formula_scale_info <- lapply(formula_output, function(output) output[["formula_scale"]])
    formula_scale_info <- formula_scale_info[!sapply(formula_scale_info, is.null)]
    if(length(formula_scale_info) == 0) formula_scale_info <- NULL
    formula_design_info <- lapply(formula_output, function(output) output[["formula_design"]])
    formula_jags_modules <- unique(unlist(lapply(formula_output, function(output) output[["jags_modules"]])), use.names = FALSE)
    formula_required_packages <- unique(unlist(lapply(formula_output, function(output) output[["required_packages"]])), use.names = FALSE)

    # add the formula syntax to the model syntax
    opening_bracket <- .JAGS_model_block_opening(model_syntax)
    syntax_start    <- substr(model_syntax, 1, opening_bracket)
    syntax_end      <- substr(model_syntax, opening_bracket + 1, nchar(model_syntax))
    model_syntax    <- paste0(syntax_start, "\n", formula_syntax, "\n", syntax_end)
  }else{
    formula_scale_info <- NULL
    formula_design_info <- NULL
    formula_add_parameters <- character()
    formula_jags_modules <- character()
    formula_required_packages <- character()
  }

  # A parameter that a formula monitors for its own needs (the scalar SD source
  # of a random-effect block) and that 'prior_list' already samples is one
  # parameter: it is monitored through its prior. A user's 'add_parameters'
  # entry that repeats a prior is still refused below.
  if(length(formula_add_parameters) > 0L){
    formula_add_parameters <- formula_add_parameters[
      !.bt_jags_parameter_base(formula_add_parameters) %in%
        .bt_jags_parameter_base(c(names(prior_list), JAGS_to_monitor(prior_list)))
    ]
  }
  add_parameters <- unique(c(add_parameters, formula_add_parameters))
  jags_modules <- unique(c(jags_modules, formula_jags_modules))
  required_packages <- unique(c(required_packages, formula_required_packages))

  if(.JAGS_prior_list_uses_BayesTools_module(prior_list)){
    jags_modules <- unique(c(jags_modules, "BayesTools"))
    required_packages <- unique(c(required_packages, "BayesTools"))
  }
  if(is.null(formula_list)){
    prior_list <- .complete_factor_metadata_prior_list(prior_list)
  }
  nested_weights <- .bt_nested_weights_bind(
    prior_list, model_syntax, c(add_parameters, .bt_jags_parameter_base(sampler_probabilities)),
    sampler_nested_weights, data = data
  )
  prior_list <- nested_weights$prior_list
  .bt_validate_jags_add_parameters(add_parameters, prior_list)

  backend_monitor <- .bt_add_backend_anchor(
    model_syntax = model_syntax,
    data = data,
    prior_list = prior_list,
    add_parameters = add_parameters,
    monitor = c(JAGS_to_monitor(prior_list), add_parameters)
  )
  model_syntax <- backend_monitor$model_syntax
  backend_anchor <- backend_monitor$backend_anchor
  if(autofit || check_readiness){
    .bt_convergence_validate_monitor_names(
      monitor = autofit_control[["monitor"]],
      monitored_names = backend_monitor$monitor
    )
  }

  # Sampler declarations: marker nodes in the final model text and their data.
  # The markers are known only here, so the module requirement is added before
  # packages and modules are loaded (workers included).
  sampler_declarations <- .bt_sampler_declarations(
    model_text = JAGS_add_priors(syntax = model_syntax, prior_list = prior_list),
    data = data,
    prior_list = prior_list,
    add_parameters = add_parameters,
    sampler_gates = sampler_gates,
    sampler_blocks = sampler_blocks
  )
  # Package and readout markers are terminal. Gate analysis above concerns
  # dependencies in the scientific model.
  sampler_declarations <- .bt_sampler_apply_plan(sampler_declarations, backend_monitor$monitor,
    prior_list, formula_design_info, data, sampler_joint, sampler_probabilities, jags_modules)
  if(!is.null(nested_weights$status)){
    sampler_declarations$declarations$nested_weights <- nested_weights$status
  }
  if(sampler_declarations$uses_module){
    jags_modules <- unique(c(jags_modules, "BayesTools"))
    required_packages <- unique(c(required_packages, "BayesTools"))
  }

  ### create the model call
  model_call <- list(
    model     = sampler_declarations$model,
    data      = sampler_declarations$data,
    inits     = .JAGS_fit_inits(prior_list, chains = chains, seed = seed),
    monitor   = backend_monitor$monitor,
    adapt     = adapt,
    burnin    = burnin,
    sample    = sample,
    thin      = thin,
    silent    = silent
  )
  probability_plan <- sampler_declarations$declarations$probabilities
  if(length(probability_plan)){
    model_call$probability_plan <- probability_plan
  }

  # Configure the actual backend topology, including automatic extensions.
  cl <- NULL
  runtime_started <- FALSE
  if(parallel){
    cl <- .JAGS_runtime_cluster_start(
      cores = min(cores, chains), chains = chains,
      packages = unique(c(required_packages, "BayesTools")), runtime_setup = runtime_setup,
      jags_modules = jags_modules, warn = !silent,
      worker_output = worker_output
    )
    on.exit(.JAGS_runtime_cluster_stop(cl), add = TRUE)
    .JAGS_run_runtime_cache(runtime_cache, "restore", chains, cl)
    model_call$cl <- cl
  }else{
    on.exit(.JAGS_finish_runtime_setup(
      if(runtime_started) runtime_setup else NULL, chains), add = TRUE)
    .JAGS_require_packages(required_packages)
    .JAGS_load_modules(jags_modules, warn = !silent)
    runtime_started <- TRUE
    .JAGS_run_runtime_setup(runtime_setup, chains)
    .JAGS_run_runtime_cache(runtime_cache, "restore", chains)
  }

  if(!is.null(seed)){
    set.seed(seed)
  }

  start_time <- .bt_autofit_clock()
  max_time_exceeded <- function(){
    elapsed <- .bt_autofit_elapsed(start_time, autofit_control[["max_time"]])
    !is.null(elapsed) && elapsed[["time"]] > autofit_control[["max_time"]][["time"]]
  }
  restart_warnings <- character()
  # incomplete-adaptation warnings of the backend runs that the fit keeps
  adaptation_warnings <- character()
  # the warning of chains that were not ready at the end of burn-in
  burnin_warnings <- character()
  adaptation <- .bt_adaptation_record(adapt)
  # the runtime reports of the backend runs that the fit keeps
  reports <- .bt_runtime_reports_empty()
  report_run <- function(result){
    reports <<- .bt_runtime_report_after_run(reports, runtime_report, result,
      chains, if(parallel) cl else NULL)
  }
  # set while the readiness check assesses a block of burn-in
  assessing <- FALSE
  run_initial <- function(call){
    assessing <<- FALSE
    run <- withRestarts(
      .bt_collect_adaptation_warnings(.JAGS_run_backend(.bt_rjags_fit, call)),
      BayesTools_jasp_progress_error = function(e, ran){
        if(ran) report_run(e)
        stop(e)
      })
    report_run(run$value)
    # An error of the readiness check itself (such as a convergence monitor
    # that the fitted columns do not contain) is not a sampler failure that
    # new initial values could repair.
    if(assessing){
      stop(run$value)
    }
    if(!inherits(run$value, "error")){
      adaptation_warnings <<- run$warnings
      # only the attempt that the fit keeps warns
      burnin_warnings <<- .bt_burnin_warn(attr(run$value, "burnin_check", exact = TRUE))
    }
    run$value
  }
  # Keep generated observations fixed when a new set of chain inits is tried.
  prepared <- tryCatch(.bt_rjags_realize_data(model_call$model, model_call$data,
    model_call$inits[[1L]][[".RNG.seed"]], model_call$monitor), error = identity)
  if(!inherits(prepared, "error")){
    model_call[c("model", "data")] <- prepared
  }
  # A BUGS data block can generate observations. Determine convergence roles
  # from the actual sampled model and its realized, fully observed inputs.
  sampling_model <- model_call$model
  data_names <- names(model_call$data)[vapply(model_call$data, function(x) !anyNA(x), logical(1))]
  # The convergence roles of the fitted coordinates, as the parameter map
  # stores them, give the selection of the readiness check and of automatic
  # fitting. The coordinates are built once, from the first columns the
  # backend returns, and again only for different columns.
  coordinates_cache <- NULL
  fit_coordinates <- function(columns){
    if(is.null(coordinates_cache) || !identical(coordinates_cache$columns, columns)){
      coordinates_cache <<- list(
        columns = columns,
        coordinates = .bt_build_parameter_coordinates(
          columns = columns,
          monitor_names = model_call$monitor,
          prior_list = prior_list,
          formula_design = formula_design_info,
          formula_scale = formula_scale_info,
          backend_anchor = backend_anchor,
          add_parameters = add_parameters,
          model_syntax = sampling_model,
          data_names = data_names
        )
      )
    }
    coordinates_cache$coordinates
  }
  if(check_readiness){
    model_call$readiness <- list(
      block           = as.integer(readiness_block),
      burnin_max_Rhat = burnin_settings[["burnin_max_Rhat"]],
      burnin_extend   = burnin_settings[["burnin_extend"]],
      max_extensions  = if(autofit) burnin_settings[["burnin_max_extend"]] else 0,
      # the burn-in sweeps that keep the public iteration number of the last
      # retained draw within the integer range of the draw geometry
      extension_room  = .Machine$integer.max -
        (as.double(adapt) + burnin + 1 + (as.double(sample) - 1) * thin),
      autofit         = autofit,
      time_exceeded   = max_time_exceeded,
      assess          = function(draws){
        assessing <<- TRUE
        values <- .bt_burnin_block_Rhat(
          block = draws,
          coordinates = fit_coordinates(colnames(draws[[1L]])),
          prior_list = prior_list,
          monitor = autofit_control[["monitor"]],
          check_indicators = autofit_control[["check_indicators"]]
        )
        assessing <<- FALSE
        values
      }
    )
  }
  if(inherits(prepared, "error")){
    fit <- prepared
    report_run(fit)
  }else if(isTRUE(dots[["is_JASP"]])){
    prefix <- dots[["is_JASP_prefix"]]
    model_call$jasp_prefix <- if(is.null(prefix)) "" else paste0(prefix, ": ")
    model_call$retain_workers <- autofit
    fit <- run_initial(model_call)
  }else{
    if(is.null(autofit_control[["restarts"]])){
      fit <- run_initial(model_call)
    }else{
      for(i in 1:autofit_control[["restarts"]]){
        fit <- run_initial(model_call)
        if(!inherits(fit, "error")){
          break
        }else if(inherits(fit, c("BayesTools_JAGS_worker_connection_error",
                                "BayesTools_sampler_probabilities_error"))){
          break
        }else if(i < autofit_control[["restarts"]]){
          restart_warnings <- c(
            restart_warnings,
            paste0(
              "JAGS fitting attempt ", i, " failed and was restarted: ",
              conditionMessage(fit), "."
            )
          )
          if(!silent){
            warning(
              restart_warnings[[length(restart_warnings)]],
              call. = FALSE,
              immediate. = TRUE
            )
          }
          # restart with different inits, seeded from the restart-seed stream
          restart_seed <- if(!is.null(seed)) .JAGS_restart_seeds(seed, i)[[i]]
          model_call$inits <- .JAGS_fit_inits(prior_list, chains = chains, seed = restart_seed)
        }
      }
    }
  }



  if(inherits(fit, "error") & !silent)
    warning(paste0("The model estimation failed with the following error: ", fit$message), immediate. = TRUE)

  extension_failed <- FALSE
  last_valid_runtime_state <- NULL
  last_valid_capture_warnings <- character()
  captured_after_success <- FALSE

  if(autofit && !inherits(fit, "error")){

    # An unsuccessful first extension must retain the cache belonging to the
    # initial valid fit, before the backend has a chance to mutate it.
    captured <- .JAGS_capture_runtime_cache(
      runtime_cache, chains, if(parallel) cl else NULL)
    last_valid_runtime_state <- captured$state
    last_valid_capture_warnings <- captured$warnings
    captured_after_success <- TRUE

    autofit_coordinates <- function(fit){
      fit_coordinates(colnames(.extract_posterior_samples(fit, as_list = FALSE)))
    }
    # sample_extend = "auto" sizes each extension from the diagnostics of the
    # failed check, which then assesses every selected quantity
    auto_extend <- identical(autofit_control[["sample_extend"]], "auto")
    autofit_converged <- function(fit){
      convergence_fit <- fit
      attr(convergence_fit, "prior_list") <- prior_list
      .bt_check_convergence(
        fit = convergence_fit,
        coordinates = autofit_coordinates(fit),
        max_Rhat = autofit_control[["max_Rhat"]],
        min_ESS = autofit_control[["min_ESS"]],
        max_error = autofit_control[["max_error"]],
        max_SD_error = autofit_control[["max_SD_error"]],
        fail_fast = !auto_extend,
        check_indicators = autofit_control[["check_indicators"]],
        monitor = autofit_control[["monitor"]],
        allow_not_assessable = autofit_control[["allow_not_assessable"]]
      )
    }
    converged <- autofit_converged(fit)
    itteration <- 1
    last_valid_fit <- fit

    if(!converged && isTRUE(dots[["is_JASP"]]))
      .JASP_progress_bar_start(n = if (!is.null(autofit_control[["max_extend"]])) autofit_control[["max_extend"]] else 10, label = paste0(if(!is.null(dots[["is_JASP_prefix"]])) paste0(dots[["is_JASP_prefix"]], ": "), "Extending the model (autofit)"))

    while(!converged){

      elapsed <- .bt_autofit_elapsed(start_time, autofit_control[["max_time"]])
      stop_constraint <- if(!is.null(elapsed) && elapsed[["time"]] > autofit_control[["max_time"]][["time"]]){
        "max_time"
      }else if(!is.null(autofit_control[["max_extend"]]) && itteration > autofit_control[["max_extend"]]){
        "max_extend"
      }
      if(is.null(stop_constraint) && auto_extend){
        # at most what the time left allows at the rate of the draws so far;
        # the clock started before the initial run, so the rate includes its
        # adaptation and burn-in
        time_cap <- .bt_autofit_time_cap(autofit_control[["max_time"]], elapsed[["secs"]],
          draws = NROW(fit[["mcmc"]][[1L]]))
        if(time_cap < 100){
          stop_constraint <- "max_time"
        }
      }
      if(!is.null(stop_constraint)){
        warning_message <- paste0("The automatic model fitting was terminated due to the '", stop_constraint, "' constraint.")
        fit <- .bt_append_fit_warnings(fit, warning_message)
        if(!silent){
          warning(warning_message, immediate. = TRUE)
        }
        break
      }

      sample_extend <- if(auto_extend){
        min(.bt_autofit_extension_size_from(converged, fit, autofit_control), time_cap,
          .bt_autofit_iteration_cap(thin))
      }else{
        autofit_control[["sample_extend"]]
      }
      if(!.bt_extension_fits_iterations(fit, sample_extend)){
        warning_message <- .bt_iteration_limit_message(sample_extend)
        fit <- .bt_append_fit_warnings(fit, warning_message)
        if(!silent){
          warning(warning_message, immediate. = TRUE)
        }
        break
      }
      run <- .bt_run_extension(fit, sample = sample_extend,
        adapt = adapt, cl = cl, silent = silent)
      extension <- run$extension
      report_run(extension)

      if(inherits(extension, "error")){
        warning_message <- paste0(
          "The model extension failed; returning the last valid fit. ",
          "Backend error: ",
          conditionMessage(extension)
        )
        fit <- .bt_append_fit_warnings(last_valid_fit, warning_message)
        extension_failed <- TRUE
        warning(
          warning_message,
          call. = FALSE,
          immediate. = TRUE
        )
        break
      }

      fit <- extension
      last_valid_fit <- fit
      adaptation <- .bt_adaptation_add_extension(adaptation, adapt, run$rebuilt)
      adaptation_warnings <- c(adaptation_warnings, run$warnings)
      captured <- .JAGS_capture_runtime_cache(
        runtime_cache, chains, if(parallel) cl else NULL)
      last_valid_runtime_state <- captured$state
      last_valid_capture_warnings <- captured$warnings
      captured_after_success <- TRUE

      converged <- autofit_converged(fit)
      itteration <- itteration + 1

      if(isTRUE(dots[["is_JASP"]]))
        .JASP_progress_bar_tick()
    }
  }

  # add information to the fitted object
  if(!inherits(fit, "error")){
    fit[["live_model"]] <- NULL
  }
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "add_parameters") <- add_parameters
  attr(fit, "required_packages") <- required_packages
  attr(fit, "jags_modules") <- jags_modules
  attr(fit, "runtime_setup") <- runtime_setup
  attr(fit, "backend_anchor") <- backend_anchor
  if(!inherits(fit, "error") && !is.null(fit$samplers) &&
     NROW(sampler_declarations$declarations$gates) > 0L){
    sampler_declarations$declarations$gate_claims <- .bt_sampler_fit_gate_claims(
      sampler_declarations$declarations$gates, prior_list, fit)
  }
  attr(fit, "sampler_declarations") <- sampler_declarations$declarations
  if(!inherits(fit, "error")){
    attr(fit, "adaptation") <- adaptation
    fit <- .bt_append_fit_warnings(fit, c(restart_warnings, adaptation_warnings, burnin_warnings))
  }
  if(!is.null(formula_scale_info)){
    # Keep formula_scale as a nested list keyed by parameter name
    # Each element contains the scaling info for that parameter's predictors
    attr(fit, "formula_scale") <- formula_scale_info
  }
  if(!is.null(formula_design_info)){
    attr(fit, "formula_design") <- formula_design_info
  }

  fit <- .bt_attach_parameter_map(
    fit,
    monitor_names = model_call$monitor,
    model_syntax = sampling_model,
    data_names = data_names
  )
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)

  attr(fit, "runtime_cache") <- runtime_cache
  attr(fit, "runtime_report") <- runtime_report
  if(!inherits(fit, "error")){
    if(isTRUE(extension_failed) || isTRUE(captured_after_success)){
      attr(fit, "runtime_state") <- last_valid_runtime_state
      fit <- .bt_append_fit_warnings(fit, last_valid_capture_warnings)
    }else{
      captured <- .JAGS_capture_runtime_cache(
        runtime_cache, chains, if(parallel) cl else NULL)
      attr(fit, "runtime_state") <- captured$state
      fit <- .bt_append_fit_warnings(fit, captured$warnings)
    }
    attr(fit, "runtime_report_values") <- .JAGS_runtime_report_value(reports$runs)
    fit <- .bt_append_fit_warnings(fit, reports$warnings)
  }
  return(fit)
}

.bt_validate_jags_formula_list_names <- function(x, name){

  if(length(x) == 0L){
    return(invisible(TRUE))
  }

  x_names <- names(x)
  if(is.null(x_names) || length(x_names) != length(x) ||
     anyNA(x_names) || any(!nzchar(x_names))){
    stop(
      "The '", name, "' argument must be a fully named list.",
      call. = FALSE
    )
  }
  if(anyDuplicated(x_names)){
    duplicate_names <- unique(x_names[duplicated(x_names)])
    stop(
      "The '", name, "' argument must not contain duplicate names ('",
      paste(duplicate_names, collapse = "', '"),
      "').",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_validate_jags_formula_lists <- function(
    formula_list,
    formula_data_list,
    formula_prior_list,
    formula_random_prior_list,
    formula_random_effects_compile_list,
    formula_scale_list){

  .bt_validate_jags_formula_list_names(formula_list, "formula_list")
  formula_names <- names(formula_list)
  companion_lists <- list(
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list,
    formula_scale_list = formula_scale_list
  )

  for(name in names(companion_lists)){
    x <- companion_lists[[name]]
    if(length(x) == 0L){
      next
    }
    if(length(formula_names) == 0L){
      stop(
        "The '", name, "' argument cannot be supplied without 'formula_list'.",
        call. = FALSE
      )
    }
    .bt_validate_jags_formula_list_names(x, name)
    unrecognized <- setdiff(names(x), formula_names)
    if(length(unrecognized) > 0L){
      stop(
        "The '", paste(unrecognized, collapse = "', '"),
        "' objects are not recognized by the '", name, "' argument.",
        call. = FALSE
      )
    }
  }

  for(name in c("formula_data_list", "formula_prior_list")){
    missing <- setdiff(formula_names, names(companion_lists[[name]]))
    if(length(missing) > 0L){
      stop(
        "The '", paste(missing, collapse = "', '"),
        "' objects are missing in the '", name, "' argument.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# JAGS node names of the formula coefficients must be unique across the
# formulas and the model's own 'prior_list'. Different formula parameters and
# terms can generate the same node name (formula 'a' with term 'b_x' and
# formula 'a_b' with term 'x' both define 'a_b_x'); the clash is reported by
# the formula parameters and terms that generate it.
.bt_check_formula_node_clashes <- function(formula_output, prior_list = NULL){

  owners <- do.call(rbind, c(
    lapply(names(formula_output), function(parameter){
      nodes <- names(formula_output[[parameter]][["prior_list"]])
      if(length(nodes) == 0L){
        return(NULL)
      }
      name_map <- formula_output[[parameter]][["formula_design"]][["name_map"]]
      terms <- name_map$term[match(nodes, name_map$jags_name)]
      terms[is.na(terms)] <- nodes[is.na(terms)]
      data.frame(
        node  = nodes,
        owner = paste0("formula parameter '", parameter, "' (term '", terms, "')"),
        stringsAsFactors = FALSE
      )
    }),
    list(data.frame(
      node  = as.character(names(prior_list)),
      owner = rep("'prior_list'", length(prior_list)),
      stringsAsFactors = FALSE
    ))
  ))
  clashes <- unique(owners$node[duplicated(owners$node)])
  if(length(clashes) == 0L){
    return(invisible(TRUE))
  }
  messages <- vapply(clashes, function(node){
    node_owners <- owners$owner[owners$node == node]
    paste0(
      paste(node_owners, collapse = " and "),
      " define the same JAGS node '", node, "'"
    )
  }, character(1))
  stop(
    paste(messages, collapse = "; "),
    ". Rename a formula parameter or predictor so that the node names differ.",
    call. = FALSE
  )
}

# The parameter name without its index: "tau[1]" is a coordinate of "tau".
.bt_jags_parameter_base <- function(x){

  sub("\\[.*$", "", x)
}

.bt_validate_jags_add_parameters <- function(add_parameters, prior_list){

  if(length(add_parameters) == 0L){
    return(invisible(TRUE))
  }
  if(any(!nzchar(add_parameters))){
    stop(
      "The 'add_parameters' argument cannot contain empty parameter names.",
      call. = FALSE
    )
  }
  if(length(prior_list) == 0L){
    return(invisible(TRUE))
  }

  prior_parameters <- unique(c(names(prior_list), JAGS_to_monitor(prior_list)))
  overlap <- add_parameters[
    .bt_jags_parameter_base(add_parameters) %in%
      .bt_jags_parameter_base(prior_parameters)
  ]
  if(length(overlap) > 0L){
    stop(
      "The 'add_parameters' argument must not include parameters already ",
      "monitored through 'prior_list': ",
      paste(unique(overlap), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @rdname JAGS_fit
JAGS_extend <- function(fit, autofit_control = list(max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05, max_time = list(time = 60, unit = "mins"), sample_extend = 1000, restarts = 10, max_extend = 10, check_indicators = FALSE, monitor = NULL, allow_not_assessable = FALSE),
                        parallel = FALSE, cores = NULL, silent = TRUE,
                        runtime_setup = attr(fit, "runtime_setup", exact = TRUE),
                        runtime_cache = attr(fit, "runtime_cache", exact = TRUE),
                        worker_output = NULL, adapt = NULL,
                        runtime_report = attr(fit, "runtime_report", exact = TRUE)){

  if(inherits(fit, "BayesTools_draws_view")){
    stop(.bt_draws_view_condition(
      c("BayesTools_draws_view_sampling_unavailable", "BayesTools_sampling_unavailable"),
      "Sampling extension is unavailable for a derived-draw view. Extend 'fit$original_fit' with 'JAGS_extend()' and regenerate the view with 'JAGS_with_draws()'.",
      operation = "extension"
    ))
  }
  if(inherits(fit, "runjags")){
    .bt_stop_refit_required(
      "This fit was produced by runjags. Refit it with JAGS_fit() before extending its chains."
    )
  }
  if(!inherits(fit, "BayesTools_fit") || inherits(fit, "error")){
    stop("'fit' must be a fitted 'BayesTools_fit'.", call. = FALSE)
  }
  backend <- fit[["backend"]]
  fitted_jags <- if(is.list(backend)) backend[["jags"]]
  jags_major <- rjags::jags.version()[, 1L]
  if(!is.character(fitted_jags) || length(fitted_jags) != 1L || is.na(fitted_jags) ||
     !startsWith(fitted_jags, paste0(jags_major, "."))){
    .bt_stop_refit_required(
      "The fit has no compatible JAGS version record. Refit it with JAGS_fit() before extending its chains."
    )
  }
  state <- fit[["end_state"]]
  chains <- length(fit[["mcmc"]])
  if(chains == 0L || !is.list(state) || length(state) != chains ||
     !all(vapply(state, function(chain){
       is.list(chain) && is.character(chain[[".RNG.name"]]) &&
         length(chain[[".RNG.name"]]) == 1L && !is.na(chain[[".RNG.name"]]) &&
         nzchar(chain[[".RNG.name"]]) && is.numeric(chain[[".RNG.state"]]) &&
         length(chain[[".RNG.state"]]) > 0L
     }, logical(1)))){
    .bt_stop_refit_required(
      "The fit does not contain the final state of every chain. Refit it with JAGS_fit() before extending."
    )
  }

  check_bool(parallel, "parallel", allow_NA = FALSE)
  check_int(cores, "cores", lower = 1, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(silent, "silent", allow_NA = FALSE)
  # Only NULL inherits; any other value is validated as for JAGS_fit().
  if(!is.null(adapt)){
    .JAGS_check_adapt(adapt)
  }
  worker_output <- .JAGS_validate_worker_output(worker_output)
  .JAGS_validate_runtime_setup(runtime_setup)
  .JAGS_validate_runtime_cache(runtime_cache)
  .JAGS_validate_runtime_report(runtime_report)
  runtime_state <- attr(fit, "runtime_state", exact = TRUE)
  last_valid_runtime_state <- runtime_state
  # Only the local input copy loses these attributes. Retained snapshots must
  # never travel inside the fit to every extension worker. The input's
  # runtime report describes earlier runs, not those of this call.
  attr(fit, "runtime_state") <- NULL
  attr(fit, "runtime_cache") <- NULL
  attr(fit, "runtime_report_values") <- NULL
  attr(fit, "runtime_report") <- NULL
  # The backend would continue a compiled model still alive in this session,
  # so that an extension depended on earlier in-memory extensions of the same
  # object. Extending from the stored chain states alone makes it depend only
  # on the fitted object, as after saving and reloading it.
  fit[["live_model"]] <- NULL
  fitted_parameter_map <- parameter_map(fit)
  backend_anchor     <- attr(fit, "backend_anchor", exact = TRUE)
  sampler_declarations <- attr(fit, "sampler_declarations", exact = TRUE)
  formula_design     <- attr(fit, "formula_design", exact = TRUE)

  coordinates <- .bt_require_fit_contract(fit)
  .bt_R2_validate_fitted_owners(fit)
  if(!is.null(formula_design)){
    if(!is.list(formula_design) || is.null(names(formula_design)) ||
       any(!nzchar(names(formula_design))) || anyDuplicated(names(formula_design))){
      .bt_stop_refit_required(
        "JAGS_extend() cannot preserve malformed formula-design metadata. Refit the model with this version of BayesTools."
      )
    }
    JAGS_formula_name_map(fit)
    for(parameter in names(formula_design)){
      .bt_regularized_check_fitted_priors(formula_design[[parameter]], attr(fit, "prior_list", exact = TRUE))
      .bt_validate_formula_design_replay_schema(
        formula_design[[parameter]],
        context = paste0("JAGS_extend() formula '", parameter, "'")
      )
      .bt_formula_scale_finalized_check(formula_design[[parameter]]$formula_scale,
        parameter, require_owner = TRUE)
      if(!identical(attr(formula_design[[parameter]]$formula_scale, "unscale_design", exact = TRUE),
                    attr(attr(fit, "formula_scale", exact = TRUE)[[parameter]], "unscale_design", exact = TRUE))){
        .bt_stop_refit_required("JAGS_extend() formula declaration carriers disagree. Refit the model with this version of BayesTools.")
      }
    }
  }
  JAGS_draw_geometry(fit)

  # extract fitting information
  prior_list        <- attr(fit, "prior_list")
  model_syntax      <- attr(fit, "model_syntax")
  required_packages <- attr(fit, "required_packages")
  jags_modules      <- attr(fit, "jags_modules")
  add_parameters    <- attr(fit, "add_parameters")
  formula_scale     <- attr(fit, "formula_scale")
  prior_list        <- .complete_factor_metadata_prior_list(prior_list)
  if(is.null(add_parameters)){
    add_parameters <- character()
  }
  .bt_validate_jags_add_parameters(add_parameters, prior_list)
  autofit_control <- JAGS_check_and_list_autofit_settings(autofit_control)
  # The stored convergence roles give the selection of the public check.
  .bt_convergence_validate_monitor(
    fit = fit,
    coordinates = coordinates,
    monitor = autofit_control[["monitor"]]
  )
  # An explicit 'adapt' applies to this call; the record keeps the initial one.
  adaptation <- .bt_fit_adaptation_record(fit)
  if(is.null(adapt)){
    adapt <- adaptation$initial
  }

  cl <- NULL
  runtime_started <- FALSE
  if(parallel){
    if(is.null(cores)){
      cores <- chains
    }
    cl <- .JAGS_runtime_cluster_start(
      cores = min(cores, chains), chains = chains,
      packages = unique(c(required_packages, "BayesTools")), runtime_setup = runtime_setup,
      jags_modules = jags_modules, warn = !silent,
      worker_output = worker_output
    )
    on.exit(.JAGS_runtime_cluster_stop(cl), add = TRUE)
    .JAGS_run_runtime_cache(runtime_cache, "restore", chains, cl, runtime_state)
    runtime_state <- NULL
  }else{
    on.exit(.JAGS_finish_runtime_setup(
      if(runtime_started) runtime_setup else NULL, chains), add = TRUE)
    .JAGS_require_packages(required_packages)
    .JAGS_load_modules(jags_modules, warn = !silent)
    runtime_started <- TRUE
    .JAGS_run_runtime_setup(runtime_setup, chains)
    .JAGS_run_runtime_cache(runtime_cache, "restore", chains, state = runtime_state)
    runtime_state <- NULL
  }

  start_time <- .bt_autofit_clock()
  supplied_draws <- NROW(fit[["mcmc"]][[1L]])
  iteration  <- 0
  converged  <- FALSE
  last_valid_fit <- fit
  extension_failed <- FALSE
  last_valid_capture_warnings <- character()
  captured_after_success <- FALSE
  reports <- .bt_runtime_reports_empty()
  # sample_extend = "auto" sizes each extension when the loop is about to run
  # it, from the diagnostics of the last failed check (assessing every
  # selected quantity); the first one from a check of the supplied fit
  auto_extend <- identical(autofit_control[["sample_extend"]], "auto")
  extend_converged <- function(fit, max_Rhat = autofit_control[["max_Rhat"]]){
    convergence_fit <- fit
    attr(convergence_fit, "prior_list") <- prior_list
    .bt_check_convergence(
      fit = convergence_fit,
      coordinates = coordinates,
      max_Rhat = max_Rhat,
      min_ESS = autofit_control[["min_ESS"]],
      max_error = autofit_control[["max_error"]],
      max_SD_error = autofit_control[["max_SD_error"]],
      fail_fast = !auto_extend,
      check_indicators = autofit_control[["check_indicators"]],
      monitor = autofit_control[["monitor"]],
      allow_not_assessable = autofit_control[["allow_not_assessable"]]
    )
  }
  failed_check <- NULL

  while(!converged){

    elapsed <- .bt_autofit_elapsed(start_time, autofit_control[["max_time"]])
    stop_constraint <- if(!is.null(elapsed) && elapsed[["time"]] > autofit_control[["max_time"]][["time"]]){
      "max_time"
    }else if(!is.null(autofit_control[["max_extend"]]) && iteration >= autofit_control[["max_extend"]]){
      "max_extend"
    }
    if(is.null(stop_constraint) && auto_extend){
      # at most what the time left allows at the rate of the draws of this
      # call; its first extension has no rate yet and is not capped
      time_cap <- .bt_autofit_time_cap(autofit_control[["max_time"]], elapsed[["secs"]],
        draws = NROW(fit[["mcmc"]][[1L]]) - supplied_draws)
      if(time_cap < 100){
        stop_constraint <- "max_time"
      }
    }
    if(!is.null(stop_constraint)){
      warning_message <- paste0("The automatic model fitting was terminated due to the '", stop_constraint, "' constraint.")
      fit <- .bt_append_fit_warnings(fit, warning_message)
      if(!silent){
        warning(warning_message, immediate. = TRUE)
      }
      break
    }

    sample_extend <- if(auto_extend){
      if(is.null(failed_check)){
        # With one chain, the check skips R-hat with a warning that the check
        # after this extension gives; this one only sizes the extension.
        failed_check <- extend_converged(fit,
          max_Rhat = if(chains > 1L) autofit_control[["max_Rhat"]])
      }
      min(.bt_autofit_extension_size_from(failed_check, fit, autofit_control), time_cap,
        .bt_autofit_iteration_cap(fit[["thin"]]))
    }else{
      autofit_control[["sample_extend"]]
    }
    if(!.bt_extension_fits_iterations(fit, sample_extend)){
      warning_message <- .bt_iteration_limit_message(sample_extend)
      fit <- .bt_append_fit_warnings(fit, warning_message)
      if(!silent){
        warning(warning_message, immediate. = TRUE)
      }
      break
    }
    run <- .bt_run_extension(fit, sample = sample_extend,
      adapt = adapt, cl = cl, silent = silent)
    extension <- run$extension
    reports <- .bt_runtime_report_after_run(reports, runtime_report, extension,
      chains, if(parallel) cl else NULL)

    if(inherits(extension, "error")){
      warning_message <- paste0(
        "The model extension failed; returning the last valid fit. ",
        "Backend error: ",
        conditionMessage(extension)
      )
      fit <- .bt_append_fit_warnings(last_valid_fit, warning_message)
      extension_failed <- TRUE
      warning(
        warning_message,
        call. = FALSE,
        immediate. = TRUE
      )
      break
    }

    # The backend returns a new object; carry the fit's warning history.
    fit <- .bt_append_fit_warnings(
      extension,
      c(attr(last_valid_fit, "warnings", exact = TRUE), run$warnings)
    )
    adaptation <- .bt_adaptation_add_extension(adaptation, adapt, run$rebuilt)
    attr(fit, "adaptation") <- adaptation
    last_valid_fit <- fit
    captured <- .JAGS_capture_runtime_cache(
      runtime_cache, chains, if(parallel) cl else NULL)
    last_valid_runtime_state <- captured$state
    last_valid_capture_warnings <- captured$warnings
    captured_after_success <- TRUE
    converged <- extend_converged(fit)

    if(!converged){
      iteration <- iteration + 1
      failed_check <- converged
    }
  }

  # add information to the fitted object
  fit[["live_model"]] <- NULL
  attr(fit, "prior_list")   <- prior_list
  attr(fit, "model_syntax") <- model_syntax
  attr(fit, "add_parameters") <- add_parameters
  attr(fit, "required_packages") <- required_packages
  attr(fit, "jags_modules") <- jags_modules
  attr(fit, "runtime_setup") <- runtime_setup
  attr(fit, "backend_anchor") <- backend_anchor
  if(!is.null(fit$samplers) && NROW(sampler_declarations$gates) > 0L){
    sampler_declarations$gate_claims <- .bt_sampler_fit_gate_claims(
      sampler_declarations$gates, prior_list, fit)
  }
  attr(fit, "sampler_declarations") <- sampler_declarations
  if(!is.null(formula_scale)){
    attr(fit, "formula_scale") <- formula_scale
  }
  if(!is.null(formula_design)){
    attr(fit, "formula_design") <- formula_design
  }

  attr(fit, "parameter_map") <- fitted_parameter_map
  fit <- .bt_attach_draw_geometry(fit)
  fit <- .bt_attach_fit_contract(fit)

  attr(fit, "runtime_cache") <- runtime_cache
  attr(fit, "runtime_report") <- runtime_report
  if(isTRUE(extension_failed) || isTRUE(captured_after_success)){
    attr(fit, "runtime_state") <- last_valid_runtime_state
    fit <- .bt_append_fit_warnings(fit, last_valid_capture_warnings)
  }else{
    captured <- .JAGS_capture_runtime_cache(
      runtime_cache, chains, if(parallel) cl else NULL)
    attr(fit, "runtime_state") <- captured$state
    fit <- .bt_append_fit_warnings(fit, captured$warnings)
  }
  attr(fit, "runtime_report_values") <- .JAGS_runtime_report_value(reports$runs)
  fit <- .bt_append_fit_warnings(fit, reports$warnings)
  return(fit)
}

# The clock of the 'max_time' limit of automatic fitting.
.bt_autofit_clock <- function(){
  Sys.time()
}

# The time on that clock since 'start_time', from one reading: 'time' in the
# unit of 'max_time', which its limit is compared with, and 'secs' in
# seconds, which the time cap of automatic extensions uses; NULL without a
# limit, and the clock is then not read.
.bt_autofit_elapsed <- function(start_time, max_time){

  if(is.null(max_time)){
    return(NULL)
  }
  now <- .bt_autofit_clock()
  c(
    time = as.numeric(difftime(now, start_time, units = max_time[["unit"]])),
    secs = as.numeric(difftime(now, start_time, units = "secs"))
  )
}
