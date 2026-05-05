# Research Log

## 2026-04-26T10:18:20.230999+00:00 Search Query

Query: `pushforward distribution change of variables additive model prior predictive density Bayesian regression priors`

Task: `Prior density transformation for standardized regression coefficients`

Top results: kneib2011jrssscas_high, brown2010ba_inference, kobayashi2008jmva_bayesian, dawid2015ba_bayesian, egidi2022cjs_avoiding, a2008misc_weakly, ley2012je_mixtures, vehtari2016misc_practical

## 2026-04-26T10:19:01.828101+00:00 Search Query

Query: `linear transformation induced prior density convolution arbitrary priors Bayesian regression standardization`

Task: `Prior density transformation for standardized regression coefficients`

Top results: castillo2015as_bayesian, brown2010ba_inference, zhang2020jasa_bayesian, kobayashi2008jmva_bayesian, kim2008jmva_class, boxgep1964misc_analysis, taddy2011jasa_dynamic, graham2017ejs_asymptotically

## 2026-04-26T10:19:43.767505+00:00 Search Query

Query: `induced priors reparameterization Bayesian model coefficients standardization interactions`

Task: `Prior density transformation for standardized regression coefficients`

Top results: a2008misc_weakly, klugkist2007csada_bayes, allen2003technometrics_supersaturated, anonndmisc_untitled2, mulder2014csada_prior, gvanvalkenhoef2015rsm_automated, consonni2008test_compatible, schuurman2016psychmethods_how

## 2026-04-26T10:21:19.980132+00:00 Coverage Report: Prior density transformation for standardized regression coefficients

Coverage snapshot:
- Task sources: 0 (citation-ready: 0)
- Discovery candidates considered: 31 from 3 recent query snapshots
- Judgment summary: pinned=0, relevant=0, supports_claim=0, irrelevant=0, unreviewed=0
- Claim role coverage: support=0, foundational=0, contrastive=0
- Method spread: bayesian (8), regression (5), nonparametric (3)
- Relation coverage: references=4, same_method=3, same_author=2, same_journal=2
- Foundational follow-ups still missing or not ready: carvalho2010biometrika_horseshoe, roberttibshirani1996jrsssbsm_regression, ghosal2000as_convergence, castillo2012as_needles
- Comparison/robustness leads: bornn2010cjs_efficient, carvalho2010biometrika_horseshoe, roberttibshirani1996jrsssbsm_regression
- Unreviewed high-signal hits: castillo2015as_bayesian, lemoine2019o_moving, hans2010statcomp_model, kobayashi2008jmva_bayesian
- Likely blind spots: top candidates are still unreviewed; no clear foundational claim support surfaced yet

## 2026-04-26T10:24:50.272108+00:00 Coverage Report: Read-only search for literature relevant to computational transformation of arbitrary custom priors under standardization, induced priors, linear transformations, additive models, interactions, and density pushforwards.

Coverage snapshot:
- Task sources: 0 (citation-ready: 0)
- Discovery candidates considered: 43 from 3 recent query snapshots
- Judgment summary: pinned=0, relevant=0, supports_claim=0, irrelevant=0, unreviewed=0
- Claim role coverage: support=10, foundational=1, contrastive=2
- Method spread: bayesian (11), computational (9), multivariate (7)
- Relation coverage: references=4, same_method=4, same_journal=2, same_author=1
- Foundational follow-ups still missing or not ready: heckerman1995ml_learning, dawid1993as_hyper, andersson1997as_characterization, cooper1992ml_bayesian
- Comparison/robustness leads: baringhaus1989cssc_recent, gelman2007sm_2, gvanvalkenhoef2015rsm_automated
- Unreviewed high-signal hits: zwet2022tas_proposal, geiger2002as_parameter, waldmann2008g_efficient, royston1993s_toolkit
- Likely blind spots: top candidates are still unreviewed

## 2026-04-26T10:26:32.800092+00:00 Coverage Report: Read-only search for computational transformation of arbitrary custom priors under standardization.

Coverage snapshot:
- Task sources: 0 (citation-ready: 0)
- Discovery candidates considered: 30 from 3 recent query snapshots
- Judgment summary: pinned=0, relevant=0, supports_claim=0, irrelevant=0, unreviewed=0
- Claim role coverage: support=5, foundational=1, contrastive=0
- Method spread: bayesian (4), multivariate (3), computational (2)
- Comparison/robustness leads: gvanvalkenhoef2015rsm_automated, consonni2008test_compatible, schuurman2016psychmethods_how
- Unreviewed high-signal hits: zwet2022tas_proposal, wang2007jasa_bayesian, anonndmisc_untitled, menard2011sf_standards
- Likely blind spots: top candidates are still unreviewed; foundational/reference-led coverage is thin

## 2026-04-26T10:29:03.110443+00:00 Note: Prior density transformation for formula scaling

Investigated original-scale prior density overlays for formula_scale regression terms. Current BayesTools path generates standardized prior samples, applies the existing coefficient unscale matrix, and plots a KDE of transformed samples. A deterministic alternative is available for scalar independent coefficient priors: for each original-scale coefficient row beta_j = sum_i M[j,i] theta_i, compute the marginal density as the convolution of scaled standardized prior marginals f_{M[j,i] theta_i}. This handles arbitrary priors with evaluable densities and ranges, including mixtures and point masses by decomposing into continuous and discrete components. Other packages mostly avoid a general arbitrary-prior density transform: rstanarm/brms use adjusted priors or transformed/prior-predictive draws; BayesFactor uses scale-adaptive g-priors; MCMCglmm gelman.prior solves only the multivariate-normal covariance transform case. Added local mockup script dev/prior-unscale-density-mockup.R for a two-predictor interaction case.

## 2026-04-30T05:31:31.121321+00:00 Search Query

Query: `bridge sampling marginal likelihood error estimate log marginal likelihood standard error Bayes factor relative error`

Task: `Decide representation of bridge-sampling Bayes factor numerical uncertainty`

Top results: overstall2010csda_default, anon2022bmcmrm_informed, fernandez2001jae_model, carter2018amapps_simple, scheibehenne2013pr_testing, gronau2017jmp_tutorial, oberauerndmisc_variance, wang2003biometrika_matching

## 2026-04-30T05:31:58.881615+00:00 Search Query

Query: `discover_claim_support(claim=Bridge-sampling numerical uncertainty for Bayes factors is more naturally represented as relative or percentage error derived from log marginal likelihood uncertainty than as lower and upper Bayes factor columns., seeds=boos1986biometrika_bootstrap,gronau2019psychometrika_simple)`

Task: `Decide representation of bridge-sampling Bayes factor numerical uncertainty`

Top results: gronau2021jss_informed, robert2014arsia_bayesian, kruschke2021nhb_bayesian, tendeiro2019psychmethods_review, amzal2006jasa_bayesian-optimal, wagenmakers2016cdps_bayesian, javiergiron2006sjs_objective, overstall2010csda_default

## 2026-04-30T05:31:59.041407+00:00 Search Query

Query: `bridge sampling marginal likelihood Bayes factor standard error log marginal likelihood relative error percentage Monte Carlo error`

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`

Top results: overstall2010csda_default, anon2022bmcmrm_informed, scheibehenne2013pr_testing, fernandez2001jae_model, robert2014arsia_bayesian, ardia2012csda_comparative, oberauerndmisc_variance, gronau2017jmp_tutorial

## 2026-04-30T05:32:21.885240+00:00 Search Query

Query: `discover_claim_support(claim=Monte Carlo uncertainty for Bayes factors is naturally propagated on the log Bayes-factor scale, and a one-column relative/percentage error is a reasonable scale-invariant display for multiplicative Bayes-factor uncertainty when the underlying standard error is for log probabilities or log marginal likelihoods., seeds=sellke2001tas_calibration,lee2003psychometrika_model)`

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`

Top results: johnstone2009aoas_approximate, gronau2020tas_informed, sellke2001tas_calibration, gronau2021jss_informed, anonndmisc_1023072291091, lindley1977biometrika_problem, butler1993as_asymptotics, shmueli2005jrssscas_useful

## 2026-04-30T05:32:29.231583+00:00 Search Query

Query: `discover_related_work(seeds=overstall2010csda_default,anon2022bmcmrm_informed,scheibehenne2013pr_testing, modes=references,cited_by,same_author,same_journal,same_method)`

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`

Top results: anonndmisc_1023072291091, wagenmakers2010cp_bayesian, smith1980jrsssbsm_bayes

## 2026-04-30T05:32:57.600327+00:00 Search Query

Query: `bridge sampling marginal likelihood error estimate standard error log marginal likelihood Bayes factor`

Task: `Decide representation of bridge-sampling Bayes factor numerical uncertainty`

Top results: overstall2010csda_default, anon2022bmcmrm_informed, fernandez2001jae_model, scheibehenne2013pr_testing, carter2018amapps_simple, vanpaemel2012pbr_using, gronau2017jmp_tutorial, wang2003biometrika_matching

## 2026-04-30T05:33:17.957373+00:00 Search Query

Query: `Meng Wong bridge sampling warp bridge sampling relative error normalizing constants ratio`

Task: `Decide representation of bridge-sampling Bayes factor numerical uncertainty`

Top results: sinharay2005jcgs_empirical, robert2009acp_computational, xu2011aoas_bayesian, llorente2023sr_marginal, ardia2012csda_comparative, fan2006jcgs_output, lin2010jasa_generating, roverato2002sjs_hyper

## 2026-04-30T05:33:39.456017+00:00 Search Query

Query: `Gronau Singmann Wagenmakers bridgesampling marginal likelihood error Bayes factor repeated runs percentage error`

Task: `Decide representation of bridge-sampling Bayes factor numerical uncertainty`

Top results: jasra2011sjs_inference, held2017biometrics_adaptive, tarantola2004statmod_mcmc, casella1987as_conditionally, dalalyan2008ml_aggregation, masson2011brm_tutorial, gramacy2010ba_shrinkage, atkinson2006jcgs_distribution

## 2026-04-30T05:33:41.262943+00:00 Workspace Materialization

- `overstall2010csda_default`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `anon2022bmcmrm_informed`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `scheibehenne2013pr_testing`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `fernandez2001jae_model`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `sellke2001tas_calibration`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `johnstone2009aoas_approximate`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `gronau2020tas_informed`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `anonndmisc_1023072291091`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `wagenmakers2010cp_bayesian`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `smith1980jrsssbsm_bayes`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-04-30T05:33:41.412456+00:00 Writing Context: Bayes factor Monte Carlo uncertainty representation in BayesTools tables

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`
Query: `bridge sampling marginal likelihood Bayes factor standard error log marginal likelihood relative error percentage Monte Carlo error`
Claim: `Monte Carlo uncertainty for Bayes factors is naturally propagated on the log Bayes-factor scale, and a one-column relative/percentage error is a reasonable scale-invariant display for multiplicative Bayes-factor uncertainty when the underlying standard error is for log probabilities or log marginal likelihoods.`
Claim ID: `claim_6bf76f2c5fca`
Prepared candidates: 10/10
Ready keys: anon2022bmcmrm_informed, anonndmisc_1023072291091, fernandez2001jae_model, gronau2020tas_informed, johnstone2009aoas_approximate, overstall2010csda_default, scheibehenne2013pr_testing, sellke2001tas_calibration, smith1980jrsssbsm_bayes, wagenmakers2010cp_bayesian
PDF-only: none
Human action: none

## 2026-04-30T05:34:06.993659+00:00 Search Query

Query: `bridge sampling marginal likelihood error estimate percentage error log marginal likelihood Bayes factor Gronau bridge_sampler error_measures`

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`

Top results: overstall2010csda_default, anon2022bmcmrm_informed, fernandez2001jae_model, oberauerndmisc_variance2, gronau2017jmp_tutorial, wang2003biometrika_matching, ardia2012csda_comparative, sinharay2005jcgs_empirical

## 2026-04-30T05:34:32.975827+00:00 Search Query

Query: `Meng Wong 1996 simulating ratios normalizing constants relative mean squared error bridge sampling`

Task: `Foundational bridge sampling uncertainty formulas`

Top results: llorente2023sr_marginal, robert2009acp_computational, lartillot2006sb_computing, sinharay2005jcgs_empirical, fan2006jcgs_output, song2006sim_model, xu2011aoas_bayesian, gronau2021jss_informed

## 2026-04-30T05:34:53.211884+00:00 Search Query

Query: `bridgesampling R package error_measures percentage error logml Bayes factor`

Task: `Bridge sampling package uncertainty output`

Top results: wang2003biometrika_matching, gronau2017jmp_tutorial, vasishth2023cbb_sample, overstall2010csda_default, tsilidis2012jjnci_evaluation, ren2018mdm_incorporating, ripatti2000biometrics_estimation, collins2016sim_sample

## 2026-04-30T05:36:26.816866+00:00 Workspace Materialization

- `gronau2017jmp_tutorial`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `sinharay2005jcgs_empirical`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `llorente2023sr_marginal`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `robert2009acp_computational`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `overstall2010csda_default`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `gelman1998statscience_simulating`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `anonndmisc_1023072291091`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-04-30T05:36:49.420604+00:00 Search Query

Query: `Meng Wong 1996 simulating ratios of normalizing constants optimal bridge sampling relative mean square error coefficient of variation`

Task: `Bayes factor Monte Carlo uncertainty representation in BayesTools tables`

Top results: robert2009acp_computational, llorente2023sr_marginal, gronau2021jss_informed, sinharay2005jcgs_empirical, overstall2010csda_default, lin2010jasa_generating, marin2010ejs_resolving, thompson1991biometrika_noise

## 2026-04-30T05:36:58.316083+00:00 Workspace Materialization

- `meng1996sinica_simulating`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none
- `fruhwirthschnatter2004ej_estimating`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=n/a, actions=none

## 2026-04-30T05:38:46.602924+00:00 Note: Bridge-sampling BF uncertainty representation

Investigated whether bridge-sampling numerical uncertainty for Bayes factors should be represented as relative/percentage error rather than lower/upper BF columns. Materialized bridge-sampling and marginal-likelihood sources: Meng & Wong (1996), Frühwirth-Schnatter (2004), Gronau et al. (2017), Sinharay & Stern (2005), Gelman & Meng (1998), Overstall & Forster (2010), Robert et al. (2009), Llorente et al. (2023), Kass & Raftery (1995). Core rationale: bridge-sampling theory defines estimator accuracy through asymptotic relative mean-squared error / coefficient of variation for ratios of normalizing constants; on the log scale, BF uncertainty is additive and can be propagated from log marginal likelihood errors, with relative BF error approximately sqrt(se_logml_1^2 + se_logml_0^2) when independent. BF-scale intervals are transforms of this quantity and can be skewed or visually over-interpreted as inferential intervals.

## 2026-04-30T05:45:50.462823+00:00 Coverage Report: Bayes factor Monte Carlo uncertainty representation in BayesTools tables

Coverage snapshot:
- Task sources: 10 (citation-ready: 10)
- Discovery candidates considered: 27 from 3 recent query snapshots
- Judgment summary: pinned=0, relevant=0, supports_claim=0, irrelevant=0, unreviewed=10
- Claim role coverage: support=8, foundational=1, contrastive=1
- Method spread: bayesian (14), model_selection (10), hypothesis_testing (6)
- Author spread: Bartoš (1), Kass (1), Fernández (1)
- Relation coverage: references=4, same_method=3, same_author=2, cited_by=1
- Foundational follow-ups still missing or not ready: fernandez2001je_benchmark
- Comparison/robustness leads: gronau2017jmp_tutorial, smith1980jrsssbsm_bayes, sinharay2005jcgs_empirical
- Unreviewed high-signal hits: overstall2010csda_default, anon2022bmcmrm_informed, scheibehenne2013pr_testing, fernandez2001jae_model
- Likely blind spots: top candidates are still unreviewed

## 2026-04-30T05:46:03.206854+00:00 Note: BF Monte Carlo error display in BayesTools tables

For BayesTools BF diagnostics, reviewed bridge-sampling and marginal-likelihood error literature. The relevant design point is that numerical uncertainty for marginal likelihood ratios is naturally handled on the log scale, while bridge-sampling accuracy is commonly summarized by relative mean-squared error/coefficient of variation. For inclusion BFs computed from model indicator frequencies, the implementation analog is the delta-method SE(log BF) = SE(p) / (p(1-p)); reporting 100 * SE(log BF) as a relative BF MC error percentage gives a one-column, scale-invariant diagnostic. Caveat: this is numerical/Monte Carlo error only and does not address MCMC non-convergence or poor bridge/proposal overlap.

## 2026-05-04T20:02:12.266028+00:00 Search Query

Query: `product space Bayes factor posterior model probability indicator MCMC all samples zero one`

Task: `Find product-space BF finite-bound paper`

Top results: scheibehenne2013pr_testing, sebastiani2001jasa_bayesian, gonen2005tas_bayesian, liu2001jasa_theory, malikov2015e_bayesian, guo2010biostatistics_multiplicity-calibrated, chen2007jasa_sieve, lodewyckx2011jmp_tutorial

## 2026-05-04T20:04:56.183185+00:00 Workspace Materialization

- `lodewyckx2011jmp_tutorial`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-05-04T20:06:00.931595+00:00 Search Query

Query: `posterior model probabilities zero one Bayes factor MCMC model index`

Task: `Find boundary handling for product-space Bayes factor estimates`

Top results: denham2003biostatistics_bayesian, tendeiro2019psychmethods_review, sabanesbove2011statcomp_bayesian, villa2017cstm_mathematics, held2010bmcmrm_nomogram, lyne2015statscience_russian, wagenmakers2022ahes_history, berkhout2023brm_tutorial

## 2026-05-04T20:07:44.761001+00:00 Search Query

Query: `Han Carlin 2001 Markov Chain Monte Carlo Methods Computing Bayes Factors product space model indicator`

Task: `Find Han Carlin product-space Bayes factor all samples boundary discussion`

Top results: scheibehenne2013pr_testing, ramoni2002ml_bayesian, he2019psychometrika_optimal, brown2003pre_statistical, bornn2011ba_bayesian, schorfheide2011wpfrbp_estimation, lodewyckx2011jmp_tutorial, vanravenzwaaij2018pbr_simple

## 2026-05-04T20:10:40.266995+00:00 Workspace Materialization

- `han2001jasa_markov`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-05-04T20:18:06.578581+00:00 Note: Product-space BF boundary/all-indicator-samples search

Searched for a remembered paper on minimum/maximum/adjusted Bayes factors from product-space/model-indicator estimation, especially all-0/all-1 indicator samples. Materialized `lodewyckx2011jmp_tutorial` and `han2001jasa_markov`; `scheibehenne2013pr_testing`, `sinharay2005jcgs_empirical`, `wagenmakers2010cp_bayesian`, `sellke2001tas_calibration`, and Kass/Raftery (`anonndmisc_1023072291091`) were already available. Strongest target is Lodewyckx et al. (2011), which defines posterior model probability as model-index visit frequency, explains that with strong evidence the disfavored model is almost never activated, and recommends bisection/calibration of prior model probabilities plus transition-matrix monitoring. Han & Carlin (2001) provide the foundational comparative review and explicitly recommend choosing prior model probabilities to balance model visit frequencies. I did not find a StatsVault-ready source recommending an add-one/pseudocount/continuity-corrected finite BF when all indicator samples are 0 or 1. Searches for minimum/maximum Bayes factor mainly led to p-value/prior-class bounds (Sellke et al.; Kass & Raftery/Edwards-Lindman-Savage), not product-space count estimation. Attempted `report_task_coverage` twice; it timed out/errored, so coverage was assessed from direct searches and materialized summaries.

## 2026-05-04T20:52:43.005100+00:00 Search Query

Query: `Quantifying uncertainty transdimensional Markov chain Monte Carlo discrete Markov models Bayes factor product space MCMC model indicator Heck Overstall Gronau Wagenmakers`

Task: `general`

Top results: rueda2007pcb_flexible, chay1972jmva_quasi-markov, lee1968jasa_maximum, heck2019statcomp_quantifying, bolt2003apm_estimation, lodewyckx2011jmp_tutorial, franck2020tas_assessing, andrieu2003ml_introduction

## 2026-05-04T20:53:40.718688+00:00 Search Query

Query: `Testing order constraints qualitative differences Bayes factors normalized maximum likelihood maximum possible Bayes factor Heck Wagenmakers Morey binomial order constrained`

Task: `general`

Top results: dunson2003biometrics_bayesian, good1967jrsssbm_bayesian, dunson2008biometrika_bayesian, tokdar2021sb_bayesian, rojo1992as_pure-tail, spjotvoll1977biometrika_ordering, heck2015spl_testing, mulder2013bjmasp_bayes

## 2026-05-04T20:55:37.203793+00:00 Search Query

Query: `Approximated adjusted fractional Bayes factors general method testing informative hypotheses Gu Mulder Hoijtink adjusted fractional Bayes factor`

Task: `general`

Top results: young1996jrsssbsm_measuring, carvalho2009biometrika_objective, bertolino2000jrsssds_bayesian, moreno1999jspi_default, johnson2013pnas_revised, gu2018bjmsp_approximated, mulder2014csada_prior, vanlissandmisc_teachers

## 2026-05-04T20:57:33.270336+00:00 Search Query

Query: `Bayesian model based sampling posterior model probabilities Horvitz Thompson Clyde model uncertainty`

Task: `general`

Top results: anonndmisc_untitled, anonndmisc_untitled2, bourne2020jgrse_stressdependent, clyde2004statscience_model, dellaportas2002statcomp_bayesian, karabatsos2006jmp_bayesian, gronau2018psychometrika_simple, hennig2015prsmpaes_probabilistic

## 2026-05-04T21:01:57.180897+00:00 Coverage Report: Identify conference-linked non-Marsman/Sekulovski sources for minimum/maximum/adjusted Bayes factor computations from product-space or Bernoulli indicator samples, especially all-0/all-1 cases.

Coverage snapshot:
- Task sources: 0 (citation-ready: 0)
- Discovery candidates considered: 44 from 3 recent query snapshots
- Judgment summary: pinned=0, relevant=0, supports_claim=0, irrelevant=0, unreviewed=0
- Claim role coverage: support=0, foundational=0, contrastive=0
- Method spread: bayesian (11), latent_variable (9), computational (8)
- Relation coverage: same_method=8, references=6, cited_by=2, same_author=2
- Foundational follow-ups still missing or not ready: nobile2004as_posterior, rousseau2011jrsssbsm_asymptotic, diebolt1994jrsssbsm_estimation, ishwaran2001jasa_bayesian
- Comparison/robustness leads: donvandenbergh2025tas_flexible, rosenbaum2018aoas_sensitivity, fruhwirth-schnatter2019adaac_from
- Unreviewed high-signal hits: anon2009e_inference, fiorentini2004e_likelihood-based, malsiner-walli2014statcomp_model-based, adamczak2009jams_quantitative
- Likely blind spots: top candidates are still unreviewed; no clear foundational claim support surfaced yet

## 2026-05-05T06:56:11.703631+00:00 Search Query

Query: `product-space model-index MCMC inclusion Bayes factor posterior inclusion probability all samples zero one uncertainty lower bound smoothing encompassing prior informative hypotheses`

Task: `Assess reporting conventions for inclusion Bayes factors from product-space/model-index MCMC when all indicator samples are 0 or 1`

Top results: sekulovski2024mbr_testing, malikov2015e_bayesian, malsiner-walli2016ajs_comparing, koslovsky2020bb_microbvs, gewekendmisc_variable, lykou2013statcomp_bayesian, jones2005statscience_experiments, guedj2013ejs_pac-bayesian

## 2026-05-05T06:59:03.214470+00:00 Search Query

Query: `Heck Davis-Stober multinomineq encompassing Bayes factors beta Jeffreys smoothing zero counts all constraints posterior samples`

Task: `Assess reporting conventions for inclusion Bayes factors from product-space/model-index MCMC when all indicator samples are 0 or 1`

Top results: hooper2008ba_exact, vandewiel2013biostatistics_bayesian, lockwood2015aoas_inferring, johnson2013pnas_revised, kolaczyk2004as_multiscale, mohan2021jasa_graphical, kaciroti2012biostatistics_bayesian, gu2014psychmethods_bayesian

## 2026-05-05T07:00:01.734180+00:00 Search Query

Query: `Klugkist Laudy Hoijtink inequality constrained analysis of variance Bayesian encompassing prior posterior prior proportion all samples satisfy constraint`

Task: `Assess reporting conventions for inclusion Bayes factors from product-space/model-index MCMC when all indicator samples are 0 or 1`

Top results: mulder2016jmp_bayes, brunk1958ams_estimation, heck2016jmp_adjusted, laudy2007smmr_bayesian, klugkist2007csada_bayes, berger1990jspi_robust, wetzels2010csda_encompassing, vanwesel2013qq_use

## 2026-05-05T07:13:34.129277+00:00 Note: Inclusion BF boundary review progress - started

Started review of local summaries for product-space/model-index inclusion Bayes factor reporting when the sampled inclusion indicator is all 0 or all 1. Delegated four focused evidence extractions: Heck et al. (2019) / MCMCprecision; Gámbaro et al. (2025) and BEAST/BSSVS adjusted BF convention; Heck & Davis-Stober / multinomineq smoothing; and Klugkist/Gu/Mulder/Hoijtink encompassing-prior informative-hypothesis BFs. Local critical-path work will inspect BayesTools inclusion-BF reporting and synthesize whether to use diagnostics, bounds, smoothing intervals, or changes to finite BF display.

## 2026-05-05T07:14:02.994471+00:00 Search Query

Query: `product-space model-index MCMC inclusion Bayes factor posterior inclusion probability all samples zero one uncertainty lower bound smoothing encompassing prior informative hypotheses`

Task: `Assess reporting conventions for inclusion Bayes factors from product-space/model-index MCMC when all indicator samples are 0 or 1`

Top results: sekulovski2024mbr_testing, malikov2015e_bayesian, malsiner-walli2016ajs_comparing, koslovsky2020bb_microbvs, gewekendmisc_variable, lykou2013statcomp_bayesian, jones2005statscience_experiments, guedj2013ejs_pac-bayesian

## 2026-05-05T07:22:47.430114+00:00 Search Query

Query: `Daniel Heck`

Task: `Find Heck Davis-Stober multinomineq`

Top results: gronaundmisc_bayesian, vandoornndmisc_bayes, heck2018jss_rrreg, kloft2025eapm_discriminant, gronau2017crsp_bayesian, heck2016jmp_adjusted, gronau2021amapps_primer, heck2015spl_testing

## 2026-05-05T07:23:41.588433+00:00 Workspace Materialization

- `klugkist2007csada_bayes`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `laudy2007smmr_bayesian`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `gu2014psychmethods_bayesian`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `gu2018bjmsp_approximated`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `mulder2014csada_prior`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `mulder2016jmp_bayes`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `heck2016jmp_adjusted`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none
- `wetzels2010csda_encompassing`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-05-05T07:24:08.425162+00:00 Workspace Materialization

- `heck2019jmp_multinomial`: readiness=ready_canonical, pdf=copied, summary=copied, local_summary=not_created, bib=written, actions=none

## 2026-05-05T07:24:35.079372+00:00 Note: Inclusion BF boundary reporting from product-space samples

Reviewed local summaries and subagent notes for reporting inclusion Bayes factors from product-space/model-index MCMC when the sampled inclusion indicator is all 0 or all 1. Current BayesTools behavior: `inclusion_BF()` returns the exact odds-ratio boundary (`0` or `Inf`) from posterior model probabilities; `runjags_inference_table(..., BF_diagnostics=TRUE)` computes delta-method `error%(BF)` for interior probabilities and returns `NA` plus a warning when one side has zero posterior visits.

Synthesis:

1. Heck et al. (2019) / MCMCprecision supports diagnostics, not finite boundary correction. Its Markov-model method estimates uncertainty of posterior model probabilities from the full model-index chain by sampling transition matrices and stationary distributions. For inclusion BFs, the principled extension is to fit the full model-index transition matrix, then compute inclusion subset probabilities from stationary-distribution draws. It explicitly does not justify adding pseudocounts to an aggregated all-0/all-1 inclusion bit; never/rarely visited models make BFs unreliable and indicate insufficient product-space mixing or resolution.

2. BEAST/BSSVS and Gámbaro et al. (2025) support an inequality-style reporting convention for saturated posterior inclusion probabilities. When posterior inclusion `p = 1`, replace only for display with the closest resolvable value `(N - 1) / N`, yielding a lower bound `BF > (N - 1) / prior_odds`, where `prior_odds = q / (1 - q)`. A symmetric BayesTools extension for all-0 indicators would use `p = 1 / N`, yielding `BF < 1 / ((N - 1) * prior_odds)`. This should be labelled as a Monte Carlo resolution bound, not a corrected/unbiased BF and not `BFadj`; Gámbaro's BFadj is a domain-specific empirical-prior adjustment for phylogeographic sampling bias.

3. Heck & Davis-Stober / multinomineq uses the encompassing-prior BF framework for inequality-constrained multinomial models. Its stepwise/automatic algorithms use beta uncertainty for estimated constraint-hit proportions and ensure positive hits for components where ratios are needed. This can inspire optional MC interval/bound displays for count-based proportions, but not a finite product-space point BF. If used, smoothing should be labelled as uncertainty regularization and should preferably use ESS-adjusted effective counts because product-space indicators are autocorrelated.

4. Klugkist/Laudy/Hoijtink/Gu informative-hypothesis sources show why an encompassing-prior BF remains finite when all posterior draws satisfy a constraint: `BF = posterior_constraint_mass / prior_constraint_mass`, so plug-in posterior mass 1 gives `1 / prior_constraint_mass`. This finite result is a property of the constrained-vs-encompassing estimand, not a correction for product-space model odds. It should be borrowed for explanation only: product-space all-1/all-0 indicators mean the complement model probability is unresolved by the chain.

Recommended BayesTools reporting direction:

- Keep the internal numeric BF estimator exact: `Inf` for all included and `0` for all excluded.
- Improve table/reporting display by optionally replacing boundary displays with one-sided resolution bounds: `> bound` for all-1 inclusion BFs and `< bound` for all-0 inclusion BFs, with formulas above and with `N` or preferably indicator ESS stated.
- Keep `error%(BF) = NA` at boundaries, but make the warning/footnote visible whenever a boundary BF is printed, not only when BF diagnostics are requested.
- Do not silently use Jeffreys/add-one smoothing as the main BF. If an interval is added, present it as a Monte Carlo uncertainty interval over the indicator probability transformed to BF scale, not as the estimate.
- Longer-term: add a full model-index diagnostic path following MCMCprecision for product-space model indicators, then derive inclusion-BF uncertainty from stationary model-probability draws. This is the principled upgrade for autocorrelated product-space chains.

## 2026-05-05T07:45:03.087888+00:00 Note: Implemented finite-sample bounds for product-space inclusion BF table reporting

Implemented BEAST/BSSVS-style reporting convention for product-space inclusion Bayes factors in `R/summary-tables.R`: when posterior inclusion probability is exactly 0 or 1 and prior inclusion probability is in (0, 1), table values now store/display finite-sample lower/upper bounds using 1/S or (S-1)/S, with '<'/'>' operator metadata and a table footnote. Bound metadata is preserved through `format_BF()` and `update.BayesTools_table()`, including inequality inversion for BF01. Updated BF Monte Carlo error column label to `error%(Inclusion BF)` and adjusted summary-table tests/reference outputs. Verification: unit profile passed; targeted summary/JAGS table fixture tests passed; full fixture profile stopped only in fixture-integrity due stale cached fit metadata.

## 2026-05-05T08:26:51.032495+00:00 Note: Adjusted inclusion BF bound reporting and interpretation

Follow-up implementation for product-space inclusion BF boundary reporting: removed the automatic finite-sample-bound explanatory footnote from default `runjags_inference_table()` output and moved the explanation into documentation (`R/summary-tables.R` roxygen and `man/BayesTools_model_tables.Rd`). Fixed table printing to format columns from the original vector so `<`/`>` bound metadata is visible in printed summaries. Added a small `BayesTools_BF` S3 vector subset method so single-row extraction preserves `bound_operator` metadata. Extended `.interpret.BF()`, `interpret()`, and `interpret2()` to report finite-sample inequalities when a BF carries `bound_operator` metadata or an explicit `BF_bound_operator` / `inference_BF_bound_operator` field. Added tests for printed signs, absence of default footnote, scalar extraction metadata, and interpretation output. Verification: unit profile passed; targeted fixture tests for summary/JAGS tables passed; full fixture profile still fails only in fixture-integrity/catalog metadata expectations unrelated to this change.

## 2026-05-05T09:23:37.267692+00:00 Note: Added records-first structured interpretation normalizer

Implemented a new records-first interpretation API in `R/interpret.R`: exported `interpret_records()` and `interpret_tables()`. The API accepts named sources (BayesTools tables, plain data frames with per-source schemas, direct record sources, or record collections) plus an explicit ordered plan. It emits traceable `BayesTools_interpret_records` data frames containing header/note/prior/footnote/evidence/estimate records with section, item_id, source, row, order, record_id, central estimate metadata, interval metadata, free-text conditioning, raw BF metadata (scale, orientation, name, bound), and canonical target-over-complement BF metadata for safe interpretation. Added dynamic `for_each` plan expansion for row-driven moderator/scale/interaction sections and a generic text renderer for simple use. Added `interpret_tables()` as an alias/wrapper for table-oriented callers. Added tests covering mixed table/direct sources, log-BF01 bounded inclusion BF normalization, explicit order, dynamic row generation, summary estimate records, generic text, alias behavior, and optional missing sources. Verification: unit profile passed; targeted fixture tests passed; full fixture profile still fails only in existing fixture-integrity catalog/cache mismatch.
