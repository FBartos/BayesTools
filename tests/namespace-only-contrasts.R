# This script is intentionally outside tests/testthat: R CMD check runs it in a
# fresh process where BayesTools is installed but not attached.

make_prior <- function(contrast){

  switch(
    contrast,
    treatment = BayesTools::prior_factor(
      "normal", list(0, 1), contrast = "treatment"
    ),
    independent = BayesTools::prior_factor(
      "normal", list(0, 1), contrast = "independent"
    ),
    orthonormal = BayesTools::prior_factor(
      "mnormal", list(0, 1), contrast = "orthonormal"
    ),
    meandif = BayesTools::prior_factor(
      "mnormal", list(0, 1), contrast = "meandif"
    ),
    cumulative = BayesTools::prior_ordered(
      BayesTools::prior("normal", list(0, 1)),
      contrast = "cumulative"
    ),
    cumulative_levels = BayesTools::prior_ordered(
      BayesTools::prior("normal", list(0, 1)),
      contrast = "cumulative_levels"
    )
  )
}

contrast_names <- c(
  "treatment",
  "independent",
  "orthonormal",
  "meandif",
  "cumulative",
  "cumulative_levels"
)

for(contrast_name in contrast_names){
  factor_data <- factor(
    rep(c("a", "b", "c"), 2),
    levels = c("a", "b", "c"),
    ordered = grepl("cumulative", contrast_name, fixed = TRUE)
  )
  result <- BayesTools::JAGS_formula(
    formula   = ~ group,
    parameter = "mu",
    data      = data.frame(group = factor_data),
    prior_list = list(
      intercept = BayesTools::prior("normal", list(0, 1)),
      group     = make_prior(contrast_name)
    )
  )

  stopifnot(
    inherits(result$formula_design, "BayesTools_formula_design"),
    nrow(result$formula_design$model_matrix) == length(factor_data)
  )
}

random_data <- data.frame(
  group = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
  id    = factor(rep(c("g1", "g2", "g3"), each = 2))
)
random_result <- suppressWarnings(BayesTools::JAGS_formula(
  formula   = ~ 1 + diag(0 + group | id),
  parameter = "mu",
  data      = random_data,
  prior_list = list(
    intercept = BayesTools::prior("normal", list(0, 1))
  ),
  prior_random = BayesTools::prior_random(
    id = BayesTools::random_block(
      sd = BayesTools::prior_factor(
        "mnormal", list(0, 1), contrast = "orthonormal"
      )
    )
  )
))
random_term <- random_result$formula_design$random_effects[[1L]]
random_prediction <- BayesTools:::.bt_random_effect_prediction_data(
  random_term,
  data = random_data
)

stopifnot(
  identical(random_term$contrasts$group, "contr.orthonormal"),
  identical(random_prediction$model_matrix, random_term$model_matrix)
)
