nested_weight_priors <- function(alpha = rep(1, 4), fine_steps = c(.025, .05, .50)){

  fine <- prior_weightfunction("one-sided", fine_steps, wf_cumulative(alpha), prior_weights = 2)
  coarse <- prior_weightfunction("one-sided", c(.05, .50), wf_cumulative(rep(1, 3)))
  other <- prior_weightfunction("one-sided", .30, wf_fixed(c(1, .35)))
  # Put the finer component first and give an unrelated branch an extra cut.
  list(bias = prior_mixture(list(fine, prior_none(), other, coarse)))
}
