% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-averaging.R
\name{mix_posteriors}
\alias{mix_posteriors}
\title{Model-average posterior distributions}
\usage{
mix_posteriors(
  model_list,
  parameters,
  is_null_list,
  conditional = FALSE,
  seed = NULL,
  n_samples = 10000
)
}
\arguments{
\item{model_list}{list of models, each of which contains marginal
likelihood estimated with bridge sampling \code{marglik} and prior model
odds \code{prior_weights}}

\item{parameters}{vector of parameters names for which inference should
be drawn}

\item{is_null_list}{list with entries for each parameter carrying either
logical vector of indicators specifying whether the model corresponds
to the null or alternative hypothesis (or an integer vector indexing models
corresponding to the null hypothesis)}

\item{conditional}{whether prior and posterior model probabilities should
be returned only for the conditional model. Defaults to \code{FALSE}}

\item{seed}{integer specifying seed for sampling posteriors for
model averaging. Defaults to \code{NULL}.}

\item{n_samples}{number of samples to be drawn for the model-averaged
posterior distribution}
}
\value{
\code{mix_posteriors} returns a named list of mixed posterior
distributions (either a vector of matrix).
}
\description{
Model-averages posterior distributions based on
a list of models, vector of parameters, and a list of
indicators of the null or alternative hypothesis models
for each parameter.
}
\seealso{
\link{ensemble_inference} \link{BayesTools_ensemble_tables}
}
