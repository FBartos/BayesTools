% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-averaging-plots.R
\name{plot_prior_list}
\alias{plot_prior_list}
\title{Plot a list of prior distributions}
\usage{
plot_prior_list(
  prior_list,
  plot_type = "base",
  x_seq = NULL,
  xlim = NULL,
  x_range_quant = NULL,
  n_points = 500,
  n_samples = 10000,
  force_samples = FALSE,
  transformation = NULL,
  transformation_arguments = NULL,
  transformation_settings = FALSE,
  rescale_x = FALSE,
  par_name = NULL,
  prior_list_mu = NULL,
  ...
)
}
\arguments{
\item{prior_list}{list of prior distributions}

\item{plot_type}{whether to use a base plot \code{"base"}
or ggplot2 \code{"ggplot"} for plotting.}

\item{x_seq}{sequence of x coordinates}

\item{xlim}{x plotting range}

\item{x_range_quant}{quantile used for
automatically obtaining \code{x_range}
if both \code{x_range} and \code{x_seq}
are unspecified. Defaults to \code{0.005}
for all but Cauchy, Student-t, Gamma, and
Inverse-gamme distributions that use
\code{0.010}.}

\item{n_points}{number of equally spaced points
in the \code{x_range} if \code{x_seq} is unspecified}

\item{n_samples}{number of samples from the prior
distribution if the density cannot be obtained
analytically (or if samples are forced with
\code{force_samples = TRUE})}

\item{force_samples}{should prior be sampled instead
of obtaining analytic solution whenever possible}

\item{transformation}{transformation to be applied
to the prior distribution. Either a character
specifying one of the prepared transformations:
\describe{
\item{lin}{linear transformation in form of \code{a + b*x}}
\item{tanh}{also known as Fisher's z transformation}
\item{exp}{exponential transformation}
}, or a list containing the transformation function \code{fun},
inverse transformation function \code{inv}, and the Jacobian of
the transformation \code{jac}. See examples for details.}

\item{transformation_arguments}{a list with named arguments for
the \code{transformation}}

\item{transformation_settings}{boolean indicating whether the
settings the \code{x_seq} or \code{x_range} was specified on
the transformed support}

\item{rescale_x}{allows to rescale x-axis in case a
weightfunction is plotted.}

\item{par_name}{a type of parameter for which the prior is
specified. Only relevant if the prior corresponds to a mu
parameter that needs to be transformed.}

\item{prior_list_mu}{list of priors for the mu parameter
required when plotting PET-PEESE}

\item{...}{additional arguments}
}
\value{
\code{plot_prior_list} returns either \code{NULL} or
an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
}
\description{
Plot a list of prior distributions
}
\seealso{
\code{\link[=prior]{prior()}} \code{\link[=lines_prior_list]{lines_prior_list()}}  \code{\link[=geom_prior_list]{geom_prior_list()}}
}
