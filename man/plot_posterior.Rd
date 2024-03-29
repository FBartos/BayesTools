% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-averaging-plots.R
\name{plot_posterior}
\alias{plot_posterior}
\title{Plot samples from the mixed posterior distributions}
\usage{
plot_posterior(
  samples,
  parameter,
  plot_type = "base",
  prior = FALSE,
  n_points = 1000,
  n_samples = 10000,
  force_samples = FALSE,
  transformation = NULL,
  transformation_arguments = NULL,
  transformation_settings = FALSE,
  rescale_x = FALSE,
  par_name = NULL,
  dots_prior = list(),
  ...
)
}
\arguments{
\item{samples}{samples from a posterior distribution for a
parameter generated by \link{mix_posteriors}.}

\item{parameter}{parameter name to be plotted. Use \code{"PETPEESE"}
for PET-PEESE plot with parameters \code{"PET"} and \code{"PEESE"},
and \code{"weightfunction"} for plotting a weightfunction with
parameters \code{"omega"}.}

\item{plot_type}{whether to use a base plot \code{"base"}
or ggplot2 \code{"ggplot"} for plotting.}

\item{prior}{whether prior distribution should be added to the figure}

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

\item{dots_prior}{additional arguments for the prior distribution plot}

\item{...}{additional arguments}
}
\value{
\code{plot_posterior} returns either \code{NULL} or
an object of class 'ggplot' if plot_type is \code{plot_type = "ggplot"}.
}
\description{
Plot samples from the mixed posterior distributions
}
\seealso{
\code{\link[=prior]{prior()}} \code{\link[=lines_prior_list]{lines_prior_list()}}  \code{\link[=geom_prior_list]{geom_prior_list()}}
}
