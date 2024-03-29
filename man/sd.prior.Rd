% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{sd.prior}
\alias{sd.prior}
\title{Prior sd}
\usage{
\method{sd}{prior}(x, ...)
}
\arguments{
\item{x}{a prior}

\item{...}{unused arguments}
}
\value{
a standard deviation of an object of class 'prior'.
}
\description{
Computes standard deviation
of a prior distribution.
}
\examples{
# create a standard normal prior distribution
p1 <- prior(distribution = "normal", parameters = list(mean = 1, sd = 1))

# compute sd of the prior distribution
sd(p1)

}
\seealso{
\code{\link[=prior]{prior()}}
}
