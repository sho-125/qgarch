\name{qgarch_utils}
\alias{qgarch_moments}
\alias{qgarch_lr_test}
\title{Utility Functions for QGARCH Workflows}
\description{Moment diagnostics and likelihood-ratio testing helpers.}
\usage{
qgarch_moments(x, scale_mean_variance = 1)

qgarch_lr_test(restricted, unrestricted, df = NULL)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{scale_mean_variance}{Optional scaling applied to the mean and variance estimates.}
  \item{restricted}{A fitted \code{"qgarch"} object representing the restricted model.}
  \item{unrestricted}{A fitted \code{"qgarch"} object representing the unrestricted model.}
  \item{df}{Optional degrees of freedom for the chi-squared test.}
}
\value{
\code{qgarch_moments()} returns a list with estimates and standard errors.
\code{qgarch_lr_test()} returns a data frame with the test statistic,
degrees of freedom, and p-value.
}
\examples{
qgarch_moments(rnorm(200))
}
