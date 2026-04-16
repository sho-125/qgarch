\name{qgarch_fit}
\alias{qgarch_fit}
\title{Fit a QGARCH-in-Mean Model}
\description{
Fits a QGARCH-M(1,1) specification using the likelihood recursions from the
supplied replication code.
}
\usage{
qgarch_fit(
  x,
  model = c("zero", "restricted", "free", "threshold"),
  threshold_indicator = NULL,
  start = NULL,
  rho = 1,
  steptol = 1e-10,
  typsize = 0.1,
  print.level = 0,
  hessian = TRUE
)
}
\arguments{
  \item{x}{A numeric vector of returns.}
  \item{model}{One of \code{"zero"}, \code{"restricted"}, \code{"free"}, or \code{"threshold"}.}
  \item{threshold_indicator}{A numeric or logical indicator used only when \code{model = "threshold"}.}
  \item{start}{Optional starting values. Supply either a numeric vector or a list of numeric vectors.}
  \item{rho}{Restriction parameter used in the mapping from \code{(alpha, beta, gamma)} to \code{lambda} for the restricted model. The default \code{rho = 1} matches the supplied replication scripts.}
  \item{steptol}{Passed to \code{stats::nlm()}.}
  \item{typsize}{Passed to \code{stats::nlm()}.}
  \item{print.level}{Passed to \code{stats::nlm()}.}
  \item{hessian}{Logical; should the Hessian be computed?}
}
\value{
An object of class \code{"qgarch"}.
}
\details{
The package implements the QGARCH-M(1,1) family used in the uploaded
replication scripts. The threshold specification follows the supplied code
closely and uses a constant \code{kappa} based on the sample mean of the
indicator.
}
\examples{
set.seed(123)
x <- rnorm(300, sd = 0.03)
fit <- qgarch_fit(x, model = "free")
fit
summary(fit)
}
