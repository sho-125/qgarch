\name{qgarch_default_starts}
\alias{qgarch_default_starts}
\title{Default Starting Values for QGARCH Estimation}
\description{Constructs a small set of default starting vectors for \code{qgarch_fit()}.}
\usage{
qgarch_default_starts(x, model = c("zero", "restricted", "free", "threshold"))
}
\arguments{
  \item{x}{A numeric vector of returns.}
  \item{model}{Model variant.}
}
\value{A list of named numeric vectors.}
\examples{
qgarch_default_starts(rnorm(100), model = "restricted")
}
