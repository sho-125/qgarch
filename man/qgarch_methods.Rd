\name{qgarch_methods}
\alias{print.qgarch}
\alias{summary.qgarch}
\alias{print.summary.qgarch}
\alias{coef.qgarch}
\alias{vcov.qgarch}
\alias{fitted.qgarch}
\alias{residuals.qgarch}
\alias{logLik.qgarch}
\alias{predict.qgarch}
\alias{plot.qgarch}
\title{Methods for \code{"qgarch"} Objects}
\description{S3 methods for fitted QGARCH objects.}
\usage{
\method{print}{qgarch}(x, digits = max(3L, getOption("digits") - 2L), ...)

\method{summary}{qgarch}(object, ...)

\method{print}{summary.qgarch}(x, digits = max(3L, getOption("digits") - 2L), ...)

\method{coef}{qgarch}(object, ...)

\method{vcov}{qgarch}(object, ...)

\method{fitted}{qgarch}(object, ...)

\method{residuals}{qgarch}(object, type = c("raw", "standardized", "eta"), ...)

\method{logLik}{qgarch}(object, ...)

\method{predict}{qgarch}(object, n.ahead = 1L, ...)

\method{plot}{qgarch}(x, which = c("sigma2", "standardized"), ...)
}
\arguments{
  \item{object, x}{A fitted \code{"qgarch"} object.}
  \item{digits}{Number of digits to print.}
  \item{type}{Residual type.}
  \item{n.ahead}{Forecast horizon.}
  \item{which}{Plot type.}
  \item{...}{Unused additional arguments passed through the generic.}
}
\value{
Standard S3 outputs. In particular, \code{predict()} returns a data frame with
forecasted conditional means and variances.
}
