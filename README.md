# qgarch

`qgarch` is an R package for estimating and working with quadratic
GARCH-in-mean models in the Campbell and Hentschel (1992)
volatility-feedback framework. 

We replicated and extended the original paper in: 
Jedrzej Bialkowski, Sanghyun Hong and Moritz Wagner (2025).
"Is no news still good news? Volatility feedback revisited." Pacific-Basin Finance Journal, Volume 91.
https://doi.org/10.1016/j.pacfin.2025.102708.

In the Campbell and Hentschel model, expected returns vary with
conditional variance, and stock returns include a volatility-feedback
term that helps capture asymmetric return dynamics. The framework is
designed to accommodate features such as predictive asymmetry, negative
skewness, and excess kurtosis in stock returns.

## What the package does

`qgarch` provides tools to:

- fit generalized QGARCH-in-mean models with nonlinear maximum likelihood,
- work with several model specifications through a common interface,
- extract coefficients, fitted values, residuals, variance-covariance matrices,
  and log-likelihood values,
- generate forecasts of the conditional mean and conditional variance,
- compute moment diagnostics, and
- compare nested specifications with likelihood-ratio tests.

The main estimation function is `qgarch_fit()`. It supports four model
specifications:

- `model = "zero"`: sets `lambda = 0`,
- `model = "restricted"`: links `lambda` to the remaining parameters through
  the restricted mapping,
- `model = "free"`: estimates `lambda` freely, and
- `model = "threshold"`: allows a state-dependent `lambda_t`.

## Installation

You can install the development version of `qgarch` directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("sho-125/qgarch")
library(qgarch)
```

## Example

The example below installs the package from GitHub, loads the bundled
monthly U.S. dataset, fits a QGARCH(1,1) model with freely estimated
volatility-feedback parameter `lambda`, and then plots the fitted
conditional variance and standardized residuals.

```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("sho-125/qgarch")

# Load package
library(qgarch)

# Load bundled example data
data("us_monthly", package = "qgarch")

# Extract excess returns
x <- us_monthly$ER

# Fit a QGARCH(1,1) model
fit11 <- qgarch_fit(
  x,
  model = "free",
  arch_order = 1L,
  garch_order = 1L
)

# Review estimation results
print(fit11)
summary(fit11)
coef(fit11)
logLik(fit11)

# Plot fitted conditional variance
plot(fit11, which = "sigma2")

# Plot standardized residuals
plot(fit11, which = "standardized")
```












