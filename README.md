# qgarch

`qgarch` is a small R package for the QGARCH-in-mean model family used in the
Campbell and Hentschel (1992) volatility-feedback framework and in the
replication scripts you supplied.

## Implemented models

- `model = "zero"`: QGARCH-M(1,1) with `lambda = 0`
- `model = "restricted"`: `lambda = rho * gamma * alpha / (1 - rho * (alpha + beta))`
- `model = "free"`: `lambda` estimated freely
- `model = "threshold"`: `lambda_t = lambda1 + lambda2 * I_t`

The package intentionally **does not** include the Markov-switching extension.

## Installation

From the package source directory:

```r
install.packages("qgarch_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Minimal example

```r
library(qgarch)
set.seed(123)
x <- rnorm(300, sd = 0.04)

fit <- qgarch_fit(x, model = "free")
print(fit)
summary(fit)
plot(fit)
predict(fit, n.ahead = 5)
```

## Notes

- The replication scripts use data-specific starting values. To make the package
  more convenient, `qgarch_fit()` tries several default starting vectors and
  keeps the best converged solution.
- The restricted-lambda model uses `rho = 1` by default so that the package
  matches the supplied replication code. You can override this with the `rho`
  argument.
