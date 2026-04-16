# Replication notes

This package was built from the uploaded Campbell and Hentschel (1992) paper
and the supplied QGARCH replication scripts.

Implemented directly from the supplied code:

- QGARCH-M(1,1) with `lambda = 0`
- QGARCH-M(1,1) with restricted `lambda`
- QGARCH-M(1,1) with free `lambda`
- Threshold extension with `lambda_t = lambda1 + lambda2 * I_t`
- Hessian-based standard errors with pseudoinverse fallback
- Likelihood-ratio tests between nested specifications
- Table 1 style moment diagnostics

Excluded on purpose:

- Markov-switching variants
- Daily QGARCH(1,2) implementation
- Data files from the replication archive

The threshold specification follows the supplied code closely, including the
constant `kappa` definition based on the sample mean of the indicator.
