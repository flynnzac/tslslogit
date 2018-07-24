# tslslogit
Implements a numerically stable method of estimating logistic regressions which can be useful when convergence issues are likely or when a system has to be robust to such problems arising unexpectedly.  It also works when there is endogeneity as well.  The estimator is based on a contraction mapping result (see included PDF, rlr_letter_zf.pdf, for details). 

Currently, an R implementation is included in the repository, and I plan to add a Stata version as well.

## R Readme

Includes a function, `tslslogit`. It takes as an argument "y" (a binary outcome), "r" (a vector of covariates), and "w" (a vector of instruments, which can be omitted if there is no endogeneity).  It includes one option (`noconstant`) which determines whether to add a constant to the instrument and covariate vectors (by default, it does).  The function returns the estimated parameter vector.
