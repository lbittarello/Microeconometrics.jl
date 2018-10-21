# Methods

This package supports the methods for regression models of [*StatsBase.jl*](http://juliastats.github.io/StatsBase.jl/stable/statmodels.html).

The following functions are available for all models: `nobs` and `response`.

The following functions are available for parametric models: `dof`, `dof_residual`, `coef`, `stderr`, `vcov`, `confint` and `coefnames`. Note that all methods refer to the second stage of two-stage models.

The following functions are available for maximum-likelihood estimators: `deviance`, `nulldeviance`, `loglikelihood`, `nullloglikelihood`, `r2` and `adjr2`. There are also R² methods for OLS and IV. The following functions are available from *StatsBase.jl*: `aic`, `aicc` and `bic`.

Most models support `predict` and `fitted`. `predict` estimates the index of single-index models. `fitted` computes the conditional outcome expectation. For example, `predict` estimates the Xβ of a logit model and `fitted` computes logistic(Xβ). Support for `residuals` depends on the availability of `fitted`. Out-of-sample forecast is supported.
