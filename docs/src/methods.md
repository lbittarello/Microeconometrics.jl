# Methods

## Methods from *StatsBase*

This package supports the
[methods for regression models](http://juliastats.github.io/StatsBase.jl/stable/statmodels.html)
of *StatsBase*.

The following functions are available for all models: `nobs` and `model_response`.

The following functions are available for parametric models: , `dof`, `dof_residual`,
`coef`, `stderr`, `vcov`, `confint`, `coefnames` and `coeftable`.
Three keywords customize the behavior of `coeftable`:
`verbose` (Boolean: suppresses printing),
`digits` (integer: controls rounding) and
`level` (float: controls the level of the confidence interval – 0.0 for none).
Note that all methods refer to the second stage of two-stage models.

The following functions are available for maximum-likelihood estimators:
`aic`, `aicc`, `bic`, `deviance`, `nulldeviance`, `loglikelihood`, `nullloglikelihood`,
`r2` and `adjr2`. There are also R² methods for OLS and IV.

Some models support `predict` and `fitted` (see the documentation).
`predict` estimates the index of single-index models.
`fitted` estimates the conditional outcome expectation.
For example, `predict` estimates the Xβ of a logit model,
whereas `fitted` estimates logistic(Xβ).
Support for `residuals` depends on the availability of `fitted`.
Out-of-sample forecast is planned for a future release.

## Additional methods

- `tstat`: the *t*-statistic (i.e., the ratio of coefficients to standard error);
- `pval`: the *p*-value of a two-sided significance test;
- `first_stage`: the first-stage estimates of a two-stage model.

## Hausman test

This function computes the difference in coefficients between two parametric models.
They return a `ParObject`, which contains the vector of differences, their covariance matrix
and labels. Our implementation is based on the GMM representation of the joint estimation
problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)).

```julia
hausman_1s(
    model₁::Union{ParModel, TwoStageModel},
    model₂::Union{ParModel, TwoStageModel},
    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))
```

This function is appropriate when both `model₁` and `model₂` were based on
the same estimation sample.

```julia
hausman_2s(
    model₁::Union{ParModel, TwoStageModel},
    model₂::Union{ParModel, TwoStageModel},
    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))
```

This function is appropriate when `model₁` and `model₂` were based on independent samples.
For example, the samples might consist of independent observations with no overlap.

```julia
hausman_2s(
    model₁::Union{ParModel, TwoStageModel},
    model₂::Union{ParModel, TwoStageModel},
    corr::CorrStructure,
    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))
```

This function is appropriate when `model₁` and `model₂` were based on dependent samples.
For example, the samples might consist of independent observations with some overlap
or clustered observations with common clusters.
The correlation structure `corr` must specify the correlation between all observations of
both estimation samples. For example, you could construct `corr` for the entire dataset
and construct the samples via the `subset` keyword to `Microdata`.
