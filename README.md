# Microeconometrics

## An example

Ordinary least squares:

```julia
julia> using DataFrames, RDatasets, Microeconometrics

julia> dset  = dataset("datasets", "Formaldehyde")
julia> d_ols = Microdata(dset, corr = Homoscedastic(), response = "OptDen", control = "Carb + 1")
julia> e_ols = fit(OLS, d_ols)

julia> coeftable(e_ols)
              Estimate  St. Err.   t-stat.   p-value      C.I. (0.95%)  
Carb            0.8763    0.0135   64.7444    <1e-99    0.8498    0.9028
(Intercept)     0.0051    0.0078    0.6492    0.5162   -0.0103    0.0204
```

To fit a model, we must first load the data (second line)
and transform it into `Microdata` (third line).
This intermediary step lets us specify the response and controls of interest,
as well as the correlation structure of the data.
Note that we must explicitly include an intercept in the control set if we want one.
We can then fit the model (fourth line) and display the results (fifth line).

## CorrStructure

This structure specifies the correlation between observations and
determines the estimation of standard errors.
Four subtypes are currently available:

#### Homoscedastic

The data are independent and identically distributed.
Only OLS and maximum-likelihood estimators support homoscedastic errors.
The constructor takes one optional argument:
the method for the estimation of the covariance matrix.
The default, `"OIM"`, uses the observed information matrix,
whereas `"OPG"` uses the outer product of the gradient.
(They are equivalent for OLS.)

#### Heteroscedastic

Observations are independent, but they may differ in distribution.
This structure leads to sandwich covariance matrices.

#### Clustered

Data are independent across clusters,
but heteroscedasticity and correlation are arbitrary within clusters.
The constructor takes two mandatory arguments:
a `DataFrame` and a `Symbol`, which specifies the column of the dataset to cluster on.

#### CrossCorrelated

This structure allows for arbitrary correlation and heteroscedasticity across observations.
For example, HAC errors fall into this category. Constructors are prefixed by `cc_`.
For the moment, only one constructor is available:

```julia
cc_timespace(
    dset::DataFrame,
    date::Symbol,
    date_bandwidth::Real,
    latitude::Symbol,
    longitude::Symbol,
    spatial_bandwidth::Real;
    date_kernel::Function = parzen,
    spatial_kernel::Function = parzen
)
```

The resulting structure allows for arbitrary correlation between observations,
as long as it declines over time and space.
`date`, `latitude` and `longitude` specify columns of `dset`.
`dset[date]` must be a `DataVector{Date}`.
Available kernels include: Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`),
and Tukey-Hanning (`tukeyhanning`). See Andrews (1991, ECTA) for formulae.
(Note that you may pass your own kernel function.)

## Microdata

This structure combines the
functionalities of `Formula`, `ModelFrame` and `ModelMatrix` from
[DataFrames](https://github.com/JuliaStats/DataFrames.jl) and
[StatsModels](https://github.com/JuliaStats/StatsModels.jl).
It comprises a `Matrix{Float64}`, a label for each column,
a map from variable sets to columns and a correlation structure (inter alia).

The main constructor takes one mandatory argument: a `DataFrame`.
The second argument is optional and specifies the estimation sample.
You may pass a `BitVector`, a `Vector{Bool}` or a `DataVector{Bool}`,
whose elements are `true` if the corresponding row of the `DataFrame` should be included
and `false` if it should be excluded.
You may then pass variable sets as keywords, following the
[syntax of `Formula`](http://juliastats.github.io/DataFrames.jl/stable/man/formulas/).
The documentation of each model described its requirements.
Conventional sets include:

- `response`: the response, outcome or dependent variable.
- `control`: exogenous explanatory variables (one must explicitly include intercepts, "+ 1").
- `treatment`: a treatment or endogenous explanatory variables.
- `instrument`: instrumental variables or excluded exogenous variables.
- `weight`: a weight vector.

The following keywords control the behavior of the constructor:

- `corr`: a `CorrStructure` (see the documentation above).
The default is `Heteroscedastic()`.
- `makecopy`: whether the correlation structure should be copied (`true`)
or a reference suffices (`false`). The default is `true`.
- `checkrank`: whether the column rank of the model matrix should be checked.
A failure indicates multicollinearity. The default is `true`.
- `normalize`: whether the weight vector should be normalized,
such that the sum is equal to the number of observations.
The default is `true`. Deactivate if you have frequency weights.

It is also possible to base new `Microdata` on existing `Microdata`.
This alternative constructor takes the parent as a mandatory argument and new variable sets as keywords.
(If you don't redefine a set, it is preserved. To suppress a set, define it as `""`.)
The new object will contain its own variable map and references to other fields of the parent,
unless you turn on `makecopy` (which defaults to `false` for this method);
in that case, other fields are copied, including the model matrix.
This functionality is useful if you wish to compare specifications.
Rather than creating separate model matrices for each specification,
you can build a master `Microdata`, holding all variables of interest,
and adjust its map as you need to.

## Estimation

To estimate a model, use `fit(model, args...)`.
The precise arguments depend on the model in question.
Read on for further information.

`fit` returns a model structure.
In the example above, the output of `fit(OLS, d_ols)` has type `OLS`.
Parametric one-stage models populate these structures with a pointer to
the estimation sample, the estimated coefficients and their covariance matrix.
Some models include additional information
– e.g., GMM estimators record the choice of weight matrix.
Two-stage models wrap the output of the first and second stages instead,
along with auxiliary objects for the second stage.

Only analytical covariance matrices are currently available.
Bootstrapping is planned for a future release. If you only need coefficients,
you may pass `novar = true` to `fit` to suppress covariance estimation.

Some functions return parametric objects (`ParObject`).
These structures contain estimates, their covariance matrix and labels.
For an example, see the documentation of `hausman` (under "Inference").

## Methods

This package extends many methods of
[StatsBase](https://github.com/JuliaStats/StatsBase.jl),
[DataFrames](https://github.com/JuliaStats/DataFrames.jl) and
[StatsModels](https://github.com/JuliaStats/StatsModels.jl), such as:

- `nobs` and `model_response`: Available for all models.
- `coef`, `stderr`, `vcov`, `confint` and `coefnames`:
Available for all parametric and two-stage models, as well as parametric objects.
- `coeftable`:
Available for all parametric and two-stage models, as well as parametric objects.
The following keywords customize behavior:
`verbose` (boolean: suppresses printing),
`digits` (integer: controls rounding) and
`level` (float: controls the level of the confidence interval – 0.0 for none).
- `dof` and `dof_residual`: Available for all parametric and two-stage models.
- `fitted`: Available for most models.
This function estimates the conditional expectation of the response.
- `predict`: Available for single-index models.
This function estimates the linear predictor (xᵀβ).
- `r2` and `adjr2`: Available if the model supports `fitted`.
By default, `r2` computes Efron's R². Maximum-likelihood estimators support variants.
- `aic`, `aicc`, `bic`, `deviance`, `nulldeviance`, `loglikelihood`, `nullloglikelihood`:
Available for maximum-likelihood estimators.

We refer the reader to the
[documentation of StatsBase](http://statsbasejl.readthedocs.io/en/latest/statmodels.html)
for further information.

We also implement additional methods:

- `tstat`: the _t_-statistic (i.e., the ratio of coefficients to standard error).
- `pval`: the _p_-value of a two-sided significance test.
- `reduced_form`: the reduced form of a linear IV model.
- `first_stage`: the first stage of a two-stage or IV model.
- `second_stage`: the second stage of a two-stage model.

Estimation of marginal effects are planned for a future release.

## Models

More models are planned. Feel free to contribute or suggest new models!

#### Ordinary least squares

```julia
fit(OLS, MD::Microdata)
```

The `Microdata` must contain: `response` and `control`.

#### Maximum-likelihood estimators

```julia
fit(Logit, MD::Microdata)
```

The `Microdata` must contain: `response` and `control`.

```julia
fit(Probit, MD::Microdata)
```

The `Microdata` must contain: `response` and `control`.

#### GMM

```julia
fit(IV, MD::Microdata; method::String = "Optimal")
```

Linear IV regression.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.
Multiple endogenous variables are supported (except for `first_stage`).
The default method is two-step optimal GMM. For two-stage least squares, use `"Unadjusted"`.

#### Treatment evaluation

```julia
fit(Abadie, FS::Micromodel, SS::ParModel, MD::Microdata; trim::AbstractFloat = 0.01)
```

This model estimates a local average response function according to Abadie (2003, JE).
In a first stage, we use model `FS` to forecast the probability of treatment take-up
and construct weights. Because of outliers, we give zero weight to observations
such that the fitted probability is below `trim` or above `1 - trim`.
In a second stage, we fit model `SS` with the weights from the first stage.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.

```julia
fit(FrölichMelly, FS::Micromodel, MD::Microdata; trim::AbstractFloat = 0.01)
```

This model estimates unconditional local average treatment effects according to
Frölich and Melly (2013, JBES).
In a first stage, we use model `FS` to forecast the probability of treatment take-up
and construct weights. Because of outliers, we give zero weight to observations
such that the fitted probability is below `trim` or above `1 - trim`.
In the second stage, we run weighted OLS on the treatment and an intercept.
The intercept gives the mean outcome in the absence of treatment for compliers.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.

## Inference

More tests are planned. Feel free to contribute or suggest new tests!

```julia
hausman_test(M1::Micromodel, M2::Micromodel [, names::Vector{String}, corr::CorrStructure])
```

This function computes the difference between the coefficients of two parametric models
and their covariance matrix. Its output is a parametric object.
By default, it computes the difference between all shared coefficients.
If you are only interested in a subset,
you may pass the labels of the relevant variables in the vector `names`.
You may also provide a correlation structure. If you don't,
`hausman_test` will attempt to use to the correlation structure of `M1`.
This default always works for homoscedastic or heteroscedastic data.
For dependent data (e.g., clustered), you must have fit both models to the same data.
Specifying a correlation structure is useful
if you wish to compare coefficients across subgroups of dependent data.
You may compute a correlation structure for the entire dataset and
feed it to `hausman_test`.

## Estimation table

The function `etable` prints a table of estimates.
It is useful if you wish to compare specifications.
It takes any number of fitted models or parametric objects as arguments.
Significance stars and standard errors accompany estimates coefficients.

It accepts the following keyword arguments:

- `digits` (integer): controls rounding.
- `compact` (boolean): suppresses the display of standard errors.
- `f` (function): replaces standard errors with the output of `f` (e.g., `tstat`).
- `stars` (matrix): a two-column matrix,
whose first column species significance levels and
whose second columns specifies corresponding adornments (strings).
The default is `[0.1 "*"; 0.05 "**"; 0.01 "***"]`.
Use `[]` to suppress the display of significance stars.
- `titles` (string vector): a vector of model labels. Defaults to `[(1), (2)...]`.

Display of summary statistics (e.g., the R² or the number of observations)
is planned for a future release.
