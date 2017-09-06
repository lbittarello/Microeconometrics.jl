# Microeconometrics

## An example

Ordinary least squares:

```julia
julia> using DataFrames, RDatasets, Microeconometrics

julia> DF    = dataset("datasets", "Formaldehyde")
julia> CS    = Homoscedastic()
julia> d_ols = Microdata(DF, response = "OptDen", control = "Carb + 1", corr = CS)
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

## Microdata

This structure combines the
functionalities of `Formula`, `ModelFrame` and `ModelMatrix` from
[DataFrames](https://github.com/JuliaStats/DataFrames.jl) and
[StatsModels](https://github.com/JuliaStats/StatsModels.jl).
It contains a `Matrix{Float64}`, a label for each column,
a map from variable sets to columns
and a [correlation structure](#corrstructure) (inter alia).

```julia
Microdata(DF::DataFrame [, subset::AbstractVector{Bool} = trues(size(df, 1))]; kwargs...)
```

The main constructor takes one mandatory argument: a `DataFrame`.
The second argument is optional and specifies the estimation sample.
Set each row to `true` if the corresponding row of the `DataFrame` should be included
and `false` if it should be excluded. By default, `Microdata` uses all complete rows.

You may then pass variable sets as keywords, following the
[syntax of `Formula`](http://juliastats.github.io/DataFrames.jl/stable/man/formulas/).
All regression models need a `response`, but other requirements may vary.
For example, `OLS` asks for `response` and `control`. Conventional sets include:

- `response`: the response, outcome or dependent variable.
- `control`: exogenous explanatory variables (n.b.: one must explicitly include intercepts, "+ 1").
- `treatment`: endogenous explanatory variables.
- `instrument`: instrumental variables (i.e., excluded exogenous variables).
- `weight`: a weight vector.

The following reserved keywords control the behavior of the constructor:

- `corr`: a [`CorrStructure`]((#corrstructure). The default is `Heteroscedastic()`.
- `makecopy`: whether the `CorrStructure` should be copied (`true`)
or a reference suffices (`false`). The default is `true`.
`Microdata` may need to modify the correlation structure,
so copying is recommended.
- `checkrank`: whether the column rank of the model matrix should be checked.
A failure indicates multicollinearity. The default is `true`.
- `normalize`: whether the weight vector should be normalized,
such that the sum is equal to the number of observations.
The default is `true`. Deactivate if you have frequency weights.

```julia
Microdata(MD::Microdata; kwargs...)
```

It is also possible to base new `Microdata` on existing `Microdata`.
This alternative constructor takes the parent as a mandatory argument
and new variable sets as keywords.
(If you don't redefine a set, it is preserved. To suppress a set, redefine it to be `""`.)
This functionality is useful if you wish to compare specifications.
Rather than building separate data matrices for each one of them,
you can build a master `Microdata`, holding all variables of interest,
and adjust its map as you go through specifications.
It is also useful if you want to switch from weighted to unweighted estimators:
use ``

The only reserved keyword is `makecopy`.
If `false` (the default), the new object will contain its own variable map
and references to other fields of the parent.
If `true`, all fields are copied, including the model matrix.

## CorrStructure

This structure specifies the correlation between observations and
determines the estimation of standard errors.
Four subtypes are currently available:

#### Homoscedastic

```julia
Homoscedastic([method::String = "OIM"])
```

Observations are independent and identically distributed.
The constructor takes the method for the estimation of the covariance matrix
as an optional argument. (Note that it is not a keyword.)
`"OIM"` uses the observed information matrix,
whereas `"OPG"` uses the outer product of the gradient.
(They are equivalent for OLS.)
Only OLS and maximum-likelihood estimators support homoscedastic errors.

#### Heteroscedastic

```julia
Heteroscedastic()
```

Observations are independent, but they may differ in distribution.
This structure leads to sandwich covariance matrices (a.k.a. Huber-Eicker-White).

#### Clustered

```julia
Clustered(DF::Microdata, cluster::Symbol)
```

Observations are independent across clusters,
but heteroscedasticity and correlation may be arbitrary within clusters.
`cluster `specifies the column of `DF` to cluster on.

#### CrossCorrelated

This structure allows for arbitrary correlation and heteroscedasticity across observations.
For example, HAC errors fall into this category. Constructors are prefixed by `cc_`.
For the moment, only one constructor is available:

```julia
cc_timespace(
    DF::DataFrame,
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
`date`, `latitude` and `longitude` specify columns of `DF`.
`DF[date]` must be a `DataVector{Date}`.
Predefined kernels include: Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`)
and Tukey-Hanning (`tukeyhanning`). See Andrews (1991, ECTA) for formulae.

## Estimation

Only analytical covariance matrices are currently available.
Bootstrapping is planned for a future release. If you only need coefficients,
you may set `novar = true` to suppress covariance estimation.

If `Microdata` contains a weight vector (i.e., a `weight` set),
estimates will be weighted. To turn it off,
you need to [adjust the map](#microdata) of your `Microdata`.

#### Type hierarchy

When you fit a model, you obtain a model structure.
For example, the output of `fit(OLS, MD)` has type `OLS`.
These structure are instances of broader abstract types, such as
`ParModel` or `TwoStageModel`, which belong in turn to the supertype `Micromodel`.

Some functions yield parametric objects (`ParObject`).
For example, see the documentation of `hausman_test` (under [Inference](#inference)).

#### One-stage models

To estimate one-stage models, use `fit`.
The output generally contains the estimated coefficients, their covariance matrix,
the estimation method and a pointer to the estimation sample.

#### Two-stage models

You can estimate two-stage models in two ways.
The simpler method consists of passing the first-stage model to `fit`.
As an alternative, you may call `first_stage` to estimate the first stage
and pass its output to `fit`.
This second strategy is useful if you wish to customize the first-step estimator.
It may also help you avoid repeated computations.
For example, `Abadie` and `FrölichMelly` have the same first stage,
so you needn't compute it twice.

The output wraps the results from the first and second stages,
along with auxiliary objects for the second stage.

## Models

More models are planned. Feel free to contribute or suggest new models!

#### Ordinary least squares

```julia
fit(OLS, MD::Microdata [; novar::Bool = false])
```

The `Microdata` must contain: `response` and `control`.

#### Maximum-likelihood estimators

```julia
fit(Logit, MD::Microdata [; novar::Bool = false])
```

The `Microdata` must contain: `response` and `control`.

```julia
fit(Probit, MD::Microdata [; novar::Bool = false])
```

The `Microdata` must contain: `response` and `control`.

#### General method of moments

```julia
fit(IV, MD::Microdata [; novar::Bool = false, method::String = "One-step GMM"])
```

Linear IV regression.
Multiple endogenous variables are supported.
For now, only two-stage least squares is implemented (`method = "One-step GMM"`).
Two-step GMM is planned for a future release.
For convenience, you may also set `method = "OLS"`.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.

Additional functionality:

- `first_stage(IV, MD::Microdata [; novar::Bool = false])`:
linear regression of the treatment on the instrument(s) and controls.
The treatment must be univariate.
- `reduced_form(IV, MD::Microdata [; novar::Bool = false])`:
linear regression of the response on the instrument(s) and controls.

#### Treatment evaluation

```julia
fit(Abadie, M₂::Type{ParModel}, M₁::Type{Micromodel}, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
fit(Abadie, M₂::Type{ParModel}, m₁::Micromodel, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
```

This model estimates a local average response function according to Abadie (2003, JE).
It assumes that both the treatment and the instrument are binary.
In a first stage, we use model `M₁` to forecast the probability of treatment take-up
and construct weights. Because of outliers, we give zero weight to observations
whose fitted score is below `trim` or above `1 - trim`.
In a second stage, we fit model `M₂` with the weights from the first stage.
Optional keywords customize the second-stage estimator.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.

Additional functionality:

- `first_stage(Abadie, M₁::Type{Micromodel}, MD::Microdata [; novar::Bool = false])`
fits the first-stage model.

```julia
fit(FrölichMelly, M₁::Type{Micromodel}, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
fit(FrölichMelly, m₁::Micromodel, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
```

This model estimates unconditional local average treatment effects according to
Frölich and Melly (2013, JBES).
It assumes that both the treatment and the instrument are binary.
In a first stage, we use model `M₁` to forecast the probability of treatment take-up
and construct weights. Because of outliers, we give zero weight to observations
such that the fitted probability is below `trim` or above `1 - trim`.
In the second stage, we run weighted OLS on the treatment and an intercept.
The intercept gives the mean outcome in the absence of treatment for compliers.
Optional keywords customize the second-stage estimator.
The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.

Additional functionality:

- `first_stage(FrölichMelly, M₁::Type{Micromodel}, MD::Microdata [; novar::Bool = false])`
fits the first-stage model.

```julia
fit(Tan, M₁::Type{Micromodel}, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
fit(Tan, m₁::Micromodel, MD::Microdata
    [; novar::Bool = false, trim::AbstractFloat = 0.0, kwargs...])
```

This model estimates unconditional local average treatment effects according to
the reweighting method of Tan (2013, JASA). The syntax is similar to `FrölichMelly`.

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

Estimation of marginal effects are planned for a future release.

## Inference

More tests are planned. Feel free to contribute or suggest new tests!

```julia
hausman_test(M₁::Micromodel, M₂::Micromodel [, names::Vector{String}, corr::CorrStructure])
```

This function computes the difference between the coefficients of two parametric models
and their covariance matrix. Its output is a parametric object.
By default, it computes the difference between all shared coefficients.
If you are only interested in a subset,
you may pass the labels of the relevant variables in the vector `names`.

You may also provide a correlation structure. If you don't,
`hausman_test` will attempt to use to the correlation structure of `M₁`.
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
- `aux` (function): replaces standard errors with the output of `aux` (e.g., `tstat`).
- `stars` (matrix): a two-column matrix,
whose first column species significance levels and
whose second columns specifies corresponding adornments (strings).
The default is `[0.1 "*"; 0.05 "**"; 0.01 "***"]`.
Use `[]` to suppress the display of significance stars.
- `titles` (string vector): a vector of model labels. Defaults to `[(1), (2)...]`.

Display of summary statistics (e.g., the R² or the number of observations)
is planned for a future release.
