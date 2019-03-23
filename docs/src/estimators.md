# Estimators

The function `fit` estimates models. It returns a model structure, which contains the estimation sample, the coefficients and their covariance matrix. For example, the output of `fit(OLS, MD)` has type `OLS`. Some have additional fields: e.g., two-stage models carry estimates from the first stage and GMM models carry the inverse of the weight matrix.

!!! note

    If you only need coefficients, pass `novar = true` to `fit`.

Model structures are subtypes of broader abstract types, such as `MLE` or `GMM`, which are ultimately instances of [`RegressionModel`](http://juliastats.github.io/StatsBase.jl/stable/statmodels.html). The type hierarchy is:

```
RegressionModel
    ParModel
        GMM
        MLE
    TwoStageModel
```

## Linear regression

### Ordinary least squares

```julia
fit(OLS, MD::Microdata)
```

The `Microdata` must contain: `response` and `control`. See the documentation for linear IV if `Microdata` includes a treatment. `OLS` is a subtype of `MLE`.

### Linear IV

```julia
fit(IV, MD::Microdata; method::String = "TSLS")
```

The following methods are currently implemented:

- `method = "TSLS"`: two-stage least squares;
- `method = "Two-step GMM"`: optimally weighted two-stage GMM with robust covariance matrix;
- `method = "Optimal GMM"`: optimally weighted two-stage GMM with simplified covariance matrix.

Additional methods are available for convenience:

- `method = "OLS"`: linear regression of the outcome on the treatment and controls;
- `method = "First stage"`: linear regression of the treatment on the instruments and controls;
- `method = "Reduced form"`: linear regression of the outcome on the instruments and controls.

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`. `IV` is a subtype of `GMM`.

## Binary choice

```julia
fit(Logit, MD::Microdata)
fit(Probit, MD::Microdata)
fit(Cloglog, MD::Microdata)
fit(Gompit, MD::Microdata)
```

The `Microdata` must contain: `response` and `control`. The outcome should be binary. The model structures are subtypes of `MLE`.

## Count data

```julia
fit(Poisson, MD::Microdata; novar::Bool = false)
```

The `Microdata` must contain: `response` and `control`. The `Microdata` may contain an `offset`. See the documentation for linear IV if `Microdata` includes a treatment. The outcome must be weakly positive. `Poisson` is a subtype of `MLE`.

```julia
fit(IVPoisson, MD::Microdata; novar::Bool = false, method::String = "One-step GMM")
fit(Mullahy, MD::Microdata; novar::Bool = false, method::String = "One-step GMM")
```

`IVPoisson` fits the exponential conditional mean model with additive errors. `Mullahy` fits the exponential conditional mean model with multiplicative errors ([Mullahy, 1997](http://www.jstor.org/stable/2951410)).

The following methods are currently implemented:

- `method = "One-step GMM"`: unweighted one-stage GMM;
- `method = "TSLS"`: one-stage GMM with the average outer product of the instrument vector as weight matrix;
- `method = "Two-step GMM"`: optimally weighted two-stage GMM with robust covariance matrix;
- `method = "Optimal GMM"`: optimally weighted two-stage GMM with simplified covariance matrix.

The models are estimated with the Gauss–Newton algorithm. The first stage of the two-stage specifications is estimated with the average outer product of the instrument vector as weight matrix.

Additional methods are available for convenience:

- `method = "Poisson"`: Poisson regression of the outcome on the treatment and controls;
- `method = "Reduced form"`: Poisson regression of the outcome on the instruments and controls.

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`. The `Microdata` may contain an `offset`. The outcome must be weakly positive. `IVPoisson` and `Mullahy` are subtypes of `GMM`.

## Reweighting methods

!!! note

    All reweighting models require the specification of a first stage. They come in two flavors. In the first, you specify the first-stage model. In the second, you pass a previously fitted model. The latter is more verbose, but it allows you to customize and reuse the first stage.

```julia
fit(IPW, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)
fit(IPW, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)
```

`IPW` estimates average treatment effects by inverse probability weighting. In a first stage, we use model `M₁` to forecast the conditional probability of treatment take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted least squares. The intercept gives the mean outcome of the untreated. We ignore observations whose score is below `trim` or above `1 - trim` (see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The `Microdata` must contain: `response`, `treatment` and `control`. The treatment must be binary. `IPW` is a subtype of `TwoStageModel`.

```julia
fit(Abadie, M₂::Type{ParModel}, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0, kwargs...)
fit(Abadie, M₂::Type{ParModel}, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0 kwargs...)
```

`Abadie` estimates local average response functions according to [Abadie (2003)](https://www.sciencedirect.com/science/article/pii/S0304407602002014). In a first stage, we use model `M₁` to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we fit `M₂` with the weights from the first stage. Keywords customize the second-stage estimator. We ignore observations whose score is below `trim` or above `1 - trim` (see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`. The treatment and the instrument must be binary. `Abadie` is a subtype of `TwoStageModel`.

```julia
fit(FrölichMelly, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)
fit(FrölichMelly, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)
```

This model estimates unconditional local average effects according to [Frölich and Melly (2013)](http://doi.org/10.1080/07350015.2013.803869). In a first stage, we use model `M₁` to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted least squares. The intercept gives mean outcome of untreated compliers. We ignore observations whose score is below `trim` or above `1 - trim` (see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`. The treatment and the instrument must be binary. `FrölichMelly` is a subtype of `TwoStageModel`.

```julia
fit(Tan, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)
fit(Tan, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)
```

This model estimates the unconditional local average treatment effects according to
[Tan (2006)](http://doi.org/10.1198/016214505000001366). In a first stage, we use model `M₁` to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted two-stage least squares. The intercept gives mean outcome of untreated compliers. We ignore observations whose score is below `trim` or above `1 - trim` (see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`. The treatment and the instrument must be binary. `Tan` is a subtype of `TwoStageModel`.
