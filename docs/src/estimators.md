# Estimators

The function `fit` estimates models. It returns a model structure, which contains
the estimation sample, coefficients and their variance matrix (inter alia).
For example, the output of `fit(OLS, MD)` has type `OLS`.
Some have additional fields: for instance, the structures of
two-stage models carry estimates from the first stage.
These structure are instances of broader abstract types, such as
`MLE` or `TwoStageModel`, which belong in turn to the supertype `Micromodel`.

If you only need coefficients, set `novar = true`.

## Models for exogenous regressors

### Ordinary least squares

```julia
fit(OLS, MD::Microdata; novar::Bool = false)
```

The `Microdata` must contain: `response` and `control`.
`OLS` is a subtype of `ParModel`.

### Binary choice

```julia
fit(Logit, MD::Microdata; novar::Bool = false)
```

The `Microdata` must contain: `response` and `control`.
`Logit` is a subtype of `MLE` and `ParModel`.

```julia
fit(Probit, MD::Microdata; novar::Bool = false)
```

The `Microdata` must contain: `response` and `control`.
`Probit` is a subtype of `MLE` and `ParModel`.

### Treatment evaluation

```julia
fit(IPW,
    M₁::Type{Micromodel},
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0,
    kwargs...)
fit(IPW,
    m₁::Micromodel,
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0)
```

The `Microdata` must contain: `response`, `treatment` and `control`.
The treatment must be binary. `IPW` is a subtype of `TwoStageModel`.

This model estimates the average treatment effect by inverse probability weighting.
In a first stage, we use model `M₁` to forecast the conditional probability of
treatment take-up (the propensity score) and construct weights, so that
the weighted covariate distribution is similar across treatment subsamples.
In the second stage, we run weighted OLS on the treatment and an intercept.
The intercept gives the mean outcome in the absence of treatment.
We ignore observations whose score is below `trim` or above `1 - trim`
(see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The first method fits the first-stage model. Keyword arguments customize this step.
The second method uses a previously estimated model instead.

## Models for endogenous regressors

### Linear IV

```julia
fit(IV, MD::Microdata; novar::Bool = false, method::String = "TSLS")
```

The `Microdata` must contain: `response`, `treatment`, `control` and `instrument`.
`IV` is a subtype of `ParModel`.

The following variants are currently implemented:

- `method = "TSLS"`: two-stage least squares;
- `method = "OLS"`: linear regression of the outcome on the treatment and controls;
- `method = "First stage"`: linear regression of the treatment on the instruments and controls;
- `method = "Reduced form"`: linear regression of the outcome on the instruments and controls.

### Reweighting models

```julia
fit(Abadie,
    M₂::Type{ParModel},
    M₁::Type{Micromodel},
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0,
    kwargs...)
fit(Abadie,
    M₂::Type{ParModel},
    m₁::Micromodel,
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0
    kwargs...)
```

The `Microdata` must contain: `response`, `treatment` and `control`.
The treatment and the instrument must be binary. `Abadie` is a subtype of `TwoStageModel`.

This model estimates a local average response function according to
[Abadie (2003)](http://doi.org/10.1016/S0304-4076(02)00201-4).
In a first stage, we use model `M₁` to forecast the conditional probability of
instrument take-up and construct weights, so that the weighted covariate distribution
is similar across treatment subsamples of the compliers.
In a second stage, we fit model `M₂` with the weights from the first stage.
Keywords customize the second-stage estimator.
We ignore observations whose score is below `trim` or above `1 - trim`
(see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The first method fits the first-stage model.
The second method uses a previously estimated model instead.

```julia
fit(FrölichMelly,
    M₁::Type{Micromodel},
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0,
    kwargs...)
fit(FrölichMelly,
    m₁::Micromodel,
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0)
```

The `Microdata` must contain: `response`, `treatment` and `control`.
The treatment and the instrument must be binary.
`FrölichMelly` is a subtype of `TwoStageModel`.

This model estimates the unconditional local average treatment effects according to
[Frölich and Melly (2013)](http://doi.org/10.1080/07350015.2013.803869).
In a first stage, we use model `M₁` to forecast the conditional probability of
instrument take-up and construct weights, so that the weighted covariate distribution
is similar across treatment subsamples of the compliers.
In the second stage, we run weighted OLS on the treatment and an intercept.
The intercept gives compliers' mean outcome in the absence of treatment.
We ignore observations whose score is below `trim` or above `1 - trim`
(see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The first method fits the first-stage model. Keyword arguments customize this step.
The second method uses a previously estimated model instead.

```julia
fit(Tan,
    M₁::Type{Micromodel},
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0,
    kwargs...)
fit(Tan,
    m₁::Micromodel,
    MD::Microdata;
    novar::Bool = false,
    trim::AbstractFloat = 0.0)
```

The `Microdata` must contain: `response`, `treatment` and `control`.
The treatment and the instrument must be binary. `Tan` is a subtype of `TwoStageModel`.

This model estimates the unconditional local average treatment effects according to
[Tan (2006)](http://doi.org/10.1198/016214505000001366).
In a first stage, we use model `M₁` to forecast the conditional probability of
instrument take-up and construct weights, so that the weighted covariate distribution
is similar across treatment subsamples of the compliers.
In the second stage, we run weighted OLS on the treatment and an intercept.
We ignore observations whose score is below `trim` or above `1 - trim`
(see [Crump et al. (2009)](http://jstor.org/stable/27798811)).

The first method fits the first-stage model. Keyword arguments customize this step.
The second method uses a previously estimated model instead.
