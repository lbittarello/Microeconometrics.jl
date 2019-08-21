# Correlation structures

Before fitting the model, you must specify the correlation between observations (a `CorrStructure`). It determines the calculation of covariance matrices. The default is always `Heteroscedastic`, i.e. independent but not identically distributed observations.

All constructors accept the Boolean keyword `adj` (omitted in the following), which defaults to `true`. If `true`, a finite-sample adjustment is applied to the covariance matrix. The adjustment factor is n / (n - 1), where n is the number of clusters for clustered data and the number of observations otherwise.

Four subtypes are currently available: `Homoscedastic`, `Heteroscedastic`, `Clustered` and `CrossCorrelated`.

## `Homoscedastic`

```julia
Homoscedastic(method::String = "OIM")
```
Observations are independent and identically distributed. The optional argument `method` is only relevant for maximum-likelihood estimators. It controls the estimation of the covariance matrix: `"OIM"` uses the observed information matrix, whereas `"OPG"` uses the outer product of the gradient. Only linear and maximum-likelihood estimators support homoscedastic errors.

## `Heteroscedastic`

```julia
Heteroscedastic()
```
Observations are independent, but they may differ in distribution. This structure leads to sandwich covariance matrices (a.k.a. Huber-Eicker-White).

## `Clustered`

```julia
Clustered(DF::DataFrame, cluster::Symbol)
```

Observations are independent across clusters, but they may differ in their joint distribution within clusters. `cluster` specifies the column of the `DataFrame` to cluster on.

## `CrossCorrelated`

This structure accommodates other correlation structures. The first argument determines the precise pattern.

### Two-way clustering

```julia
CrossCorrelated("Two-way clustering", DF::DataFrame, c₁::Symbol, c₂::Symbol)
```
if two observations share any cluster, they may be arbitrarily correlated.

### Correlation across time

```julia
CrossCorrelated("Time",
        DF::DataFrame,
        time::Symbol,
        bandwidth::Real,
        kernel::Function = parzen
    )
```

The maximum possible correlation between two observations declines with the time difference between them. The actual correlation is arbitrary below that limit. (See [Conley (1999)](https://www.sciencedirect.com/science/article/pii/S0304407698000840).) The bandwidth and the kernel function control the upper bound. `time` specifies the column of `DF` that contains the date of each observation (of type `Date`).

The following kernels are predefined for convenience: Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`) and Tukey-Hanning (`tukeyhanning`). See [Andrews (1991)](http://jstor.org/stable/2938229) for formulae.

!!! warning

    The resulting covariance matrices differ from the Newey-West estimator, which assumes independence across units (though observations for the same unit may correlate across time).

### Correlation across space

```julia
CrossCorrelated("Space",
        DF::DataFrame,
        latitude::Symbol,
        longitude::Symbol,
        bandwidth::Real,
        kernel::Function = parzen
    )
```

The maximum possible correlation between two observations declines with the spatial distance between them. The actual correlation is arbitrary below that limit. (See [Conley (1999)](https://www.sciencedirect.com/science/article/pii/S0304407698000840).) The bandwidth and the kernel function control the upper bound. `latitude` and `longitude` specify the columns of `DF` that contain the coordinates of each observation in decimal degrees (of type `Float64`).

The following kernels are predefined for convenience: Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`) and Tukey-Hanning (`tukeyhanning`). See [Andrews (1991)](http://jstor.org/stable/2938229) for formulae.

### Correlation across time and space

```julia
CrossCorrelated("Time and space",
        DF::DataFrame,
        time::Symbol,
        bandwidth_time::Real,
        latitude::Symbol,
        longitude::Symbol,
        bandwidth_space::Real,
        kernel::Function = parzen
    )
```

The maximum possible correlation between two observations declines with the time difference and the spatial distance between them. The actual correlation is arbitrary below that limit. (See [Conley (1999)](https://www.sciencedirect.com/science/article/pii/S0304407698000840).) The bandwidths and the kernel function control the upper bound. `time` specifies the column of `DF` that contains the date of each observation. `latitude` and `longitude` specify the columns of `DF` that contain the coordinates of each observation in decimal degrees (`Float64`).

The following kernels are predefined for convenience: Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`) and Tukey-Hanning (`tukeyhanning`). See [Andrews (1991)](http://jstor.org/stable/2938229) for formulae.
