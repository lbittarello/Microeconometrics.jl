# Model specification

Before fitting the model, one must specify the correlation pattern of the error term
(a `CorrStructure`) and the components of the model (a `Microdata`).

## `CorrStructure`

This structure specifies the correlation between observations.
It determines the calculation of standard errors.

All constructors accept the Boolean keyword `adj`, which defaults to `true`.
If `true`, a finite-sample adjustment is applied to the variance matrix.
The adjustment factor is  / (n - 1),
where n is the number of clusters for clustered data
and the number of observations otherwise.

Four subtypes are currently available:

### Homoscedastic

```julia
Homoscedastic(method::String = "OIM"; adj::Bool = true)
```

Observations are independent and identically distributed.
The optional `method` argument controls the estimation of the variance matrix:
`"OIM"` uses the observed information matrix,
whereas `"OPG"` uses the outer product of the gradient.
(They are equivalent for OLS.)
Only OLS and maximum-likelihood estimators support homoscedastic errors.

### Heteroscedastic

```julia
Heteroscedastic(; adj::Bool = true)
```

Observations are independent, but they may differ in distribution.
This structure leads to sandwich covariance matrices (a.k.a. Huber-Eicker-White).

### Clustered

```julia
Clustered(DF::Microdata, cluster::Symbol; adj::Bool = true)
```

Observations are independent across clusters,
but their joint distribution within clusters is arbitrary.
`cluster` specifies the column of `DF` to cluster on.

### CrossCorrelated

This structure accommodates other correlation structures.
The first argument determines the precise pattern.
The following methods are available:

```julia
CrossCorrelated("Two-way clustering", DF::DataFrame, cluster₁::Symbol, cluster₂::Symbol; adj::Bool = true)
```

Observations may be arbitrarily correlated if they share any cluster.

```julia
CrossCorrelated("Time", DF::DataFrame, time::Symbol, bandwidth::Real; adj::Bool = true)
CrossCorrelated("Time", DF::DataFrame, time::Symbol, kernel::function; adj::Bool = true)
```

The maximum possible correlation between two observations declines
with the time difference between them. Correlation is arbitrary below that limit.
The bandwidth and the kernel function control the upper bound.
`time` specifies the column of `DF` that contains the date of each observation (`Date`).

The first method takes a bandwidth and uses the Parzen kernel.
The second method takes a kernel function instead, which must incorporate the bandwidth.
The first method is equivalent to setting `x -> parzen(x / bandwidth)`.
The following kernels are predefined for convenience:
Bartlett (`bartlett`), Parzen (`parzen`), Truncated (`truncated`)
and Tukey-Hanning (`tukeyhanning`).
See [Andrews (1991)](http://jstor.org/stable/2938229) for formulae.

```julia
CrossCorrelated("Space", DF::DataFrame, latitude::Symbol, longitude::Symbol, bandwidth::Real; adj::Bool = true)
CrossCorrelated("Space", DF::DataFrame, latitude::Symbol, longitude::Symbol, kernel::function; adj::Bool = true)
```

The maximum possible correlation between two observations declines
with the spatial distance between them. Correlation is arbitrary below that limit.
The bandwidth and the kernel function control the upper bound.
`latitude` and `longitude` specify the columns of `DF`
that contain the coordinates of each observation (`Float64`).

For an explanation of the difference between the two methods,
see `CrossCorrelated("time", args...)` above.

```julia
CrossCorrelated("Time and space",
    DF::DataFrame,
    time::Symbol,
    bandwidth_time::Real,
    latitude::Symbol,
    longitude::Symbol,
    bandwidth_space::Real;
    adj::Bool = true)
CrossCorrelated("Time and space",
    DF::DataFrame,
    time::Symbol,
    kernel_time::Function,
    latitude::Symbol,
    longitude::Symbol,
    kernel_space::Function;
    adj::Bool = true)
```

The maximum possible correlation between two observations declines
with the time difference and the spatial distance between them.
Correlation is arbitrary below that limit.
The bandwidths and the kernel functions control the upper bound.
`time` specifies the column of `DF` that contains the date of each observation.
`latitude` and `longitude` specify the columns of `DF`
that contain the coordinates of each observation (`Float64`).

For an explanation of the difference between the two methods,
see `CrossCorrelated("time", args...)` above.

## `Microdata`

This structure combines the functionalities of `Formula`, `ModelFrame` and `ModelMatrix` from
[*StatsModels*](https://github.com/JuliaStats/StatsModels.jl).
It contains a `Matrix{Float64}` (the data matrix),
a map from model components to matrix columns,
a correlation structure and weights (inter alia).

```julia
Microdata(
    DF::DataFrame;
    vcov::CorrStructure = Heteroscedastic(),
    weights::AbstractWeights = UnitWeights(size(DF, 1)),
    subset::AbstractVector{Bool} = trues(size(DF, 1)),
    kwargs...)
```

`subset` determines the estimation sample.
Set a row to `true` if the corresponding row of `DF` should be included
and `false` if it should be excluded.
This keyword is useful in two situations.
First, you have precomputed `vcov` based on the entire sample.
`Microdata` will copy the correlation structure and restrict it to relevant observations.
Second, you are comparing subgroup effects
and observations in different subgroups may correlate
(e.g., they may belong to the same cluster).
`hausman_2s` will account for that correlation
if the `Microdata`s were constructed with `subset`.

`weights` is a [weight vector](http://juliastats.github.io/StatsBase.jl/stable/weights.html).
Except for frequency weights, the weight vector is normalized
to sum up to the number of observations in the sample.

Additional keywords determine the model components.
All regression models need a `response`, but other requirements may vary.
For example, `OLS` asks for `response` and `control`.
See the [introduction](#getting-started) for examples. Conventional sets include:

- `response`: the response, outcome or dependent variable;
- `control`: exogenous explanatory variables (n.b.: you must explicitly include intercepts);
- `treatment`: endogenous explanatory variables;
- `instrument`: instrumental variables (i.e., excluded exogenous variables).

You pass these sets as strings, following
[syntax of `Formula`](http://juliastats.github.io/StatsModels.jl/latest/formula.html).

```julia
Microdata(MD::Microdata; kwargs...)
```

It is also possible to base new `Microdata` on existing `Microdata`.
This constructor allows you to reassign variables to new sets.
You can create new variable sets. If you do not redefine a set, it is preserved.
To suppress a set, redefine it to `""`.
You cannot add new variables, modify the correlation structure, restrict the sample
or reweight observations.

This functionality is useful if you wish to compare specifications.
Rather than building separate data matrices for each one of them,
you can build a master `Microdata`, holding all variables of interest,
and adjust its map as you go through specifications.
