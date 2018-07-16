# Model specification

`Microdata` combines the functionalities of `ModelFrame` and `ModelMatrix` from [*StatsModels.jl*](https://github.com/JuliaStats/StatsModels.jl). It contains a `Matrix{Float64}` (the data matrix), a map from model components to matrix columns, a correlation structure and weights (inter alia).

```julia
Microdata(
        DF::DataFrame,
        model::Dict{Symbol, String};
        contrasts::Dict,
        subset::AbstractVector{Bool} = trues(size(DF, 1)),
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights(size(DF, 1))
    )
```

To construct a `Microdata`, two arguments are compulsory: a [`DataFrame`](http://juliadata.github.io/DataFrames.jl/stable/), which contains the data, and a dictionary, which specifies the components of the model of interest.

All regression models need a `response`, but other requirements may vary. For example, `OLS` asks for `response` and `control`. See the [tutorial](getting_started.md) for examples. You pass these sets as strings, following [syntax of `Formula`](http://juliastats.github.io/StatsModels.jl/latest/formula.html). Conventional sets include:

- `response`: the response (a.k.a. outcome or dependent variable);
- `control`: exogenous explanatory variables (n.b.: you must explicitly include intercepts, `+ 1`);
- `treatment`: endogenous explanatory variables;
- `instrument`: instrumental variables (i.e. excluded exogenous variables).

As for the keywords:
- `contrasts`: a dictionary from column labels to [contrast schemes](https://juliastats.github.io/StatsModels.jl/latest/contrasts.html).
- `subset` determines the estimation sample. Set a row to `true` if the corresponding row of `DF` should be included and `false` if it should be excluded. This keyword is useful if you are comparing subgroup effects and observations in different subgroups may correlate
(e.g., they may belong to the same cluster). [`hausman_2s`](methods.md#hausman-test) will account for that correlation if the `Microdata` were constructed with `subset`.
- `weights` is a [weight vector](http://juliastats.github.io/StatsBase.jl/stable/weights.html). Except for frequency weights, the weight vector is normalized to sum up to the number of observations in the sample.
- `vcov` is a [correlation structure](correlation_structures.md).

It is also possible to base new `Microdata` on existing `Microdata`:
```julia
Microdata(MD::Microdata; kwargs...)
```
This constructor allows you to reassign variables to new sets. You can also create new variable sets. If you do not redefine a set, it is preserved. To suppress a set, redefine it to `""`. You cannot add new variables to the data matrix, modify the correlation structure, restrict the sample or reweight observations.

This functionality is useful if you wish to compare specifications. Rather than building separate data matrices for each one of them, you can build a master `Microdata`, holding all variables of interest, and adjust its map as you go through specifications.
