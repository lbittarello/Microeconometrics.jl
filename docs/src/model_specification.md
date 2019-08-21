# Model specification

Before you estimate a model, you must load the data and specify the role of each variable in the model. This packages defines `Microdata` for that purpose.

```julia
Microdata(
        DF::DataFrame,
        model::Dict{Symbol, String};
        hints::Dict{Symbol, TermOrTerms},
        subset::AbstractVector{Bool} = trues(size(DF, 1)),
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights(size(DF, 1))
    )
```

To construct a `Microdata`, two arguments are compulsory: a [`DataFrame`](http://juliadata.github.io/DataFrames.jl/stable/), which contains the data, and a dictionary, which specifies the components of the model of interest. To construct this dictionary, use the macro `@micromodel`.

All regression models need a `response`, but other requirements may vary. (Check the [documentation](estimators.md)!) For example, `OLS` asks for `response` and `control`. In defining these sets, follow the [syntax of `Formula`](http://juliastats.github.io/StatsModels.jl/latest/formula.html). (See the [tutorial](getting_started.md) for examples.) Conventional sets include:

- `response`: the response (a.k.a. outcome or dependent variable);
- `control`: exogenous explanatory variables (n.b.: you must explicitly include intercepts, `+ 1`);
- `offset`: an exogenous variable whose coefficient is constrained to unity;
- `treatment`: endogenous explanatory variables;
- `instrument`: instrumental variables (i.e. excluded exogenous variables).

As for the keywords:
- `hints`: a dictionary from column labels to [schemas](https://juliastats.github.io/StatsModels.jl/latest/internals.html) or [contrasts](https://juliastats.github.io/StatsModels.jl/latest/contrasts/).
- `subset` determines the estimation sample. Set an entry to `true` if the corresponding row of `DF` should be included and `false` if it should be excluded. This keyword is useful if you are comparing subgroups and observations in different subgroups may correlate (e.g., they may belong to the same cluster). [`hausman_2s`](hypothesis_tests.md#hausman-test) will account for that correlation if the `Microdata` were constructed with `subset`.
- `weights` is a [weight vector](http://juliastats.github.io/StatsBase.jl/stable/weights.html). Except for frequency weights, the weight vector is normalized to sum up to the number of observations in the sample.
- `vcov` is a [correlation structure](correlation_structures.md).
