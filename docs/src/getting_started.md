# Getting Started

## Installation

To install the latest release, run:
```julia
julia> Pkg.add("Microeconometrics")
```
To update to the latest release, run:
```julia
julia> Pkg.update("Microeconometrics")
```
To install the last version, run:
```julia
julia> Pkg.add(PackageSpec(name = "Microeconometrics", rev = "master"))
```
To update to the last version, run:
```julia
julia> Pkg.update(PackageSpec(name = "Microeconometrics", rev = "master"))
```

## Example I: OLS

For an introduction to the package, let's consider an example: ordinary least squares.

We first load some modules:
```julia
julia> using CSV
julia> using DataFrames
julia> using Microeconometrics
```

Next, we load the data:
```julia
julia> S = CSV.read("admit.csv") ; categorical!(S, :rank) ;
```
This sample is available [here](http://github.com/lbittarello/Microeconometrics.jl/tree/master/data), if you wish to replicate this exercise. It comprises 400 observations and four variables (`admit`, `gre`, `gpa` and `rank`). The second command converts `S[:rank]` into a [categorical array](http://juliadata.github.io/DataFrames.jl/stable/man/categorical.html) (a.k.a. a factor variable).

We wish to regress `admit` on `gre`, `gpa`, `rank` and an intercept. We first specify the model:
```julia
julia> M = @micromodel(response => admit, control => gre + gpa + rank + 1) ;
```
The syntax is due to [*StatsModels.jl*](http://juliastats.github.io/StatsModels.jl/latest/formula.html). The `@micromodel` macro generalizes `@formula`.

We then define the [correlation structure](correlation_structures.md) of the errors. Let's assume that observations are independent and identically distributed, so that errors are homoscedastic:
```julia
julia> C = Homoscedastic() ;
```

We now construct the [estimation sample](model_specification.md):
```julia
julia> D = Microdata(S, M, vcov = C) ;
```

We can finally fit the model and visualize the results:
```julia
julia> E = fit(OLS, D)

OLS

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0002    2.0384    0.0415    0.0000  0.0008
gpa             0.1555    0.0640    2.4317    0.0150    0.0302  0.2809
rank: 2        -0.1624    0.0677   -2.3978    0.0165   -0.2951 -0.0296
rank: 3        -0.2906    0.0702   -4.1365    0.0000   -0.4282 -0.1529
rank: 4        -0.3230    0.0793   -4.0726    0.0000   -0.4785 -0.1676
(Intercept)    -0.2589    0.2160   -1.1987    0.2306   -0.6822  0.1644
```

The corresponding code in Stata is:
```stata
. import delimited "admit.csv"
. regress admit gre gpa i.rank
```

To obtain summary statistics, we can apply the methods for regression models of [*StatsBase.jl*](http://juliastats.github.io/StatsBase.jl/stable/statmodels.html):
```julia
julia> nobs(e_ols)

400

julia> r2(e_ols)

0.10040062851886422
```

## Example II: Comparing models

To illustrate more advanced features, suppose that we want to compare specifications. As a first exercise, we wonder if we could drop the rank fixed effects.

We start with the full model. We now assume that the data are heteroscedastic. Because it is the default, we need not specify it.
```julia
julia> M₁ = @micromodel(response => admit, control => gre + gpa + rank + 1) ;
julia> D₁ = Microdata(S, M₁) ;
julia> E₁ = fit(OLS, D₁)

OLS

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0002    2.0501    0.0404    0.0000  0.0008
gpa             0.1555    0.0653    2.3833    0.0172    0.0276  0.2834
rank: 2        -0.1624    0.0729   -2.2266    0.0260   -0.3053 -0.0194
rank: 3        -0.2906    0.0730   -3.9827    0.0001   -0.4336 -0.1476
rank: 4        -0.3230    0.0780   -4.1408    0.0000   -0.4759 -0.1701
(Intercept)    -0.2589    0.2110   -1.2268    0.2199   -0.6725  0.1547
```

Before we estimate the reduced model, we must redefine the control set.:
```julia
julia> M₂ = @micromodel(response => admit, :control => gre + gpa + 1) ;
julia> D₂ = Microdata(S, M₂) ;
```
We can now fit the reduced model:
```julia
julia> E₂ = fit(OLS, D₂)

OLS

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0005    0.0002    2.5642    0.0103    0.0001  0.0010
gpa             0.1542    0.0650    2.3737    0.0176    0.0269  0.2816
(Intercept)    -0.5279    0.2087   -2.5293    0.0114   -0.9370 -0.1188
```

The coefficients on `gre` and `gpa` seem to be robust. For a formal equality test, we use a Hausman test:
```julia
julia> H = hausman_1s(E₁, E₂, ["gre", "gpa"]) ;
julia> tstat(H)

2-element Array{Float64,1}:
  -2.07838
  0.0749552
```
As it turns out, the difference between the coefficients on `gre` is statistically significant.

To further investigate this result, we wish to estimate separate effects by rank. The keyword `subset` helps us construct the appropriate `Microdata`:
```julia
julia> M = @micromodel(response => admit, control => gre + gpa + 1) ;

julia> I₁ = (S[:rank] .== 1) ;
julia> D₁ = Microdata(S, M, subset = I₁) ;
julia> E₁ = fit(OLS, D₁)

OLS

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0006    0.0006    1.0386    0.2990   -0.0005  0.0017
gpa             0.2508    0.1929    1.3000    0.1936   -0.1273  0.6290
(Intercept)    -0.6806    0.5662   -1.2021    0.2293   -1.7903  0.4291

julia> I₂ = (S[:rank] .== 2) ;
julia> D₂ = Microdata(S, M, subset = I₂) ;
julia> E₂ = fit(OLS, D₂)

OLS

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0004    0.8794    0.3792   -0.0004  0.0011
gpa             0.1771    0.1028    1.7219    0.0851   -0.0245  0.3787
(Intercept)    -0.4470    0.3641   -1.2277    0.2196   -1.1607  0.2666

julia> H = hausman_2s(E₁, E₂, ["gre", "gpa"]) ;
julia> tstat(H)

2-element Array{Float64,1}:
  0.334261
  0.337304
```
We used the function `hausman_2s` because these estimates are based on different samples. The difference in the effect of `gre` between ranks 1 and 2 is not significant.
