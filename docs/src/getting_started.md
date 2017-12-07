# Getting Started

## Installation

To install the package, run:
```julia
Pkg.add("Microeconometrics")
```
To update to the latest release, run:
```julia
Pkg.update("Microeconometrics")
```
To obtain the last version, run:
```julia
Pkg.checkout("Microeconometrics")
```

You will also need [*DataFrames*](http://juliadata.github.io/DataFrames.jl/stable/),
as well as a package to import your data
(e.g., [*CSV*](http://juliadata.github.io/CSV.jl/stable/)).

## Example I: OLS

For an introduction to the package, let's consider an example: ordinary least squares.

We first load the modules:
```julia
julia> using CSV
julia> using DataFrames
julia> using Microeconometrics
```
We then load the data:
```julia
julia> dst = CSV.read("admit.csv")
julia> categorical!(dst, :rank)
```
This sample is available
[here](http://github.com/lbittarello/Microeconometrics.jl/tree/master/data),
if you wish to replicate this exercise.
The last line converts `dst[:rank]` into a
[categorical array](http://juliadata.github.io/DataFrames.jl/stable/man/categorical.html).

We now define the [correlation structure](#corrstructure) of the error term.
Let's assume that observations are independent and identically distributed,
so that errors are homoscedastic:
```julia
julia> corr = Homoscedastic();
```
We next specify the model by constructing a [`Microdata`](#microdata):
```julia
julia> dta = Microdata(dst,
                vcov = corr,
                response = "admit",
                control = "gre + gpa + rank + 1");
```

We can finally fit the model and visualize the results:
```julia
julia> e_ols = fit(OLS, dta);
julia> coeftable(e_ols);

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0002    2.0384    0.0415       0.0  0.0008
gpa             0.1555     0.064    2.4317    0.0150    0.0302  0.2809
rank: 2        -0.1624    0.0677   -2.3978    0.0165   -0.2951 -0.0296
rank: 3        -0.2906    0.0702   -4.1365    <1e-99   -0.4282 -0.1529
rank: 4         -0.323    0.0793   -4.0726    <1e-99   -0.4785 -0.1676
(Intercept)    -0.2589     0.216   -1.1987    0.2306   -0.6822  0.1644
```

The corresponding code in Stata is:
```stata
. import delimited "admit.csv"
. regress admit gre gpa i.rank
```

We can apply the methods for regression models of
[*StatsBase*](http://juliastats.github.io/StatsBase.jl/stable/statmodels.html)
to obtain statistics:
```julia
julia> nobs(e_ols)
400
julia> r2(e_ols)
0.10040062851886422
```

## Example II: Comparing models

To illustrate more advanced features, suppose that we want to compare specifications.
As a first exercise, we wonder if we could drop the rank fixed effects.

We start with the full model. We now assume that the data are heteroscedastic:
```julia
julia> corr = Heteroscedastic();
julia> dta₁ = Microdata(dst, vcov = corr, response = "admit", control = "gre + gpa + rank + 1");
julia> e₁ = fit(OLS, dta₁);
julia> coeftable(e₁);

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0002    2.0501    0.0404       0.0  0.0008
gpa             0.1555    0.0653    2.3833    0.0172    0.0276  0.2834
rank: 2        -0.1624    0.0729   -2.2266    0.0260   -0.3053 -0.0194
rank: 3        -0.2906     0.073   -3.9827    0.0001   -0.4336 -0.1476
rank: 4         -0.323     0.078   -4.1408    <1e-99   -0.4759 -0.1701
(Intercept)    -0.2589     0.211   -1.2268    0.2199   -0.6725  0.1547
```

Before we estimate the reduced model, we need to redefine the control set.
There are two approaches to this task.
We can construct a new `Microdata` from scratch:

```julia
julia> dta₂ = Microdata(dst, vcov = corr, response = "admit", control = "gre + gpa + 1");
```

Or we can modify the control set of our existing `Microdata`:
```julia
julia> dta₂ = Microdata(dta₁, control = "gre + gpa + 1");
```
The second approach does not reconstruct the underlying data matrix.
Therefore, it is faster and uses less memory.
On the other hand, it only allows us to reassign existing variables across variable sets.
We cannot add new variables or modify the correlation structure of the error term.

We next fit the reduced model:
```julia
julia> e₂ = fit(OLS, dta₂);
julia> coeftable(e₂);

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0005    0.0002    2.5642    0.0103    0.0001   0.001
gpa             0.1542     0.065    2.3737    0.0176    0.0269  0.2816
(Intercept)    -0.5279    0.2087   -2.5293    0.0114    -0.937 -0.1188
```
The coeffients on `gre` and `gpa` seem to be robust.

For a formal equality test, we use a Hausman test:
```julia
julia> ht = hausman_1s(e₁, e₂, ["gre", "gpa"]);
julia> tstat(ht)

2-element Array{Float64,1}:
 -2.07838
  0.0749552
```
As it turns out, the difference between the coefficients on `gre` is statistically significant.

To further investigate this result, we wish to estimate separate effects by rank.
The keyword `subset` helps us construct the appropriate `Microdata`s:
```julia
julia> idx₁ = (dst[:rank] .== 1);
julia> dta₁ = Microdata(dst, subset = idx₁, vcov = corr, response = "admit", control = "gre + gpa + 1");
julia> e₁   = fit(OLS, dta₁);
julia> coeftable(e₁);

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0006    0.0006    1.0386    0.2990   -0.0005  0.0017
gpa             0.2508    0.1929       1.3    0.1936   -0.1273   0.629
(Intercept)    -0.6806    0.5662   -1.2021    0.2293   -1.7903  0.4291

julia> idx₂ = (dst[:rank] .== 2);
julia> dta₂ = Microdata(dst, subset = idx₂, vcov = corr, response = "admit", control = "gre + gpa + 1");
julia> e₂   = fit(OLS, dta₂);
julia> coeftable(e₂);

              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  
gre             0.0004    0.0004    0.8794    0.3792   -0.0004  0.0011
gpa             0.1771    0.1028    1.7219    0.0851   -0.0245  0.3787
(Intercept)     -0.447    0.3641   -1.2277    0.2196   -1.1607  0.2666

julia> ht = hausman_2s(e₁, e₂, ["gre", "gpa"]);
julia> tstat(ht)

2-element Array{Float64,1}:
 0.334261
 0.337304
```
We used the function `hausman_2s` because these estimates are based on different samples.
The difference in the effect of `gre` between ranks 1 and 2 is not significant.
