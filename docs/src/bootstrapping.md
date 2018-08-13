# Bootstrapping

This package does not provide support for bootstrap standard errors at the moment. Nonetheless, it is possible to bootstrap with the existing tools. This tutorial provides some sample code.

We first load some packages:
```julia
using StatsBase
using DataFrames
using CSV
using Microeconometrics
```

We then set up the problem:
```julia
S        = CSV.read(joinpath(datadir, "auto.csv")) ;
S[:gpmw] = ((1.0 ./ S[:mpg]) ./ S[:weight]) * 100000 ;
M        = Dict(:response => "gpmw", :control => "foreign + 1") ;
D        = Microdata(S, M) ;
```

Next, we obtain the coefficient estimates:
```julia
E = fit(OLS, D, novar = true) ;
```

We can now set up the bootstrap:
```julia
srand(0101)

reps = 1000 ;
n    = nobs(E) ;
wgts = fill(0, n) ;
B    = Array{Float64}(reps, dof(E)) ;
```
The vector `wgts` will translate the draw of a bootstrap sample into an input for `Microdata`. The matrix `B` will contain the sample of coefficient estimates. Don't forget to set the seed for the sake of reproducibility!

The algorithm is:
```julia
for b = 1:reps

    wgts .= 0
    draw  = rand(1:n, n)

    for d in draw
        wgts[d] += 1
    end

    Db      = Microdata(S, M, weights = fweights(wgts))
    Eb      = fit(OLS, Db, novar = true)
    B[b, :] = coef(Eb)'
end
```
Note that we do not compute the covariance matrix at each step, which saves us some time.

We can finally see the results:
```julia
E.V = cov(B) ;
coeftable(E)
```
The output is:
```julia
                   Estimate  St. Err.   t-stat.   p-value      C.I. (95%)
foreign: Foreign     0.2462    0.0682    3.6072    0.0003    0.1124  0.3799
(Intercept)           1.609    0.0237   67.9372    <1e-99    1.5626  1.6554
```

You can easily adapt this code to more complex problems (e.g., critical values) or parallelize it for additional speed!