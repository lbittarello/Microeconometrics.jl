# Contributing

This section helps you specify your own estimator. You can then submit a pull request to add it to *Microeconometrics.jl*!

We will analyze the implementation of OLS. Although it is a simple model, others follow the same steps.

The first step is defining the output `struct`:
```julia
mutable struct OLS <: MLE

    sample::Microdata     # estimation sample
    β::Vector{Float64}    # coefficient vector
    V::Matrix{Float64}    # variance matrix

    OLS() = new()
end
```
Every estimator has these fields. Internal utilities rely on them. Some have additional fields (e.g., `IV` stores the estimation method and the weight matrix). Two-stage models store the estimators for each stage instead.

We then define an uninitialized constructor:
```julia
function OLS(MD::Microdata)
    obj        = OLS()
    obj.sample = MD
    return obj
end
```
The functions `_fit!` and `_vcov!` will later set `β` and `V`.

The next step is overloading the function `fit` from [*StatsBase.jl*](https://github.com/JuliaStats/StatsBase.jl):
```julia
function fit(::OLS, MD::Microdata; novar::Bool = false)

    obj = OLS(MD)

    _fit!(obj, getweights(obj))
    novar || _vcov!(obj, getcorr(obj), getweights(obj))

    return obj
end
```
For OLS, we only need to initialize the output object and pass it to `_fit!` and `_vcov!`. (It is not actually necessary to extend `fit` unless you need to perform additional steps before the estimation, as the fallback will suffice.) Note the utilities `getcorr` and `getweights`.

We can now estimate the coefficients. For efficiency, we write separate functions for unweighted and weighted data.
```julia
function _fit!(obj::OLS, w::UnitWeights)
    y     = getvector(obj, :response)
    x     = getmatrix(obj, :control)
    obj.β = x \ y
end

function _fit!(obj::OLS, w::AbstractWeights)
    y     = getvector(obj, :response)
    x     = getmatrix(obj, :control)
    v     = Diagonal(w) * x
    obj.β =  (v' * x) \ (v' * y)
end
```
Notice the internal utilities `getvector` and `getmatrix`. Their first argument is a `Microdata` or a model structure. The following arguments are the model components of interest. You can request several components from `getmatrix` at once. For example, IV needs `x = getmatrix(obj, :treatment, :control)` and `z = getmatrix(obj, :instrument, :control)`. A single matrix is returned in all cases.

!!! warning

    `getvector` and `getmatrix` return views into the underlying data matrix. You should never modify their output, as you would irremediably alter the data. If you need to perform an in-place operation, make a copy beforehand.

OLS does not require nonlinear optimization. If your estimator needs it, you can use the tools of [*Optim.jl*](http://julianlsolvers.github.io/Optim.jl/stable/). See the implementation of `Logit` for an example.

We must now define `score` and `jacobian`. These functions are the building blocks of the variance estimator. The score is the vector of moment conditions. For OLS, it is −xᵢ (yᵢ − xᵢ'β) (the derivative of the objective function). `score` should return the matrix of score vectors in row form. The Jacobian matrix is the derivative of the moment conditions. For OLS, it is xᵢ xᵢ'. `jacobian` should return the weighted sum of Jacobians (i.e., the expected Jacobian × the number of observations).
```julia
function score(obj::OLS)
    x = copy(getmatrix(obj, :control))
    û = residuals(obj)
    return - x .* û
end

function jacobian(obj::OLS, w::UnitWeights)
    x = getmatrix(obj, :control)
    return x' * x
end

function jacobian(obj::OLS, w::AbstractWeights)
    x = getmatrix(obj, :control)
    v = copy(x) .* w
    return x' * v
end
```
`score` returns the score for each observation, so it ignores weights. `jacobian` returns an expectation; therefore, it must account for weights.

We do not need to extend `_vcov!`. The default method will call `score` and `jacobian` and construct the appropriate estimator, accounting for the correlation structure of the data and the type of weights.

We now overload `predict` and `fitted`. For OLS, these functions are equivalent.
```julia
predict(obj::OLS) = getmatrix(obj, :control) * obj.β
fitted(obj::OLS)  = predict(obj)
```
The next step is optional. We extend `jacobexp`, which computes the derivative of fitted values.
```julia
jacobexp(obj::OLS) = copy(getmatrix(obj, :control))
```
`jacobexp` is only necessary when the estimator serves as the first stage of a two-stage estimator. By extending it, you make your estimator available to two-stage estimators.

We conclude with a function to retrieve coefficient labels:
```julia
coefnames(obj::OLS) = getnames(obj, :control)
```
The syntax of `getnames` is similar to that of `getmatrix`.

You can implement additional methods. For example, *Microeconometrics.jl* extends `r2`, `adjr2`, `loglikelihood`, etc. to OLS.
