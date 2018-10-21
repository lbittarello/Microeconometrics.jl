#==========================================================================================#

# TYPE

mutable struct OLS <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    OLS() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function OLS(MD::Microdata)
    obj        = OLS()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::OLS, ::UnitWeights)
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

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::OLS) = Diagonal(residuals(obj)) * getmatrix(obj, :control)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::OLS, ::UnitWeights)      = - crossprod(getmatrix(obj, :control))
jacobian(obj::OLS, w::AbstractWeights) = - crossprod(getmatrix(obj, :control), w)

# HOMOSCEDASTIC VARIANCE MATRIX

function _vcov!(obj::OLS, corr::Homoscedastic, w::UnitWeights)
    σ²    = sum(abs2, residuals(obj)) / dof_residual(obj)
    obj.V = lmul!(-σ², inv(jacobian(obj, w)))
end

function _vcov!(obj::OLS, corr::Homoscedastic, w::AbstractWeights)
    σ²    = sum(abs2.(residuals(obj)), w) / dof_residual(obj)
    obj.V = lmul!(-σ², inv(jacobian(obj, w)))
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::OLS, MD::Microdata)
    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")
    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::OLS, MD::Microdata) = predict(obj, MD)

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::OLS) = copy(getmatrix(obj, :control))

#==========================================================================================#

# UTILITIES

coefnames(obj::OLS) = getnames(obj, :control)
mtitle(obj::OLS)    = "OLS"

# COEFFICIENT OF DETERMINATION

adjr2(obj::OLS) = 1.0 - (1.0 - r2(obj)) * (nobs(obj) - 1) / dof_residual(obj)
r2(obj::OLS)    = _r2(obj, getweights(obj))

function _r2(obj::OLS, ::UnitWeights)
    y   = response(obj)
    ŷ   = fitted(obj)
    rss = sum(abs2, y .- ŷ)
    tss = sum(abs2, y .- mean(y))
    return 1.0 - rss / tss
end

function _r2(obj::OLS, w::AbstractWeights)
    y   = response(obj)
    μ   = mean(y, w)
    ŷ   = fitted(obj)
    rss = sum(abs2.(y .- ŷ), w)
    tss = sum(abs2.(y .- μ), w)
    return 1.0 - rss / tss
end

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::OLS, w::AbstractWeights)
    n  = sum(w)
    r  = residuals(obj)
    r .= abs2.(r)
    σ² = sum(r, w) / n
    return - (log(σ²) + log2π + 1.0) * n / 2.0
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::OLS, w::AbstractWeights)
    n  = sum(w)
    r  = copy(getvector(obj, :response))
    μ  = mean(r, w)
    r .= abs2.(r .- μ)
    σ² = sum(r, w) / n
    return - (log(σ²) + log2π + 1.0) * n / 2.0
end

# DEVIANCE

function _deviance(obj::OLS, w::AbstractWeights)
    r  = residuals(obj)
    r .= abs2.(r)
    return sum(r, w)
end

# DEVIANCE UNDER NULL MODEL

function _nulldeviance(obj::OLS, w::AbstractWeights)
    r  = copy(getvector(obj, :response))
    μ  = mean(r, w)
    r .= abs2.(r .- μ)
    return sum(r)
end
