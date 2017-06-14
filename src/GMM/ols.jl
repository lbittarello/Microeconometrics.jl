#==========================================================================================#

# TYPE

mutable struct OLS <: ParModel

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    OLS() = new()
end

#==========================================================================================#

# ESTIMATION

function _fit(obj::OLS)
    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    return x \ y
end

function _fit(obj::OLS, w::AbstractVector)
    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    z = transpose(x .* w)
    return (z * x) \ (z * y)
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::OLS) = getmatrix(obj, :control) .* residuals(obj)

function score(obj::OLS, w::AbstractVector)
    x  = getmatrix(obj, :control)
    r  = residuals(obj)
    r .= r .* w
    return x .* r
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::OLS)
    x = getmatrix(obj, :control)
    return - x' * x
end

function jacobian(obj::OLS, w::AbstractVector)
     x  = getmatrix(obj, :control)
     return - x' * (x .* w)
end

# HOMOSCEDASTIC VARIANCE MATRIX

function _vcov(obj::OLS, corr::Homoscedastic)
    σ = std(residuals(obj)) * (nobs(obj) - 1) / dof_residual(obj)
    return σ * inv(jacobian(obj))
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::OLS) = getmatrix(obj, :control) * obj.β

# FITTED VALUES

fitted(obj::OLS) = predict(obj)

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::OLS) = getmatrix(obj, :control)

#==========================================================================================#

# UTILITIES

coefnames(obj::OLS) = getnames(obj, :control)
