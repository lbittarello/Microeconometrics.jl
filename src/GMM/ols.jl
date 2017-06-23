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
    z = scale!(transpose(x), w)
    return (z * x) \ (z * y)
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::OLS) = scale!(residuals(obj), copy(getmatrix(obj, :control)))

function score(obj::OLS, w::AbstractVector)
    return scale!(w, scale!(residuals(obj), copy(getmatrix(obj, :control))))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::OLS)
    x = getmatrix(obj, :control)
    return crossprod(x, neg = true)
end

function jacobian(obj::OLS, w::AbstractVector)
     x = getmatrix(obj, :control)
     return crossprod(x, w, neg = true)
end

# HOMOSCEDASTIC VARIANCE MATRIX

function _vcov(obj::OLS, corr::Homoscedastic)
    return scale!(- sum(abs2, residuals(obj)) / dof_residual(obj), inv(jacobian(obj)))
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
