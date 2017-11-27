#==========================================================================================#

# TYPE

mutable struct OLS{T} <: ParModel{T}

    sample::Microdata{T}
    β::Vector{Float64}
    V::Matrix{Float64}

    OLS{T}() where {T} = new()
end

#==========================================================================================#

# CONSTRUCTOR

function OLS(MD::Microdata{T}) where {T}
    obj        = OLS{T}()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::OLS)
    y     = getvector(obj, :response)
    x     = getmatrix(obj, :control)
    obj.β =  x \ y
end

function _fit!(obj::OLS, w::AbstractVector)
    y     = getvector(obj, :response)
    x     = getmatrix(obj, :control)
    v     = scale!(w, copy(x))
    obj.β =  (v' * x) \ (v' * y)
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::OLS) = scale!(residuals(obj), copy(getmatrix(obj, :control)))

function score(obj::OLS, w::AbstractVector)
    return scale!(w, scale!(residuals(obj), copy(getmatrix(obj, :control))))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::OLS)                    = crossprod(getmatrix(obj, :control), neg = true)
jacobian(obj::OLS, w::AbstractVector) = crossprod(getmatrix(obj, :control), w, neg = true)

# HOMOSCEDASTIC VARIANCE MATRIX

function _vcov!(obj::OLS{Homoscedastic})
    obj.V = scale!(- sum(abs2, residuals(obj)) / dof_residual(obj), inv(jacobian(obj)))
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
