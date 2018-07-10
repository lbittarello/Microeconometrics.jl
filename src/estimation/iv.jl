#==========================================================================================#

# TYPE

mutable struct IV <: ParModel

    method::String
    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    IV() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function IV(MD::Microdata; method::String = "TSLS")

    obj        = IV()
    obj.sample = MD

    if length(MD.map[:treatment]) == length(MD.map[:instrument])
        obj.method = "Method of moments"
    elseif (method == "TSLS") | (method == "2SLS")
        obj.method = "TSLS"
    else
        throw("unknown method")
    end

    return obj
end

#==========================================================================================#

# INTERFACE

function fit(::Type{IV}, MD::Microdata; novar::Bool = false, method::String = "TSLS")

    if method == "OLS"
        FSD               = Microdata(MD)
        FSD.map[:control] = vcat(FSD.map[:treatment], FSD.map[:control])
        pop!(FSD.map, :treatment)
        pop!(FSD.map, :instrument)
        obj = OLS(FSD)
    elseif method == "Reduced form"
        FSD               = Microdata(MD)
        FSD.map[:control] = vcat(FSD.map[:instrument], FSD.map[:control])
        pop!(FSD.map, :treatment)
        pop!(FSD.map, :instrument)
        obj = OLS(FSD)
    elseif method == "First stage"
        FSD                = Microdata(MD)
        FSD.map[:response] = FSD.map[:treatment]
        FSD.map[:control]  = vcat(FSD.map[:instrument], FSD.map[:control])
        pop!(FSD.map, :treatment)
        pop!(FSD.map, :instrument)
        obj = OLS(FSD)
    else
        obj = IV(MD, method = method)
    end

    _fit!(obj, getweights(obj))
    novar || _vcov!(obj, getcorr(obj), getweights(obj))

    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::IV, w::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "TSLS"
        γ     = z \ x
        zγ    = z * γ
        obj.β = zγ \ y
    end
end

function _fit!(obj::IV, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = scale!(values(w), copy(z))

    if obj.method == "Method of moments"
        obj.β = (v' * x) \ (v' * y)
    elseif obj.method == "TSLS"
        γ     = (v' * x) \ (v' * y)
        vγ    = v * γ
        obj.β = (vγ' * x) \ (vγ' * y)
    end
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

function score(obj::IV)

    z = getmatrix(obj, :instrument, :control)
    s = scale!(residuals(obj), copy(z))

    if obj.method == "Method of moments"
        return s
    elseif obj.method == "TSLS"
        x = getmatrix(obj, :treatment, :control)
        γ = z \ x
        return s * γ
    end
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV, w::UnitWeights)

    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        return scale!(- 1.0, z' * x)
    elseif obj.method == "TSLS"
        γ  = z \ x
        zγ = z * γ
        return scale!(- 1.0, zγ' * x)
    end
end

function jacobian(obj::IV, w::AbstractWeights)

    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = scale!(values(w), copy(z))

    if obj.method == "Method of moments"
        return scale!(- 1.0, v' * x)
    elseif obj.method == "TSLS"
        γ  = (v' * z) \ (v' * x)
        vγ = v * γ
        return scale!(- 1.0, vγ' * x)
    end
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::IV, MD::Microdata)

    if getnames(obj, :treatment, :control) != getnames(MD, :treatment, :control)
        throw("some variables are missing")
    end

    getmatrix(MD, :treatment, :control) * obj.β
end

# FITTED VALUES

fitted(obj::IV, MD::Microdata) = predict(obj, MD)

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::IV) = getmatrix(obj, :treatment, :control)

#==========================================================================================#

# UTILITIES

coefnames(obj::IV) = getnames(obj, :treatment, :control)
adjr2(obj::IV)     = 1.0 - (1.0 - r2(obj)) * (nobs(obj) - 1) / dof_residual(obj)
r2(obj::IV)        = _r2(obj, getweights(obj))

function _r2(obj::IV, w::UnitWeights)
    y   = model_response(obj)
    ŷ   = fitted(obj)
    rss = sum(abs2, y .- ŷ)
    tss = sum(abs2, y .- mean(y))
    return 1.0 - rss / tss
end

function _r2(obj::IV, w::AbstractWeights)
    y   = model_response(obj)
    μ   = mean(y, w)
    ŷ   = fitted(obj)
    rss = sum(abs2.(y .- ŷ), w)
    tss = sum(abs2.(y .- μ), w)
    return 1.0 - rss / tss
end
