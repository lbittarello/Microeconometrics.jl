#==========================================================================================#

# TYPE

mutable struct IV <: GMM

    method::String
    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}
    W::Matrix{Float64}

    IV() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function IV(MD::Microdata, method::String)
    obj        = IV()
    obj.sample = MD
    obj.method = method
    return obj
end

#==========================================================================================#

# INTERFACE

function fit(::Type{IV}, MD::Microdata; novar::Bool = false, method::String = "TSLS")

    if method == "OLS"

        FSM               = Dict(:treatment => "", :instrument => "")
        FSD               = Microdata(MD, FSM)
        FSD.map[:control] = vcat(MD.map[:treatment], MD.map[:control])
        obj               = OLS(FSD)

        _fit!(obj, getweights(obj))

    elseif method == "Reduced form"

        FSM               = Dict(:treatment => "", :instrument => "")
        FSD               = Microdata(MD, FSM)
        FSD.map[:control] = vcat(MD.map[:instrument], MD.map[:control])
        obj               = OLS(FSD)

        _fit!(obj, getweights(obj))

    elseif method == "First stage"

        FSM                = Dict(:treatment => "", :instrument => "")
        FSD                = Microdata(MD, FSM)
        FSD.map[:response] = MD.map[:treatment]
        FSD.map[:control]  = vcat(MD.map[:instrument], MD.map[:control])
        obj                = OLS(FSD)

        _fit!(obj, getweights(obj))

    elseif length(MD.map[:treatment]) == length(MD.map[:instrument])

        obj   = IV(MD, "Method of moments")
        k     = length(MD.map[:instrument]) + length(MD.map[:control])
        obj.W = Matrix{Float64}(I, k, k)

        _fit!(obj, getweights(obj))

    elseif (method == "TSLS") | (method == "2SLS")

        obj   = IV(MD, "Two-step GMM")
        obj.W = crossprod(getmatrix(obj, :instrument, :control), getweights(obj))

        _fit!(obj, getweights(obj))

    elseif (method == "Two-step GMM") | (method == "Optimal GMM")

        obj   = IV(MD, method) ; _fit!(obj, getweights(obj))
        obj.W = wmatrix(obj, getcorr(obj), getweights(obj))

        _fit!(obj, obj.W, getweights(obj))

    else
        throw("unknown method")
    end

    novar || _vcov!(obj, getcorr(obj), getweights(obj))

    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::IV, ::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        obj.β = (z' * x) \ (z' * y)
    else
        zγ    = z * (z \ x)
        obj.β = (zγ' * x) \ (zγ' * y)
    end
end

function _fit!(obj::IV, W::Matrix{Float64}, ::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    zγ    = z * (W \ (z' * x))
    obj.β = (zγ' * x) \ (zγ' * y)
end

function _fit!(obj::IV, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = Diagonal(w) * z

    if obj.method == "Method of moments"
        obj.β = (v' * x) \ (v' * y)
    else
        γ     = (v' * z) \ (v' * x)
        vγ    = v * γ
        obj.β = (vγ' * x) \ (vγ' * y)
    end
end

function _fit!(obj::IV, W::Matrix{Float64}, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = Diagonal(w) * z

    vγ    = v * (W \ (v' * x))
    obj.β = (vγ' * x) \ (vγ' * y)
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::IV) = Diagonal(residuals(obj)) * getmatrix(obj, :instrument, :control)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV, ::UnitWeights)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    return - z' * x
end

function jacobian(obj::IV, w::AbstractWeights)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = Diagonal(w) * z
    return - v' * x
end

# VARIANCE MATRIX

function _vcov!(obj::IV, corr::Homoscedastic, ::UnitWeights)

    x  = getmatrix(obj, :treatment, :control)
    z  = getmatrix(obj, :instrument, :control)
    r  = residuals(obj)
    σ² = sum(abs2, r) / dof_residual(obj)
    γ  = z \ x
    V  = x' * z * γ

    obj.V = lmul!(σ², inv(V))
end

function _vcov!(obj::IV, corr::Homoscedastic, w::AbstractWeights)

    x  = getmatrix(obj, :treatment, :control)
    z  = getmatrix(obj, :instrument, :control)
    v  = Diagonal(w) * z
    r  = residuals(obj)
    σ² = sum(abs2.(r), w) / dof_residual(obj)
    γ  = (v' * z) \ (v' * x)
    V  = x' * v * γ

    obj.V = lmul!(σ², inv(V))
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::IV, MD::Microdata)
    if getnames(obj, :treatment, :control) != getnames(MD, :treatment, :control)
        throw("missing variables")
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
mtitle(obj::IV)    =  "Linear GMM"

# COEFFICIENT OF DETERMINATION

adjr2(obj::IV) = 1.0 - (1.0 - r2(obj)) * (nobs(obj) - 1) / dof_residual(obj)
r2(obj::IV)    = _r2(obj, getweights(obj))

function _r2(obj::IV, ::UnitWeights)
    y   = response(obj)
    ŷ   = fitted(obj)
    rss = sum(abs2, y .- ŷ)
    tss = sum(abs2, y .- mean(y))
    return 1.0 - rss / tss
end

function _r2(obj::IV, w::AbstractWeights)
    y   = response(obj)
    μ   = mean(y, w)
    ŷ   = fitted(obj)
    rss = sum(abs2.(y .- ŷ), w)
    tss = sum(abs2.(y .- μ), w)
    return 1.0 - rss / tss
end
