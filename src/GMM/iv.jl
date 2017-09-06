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

function IV(MD::Microdata; method::String = "One-step")

    obj        = IV()
    obj.sample = MD

    if length(MD.map[:treatment]) == length(MD.map[:instrument])
        obj.method = "Method of moments"
    elseif (method == "One-step") | (method == "One-step GMM") | (method == "2SLS")
        obj.method = "One-step GMM"
    else
        throw("unknown method")
    end

    return obj
end

#==========================================================================================#

# INTERFACE

function fit(::Type{IV}, MD::Microdata; novar::Bool = false, method::String = "One-step")

    if method == "OLS"
        FSD               = Microdata(MD)
        FSD.map[:control] = vcat(FSD.map[:treatment], FSD.map[:control])
        pop!(FSD.map, :treatment)
        pop!(FSD.map, :instrument)
        obj = OLS(FSD)
    else
        obj = IV(MD, method = method)
    end

    if checkweight(MD)
        w = getvector(MD, :weight)
        _fit!(obj, w)
        novar || (obj.V = _vcov(obj, MD.corr, w))
    else
        _fit!(obj)
        novar || (obj.V = _vcov(obj, MD.corr))
    end

    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::IV)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :treatment, :control)
    z  = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "One-step GMM"
        zγ    = z * (z \ x)
        obj.β = zγ \ y
    end
end

function _fit!(obj::IV, w::AbstractVector)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = scale!(transpose(z), w)

    if obj.method == "Method of moments"
        obj.β = (v * x) \ (v * y)
    elseif obj.method == "One-step GMM"
        vγ    = transpose((v * z) \ (v * x)) * v
        obj.β = (vγ * x) \ (vγ * y)
    end
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

function score(obj::IV)

    z = getmatrix(obj, :instrument, :control)
    s = scale!(residuals(obj), copy(z))

    if obj.method == "Method of moments"
        return s
    elseif obj.method == "One-step GMM"
        x = getmatrix(obj, :treatment, :control)
        return s * (z \ x)
    end
end

function score(obj::IV, w::AbstractVector)

    z = getmatrix(obj, :instrument, :control)
    v = scale!(w, copy(z))
    s = scale!(residuals(obj), copy(v))

    if obj.method == "Method of moments"
        return s
    elseif obj.method == "One-step GMM"
        x = getmatrix(obj, :treatment, :control)
        return s * ((v' * z) \ (v' * x))
    end
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV)

    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        return scale!(- 1.0, z' * x)
    elseif obj.method == "One-step GMM"
        zγ = z * (z \ x)
        return scale!(- 1.0, zγ' * x)
    end
end

function jacobian(obj::IV, w::AbstractVector)

    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = scale!(transpose(z), w)

    if obj.method == "Method of moments"
        return scale!(- 1.0, v * x)
    elseif obj.method == "One-step GMM"
        vγ = transpose((v * z) \ (v * x)) * v
        return scale!(- 1.0, vγ * x)
    end
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::IV) = getmatrix(obj, :treatment, :control) * obj.β

# FITTED VALUES

fitted(obj::IV) = predict(obj)

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::IV) = getmatrix(obj, :treatment, :control)

#==========================================================================================#

# UTILITIES

coefnames(obj::IV) = getnames(obj, :treatment, :control)

# FIRST-STAGE (OLS, OUTCOME IS TREATMENT AND INSTRUMENTS ENTER AS CONTROLS)

function first_stage(::Type{IV}, MD::Microdata)
    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:treatment]
    FSD.map[:control]  = vcat(FSD.map[:instrument], FSD.map[:control])
    pop!(FSD.map, :treatment)
    pop!(FSD.map, :instrument)
    return fit(OLS, FSD)
end

first_stage(MM::IV) = first_stage(IV, MM.sample)

# REDUCED FORM (OLS, INSTRUMENTS REPLACE TREATMENTS)

function reduced_form(::Type{IV}, MD::Microdata)
    RF               = Microdata(MD)
    RF.map[:control] = vcat(FSD.map[:instrument], FSD.map[:control])
    pop!(RF.map, :treatment)
    pop!(RF.map, :instrument)
    return fit(OLS, RF)
end

reduced_form(MM::IV) = reduced_form(IV, MM.sample)
