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

function IV(MD::Microdata; method::String = "Unadjusted")

    obj        = IV()
    obj.sample = MD

    if length(MD.map[:treatment]) == length(MD.map[:instrument])
        obj.method = "Method of moments"
    elseif (method == "Unadjusted") | (method == "Unadjusted GMM") | (method == "2SLS")
        obj.method = "Unadjusted GMM"
    elseif (method == "Optimal") | (method == "Optimal GMM")
        obj.method = "Unadjusted GMM"
    else
        throw("unknown method; choose between Unadjusted and Optimal")
    end

    return obj
end

#==========================================================================================#

# INTERFACE

function fit(::Type{IV}, MD::Microdata; novar::Bool = false, method::String = "Unadjusted")

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
        obj.β = _fit(obj, w)
        novar || (obj.V = _vcov(obj, MD.corr, w))
    else
        obj.β = _fit(obj)
        novar || (obj.V = _vcov(obj, MD.corr))
    end

    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit(obj::IV)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "Unadjusted GMM"
        obj.β = (z * (z \ x)) \ y
    elseif obj.method == "Optimal GMM"
        obj.β = (z * (z \ x)) \ y
        W     = _opg(obj, getcorr(obj))
        zz    = W \ z
        obj.β = (zz' * x) \ (zz' * y)
    end
end

function _fit(obj::IV, w::AbstractVector)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    v = scale!(transpose(z), w)

    if obj.method == "Method of moments"
        obj.β = (v * x) \ (v * y)
    elseif obj.method == "Unadjusted GMM"
        obj.β = (z * (v * z) \ (v * x)) \ y
    elseif obj.method == "Optimal GMM"
        obj.β = (z * (v * z) \ (v * x)) \ y
        W     = _opg(obj, getcorr(obj), w)
        vv    = v / W
        obj.β = (vv * x) \ (vv * y)
    end
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::IV) = scale!(residuals(obj), copy(getmatrix(obj, :instrument, :control)))

function score(obj::IV, w::AbstractVector)
    return scale!(w, scale!(residuals(obj), copy(getmatrix(obj, :instrument, :control))))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    return scale!(- 1.0, z' * x)
end

function jacobian(obj::IV, w::AbstractVector)
    x = getmatrix(obj, :treatment, :control)
    v = scale!(transpose(getmatrix(obj, :instrument, :control)), w)
    return scale!(- 1.0, v * x)
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
