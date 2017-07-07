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

    if method == "OLS"
        obj.method = "OLS"
    elseif length(MD.map[:treatment]) == length(MD.map[:instrument])
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
        wmat  = _opg(obj, getcorr(obj))
        zw    = wmat \ z
        obj.β = (zw' * x) \ (zw' * y)
    else
        throw("unknown method")
    end
end

function _fit(obj::IV, w::AbstractVector)

    y = scale!(w, copy(getvector(obj, :response)))
    x = scale!(w, copy(getmatrix(obj, :treatment, :control)))
    z = getmatrix(obj, :instrument, :control)

    if obj.method == "Method of moments"
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "Unadjusted GMM"
        obj.β = (z * (z \ x)) \ y
    elseif obj.method == "Optimal GMM"
        obj.β = (z * (z \ x)) \ y
        wmat  = _opg(obj, getcorr(obj), w)
        zw    = wmat \ z
        obj.β = (zw' * x) \ (zw' * y)
    else
        throw("unknown method")
    end
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::IV) = scale!(residuals(obj), copy(getmatrix(obj, :instrument, :control)))

function score(obj::IV, w::AbstractVector)
    z  = copy(getmatrix(obj, :instrument, :control))
    r  = residuals(obj)
    r .= r .* w
    return scale!(w, scale!(r, z))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    return scale!(- 1.0, z' * x)
end

function jacobian(obj::IV, w::AbstractVector)
    x = scale!(w, copy(getmatrix(obj, :treatment, :control)))
    z = getmatrix(obj, :instrument, :control)
    return scale!(- 1.0, z' * x)
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
    FSD                  = Microdata(MD)
    FSD.map[:response]   = FSD.map[:treatment]
    FSD.map[:control]    = vcat(FSD.map[:instrument], FSD.map[:control])
    FSD.map[:treatment]  = [0]
    FSD.map[:instrument] = [0]
    return fit(OLS, FSD)
end

first_stage(MM::IV) = first_stage(IV, MM.sample)

# REDUCED FORM (OLS, INSTRUMENTS REPLACE TREATMENTS)

function reduced_form(::Type{IV}, MD::Microdata)
    RF                  = Microdata(MD)
    RF.map[:control]    = vcat(FSD.map[:instrument], FSD.map[:control])
    RF.map[:treatment]  = [0]
    RF.map[:instrument] = [0]
    return fit(OLS, RF)
end

reduced_form(MM::IV) = reduced_form(IV, MM.sample)
