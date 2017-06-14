#==========================================================================================#

# TYPE

mutable struct IV <: GMM

    method::String
    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    IV() = new()
end

#==========================================================================================#

# INTERFACE

function fit(::Type{IV}, MD::Microdata; method::String = "Unadjusted", novar::Bool = false)

    output        = IV()
    output.method = method
    output.sample = MD

    if checkweight(MD)
        w = getvector(MD, :weight)
        output.β = _fit(output, w)
        novar || (output.V = _vcov(output, MD.corr, w))
    else
        output.β = _fit(output)
        novar || (output.V = _vcov(output, MD.corr))
    end

    return output
end

#==========================================================================================#

# ESTIMATION

function _fit(obj::IV)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if size(x, 2) == size(z, 2)
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "Unadjusted"
        obj.β = (z * (z \ x)) \ y
    elseif obj.method == "Optimal"
        obj.β = (z * (z \ x)) \ y
        wmat  = _opg(obj, getcorr(obj))
        zw    = wmat \ z
        obj.β = (zw' * x) \ (zw' * y)
    else
        throw("unknown method")
    end
end

function _fit(obj::IV, w::AbstractVector)

    y = getvector(obj, :response) .* w
    x = getmatrix(obj, :treatment, :control) .* w
    z = getmatrix(obj, :instrument, :control)

    if size(x, 2) == size(z, 2)
        obj.β = (z' * x) \ (z' * y)
    elseif obj.method == "Unadjusted"
        obj.β = (z * (z \ x)) \ y
    elseif obj.method == "Optimal"
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

score(obj::IV) = getmatrix(obj, :instrument, :control) .* residuals(obj)

function score(obj::IV, w::AbstractVector)
    z  = getmatrix(obj, :instrument, :control)
    r  = residuals(obj)
    r .= r .* w
    return z .* r
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::IV)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)
    return - z' * x
end

function jacobian(obj::IV, w::AbstractVector)
    x = getmatrix(obj, :treatment, :control) .* w
    z = getmatrix(obj, :instrument, :control)
    return - z' * x
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
