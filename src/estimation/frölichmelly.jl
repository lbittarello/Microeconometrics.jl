#==========================================================================================#

# TYPE

mutable struct FrölichMelly <: TwoStageModel

    first_stage::Micromodel
    second_stage::OLS
    pscore::Vector{Float64}
    eweights::ProbabilityWeights{Float64, Float64, Vector{Float64}}

    FrölichMelly() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(::Type{FrölichMelly}, MM::Type{<:Micromodel}, MD::Microdata; kwargs...)

    FSM                = Dict(:treatment => "", :instrument => "")
    FSD                = Microdata(MD, FSM)
    FSD.map[:response] = MD.map[:instrument]

    return fit(MM, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit(
        ::Type{FrölichMelly},
        MM::Type{<:Micromodel},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    )

    m = first_stage(FrölichMelly, MM, MD; novar = novar)
    return fit(FrölichMelly, m, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{FrölichMelly},
        MM::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
    )

    w = getweights(MD)
    d = getvector(MD, :treatment)
    z = getvector(MD, :instrument)
    π = fitted(MM)
    v = [(2.0 * di - 1.0) * (zi - πi) / (πi * (1.0 - πi)) for (di, zi, πi) in zip(d, z, π)]

    v[((trim .> π) .| (1.0 - trim .< π))] .= 0.0

    SSD               = Microdata(MD, Dict{Symbol,String}())
    SSD.map[:control] = vcat(SSD.map[:treatment], 1)
    obj               = FrölichMelly()
    obj.first_stage   = MM
    obj.second_stage  = OLS(SSD)
    obj.pscore        = π
    obj.eweights       = pweights(v)

    _fit!(second_stage(obj), reweight(w, obj.eweights))
    novar || _vcov!(obj, getcorr(obj), w)

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

function score(obj::FrölichMelly)
    return lmul!(Diagonal(obj.eweights), score(second_stage(obj)))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::FrölichMelly, ::UnitWeights) = jacobian(second_stage(obj), obj.eweights)

function jacobian(obj::FrölichMelly, w::AbstractWeights)
    return jacobian(second_stage(obj), reweight(w, obj.eweights))
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::FrölichMelly, ::UnitWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    D = [- (2.0 * di - 1.0) * (zi / abs2(πi) + (1.0 - zi) / abs2(1.0 - πi))
         for (di, zi, πi) in zip(d, z, π)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

function crossjacobian(obj::FrölichMelly, w::AbstractWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    D = [- wi * (2.0 * di - 1.0) * (zi / abs2(πi) + (1.0 - zi) / abs2(1.0 - πi))
         for (di, zi, πi, wi) in zip(d, z, π, w)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::FrölichMelly) = predict(second_stage(obj))

# FITTED VALUES

fitted(obj::FrölichMelly) = fitted(second_stage(obj))

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::FrölichMelly) = jacobexp(second_stage(obj))

#==========================================================================================#

# UTILITIES

coefnames(obj::FrölichMelly) = coefnames(second_stage(obj))
mtitle(obj::FrölichMelly)    = "Frölich and Melly (2013)"
