#==========================================================================================#

# TYPE

mutable struct Abadie <: TwoStageModel

    first_stage::Micromodel
    second_stage::ParModel
    pscore::Vector{Float64}
    eweights::ProbabilityWeights{Float64, Float64, Vector{Float64}}

    Abadie() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(::Type{Abadie}, M₁::Type{<:Micromodel}, MD::Microdata; kwargs...)

    FSM                = Dict(:treatment => "", :instrument => "")
    FSD                = Microdata(MD, FSM)
    FSD.map[:response] = MD.map[:instrument]

    return fit(M₁, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit(
        ::Type{Abadie},
        M₂::Type{<:Micromodel},
        M₁::Type{<:ParModel},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    )

    m₁ = first_stage(Abadie, M₁, MD, novar = novar)
    return fit(Abadie, M₂, m₁, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{Abadie},
        M₂::Type{<:Micromodel},
        M₁::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    )

    w = getweights(MD)
    d = getvector(MD, :treatment)
    z = getvector(MD, :instrument)
    π = fitted(M₁)
    v = [1.0 - (1.0 - di) * zi / πi - di * (1.0 - zi) / (1.0 - πi)
         for (di, zi, πi) in zip(d, z, π)]

    v[((trim .> π) .| (1.0 - trim .< π))] .= 0.0

    SSD               = Microdata(MD, Dict{Symbol,String}())
    SSD.map[:control] = vcat(SSD.map[:treatment], SSD.map[:control])
    obj               = Abadie()
    obj.first_stage   = M₁
    obj.second_stage  = M₂(SSD; kwargs...)
    obj.pscore        = π
    obj.eweights       = pweights(v)

    _fit!(second_stage(obj), reweight(w, obj.eweights))
    novar || _vcov!(obj, getcorr(obj), w)

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Abadie) = lmul!(Diagonal(obj.eweights), score(second_stage(obj)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::Abadie, ::UnitWeights) = jacobian(second_stage(obj), obj.eweights)

function jacobian(obj::Abadie, w::AbstractWeights)
    return jacobian(second_stage(obj), reweight(w, obj.eweights))
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::Abadie, ::UnitWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    D = [(1.0 - di) * zi / abs2(πi) - di * (1.0 - zi) / abs2(1.0 - πi)
         for (di, zi, πi) in zip(d, z, π)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

function crossjacobian(obj::Abadie, w::AbstractWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    D = [wi * ((1.0 - di) * zi / abs2(πi) - di * (1.0 - zi) / abs2(1.0 - πi))
         for (di, zi, πi, wi) in zip(d, z, π, w)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::Abadie) = predict(second_stage(obj))

# FITTED VALUES

fitted(obj::Abadie) = fitted(second_stage(obj))

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::Abadie) = jacobexp(second_stage(obj))

#==========================================================================================#

# UTILITIES

coefnames(obj::Abadie) = coefnames(second_stage(obj))
mtitle(obj::Abadie)    = "Abadie (2003)"
