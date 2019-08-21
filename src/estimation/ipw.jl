#==========================================================================================#

# TYPE

mutable struct IPW <: TwoStageModel

    first_stage::OneStageModel
    second_stage::OLS
    pscore::Vector{Float64}
    eweights::ProbabilityWeights{Float64, Float64, Vector{Float64}}

    IPW() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(::Type{IPW}, MM::Type{<:OneStageModel}, MD::Microdata; kwargs...)

    FSD                    = Microdata(MD)
    FSD.mapping[:response] = MD.mapping[:treatment]

    pop!(FSD.mapping, :treatment)

    return fit(MM, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit(
        ::Type{IPW}, MM::Type{<:OneStageModel}, MD::Microdata; novar::Bool = false, kwargs...
    )

    m = first_stage(IPW, MM, MD; novar = novar)
    return fit(IPW, m, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{IPW},
        MM::OneStageModel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
    )

    w = getweights(MD)
    d = getvector(MD, :treatment)
    π = fitted(MM)
    v = [(1.0 - di) / (1.0 - πi) + di / πi for (di, πi) in zip(d, π)]

    v[((trim .> π) .| (1.0 - trim .< π))] .= 0.0

    SSD                   = Microdata(MD)
    SSD.mapping[:control] = asgn(MD.model, (MD.model[:treatment], InterceptTerm{true}()))
    
    obj              = IPW()
    obj.first_stage  = MM
    obj.second_stage = OLS(SSD)
    obj.pscore       = π
    obj.eweights     = pweights(v)

    _fit!(second_stage(obj), reweight(w, obj.eweights))
    novar || _vcov!(obj, getcorr(obj), w)

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::IPW) = lmul!(Diagonal(obj.eweights), score(second_stage(obj)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::IPW, ::UnitWeights) = jacobian(second_stage(obj), obj.eweights)

function jacobian(obj::IPW, w::AbstractWeights)
    return jacobian(second_stage(obj), reweight(w, obj.eweights))
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::IPW, ::UnitWeights)

    d = getvector(obj, :treatment)
    π = obj.pscore
    D = [(1.0 - di) / abs2(1.0 - πi) - di / abs2(πi) for (di, πi) in zip(d, π)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

function crossjacobian(obj::IPW, w::AbstractWeights)

    d = getvector(obj, :treatment)
    π = obj.pscore
    D = [wi * ((1.0 - di) / abs2(1.0 - πi) - di / abs2(πi)) for (di, πi, wi)
         in zip(d, π, w)]

    D[iszero.(obj.eweights)] .= 0.0

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * lmul!(Diagonal(D), g₁)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::IPW) = predict(second_stage(obj))

# FITTED VALUES

fitted(obj::IPW) = fitted(second_stage(obj))

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::IPW) = jacobexp(second_stage(obj))

#==========================================================================================#

# UTILITIES

coefnames(obj::IPW) = coefnames(second_stage(obj))
mtitle(obj::IPW)    = "Inverse probability weighting"
