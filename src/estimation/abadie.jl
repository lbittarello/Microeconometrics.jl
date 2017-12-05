#==========================================================================================#

# TYPE

mutable struct Abadie <: TwoStageModel

    first_stage::Micromodel
    second_stage::ParModel
    pscore::Vector{Float64}
    weights::PWeights

    Abadie() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(
        ::Type{Abadie}, ::Type{M}, MD::Microdata; kwargs...
    ) where {M <: Micromodel}

    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:instrument]
    pop!(FSD.map, :treatment)
    pop!(FSD.map, :instrument)
    return fit(M, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit(
        ::Type{Abadie},
        ::Type{M₂},
        ::Type{M₁},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    ) where {M₂ <: Micromodel, M₁ <: ParModel}

    m₁ = first_stage(Abadie, M₁, MD, novar = novar)
    return fit(Abadie, M₂, m₁, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{Abadie},
        ::Type{M},
        MM::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    ) where {M <: Micromodel}

    invtrim = one(trim) - trim
    d       = getvector(MD, :treatment)
    z       = getvector(MD, :instrument)
    π       = fitted(MM)
    v       = fill(0.0, size(MD, 1))

    @inbounds for (i, (di, zi, πi)) in enumerate(zip(d, z, π))
        d0 = iszero(di)
        z0 = iszero(zi)
        if trim <= πi <= invtrim
            if d0 & z0
                v[i] = 1.0
            elseif d0 & !z0
                v[i] = 1.0 - 1.0 / πi
            elseif !d0 & z0
                v[i] = 1.0 - 1.0 / (1.0 - πi)
            elseif !d0 & !z0
                v[i] = 1.0
            end
        end
    end

    SSD               = Microdata(MD)
    SSD.map[:control] = vcat(SSD.map[:treatment], SSD.map[:control])
    obj               = Abadie()
    obj.first_stage   = MM
    obj.second_stage  = M(SSD; kwargs...)
    obj.pscore        = π
    obj.weights       = pweights(v)

    _fit!(second_stage(obj), reweight(w, obj.weights))
    novar || _vcov!(obj, getcorr(obj), w)

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Abadie) = scale!(obj.weights, score(second_stage(obj)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::Abadie, w::UnitWeights) = jacobian(second_stage(obj), obj.weights)

function jacobian(obj::Abadie, w::AbstractWeights)
    return jacobian(second_stage(obj), reweight(w, obj.weights))
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::Abadie, w::UnitWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    v = obj.weights
    D = fill(0.0, length(v))

    @inbounds for (i, (di, zi, πi, vi)) in enumerate(zip(d, z, π, v))
        if !iszero(vi)
            d0 = iszero(di)
            z0 = iszero(zi)
            if d0 & z0
                D[i] = 0.0
            elseif d0 & !z0
                D[i] = 1.0 / abs2(πi)
            elseif !d0 & z0
                D[i] = - 1.0 / abs2(1.0 - πi)
            elseif !d0 & !z0
                D[i] = 0.0
            end
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
end

function crossjacobian(obj::Abadie, w::AbstractWeights)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.pscore
    v = obj.weights
    D = fill(0.0, length(v))

    @inbounds for (i, (di, zi, πi, vi, wi)) in enumerate(zip(d, z, π, v, values(w)))
        if !iszero(vi)
            d0 = iszero(di)
            z0 = iszero(zi)
            if d0 & z0
                D[i] = 0.0
            elseif d0 & !z0
                D[i] = wi / abs2(πi)
            elseif !d0 & z0
                D[i] = - wi / abs2(1.0 - πi)
            elseif !d0 & !z0
                D[i] = 0.0
            end
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
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
