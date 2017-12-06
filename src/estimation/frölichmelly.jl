#==========================================================================================#

# TYPE

mutable struct FrölichMelly <: TwoStageModel

    first_stage::Micromodel
    second_stage::OLS
    pscore::Vector{Float64}
    weights::PWeights

    FrölichMelly() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(
        ::Type{FrölichMelly}, ::Type{M}, MD::Microdata; kwargs...
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
        ::Type{FrölichMelly},
        ::Type{M},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    ) where {M <: Micromodel}

    m = first_stage(FrölichMelly, M, MD; novar = novar, kwargs...)
    return fit(FrölichMelly, m, MD; novar = novar)
end

function fit(
        ::Type{FrölichMelly},
        MM::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
    )

    invtrim = one(trim) - trim
    w       = getweights(MD)
    d       = getvector(MD, :treatment)
    z       = getvector(MD, :instrument)
    π       = fitted(MM)
    v       = fill(0.0, size(MD, 1))

    @inbounds for (i, (di, zi, πi)) in enumerate(zip(d, z, π))
        d0 = iszero(di)
        z0 = iszero(zi)
        if trim <= πi <= invtrim
            if d0 & z0
                v[i] = 1.0 / (1.0 - πi)
            elseif d0 & !z0
                v[i] = - 1.0 / πi
            elseif !d0 & z0
                v[i] = - 1.0 / (1.0 - πi)
            elseif !d0 & !z0
                v[i] = 1.0 / πi
            end
        end
    end

    SSD               = Microdata(MD)
    SSD.map[:control] = vcat(SSD.map[:treatment], 1)
    obj               = FrölichMelly()
    obj.first_stage   = MM
    obj.second_stage  = OLS(SSD)
    obj.pscore        = π
    obj.weights       = pweights(v)

    _fit!(second_stage(obj), reweight(w, obj.weights))
    novar || _vcov!(obj, getcorr(obj), w)

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::FrölichMelly) = scale!(obj.weights, score(second_stage(obj)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::FrölichMelly, w::UnitWeights) = jacobian(second_stage(obj), obj.weights)

function jacobian(obj::FrölichMelly, w::AbstractWeights)
    return jacobian(second_stage(obj), reweight(w, obj.weights))
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::FrölichMelly, w::UnitWeights)

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
                D[i] = 1.0 / abs2(1.0 - πi)
            elseif d0 & !z0
                D[i] = 1.0 / abs2(πi)
            elseif !d0 & z0
                D[i] = - 1.0 / abs2(1.0 - πi)
            elseif !d0 & !z0
                D[i] = - 1.0 / abs2(πi)
            end
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
end

function crossjacobian(obj::FrölichMelly, w::AbstractWeights)

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
                D[i] = wi / abs2(1.0 - πi)
            elseif d0 & !z0
                D[i] = wi / abs2(πi)
            elseif !d0 & z0
                D[i] = - wi / abs2(1.0 - πi)
            elseif !d0 & !z0
                D[i] = - wi / abs2(πi)
            end
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
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
