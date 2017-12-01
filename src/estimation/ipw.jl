#==========================================================================================#

# TYPE

mutable struct IPW{T} <: TwoStageModel{T}

    first_stage::Micromodel{T}
    second_stage::OLS{T}
    mat::Matrix{Float64}

    IPW{T}() where {T} = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage(
        ::Type{IPW}, ::Type{M}, MD::Microdata; kwargs...
    ) where {M <: Micromodel{T} where T}

    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:treatment]
    pop!(FSD.map, :treatment)

    return fit(M, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit(
        ::Type{IPW},
        ::Type{M},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    ) where {M <: Micromodel{T} where T}

    m = first_stage(IPW, M, MD, novar = novar)
    return fit(IPW, m, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{IPW},
        MM::Micromodel{T},
        MD::Microdata{T};
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    ) where {T}

    SSD               = Microdata(MD)
    SSD.map[:control] = vcat(SSD.map[:treatment], 1)
    obj               = IPW{T}()
    obj.first_stage   = MM
    obj.second_stage  = OLS(SSD)

    invtrim = one(trim) - trim
    t       = getvector(SSD, :treatment)
    π       = fitted(MM)

    obj.mat       = Matrix{Float64}(nobs(SSD), 2)
    obj.mat[:, 1] = π
    obj.mat[:, 2] = 0.0

    @inbounds for (i, (ti, πi)) in enumerate(zip(t, π))
        if trim <= πi <= invtrim
            obj.mat[i, 2] = (iszero(ti) ? (1.0 / (1.0 - πi)) : (1.0 / πi))
        end
    end

    if checkweight(SSD)
        w = getvector(SSD, :weight)
        _fit!(second_stage(obj), w .* obj.mat[:, 2])
        novar || _vcov!(obj, w)
    else
        _fit!(second_stage(obj), obj.mat[:, 2])
        novar || _vcov!(obj)
    end

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::IPW)                    = score(second_stage(obj), obj.mat[:, 2])
score(obj::IPW, w::AbstractVector) = score(second_stage(obj), w .* obj.mat[:, 2])

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::IPW)                    = jacobian(second_stage(obj), obj.mat[:, 2])
jacobian(obj::IPW, w::AbstractVector) = jacobian(second_stage(obj), w .* obj.mat[:, 2])

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::IPW)

    t = getvector(obj, :treatment)
    π = obj.mat[:, 1]
    v = obj.mat[:, 2]
    D = fill(0.0, nobs(obj))

    @inbounds for (i, (ti, πi, vi)) in enumerate(zip(t, π, v))
        if !iszero(vi)
            D[i] = (iszero(ti) ? (1.0 / abs2(1.0 - πi)) : (- 1.0 / abs2(πi)))
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
end

function crossjacobian(obj::IPW, w::AbstractVector)

    t = getvector(obj, :treatment)
    π = obj.mat[:, 1]
    v = obj.mat[:, 2]
    D = fill(0.0, nobs(obj))

    @inbounds for (i, (ti, πi, vi, wi)) in enumerate(zip(t, π, v, w))
        if !iszero(vi)
            D[i] = (iszero(ti) ? (wi / abs2(1.0 - πi)) : (- wi / abs2(πi)))
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
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
