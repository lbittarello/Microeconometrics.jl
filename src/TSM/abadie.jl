#==========================================================================================#

# TYPE

mutable struct Abadie <: TwoStageModel

    first_stage::Micromodel
    second_stage::ParModel
    mat::Matrix{Float64}

    Abadie() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage{M <: Micromodel}(::Type{Abadie}, ::Type{M}, MD::Microdata; kwargs...)
    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:instrument]
    pop!(FSD.map, :treatment)
    pop!(FSD.map, :instrument)
    return fit(M, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit{M₁ <: Micromodel, M₂ <: ParModel}(
        ::Type{Abadie},
        ::Type{M₂},
        ::Type{M₁},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    )

    m₁ = first_stage(Abadie, M₁, MD, novar = novar)
    return fit(Abadie, M₂, m₁, MD; novar = novar, kwargs...)
end

function fit{M <: Micromodel}(
        ::Type{Abadie},
        ::Type{M},
        MM::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    )

    SSD               = Microdata(MD)
    SSD.map[:control] = vcat(SSD.map[:treatment], SSD.map[:control])
    obj               = Abadie()
    obj.first_stage   = MM
    obj.second_stage  = M(SSD; kwargs...)

    invtrim = one(trim) - trim
    d       = getvector(SSD, :treatment)
    z       = getvector(SSD, :instrument)
    π       = fitted(MM)

    obj.mat       = Matrix{Float64}(nobs(SSD), 2)
    obj.mat[:, 1] = π
    obj.mat[:, 2] = 0.0

    @inbounds for (i, (di, zi, πi)) in enumerate(zip(d, z, π))
        d0 = iszero(di)
        z0 = iszero(zi)
        if trim <= πi <= invtrim
            if d0 & z0
                obj.mat[i, 2] = 1.0
            elseif d0 & !z0
                obj.mat[i, 2] = 1.0 - 1.0 / πi
            elseif !d0 & z0
                obj.mat[i, 2] = 1.0 - 1.0 / (1.0 - πi)
            elseif !d0 & !z0
                obj.mat[i, 2] = 1.0
            end
        end
    end

    if checkweight(SSD)
        w = getvector(SSD, :weight)
        _fit!(second_stage(obj), w .* obj.mat[:, 2])
        novar || (obj.second_stage.V = _vcov(obj, SSD.corr, w))
    else
        _fit!(second_stage(obj), obj.mat[:, 2])
        novar || (obj.second_stage.V = _vcov(obj, SSD.corr))
    end

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Abadie)                    = score(second_stage(obj), obj.mat[:, 2])
score(obj::Abadie, w::AbstractVector) = score(second_stage(obj), w .* obj.mat[:, 2])

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::Abadie)                    = jacobian(second_stage(obj), obj.mat[:, 2])
jacobian(obj::Abadie, w::AbstractVector) = jacobian(second_stage(obj), w .* obj.mat[:, 2])

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::Abadie)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.mat[:, 1]
    v = obj.mat[:, 2]
    D = zeros(nobs(obj))

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

function crossjacobian(obj::Abadie, w::AbstractVector)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.mat[:, 1]
    v = obj.mat[:, 2]
    D = zeros(π)

    @inbounds for (i, (di, zi, πi, vi, wi)) in enumerate(zip(d, z, π, v, w))
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
