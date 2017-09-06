#==========================================================================================#

# TYPE

mutable struct Tan <: TwoStageModel

    first_stage::Micromodel
    second_stage::IV
    mat::AbstractVector

    Tan() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage{M <: Micromodel}(::Type{Tan}, ::Type{M}, MD::Microdata; kwargs...)
    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:instrument]
    pop!(FSD.map, :treatment)
    pop!(FSD.map, :instrument)
    return fit(M, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit{M <: Micromodel}(
        ::Type{Tan},
        ::Type{M},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    )

    m = first_stage(Tan, M, MD, novar = novar)
    return fit(Tan, m, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{Tan},
        m::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    )

    SSD              = Microdata(MD, control = "1")
    obj              = Tan()
    obj.first_stage  = m
    obj.second_stage = IV(SSD)

    invtrim = one(trim) - trim
    z       = getvector(SSD, :instrument)
    p       = mean(z)
    π       = fitted(obj.first_stage)
    obj.mat = zeros(π)

    @inbounds for (i, (zi, πi)) in enumerate(zip(z, π))
        if trim < πi < invtrim
            obj.mat[i] = (iszero(zi) ? ((1.0 - p) / (1.0 - πi)) : (p / πi))
        end
    end

    if checkweight(SSD)
        w = getvector(SSD, :weight)
        _fit!(second_stage(obj), w .* obj.mat)
        novar || (obj.second_stage.V = _vcov(obj, SSD.corr, w))
    else
        _fit!(second_stage(obj), obj.mat)
        novar || (obj.second_stage.V = _vcov(obj, SSD.corr))
    end

    return obj
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Tan)                    = score(second_stage(obj), obj.mat)
score(obj::Tan, w::AbstractVector) = score(second_stage(obj), w .* obj.mat)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::Tan)                    = jacobian(second_stage(obj), obj.mat)
jacobian(obj::Tan, w::AbstractVector) = jacobian(second_stage(obj), w .* obj.mat)

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::Tan)

    z = getvector(obj, :instrument)
    p = mean(z)
    π = fitted(obj.first_stage)
    D = zeros(π)

    @inbounds for (i, (zi, πi, vi)) in enumerate(zip(z, π, obj.mat))
        if !iszero(vi)
            D[i] = (iszero(zi) ? ((1.0 - p) / abs2(1.0 - πi)) : (- p / abs2(πi)))
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
end

function crossjacobian(obj::Tan, w::AbstractVector)

    z = getvector(obj, :instrument)
    π = fitted(obj.first_stage)
    D = zeros(π)

    @inbounds for (i, (zi, πi, vi, wi)) in enumerate(zip(z, π, obj.mat, w))
        if !iszero(vi)
            D[i] = (iszero(zi) ? (wi * (1.0 - p) / abs2(1.0 - πi)) : (- wi * p / abs2(πi)))
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * scale!(D, g₁)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::Tan) = predict(second_stage(obj))

# FITTED VALUES

fitted(obj::Tan) = fitted(second_stage(obj))

# DERIVATIVE OF FITTED VALUES

jacobexp(obj::Tan) = jacobexp(second_stage(obj))

#==========================================================================================#

# UTILITIES

coefnames(obj::Tan) = coefnames(second_stage(obj))
