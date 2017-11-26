#==========================================================================================#

# TYPE

mutable struct FrölichMelly <: TwoStageModel

    first_stage::Micromodel
    second_stage::OLS
    mat::Matrix{Float64}

    FrölichMelly() = new()
end

#==========================================================================================#

# FIRST STAGE

function first_stage{M <: Micromodel}(
        ::Type{FrölichMelly}, ::Type{M}, MD::Microdata; kwargs...
    )

    FSD                = Microdata(MD)
    FSD.map[:response] = FSD.map[:instrument]
    pop!(FSD.map, :treatment)
    pop!(FSD.map, :instrument)

    return fit(M, FSD; kwargs...)
end

#==========================================================================================#

# ESTIMATION

function fit{M <: Micromodel}(
        ::Type{FrölichMelly},
        ::Type{M},
        MD::Microdata;
        novar::Bool = false,
        kwargs...
    )

    m = first_stage(FrölichMelly, M, MD, novar = novar)
    return fit(FrölichMelly, m, MD; novar = novar, kwargs...)
end

function fit(
        ::Type{FrölichMelly},
        MM::Micromodel,
        MD::Microdata;
        novar::Bool = false,
        trim::AbstractFloat = 0.0,
        kwargs...
    )

    SSD               = Microdata(MD)
    SSD.map[:control] = vcat(SSD.map[:treatment], 1)
    obj               = FrölichMelly()
    obj.first_stage   = MM
    obj.second_stage  = OLS(SSD)

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
                obj.mat[i, 2] = 1.0 / (1.0 - πi)
            elseif d0 & !z0
                obj.mat[i, 2] = - 1.0 / πi
            elseif !d0 & z0
                obj.mat[i, 2] = - 1.0 / (1.0 - πi)
            elseif !d0 & !z0
                obj.mat[i, 2] = 1.0 / πi
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

score(obj::FrölichMelly)                    = score(second_stage(obj), obj.mat[:, 2])
score(obj::FrölichMelly, w::AbstractVector) = score(second_stage(obj), w .* obj.mat[:, 2])

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::FrölichMelly) = jacobian(second_stage(obj), obj.mat[:, 2])

function jacobian(obj::FrölichMelly, w::AbstractVector)
    jacobian(second_stage(obj), w .* obj.mat[:, 2])
end

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::FrölichMelly)

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

function crossjacobian(obj::FrölichMelly, w::AbstractVector)

    d = getvector(obj, :treatment)
    z = getvector(obj, :instrument)
    π = obj.mat[:, 1]
    v = obj.mat[:, 2]
    D = zeros(nobs(obj))

    @inbounds for (i, (di, zi, πi, vi, wi)) in enumerate(zip(d, z, π, v, w))
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
