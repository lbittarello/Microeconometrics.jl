#==========================================================================================#

# TYPE

mutable struct Abadie <: TwoStageModel

    first_stage::Micromodel
    second_stage::ParModel
    mat::Vector{Float64}

    Abadie() = new()
end

#==========================================================================================#

# ESTIMATION

function fit{M₁ <: Micromodel, M₂ <: ParModel}(
        ::Type{Abadie},
        ::Type{M₁},
        ::Type{M₂},
        MD::Microdata;
        trim::AbstractFloat = 0.01,
        novar::Bool = false
    )

    invtrim = one(trim) - trim

    FSD                  = Microdata(MD)
    FSD.map[:response]   = FSD.map[:instrument]
    FSD.map[:instrument] = [0]
    FSD.map[:treatment]  = [0]
    SSD                  = Microdata(MD)
    SSD.map[:control]    = vcat(SSD.map[:treatment], SSD.map[:control])

    output                     = Abadie()
    output.first_stage         = fit(M₁, FSD)
    output.second_stage        = M₂()
    output.second_stage.sample = SSD

    d          = getvector(SSD, :treatment)
    z          = getvector(SSD, :instrument)
    π          = fitted(output.first_stage)
    output.mat = zeros(π)

    @inbounds for (i, (di, zi, πi)) in enumerate(zip(d, z, π))
        d0 = iszero(di)
        z0 = iszero(zi)
        if trim < πi < invtrim
            if d0 & z0
                output.mat[i] = 1.0
            elseif d0 & !z0
                output.mat[i] = 1.0 - 1.0 / πi
            elseif !d0 & z0
                output.mat[i] = 1.0 - 1.0 / (1.0 - πi)
            elseif !d0 & !z0
                output.mat[i] = 1.0
            end
        end
    end

    if checkweight(SSD)
        w = getvector(SSD, :weight)
        output.second_stage.β = _fit(output.second_stage, w .* output.mat)
        novar || (output.second_stage.V = _vcov(output, SSD.corr, w))
    else
        output.second_stage.β = _fit(output.second_stage, output.mat)
        novar || (output.second_stage.V = _vcov(output, SSD.corr))
    end

    return output
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Abadie)                    = score(second_stage(obj), obj.mat)
score(obj::Abadie, w::AbstractVector) = score(second_stage(obj), w .* obj.mat)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

jacobian(obj::Abadie)                    = jacobian(second_stage(obj), obj.mat)
jacobian(obj::Abadie, w::AbstractVector) = jacobian(second_stage(obj), w .* obj.mat)

# EXPECTED JACOBIAN OF SCORE W.R.T. FIRST-STAGE PARAMETERS × NUMBER OF OBSERVATIONS

function crossjacobian(obj::Abadie)

    d = getvector(obj.second_stage, :treatment)
    z = getvector(obj.second_stage, :instrument)
    π = fitted(obj.first_stage)
    D = zeros(π)

    @inbounds for (i, (di, zi, πi, wi)) in enumerate(zip(d, z, π, obj.mat))
        if !iszero(wi)
            d0 = iszero(di)
            z0 = iszero(zi)
            if d0 & z0
                D[i] = 0.0
            elseif d0 & !z0
                D[i] = 1.0 / πi^2
            elseif !d0 & z0
                D[i] = - 1.0 / (1.0 - πi)^2
            elseif !d0 & !z0
                D[i] = 0.0
            end
        end
    end

    g₁ = jacobexp(obj.first_stage)
    g₂ = score(obj.second_stage)

    return g₂' * (g₁ .* D)
end

function crossjacobian(obj::Abadie, w::AbstractVector)

    d = getvector(obj.second_stage, :treatment)
    z = getvector(obj.second_stage, :instrument)
    π = fitted(obj.first_stage)
    D = zeros(π)

    @inbounds for (i, (di, zi, πi, wi)) in enumerate(zip(d, z, π, obj.mat))
        if !iszero(wi)
            d0 = iszero(di)
            z0 = iszero(zi)
            if d0 & z0
                D[i] = 0.0
            elseif d0 & !z0
                D[i] = 1.0 / πi^2
            elseif !d0 & z0
                D[i] = - 1.0 / (1.0 - πi)^2
            elseif !d0 & !z0
                D[i] = 0.0
            end
        end
    end

    g₁ = jacobexp(obj.first_stage, w)
    g₂ = score(obj.second_stage, w)

    return g₂' * (g₁ .* D)
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
