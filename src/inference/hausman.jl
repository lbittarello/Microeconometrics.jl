#==========================================================================================#

# ONE SAMPLE: INTERFACE

function hausman_1s(obj₁::ParOr2Stage, obj₂::ParOr2Stage)
    return hausman_1s(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_1s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, name::String)
    return hausman_1s(obj₁, obj₂, [name])
end

function hausman_1s(
        obj₁::ParOr2Stage{T},
        obj₂::ParOr2Stage{T},
        names::Vector{String}
    ) where {T <: CorrStructure}

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))
    W₁ = checkweight(obj₁)
    W₂ = checkweight(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    output       = ParObject()
    output.names = copy(names)

    if W₁ | W₂
        w₁ = (W₁ ? getvector(obj₁, :weight) : fill(1.0, nobs(obj₁)))
        w₂ = (W₂ ? getvector(obj₂, :weight) : fill(1.0, nobs(obj₂)))
        _hausman_1s!(output, obj₁, i₁, obj₂, i₂, w₁, w₂)
    else
        _hausman_1s!(output, obj₁, i₁, obj₂, i₂)
    end

    return output
end

# ONE SAMPLE: ESTIMATION

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int}
    ) where {T <: Heteroscedastic}

    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        w₁::AbstractVector,
        w₂::AbstractVector
    ) where {T <: Heteroscedastic}

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w₁)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int}
    ) where {T <: ClusterOrCross}

    Ω   = getcorr(obj₁).mat
    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * Ω * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        w₁::AbstractVector,
        w₂::AbstractVector
    ) where {T <: ClusterOrCross}

    Ω   = getcorr(obj₁).mat
    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w₁)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * Ω * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

#==========================================================================================#

# TWO INDEPENDENT SAMPLES

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage)
    return hausman(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, name::String)
    return hausman(obj₁, obj₂, [name])
end

function hausman_2s(
        obj₁::ParOr2Stage{T},
        obj₂::ParOr2Stage{T},
        names::Vector{String}
    ) where {T <: CorrStructure}

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    output       = ParObject()
    output.names = copy(names)
    output.β     = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.names = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂)

    return output
end

#==========================================================================================#

# TWO DEPENDENT SAMPLES: INTERFACE

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, corr::CorrStructure)
    return hausman_2s(obj₁, obj₂, corr, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, corr::CorrStructure, name::String)
    return hausman_2s(obj₁, obj₂, corr, [name])
end

function hausman_2s(
        obj₁::ParOr2Stage{T},
        obj₂::ParOr2Stage{T},
        corr::T,
        names::Vector{String}
    ) where {T <: CorrStructure}

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))
    W₁ = checkweight(obj₁)
    W₂ = checkweight(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    output       = ParObject()
    output.names = copy(names)

    if W₁ | W₂
        w₁ = (W₁ ? getvector(obj₁, :weight) : fill(1.0, nobs(obj₁)))
        w₂ = (W₂ ? getvector(obj₂, :weight) : fill(1.0, nobs(obj₂)))
        _hausman_2s!(output, obj₁, i₁, obj₂, i₂, w₁, w₂)
    else
        _hausman_2s!(output, obj₁, i₁, obj₂, i₂, corr)
    end

    return output
end

# TWO DEPENDENT SAMPLES: ESTIMATION

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        corr::T
    ) where {T <: Heteroscedastic}

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = msng₁[touse]
    touse₂ = msng₂[touse]

    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        corr::T,
        w₁::AbstractVector,
        w₂::AbstractVector
    ) where {T <: Heteroscedastic}

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = msng₁[touse]
    touse₂ = msng₂[touse]

    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        corr::T
    ) where {T <: ClusterOrCross}

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse₁ = corr.msng .* msng₁
    touse₁ = touse₁[corr.msng]
    touse₂ = corr.msng .* msng₂
    touse₂ = touse₂[corr.msng]

    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * corr.mat[touse₁, touse₂] * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage{T},
        i₁::Vector{Int},
        obj₂::ParOr2Stage{T},
        i₂::Vector{Int},
        corr::T,
        w₁::AbstractVector,
        w₂::AbstractVector
    ) where {T <: ClusterOrCross}

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse₁ = corr.msng .* msng₁
    touse₁ = touse₁[corr.msng]
    touse₂ = corr.msng .* msng₂
    touse₂ = touse₂[corr.msng]

    Ψ₁  = influence(obj₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * corr.mat[touse₁, touse₂] * ψ₂

    adjfactor!(V₁₂, obj₁)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end
