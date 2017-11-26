#==========================================================================================#

# INTERFACE

function hausman(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage)
    return hausman(obj₁, obj₂, getcorr(obj₁), intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage, corr::CorrStructure)
    return hausman(obj₁, obj₂, corr, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage, names::String)
    return hausman(obj₁, obj₂, getcorr(obj₁), [names])
end

function hausman(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage, names::Vector{String})
    return hausman(obj₁, obj₂, getcorr(obj₁), names)
end

function hausman(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage, corr::CorrStructure, names::String)
    return hausman(obj₁, obj₂, corr, [names])
end

function hausman(
        obj₁::ParOrTwoStage,
        obj₂::ParOrTwoStage,
        corr::CorrStructure,
        names::Vector{String}
    )

    if (typeof(getcorr(obj₁)) != typeof(corr)) | (typeof(getcorr(obj₂)) != typeof(corr))
        throw("different correlation structures")
    end

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    W₁ = checkweight(obj₁)
    W₂ = checkweight(obj₂)

    output       = ParObject()
    output.names = copy(names)

    if W₁ & W₂
        w₁ = getvector(obj₁, :weight)
        w₂ = getvector(obj₂, :weight)
        _hausman!(output, obj₁, i₁, obj₂, i₂, corr, w₁, w₂)
    elseif W₁ & !W₂
        w₁ = getvector(obj₁, :weight)
        w₂ = fill(1.0, nobs(obj₂))
        _hausman!(output, obj₁, i₁, obj₂, i₂, corr, w₁, w₂)
    elseif !W₁ & W₂
        w₁ = fill(1.0, nobs(obj₁))
        w₂ = getvector(obj₂, :weight)
        _hausman!(output, obj₁, i₁, obj₂, i₂, corr, w₁, w₂)
    else
        _hausman!(output, obj₁, i₁, obj₂, i₂, corr)
    end

    return output
end

#==========================================================================================#

# TWO-SAMPLE HAUSMAN TEST

function _hausman!(
        output::ParObject,
        obj₁::ParOrTwoStage,
        i₁::Vector{Int},
        obj₂::ParOrTwoStage,
        i₂::Vector{Int},
        corr::Heteroscedastic
    )

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = msng₁[touse]
    touse₂ = msng₂[touse]

    β₁  = coef(obj₁)
    β₂  = coef(obj₂)
    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = view(ψ₁, touse₁, i₁)' * view(ψ₂, touse₂, i₂)
    V₁  = view(vcov(obj₁), i₁, i₁)
    V₂  = view(vcov(obj₂), i₂, i₂)

    V  = V₁₂' + V₁₂
    V .= V₁ .+ V₂ .- V

    output.β = view(β₁, i₁) - view(β₂, i₂)
    output.V = V
end

function _hausman!(
        output::ParObject,
        obj₁::ParOrTwoStage,
        i₁::Vector{Int},
        obj₂::ParOrTwoStage,
        i₂::Vector{Int},
        corr::Heteroscedastic,
        w₁::AbstractVector,
        w₂::AbstractVector
    )

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = msng₁[touse]
    touse₂ = msng₂[touse]

    β₁  = coef(obj₁)
    β₂  = coef(obj₂)
    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = view(ψ₁, touse₁, i₁)' * view(ψ₂, touse₂, i₂)
    V₁  = view(vcov(obj₁), i₁, i₁)
    V₂  = view(vcov(obj₂), i₂, i₂)

    V  = V₁₂' + V₁₂
    V .= V₁ .+ V₂ .- V

    output.β = view(β₁, i₁) - view(β₂, i₂)
    output.V = V
end

function _hausman!(
        output::ParObject,
        obj₁::ParOrTwoStage,
        i₁::Vector{Int},
        obj₂::ParOrTwoStage,
        i₂::Vector{Int},
        corr::ClusterOrCross
    )

    corr₁ = getcorr(obj₁)
    corr₂ = getcorr(obj₂)

    if any(corr.msng[corr₁.msng] .== false) | any(corr.msng[corr₂.msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    end

    touse₁ = corr.msng .* corr₁.msng
    touse₁ = touse₁[corr.msng]
    touse₂ = corr.msng .* corr₂.msng
    touse₂ = touse₂[corr.msng]

    β₁  = coef(obj₁)
    β₂  = coef(obj₂)
    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁₂ = view(V₁₂, i₁, i₂)
    V₁  = view(vcov(obj₁), i₁, i₁)
    V₂  = view(vcov(obj₂), i₂, i₂)

    V  = adjfactor!(V₁₂' + V₁₂, corr)
    V .= V₁ .+ V₂ .- V

    output.β = view(β₁, i₁) - view(β₂, i₂)
    output.V = V
end

function _hausman!(
        output::ParObject,
        obj₁::ParOrTwoStage,
        i₁::Vector{Int},
        obj₂::ParOrTwoStage,
        i₂::Vector{Int},
        corr::ClusterOrCross,
        w₁::AbstractVector,
        w₂::AbstractVector
    )

    corr₁ = getcorr(obj₁)
    corr₂ = getcorr(obj₂)

    if any(corr.msng[corr₁.msng] .== false) | any(corr.msng[corr₂.msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    end

    touse₁ = corr.msng .* corr₁.msng
    touse₁ = touse₁[corr.msng]
    touse₂ = corr.msng .* corr₂.msng
    touse₂ = touse₂[corr.msng]

    β₁  = coef(obj₁)
    β₂  = coef(obj₂)
    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁₂ = view(V₁₂, i₁, i₂)
    V₁  = view(vcov(obj₁), i₁, i₁)
    V₂  = view(vcov(obj₂), i₂, i₂)

    V  = adjfactor!(V₁₂' + V₁₂, corr)
    V .= V₁ .+ V₂ .- V

    output.β = view(β₁, i₁) - view(β₂, i₂)
    output.V = V
end
