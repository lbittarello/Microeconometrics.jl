#==========================================================================================#

# INTERFACE

function hausman_test(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage)
    names = intersect(coefnames(obj₁), coefnames(obj₂))
    return hausman_test(obj₁, obj₂, names, getcorr(obj₁))
end

function hausman_test(obj₁::ParOrTwoStage, obj₂::ParOrTwoStage, corr::CorrStructure)
    names = intersect(coefnames(obj₁), coefnames(obj₂))
    return hausman_test(obj₁, obj₂, names, corr)
end

function hausman_test(
        obj₁::ParOrTwoStage,
        obj₂::ParOrTwoStage,
        names::String,
        corr::CorrStructure = getcorr(obj₁)
    )

    return hausman_test(obj₁, obj₂, [names], corr)
end

function hausman_test(
        obj₁::ParOrTwoStage,
        obj₂::ParOrTwoStage,
        names::Vector{String},
        corr::CorrStructure = getcorr(obj₁)
    )

    if typeof(getcorr(obj₁)) != typeof(corr)
        throw("different correlation structures")
    end
    if typeof(getcorr(obj₂)) != typeof(corr)
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
        w₂ = ones(Float64, nobs(obj₂))
        _hausman!(output, obj₁, i₁, obj₂, i₂, corr, w₁, w₂)
    elseif !W₁ & W₂
        w₁ = ones(Float64, nobs(obj₁))
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

    touse₁ = view(getcorr(obj₂).msng, getcorr(obj₁).msng)
    touse₂ = view(getcorr(obj₁).msng, getcorr(obj₂).msng)

    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = view(ψ₁, touse₁, :)' * view(ψ₂, touse₂, :)
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    output.β = β₁[i₁] - β₂[i₂]
    output.V = V₁[i₁, i₁] + V₂[i₂, i₂] - (V₁₂[i₁, i₂] + V₁₂[i₁, i₂]')
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

    touse₁ = find(corr.msng .* getcorr(obj₁).msng)
    touse₂ = find(corr.msng .* getcorr(obj₂).msng)

    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = view(ψ₁, touse₁, :)' * view(ψ₂, touse₂, :)
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    output.β = β₁[i₁] - β₂[i₂]
    output.V = V₁[i₁, i₁] + V₂[i₂, i₂] - (V₁₂[i₁, i₂] + V₁₂[i₁, i₂]')
end

function _hausman!(
        output::ParObject,
        obj₁::ParOrTwoStage,
        i₁::Vector{Int},
        obj₂::ParOrTwoStage,
        i₂::Vector{Int},
        corr::ClusterOrCross
    )

    if any(corr.msng[getcorr(obj₁).msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    elseif any(corr.msng[getcorr(obj₂).msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    end

    touse₁ = find(corr.msng .* getcorr(obj₁).msng)
    touse₂ = find(corr.msng .* getcorr(obj₂).msng)

    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    output.β = β₁[i₁] - β₂[i₂]
    output.V = V₁[i₁, i₁] + V₂[i₂, i₂] - _adjcluster(corr) * (V₁₂[i₁, i₂] + V₁₂[i₁, i₂]')
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

    if any(corr.msng[getcorr(obj₁).msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    elseif any(corr.msng[getcorr(obj₂).msng] .== false)
        throw("joint correlation structure does not overlap with both estimation samples")
    end

    touse₁ = find(corr.msng .* getcorr(obj₁).msng)
    touse₂ = find(corr.msng .* getcorr(obj₂).msng)

    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    output.β = β₁[i₁] - β₂[i₂]
    output.V = V₁[i₁, i₁] + V₂[i₂, i₂] - _adjcluster(corr) * (V₁₂[i₁, i₂] + V₁₂[i₁, i₂]')
end
