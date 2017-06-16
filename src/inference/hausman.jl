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

    touse₁ = view(obj₂.sample.msng, obj₁.sample.msng)
    touse₂ = view(obj₁.sample.msng, obj₂.sample.msng)

    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = view(ψ₁, touse₁, :)' * view(ψ₂, touse₂, :)
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    V  = transpose(V₁₂[i₁, i₂])
    V .= V₁[i₁, i₁] .+ V₂[i₂, i₂] .- V₁₂[i₁, i₂] .- V

    output.β = β₁[i₁] - β₂[i₂]
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

    touse₁ = view(obj₂.sample.msng, obj₁.sample.msng)
    touse₂ = view(obj₁.sample.msng, obj₂.sample.msng)

    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = view(ψ₁, touse₁, :)' * view(ψ₂, touse₂, :)
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    V  = transpose(V₁₂[i₁, i₂])
    V .= V₁[i₁, i₁] .+ V₂[i₂, i₂] .- V₁₂[i₁, i₂] .- V

    output.β = β₁[i₁] - β₂[i₂]
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

    touse₁ = find(corr.msng .* corr₁.msng)
    touse₂ = find(corr.msng .* corr₂.msng)

    ψ₁  = influence(obj₁)
    ψ₂  = influence(obj₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    V   = transpose(V₁₂[i₁, i₂])
    V .+= V₁₂[i₁, i₂]
    _adjcluster!(V, corr)
    V .= V₁[i₁, i₁] .+ V₂[i₂, i₂] .- V

    output.β = β₁[i₁] - β₂[i₂]
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

    touse₁ = find(corr.msng .* corr₁.msng)
    touse₂ = find(corr.msng .* corr₂.msng)

    ψ₁  = influence(obj₁, w₁)
    ψ₂  = influence(obj₂, w₂)
    V₁₂ = ψ₁' * view(corr.mat, touse₁, touse₂) * ψ₂
    V₁  = vcov(obj₁)
    V₂  = vcov(obj₂)
    β₁  = coef(obj₁)
    β₂  = coef(obj₂)

    V   = transpose(V₁₂[i₁, i₂])
    V .+= V₁₂[i₁, i₂]
    _adjcluster!(V, corr)
    V .= V₁[i₁, i₁] .+ V₂[i₂, i₂] .- V

    output.β = β₁[i₁] - β₂[i₂]
    output.V = V
end
