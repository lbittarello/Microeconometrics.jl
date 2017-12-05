#==========================================================================================#

# ONE SAMPLE: INTERFACE

function hausman_1s(obj₁::ParOr2Stage, obj₂::ParOr2Stage)
    return hausman_1s(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_1s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, name::String)
    return hausman_1s(obj₁, obj₂, [name])
end

function hausman_1s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, names::Vector{String})

    i₁    = indexin(names, coefnames(obj₁))
    i₂    = indexin(names, coefnames(obj₂))
    w₁    = getweights(obj₁)
    w₂    = getweights(obj₂)
    corr₁ = getcorr(obj₁)
    corr₂ = getcorr(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")
    isequal(corr₁, corr₂) || throw("different correlation structures")
    isequal(getweights(obj₁), getweights(obj₂)) || throw("different weighting schemes")

    output       = ParObject()
    output.names = copy(names)

    _hausman_1s!(output, obj₁, i₁, obj₂, i₂, w₁, corr₁)

    return output
end

# ONE SAMPLE: ESTIMATION

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w::UnitWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w::FrequencyWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; scale!(w, ψ₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w::AbstractWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁) ; scale!(w, ψ₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; scale!(w, ψ₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w::UnitWeights,
        corr::ClusterOrCross
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * corr.mat * ψ₂

    adjfactor!(V₁₂, obj₁, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_1s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w::AbstractWeights,
        corr::ClusterOrCross
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁) ; scale!(w, ψ₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; scale!(w, ψ₂)
    V₁₂ = ψ₁' * corr.mat * ψ₂

    adjfactor!(V₁₂, obj₁, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

#==========================================================================================#

# TWO INDEPENDENT SAMPLES

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage)
    return hausman_2s(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, name::String)
    return hausman_2s(obj₁, obj₂, [name])
end

function hausman_2s(obj₁::ParOr2Stage, obj₂::ParOr2Stage, names::Vector{String})

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    output       = ParObject()
    output.names = copy(names)
    output.β     = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V     = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂)

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
        obj₁::ParOr2Stage, obj₂::ParOr2Stage, corr::CorrStructure, names::Vector{String}
    )

    i₁ = indexin(names, coefnames(obj₁))
    i₂ = indexin(names, coefnames(obj₂))
    w₁ = getweights(obj₁)
    w₂ = getweights(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    if !(typeof(corr) == typeof(getcorr(obj₁)) == typeof(getcorr(obj₂)))
        throw("incompatible correlation structures")
    end

    output       = ParObject()
    output.names = copy(names)

    _hausman_2s!(output, obj₁, i₁, w₁, obj₂, i₂, w₂, corr)

    return output
end

# TWO DEPENDENT SAMPLES: ESTIMATION

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        w₁::UnitWeights,
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w₂::UnitWeights,
        corr::Heteroscedastic
    )

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = touse[msng₁]
    touse₂ = touse[msng₂]

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, touse₁, i₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, touse₂, i₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁, obj₂, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        w₁::AbstractWeights,
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w₂::AbstractWeights,
        corr::Heteroscedastic
    )

    msng₁  = getmsng(obj₁)
    msng₂  = getmsng(obj₂)
    touse  = msng₁ .* msng₂
    touse₁ = touse[msng₁]
    touse₂ = touse[msng₂]

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, touse₁, i₁) ; scale!(w₁, ψ₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, touse₂, i₂) ; scale!(w₂, ψ₂)
    V₁₂ = ψ₁' * ψ₂

    adjfactor!(V₁₂, obj₁, obj₂, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        w₁::UnitWeights,
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w₂::UnitWeights,
        corr::ClusterOrCross
    )

    msng₁ = getmsng(obj₁)
    msng₂ = getmsng(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * corr.mat[msng₁, msng₂] * ψ₂

    adjfactor!(V₁₂, obj₁, obj₂, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_2s!(
        output::ParObject,
        obj₁::ParOr2Stage,
        i₁::Vector{Int},
        w₁::AbstractWeights,
        obj₂::ParOr2Stage,
        i₂::Vector{Int},
        w₂::AbstractWeights,
        corr::ClusterOrCross
    )

    msng₁ = getmsng(obj₁)
    msng₂ = getmsng(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁) ; scale!(w₁, ψ₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, :, i₂) ; scale!(w₂, ψ₂)
    V₁₂ = ψ₁' * corr.mat[msng₁, msng₂] * ψ₂

    adjfactor!(V₁₂, obj₁, obj₂, corr)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end
