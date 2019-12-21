#==========================================================================================#

# ONE SAMPLE: INTERFACE

function hausman_test(obj₁::ParModel, obj₂::ParModel)
    hausman_test(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function hausman_test(obj₁::ParModel, obj₂::ParModel, name::String)
    hausman_test(obj₁, obj₂, [name])
end

function hausman_test(obj₁::ParModel, obj₂::ParModel, names::Vector{String})

    i₁    = findall((in)(names), coefnames(obj₁))
    i₂    = findall((in)(names), coefnames(obj₂))
    w₁    = getweights(obj₁)
    w₂    = getweights(obj₂)
    corr₁ = getcorr(obj₁)
    corr₂ = getcorr(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")
    isequal(corr₁, corr₂) || throw("different correlation structures")
    isequal(w₁, w₂) || throw("different weighting schemes")

    output       = ParEstimate()
    output.names = copy(names)

    _hausman_test!(output, obj₁, i₁, obj₂, i₂, w₁, corr₁)

    return output
end

# ONE SAMPLE: ESTIMATION

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::UnitWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::FrequencyWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; lmul!(Diagonal(w), ψ₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::AbstractWeights,
        corr::Heteroscedastic
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁) ; lmul!(Diagonal(w), ψ₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; lmul!(Diagonal(w), ψ₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::UnitWeights,
        corr::Clustered
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = corr.mat * view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = corr.mat * view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::AbstractWeights,
        corr::Clustered
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = corr.mat * view(Ψ₁, :, i₁) ; lmul!(Diagonal(w), ψ₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = corr.mat * view(Ψ₂, :, i₂) ; lmul!(Diagonal(w), ψ₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::UnitWeights,
        corr::CrossCorrelated
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * corr.mat * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _hausman_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        obj₂::ParModel,
        i₂::Vector{Int},
        w::AbstractWeights,
        corr::CrossCorrelated
    )

    Ψ₁  = influence(obj₁, w)
    ψ₁  = view(Ψ₁, :, i₁) ; lmul!(Diagonal(w), ψ₁)
    Ψ₂  = influence(obj₂, w)
    ψ₂  = view(Ψ₂, :, i₂) ; lmul!(Diagonal(w), ψ₂)
    V₁₂ = ψ₁' * corr.mat * ψ₂

    corr.corrected && lmul!(varcorrection(corr, w), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

#==========================================================================================#

# TWO INDEPENDENT SAMPLES

function chow_test(obj₁::ParModel, obj₂::ParModel)
    chow_test(obj₁, obj₂, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function chow_test(obj₁::ParModel, obj₂::ParModel, name::String)
    chow_test(obj₁, obj₂, [name])
end

function chow_test(obj₁::ParModel, obj₂::ParModel, names::Vector{String})

    i₁ = findall((in)(names), coefnames(obj₁))
    i₂ = findall((in)(names), coefnames(obj₂))

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    output       = ParEstimate()
    output.names = copy(names)
    output.β     = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V     = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂)

    return output
end

#==========================================================================================#

# TWO DEPENDENT SAMPLES: INTERFACE

function chow_test(obj₁::ParModel, obj₂::ParModel, corr::CorrStructure)
    chow_test(obj₁, obj₂, corr, intersect(coefnames(obj₁), coefnames(obj₂)))
end

function chow_test(obj₁::ParModel, obj₂::ParModel, corr::CorrStructure, name::String)
    chow_test(obj₁, obj₂, corr, [name])
end

function chow_test(
        obj₁::ParModel,
        obj₂::ParModel,
        corr::CorrStructure,
        names::Vector{String}
    )

    i₁ = findall((in)(names), coefnames(obj₁))
    i₂ = findall((in)(names), coefnames(obj₂))
    w₁ = getweights(obj₁)
    w₂ = getweights(obj₂)

    (iszero(i₁) | iszero(i₂)) && throw("missing coefficients in at least one model")

    if !(typeof(corr) == typeof(getcorr(obj₁)) == typeof(getcorr(obj₂)))
        throw("incompatible correlation structures")
    end

    output       = ParEstimate()
    output.names = copy(names)

    _chow_test!(output, obj₁, i₁, w₁, obj₂, i₂, w₂, corr)

    return output
end

# TWO DEPENDENT SAMPLES: ESTIMATION

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::UnitWeights,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::UnitWeights,
        corr::Heteroscedastic
    )

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)
    touse       = nonmissing₁ .* nonmissing₂
    touse₁      = touse[nonmissing₁]
    touse₂      = touse[nonmissing₂]

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, touse₁, i₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, touse₂, i₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::W,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::W,
        corr::Heteroscedastic
    ) where W <: AbstractWeights

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)
    touse       = nonmissing₁ .* nonmissing₂
    touse₁      = touse[nonmissing₁]
    touse₂      = touse[nonmissing₂]

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, touse₁, i₁) ; lmul!(Diagonal(w₁), ψ₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, touse₂, i₂) ; lmul!(Diagonal(w₂), ψ₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::UnitWeights,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::UnitWeights,
        corr::Clustered
    )

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(corr.mat, :, nonmissing₁) * view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(corr.mat, :, nonmissing₂) * view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::W,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::W,
        corr::Clustered
    ) where W <: AbstractWeights

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(corr.mat, :, nonmissing₁) * view(Ψ₁, :, i₁) ; lmul!(Diagonal(w₁), ψ₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(corr.mat, :, nonmissing₂) * view(Ψ₂, :, i₂) ; lmul!(Diagonal(w₂), ψ₂)
    V₁₂ = ψ₁' * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::UnitWeights,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::UnitWeights,
        corr::CorrStructure
    )

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, :, i₂)
    V₁₂ = ψ₁' * view(corr.mat, nonmissing₁, nonmissing₂) * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end

function _chow_test!(
        output::ParEstimate,
        obj₁::ParModel,
        i₁::Vector{Int},
        w₁::W,
        obj₂::ParModel,
        i₂::Vector{Int},
        w₂::W,
        corr::CorrStructure
    ) where W <: AbstractWeights

    nonmissing₁ = getnonmissing(obj₁)
    nonmissing₂ = getnonmissing(obj₂)

    Ψ₁  = influence(obj₁, w₁)
    ψ₁  = view(Ψ₁, :, i₁) ; lmul!(Diagonal(w₁), ψ₁)
    Ψ₂  = influence(obj₂, w₂)
    ψ₂  = view(Ψ₂, :, i₂) ; lmul!(Diagonal(w₂), ψ₂)
    V₁₂ = ψ₁' * view(corr.mat, nonmissing₁, nonmissing₂) * ψ₂

    corr.corrected && lmul!(varcorrection(corr, UnitWeights(sum(touse))), V₁₂)

    output.β = view(coef(obj₁), i₁) - view(coef(obj₂), i₂)
    output.V = view(vcov(obj₁), i₁, i₁) + view(vcov(obj₂), i₂, i₂) - (V₁₂' + V₁₂)
end
