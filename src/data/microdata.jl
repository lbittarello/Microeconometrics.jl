mutable struct Microdata
    nonmissing::BitVector
    model::Dict{Symbol, SM.TermOrTerms}
    mapping::Dict{Symbol, Vector{Int}}
    matrix::AbstractMatrix{Float64}
    weights::AbstractWeights{Int, Float64, Vector{Float64}}
    corr::CorrStructure
end

# CONSTRUCTORS

function Microdata(
        df::AbstractDataFrame,
        model::Dict{Symbol, SM.TermOrTerms};
        hints::Dict{Symbol} = Dict{Symbol, Any}(),
        subset::AbstractVector{Bool} = [true],
        corr::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights{Float64}(size(df, 1))
    )

    columns = termvars(model)

    nonmissing = .!iszero.(weights) .& corr.nonmissing .& subset
    nonmissing[nonmissing] .&= completecases(view(df, nonmissing, columns))

    corr    = update_corr(corr, nonmissing)
    weights = update_weights(weights, nonmissing)
    data    = view(df, nonmissing, columns)
    model   = apply_schema(model, schema(model, data, hints))
    mapping = asgn(model)
    matrix  = modelcols(model, data)

    return Microdata(nonmissing, model, mapping, matrix, weights, corr)
end

function Microdata(MD::Microdata)
    Microdata(
        copy(MD.nonmissing),
        copy(MD.model),
        copy(MD.mapping),
        MD.matrix,
        MD.weights,
        MD.corr
    )
end
