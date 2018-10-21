#==========================================================================================#

# MICRODATA

mutable struct Microdata
    nonmissing::BitVector
    corr::CorrStructure
    weights::AbstractWeights{Int, Float64, Vector{Float64}}
    mat::Matrix{Float64}
    names::Vector{String}
    map::Dict{Symbol, Vector{Int}}
    terms::Vector
end

# CONSTRUCTOR

function Microdata(
        df::DataFrame, model::Dict{Symbol, String};
        contrasts::Dict = Dict(),
        subset::AbstractVector{Bool} = trues(size(df, 1)),
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = uweights(fill(1.0, size(df, 1)))
    )

    input    = reduce((x, y) -> x * " + " * model[y], keys(model), init = "")
    formula  = @eval @formula $nothing ~ $(Meta.parse(input))
    terms    = StatsModels.Terms(formula) ; terms.intercept = true
    eterms   = terms.eterms
    msng     = BitVector(completecases(df[terms.eterms]))
    msng    .= msng .* .!iszero.(weights) .* BitVector(subset)
    new_corr = adjmissing!(msng, vcov)
    new_wts  = parse_weights(weights, msng)
    new_df   = DataFrame(map(x -> disallowmissing(df[x][msng]), eterms), Symbol.(eterms))
    frame    = ModelFrame(new_df, terms, msng, StatsModels.evalcontrasts(new_df, contrasts))
    names    = StatsModels.coefnames(frame)
    mat      = ModelMatrix(frame)
    eterms   = eterms[mat.assign[2:end]]
    new_map  = Dict{Symbol, Vector{Int}}()

    for (i, j) in model
        new_map[i] = assign_columns(j, eterms)
    end

    return Microdata(msng, new_corr, new_wts, mat.m, names, new_map, eterms)
end

# REASSIGN VARIABLE SETS

function Microdata(MD::Microdata, model::Dict{Symbol, String})

    new_map = copy(MD.map)

    for (i, j) in model
        (j == "") ? pop!(new_map, i) : (new_map[i] = assign_columns(j, MD.terms))
    end

    Microdata(MD.nonmissing, MD.corr, MD.weights, MD.mat, MD.names, new_map, MD.terms)
end
