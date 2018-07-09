#==========================================================================================#

# MICRODATA

mutable struct Microdata
    msng::BitVector
    mat::ModelMatrix{Matrix{Float64}}
    corr::CorrStructure
    weights::AbstractWeights{Int, Float64, Vector{Float64}}
    names::Vector{String}
    map::Dict{Symbol, Vector{Int}}
    terms::Terms
end

# CONSTRUCTOR

function Microdata(
        df::DataFrame, model::Dict{Symbol, String};
        contrasts::Dict = Dict(),
        subset::AbstractVector{Bool} = trues(size(df, 1)),
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights(fill(1.0, size(df, 1)))
    )

    input    = reduce((x, y) -> x * " + " * model[y], "", keys(model))
    formula  = @eval @formula $nothing ~ $(parse(input))
    terms    = Terms(formula)
    eterms   = terms.eterms
    msng     = BitVector(completecases(df[:, terms.eterms]))
    msng    .= msng .* (.!iszero.(weights)) .* BitVector(subset)
    new_corr = adjmsng!(msng, vcov)
    new_wts  = parse_weights(weights, msng)
    new_df   = DataFrame(map(x -> disallowmissing(df[x][msng]), eterms), Symbol.(eterms))
    frame    = ModelFrame(new_df, terms, msng, evalcontrasts(df, contrasts))
    names    = StatsModels.coefnames(frame)
    mat      = ModelMatrix(frame)
    new_map  = Dict{Symbol, Vector{Int}}()

    for (i, j) in model
        new_map[i] = assign_columns(j, terms, mat.assign)
    end

    return Microdata(msng, mat, new_corr, new_wts, names, new_map, terms)
end

# REASSIGN VARIABLE SETS

function Microdata(MD::Microdata, model::Dict{Symbol, String})

    new_map = copy(MD.map)

    for (i, j) in model
        if j == ""
            pop!(new_map, i)
        else
            new_map[i] = assign_columns(j, MD.terms, MD.mat.assign)
        end
    end

    Microdata(MD.msng, MD.mat, MD.corr, MD.weights, MD.names, new_map, MD.terms)
end
