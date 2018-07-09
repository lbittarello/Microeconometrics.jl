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
        df::DataFrame;
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights(fill(1.0, size(df, 1))),
        subset::AbstractVector{Bool} = trues(size(df, 1)),
        kwargs...
    )

    input    = reduce((x, y) -> x * " + " * y[2], "", kwargs)
    formula  = @eval @formula $nothing ~ $(parse(input))
    terms    = Terms(formula)
    eterms   = terms.eterms
    msng     = BitVector(completecases(df[:, terms.eterms]))
    msng    .= msng .* (.!iszero.(weights)) .* BitVector(subset)
    new_corr = adjmsng!(msng, vcov)
    new_wts  = parse_weights(weights, msng)
    new_df   = DataFrame(map(x -> disallowmissing(df[x][msng]), eterms), Symbol.(eterms))
    frame    = ModelFrame(new_df, terms, msng)
    names    = StatsModels.coefnames(frame)
    mat      = ModelMatrix(frame)

    new_map = Dict{Symbol, Vector{Int}}()

    for (i, j) in kwargs
        new_map[i] = assign_columns(j, terms, mat.assign)
    end

    return Microdata(msng, mat, new_corr, new_wts, names, new_map, terms)
end

# REASSIGN VARIABLE SETS

function Microdata(MD::Microdata; kwargs...)

    new_map = copy(MD.map)

    for (i, j) in kwargs
        if j == ""
            pop!(new_map, i)
        else
            new_map[i] = assign_columns(j, MD.terms, MD.mat.assign)
        end
    end

    Microdata(MD.msng, MD.mat, MD.corr, MD.weights, MD.names, new_map, MD.terms)
end
