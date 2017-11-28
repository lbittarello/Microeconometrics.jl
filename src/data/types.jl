#==========================================================================================#

# TYPES

mutable struct Microdata{T <: CorrStructure}
    msng::BitVector
    mat::Matrix{Float64}
    names::Vector{String}
    map::Dict{Symbol, Vector{Int}}
    corr::T
    terms::Terms
    assign::Vector{Int}
end

#==========================================================================================#

# COPY

function copy(MD::Microdata{T}) where {T}

    msng   = copy(MD.msng)
    mat    = copy(MD.mat)
    names  = copy(MD.names)
    map    = copy(MD.map)
    corr   = copy(MD.corr)
    terms  = copy(MD.terms)
    assign = copy(MD.assign)

    return Microdata{T}(msng, mat, names, map, newcorr, terms, assign)
end

#==========================================================================================#

# CONSTRUCTOR

function Microdata(df::DataFrame; kwargs...)
    return Microdata(df, Heteroscedastic(), trues(size(df, 1)); kwargs...)
end

function Microdata(df::DataFrame, corr::CorrStructure; kwargs...)
    return Microdata(df, corr, trues(size(df, 1)); kwargs...)
end

function Microdata(df::DataFrame, subset::AbstractVector{Bool}; kwargs...)
    return Microdata(df, Heteroscedastic(), subset; kwargs...)
end

function Microdata(
        df::DataFrame,
        corr::T,
        subset::AbstractVector{Bool};
        kwargs...
    ) where {T <: CorrStructure}

    input   = reduce((x, y) -> x * " + " * y[2], "", kwargs)
    formula = StatsModels.Formula(nothing, parse(input))
    terms   = StatsModels.Terms(formula)
    msng    = BitVector(completecases(df[:, terms.eterms]))
    msng   .= msng .* BitVector(subset)
    newcorr = adjmsng!(msng, corr)
    frame   = ModelFrame(terms, df[msng, :])
    names   = coefnames(frame)
    mat     = ModelMatrix(frame)

    map = Dict{Symbol, Vector{Int}}()

    for (i, j) in kwargs
        map[i] = assign_columns(j, terms, mat.assign)
    end

    return Microdata{T}(msng, mat.m, names, map, newcorr, terms, mat.assign)
end

#==========================================================================================#

# REASSIGN VARIABLE SETS

function Microdata(MD::Microdata{T}; kwargs...) where {T}

    map = copy(MD.map)

    for (i, j) in kwargs
        (j == "") ? pop!(map, i) : (map[i] = assign_columns(j, MD.terms, MD.assign))
    end

    Microdata{T}(MD.msng, MD.mat, MD.names, map, MD.corr, MD.terms, MD.assign)
end
