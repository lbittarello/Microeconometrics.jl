#==========================================================================================#

# TYPES

mutable struct Microdata
    msng::BitVector
    mat::Matrix{Float64}
    names::Vector{String}
    map::Dict{Symbol, Vector{Int}}
    corr::CorrStructure
    terms::Terms
    assign::Vector{Int}
end

#==========================================================================================#

# COPY

function copy(MD::Microdata)

    msng   = copy(MD.msng)
    mat    = copy(MD.mat)
    names  = copy(MD.names)
    map    = copy(MD.map)
    corr   = copy(MD.corr)
    terms  = copy(MD.terms)
    assign = copy(MD.assign)

    return Microdata(msng, mat, names, map, newcorr, terms, assign)
end

#==========================================================================================#

# CONSTRUCTOR

function Microdata(
        df::DataFrame,
        subset::AbstractVector{Bool} = trues(size(df, 1));
        corr::CorrStructure = Heteroscedastic(),
        makecopy::Bool = true,
        checkrank::Bool = true,
        normalize::Bool = true,
        kwargs...
    )

    input = ""

    for (i, j) in kwargs
        input = input * " + " * j
    end

    formula = DataFrames.Formula(nothing, parse(input))
    terms   = DataFrames.Terms(formula)
    msng    = convert(BitVector, completecases(df[:, terms.eterms]))
    msng   .= msng .* BitVector(subset)
    newcorr = adjmsng!(msng, corr, makecopy)
    frame   = ModelFrame(terms, df[msng, :])
    names   = coefnames(frame)
    mat     = ModelMatrix(frame)

    if checkrank
        (rank(mat.m) == size(mat.m, 2)) || throw("model matrix does not have full rank")
    end

    map = Dict{Symbol, Vector{Int}}()

    for (i, j) in kwargs
        map[i] = assign_columns(j, terms, mat.assign)
    end

    if haskey(map, :weight) & normalize
        w     = view(mat.m, :, map[:weight])
        nrm_w = length(w) / sum(w)
        w    .= w .* nrm_w
    end

    return Microdata(msng, mat.m, names, map, newcorr, terms, mat.assign)
end

#==========================================================================================#

# REASSIGN VARIABLE SETS

function Microdata(MD::Microdata; makecopy::Bool = false, kwargs...)

    map = copy(MD.map)

    for (i, j) in kwargs
        (j == "") ? pop!(map, i) : (map[i] = assign_columns(j, MD.terms, MD.assign))
    end

    newmd = Microdata(MD.msng, MD.mat, MD.names, map, MD.corr, MD.terms, MD.assign)

    makecopy ? (return copy(newmd)) : (return newmd)
end
