#==========================================================================================#

# AUXILIARY FUNCTIONS

macro component!(x)
    SM.terms!(SM.sort_terms!(SM.parse!(x)))
end

function component(x::Expr)

    (x.args[1] == :(=>)) || throw(ArgumentError("invalid expression encountered: $x"))

    y = @eval @component! $(x.args[3])

    if (typeof(y) <: Tuple) && (typeof(y[1]) <: SM.ConstantTerm)
        y = (y[2:end]..., y[1])
    end

    return y
end

macro micromodel(x...)
    Dict{Symbol, SM.TermOrTerms}(xi.args[2] => component(xi) for xi in x)
end

terms(y::Dict{Symbol, SM.TermOrTerms})    = mapreduce(terms, union, values(y))
termvars(y::Dict{Symbol, SM.TermOrTerms}) = mapreduce(termvars, union, values(y))
eterms(y::Dict{Symbol, SM.TermOrTerms})   = Tuple(mapreduce(SM.vectorize, union, values(y)))

function schema(y::Dict{Symbol, SM.TermOrTerms}, df::DataFrame, hints::Dict)
    schema(filter(SM.needs_schema, terms(y)), df, hints)
end

function apply_schema(y::Dict{Symbol, SM.TermOrTerms}, schema::SM.Schema)
    Dict{Symbol, SM.TermOrTerms}(i => apply_schema(v, schema) for (i, v) in y)
end

has_schema(y::Dict{Symbol, SM.TermOrTerms}) = all(has_schema, values(y))

function modelcols(y::Dict{Symbol, SM.TermOrTerms}, df::DataFrame)
    modelcols(SM.collect_matrix_terms(eterms(y)), df)
end

function asgn(y::Dict{Symbol, SM.TermOrTerms})

    Y = eterms(y)
    Y = Y[SM.asgn(Y)]

    return Dict{Symbol, Vector{Int}}(i => findall((in)(SM.vectorize(v)), Y) for (i, v) in y)
end

function asgn(y::Dict{Symbol, SM.TermOrTerms}, t::SM.TermOrTerms)

    Y = eterms(y)
    Y = Y[SM.asgn(Y)]

    return [findfirst(isequal(ti), Y) for ti in t]
end

#==========================================================================================#

# TYPE

mutable struct Microdata
    nonmissing::BitVector
    model::Dict{Symbol, SM.TermOrTerms}
    mapping::Dict{Symbol, Vector{Int}}
    matrix::Matrix{Float64}
    weights::AbstractWeights{Int, Float64, Vector{Float64}}
    corr::CorrStructure    
end

# CONSTRUCTORS

function Microdata(
        df::DataFrame, model::Dict{Symbol, SM.TermOrTerms};
        hints::Dict{Symbol} = Dict{Symbol, Any}(),
        subset::AbstractVector{Bool} = trues(size(df, 1)),
        vcov::CorrStructure = Heteroscedastic(),
        weights::AbstractWeights = UnitWeights(Float64, size(df, 1))
    )

    nonmissing  = completecases(df, termvars(model))
    nonmissing .= nonmissing .* .!iszero.(weights) .* subset
    corr        = parse_corr!(nonmissing, vcov)
    weights     = parse_weights(nonmissing, weights)
    columns     = termvars(model)
    data        = DataFrame(map(x -> disallowmissing(df[nonmissing, x]), columns), columns)
    model       = apply_schema(model, schema(model, data, hints))
    mapping     = asgn(model)
    matrix      = modelcols(model, data)
    mapping     = asgn(model)

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

#==========================================================================================#

# COPY

copy(m::Microdata) = Microdata([copy(getfield(m, k))     for k in fieldnames(m)]...)
