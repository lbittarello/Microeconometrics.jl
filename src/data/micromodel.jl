macro component!(x)
    SM.terms!(SM.sort_terms!(SM.parse!(x)))
end

function component(x::Expr)

    (x.args[1] != :(=>)) && throw(ArgumentError("invalid expression encountered: $x"))

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

function schema(y::Dict{Symbol, SM.TermOrTerms}, df::AbstractDataFrame, hints::Dict)
    schema(filter(SM.needs_schema, terms(y)), df, hints)
end

function apply_schema(y::Dict{Symbol, SM.TermOrTerms}, schema::SM.Schema)
    Dict{Symbol, SM.TermOrTerms}(i => apply_schema(v, schema) for (i, v) in y)
end

has_schema(y::Dict{Symbol, SM.TermOrTerms}) = all(has_schema, values(y))

function modelcols(y::Dict{Symbol, SM.TermOrTerms}, df::AbstractDataFrame)
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
