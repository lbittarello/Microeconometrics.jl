#==========================================================================================#

# PARSE WEIGHTS

function parse_weights(nonmissing::BitVector, w::AnalyticWeights)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    AnalyticWeights(v, n)
end

function parse_weights(nonmissing::BitVector, w::FrequencyWeights)
    v = w[nonmissing]
    s = Int(sum(v))
    FrequencyWeights(float(v), s)
end

function parse_weights(nonmissing::BitVector, w::ProbabilityWeights)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    ProbabilityWeights(v, n)
end

parse_weights(nonmissing::BitVector, w::UnitWeights) = UnitWeights(Float64, sum(nonmissing))

function parse_weights(nonmissing::BitVector, w::Weights)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    Weights(v, n)
end

#==========================================================================================#

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

parse_corr!(nonmissing::BitVector, corr::Homoscedastic)   = copy(corr)
parse_corr!(nonmissing::BitVector, corr::Heteroscedastic) = copy(corr)

function parse_corr!(nonmissing::BitVector, corr::Clustered)

    (nonmissing == corr.nonmissing) && (return copy(corr))

    nonmissing .*= corr.nonmissing
    touse        = nonmissing[corr.nonmissing]
    new_ic       = corr.ic[touse]
    iter         = sort(unique(corr.ic))
    new_iter     = sort(unique(new_ic))
    new_nc       = length(new_iter)
    new_mat      = corr.mat[findall((in)(new_iter), iter), touse]

    return Clustered(corr.adj, nonmissing, new_mat, new_ic, new_nc)
end

function parse_corr!(nonmissing::BitVector, corr::CrossCorrelated)

    (nonmissing == corr.nonmissing) && (return copy(corr))

    nonmissing .*= corr.nonmissing
    touse        = nonmissing[corr.nonmissing]
    new_mat      = Symmetric(corr.mat[touse, touse])

    return CrossCorrelated(corr.adj, nonmissing, SuiteSparse.CHOLMOD.Sparse(new_mat))
end

#==========================================================================================#

# RETRIEVAL

function getvector(obj::Microdata, x::Symbol)
    (length(obj.mapping[x]) == 1) || throw(ArgumentError(string(x) * " is not a vector"))
    return view(obj.matrix, :, obj.mapping[x]...)
end

function getmatrix(obj::Microdata, args...)
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    return view(obj.matrix, :, x)
end

function getnames(obj::Microdata, args...)
    n = coefnames(eterms(obj.model))
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    return n[x]
end

getcorr(obj::Microdata)       = obj.corr
getnonmissing(obj::Microdata) = obj.nonmissing
getweights(obj::Microdata)    = obj.weights