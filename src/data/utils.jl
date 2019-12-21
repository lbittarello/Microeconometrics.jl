#==========================================================================================#

# PARSE WEIGHTS

for W in (AnalyticWeights, FrequencyWeights, ProbabilityWeights, Weights)
    @eval begin
        function update_weights(w::$W{T, Bool}, nonmissing::BitVector) where T
            UnitWeights{Float64}(sum(nonmissing))
        end
    end
end

function update_weights(w::FrequencyWeights, nonmissing::BitVector)
    v = w[nonmissing]
    s = Int(sum(v))
    FrequencyWeights(float(v), s)
end

function update_weights(w::UnitWeights, nonmissing::BitVector)
    UnitWeights{Float64}(sum(nonmissing))
end

function update_weights(w::AbstractWeights, nonmissing::BitVector)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    Weights(v, n)
end

#==========================================================================================#

# COPY

function Base.copy(corr::C) where {C <: CorrStructure}
    C((copy(getfield(corr, f)) for f in fieldnames(C))...)
end

# EQUALITY

function Base.isequal(corr₁::C, corr₂::C) where {C <: CorrStructure}
    mapreduce(i -> isequal(getfield.([corr₁, corr₂], i)...), *, fieldnames(C), init = true)
end

Base.isequal(corr₁::CorrStructure, corr₂::CorrStructure) = false

# INTERSECTION BETWEEN A FRAME AND CORRELATION STRUCTURE

update_corr(corr::Homoscedastic, nonmissing::BitVector)   = copy(corr)
update_corr(corr::Heteroscedastic, nonmissing::BitVector) = copy(corr)

function update_corr(corr::Clustered, nonmissing::BitVector)

    isequal(nonmissing, corr.nonmissing) && (return copy(corr))

    nonmissing    .*= corr.nonmissing
    touse           = nonmissing[corr.nonmissing]
    new_id_clusters = corr.id_clusters[touse]
    iter            = sort(unique(corr.id_clusters))
    new_iter        = sort(unique(new_id_clusters))
    new_n_clusters  = length(new_iter)
    new_matrix      = corr.matrix[findall((in)(new_iter), iter), touse]

    Clustered(corr.corrected, nonmissing, new_matrix, new_id_clusters, new_n_clusters)
end

function update_corr(corr::CrossCorrelated, nonmissing::BitVector)

    isequal(nonmissing, corr.nonmissing) && (return copy(corr))

    nonmissing .*= corr.nonmissing
    touse        = nonmissing[corr.nonmissing]
    new_matrix   = SuiteSparse.CHOLMOD.Sparse(Symmetric(corr.matrix[touse, touse]))

    CrossCorrelated(corr.corrected, nonmissing, new_matrix)
end

#==========================================================================================#

# RETRIEVAL

function getvector(obj::Microdata, x::Symbol)
    (length(obj.mapping[x]) == 1) || throw(ArgumentError(string(x) * " is not a vector"))
    view(obj.matrix, :, obj.mapping[x]...)
end

function getmatrix(obj::Microdata, args...)
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    view(obj.matrix, :, x)
end

function getnames(obj::Microdata, args...)
    n = coefnames(eterms(obj.model))
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    n[x]
end

getcorr(obj::Microdata)       = obj.corr
getnonmissing(obj::Microdata) = obj.nonmissing
getweights(obj::Microdata)    = obj.weights