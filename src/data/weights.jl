#==========================================================================================#

# UNIT WEIGHTS

struct UnitWeights{S <: Real, T <: Real, V <: AbstractVector{T}} <: AbstractWeights{S, T, V}
    values::V
    sum::S
end

uweights(vs::AbstractVector) = UnitWeights(vs, sum(vs))

#==========================================================================================#

# CONVERT WEIGHTS

function parse_weights(w::AnalyticWeights, nonmissing::BitVector)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    AnalyticWeights(v, n)
end

function parse_weights(w::FrequencyWeights, nonmissing::BitVector)
    v = w[nonmissing]
    s = Int(sum(v))
    FrequencyWeights(float(v), s)
end

function parse_weights(w::ProbabilityWeights, nonmissing::BitVector)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    ProbabilityWeights(v, n)
end

function parse_weights(w::UnitWeights, nonmissing::BitVector)
    UnitWeights(float(w[nonmissing]), sum(nonmissing))
end

function parse_weights(w::Weights, nonmissing::BitVector)
    n  = sum(nonmissing)
    v  = float(w[nonmissing])
    s  = n / sum(v)
    v .= v .* s
    Weights(v, n)
end

#==========================================================================================#

# REWEIGHTING

reweight(w::UnitWeights, v::ProbabilityWeights)     = v
reweight(w::AbstractWeights, v::ProbabilityWeights) = pweights(w .* v)

#==========================================================================================#

# OPERATIONS

sum(v::AbstractArray{T,N}  where N where T<:Real, ::UnitWeights) = sum(v)
mean(v::AbstractArray{T,N} where N where T<:Real, ::UnitWeights) = mean(v)

function isequal(x::UnitWeights, y::UnitWeights)
    return isequal(x.sum, y.sum) && isequal(x.values, y.values)
end
