#==========================================================================================#

# UNIT WEIGHTS

struct UnitWeights{S <: Real, T <: Real, V <: AbstractVector{T}} <: AbstractWeights{S, T, V}
    values::V
    sum::S
end

uweights(vs::AbstractVector) = UnitWeights(vs)
uweights(vs::AbstractArray)  = UnitWeights(vec(vs))

function UnitWeights(vs::AbstractVector{T}, s::S = sum(vs)) where {S <: Real, T <: Real}
    UnitWeights{S, T, typeof(vs)}(vs, s)
end

#==========================================================================================#

# EQUALITY TEST

function Base.isequal(w₁::W, w₂::W) where {W <: AbstractWeights}
    return isequal(values(w₁), values(w₂)) & isequal(sum(w₁), sum(w₂))
end

Base.isequal(w₁::AbstractWeights, w₂::AbstractWeights) = false

#==========================================================================================#

# CONVERT WEIGHTS

function parse_weights(w::AnalyticWeights, msng::BitVector)
    n  = sum(msng)
    v  = float(w[msng])
    s  = n / sum(v)
    v .= v .* s
    AnalyticWeights(v, n)
end

function parse_weights(w::FrequencyWeights, msng::BitVector)
    v = w[msng]
    s = Int(sum(v))
    FrequencyWeights(float(v), s)
end

function parse_weights(w::ProbabilityWeights, msng::BitVector)
    n  = sum(msng)
    v  = float(w[msng])
    s  = n / sum(v)
    v .= v .* s
    ProbabilityWeights(v, n)
end

function parse_weights(w::UnitWeights, msng::BitVector)
    UnitWeights(float(w[msng]), sum(msng))
end

function parse_weights(w::Weights, msng::BitVector)
    n  = sum(msng)
    v  = float(w[msng])
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

Base.sum(v::AbstractArray, w::UnitWeights)  = sum(v)
Base.mean(v::AbstractArray, w::UnitWeights) = mean(v)
