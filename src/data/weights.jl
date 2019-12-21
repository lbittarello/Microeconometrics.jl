#==========================================================================================#

# UNIT WEIGHTS

struct UnitWeights{T<:Real} <: AbstractWeights{Int, T, V where V<:Vector{T}}
    len::Int
end

Base.eltype(wv::UnitWeights{T}) where T = T
Base.sum(wv::UnitWeights{T}) where T = convert(T, length(wv))
Base.isempty(wv::UnitWeights) = iszero(wv.len)
Base.length(wv::UnitWeights) = wv.len
Base.size(wv::UnitWeights) = Tuple(length(wv))

Base.convert(::Type{Vector}, wv::UnitWeights{T}) where {T} = ones(T, length(wv))

@propagate_inbounds function Base.getindex(wv::UnitWeights{T}, i::Integer) where T
    @boundscheck checkbounds(wv, i)
    one(T)
end

@propagate_inbounds function Base.getindex(wv::UnitWeights{T}, i::AbstractArray{<:Int}) where T
    @boundscheck checkbounds(wv, i)
    UnitWeights{T}(length(i))
end

Base.getindex(wv::UnitWeights{T}, ::Colon) where {T} = UnitWeights{T}(wv.len)

@inline function varcorrection(w::UnitWeights, corrected::Bool=false)
    corrected ? (1 / (w.len - 1)) : (1 / w.len)
end

function Base.sum(A::AbstractArray, w::UnitWeights; dims::Union{Colon, Int}=:)
    a = (dims === :) ? length(A) : size(A, dims)
    a != length(w) && throw(DimensionMismatch("Inconsistent array dimension."))
    return sum(A, dims=dims)
end

function mean(A::AbstractArray, w::UnitWeights; dims::Union{Colon, Int}=:)
    a = (dims === :) ? length(A) : size(A, dims)
    a != length(w) && throw(DimensionMismatch("Inconsistent array dimension."))
    return mean(A, dims=dims)
end

Base.isequal(x::UnitWeights, y::UnitWeights) = isequal(x.len, y.len)
Base.:(==)(x::UnitWeights, y::UnitWeights)   = (x.len == y.len)