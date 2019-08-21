#==========================================================================================#

# UNIT WEIGHTS

struct UnitWeights{S<:Real, T<:Real} <: AbstractWeights{S, T, V where V<:Vector{T}}

    el::T
    sum::S

    UnitWeights{S, T}(el::Type{<:Real}, s::Real) where {S, T} = new(one(el), s)
end

UnitWeights(::T, s::S)       where {S, T} = UnitWeights{S, T}(T, s)
UnitWeights(::Type{T}, s::S) where {S, T} = UnitWeights{S, T}(T, s)

Base.eltype(wv::UnitWeights{S, T})  where {S, T} = T
Base.values(wv::UnitWeights{S, T})  where {S, T} = fill(wv.el, length(wv))
Base.sum(wv::UnitWeights{S, T})     where {S, T} = wv.sum
Base.isempty(wv::UnitWeights{S, T}) where {S, T} = iszero(wv.sum)
Base.length(wv::UnitWeights{S, T})  where {S, T} = Int(wv.sum)
Base.size(wv::UnitWeights{S, T})    where {S, T} = Tuple(length(wv))

function Base.getindex(wv::UnitWeights{S, T}, i::Int) where {S, T}
    (i > 0) & (i <= length(wv)) ? wv.el : Base.throw_boundserror(wv, i)
end

function Base.getindex(wv::UnitWeights{S, T}, i::AbstractRange{<:Int}) where {S, T}
    if (minimum(i) > 0) & (maximum(i) <= length(wv))
        return fill(wv.el, length(i))
    else
        Base.throw_boundserror(wv, i)
    end
end

function Base.getindex(wv::UnitWeights{S, T}, i::AbstractArray{<:Int}) where {S, T}
    if (minimum(i) > 0) & (maximum(i) <= length(wv))
        return fill(wv.el, size(i))
    else
        Base.throw_boundserror(wv, i)
    end
end

Base.getindex(wv::UnitWeights{S, T}, ::Colon) where {S, T} = fill(wv.el, length(wv))

sum(v::AbstractArray, ::UnitWeights)  = sum(v)
mean(v::AbstractArray, ::UnitWeights) = mean(v)

isequal(x::UnitWeights, y::UnitWeights) = (isequal(x.sum, y.sum) && isequal(x.el, y.el))
