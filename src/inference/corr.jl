#==========================================================================================#

# VCE TYPES

abstract type CorrStructure
end

mutable struct Homoscedastic <: CorrStructure
    adj::Bool
    method::String
end

mutable struct Heteroscedastic <: CorrStructure
    adj::Bool
end

mutable struct Clustered <: CorrStructure
    adj::Bool
    msng::BitVector
    mat::AbstractMatrix
    ic::AbstractVector
    nc::Int
end

mutable struct CrossCorrelated <: CorrStructure
    adj::Bool
    msng::BitVector
    mat::AbstractMatrix
end

ClusterOrCross = Union{Clustered, CrossCorrelated}

#==========================================================================================#

# HOMOSCEDASTIC

Homoscedastic(method::String = "OIM"; adj::Bool = true) = Homoscedastic(adj, method)

# HETEROSCEDASTIC

Heteroscedastic(; adj::Bool = true) = Heteroscedastic(adj)

#==========================================================================================#

# CLUSTERED

function Clustered(df::DataFrame, x::Symbol; adj::Bool = true)

    msng  = BitVector(size(df, 1))
    msng .= .!ismissing.(df[x])
    ic    = Array(df[x][msng])
    n     = sum(msng)
    iter  = unique(ic)
    nc    = length(iter)
    idx₁  = Vector{Int}(1:n)
    idx₂  = Vector{Int}(1:n)

    @inbounds for i in iter
        idx = findin(ic, [i])
        for j = 1:length(idx)
            for k = 1:(j - 1)
                push!(idx₁, idx[k])
                push!(idx₂, idx[j])
            end
        end
    end

    val = fill(1.0, length(idx₁))

    return Clustered(adj, msng, Symmetric(sparse(idx₁, idx₂, val)), ic, nc)
end

#==========================================================================================#

# CROSS CORRELATED: GENERAL CONSTRUCTOR

function CrossCorrelated(variant::String, args...; kwargs...)
    if variant == "Two-way clustering"
        return cc_twowayclustering(args...; kwargs...)
    elseif variant == "Time"
        return cc_time(args...; kwargs...)
    elseif variant == "Space"
        return cc_space(args...; kwargs...)
    elseif variant == "Time and space"
        return cc_timespace(args...; kwargs...)
    else
        throw("unknown variant")
    end
end

# TWOWAY CLUSTERING

function cc_twowayclustering(df::DataFrame, x₁::Symbol, x₂::Symbol; adj::Bool = true)

    msng   = BitVector(size(df, 1))
    msng  .= .!(ismissing.(df[x₁]) .& ismissing.(df[x₂]))
    ic₁    = df[x₁][msng]
    ic₂    = df[x₂][msng]
    n      = sum(msng)
    iter₁  = unique(ic₁)
    iter₂  = unique(ic₂)
    idx₁   = Vector{Int}(1:n)
    idx₂   = Vector{Int}(1:n)

    @inbounds for i in iter₁
        idx = findin(ic₁, [i])
        nix  = length(idx)
        for j = 1:nix
            for k = 1:(j - 1)
                push!(idx₁, idx[k])
                push!(idx₂, idx[j])
            end
        end
    end
    @inbounds for i in iter₂
        idx = findin(ic₂, [i])
        nix  = length(idx)
        for j = 1:nix
            for k = 1:(j - 1)
                push!(idx₁, idx[k])
                push!(idx₂, idx[j])
            end
        end
    end

    val = fill(1.0, length(idx₁))

    return CrossCorrelated(adj, msng, Symmetric(sparse(idx₁, idx₂, val, n, n, max)))
end

# CORRELATION ACROSS TIME

function cc_time(df::DataFrame, x::Symbol, b::Real; kwargs...)
    k(z) = parzen(z / float(b))
    return cc_time(df, x, k; kwargs...)
end

function cc_time(df::DataFrame, x::Symbol, k::Function; adj::Bool = true)

    msng  = BitVector(size(df, 1))
    msng .= .!ismissing.(df[x])
    xx    = df[x][msng] :: Vector{Date}
    n     = sum(msng)
    mat   = speye(Float64, n, n)

    for i = 1:n
        for j = 1:(i - 1)
            w = Dates.value(xx[i] - xx[j])
            w = k(w)
            (w > 0.0) && (mat[j, i] = w)
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(mat))
end

# CORRELATION ACROSS SPACE

function cc_space(df::DataFrame, y::Symbol, x::Symbol, b::Real; kwargs...)
    k(z) = parzen(z / float(b))
    return cc_space(df, y, x, k; kwargs...)
end

function cc_space(df::DataFrame, y::Symbol, x::Symbol, k::Function; adj::Bool = true)

    msng  = BitVector(size(df, 1))
    msng .= .!(ismissing.(df[y]) .& ismissing.(df[x]))
    yy    = df[y][msng] :: Vector{Float64}
    xx    = df[x][msng] :: Vector{Float64}
    n     = sum(msng)
    mat   = speye(Float64, n, n)

    for i = 1:n
        for j = 1:(i - 1)
            w = geodistance(yy[i], xx[i], yy[j], xx[j])
            w = k(w)
            (w > 0.0) && (mat[j, i] = w)
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(mat))
end

# CORRELATION ACROSS TIME AND SPACE

function cc_timespace(
        df::DataFrame, x₁::Symbol, b₁::Real, y₂::Symbol, x₂::Symbol, b₂::Real; kwargs...
    )

    k₁(z) = parzen(z / float(b₁))
    k₂(z) = parzen(z / float(b₂))

    cc_timespace(df, x₁, k₁, y₂, x₂, k₂; kwargs...)
end

function cc_timespace(
        df::DataFrame,
        x₁::Symbol,
        k₁::Function,
        y₂::Symbol,
        x₂::Symbol,
        k₂::Function;
        adj::Bool = true
    )

    msng  = BitVector(size(df, 1))
    msng .= .!(ismissing.(df[x₁]) .& ismissing.(df[y₂]) .& ismissing.(df[x₂]))
    xx₁   = df[x₁][msng] :: Vector{Date}
    yy₂   = df[y₂][msng] :: Vector{Float64}
    xx₂   = df[x₂][msng] :: Vector{Float64}
    n     = sum(msng)
    mat   = speye(Float64, n, n)

    for i = 1:n
        for j = 1:(i - 1)
            w₁ = Dates.value(xx₁[i] - xx₁[j])
            w₁ = k₁(w₁)
            if w₁ > 0.0
                w₂ = geodistance(yy₂[i], xx₂[i], yy₂[j], xx₂[j])
                w₂ = k₂(w₂)
                if w₂ > 0.0
                    mat[j, i] = w₁ * w₂
                end
            end
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(mat))
end
