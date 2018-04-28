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

function cc_time(
        df::DataFrame,
        x::Symbol,
        b::Real,
        k::Function;
        adj::Bool = true
    )

    msng  = BitVector(size(df, 1))
    msng .= .!ismissing.(df[x])
    xx    = Vector{Date}(df[x][msng])
    n     = sum(msng)
    idx₁  = Vector{Int}(1:n)
    idx₂  = Vector{Int}(1:n)
    val   = fill(1.0, n)

    kernel(z) = k(z / float(b))

    for i = 1:n
        for j = 1:(i - 1)
            w = Dates.value(xx[i] - xx[j])
            w = k(w)
            if w > 0.0
                push!(idx₁, j)
                push!(idx₂, i)
                push!(val, w)
            end
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(sparse(idx₁, idx₂, val)))
end

# CORRELATION ACROSS SPACE

function cc_space(
        df::DataFrame,
        y::Symbol,
        x::Symbol,
        b::Real,
        k::Function = parzen;
        adj::Bool = true
    )

    msng  = BitVector(size(df, 1))
    msng .= .!(ismissing.(df[y]) .& ismissing.(df[x]))
    yy    = Vector{Float64}(df[y][msng])
    xx    = Vector{Float64}(df[x][msng])
    n     = sum(msng)
    idx₁  = Vector{Int}(1:n)
    idx₂  = Vector{Int}(1:n)
    val   = fill(1.0, n)

    kernel(z) = k(z / float(b))

    for i = 1:n
        for j = 1:(i - 1)
            w = geodistance(yy[i], xx[i], yy[j], xx[j])
            w = kernel(w)
            if w > 0.0
                push!(idx₁, j)
                push!(idx₂, i)
                push!(val, w)
            end
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(sparse(idx₁, idx₂, val)))
end

# CORRELATION ACROSS TIME AND SPACE

function cc_timespace(
        df::DataFrame,
        x₁::Symbol,
        b₁::Real,
        y₂::Symbol,
        x₂::Symbol,
        b₂::Real,
        k::Function = parzen;
        adj::Bool = true
    )

    msng  = BitVector(size(df, 1))
    msng .= .!(ismissing.(df[x₁]) .& ismissing.(df[y₂]) .& ismissing.(df[x₂]))
    xx₁   = Vector{Date}(df[x₁][msng])
    yy₂   = Vector{Float64}(df[y₂][msng])
    xx₂   = Vector{Float64}(df[x₂][msng])
    n     = sum(msng)
    idx₁  = Vector{Int}(1:n)
    idx₂  = Vector{Int}(1:n)
    val   = fill(1.0, n)

    b₁ = float(b₁)
    b₂ = float(b₂)

    kernel(z₁, z₂) = k(sqrt(abs2(z₁) + abs2(z₂)))

    for i = 1:n
        for j = 1:(i - 1)
            w₁ = Dates.value(xx₁[i] - xx₁[j]) / b₁
            if w₁ <= 1.0
                w₂ = geodistance(yy₂[i], xx₂[i], yy₂[j], xx₂[j]) / b₂
                if w₂ <= 1.0
                    push!(idx₁, j)
                    push!(idx₂, i)
                    push!(val, kernel(w₁, w₂))
                end
            end
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(sparse(idx₁, idx₂, val)))
end
