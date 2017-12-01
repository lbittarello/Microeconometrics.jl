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

# COPY

copy(corr::Homoscedastic)   = Homoscedastic(corr.adj, corr.method)
copy(corr::Heteroscedastic) = Heteroscedastic(corr.adj)

function copy(corr::CrossCorrelated)
    CrossCorrelated(corr.adj, copy(corr.msng), copy(corr.mat))
end

function copy(corr::Clustered)
    Clustered(corr.adj, copy(corr.msng), copy(corr.mat), copy(corr.ic), copy(corr.nc))
end

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
    iter  = unique(id)
    nc    = length(iter)
    mat   = spzeros(Float64, n, n)

    for i in iter
        ii = findin(id, [i])
        mat[ii, ii] .= 1.0
    end

    return Clustered(adj, msng, mat, ic, nc)
end

#==========================================================================================#

# CROSS CORRELATED: GENERAL CONSTRUCTOR

function CrossCorrelated(variant::String, args...; kwargs...)
    if variant == "twoway clustering"
        return cc_twowayclustering(args...; kwargs...)
    elseif variant == "time"
        return cc_time(args...; kwargs...)
    elseif variant == "space"
        return cc_space(args...; kwargs...)
    elseif variant == "time and space"
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
    mat    = spzeros(Float64, n, n)

    for i in iter₁
        ii = findin(id₁, [i])
        mat[ii, ii] .= 1.0
    end
    for i in iter₂
        ii = findin(id₂, [i])
        mat[ii, ii] .= 1.0
    end

    return CrossCorrelated(adj, msng, mat)
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
