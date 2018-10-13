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
    nonmissing::BitVector
    mat::AbstractMatrix
    ic::AbstractVector
    nc::Int
end

mutable struct CrossCorrelated <: CorrStructure
    adj::Bool
    nonmissing::BitVector
    mat::AbstractMatrix
end

#==========================================================================================#

# HOMOSCEDASTIC

Homoscedastic(method::String = "OIM"; adj::Bool = true) = Homoscedastic(adj, method)

# HETEROSCEDASTIC

Heteroscedastic(; adj::Bool = true) = Heteroscedastic(adj)

#==========================================================================================#

# CLUSTERED

function Clustered(df::DataFrame, x::Symbol; adj::Bool = true)

    nonmissing  = BitVector(undef, size(df, 1))
    nonmissing .= .!ismissing.(df[x])
    ic          = disallowmissing(df[x][nonmissing])
    n           = sum(nonmissing)
    iter        = sort(unique(ic))
    nc          = length(iter)
    idx₁        = Vector{Int}()
    idx₂        = Vector{Int}()

    @inbounds for (i, ci) in enumerate(iter)
        idx = findall((in)([ci]), ic)
        append!(idx₁, fill(i, length(idx)))
        append!(idx₂, idx)
    end

    val = fill(1.0, length(idx₁))

    return Clustered(adj, nonmissing, sparse(idx₁, idx₂, val), ic, nc)
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

    nonmissing  = BitVector(undef, size(df, 1))
    nonmissing .= .!(ismissing.(df[x₁]) .& ismissing.(df[x₂]))
    ic₁         = disallowmissing(df[x₁][nonmissing])
    ic₂         = disallowmissing(df[x₂][nonmissing])
    n           = sum(nonmissing)
    iter₁       = unique(ic₁)
    iter₂       = unique(ic₂)
    idx₁        = Vector{Int}(1:n)
    idx₂        = Vector{Int}(1:n)
    nn          = n

    @inbounds for i in iter₁

        idx = findall((in)([i]), ic₁)
        nix = length(idx)
        nel = Int(nix * (nix - 1) / 2)

        append!(idx₁, Vector{Int}(nel))
        append!(idx₂, Vector{Int}(nel))

        for j = 1:nix
            for k = 1:(j - 1)
                nn += 1
                idx₁[nn] = idx[k]
                idx₂[nn] = idx[j]
            end
        end
    end

    @inbounds for i in iter₂

        idx = findall((in)([i]), ic₂)
        nix = length(idx)
        nel = Int(nix * (nix - 1) / 2)

        append!(idx₁, Vector{Int}(nel))
        append!(idx₂, Vector{Int}(nel))

        for j = 1:nix
            for k = 1:(j - 1)
                nn += 1
                idx₁[nn] = idx[k]
                idx₂[nn] = idx[j]
            end
        end
    end

    val = fill(1.0, length(idx₁))
    V   = Symmetric(sparse(idx₁, idx₂, val, n, n, max))

    return CrossCorrelated(adj, nonmissing, SuiteSparse.CHOLMOD.Sparse(V))
end

# CORRELATION ACROSS TIME

function cc_time(
        df::DataFrame,
        x::Symbol,
        b::Real,
        k::Function = parzen;
        adj::Bool = true
    )

    nonmissing  = BitVector(undef, size(df, 1))
    nonmissing .= .!ismissing.(df[x])
    xx          = Vector{Date}(df[x][nonmissing])
    n           = sum(nonmissing)
    idx₁        = Vector{Int}(1:n)
    idx₂        = Vector{Int}(1:n)
    val         = fill(1.0, n)

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

    V = Symmetric(sparse(idx₁, idx₂, val))

    return CrossCorrelated(adj, nonmissing, SuiteSparse.CHOLMOD.Sparse(V))
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

    nonmissing  = BitVector(undef, size(df, 1))
    nonmissing .= .!(ismissing.(df[y]) .& ismissing.(df[x]))
    yy          = Vector{Float64}(df[y][nonmissing])
    xx          = Vector{Float64}(df[x][nonmissing])
    n           = sum(nonmissing)
    idx₁        = Vector{Int}(1:n)
    idx₂        = Vector{Int}(1:n)
    val         = fill(1.0, n)

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

    V = Symmetric(sparse(idx₁, idx₂, val))

    return CrossCorrelated(adj, nonmissing, SuiteSparse.CHOLMOD.Sparse(V))
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

    nonmissing  = BitVector(undef, size(df, 1))
    nonmissing .= .!(ismissing.(df[x₁]) .& ismissing.(df[y₂]) .& ismissing.(df[x₂]))
    xx₁         = Vector{Date}(df[x₁][nonmissing])
    yy₂         = Vector{Float64}(df[y₂][nonmissing])
    xx₂         = Vector{Float64}(df[x₂][nonmissing])
    n           = sum(nonmissing)
    idx₁        = Vector{Int}(1:n)
    idx₂        = Vector{Int}(1:n)
    val         = fill(1.0, n)

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

    V = Symmetric(sparse(idx₁, idx₂, val))

    return CrossCorrelated(adj, nonmissing, SuiteSparse.CHOLMOD.Sparse(V))
end
