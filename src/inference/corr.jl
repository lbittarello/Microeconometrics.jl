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

    msng  = BitVector(length(df[x]))
    msng .= (isna.(df[x]) .== false)
    id    = Array(df[x][msng])
    n     = length(id)
    iter  = unique(id)
    nc    = length(iter)
    mat   = spzeros(Float64, n, n)

    for i in iter
        ii = findin(id, [i])
        mat[ii, ii] = 1.0
    end

    return Clustered(adj, msng, mat, id, nc)
end

#==========================================================================================#

# CORRELATION ACROSS TIME AND SPACE

function cc_timespace(
        df::DataFrame,
        x1::Symbol,
        b1::Real,
        y2::Symbol,
        x2::Symbol,
        b2::Real;
        k1::Function = parzen,
        k2::Function = parzen,
        adj::Bool = true
    )

    _timespace(df[x1], float(b1), df[y2], df[x2], float(b2), k1, k2, adj)
end

function _timespace(
        x1::Vector{Date},
        b1::Float64,
        y2::Vector{Float64},
        x2::Vector{Float64},
        b2::Float64,
        k1::Function,
        k2::Function,
        adj::Bool
    )

    msng  = Array{Bool}(length(x1))
    msng .= (ismissing.(x1) .* ismissing.(y2) .* ismissing.(x2) .== false)
    n     = sum(msng)
    mat   = speye(Float64, n, n)
    idx   = findin(msng, [true])

    for (i,ii) in enumerate(idx)
        for (j,jj) in enumerate(idx[1:(i - 1)])
            w1 = Dates.value(x1[ii] - x1[jj])
            w1 = k1(float(w1) / b1)
            if w1 > 0.0
                w2 = geodistance(y2[ii], x2[ii], y2[jj], x2[jj])
                w2 = k2(w2 / b2)
                if w2 > 0.0
                    mat[j, i] = w1 * w2
                end
            end
        end
    end

    return CrossCorrelated(adj, msng, Symmetric(mat))
end
