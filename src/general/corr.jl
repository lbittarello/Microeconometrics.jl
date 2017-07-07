#==========================================================================================#

# VCE TYPES

abstract type CorrStructure
end

mutable struct Homoscedastic <: CorrStructure
    method::String
end

mutable struct Heteroscedastic <: CorrStructure
end

mutable struct Clustered <: CorrStructure
    msng::BitVector
    mat::AbstractMatrix
    ic::AbstractVector
    nc::Int
end

mutable struct CrossCorrelated <: CorrStructure
    msng::BitVector
    mat::AbstractMatrix
end

const ClusterOrCross = Union{Clustered, CrossCorrelated}

#==========================================================================================#

# COPY

copy(corr::Homoscedastic)   = Homoscedastic()
copy(corr::Heteroscedastic) = Heteroscedastic()
copy(corr::CrossCorrelated) = CrossCorrelated(copy(corr.msng), copy(corr.mat))

function copy(corr::Clustered)
    return Clustered(copy(corr.msng), copy(corr.mat), copy(corr.ic), copy(corr.nc))
end

#==========================================================================================#

# HOMOSCEDASTIC

Homoscedastic() = Homoscedastic("OIM")

#==========================================================================================#

# CLUSTERED

function Clustered(df::DataFrame, x::Symbol)

    msng  = Array{Bool}(length(df[x]))
    msng .= (isna.(df[x]) .== false)
    n     = sum(msng)
    dfx   = Array(df[x][msng])
    iter  = unique(dfx)
    mat   = spzeros(Float64, n, n)

    for i in iter
        ii = findin(dfx, [i])
        mat[ii, ii] = 1.0
    end

    return Clustered(msng, mat, dfx, length(iter))
end

# ADJUSTMENT FOR CLUSTERED COVARIANCE MATRICES

_adjcluster!(V::Matrix, corr::Clustered)     = (c  = corr.nc; scale!(c / (c - 1), V))
_adjcluster!(V::Matrix, corr::CorrStructure) = V

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
        k2::Function = parzen
    )

    _timespace(df[x1], float(b1), df[y2], df[x2], float(b2), k1, k2)
end

function _timespace(
        x1::DataVector{Date},
        b1::Float64,
        y2::DataVector{Float64},
        x2::DataVector{Float64},
        b2::Float64,
        k1::Function,
        k2::Function
    )

    msng  = Array{Bool}(length(x1))
    msng .= (isna.(x1) .* isna.(y2) .* isna.(x2) .== false)
    n     = sum(msng)
    mat   = speye(Float64, n, n)
    idx   = findin(msng, [true])

    for (i,ii) in enumerate(idx)
        for (j,jj) in enumerate(idx[1:(i - 1)])
            w1 = Dates.value(x1[ii] - x1[jj])
            w1 = k1(float(w1) / b1)
            if w1 > 0.0
                w2 = distance(y2[ii], x2[ii], y2[jj], x2[jj])
                w2 = k2(w2 / b2)
                if w2 > 0.0
                    mat[j, i] = w1 * w2
                end
            end
        end
    end

    return CrossCorrelated(msng, Symmetric(mat))
end
