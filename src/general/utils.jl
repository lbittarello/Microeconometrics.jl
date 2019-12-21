#==========================================================================================#

# CROSS PRODUCTS

crossprod(x::AbstractMatrix)                          = Symmetric(x' * x)
crossprod(x::AbstractMatrix, w::UnitWeights)          = Symmetric(x' * x)
crossprod(x::AbstractMatrix, w::AbstractVector)       = Symmetric(x' * (Diagonal(w) * x))
crossprod(x::AbstractMatrix, w::AbstractSparseMatrix) = Symmetric(x' * Matrix(w * x))
crossprod(x::AbstractMatrix, w::AbstractMatrix)       = Symmetric(x' * (w * x))

# SCALING OF SYMMETRIC MATRICES

function LinearAlgebra.lmul!(s::Number, X::Symmetric)
    n = size(X, 1)
    if X.uplo == "U"
        @inbounds @simd for i in 1:n
            X.data[i, i] *= s
            for j in 1:(i - 1)
                X.data[j, i] *= s
                X.data[i, j] = X.data[j, i]
            end
        end
    else
        @inbounds @simd for i in 1:n
            X.data[i, i] *= s
            for j in (i + 1):n
                X.data[j, i] *= s
                X.data[i, j] = X.data[j, i]
            end
        end
    end
    X
end

#==========================================================================================#

# DISTANCE ON THE GLOBE

havd(x::Float64) = 0.5 * (1.0 - cosd(x))

function geodistance(y0::Float64, x0::Float64, y1::Float64, x1::Float64)
    h = havd(y1 - y0) + havd(x1 - x0) * cosd(y0) * cosd(y1)
    return 12742.0 * asin(min(1.0, sqrt(h)))
end

#==========================================================================================#

# KERNELS

bartlett(x::Real)        = bartlett(float(x))
bartlett(x::Float64)     = (abs(x) >= 1.0) ? 0.0 : (1.0 - abs(x))
truncated(x::Real)       = truncated(float(x))
truncated(x::Float64)    = (abs(x) >= 1.0) ? 0.0 : 1.0
tukeyhanning(x::Real)    = tukeyhanning(float(x))
tukeyhanning(x::Float64) = (abs(x) >= 1.0) ? 0.0 : 0.5 * (1.0 + cospi(x))

parzen(x::Real) = parzen(float(x))

function parzen(x::Float64)
    if abs(x) <= 0.5
        return 1.0 - 6.0 * abs2(x) + 6 * abs(x)^3
    elseif abs(x) <= 1.0
        return 2.0 * (1.0 - abs(x))^3
    else
        return 0.0
    end
end

gallant(x::Real)    = parzen(float(x))
gallant(x::Float64) = parzen(x)

#==========================================================================================#

# FIRST AND SECOND STAGES

first_stage(obj::TwoStageModel)  = obj.first_stage
second_stage(obj::TwoStageModel) = obj.second_stage

switch_stage(obj::OneStageModel) = obj
switch_stage(obj::TwoStageModel) = obj.second_stage
switch_stage(obj::ParEstimate)   = obj

#==========================================================================================#

# RETRIEVAL

getvector(obj::ParModel, x::Symbol) = getvector(switch_stage(obj).sample, x)
getmatrix(obj::ParModel, args...)   = getmatrix(switch_stage(obj).sample, args...)
getnames(obj::ParModel, args...)    = getnames(switch_stage(obj).sample, args...)
getcorr(obj::ParModel)              = getcorr(switch_stage(obj).sample)
getnonmissing(obj::AnyModel)        = getnonmissing(switch_stage(obj).sample)
getweights(obj::AnyModel)           = getweights(switch_stage(obj).sample)

#==========================================================================================#

# REWEIGHTING

reweight(w::UnitWeights, v::ProbabilityWeights)     = v
reweight(w::AbstractWeights, v::ProbabilityWeights) = ProbabilityWeights(w .* v)

#==========================================================================================#

# FORMATTER

function formatter(X::Array{<: Real}, d::Int)

    table = format.("{:.$(d)f}", X)

    for c = 1:size(table, 2)
        l = maximum(length.(table[:, c]))
        table[:, c] .= format.("{:>$(l)s}", table[:, c])
    end

    return table
end
