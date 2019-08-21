#==========================================================================================#

# CROSS PRODUCTS

crossprod(x::AbstractMatrix)                          = x' * x
crossprod(x::AbstractMatrix, w::UnitWeights)          = x' * x
crossprod(x::AbstractMatrix, w::AbstractVector)       = x' * (Diagonal(w) * x)
crossprod(x::AbstractMatrix, w::AbstractSparseMatrix) = x' * Matrix(w * x)
crossprod(x::AbstractMatrix, w::AbstractMatrix)       = x' * (w * x)

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
reweight(w::AbstractWeights, v::ProbabilityWeights) = pweights(w .* v)

#==========================================================================================#

# FORMATTER

function frmtr(X::Array{<: Real}, d::Int)

    table = format.("{:.$(d)f}", X)

    for c = 1:size(table, 2)
        l = maximum(length.(table[:, c]))
        table[:, c] .= format.("{:>$(l)s}", table[:, c])
    end

    return table
end
