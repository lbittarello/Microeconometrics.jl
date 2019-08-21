abstract type MLE <: RegressionModel end
abstract type GMM <: RegressionModel end
abstract type TwoStageModel <: RegressionModel end

mutable struct ParEstimate

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParEstimate() = new()
end

const OneStageModel = Union{MLE, GMM}
const ParModel      = Union{MLE, GMM, TwoStageModel}
const AnyModel      = Union{MLE, GMM, TwoStageModel}
const ParObject     = Union{MLE, GMM, TwoStageModel, ParEstimate}
