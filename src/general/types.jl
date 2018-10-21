abstract type Micromodel <: RegressionModel end

abstract type ParModel <: Micromodel end

abstract type MLE <: ParModel end

abstract type GMM <: ParModel end

abstract type TwoStageModel <: Micromodel end

mutable struct ParEstimate

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParEstimate() = new()
end

const Par1S  = Union{ParModel, ParEstimate}
const Par2S  = Union{ParModel, TwoStageModel, ParEstimate}
const ParM2S = Union{ParModel, TwoStageModel}
