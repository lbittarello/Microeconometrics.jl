abstract type Micromodel <: StatsBase.RegressionModel end

abstract type ParModel <: Micromodel end

abstract type MLE <: ParModel end

abstract type GMM <: Micromodel end

abstract type TwoStageModel <: Micromodel end

mutable struct ParEstimate

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParObject() = new()
end

const Par1S  = Union{ParModel, GMM, ParEstimate}
const Par2S  = Union{ParModel, GMM, TwoStageModel, ParEstimate}
const ParM1S = Union{ParModel, GMM}
const ParM2S = Union{ParModel, GMM, TwoStageModel}