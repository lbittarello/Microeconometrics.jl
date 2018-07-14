abstract type Micromodel <: RegressionModel end

abstract type ParModel <: Micromodel end

abstract type MLE <: ParModel end

abstract type GMM <: Micromodel end

abstract type TwoStageModel <: Micromodel end

mutable struct ParObject

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParObject() = new()
end

const ParOr2Stage  = Union{ParModel, TwoStageModel, GMM}
const ParOrGMM     = Union{ParModel, GMM}
const ParObjects   = Union{ParModel, GMM, ParObject}
const MicroObjects = Union{Micromodel, ParObject}