abstract type Micromodel <: RegressionModel
end

abstract type ParModel <: Micromodel
end

mutable struct ParObject

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParObject() = new()
end

abstract type TwoStageModel <: Micromodel
end

const ParOrTwoStage = Union{ParModel, TwoStageModel}
