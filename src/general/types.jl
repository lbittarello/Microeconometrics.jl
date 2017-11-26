abstract type Micromodel <: RegressionModel
end

abstract type ParModel <: Micromodel
end

abstract type MLE <: ParModel
end

abstract type TwoStageModel <: Micromodel
end

mutable struct ParObject

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParObject() = new()
end

const ParOrTwoStage = Union{ParModel, TwoStageModel}
