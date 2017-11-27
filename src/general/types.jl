abstract type Micromodel{T <: CorrStructure} <: RegressionModel end

abstract type ParModel{T} <: Micromodel{T} end

abstract type MLE{T} <: ParModel{T} end

abstract type TwoStageModel{T} <: Micromodel{T} end

mutable struct ParObject

    names::Vector{String}
    β::AbstractArray
    V::AbstractArray

    ParObject() = new()
end

ParOr2Stage{T} = Union{ParModel{T}, TwoStageModel{T} where T}
