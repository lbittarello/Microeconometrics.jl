#==========================================================================================#

# CROSS PRODUCTS

crossprod(x::AbstractMatrix)                    = x' * x
crossprod(x::AbstractMatrix, w::UnitWeights)    = x' * x
crossprod(x::AbstractMatrix, w::AbstractVector) = x' * Diagonal(w) * x
crossprod(x::AbstractMatrix, w::AbstractMatrix) = x' * w * x