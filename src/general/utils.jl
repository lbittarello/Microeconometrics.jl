#==========================================================================================#

# CROSS PRODUCTS

crossprod(x::AbstractMatrix)                          = x' * x
crossprod(x::AbstractMatrix, w::UnitWeights)          = x' * x
crossprod(x::AbstractMatrix, w::AbstractVector)       = x' * (Diagonal(w) * x)
crossprod(x::AbstractMatrix, w::AbstractSparseMatrix) = x' * Matrix(w * x)
crossprod(x::AbstractMatrix, w::AbstractMatrix)       = x' * (w * x)

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
