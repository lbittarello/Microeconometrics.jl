#==========================================================================================#

# CROSS PRODUCTS

function crossprod(x::AbstractMatrix; neg::Bool = false)
    neg ? scale!(- 1.0, x' * x) : (x' * x)
end

function crossprod(x::AbstractMatrix, w::AbstractVector; neg::Bool = false)
    neg ? scale!(- 1.0, x' * scale!(w, copy(x))) : (x' * scale!(w, copy(x)))
end

function crossprod(x::AbstractMatrix, w::AbstractMatrix; neg::Bool = false)
    neg ? scale!(- 1.0, x' * w * x) : (x' * w * x)
end
