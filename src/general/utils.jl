function crossprod(x::AbstractMatrix; neg = false)
    neg ? scale!(- 1.0, x' * x) : (x' * x)
end

function crossprod(x::AbstractMatrix, w::AbstractVector; neg = false)
    neg ? scale!(- 1.0, x' * scale!(w, copy(x))) : (x' * scale!(w, copy(x)))
end
