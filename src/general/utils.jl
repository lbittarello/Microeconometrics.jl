#==========================================================================================#

# CROSS PRODUCTS

function crossprod(x::AbstractMatrix; neg::Bool = false)
    neg ? scale!(- 1.0, x' * x) : (x' * x)
end

function crossprod(x::AbstractMatrix, w::AbstractVector; neg::Bool = false)
    neg ? scale!(- 1.0, x' * scale!(w, copy(x))) : (x' * scale!(w, copy(x)))
end

#==========================================================================================#

# DISTANCE ON THE GLOBE

havd(x::Float64) = 0.5 * (1.0 - cosd(x))

function distance(y0::Float64, x0::Float64, y1::Float64, x1::Float64)
    h = havd(y1 - y0) + havd(x1 - x0) * cosd(y0) * cosd(y1)
    return 12742.0 * asin(min(1.0, sqrt(h)))
end

#==========================================================================================#

# KERNELS

bartlett(x::Float64)     = (abs(x) >= 1.0) ? 0.0 : (1.0 - abs(x))
truncated(x::Float64)    = (abs(x) >= 1.0) ? 0.0 : 1.0
tukeyhanning(x::Float64) = (abs(x) >= 1.0) ? 0.0 : 0.5 * (1.0 + cospi(x))

function parzen(x::Float64)
    if abs(x) <= 0.5
        return 1.0 - 6.0 * abs2(x) + 6 * abs(x)^3
    elseif abs(x) <= 1.0
        return 2.0 * (1.0 - abs(x))^3
    else
        return 0.0
    end
end
