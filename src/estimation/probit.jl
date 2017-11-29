#==========================================================================================#

# TYPE

mutable struct Probit{T} <: MLE{T}

    sample::Microdata{T}
    β::Vector{Float64}
    V::Matrix{Float64}

    Probit{T}() where {T} = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Probit(MD::Microdata{T}) where {T}
    obj        = Probit{T}()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Probit)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    p  = mean(y)
    p  = 1.0 / normpdf(norminvcdf(p))
    β₀ = scale!(p, x \ y)

    μ  = x * β₀
    r  = similar(y)
    ω  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ll += (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
        end

        return - ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            r[i] = - normpdf(μi) / ηi
        end

        g[:] = x' * r

    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            r[i] = - normpdf(μi) / ηi
            ll  += (iszero(yi) ? log(- ηi) : log(ηi))
        end

        g[:] = x' * r

        return - ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            ηi   = normpdf(μi) / ηi
            ω[i] = abs2(ηi) + μi * ηi
        end

        h[:, :] = crossprod(x, ω)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

function _fit!(obj::Probit, w::AbstractVector)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    p  = mean(y)
    p  = 1.0 / normpdf(norminvcdf(p))
    β₀ = scale!(p, x \ y)

    μ  = x * β₀
    r  = similar(y)
    ω  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ll += wi * (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
        end

        return - ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            r[i] = - wi * normpdf(μi) / ηi
        end

        g[:] = x' * r

    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            r[i] = - wi * normpdf(μi) / ηi
            ll  += wi * (iszero(yi) ? log(- ηi) : log(ηi))
        end

        g[:] = x' * r

        return - ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (iszero(yi) ? (normcdf(μi) - 1.0) : normcdf(μi))
            ηi   = normpdf(μi) / ηi
            ω[i] = wi * (abs2(ηi) + μi * ηi)
        end

        h[:, :] = crossprod(x, ω)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

#==========================================================================================#

# SCORE (DERIVATIVE OF THE LIKELIHOOD FUNCTION)

function score(obj::Probit)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    ω = x * obj.β

    @inbounds for (i, (yi, ωi)) in enumerate(zip(y, ω))
        ηi   = (iszero(yi) ? (normcdf(ωi) - 1.0) : normcdf(ωi))
        ω[i] = normpdf(ωi) / ηi
    end

    return scale!(ω, copy(x))
end

function score(obj::Probit, w::AbstractVector)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    ω = x * obj.β

    @inbounds for (i, (yi, ωi, wi)) in enumerate(zip(y, ω, w))
        ηi   = (iszero(yi) ? (normcdf(ωi) - 1.0) : normcdf(ωi))
        ω[i] = wi * normpdf(ωi) / ηi
    end

    return scale!(ω, copy(x))
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Probit)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    ω = x * obj.β

    @inbounds for (i, ωi) in enumerate(ω)
        ηi   = (iszero(yi) ? (normcdf(ωi) - 1.0) : normcdf(ωi))
        ηi   = normpdf(ωi) / ηi
        ω[i] = abs2(ηi) + ωi * ηi
    end

    return crossprod(x, ω, neg = true)
end

function jacobian(obj::Probit, w::AbstractVector)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    ω = x * obj.β

    @inbounds for (i, (ωi, wi)) in enumerate(zip(ω, w))
        ηi   = (iszero(yi) ? (normcdf(ωi) - 1.0) : normcdf(ωi))
        ηi   = normpdf(ωi) / ηi
        ω[i] = wi * (abs2(ηi) + ωi * ηi)
    end

    return crossprod(x, ω, neg = true)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::Probit) = getmatrix(obj, :control) * obj.β

# FITTED VALUES

fitted(obj::Probit) = normcdf.(predict(obj))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Probit)
    x  = copy(getmatrix(obj, :control))
    ϕ  = x * obj.β
    ϕ .= normpdf.(ϕ)
    return scale!(ϕ, x)
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Probit) = getnames(obj, :control)

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Probit)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi) in zip(y, μ)
        ll += (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
    end

    return ll
end

function _loglikelihood(obj::Probit, w::AbstractVector)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi, wi) in zip(y, μ, w)
        ll += wi * (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
    end

    return ll
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Probit)
    y = getvector(obj, :response)
    μ = mean(y)
    return length(y) * (μ * log(μ) + (1.0 - μ) * log(1.0 - μ))
end

function _nullloglikelihood(obj::Probit, w::AbstractVector)
    y  = getvector(obj, :response)
    w₀ = view(w, y .== 0.0)
    w₁ = view(w, y .== 1.0)
    s₁ = sum(w₁)
    μ₁ = mean(w₁) / s₁
    return sum(w₁) * log(μ₁) + sum(w₀) * log(1.0 - μ₁)
end

# DEVIANCE

_deviance(obj::Probit)                    = - 2.0 * _loglikelihood(obj)
_deviance(obj::Probit, w::AbstractVector) = - 2.0 * _loglikelihood(obj, w)

# DEVIANCE UNDER NULL MODEL

_nulldeviance(obj::Probit)                    = - 2.0 * _nullloglikelihood(obj)
_nulldeviance(obj::Probit, w::AbstractVector) = - 2.0 * _nullloglikelihood(obj, w)
