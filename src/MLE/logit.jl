#==========================================================================================#

# TYPE

mutable struct Logit <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    Logit() = new()
end

#==========================================================================================#

# ESTIMATION

function _fit(obj::Logit)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    p  = mean(y)
    β₀ = (x \ y) / (p * (1.0 - p))

    μ  = x * β₀
    r  = similar(y)
    ω  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ll -= (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
        end

        return - ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = logistic(μi)
            r[i] = yi - ηi
        end

        g[:] = - x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = logistic(μi)
            r[i] = yi - ηi
            ll  += (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        g[:] = - x' * r

        return - ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, μi) in enumerate(μ)
            ηi   = logistic(μi)
            ω[i] = ηi * (1.0 - ηi)
        end

        h[:, :] = x' * (x .* ω)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!), β₀, Newton())

    if Optim.converged(res)
        return Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

function _fit(obj::Logit, w::AbstractVector)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    p  = mean(y)
    β₀ = (x \ y) / (p * (1.0 - p))

    μ  = x * β₀
    r  = similar(y)
    ω  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ll -= wi * (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        return - ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = logistic(μi)
            r[i] = wi * (yi - ηi)
        end

        g[:] = - x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = logistic(μi)
            r[i] = wi * (yi - ηi)
            ll  += wi * (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        g[:] = - x' * r

        return - ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (μi, wi)) in enumerate(zip(μ, w))
            ηi   = logistic(μi)
            ω[i] = wi * ηi * (1.0 - ηi)
        end

        h[:, :] = x' * (x .* ω)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!), β₀, Newton())

    if Optim.converged(res)
        return Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

#==========================================================================================#

# SCORE (DERIVATIVE OF THE LIKELIHOOD FUNCTION)

score(obj::Logit) = getmatrix(obj, :control) .* residuals(obj)

function score(obj::Logit, w::AbstractVector)
    x  = getmatrix(obj, :control)
    r  = residuals(obj)
    r .= r .* w
    return x .* r
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Logit)

    x  = getmatrix(obj, :control)
    η  = x * obj.β

    @inbounds for (i, ηi) in enumerate(η)
        ψ     = logistic(ηi)
        η[i] .= ψ * (1.0 - ψ)
    end

    return - x' * (x .* η)
end

function jacobian(obj::Logit, w::AbstractVector)

    x  = getmatrix(obj, :control)
    η  = x * obj.β

    @inbounds for (i, (ηi, wi)) in enumerate(zip(η, w))
        ψ     = logistic(ηi)
        η[i] .= wi * ψ * (1.0 - ψ)
    end

    return - x' * (x .* η)
end

#==========================================================================================#

# LINEAR PREDICTOR

predict(obj::Logit) = getmatrix(obj, :control) * obj.β

# FITTED VALUES

fitted(obj::Logit)  = logistic.(predict(obj))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Logit)

    x  = getmatrix(obj, :control)
    η  = x * obj.β

    @inbounds for (i, ηi) in enumerate(η)
        ψ     = logistic(ηi)
        η[i] .= ψ * (1.0 - ψ)
    end

    return x .* η
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Logit) = getnames(obj, :control)

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Logit)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi) in zip(y, μ)
        ll -= (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
    end

    return ll
end

function _loglikelihood(obj::Logit, w::AbstractVector)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    for (yi, μi, wi) in zip(y, μ, w)
        ll -= wi * (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
    end

    return ll
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Logit)
    y = getvector(obj, :response)
    μ = mean(y)
    return length(y) * (μ * log(μ) + (1.0 - μ) * log(1.0 - μ))
end

function _nullloglikelihood(obj::Logit, w::AbstractVector)
    y  = getvector(obj, :response)
    w0 = view(w, y .== 0.0)
    w1 = view(w, y .== 1.0)
    μ  = mean(w1)
    return sum(w1) * log(μ) + sum(w0) * log(1.0 - μ)
end

# DEVIANCE

_deviance(obj::Logit)                        = - 2.0 * _loglikelihood(obj)
_deviance(obj::Logit, w::AbstractVector)     = - 2.0 * _loglikelihood(obj, w)

# DEVIANCE UNDER NULL MODEL

_nulldeviance(obj::Logit)                    = - 2.0 * _nullloglikelihood(obj)
_nulldeviance(obj::Logit, w::AbstractVector) = - 2.0 * _nullloglikelihood(obj, w)
