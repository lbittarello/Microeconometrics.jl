#==========================================================================================#

# TYPE

mutable struct Logit <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    Logit() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Logit(MD::Microdata)
    obj        = Logit()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Logit, w::UnitWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    p  = mean(y)
    p  = 1.0 / (p * (1.0 - p))
    β₀ = scale!(p, x \ y)

    μ  = x * β₀
    r  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ll += (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            r[i] = logistic(μi) - yi
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = logistic(μi)
            r[i] = ηi - yi
            ll  -= (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        g[:] = x' * r

        return ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, μi) in enumerate(μ)
            ηi   = logistic(μi)
            r[i] = ηi * (1.0 - ηi)
        end

        h[:, :] = crossprod(x, r)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

function _fit!(obj::Logit, w::AbstractWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    w  = values(w)

    p  = mean(y)
    p  = 1.0 / (p * (1.0 - p))
    β₀ = scale!(p, x \ y)

    μ  = x * β₀
    r  = similar(μ)

    function L(β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ll += wi * (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            r[i] = wi * (logistic(μi) - yi)
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = logistic(μi)
            r[i] = wi * (ηi - yi)
            ll  -= wi * (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        g[:] = x' * r

        return ll
    end

    function H!(h::Matrix, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (μi, wi)) in enumerate(zip(μ, w))
            ηi   = logistic(μi)
            r[i] = wi * ηi * (1.0 - ηi)
        end

        h[:, :] = crossprod(x, r)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

#==========================================================================================#

# SCORE (DERIVATIVE OF THE LIKELIHOOD FUNCTION)

score(obj::Logit) = scale!(residuals(obj), copy(getmatrix(obj, :control)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Logit, w::UnitWeights)

    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, vi) in enumerate(v)
        λ     = logistic(vi)
        v[i] .= λ * (1.0 - λ)
    end

    return crossprod(x, v, neg = true)
end

function jacobian(obj::Logit, w::AbstractWeights)

    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (vi, wi)) in enumerate(zip(v, values(w)))
        λ     = logistic(vi)
        v[i] .= wi * λ * (1.0 - λ)
    end

    return crossprod(x, v, neg = true)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Logit, MD::Microdata)

    if getnames(obj, :control) != getnames(MD, :control)
        throw("some variables are missing")
    end

    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::Logit, MD::Microdata) = logistic.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Logit)

    x = getmatrix(obj, :control)
    η = x * obj.β

    @inbounds for (i, ηi) in enumerate(η)
        λ     = logistic(ηi)
        η[i] .= λ * (1.0 - λ)
    end

    return scale!(η, copy(x))
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Logit) = getnames(obj, :control)

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Logit, w::UnitWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi) in zip(y, μ)
        ll -= (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
    end

    return ll
end

function _loglikelihood(obj::Logit, w::AbstractWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    for (yi, μi, wi) in zip(y, μ, values(w))
        ll -= wi * (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
    end

    return ll
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Logit, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    return nobs(obj) * (μ * log(μ) + (1.0 - μ) * log(1.0 - μ))
end

# DEVIANCE

deviance(obj::Logit) = - 2.0 * _loglikelihood(obj, getweights(obj))

# DEVIANCE UNDER NULL MODEL

nulldeviance(obj::Logit) = - 2.0 * _nullloglikelihood(obj, getweights(obj))
