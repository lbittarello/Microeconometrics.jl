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

function _fit!(obj::Logit, ::UnitWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    p  = mean(y)
    p  = 1.0 / (p * (1.0 - p))
    β₀ = lmul!(p, x \ y)

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ll += (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            μ[i] = logistic(μi) - yi
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = logistic(μi)
            μ[i] = ηi - yi
            ll  -= (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, μi) in enumerate(μ)
            ηi   = logistic(μi)
            μ[i] = ηi * (1.0 - ηi)
        end

        xx .= x .* μ
        mul!(h, transpose(x), xx)
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
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    p  = mean(y)
    p  = 1.0 / (p * (1.0 - p))
    β₀ = lmul!(p, x \ y)

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ll += wi * (iszero(yi) ? log1pexp(μi) : log1pexp(- μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            μ[i] = wi * (logistic(μi) - yi)
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = logistic(μi)
            μ[i] = wi * (ηi - yi)
            ll  -= wi * (iszero(yi) ? log(1.0 - ηi) : log(ηi))
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (μi, wi)) in enumerate(zip(μ, w))
            ηi   = logistic(μi)
            μ[i] = wi * ηi * (1.0 - ηi)
        end

        xx .= x .* μ
        mul!(h, transpose(x), xx)
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

score(obj::Logit) = Diagonal(residuals(obj)) * getmatrix(obj, :control)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Logit, ::UnitWeights)

    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, vi) in enumerate(v)
        ηi   = logistic(vi)
        v[i] = ηi * (1.0 - ηi)
    end

    return - crossprod(x, v)
end

function jacobian(obj::Logit, w::AbstractWeights)

    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (vi, wi)) in enumerate(zip(v, w))
        ηi   = logistic(vi)
        v[i] = wi * ηi * (1.0 - ηi)
    end

    return - crossprod(x, v)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Logit, MD::Microdata)
    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")
    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::Logit, MD::Microdata) = logistic.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Logit)

    x = copy(getmatrix(obj, :control))
    v = x * obj.β

    @inbounds for (i, vi) in enumerate(v)
        ηi   = logistic(vi)
        v[i] = ηi * (1.0 - ηi)
    end

    return lmul!(Diagonal(v), x)
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Logit) = getnames(obj, :control)
mtitle(obj::Logit)    = "Logit MLE"

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Logit, ::UnitWeights)

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

    for (yi, μi, wi) in zip(y, μ, w)
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
