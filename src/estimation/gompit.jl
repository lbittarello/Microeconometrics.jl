#==========================================================================================#

# TYPE

mutable struct Gompit <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    Gompit() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Gompit(MD::Microdata)
    obj        = Gompit()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Gompit, ::UnitWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(- log(mean(y))))

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ηi  = - exp(μi)
            ll -= (iszero(yi) ? log(1.0 - exp(ηi)) : ηi)
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = ηi * (eηi - yi) / (1.0 - eηi)
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = ηi * (eηi - yi) / (1.0 - eηi)
            ll  -= (iszero(yi) ? log(1.0 - eηi) : ηi)
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = ηi * ((1.0 - yi) * ((ηi - 1.0) * eηi + 1.0) / (1.0 - eηi)^2 - 1.0)
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

function _fit!(obj::Gompit, w::AbstractWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(- log(mean(y, w))))

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ηi  = - exp(μi)
            ll -= wi * (iszero(yi) ? log(1.0 - eηi) : ηi)
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = wi * ηi * (eηi - yi) / (1.0 - eηi)
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = wi * ηi * (eηi - yi) / (1.0 - eηi)
            ll  -= wi * (iszero(yi) ? log(1.0 - eηi) : ηi)
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = - exp(μi)
            eηi  = exp(ηi)
            μ[i] = wi * ηi * ((1.0 - yi) * ((ηi - 1.0) * eηi + 1.0) / (1.0 - eηi)^2 - 1.0)
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

function score(obj::Gompit)

    y = getvector(obj, :response)
    x = copy(getmatrix(obj, :control))
    v = x * obj.β

    @inbounds for (i, (yi, vi)) in enumerate(zip(y, v))
        ηi   = - exp(vi)
        eηi  = exp(ηi)
        v[i] = ηi * (yi - eηi) / (1.0 - eηi)
    end

    return lmul!(Diagonal(v), x)
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Gompit, ::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (yi, vi)) in enumerate(zip(y, v))
        ηi   = - exp(vi)
        eηi  = exp(ηi)
        v[i] = ηi * (1.0 - (1.0 - yi) * ((ηi - 1.0) * eηi + 1.0) / (1.0 - eηi)^2)
    end

    return crossprod(x, v)
end

function jacobian(obj::Gompit, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (yi, vi, wi)) in enumerate(zip(y, v, w))
        ηi   = - exp(vi)
        eηi  = exp(ηi)
        v[i] = wi * ηi * (1.0 - (1.0 - yi) * ((ηi - 1.0) * eηi + 1.0) / (1.0 - eηi)^2)
    end

    return crossprod(x, v)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Gompit, MD::Microdata)
    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")
    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::Gompit, MD::Microdata) = exp.(.- exp.(predict(obj, MD)))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Gompit)
    x  = copy(getmatrix(obj, :control))
    v  = x * obj.β
    v .= .- exp.(v .- exp.(v))
    return lmul!(Diagonal(v), x)
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Gompit) = getnames(obj, :control)
mtitle(obj::Gompit)    = "Gompit MLE"

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Gompit, ::UnitWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi) in zip(y, μ)
        ηi  = - exp(μi)
        ll += (iszero(yi) ? log(1.0 - exp(ηi)) : ηi)
    end

    return ll
end

function _loglikelihood(obj::Gompit, w::AbstractWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    for (yi, μi, wi) in zip(y, μ, w)
        ηi  = - exp(μi)
        ll += wi * (iszero(yi) ? log(1.0 - exp(ηi)) : ηi)
    end

    return ll
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Gompit, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    return nobs(obj) * (μ * log(μ) + (1.0 - μ) * log(1.0 - μ))
end

# DEVIANCE

deviance(obj::Gompit) = - 2.0 * _loglikelihood(obj, getweights(obj))

# DEVIANCE UNDER NULL MODEL

nulldeviance(obj::Gompit) = - 2.0 * _nullloglikelihood(obj, getweights(obj))
