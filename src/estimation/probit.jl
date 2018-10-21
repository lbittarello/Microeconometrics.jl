#==========================================================================================#

# TYPE

mutable struct Probit <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    Probit() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Probit(MD::Microdata)
    obj        = Probit()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Probit, ::UnitWeights)

    y  = iszero.(getvector(obj, :response))
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    p  = mean(y)
    p  = 1.0 / normpdf(norminvcdf(p))
    β₀ = lmul!(p, x \ y)

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi) in zip(y, μ)
            ll -= (yi ? normlogccdf(μi) : normlogcdf(μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            μ[i] = - normpdf(μi) / ηi
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            μ[i] = - normpdf(μi) / ηi
            ll  -= (yi ? log(- ηi) : log(ηi))
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            ηi   = normpdf(μi) / ηi
            μ[i] = muladd(μi, ηi, abs2(ηi))
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

function _fit!(obj::Probit, w::AbstractWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    p  = mean(y)
    p  = 1.0 / normpdf(norminvcdf(p))
    β₀ = lmul!(p, x \ y)

    function L(β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (yi, μi, wi) in zip(y, μ, w)
            ll -= wi * (yi ? normlogccdf(μi) : normlogcdf(μi))
        end

        return ll
    end

    function G!(g::Vector, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            μ[i] = - wi * normpdf(μi) / ηi
        end

        mul!(g, transpose(x), μ)
    end

    function LG!(g::Vector, β::Vector)

        mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            μ[i] = - wi * normpdf(μi) / ηi
            ll  -= wi * (yi ? log(- ηi) : log(ηi))
        end

        mul!(g, transpose(x), μ)

        return ll
    end

    function H!(h::Matrix, β::Vector)

        mul!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = (yi ? (normcdf(μi) - 1.0) : normcdf(μi))
            ηi   = normpdf(μi) / ηi
            μ[i] = wi * muladd(μi, ηi, abs2(ηi))
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

function score(obj::Probit)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (yi, vi)) in enumerate(zip(y, v))
        ηi   = (iszero(yi) ? (normcdf(vi) - 1.0) : normcdf(vi))
        v[i] = normpdf(vi) / ηi
    end

    return Diagonal(v) * x
end

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Probit, ::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (yi, vi)) in enumerate(zip(y, v))
        ηi   = (iszero(yi) ? (normcdf(vi) - 1.0) : normcdf(vi))
        ηi   = normpdf(vi) / ηi
        v[i] = muladd(vi, ηi, abs2(ηi))
    end

    return - crossprod(x, v)
end

function jacobian(obj::Probit, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :control)
    v = x * obj.β

    @inbounds for (i, (yi, vi, wi)) in enumerate(zip(y, v, w))
        ηi   = (iszero(yi) ? (normcdf(vi) - 1.0) : normcdf(vi))
        ηi   = normpdf(vi) / ηi
        v[i] = wi * muladd(vi, ηi, abs2(ηi))
    end

    return - crossprod(x, v)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Probit, MD::Microdata)
    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")
    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::Probit, MD::Microdata) = normcdf.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Probit)
    x  = getmatrix(obj, :control)
    v  = x * obj.β
    v .= normpdf.(v)
    return Diagonal(v) * x
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Probit) = getnames(obj, :control)
mtitle(obj::Probit)    = "Probit MLE"

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Probit, ::UnitWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi) in zip(y, μ)
        ll += (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
    end

    return ll
end

function _loglikelihood(obj::Probit, w::AbstractWeights)

    y  = getvector(obj, :response)
    μ  = predict(obj)
    ll = 0.0

    @inbounds for (yi, μi, wi) in zip(y, μ, w)
        ll += wi * (iszero(yi) ? normlogccdf(μi) : normlogcdf(μi))
    end

    return ll
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Probit, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    return nobs(obj) * (μ * log(μ) + (1.0 - μ) * log(1.0 - μ))
end

# DEVIANCE

deviance(obj::Probit) = - 2.0 * _loglikelihood(obj, getweights(obj))

# DEVIANCE UNDER NULL MODEL

nulldeviance(obj::Probit) = - 2.0 * _nullloglikelihood(obj, getweights(obj))
