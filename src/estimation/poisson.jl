#==========================================================================================#

# TYPE

mutable struct Poisson <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}

    Poisson() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Poisson(MD::Microdata)
    obj        = Poisson()
    obj.sample = MD
    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Poisson, w::UnitWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))

    μ  = x * β₀
    r  = similar(μ)

    function L(β::Vector)
        A_mul_B!(μ, x, β)
        return sum(exp(μ) .- y .* μ)
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            r[i] = exp(μi) - yi
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = exp(μi)
            r[i] = ηi - yi
            ll  += ηi - yi * μi
        end

        g[:] = x' * r

        return ll
    end

    function H!(h::Matrix, β::Vector)
        A_mul_B!(μ, x, β)
        r .= exp.(μ)
        h[:, :] = crossprod(x, r)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

function _fit!(obj::Poisson, w::AbstractWeights)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y, w)))

    μ  = x * β₀
    r  = similar(μ)

    function L(β::Vector)
        A_mul_B!(μ, x, β)
        return sum(exp(μ) .- y .* μ, w)
    end

    function G!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            r[i] = wi * (exp(μi) - yi)
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        A_mul_B!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = exp(μi)
            r[i] = wi * (ηi - yi)
            ll  += wi * (ηi - yi * μi)
        end

        g[:] = x' * r

        return ll
    end

    function H!(h::Matrix, β::Vector)
        A_mul_B!(μ, x, β)
        r .= w .* exp.(μi)
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

score(obj::Poisson) = scale!(residuals(obj), copy(getmatrix(obj, :control)))

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Poisson, w::UnitWeights)
    x  = getmatrix(obj, :control)
    v  = x * obj.β
    v .= exp.(v)
    return crossprod(x, v, neg = true)
end

function jacobian(obj::Poisson, w::AbstractWeights)
    x  = getmatrix(obj, :control)
    v  = x * obj.β
    v .= w .* exp.(v)
    return crossprod(x, v, neg = true)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Poisson, MD::Microdata)

    if getnames(obj, :control) != getnames(MD, :control)
        throw("some variables are missing")
    end

    getmatrix(MD, :control) * obj.β
end

# FITTED VALUES

fitted(obj::Poisson, MD::Microdata) = exp.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Poisson)
    x = getmatrix(obj, :control)
    v = x * obj.β
    v = exp.(v)
    return scale!(v, copy(x))
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Poisson) = getnames(obj, :control)

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = predict(obj)
    return sum(y .* μ .- exp.(μ) .- lgamma.(1.0 .+ y), w)
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    λ = log(μ)
    return sum(y .* λ .- μ .- lgamma.(1.0 .+ y), w)
end

# DEVIANCE

function _deviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = predict(obj)
    return 2 * sum(xlogx.(y) .- y .* μ .+ exp.(μ), w)
end

# DEVIANCE UNDER NULL MODEL

function _nulldeviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    λ = log(μ)
    return 2 * sum(xlogx.(y) .- y .* λ .+ μ, w)
end