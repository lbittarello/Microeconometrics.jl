#==========================================================================================#

# TYPE

mutable struct Poisson <: MLE

    sample::Microdata
    β::Vector{Float64}
    V::AbstractMatrix{Float64}

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

function _fit!(obj::Poisson, ::UnitWeights)

    has_offset = haskey(obj.sample.mapping, :offset)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    r  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))

    has_offset && (xo = getmatrix(obj, :offset, :control))

    function L(β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        return sum(exp.(μ) .- y .* μ)
    end

    function G!(g::Vector, β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        r .= exp.(μ) .- y
        mul!(g, transpose(x), r)
    end

    function LG!(g::Vector, β::Vector)

        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = exp(μi)
            r[i] = ηi - yi
            ll  += ηi - yi * μi
        end

        mul!(g, transpose(x), r)

        return ll
    end

    function H!(h::AbstractMatrix, β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        xx .= x .* exp.(μ)
        mul!(h, transpose(x), xx)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("likelihood maximization did not converge")
    end
end

function _fit!(obj::Poisson, w::AbstractWeights)

    has_offset = haskey(obj.sample.mapping, :offset)

    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    r  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y, w)))

    has_offset && (xo = getmatrix(obj, :offset, :control))

    function L(β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        return sum(exp.(μ) .- y .* μ, w)
    end

    function G!(g::Vector, β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        r .= w .* (exp.(μ) .- y)
        mul!(g, transpose(x), r)
    end

    function LG!(g::Vector, β::Vector)

        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = exp(μi)
            r[i] = wi * (ηi - yi)
            ll  += wi * (ηi - yi * μi)
        end

        mul!(g, transpose(x), r)

        return ll
    end

    function H!(h::AbstractMatrix, β::Vector)
        has_offset ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        r  .= w .* exp.(μi)
        xx .= x .* r
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

score(obj::Poisson) = Diagonal(residuals(obj)) * getmatrix(obj, :control)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Poisson, ::UnitWeights)

    x = getmatrix(obj, :control)

    if haskey(obj.sample.mapping, :offset)
        xo = getmatrix(obj, :offset, :control)
        v  = xo * vcat(1.0, obj.β)
    else
        v  = x * obj.β
    end

    v .= exp.(v)

    return - crossprod(x, v)
end

function jacobian(obj::Poisson, w::AbstractWeights)

    x = getmatrix(obj, :control)

    if haskey(obj.sample.mapping, :offset)
        xo = getmatrix(obj, :offset, :control)
        v  = xo * vcat(1.0, obj.β)
    else
        v  = x * obj.β
    end

    v .= w .* exp.(v)

    return - crossprod(x, v)
end

#==========================================================================================#

# LINEAR PREDICTOR

function linear_predictor(obj::Poisson, MD::Microdata)

    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")

    if haskey(obj.sample.mapping, :offset)
        return getmatrix(MD, :offset, :control) * vcat(1.0, obj.β)
    else
        return getmatrix(MD, :control) * obj.β
    end
end

# FITTED VALUES

predict(obj::Poisson, MD::Microdata) = exp.(linear_predictor(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Poisson)

    if haskey(obj.sample.mapping, :offset)
        v = getmatrix(obj, :offset, :control) * vcat(1.0, obj.β)
    else
        v = getmatrix(obj, :control) * obj.β
    end

    v .= exp.(v)

    return Diag(v) * getmatrix(obj, :control)
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Poisson) = getnames(obj, :control)

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = linear_predictor(obj)
    return sum(y .* μ .- exp.(μ) .- loggamma.(1.0 .+ y), w)
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Poisson, w::AbstractWeights)
    if haskey(obj.sample.mapping, :offset)
        y = getvector(obj, :response)
        o = getmatrix(obj, :offset)
        β = log(mean(y, w) / mean(exp.(o), w))
        return sum(y .* (o .+ β) .- exp.(o .+ β) .- loggamma.(1.0 .+ y), w)
    else
        y = getvector(obj, :response)
        η = mean(y, w)
        μ = log(η)
        return sum(y .* μ .- η .- loggamma.(1.0 .+ y), w)
    end
end

# DEVIANCE

function _deviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = linear_predictor(obj)
    return 2 * sum(xlogx.(y) .- y .* (1.0 .+ μ) .+ exp.(μ), w)
end

# DEVIANCE UNDER NULL MODEL

function _nulldeviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    return 2 * sum(xlogx.(y) .- y .* (1.0 .+ log(μ)) .+ μ, w)
end
