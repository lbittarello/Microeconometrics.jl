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

function _fit!(obj::Poisson, ::UnitWeights)

    O  = haskey(obj.sample.map, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    r  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))

    O && (xo = getmatrix(obj, :offset, :control))

    function L(β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        return sum(exp.(μ) .- y .* μ)
    end

    function G!(g::Vector, β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        r .= exp.(μ) .- y
        mul!(g, transpose(x), r)
    end

    function LG!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            ηi   = exp(μi)
            r[i] = ηi - yi
            ll  += ηi - yi * μi
        end

        mul!(g, transpose(x), r)

        return ll
    end

    function H!(h::Matrix, β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
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

    O  = haskey(obj.sample.map, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    μ  = Array{Float64}(undef, length(y))
    r  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y, w)))

    O && (xo = getmatrix(obj, :offset, :control))

    function L(β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        return sum(exp.(μ) .- y .* μ, w)
    end

    function G!(g::Vector, β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        r .= w .* (exp.(μ) .- y)
        mul!(g, transpose(x), r)
    end

    function LG!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
        ll = 0.0

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            ηi   = exp(μi)
            r[i] = wi * (ηi - yi)
            ll  += wi * (ηi - yi * μi)
        end

        mul!(g, transpose(x), r)

        return ll
    end

    function H!(h::Matrix, β::Vector)
        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)
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

    if haskey(obj.sample.map, :offset)
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

    if haskey(obj.sample.map, :offset)
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

function predict(obj::Poisson, MD::Microdata)

    (getnames(obj, :control) != getnames(MD, :control)) && throw("missing variables")

    if haskey(obj.sample.map, :offset)
        return getmatrix(MD, :offset, :control) * vcat(1.0, obj.β)
    else
        return getmatrix(MD, :control) * obj.β
    end
end

# FITTED VALUES

fitted(obj::Poisson, MD::Microdata) = exp.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Poisson)

    if haskey(obj.sample.map, :offset)
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
mtitle(obj::Poisson)    = "Poisson MLE"

# LIKELIHOOD FUNCTION

function _loglikelihood(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = predict(obj)
    return sum(y .* μ .- exp.(μ) .- lgamma.(1.0 .+ y), w)
end

# LIKELIHOOD FUNCTION UNDER NULL MODEL

function _nullloglikelihood(obj::Poisson, ::UnitWeights)

    if haskey(obj.sample.map, :offset)

        y  = getvector(obj, :response)
        o  = getmatrix(obj, :offset)
        β₀ = [log(mean(y))]

        function L(β::Vector)
            β = β[1]
            return sum(exp.(o .+ β) .- y .* (o .+ β))
        end

        function G!(g::Vector, β::Vector)
            β    = β[1]
            g[1] = sum(exp.(o .+ β) .- y)
        end

        function LG!(g::Vector, β::Vector)

            β    = β[1]
            ll   = 0.0
            g[1] = 0.0

            @inbounds for (i, (yi, oi)) in enumerate(zip(y, o))
                μi    = oi + β
                ηi    = exp(μi)
                g[1] += ηi - yi
                ll   += ηi - yi * μi
            end

            return ll
        end

        function H!(h::Matrix, β::Vector)
            β       = β[1]
            h[1, 1] = sum(exp.(o .+ β))
        end

        res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())
        β   = Optim.minimizer(res)

        return sum(y .* (o .+ β) .- exp.(o .+ β) .- lgamma.(1.0 .+ y))

    else

        y = getvector(obj, :response)
        η = mean(y, w)
        μ = log(η)
        return sum(y .* μ .- η .- lgamma.(1.0 .+ y))
    end
end

function _nullloglikelihood(obj::Poisson, w::AbstractWeights)

    if haskey(obj.sample.map, :offset)

        y  = getvector(obj, :response)
        o  = getmatrix(obj, :offset)
        β₀ = [log(mean(y, w))]

        function L(β::Vector)
            β = β[1]
            return sum(exp.(o .+ β) .- y .* (o .+ β), w)
        end

        function G!(g::Vector, β::Vector)
            β    = β[1]
            g[1] = sum(exp.(o .+ β) .- y, w)
        end

        function LG!(g::Vector, β::Vector)

            β    = β[1]
            ll   = 0.0
            g[1] = 0.0

            @inbounds for (i, (yi, oi, wi)) in enumerate(zip(y, o, w))
                μi    = oi + β
                ηi    = exp(μi)
                g[1] += wi * (ηi - yi)
                ll   += wi * (ηi - yi * μi)
            end

            return ll
        end

        function H!(h::Matrix, β::Vector)
            β       = β[1]
            h[1, 1] = sum(exp.(o .+ β), w)
        end

        res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())
        β   = Optim.minimizer(res)

        return sum(y .* (o .+ β) .- exp.(o .+ β) .- lgamma.(1.0 .+ y), w)

    else

        y = getvector(obj, :response)
        η = mean(y, w)
        μ = log(η)
        return sum(y .* μ .- η .- lgamma.(1.0 .+ y), w)
    end
end

# DEVIANCE

function _deviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = predict(obj)
    return 2 * sum(xlogx.(y) .- y .* (1.0 .+ μ) .+ exp.(μ), w)
end

# DEVIANCE UNDER NULL MODEL

function _nulldeviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    μ = mean(y, w)
    return 2 * sum(xlogx.(y) .- y .* (1.0 .+ log(μ)) .+ μ, w)
end
