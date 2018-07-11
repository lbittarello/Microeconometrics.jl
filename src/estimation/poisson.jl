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

    O  = haskey(obj.sample.map, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    r  = similar(y)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))

    if O
        xo = getmatrix(obj, :offset, :control)
        μ  = xo * vcat(0.0, β₀)
    else
        μ = x * β₀
    end

    function L(β::Vector)
        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
        return sum(exp.(μ) .- y .* μ)
    end

    function G!(g::Vector, β::Vector)

        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi)) in enumerate(zip(y, μ))
            r[i] = exp(μi) - yi
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
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
        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
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

    O  = haskey(obj.sample.map, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :control)
    r  = similar(y)

    β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y, w)))

    if O
        xo = getmatrix(obj, :offset, :control)
        μ  = xo * vcat(0.0, β₀)
    else
        μ = x * β₀
    end

    function L(β::Vector)
        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
        return sum(exp.(μ) .- y .* μ, w)
    end

    function G!(g::Vector, β::Vector)

        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)

        @inbounds for (i, (yi, μi, wi)) in enumerate(zip(y, μ, w))
            r[i] = wi * (exp(μi) - yi)
        end

        g[:] = x' * r
    end

    function LG!(g::Vector, β::Vector)

        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
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
        O ? A_mul_B!(μ, xo, vcat(1.0, β)) : A_mul_B!(μ, x, β)
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
    if haskey(obj.sample.map, :offset)
        x  = getmatrix(obj, :control)
        xo = getmatrix(obj, :offset, :control)
        v  = xo * vcat(1.0, obj.β)
        v .= exp.(v)
        return crossprod(x, v, neg = true)
    else
        x  = getmatrix(obj, :control)
        v  = x * obj.β
        v .= exp.(v)
        return crossprod(x, v, neg = true)
    end
end

function jacobian(obj::Poisson, w::AbstractWeights)
    if haskey(obj.sample.map, :offset)
        x  = getmatrix(obj, :control)
        xo = getmatrix(obj, :offset, :control)
        v  = xo * vcat(1.0, obj.β)
        v .= w .* exp.(v)
        return crossprod(x, v, neg = true)
    else
        x  = getmatrix(obj, :control)
        v  = x * obj.β
        v .= w .* exp.(v)
        return crossprod(x, v, neg = true)
    end
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Poisson, MD::Microdata)
    if getnames(obj, :control) != getnames(MD, :control)
        throw("some variables are missing")
    end
    if haskey(MD.map, :offset)
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
        x  = getmatrix(obj, :control)
        xo = getmatrix(obj, :offset, :control)
        v  = xo * vcat(1.0, obj.β)
        v .= exp.(v)
        return scale!(v, copy(x))
    else
        x  = getmatrix(obj, :control)
        v  = x * obj.β
        v .= exp.(v)
        return scale!(v, copy(x))
    end
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

function _nullloglikelihood(obj::Poisson, w::UnitWeights)
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
    return 2 * sum(xlogx.(y) .- y .* (1.0 + μ) .+ exp.(μ), w)
end

# DEVIANCE UNDER NULL MODEL

function _nulldeviance(obj::Poisson, w::AbstractWeights)
    y = getvector(obj, :response)
    η = mean(y, w)
    λ = 1.0 + log(μ)
    return 2 * sum(xlogx.(y) .- y .* λ .+ η, w)
end
