#==========================================================================================#

# TYPE

mutable struct Mullahy <: GMM

    method::String
    sample::Microdata
    β::Vector{Float64}
    V::Matrix{Float64}
    W::Matrix{Float64}

    Mullahy() = new()
end

#==========================================================================================#

# CONSTRUCTOR

function Mullahy(MD::Microdata, W::Matrix{Float64}, method::String)
    obj        = Mullahy()
    obj.sample = MD
    obj.W      = W
    obj.method = method
    return obj
end

#==========================================================================================#

# INTERFACE

function fit(
        ::Type{Mullahy},
        MD::Microdata;
        novar::Bool = false,
        method::String = "One-step GMM"
    )

    if method == "Poisson"

        FSD                   = Microdata(MD)
        FSD.mapping[:control] = vcat(MD.mapping[:treatment], MD.mapping[:control])

        pop!(FSD.mapping, :treatment)
        pop!(FSD.mapping, :instrument)

        obj = Poisson(FSD)

        _fit!(obj, getweights(obj))

    elseif method == "Reduced form"

        FSD                   = Microdata(MD)
        FSD.mapping[:control] = vcat(MD.mapping[:instrument], MD.mapping[:control])

        pop!(FSD.mapping, :treatment)

        obj = Poisson(FSD)

        _fit!(obj, getweights(obj))

    else

        if length(MD.mapping[:treatment]) == length(MD.mapping[:instrument])
            k   = length(MD.mapping[:instrument]) + length(MD.mapping[:control])
            obj = Mullahy(MD, Matrix{Float64}(I, k, k), "Method of moments")
        elseif method == "One-step GMM"
            k   = length(MD.mapping[:instrument]) + length(MD.mapping[:control])
            obj = Mullahy(MD, Matrix{Float64}(I, k, k), "One-step GMM")
        elseif (method == "TSLS") | (method == "2SLS")
            W   = crossprod(getmatrix(MD, :instrument, :control), getweights(MD))
            obj = Mullahy(MD, W, "Two-step GMM")
        elseif (method == "Two-step GMM") | (method == "Optimal GMM")
            W     = crossprod(getmatrix(MD, :instrument, :control), getweights(MD))
            obj   = Mullahy(MD, W, method) ; _fit!(obj, getweights(obj))
            obj.W = wmatrix(obj, getcorr(obj), getweights(obj))
        else
            throw("unknown method")
        end

        _fit!(obj, getweights(obj))
    end

    novar || _vcov!(obj, getcorr(obj), getweights(obj))

    return obj
end

#==========================================================================================#

# ESTIMATION

function _fit!(obj::Mullahy, ::UnitWeights)

    O  = haskey(obj.sample.mapping, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :treatment, :control)
    z  = getmatrix(obj, :instrument, :control)
    W  = obj.W * nobs(obj)^2
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    if isdefined(obj, :β)
        β₀ = obj.β
    else
        β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))
    end

    O && (xo = getmatrix(obj, :offset, :treatment, :control))

    function L(β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ .= muladd.(y, exp.(- μ), - 1.0)
        m  = z' * μ
        mw = W \ m

        return 0.5 * dot(m, mw)
    end

    function G!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= y .* exp.(- μ)
        xx .= x .* μ
        d   = - xx' * z
        μ  .= μ .- 1.0
        mw  = W \ (z' * μ)

        mul!(g, d, mw)
    end

    function LG!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= y .* exp.(- μ)
        xx .= x .* μ
        d   = - xx' * z
        μ  .= μ .- 1.0
        m   = z' * μ
        mw  = W \ m

        mul!(g, d, mw)

        return 0.5 * dot(m, mw)
    end

    function H!(h::Matrix, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= y .* exp.(- μ)
        xx .= x .* μ
        d   = z' * xx
        dw  = W \ d

        mul!(h, transpose(d), dw)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("minimization did not converge")
    end
end

function _fit!(obj::Mullahy, w::AbstractWeights)

    O  = haskey(obj.sample.mapping, :offset)
    y  = getvector(obj, :response)
    x  = getmatrix(obj, :treatment, :control)
    z  = getmatrix(obj, :instrument, :control)
    W  = obj.W * nobs(obj)^2
    μ  = Array{Float64}(undef, length(y))
    xx = Array{Float64}(undef, size(x)...)

    if isdefined(obj, :β)
        β₀ = obj.β
    else
        β₀ = vcat(fill(0.0, size(x, 2) - 1), log(mean(y)))
    end

    O && (xo = getmatrix(obj, :offset, :treatment, :control))

    function L(β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ .= w .* muladd.(y, exp.(- μ), - 1.0)
        m  = z' * μ
        mw = W \ m

        return 0.5 * dot(m, mw)
    end

    function G!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= y .* exp.(- μ)
        xx .= x .* μ .* w
        d   = - xx' * z
        μ  .= w .* (μ .- 1.0)
        mw  = W \ (z' * μ)

        mul!(g, d, mw)
    end

    function LG!(g::Vector, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= y .* exp.(- μ)
        xx .= x .* μ .* w
        d   = - xx' * z
        μ  .= w .* (μ .- 1.0)
        m   = z' * μ
        mw  = W \ m

        mul!(g, d, mw)

        return 0.5 * dot(m, mw)
    end

    function H!(h::Matrix, β::Vector)

        O ? mul!(μ, xo, vcat(1.0, β)) : mul!(μ, x, β)

        μ  .= w .* y .* exp.(- μ)
        xx .= x .* μ
        d   = z' * xx
        dw  = W \ d

        mul!(h, transpose(d), dw)
    end

    res = optimize(TwiceDifferentiable(L, G!, LG!, H!, β₀), β₀, Newton())

    if Optim.converged(res)
        obj.β = Optim.minimizer(res)
    else
        throw("minimization did not converge")
    end
end

#==========================================================================================#

# SCORE (MOMENT CONDITIONS)

score(obj::Mullahy) = Diagonal(residuals(obj)) * getmatrix(obj, :instrument, :control)

# EXPECTED JACOBIAN OF SCORE × NUMBER OF OBSERVATIONS

function jacobian(obj::Mullahy, ::UnitWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if haskey(obj.sample.mapping, :offset)
        x = getmatrix(obj, :offset, :treatment, :control)
        v = xo * vcat(1.0, obj.β)
    else
        x = getmatrix(obj, :treatment, :control)
        v = x * obj.β
    end

    v .= y .* exp.(- v)

    return - z' * (x .* v)
end

function jacobian(obj::Mullahy, w::AbstractWeights)

    y = getvector(obj, :response)
    x = getmatrix(obj, :treatment, :control)
    z = getmatrix(obj, :instrument, :control)

    if haskey(obj.sample.mapping, :offset)
        x = getmatrix(obj, :offset, :treatment, :control)
        v = xo * vcat(1.0, obj.β)
    else
        x = getmatrix(obj, :treatment, :control)
        v = x * obj.β
    end

    v .= w .* y .* exp.(- v)

    return - z' * (x .* v)
end

#==========================================================================================#

# LINEAR PREDICTOR

function predict(obj::Mullahy, MD::Microdata)
    if getnames(obj, :treatment, :control) != getnames(MD, :treatment, :control)
        throw("missing variables")
    end
    if haskey(obj.sample.mapping, :offset)
        return getmatrix(MD, :offset, :treatment, :control) * vcat(1.0, obj.β)
    else
        return getmatrix(MD, :treatment, :control) * obj.β
    end
end

# FITTED VALUES

fitted(obj::Mullahy, MD::Microdata) = exp.(predict(obj, MD))

# DERIVATIVE OF FITTED VALUES

function jacobexp(obj::Mullahy)

    if haskey(obj.sample.mapping, :offset)
        v = getmatrix(obj, :offset, :treatment, :control) * vcat(1.0, obj.β)
    else
        v = getmatrix(obj, :treatment, :control) * obj.β
    end

    v .= exp.(v)

    return Diagonal(v) * getmatrix(obj, :treatment, :control)
end

# RESIDUALS

function residuals(obj::Mullahy, MD::Microdata)
    r  = fitted(obj, MD)
    r .= response(obj) ./ r .- 1.0
    return r
end

#==========================================================================================#

# UTILITIES

coefnames(obj::Mullahy) = getnames(obj, :treatment, :control)
mtitle(obj::Mullahy)    = "IV Poisson with multiplicative errors"
