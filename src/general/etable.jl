#==========================================================================================#

# TABLE OF ESTIMATES

function etable(args...;
        digits::Int = 4,
        aux::Union{Function, Nothing} = nothing,
        stars::Matrix{Any} = [0.1 "*"; 0.05 "**"; 0.01 "***"],
        titles::Vector{String} = [""]
    )

    N                  = length(args)
    cutpoints, symbols = _parsestars(stars)

    (titles == [""]) && (titles = String["(" * string(i) * ")" for i = 1:N])

    μ = Vector{Vector{String}}(undef, N)
    β = Vector{Vector{String}}(undef, N)
    σ = Vector{Vector{String}}(undef, N)

    for (i, ai) in enumerate(args)

        μ[i] = coefnames(ai)
        β[i] = format.("{:.$(digits)f}", coef(ai))

        if typeof(aux) != Nothing
            σ[i] = "(" .* format.("{:.$(digits)f}", aux(ai)) .* ")"
        end

        if cutpoints != []
            for (j, pj) in enumerate(pval(ai))
                s        = searchsortedlast(cutpoints, pj)
                β[i][j] *= symbols[s]
            end
        end
    end

    inter = union(μ...)

    if typeof(aux) == Nothing
        output = vcat("", inter)
    else
        output          = fill("", 2 * length(inter) + 1)
        output[2:2:end] = inter
    end

    w1 = similar(output)
    w2 = similar(output[2:end])

    _alignatchar!(output)

    for (i, (ni, ti, βi)) in enumerate(zip(μ, titles, β))
        idx = findall((in)(ni), inter)
        (typeof(aux) == Nothing) || (idx .= 2 .* idx .- 1)
        w2      .= ""
        w2[idx] .= βi
        (typeof(aux) == Nothing) || (w2[(idx .+ 1)] .= σ[i])
        _alignatchar!(w2, '.')
        w1[1]      = ti
        w1[2:end] .= w2
        _alignatchar!(w1)
        output .*= "    " .* w1
    end

    for i in output
        println(i)
    end
end

#==========================================================================================#

# PARSING STAR SCHEME

function _parsestars(x::Matrix{Any})
    if size(x, 1) != 0
        s         = sort(x, dims = 1, rev = true)
        cutpoints = vcat(0.0, convert(Vector{Float64}, s[:, 1]), 1.0)
        symbols   = vcat(convert(Vector{String}, s[:, 2]), "")
    else
        cutpoints = Vector{Float64}()
        symbols   = Vector{String}()
    end
    return cutpoints, symbols
end

#==========================================================================================#

# ALIGN STRING VECTOR AT FIRST CHARACTER

function _alignatchar!(x::Vector{String})
    len    = length.(x)
    maxlen = maximum(len)
    len    = maxlen .- len
    add    = " " .^ len
    x     .= x .* add
end

# ALIGN STRING VECTOR AT FIRST OCCURRENCE OF CHARACTER

function _alignatchar!(x::Vector{String}, y::Char)

    pos = fill(0, length(x))

    for (xi, xx) in enumerate(x)
        (xx == "") || (pos[xi] = findlast(z -> z == y, xx))
    end
    
    maxpos = maximum(pos)
    add    = " " .^ (maxpos .- pos)
    x     .= add .* x
end
