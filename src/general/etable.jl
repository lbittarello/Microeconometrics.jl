#==========================================================================================#

# TABLE OF ESTIMATES

function etable(args...;
        digits::Int = 4,
        f::Function = stderr,
        compact::Bool = false,
        stars::Matrix{Any} = [0.1 "*"; 0.05 "**"; 0.01 "***"],
        titles::Vector{String} = []
    )

    N                  = length(args)
    fspec              = FormatSpec("0.$(digits)f")
    cutpoints, symbols = _parsestars(stars)

    (titles == []) && (titles = String["(" * string(i) * ")" for i = 1:N])

    β     = Vector{Vector{String}}(N)
    σ     = Vector{Vector{String}}(N)
    names = Vector{Vector{String}}(N)

    for (i, ai) in enumerate(args)

        β[i]     = fmt.(fspec, coef(ai))
        names[i] = coefnames(ai)

        compact || (σ[i] = "(" .* fmt.(fspec, f(ai)) .* ")")

        if cutpoints != []
            for (j, pj) in enumerate(pval(ai))
                s        = searchsortedlast(cutpoints, pj)
                β[i][j] *= symbols[s]
            end
        end
    end

    inter = union(names...)

    if compact
        output = vcat("", inter)
    else
        output          = fill("", 2 * length(inter) + 1)
        output[2:2:end] = inter
    end

    w1 = similar(output)
    w2 = similar(output[2:end])

    _fixwidth!(output)

    for (i, (ni, ti, βi)) in enumerate(zip(names, titles, β))
        idx = findin(inter, ni)
        compact || (idx .= 2 * idx - 1)
        w2      .= ""
        w2[idx] .= βi
        compact || (w2[idx + 1] .= σ[i])
        _alignatchar!(w2, '.')
        w1[1]     .= ti
        w1[2:end] .= w2
        _fixwidth!(w1)
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
        s         = sortrows(x)
        cutpoints = vcat(0.0, convert(Vector{Float64}, s[:, 1]), 1.0)
        symbols   = vcat(convert(Vector{String}, s[:, 2]), "")
    else
        cutpoints = Vector{Float64}()
        symbols   = Vector{String}()
    end
    return cutpoints, symbols
end

#==========================================================================================#

# HARMONIZE WIDTH OF STRING VECTOR

function _fixwidth!(x::Vector{String}; rev::Bool = true)
    len    = length.(x)
    maxlen = maximum(len)
    len    = maxlen .- len
    add    = " " .^ len
    rev ? (x .= x .* add) : (x .= add .* x)
end

#==========================================================================================#

# ALIGN STRING VECTOR AT FIRST OCCURRENCE OF CHARACTER

function _alignatchar!(x::Vector{String}, y::Char)
    pos    = search.(x, y)
    maxpos = maximum(pos)
    add    = " " .^ (maxpos .- pos)
    x     .= add .* x
end
