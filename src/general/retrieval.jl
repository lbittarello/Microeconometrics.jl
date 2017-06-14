#==========================================================================================#

# GET VARIABLE SETS

function getvector(MD::Microdata, x::Symbol)
    haskey(MD.map, x)       || throw(string(x) * " not found")
    iszero(MD.map[x])       && throw(string(x) * " not found")
    (length(MD.map[x]) > 1) && throw(string(x) * "is not a vector")
    return view(MD.mat, :, MD.map[x]...)
end

getvector(MM::Micromodel, x::Symbol) = getvector(MM.sample, x)

function getmatrix(MD::Microdata, args...)

    n = length(args)

    for i in args
        haskey(MD.map, i) || throw(string(x) * " not found")
        iszero(MD.map[i]) && throw(string(i) * " not found")
    end

    x = Vector{Int64}()

    for i in args
        x = vcat(x, MD.map[i])
    end

    return view(MD.mat, :, x)
end

getmatrix(MM::Micromodel, args...) = getmatrix(MM.sample, args...)

#==========================================================================================#

# GET VARIABLE NAMES

function getnames(MD::Microdata, args...)

    n = length(args)

    for i in args
        haskey(MD.map, i) || throw(string(x) * " not found")
        iszero(MD.map[i]) && throw(string(i) * " not found")
    end

    x = Vector{Int64}()

    for i in args
        x = vcat(x, MD.map[i])
    end

    return MD.names[x]
end

getnames(MM::Micromodel, args...) = getnames(MM.sample, args...)

#==========================================================================================#

# GET CORRELATION STRUCTURE

getcorr(obj::ParModel)      = obj.sample.corr
getcorr(obj::TwoStageModel) = second_stage(obj).sample.corr

#==========================================================================================#

# CHECK WEIGHT

checkweight(MD::Microdata)     = (haskey(MD.map, :weight) && !iszero(MD.map[:weight]))
checkweight(MM::Micromodel)    = checkweight(MM.sample)
checkweight(MM::TwoStageModel) = checkweight(second_stage(MM))
