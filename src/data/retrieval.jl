#==========================================================================================#

# GET VARIABLE SETS

function getvector(MD::Microdata, x::Symbol)
    haskey(MD.map, x)       || throw(string(x) * " not found")
    iszero(MD.map[x])       && throw(string(x) * " not found")
    (length(MD.map[x]) > 1) && throw(string(x) * " is not a vector")
    return view(MD.mat, :, MD.map[x]...)
end

getvector(MM::ParModel, x::Symbol)      = getvector(MM.sample, x)
getvector(MM::TwoStageModel, x::Symbol) = getvector(MM.second_stage.sample, x)

function getmatrix(MD::Microdata, args...)

    n = length(args)

    for i in args
        haskey(MD.map, i) || throw(string(i) * " not found")
        iszero(MD.map[i]) && throw(string(i) * " not found")
    end

    x = Vector{Int64}()

    for i in args
        x = vcat(x, MD.map[i])
    end

    return view(MD.mat, :, x)
end

getmatrix(MM::ParModel, args...)        = getmatrix(MM.sample, args...)
getmatrix(MM::TwoStageModel, x::Symbol) = getmatrix(MM.second_stage.sample, x)

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

getnames(MM::ParModel, args...)      = getnames(MM.sample, args...)
getnames(MM::TwoStageModel, args...) = getnames(MM.second_stage.sample, args...)

#==========================================================================================#

# GET CORRELATION STRUCTURE

getcorr(obj::Microdata)     = obj.corr
getcorr(obj::ParModel)      = getcorr(obj.sample)
getcorr(obj::TwoStageModel) = getcorr(second_stage(obj).sample)

#==========================================================================================#

# GET INDICATOR OF MISSING DATA

getnonmissing(obj::Microdata)     = obj.nonmissing
getnonmissing(obj::ParModel)      = getnonmissing(obj.sample)
getnonmissing(obj::TwoStageModel) = getnonmissing(second_stage(obj).sample)

#==========================================================================================#

# GET WEIGHT

getweights(MD::Microdata)     = MD.weights
getweights(MM::Micromodel)    = getweights(MM.sample)
getweights(MM::TwoStageModel) = getweights(MM.second_stage.sample)
