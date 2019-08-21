#==========================================================================================#

# GET VARIABLE SETS

function getvector(obj::Microdata, x::Symbol)
    (length(obj.mapping[x]) == 1) || throw(ArgumentError(string(x) * " is not a vector"))
    return view(obj.matrix, :, obj.mapping[x]...)
end

getvector(obj::ParModel, x::Symbol)      = getvector(obj.sample, x)
getvector(obj::TwoStageModel, x::Symbol) = getvector(obj.second_stage.sample, x)

function getmatrix(obj::Microdata, args...)
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    return view(obj.matrix, :, x)
end

getmatrix(obj::ParModel, args...)      = getmatrix(obj.sample, args...)
getmatrix(obj::TwoStageModel, args...) = getmatrix(obj.second_stage.sample, args...)

#==========================================================================================#

# GET VARIABLE NAMES

function getnames(obj::Microdata, args...)
    n = coefnames(eterms(obj.model))
    x = mapreduce(i -> obj.mapping[i], vcat, args, init = Vector{Int64}())
    return n[x]
end

getnames(obj::ParModel, args...)      = getnames(obj.sample, args...)
getnames(obj::TwoStageModel, args...) = getnames(obj.second_stage.sample, args...)

#==========================================================================================#

# GET CORRELATION STRUCTURE

getcorr(obj::Microdata)     = obj.corr
getcorr(obj::ParModel)      = getcorr(obj.sample)
getcorr(obj::TwoStageModel) = getcorr(obj.second_stage.sample)

#==========================================================================================#

# GET INDICATOR OF MISSING DATA

getnonmissing(obj::Microdata)     = obj.nonmissing
getnonmissing(obj::ParModel)      = getnonmissing(obj.sample)
getnonmissing(obj::TwoStageModel) = getnonmissing(obj.second_stage.sample)

#==========================================================================================#

# GET WEIGHT

getweights(obj::Microdata)     = obj.weights
getweights(obj::Micromodel)    = getweights(obj.sample)
getweights(obj::TwoStageModel) = getweights(obj.second_stage.sample)
