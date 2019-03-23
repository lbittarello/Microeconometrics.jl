#==========================================================================================#

# INTERFACE

function fit(M::Type{<: ParModel}, MD::Microdata; novar::Bool = false, kwargs...)

    obj = M(MD; kwargs...)

    _fit!(obj, MD.weights)
    novar || _vcov!(obj, getcorr(obj), MD.weights)

    return obj
end

#==========================================================================================#

# FIRST AND SECOND STAGES

first_stage(obj::TwoStageModel)  = obj.first_stage
second_stage(obj::TwoStageModel) = obj.second_stage

#==========================================================================================#

# ESTIMATES

coef(obj::Par1S)         = obj.β
coef(obj::TwoStageModel) = coef(second_stage(obj))
vcov(obj::Par1S)         = obj.V
vcov(obj::TwoStageModel) = vcov(second_stage(obj))
stderror(obj::Par2S)     = sqrt.(diag(vcov(obj)))
tstat(obj::Par2S)        = coef(obj) ./ stderror(obj)
pval(obj::Par2S)         = 2.0 * normccdf.(abs.(tstat(obj)))

function confint(obj::Par2S, level::Real = 0.95)
    return coef(obj) .+ norminvcdf((1.0 - level) / 2.0) * stderror(obj) .* [1.0 -1.0]
end

#==========================================================================================#

# SUMMARY STATISTICS

nobs(obj::Microdata)     = sum(obj.weights)
nobs(obj::Micromodel)    = nobs(obj.sample)
nobs(obj::TwoStageModel) = nobs(second_stage(obj))

dof(obj::ParModel)        = length(coef(obj))
dof(obj::TwoStageModel)   = dof(second_stage(obj))
dof_residual(obj::ParM2S) = nobs(obj) - dof(obj)

loglikelihood(obj::MLE)     = _loglikelihood(obj, getweights(obj))
nullloglikelihood(obj::MLE) = _nullloglikelihood(obj::MLE, getweights(obj))
deviance(obj::MLE)          = _deviance(obj::MLE, getweights(obj))
nulldeviance(obj::MLE)      = _nulldeviance(obj::MLE, getweights(obj))

#==========================================================================================#

# PREDICTION

predict(obj::ParM2S)      = predict(obj, obj.sample)
fitted(obj::ParM2S)       = fitted(obj, obj.sample)
residuals(obj::ParM2S)    = residuals(obj, obj.sample)
response(obj::Micromodel) = getvector(obj, :response)

function residuals(obj::ParM2S, MD::Microdata)
    r  = fitted(obj, MD)
    r .= response(obj) .- r
    return r
end

#==========================================================================================#

# COEFFICIENT LABELS

coefnames(obj::ParEstimate) = obj.names

# OUTPUT

function coeftable(obj::Par2S; level::Float64 = 0.95, digits::Int = 4)

    table = frmtr(hcat(coef(obj), stderror(obj), tstat(obj), pval(obj)), digits)
    label = [" Estimate", " St. Err.", "  t-stat.", "  p-value"]

    if level > 0.0
        lprint = format("{:.0d}", 100 * 0.95)
        table  = hcat(table, frmtr(confint(obj, level), digits))
        label  = vcat(label, ["     C.I.", "($(lprint)%)  "])
    end

    CT = CoefTable(table, label, coefnames(obj))

    return CT
end

function show(io::IO, obj::Par1S)
    if isdefined(obj, :V)
        println(io, "$(mtitle(obj))\n\n", coeftable(obj))
    else
        println(io, "$(mtitle(obj))\n\n", coef(obj))
    end
end

function show(io::IO, obj::TwoStageModel)
    if isdefined(obj.second_stage, :V)
        println(io, "$(mtitle(obj))\n\n", coeftable(obj))
    else
        println(io, "$(mtitle(obj))\n\n", coef(obj))
    end
end
