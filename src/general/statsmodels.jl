#==========================================================================================#

# INTERFACE

function fit(M::Type{<: OneStageModel}, MD::Microdata; novar::Bool = false, kwargs...)

    obj = M(MD; kwargs...)

    _fit!(obj, MD.weights)
    novar || _vcov!(obj, getcorr(obj), MD.weights)

    return obj
end

#==========================================================================================#

# ESTIMATES

coef(obj::ParObject)     = switch_stage(obj).β
vcov(obj::ParObject)     = switch_stage(obj).V
stderror(obj::ParObject) = sqrt.(diag(vcov(obj)))
tstat(obj::ParObject)    = coef(obj) ./ stderror(obj)
pval(obj::ParObject)     = 2.0 * normccdf.(abs.(tstat(obj)))

function confint(obj::ParObject, level::Real = 0.95)
    return coef(obj) .+ norminvcdf((1.0 - level) / 2.0) * stderror(obj) .* [1.0 -1.0]
end

#==========================================================================================#

# SUMMARY STATISTICS

nobs(obj::Microdata) = sum(obj.weights)
nobs(obj::AnyModel)  = nobs(switch_stage(obj).sample)

dof(obj::ParModel)          = length(coef(obj))
dof_residual(obj::ParModel) = nobs(obj) - dof(obj)

loglikelihood(obj::MLE)     = _loglikelihood(obj, getweights(obj))
nullloglikelihood(obj::MLE) = _nullloglikelihood(obj::MLE, getweights(obj))
deviance(obj::MLE)          = _deviance(obj::MLE, getweights(obj))
nulldeviance(obj::MLE)      = _nulldeviance(obj::MLE, getweights(obj))

#==========================================================================================#

# PREDICTION

linear_predictor(obj::AnyModel) = linear_predictor(obj, switch_stage(obj).sample)
predict(obj::AnyModel)          = predict(obj, switch_stage(obj).sample)
fitted(obj::AnyModel)           = predict(obj)
residuals(obj::AnyModel)        = residuals(obj, switch_stage(obj).sample)
response(obj::AnyModel)         = getvector(obj, :response)

function residuals(obj::AnyModel, MD::Microdata)
    y  = response(obj)
    r  = predict(obj, MD)
    r .= y .- r
    return r
end

#==========================================================================================#

# COEFFICIENT LABELS

coefnames(obj::ParEstimate) = obj.names

# OUTPUT

function coeftable(obj::ParObject; level::Float64 = 0.95, digits::Int = 4)

    table = formatter(hcat(coef(obj), stderror(obj), tstat(obj), pval(obj)), digits)
    label = [" Estimate", " St. Err.", "  t-stat.", "  p-value"]

    if level > 0.0
        table  = hcat(table, formatter(confint(obj, level), digits))
        label  = vcat(label, ["     C.I.", "($(format("{:.0d}", 100 * level))%)  "])
    end

    return CoefTable(table, label, coefnames(obj))
end

function Base.show(io::IO, obj::ParObject)
    isdefined(switch_stage(obj), :V) ? println(io, coeftable(obj)) : println(io, coef(obj))
end
