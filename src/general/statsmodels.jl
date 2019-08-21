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
nobs(obj::ParModel)  = nobs(switch_stage(obj).sample)

dof(obj::ParModel)          = length(coef(obj))
dof_residual(obj::ParModel) = nobs(obj) - dof(obj)

loglikelihood(obj::MLE)     = _loglikelihood(obj, getweights(obj))
nullloglikelihood(obj::MLE) = _nullloglikelihood(obj::MLE, getweights(obj))
deviance(obj::MLE)          = _deviance(obj::MLE, getweights(obj))
nulldeviance(obj::MLE)      = _nulldeviance(obj::MLE, getweights(obj))

#==========================================================================================#

# PREDICTION

predict(obj::AnyModel)   = predict(obj, obj.sample)
fitted(obj::ParModel)    = fitted(obj, obj.sample)
residuals(obj::AnyModel) = residuals(obj, obj.sample)
response(obj::AnyModel)  = getvector(obj, :response)

function residuals(obj::AnyModel, MD::Microdata)
    y  = response(obj)
    r  = fitted(obj, MD)
    r .= y .- r
    return r
end

#==========================================================================================#

# COEFFICIENT LABELS

coefnames(obj::ParEstimate) = obj.names

# OUTPUT

function coeftable(obj::ParObject; level::Float64 = 0.95, digits::Int = 4)

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

function show(io::IO, obj::ParObject)
    if isdefined(switch_stage(obj), :V)
        println(io, "$(mtitle(obj))\n\n", coeftable(obj))
    else
        println(io, "$(mtitle(obj))\n\n", coef(obj))
    end
end
