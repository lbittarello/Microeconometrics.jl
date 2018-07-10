#==========================================================================================#

# INTERFACE

function fit(::Type{M}, MD::Microdata; novar::Bool = false, kwargs...) where {M <: ParModel}

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

coef(obj::Union{ParModel, ParObject})     = obj.β
coef(obj::TwoStageModel)                  = coef(second_stage(obj))
vcov(obj::Union{ParModel, ParObject})     = obj.V
vcov(obj::TwoStageModel)                  = vcov(second_stage(obj))
stderr(obj::Union{Micromodel, ParObject}) = sqrt.(diag(vcov(obj)))
tstat(obj::Union{Micromodel, ParObject})  = coef(obj) ./ stderr(obj)
pval(obj::Union{Micromodel, ParObject})   = 2.0 * normccdf.(abs.(tstat(obj)))

function confint(obj::Union{Micromodel, ParObject}, level::Float64 = 0.95)
    return coef(obj) .+ norminvcdf((1.0 - level) / 2.0) * stderr(obj) .* [1.0 -1.0]
end

#==========================================================================================#

# SUMMARY STATISTICS

nobs(obj::Microdata)     = sum(obj.weights)
nobs(obj::Micromodel)    = nobs(obj.sample)
nobs(obj::TwoStageModel) = nobs(second_stage(obj))

dof(obj::ParModel)             = length(coef(obj))
dof(obj::TwoStageModel)        = dof(second_stage(obj))
dof_residual(obj::ParOr2Stage) = nobs(obj) - dof(obj)

loglikelihood(obj::MLE)     = _loglikelihood(obj, getweights(obj))
nullloglikelihood(obj::MLE) = _nullloglikelihood(obj::MLE, getweights(obj))
deviance(obj::MLE)          = _deviance(obj::MLE, getweights(obj))
nulldeviance(obj::MLE)      = _nulldeviance(obj::MLE, getweights(obj))

#==========================================================================================#

# PREDICTION

predict(obj::ParModel) = predict(obj, obj.sample)
fitted(obj::ParModel)  = fitted(obj, obj.sample)

model_response(obj::Micromodel) = getvector(obj, :response)

function residuals(obj::Micromodel)
    r  = fitted(obj)
    r .= model_response(obj) .- r
    return r
end

#==========================================================================================#

# COEFFICIENT NAMES

coefnames(obj::ParObject) = obj.names

# OUTPUT

function coeftable(
        obj::Union{ParOr2Stage, ParObject};
        level::Float64 = 0.95,
        digits::Int = 4
    )

    table = round.(hcat(coef(obj), stderr(obj), tstat(obj), pval(obj)), digits)
    label = [" Estimate", " St. Err.", "  t-stat.", "  p-value"]

    if level > 0.0
        lprint = fmt(FormatSpec("0.0d"), 100 * level)
        table  = hcat(table, round.(confint(obj, level), digits))
        label  = vcat(label, ["     C.I.", "($(lprint)%)  "])
    end

    CT = CoefTable(table, label, coefnames(obj), 4)

    return CT
end
