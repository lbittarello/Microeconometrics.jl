#==========================================================================================#

# INTERFACE

function fit{M <: ParModel}(::Type{M}, MD::Microdata; novar::Bool = false)

    output        = M()
    output.sample = MD

    if checkweight(MD)
        w = getvector(MD, :weight)
        output.β = _fit(output, w)
        novar || (output.V = _vcov(output, MD.corr, w))
    else
        output.β = _fit(output)
        novar || (output.V = _vcov(output, MD.corr))
    end

    return output
end

#==========================================================================================#

# FIRST-STAGE AND REDUCED FORM

first_stage(obj::TwoStageModel)  = obj.first_stage
second_stage(obj::TwoStageModel) = obj.second_stage
reduced_form(obj::ParModel)      = reduced_form(typeof(obj), obj.sample)

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

nobs(obj::Micromodel)            = size(obj.sample.mat, 1)
nobs(obj::TwoStageModel)         = nobs(second_stage(obj))
dof(obj::ParModel)               = length(coef(obj))
dof(obj::TwoStageModel)          = dof(second_stage(obj))
dof_residual(obj::ParOrTwoStage) = nobs(obj) - dof(obj)

r2(obj::Micromodel) = (checkweight(obj) ? _r2(obj, getvector(obj, :weight)) : _r2(obj))
adjr2(obj::ParOrTwoStage) = 1.0 - ((nobs(obj) - 1) / dof_residual(obj)) * (1.0 - r2(obj))

function _r2(obj::Micromodel)
    y   = model_response(obj)
    fit = fitted(obj)
    rss = sum(abs2, y .- fit)
    tss = sum(abs2, y .- mean(y))
    return 1.0 - rss / tss
end

function _r2(obj::Micromodel, w::AbstractVector)
    y   = model_response(obj)
    μ   = sum(w .* y) / sum(w)
    fit = fitted(obj)
    rss = sum(w .* abs2.(y .- fit))
    tss = sum(w .* abs2.(y .- μ))
    return 1.0 - rss / tss
end

function loglikelihood(obj::MLE)
    checkweight(obj) ? _loglikelihood(obj, getvector(obj, :weight)) : _loglikelihood(obj)
end

function nullloglikelihood(obj::MLE)
    if checkweight(obj)
        _nullloglikelihood(obj, getvector(obj, :weight))
    else
        _nullloglikelihood(obj)
    end
end

function deviance(obj::MLE)
    checkweight(obj) ? _deviance(obj, getvector(obj, :weight)) : _deviance(obj)
end

function nulldeviance(obj::MLE)
    checkweight(obj) ? _nulldeviance(obj, getvector(obj, :weight)) : _nulldeviance(obj)
end

#==========================================================================================#

# PREDICTION AND RELATED

model_response(obj::Micromodel)    = getvector(obj, :response)
model_response(obj::TwoStageModel) = getvector(second_stage(obj), :response)

function residuals(obj::Micromodel)
    r  = fitted(obj)
    r .= model_response(obj) .- r
    return r
end

#==========================================================================================#

# OUTPUT

coefnames(obj::ParObject) = obj.names

function coeftable(
        obj::Union{ParModel, ParObject, TwoStageModel};
        level::Float64 = 0.05,
        digits::Int = 4,
        verbose::Bool = true
    )

    table = round.(hcat(coef(obj), stderr(obj), tstat(obj), pval(obj)), digits)
    label = [" Estimate", " St. Err.", "  t-stat.", "  p-value"]

    if level > 0.0
        lprint = fmt.(FormatSpec("0.2f"), level)
        table  = hcat(table, confint(obj, level))
        label  = vcat(label, ["     C.I.", "($(lprint)%)  "])
    end

    CT = CoefTable(table, label, coefnames(obj), 4)

    verbose && println(CT)

    return CT
end
