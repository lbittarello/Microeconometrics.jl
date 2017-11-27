#==========================================================================================#

# INTERFACE

function fit(::Type{M}, MD::Microdata; novar::Bool = false, kwargs...) where {M <: ParModel}

    obj = M(MD; kwargs...)

    if checkweight(MD)
        w = getvector(MD, :weight)
        _fit!(obj, w)
        novar || _vcov!(obj, w)
    else
        _fit!(obj)
        novar || _vcov!(obj)
    end

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

nobs(obj::Microdata)           = size(obj.mat, 1)
nobs(obj::Micromodel)          = nobs(obj.sample)
nobs(obj::TwoStageModel)       = nobs(second_stage(obj))
dof(obj::ParModel)             = length(coef(obj))
dof(obj::TwoStageModel)        = dof(second_stage(obj))
dof_residual(obj::ParOr2Stage) = nobs(obj) - dof(obj)

adjr2(obj::ParOr2Stage) = 1.0 - ((nobs(obj) - 1) / dof_residual(obj)) * (1.0 - r2(obj))
r2(obj::Micromodel)     = (checkweight(obj) ? _r2(obj, getvector(obj, :weight)) : _r2(obj))

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
    if checkweight(obj)
        _loglikelihood(obj, getvector(obj, :weight))
    else
        _loglikelihood(obj)
    end
end

function nullloglikelihood(obj::MLE)
    if checkweight(obj)
        _nullloglikelihood(obj, getvector(obj, :weight))
    else
        _nullloglikelihood(obj)
    end
end

function deviance(obj::MLE)
    if checkweight(obj)
        _deviance(obj, getvector(obj, :weight))
    else
        _deviance(obj)
    end
end

function nulldeviance(obj::MLE)
    if checkweight(obj)
        _nulldeviance(obj, getvector(obj, :weight))
    else
        _nulldeviance(obj)
    end
end

#==========================================================================================#

# PREDICTION

model_response(obj::Micromodel)    = getvector(obj, :response)
model_response(obj::TwoStageModel) = getvector(second_stage(obj), :response)

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
        digits::Int = 4,
        verbose::Bool = true
    )

    table = round.(hcat(coef(obj), stderr(obj), tstat(obj), pval(obj)), digits)
    label = [" Estimate", " St. Err.", "  t-stat.", "  p-value"]

    if level > 0.0
        lprint = fmt(FormatSpec("0.0d"), 100 * level)
        table  = hcat(table, round.(confint(obj, level), digits))
        label  = vcat(label, ["     C.I.", "($(lprint)%)  "])
    end

    CT = CoefTable(table, label, coefnames(obj), 4)

    verbose && println(CT)

    return CT
end
