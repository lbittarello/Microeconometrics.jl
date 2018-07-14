__precompile__(true)

module Microeconometrics

using StatsBase
using StatsFuns
using DataFrames
using StatsModels
using Optim
using Formatting

import Base:        copy, deepcopy, size
import StatsBase:   RegressionModel, CoefTable, coefnames, coeftable
import StatsBase:   fit, coef, confint, stderr, vcov
import StatsBase:   loglikelihood, nullloglikelihood, deviance, nulldeviance
import StatsBase:   aic, aicc, bic, dof, dof_residual, nobs, r2, r², adjr2, adjr²
import StatsBase:   model_response, predict, fitted, residuals
import StatsModels: Terms, evalcontrasts

const PWeights = ProbabilityWeights{Float64, Float64, Vector{Float64}}

include("./inference/corr.jl")
include("./data/types.jl")
include("./data/weights.jl")
include("./general/types.jl")
include("./data/retrieval.jl")
include("./data/utils.jl")
include("./general/etable.jl")
include("./general/statsmodels.jl")
include("./general/utils.jl")
include("./inference/adjfactor.jl")
include("./inference/hausman.jl")
include("./inference/utils.jl")
include("./inference/vcov.jl")
include("./estimation/ols.jl")
include("./estimation/iv.jl")
include("./estimation/logit.jl")
include("./estimation/probit.jl")
include("./estimation/cloglog.jl")
include("./estimation/poisson.jl")
include("./estimation/ipw.jl")
include("./estimation/abadie.jl")
include("./estimation/frölichmelly.jl")
include("./estimation/tan.jl")

export

    CorrStructure,
        Homoscedastic,
        Heteroscedastic,
        Clustered,
        CrossCorrelated,
    Microdata, Microdata!,
    OLS, IV, Logit, Probit, Cloglog, Poisson,
    IPW, Abadie, FrölichMelly, Tan,
    fit, first_stage, second_stage,
    coef, confint, pval, stderr, tstat, vcov, hausman_1s, hausman_2s,
    loglikelihood, nullloglikelihood, deviance, nulldeviance,
    aic, aicc, bic, dof, dof_residual, nobs, r2, r², adjr2, adjr²,
    model_response, predict, fitted, residuals,
    coefnames, coeftable, etable

end # module
