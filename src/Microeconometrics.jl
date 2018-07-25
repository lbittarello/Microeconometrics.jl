__precompile__(true)

module Microeconometrics

using LinearAlgebra
using SparseArrays
using Statistics
using StatsBase
using StatsFuns
using DataFrames
using StatsModels
using Optim
using Formatting

import Base:        copy, deepcopy, size
import StatsBase:   RegressionModel, CoefTable, coefnames, coeftable, islinear
import StatsBase:   fit, coef, confint, informationmatrix, score, stderror, vcov
import StatsBase:   deviance, loglikelihood, nulldeviance, nullloglikelihood
import StatsBase:   adjr2, adjr², aic, aicc, bic, dof, dof_residual, mss, nobs, r2, r², rss
import StatsBase:   fitted, meanresponse, predict, residuals, response
import StatsModels: Terms

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
include("./estimation/ivpoisson.jl")
include("./estimation/mullahy.jl")
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
    Microdata,
    OLS, IV, IVPoisson, Mullahy,
    Logit, Probit, Cloglog, Poisson,
    IPW, Abadie, FrölichMelly, Tan,
    fit, coef, confint, hausman_1s, hausman_2s, informationmatrix, score, stderror, vcov,
    deviance, loglikelihood, nulldeviance, nullloglikelihood,
    adjr2, adjr², aic, aicc, bic, dof, dof_residual, mss, nobs, r2, r², rss,
    fitted, meanresponse, predict, residuals, response,
    coefnames, coeftable, etable, islinear

end # module
