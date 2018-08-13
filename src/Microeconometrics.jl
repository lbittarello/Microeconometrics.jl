__precompile__(true)

module Microeconometrics

using DataFrames
using Formatting
using LinearAlgebra
using Optim
using SparseArrays
using SpecialFunctions: lgamma
using StatsBase
using StatsFuns
using StatsModels

import Base:       copy, deepcopy, isequal, sum
import Statistics: mean
import StatsBase:  fit, coef, coefnames, coeftable, confint, stderror, vcov
import StatsBase:  deviance, loglikelihood, nulldeviance, nullloglikelihood
import StatsBase:  adjr2, aic, aicc, bic, dof, dof_residual, nobs, r2
import StatsBase:  fitted, predict, residuals, response

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
    Microdata, Homoscedastic, Heteroscedastic, Clustered, CrossCorrelated,
    OLS, Logit, Probit, Cloglog, Poisson,
    IV, IVPoisson, Mullahy,
    IPW, Abadie, FrölichMelly, Tan

export
    fit, coef, confint, stderror, vcov,
    deviance, loglikelihood, nulldeviance, nullloglikelihood,
    adjr2, aic, aicc, bic, dof, dof_residual, nobs, r2,
    fitted, predict, residuals, response,
    coefnames, coeftable, etable,
    hausman_1s, hausman_2s, pval, tstat

end # module
