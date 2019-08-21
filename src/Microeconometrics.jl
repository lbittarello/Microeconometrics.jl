__precompile__(true)

module Microeconometrics

using DataFrames
using Dates
using Format: format
using LinearAlgebra
using Optim
using SparseArrays
using SpecialFunctions: lgamma
using StatsFuns
using StatsModels
    const SM = StatsModels
using SuiteSparse

using StatsBase: AbstractWeights, CoefTable, RegressionModel
using StatsBase: AnalyticWeights, FrequencyWeights, ProbabilityWeights, Weights
using StatsBase: mean, sum, pweights

import Base:        copy, isequal, show, sum
import Statistics:  mean
import StatsBase:   fit, coef, coefnames, coeftable, confint, stderror, vcov
import StatsBase:   deviance, loglikelihood, nulldeviance, nullloglikelihood
import StatsBase:   adjr2, dof, dof_residual, nobs, r2
import StatsBase:   fitted, predict, residuals, response
import StatsModels: terms, termvars, schema, apply_schema, has_schema, modelcols

include("./inference/corr.jl")
include("./data/types.jl")
include("./data/weights.jl")
include("./general/types.jl")
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
include("./estimation/gompit.jl")
include("./estimation/poisson.jl")
include("./estimation/ivpoisson.jl")
include("./estimation/mullahy.jl")
include("./estimation/ipw.jl")
include("./estimation/abadie.jl")
include("./estimation/frölichmelly.jl")
include("./estimation/tan.jl")

export
    @micromodel, Microdata,
    Homoscedastic, Heteroscedastic, Clustered, CrossCorrelated,
    OLS, Logit, Probit, Cloglog, Gompit, Poisson,
    IV, IVPoisson, Mullahy,
    IPW, Abadie, FrölichMelly, Tan

export
    fit, coef, confint, stderror, vcov,
    deviance, loglikelihood, nulldeviance, nullloglikelihood,
    adjr2, dof, dof_residual, nobs, r2,
    fitted, predict, residuals, response,
    coefnames, coeftable, etable,
    hausman_1s, hausman_2s, pval, tstat

end # module
