__precompile__(true)

module Microeconometrics

using StatsFuns
using DataFrames
using StatsModels
using Optim
using Formatting

import Base:        copy
import StatsBase:   RegressionModel
import StatsBase:   fit, coef, confint, stderr, vcov
import StatsBase:   loglikelihood, nullloglikelihood, deviance, nulldeviance
import StatsBase:   aic, aicc, bic, dof, dof_residual, nobs, r2, r², adjr2, adjr², coeftable
import StatsBase:   model_response, predict, fitted, residuals
import StatsModels: coefnames, Terms

include("./inference/corr.jl")
include("./inference/utils.jl")
include("./general/types.jl")
include("./general/utils.jl")
include("./data/types.jl")
include("./data/utils.jl")
include("./data/retrieval.jl")
include("./general/statsmodel.jl")
include("./general/etable.jl")
include("./inference/vcov.jl")
include("./inference/hausman.jl")
include("./estimation/logit.jl")
include("./estimation/probit.jl")
include("./estimation/ols.jl")
include("./estimation/iv.jl")
include("./estimation/abadie.jl")
include("./estimation/frölichmelly.jl")
include("./estimation/tan.jl")

export

    CorrStructure,
        Homoscedastic,
        Heteroscedastic,
        Clustered,
        CrossCorrelated,
            cc_timespace,
    Kernel,
        Bartlett,
        Truncated,
        TukeyHanning,
        Parzen, Gallant,
    Microdata,
    OLS, IV, Logit, Probit, Abadie, FrölichMelly, Tan,
    fit, first_stage, second_stage, reduced_form,
    coef, confint, pval, stderr, tstat, vcov, hausman,
    loglikelihood, nullloglikelihood, deviance, nulldeviance,
    aic, aicc, bic, dof, dof_residual, nobs, r2, r², adjr2, adjr², coeftable, etable,
    model_response, predict, fitted, residuals

end # module
