__precompile__()

module Microeconometrics

using DataFrames
using StatsFuns
using Optim
using Formatting

import Base:       copy
import StatsBase:  fit, loglikelihood, nullloglikelihood, deviance, nulldeviance
import StatsBase:  coef, coeftable, confint, vcov
import StatsBase:  dof, dof_residual, nobs, r2, adjr2
import StatsBase:  model_response, predict, fitted, residuals
import DataFrames: coefnames, Terms

include("./general/types.jl")
include("./general/corr.jl")
include("./general/corr_utils.jl")
include("./general/data.jl")
include("./general/data_utils.jl")
include("./general/retrieval.jl")
include("./general/statsmodel.jl")
include("./general/etable.jl")
include("./general/utils.jl")
include("./inference/vcov.jl")
include("./inference/hausman.jl")
include("./MLE/logit.jl")
include("./MLE/probit.jl")
include("./GMM/ols.jl")
include("./GMM/iv.jl")
include("./TSM/abadie.jl")
include("./TSM/frölichmelly.jl")

export CorrStructure, Microdata
export Homoscedastic, Heteroscedastic, Clustered, CrossCorrelated
export cc_timespace, truncated, bartlett, parzen, tukeyhanning
export ParObject, OLS, IV, Logit, Probit, Abadie, FrölichMelly
export first_stage, second_stage, reduced_form, tstat, pval, etable
export hausman_test

export fit, loglikelihood, nullloglikelihood, deviance, nulldeviance
export coefnames, coef, coeftable, confint, vcov
export dof, dof_residual, nobs, r2, adjr2
export model_response, predict, fitted, residuals

end # module
