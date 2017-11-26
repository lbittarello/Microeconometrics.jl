__precompile__(true)

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
include("./general/utils.jl")
include("./inference/corr.jl")
include("./inference/utils.jl")
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
    first_stage, second_stage, reduced_form, tstat, pval, etable, hausman

end # module
