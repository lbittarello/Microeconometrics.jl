var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Microeconometrics.jl-1",
    "page": "Introduction",
    "title": "Microeconometrics.jl",
    "category": "section",
    "text": "This package provides support for microeconometric models. It supports complex covariance structures (clustered, etc.) and weighted data. Please report bugs by opening an issue. Information on specific versions can be found on the release page."
},

{
    "location": "index.html#Supported-estimators-1",
    "page": "Introduction",
    "title": "Supported estimators",
    "category": "section",
    "text": "More models are planned. If your preferred model is not currently available, file an issue or contribute!Models for exogenous regressors\nOrdinary least squares (OLS)\nLogit\nProbit\nInverse probability weighting (IPW)\nModels for endogenous regressors\nTwo-stage least squares\nAbadie (2002)\nFrölich and Melly (2013)\nTan (2006)"
},

{
    "location": "index.html#Package-manual-1",
    "page": "Introduction",
    "title": "Package manual",
    "category": "section",
    "text": "Pages = [\n        \"getting_started.md\",\n        \"model_specification.md\",\n        \"estimators.md\",\n        \"methods.md\",\n        \"contributing.md\",\n        \"to_do.md\",\n    ]\nDepth = 2"
},

{
    "location": "getting_started.html#",
    "page": "Getting started",
    "title": "Getting started",
    "category": "page",
    "text": ""
},

{
    "location": "getting_started.html#Getting-Started-1",
    "page": "Getting started",
    "title": "Getting Started",
    "category": "section",
    "text": ""
},

{
    "location": "getting_started.html#Installation-1",
    "page": "Getting started",
    "title": "Installation",
    "category": "section",
    "text": "To install the package, run:Pkg.add(\"Microeconometrics\")To update to the latest release, run:Pkg.update(\"Microeconometrics\")To obtain the last version, run:Pkg.checkout(\"Microeconometrics\")You will also need DataFrames, as well as a package to import your data (e.g., CSV)."
},

{
    "location": "getting_started.html#Example-I:-OLS-1",
    "page": "Getting started",
    "title": "Example I: OLS",
    "category": "section",
    "text": "For an introduction to the package, let's consider an example: ordinary least squares.We first load the modules:julia> using CSV\njulia> using DataFrames\njulia> using MicroeconometricsWe then load the data:julia> dst = CSV.read(\"admit.csv\")\njulia> categorical!(dst, :rank)This sample is available here, if you wish to replicate this exercise. The last line converts dst[:rank] into a categorical array.We now define the correlation structure of the error term. Let's assume that observations are independent and identically distributed, so that errors are homoscedastic:julia> corr = Homoscedastic();We next specify the model by constructing a Microdata:julia> dta = Microdata(dst, vcov = corr, response = \"admit\", control = \"gre + gpa + rank + 1\");We can finally fit the model and visualize the results:julia> e_ols = fit(OLS, dta);\njulia> coeftable(e_ols);\n\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \ngre             0.0004    0.0002    2.0384    0.0415       0.0  0.0008\ngpa             0.1555     0.064    2.4317    0.0150    0.0302  0.2809\nrank: 2        -0.1624    0.0677   -2.3978    0.0165   -0.2951 -0.0296\nrank: 3        -0.2906    0.0702   -4.1365    <1e-99   -0.4282 -0.1529\nrank: 4         -0.323    0.0793   -4.0726    <1e-99   -0.4785 -0.1676\n(Intercept)    -0.2589     0.216   -1.1987    0.2306   -0.6822  0.1644The corresponding code in Stata is:. import delimited \"admit.csv\"\n. regress admit gre gpa i.rankWe can apply the methods for regression models of StatsBase to obtain statistics:julia> nobs(e_ols)\n400\njulia> r2(e_ols)\n0.10040062851886422"
},

{
    "location": "getting_started.html#Example-II:-Comparing-models-1",
    "page": "Getting started",
    "title": "Example II: Comparing models",
    "category": "section",
    "text": "To illustrate more advanced features, suppose that we want to compare specifications. As a first exercise, we wonder if we could drop the rank fixed effects.We start with the full model. We now assume that the data are heteroscedastic:julia> corr = Heteroscedastic();\njulia> dta₁ = Microdata(dst, vcov = corr, response = \"admit\", control = \"gre + gpa + rank + 1\");\njulia> e₁ = fit(OLS, dta₁);\njulia> coeftable(e₁);\n\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \ngre             0.0004    0.0002    2.0501    0.0404       0.0  0.0008\ngpa             0.1555    0.0653    2.3833    0.0172    0.0276  0.2834\nrank: 2        -0.1624    0.0729   -2.2266    0.0260   -0.3053 -0.0194\nrank: 3        -0.2906     0.073   -3.9827    0.0001   -0.4336 -0.1476\nrank: 4         -0.323     0.078   -4.1408    <1e-99   -0.4759 -0.1701\n(Intercept)    -0.2589     0.211   -1.2268    0.2199   -0.6725  0.1547Before we estimate the reduced model, we need to redefine the control set. There are two approaches to this task. We can construct a new Microdata from scratch:julia> dta₂ = Microdata(dst, vcov = corr, response = \"admit\", control = \"gre + gpa + 1\");Or we can modify the control set of our existing Microdata:julia> dta₂ = Microdata(dta₁, control = \"gre + gpa + 1\");The second approach does not reconstruct the underlying data matrix. Therefore, it is faster and uses less memory. On the other hand, it only allows us to reassign existing variables across variable sets. We cannot add new variables or modify the correlation structure of the error term.We next fit the reduced model:julia> e₂ = fit(OLS, dta₂);\njulia> coeftable(e₂);\n\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \ngre             0.0005    0.0002    2.5642    0.0103    0.0001   0.001\ngpa             0.1542     0.065    2.3737    0.0176    0.0269  0.2816\n(Intercept)    -0.5279    0.2087   -2.5293    0.0114    -0.937 -0.1188The coeffients on gre and gpa seem to be robust.For a formal equality test, we use a Hausman test:julia> ht = hausman_1s(e₁, e₂, [\"gre\", \"gpa\"]);\njulia> tstat(ht)\n\n2-element Array{Float64,1}:\n -2.07838\n  0.0749552The function hausman_1s estimates the difference between two estimates based on the same sample. As it turns out, the difference between the coefficients on gre is statistically significant.To further investigate this result, we wish to estimate separate effects by rank. The keyword subset helps us construct the appropriate Microdatas:julia> idx₁ = (dst[:rank] .== 1);\njulia> dta₁ = Microdata(dst, subset = idx₁, vcov = corr, response = \"admit\", control = \"gre + gpa + 1\");\njulia> e₁   = fit(OLS, dta₁);\njulia> coeftable(e₁);\n\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \ngre             0.0006    0.0006    1.0386    0.2990   -0.0005  0.0017\ngpa             0.2508    0.1929       1.3    0.1936   -0.1273   0.629\n(Intercept)    -0.6806    0.5662   -1.2021    0.2293   -1.7903  0.4291\n\njulia> idx₂ = (dst[:rank] .== 2);\njulia> dta₂ = Microdata(dst, subset = idx₂, vcov = corr, response = \"admit\", control = \"gre + gpa + 1\");\njulia> e₂   = fit(OLS, dta₂);\njulia> coeftable(e₂);\n\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \ngre             0.0004    0.0004    0.8794    0.3792   -0.0004  0.0011\ngpa             0.1771    0.1028    1.7219    0.0851   -0.0245  0.3787\n(Intercept)     -0.447    0.3641   -1.2277    0.2196   -1.1607  0.2666\n\njulia> ht = hausman_2s(e₁, e₂, [\"gre\", \"gpa\"]);\njulia> tstat(ht)\n\n2-element Array{Float64,1}:\n 0.334261\n 0.337304We used the function hausman_2s because these estimates are based on different samples. The difference in the effect of gre between ranks 1 and 2 is not significant."
},

{
    "location": "model_specification.html#",
    "page": "Model specification",
    "title": "Model specification",
    "category": "page",
    "text": ""
},

{
    "location": "model_specification.html#Model-specification-1",
    "page": "Model specification",
    "title": "Model specification",
    "category": "section",
    "text": "Before fitting the model, one must specify the correlation pattern of the error term (a CorrStructure) and the components of the model (a Microdata)."
},

{
    "location": "model_specification.html#CorrStructure-1",
    "page": "Model specification",
    "title": "CorrStructure",
    "category": "section",
    "text": "This structure specifies the correlation between observations. It determines the calculation of standard errors.All constructors accept the Boolean keyword adj, which defaults to true. If true, a finite-sample adjustment is applied to the variance matrix. The adjustment factor is  / (n - 1), where n is the number of clusters for clustered data and the number of observations otherwise.Four subtypes are currently available:"
},

{
    "location": "model_specification.html#Homoscedastic-1",
    "page": "Model specification",
    "title": "Homoscedastic",
    "category": "section",
    "text": "Homoscedastic(method::String = \"OIM\"; adj::Bool = true)Observations are independent and identically distributed. The optional method argument controls the estimation of the variance matrix: \"OIM\" uses the observed information matrix, whereas \"OPG\" uses the outer product of the gradient. (They are equivalent for OLS.) Only OLS and maximum-likelihood estimators support homoscedastic errors."
},

{
    "location": "model_specification.html#Heteroscedastic-1",
    "page": "Model specification",
    "title": "Heteroscedastic",
    "category": "section",
    "text": "Heteroscedastic(; adj::Bool = true)Observations are independent, but they may differ in distribution. This structure leads to sandwich covariance matrices (a.k.a. Huber-Eicker-White)."
},

{
    "location": "model_specification.html#Clustered-1",
    "page": "Model specification",
    "title": "Clustered",
    "category": "section",
    "text": "Clustered(DF::Microdata, cluster::Symbol; adj::Bool = true)Observations are independent across clusters, but their joint distribution within clusters is arbitrary. cluster specifies the column of DF to cluster on."
},

{
    "location": "model_specification.html#CrossCorrelated-1",
    "page": "Model specification",
    "title": "CrossCorrelated",
    "category": "section",
    "text": "This structure accommodates other correlation structures. The first argument determines the precise pattern. The following methods are available:CrossCorrelated(\"Two-way clustering\", DF::DataFrame, cluster₁::Symbol, cluster₂::Symbol; adj::Bool = true)Observations may be arbitrarily correlated if they share any cluster.CrossCorrelated(\"Time\", DF::DataFrame, time::Symbol, bandwidth::Real; adj::Bool = true)\nCrossCorrelated(\"Time\", DF::DataFrame, time::Symbol, kernel::function; adj::Bool = true)The maximum possible correlation between two observations declines with the time difference between them. Correlation is arbitrary below that limit. The bandwidth and the kernel function control the upper bound. time specifies the column of DF that contains the date of each observation (Date).The first method takes a bandwidth and uses the Parzen kernel. The second method takes a kernel function instead, which must incorporate the bandwidth. The first method is equivalent to setting x -> parzen(x / bandwidth). The following kernels are predefined for convenience: Bartlett (bartlett), Parzen (parzen), Truncated (truncated) and Tukey-Hanning (tukeyhanning). See Andrews (1991) for formulae.CrossCorrelated(\"Space\", DF::DataFrame, latitude::Symbol, longitude::Symbol, bandwidth::Real; adj::Bool = true)\nCrossCorrelated(\"Space\", DF::DataFrame, latitude::Symbol, longitude::Symbol, kernel::function; adj::Bool = true)The maximum possible correlation between two observations declines with the spatial distance between them. Correlation is arbitrary below that limit. The bandwidth and the kernel function control the upper bound. latitude and longitude specify the columns of DF that contain the coordinates of each observation (Float64).For an explanation of the difference between the two methods, see CrossCorrelated(\"time\", args...) above.CrossCorrelated(\"Time and space\",\n    DF::DataFrame,\n    time::Symbol,\n    bandwidth_time::Real,\n    latitude::Symbol,\n    longitude::Symbol,\n    bandwidth_space::Real;\n    adj::Bool = true)\nCrossCorrelated(\"Time and space\",\n    DF::DataFrame,\n    time::Symbol,\n    kernel_time::Function,\n    latitude::Symbol,\n    longitude::Symbol,\n    kernel_space::Function;\n    adj::Bool = true)The maximum possible correlation between two observations declines with the time difference and the spatial distance between them. Correlation is arbitrary below that limit. The bandwidths and the kernel functions control the upper bound. time specifies the column of DF that contains the date of each observation. latitude and longitude specify the columns of DF that contain the coordinates of each observation (Float64).For an explanation of the difference between the two methods, see CrossCorrelated(\"time\", args...) above."
},

{
    "location": "model_specification.html#Microdata-1",
    "page": "Model specification",
    "title": "Microdata",
    "category": "section",
    "text": "This structure combines the functionalities of Formula, ModelFrame and ModelMatrix from StatsModels. It contains a Matrix{Float64} (the data matrix), a map from model components to matrix columns, a correlation structure and weights (inter alia).Microdata(\n    DF::DataFrame;\n    vcov::CorrStructure = Heteroscedastic(),\n    weights::AbstractWeights = UnitWeights(size(DF, 1)),\n    subset::AbstractVector{Bool} = trues(size(DF, 1)),\n    kwargs...)subset determines the estimation sample. Set a row to true if the corresponding row of DF should be included and false if it should be excluded. This keyword is useful in two situations. First, you have precomputed vcov based on the entire sample. Microdata will copy the correlation structure and restrict it to relevant observations. Second, you are comparing subgroup effects and observations in different subgroups may correlate (e.g., they may belong to the same cluster). hausman_2s will account for that correlation if the Microdatas were constructed with subset.weights is a weight vectors. Except for frequency weights, the weight vector is normalized to sum up to the number of observations in the sample.Additional keywords determine the model components. All regression models need a response, but other requirements may vary. For example, OLS asks for response and control. See the introduction for examples. Conventional sets include:response: the response, outcome or dependent variable;\ncontrol: exogenous explanatory variables (n.b.: you must explicitly include intercepts);\ntreatment: endogenous explanatory variables;\ninstrument: instrumental variables (i.e., excluded exogenous variables).You pass these sets as strings, following syntax of Formula.Microdata(MD::Microdata; kwargs...)It is also possible to base new Microdata on existing Microdata. This constructor allows you to reassign variables to new sets. You can create new variable sets. If you do not redefine a set, it is preserved. To suppress a set, redefine it to \"\". You cannot add new variables, modify the correlation structure, restrict the sample or reweight observations.This functionality is useful if you wish to compare specifications. Rather than building separate data matrices for each one of them, you can build a master Microdata, holding all variables of interest, and adjust its map as you go through specifications."
},

{
    "location": "estimators.html#",
    "page": "Estimators",
    "title": "Estimators",
    "category": "page",
    "text": ""
},

{
    "location": "estimators.html#Estimators-1",
    "page": "Estimators",
    "title": "Estimators",
    "category": "section",
    "text": "The function fit estimates models. It returns a model structure, which contains the estimation sample, coefficients and their variance matrix (inter alia). For example, the output of fit(OLS, MD) has type OLS. Some have additional fields: for instance, the structures of two-stage models carry estimates from the first stage. These structure are instances of broader abstract types, such as MLE or TwoStageModel, which belong in turn to the supertype Micromodel.If you only need coefficients, set novar = true."
},

{
    "location": "estimators.html#Models-for-exogenous-regressors-1",
    "page": "Estimators",
    "title": "Models for exogenous regressors",
    "category": "section",
    "text": ""
},

{
    "location": "estimators.html#Ordinary-least-squares-1",
    "page": "Estimators",
    "title": "Ordinary least squares",
    "category": "section",
    "text": "fit(OLS, MD::Microdata; novar::Bool = false)The Microdata must contain: response and control. OLS is a subtype of ParModel."
},

{
    "location": "estimators.html#Binary-choice-1",
    "page": "Estimators",
    "title": "Binary choice",
    "category": "section",
    "text": "fit(Logit, MD::Microdata; novar::Bool = false)The Microdata must contain: response and control. Logit is a subtype of MLE and ParModel.fit(Probit, MD::Microdata; novar::Bool = false)The Microdata must contain: response and control. Probit is a subtype of MLE and ParModel."
},

{
    "location": "estimators.html#Treatment-evaluation-1",
    "page": "Estimators",
    "title": "Treatment evaluation",
    "category": "section",
    "text": "fit(IPW,\n    M₁::Type{Micromodel},\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0,\n    kwargs...)\nfit(IPW,\n    m₁::Micromodel,\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0)The Microdata must contain: response, treatment and control. The treatment must be binary. IPW is a subtype of TwoStageModel.This model estimates the average treatment effect by inverse probability weighting. In a first stage, we use model M₁ to forecast the conditional probability of treatment take-up (the propensity score) and construct weights, so that the weighted covariate distribution is similar across treatment subsamples. In the second stage, we run weighted OLS on the treatment and an intercept. The intercept gives the mean outcome in the absence of treatment. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The first method fits the first-stage model. Keyword arguments customize this step. The second method uses a previously estimated model instead."
},

{
    "location": "estimators.html#Models-for-endogenous-regressors-1",
    "page": "Estimators",
    "title": "Models for endogenous regressors",
    "category": "section",
    "text": ""
},

{
    "location": "estimators.html#Linear-IV-1",
    "page": "Estimators",
    "title": "Linear IV",
    "category": "section",
    "text": "fit(IV, MD::Microdata; novar::Bool = false, method::String = \"TSLS\")The Microdata must contain: response, treatment, control and instrument. IV is a subtype of ParModel.The following variants are currently implemented:method = \"TSLS\": two-stage least squares;\nmethod = \"OLS\": linear regression of the outcome on the treatment and controls;\nmethod = \"First stage\": linear regression of the treatment on the instruments and controls;\nmethod = \"Reduced form\": linear regression of the outcome on the instruments and controls."
},

{
    "location": "estimators.html#Reweighting-models-1",
    "page": "Estimators",
    "title": "Reweighting models",
    "category": "section",
    "text": "fit(Abadie,\n    M₂::Type{ParModel},\n    M₁::Type{Micromodel},\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0,\n    kwargs...)\nfit(Abadie,\n    M₂::Type{ParModel},\n    m₁::Micromodel,\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0\n    kwargs...)The Microdata must contain: response, treatment and control. The treatment and the instrument must be binary. Abadie is a subtype of TwoStageModel.This model estimates a local average response function according to Abadie (2003). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct weights, so that the weighted covariate distribution is similar across treatment subsamples of the compliers. In a second stage, we fit model M₂ with the weights from the first stage. Keywords customize the second-stage estimator. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The first method fits the first-stage model. The second method uses a previously estimated model instead.fit(FrölichMelly,\n    M₁::Type{Micromodel},\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0,\n    kwargs...)\nfit(FrölichMelly,\n    m₁::Micromodel,\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0)The Microdata must contain: response, treatment and control. The treatment and the instrument must be binary. FrölichMelly is a subtype of TwoStageModel.This model estimates the unconditional local average treatment effects according to Frölich and Melly (2013). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct weights, so that the weighted covariate distribution is similar across treatment subsamples of the compliers. In the second stage, we run weighted OLS on the treatment and an intercept. The intercept gives compliers' mean outcome in the absence of treatment. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The first method fits the first-stage model. Keyword arguments customize this step. The second method uses a previously estimated model instead.fit(Tan,\n    M₁::Type{Micromodel},\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0,\n    kwargs...)\nfit(Tan,\n    m₁::Micromodel,\n    MD::Microdata;\n    novar::Bool = false,\n    trim::AbstractFloat = 0.0)The Microdata must contain: response, treatment and control. The treatment and the instrument must be binary. Tan is a subtype of TwoStageModel.This model estimates the unconditional local average treatment effects according to Tan (2006). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct weights, so that the weighted covariate distribution is similar across treatment subsamples of the compliers. In the second stage, we run weighted OLS on the treatment and an intercept. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The first method fits the first-stage model. Keyword arguments customize this step. The second method uses a previously estimated model instead."
},

{
    "location": "methods.html#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "methods.html#Methods-1",
    "page": "Methods",
    "title": "Methods",
    "category": "section",
    "text": ""
},

{
    "location": "methods.html#Methods-from-*StatsBase*-1",
    "page": "Methods",
    "title": "Methods from StatsBase",
    "category": "section",
    "text": "This package supports the methods for regression models of StatsBase.These following functions are available for all parametric models: nobs, dof, dof_residual, model_response, coef, stderr, vcov, confint, coefnames and coeftable. Three keywords customize the behavior of coeftable: verbose (Boolean: suppresses printing), digits (integer: controls rounding) and level (float: controls the level of the confidence interval – 0.0 for none). Note that all methods refer to the second stage of two-stage models.The following functions are available for maximum-likelihood estimators: aic, aicc, bic, deviance, nulldeviance, loglikelihood, nullloglikelihood, r2 and adjr2. Both R² functions are also available for OLS and IV.Some models support predict, fitted and residuals. predict estimates the index of single-index models. fitted estimates the conditional outcome expectation. For example, predict estimates the Xβ of a logit model, whereas fitted estimates logistic(Xβ).We also implement additional methods:tstat: the t-statistic (i.e., the ratio of coefficients to standard error);\npval: the p-value of a two-sided significance test;\nfirst_stage: the first-stage estimates of a two-stage model."
},

{
    "location": "methods.html#Hausman-test-1",
    "page": "Methods",
    "title": "Hausman test",
    "category": "section",
    "text": "These function computes the difference in coefficients between two parametric models. They return a ParObject, which contains the vector of differences, their variance matrix and labels. Our implementation is based on the GMM representation of the joint estimation problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)).hausman_1s(\n    model₁::Union{ParModel, TwoStageModel},\n    model₂::Union{ParModel, TwoStageModel},\n    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))This function is appropriate when both model₁ and model₂ were based on the same estimation sample.hausman_2s(\n    model₁::Union{ParModel, TwoStageModel},\n    model₂::Union{ParModel, TwoStageModel},\n    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))This function is appropriate when both model₁ and model₂ were based on independent samples. For example, the samples might consist of independent observations with no overlap.hausman_2s(\n    model₁::Union{ParModel, TwoStageModel},\n    model₂::Union{ParModel, TwoStageModel},\n    corr::CorrStructure,\n    names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂)))This function is appropriate when both model₁ and model₂ were based on dependent samples. For example, the samples might consist of independent observations with some overlap or clustered observations with common clusters. The correlation structure corr must specify the correlation between all observations of both estimation samples. For example, you could precompute corr for the entire dataset and construct the two estimation samples via the subset keyword to Microdata."
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "This section explains how to specify your own estimator. You can then submit a pull request to add it to Microeconometrics!We will analyze the implementation of OLS. Although it is a simple model, other models follow the same steps.The first step is defining the output struct:mutable struct OLS <: ParModel\n\n    sample::Microdata     # estimation sample\n    β::Vector{Float64}    # coefficient vector\n    V::Matrix{Float64}    # variance matrix\n\n    OLS() = new()\nendEvery estimator has these fields. Internal utilities rely on them. Some have additional fields (e.g., IV stores the estimation method). Two-stage models store the estimators for each stage instead.We the define an uninitialized constructor:function OLS(MD::Microdata)\n    obj        = OLS()\n    obj.sample = MD\n    return obj\nendThe functions _fit! and _vcov! will later set β and V.The next step is overloading the function fit from StatsBase:function fit(::OLS, MD::Microdata; novar::Bool = false)\n\n    obj = OLS(MD)\n\n    _fit!(obj, getweights(obj))\n    novar || _vcov!(obj, getcorr(obj), getweights(obj))\n\n    return obj\nendFor OLS, we only need to initialize the output object and pass it to _fit! and _vcov!. (It is not actually necessary to extend fit unless you need to perform additional steps before the estimation, as the default implementation will suffice.) Note the utilities getcorr and getweights.We can now estimate the coefficients. For efficiency, we write separate functions for unweighted and weighted data.function _fit!(obj::OLS, w::UnitWeights)\n    y     = getvector(obj, :response)\n    x     = getmatrix(obj, :control)\n    obj.β = x \\ y\nend\n\nfunction _fit!(obj::OLS, w::AbstractWeights)\n    y     = getvector(obj, :response)\n    x     = getmatrix(obj, :control)\n    v     = copy(x) .* values(w)\n    obj.β =  (v' * x) \\ (v' * y)\nendNotice the internal utilities getvector and getmatrix. Their first argument is a Microdata or a model structure. The following arguments are the model components of interest. You can request several components from getmatrix at once. For example, IV needs x = getmatrix(obj, :treatment, :control) and z = getmatrix(obj, :instrument, :control). A single matrix is returned in all cases.warning: Warning\ngetvector and getmatrix return views into the underlying data matrix. You should never modify their output, as you would irremediably alter the data. If you need to perform some operation, make a copy beforehand.OLS does not require nonlinear optimization. If your estimator needs it, you can use the tools of Optim. See the implementation of Logit for an example.We must now define score and jacobian. These functions are the building blocks of the variance estimator. The score is the vector of moment conditions. For OLS, it is −xᵢ (yᵢ - xᵢ'β) (the derivative of the objective function). score should return the matrix of score vectors (in row form). The Jacobian matrix is the derivative of the moment conditions. For OLS, it is xᵢ'xᵢ. jacobian should return the weighted sum of Jacobians (i.e., the expected Jacobian × the number of observations).function score(obj::OLS)\n    x = copy(getmatrix(obj, :control))\n    û = residuals(obj)\n    return - x .* û\nend\n\nfunction jacobian(obj::OLS, w::UnitWeights)\n    x = getmatrix(obj, :control)\n    return x' * x\nend\n\nfunction jacobian(obj::OLS, w::AbstractWeights)\n    x = getmatrix(obj, :control)\n    v = copy(x) .* values(w)\n    return x' * v\nendscore returns the score for each observation, so it ignores weights. jacobian returns an expectation; therefore, it must account for weights.We do not need to extend _vcov!. The default method will call score and jacobian and construct the appropriate estimator, accounting for the correlation structure of the data and the type of weights.We now overload predict and fitted. For OLS, these functions are equivalent.predict(obj::OLS) = getmatrix(obj, :control) * obj.β\nfitted(obj::OLS)  = predict(obj)The next step is optional. We extend jacobexp, which computes the derivative of fitted values.jacobexp(obj::OLS) = getmatrix(obj, :control)jacobexp is only necessary when the estimator in question serves as the first stage of a two-stage estimator. By extending it, you make your estimator available to two-stage estimators.We conclude with a function to retrieve coefficient labels:coefnames(obj::OLS) = getnames(obj, :control)You can implement additional methods. For example, Microeconometrics extends r2 and adjr2 to OLS (the default method only cover subtypes of MLE)."
},

{
    "location": "to_do.html#",
    "page": "To do",
    "title": "To do",
    "category": "page",
    "text": ""
},

{
    "location": "to_do.html#To-do-1",
    "page": "To do",
    "title": "To do",
    "category": "section",
    "text": "More models;\nMarginal effects for nonlinear models;\nA framework for panel data;\nA framework for GMM;\nA framework for semiparametric and nonparametric models;\nWald test;\nLikelihood-based tests;\nBootstrapping."
},

]}
