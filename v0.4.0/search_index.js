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
    "text": "This package provides support for microeconometric estimation. It supports complex weighted data and covariance structures (e.g., clustered). Please report bugs by opening an issue. Information about specific versions can be found on the release page."
},

{
    "location": "index.html#Supported-estimators-1",
    "page": "Introduction",
    "title": "Supported estimators",
    "category": "section",
    "text": "More models are planned. If your preferred model is not currently available, file an issue or contribute!Linear regression\nOrdinary least squares\nTwo-stage least squares\nLinear GMM\nBinary choice\nLogit\nProbit\nComplementary log-log\nCount data\nPoisson\nIV Poisson with additive errors\nIV Poisson with multiplicative errors\nReweighting methods\nInverse probability weighting\nAbadie (2003)\nFrölich and Melly (2013)\nTan (2006)"
},

{
    "location": "index.html#Package-manual-1",
    "page": "Introduction",
    "title": "Package manual",
    "category": "section",
    "text": "Pages = [\r\n        \"getting_started.md\",\r\n        \"model_specification.md\",\r\n        \"correlation_structures.md\",\r\n        \"estimators.md\",\r\n        \"methods.md\",\r\n        \"hypothesis_tests.md\",\r\n        \"estimation_tables.md\",\r\n        \"bootstrapping.md\",\r\n        \"contributing.md\",\r\n        \"to_do.md\",\r\n    ]\r\nDepth = 2"
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
    "text": "To install the latest release, run:julia> Pkg.add(\"Microeconometrics\")To update to the latest release, run:julia> Pkg.update(\"Microeconometrics\")To install the last version, run:julia> Pkg.add(PackageSpec(name = \"Microeconometrics\", rev = \"master\"))To update to the last version, run:julia> Pkg.update(PackageSpec(name = \"Microeconometrics\", rev = \"master\"))"
},

{
    "location": "getting_started.html#Example-I:-OLS-1",
    "page": "Getting started",
    "title": "Example I: OLS",
    "category": "section",
    "text": "For an introduction to the package, let\'s consider an example: ordinary least squares.We first load some modules:julia> using CSV\r\njulia> using DataFrames\r\njulia> using MicroeconometricsNext, we load the data:julia> S = CSV.read(\"admit.csv\") ; categorical!(S, :rank) ;This sample is available here, if you wish to replicate this exercise. It comprises 400 observations and four variables (admit, gre, gpa and rank). The second command converts S[:rank] into a categorical array (a.k.a. factor variable).We wish to regress admit on gre, gpa, rank and an intercept. We specify this model via a dictionary:julia> M = Dict(:response => \"admit\", :control => \"gre + gpa + rank + 1\") ;The syntax is due to StatsModels.jl.We then define the correlation structure of the errors. Let\'s assume that observations are independent and identically distributed, so that errors are homoscedastic:julia> C = Homoscedastic() ;We now construct the estimation sample:julia> D = Microdata(S, M, vcov = C) ;We can finally fit the model and visualize the results:julia> E = fit(OLS, D) ;\r\njulia> coeftable(E)\r\n\r\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \r\ngre             0.0004    0.0002    2.0384    0.0415       0.0  0.0008\r\ngpa             0.1555     0.064    2.4317    0.0150    0.0302  0.2809\r\nrank: 2        -0.1624    0.0677   -2.3978    0.0165   -0.2951 -0.0296\r\nrank: 3        -0.2906    0.0702   -4.1365    <1e-99   -0.4282 -0.1529\r\nrank: 4         -0.323    0.0793   -4.0726    <1e-99   -0.4785 -0.1676\r\n(Intercept)    -0.2589     0.216   -1.1987    0.2306   -0.6822  0.1644The corresponding code in Stata is:. import delimited \"admit.csv\"\r\n. regress admit gre gpa i.rankTo obtain summary statistics, we can apply the methods for regression models of StatsBase.jl:julia> nobs(e_ols)\r\n\r\n400\r\n\r\njulia> r2(e_ols)\r\n\r\n0.10040062851886422"
},

{
    "location": "getting_started.html#Example-II:-Comparing-models-1",
    "page": "Getting started",
    "title": "Example II: Comparing models",
    "category": "section",
    "text": "To illustrate more advanced features, suppose that we want to compare specifications. As a first exercise, we wonder if we could drop the rank fixed effects.We start with the full model. We now assume that the data are heteroscedastic. Because it is the default, we need not specify it.julia> M₁ = Dict(:response => \"admit\", :control => \"gre + gpa + rank + 1\") ;\r\njulia> D₁ = Microdata(S, M₁);\r\njulia> E₁ = fit(OLS, D₁);\r\njulia> coeftable(E₁)\r\n\r\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \r\ngre             0.0004    0.0002    2.0501    0.0404       0.0  0.0008\r\ngpa             0.1555    0.0653    2.3833    0.0172    0.0276  0.2834\r\nrank: 2        -0.1624    0.0729   -2.2266    0.0260   -0.3053 -0.0194\r\nrank: 3        -0.2906     0.073   -3.9827    0.0001   -0.4336 -0.1476\r\nrank: 4         -0.323     0.078   -4.1408    <1e-99   -0.4759 -0.1701\r\n(Intercept)    -0.2589     0.211   -1.2268    0.2199   -0.6725  0.1547Before we estimate the reduced model, we must redefine the control set. There are two approaches to this task. We can construct a new Microdata from scratch:julia> M₂ = Dict(:response => \"admit\", :control => \"gre + gpa + 1\") ;\r\njulia> D₂ = Microdata(S, M₂) ;Or we can modify the control set of our existing Microdata:julia> M₂ = Dict(:control => \"gre + gpa + 1\") ;\r\njulia> D₂ = Microdata(D₁, M₂) ;The second approach does not reconstruct the underlying data matrix. Therefore, it is faster and uses less memory. On the other hand, it only allows us to reassign existing variables across variable sets. We cannot add new variables to the data matrix or modify the correlation structure of the error term.We can now fit the reduced model:julia> E₂ = fit(OLS, D₂) ;\r\njulia> coeftable(E₂)\r\n\r\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \r\ngre             0.0005    0.0002    2.5642    0.0103    0.0001   0.001\r\ngpa             0.1542     0.065    2.3737    0.0176    0.0269  0.2816\r\n(Intercept)    -0.5279    0.2087   -2.5293    0.0114    -0.937 -0.1188The coefficients on gre and gpa seem to be robust. For a formal equality test, we use a Hausman test:julia> H = hausman_1s(E₁, E₂, [\"gre\", \"gpa\"]) ;\r\njulia> tstat(H)\r\n\r\n2-element Array{Float64,1}:\r\n  -2.07838\r\n  0.0749552As it turns out, the difference between the coefficients on gre is statistically significant.To further investigate this result, we wish to estimate separate effects by rank. The keyword subset helps us construct the appropriate Microdata:julia> M = Dict(:response => \"admit\", :control => \"gre + gpa + 1\") ;\r\n\r\njulia> I₁ = (S[:rank] .== 1);\r\njulia> D₁ = Microdata(S, M, subset = I₁) ;\r\njulia> E₁ = fit(OLS, D₁) ;\r\njulia> coeftable(E₁)\r\n\r\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \r\ngre             0.0006    0.0006    1.0386    0.2990   -0.0005  0.0017\r\ngpa             0.2508    0.1929       1.3    0.1936   -0.1273   0.629\r\n(Intercept)    -0.6806    0.5662   -1.2021    0.2293   -1.7903  0.4291\r\n\r\njulia> I₂ = (S[:rank] .== 2) ;\r\njulia> D₂ = Microdata(S, M, subset = I₂) ;\r\njulia> E₂ = fit(OLS, D₂) ;\r\njulia> coeftable(E₂)\r\n\r\n              Estimate  St. Err.   t-stat.   p-value      C.I. (95%)  \r\ngre             0.0004    0.0004    0.8794    0.3792   -0.0004  0.0011\r\ngpa             0.1771    0.1028    1.7219    0.0851   -0.0245  0.3787\r\n(Intercept)     -0.447    0.3641   -1.2277    0.2196   -1.1607  0.2666\r\n\r\njulia> H = hausman_2s(E₁, E₂, [\"gre\", \"gpa\"]) ;\r\njulia> tstat(H)\r\n\r\n2-element Array{Float64,1}:\r\n  0.334261\r\n  0.337304We used the function hausman_2s because these estimates are based on different samples. The difference in the effect of gre between ranks 1 and 2 is not significant."
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
    "text": "Before you estimate a model, you must load the data and specify the role of each variable in the model. This packages defines Microdata for that purpose. This structure combines the functionalities of ModelFrame and ModelMatrix from StatsModels.jl.Microdata(\r\n        DF::DataFrame,\r\n        model::Dict{Symbol, String};\r\n        contrasts::Dict,\r\n        subset::AbstractVector{Bool} = trues(size(DF, 1)),\r\n        vcov::CorrStructure = Heteroscedastic(),\r\n        weights::AbstractWeights = UnitWeights(size(DF, 1))\r\n    )To construct a Microdata, two arguments are compulsory: a DataFrame, which contains the data, and a dictionary, which specifies the components of the model of interest.All regression models need a response, but other requirements may vary. (Check the documentation!) For example, OLS asks for response and control. You pass these sets as strings, following the syntax of Formula. (See the tutorial for examples.) Conventional sets include:response: the response (a.k.a. outcome or dependent variable);\ncontrol: exogenous explanatory variables (n.b.: you must explicitly include intercepts, + 1);\noffset: an exogenous variable whose coefficient is constrained to unity;\ntreatment: endogenous explanatory variables;\ninstrument: instrumental variables (i.e. excluded exogenous variables).As for the keywords:contrasts: a dictionary from column labels to contrast schemes.\nsubset determines the estimation sample. Set an entry to true if the corresponding row of DF should be included and false if it should be excluded. This keyword is useful if you are comparing subgroups and observations in different subgroups may correlate (e.g., they may belong to the same cluster). hausman_2s will account for that correlation if the Microdata were constructed with subset.\nweights is a weight vector. Except for frequency weights, the weight vector is normalized to sum up to the number of observations in the sample.\nvcov is a correlation structure.It is also possible to base new Microdata on existing Microdata:Microdata(MD::Microdata, model::Dict{Symbol, String})This constructor allows you to reassign variables to new sets. You can also create new variable sets. If you do not redefine a set, it is preserved. To suppress a set, redefine it to \"\". You cannot add new variables to the data matrix, modify the correlation structure, restrict the sample or reweight observations.This functionality is useful if you wish to compare specifications. Rather than building separate data matrices for each one of them, you can build a master Microdata, holding all variables of interest, and adjust its map as you go through specifications."
},

{
    "location": "correlation_structures.html#",
    "page": "Correlation structures",
    "title": "Correlation structures",
    "category": "page",
    "text": ""
},

{
    "location": "correlation_structures.html#Correlation-structures-1",
    "page": "Correlation structures",
    "title": "Correlation structures",
    "category": "section",
    "text": "Before fitting the model, you must specify the correlation between observations (a CorrStructure). It determines the calculation of covariance matrices. The default is always Heteroscedastic, i.e. independent but not identically distributed observations.All constructors accept the Boolean keyword adj (omitted in the following), which defaults to true. If true, a finite-sample adjustment is applied to the covariance matrix. The adjustment factor is n / (n - 1), where n is the number of clusters for clustered data and the number of observations otherwise.Four subtypes are currently available: Homoscedastic, Heteroscedastic, Clustered and CrossCorrelated."
},

{
    "location": "correlation_structures.html#Homoscedastic-1",
    "page": "Correlation structures",
    "title": "Homoscedastic",
    "category": "section",
    "text": "Homoscedastic(method::String = \"OIM\")Observations are independent and identically distributed. The optional argument method is only relevant for maximum-likelihood estimators. It controls the estimation of the covariance matrix: \"OIM\" uses the observed information matrix, whereas \"OPG\" uses the outer product of the gradient. Only linear and maximum-likelihood estimators support homoscedastic errors."
},

{
    "location": "correlation_structures.html#Heteroscedastic-1",
    "page": "Correlation structures",
    "title": "Heteroscedastic",
    "category": "section",
    "text": "Heteroscedastic()Observations are independent, but they may differ in distribution. This structure leads to sandwich covariance matrices (a.k.a. Huber-Eicker-White)."
},

{
    "location": "correlation_structures.html#Clustered-1",
    "page": "Correlation structures",
    "title": "Clustered",
    "category": "section",
    "text": "Clustered(DF::DataFrame, cluster::Symbol)Observations are independent across clusters, but they may differ in their joint distribution within clusters. cluster specifies the column of the DataFrame to cluster on."
},

{
    "location": "correlation_structures.html#CrossCorrelated-1",
    "page": "Correlation structures",
    "title": "CrossCorrelated",
    "category": "section",
    "text": "This structure accommodates other correlation structures. The first argument determines the precise pattern."
},

{
    "location": "correlation_structures.html#Two-way-clustering-1",
    "page": "Correlation structures",
    "title": "Two-way clustering",
    "category": "section",
    "text": "CrossCorrelated(\"Two-way clustering\", DF::DataFrame, c₁::Symbol, c₂::Symbol)if two observations share any cluster, they may be arbitrarily correlated."
},

{
    "location": "correlation_structures.html#Correlation-across-time-1",
    "page": "Correlation structures",
    "title": "Correlation across time",
    "category": "section",
    "text": "CrossCorrelated(\"Time\",\r\n        DF::DataFrame,\r\n        time::Symbol,\r\n        bandwidth::Real,\r\n        kernel::Function = parzen\r\n    )The maximum possible correlation between two observations declines with the time difference between them. The actual correlation is arbitrary below that limit. (See Conley (1999).) The bandwidth and the kernel function control the upper bound. time specifies the column of DF that contains the date of each observation (of type Date).The following kernels are predefined for convenience: Bartlett (bartlett), Parzen (parzen), Truncated (truncated) and Tukey-Hanning (tukeyhanning). See Andrews (1991) for formulae.warning: Warning\nThe resulting covariance matrices differ from the Newey-West estimator, which assumes independence across units (though observations for the same unit may correlate across time)."
},

{
    "location": "correlation_structures.html#Correlation-across-space-1",
    "page": "Correlation structures",
    "title": "Correlation across space",
    "category": "section",
    "text": "CrossCorrelated(\"Space\",\r\n        DF::DataFrame,\r\n        latitude::Symbol,\r\n        longitude::Symbol,\r\n        bandwidth::Real,\r\n        kernel::Function = parzen\r\n    )The maximum possible correlation between two observations declines with the spatial distance between them. The actual correlation is arbitrary below that limit. (See Conley (1999).) The bandwidth and the kernel function control the upper bound. latitude and longitude specify the columns of DF that contain the coordinates of each observation in radians (of type Float64).The following kernels are predefined for convenience: Bartlett (bartlett), Parzen (parzen), Truncated (truncated) and Tukey-Hanning (tukeyhanning). See Andrews (1991) for formulae."
},

{
    "location": "correlation_structures.html#Correlation-across-time-and-space-1",
    "page": "Correlation structures",
    "title": "Correlation across time and space",
    "category": "section",
    "text": "CrossCorrelated(\"Time and space\",\r\n        DF::DataFrame,\r\n        time::Symbol,\r\n        bandwidth_time::Real,\r\n        latitude::Symbol,\r\n        longitude::Symbol,\r\n        bandwidth_space::Real,\r\n        kernel::Function = parzen\r\n    )The maximum possible correlation between two observations declines with the time difference and the spatial distance between them. The actual correlation is arbitrary below that limit. (See Conley (1999).) The bandwidths and the kernel function control the upper bound. time specifies the column of DF that contains the date of each observation. latitude and longitude specify the columns of DF that contain the coordinates of each observation in radians (Float64).The following kernels are predefined for convenience: Bartlett (bartlett), Parzen (parzen), Truncated (truncated) and Tukey-Hanning (tukeyhanning). See Andrews (1991) for formulae."
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
    "text": "The function fit estimates models. It returns a model structure, which contains the estimation sample, the coefficients and their covariance matrix. For example, the output of fit(OLS, MD) has type OLS. Some have additional fields: e.g., two-stage models carry estimates from the first stage and GMM models carry the inverse of the weight matrix.note: Note\nIf you only need coefficients, pass novar = true to fit.Model structures are subtypes of broader abstract types, such as MLE or GMM, which are ultimately instances of RegressionModel. The type hierarchy is:RegressionModel\r\n    GMM\r\n    ParModel\r\n        MLE\r\n    TwoStageModel"
},

{
    "location": "estimators.html#Linear-regression-1",
    "page": "Estimators",
    "title": "Linear regression",
    "category": "section",
    "text": ""
},

{
    "location": "estimators.html#Ordinary-least-squares-1",
    "page": "Estimators",
    "title": "Ordinary least squares",
    "category": "section",
    "text": "fit(OLS, MD::Microdata)The Microdata must contain: response and control. See the documentation for linear IV if Microdata includes a treatment. OLS is a subtype of ParModel."
},

{
    "location": "estimators.html#Linear-IV-1",
    "page": "Estimators",
    "title": "Linear IV",
    "category": "section",
    "text": "fit(IV, MD::Microdata; method::String = \"TSLS\")The following methods are currently implemented:method = \"TSLS\": two-stage least squares;\nmethod = \"Two-step GMM\": optimally weighted two-stage GMM with robust covariance matrix;\nmethod = \"Optimal GMM\": optimally weighted two-stage GMM with simplified covariance matrix.Additional methods are available for convenience:method = \"OLS\": linear regression of the outcome on the treatment and controls;\nmethod = \"First stage\": linear regression of the treatment on the instruments and controls;\nmethod = \"Reduced form\": linear regression of the outcome on the instruments and controls.The Microdata must contain: response, treatment, control and instrument. IV is a subtype of GMM."
},

{
    "location": "estimators.html#Binary-choice-1",
    "page": "Estimators",
    "title": "Binary choice",
    "category": "section",
    "text": "fit(Logit, MD::Microdata)\r\nfit(Probit, MD::Microdata)\r\nfit(Cloglog, MD::Microdata)The Microdata must contain: response and control. The outcome should be binary. The model structures are subtypes of MLE."
},

{
    "location": "estimators.html#Count-data-1",
    "page": "Estimators",
    "title": "Count data",
    "category": "section",
    "text": "fit(Poisson, MD::Microdata; novar::Bool = false)The Microdata must contain: response and control. The Microdata may contain an offset. See the documentation for linear IV if Microdata includes a treatment. The outcome must be weakly positive. Poisson is a subtype of MLE.fit(IVPoisson, MD::Microdata; novar::Bool = false, method::String = \"One-step GMM\")\r\nfit(Mullahy, MD::Microdata; novar::Bool = false, method::String = \"One-step GMM\")IVPoisson fits the exponential conditional mean model with additive errors. Mullahy fits the exponential conditional mean model with multiplicative errors (Mullahy, 1997).The following methods are currently implemented:method = \"One-step GMM\": unweighted one-stage GMM;\nmethod = \"TSLS\": one-stage GMM with the average outer product of the instrument vector as weight matrix;\nmethod = \"Two-step GMM\": optimally weighted two-stage GMM with robust covariance matrix;\nmethod = \"Optimal GMM\": optimally weighted two-stage GMM with simplified covariance matrix.The models are estimated with the Gauss–Newton algorithm. The first stage of the two-stage specifications is estimated with the average outer product of the instrument vector as weight matrix.Additional methods are available for convenience:method = \"Poisson\": Poisson regression of the outcome on the treatment and controls;\nmethod = \"Reduced form\": Poisson regression of the outcome on the instruments and controls.The Microdata must contain: response, treatment, control and instrument. The Microdata may contain an offset. The outcome must be weakly positive. IVPoisson and Mullahy are subtypes of GMM."
},

{
    "location": "estimators.html#Reweighting-methods-1",
    "page": "Estimators",
    "title": "Reweighting methods",
    "category": "section",
    "text": "note: Note\nAll reweighting models require the specification of a first stage. They come in two flavors. In the first, you specify the first-stage model. In the second, you pass a previously fitted model. The latter is more verbose, but it allows you to customize and reuse the first stage.fit(IPW, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)\r\nfit(IPW, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)IPW estimates average treatment effects by inverse probability weighting. In a first stage, we use model M₁ to forecast the conditional probability of treatment take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted least squares. The intercept gives the mean outcome of the untreated. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The Microdata must contain: response, treatment and control. The treatment must be binary. IPW is a subtype of TwoStageModel.fit(Abadie, M₂::Type{ParModel}, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0, kwargs...)\r\nfit(Abadie, M₂::Type{ParModel}, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0 kwargs...)Abadie estimates local average response functions according to Abadie (2003). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we fit M₂ with the weights from the first stage. Keywords customize the second-stage estimator. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The Microdata must contain: response, treatment, control and instrument. The treatment and the instrument must be binary. Abadie is a subtype of TwoStageModel.fit(FrölichMelly, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)\r\nfit(FrölichMelly, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)This model estimates unconditional local average effects according to Frölich and Melly (2013). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted least squares. The intercept gives mean outcome of untreated compliers. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The Microdata must contain: response, treatment, control and instrument. The treatment and the instrument must be binary. FrölichMelly is a subtype of TwoStageModel.fit(Tan, M₁::Type{Micromodel}, MD::Microdata; trim::AbstractFloat = 0.0)\r\nfit(Tan, M₁::Micromodel, MD::Microdata; trim::AbstractFloat = 0.0)This model estimates the unconditional local average treatment effects according to Tan (2006). In a first stage, we use model M₁ to forecast the conditional probability of instrument take-up and construct estimation weights. In the second stage, we regress the outcome on the treatment and an intercept by weighted two-stage least squares. The intercept gives mean outcome of untreated compliers. We ignore observations whose score is below trim or above 1 - trim (see Crump et al. (2009)).The Microdata must contain: response, treatment, control and instrument. The treatment and the instrument must be binary. Tan is a subtype of TwoStageModel."
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
    "text": "This package supports the methods for regression models of StatsBase.jl.The following functions are available for all models: nobs and response.The following functions are available for parametric models: dof, dof_residual, coef, stderr, vcov, confint and coefnames. Note that all methods refer to the second stage of two-stage models.The following functions are available for maximum-likelihood estimators: deviance, nulldeviance, loglikelihood, nullloglikelihood, r2 and adjr2. There are also R² methods for OLS and IV. The following functions are available from StatsBase.jl: aic, aicc and bic.Most models support predict and fitted. predict estimates the index of single-index models. fitted computes the conditional outcome expectation. For example, predict estimates the Xβ of a logit model and fitted computes logistic(Xβ). Support for residuals depends on the availability of fitted. Out-of-sample forecast is supported."
},

{
    "location": "hypothesis_tests.html#",
    "page": "Hypothesis tests",
    "title": "Hypothesis tests",
    "category": "page",
    "text": ""
},

{
    "location": "hypothesis_tests.html#Hypothesis-tests-1",
    "page": "Hypothesis tests",
    "title": "Hypothesis tests",
    "category": "section",
    "text": ""
},

{
    "location": "hypothesis_tests.html#Significance-tests-1",
    "page": "Hypothesis tests",
    "title": "Significance tests",
    "category": "section",
    "text": "tstat(model::Union{GMM, ParModel, TwoStageModel})This function returns the t-statistic (i.e. the ratio of coefficients to standard error).pval(model::Union{GMM, ParModel, TwoStageModel})This function returns the p-value of a two-sided significance test."
},

{
    "location": "hypothesis_tests.html#Hausman-test-1",
    "page": "Hypothesis tests",
    "title": "Hausman test",
    "category": "section",
    "text": "This procedure tests the difference in coefficients between two parametric models (Hausman, 1984). Our implementation is based on the GMM representation of the joint estimation problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)). The optional argument names specifies the coefficients of interest as they appear on regression tables (be careful with categorical variables!). The output is a ParObject, which contains the vector of differences, their covariance matrix and labels.hausman_1s(\r\n        model₁::Union{GMM, ParModel, TwoStageModel},\r\n        model₂::Union{GMM, ParModel, TwoStageModel},\r\n        names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))\r\n    )This function is appropriate when model₁ and model₂ were based on a single estimation sample.hausman_2s(\r\n        model₁::Union{GMM, ParModel, TwoStageModel},\r\n        model₂::Union{GMM, ParModel, TwoStageModel},\r\n        names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))\r\n    )This function is appropriate when model₁ and model₂ were based on independent samples. For example, the samples might consist of independent observations with no overlap.hausman_2s(\r\n        model₁::Union{GMM, ParModel, TwoStageModel},\r\n        model₂::Union{GMM, ParModel, TwoStageModel},\r\n        corr::CorrStructure,\r\n        names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))\r\n    )This function is appropriate when model₁ and model₂ were based on dependent samples. For example, the samples might consist of independent observations with some overlap or clustered observations with shared clusters. The correlation structure corr must specify the correlation between all observations of both estimation samples. For example, you could construct corr for the entire dataset and construct the samples via the subset keyword to Microdata."
},

{
    "location": "estimation_tables.html#",
    "page": "Estimation tables",
    "title": "Estimation tables",
    "category": "page",
    "text": ""
},

{
    "location": "estimation_tables.html#Estimation-tables-1",
    "page": "Estimation tables",
    "title": "Estimation tables",
    "category": "section",
    "text": ""
},

{
    "location": "estimation_tables.html#Single-model-1",
    "page": "Estimation tables",
    "title": "Single model",
    "category": "section",
    "text": "coeftable(\r\n        model::Union{GMM, ParModel, TwoStageModel};\r\n        digits::Int = 4,\r\n        level::AbstractFloat = 0.95\r\n    )This function displays coefficients, standard errors and other information. The keyword digits controls rounding. The keyword level controls the level of the confidence interval (0.0 for none)."
},

{
    "location": "estimation_tables.html#Multiple-models-1",
    "page": "Estimation tables",
    "title": "Multiple models",
    "category": "section",
    "text": "etable(\r\n        args...;\r\n        digits::Int = 4,\r\n        aux::Union{Function, Nothing} = nothing,\r\n        stars::Matrix{Any} = [0.1 \"*\"; 0.05 \"**\"; 0.01 \"***\"],\r\n        titles::Vector{String} = [\"\"]\r\n    )This function displays a simple regression table. The keyword arguments are:digits: the number of digits on display.\naux: an auxiliary statistic (e.g., stderr), displayed below each coefficient.\nstars: the star scheme.\ntitles: the title of each regression (defaults to numbers)."
},

{
    "location": "bootstrapping.html#",
    "page": "Bootstrapping",
    "title": "Bootstrapping",
    "category": "page",
    "text": ""
},

{
    "location": "bootstrapping.html#Bootstrapping-1",
    "page": "Bootstrapping",
    "title": "Bootstrapping",
    "category": "section",
    "text": "This package does not provide support for bootstrap standard errors at the moment. Nonetheless, it is possible to bootstrap with the existing tools. This tutorial provides some sample code.We first load some packages:using StatsBase\r\nusing DataFrames\r\nusing CSV\r\nusing MicroeconometricsWe then set up the problem:S        = CSV.read(joinpath(datadir, \"auto.csv\")) ;\r\nS[:gpmw] = ((1.0 ./ S[:mpg]) ./ S[:weight]) * 100000 ;\r\nM        = Dict(:response => \"gpmw\", :control => \"foreign + 1\") ;\r\nD        = Microdata(S, M) ;Next, we obtain the coefficient estimates:E = fit(OLS, D, novar = true) ;We can now set up the bootstrap:srand(0101)\r\n\r\nreps = 1000 ;\r\nn    = nobs(E) ;\r\nwgts = fill(0, n) ;\r\nB    = Array{Float64}(reps, dof(E)) ;The vector wgts will translate the draw of a bootstrap sample into an input for Microdata. The matrix B will contain the sample of coefficient estimates. Don\'t forget to set the seed for the sake of reproducibility!The algorithm is:for b = 1:reps\r\n\r\n    wgts .= 0\r\n    draw  = rand(1:n, n)\r\n\r\n    for d in draw\r\n        wgts[d] += 1\r\n    end\r\n\r\n    Db      = Microdata(S, M, weights = fweights(wgts))\r\n    Eb      = fit(OLS, Db, novar = true)\r\n    B[b, :] = coef(Eb)\'\r\nendNote that we do not compute the covariance matrix at each step, which saves us some time.We can finally see the results:E.V = cov(B) ;\r\ncoeftable(E)The output is:                   Estimate  St. Err.   t-stat.   p-value      C.I. (95%)\r\nforeign: Foreign     0.2462    0.0682    3.6072    0.0003    0.1124  0.3799\r\n(Intercept)           1.609    0.0237   67.9372    <1e-99    1.5626  1.6554You can easily adapt this code to more complex problems (e.g., critical values) or parallelize it for additional speed!"
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
    "text": "This section helps you specify your own estimator. You can then submit a pull request to add it to Microeconometrics.jl!We will analyze the implementation of OLS. Although it is a simple model, others follow the same steps.The first step is defining the output struct:mutable struct OLS <: ParModel\r\n\r\n    sample::Microdata     # estimation sample\r\n    β::Vector{Float64}    # coefficient vector\r\n    V::Matrix{Float64}    # variance matrix\r\n\r\n    OLS() = new()\r\nendEvery estimator has these fields. Internal utilities rely on them. Some have additional fields (e.g., IV stores the estimation method and the weight matrix). Two-stage models store the estimators for each stage instead.We then define an uninitialized constructor:function OLS(MD::Microdata)\r\n    obj        = OLS()\r\n    obj.sample = MD\r\n    return obj\r\nendThe functions _fit! and _vcov! will later set β and V.The next step is overloading the function fit from StatsBase.jl:function fit(::OLS, MD::Microdata; novar::Bool = false)\r\n\r\n    obj = OLS(MD)\r\n\r\n    _fit!(obj, getweights(obj))\r\n    novar || _vcov!(obj, getcorr(obj), getweights(obj))\r\n\r\n    return obj\r\nendFor OLS, we only need to initialize the output object and pass it to _fit! and _vcov!. (It is not actually necessary to extend fit unless you need to perform additional steps before the estimation, as the fallback will suffice.) Note the utilities getcorr and getweights.We can now estimate the coefficients. For efficiency, we write separate functions for unweighted and weighted data.function _fit!(obj::OLS, w::UnitWeights)\r\n    y     = getvector(obj, :response)\r\n    x     = getmatrix(obj, :control)\r\n    obj.β = x \\ y\r\nend\r\n\r\nfunction _fit!(obj::OLS, w::AbstractWeights)\r\n    y     = getvector(obj, :response)\r\n    x     = getmatrix(obj, :control)\r\n    v     = Diagonal(w) * x\r\n    obj.β =  (v\' * x) \\ (v\' * y)\r\nendNotice the internal utilities getvector and getmatrix. Their first argument is a Microdata or a model structure. The following arguments are the model components of interest. You can request several components from getmatrix at once. For example, IV needs x = getmatrix(obj, :treatment, :control) and z = getmatrix(obj, :instrument, :control). A single matrix is returned in all cases.warning: Warning\ngetvector and getmatrix return views into the underlying data matrix. You should never modify their output, as you would irremediably alter the data. If you need to perform an in-place operation, make a copy beforehand.OLS does not require nonlinear optimization. If your estimator needs it, you can use the tools of Optim.jl. See the implementation of Logit for an example.We must now define score and jacobian. These functions are the building blocks of the variance estimator. The score is the vector of moment conditions. For OLS, it is −xᵢ (yᵢ − xᵢ\'β) (the derivative of the objective function). score should return the matrix of score vectors in row form. The Jacobian matrix is the derivative of the moment conditions. For OLS, it is xᵢ xᵢ\'. jacobian should return the weighted sum of Jacobians (i.e., the expected Jacobian × the number of observations).function score(obj::OLS)\r\n    x = copy(getmatrix(obj, :control))\r\n    û = residuals(obj)\r\n    return - x .* û\r\nend\r\n\r\nfunction jacobian(obj::OLS, w::UnitWeights)\r\n    x = getmatrix(obj, :control)\r\n    return x\' * x\r\nend\r\n\r\nfunction jacobian(obj::OLS, w::AbstractWeights)\r\n    x = getmatrix(obj, :control)\r\n    v = copy(x) .* w\r\n    return x\' * v\r\nendscore returns the score for each observation, so it ignores weights. jacobian returns an expectation; therefore, it must account for weights.We do not need to extend _vcov!. The default method will call score and jacobian and construct the appropriate estimator, accounting for the correlation structure of the data and the type of weights.We now overload predict and fitted. For OLS, these functions are equivalent.predict(obj::OLS) = getmatrix(obj, :control) * obj.β\r\nfitted(obj::OLS)  = predict(obj)The next step is optional. We extend jacobexp, which computes the derivative of fitted values.jacobexp(obj::OLS) = copy(getmatrix(obj, :control))jacobexp is only necessary when the estimator serves as the first stage of a two-stage estimator. By extending it, you make your estimator available to two-stage estimators.We conclude with a function to retrieve coefficient labels:coefnames(obj::OLS) = getnames(obj, :control)The syntax of getnames is similar to that of getmatrix.You can implement additional methods. For example, Microeconometrics.jl extends r2 and adjr2 to OLS (the fallback method only cover subtypes of MLE)."
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
    "text": "More models;\nMarginal effects for nonlinear models;\nA framework for panel data;\nA framework for semiparametric and nonparametric models;\nWald test;\nLikelihood-based tests."
},

]}
