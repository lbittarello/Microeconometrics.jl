# Hypothesis tests

## Significance tests

```julia
tstat(model::Union{GMM, ParModel, TwoStageModel, ParEstimate})
```

This function returns the *t*-statistic (i.e. the ratio of coefficients to standard error).

```julia
pval(model::Union{GMM, ParModel, TwoStageModel, ParEstimate})
```

This function returns the *p*-value of a two-sided significance test.

## Hausman test

This procedure tests the difference in coefficients between two parametric models ([Hausman, 1984](http://jstor.org/stable/1913827)). It can be used in replacement of the Chow test and the Sobel test. Our implementation is based on the GMM representation of the joint estimation problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)). Neither model need be efficient.

The optional argument `names` specifies the coefficients of interest as they appear on regression tables (be careful with categorical variables!). The output is a `ParEstimate`, which contains the vector of differences, their covariance matrix and labels.

```julia
hausman_1s(
        model₁::Union{GMM, ParModel, TwoStageModel},
        model₂::Union{GMM, ParModel, TwoStageModel},
        names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
)
```

This function is appropriate when `model₁` and `model₂` were based on a single estimation sample.

```julia
 hausman_2s(
         model₁::Union{GMM, ParModel, TwoStageModel},
         model₂::Union{GMM, ParModel, TwoStageModel},
         names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
 )
```

This function is appropriate when `model₁` and `model₂` were based on independent samples. For example, the samples might consist of independent observations with no overlap.

```julia
 hausman_2s(
         model₁::Union{GMM, ParModel, TwoStageModel},
         model₂::Union{GMM, ParModel, TwoStageModel},
         corr::CorrStructure,
         names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
 )
```

This function is appropriate when `model₁` and `model₂` were based on dependent samples. For example, the samples might consist of independent observations with some overlap or clustered observations with shared clusters. The correlation structure `corr` must specify the correlation between all observations of both estimation samples. For example, you could construct `corr` for the entire dataset and construct the samples via the `subset` keyword to `Microdata`.
