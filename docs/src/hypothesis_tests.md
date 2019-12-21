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

```julia
hausman_test(
        model₁::Union{GMM, ParModel, TwoStageModel},
        model₂::Union{GMM, ParModel, TwoStageModel},
        names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
)
```

This procedure tests the difference in coefficients between two parametric models ([Hausman, 1984](http://jstor.org/stable/1913827)). Both models must have been estimated on the sample. It can be used in replacement of the Sobel test. This implementation is based on the GMM representation of the joint estimation problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)). Neither model need be efficient.

The optional argument `names` specifies the coefficients of interest as they appear on regression tables (be careful with categorical variables!). The output is a `ParEstimate`, which contains the vector of differences, their covariance matrix and labels.

## Chow test

This procedure tests the difference in coefficients between two parametric models ([Hausman, 1984](http://jstor.org/stable/1913827)). The models need not have been estimated on the sample. This implementation is based on the GMM representation of the joint estimation problem (see Subsection 8.3.2 of Cameron and Trivedi (2005)).

The optional argument `names` specifies the coefficients of interest as they appear on regression tables (be careful with categorical variables!). The output is a `ParEstimate`, which contains the vector of differences, their covariance matrix and labels.

```julia
 chow_test(
         model₁::Union{GMM, ParModel, TwoStageModel},
         model₂::Union{GMM, ParModel, TwoStageModel},
         names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
 )
```

This function is appropriate when `model₁` and `model₂` were based on independent samples. For example, the samples might consist of independent observations with no overlap.

```julia
 chow_test(
         model₁::Union{GMM, ParModel, TwoStageModel},
         model₂::Union{GMM, ParModel, TwoStageModel},
         corr::CorrStructure,
         names::Vector{String} = intersect(coefnames(obj₁), coefnames(obj₂))
 )
```

This function is appropriate when `model₁` and `model₂` were based on dependent samples. For example, two clustered samples might have one or more clusters in common. We must then take the correlation across samples within shared clusters into account.

The correlation structure `corr` must specify the correlation between all observations of both estimation samples. For instance, you could construct `corr` for the entire dataset and construct the subsamples via the `subset` keyword of `Microdata`.
