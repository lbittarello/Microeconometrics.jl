# Microeconometrics.jl

This package provides support for microeconometric estimation.
It supports complex covariance structures (clustered, etc.) and weighted data.
Please report bugs by
[opening an issue](https://github.com/lbittarello/Microeconometrics.jl/issues/new).
Information on specific versions can be found on the
[release page](https://github.com/lbittarello/Microeconometrics.jl/releases).

## Supported estimators

More models are planned.
If your preferred model is not currently available,
[file an issue](https://github.com/lbittarello/Microeconometrics.jl/issues/new) or
[contribute](#contributing)!

- [Models for exogenous regressors](estimators.md#Models-for-exogenous-regressors-1)
    - Ordinary least squares
    - Logit
    - Probit
    - Inverse probability weighting
- [Models for endogenous regressors](estimators.md#Models-for-endogenous-regressors-1)
    - Two-stage least squares
    - Abadie (2002)
    - Frölich and Melly (2013)
    - Tan (2006)

## Package manual

```@contents
Pages = [
        "getting_started.md",
        "model_specification.md",
        "estimators.md",
        "methods.md",
        "contributing.md",
        "to_do.md",
    ]
Depth = 2
```
