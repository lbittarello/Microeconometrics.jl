# Microeconometrics.jl

| **Documentation** | **Master status** |
|:-----------------:|:-----------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] |

Microeconometric estimation in Julia.

Currently supported estimators:
- Linear regression
    - Ordinary least squares
    - Two-stage least squares
    - Linear GMM
- Binary choice
    - Logit
    - Probit
    - Complementary log-log
    - Gompit
- Count data
    - Poisson
    - IV Poisson with additive errors
    - IV Poisson with multiplicative errors
- Reweighting methods
    - Inverse probability weighting
    - Abadie (2003)
    - Frölich and Melly (2013)
    - Tan (2006)

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://lbittarello.github.io/Microeconometrics.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://lbittarello.github.io/Microeconometrics.jl/stable/

[travis-img]: https://travis-ci.org/lbittarello/Microeconometrics.jl.svg?branch=master
[travis-url]: https://travis-ci.org/lbittarello/Microeconometrics.jl
