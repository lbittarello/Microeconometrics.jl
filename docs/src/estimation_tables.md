# Estimation tables

## Single model

```julia
coeftable(
    model::Union{GMM, ParModel, TwoStageModel};
    digits::Int = 4,
    level::AbstractFloat = 0.95
)
```

This function displays coefficients, standard errors and other information. The keyword `digits` controls rounding. The keyword `level` controls the level of the confidence interval (0.0 for none).

## Multiple models

```julia
etable(
    args...;
    digits::Int = 4,
    aux::Union{Function, Nothing} = nothing,
    stars::Matrix{Any} = [0.1 "*"; 0.05 "**"; 0.01 "***"],
    titles::Vector{String} = [""]
)
```

This function displays a simple regression table. The keyword arguments are:

- `digits`: the number of digits on display.
- `aux`: an auxiliary statistic (e.g., stderr), displayed below each coefficient.
- `stars`: the star scheme.
- `titles`: the title of each regression (defaults to numbers).
