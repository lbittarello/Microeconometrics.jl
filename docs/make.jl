using Documenter
using Microeconometrics

makedocs(
    modules  = [Microeconometrics],
    doctest  = false,
    clean    = false,
    sitename = "Microeconometrics.jl",
    format   = :html,
    pages    = Any[
        "Introduction"           => "index.md",
        "Getting started"        => "getting_started.md",
        "Model specification"    => "model_specification.md",
        "Correlation structures" => "correlation_structures.md",
        "Estimators"             => "estimators.md",
        "Methods"                => "methods.md",
        "Hypothesis tests"       => "hypothesis_tests.md",
        "Estimation tables"      => "estimation_tables.md",
        "Bootstrapping"          => "bootstrapping.md",
        "Contributing"           => "contributing.md",
        "To do"                  => "to_do.md",
    ]
)

deploydocs(
    repo   = "github.com/lbittarello/Microeconometrics.jl.git",
    target = "build",
    julia  = "0.7",
    deps   = nothing,
    make   = nothing,
)
