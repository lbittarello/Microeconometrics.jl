using Documenter
using Microeconometrics

makedocs(
    modules  = [Microeconometrics],
    sitename = "Microeconometrics.jl",
    format   = Documenter.HTML(),
    doctest  = false,
    pages    = [
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

deploydocs(repo = "github.com/lbittarello/Microeconometrics.jl.git")