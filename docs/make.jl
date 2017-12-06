using Documenter
using Microeconometrics

makedocs(
    modules  = [Microeconometrics],
    doctest  = false,
    clean    = false,
    sitename = "Microeconometrics.jl",
    format   = :html,
    pages    = Any[
        "Introduction"        => "index.md",
        "Getting started"     => "getting_started.md",
        "Model specification" => "model_specification.md",
        "Estimators"          => "estimators.md",
        "Methods"             => "methods.md",
        "Contributing"        => "contributing.md",
        "To do"               => "to_do.md",
    ]
)

deploydocs(
    repo   = "github.com/lbittarello/Microeconometrics.jl.git",
    target = "build",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
