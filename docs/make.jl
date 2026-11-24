import Pkg
Pkg.instantiate()

using Documenter
using RamboOnDiet

makedocs(
    sitename = "RamboOnDiet",
    modules = [RamboOnDiet],
    source = "src",
    build = "build",
    format = Documenter.HTML(repolink = "https://github.com/mmikhasenko/RamboOnDiet.jl"),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
        "Algorithm" => "algorithm.md",
        "Implementation" => "implementation.md",
        "Edge cases" => "edge-cases.md",
        "Three-Body Tutorial" => "generated/three-body.md",
        "B Decay Tutorial" => "generated/b-decay.md",
        "Four-Body Mass Window" => "generated/four-body-mass-window.md",
        "For maintainers" => "maintainers.md",
    ],
)
