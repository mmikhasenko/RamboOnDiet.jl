import Pkg
Pkg.instantiate()

using Documenter
using Literate
using RemboOnDiet

const DOCS_ROOT = @__DIR__
const LITERATE_OUTPUT = joinpath(DOCS_ROOT, "src", "generated")

mkpath(LITERATE_OUTPUT)

for tutorial in ("three-body.jl", "b-decay.jl")
    Literate.markdown(
        joinpath(DOCS_ROOT, tutorial),
        LITERATE_OUTPUT;
        documenter = true,
        flavor = Literate.DocumenterFlavor(),
    )
end

makedocs(
    sitename = "RemboOnDiet",
    modules = [RemboOnDiet],
    source = "src",
    build = "build",
    format = Documenter.HTML(repolink = "https://github.com/mmikhasenko/RamboOnDiet.jl"),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Algorithm" => "algorithm.md",
        "Implementation" => "implementation.md",
        "Three-Body Tutorial" => "generated/three-body.md",
        "B Decay Tutorial" => "generated/b-decay.md",
    ],
)
