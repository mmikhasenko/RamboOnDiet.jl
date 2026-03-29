import Pkg
Pkg.instantiate()

using Documenter
using Literate
using RemboOnDiet

const DOCS_ROOT = @__DIR__
const LITERATE_OUTPUT = joinpath(DOCS_ROOT, "src", "generated")

mkpath(LITERATE_OUTPUT)

Literate.markdown(
    joinpath(DOCS_ROOT, "three-body.jl"),
    LITERATE_OUTPUT;
    documenter = true,
    flavor = Literate.DocumenterFlavor(),
)

makedocs(
    sitename = "RemboOnDiet",
    modules = [RemboOnDiet],
    source = "src",
    build = "build",
    format = Documenter.HTML(),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Three-Body Tutorial" => "generated/three-body.md",
    ],
)
