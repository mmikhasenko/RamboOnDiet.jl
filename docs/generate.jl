import Pkg
Pkg.instantiate()

using Literate

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

# Lines between `@hide-start` and `@hide-end` are executed but omitted from the
# page, so the rendered summary is not preceded by the string that builds it.
function hide_marked_region(content)
    hiding = false
    lines = String[]
    for line in eachsplit(content, '\n', keepempty = true)
        if occursin("@hide-start", line)
            hiding = true
            continue
        elseif occursin("@hide-end", line)
            hiding = false
            continue
        end
        push!(lines, hiding ? line * " #hide" : line)
    end
    return join(lines, '\n')
end

# Executed locally. The committed Markdown uses plain Julia fences, so
# `docs/make.jl` renders the page without resampling.
Literate.markdown(
    joinpath(DOCS_ROOT, "four-body-mass-window.jl"),
    LITERATE_OUTPUT;
    execute = true,
    preprocess = hide_marked_region,
)
