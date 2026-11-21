# Audit public docstrings and the README examples.
#
# From the package root:
#   julia --project scripts/audit_documentation.jl

import Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT; io = devnull)

using RamboOnDiet
import FourVectors
import LorentzVectorBase
import LinearAlgebra
import Random
import StaticArrays

# `?RamboOnDiet` shows README.md. The help-page graph starts at the generator.
const ENTRY_POINT = :PhaseSpaceGenerator

const DEPENDENCIES =
    (FourVectors, LorentzVectorBase, StaticArrays, LinearAlgebra, Random, Base)

const REF_RE = r"\[`([^`]+)`\]\(@ref(?:\s[^)]*)?\)"
const JULIA_FENCE_RE = r"^[ \t]*```julia[ \t]*\r?\n(.*?)^[ \t]*```"sm

struct Link
    raw::String
    name::Symbol
    qualified::Bool
    module_name::Union{Nothing,Symbol}
end

function docstr_text(docstr)
    text = docstr.text
    if text isa AbstractString
        return String(text)
    end
    return join(string.(text))
end

function raw_doc_text(mod::Module, name::Symbol)
    binding = Base.Docs.Binding(mod, name)
    meta = Base.Docs.meta(mod)
    haskey(meta, binding) || return ""
    return join((docstr_text(docstr) for docstr in values(meta[binding].docs)), "\n\n")
end

function dedent(code::AbstractString)
    lines = collect(eachline(IOBuffer(String(code))))
    while !isempty(lines) && isempty(strip(lines[1]))
        popfirst!(lines)
    end
    while !isempty(lines) && isempty(strip(lines[end]))
        pop!(lines)
    end
    isempty(lines) && return ""
    width = minimum(lines) do line
        isempty(strip(line)) && return typemax(Int)
        m = match(r"^[ \t]*", line)
        return length(m.match)
    end
    width == typemax(Int) && return ""
    stripped = map(lines) do line
        return length(line) ≤ width ? "" : line[(width+1):end]
    end
    return join(stripped, "\n")
end

function julia_blocks(text::AbstractString)
    return [dedent(m.captures[1]) for m in eachmatch(JULIA_FENCE_RE, text)]
end

function parse_link(raw::AbstractString)
    parts = split(strip(raw), '.')
    if length(parts) == 1
        return Link(String(raw), Symbol(parts[1]), false, nothing)
    end
    return Link(String(raw), Symbol(parts[end]), true, Symbol(parts[1]))
end

function links_in(text::AbstractString)
    return [parse_link(m.captures[1]) for m in eachmatch(REF_RE, text)]
end

function exported_names()
    # The module binding is not a callable help page.
    return sort!(filter(n -> n !== nameof(RamboOnDiet), collect(names(RamboOnDiet))))
end

function is_local_package_link(link::Link)
    link.module_name === nothing && return true
    return link.module_name === nameof(RamboOnDiet) || link.module_name === :RamboOnDiet
end

function has_docs(mod::Module, name::Symbol)
    !isempty(strip(raw_doc_text(mod, name))) && return true
    isdefined(mod, name) || return false
    owner = parentmodule(getfield(mod, name))
    owner === mod && return false
    return !isempty(strip(raw_doc_text(owner, name)))
end

# A dependency export counts when that package is the call site and the
# function carries a docstring, including one written in the module that defines it.
function documented_export(mod::Module, name::Symbol)
    name in names(mod) || return false
    return has_docs(mod, name)
end

function foreign_home(name::Symbol)
    for mod in DEPENDENCIES
        documented_export(mod, name) || continue
        return mod
    end
    return nothing
end

function example_calls(examples, name::Symbol)
    pattern = Regex("\\b$(name)\\s*\\(")
    return any(example -> occursin(pattern, example), examples)
end

function execute_block(code::AbstractString, label::AbstractString)
    if occursin(r"\bPkg\.add\b", code)
        execute_install_example(code, label)
        return
    end
    mod = Module()
    include_string(mod, code, label)
    return nothing
end

function execute_install_example(code::AbstractString, label::AbstractString)
    mktempdir() do dir
        previous = Base.active_project()
        Pkg.activate(dir; io = devnull)
        try
            withenv("CI" => "true", "JULIA_PKG_PRECOMPILE_AUTO" => "0") do
                mod = Module()
                include_string(mod, code, label)
            end
        finally
            Pkg.activate(previous; io = devnull)
        end
    end
    return nothing
end

function related_body(text::AbstractString)
    matches = collect(eachmatch(r"(?m)^## Related\s*$", text))
    length(matches) == 1 || return nothing
    tail = text[matches[1].offset:end]
    occursin(r"(?m)^## (?!Related\b)", tail[length(matches[1].match):end]) && return nothing
    return tail
end

function mermaid_lines(edges)
    lines = String["flowchart LR"]
    for source in sort!(collect(keys(edges)))
        for target in sort!(collect(edges[source]))
            push!(lines, "    $(source) --> $(target)")
        end
    end
    return lines
end

function counts_table(exports, edges)
    incoming = Dict(name => 0 for name in exports)
    for source in exports
        for target in edges[source]
            incoming[target] += 1
        end
    end
    lines = String["| Name | Incoming | Outgoing |", "| --- | ---: | ---: |"]
    for name in exports
        label = name === ENTRY_POINT ? "$(name) (entry)" : string(name)
        push!(lines, "| $(label) | $(incoming[name]) | $(length(edges[name])) |")
    end
    return lines, incoming
end

function reachable_from(entry, edges)
    seen = Set{Symbol}()
    stack = Symbol[entry]
    while !isempty(stack)
        node = pop!(stack)
        node in seen && continue
        push!(seen, node)
        for target in get(edges, node, ())
            push!(stack, target)
        end
    end
    return seen
end

function graph_markdown(exports, edges)
    table, _ = counts_table(exports, edges)
    mermaid = mermaid_lines(edges)
    return join(
        [
            "# Discoverability graph",
            "",
            "Maintainer map of help links between exported names. An arrow points from a help page to the page a reader opens next. The entry point is `$(ENTRY_POINT)`. Any other name with no incoming link is hard to find. A name with no outgoing link is a dead end. Every exported name needs an outgoing link and a path from the entry point. This file is not part of the user guide.",
            "",
            "```mermaid",
            mermaid...,
            "```",
            "",
            "## Link counts",
            "",
            table...,
            "",
        ],
        "\n",
    )
end

function write_graph(exports, edges)
    path = joinpath(pkgdir(RamboOnDiet), "docs", "discoverability_graph.md")
    write(path, graph_markdown(exports, edges))
    return path
end

function print_graph(exports, edges)
    table, _ = counts_table(exports, edges)
    println("Discoverability graph (entry: $(ENTRY_POINT))")
    println()
    println("```mermaid")
    foreach(println, mermaid_lines(edges))
    println("```")
    println()
    foreach(println, table)
    println()
end

function audit()
    failures = String[]
    exports = exported_names()
    if !(ENTRY_POINT in exports)
        push!(failures, "entry point $(ENTRY_POINT) is not exported")
    end

    docs = Dict{Symbol,String}()
    examples = Dict{Symbol,Vector{String}}()
    edges = Dict{Symbol,Set{Symbol}}(name => Set{Symbol}() for name in exports)

    for name in exports
        text = raw_doc_text(RamboOnDiet, name)
        docs[name] = text
        if isempty(strip(text))
            push!(failures, "$(name) has no docstring")
            examples[name] = String[]
            continue
        end
        blocks = julia_blocks(text)
        examples[name] = blocks
        if isempty(blocks)
            push!(failures, "$(name) has no ```julia example")
        end
        related = related_body(text)
        if related === nothing
            push!(
                failures,
                "$(name) needs one ## Related section at the end of the docstring",
            )
        elseif isempty(links_in(related))
            push!(failures, "$(name) ## Related section has no [`name`](@ref) link")
        end

        for link in links_in(text)
            if is_local_package_link(link)
                if isdefined(RamboOnDiet, link.name)
                    if !has_docs(RamboOnDiet, link.name)
                        push!(
                            failures,
                            "$(name) links to [`$(link.raw)`](@ref), which has no docstring",
                        )
                    end
                else
                    home = foreign_home(link.name)
                    if home === nothing
                        push!(
                            failures,
                            "$(name) links to [`$(link.raw)`](@ref), which has no docstring in RamboOnDiet or its dependencies",
                        )
                    elseif !example_calls(blocks, link.name)
                        push!(
                            failures,
                            "$(name) links to [`$(link.name)`](@ref) in $(nameof(home)), but no example calls it",
                        )
                    end
                end
                if link.name in exports && link.name !== name
                    push!(edges[name], link.name)
                end
            else
                modname = link.module_name
                home = nothing
                for mod in DEPENDENCIES
                    nameof(mod) === modname || continue
                    home = mod
                    break
                end
                if home === nothing ||
                   !isdefined(home, link.name) ||
                   !has_docs(home, link.name)
                    push!(
                        failures,
                        "$(name) links to [`$(link.raw)`](@ref), which is not a documented name in that package",
                    )
                elseif !example_calls(blocks, link.name)
                    push!(
                        failures,
                        "$(name) links to [`$(link.raw)`](@ref), but no example calls $(link.name)",
                    )
                end
            end
        end
    end

    for name in exports
        for (index, block) in enumerate(examples[name])
            label = "docstring $(name) example $(index)"
            try
                execute_block(block, label)
            catch err
                push!(failures, "$(label) failed: $(sprint(showerror, err))")
            end
        end
    end

    readme_path = joinpath(pkgdir(RamboOnDiet), "README.md")
    readme = read(readme_path, String)
    for (index, block) in enumerate(julia_blocks(readme))
        label = "README.md example $(index)"
        try
            execute_block(block, label)
        catch err
            push!(failures, "$(label) failed: $(sprint(showerror, err))")
        end
    end

    for name in exports
        if isempty(edges[name])
            push!(failures, "$(name) has no outgoing link to another exported name")
        end
    end
    if ENTRY_POINT in exports
        seen = reachable_from(ENTRY_POINT, edges)
        for name in exports
            if name ∉ seen
                push!(failures, "$(name) cannot be reached from $(ENTRY_POINT)")
            end
        end
    end

    graph_path = write_graph(exports, edges)
    print_graph(exports, edges)

    if isempty(failures)
        println(
            "Audit passed for $(length(exports)) exports. Graph written to $(graph_path).",
        )
        return 0
    end
    println(stderr, "Audit failed ($(length(failures))):")
    for failure in failures
        println(stderr, "  ", failure)
    end
    return 1
end

exit(audit())
