# For maintainers

This page is for people changing the package or the manual.

## Building the manual

From the repository root:

```bash
julia -O1 --project=docs docs/generate.jl
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

`docs/generate.jl` turns the Literate sources in `docs/` into Markdown under
`docs/src/generated/`. `docs/make.jl` renders that Markdown. After editing a
tutorial `.jl` file, run `docs/generate.jl` and commit the Markdown in
`docs/src/generated/`.

CI runs only `docs/make.jl`. It does not regenerate the Literate pages.

`-O1` avoids an LLVM optimizer crash seen with Julia 1.11.5 on Apple Silicon
while `docs/generate.jl` executes the four-body mass-window sample.

## What each build step runs

The three-body and B-decay tutorials become `@example` blocks. Documenter runs
those blocks during `docs/make.jl`.

`docs/four-body-mass-window.jl` is executed by `docs/generate.jl`, and the
resulting [page](generated/four-body-mass-window.md) is committed. Its code
fences are plain Julia, so `docs/make.jl` renders the recorded sample and does
not draw the events again. Lines between `@hide-start` and `@hide-end` still run; they
are omitted from the page. That region holds the angle cross-check and the
text of the result table.
