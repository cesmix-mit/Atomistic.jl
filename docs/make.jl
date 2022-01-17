push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Atomistic
using Documenter
using Literate

DocMeta.setdocmeta!(Atomistic, :DocTestSetup, :(using Atomistic); recursive = true)

Literate.markdown(joinpath(@__DIR__, "src", "usage.jl"), joinpath(@__DIR__, "src"); documenter = true)

makedocs(;
    modules = [Atomistic],
    authors = "CESMIX-MIT",
    repo = "https://github.com/cesmix-mit/Atomistic.jl/blob/{commit}{path}#{line}",
    sitename = "Atomistic.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://cesmix-mit.github.io/Atomistic.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Implementing the Atomistic API" => "extension.md",
        "Using Atomistic-Compatible Packages" => "usage.md",
        "API Reference" => "api.md",
    ],
    doctest = true,
    linkcheck = true,
    strict = true
)

deploydocs(;
    repo = "github.com/cesmix-mit/Atomistic.jl.git",
    devbranch = "main",
    push_preview = true
)