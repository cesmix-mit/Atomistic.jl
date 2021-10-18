push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using Atomistic
using Documenter

DocMeta.setdocmeta!(Atomistic, :DocTestSetup, :(using Atomistic); recursive=true)

makedocs(;
    modules=[Atomistic],
    authors="CESMIX-MIT",
    repo="https://github.com/cesmix-mit/Atomistic.jl/blob/{commit}{path}#{line}",
    sitename="Atomistic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cesmix-mit.github.io/Atomistic.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md"
    ],
    doctest=true,
    linkcheck=true,
    strict=true
)

deploydocs(;
    repo="github.com/cesmix-mit/Atomistic.jl.git",
    devbranch="main",
    push_preview=true,
)