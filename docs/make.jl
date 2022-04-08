using MuTuner
using Documenter

DocMeta.setdocmeta!(MuTuner, :DocTestSetup, :(using MuTuner); recursive=true)

makedocs(;
    modules=[MuTuner],
    authors="Benjamin Cohen-Stead <benwcs@gmail.com>",
    repo="https://github.com/cohensbw/MuTuner.jl/blob/{commit}{path}#{line}",
    sitename="MuTuner.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cohensbw.github.io/MuTuner.jl",
        assets=String[],
    ),
    pages=[
        "Home"     => "index.md",
        "Examples" => "examples.md",
        "API"      => "api.md"
    ],
)

deploydocs(;
    repo="github.com/cohensbw/MuTuner.jl",
    devbranch="master",
)
