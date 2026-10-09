using Documenter
using ModVisc

makedocs(
    sitename="ModVisc.jl",
    remotes=nothing,
    format=Documenter.HTML(
        edit_link=nothing,
        repolink="https://github.com/albert-de-montserrat/ModVisc.jl",
    ),
    modules=[ModVisc],
    pages=[
        "Home" => "index.md",
        "Mathematical model" => "model.md",
        "Paper figures" => "paper_figures.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(
    repo="github.com/albert-de-montserrat/ModVisc.jl.git",
    devbranch="main",
)
