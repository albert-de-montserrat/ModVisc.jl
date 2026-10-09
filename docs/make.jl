using Documenter
using ModVisc

makedocs(
    sitename="ModVisc.jl",
    remotes=nothing,
    format=Documenter.HTML(edit_link=nothing, repolink=nothing),
    modules=[ModVisc],
    pages=[
        "Home" => "index.md",
        "Mathematical model" => "model.md",
        "Mathematical validation" => "validation.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(
    repo="github.com/albert-de-montserrat/ModVisc.jl.git",
    devbranch="main",
)
