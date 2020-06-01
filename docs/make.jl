using PoissonFE
using Documenter

makedocs(;
    modules=[PoissonFE],
    authors="Evan Wright <evan.btbucket@gmail.com> and contributors",
    repo="https://github.com/ew-git/PoissonFE.jl/blob/{commit}{path}#L{line}",
    sitename="PoissonFE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ew-git.github.io/PoissonFE.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ew-git/PoissonFE.jl",
)
