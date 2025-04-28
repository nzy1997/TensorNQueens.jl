using TensorNQueens
using Documenter

DocMeta.setdocmeta!(TensorNQueens, :DocTestSetup, :(using TensorNQueens); recursive=true)

makedocs(;
    modules=[TensorNQueens],
    authors="nzy1997",
    sitename="TensorNQueens.jl",
    format=Documenter.HTML(;
        canonical="https://nzy1997.github.io/TensorNQueens.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nzy1997/TensorNQueens.jl",
    devbranch="main",
)
