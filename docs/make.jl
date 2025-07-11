cd(@__DIR__)
using SignalDecomposition
using PredefinedDynamicalSystems

pages = [
    "index.md",
    "examples.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

build_docs_with_style(pages, SignalDecomposition;
    authors = "George Datseris <datseris.george@gmail.com>",
    expandfirst = ["index.md"], #  this is the first script that loads colorscheme
)
