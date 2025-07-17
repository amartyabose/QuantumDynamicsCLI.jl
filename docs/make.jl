using Documenter
using DocumenterCitations
using QuantumDynamicsCLI

bib = CitationBibliography(joinpath(@__DIR__, "src", "library.bib"))
makedocs(;
    plugins=[bib],
    modules=[QuantumDynamicsCLI],
    sitename="QuantumDynamicsCLI.jl",
    pages=[
        "Introduction" => "index.md",
        "Documentation" => [
            "qdsim Inputs" => "./documentation/ParseInput.md",
            "Simulate Module" => "./documentation/Simulate.md",
            "Post Module" => "./documentation/Post.md",
            "Comonicon" => "./documentation/ComonIcon.md",
        ],
        "Examples" => [
            "Spin-Boson Simulations" => "./documentation/Spin-Boson.md",
        ],
        "References" => "./documentation/References.md"
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamicsCLI.jl.git"
)