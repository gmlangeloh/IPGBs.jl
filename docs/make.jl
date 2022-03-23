using Documenter
using IPGBs

makedocs(
    sitename = "IPGBs",
    format = Documenter.HTML(),
    modules = [IPGBs],
    pages = [
        "index.md",
        "IPGBs" => "IPGBs.md",
        "BinomialSets" => "BinomialSets.md",
        "Binomials" => "Binomials.md",
        "Buchberger" => "Buchberger.md",
        "FGLM" => "FGLM.md",
        "FastBitSets" => "FastBitSets.md",
        "FastComparator" => "FastComparator.md",
        "FourTi2" => "FourTi2.md",
        "GBAlgorithms" => "GBAlgorithms.md",
        "GBElements" => "GBElements.md",
        "GBTools" => "GBTools.md",
        "GradedBinomials" => "GradedBinomials.md",
        "IPInstances" => "IPInstances.md",
        "Markov" => "Markov.md",
        "MonomialHeaps" => "MonomialHeaps.md",
        "Orders" => "Orders.md",
        "SignatureAlgorithms" => "SignatureAlgorithms.md",
        "SignaturePolynomials" => "SignaturePolynomials.md",
        "SolverTools" => "SolverTools.md",
        "StandardDecomposition" => "StandardDecomposition.md",
        "Statistics" => "Statistics.md",
        "SupportTrees" => "SupportTrees.md",
        "TriangleHeaps" => "TriangleHeaps.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
