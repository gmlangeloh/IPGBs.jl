using Documenter
using IPGBs

makedocs(
  sitename = "IPGBs",
  authors = "Gabriel Mattos Langeloh",
  format = Documenter.HTML(
    # prettyurls is necessary to get direct links from each page
    prettyurls = get(ENV, "CI", nothing) == true,
    mathengine = Documenter.MathJax2(),
    collapselevel = 1,
  ),
  modules = [IPGBs],
  doctest = false,
  pages = [
    "index.md",
    "IPGBs" => "IPGBs.md",
    "BinomialSets" => "BinomialSets.md",
    "Binomials" => "Binomials.md",
    "Buchberger" => "Buchberger.md",
    "FastBitSets" => "FastBitSets.md",
    "FastComparator" => "FastComparator.md",
    "FeasibleGraphs" => "FeasibleGraphs.md",
    "FGLM" => "FGLM.md",
    "FourTi2" => "FourTi2.md",
    "GBAlgorithms" => "GBAlgorithms.md",
    "GBElements" => "GBElements.md",
    "GBTools" => "GBTools.md",
    "GradedBinomials" => "GradedBinomials.md",
    "IPInstances" => "IPInstances.md",
    "Markov" => "Markov.md",
    "MatrixTools" => "MatrixTools.md",
    "MonomialHeaps" => "MonomialHeaps.md",
    "MultiObjectiveAlgorithms" => "MultiObjectiveAlgorithms.md",
    "MultiObjectiveStats" => "MultiObjectiveStats.md",
    "MultiObjectiveTools" => "MultiObjectiveTools.md",
    "Optimize" => "Optimize.md",
    "Orders" => "Orders.md",
    "SignatureAlgorithms" => "SignatureAlgorithms.md",
    "SignaturePolynomials" => "SignaturePolynomials.md",
    "SingleObjective" => "SingleObjective.md",
    "SolverTools" => "SolverTools.md",
    "StandardDecomposition" => "StandardDecomposition.md",
    "Statistics" => "Statistics.md",
    "SupportMatrixTrees" => "SupportMatrixTrees.md",
    "SupportTrees" => "SupportTrees.md",
    "TriangleHeaps" => "TriangleHeaps.md",
    "Walkback" => "Walkback.md"
  ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
  repo = "github.com/gmlangeloh/IPGBs.jl.git"
)
