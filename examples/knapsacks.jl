# In this example, we explore GBs of unbounded / binary knapsacks
# We also generate some figures of their digraphs to illustrate the results
# The reason this is not a notebook is that plotting interacts weirdly with VSCode
# notebooks in Julia. It appears to be a bug in Jupyter or something
# https://github.com/JuliaPlots/Plots.jl/issues/4827

using IPGBs.FeasibleGraphs
using MIPMatrixTools.IPInstances

function basic_knapsack(binary = false)
    A = [15 18 19 17 20]
    b = [45]
    C = [1 1 1 1 1]
    u = [binary ? 1 : nothing for i in 1:5]
    return IPInstance(A, b, C, u)
end

plot_feasible_graph(basic_knapsack(true))
