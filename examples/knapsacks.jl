# In this example, we explore GBs of unbounded / binary knapsacks
# We can use them to generate some figures of their digraphs to illustrate the results
# The reason this is not a notebook is that plotting interacts weirdly with VSCode
# notebooks in Julia. It appears to be a bug in Jupyter or something
# https://github.com/JuliaPlots/Plots.jl/issues/4827

using IPGBs.IPInstances

function basic_knapsack(binary = false)
    A = [15 18 19 17]
    b = [34]
    C = [1 2 3 4]
    u = [binary ? 1 : nothing for i in 1:4]
    return IPInstance(A, b, C, u)
end
