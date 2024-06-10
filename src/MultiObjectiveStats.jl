module MultiObjectiveStats

export Stats, new_gb_size_and_time

"""
Struct storing the data of a MultiObjectiveAlgorithms.moip_gb_solve run to analyze
and include as a table in papers.
"""
mutable struct Stats
    gb_total_time::AbstractFloat
    normalform_time::AbstractFloat
    total_time::AbstractFloat
    solver_time::AbstractFloat
    initial_time::AbstractFloat
    pareto_size::Int
    gb_times::Vector{AbstractFloat}
    gb_sizes::Vector{Int}
    num_ips::Int

    Stats() = new(0.0, 0.0, 0.0, 0.0, time(), 0, Float64[], Int[], 0)
end

function header(num_objectives :: Int)
    print("pareto_size total_time solver_time num_ips gb_total_time ")
    for i in 1:num_objectives
        print("gb_sizes$i gb_times$i ")
    end
    println("normalform_time")
end

function terminate(stats::Stats, pareto)
    stats.total_time = time() - stats.initial_time
    stats.pareto_size = length(pareto)
end

function new_gb_size_and_time(stats:: Stats, size :: Int, time :: Float64)
    push!(stats.gb_times, time)
    push!(stats.gb_sizes, size)
end

"""
Shows a `MultiObjecitveStats` struct with all data in the same row to easily
make tables for papers.
"""
function Base.show(
    io::IO,
    stats::Stats
)
    print(io, stats.pareto_size, " ", stats.total_time, " ", stats.solver_time, " ")
    print(io, stats.num_ips, " ", stats.gb_total_time, " ")
    for i in eachindex(stats.gb_times)
        print(io, stats.gb_sizes[i], " ", stats.gb_times[i], " ")
    end
    print(io, stats.normalform_time)
end

end
