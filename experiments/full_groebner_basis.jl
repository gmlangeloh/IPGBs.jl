#Run experiments
#Relevant data: running time, GB size, number of reductions, memory use
#Variations of the algorithm: with and without binary truncation

using IPGBs
using IPGBs.FourTi2

function run_instance(path)
    try
        println("IPGBs, no binary truncation")
        ipgbs_stats = @timed groebner_basis(path, quiet=false, use_binary_truncation=false)
        println("time => ", ipgbs_stats.time)
        println("memory => ", ipgbs_stats.bytes / 1024^3) # in GB
        println("IPGBs, with binary truncation")
        ipgbs_stats = @timed groebner_basis(path, quiet=false, use_binary_truncation=true)
        println("time => ", ipgbs_stats.time)
        println("memory => ", ipgbs_stats.bytes / 1024^3) # in GB
        fti2_stats = @timed groebner(path)
        println("4ti2 results")
        println("4ti2 time => ", fti2_stats.time) #memory isn't easily available
    catch e
        println("Error: ", e)
    end
end

function all_instances()
    #Precompile stuff to avoid timing issues
    run_instance("knapsack_binary/knapsack_binary_10_1.mps")

    #Now run the experiments for all instances
    for dirname in readdir()
        if isdir(dirname)
            for filename in readdir(dirname, join=true)
                if endswith(filename, ".mps")
                    println("Running ", filename)
                    #Redirect stdout to a filename.log file
                    fn() = run_instance(filename)
                    open(filename * ".log", "w") do f
                        redirect_stdout(fn, f)
                    end
                end
            end
        end
    end
end

all_instances()
