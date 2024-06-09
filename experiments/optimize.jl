using IPGBs
using JuMP
using IPGBs.IPInstances

function test_optimize(path; init_sol = nothing)
    #model = read_from_file(path)
    #set_silent(model)
    #optimize!(model)
    #value = round(Int, objective_value(model))
    if isnothing(init_sol)
        ret, t, _, _, _ = @timed IPGBs.optimize(path, quiet=false)
    else
        ret, t, _, _, _ = @timed IPGBs.optimize(path, quiet=false, solution=init_sol)
    end
    #val = ret[2]
    #println("IPGBs opt: ", val, " time: ", t)
    #return ret[1]
end

function optimize_all()
    for dirname in readdir()
        if isdir(dirname)
            for filename in readdir(dirname, join=true)
                if endswith(filename, ".mps")
                    println("Running ", filename)
                    #Redirect stdout to a filename.log file
                    fn() = test_optimize(filename)
                    open(filename * ".optlog", "w") do f
                        redirect_stdout(fn, f)
                    end
                end
            end
        end
    end
end

optimize_all()
