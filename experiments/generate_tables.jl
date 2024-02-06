using StatsKit

function get_line_value(line :: String)
    return parse(Float64, split(line, "=>")[2])
end

function instance_name(filename)
    #take filename, remove path and extension
    return split(split(filename, "/")[end], ".")[1]
end

function read_log(filename :: String, log :: Vector{String})
    no_trunc = false
    with_trunc = false
    fti2 = false
    no_trunc_time = 0.0
    no_trunc_mem = 0.0
    no_trunc_reduced_pairs = 0
    with_trunc_time = 0.0
    with_trunc_mem = 0.0
    with_trunc_reduced_pairs = 0
    fti2_time = 0.0
    #Parse log file
    for line in log
        #Figure out which algorithm we're dealing with
        if startswith(line, "IPGBs, no binary truncation")
            no_trunc = true
            with_trunc = false
            fti2 = false
            continue
        elseif startswith(line, "IPGBs, with binary truncation")
            no_trunc = false
            with_trunc = true
            fti2 = false
            continue
        elseif startswith(line, "4ti2 results")
            no_trunc = false
            with_trunc = false
            fti2 = true
            continue
        end
        #Get data from algorithm
        if startswith(line, "time =>")
            time = get_line_value(line)
            if no_trunc
                no_trunc_time += time
            elseif with_trunc
                with_trunc_time += time
            end
        elseif startswith(line, "memory =>")
            memory = get_line_value(line)
            if no_trunc
                no_trunc_mem += memory
            elseif with_trunc
                with_trunc_mem += memory
            end
        elseif startswith(line, "reduced_pairs =>")
            reduced_pairs = get_line_value(line)
            if no_trunc
                no_trunc_reduced_pairs = reduced_pairs
            elseif with_trunc
                with_trunc_reduced_pairs = reduced_pairs
            end
        elseif startswith(line, "4ti2 time =>")
            fti2_time = get_line_value(line)
        end
    end
    #Finished reading log file, make a vector with the data
    #This vector will later be turned into a dataframe
    #Vector components: instance_name, n, m, no_trunc_time, no_trunc_mem,
    #no_trunc_reduced_pairs, with_trunc_time, with_trunc_mem, with_trunc_reduced_pairs,
    #fti2_time
    v = [instance_name(filename), no_trunc_time, no_trunc_mem, no_trunc_reduced_pairs,
        with_trunc_time, with_trunc_mem, with_trunc_reduced_pairs, fti2_time]
    return v
end

function read_output_data()
    data = []
    for dirname in readdir()
        if isdir(dirname)
            for filename in readdir(dirname, join = true)
                if endswith(filename, ".log")
                    open(filename) do f
                        log = readlines(f)
                        v = read_log(filename, log)
                        push!(data, v)
                    end
                end
            end
        end
    end
    return data
end

function make_dataframe(data)
    col_names = [:instance, :no_trunc_time, :no_trunc_mem, :no_trunc_reduced_pairs,
        :with_trunc_time, :with_trunc_mem, :with_trunc_reduced_pairs, :fti2_time]
    cols = []
    for i in 1:length(col_names)
        if i == 1
            col = String[]
        elseif i == 4 || i == 7
            col = Int[]
        else
            col = Float64[]
        end
        for j in eachindex(data)
            push!(col, data[j][i])
        end
        push!(cols, col)
    end
    df = DataFrame(cols, col_names)
    return df
end

function main()
    data = read_output_data()
    df = make_dataframe(data)
    return df
end
