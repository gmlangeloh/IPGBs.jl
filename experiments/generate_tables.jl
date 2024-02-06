using StatsKit

function get_line_value(line :: String)
    return parse(Float64, split(line, "=>")[2])
end

function instance_name(filename)
    #take filename, remove path and extension
    no_path = split(filename, "/")[end]
    #To remove the extension, simply remove the last two elements of the array
    #This is because the extension is .mps.log
    name_with_reps = join(split(no_path, ".")[1:end-2], ".")
    #Delete the number of reps, the number after last _
    name = join(split(name_with_reps, "_")[1:end-1], "_")
    return name
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
    max_basis_size = 0
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
        elseif startswith(line, "max_basis_size =>")
            max_basis_size = parse(Int, split(line, "=>")[2])
        end
    end
    #Finished reading log file, make a vector with the data
    #This vector will later be turned into a dataframe
    #Vector components: instance_name, n, m, no_trunc_time, no_trunc_mem,
    #no_trunc_reduced_pairs, with_trunc_time, with_trunc_mem, with_trunc_reduced_pairs,
    #fti2_time, max_basis_size
    v = [instance_name(filename), no_trunc_time, no_trunc_mem, no_trunc_reduced_pairs,
        with_trunc_time, with_trunc_mem, with_trunc_reduced_pairs, fti2_time, max_basis_size]
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
        :with_trunc_time, :with_trunc_mem, :with_trunc_reduced_pairs, :fti2_time,
        :gb_size]
    cols = []
    for i in 1:length(col_names)
        if i == 1
            col = String[]
        elseif i == 4 || i == 7 || i == 9
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

function problem_from_instance_name(name)
    if startswith(name, "knapsack_binary_multidimensional")
        return "0-1MDKP"
    elseif startswith(name, "knapsack_binary")
        return "0-1KP"
    elseif startswith(name, "knapsack_multidimensional")
        return "MDKP"
    elseif startswith(name, "knapsack_unbounded")
        return "UKP"
    elseif startswith(name, "set_cover")
        return "SCP"
    elseif startswith(name, "set_packing")
        return "SPP"
    end
    return ""
end

function parameters_from_instance_name(name)
    #Instance name formats:
    # knapsack_binary_n
    # knapsack_binary_multidimensional_n_m
    # knapsack_multidimensional_n_m
    # knapsack_unbounded_n
    # knapsack_unbounded_n_corr
    # set_cover_n_m_p
    # set_packing_n_m_p
    #Extract n, m from instance name
    n = 0
    m = 1
    if startswith(name, "knapsack_binary")
        if occursin("multidimensional", name)
            n = parse(Int, split(name, "_")[4])
            m = parse(Int, split(name, "_")[5])
        else
            n = parse(Int, split(name, "_")[3])
        end
    elseif startswith(name, "knapsack_multidimensional")
        n = parse(Int, split(name, "_")[3])
        m = parse(Int, split(name, "_")[4])
    elseif startswith(name, "knapsack_unbounded")
        if occursin("corr", name)
            n = parse(Int, split(name, "_")[3])
        else
            n = parse(Int, split(name, "_")[3])
        end
    elseif startswith(name, "set_cover")
        n = parse(Int, split(name, "_")[3])
        m = parse(Int, split(name, "_")[4])
    elseif startswith(name, "set_packing")
        n = parse(Int, split(name, "_")[3])
        m = parse(Int, split(name, "_")[4])
    end
    return n, m
end

function summarize_dataframe(df)
    instance_names = unique(df.instance)
    data_cols = names(df)
    filter!(name -> name != "instance", data_cols)
    problem_names = String[]
    ns = Int[]
    ms = Int[]
    gb_sizes = Float64[]
    means_no_trunc_time = Float64[]
    means_no_trunc_mem = Float64[]
    means_no_trunc_reduced_pairs = Float64[]
    means_with_trunc_time = Float64[]
    means_with_trunc_mem = Float64[]
    means_with_trunc_reduced_pairs = Float64[]
    means_percent_reduced_pairs = Float64[]
    means_percent_time = Float64[]
    means_percent_mem = Float64[]
    means_fti2_time = Float64[]
    for name in instance_names
        println("Processing: ", name)
        inst_data = df[df.instance .== name, :]
        prob_name = problem_from_instance_name(name)
        push!(problem_names, prob_name)
        n, m = parameters_from_instance_name(name)
        push!(ns, n)
        push!(ms, m)
        m_gb_size = mean(inst_data.gb_size)
        push!(gb_sizes, m_gb_size)
        m_no_trunc_time = mean(inst_data.no_trunc_time)
        push!(means_no_trunc_time, m_no_trunc_time)
        m_no_trunc_mem = mean(inst_data.no_trunc_mem)
        push!(means_no_trunc_mem, m_no_trunc_mem)
        m_no_trunc_reduced_pairs = mean(inst_data.no_trunc_reduced_pairs)
        push!(means_no_trunc_reduced_pairs, m_no_trunc_reduced_pairs)
        m_with_trunc_time = mean(inst_data.with_trunc_time)
        push!(means_with_trunc_time, m_with_trunc_time)
        m_with_trunc_mem = mean(inst_data.with_trunc_mem)
        push!(means_with_trunc_mem, m_with_trunc_mem)
        m_with_trunc_reduced_pairs = mean(inst_data.with_trunc_reduced_pairs)
        push!(means_with_trunc_reduced_pairs, m_with_trunc_reduced_pairs)
        m_percent_reduced_pairs = 100 * m_with_trunc_reduced_pairs / m_no_trunc_reduced_pairs
        push!(means_percent_reduced_pairs, m_percent_reduced_pairs)
        m_percent_time = 100 * m_with_trunc_time / m_no_trunc_time
        push!(means_percent_time, m_percent_time)
        m_percent_mem = 100 * m_with_trunc_mem / m_no_trunc_mem
        push!(means_percent_mem, m_percent_mem)
        m_fti2_time = mean(inst_data.fti2_time)
        push!(means_fti2_time, m_fti2_time)
    end
    col_names = ["Problem", "n", "m", "|G|", "t (no binary)", "m (no binary)", "reduced pairs (no binary)",
        "t (with binary)", "m (with binary)", "reduced pairs (with binary)", "t (4ti2)"]
    summarized_df = DataFrame([problem_names, ns, ms, gb_sizes, means_no_trunc_time, means_no_trunc_mem,
        means_no_trunc_reduced_pairs, means_with_trunc_time, means_with_trunc_mem,
        means_with_trunc_reduced_pairs, means_fti2_time], col_names)
    percent_col_names = ["Problem", "n", "m", "% time", "% memory", "% reduced pairs"]
    percent_df = DataFrame([problem_names, ns, ms, means_percent_time, means_percent_mem,
        means_percent_reduced_pairs], percent_col_names)
    return summarized_df, percent_df
end

function main()
    data = read_output_data()
    df = make_dataframe(data)
    sum_df, percent_df = summarize_dataframe(df)
    return sum_df, percent_df
end
