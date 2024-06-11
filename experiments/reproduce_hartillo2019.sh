#!/usr/bin/bash

julia run_moip_knapsack.jl grobnecon 1 > table1.dat
julia run_moip_knapsack.jl grobnecon 2 > table2.dat
julia run_moip_knapsack.jl grobnecon 3 > table3.dat
julia run_moip_knapsack.jl grobnecon 4 > table4.dat
julia run_moip_knapsack.jl grobnecon 5 > table5.dat
julia run_moip_knapsack.jl grobnecon 6 > table6.dat

#julia run_knapsack.jl augmecon 1 > table1augmecon.dat
#julia run_knapsack.jl augmecon 2 > table2augmecon.dat
#julia run_knapsack.jl augmecon 3 > table3augmecon.dat
#julia run_knapsack.jl augmecon 4 > table4augmecon.dat
#julia run_knapsack.jl augmecon 5 > table5augmecon.dat
#julia run_knapsack.jl augmecon 6 > table6augmecon.dat
