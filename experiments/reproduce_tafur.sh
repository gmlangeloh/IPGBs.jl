#!/usr/bin/bash

julia run_3objective.jl 4.15 > table415.dat &
julia run_3objective.jl 4.16 > table416.dat &
julia run_3objective.jl 4.17 > table417.dat &
julia run_3objective.jl 4.18 > table418.dat

./clean_4ti2_files.sh
