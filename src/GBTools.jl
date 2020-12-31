#TODO include more conversion functions from my formats to 4ti2 formats
#This will help automate all the tests I need
module GBTools

function isincluded(
    gb1 :: Vector{Vector{Int}},
    gb2 :: Vector{Vector{Int}}
) :: Bool
    for g in gb1
        if !(g in gb2)
            return false
        end
    end
    return true
end

function isequal(
    gb1 :: Vector{Vector{Int}},
    gb2 :: Vector{Vector{Int}}
) :: Bool
    return isincluded(gb1, gb2) && isincluded(gb2, gb1)
end

function tomatrix(
    gb :: Vector{Vector{Int}}
) :: Array{Int, 2}
    M = foldl(hcat, gb)
    return M'
end

function tovector(
    gb :: Array{Int, 2}
) :: Vector{Vector{Int}}
    return [ gb[i, :] for i in 1:size(gb, 1) ]
end

end
