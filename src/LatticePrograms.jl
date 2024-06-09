module LatticePrograms

"""
    Lattice program in the form:

    min {c^T x : x - u in L}

    where
        c : objective
        u : fiber_solution
        L : the lattice generated by the columns of basis
"""
struct LatticeProgram
    basis :: Matrix{Int} # n x (n - d) matrix
    fiber_solution :: Vector{Int}
    objective :: Vector{Int}
    n :: Int
    d :: Int
end

complement(inds :: Vector{Int}, n :: Int) = setdiff(1:n, inds)

function relaxation(lat :: LatticeProgram, away_from :: Vector{Int})
    complement = complement(away_from, lat.n)
    new_fiber = lat.fiber_solution[complement]
    new_objective = lat.objective[complement]
    new_basis = lat.basis[:, complement]
    return LatticeProgram(new_basis, new_fiber, new_objective, length(complement), lat.d)
end

end