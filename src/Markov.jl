"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis
"""
module Markov

using AbstractAlgebra

"""
Given a saturated lattice in the form L = {x | Ax = 0} find a lattice basis
for it using the AbstractAlgebra library.
"""
function lattice_basis(
    A :: Array{Int, 2}
) :: Vector{Vector{BigInt}}
    m, n = size(A)
    space = MatrixSpace(ZZ, m, n)
    rank, basis = kernel(A)
    julia_form = Array(basis)
    return [ julia_form[:, j] for j in 1:size(julia_form, 2)]
end

end
