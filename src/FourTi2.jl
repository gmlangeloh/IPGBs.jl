"""
Basic Julia-4ti2 interface using system calls.
Currently includes the minimize, groebner, normalform, markov and graver 4ti2
commands.
"""
module FourTi2
export minimize, groebner, normalform, markov, groebnernf, graver

using DelimitedFiles
using LinearAlgebra

"""
Internal use.

Writes a Julia array to a file in a format that can be read by 4ti2.
"""
function _4ti2_write(
    v :: Array{T},
    filename :: String
) where {T <: Real}
    open(filename, "w") do io
        if ndims(v) == 1
            writedlm(io, collect(size(v'))')
            writedlm(io, v')
        else
            writedlm(io, collect(size(v))')
            writedlm(io, v)
        end
    end
end

"""
Internal use.

Reads a Julia array from a given file written by 4ti2.
"""
function _4ti2_read(
    filename :: String
)
    x = readdlm(filename, header=true)[1]
    return x
end

"""
Internal use.

Removes the temporary files used to communicate between Julia and 4ti2.
"""
function _4ti2_clear(
    filename :: String
)
    extensions = [
        ".zsol", ".lat", ".mat", ".cost", ".gro", "gra", ".feas", ".nf",
        ".sign", ".min"
    ]
    for ext in extensions
        path = filename * ext
        rm = `rm -f $path`
        run(rm)
    end
end

"""
Writes an integer programming problem in 4ti2-readable format.
"""
function write_4ti2_input(
    A :: Array{Int, 2},
    c :: Array{T},
    xinit :: Vector{Int},
    project_name = "tmp" :: String,
    nonnegative = true :: String
) where {T <: Real}
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    xinit_file = project_name * ".zsol"
    _4ti2_write(xinit, xinit_file)
    write_4ti2_sign(nonnegative, size(A, 2), project_name)
end

"""
Writes sign file for 4ti2.
"""
function write_4ti2_sign(
    nonnegative :: Bool,
    nvars :: Int,
    project_name = "tmp" :: String
)
    sign_file = project_name * ".sign"
    if nonnegative
        _4ti2_write(ones(Int, nvars), sign_file)
    else
        _4ti2_write(zeros(Int, nvars), sign_file)
    end
end

"""
Interface with the 4ti2 minimize command.
Solves the Integer Programming problem:

min c' * x
s.t. A * x = b
x >= 0

where b = A * xinit using Gröbner Basis methods / 4ti2. This is essentially
equivalent to calling `groebner` followed by `normalform`, but a few
optimizations are applied (group relaxations + upper bounds)

I found that --precision=arb is often necessary here to avoid errors.
"""
function minimize(
    A :: Array{Int, 2},
    c :: Array{T},
    xinit :: Vector{Int},
    project_name = "tmp" :: String,
    nonnegative = true :: Bool;
    timeout :: Union{Int, Nothing} = nothing
) :: Tuple{Vector{Int}, Int} where {T <: Real}
    _4ti2_clear(project_name)
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    xinit_file = project_name * ".zsol"
    _4ti2_write(xinit, xinit_file)
    write_4ti2_sign(nonnegative, size(A, 2), project_name)

    #Run 4ti2
    if timeout == nothing
        cmd = `minimize -q --precision=arb $project_name`
    else
        cmd = `timeout $timeout minimize -q --algorithm=weighted --precision=arb $project_name`
    end
    try
        run(cmd)
    catch e #Deal with timeouts by just returning an empty solution
        if isa(e, ProcessFailedException)
            return Int[], 0
        else
            throw(e)
        end
    end

    #Retrieve results
    out_file = project_name * ".min"
    x = _4ti2_read(out_file)
    x = reshape(Int.(x), size(x, 2))

    if ndims(c) > 1
        obj = c[1, :]
    else
        obj = c
    end

    return x, obj' * x
end

"""
Interface with the 4ti2 groebner command.
Computes the Gröbner Basis of the toric ideal given by `A` using `c` as the
ordering.

If `lattice` is true, A is a lattice basis instead.

Returns the Gröbner Basis as a matrix where each row corresponds to a pure
binomial in it.
"""
function groebner(
    A :: Array{Int, 2},
    c :: Union{Array{T}, Nothing} = nothing,
    project_name :: String = "tmp",
    nonnegative :: Bool = true;
    truncation_sol = [],
    lattice :: Bool = false,
    quiet :: Bool = true
) :: Array{Int, 2} where {T <: Real}
    _4ti2_clear(project_name)
    #Write the project files
    if !lattice
        matrix_file = project_name * ".mat"
        _4ti2_write(A, matrix_file)
    else
        lattice_file = project_name * ".lat"
        _4ti2_write(A, lattice_file)
    end
    if !isnothing(c)
        obj_file = project_name * ".cost"
        _4ti2_write(c, obj_file)
    end
    write_4ti2_sign(nonnegative, size(A, 2), project_name)
    #Set options for truncated bases
    truncation_opt = ""
    if length(truncation_sol) > 0
        truncation_file = project_name * ".zsol"
        _4ti2_write(truncation_sol, truncation_file)
        truncation_opt = "--truncation=ip"
    end
    #Run 4ti2
    quiet_opt = quiet ? "-q" : ""
    cmd = `groebner $quiet_opt $truncation_opt $project_name`
    run(cmd)

    out_file = project_name * ".gro"
    gb = _4ti2_read(out_file)
    gb = Int.(gb)

    #Returns the gb as a matrix where rows are elements of the test set
    return gb
end

"""
Interface with the 4ti2 normalform command.
Computes the normal form of `xinit` with respect to the Gröbner Basis of the toric
ideal given by `A` using `c` as an ordering. This can be called multiple times
for various xinit without recomputing the Gröbner Basis.

Returns the normal form and its value using `c` as objective function.
"""
function normalform(
    A :: Array{Int, 2},
    c :: Array{T},
    xinit :: Vector{Int},
    project_name = "tmp" :: String,
    nonnegative = true :: Bool
) :: Tuple{Vector{Int}, Int} where {T <: Real}
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    xinit_file = project_name * ".feas"
    _4ti2_write(xinit, xinit_file)
    write_4ti2_sign(nonnegative, size(A, 2), project_name)

    #Run 4ti2
    cmd = `normalform -q $project_name`
    run(cmd)

    #Retrieve results
    out_file = project_name * ".nf"
    x = _4ti2_read(out_file)
    x = reshape(Int.(x), size(x, 2))

    if ndims(c) > 1
        obj = c[1, :]
    else
        obj = c
    end

    return x, obj' * x
end

"""
Calls 4ti2's `groebner` command followed by `normalform`. This is equivalent to
`minimize`, except it generates a Gröbner basis in as output and doesn't apply
the optimizations of the latter command.

Returns the optimal solution to the given IP along with its objective value.
"""
function groebnernf(
    A :: Array{Int, 2},
    c :: Array{T},
    xinit :: Vector{Int},
    project_name = "tmp" :: String,
    nonnegative = true :: Bool
) :: Tuple{Vector{Int}, Int} where {T <: Real}
    _ = groebner(A, c, project_name, nonnegative, truncation_sol=xinit)
    return normalform(A, c, xinit, project_name, nonnegative)
end

"""
Returns the minimal Markov basis of A with respect to c as a matrix,
where the rows are the elements of the basis.
"""
function markov(
    A :: Array{Int, 2},
    c :: Array{T},
    project_name = "tmp" :: String,
    nonnegative = true :: Bool;
    truncation_sol = []
) :: Array{Int, 2} where {T <: Real}
    _4ti2_clear(project_name)
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    write_4ti2_sign(nonnegative, size(A, 2), project_name)
    #Set options for truncated bases
    truncation_opt = ""
    if length(truncation_sol) > 0
        truncation_file = project_name * ".zsol"
        _4ti2_write(truncation_sol, truncation_file)
        truncation_opt = "--truncation=ip"
    end
    #Run 4ti2
    cmd = `markov -q $truncation_opt $project_name`
    run(cmd)

    out_file = project_name * ".mar"
    markov_basis = Int.(_4ti2_read(out_file))

    #Returns the markov basis as a matrix where rows are elements of the test set
    return markov_basis
end

"""
Calls 4ti2's graver command.

TODO Implement other options (lower and upper bound, arbitrary signs,
support to lattice bases instead of matrices...)
"""
function graver(
    A :: Array{Int, 2},
    project_name :: String = "tmp",
    nonnegative :: Bool = true
) :: Array{Int, 2}
    _4ti2_clear(project_name)
    #Write project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    #Run 4ti2
    cmd = `graver -q $project_name`
    run(cmd)

    out_file = project_name * ".gra"
    graver_basis = Int.(_4ti2_read(out_file))
    return graver_basis
end

end
