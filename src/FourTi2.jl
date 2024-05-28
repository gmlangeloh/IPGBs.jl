"""
Basic Julia-4ti2 interface using system calls.
Currently includes the minimize, groebner, normalform, markov and graver 4ti2
commands.
"""
module FourTi2
export minimize, groebner, normalform, markov, groebnernf, graver

using DelimitedFiles
using JuMP
using LinearAlgebra

using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances
using MIPMatrixTools.MatrixTools

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
        ".sign", ".min", ".mar", ".err"
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
    nonnegative :: Vector{Bool},
    project_name = "tmp" :: String,
) where {T <: Real}
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    xinit_file = project_name * ".zsol"
    _4ti2_write(xinit, xinit_file)
    write_4ti2_sign(nonnegative, project_name)
end

"""
Writes sign file for 4ti2.
"""
function write_4ti2_sign(
    nonnegative :: Vector{Bool},
    project_name = "tmp" :: String
)
    sign_file = project_name * ".sign"
    _4ti2_write(Int.(nonnegative), sign_file)
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
    xinit :: Vector{Int};
    nonnegative :: Vector{Bool} = Bool[],
    project_name :: String = "tmp",
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
    if isempty(nonnegative)
        nonnegative = fill(true, size(A, 2))
    end
    write_4ti2_sign(nonnegative, project_name)

    #Run 4ti2
    if isnothing(timeout)
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

function minimize(
    instance :: IPInstance;
    solution :: Vector{Int} = Int[],
    project_name :: String = "tmp"
)
    nonnegative = IPInstances.nonnegative_variables(instance)
    int_objective = IPInstances.integer_objective(instance)
    if isempty(solution)
        init_sol = MatrixTools.initial_solution(instance.A, instance.b)
    else
        init_sol = solution
    end
    return minimize(
        instance.A, int_objective, init_sol, nonnegative=nonnegative,
        project_name=project_name
    )
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
    A::Array{Int,2},
    c::Array{T};
    nonnegative::Vector{Bool} = Bool[],
    project_name::String = "tmp",
    truncation_sol = Int[],
    lattice::Bool = false,
    quiet::Bool = true,
    markov::Union{Nothing,Vector{Vector{Int}}} = nothing
)::Array{Int,2} where {T<:Real}
    _4ti2_clear(project_name)
    #Write the project files
    if !lattice
        matrix_file = project_name * ".mat"
        _4ti2_write(A, matrix_file)
    else
        lattice_file = project_name * ".lat"
        _4ti2_write(A, lattice_file)
    end
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    if isempty(nonnegative)
        nonnegative = fill(true, size(A, 2))
    end
    write_4ti2_sign(nonnegative, project_name)
    if !isnothing(markov)
        markov_file = project_name * ".mar"
        M = zeros(Int, length(markov), length(markov[1]))
        for i in eachindex(markov)
            for j in eachindex(markov[i])
                M[i, j] = markov[i][j]
            end
        end
        _4ti2_write(M, markov_file)
    end
    #Set options for truncated bases
    truncation_opt = ""
    if length(truncation_sol) > 0
        truncation_file = project_name * ".zsol"
        _4ti2_write(truncation_sol, truncation_file)
        truncation_opt = "--truncation=lp"
    end
    #Run 4ti2
    quiet_opt = quiet ? "-q" : ""
    cmd = `groebner -parb $quiet_opt $truncation_opt $project_name`
    error_file = project_name * ".err"

    unbounded = false
    try
        run(pipeline(cmd, stderr=error_file))
    catch
        #If 4ti2 fails, check why by reading the error file
        #If the cost function is unbounded, deal with it.
        #default to returning a basis with the zero vector in these cases.
        err = open(error_file, "r")
        error_message = read(err, String)
        if !occursin("Cost function is not bounded.", error_message)
            rethrow()
        end
        close(err)
        unbounded = true
    end
    if unbounded
        return zeros(Int, 1, size(A, 2))
    end

    out_file = project_name * ".gro"
    gb = _4ti2_read(out_file)
    gb = Int.(gb)

    #Returns the gb as a matrix where rows are elements of the test set
    return gb
end

function make_project(
    instance :: IPInstance;
    markov :: Union{Nothing,Vector{Vector{Int}}} = nothing,
    project_name :: String = "tmp"
)
    nonnegative = IPInstances.nonnegative_variables(instance)
    int_objective = IPInstances.integer_objective(instance)
    truncation_sol = MatrixTools.initial_solution(instance.A, instance.b)
    _4ti2_clear(project_name)
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(instance.A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(int_objective, obj_file)
    if isempty(nonnegative)
        nonnegative = fill(true, size(instance.A, 2))
    end
    write_4ti2_sign(nonnegative, project_name)
    if !isnothing(markov)
        markov_file = project_name * ".mar"
        M = zeros(Int, length(markov), length(markov[1]))
        for i in eachindex(markov)
            for j in eachindex(markov[i])
                M[i, j] = markov[i][j]
            end
        end
        _4ti2_write(M, markov_file)
    end
    #Set options for truncated bases
    truncation_opt = ""
    if length(truncation_sol) > 0
        truncation_file = project_name * ".zsol"
        _4ti2_write(truncation_sol, truncation_file)
        truncation_opt = "--truncation=lp"
    end
end

function groebner(
    instance :: IPInstance;
    markov :: Union{Nothing,Vector{Vector{Int}}} = nothing,
    project_name :: String = "tmp"
)
    nonnegative = IPInstances.nonnegative_variables(instance)
    int_objective = IPInstances.integer_objective(instance)
    init_sol = MatrixTools.initial_solution(instance.A, instance.b)
    return groebner(
        instance.A, int_objective, nonnegative=nonnegative, markov=markov,
        truncation_sol=init_sol, project_name=project_name
    )
end

function groebner(
    model :: JuMP.Model;
    project_name :: String = "tmp"
)
    return groebner(IPInstance(model), project_name=project_name)
end

function groebner(
    filename :: String;
    project_name :: String = "tmp"
)
    return groebner(IPInstance(filename), project_name=project_name)
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
    c :: Array{Int},
    xinit :: Vector{Int};
    nonnegative :: Vector{Bool} = Bool[],
    project_name = "tmp" :: String,
) :: Tuple{Vector{Int}, Int}
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    xinit_file = project_name * ".feas"
    _4ti2_write(xinit, xinit_file)
    if isempty(nonnegative)
        nonnegative = fill(true, size(A, 2))
    end
    write_4ti2_sign(nonnegative, project_name)

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

function normalform(
    instance :: IPInstance,
    xinit :: Vector{Int};
    project_name :: String = "tmp"
) :: Tuple{Vector{Int}, Int}
    return normalform(
        instance.A, round.(Int, instance.C), xinit,
        nonnegative=IPInstances.nonnegative_variables(instance),
        project_name=project_name
    )
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
    xinit :: Vector{Int};
    nonnegative :: Vector{Bool},
    project_name = "tmp" :: String,
) :: Tuple{Matrix{Int}, Vector{Int}, Int} where {T <: Real}
    gb = groebner(
        A, c, nonnegative=nonnegative, project_name=project_name,
        truncation_sol=xinit
    )
    sol, val = normalform(
        A, round(Int, c), xinit,
        project_name=project_name, nonnegative=nonnegative
    )
    return gb, sol, val
end

function groebnernf(
    instance :: IPInstance,
    markov :: Vector{Vector{Int}},
    xinit :: Vector{Int};
    project_name :: String = "tmp"
) :: Tuple{Matrix{Int}, Vector{Int}, Int}
    gb = groebner(instance, markov=markov, project_name=project_name)
    sol, val = normalform(instance, xinit, project_name=project_name)
    return gb, sol, val
end

"""
Returns the minimal Markov basis of A with respect to c as a matrix,
where the rows are the elements of the basis.
"""
function markov(
    A :: Array{Int, 2},
    c :: Array{T};
    nonnegative :: Vector{Bool} = Bool[],
    project_name :: String = "tmp",
    truncation_sol :: Vector{Int} = Int[]
) :: Array{Int, 2} where {T <: Real}
    _4ti2_clear(project_name)
    #Write the project files
    matrix_file = project_name * ".mat"
    _4ti2_write(A, matrix_file)
    obj_file = project_name * ".cost"
    _4ti2_write(c, obj_file)
    if isempty(nonnegative)
        nonnegative = fill(true, size(A, 2))
    end
    write_4ti2_sign(nonnegative, project_name)
    #Set options for truncated bases
    truncation_opt = ""
    if length(truncation_sol) > 0
        truncation_file = project_name * ".zsol"
        _4ti2_write(truncation_sol, truncation_file)
        truncation_opt = "--truncation=ip"
    end
    #Run 4ti2
    cmd = `markov -parb -q $truncation_opt $project_name`
    run(cmd)

    out_file = project_name * ".mar"
    markov_basis = Int.(_4ti2_read(out_file))

    #Returns the markov basis as a matrix where rows are elements of the test set
    return markov_basis
end

function markov(
    instance :: IPInstance;
    project_name :: String = "tmp"
)
    nonnegative = IPInstances.nonnegative_variables(instance)
    int_objective = IPInstances.integer_objective(instance)
    init_sol = MatrixTools.initial_solution(instance.A, instance.b)
    return markov(instance.A, int_objective, nonnegative=nonnegative,
        truncation_sol=init_sol, project_name=project_name)
end

"""
Calls 4ti2's graver command.

TODO: Implement other options (lower and upper bound, arbitrary signs,
support to lattice bases instead of matrices...)
"""
function graver(
    A :: Array{Int, 2},
    nonnegative :: Vector{Bool},
    project_name :: String = "tmp",
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
