# IPGBs

The IPGBs package provides tools for the solution of Integer Programs with Gröbner Bases. Currently, it provides
a state-of-the-art implementation of the Geometric Buchberger algorithm and a work in progress implementation
of Project-and-lift.

It also provides a Julia interface to [4ti2](https://github.com/4ti2/4ti2), another state-of-the-art
Gröbner basis solver. In particular, it allows using the `groebner`, `normalform` and `minimize`
commands with either input matrices, as in the original 4ti2 interface, or using a JuMP model.

## Installation

IPGBs 4ti2 interface depends on a local installation of 4ti2 and on the desired 4ti2 commands being available in the PATH.

## Using IPGBs

The easiest way to use IPGBs is by building an IP with JuMP and then calling IPGBs `groebner_basis` command.
For example:

```julia
    using IPGBs
    using JuMP
    using IPGBs.FourTi2 #To access the 4ti2 interface

    c = [8, 5, 10, 6, 4]
    A = [10 6 2 5 3]
    b = [13]
    model = Model()
    @variable(model, x[1:5] >= 0, Int)
    @objective(model, Max, c' * x)
    @constraint(model, A * x <= b)

    #Computing the Gröbner basis with IPGBs
    gb_ipgbs = groebner_basis(model)

    #Computing the same Gröbner basis with the 4ti2 interface
    gb_4ti2 = FourTi2.groebner(model)
```

Alternatively, if the IP is stored in any file format that is readable by JuMP, IPGBs can be run
directly.

```julia
    using IPGBs
    using JuMP
    using IPGBs.FourTi2 #To access the 4ti2 interface

    gb_ipgbs = groebner_basis("integer_programming_instance.mps")
    gb_4ti2 = FourTi2.groebner("integer_programming_instance.mps")
```

## Experimental results

Current experimental results for IPGBs can be found in the `experiments` directory.
It provides .mps files for each instance and logs from which results can be extracted with
the `generate_tables.jl` script. The script `full_groebner_basis.jl` can be used to
rerun IPGBs and 4ti2 over all instances in this directory.
