# Reflector_1D_julia

Reflector_1D_julia is a julia program aimed to solve the inverse reflector problem in 1D with a plane source emitting toward its normal and a a target plane in the far field.

# Dependencies

This program depends on the following julia packages:
+ PyPlot
+ Grid

You can install then with the following command in julia:

``` sh
Pkg.add("packageName")
```

# Parameters
For instance, all the parameters have to be set directly in the reflector_1D_v0.1.jl file.
Their denomination is the following:
+ Nx : Number of samples of the source
+ source_box : Source support
+ Np : Number of samples of the target
+ target_plane_box : Target plane support
+ e_xi : Unit vector of the target plane
+ n_plan : Orthogonal vector of the target plane. Its norm is the distance from the reflector
