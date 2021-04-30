# 🚀 Sparse Positive Definite Linear System Solver Benchmark 📏

This is an informal benchmark for _k_-harmonic diffusion problems on triangle
meshes found commonly in geometry processing.

> Current upshot: 
> 
> - CholMod is **very good** 🐐
> - Pardiso **holds up alright** 🏆
> - Eigen LLT/LDLT are **OK** but suffer for less sparse systems 🏅
> - catamari is **OK** for medium but inaccurate for big systems (bonus: it's MPL2 🆓)
> - Eigen LU¹ is significantly **slower** 🐌
>
> ¹These systems are
> [SPD](https://en.wikipedia.org/wiki/Definite_symmetric_matrix) so LU is not a
> good choice, but provides a reference.

## Clone

    git clone --recursive https://github.com/alecjacobson/sparse-solver-benchmark

## Build

    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=Release
    make

## Run

    ./sparse_solver_benchmark [path to triangle mesh]

## Example

Running

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply 

on my MacBook 2.3 GHz Quad-Core Intel Core i7 with 32 GB Ram will produce:

> # harmonic
> 
> |                         Method |      Factor |       Solve |     L∞ norm |
> |-------------------------------:|------------:|------------:|------------:|
> |    Eigen::CholmodSupernodalLLT |   0.86 secs |    0.1 secs | 1.13266e-10 |
> |           Eigen::SimplicialLLT |    1.5 secs |   0.12 secs | 1.19565e-10 |
> |          Eigen::SimplicialLDLT |    1.5 secs |   0.14 secs | 1.48007e-10 |
> |            catamari::SparseLDL |    1.7 secs |   0.11 secs | 1.27168e-10 |
> |              Eigen::PardisoLLT |    1.9 secs |   0.59 secs | 1.05662e-10 |
> |                Eigen::SparseLU |    4.8 secs |   0.18 secs | 6.85971e-11 |
> 
> # biharmonic
> 
> |                         Method |      Factor |       Solve |     L∞ norm |
> |-------------------------------:|------------:|------------:|------------:|
> |    Eigen::CholmodSupernodalLLT |    1.7 secs |   0.13 secs | 8.33117e-05 |
> |           Eigen::SimplicialLLT |     12 secs |   0.37 secs | 5.58785e-05 |
> |          Eigen::SimplicialLDLT |     12 secs |   0.41 secs | 6.92762e-05 |
> |            catamari::SparseLDL |     13 secs |   0.33 secs | 0.0002359 |
> |              Eigen::PardisoLLT |      4 secs |   0.69 secs | 6.34374e-05 |
> |                Eigen::SparseLU |     31 secs |    0.6 secs | 4.10639e-05 |
> 
> # triharmonic
> 
> |                         Method |      Factor |       Solve |     L∞ norm |
> |-------------------------------:|------------:|------------:|------------:|
> |    Eigen::CholmodSupernodalLLT |    3.3 secs |   0.19 secs | 42.2496 |
> |           Eigen::SimplicialLLT |     41 secs |    0.8 secs | 36.0525 |
> |          Eigen::SimplicialLDLT |     41 secs |   0.89 secs | 32.1019 |
> |            catamari::SparseLDL |     50 secs |   0.78 secs | 150.97 |
> |              Eigen::PardisoLLT |    6.9 secs |   0.81 secs | 96.7205 |
> |                Eigen::SparseLU | 1.3e+02 secs |    1.2 secs | 22.0579 |

Obviously [YMMV](https://www.google.com/search?q=YMMV), if you find something
interesting [let me know!](https://github.com/alecjacobson/sparse-solver-benchmark/issues).

## What about this other solver XYZ?

Please [submit a pull
request](https://github.com/alecjacobson/sparse-solver-benchmark/pulls) with a
wrapper for solver XYZ. The more the merrier. 


## What are the systems being solved?

This code will build a discretization of the ∆ᵏ operator and solve a system of
the form:

    ∆ᵏ u + u = x

where x is the surface's embedding. In matrix form this is:


    (Wᵏ + M) u = M x

where Wᵏ is defined recursively as:

    W¹ = L
    Wᵏ⁺¹ = Wᵏ M⁻¹ L

and `L` is the discrete Laplacian and `M` is the discrete mass matrix.

This is a form of smoothing (k=1 is implicit mean curvature flow "Implicit
Fairing of Irregular Meshes" Desbrun et al. 1999, k≥2 is higher order, e.g., "Mixed
Finite Elements for Variational Surface Modeling" Jacobson et al. 2010)

For k=1, the system is generally OK w.r.t. conditioning and the sparsity for a
minifold mesh will be 7 non-zeros per row (on average).

For k=3, the system can get really badly scaled and starts to become more dense
(~40 non-zeros per row).
