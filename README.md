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
you might need to emit the following as well:
    git submodule update --init --recursive

## Build

    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/10.2.0_4/bin/g++-10
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

I got the following in my Mac (16GB RAM)
# Harmonic

|                         Method |      Factor |       Solve |     L∞ norm |
|-------------------------------:|------------:|------------:|------------:|
|    Eigen::CholmodSupernodalLLT |    1.4 secs |   0.11 secs | 1.33919e-10 |
|           Eigen::SimplicialLLT |    1.8 secs |   0.14 secs | 9.85347e-11 |
|          Eigen::SimplicialLDLT |    1.8 secs |   0.14 secs | 6.01932e-11 |
|            catamari::SparseLDL |      2 secs |   0.12 secs | 1.15886e-10 |
|              Eigen::PardisoLLT |    2.2 secs |   0.52 secs | 3.30178e-11 |
|                Eigen::SparseLU |    5.8 secs |   0.22 secs | 8.13434e-11 |
|             Sympiler::Cholesky |    1.3 secs |  0.039 secs | 2.45269e-11 |

# Biharmonic

|                         Method |      Factor |       Solve |     L∞ norm |
|-------------------------------:|------------:|------------:|------------:|
|    Eigen::CholmodSupernodalLLT |    2.5 secs |   0.17 secs | 3.10268e-05 |
|           Eigen::SimplicialLLT |     14 secs |   0.41 secs | 4.64678e-05 |
|          Eigen::SimplicialLDLT |     14 secs |   0.44 secs | 2.12954e-05 |
|            catamari::SparseLDL |     12 secs |   0.43 secs | 8.58138e-05 |
|              Eigen::PardisoLLT |    4.3 secs |   0.66 secs | 3.91748e-05 |
|                Eigen::SparseLU |     35 secs |   0.61 secs | 1.20096e-05 |
|             Sympiler::Cholesky |    2.7 secs |  0.084 secs | 4.04411e-05 |

# Triharmonic

|                         Method |      Factor |       Solve |     L∞ norm |
|-------------------------------:|------------:|------------:|------------:|
|    Eigen::CholmodSupernodalLLT |    4.8 secs |    0.3 secs | 56.5885 |
|           Eigen::SimplicialLLT |     47 secs |   0.86 secs | 61.6743 |
|          Eigen::SimplicialLDLT |     47 secs |   0.98 secs | 30.9774 |
|            catamari::SparseLDL |     49 secs |   0.86 secs | 94.3167 |
|              Eigen::PardisoLLT |      8 secs |   0.94 secs | 25.9679 |
|                Eigen::SparseLU | 1.5e+02 secs |    1.4 secs | 16.2049 |
|             Sympiler::Cholesky |    4.6 secs |   0.14 secs | 54.6663 |
