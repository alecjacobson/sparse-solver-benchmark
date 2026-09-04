# 🚀 Sparse Positive Definite Linear System Solver Leaderboard 📏

[![CI](https://github.com/alecjacobson/sparse-solver-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/alecjacobson/sparse-solver-benchmark/actions/workflows/ci.yml)

This is an informal leaderboard/benchmark for _k_-harmonic diffusion problems
on triangle meshes found commonly in geometry processing. Solvers are ranked
by total (factor + solve) time on each system; 🥇🥈🥉 mark the podium.

> Current upshot (see full results below, from an NVIDIA L40 / dual Intel
> Xeon Platinum 8362 machine):
>
> - **NVIDIA cuDSS** 🚀 wins every system here once the mesh is large (720K
>   vertices) — and its *solve* step is essentially free (single-digit
>   milliseconds) once factored, since it stays resident on the GPU.
> - **MKL Pardiso** 🏆 is the strongest CPU-only runner-up throughout.
> - **CHOLMOD / Eigen SimplicialLLT·LDLT / catamari** 🥈 remain solid, accurate,
>   dependency-light choices when no GPU is available.
> - **Eigen's iterative solvers** (BiCGSTAB/CG + IncompleteLUT) get
>   unreliable as k grows and can diverge outright at k=3 — this benchmark
>   caps them at 200 iterations so a hard system fails fast with an honest
>   (possibly huge, even `nan`) residual instead of hanging.
> - **Eigen SparseLU**¹ is still the slowest by a wide margin, as expected —
>   but see the ⚠️ note below if you see it apparently taking *minutes*
>   instead of seconds, that's a build misconfiguration, not real solver cost.
>
> ¹These systems are
> [SPD](https://en.wikipedia.org/wiki/Definite_symmetric_matrix) so LU is not
> a good choice, but provides a reference.

> [!WARNING]
> **Build gotcha: don't define `EIGEN_USE_MKL_ALL`.** It routes *every* Eigen
> dense matrix product — including the many tiny internal panel updates
> `Eigen::SparseLU` performs during factorization — through MKL's threaded
> GEMM. MKL's `dgemm_direct` path calls `pthread_create()` per small update in
> that hot loop, which turned a ~5 second `SparseLU` factorization into an
> ~850 second one (verified with `gdb`) on this machine. `Eigen::PardisoLLT`
> doesn't need this define — it links MKL directly, independent of Eigen's
> own dense dispatch. `CMakeLists.txt` deliberately only defines
> `IGL_WITH_MKL`, not `EIGEN_USE_MKL_ALL`. Similarly, if `SuiteSparse`'s own
> `find_package(BLAS)` picks a different MKL threading layer (e.g.
> `libiomp5`) than the one this project links (`libgomp`), you get two
> competing OpenMP runtimes spinning against each other in one process —
> `CMakeLists.txt` forces `BLA_VENDOR Intel10_64lp_seq` (sequential MKL) for
> SuiteSparse's own BLAS/LAPACK to avoid this.

## Clone

    git clone --recursive https://github.com/alecjacobson/sparse-solver-benchmark

## Build

    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=Release
    make -j8

CMake options (all default `ON`, gracefully disable themselves with a
warning if their dependency isn't found):

| Option              | Enables                                             |
|----------------------|------------------------------------------------------|
| `IGL_WITH_CHOLMOD`  | SuiteSparse CHOLMOD + UMFPACK (built from source)   |
| `IGL_WITH_MKL`      | Intel MKL Pardiso                                   |
| `IGL_WITH_GPL`      | GPL-licensed Eigen sparse solver code                |
| `IGL_WITH_CUDSS`    | NVIDIA cuDSS + cuSOLVER (requires the CUDA toolkit) |

## Run

    ./sparse_solver_benchmark [path to triangle mesh]
    ./sparse_solver_benchmark --grid 20 --check   # fast synthetic-mesh correctness check
    ./sparse_solver_benchmark mesh.ply --csv results.csv   # also dump raw results

## Testing

    ctest --test-dir build --output-on-failure

Runs the benchmark on a small synthetic grid mesh (`igl::triangulated_grid`,
not the 720K-vertex dragon) and fails if any compiled-in solver's L∞ residual
exceeds a per-_k_ tolerance — a correctness regression test, not a
performance one. GitHub Actions runs this on Linux (full solver matrix),
Windows, and macOS (Eigen + CHOLMOD/UMFPACK) on every push/PR; see
`.github/workflows/ci.yml`. Hosted runners have no GPU, so cuDSS/cuSOLVER
still compile there but detect the missing CUDA device at runtime and report
`skipped` rather than fail.

## Example

Running

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

on an NVIDIA L40 / dual Intel Xeon Platinum 8362 (128 threads) machine
produces:

# Harmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |            Eigen::SimplicialLDLT |      1.3 secs |     0.13 secs | 1.12086e-10 |
| 🥈 2 |             Eigen::SimplicialLLT |      1.3 secs |     0.12 secs | 4.55334e-11 |
| 🥉 3 |              catamari::SparseLDL |      1.5 secs |     0.12 secs | 3.82439e-11 |
|    4 |                     NVIDIA cuDSS |      2.2 secs |   0.0027 secs | 1.59312e-10 |
|    5 |   Eigen::BiCGSTAB\<IncompleteLUT\> |      1.6 secs |      1.4 secs | 1.23985e-10 |
|    6 |                Eigen::PardisoLLT |      3.1 secs |     0.88 secs | 7.58549e-11 |
|    7 |         Eigen::CG\<IncompleteLUT\> |      1.7 secs |        3 secs | 8.96274e-11 |
|    8 |                  Eigen::SparseLU |      5.3 secs |     0.21 secs | 2.37845e-11 |
|    9 |      Eigen::CholmodSupernodalLLT |      8.2 secs |     0.73 secs | 6.63736e-11 |
|   10 |        NVIDIA cuSOLVER (Sp Chol) |        0 secs |      9.7 secs | 5.25522e-11 |
|   11 |                 Eigen::UmfPackLU |       28 secs |     0.78 secs | 4.20999e-11 |

# Biharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |        4 secs |   0.0084 secs | 9.26244e-05 |
| 🥈 2 |                Eigen::PardisoLLT |        5 secs |     0.91 secs | 6.93083e-05 |
| 🥉 3 |            Eigen::SimplicialLDLT |      9.4 secs |     0.43 secs | 4.80425e-05 |
|    4 |             Eigen::SimplicialLLT |      9.5 secs |     0.44 secs | 2.60041e-05 |
|    5 |              catamari::SparseLDL |       11 secs |     0.45 secs | 3.05382e-05 |
|    6 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       11 secs |      4.6 secs | 5.56194e-05 |
|    7 |         Eigen::CG\<IncompleteLUT\> |       11 secs |      6.2 secs | 4.73183e-05 |
|    8 |        NVIDIA cuSOLVER (Sp Chol) |        0 secs |       20 secs | 3.40111e-05 |
|    9 |                  Eigen::SparseLU |       35 secs |      0.8 secs | 2.46911e-05 |
|   10 |      Eigen::CholmodSupernodalLLT |       46 secs |      6.7 secs | 9.99686e-05 |
|   11 |                 Eigen::UmfPackLU |  1.7e+02 secs |      2.3 secs | 8.19072e-05 |

# Triharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      6.7 secs |   0.0054 secs |     15.7028 |
| 🥈 2 |                Eigen::PardisoLLT |      8.7 secs |      1.8 secs |       10.86 |
| 🥉 3 |        NVIDIA cuSOLVER (Sp Chol) |        0 secs |       37 secs |     38.1472 |
|    4 |            Eigen::SimplicialLDLT |       37 secs |     0.95 secs |     93.8209 |
|    5 |             Eigen::SimplicialLLT |       38 secs |        1 secs |     39.0697 |
|    6 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       41 secs |      1.6 secs |         nan |
|    7 |              catamari::SparseLDL |       45 secs |     0.93 secs |     25.3056 |
|    8 |      Eigen::CholmodSupernodalLLT |       96 secs |       16 secs |     46.7459 |
|    9 |                  Eigen::SparseLU |  1.5e+02 secs |      1.9 secs |     37.5559 |
|   10 |         Eigen::CG\<IncompleteLUT\> |       41 secs |  1.4e+02 secs |         nan |
|   11 |                 Eigen::UmfPackLU |  2.5e+02 secs |  1.4e-06 secs |     46.7459 |

(BiCGSTAB/CG are capped at 200 iterations — the `nan`/large-residual rows at
k=3 reflect genuine non-convergence/divergence on this badly-scaled system,
not a bug; see the ⚠️ note above and "What are the systems being solved?"
below for why k=3 is inherently harder.)

Obviously [YMMV](https://www.google.com/search?q=YMMV), if you find something
interesting [let me know!](https://github.com/alecjacobson/sparse-solver-benchmark/issues).

## What about this other solver XYZ?

Please [submit a pull
request](https://github.com/alecjacobson/sparse-solver-benchmark/pulls) with a
wrapper for solver XYZ. The more the merrier. Ideas not yet covered here:
[Ginkgo](https://github.com/ginkgo-project/ginkgo) (GPU iterative/direct,
including a cuDSS backend), [AMGCL](https://github.com/ddemidov/amgcl),
[MUMPS](https://mumps-solver.org/), [PaStiX](https://gitlab.inria.fr/solverstack/pastix).

To add a solver: if it exposes an Eigen-compatible `Factor(Q)` /
`factor.solve(rhs)` interface, instantiate `solve<YourSolver>(...)` in the
`k` loop in `main.cpp`. Otherwise (like `catamari::SparseLDL` or the cuDSS/
cuSOLVER wrappers), add a dedicated function that calls `record(...)` with
the timing/residual, following the existing examples.

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
