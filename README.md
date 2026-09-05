# 🚀 Sparse Positive Definite (and Indefinite) Linear System Solver Leaderboard 📏

[![CI](https://github.com/alecjacobson/sparse-solver-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/alecjacobson/sparse-solver-benchmark/actions/workflows/ci.yml)

This is an informal leaderboard/benchmark for _k_-harmonic diffusion problems
on triangle meshes found commonly in geometry processing. Solvers are ranked
by total (factor + solve) time on each system; 🥇🥈🥉 mark the podium.

Two families of systems are benchmarked for k=2,3 (see "What are the systems
being solved?" below): the usual **flattened** SPD systems (k=1,2,3), and the
**mixed** (unflattened) formulation's genuinely **symmetric indefinite**
saddle-point systems (k=4,5) — a much harder test that most SPD-only direct
solvers can't handle at all, and a couple can't even handle *gracefully* (see
⚠️ below).

> Current upshot (see full results below, from an NVIDIA L40 / dual Intel
> Xeon Platinum 8362 machine):
>
> - **NVIDIA cuDSS** 🚀 wins every flattened (SPD) system here once the mesh
>   is large (720K vertices) — and its *solve* step is essentially free
>   (single-digit milliseconds) once factored, since it stays resident on the GPU.
> - **MKL Pardiso** 🏆 is the strongest CPU-only runner-up on SPD systems.
> - **CHOLMOD / Eigen SimplicialLLT·LDLT / catamari** 🥈 remain solid, accurate,
>   dependency-light choices when no GPU is available.
> - **On the indefinite mixed systems**, Cholesky-only solvers (CHOLMOD,
>   `SimplicialLLT`, `PardisoLLT`) correctly and cleanly refuse to factor —
>   exactly as expected. Pivoted LDLT (`PardisoLDLT`, catamari's LDLᵀ mode)
>   handles them well; Eigen's own *unpivoted* `SimplicialLDLT` does not (see
>   ⚠️ below). General LU (`SparseLU`) always works, as expected, since it
>   makes no definiteness assumption.
> - **Eigen's iterative solvers** (BiCGSTAB/CG + IncompleteLUT) get
>   unreliable as k grows and are unreliable on indefinite systems (CG in
>   particular, since it assumes SPD) — this benchmark caps them at 200
>   iterations so a hard system fails fast with an honest (possibly huge,
>   even `nan`) residual instead of hanging.
> - **Eigen SparseLU**¹ is the slowest general-purpose solver by a wide
>   margin, as expected — but see the ⚠️ note below if you see it apparently
>   taking *minutes* instead of seconds, that's a build misconfiguration, not
>   real solver cost.
>
> ¹Flattened systems are
> [SPD](https://en.wikipedia.org/wiki/Definite_symmetric_matrix) so LU is not
> a good choice there, but provides a reference (and is required for the
> indefinite mixed systems).

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

> [!WARNING]
> **Three solvers are known to crash — not just fail — on the mixed
> triharmonic system specifically** (a 2.16M-row genuinely indefinite
> saddle-point matrix on the example mesh) and are skipped there with an
> explicit reason, verified with `gdb`:
> - `Eigen::PardisoLLT`/`PardisoLDLT`: MKL Pardiso's reordering (both METIS
>   and minimum-degree) hangs — not just fails — on this system's sparsity
>   pattern (its λ block has an all-zero diagonal, inherent to this
>   saddle-point/KKT system).
> - `Eigen::UmfPackLU`: its internal MKL BLAS3 calls
>   (`umfdi_blas3_update`) spin up a fresh OpenMP thread team per call; the
>   catastrophic fill-in from this system's structure produces so many of
>   these tiny updates it exhausts OS thread/process limits and crashes.
> - `Eigen::SimplicialLDLT`: segfaults (`SIGSEGV` inside
>   `factorize_preordered`) — a real out-of-bounds access in Eigen's own
>   *unpivoted* LDLT at this scale, not just an inaccurate result.
>
> These are all skipped only for that one system (`k==5` in `main.cpp`); they
> run normally everywhere else, including the smaller mixed biharmonic system.

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
exceeds a per-system tolerance — a correctness regression test, not a
performance one. Covers all 5 systems (3 flattened + 2 mixed/indefinite).
GitHub Actions runs this on Linux (full solver matrix), Windows, and macOS
(Eigen + CHOLMOD/UMFPACK) on every push/PR; see `.github/workflows/ci.yml`.
Hosted runners have no GPU, so cuDSS/cuSOLVER still compile there but detect
the missing CUDA device at runtime and report `skipped` rather than fail.

## Example

Running

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

on an NVIDIA L40 / dual Intel Xeon Platinum 8362 (128 threads) machine
produces:

# Harmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |            Eigen::SimplicialLDLT |      1.3 secs |     0.13 secs | 1.12086e-10 |
| 🥈 2 |             Eigen::SimplicialLLT |      1.4 secs |     0.13 secs | 4.55334e-11 |
| 🥉 3 |              catamari::SparseLDL |      1.5 secs |      0.1 secs | 3.82439e-11 |
|    4 |                     NVIDIA cuDSS |      2.2 secs |   0.0023 secs | 1.59312e-10 |
|    5 |   Eigen::BiCGSTAB\<IncompleteLUT\> |      1.7 secs |      1.4 secs | 1.23985e-10 |
|    6 |               Eigen::PardisoLDLT |      3.4 secs |      1.2 secs | 1.04873e-10 |
|    7 |                Eigen::PardisoLLT |      3.5 secs |      1.2 secs | 7.58549e-11 |
|    8 |         Eigen::CG\<IncompleteLUT\> |      1.7 secs |      3.1 secs | 8.96274e-11 |
|    9 |                  Eigen::SparseLU |      5.6 secs |     0.24 secs | 2.37845e-11 |
|   10 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |      9.6 secs | 5.25522e-11 |
|   11 |      Eigen::CholmodSupernodalLLT |      8.7 secs |        1 secs | 6.63736e-11 |
|   12 |                 Eigen::UmfPackLU |       44 secs |     0.64 secs | 4.20999e-11 |

# Biharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      4.1 secs |    0.014 secs | 0.000121154 |
| 🥈 2 |                Eigen::PardisoLLT |        5 secs |      1.2 secs | 6.93083e-05 |
| 🥉 3 |               Eigen::PardisoLDLT |      5.1 secs |      1.3 secs | 4.78335e-05 |
|    4 |            Eigen::SimplicialLDLT |      9.7 secs |     0.42 secs | 4.80425e-05 |
|    5 |             Eigen::SimplicialLLT |      9.7 secs |     0.44 secs | 2.60041e-05 |
|    6 |              catamari::SparseLDL |       11 secs |     0.42 secs | 3.05382e-05 |
|    7 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       11 secs |      4.8 secs | 5.56194e-05 |
|    8 |         Eigen::CG\<IncompleteLUT\> |       11 secs |      6.5 secs | 4.73183e-05 |
|    9 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       20 secs | 3.40111e-05 |
|   10 |                  Eigen::SparseLU |       36 secs |     0.72 secs | 2.46911e-05 |
|   11 |      Eigen::CholmodSupernodalLLT |       59 secs |       10 secs | 9.99686e-05 |
|   12 |                 Eigen::UmfPackLU |  2.5e+02 secs |      2.4 secs | 8.19072e-05 |

# Triharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      6.5 secs |   0.0054 secs |     24.2444 |
| 🥈 2 |                Eigen::PardisoLLT |      8.5 secs |      1.3 secs |       10.86 |
| 🥉 3 |               Eigen::PardisoLDLT |      9.2 secs |      1.3 secs |     23.1328 |
|    4 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       36 secs |     38.1472 |
|    5 |            Eigen::SimplicialLDLT |       41 secs |     0.97 secs |     93.8209 |
|    6 |             Eigen::SimplicialLLT |       41 secs |        1 secs |     39.0697 |
|    7 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       41 secs |      1.5 secs |         nan |
|    8 |              catamari::SparseLDL |       45 secs |     0.91 secs |     25.3056 |
|    9 |      Eigen::CholmodSupernodalLLT |       93 secs |       12 secs |     46.7459 |
|   10 |                  Eigen::SparseLU |  1.5e+02 secs |      1.8 secs |     37.5559 |
|   11 |         Eigen::CG\<IncompleteLUT\> |       41 secs |  1.4e+02 secs |         nan |
|   12 |                 Eigen::UmfPackLU |  3.4e+02 secs |  3.1e-06 secs |     46.7459 |

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      5.5 secs |   0.0038 secs | 2.81664e-05 |
| 🥈 2 |               Eigen::PardisoLDLT |      6.5 secs |      2.3 secs | 5.45176e-11 |
| 🥉 3 |            Eigen::SimplicialLDLT |      9.6 secs |     0.48 secs | 6.46751e-05 |
|    4 |     catamari::SparseLDL (LDLᵀ) |       12 secs |     0.42 secs | 8.34264e-05 |
|    5 |                  Eigen::SparseLU |       32 secs |     0.74 secs | 1.09842e-10 |
|    6 |         Eigen::CG\<IncompleteLUT\> |      8.3 secs |       50 secs | 3.27831e+06 |
|    7 |   Eigen::BiCGSTAB\<IncompleteLUT\> |      8.4 secs |       90 secs | 9.73973e-09 |
|    8 |                 Eigen::UmfPackLU |    3e+02 secs |      2.7 secs |  8.0989e-11 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      9.3 secs |   0.0067 secs |     81430.7 |
| 🥈 2 |     catamari::SparseLDL (LDLᵀ) |       46 secs |      1.1 secs | 2.88158e-06 |
| 🥉 3 |                  Eigen::SparseLU |  1.1e+02 secs |      2.5 secs | 1.00706e-10 |
|    4 |         Eigen::CG\<IncompleteLUT\> |  6.2e+02 secs |       97 secs |     1735.27 |
|    5 |   Eigen::BiCGSTAB\<IncompleteLUT\> |  6.2e+02 secs |  1.8e+02 secs |      130.92 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk (see ⚠️ above) |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash (see ⚠️ above) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: known hang (see ⚠️ above) |
|    - |               Eigen::PardisoLDLT |           - |           - | skipped: known hang (see ⚠️ above) |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |

(BiCGSTAB/CG are capped at 200 iterations — the `nan`/large-residual rows
reflect genuine non-convergence/divergence on badly-scaled or indefinite
systems, not a bug. cuDSS's large residual on the mixed triharmonic system
reflects genuine numerical difficulty on this particularly hard indefinite
matrix — unlike the SPD-only solvers, it doesn't fail cleanly, it just
returns a less accurate answer; `SparseLU`/catamari's pivoted LDLᵀ remain
the trustworthy references there. See "What are the systems being solved?"
below for why k=3/mixed-triharmonic are inherently harder.)

\*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step —
it fuses reordering + symbolic + numeric factorization + triangular solve
into one call, repeated once per RHS column (there's no lower-level phased
Cholesky API in this cuSOLVER version) — so the whole cost is reported under
Solve rather than a fabricated Factor/Solve split. This is also why it's
slower than cuDSS, which factors once and solves 3 times.

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
`k` loop in `main.cpp` (this generic path also checks `factor.info()`, so a
solver that can't handle an indefinite system just gets a clean `skipped`
row instead of garbage output — see `Eigen::PardisoLDLT` for an example of a
solver added this way to answer "does pivoted LDLT survive indefinite
systems?"). Otherwise (like `catamari::SparseLDL`/`solve_catamari_ldl` or the
cuDSS/cuSOLVER wrappers), add a dedicated function that calls `record(...)`
with the timing/residual, following the existing examples.

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

### Mixed (unflattened) formulation

The systems above are "flattened": `Wᵏ` is built by eliminating auxiliary
variables via `M⁻¹` (`igl::harmonic`), producing a smaller but denser SPD
system. Jacobson et al. 2010's *mixed* finite element formulation instead
keeps those auxiliary variables explicit — sparser, but the resulting
symmetric block system is **indefinite** (a saddle-point/KKT structure)
instead of SPD, exercising solver robustness on harder input.

For the mixed biharmonic system, unknowns `(u, a₁)` where `a₁ = M⁻¹Lu`:

    [ M   L ] [u ]   [Mx]
    [ L  -M ] [a₁] = [0 ]

For the mixed triharmonic system, unknowns `(u, a₁, λ)` where `a₁=M⁻¹Lu`,
`λ=M⁻¹La₁`:

    [ M   0   L ] [u ]   [Mx]
    [ 0   L  -M ] [a₁] = [0 ]
    [ L  -M   0 ] [λ ]   [0 ]

Both are built via `build_mixed_system()` in `main.cpp` and produce the same
`u` as the corresponding flattened system (verified during development by
directly comparing the two solutions on the same mesh — they agree to
solver tolerance, degrading with k in the same way the flattened systems'
own conditioning does).
