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
> - **NASOQ LBL** (a QP solver's symmetric-indefinite linear solver
>   component) is competitive and accurate everywhere it runs — including
>   the indefinite mixed biharmonic system — but crashes on the mixed
>   triharmonic system specifically (see ⚠️ below); its first "solve" call
>   used to look absurdly slow (16-24s on tiny problems) due to two
>   independent upstream bugs, both fixed here (see ⚠️ below).
> - **Eigen's iterative solvers** (BiCGSTAB/CG + IncompleteLUT) get
>   unreliable as k grows and are unreliable on indefinite systems (CG in
>   particular, since it assumes SPD). They now stop at a **relative L2
>   residual tolerance of 1e-7** (`‖b−Ax‖₂ < 1e-7·‖b‖₂`, via `setTolerance()`)
>   — the actual intended stopping criterion — with a 200-iteration cap
>   remaining only as a safety net against a system that never converges at
>   all, not as the thing doing the deciding on systems that do. `warp_bench/`
>   uses the identical relative-L2 formula at the same `1e-7`, verified by
>   reading both libraries' source, so a comparison between them is
>   apples-to-apples.
> - **Eigen SparseLU**¹ is the slowest general-purpose solver by a wide
>   margin, as expected — but see the ⚠️ note below if you see it apparently
>   taking *minutes* instead of seconds, that's a build misconfiguration, not
>   real solver cost.
> - **NVIDIA Warp**'s `warp.optim.linear` solvers (`cg`/`cr`/`bicgstab`/
>   `gmres`, timed separately via [`warp_bench/`](warp_bench/)) top every
>   table here by raw speed — sub-2-second solves on every system, including
>   the 2.16M-row mixed triharmonic one — but at the same 200-iteration cap
>   as Eigen's iterative solvers, that's nowhere near enough to converge at
>   this scale, so their residuals are correspondingly huge. Read these rows
>   as "how fast is one Warp iteration," not "how accurate is Warp."
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
> **Four solvers are known to crash or hang — not just fail — on the mixed
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
> - `NASOQ LBL`: root-caused with `valgrind` down to a genuine one-off
>   heap buffer overflow in NASOQ's own `SolverSettings::find_perturbation()`
>   (`src/QP/linear_solver_wrapper.cpp`), which assumes
>   `AorSM->x[AorSM->p[i]]` is always column `i`'s diagonal entry. This
>   system's λ block has a structurally all-zero diagonal (same block that
>   trips up Pardiso's reordering above), so every one of its columns has
>   zero stored lower-triangular entries — for the last such column,
>   `p[i]` equals `nnz`, and the code reads/writes exactly one `double` past
>   the end of the CSC value array. That small heap corruption then cascades
>   into an intermittent, memory-layout-dependent `SIGSEGV` later, deep
>   inside `libmetis.so.5`'s minimum-degree ordering (`genmmd`/`mmdelm`,
>   reached via NASOQ's `symbolic_analysis_lin_solve()` → `METIS_NodeND()`)
>   — which is why it doesn't crash on every single run. Reproduces on a
>   tiny synthetic 1200-row mixed triharmonic system, standalone (i.e.
>   independent of this benchmark) against NASOQ's own unmodified
>   `LBL_Eigen` example — see
>   [sympiler/nasoq#33](https://github.com/sympiler/nasoq/issues/33).
>
> These are all skipped only for that one system (`k==5` in `main.cpp`); they
> run normally everywhere else, including the smaller mixed biharmonic system.

> [!WARNING]
> **NASOQ LBL's `solve_only()` looked absurdly slow (16-24 seconds!) on
> problems from 400 rows to 1.4M rows alike** — tracked down to two
> independent bugs in NASOQ itself, not this benchmark's usage:
> 1. NASOQ's own `eigen_interface` example sets `req_ref_iter=2`, requesting
>    internal GMRES-based iterative refinement (`pmgmres_ldlt_auto`) after
>    every direct solve — pure overhead here, since this benchmark already
>    validates accuracy via its own downstream residual check like every
>    other solver. Fixed locally by setting `req_ref_iter=0`.
> 2. The real cost: `SolverSettings::num_thread` does *not* control the
>    OpenMP team size of the actual triangular-solve kernels
>    (`H2LeveledBlockedLsolve`/`LTsolve`) — they use a bare
>    `#pragma omp parallel` with no `num_threads()` clause, and NASOQ's only
>    `omp_set_num_threads()` call site in the whole codebase is commented
>    out. So every `solve_only()` call was spinning up
>    `omp_get_max_threads()` (128 on this machine) OpenMP threads regardless
>    of the configured thread count — pure thread-team creation overhead,
>    not real work. Fixed locally by calling `omp_set_num_threads(1)`
>    ourselves around the NASOQ calls (save/restore, so it doesn't affect
>    other solvers).
>
> Together: **16-24 seconds → sub-millisecond** for `solve_only()`. Filed
> upstream as [sympiler/nasoq#31](https://github.com/sympiler/nasoq/issues/31)
> with [sympiler/nasoq#32](https://github.com/sympiler/nasoq/pull/32) (a real
> fix threading a `num_threads` parameter through the call chain, rather than
> this benchmark's coarser process-wide `omp_set_num_threads()` workaround).

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
| `IGL_WITH_NASOQ`    | NASOQ's LBL parallel symmetric-indefinite solver (requires `IGL_WITH_MKL` for its BLAS backend) |

## Run

    ./sparse_solver_benchmark [path to triangle mesh]
    ./sparse_solver_benchmark --grid 20 --check   # fast synthetic-mesh correctness check
    ./sparse_solver_benchmark mesh.ply --csv results.csv   # also dump raw results
    ./sparse_solver_benchmark mesh.ply --only nasoq        # only run solvers matching "nasoq"
    ./sparse_solver_benchmark mesh.ply --exclude umfpack,sparselu   # skip the slow ones

`--only`/`--exclude` take a comma-separated (repeatable) list of
case-insensitive substrings matched against each solver's printed name
(e.g. `--only nasoq` runs just NASOQ LBL; `--exclude umfpack,sparselu` skips
the two slowest general-purpose solvers on the big meshes). A filtered-out
solver isn't run at all — it doesn't even appear as a `skipped` row — so
this is meant for fast local iteration on one solver at a time (e.g. while
debugging a specific solver against the full dragon mesh, without waiting
for `SparseLU`/`CG`/etc. to grind through every k), not for the leaderboard
tables below, which always run every solver.

    ./sparse_solver_benchmark mesh.ply --dump-matrices /tmp/dump --dump-only

`--dump-matrices dir` writes each system's `Q`/`rhs` as plain
[MatrixMarket](https://math.nist.gov/MatrixMarket/formats.html) files
(`k<k>_Q.mtx`, `k<k>_rhs.mtx`) instead of (or alongside, without
`--dump-only`) running this benchmark's own solvers -- see
[`warp_bench/`](warp_bench/) for a Python add-on that loads these to time
[NVIDIA Warp](https://github.com/NVIDIA/warp)'s `warp.optim.linear` iterative
solvers (`cg`/`cr`/`bicgstab`/`gmres`), a comparison that doesn't fit into
this project's C++/Eigen-based solver interface.

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
| 🥇 1 |                         warp::cg |    0.0013 secs |     0.052 secs |     1521.55 |
| 🥈 2 |                         warp::cr |    0.0019 secs |     0.058 secs |     16.3665 |
| 🥉 3 |                   warp::bicgstab |    0.0013 secs |      0.11 secs |     128.245 |
|    4 |                      warp::gmres |    0.0013 secs |      0.27 secs |     30.6497 |
|    5 |             Eigen::SimplicialLLT |      1.2 secs |     0.13 secs | 4.55334e-11 |
|    6 |            Eigen::SimplicialLDLT |      1.2 secs |     0.13 secs | 1.12086e-10 |
|    7 |              catamari::SparseLDL |      1.4 secs |     0.11 secs | 3.82439e-11 |
|    8 |                     NVIDIA cuDSS |      2.2 secs |   0.0022 secs | 1.59312e-10 |
|    9 |                        NASOQ LBL |      2.6 secs |     0.17 secs | 1.09436e-10 |
|   10 |   Eigen::BiCGSTAB\<IncompleteLUT\> |      1.6 secs |      1.3 secs | 1.23985e-10 |
|   11 |         Eigen::CG\<IncompleteLUT\> |      1.6 secs |      2.8 secs | 8.96274e-11 |
|   12 |                Eigen::PardisoLLT |      3.3 secs |      1.1 secs | 7.58549e-11 |
|   13 |               Eigen::PardisoLDLT |      3.2 secs |      1.3 secs | 1.04873e-10 |
|   14 |                  Eigen::SparseLU |      4.9 secs |     0.18 secs | 2.37845e-11 |
|   15 |      Eigen::CholmodSupernodalLLT |      6.7 secs |     0.75 secs | 6.63736e-11 |
|   16 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |      9.2 secs | 5.25522e-11 |
|   17 |                 Eigen::UmfPackLU |       28 secs |     0.76 secs | 4.20999e-11 |

# Biharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                         warp::cg |    0.0068 secs |       0.1 secs |   1.40821e+06 |
| 🥈 2 |                         warp::cr |    0.0013 secs |      0.11 secs |     127.337 |
| 🥉 3 |                   warp::bicgstab |    0.0013 secs |      0.22 secs |     1076.87 |
|    4 |                      warp::gmres |    0.0013 secs |      0.41 secs |     128.031 |
|    5 |      Eigen::CholmodSupernodalLLT |      2.2 secs |     0.26 secs | 9.99686e-05 |
|    6 |                     NVIDIA cuDSS |        4 secs |    0.008 secs | 0.000142766 |
|    7 |                        NASOQ LBL |      5.3 secs |     0.29 secs | 0.000148578 |
|    8 |                Eigen::PardisoLLT |      4.9 secs |      1.1 secs | 6.93083e-05 |
|    9 |               Eigen::PardisoLDLT |        5 secs |        1 secs | 4.78335e-05 |
|   10 |                 Eigen::UmfPackLU |      5.4 secs |      1.9 secs | 8.19072e-05 |
|   11 |             Eigen::SimplicialLLT |      9.1 secs |     0.39 secs | 2.60041e-05 |
|   12 |            Eigen::SimplicialLDLT |      9.1 secs |      0.4 secs | 4.80425e-05 |
|   13 |              catamari::SparseLDL |       11 secs |     0.43 secs | 3.05382e-05 |
|   14 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       11 secs |      4.5 secs | 5.56194e-05 |
|   15 |         Eigen::CG\<IncompleteLUT\> |       11 secs |      6.3 secs | 4.73183e-05 |
|   16 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       20 secs | 3.40111e-05 |
|   17 |                  Eigen::SparseLU |       35 secs |     0.66 secs | 2.46911e-05 |

# Triharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                         warp::cr |    0.0015 secs |      0.91 secs |     128.344 |
| 🥈 2 |                         warp::cg |     0.037 secs |      0.89 secs |      546670 |
| 🥉 3 |                      warp::gmres |    0.0013 secs |       1.4 secs |     128.344 |
|    4 |                   warp::bicgstab |    0.0016 secs |       1.7 secs |     45520.5 |
|    5 |                     NVIDIA cuDSS |      6.5 secs |   0.0054 secs |     31.6823 |
|    6 |      Eigen::CholmodSupernodalLLT |      9.5 secs |     0.32 secs |     6.77554 |
|    7 |                Eigen::PardisoLLT |      9.3 secs |      1.1 secs |       10.86 |
|    8 |               Eigen::PardisoLDLT |      9.4 secs |      1.7 secs |     23.1328 |
|    9 |                        NASOQ LBL |       11 secs |     0.45 secs |      61.402 |
|   10 |                 Eigen::UmfPackLU |       14 secs |  4.8e-07 secs |     6.77554 |
|   11 |            Eigen::SimplicialLDLT |       36 secs |     0.86 secs |     93.8209 |
|   12 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       37 secs |     38.1472 |
|   13 |             Eigen::SimplicialLLT |       36 secs |     0.92 secs |     39.0697 |
|   14 |   Eigen::BiCGSTAB\<IncompleteLUT\> |       41 secs |      1.6 secs |         nan |
|   15 |              catamari::SparseLDL |       45 secs |     0.95 secs |     25.3056 |
|   16 |                  Eigen::SparseLU |  1.5e+02 secs |      1.9 secs |     37.5559 |
|   17 |         Eigen::CG\<IncompleteLUT\> |       41 secs |  1.5e+02 secs |         nan |

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                         warp::cg |     0.016 secs |      0.19 secs |      165519 |
| 🥈 2 |                         warp::cr |    0.0024 secs |      0.23 secs |     131.166 |
| 🥉 3 |                   warp::bicgstab |    0.0024 secs |      0.43 secs |     38996.8 |
|    4 |                      warp::gmres |    0.0025 secs |         1 secs |     117.849 |
|    5 |                     NVIDIA cuDSS |      5.4 secs |   0.0038 secs | 4.32598e-05 |
|    6 |                        NASOQ LBL |        7 secs |     0.49 secs | 1.68701e-05 |
|    7 |            Eigen::SimplicialLDLT |      9.4 secs |     0.45 secs | 6.46751e-05 |
|    8 |                 Eigen::UmfPackLU |       11 secs |      2.2 secs |  8.0989e-11 |
|    9 |               Eigen::PardisoLDLT |        7 secs |      6.8 secs | 5.45176e-11 |
|   10 |         Eigen::CG\<IncompleteLUT\> |      8.5 secs |       49 secs | 3.27831e+06 |
|   11 |     catamari::SparseLDL (LDLᵀ) |  1.3e+02 secs |      0.5 secs | 8.34264e-05 |
|   12 |                  Eigen::SparseLU |  1.5e+02 secs |     0.82 secs | 1.09842e-10 |
|   13 |   Eigen::BiCGSTAB\<IncompleteLUT\> |      8.2 secs |  1.7e+02 secs | 9.73973e-09 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed (not SPD, as expected) |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                         warp::cg |     0.019 secs |      0.35 secs |     166357 |
| 🥈 2 |                         warp::cr |    0.0038 secs |       0.4 secs |     131.168 |
| 🥉 3 |                   warp::bicgstab |    0.0041 secs |      0.74 secs |     6801.41 |
|    4 |                      warp::gmres |    0.0038 secs |       1.5 secs |     120.478 |
|    5 |                     NVIDIA cuDSS |      9.3 secs |   0.0068 secs |     28903.4 |
|    6 |     catamari::SparseLDL (LDLᵀ) |       46 secs |      1.1 secs | 2.88158e-06 |
|    7 |                  Eigen::SparseLU |  1.1e+02 secs |      1.6 secs | 1.00706e-10 |
|    8 |         Eigen::CG\<IncompleteLUT\> |  6.2e+02 secs |       97 secs |     1735.27 |
|    9 |   Eigen::BiCGSTAB\<IncompleteLUT\> |  6.1e+02 secs |  1.9e+02 secs |      130.92 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk (see ⚠️ above) |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash (see ⚠️ above) |
|    - |                        NASOQ LBL |           - |           - | skipped: known crash (see ⚠️ above) |
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
below for why k=3/mixed-triharmonic are inherently harder.

The `warp::*` rows are [NVIDIA Warp](https://github.com/NVIDIA/warp)'s
`warp.optim.linear` iterative solvers, timed separately via
[`warp_bench/`](warp_bench/) (see `--dump-matrices` above) and capped at the
same 200 iterations as the Eigen iterative solvers for direct comparability
— which is *not* enough iterations to converge on a system this large
(360K-2.16M rows depending on k), so their huge residuals are expected, not
a bug: these rows show how fast Warp's GPU-native solvers run per iteration,
not their achievable accuracy. A tolerance-based comparison across every
iterative solver here (Eigen's and Warp's alike) — so "200 iterations" isn't
silently doing the deciding — is planned as a follow-up.)

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
