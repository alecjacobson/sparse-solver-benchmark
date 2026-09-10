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
>   On the indefinite mixed systems it's more of a mixed bag (pun intended):
>   with iterative refinement enabled it's essentially exact on the mixed
>   biharmonic system, but still silently wrong on the harder mixed
>   triharmonic one — see ⚠️ below, this cuDSS version's pivoting for
>   symmetric indefinite matrices has a real gap (Bunch-Kaufman pivoting
>   isn't supported yet), confirmed directly with NVIDIA.
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
>   residual tolerance of 1e-11** (`‖b−Ax‖₂ < 1e-11·‖b‖₂`, via `setTolerance()`)
>   — the actual intended stopping criterion — with a deliberately generous
>   20000-iteration cap and a **10-minute wall-clock time limit** as two
>   independent safety nets against a system that never converges at all
>   (e.g. CG on indefinite input), not as the thing doing the deciding on
>   systems that do. The time limit exists because the iteration cap alone
>   isn't enough to bound wall time in practice: on the real dragon mesh, an
>   iterative solver on the badly-scaled flattened triharmonic system was
>   observed (via `gdb`, confirming it was genuinely still computing, not
>   hung) to run for multiple *hours* without reaching either the tolerance
>   or the iteration cap. Solving happens in time-boxed chunks (via
>   `solveWithGuess()`, continuing from the previous chunk's partial answer)
>   so a solver still short of convergence at the 10-minute mark is stopped
>   and reports its best-effort accuracy at that point — marked with a `†`
>   in the leaderboard tables — rather than either running unboundedly or
>   being cut off with no time guarantee at all. `warp_bench/` uses the
>   identical relative-L2 tolerance formula and the same 10-minute limit
>   (chunked the same way), verified by reading both libraries' source, so a
>   comparison between them is apples-to-apples. `1e-11` (not the more
>   obvious-looking `1e-7`) because the relative-L2 criterion and the L∞
>   absolute residual this benchmark displays aren't the same quantity, and
>   on a large vector with unevenly-distributed error the gap between them
>   can be substantial — see "How is the L∞ norm computed?" below for the
>   measured numbers that drove this choice.
> - **Eigen SparseLU**¹ is the slowest general-purpose solver by a wide
>   margin, as expected — but see the ⚠️ note below if you see it apparently
>   taking *minutes* instead of seconds, that's a build misconfiguration, not
>   real solver cost.
> - **NVIDIA Warp**'s `warp.optim.linear` solvers (`cg`/`cr`/`bicgstab`/
>   `gmres`, timed separately via [`warp_bench/`](warp_bench/)) run to the
>   same `1e-11` relative tolerance as Eigen's iterative solvers, but with
>   only Jacobi preconditioning (the strongest Warp currently offers — Eigen
>   uses incomplete-LU, meaningfully stronger) and the same 20000-iteration/
>   10-minute safety nets. At this mesh's real scale (360K-2.16M rows), that
>   combination usually isn't enough to actually converge, so most `warp::*`
>   rows are correctly caught by the did-not-actually-succeed check (see
>   below) and don't appear ranked at all — a genuine, if slightly deflating,
>   result: naive Warp iterative solvers aren't a free win over tuned CPU/GPU
>   direct solvers at this problem size without a better preconditioner.
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
> **cuDSS "succeeds" (`CUDSS_STATUS_SUCCESS`, no error) on the mixed
> triharmonic system while silently returning a numerically useless
> answer** (L∞ residual in the tens of thousands — see the table above,
> where it's correctly caught by this benchmark's did-not-actually-succeed
> check rather than shown as a real result). Root cause, confirmed directly
> with NVIDIA's cuDSS team: this cuDSS version's only pivoting strategy for
> `CUDSS_MTYPE_SYMMETRIC` matrices is `CUDSS_PIVOT_DIAGONAL` (diagonal-only
> search) — `CUDSS_PIVOT_BUNCH_KAUFMAN`, the safe block-pivoting strategy
> NASOQ uses (see above), is explicitly "reserved for future, not supported
> yet" in `cudss_data_types.h`. Diagonal-only pivoting can't find a safe
> pivot in a block that's structurally all-zero on the diagonal (this
> system's λ block, same root cause as the Pardiso/NASOQ issues above), so
> the factorization silently produces garbage rather than erroring out.
> NVIDIA confirmed Bunch-Kaufman pivoting is coming in cuDSS's next release,
> and suggested enabling iterative refinement (`CUDSS_CONFIG_IR_N_STEPS`) as
> a workaround in the meantime — this benchmark now sets it to `2` for
> `CUDSS_MTYPE_SYMMETRIC` solves (see `solve_cudss()`). It genuinely fixes
> the **mixed biharmonic** system (residual improved from `2.6e-05` to
> `5.5e-11`, essentially machine precision) but does **not** fix the mixed
> triharmonic one — empirically (`IR_N_STEPS` from 2 up to 20 tested on a
> small synthetic case) the result is at best inconsistent, ranging from
> ~1e-6 to ~5.9 residual run-to-run on the *identical* input, matching
> NVIDIA's own statement that enabling their deterministic mode makes a
> wrong answer *consistently* wrong rather than fixing its accuracy — so
> this benchmark doesn't rely on IR alone to detect success here; the
> did-not-actually-succeed check above is the real safety net.

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
produces (tables below predate `kMaxIterativeIterations` being raised from
200 to 20000, and predate the 10-minute wall-clock time limit described
above -- a full rerun at the current settings is a pending follow-up, since
it substantially increases total run time (an iterative solver was observed
to legitimately run for hours on the hardest system before the time limit
was added); the numbers are still directionally correct, but some
iterative-solver rows may now converge further/differently, or show a `†`
marker if they hit the new time limit instead):

# Harmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |            Eigen::SimplicialLDLT |      1.2 secs |     0.13 secs | 1.12086e-10 |
| 🥈 2 |             Eigen::SimplicialLLT |      1.2 secs |     0.13 secs | 4.55334e-11 |
| 🥉 3 |              catamari::SparseLDL |      1.4 secs |     0.11 secs | 3.82439e-11 |
|    4 |   Eigen::BiCGSTAB<IncompleteLUT> |      1.6 secs |     0.76 secs | 1.37219e-10 |
|    5 |                     NVIDIA cuDSS |      2.4 secs |   0.0026 secs | 1.59312e-10 |
|    6 |         Eigen::CG<IncompleteLUT> |      1.6 secs |     0.84 secs |  1.4579e-05 |
|    7 |                        NASOQ LBL |      2.6 secs |     0.18 secs | 1.09436e-10 |
|    8 |                Eigen::PardisoLLT |        3 secs |        1 secs | 7.58549e-11 |
|    9 |               Eigen::PardisoLDLT |      3.2 secs |      1.1 secs | 1.04873e-10 |
|   10 |                  Eigen::SparseLU |        5 secs |     0.21 secs | 2.37845e-11 |
|   11 |      Eigen::CholmodSupernodalLLT |      6.3 secs |     0.53 secs | 6.63736e-11 |
|   12 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |      9.6 secs | 5.25522e-11 |
|   13 |                 Eigen::UmfPackLU |       30 secs |     0.65 secs | 4.20999e-11 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: L∞ residual 14.51 is 6.1e+11 x the best solver's (2.378e-11) on this system |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: L∞ residual 33.35 is 1.4e+12 x the best solver's (2.378e-11) on this system |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: L∞ residual 57.67 is 2.42e+12 x the best solver's (2.378e-11) on this system |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: L∞ residual 31.49 is 1.32e+12 x the best solver's (2.378e-11) on this system |

*(fused): this solver's API has no separate factor step; the whole
 analysis+factor+solve cost is reported under Solve instead.

# Biharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |      Eigen::CholmodSupernodalLLT |      2.2 secs |     0.25 secs | 9.99686e-05 |
| 🥈 2 |                     NVIDIA cuDSS |        4 secs |    0.018 secs | 6.79479e-05 |
| 🥉 3 |                        NASOQ LBL |      5.3 secs |     0.29 secs | 0.000148578 |
|    4 |                Eigen::PardisoLLT |      4.9 secs |     0.83 secs | 6.93083e-05 |
|    5 |               Eigen::PardisoLDLT |        5 secs |      1.1 secs | 4.78335e-05 |
|    6 |                 Eigen::UmfPackLU |      5.4 secs |        2 secs | 8.19072e-05 |
|    7 |             Eigen::SimplicialLLT |      9.1 secs |     0.39 secs | 2.60041e-05 |
|    8 |            Eigen::SimplicialLDLT |      9.2 secs |     0.53 secs | 4.80425e-05 |
|    9 |              catamari::SparseLDL |       11 secs |     0.42 secs | 3.05382e-05 |
|   10 |         Eigen::CG<IncompleteLUT> |       12 secs |      2.3 secs |  5.3789e-05 |
|   11 |   Eigen::BiCGSTAB<IncompleteLUT> |       12 secs |      2.5 secs | 7.33713e-05 |
|   12 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       20 secs | 3.40111e-05 |
|   13 |                  Eigen::SparseLU |       34 secs |     0.74 secs | 2.46911e-05 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: L∞ residual 1.963e+04 is 7.95e+08 x the best solver's (2.469e-05) on this system |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: L∞ residual 529.9 is 2.15e+07 x the best solver's (2.469e-05) on this system |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: L∞ residual 552.7 is 2.24e+07 x the best solver's (2.469e-05) on this system |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: L∞ residual 98.38 is 3.98e+06 x the best solver's (2.469e-05) on this system |

*(fused): this solver's API has no separate factor step; the whole
 analysis+factor+solve cost is reported under Solve instead.

# Triharmonic

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                      warp::gmres |  0.00046 secs |     0.68 secs |     124.288 |
| 🥈 2 |                     NVIDIA cuDSS |      6.5 secs |   0.0054 secs |     38.0908 |
| 🥉 3 |      Eigen::CholmodSupernodalLLT |      9.6 secs |     0.32 secs |     6.77554 |
|    4 |                Eigen::PardisoLLT |      8.9 secs |      1.3 secs |       10.86 |
|    5 |               Eigen::PardisoLDLT |        9 secs |      1.2 secs |     23.1328 |
|    6 |                        NASOQ LBL |       11 secs |     0.44 secs |      61.402 |
|    7 |                 Eigen::UmfPackLU |       14 secs |  7.2e-07 secs |     6.77554 |
|    8 |             Eigen::SimplicialLLT |       35 secs |     0.88 secs |     39.0697 |
|    9 |            Eigen::SimplicialLDLT |       35 secs |     0.91 secs |     93.8209 |
|   10 |        NVIDIA cuSOLVER (Sp Chol) |     (fused)* |       36 secs |     38.1472 |
|   11 |              catamari::SparseLDL |       45 secs |     0.84 secs |     25.3056 |
|   12 |                  Eigen::SparseLU |  1.5e+02 secs |      1.8 secs |     37.5559 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: L∞ residual 1.243e+09 is 1.83e+08 x the best solver's (6.776) on this system |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: L∞ residual 1.748e+06 is 2.58e+05 x the best solver's (6.776) on this system |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: L∞ residual 1.988e+07 is 2.93e+06 x the best solver's (6.776) on this system |

*(fused): this solver's API has no separate factor step; the whole
 analysis+factor+solve cost is reported under Solve instead.

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |                     NVIDIA cuDSS |      5.4 secs |    0.014 secs | 5.47711e-11 |
| 🥈 2 |                        NASOQ LBL |      6.9 secs |      0.5 secs | 1.68701e-05 |
| 🥉 3 |               Eigen::PardisoLDLT |      6.1 secs |      2.3 secs | 5.45176e-11 |
|    4 |            Eigen::SimplicialLDLT |      8.9 secs |     0.45 secs | 6.46751e-05 |
|    5 |       catamari::SparseLDL (LDLᵀ) |       12 secs |     0.43 secs | 8.34264e-05 |
|    6 |                 Eigen::UmfPackLU |       11 secs |      2.2 secs |  8.0989e-11 |
|    7 |                  Eigen::SparseLU |       30 secs |      0.8 secs | 1.09842e-10 |
|    8 |   Eigen::BiCGSTAB<IncompleteLUT> |      8.3 secs |       67 secs | 3.33011e-08 |
|    - |         Eigen::CG<IncompleteLUT> |           - |           - | skipped: did not actually succeed: L∞ residual 3.278e+06 is 6.01e+16 x the best solver's (5.452e-11) on this system |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: L∞ residual 5067 is 9.29e+13 x the best solver's (5.452e-11) on this system |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: L∞ residual 4.021e+14 is 7.38e+24 x the best solver's (5.452e-11) on this system |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: L∞ residual 5889 is 1.08e+14 x the best solver's (5.452e-11) on this system |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: L∞ residual 1.511 is 2.77e+10 x the best solver's (5.452e-11) on this system |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve |     L∞ norm |
|-----:|--------------------------------:|------------:|------------:|------------:|
| 🥇 1 |       catamari::SparseLDL (LDLᵀ) |       46 secs |        1 secs | 2.88158e-06 |
| 🥈 2 |                  Eigen::SparseLU |    1e+02 secs |      1.7 secs | 1.00706e-10 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: L∞ residual 130.9 is 1.3e+12 x the best solver's (1.007e-10) on this system |
|    - |         Eigen::CG<IncompleteLUT> |           - |           - | skipped: did not actually succeed: L∞ residual 1735 is 1.72e+13 x the best solver's (1.007e-10) on this system |
|    - |                     NVIDIA cuDSS |           - |           - | skipped: did not actually succeed: L∞ residual 8.284e+04 is 8.23e+14 x the best solver's (1.007e-10) on this system |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: L∞ residual 8.638e+07 is 8.58e+17 x the best solver's (1.007e-10) on this system |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: L∞ residual 1.378 is 1.37e+10 x the best solver's (1.007e-10) on this system |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: L∞ residual 2.414 is 2.4e+10 x the best solver's (1.007e-10) on this system |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: L∞ residual 1.392 is 1.38e+10 x the best solver's (1.007e-10) on this system |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk: excessive MKL thread churn on this system's fill-in |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale |
|    - |                        NASOQ LBL |           - |           - | skipped: known crash: SIGSEGV in libmetis genmmd/mmdelm via NASOQ's symbolic_analysis_lin_solve on this system's sparsity pattern |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - |               Eigen::PardisoLDLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |

(Every solver here — direct or iterative, C++ or the `warp::*` rows from
[`warp_bench/`](warp_bench/), see `--dump-matrices` above — is checked
against every other solver's result on the same system: if a row's L∞
residual is enormous relative to the best actually achieved there, it's
reclassified as "did not actually succeed" and moved to the skipped section
instead of ranking with a misleadingly-real-looking number, even if the
solver itself reported success. This catches genuine silent failures like
cuDSS's on the mixed triharmonic system (see ⚠️ below) as well as iterative
solvers that hit their iteration cap without reaching tolerance — both
Eigen's BiCGSTAB/CG and Warp's `cg`/`cr`/`bicgstab`/`gmres` now run to a real
relative-L2 tolerance of `1e-11` (not just an iteration cap; see the ⚠️ note
below) but can still fail to reach it on a large-enough or hard-enough
system, which is exactly what most of the `warp::*` rows above show: at this
mesh's scale (360K-2.16M rows) and with only Jacobi preconditioning (the
strongest Warp currently offers; Eigen's rows use incomplete-LU, a
meaningfully stronger preconditioner Warp doesn't have), 200 iterations
usually isn't enough to converge, so most Warp rows are correctly caught by
this check rather than ranking on raw iteration speed alone.)

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

### How is the L∞ norm computed, and why is it so large for k=3?

The `L∞ norm` column is the residual's absolute max-component norm,
`‖rhs − Q·U‖∞`, computed with `(rhs-Q*U).array().abs().maxCoeff()` in
`main.cpp` — the single largest absolute error across *every* row and *every*
right-hand-side column at once (`rhs = M*x` has 3 columns, the mesh's own
x/y/z coordinates), not a per-column or relative quantity. It's **absolute**,
not normalized by `‖rhs‖` or `‖Q‖` — deliberately, so it's the same metric
for every solver regardless of algorithm, and comparable to the per-system
`g_check_tol[]` correctness thresholds used by `--check`/CTest. It is *not*
comparable across different k or between the flattened/mixed systems, whose
residuals live at genuinely different scales (see below) — only within the
same system, across solvers, is it apples-to-apples.

k=3's residuals (tens to hundreds, vs. ~1e-10 for k=1) look alarming at a
glance, but they're a direct, unavoidable consequence of `Wᵏ⁺¹ = Wᵏ M⁻¹ L`
being applied recursively: each application inverts the mass matrix `M`
again, and `M`'s entries scale with triangle area, so on a mesh with small
triangles `M⁻¹`'s entries (and thus `Q`'s) can become enormous. Measured
directly (via `--dump-matrices`, see above) on this README's two example
meshes' own k=3 (flattened Triharmonic) systems:

| Mesh | `Q`'s largest entry (abs) | `rhs`'s largest entry (abs) | Absolute residual | Residual ÷ `‖rhs‖∞` | Residual ÷ `‖Q‖∞` |
|---|---:|---:|---:|---:|---:|
| synthetic 20×20 grid | 1.8×10⁷ | 2.6×10⁻³ | 1.8×10⁻⁸ | 6.9×10⁻⁶ | **1.0×10⁻¹⁵** |
| `xyzrgb_dragon-720K.ply` | **3.0×10¹⁵** | 1.3×10² | 47.9 | 0.37 | **1.6×10⁻¹⁴** |

On the dragon mesh, `Q` has individual entries as large as 3×10¹⁵ — right at
the edge of what a double (≈15-16 significant decimal digits) can represent
alongside `Q`'s much smaller entries (as small as 2.8×10⁻⁸) in the same
matrix. A residual of "47.9" looks large next to `rhs`'s own scale (~131),
but relative to `Q`'s own scale it's ~1.6×10⁻¹⁴ — right at double-precision
machine epsilon (~2.2×10⁻¹⁶) times `Q`'s magnitude. In other words: this is
the floating-point roundoff floor imposed by `Q` containing such
astronomically large entries, not an inaccurate solve — a mathematically
exact solve, computed in double precision, would show essentially the same
absolute residual, because the roundoff in forming `Q·U` itself scales with
`Q`'s magnitude. `g_check_tol[3] = 1e3` (see `main.cpp`) was tuned
empirically with exactly this in mind; the other flattened systems'
comparatively tiny tolerances (`1e-4`, `1e-1`) reflect `Q` staying at a much
saner scale for k=1,2.

#### Why is `kIterativeTolerance` 1e-11, not the more obvious-looking 1e-7?

The iterative solvers' internal stopping criterion is a *relative L2*
quantity (`‖b−Ax‖₂ < tol·‖b‖₂`), but the `L∞ norm` column above is an
*absolute, max-component* quantity — not the same number, and on a large
vector they can differ substantially if the residual error isn't spread
evenly across components. Measured directly on the dragon mesh's Harmonic
(k=1) system, `‖b‖₂` (3824) is ~30x `‖b‖∞` (130); at `tol=1e-7`, Warp's
`cg` genuinely met that relative-L2 target (verified: its own
internally-tracked residual and an independently recomputed one agreed
exactly, so this isn't drift) while the resulting *absolute* L∞ residual
was only `1.5e-4` — five orders of magnitude looser than `1e-7` might
suggest. `tol=1e-10` narrowed the gap but still landed on the wrong side
of a clean `1e-7` L∞ target for some solvers/columns (`1.5e-7`-`2.8e-7`
across `cg`/`cr`/`bicgstab`/`gmres`); `tol=1e-11` gives solid margin
(measured `9.4e-9`-`3.2e-8`) at negligible extra iteration cost (CG's
iteration count grows only mildly per decade of tolerance on a
well-conditioned system like this one — 730 iterations at `1e-11` vs. 549
at `1e-7`, both well under a second).

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
