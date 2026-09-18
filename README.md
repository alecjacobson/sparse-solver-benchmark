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

> Current upshot (qualitative conclusions below hold across machines; see
> [`leaderboards/`](leaderboards/) for full per-machine results tables, each
> file named after the machine it was measured on with its specs at the top):
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
> - **Accuracy is now reported as componentwise relative backward error**
>   (LAPACK's BERR — see "How is accuracy measured?" below), not an absolute
>   residual. Every direct solver here lands within a couple of orders of
>   magnitude of machine epsilon (~1e-13 to 1e-16) *regardless of k* — the
>   old absolute-residual metric made k=3/mixed-triharmonic look like a
>   fundamentally less accurate solve than k=1, when in fact it's the same
>   quality of double-precision solve throughout, just on a matrix with much
>   larger entries. Iterative solvers (Eigen's BiCGSTAB/CG, Warp's
>   `cg`/`cr`/`bicgstab`/`gmres`) are driven toward a common external target
>   of **backward error < 1e-8**, checked ourselves after every chunk of a
>   time-boxed solve (not each library's own internal convergence test —
>   comparing solvers fairly requires one shared target, not each library's
>   private notion of "converged") — with a 1,000,000-iteration cap, stall
>   detection (stop if backward error hasn't improved meaningfully over
>   several consecutive chunks — a genuine "not making progress" signal, not
>   an arbitrary budget), and a **10-minute wall-clock time limit** as
>   independent safety nets against a system that never converges at all
>   (e.g. CG on indefinite input), not as the thing doing the deciding on
>   systems that do. The time limit exists because an iteration cap alone
>   isn't enough to bound wall time in practice: on the real dragon mesh, an
>   iterative solver on the badly-scaled flattened triharmonic system was
>   observed (via `gdb`, confirming it was genuinely still computing, not
>   hung) to run for multiple *hours* without reaching either target. A
>   solver still short of the target at the 10-minute mark is stopped and
>   reports its best-effort backward error at that point — marked with a `†`
>   in the leaderboard tables — rather than either running unboundedly or
>   being cut off with no time guarantee at all.
> - **Eigen SparseLU**¹ is the slowest general-purpose solver by a wide
>   margin, as expected — but see the ⚠️ note below if you see it apparently
>   taking *minutes* instead of seconds, that's a build misconfiguration, not
>   real solver cost.
> - **NVIDIA Warp**'s `warp.optim.linear` solvers (`cg`/`cr`/`bicgstab`/
>   `gmres`, timed separately via [`warp_bench/`](warp_bench/)) use the
>   identical backward-error target and safety nets as Eigen's iterative
>   solvers (verified by reading both libraries' source, so a comparison
>   between them is apples-to-apples), but with only Jacobi preconditioning
>   (the strongest Warp currently offers — Eigen uses incomplete-LU,
>   meaningfully stronger). `cg` (the best-conditioned method here) reaches
>   near-machine-precision backward error on most systems even at this
>   mesh's real scale (360K-2.16M rows); `cr`/`bicgstab`/`gmres` are much
>   less consistent, especially on the indefinite mixed systems, and often
>   don't reach the target within budget — a genuine result about these
>   specific solver/preconditioner combinations, not an artifact of the
>   metric.
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
> `IGL_WITH_MKL`, not `EIGEN_USE_MKL_ALL`. A related, since-reverted fix is
> worth knowing about too: `SuiteSparse`'s own `find_package(BLAS)` is
> independent of this project's own MKL linkage, and could in principle pick
> a different MKL threading layer, which was once diagnosed (via `gdb`) as a
> ~100x thread-pool-contention slowdown and "fixed" by forcing SuiteSparse's
> BLAS to sequential (`BLA_VENDOR Intel10_64lp_seq`) — but that also
> serializes CHOLMOD/UMFPACK's own internal BLAS3 calls, and later A/B timing
> (full dragon mesh, k=1..3) showed plain threaded MKL BLAS is faster across
> the board with no reappearance of the contention, so `CMakeLists.txt` now
> leaves SuiteSparse's BLAS threaded (`BLA_VENDOR Intel10_64lp`). If that
> contention bug resurfaces, this is the first place to look.

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
> answer** — correctly caught by this benchmark's did-not-actually-succeed
> check rather than shown as a real result. Root cause, confirmed directly
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
> the **mixed biharmonic** system (essentially machine precision afterward)
> but does **not** fix the mixed triharmonic one — empirically (`IR_N_STEPS`
> from 2 up to 20 tested on a small synthetic case) the result is at best
> inconsistent, sometimes fixed and sometimes not run-to-run on the
> *identical* input (measured with the benchmark's old absolute-residual
> metric, before the switch to backward error below: ~1e-6 to ~5.9), matching
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
| `IGL_WITH_CHOLMOD`  | SuiteSparse CHOLMOD + UMFPACK (built from the `SuiteSparse/` submodule) |
| `IGL_WITH_MKL`      | Intel MKL Pardiso                                   |
| `IGL_WITH_GPL`      | GPL-licensed Eigen sparse solver code                |
| `IGL_WITH_CUDSS`    | NVIDIA cuDSS + cuSOLVER (requires the CUDA toolkit) |
| `IGL_WITH_NASOQ`    | NASOQ's LBL parallel symmetric-indefinite solver (requires `IGL_WITH_MKL` for its BLAS backend; built from the `nasoq/` submodule) |
| `IGL_WITH_CATAMARI` | catamari's SparseLDL solver (built from the `catamari`/`mantis`/`quotient` submodules) |
| `IGL_WITH_MA57`     | MA57 (symla), a from-scratch header-only symmetric-indefinite solver (built from the `ma57/` submodule; also requires system Eigen3 and SuiteSparse AMD, e.g. `apt install libeigen3-dev libsuitesparse-dev`) |

Every option above is probed at configure time (submodule checked out?
dependency findable?) and quietly turns itself `OFF` with a `message(WARNING
...)` if not -- a fresh `cmake ..` with no flags always succeeds and just
benchmarks whatever's actually available on the machine, rather than hard
`FATAL_ERROR`ing on a missing optional dependency.

## Run

    ./sparse_solver_benchmark [path to triangle mesh]
    ./sparse_solver_benchmark --grid 20 --check   # fast synthetic-mesh correctness check
    ./sparse_solver_benchmark --dragon 20000      # fast, real-mesh dev-loop (see below)
    ./sparse_solver_benchmark mesh.ply --csv results.csv   # also dump raw results
    ./sparse_solver_benchmark mesh.ply --only nasoq        # only run solvers matching "nasoq"
    ./sparse_solver_benchmark mesh.ply --exclude umfpack,sparselu   # skip the slow ones

`--dragon N` decimates (via `igl::decimate`) the checked-in
`xyzrgb_dragon-720K.ply` down to ~N faces and caches the result under
`meshes/` (gitignored) next to wherever the full-resolution mesh is found, so
repeated runs at the same N are instant after the first. Unlike `--grid`'s
synthetic triangulated grid -- trivially reorderable/banded, so it doesn't
exercise fill-in the way a real mesh's cotangent Laplacian does -- this is a
genuine (if smaller) mesh, useful for iterating on solver/tolerance logic
without waiting on the full 720K-vertex mesh. Use the full mesh only for
final leaderboard numbers (see [`leaderboards/`](leaderboards/) below).

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
not the 720K-vertex dragon) and fails if any compiled-in solver's backward
error exceeds `kBackwardErrorDivergedThreshold` (`1e-6`, the same fixed bar
used for the leaderboard's own did-not-actually-succeed check, applied
uniformly across all 5 systems since backward error is scale-invariant) —
a correctness regression test, not a performance one. Covers all 5 systems
(3 flattened + 2 mixed/indefinite).
GitHub Actions runs this on Linux (full solver matrix), Windows, and macOS
(Eigen + CHOLMOD/UMFPACK) on every push/PR; see `.github/workflows/ci.yml`.
Hosted runners have no GPU, so cuDSS/cuSOLVER still compile there but detect
the missing CUDA device at runtime and report `skipped` rather than fail.

## Leaderboards

Full per-machine results tables live under [`leaderboards/`](leaderboards/),
one markdown file per machine, named after its hardware (e.g.
`l40-2x-xeon-platinum-8362.md`) with the exact CPU/GPU/RAM/OS/mesh/date at
the top of the file. Results aren't comparable across machines (different
CPU, GPU, core count, memory bandwidth), so rather than one shared leaderboard
this project keeps a report card per machine it's been run on -- add your
own by running the benchmark and dropping a new file in that folder, named
after your hardware. Put any `[!NOTE]` callouts (explaining an anomaly,
a new solver, a fix applied since the last run, etc.) *after* the result
tables, not before -- the tables are what someone's scanning for first.

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
systems?"). Otherwise (like `catamari::SparseLDL`/`solve_catamari_ldl`,
MA57/`solve_symla`, or the cuDSS/cuSOLVER wrappers), add a dedicated function
that calls `record(...)` with the timing/residual, following the existing
examples.

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

### How is accuracy measured? (componentwise relative backward error)

The `Backward error` column is the **componentwise relative backward
error** — LAPACK's BERR — computed in `backward_error()` in `main.cpp`:

    η_cw = max_ij |B - A·X|_ij / (|A|·|X| + |B|)_ij

where `|A|` means entrywise absolute value (not an induced norm) and the
max is over every row *and* every right-hand-side column at once (`B = M*x`
has 3 columns, the mesh's own x/y/z coordinates). Conceptually: the
smallest relative, entrywise perturbation of `A` and `b` that would make
the computed `X` an *exact* solution — how much would you have to distrust
the problem's own entries, not the answer, for this to be exactly right?
`|A|·|X|` is one sparse-times-dense product (`|A|` has `A`'s exact
sparsity pattern, no dense matrix ever formed), the same complexity class
as computing `A·X` itself.

This replaced a plain absolute residual (`‖B−AX‖∞`, still computable as
`(rhs-Q*U).cwiseAbs()` if you want it) as this benchmark's accuracy metric
because backward error is **scale-invariant**: it correctly reports a
near-machine-precision solve even when `A` has astronomically large entries
or wildly different row scales (as in the mixed systems' `M`/`L` blocks),
because the denominator is formed with the same arithmetic — and thus
subject to the same rounding — that produced the residual in the numerator.
Measured directly (via `--dump-matrices`) on the same k=3 (flattened
Triharmonic) system on both example meshes:

| Mesh | `Q`'s largest entry (abs) | Absolute residual | Absolute ÷ `‖rhs‖∞` | Backward error `η_cw` |
|---|---:|---:|---:|---:|
| synthetic 20×20 grid | 1.8×10⁷ | 1.8×10⁻⁸ | 6.9×10⁻⁶ | **~10⁻¹⁶** |
| `xyzrgb_dragon-720K.ply` | **3.0×10¹⁵** | 47.9 | 0.37 | **~10⁻¹⁶** |

An absolute residual of "47.9" on the dragon mesh looks like a badly wrong
answer next to `rhs`'s own scale (~131) — that's what this benchmark used
to report, and it's genuinely misleading: `Q` there has entries up to
3×10¹⁵ (`Wᵏ⁺¹ = Wᵏ M⁻¹ L` applied recursively keeps inverting the mass
matrix, whose entries scale with triangle area, so small triangles push
`Q`'s entries to astronomical magnitudes), so an absolute residual around
48 is actually *at the floating-point roundoff floor* for arithmetic of
that magnitude — essentially machine precision, once you measure it against
what actually produced it rather than against `rhs`'s unrelated scale.
Backward error reports this correctly as `~1e-16` on **both** meshes and
**every** k, direct or iterative solver alike, which is why this benchmark
no longer needs a per-k tolerance table (`g_check_tol[]` used to range from
`1e-4` to `1e3`, hand-tuned per system) — a single fixed bar
(`kBackwardErrorDivergedThreshold = 1e-6`) now works everywhere, for both
`--check`'s correctness gate and the leaderboard's did-not-actually-succeed
reclassification (see below).

Iterative solvers are driven toward `kBackwardErrorTarget = 1e-8`, checked
ourselves after every chunk of the time-boxed solve — not via Eigen's
`setTolerance()` or Warp's `tol=`, both of which are a *relative L2*
quantity (`‖b−Ax‖₂ < tol·‖b‖₂`), a different norm entirely, and (more
importantly) each library's own private notion of "converged" rather than
one shared, externally-verified target. Comparing solvers by timing them to
different accuracy levels wouldn't be a fair leaderboard.

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
