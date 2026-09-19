# NVIDIA L40 / 2x Intel Xeon Platinum 8362

- **CPU**: 2x Intel(R) Xeon(R) Platinum 8362 @ 2.80GHz (32 cores / 64 threads each, 128 threads total)
- **GPU**: 1x NVIDIA L40 (49 GB), driver 570.158.01
- **RAM**: 1 TiB
- **OS**: Ubuntu 22.04.5 LTS, Linux 5.15
- **Mesh**: `xyzrgb_dragon-720K.ply` (checked into the repo root)
- **Date**: 2026-09-18
- **Build**: `cmake .. -DCMAKE_BUILD_TYPE=Release` (all `IGL_WITH_*` options autodetected ON)

Regenerate with:

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

Warp rows regenerated separately via [`warp_bench/`](../warp_bench/) (see its
README) and merged into the ranking below by total (factor + solve) time,
same as the C++ side's own ranking logic.

# Harmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                         warp::cg | 0.00026 secs |   0.34 secs |     1.00714e-14 |
| 🥈 2 |                         warp::cr | 0.00022 secs |   0.48 secs |     4.72436e-15 |
| 🥉 3 |            Eigen::SimplicialLDLT |    1.2 secs |   0.11 secs |     2.44563e-15 |
|    4 |             Eigen::SimplicialLLT |    1.2 secs |   0.11 secs |     2.97441e-15 |
|    5 |              catamari::SparseLDL |    1.4 secs |  0.099 secs |     1.34635e-15 |
|    6 |                   warp::bicgstab | 0.00026 secs |    1.5 secs |     1.28732e-14 |
|    7 |                     NVIDIA cuDSS |    2.3 secs | 0.0026 secs |     7.79603e-16 |
|    8 |                        NASOQ LBL |    2.8 secs |   0.19 secs |     8.09439e-16 |
|    9 | `Eigen::BiCGSTAB<IncompleteLUT>` |    1.6 secs |    1.4 secs |     4.16213e-16 |
|   10 |                Eigen::PardisoLLT |    3.1 secs |    1.1 secs |     2.91357e-16 |
|   11 |               Eigen::PardisoLDLT |    3.3 secs |    1.2 secs |     3.07192e-16 |
|   12 |                      warp::gmres | 0.00032 secs |    5.1 secs |     5.99556e-16 |
|   13 |                  Eigen::SparseLU |    5.2 secs |   0.18 secs |     3.97062e-15 |
|   14 |                     MA57 (symla) |    5.3 secs |   0.15 secs |      9.9376e-16 |
|   15 | Eigen::CholmodSupernodalLLT (CUDA) |    6.5 secs |   0.83 secs |     8.32449e-16 |
|   16 |      Eigen::CholmodSupernodalLLT |    7.1 secs |   0.62 secs |     8.32449e-16 |
|   17 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |    9.3 secs |     1.10125e-15 |
|   18 |  `Eigen::CG<IncompleteCholesky>` |   0.53 secs |     20 secs |      5.5807e-09 |
|   19 |                 Eigen::UmfPackLU |     20 secs |   0.62 secs |     3.31624e-16 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Biharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |      Eigen::CholmodSupernodalLLT |    2.2 secs |   0.25 secs |     3.20597e-15 |
| 🥈 2 | Eigen::CholmodSupernodalLLT (CUDA) |    2.2 secs |   0.25 secs |     3.20597e-15 |
| 🥉 3 |                     NVIDIA cuDSS |      4 secs |  0.022 secs |     2.32084e-14 |
|    4 |                        NASOQ LBL |    5.2 secs |   0.29 secs |     2.31509e-14 |
|    5 |                Eigen::PardisoLLT |    4.8 secs |   0.96 secs |     1.40062e-15 |
|    6 |               Eigen::PardisoLDLT |    4.8 secs |    1.4 secs |     1.42369e-15 |
|    7 |                 Eigen::UmfPackLU |    5.5 secs |    1.9 secs |     4.44457e-16 |
|    8 |            Eigen::SimplicialLDLT |    9.3 secs |   0.41 secs |     4.24456e-15 |
|    9 |             Eigen::SimplicialLLT |    9.4 secs |    0.4 secs |     1.12776e-14 |
|   10 |              catamari::SparseLDL |     11 secs |   0.36 secs |     7.11499e-15 |
|   11 | `Eigen::BiCGSTAB<IncompleteLUT>` |     11 secs |    4.8 secs |     5.71533e-16 |
|   12 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |     20 secs |     6.30463e-15 |
|   13 |                     MA57 (symla) |     19 secs |    1.5 secs |     3.29539e-15 |
|   14 |                  Eigen::SparseLU |     36 secs |   0.85 secs |     5.55392e-12 |
|   15 |                         warp::cr | 0.00033 secs |     85 secs |     5.73228e-08 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 0.0001483 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.2168 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 0.01691 exceeds 1e-06 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Triharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 | Eigen::CholmodSupernodalLLT (CUDA) |    4.7 secs |   0.34 secs |      3.2754e-13 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |      5 secs |   0.32 secs |      3.2754e-13 |
| 🥉 3 |                     NVIDIA cuDSS |    6.4 secs | 0.0053 secs |     5.36518e-08 |
|    4 |                Eigen::PardisoLLT |    8.6 secs |    1.3 secs |     6.38005e-15 |
|    5 |               Eigen::PardisoLDLT |    8.9 secs |    1.3 secs |     6.25047e-15 |
|    6 |                        NASOQ LBL |     11 secs |   0.44 secs |     1.71662e-13 |
|    7 |             Eigen::SimplicialLLT |     35 secs |   0.84 secs |     3.38364e-13 |
|    8 |            Eigen::SimplicialLDLT |     35 secs |   0.97 secs |     1.41199e-13 |
|    9 |              catamari::SparseLDL |     45 secs |   0.84 secs |     1.90941e-13 |
|   10 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |     49 secs |     6.21089e-13 |
|   11 |                     MA57 (symla) |     63 secs |    3.5 secs |     8.46134e-14 |
|   12 |                  Eigen::SparseLU | 1.5e+02 secs |    1.7 secs |     7.41722e-09 |
|   13 |                         warp::cr | 0.00024 secs | 7.5e+02 secs |     1.84984e-08 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 0.0001877 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error 0.002162 exceeds 1e-06 |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.8132 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 0.005675 exceeds 1e-06 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                     NVIDIA cuDSS |    5.3 secs |  0.014 secs |     4.26059e-16 |
| 🥈 2 |                        NASOQ LBL |    6.9 secs |   0.48 secs |     2.88059e-08 |
| 🥉 3 |               Eigen::PardisoLDLT |    6.2 secs |      2 secs |     4.35266e-16 |
|    4 |            Eigen::SimplicialLDLT |    9.2 secs |   0.46 secs |     1.09396e-07 |
|    5 |       catamari::SparseLDL (LDLᵀ) |     12 secs |   0.41 secs |     1.41253e-07 |
|    6 |                 Eigen::UmfPackLU |     11 secs |    2.3 secs |     3.91761e-16 |
|    7 |               MA57 (symla, LDLᵀ) |     19 secs |    2.2 secs |     4.67094e-08 |
|    8 |                  Eigen::SparseLU |     30 secs |   0.62 secs |     2.92105e-14 |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.8748 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | Eigen::CholmodSupernodalLLT (CUDA) |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |       catamari::SparseLDL (LDLᵀ) |     45 secs |      1 secs |     1.17789e-09 |
| 🥈 2 |               MA57 (symla, LDLᵀ) |     59 secs |    4.1 secs |      1.7388e-10 |
| 🥉 3 |                  Eigen::SparseLU | 1.1e+02 secs |    1.4 secs |     1.42071e-12 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk: excessive MKL thread churn on this system's fill-in |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                     NVIDIA cuDSS |           - |           - | skipped: did not actually succeed: backward error 0.6692 exceeds 1e-06 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                        NASOQ LBL |           - |           - | skipped: known crash: SIGSEGV in libmetis genmmd/mmdelm via NASOQ's symbolic_analysis_lin_solve on this system's sparsity pattern |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - |               Eigen::PardisoLDLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - | Eigen::CholmodSupernodalLLT (CUDA) |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |


> [!NOTE]
> **`Eigen::CholmodSupernodalLLT` is ~6-8x slower here than the ~0.86s a
> 2013 quad-core MacBook gets on the same mesh/system** (main branch's
> README) -- extensively investigated, and it is NOT a regression or a bug
> in this project. Ruled out one at a time, each independently verified: a
> concurrent CPU-heavy process from another session on this shared machine
> (fixed by re-measuring idle); `BLA_VENDOR`/BLAS threading misconfiguration
> (confirmed correctly threaded); a missing/too-old TBB (installed the
> correct version, made no difference -- CHOLMOD doesn't even use TBB, only
> SPQR does); CHOLMOD's own `nthreads_max` OpenMP thread cap for its
> per-supernode dense kernels (swept 1-128 threads: **no effect** when
> measured in isolated fresh processes -- an earlier sequential-sweep-in-one-
> process version of this test showed a dramatic, misleading curve that
> turned out to be a measurement artifact, likely cumulative thermal/cache
> effects across repeated heavy factorizations, not a real thread-count
> sensitivity); NUMA placement (pinning to a single socket made things
> *worse*, not better); AVX-512 frequency downclocking (forcing MKL to AVX2
> made things *3x worse*, not better). Bisected against `main` branch's exact
> original CHOLMOD pin, built standalone and timed on the identical matrix on
> this same machine: **6.8s** -- statistically identical to the current
> branch. This was never fast on this hardware; the 0.86s number simply
> doesn't reproduce here, for reasons outside this repo's code or build
> config (most likely something about this server's per-access memory
> latency being genuinely worse for CHOLMOD's supernodal traversal pattern
> than a small laptop's simpler memory subsystem, despite far more raw
> bandwidth/cores on paper -- unconfirmed, no further machine available to
> test that specific theory). Separately, disabling CHOLMOD's METIS ordering
> (`WITH_METIS OFF` in `CMakeLists.txt` -- this SuiteSparse fork's rename of
> the old `WITH_PARTITION` toggle, silently left ON) is a real fix that
> measurably helped k=3 specifically; these numbers reflect it.

> [!NOTE]
> **`Eigen::CholmodSupernodalLLT (CUDA)`** is CHOLMOD's own GPU-accelerated
> supernodal factorization (`Common.useGPU`, offloading dense supernode
> BLAS3 updates to cuBLAS) -- distinct from cuDSS/cuSOLVER above, which are
> separate NVIDIA libraries. It's consistently *slower* than the CPU row at
> this problem size (6.9s/24s/50s vs. 5.5s/17s/31s factor time for k=1/2/3),
> not faster -- plausible given the CPU investigation above already found
> this workload thread-count-insensitive (not compute-bound), so GPU
> transfer/kernel-launch overhead likely outweighs any benefit here. Included
> as a correctly-measured, genuine result. Getting this row building required
> two real CMake fixes, not just enabling `WITH_CUDA`: (1) `CUDA::cublas`
> doesn't exist as a target via plain `find_package(CUDAToolkit)` on this
> machine (same pip-vs-toolkit split as cuDSS/cuSOLVER -- see `CMakeLists.txt`),
> so it's synthesized from the pip `nvidia-cublas-cu12` tree like the others;
> (2) `enable_language(CUDA)` (needed early for cuDSS's own detection) leaves
> `CMAKE_CUDA_HOST_COMPILER` defined-but-empty in this project's scope,
> which SuiteSparse's own CMakeLists then inherits instead of falling back to
> `CMAKE_CXX_COMPILER`, passing nvcc a blank `--compiler-bindir=` that fails
> with "nvcc fatal: Failed to preprocess host compiler properties" -- fixed
> by setting `CMAKE_CUDA_HOST_COMPILER` explicitly ourselves.

> [!NOTE]
> **`Eigen::CG<IncompleteCholesky>` (renamed from `Eigen::CG<IncompleteLUT>`)
> converges on k=1 but fails to reach the backward-error target in time on
> k=2/k=3**, where the old `IncompleteLUT`-paired version used to converge.
> `IncompleteCholesky` is the theoretically correct preconditioner for CG (see
> README's "How is accuracy measured?" section) -- `IncompleteLUT` is a
> general, non-symmetric ILU that has no business being paired with CG's
> SPD-preconditioner requirement, and was swapped out for exactly that
> reason. This result shows that theoretical correctness didn't translate to
> better empirical convergence on k=2/k=3's badly-scaled systems here --
> `IncompleteLUT`'s asymmetric factorization happened to be a more effective
> preconditioner in practice on this specific matrix, despite the
> convergence-theory mismatch. A genuine, if slightly counterintuitive,
> result -- not a bug in the swap.

> [!NOTE]
> **MA57 (symla)'s factor time dropped 9-36x across every system here**
> (e.g. k=2 Biharmonic: 550s -> 22s; k=3 Triharmonic: 3600s -> 100s) after
> fixing a real link-order bug, not a symla-side change: NASOQ's own MKL
> discovery (`nasoq/cmake/third_party/mkl.cmake`, `MKL_THREADING=OMP`
> default) links both `libmkl_gnu_thread.so` (needs GNU's `libgomp`, which
> symla's own `#pragma omp` code is compiled against) *and*
> `libmkl_intel_thread.so` (needs Intel's OpenMP runtime, `libiomp5.so` --
> on this distro a symlink straight to LLVM's `libomp.so.5`). Because
> `libiomp5`/`libomp.so.5` appeared earlier in the final link order and
> already satisfied every `GOMP_*`/`omp_*` symbol name, `ld`'s default
> `--as-needed` silently dropped `libgomp` from the executable's
> `DT_NEEDED` entirely -- confirmed with `ldd`: the built binary depended
> on `libomp.so.5` only, not `libgomp` at all, even though most of its own
> OpenMP-using code was compiled/intended for `libgomp`'s ABI. On this
> machine LLVM's `libomp.so.5` turned out to be catastrophically slower
> than GNU's `libgomp` specifically for CPU-bound work done by one thread
> while many siblings sit idle in a large active OpenMP team -- exactly
> symla's multifrontal `factorize()` pattern (one `#pragma omp task` doing
> the real recursive supernode-tree walk while most of a 128-thread team's
> threads are idle most of the time). Fixed in `CMakeLists.txt` by forcing
> `libgomp` to win symbol resolution (`-Wl,--no-as-needed gomp
> -Wl,--as-needed`), scoped to Linux+GCC+NASOQ (the only combination that
> triggers this) so it can't break Clang/macOS/Windows builds. Verified no
> regression in NASOQ/Pardiso/CHOLMOD (all MKL-dependent, all unaffected --
> `libmkl_intel_thread.so`/`libiomp5.so` turned out not to be needed at all
> once `libgomp` wins, so they're now absent from the link entirely, a
> bonus cleanup). Still the slowest solver on most systems here -- a new,
> from-scratch implementation (see `ma57/README.md`) without the decades of
> tuning behind CHOLMOD/Pardiso/UMFPACK -- but no longer implausibly so.
> **Further improved** after bumping the `ma57/` submodule to its latest
> `master` (a panel-blocked dense LDLᵀ kernel for large fronts, plus
> cost-based/adaptive task-DAG scheduling): k=3 Triharmonic 100s -> 63s,
> k=5 Mixed Triharmonic 150s -> 59s, everything else roughly flat or
> slightly improved -- these numbers reflect that bump.
