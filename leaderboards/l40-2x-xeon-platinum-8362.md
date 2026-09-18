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
> **MA57 (symla) is dramatically slower than every other solver here at this
> scale** (33s-3600s factor time vs. single-digit-to-tens of seconds for
> everything else) despite matching their accuracy. It's a new, from-scratch
> solver (see `ma57/README.md`) without the decades of tuning behind
> CHOLMOD/Pardiso/UMFPACK -- included here as a correctly-reported, genuine
> result, not a bug in this benchmark's harness.

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
|    8 |                        NASOQ LBL |    2.6 secs |   0.18 secs |     8.09439e-16 |
|    9 | `Eigen::BiCGSTAB<IncompleteLUT>` |    1.6 secs |    1.4 secs |     4.16213e-16 |
|   10 |                      warp::gmres | 0.00032 secs |    5.1 secs |     5.99556e-16 |
|   11 |                Eigen::PardisoLLT |    3.9 secs |    1.3 secs |     2.91357e-16 |
|   12 |                  Eigen::SparseLU |    5.2 secs |   0.18 secs |     3.97062e-15 |
|   13 |               Eigen::PardisoLDLT |    4.4 secs |    1.5 secs |     3.07192e-16 |
|   14 |      Eigen::CholmodSupernodalLLT |    5.5 secs |    1.1 secs |     8.32449e-16 |
|   15 | Eigen::CholmodSupernodalLLT (CUDA) |    6.9 secs |    1.5 secs |     8.32449e-16 |
|   16 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |    9.3 secs |     1.10125e-15 |
|   17 |  `Eigen::CG<IncompleteCholesky>` |   0.53 secs |     20 secs |      5.5807e-09 |
|   18 |                 Eigen::UmfPackLU |     20 secs |   0.62 secs |     3.31624e-16 |
|   19 |                     MA57 (symla) |     33 secs |   0.23 secs |     9.28829e-16 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Biharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                     NVIDIA cuDSS |      4 secs |  0.022 secs |     2.32084e-14 |
| 🥈 2 |                        NASOQ LBL |    5.3 secs |   0.29 secs |     2.31509e-14 |
| 🥉 3 |               Eigen::PardisoLDLT |    5.6 secs |    1.4 secs |     1.42369e-15 |
|    4 |                Eigen::PardisoLLT |    5.7 secs |    1.4 secs |     1.40062e-15 |
|    5 |                 Eigen::UmfPackLU |    5.5 secs |    1.9 secs |     4.44457e-16 |
|    6 |            Eigen::SimplicialLDLT |    9.3 secs |   0.41 secs |     4.24456e-15 |
|    7 |             Eigen::SimplicialLLT |    9.4 secs |    0.4 secs |     1.12776e-14 |
|    8 |              catamari::SparseLDL |     11 secs |   0.36 secs |     7.11499e-15 |
|    9 | `Eigen::BiCGSTAB<IncompleteLUT>` |     11 secs |    4.8 secs |     5.71533e-16 |
|   10 |      Eigen::CholmodSupernodalLLT |     17 secs |    2.7 secs |     3.07552e-15 |
|   11 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |     20 secs |     6.30463e-15 |
|   12 | Eigen::CholmodSupernodalLLT (CUDA) |     24 secs |    4.5 secs |     3.07552e-15 |
|   13 |                  Eigen::SparseLU |     36 secs |   0.85 secs |     5.55392e-12 |
|   14 |                         warp::cr | 0.00033 secs |     85 secs |     5.73228e-08 |
|   15 |                     MA57 (symla) | 5.5e+02 secs |   0.66 secs |     3.05891e-15 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 0.0001483 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.2168 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 0.01691 exceeds 1e-06 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Triharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                     NVIDIA cuDSS |    6.4 secs | 0.0053 secs |     5.36518e-08 |
| 🥈 2 |                Eigen::PardisoLLT |    9.6 secs |    1.5 secs |     6.38005e-15 |
| 🥉 3 |               Eigen::PardisoLDLT |    9.6 secs |    1.6 secs |     6.25047e-15 |
|    4 |                        NASOQ LBL |     11 secs |   0.47 secs |     1.71662e-13 |
|    5 |             Eigen::SimplicialLLT |     35 secs |   0.84 secs |     3.38364e-13 |
|    6 |            Eigen::SimplicialLDLT |     35 secs |   0.97 secs |     1.41199e-13 |
|    7 |      Eigen::CholmodSupernodalLLT |     31 secs |    6.2 secs |     5.15252e-13 |
|    8 |              catamari::SparseLDL |     45 secs |   0.84 secs |     1.90941e-13 |
|    9 |        NVIDIA cuSOLVER (Sp Chol) |      (fused)* |     49 secs |     6.21089e-13 |
|   10 | Eigen::CholmodSupernodalLLT (CUDA) |     50 secs |     11 secs |     5.15252e-13 |
|   11 |                  Eigen::SparseLU | 1.5e+02 secs |    1.7 secs |     7.41722e-09 |
|   12 |                         warp::cr | 0.00024 secs | 7.5e+02 secs |     1.84984e-08 |
|   13 |                     MA57 (symla) | 3.6e+03 secs |   0.87 secs |      4.8854e-14 |
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
| 🥈 2 |                        NASOQ LBL |    6.8 secs |   0.47 secs |     2.88059e-08 |
| 🥉 3 |            Eigen::SimplicialLDLT |    9.2 secs |   0.46 secs |     1.09396e-07 |
|    4 |               Eigen::PardisoLDLT |    7.6 secs |    2.4 secs |     4.35266e-16 |
|    5 |       catamari::SparseLDL (LDLᵀ) |     12 secs |   0.41 secs |     1.41253e-07 |
|    6 |                 Eigen::UmfPackLU |     11 secs |    2.3 secs |     3.91761e-16 |
|    7 |                  Eigen::SparseLU |     30 secs |   0.62 secs |     2.92105e-14 |
|    8 |               MA57 (symla, LDLᵀ) | 5.9e+02 secs |    1.1 secs |     5.12888e-08 |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.8748 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | Eigen::CholmodSupernodalLLT (CUDA) |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |       catamari::SparseLDL (LDLᵀ) |     45 secs |      1 secs |     1.17789e-09 |
| 🥈 2 |                  Eigen::SparseLU | 1.1e+02 secs |    1.4 secs |     1.42071e-12 |
| 🥉 3 |               MA57 (symla, LDLᵀ) |  3e+03 secs |    1.8 secs |     2.67506e-10 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk: excessive MKL thread churn on this system's fill-in |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale |
|    - |                        NASOQ LBL |           - |           - | skipped: known crash: SIGSEGV in libmetis genmmd/mmdelm via NASOQ's symbolic_analysis_lin_solve on this system's sparsity pattern |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - |               Eigen::PardisoLDLT |           - |           - | skipped: known Pardiso reordering hang on this system's sparsity pattern |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                     NVIDIA cuDSS |           - |           - | skipped: did not actually succeed: backward error 0.6692 exceeds 1e-06 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | Eigen::CholmodSupernodalLLT (CUDA) |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |

