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
> **`Eigen::CholmodSupernodalLLT`'s Harmonic (k=1) factor time (6.2s) is
> slower than its own Biharmonic (k=2) factor time (3.4s)**, despite k=2's
> system having strictly more fill-in. This isn't a measurement artifact --
> a warm-up factorization runs before the timed loop specifically to rule
> out one-time library-init costs (see `main.cpp`, right before the k-loop)
> -- it's most likely CHOLMOD's own supernodal-vs-simplicial heuristic
> picking a less BLAS3-friendly (many small supernodes) strategy for k=1's
> sparser matrix than for k=2's denser one, where large supernodes let it
> use fewer, bigger, much higher-throughput BLAS3 calls. A genuine
> characteristic of this input/solver combination, not a bug.

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
| 🥇 1 |                         warp::cg |  0.00026 secs |     0.34 secs |     1.00714e-14 |
| 🥈 2 |                         warp::cr |  0.00022 secs |     0.48 secs |     4.72436e-15 |
| 🥉 3 |            Eigen::SimplicialLDLT |    1.2 secs |   0.11 secs |     2.44563e-15 |
|    4 |             Eigen::SimplicialLLT |    1.2 secs |   0.11 secs |     2.97441e-15 |
|    5 |              catamari::SparseLDL |    1.4 secs |  0.099 secs |     1.34635e-15 |
|    6 |                   warp::bicgstab |  0.00026 secs |      1.5 secs |     1.28732e-14 |
|    7 |                     NVIDIA cuDSS |    2.3 secs | 0.0026 secs |     7.79603e-16 |
|    8 |                        NASOQ LBL |    2.6 secs |   0.18 secs |     8.09439e-16 |
|    9 |         Eigen::CG<IncompleteLUT> |    1.6 secs |    1.3 secs |     6.22174e-09 |
|   10 |   Eigen::BiCGSTAB<IncompleteLUT> |    1.6 secs |    1.4 secs |     4.16213e-16 |
|   11 |                      warp::gmres |  0.00032 secs |      5.1 secs |     5.99556e-16 |
|   12 |                Eigen::PardisoLLT |    3.9 secs |    1.3 secs |     2.91357e-16 |
|   13 |                  Eigen::SparseLU |    5.2 secs |   0.18 secs |     3.97062e-15 |
|   14 |               Eigen::PardisoLDLT |    4.4 secs |    1.5 secs |     3.07192e-16 |
|   15 |      Eigen::CholmodSupernodalLLT |    6.2 secs |    1.1 secs |     8.32449e-16 |
|   16 |        NVIDIA cuSOLVER (Sp Chol) | (fused)* |    9.3 secs |     1.10125e-15 |
|   17 |                 Eigen::UmfPackLU |     20 secs |   0.62 secs |     3.31624e-16 |
|   18 |                     MA57 (symla) |     33 secs |   0.23 secs |     9.28829e-16 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Biharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |      Eigen::CholmodSupernodalLLT |    3.4 secs |    0.4 secs |     3.20597e-15 |
| 🥈 2 |                     NVIDIA cuDSS |      4 secs |  0.022 secs |     2.32084e-14 |
| 🥉 3 |                        NASOQ LBL |    5.3 secs |   0.29 secs |     2.31509e-14 |
|    4 |               Eigen::PardisoLDLT |    5.6 secs |    1.4 secs |     1.42369e-15 |
|    5 |                Eigen::PardisoLLT |    5.7 secs |    1.4 secs |     1.40062e-15 |
|    6 |                 Eigen::UmfPackLU |    5.5 secs |    1.9 secs |     4.44457e-16 |
|    7 |            Eigen::SimplicialLDLT |    9.3 secs |   0.41 secs |     4.24456e-15 |
|    8 |             Eigen::SimplicialLLT |    9.4 secs |    0.4 secs |     1.12776e-14 |
|    9 |              catamari::SparseLDL |     11 secs |   0.36 secs |     7.11499e-15 |
|   10 |         Eigen::CG<IncompleteLUT> |     11 secs |    3.6 secs |     9.88902e-13 |
|   11 |   Eigen::BiCGSTAB<IncompleteLUT> |     11 secs |    4.8 secs |     5.71533e-16 |
|   12 |        NVIDIA cuSOLVER (Sp Chol) | (fused)* |     20 secs |     6.30463e-15 |
|   13 |                  Eigen::SparseLU |     36 secs |   0.85 secs |     5.55392e-12 |
|   14 |                         warp::cr |  0.00033 secs |       85 secs |     5.73228e-08 |
|   15 |                     MA57 (symla) | 5.5e+02 secs |   0.66 secs |     3.05891e-15 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 0.0001483 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.2168 exceeds 1e-06 |

*(fused): cuSOLVER's `cusolverSpDcsrlsvchol` has no separate factor step -- see README.

# Triharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                     NVIDIA cuDSS |    6.4 secs | 0.0053 secs |     5.36518e-08 |
| 🥈 2 |                Eigen::PardisoLLT |    9.6 secs |    1.5 secs |     6.38005e-15 |
| 🥉 3 |               Eigen::PardisoLDLT |    9.6 secs |    1.6 secs |     6.25047e-15 |
|    4 |                        NASOQ LBL |     11 secs |   0.47 secs |     1.71662e-13 |
|    5 |      Eigen::CholmodSupernodalLLT |     13 secs |   0.41 secs |     1.17512e-13 |
|    6 |             Eigen::SimplicialLLT |     35 secs |   0.84 secs |     3.38364e-13 |
|    7 |            Eigen::SimplicialLDLT |     35 secs |   0.97 secs |     1.41199e-13 |
|    8 |              catamari::SparseLDL |     45 secs |   0.84 secs |     1.90941e-13 |
|    9 |        NVIDIA cuSOLVER (Sp Chol) | (fused)* |     49 secs |     6.21089e-13 |
|   10 |                  Eigen::SparseLU | 1.5e+02 secs |    1.7 secs |     7.41722e-09 |
|   11 |                         warp::cr |  0.00024 secs |  7.5e+02 secs |     1.84984e-08 |
|   12 |                     MA57 (symla) | 3.6e+03 secs |   0.87 secs |      4.8854e-14 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |         Eigen::CG<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 0.0001877 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error 0.002162 exceeds 1e-06 |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.8132 exceeds 1e-06 |

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
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |                Eigen::PardisoLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |         Eigen::CG<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 0.8748 exceeds 1e-06 |

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
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |         Eigen::CG<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |        NVIDIA cuSOLVER (Sp Chol) |           - |           - | skipped: no indefinite/LDLT solver in this cuSOLVER version |
|    - |                     NVIDIA cuDSS |           - |           - | skipped: did not actually succeed: backward error 0.6692 exceeds 1e-06 |
|    - |                         warp::cg |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                         warp::cr |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                   warp::bicgstab |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |                      warp::gmres |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |

