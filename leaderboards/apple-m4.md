# Apple M4 (MacBook, macOS)

- **CPU**: Apple M4 (10 cores)
- **GPU**: none (Apple Silicon has no CUDA support -- NVIDIA cuDSS/cuSOLVER
  rows are absent from this report)
- **RAM**: 32 GB
- **OS**: macOS 26.6.2 (Darwin 25.6.0), arm64
- **Mesh**: `xyzrgb_dragon-720K.ply` (checked into the repo root)
- **Date**: 2026-09-18
- **Build**: `cmake .. -DCMAKE_BUILD_TYPE=Release` (all `IGL_WITH_*` options
  autodetected; `IGL_WITH_MKL` and everything gated behind it --
  `Eigen::PardisoLLT`/`PardisoLDLT`, catamari's BLAS kernels -- auto-disable
  since there's no Intel MKL on Apple Silicon; NASOQ LBL now builds via
  Homebrew OpenBLAS instead, see NOTE at the bottom)

Regenerate with:

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

# Harmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |        Accelerate SparseCholesky |     0.27 secs |    0.038 secs |     1.18122e-15 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |      0.3 secs |    0.032 secs |     1.52008e-15 |
| 🥉 3 |                        NASOQ LBL |     0.35 secs |    0.061 secs |     7.96362e-16 |
|    4 |            Eigen::SimplicialLDLT |      0.6 secs |    0.054 secs |     2.33654e-15 |
|    5 |             Eigen::SimplicialLLT |     0.68 secs |    0.054 secs |     3.01865e-15 |
|    6 |              catamari::SparseLDL |     0.78 secs |    0.042 secs |     1.35589e-15 |
|    7 |                 Eigen::UmfPackLU |     0.71 secs |      0.2 secs |     2.71385e-16 |
|    8 | `Eigen::BiCGSTAB<IncompleteLUT>` |     0.85 secs |     0.59 secs |     3.91637e-16 |
|    9 |                  Eigen::SparseLU |      2.3 secs |    0.084 secs |     2.38898e-15 |
|   10 |                     MA57 (symla) |      2.8 secs |    0.064 secs |      1.0304e-15 |
|   11 |  `Eigen::CG<IncompleteCholesky>` |     0.27 secs |      9.9 secs |      4.0063e-09 |

# Biharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |        Accelerate SparseCholesky |     0.59 secs |    0.082 secs |     1.34589e-14 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |     0.63 secs |    0.074 secs |     4.91755e-15 |
| 🥉 3 |                        NASOQ LBL |     0.84 secs |     0.12 secs |     3.75248e-15 |
|    4 |                 Eigen::UmfPackLU |      2.3 secs |     0.74 secs |     5.27979e-16 |
|    5 |            Eigen::SimplicialLDLT |      4.7 secs |      0.2 secs |      7.8414e-15 |
|    6 |             Eigen::SimplicialLLT |      4.7 secs |      0.2 secs |     4.61287e-15 |
|    7 |              catamari::SparseLDL |      6.6 secs |     0.17 secs |     3.66465e-15 |
|    8 | `Eigen::BiCGSTAB<IncompleteLUT>` |        6 secs |      1.8 secs |     5.55785e-16 |
|    9 |                  Eigen::SparseLU |       13 secs |     0.28 secs |     4.51828e-12 |
|   10 |                     MA57 (symla) |       20 secs |     0.12 secs |     4.95481e-15 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 0.006303 exceeds 1e-06 |

# Triharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |      Eigen::CholmodSupernodalLLT |      1.3 secs |     0.14 secs |     2.38957e-13 |
| 🥈 2 |        Accelerate SparseCholesky |      1.6 secs |     0.22 secs |     5.66792e-13 |
| 🥉 3 |                        NASOQ LBL |      1.7 secs |     0.21 secs |     5.98064e-14 |
|    4 |             Eigen::SimplicialLLT |       26 secs |     0.68 secs |     5.90259e-13 |
|    5 |            Eigen::SimplicialLDLT |       26 secs |     0.74 secs |     3.99303e-13 |
|    6 |              catamari::SparseLDL |       35 secs |     0.41 secs |     9.05897e-13 |
|    7 |                  Eigen::SparseLU |       64 secs |      1.4 secs |     5.65782e-09 |
|    8 |                     MA57 (symla) |       95 secs |     0.26 secs |     7.25797e-14 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 0.006508 exceeds 1e-06 |

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                        NASOQ LBL |      1.2 secs |     0.19 secs |     8.34662e-09 |
| 🥈 2 |            Eigen::SimplicialLDLT |      5.7 secs |     0.28 secs |      1.2478e-07 |
| 🥉 3 |                 Eigen::UmfPackLU |      5.1 secs |      1.2 secs |     3.87298e-16 |
|    4 |     catamari::SparseLDL (LDLᵀ) |      7.1 secs |     0.19 secs |     4.69765e-08 |
|    5 |                  Eigen::SparseLU |       13 secs |     0.31 secs |     3.26471e-14 |
|    6 |             MA57 (symla, LDLᵀ) |       22 secs |     0.21 secs |     2.25769e-08 |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |        Accelerate SparseCholesky |           - |           - | skipped: factorization failed: status -1 |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 0.2572 exceeds 1e-06 |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |     catamari::SparseLDL (LDLᵀ) |       26 secs |     0.52 secs |     1.55308e-09 |
| 🥈 2 |                  Eigen::SparseLU |       50 secs |      1.4 secs |     3.47129e-12 |
| 🥉 3 |             MA57 (symla, LDLᵀ) |    1e+02 secs |     0.46 secs |     2.65418e-10 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk: excessive MKL thread churn on this system's fill-in |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale |
|    - |                        NASOQ LBL |           - |           - | skipped: known crash: SIGSEGV in libmetis genmmd/mmdelm via NASOQ's symbolic_analysis_lin_solve on this system's sparsity pattern |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |        Accelerate SparseCholesky |           - |           - | skipped: factorization failed: status -1 |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |  `Eigen::CG<IncompleteCholesky>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - | `Eigen::BiCGSTAB<IncompleteLUT>` |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |

> [!NOTE]
> **`Accelerate SparseCholesky`** (Apple's Accelerate Sparse Solvers library,
> `IGL_WITH_ACCELERATE_SPARSE`, macOS-only) is new on this report card. It
> ships in the OS on every Mac -- no submodule (unlike CHOLMOD's SuiteSparse)
> or external install (unlike MKL) required -- and takes 🥇 first place on
> all three SPD systems here, edging out `Eigen::CholmodSupernodalLLT` by a
> small margin on factor time while matching its accuracy. Like every other
> SPD-only solver here, it correctly fails/skips on the indefinite mixed
> k=4,5 systems (Accelerate has no exposed LDLT/indefinite factorization in
> this API).

> [!NOTE]
> **This machine's `SuiteSparse` build previously had CHOLMOD's METIS-based
> nested-dissection ordering silently enabled** (a newer SuiteSparse fork
> renamed the old `WITH_PARTITION` toggle this project's `CMakeLists.txt` set
> to `WITH_METIS`, default `ON`, so the old toggle became a no-op). A/B
> timing showed no measurable effect on k=1,2 but a ~4x triharmonic (k=3)
> factor-time regression (4.4s vs 1.1s). Fixed by explicitly setting
> `WITH_METIS OFF`; the numbers below reflect the fix.

> [!NOTE]
> **`NASOQ LBL` is now enabled on this machine**, previously disabled/absent
> from this report because NASOQ's own `NASOQ_BLAS_BACKEND` only supported
> `MKL` or `OpenBLAS`, and there's no Intel MKL on Apple Silicon. It now
> builds against Homebrew's `openblas` (`NASOQ_BLAS_BACKEND=OpenBLAS`,
> `NASOQ_USE_CLAPACK=ON`, see `CMakeLists.txt`'s `IGL_WITH_NASOQ` block). Two
> real bugs had to be fixed to get correct (rather than all-NaN) results:
> (1) NASOQ's vendored CLAPACK (a netlib f2c translation used for
> `DSYTRF`/`DLAPMT` under `NASOQ_USE_CLAPACK=ON`) defines its Fortran
> `integer`/`logical` types as 8-byte `long int` via f2c.h, but calls out to
> OpenBLAS routines (`dsyrk_`, `dgemv_`, etc.) that use the standard 4-byte
> LP64 Fortran `INTEGER` -- a latent ABI mismatch on every 64-bit platform,
> fixed by patching the fetched `f2c.h` to use 4-byte ints (see
> `nasoq/cmake/third_party/clapack.cmake` and
> `nasoq/include/nasoq/clapacke/clapacke.h`); (2) the actual all-NaN cause on
> this machine was that `main.cpp`'s `solve_nasoq_lbl` requested
> `ldl_variant = 4` ("Parallel SBK"), whose implementation in
> `nasoq/src/QP/linear_solver_wrapper.cpp` is entirely wrapped in
> `#ifdef OPENMP` with no fallback; since CMake's `find_package(OpenMP)`
> fails for AppleClang on this machine, NASOQ never defines `OPENMP` or
> compiles the parallel kernels, so `numerical_factorization()` silently
> returned "success" having factorized nothing. Switched to the
> unconditionally-available serial variant (`ldl_variant = 2`), which is
> also the right choice given this file already forces `num_thread = 1`.
> With both fixes, NASOQ LBL lands respectably (top 3 on Harmonic/
> Biharmonic/Triharmonic, and 🥇 first on Mixed Biharmonic). It's still
> unavailable on Mixed Triharmonic, where it crashes with a SIGSEGV inside
> `libmetis` during symbolic analysis on this system's sparsity pattern --
> unrelated to the two fixes above (filed as
> [sympiler/nasoq#33](https://github.com/sympiler/nasoq/issues/33); still
> reproduces on this machine, though here as an infinite hang inside
> `libmetis__mmdelm` rather than a SIGSEGV -- same underlying heap
> corruption, different downstream symptom, see that issue's comments).

> [!NOTE]
> **`NASOQ LBL`'s numbers above use `AMD` ordering, not NASOQ's previous
> (and still default) `METIS` ordering.** NASOQ's own `symbolic_phase.cpp`
> already had a complete, working AMD-ordering code path sitting right next
> to the METIS one, but it was permanently unreachable dead code -- `METIS`
> was unconditionally `#define`'d in a core header with no way to disable
> it. Exposed as a real `NASOQ_ORDERING` CMake option in
> [sympiler/nasoq#36](https://github.com/sympiler/nasoq/pull/36). A/B timing
> on this machine showed AMD is consistently **~3-4x faster than METIS** at
> matching accuracy on every system here: Harmonic 0.35s vs. 1.4s,
> Biharmonic 0.84s vs. 2.7s, Triharmonic 1.7s vs. 4.7s, Mixed Biharmonic
> 1.2s vs. 3.7s (factor time) -- the same lesson as this report's earlier
> CHOLMOD/METIS note above. This project's `CMakeLists.txt` now sets
> `NASOQ_ORDERING AMD` explicitly for its own NASOQ integration.

> [!NOTE]
> **A separate attempt to switch NASOQ's BLAS backend from OpenBLAS to
> Apple's Accelerate framework was evaluated and not adopted.** The
> nondeterminism originally reported in the unmerged
> [sympiler/nasoq#27](https://github.com/sympiler/nasoq/pull/27) was
> root-caused (a `void`-vs-`int` LAPACK return-value mismatch in
> `nasoq::clapacke::LAPACKE_dlapmt`, garbage-reading a nonexistent return
> value from Accelerate's `dlapmt_`) and fixed in
> [sympiler/nasoq#35](https://github.com/sympiler/nasoq/pull/35), verified
> deterministic across 20 repeated test runs plus ASan/UBSan-clean stress
> tests. But at the scale this benchmark actually runs at, a second,
> distinct, unresolved problem showed up on the indefinite systems:
> backward error grew with problem size, and the full-mesh Mixed Biharmonic
> case hung for 15+ minutes instead of OpenBLAS's 3.7s. Likely a separate
> bug in Accelerate's Bunch-Kaufman pivoting or supernodal blocking at
> scale, not investigated further -- OpenBLAS remains this project's NASOQ
> backend on Apple Silicon.
> Apple Accelerate support for NASOQ's BLAS backend
> (https://github.com/sympiler/nasoq/pull/27) remains unmerged and was not
> pursued in this pass; OpenBLAS was sufficient to get correct results.
