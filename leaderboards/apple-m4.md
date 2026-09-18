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
  `Eigen::PardisoLLT`/`PardisoLDLT`, NASOQ LBL, catamari's BLAS kernels --
  auto-disable since there's no Intel MKL on Apple Silicon)

Regenerate with:

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

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

# Harmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |        Accelerate SparseCholesky |     0.24 secs |    0.035 secs |     1.23452e-15 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |     0.28 secs |     0.03 secs |     1.52008e-15 |
| 🥉 3 |             Eigen::SimplicialLLT |     0.59 secs |    0.053 secs |     3.01865e-15 |
|    4 |            Eigen::SimplicialLDLT |      0.6 secs |    0.052 secs |     2.33654e-15 |
|    5 |              catamari::SparseLDL |     0.73 secs |    0.043 secs |     1.35589e-15 |
|    6 |                 Eigen::UmfPackLU |     0.69 secs |     0.19 secs |     2.71385e-16 |
|    7 |   Eigen::BiCGSTAB<IncompleteLUT> |     0.82 secs |     0.55 secs |     3.91637e-16 |
|    8 |                  Eigen::SparseLU |        2 secs |    0.068 secs |     2.38898e-15 |
|    9 |                     MA57 (symla) |      2.8 secs |    0.051 secs |      1.0304e-15 |
|   10 |    Eigen::CG<IncompleteCholesky> |     0.22 secs |      9.9 secs |     9.15631e-10 |

# Biharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |        Accelerate SparseCholesky |     0.57 secs |    0.078 secs |     1.19627e-14 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |      0.6 secs |    0.067 secs |     4.91755e-15 |
| 🥉 3 |                 Eigen::UmfPackLU |      2.1 secs |     0.72 secs |     5.27979e-16 |
|    4 |            Eigen::SimplicialLDLT |      4.5 secs |     0.19 secs |      7.8414e-15 |
|    5 |             Eigen::SimplicialLLT |      4.5 secs |     0.19 secs |     4.61287e-15 |
|    6 |              catamari::SparseLDL |      6.3 secs |     0.16 secs |     3.66465e-15 |
|    7 |   Eigen::BiCGSTAB<IncompleteLUT> |      6.1 secs |      1.8 secs |     5.55785e-16 |
|    8 |                  Eigen::SparseLU |       13 secs |     0.25 secs |     4.51828e-12 |
|    9 |                     MA57 (symla) |       18 secs |     0.12 secs |     4.95481e-15 |
|    - |    Eigen::CG<IncompleteCholesky> |           - |           - | skipped: did not actually succeed: backward error 0.002724 exceeds 1e-06 |

# Triharmonic

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |        Accelerate SparseCholesky |     0.99 secs |     0.12 secs |     5.94182e-13 |
| 🥈 2 |      Eigen::CholmodSupernodalLLT |      1.2 secs |     0.11 secs |     2.38957e-13 |
| 🥉 3 |             Eigen::SimplicialLLT |       15 secs |     0.41 secs |     5.90259e-13 |
|    4 |            Eigen::SimplicialLDLT |       15 secs |     0.41 secs |     3.99303e-13 |
|    5 |              catamari::SparseLDL |       27 secs |     0.37 secs |     9.05897e-13 |
|    6 |                  Eigen::SparseLU |       55 secs |     0.68 secs |     5.65782e-09 |
|    7 |                     MA57 (symla) |       71 secs |     0.22 secs |     7.25797e-14 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error is NaN (solver diverged) |
|    - |    Eigen::CG<IncompleteCholesky> |           - |           - | skipped: did not actually succeed: backward error 0.000616 exceeds 1e-06 |

# Mixed Biharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |                 Eigen::UmfPackLU |      4.2 secs |      0.9 secs |     3.87298e-16 |
| 🥈 2 |            Eigen::SimplicialLDLT |      5.7 secs |     0.23 secs |      1.2478e-07 |
| 🥉 3 |     catamari::SparseLDL (LDLᵀ) |      6.4 secs |     0.18 secs |     4.69765e-08 |
|    4 |                  Eigen::SparseLU |       11 secs |     0.25 secs |     3.26471e-14 |
|    5 |             MA57 (symla, LDLᵀ) |       19 secs |     0.18 secs |     2.25769e-08 |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |        Accelerate SparseCholesky |           - |           - | skipped: factorization failed: status -1 |
|    - |    Eigen::CG<IncompleteCholesky> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 0.0005097 exceeds 1e-06 |

# Mixed Triharmonic (unflattened, indefinite)

| Rank |                          Method |      Factor |       Solve | Backward error |
|-----:|--------------------------------:|------------:|------------:|----------------:|
| 🥇 1 |     catamari::SparseLDL (LDLᵀ) |       24 secs |     0.42 secs |     1.55308e-09 |
| 🥈 2 |                  Eigen::SparseLU |       42 secs |      1.1 secs |     3.47129e-12 |
| 🥉 3 |             MA57 (symla, LDLᵀ) |       80 secs |     0.38 secs |     2.65418e-10 |
|    - |                 Eigen::UmfPackLU |           - |           - | skipped: known crash risk: excessive MKL thread churn on this system's fill-in |
|    - |            Eigen::SimplicialLDLT |           - |           - | skipped: known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale |
|    - |             Eigen::SimplicialLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |      Eigen::CholmodSupernodalLLT |           - |           - | skipped: factorization failed: NumericalIssue (not SPD/singular?) |
|    - |        Accelerate SparseCholesky |           - |           - | skipped: factorization failed: status -1 |
|    - |    Eigen::CG<IncompleteCholesky> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
|    - |   Eigen::BiCGSTAB<IncompleteLUT> |           - |           - | skipped: did not actually succeed: backward error 1 exceeds 1e-06 |
