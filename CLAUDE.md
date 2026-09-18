# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-file (`main.cpp`) informal benchmark/leaderboard comparing sparse
solvers on _k_-harmonic diffusion problems over triangle meshes. Two families
of systems: the **flattened** SPD systems `(Wᵏ + M) u = M x` for k=1,2,3
(harmonic/biharmonic/triharmonic, via `igl::harmonic`), where `L` is the
cotangent Laplacian, `M` is the mass matrix, and `Wᵏ⁺¹ = Wᵏ M⁻¹ L`; and the
**mixed** (unflattened) formulation's genuinely symmetric **indefinite**
saddle-point systems (`build_mixed_system()`, k=4,5 internally), which most
SPD-only solvers correctly refuse to factor. Each solver factors and solves
its system(s) and the benchmark prints factor time, solve time, and
componentwise relative backward error (see README's "How is accuracy
measured?" section) per solver per system.

## Clone

This repo relies on several git submodules that must be cloned recursively:

    git clone --recursive https://github.com/alecjacobson/sparse-solver-benchmark

Submodules: `libigl`, `SuiteSparse` (fork at sergiud/SuiteSparse), `catamari`,
`mantis`, `quotient`, `nasoq`, `ma57` (symla). All except `libigl` are optional
at configure time (see CMake options below) — a missing/uninitialized
submodule just disables that solver with a `message(WARNING ...)`, it doesn't
fail `cmake` configure.

## Build

    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=Release
    make -j<N>

Always pass an explicit, bounded `-j<N>` to `make`/`cmake --build` — this project
links SuiteSparse, MKL, and libigl and an unbounded parallel build can exhaust
machine resources.

Key CMake options (in `CMakeLists.txt`, all default ON, all autoconfigure —
each probes its own submodule/system dependency and turns itself OFF with a
warning if unavailable, rather than hard-failing `cmake` configure):
- `IGL_WITH_CHOLMOD` — build/link SuiteSparse (CHOLMOD + UMFPACK) as the fastest solver path.
- `IGL_WITH_MKL` — link Intel MKL and enable `Eigen::PardisoLLT`/`PardisoLDLT`.
  Deliberately does NOT define `EIGEN_USE_MKL_ALL` — see the `[!WARNING]` in
  README about why (a documented ~100x slowdown it causes in `SparseLU`).
- `IGL_WITH_GPL` — allow GPL-licensed Eigen sparse solver code; if off, falls back to slower `SparseLU`.
- `IGL_WITH_CUDSS` — NVIDIA cuDSS/cuSOLVER GPU solvers (requires CUDA toolkit + `libcudss0-dev-cuda-12`).
  On a GPU dev machine without a full system/apt CUDA install, `nvcc` and the
  CUDA math libraries (cuSOLVER/cuSPARSE/cuDSS) can come from several
  different places at once — e.g. this sort of layout has been seen on dev
  containers: `nvcc` itself assembled under `$HOME` (e.g. `~/warp-cuda/<ver>/bin`,
  for Warp's native builds — see CLAUDE.md's global "Warp / CUDA" section),
  `cusolver`/`cusparse` from real apt `libcusolver-dev`/`libcusparse-dev`
  packages in `/usr/local/cuda-*`, and `cudss` only from a pip
  `nvidia-cudss-cu12` wheel (no apt package for it exists at all). CMake's
  standard `find_package(CUDAToolkit)` only looks next to `nvcc`, so it won't
  find components that live in a different tree. `CMakeLists.txt`'s
  `IGL_WITH_CUDSS` block handles this by: (1) extending `nvcc` search hints
  beyond `/usr/local/cuda` to `$HOME/warp-cuda/*/bin`, `$HOME/projects/cuda*/bin`,
  and `$ENV{CUDA_ROOT}` (a convention some sibling projects on such machines
  already use); (2) falling back to `python3 -c "import nvidia, os; ...
  os.path.dirname(nvidia.__file__)"` to locate any pip-installed
  `nvidia-*-cu12` packages when `CUDA::cusolver`/`CUDA::cusparse` targets
  don't exist and/or cuDSS isn't found via the normal system search, using
  `file(GLOB ...)` instead of `find_library` since those wheels ship only
  version-suffixed `.so` files (e.g. `libcusolver.so.11`, no unversioned
  symlink) that `find_library`'s naming convention won't match. All of this
  is additive/best-effort — a machine with none of these has this code do
  nothing and fall through to the existing graceful-disable warning.
- `IGL_WITH_NASOQ` — NASOQ's LBL symmetric-indefinite solver (requires `IGL_WITH_MKL`).
- `IGL_WITH_CATAMARI` — catamari's SparseLDL (needs `mantis`/`quotient` submodules too).
- `IGL_WITH_MA57` — MA57 (symla), header-only symmetric-indefinite solver (needs system
  Eigen3 + SuiteSparse AMD, e.g. `apt install libeigen3-dev libsuitesparse-dev`).

Disabling `IGL_WITH_CHOLMOD` falls back to `SimplicialLLT`; disabling `IGL_WITH_GPL`
falls back to `SparseLU` (both are noted as slower in the CMake warnings).

## Run

    ./sparse_solver_benchmark [path to triangle mesh]
    ./sparse_solver_benchmark --grid N --check    # synthetic grid mesh, correctness check
    ./sparse_solver_benchmark --dragon N          # decimated real mesh, fast dev loop

Example mesh checked into the repo: `xyzrgb_dragon-720K.ply`.

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

For iterating on solver/tolerance logic, use `--dragon N` (decimates the
checked-in dragon mesh to ~N faces, caches under `meshes/`), NOT the full
720K-vertex mesh — using the full mesh mid-development is slow and has
previously led to multi-hour debugging loops (see git history around the
iterative-solver stopping-criteria commits). `--grid N` (a synthetic
triangulated grid) is even faster but trivially reorderable/banded, so it
doesn't exercise fill-in the way a real mesh does; prefer `--dragon` when
fill-in-dependent behavior (ordering, factor time) matters.

Output is a markdown table per system (Harmonic/Biharmonic/Triharmonic/Mixed
Biharmonic/Mixed Triharmonic) with columns: Method, Factor time, Solve time,
Backward error (componentwise relative, LAPACK's BERR — see README).
Full per-machine leaderboard runs live under `leaderboards/`, one file per
machine (see README's "Leaderboards" section).

## Architecture notes

- `main.cpp` is the entire program. There is a generic `solve<Factor>(...)`
  template that works for any Eigen sparse solver class implementing the usual
  `Factor(Q)` / `factor.solve(rhs)` interface (CholmodSupernodalLLT, UmfPackLU,
  SimplicialLLT, SimplicialLDLT, PardisoLLT, SparseLU, BiCGSTAB, ConjugateGradient).
- `catamari::SparseLDL<double>` and MA57 (`symla::SymLDLT<double>`) don't share
  Eigen's solver interface. Catamari gets a full template specialization of
  `solve<>` (manual `CoordinateMatrix`/`BlasMatrix` conversion, catamari's own
  `Factor`/`Solve` API); MA57 gets a standalone `solve_symla()` function (no
  `.info()`, uses `.isSingular()` instead, `.solve()` takes/returns a plain
  dense matrix) — same pattern NASOQ (`solve_nasoq_lbl`, if present) follows.
- To add a new solver: either instantiate `solve<YourEigenCompatibleSolver>(...)`
  in the `for(int k=1;k<=5;k++)` loop if it follows the Eigen factor/solve
  pattern, or add a new template specialization of `solve<>` / a standalone
  function if it has a different API (as catamari/MA57 do).
- Every optional solver's includes AND usage must both be `#ifdef IGL_WITH_*`
  guarded, not just the usage — an unguarded `#include <Eigen/CholmodSupport>`
  or `#include "catamari.hpp"` breaks compilation on a machine that disabled
  that solver, even if the call sites are properly guarded.
- `Eigen::PardisoLLT`/`PardisoLDLT` usage is gated behind `#ifdef IGL_WITH_MKL`.
- Timing uses a closured `Timer`/`tictoc()`-style helper (via
  `igl::get_seconds()`) reset before each factor/solve call. A small
  `#ifdef IGL_WITH_CHOLMOD` warm-up factorization runs before the main k-loop
  starts, so whichever solver happens to run first doesn't absorb one-time
  library-init/first-allocation costs into its own timed factor step.
