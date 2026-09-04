# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-file (`main.cpp`) informal benchmark comparing sparse positive-definite
linear solvers on _k_-harmonic diffusion problems over triangle meshes. It builds
`(Wᵏ + M) u = M x` for k=1,2,3 (harmonic/biharmonic/triharmonic), where `L` is the
cotangent Laplacian, `M` is the mass matrix, and `Wᵏ` is defined recursively as
`W¹ = L`, `Wᵏ⁺¹ = Wᵏ M⁻¹ L` (via `igl::harmonic`). Each solver factors and solves
this system and the benchmark prints factor time, solve time, and L∞ residual norm.

## Clone

This repo relies on several git submodules that must be cloned recursively:

    git clone --recursive https://github.com/alecjacobson/sparse-solver-benchmark

Submodules: `libigl`, `SuiteSparse` (fork at sergiud/SuiteSparse), `catamari`,
`mantis`, `quotient`.

## Build

    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=Release
    make -j<N>

Always pass an explicit, bounded `-j<N>` to `make`/`cmake --build` — this project
links SuiteSparse, MKL, and libigl and an unbounded parallel build can exhaust
machine resources.

Key CMake options (in `CMakeLists.txt`, all default ON):
- `IGL_WITH_CHOLMOD` — build/link SuiteSparse (CHOLMOD + UMFPACK) as the fastest solver path.
- `IGL_WITH_MKL` — link Intel MKL and enable `Eigen::PardisoLLT`, plus `EIGEN_USE_MKL_ALL`.
- `IGL_WITH_GPL` — allow GPL-licensed Eigen sparse solver code; if off, falls back to slower `SparseLU`.

Disabling `IGL_WITH_CHOLMOD` falls back to `SimplicialLLT`; disabling `IGL_WITH_GPL`
falls back to `SparseLU` (both are noted as slower in the CMake warnings).

## Run

    ./sparse_solver_benchmark [path to triangle mesh]

Example mesh checked into the repo: `xyzrgb_dragon-720K.ply`.

    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply

Output is a markdown table per k-value (Harmonic/Biharmonic/Triharmonic) with
columns: Method, Factor time, Solve time, L∞ norm (residual `‖rhs - Q*U‖∞`).

## Architecture notes

- `main.cpp` is the entire program. There is a generic `solve<Factor>(...)`
  template that works for any Eigen sparse solver class implementing the usual
  `Factor(Q)` / `factor.solve(rhs)` interface (CholmodSupernodalLLT, UmfPackLU,
  SimplicialLLT, SimplicialLDLT, PardisoLLT, SparseLU, BiCGSTAB, ConjugateGradient).
- `catamari::SparseLDL<double>` does not share Eigen's solver interface, so it has
  a full template specialization of `solve<>` that does manual
  `CoordinateMatrix`/`BlasMatrix` conversion and calls catamari's own
  `Factor`/`Solve` API.
- To add a new solver: either instantiate `solve<YourEigenCompatibleSolver>(...)`
  in the `for(int k=1;k<=3;k++)` loop if it follows the Eigen factor/solve
  pattern, or add a new template specialization of `solve<>` if it has a
  different API (as catamari does).
- `Eigen::PardisoLLT` usage is gated behind `#ifdef IGL_WITH_MKL`.
- Timing uses a closured `tictoc()` lambda (via `igl::get_seconds()`) reset
  before each factor/solve call.
