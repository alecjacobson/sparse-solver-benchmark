# warp_bench

Times [NVIDIA Warp](https://github.com/NVIDIA/warp)'s `warp.optim.linear`
iterative solvers (`cg`, `cr`, `bicgstab`, `gmres`) against the same
_k_-harmonic systems the C++ benchmark solves, without adding a Python/Warp
dependency to the main C++ build. The C++ side only gains a flag to dump its
`Q`/`rhs` matrices to plain [MatrixMarket](https://math.nist.gov/MatrixMarket/formats.html)
files; everything Warp-specific lives here.

## Usage

    cd ../build
    ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply --dump-matrices /tmp/dump --dump-only
    cd ../warp_bench
    uv run bench_warp.py --dir /tmp/dump

`--dump-only` skips this benchmark's own (C++) solvers entirely, so the dump
finishes in seconds regardless of mesh size -- combine with `--only`/
`--exclude` instead if you want the C++ side's own solvers to still run
alongside the dump.

`bench_warp.py` is a [PEP 723](https://peps.python.org/pep-0723/) script
(`uv run` resolves `warp-lang`/`scipy`/`numpy` automatically, no venv setup
needed). Add `--csv results.csv` to get output in the same schema as the C++
benchmark's own `--csv` (`k,method,factor_secs,solve_secs,linf_residual,
skipped,fused_factor`), so results from both sides can be concatenated into
one table.

## Fairness / timing methodology

See the docstring at the top of `bench_warp.py` for the full rationale. In
short: each solver gets one untimed warm-up call to absorb Warp's kernel JIT
compilation and CUDA-graph capture (both one-time-per-process costs, not
representative of steady-state solve cost), then the timed call is
bracketed by a single `wp.synchronize()` before and after -- no host syncs
inside the solve itself, matching `check_every=0` (the default here), which
runs the whole solve as one CUDA-graph replay.

## Why a separate add-on instead of wiring this into the C++ binary?

Warp is a Python/CUDA library with no C++ embedding story that fits this
project's existing Eigen-based `solve<Factor>()` pattern, and pulling in a
Python interpreter + Warp's own JIT/build system as a build-time dependency
of a small C++ CLI tool would be disproportionate. Dumping to a standard,
solver-agnostic file format instead means any other external tool (MATLAB,
a different Python library, ...) can reuse the same dumps without this
benchmark needing to know about it.
