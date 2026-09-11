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
benchmark's own `--csv` (`k,method,factor_secs,solve_secs,backward_error,
skipped,fused_factor,timed_out,iterations_used`), so results from both sides
can be concatenated into one table.

## Accuracy metric and fairness / timing methodology

Accuracy is reported as componentwise relative backward error (LAPACK's
BERR), computed the same way as the C++ side's `backward_error()` in
`main.cpp` -- see the main README's "How is accuracy measured?" section for
the full rationale (short version: it's scale-invariant, so a single fixed
threshold works for every system here, unlike an absolute residual).

Iterative solvers are driven toward `--berr-target` (default `1e-8`),
checked externally after every chunk of a time-boxed, warm-started solve --
not via Warp's own internal `tol=` convergence test, since comparable
accuracy across solvers/libraries requires one shared external target, not
each library's private notion of "converged". This does mean a host sync
(and a backward-error computation via scipy) every chunk, not zero host
syncs during the whole solve -- see the docstring at the top of
`bench_warp.py` for the rest of the timing methodology (warm-up call to
absorb JIT/CUDA-graph-capture cost, `wp.synchronize()` bracketing, etc.),
which otherwise still applies.

## Why a separate add-on instead of wiring this into the C++ binary?

Warp is a Python/CUDA library with no C++ embedding story that fits this
project's existing Eigen-based `solve<Factor>()` pattern, and pulling in a
Python interpreter + Warp's own JIT/build system as a build-time dependency
of a small C++ CLI tool would be disproportionate. Dumping to a standard,
solver-agnostic file format instead means any other external tool (MATLAB,
a different Python library, ...) can reuse the same dumps without this
benchmark needing to know about it.
