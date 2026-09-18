#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "warp-lang",
#     "scipy",
#     "numpy",
# ]
# ///
"""Times NVIDIA Warp's warp.optim.linear iterative solvers (cg, cr, bicgstab,
gmres) against the Q/rhs systems dumped by the C++ benchmark's
`--dump-matrices` flag.

Usage:
    cd ../build && ./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply \
        --dump-matrices /tmp/dump --dump-only
    uv run --with warp-lang --with scipy bench_warp.py --dir /tmp/dump

Timing/accuracy methodology (see the "fairness" note in the README): each
solver gets one untimed warm-up call first, to absorb Warp's kernel JIT
compilation and (with use_cuda_graph=True, the default) first-launch
CUDA-graph capture -- neither is representative of steady-state solver
cost, and both are one-time per-process costs a real caller would also pay
only once.

Accuracy is componentwise relative backward error (LAPACK's BERR, see
backward_error() below), not Warp's own tol= convergence test -- iterative
solvers are driven toward --berr-target by solving in chunks (each an
independent call with maxiter=chunk, x carrying over as the warm-start
guess) and checking backward_error() ourselves after each chunk, so the
decision to keep going is a fair, externally-verified one rather than each
library's private notion of "converged". This does mean a host sync (and a
scipy-side backward-error computation) once per chunk, not zero host syncs
for the whole solve -- chunk sizes are adaptive (~5s of work each) so this
overhead stays small relative to total solve time. Within a single chunk,
check_every=0 still runs it as one CUDA graph replay with no host
synchronization inside that chunk -- the fastest-possible, most GPU-native
way to run these solvers (confirmed empirically on this GPU that
check_every=0's device-side conditional-graph early exit genuinely works,
not just always running to maxiter as warp.optim.linear's own docstrings
warn can happen without that hardware support). Each chunk call passes
Warp's own tol=1e-14 -- deliberately far tighter than --berr-target, so
Warp's internal check essentially never fires before the chunk's requested
iteration count is used up; the loop's real stopping decision is always
our own external backward_error() check between chunks, not Warp's.
"""
import argparse
import math
import os
import time

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

import warp as wp
import warp.optim.linear as wol

SYSTEM_NAMES = {
    1: "Harmonic",
    2: "Biharmonic",
    3: "Triharmonic",
    4: "Mixed Biharmonic (unflattened, indefinite)",
    5: "Mixed Triharmonic (unflattened, indefinite)",
}

SOLVERS = {
    "warp::cg": wol.cg,
    "warp::cr": wol.cr,
    "warp::bicgstab": wol.bicgstab,
    "warp::gmres": wol.gmres,
}


def load_system(dump_dir, k):
    Q = sio.mmread(os.path.join(dump_dir, f"k{k}_Q.mtx")).tocsr()
    rhs = np.atleast_2d(sio.mmread(os.path.join(dump_dir, f"k{k}_rhs.mtx")))
    if rhs.shape[0] != Q.shape[0]:
        rhs = rhs.T
    return Q, rhs


def to_warp_matrix(Q, device):
    Qc = Q.tocoo()
    rows = wp.array(Qc.row.astype(np.int32), dtype=wp.int32, device=device)
    cols = wp.array(Qc.col.astype(np.int32), dtype=wp.int32, device=device)
    vals = wp.array(Qc.data.astype(np.float64), dtype=wp.float64, device=device)
    return wp.sparse.bsr_from_triplets(Q.shape[0], Q.shape[1], rows, cols, vals)


# Componentwise relative backward error (LAPACK's BERR), matching
# main.cpp's backward_error() exactly -- see its comment for the full
# rationale. eta_cw = max_i |b-Ax|_i / (|A||x|+|b|)_i. |A| is entrywise
# absolute value, not an induced norm; |A|*|x| is one sparse matrix-vector
# product (Q_abs, precomputed once per system since it doesn't depend on
# x, has identical sparsity to Q), same complexity as computing Q@x itself.
def backward_error(Q, Q_abs, b, x):
    r = np.abs(b - Q @ x)
    d = Q_abs @ np.abs(x) + np.abs(b)
    d = np.maximum(d, np.finfo(float).tiny)
    return float(np.max(r / d))


# kMaxIterativeIterations (20000) alone turned out to be impractical as the
# sole backstop on the C++ side -- see main.cpp's kIterativeTimeLimitSeconds
# comment: on the real dragon mesh, an iterative solver on the badly-scaled
# flattened triharmonic system was observed (via gdb) to genuinely still be
# computing after multiple *hours*, nowhere near either the tolerance or the
# iteration cap. Same fix here: an actual wall-clock deadline, matching the
# C++ side's value, checked between chunks of a warm-started solve rather
# than trusting iteration count alone to bound wall time.
TIME_LIMIT_SECONDS = 600.0  # matches main.cpp's kIterativeTimeLimitSeconds


# A genuinely non-convergent solve (e.g. cg on an indefinite system, which
# it isn't designed to handle) doesn't reliably improve at all -- it can
# oscillate or plateau rather than approach berr_target. Once max_total_iters
# was raised (see main.cpp's kMaxIterativeIterations comment) to give real,
# large, well-conditioned systems room to actually converge, that same
# generous budget became a liability for a small system that will NEVER
# converge. Detect the stall directly instead of just capping it: if
# backward error hasn't improved by at least STALL_IMPROVEMENT_FRACTION
# (relatively) over STALL_CHUNK_LIMIT consecutive chunks, stop -- a genuine
# "not making progress" signal, not an arbitrary budget, so it doesn't cut
# off systems that are actually still converging.
STALL_IMPROVEMENT_FRACTION = 0.01
STALL_CHUNK_LIMIT = 5
# Hard ceiling on chunk size regardless of the ~5s-of-work time estimate:
# on a cheap/tiny problem where each iteration costs microseconds, "~5s of
# work" can mean hundreds of thousands of iterations in a single chunk,
# which defeats stall detection above (it only checks BETWEEN chunks).
MAX_CHUNK_SIZE = 5000


def solve_time_limited(fn, A, Q, Q_abs, b_np, b, x, M, check_every, max_total_iters, berr_target):
    """Chunked, wall-clock-limited solve mirroring main.cpp's solve_rhs():
    repeatedly calls fn() with a growing iteration budget, x carrying over
    as the warm-start guess between chunks (per warp.optim.linear's own
    docs, x is "initial guess and solution vector"), checking elapsed time
    between chunks. The stopping criterion is backward_error(x) < berr_target
    -- an external check we compute ourselves after every chunk, not Warp's
    own internal tol= convergence test (see main.cpp's kBackwardErrorTarget
    comment for why: comparable accuracy across libraries requires a common
    external target, not each library's private notion of "converged").
    Returns (timed_out, iterations_used); x is updated in-place with the
    best-effort result either way."""
    deadline = time.perf_counter() + TIME_LIMIT_SECONDS
    chunk = 10
    iterations_used = 0
    best_berr = float("inf")
    stall_chunks = 0
    while True:
        t0 = time.perf_counter()
        # tol left at Warp's own tight default -- irrelevant to whether we
        # continue looping, since that decision is made below via our own
        # backward_error() check, but still passed so Warp's internal
        # bookkeeping (atol_sq, etc.) has a sane value.
        fn(A, b, x, tol=1e-14, maxiter=chunk, M=M, check_every=check_every, use_cuda_graph=True)
        wp.synchronize()
        chunk_dt = time.perf_counter() - t0
        iterations_used += chunk
        x_np = x.numpy()
        berr = backward_error(Q, Q_abs, b_np, x_np)
        # NaN unambiguously means the solve has diverged beyond any hope of
        # recovery -- stop immediately rather than folding it into stall
        # detection below (an earlier version treated NaN as "improved" to
        # avoid a spurious stall-counter increment, but that reset the
        # counter to 0 every chunk once NaN appeared, so it never fired,
        # and `berr < target` is also always False for NaN (IEEE 754) --
        # the loop had no way out except the full max_total_iters cap).
        if math.isnan(berr):
            return False, iterations_used
        if berr < berr_target or iterations_used >= max_total_iters:
            return False, iterations_used
        if berr < best_berr * (1.0 - STALL_IMPROVEMENT_FRACTION):
            best_berr = berr
            stall_chunks = 0
        else:
            stall_chunks += 1
            if stall_chunks >= STALL_CHUNK_LIMIT:
                return False, iterations_used
        now = time.perf_counter()
        if now >= deadline:
            return True, iterations_used
        remaining = max_total_iters - iterations_used
        # Target ~5s of work, or however much time is actually left before
        # the deadline if less, so we don't overshoot it by a wide margin
        # on a slow final chunk.
        target_secs = min(5.0, deadline - now)
        if chunk_dt > 1e-6:
            per_iter = chunk_dt / chunk
            chunk = max(10, min(remaining, int(target_secs / per_iter), MAX_CHUNK_SIZE))
        else:
            chunk = max(10, min(remaining, chunk * 4, MAX_CHUNK_SIZE))


def run_solver(fn, A, Q, Q_abs, rhs, device, maxiter, berr_target, check_every):
    """Solves every RHS column separately (warp's solvers take a single b/x
    vector at a time, same constraint NASOQ's solve_only() has on the C++
    side); returns (t_factor, t_solve, backward_error, timed_out, iterations_used).

    Uses Warp's Jacobi ("diag") preconditioner -- the closest thing Warp
    offers to the Eigen side's IncompleteLUT (Warp has no incomplete-LU
    preconditioner at all, only diag/diag_abs/id). Jacobi is a strictly
    weaker preconditioner than ILU, so this narrows but does not close the
    fairness gap between the two sides; see the README note next to the
    leaderboard tables."""
    n, nrhs = rhs.shape

    # Warm-up: absorb JIT compilation of both the preconditioner-construction
    # kernel and the solver's own kernels/CUDA-graph capture, none of which
    # are representative of steady-state cost (each only pays once per
    # process, cached across every later call/matrix size).
    wol.preconditioner(A, "diag")
    b0 = wp.array(rhs[:, 0].copy(), dtype=wp.float64, device=device)
    x0 = wp.zeros(n, dtype=wp.float64, device=device)
    fn(A, b0, x0, tol=1e-14, maxiter=10, M=wol.preconditioner(A, "diag"), check_every=check_every, use_cuda_graph=True)
    wp.synchronize()

    t0 = time.perf_counter()
    M = wol.preconditioner(A, "diag")
    wp.synchronize()
    t_factor = time.perf_counter() - t0

    U = np.zeros((n, nrhs))
    t_solve = 0.0
    timed_out = False
    iterations_used = 0
    for c in range(nrhs):
        b_np = rhs[:, c].copy()
        b = wp.array(b_np, dtype=wp.float64, device=device)
        x = wp.zeros(n, dtype=wp.float64, device=device)
        wp.synchronize()
        t0 = time.perf_counter()
        col_timed_out, col_iters = solve_time_limited(fn, A, Q, Q_abs, b_np, b, x, M, check_every, maxiter, berr_target)
        t_solve += time.perf_counter() - t0
        timed_out = timed_out or col_timed_out
        iterations_used += col_iters
        U[:, c] = x.numpy()

    berr = backward_error(Q, Q_abs, rhs, U)
    return t_factor, t_solve, berr, timed_out, iterations_used


# Same reclassification the C++ benchmark applies (see main.cpp's
# kBackwardErrorDivergedThreshold and print_leaderboard()): a solver can
# report "success" while its backward error is still far above what a real
# double-precision solve should achieve, which is a genuine failure to
# converge, not a data point worth ranking. Because backward error is
# scale-invariant (see backward_error() above), a single fixed threshold
# works uniformly across every system -- no per-k tuning, no comparison
# against what other solvers achieved needed.
DIVERGED_THRESHOLD = 1e-6  # matches main.cpp's kBackwardErrorDivergedThreshold


def print_leaderboard(rows):
    # rows: list of dicts with name/t_factor/t_solve/residual/timed_out/iterations_used
    #
    # NaN comparisons are always False in IEEE 754, so `nan > threshold` is
    # also False -- a NaN backward error needs an explicit isnan() check
    # rather than relying on the threshold comparison to catch it.
    ok = [r for r in rows if not math.isnan(r["residual"]) and r["residual"] <= DIVERGED_THRESHOLD]
    diverged = [r for r in rows if math.isnan(r["residual"]) or r["residual"] > DIVERGED_THRESHOLD]
    ok = sorted(ok, key=lambda r: r["t_factor"] + r["t_solve"])

    medals = ["\U0001F947", "\U0001F948", "\U0001F949"]
    any_timed_out = any(r["timed_out"] for r in ok)
    print("\n| Rank |          Method |      Factor |       Solve | Backward error |")
    print("|-----:|-----------------:|------------:|------------:|----------------:|")
    for i, r in enumerate(ok):
        rank = medals[i] + f" {i+1}" if i < 3 else f"   {i+1}"
        name = r["name"] + ("†" if r["timed_out"] else "")
        print(f"| {rank} | {name:>16} | {r['t_factor']:>9.2g} secs | {r['t_solve']:>9.2g} secs | {r['residual']:>15.6g} |")
    for r in diverged:
        name = r["name"]
        if math.isnan(r["residual"]):
            print(f"|    - | {name:>16} |           - |           - | "
                  f"skipped: did not actually succeed: backward error is NaN (solver diverged) |")
            continue
        print(f"|    - | {name:>16} |           - |           - | "
              f"skipped: did not actually succeed: backward error {r['residual']:.4g} exceeds {DIVERGED_THRESHOLD:.4g} |")
    if any_timed_out:
        print(f"\n† hit the {TIME_LIMIT_SECONDS/60:.0f}-minute iterative-solver time limit before reaching "
              f"--berr-target -- backward error is a snapshot at cutoff, not a converged result.")
    return diverged


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory of k<k>_Q.mtx/k<k>_rhs.mtx from --dump-matrices")
    ap.add_argument("--device", default=None, help="warp device (default: cuda:0 if available, else cpu)")
    ap.add_argument("--maxiter", type=int, default=1000000, help="safety-net cap, matches the C++ benchmark's kMaxIterativeIterations -- deliberately high (measured: 200000 iterations only took 26.6s on the dragon mesh's biharmonic system, reaching backward error ~2.85e-13) so --berr-target and the time limit, not this, are what usually decide")
    ap.add_argument("--berr-target", type=float, default=1e-8, help="componentwise relative backward error target iterative solvers are driven toward (checked externally after each chunk, not via Warp's own internal tol=), matching the C++ benchmark's kBackwardErrorTarget")
    ap.add_argument("--check-every", type=int, default=0, help="0 disables host-side convergence checks (pure CUDA-graph replay, but no early exit without device-side conditional graphs); >0 enables early exit at the cost of host syncs")
    ap.add_argument("--csv", default=None, help="write results in the same schema as the C++ benchmark's --csv")
    args = ap.parse_args()

    wp.init()
    device = args.device or ("cuda:0" if wp.get_cuda_device_count() > 0 else "cpu")
    print(f"# warp {wp.config.version} on {device}", flush=True)

    csv_f = open(args.csv, "w") if args.csv else None
    if csv_f:
        csv_f.write("k,method,factor_secs,solve_secs,backward_error,skipped,fused_factor,timed_out,iterations_used\n")

    for k in sorted(SYSTEM_NAMES):
        q_path = os.path.join(args.dir, f"k{k}_Q.mtx")
        if not os.path.exists(q_path):
            continue
        print(f"\n# {SYSTEM_NAMES[k]}")
        Q, rhs = load_system(args.dir, k)
        A = to_warp_matrix(Q, device)
        Q_abs = Q.copy()
        Q_abs.data = np.abs(Q_abs.data)

        rows = []
        for name, fn in SOLVERS.items():
            t_factor, t_solve, residual, timed_out, iterations_used = run_solver(
                fn, A, Q, Q_abs, rhs, device, args.maxiter, args.berr_target, args.check_every
            )
            rows.append({
                "name": name, "t_factor": t_factor, "t_solve": t_solve, "residual": residual,
                "timed_out": timed_out, "iterations_used": iterations_used,
            })
        diverged = print_leaderboard(rows)
        diverged_names = {r["name"] for r in diverged}
        if csv_f:
            for r in rows:
                skipped = 1 if r["name"] in diverged_names else 0
                csv_f.write(f"{k},{r['name']},{r['t_factor']:.9g},{r['t_solve']:.9g},{r['residual']:.9g},"
                             f"{skipped},0,{1 if r['timed_out'] else 0},{r['iterations_used']}\n")

    if csv_f:
        csv_f.close()


if __name__ == "__main__":
    main()
