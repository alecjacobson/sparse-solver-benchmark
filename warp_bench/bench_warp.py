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

Timing methodology (see the "fairness" note in the README): each solver gets
one untimed warm-up call first, to absorb Warp's kernel JIT compilation and
(with use_cuda_graph=True, the default) first-launch CUDA-graph capture --
neither is representative of steady-state solver cost, and both are one-time
per-process costs a real caller would also pay only once. The timed call is
bracketed by a single wp.synchronize() immediately before starting the clock
and one immediately after, matching how the C++ side's own Timer works
(wall-clock around a single blocking call) rather than summing per-iteration
host-device round trips. check_every=0 (default here) additionally disables
Warp's own host-side convergence checks so the whole solve runs as a single
CUDA graph replay with no host synchronization in the loop at all -- the
fastest-possible, most GPU-native way to run these solvers. The tradeoff:
without device-side conditional-graph support, check_every=0 means the solver
always runs to `maxiter` rather than early-exiting on convergence (see
warp.optim.linear's own docstrings) -- pass --check-every N>0 if you want
early exit and don't mind the extra host syncs that requires.
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


def run_solver(fn, A, Q_scipy, rhs, device, maxiter, tol, check_every):
    """Solves every RHS column separately (warp's solvers take a single b/x
    vector at a time, same constraint NASOQ's solve_only() has on the C++
    side); returns (t_factor, t_solve, linf_residual).

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
    fn(A, b0, x0, tol=tol, maxiter=maxiter, M=wol.preconditioner(A, "diag"), check_every=check_every, use_cuda_graph=True)
    wp.synchronize()

    t0 = time.perf_counter()
    M = wol.preconditioner(A, "diag")
    wp.synchronize()
    t_factor = time.perf_counter() - t0

    U = np.zeros((n, nrhs))
    t_solve = 0.0
    for c in range(nrhs):
        b = wp.array(rhs[:, c].copy(), dtype=wp.float64, device=device)
        x = wp.zeros(n, dtype=wp.float64, device=device)
        wp.synchronize()
        t0 = time.perf_counter()
        fn(A, b, x, tol=tol, maxiter=maxiter, M=M, check_every=check_every, use_cuda_graph=True)
        wp.synchronize()
        t_solve += time.perf_counter() - t0
        U[:, c] = x.numpy()

    residual = float(np.abs(rhs - Q_scipy @ U).max())
    return t_factor, t_solve, residual


# Same reclassification the C++ benchmark applies (see main.cpp's
# kDivergedFactor/kDivergedFloor and print_leaderboard()): a solver can
# report "success" while its residual is still enormous relative to what's
# actually achievable on this system, which is a real failure to converge,
# not a data point worth ranking. reference_residual, when given (from the
# C++ side's --csv, via --reference-csv below), is the true best across
# every solver on this k -- not just Warp's own four -- so a Warp solver
# that's uniformly bad relative to a good direct solve on the same system
# gets caught too, not just relative to its own (possibly also-bad) peers.
DIVERGED_FACTOR = 1e6
DIVERGED_FLOOR = 1e-3
# See main.cpp's kDivergedAbsoluteCap: on a hard system where even the best
# achieved residual is already large (e.g. the flattened triharmonic
# system), DIVERGED_FACTOR alone produces a threshold so loose that a truly
# nonsensical residual can slip through uncaught. Same value as the C++ side
# for consistency.
DIVERGED_ABSOLUTE_CAP = 1e4


def print_leaderboard(rows, reference_residual=None):
    # NaN comparisons are always False in IEEE 754, so `nan <= threshold` and
    # `nan > threshold` are both False -- a NaN residual would silently
    # vanish from both the "ok" and "diverged" lists below (and could also
    # corrupt min() below, which has no defined NaN-skipping behavior)
    # without this explicit filter. Route NaN rows straight to diverged.
    finite_rows = [r for r in rows if not math.isnan(r[3])]
    nan_rows = [r for r in rows if math.isnan(r[3])]

    residuals = [r[3] for r in finite_rows]
    candidates = residuals + ([reference_residual] if reference_residual is not None else [])
    min_residual = min(candidates) if candidates else float("inf")
    threshold = min(max(min_residual * DIVERGED_FACTOR, DIVERGED_FLOOR), DIVERGED_ABSOLUTE_CAP)

    ok = [r for r in finite_rows if r[3] <= threshold]
    diverged = [r for r in finite_rows if r[3] > threshold] + nan_rows
    ok = sorted(ok, key=lambda r: r[1] + r[2])

    medals = ["\U0001F947", "\U0001F948", "\U0001F949"]
    print("\n| Rank |          Method |      Factor |       Solve |     L∞ norm |")
    print("|-----:|-----------------:|------------:|------------:|------------:|")
    for i, (name, t_factor, t_solve, residual) in enumerate(ok):
        rank = medals[i] + f" {i+1}" if i < 3 else f"   {i+1}"
        print(f"| {rank} | {name:>16} | {t_factor:>9.2g} secs | {t_solve:>9.2g} secs | {residual:>10.6g} |")
    for name, t_factor, t_solve, residual in diverged:
        if math.isnan(residual):
            print(f"|    - | {name:>16} |           - |           - | "
                  f"skipped: did not actually succeed: L∞ residual is NaN (solver diverged) |")
            continue
        ratio = residual / min_residual if min_residual > 0 else float("inf")
        print(f"|    - | {name:>16} |           - |           - | "
              f"skipped: did not actually succeed: L∞ residual {residual:.4g} "
              f"is {ratio:.3g} x the best solver's ({min_residual:.4g}) on this system |")
    return diverged


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory of k<k>_Q.mtx/k<k>_rhs.mtx from --dump-matrices")
    ap.add_argument("--device", default=None, help="warp device (default: cuda:0 if available, else cpu)")
    ap.add_argument("--maxiter", type=int, default=20000, help="safety-net cap, matches the C++ benchmark's kMaxIterativeIterations -- deliberately high so --tol (the actual intended stopping criterion) is what usually decides, not this")
    ap.add_argument("--tol", type=float, default=1e-7, help="relative L2 residual tolerance (||b-Ax||_2 < tol*||b||_2), matches the C++ benchmark's kIterativeTolerance -- verified this is the same convergence formula Eigen's setTolerance() uses, so a given value means the same thing on both sides")
    ap.add_argument("--check-every", type=int, default=0, help="0 disables host-side convergence checks (pure CUDA-graph replay, but no early exit without device-side conditional graphs); >0 enables early exit at the cost of host syncs")
    ap.add_argument("--csv", default=None, help="write results in the same schema as the C++ benchmark's --csv")
    ap.add_argument("--reference-csv", default=None, help="C++ benchmark's own --csv output; when given, the diverged-vs-best check compares against the true best across ALL solvers on each k, not just Warp's own four")
    args = ap.parse_args()

    reference = {}
    if args.reference_csv:
        with open(args.reference_csv) as f:
            next(f)  # header
            for line in f:
                k, method, t_factor, t_solve, residual, skipped, fused = line.strip().split(",")
                if skipped == "0":
                    k = int(k)
                    reference[k] = min(reference.get(k, float("inf")), float(residual))

    wp.init()
    device = args.device or ("cuda:0" if wp.get_cuda_device_count() > 0 else "cpu")
    print(f"# warp {wp.config.version} on {device}", flush=True)

    csv_f = open(args.csv, "w") if args.csv else None
    if csv_f:
        csv_f.write("k,method,factor_secs,solve_secs,linf_residual,skipped,fused_factor\n")

    for k in sorted(SYSTEM_NAMES):
        q_path = os.path.join(args.dir, f"k{k}_Q.mtx")
        if not os.path.exists(q_path):
            continue
        print(f"\n# {SYSTEM_NAMES[k]}")
        Q, rhs = load_system(args.dir, k)
        A = to_warp_matrix(Q, device)

        rows = []
        for name, fn in SOLVERS.items():
            t_factor, t_solve, residual = run_solver(
                fn, A, Q, rhs, device, args.maxiter, args.tol, args.check_every
            )
            rows.append((name, t_factor, t_solve, residual))
        diverged = print_leaderboard(rows, reference.get(k))
        diverged_names = {r[0] for r in diverged}
        if csv_f:
            for name, t_factor, t_solve, residual in rows:
                skipped = 1 if name in diverged_names else 0
                csv_f.write(f"{k},{name},{t_factor:.9g},{t_solve:.9g},{residual:.9g},{skipped},0\n")

    if csv_f:
        csv_f.close()


if __name__ == "__main__":
    main()
