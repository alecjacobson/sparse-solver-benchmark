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
    side); returns (t_factor, t_solve, linf_residual)."""
    n, nrhs = rhs.shape

    t0 = time.perf_counter()
    b0 = wp.array(rhs[:, 0].copy(), dtype=wp.float64, device=device)
    x0 = wp.zeros(n, dtype=wp.float64, device=device)
    wp.synchronize()
    t_factor = time.perf_counter() - t0

    # Warm-up: absorb JIT compilation + first-launch CUDA-graph capture.
    fn(A, b0, x0, tol=tol, maxiter=maxiter, check_every=check_every, use_cuda_graph=True)
    wp.synchronize()

    U = np.zeros((n, nrhs))
    t_solve = 0.0
    for c in range(nrhs):
        b = wp.array(rhs[:, c].copy(), dtype=wp.float64, device=device)
        x = wp.zeros(n, dtype=wp.float64, device=device)
        wp.synchronize()
        t0 = time.perf_counter()
        fn(A, b, x, tol=tol, maxiter=maxiter, check_every=check_every, use_cuda_graph=True)
        wp.synchronize()
        t_solve += time.perf_counter() - t0
        U[:, c] = x.numpy()

    residual = float(np.abs(rhs - Q_scipy @ U).max())
    return t_factor, t_solve, residual


def print_leaderboard(rows):
    rows = sorted(rows, key=lambda r: r[1] + r[2])
    medals = ["\U0001F947", "\U0001F948", "\U0001F949"]
    print("\n| Rank |          Method |      Factor |       Solve |     L∞ norm |")
    print("|-----:|-----------------:|------------:|------------:|------------:|")
    for i, (name, t_factor, t_solve, residual) in enumerate(rows):
        rank = medals[i] + f" {i+1}" if i < 3 else f"   {i+1}"
        print(f"| {rank} | {name:>16} | {t_factor:>9.2g} secs | {t_solve:>9.2g} secs | {residual:>10.6g} |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="directory of k<k>_Q.mtx/k<k>_rhs.mtx from --dump-matrices")
    ap.add_argument("--device", default=None, help="warp device (default: cuda:0 if available, else cpu)")
    ap.add_argument("--maxiter", type=int, default=200, help="matches the C++ benchmark's kMaxIterativeIterations cap")
    ap.add_argument("--tol", type=float, default=None, help="relative tolerance (default: warp's own dtype-based default, ~1e-12 for float64)")
    ap.add_argument("--check-every", type=int, default=0, help="0 disables host-side convergence checks (pure CUDA-graph replay, but no early exit without device-side conditional graphs); >0 enables early exit at the cost of host syncs")
    ap.add_argument("--csv", default=None, help="write results in the same schema as the C++ benchmark's --csv")
    args = ap.parse_args()

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
            if csv_f:
                csv_f.write(f"{k},{name},{t_factor:.9g},{t_solve:.9g},{residual:.9g},0,0\n")
        print_leaderboard(rows)

    if csv_f:
        csv_f.close()


if __name__ == "__main__":
    main()
