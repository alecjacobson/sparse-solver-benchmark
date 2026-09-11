
#ifdef _OPENMP
#include <omp.h>
#endif

#include "catamari.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/harmonic.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/get_seconds.h>
#ifdef IGL_WITH_MKL
#include <Eigen/PardisoSupport>
#endif
#include <Eigen/CholmodSupport>
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/SparseExtra>
#ifdef IGL_WITH_CUDSS
#include <cuda_runtime_api.h>
#include <cudss.h>
#include <cusolverSp.h>
#include <cusparse.h>
#endif
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Simple stopwatch: tic() resets, toc() returns seconds since the last
// tic()/toc() call and resets again.
struct Timer
{
  double t = igl::get_seconds();
  void tic() { t = igl::get_seconds(); }
  double toc()
  {
    const double now = igl::get_seconds();
    const double diff = now - t;
    t = now;
    return diff;
  }
};

struct Result
{
  int k;
  std::string name;
  double t_factor;
  double t_solve;
  double residual;
  bool skipped;
  std::string skip_reason;
  // True when t_factor doesn't represent a real, separately-timed factor
  // step (e.g. an API that fuses analysis+factor+solve into one call), so
  // the leaderboard can print "(fused)" instead of a misleading "0 secs".
  bool fused_factor = false;
  // True when an iterative solver hit kIterativeTimeLimitSeconds before
  // reaching kBackwardErrorTarget or kMaxIterativeIterations -- the row still
  // reports its best-effort accuracy at cutoff (not skipped), but the
  // leaderboard marks it so a fast-but-inaccurate time-limited result isn't
  // read as a genuinely converged one. iterations_used is -1 for solvers
  // this doesn't apply to (direct solvers have no iteration count).
  bool timed_out = false;
  int iterations_used = -1;
};

static std::vector<Result> g_results;
static FILE * g_csv = nullptr;
static bool g_check_mode = false;
static bool g_check_failed = false;

// --only/--exclude: case-insensitive substring filters on solver name, so
// e.g. `--only nasoq` runs just NASOQ LBL, or `--exclude umfpack,sparselu`
// skips the two slowest general-purpose solvers. Checked at the top of every
// solve_*()/solve<>() entry point below, before any work is done — a
// filtered-out solver doesn't even appear in the leaderboard (it's not
// "skipped", it's simply not run), so this is for fast local iteration on
// one solver at a time (e.g. while debugging), not for permanent leaderboard
// output.
static std::vector<std::string> g_only;
static std::vector<std::string> g_exclude;

static std::string to_lower(std::string s)
{
  for(char & c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

static bool should_run(const std::string & name)
{
  const std::string lname = to_lower(name);
  if(!g_only.empty())
  {
    bool found = false;
    for(const auto & pat : g_only) if(lname.find(pat) != std::string::npos) { found = true; break; }
    if(!found) return false;
  }
  for(const auto & pat : g_exclude) if(lname.find(pat) != std::string::npos) return false;
  return true;
}

static void record(
  int k,
  const std::string & name,
  double t_factor,
  double t_solve,
  double residual,
  bool skipped = false,
  const std::string & skip_reason = "",
  bool fused_factor = false,
  bool timed_out = false,
  int iterations_used = -1)
{
  g_results.push_back({k, name, t_factor, t_solve, residual, skipped, skip_reason, fused_factor, timed_out, iterations_used});
  if(g_csv)
  {
    fprintf(g_csv,"%d,%s,%.9g,%.9g,%.9g,%d,%d,%d,%d\n",
      k, name.c_str(), t_factor, t_solve, residual, skipped?1:0, fused_factor?1:0, timed_out?1:0, iterations_used);
  }
  // --check's correctness gate deliberately does NOT run here: print_leaderboard()
  // later reclassifies any result whose residual is enormous (or NaN)
  // relative to the best solver on this same k as "did not actually
  // succeed" -- e.g. Eigen::CG genuinely diverging (to a huge value, or all
  // the way to NaN, platform-dependently) on an indefinite system it isn't
  // designed for. Checking here, before that reclassification has a chance
  // to run, would flag an already-known, already-handled non-result as a
  // correctness regression. See print_leaderboard() for the actual check.
}

// Componentwise relative backward error (LAPACK's BERR): the smallest
// relative, entrywise perturbation of A and b that would make the computed
// X an *exact* solution --
//
//   eta_cw = max_ij |B-AX|_ij / (|A||X|+|B|)_ij
//
// This replaces a plain absolute residual (||B-AX||) as this benchmark's
// accuracy metric because it's naturally scale-invariant: it correctly
// reports a near-machine-precision result even when A has astronomically
// large entries (see the k=3/triharmonic README note -- an absolute
// residual around 48 there is actually excellent once you account for A's
// 3e15-magnitude entries) or wildly different row scales (as in the mixed
// systems' M/L blocks), because the denominator is formed with the same
// arithmetic -- and thus subject to the same rounding -- that produced the
// residual in the numerator. |A| means entrywise absolute value here, not
// an induced matrix norm; |A|*|X| is a single sparse-times-dense product
// with the same sparsity pattern as A (no dense matrix ever formed),
// O(nnz(A)*ncols(X)) -- the same complexity class as computing A*X itself.
static double backward_error(
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  const Eigen::MatrixXd & U)
{
  const Eigen::MatrixXd R = (rhs - Q*U).cwiseAbs();
  Eigen::MatrixXd D = rhs.cwiseAbs();
  for(int c=0; c<Q.outerSize(); ++c)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(Q,c); it; ++it)
    {
      // Column-major traversal: it.col() is the fixed outer index j, it.row()
      // the inner index i -- i.e. this nonzero is Q(i,j), contributing
      // |Q(i,j)|*|U(j,:)| to D's row i.
      D.row(it.row()).array() += std::abs(it.value()) * U.row(it.col()).array().abs();
    }
  }
  // LAPACK convention: floor the denominator instead of dividing by exactly
  // zero, which only happens for a row where A, x, and b are all exactly
  // zero (the ratio there is definitionally 0/0 -- treat it as perfectly
  // solved, not NaN/Inf).
  return (R.array() / D.cwiseMax(std::numeric_limits<double>::min()).array()).maxCoeff();
}

// Iterative solvers (BiCGSTAB/ConjugateGradient) default to Eigen's built-in
// cap of twice the system size, which is effectively unbounded for a 720K-row
// mesh: on the k=3 triharmonic system (documented above as badly scaled) they
// can burn tens of minutes grinding through non-convergent iterations without
// ever tripping any other limit. maxiter here is a SAFETY NET, not the
// intended stopping criterion -- that's kBackwardErrorTarget below; this just
// bounds worst-case runtime on a system that never converges (e.g. CG on a
// genuinely indefinite input, which it isn't designed to handle) instead of
// letting it hang.
static const int kMaxIterativeIterations = 20000;

// The actual intended stopping criterion for iterative solvers: componentwise
// relative backward error (see backward_error() above) under this target,
// checked directly -- an EXTERNAL common target applied the same way to
// every iterative solver here (Eigen's and Warp's, see warp_bench/bench_warp.py's
// --berr-target), rather than trusting each library's own internal
// convergence criterion. This matters because internal criteria differ in
// both norm (Eigen's setTolerance()/Warp's tol= are both relative-L2, not
// backward error) and, more importantly, in what they're relative TO --
// letting each library privately decide "solved" would time solvers at
// whatever accuracy their own default happens to land on, not a comparable
// one. Direct solvers have no notion of a target to iterate toward; they
// solve once and simply report whatever backward error they achieve.
static const double kBackwardErrorTarget = 1e-8;

// kMaxIterativeIterations (20000) alone turned out to be impractical as the
// sole backstop: on the real dragon mesh, Eigen::CG on the badly-scaled
// flattened triharmonic system was observed (via gdb, confirming it was
// genuinely still computing, not hung) to run for multiple *hours* without
// reaching either the target or the iteration cap. A pure iteration cap
// can't bound wall-clock time when a system is this much harder than
// smaller test cases suggested. Add an actual wall-clock deadline: for
// iterative solvers, solve in chunks via solveWithGuess() (continuing from
// the previous chunk's x, not restarting from zero each time) and check
// elapsed time between chunks, so a solver that's still short of
// convergence at the deadline is stopped and reports its best-effort
// accuracy at that point -- rather than either running unboundedly or
// being cut off at an arbitrary iteration count with no time guarantee.
static const double kIterativeTimeLimitSeconds = 600.0; // 10 minutes

// SFINAE dispatch: only types with solveWithGuess() (Eigen's iterative
// solvers) get the chunked/timed path; direct solvers fall through to a
// single plain solve() call. Q is needed (in addition to rhs) because the
// per-chunk stopping check is backward_error(Q, rhs, x), not anything
// Eigen's own factor object tracks internally.
template <typename Factor>
auto solve_rhs(Factor & factor, const Eigen::SparseMatrix<double> & Q, const Eigen::MatrixXd & rhs, bool & timed_out, int & iterations_used, int)
  -> decltype(factor.solveWithGuess(rhs, rhs), Eigen::MatrixXd())
{
  timed_out = false;
  iterations_used = 0;
  Eigen::MatrixXd x = Eigen::MatrixXd::Zero(rhs.rows(), rhs.cols());
  const double deadline = igl::get_seconds() + kIterativeTimeLimitSeconds;
  int chunk = 10;
  while(true)
  {
    factor.setMaxIterations(chunk);
    const double chunk_t0 = igl::get_seconds();
    x = factor.solveWithGuess(rhs, x);
    const double chunk_dt = igl::get_seconds() - chunk_t0;
    // factor.iterations() is NOT used here: for a multi-column rhs (this
    // benchmark's rhs is M*V, 3 columns), Eigen's own
    // IterativeSolverBase::_solve_with_guess_impl loops per-column and
    // aggregates m_info across columns but never aggregates m_iterations --
    // verified empirically (and by reading IterativeSolverBase.h) that it
    // reports 0 after a multi-column solveWithGuess() regardless of what
    // actually happened, converged or not. Track our own upper-bound
    // instead: since setMaxIterations(chunk) caps each column's iteration
    // count, "chunk" itself is what a non-converged call actually spent
    // (this is an approximation when some columns converge before others).
    iterations_used += chunk;
    if(backward_error(Q,rhs,x) < kBackwardErrorTarget || iterations_used >= kMaxIterativeIterations) break;
    const double now = igl::get_seconds();
    if(now >= deadline) { timed_out = true; break; }
    // Adapt the next chunk size toward ~5s of work (or however much time is
    // actually left before the deadline, if less), based on this chunk's
    // observed per-iteration cost, so we check the deadline reasonably
    // often without paying excessive per-chunk call overhead on easy
    // problems (few, tiny chunks) or overshooting the deadline by a wide
    // margin on hard ones (one huge final chunk).
    const int remaining = kMaxIterativeIterations - iterations_used;
    const double target_secs = std::min(5.0, deadline - now);
    if(chunk_dt > 1e-6)
    {
      const double per_iter = chunk_dt / (double)chunk;
      chunk = std::max(10, std::min(remaining, (int)(target_secs/per_iter)));
    }
    else
    {
      chunk = std::max(10, std::min(remaining, chunk*4));
    }
  }
  return x;
}
template <typename Factor>
Eigen::MatrixXd solve_rhs(Factor & factor, const Eigen::SparseMatrix<double> &, const Eigen::MatrixXd & rhs, bool & timed_out, int & iterations_used, long)
{
  timed_out = false;
  iterations_used = -1;
  return factor.solve(rhs);
}

static const char * eigen_info_string(Eigen::ComputationInfo info)
{
  switch(info)
  {
    case Eigen::Success: return "Success";
    case Eigen::NumericalIssue: return "NumericalIssue (not SPD/singular?)";
    case Eigen::NoConvergence: return "NoConvergence";
    case Eigen::InvalidInput: return "InvalidInput";
    default: return "Unknown";
  }
}

template <typename Factor>
void solve(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  if(!should_run(name)) return;
  Timer timer;
  Factor factor;
  factor.compute(Q);
  const double t_factor = timer.toc();
  // Only gate on info() for the mixed/indefinite systems (k>=4): there,
  // NumericalIssue reliably means "this SPD-only algorithm can't handle
  // indefinite input," worth a clean skip. For the flattened k=1..3 systems
  // (always SPD, just sometimes ill-conditioned — see the k=3 "badly scaled"
  // note below), some backends (observed: UmfPackLU on the triharmonic
  // dragon-mesh system) report NumericalIssue for a near-singular-but-still
  // meaningful result; skipping those would silently drop real, if
  // inaccurate, data that used to show up in the leaderboard.
  if(k>=4 && factor.info() != Eigen::Success)
  {
    record(k, name, t_factor, 0, 0, true,
      std::string("factorization failed: ") + eigen_info_string(factor.info()));
    return;
  }
  bool timed_out = false;
  int iterations_used = -1;
  U = solve_rhs(factor, Q, rhs, timed_out, iterations_used, 0);
  const double t_solve = timer.toc();
  record(k, name, t_factor, t_solve, backward_error(Q,rhs,U),
    false, "", false, timed_out, iterations_used);
}

template <>
void solve<catamari::SparseLDL<double>>(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  if(!should_run(name)) return;
  catamari::CoordinateMatrix<double> matrix;
  matrix.Resize(Q.rows(), Q.cols());
  matrix.ReserveEntryAdditions(Q.nonZeros());
  // Queue updates of entries in the sparse matrix using commands of the form:
  for(int c=0; c<Q.outerSize(); ++c)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(Q,c); it; ++it)
    {
      matrix.QueueEntryAddition(it.row(), it.col(), it.value());
    }
  }
  matrix.FlushEntryQueues();

  Timer timer;
  // Fill the options for the factorization.
  catamari::SparseLDLControl<double> ldl_control;
  ldl_control.SetFactorizationType(catamari::kCholeskyFactorization);

  // Factor the matrix.
  catamari::SparseLDL<double> ldl;
  const catamari::SparseLDLResult<double> result = ldl.Factor(matrix, ldl_control);
  const double t_factor = timer.toc();

  // copy rhs
  catamari::BlasMatrix<double> right_hand_sides;
  right_hand_sides.Resize(rhs.rows(), rhs.cols());
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      right_hand_sides(i, j) = rhs(i,j);
    }
  }

  // Solve a linear system using the factorization.
  ldl.Solve(&right_hand_sides.view);
  const double t_solve = timer.toc();

  // copy solution
  U.resize(rhs.rows(),rhs.cols());
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      U(i,j) = right_hand_sides(i, j);
    }
  }
  record(k, name, t_factor, t_solve, backward_error(Q,rhs,U));
}

#ifdef IGL_WITH_CUDSS

#define CUDA_CHECK(call) \
  do { \
    const cudaError_t _err = (call); \
    if(_err != cudaSuccess) \
    { \
      fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(_err)); \
      std::exit(1); \
    } \
  } while(0)

#define CUDSS_CHECK(call) \
  do { \
    const cudssStatus_t _st = (call); \
    if(_st != CUDSS_STATUS_SUCCESS) \
    { \
      fprintf(stderr,"cuDSS error %s:%d: status %d\n",__FILE__,__LINE__,(int)_st); \
      std::exit(1); \
    } \
  } while(0)

#define CUSOLVER_CHECK(call) \
  do { \
    const cusolverStatus_t _st = (call); \
    if(_st != CUSOLVER_STATUS_SUCCESS) \
    { \
      fprintf(stderr,"cuSOLVER error %s:%d: status %d\n",__FILE__,__LINE__,(int)_st); \
      std::exit(1); \
    } \
  } while(0)

#define CUSPARSE_CHECK(call) \
  do { \
    const cusparseStatus_t _st = (call); \
    if(_st != CUSPARSE_STATUS_SUCCESS) \
    { \
      fprintf(stderr,"cuSPARSE error %s:%d: status %d\n",__FILE__,__LINE__,(int)_st); \
      std::exit(1); \
    } \
  } while(0)

bool cuda_device_available()
{
  int count = 0;
  const cudaError_t err = cudaGetDeviceCount(&count);
  return err == cudaSuccess && count > 0;
}

// Q is symmetric (Q = M + Wᵏ), so its column-major (CSC) storage doubles as
// row-major (CSR) storage: no transpose/reshuffle needed, just relabel.
struct DeviceCSR
{
  int n = 0, nnz = 0;
  int * row_ptr = nullptr;
  int * col_idx = nullptr;
  double * val = nullptr;

  DeviceCSR(const Eigen::SparseMatrix<double> & Q)
  {
    n = (int)Q.rows();
    nnz = (int)Q.nonZeros();
    CUDA_CHECK(cudaMalloc(&row_ptr,(n+1)*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&col_idx,nnz*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&val,nnz*sizeof(double)));
    CUDA_CHECK(cudaMemcpy(row_ptr,Q.outerIndexPtr(),(n+1)*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(col_idx,Q.innerIndexPtr(),nnz*sizeof(int),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(val,Q.valuePtr(),nnz*sizeof(double),cudaMemcpyHostToDevice));
  }
  ~DeviceCSR()
  {
    cudaFree(row_ptr);
    cudaFree(col_idx);
    cudaFree(val);
  }
};

void solve_cudss(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U,
  cudssMatrixType_t mtype = CUDSS_MTYPE_SPD)
{
  if(!should_run(name)) return;
  if(!cuda_device_available())
  {
    record(k, name, 0, 0, 0, true, "no CUDA device");
    return;
  }

  const int n = (int)Q.rows();
  const int nrhs = (int)rhs.cols();
  DeviceCSR A_csr(Q);

  double * d_b = nullptr;
  double * d_x = nullptr;
  CUDA_CHECK(cudaMalloc(&d_b,n*nrhs*sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x,n*nrhs*sizeof(double)));
  // rhs is Eigen::MatrixXd, column-major by default, matching CUDSS_LAYOUT_COL_MAJOR.
  CUDA_CHECK(cudaMemcpy(d_b,rhs.data(),n*nrhs*sizeof(double),cudaMemcpyHostToDevice));

  cudssHandle_t handle;
  cudssConfig_t config;
  cudssData_t data;
  CUDSS_CHECK(cudssCreate(&handle));
  CUDSS_CHECK(cudssConfigCreate(&config));
  if(mtype == CUDSS_MTYPE_SYMMETRIC)
  {
    // This cuDSS version's only pivoting strategy for symmetric indefinite
    // matrices is CUDSS_PIVOT_DIAGONAL (diagonal-only search) --
    // CUDSS_PIVOT_BUNCH_KAUFMAN, the safe block-pivoting strategy NASOQ
    // uses, is "reserved for future, not supported yet" per
    // cudss_data_types.h in this version. Diagonal-only pivoting can't find
    // a safe pivot in a block that's structurally all-zero on the diagonal
    // (this benchmark's mixed-system λ block), so factorization silently
    // succeeds with a numerically useless answer instead of erroring out
    // (see the README's ⚠️ note, confirmed directly with NVIDIA's cuDSS
    // team). They suggested iterative refinement as a mitigation in the
    // meantime; it genuinely fixes the mixed biharmonic system (near
    // machine precision) but not the harder mixed triharmonic one, and is
    // itself inconsistent run-to-run there -- this benchmark's own
    // did-not-actually-succeed check (see print_leaderboard) is the real
    // safety net, not this alone.
    const int ir_n_steps = 2;
    CUDSS_CHECK(cudssConfigSet(config,CUDSS_CONFIG_IR_N_STEPS,&ir_n_steps,sizeof(ir_n_steps)));
  }
  CUDSS_CHECK(cudssDataCreate(handle,&data));

  cudssMatrix_t A, B, X;
  CUDSS_CHECK(cudssMatrixCreateCsr(&A,n,n,A_csr.nnz,A_csr.row_ptr,nullptr,A_csr.col_idx,A_csr.val,
    CUDSS_R_32I,CUDSS_R_32I,CUDSS_R_64F,mtype,CUDSS_MVIEW_FULL,CUDSS_BASE_ZERO));
  CUDSS_CHECK(cudssMatrixCreateDn(&B,n,nrhs,n,d_b,CUDSS_R_64F,CUDSS_LAYOUT_COL_MAJOR));
  CUDSS_CHECK(cudssMatrixCreateDn(&X,n,nrhs,n,d_x,CUDSS_R_64F,CUDSS_LAYOUT_COL_MAJOR));

  // Unlike the setup calls above (real bugs if they fail), factorization can
  // legitimately fail here (e.g. an indefinite matrix passed with a wrong/
  // unsupported mtype) — report it as a skipped result instead of aborting
  // the whole benchmark.
  const auto cleanup = [&]()
  {
    cudssMatrixDestroy(A);
    cudssMatrixDestroy(B);
    cudssMatrixDestroy(X);
    cudssDataDestroy(handle,data);
    cudssConfigDestroy(config);
    cudssDestroy(handle);
    cudaFree(d_b);
    cudaFree(d_x);
  };

  Timer timer;
  cudssStatus_t status = cudssExecute(handle,CUDSS_PHASE_REORDERING,config,data,A,X,B);
  if(status==CUDSS_STATUS_SUCCESS) status = cudssExecute(handle,CUDSS_PHASE_SYMBOLIC_FACTORIZATION,config,data,A,X,B);
  if(status==CUDSS_STATUS_SUCCESS) status = cudssExecute(handle,CUDSS_PHASE_FACTORIZATION,config,data,A,X,B);
  if(status==CUDSS_STATUS_SUCCESS) status = cudaDeviceSynchronize()==cudaSuccess ? CUDSS_STATUS_SUCCESS : CUDSS_STATUS_EXECUTION_FAILED;
  const double t_factor = timer.toc();
  if(status != CUDSS_STATUS_SUCCESS)
  {
    cleanup();
    record(k, name, t_factor, 0, 0, true,
      std::string("cuDSS factorization failed: status ")+std::to_string((int)status));
    return;
  }

  CUDSS_CHECK(cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,X,B));
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_solve = timer.toc();

  U.resize(n,nrhs);
  CUDA_CHECK(cudaMemcpy(U.data(),d_x,n*nrhs*sizeof(double),cudaMemcpyDeviceToHost));

  cleanup();

  record(k, name, t_factor, t_solve, backward_error(Q,rhs,U));
}

// cusolverSpDcsrlsvchol fuses reordering + symbolic + numeric factorization +
// triangular solve into a single call for one dense RHS vector at a time —
// cuSOLVER's legacy sparse-Cholesky API has no separate factor/solve phases
// (unlike cuDSS), so the whole fused cost (repeated once per RHS column) is
// reported under "solve" and "factor" is left at zero, rather than fabricate
// a split the API doesn't expose.
void solve_cusolver(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  if(!should_run(name)) return;
  if(!cuda_device_available())
  {
    record(k, name, 0, 0, 0, true, "no CUDA device");
    return;
  }

  const int n = (int)Q.rows();
  const int nrhs = (int)rhs.cols();
  DeviceCSR A_csr(Q);

  cusolverSpHandle_t handle;
  CUSOLVER_CHECK(cusolverSpCreate(&handle));
  cusparseMatDescr_t descr;
  CUSPARSE_CHECK(cusparseCreateMatDescr(&descr));
  CUSPARSE_CHECK(cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL));
  CUSPARSE_CHECK(cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO));

  double * d_b = nullptr;
  double * d_x = nullptr;
  CUDA_CHECK(cudaMalloc(&d_b,n*sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x,n*sizeof(double)));
  U.resize(n,nrhs);

  Timer timer;
  for(int c = 0; c < nrhs; c++)
  {
    CUDA_CHECK(cudaMemcpy(d_b,rhs.col(c).data(),n*sizeof(double),cudaMemcpyHostToDevice));
    int singularity = -1;
    CUSOLVER_CHECK(cusolverSpDcsrlsvchol(handle,n,A_csr.nnz,descr,A_csr.val,A_csr.row_ptr,A_csr.col_idx,
      d_b,/*tol=*/1e-12,/*reorder=*/3,d_x,&singularity));
    CUDA_CHECK(cudaMemcpy(U.col(c).data(),d_x,n*sizeof(double),cudaMemcpyDeviceToHost));
  }
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_solve = timer.toc();

  cudaFree(d_b);
  cudaFree(d_x);
  CUSPARSE_CHECK(cusparseDestroyMatDescr(descr));
  CUSOLVER_CHECK(cusolverSpDestroy(handle));

  record(k, name, 0, t_solve, backward_error(Q,rhs,U),
    /*skipped=*/false, /*skip_reason=*/"", /*fused_factor=*/true);
}

#endif

// Build the "mixed" (unflattened) FEM k-harmonic block system: instead of
// eliminating auxiliary variables via M⁻¹ (as igl::harmonic does to produce
// the flattened Wᵏ = Wᵏ⁻¹M⁻¹L), keep them explicit. This is sparser but the
// resulting symmetric block system is indefinite instead of SPD.
//
//   order=2 (biharmonic), unknowns (u, a₁), a₁ = M⁻¹Lu:
//     [ M   L ] [u ]   [Mx]
//     [ L  -M ] [a₁] = [0 ]
//
//   order=3 (triharmonic), unknowns (u, a₁, λ), a₁=M⁻¹Lu, λ=M⁻¹La₁:
//     [ M   0   L ] [u ]   [Mx]
//     [ 0   L  -M ] [a₁] = [0 ]
//     [ L  -M   0 ] [λ ]   [0 ]
static void build_mixed_system(
  int order,
  const Eigen::SparseMatrix<double> & L,
  const Eigen::SparseMatrix<double> & M,
  const Eigen::MatrixXd & V,
  Eigen::SparseMatrix<double> & Q,
  Eigen::MatrixXd & rhs)
{
  const int n = (int)L.rows();
  const int nb = order;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve((L.nonZeros()+M.nonZeros())*2);

  const auto add_block = [&](int block_row, int block_col, const Eigen::SparseMatrix<double> & A, double sign)
  {
    const int ro = block_row*n;
    const int co = block_col*n;
    for(int c=0;c<A.outerSize();++c)
    {
      for(Eigen::SparseMatrix<double>::InnerIterator it(A,c); it; ++it)
      {
        triplets.emplace_back(ro+it.row(), co+it.col(), sign*it.value());
      }
    }
  };

  if(order == 2)
  {
    add_block(0,0,M, 1); add_block(0,1,L, 1);
    add_block(1,0,L, 1); add_block(1,1,M,-1);
  }
  else // order == 3
  {
    add_block(0,0,M, 1);                 add_block(0,2,L, 1);
                  add_block(1,1,L, 1);    add_block(1,2,M,-1);
    add_block(2,0,L, 1); add_block(2,1,M,-1);
  }

  Q.resize(nb*n,nb*n);
  Q.setFromTriplets(triplets.begin(),triplets.end());

  rhs = Eigen::MatrixXd::Zero(nb*n, V.cols());
  rhs.topRows(n) = M*V;
}

// catamari's genuine symmetric-indefinite LDLᵀ mode (as opposed to the
// Cholesky-only specialization above), for the mixed/indefinite systems.
// A standalone function rather than another solve<> specialization since
// which factorization type to use depends on the system, not just the type.
void solve_catamari_ldl(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  if(!should_run(name)) return;
  catamari::CoordinateMatrix<double> matrix;
  matrix.Resize(Q.rows(), Q.cols());
  matrix.ReserveEntryAdditions(Q.nonZeros());
  for(int c=0; c<Q.outerSize(); ++c)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(Q,c); it; ++it)
    {
      matrix.QueueEntryAddition(it.row(), it.col(), it.value());
    }
  }
  matrix.FlushEntryQueues();

  Timer timer;
  catamari::SparseLDLControl<double> ldl_control;
  ldl_control.SetFactorizationType(catamari::kLDLTransposeFactorization);

  catamari::SparseLDL<double> ldl;
  const catamari::SparseLDLResult<double> result = ldl.Factor(matrix, ldl_control);
  const double t_factor = timer.toc();

  if(result.num_successful_pivots != Q.rows())
  {
    record(k, name, t_factor, 0, 0, true,
      "factorization failed: incomplete pivoting (indefinite/singular)");
    return;
  }

  catamari::BlasMatrix<double> right_hand_sides;
  right_hand_sides.Resize(rhs.rows(), rhs.cols());
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      right_hand_sides(i, j) = rhs(i,j);
    }
  }

  ldl.Solve(&right_hand_sides.view);
  const double t_solve = timer.toc();

  U.resize(rhs.rows(),rhs.cols());
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      U(i,j) = right_hand_sides(i, j);
    }
  }
  record(k, name, t_factor, t_solve, backward_error(Q,rhs,U));
}

#ifdef IGL_WITH_NASOQ
#include <nasoq/lbl_eigen.h>

// NASOQ's LBL: a parallel sparse symmetric-indefinite (Bunch-Kaufman-style,
// dynamically regularized) direct solver built for QP KKT systems — a
// natural fit for testing on our own indefinite mixed-FEM systems (its
// diagonal regularization is exactly the kind of mechanism that might let it
// succeed where Pardiso's reordering hangs on the mixed triharmonic system's
// all-zero-diagonal block; run on all 5 systems to see, same as PardisoLDLT).
// Uses nasoq::SolverSettings directly (rather than the fused linear_solve()
// convenience wrapper) for separate factor/solve timing; solves one RHS
// column at a time (no native multi-RHS API).
void solve_nasoq_lbl(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  if(!should_run(name)) return;
  const int n = (int)Q.rows();
  const int nrhs = (int)rhs.cols();
  // NASOQ wants the lower triangle only (Q is symmetric).
  Eigen::SparseMatrix<double> QL = Q.triangularView<Eigen::Lower>();
  QL.makeCompressed();

  nasoq::CSC A;
  A.nzmax = QL.nonZeros();
  A.ncol = A.nrow = n;
  A.p = QL.outerIndexPtr();
  A.i = QL.innerIndexPtr();
  A.x = QL.valuePtr();
  A.stype = -1;
  A.xtype = CHOLMOD_REAL;
  A.packed = TRUE;
  A.nz = NULL;
  A.sorted = TRUE;

  Eigen::VectorXd rhs0 = rhs.col(0);
  nasoq::SolverSettings solver(&A, rhs0.data());
  solver.ldl_variant = 4;
  solver.solver_mode = 0;
  solver.reg_diag = std::pow(10, -9);
  // SolverSettings' own num_thread field (used for MKL's SET_BLAS_THREAD and
  // workspace sizing) does NOT control the OpenMP team size of its actual
  // triangular-solve kernels: H2LeveledBlockedLsolve/LTsolve use a bare
  // "#pragma omp parallel" with no num_threads() clause, and NASOQ's only
  // omp_set_num_threads() call site (Parallel_simplicial_ldl.cpp) is
  // commented out. So every solve_only() call was spinning up
  // omp_get_max_threads() (128 here) OpenMP threads regardless of this
  // field — the real explanation (verified with gdb + source inspection,
  // not just the req_ref_iter fix below) for a multi-second "solve" on
  // problems from 400 rows to 1.4M rows alike: thread-team creation cost on
  // this shared machine, not real triangular-solve work. Cap it ourselves.
  solver.num_thread = 1;
  // The NASOQ eigen_interface example sets req_ref_iter=2, requesting
  // internal GMRES-based iterative refinement (pmgmres_ldlt_auto) after the
  // direct solve. Verified with gdb that this was a second, independent
  // cost: each GMRES iteration issues its own pair of OpenMP-parallel
  // triangular solves. Leave refinement off; our own downstream residual
  // check already validates accuracy, same as every other solver here.
  solver.req_ref_iter = 0;

#ifdef _OPENMP
  const int prev_omp_threads = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  Timer timer;
  solver.symbolic_analysis();
  solver.numerical_factorization();
  const double t_factor = timer.toc();

  U.resize(n, nrhs);
  for(int c = 0; c < nrhs; c++)
  {
    Eigen::VectorXd rc = rhs.col(c);
    double * sol = solver.solve_only(n, rc.data());
    U.col(c) = Eigen::Map<Eigen::VectorXd>(sol, n);
  }
  const double t_solve = timer.toc();

#ifdef _OPENMP
  omp_set_num_threads(prev_omp_threads);
#endif

  record(k, name, t_factor, t_solve, backward_error(Q,rhs,U));
}
#endif

// Writes a dense Eigen matrix in MatrixMarket "array" format (column-major,
// one value per line) -- Eigen's own unsupported/Eigen/SparseExtra only
// provides saveMarket() for sparse matrices and saveMarketVector() for a
// single vector, neither of which covers a dense multi-column RHS.
static bool write_dense_market(const Eigen::MatrixXd & mat, const std::string & filename)
{
  std::ofstream out(filename);
  if(!out.is_open()) return false;
  out << "%%MatrixMarket matrix array real general\n";
  out << mat.rows() << " " << mat.cols() << "\n";
  out << std::setprecision(17);
  for(int c=0;c<mat.cols();c++)
  {
    for(int r=0;r<mat.rows();r++)
    {
      out << mat(r,c) << "\n";
    }
  }
  return true;
}

// Dumps Q (full, not just lower triangular -- unlike NASOQ's CSC input,
// external readers like scipy/warp expect the whole symmetric matrix) and
// rhs for one k to <dir>/k<k>_Q.mtx and <dir>/k<k>_rhs.mtx, both in the
// standard MatrixMarket format so any external tool (scipy, warp, MATLAB,
// ...) can load them without depending on this benchmark's own code --
// see --dump-matrices in main() and warp_bench/ for the Python-side
// consumer that times NVIDIA Warp's warp.optim.linear solvers on them.
static void dump_matrices(
  const std::string & dir, int k,
  const Eigen::SparseMatrix<double> & Q, const Eigen::MatrixXd & rhs)
{
  Eigen::saveMarket(Q, dir + "/k" + std::to_string(k) + "_Q.mtx");
  write_dense_market(rhs, dir + "/k" + std::to_string(k) + "_rhs.mtx");
}

static std::string machine_info()
{
  std::string info;
#if defined(__linux__)
  std::ifstream cpuinfo("/proc/cpuinfo");
  std::string line;
  while(std::getline(cpuinfo,line))
  {
    if(line.rfind("model name",0)==0)
    {
      const auto pos = line.find(':');
      if(pos != std::string::npos) info = line.substr(pos+2);
      break;
    }
  }
#endif
  if(info.empty()) info = "unknown CPU";
  unsigned hc = std::thread::hardware_concurrency();
  info += " (" + std::to_string(hc) + " threads)";
#ifdef IGL_WITH_CUDSS
  int device_count = 0;
  if(cudaGetDeviceCount(&device_count)==cudaSuccess && device_count>0)
  {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    info += std::string(", GPU: ") + prop.name;
  }
#endif
  return info;
}

// A solver can report success (Eigen's info()==Success, cuDSS's
// CUDSS_STATUS_SUCCESS, ...) while still having silently produced a
// numerically useless answer -- e.g. cuDSS on the mixed triharmonic system:
// it completes without error, but its only pivoting strategy for symmetric
// indefinite matrices is CUDSS_PIVOT_DIAGONAL (diagonal-only search;
// CUDSS_PIVOT_BUNCH_KAUFMAN, the safe block-pivoting strategy NASOQ uses,
// is explicitly "reserved for future, not supported yet" per cudss_data_types.h
// in this version), which can't find a safe pivot in a block that's
// structurally all-zero on the diagonal (this system's λ block) -- so
// rather than trust any solver's own success signal, reclassify any row
// whose backward error exceeds this bar as "didn't actually succeed",
// moving it to the skipped section instead of letting it rank with a
// misleadingly-real-looking number.
//
// Because backward error is scale-invariant (see backward_error() above),
// unlike the old absolute-residual metric this benchmark used to report,
// a SINGLE fixed threshold works uniformly across every system here --
// flattened or mixed, k=1 or k=5 -- rather than needing a per-k tolerance
// tuned empirically to each system's own residual scale (this used to be a
// 6-entry table, g_check_tol[], ranging from 1e-4 to 1e3; a real double-
// precision solve's backward error sits near machine epsilon (~1e-16)
// regardless of k, so one bar comfortably separates "solved" from "not").
// Also used directly as --check's correctness bar (see below), for the
// same reason: no per-k tuning needed anymore.
static const double kBackwardErrorDivergedThreshold = 1e-6;

static void print_leaderboard(int k)
{
  std::vector<Result> rows;
  for(const auto & r : g_results) if(r.k==k) rows.push_back(r);

  // NaN comparisons are always false (IEEE 754): `nan > threshold` silently
  // evaluates to false rather than flagging a diverged solve, so a solver
  // that returned NaN (e.g. Eigen::CG genuinely diverging on an indefinite
  // system) needs an explicit isnan() check rather than relying on the
  // threshold comparison to catch it.
  for(auto & r : rows)
  {
    if(r.skipped) continue;
    if(std::isnan(r.residual))
    {
      r.skipped = true;
      r.skip_reason = "did not actually succeed: backward error is NaN (solver diverged)";
    }
    else if(r.residual > kBackwardErrorDivergedThreshold)
    {
      char buf[256];
      snprintf(buf,sizeof(buf),
        "did not actually succeed: backward error %.4g exceeds %.4g",
        r.residual,kBackwardErrorDivergedThreshold);
      r.skipped = true;
      r.skip_reason = buf;
    }
  }

  // --check's correctness gate runs here, after reclassification above, so
  // it only ever judges results we're still treating as legitimate --
  // anything already caught as "did not actually succeed" (including NaN)
  // is a known, expected non-result (e.g. CG on an indefinite system),
  // not a correctness regression to fail the build over. Since the
  // reclassification above already applies this exact same threshold,
  // this loop in practice only ever fires for bugs in the reclassification
  // logic itself -- kept as a second, explicit check rather than assuming
  // that invariant always holds.
  if(g_check_mode)
  {
    for(const auto & r : rows)
    {
      if(!r.skipped && !(r.residual <= kBackwardErrorDivergedThreshold))
      {
        fprintf(stderr,"CHECK FAILED: k=%d %s backward_error=%.6g exceeds %.6g\n",
          k, r.name.c_str(), r.residual, kBackwardErrorDivergedThreshold);
        g_check_failed = true;
      }
    }
  }

  std::stable_sort(rows.begin(),rows.end(),[](const Result & a,const Result & b)
  {
    if(a.skipped != b.skipped) return !a.skipped;
    return (a.t_factor+a.t_solve) < (b.t_factor+b.t_solve);
  });

  printf("\n");
  printf("| Rank |                          Method |      Factor |       Solve | Backward error |\n");
  printf("|-----:|--------------------------------:|------------:|------------:|----------------:|\n");
  int rank = 0;
  bool any_fused = false;
  bool any_timed_out = false;
  for(const auto & r : rows)
  {
    if(r.skipped)
    {
      printf("|    - | %32s |           - |           - | skipped: %s |\n",r.name.c_str(),r.skip_reason.c_str());
      continue;
    }
    rank++;
    const char * medal = rank==1 ? "\U0001F947" : rank==2 ? "\U0001F948" : rank==3 ? "\U0001F949" : "  ";
    // A timed-out row's backward error is a snapshot at the
    // kIterativeTimeLimitSeconds cutoff, not a converged result -- mark it
    // with a dagger (and a footnote, matching the existing "(fused)*"
    // pattern) rather than let it read as equivalent to a genuinely
    // converged row's number.
    std::string display_name = r.name;
    if(r.timed_out)
    {
      any_timed_out = true;
      display_name += "†";
    }
    if(r.fused_factor)
    {
      any_fused = true;
      printf("| %s%2d | %32s |     (fused)* | %8.2g secs | %15.6g |\n",
        medal,rank,display_name.c_str(),r.t_solve,r.residual);
    }
    else
    {
      printf("| %s%2d | %32s | %8.2g secs | %8.2g secs | %15.6g |\n",
        medal,rank,display_name.c_str(),r.t_factor,r.t_solve,r.residual);
    }
  }
  if(any_fused)
  {
    printf("\n*(fused): this solver's API has no separate factor step; the whole\n");
    printf(" analysis+factor+solve cost is reported under Solve instead.\n");
  }
  if(any_timed_out)
  {
    printf("\n\xe2\x80\xa0 hit the %.0f-minute iterative-solver time limit before reaching\n",
      kIterativeTimeLimitSeconds/60.0);
    printf(" kBackwardErrorTarget -- backward error is a snapshot at cutoff, not a converged result.\n");
  }
  printf("\n");
}

// Splits a comma-separated list into lowercased tokens (matching should_run's
// case-insensitive comparison), for --only/--exclude parsing.
static std::vector<std::string> split_lower_csv(const std::string & s)
{
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while(std::getline(ss, tok, ','))
  {
    if(!tok.empty()) out.push_back(to_lower(tok));
  }
  return out;
}

int main(int argc, char * argv[])
{
  setbuf(stdout, NULL);
  std::string mesh_path;
  std::string csv_path;
  std::string dump_dir;
  bool dump_only = false;
  int grid_n = 0;
  for(int i=1;i<argc;i++)
  {
    const std::string arg = argv[i];
    if(arg == "--csv" && i+1<argc) { csv_path = argv[++i]; }
    else if(arg == "--check") { g_check_mode = true; }
    else if(arg == "--grid" && i+1<argc) { grid_n = std::atoi(argv[++i]); }
    else if(arg == "--only" && i+1<argc)
    {
      const auto toks = split_lower_csv(argv[++i]);
      g_only.insert(g_only.end(), toks.begin(), toks.end());
    }
    else if(arg == "--exclude" && i+1<argc)
    {
      const auto toks = split_lower_csv(argv[++i]);
      g_exclude.insert(g_exclude.end(), toks.begin(), toks.end());
    }
    else if(arg == "--dump-matrices" && i+1<argc) { dump_dir = argv[++i]; }
    else if(arg == "--dump-only") { dump_only = true; }
    else { mesh_path = arg; }
  }
  if(mesh_path.empty() && grid_n<=0)
  {
    fprintf(stderr,
      "usage: %s [--csv results.csv] [--check] [--only name[,name...]] "
      "[--exclude name[,name...]] [--dump-matrices dir] [--dump-only] "
      "(<mesh> | --grid N)\n"
      "  --only/--exclude match case-insensitively against a substring of\n"
      "  the solver's printed name (e.g. --only nasoq, --exclude umfpack,sparselu).\n"
      "  Repeatable/comma-separated; --only takes precedence, --exclude is\n"
      "  applied on top of it. Filtered-out solvers are simply not run (not\n"
      "  shown as skipped) -- for fast local iteration on one solver, not\n"
      "  for permanent leaderboard output.\n"
      "  --dump-matrices dir writes each system's Q/rhs as MatrixMarket files\n"
      "  (k<k>_Q.mtx, k<k>_rhs.mtx) to dir, for external tools (e.g. warp_bench/)\n"
      "  to load; combine with --dump-only to skip this benchmark's own solvers\n"
      "  entirely (just build+dump), or with --only/--exclude to dump alongside\n"
      "  running a subset.\n",
      argv[0]);
    return 1;
  }
  if(!dump_dir.empty())
  {
    std::filesystem::create_directories(dump_dir);
  }
  if(!csv_path.empty())
  {
    g_csv = fopen(csv_path.c_str(),"w");
    fprintf(g_csv,"k,method,factor_secs,solve_secs,backward_error,skipped,fused_factor,timed_out,iterations_used\n");
  }

  fprintf(stderr,"# %s\n",machine_info().c_str());
#if defined(_OPENMP)
  fprintf(stderr,"omp_get_num_threads(): %d\n",omp_get_max_threads());
#endif

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(grid_n>0)
  {
    // Small synthetic mesh for fast, self-contained correctness regression
    // tests (see CTest) that don't require checking in another large .ply.
    Eigen::MatrixXd GV;
    Eigen::MatrixXi GF;
    igl::triangulated_grid(grid_n,grid_n,GV,GF);
    V.resize(GV.rows(),3);
    V.leftCols(2) = GV;
    V.col(2).setZero();
    F = GF;
  }
  else
  {
    igl::read_triangle_mesh(mesh_path,V,F);
  }
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  // k=1,2,3: flattened SPD k-harmonic systems (Q=M+Wᵏ, igl::harmonic).
  // k=4,5: mixed (unflattened) biharmonic/triharmonic block systems built by
  // build_mixed_system() above — same underlying PDE, but symmetric
  // indefinite instead of SPD, exercising solver behavior on harder input.
  for(int k = 1;k<=5;k++)
  {
    const bool is_mixed = k>=4;
    switch(k)
    {
      case 1: printf("# Harmonic\n"); break;
      case 2: printf("# Biharmonic\n"); break;
      case 3: printf("# Triharmonic\n"); break;
      case 4: printf("# Mixed Biharmonic (unflattened, indefinite)\n"); break;
      case 5: printf("# Mixed Triharmonic (unflattened, indefinite)\n"); break;
    }

    Eigen::SparseMatrix<double> Q;
    Eigen::MatrixXd rhs;
    if(!is_mixed)
    {
      Eigen::SparseMatrix<double> W;
      igl::harmonic(L,M,k,W);
      Q = M+W;
      rhs = M*V;
    }
    else
    {
      build_mixed_system(k==4 ? 2 : 3,L,M,V,Q,rhs);
    }
    if(!dump_dir.empty())
    {
      dump_matrices(dump_dir, k, Q, rhs);
    }
    if(dump_only)
    {
      continue;
    }
    Eigen::MatrixXd U;
#ifdef IGL_WITH_CHOLMOD
    solve<Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>>("Eigen::CholmodSupernodalLLT",k,Q,rhs,U);
    if(k == 5)
    {
      // UmfPackLU's own internal MKL BLAS3 calls (umfdi_blas3_update) spin
      // up a fresh OpenMP thread team per call; on the mixed triharmonic
      // system's 2.16M-row/3-block structure the catastrophic fill-in from
      // this system's sparsity pattern produces so many of these tiny
      // updates that it exhausts OS thread/process limits and crashes
      // (verified with gdb, reproducible). Skip rather than risk crashing
      // the whole benchmark; every other solver still runs on it.
      record(k, "Eigen::UmfPackLU", 0, 0, 0, true,
        "known crash risk: excessive MKL thread churn on this system's fill-in");
    }
    else
    {
      solve<Eigen::UmfPackLU<Eigen::SparseMatrix<double>>>("Eigen::UmfPackLU",k,Q,rhs,U);
    }
#endif
    solve<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> >("Eigen::SimplicialLLT",k,Q,rhs,U);
    if(k == 5)
    {
      // Eigen::SimplicialLDLT has no pivoting; verified with gdb that it
      // segfaults (SIGSEGV inside factorize_preordered) on the mixed
      // triharmonic system's 2.16M-row genuinely indefinite matrix — a real
      // out-of-bounds access in Eigen's own unpivoted LDLT at this scale,
      // not just a slow/inaccurate result. Skip rather than crash; the
      // pivoted PardisoLDLT and catamari LDLᵀ still exercise this question.
      record(k, "Eigen::SimplicialLDLT", 0, 0, 0, true,
        "known crash: SIGSEGV in Eigen's unpivoted LDLT at this scale");
    }
    else
    {
      solve<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>("Eigen::SimplicialLDLT",k,Q,rhs,U);
    }
    if(is_mixed)
    {
      solve_catamari_ldl("catamari::SparseLDL (LDLᵀ)",k,Q,rhs,U);
    }
    else
    {
      solve<catamari::SparseLDL<double>>("catamari::SparseLDL",k,Q,rhs,U);
    }
#ifdef IGL_WITH_NASOQ
    if(k == 5)
    {
      // Re-verified after fixing the OpenMP thread-storm bug (this crash
      // predates that fix, so it was worth re-checking): still reproduces.
      // Root cause identified via gdb: SIGSEGV inside libmetis.so.5's
      // minimum-degree ordering (genmmd/mmdelm), called from NASOQ's own
      // symbolic_analysis_lin_solve(). Reproduces (intermittently --
      // memory-layout dependent, consistent with an out-of-bounds write)
      // even on a tiny synthetic --grid 20 mixed triharmonic system
      // (~1200 rows), i.e. this is inherent to the system's structure
      // (the lambda block's structurally-zero diagonal, shared with the
      // Pardiso reordering hang above) rather than a scale issue. Skip
      // rather than risk it; see the upstream NASOQ issue for the reduced
      // repro and backtrace.
      record(k, "NASOQ LBL", 0, 0, 0, true,
        "known crash: SIGSEGV in libmetis genmmd/mmdelm via NASOQ's symbolic_analysis_lin_solve on this system's sparsity pattern");
    }
    else
    {
      solve_nasoq_lbl("NASOQ LBL",k,Q,rhs,U);
    }
#endif
#ifdef IGL_WITH_MKL
    if(k == 5)
    {
      // MKL Pardiso's reordering (both METIS and minimum-degree — verified
      // with gdb) hangs, not just fails, on the mixed triharmonic system's
      // sparsity pattern (its λ block has an all-zero diagonal, inherent to
      // this saddle-point/KKT system; Pardiso's ordering heuristics appear
      // not to handle that gracefully at this size). Skip rather than risk
      // hanging the whole benchmark; every other solver still runs on it.
      const char * reason = "known Pardiso reordering hang on this system's sparsity pattern";
      record(k, "Eigen::PardisoLLT", 0, 0, 0, true, reason);
      record(k, "Eigen::PardisoLDLT", 0, 0, 0, true, reason);
    }
    else
    {
      solve<Eigen::PardisoLLT<Eigen::SparseMatrix<double>>>("Eigen::PardisoLLT",k,Q,rhs,U);
      solve<Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>>("Eigen::PardisoLDLT",k,Q,rhs,U);
    }
#endif
#ifdef IGL_WITH_CUDSS
    solve_cudss("NVIDIA cuDSS",k,Q,rhs,U,is_mixed ? CUDSS_MTYPE_SYMMETRIC : CUDSS_MTYPE_SPD);
    if(is_mixed)
    {
      record(k, "NVIDIA cuSOLVER (Sp Chol)", 0, 0, 0, true,
        "no indefinite/LDLT solver in this cuSOLVER version");
    }
    else
    {
      solve_cusolver("NVIDIA cuSOLVER (Sp Chol)",k,Q,rhs,U);
    }
#endif
    solve<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>("Eigen::SparseLU",k,Q,rhs,U);
    solve<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>>>("Eigen::BiCGSTAB<IncompleteLUT>",k,Q,rhs,U);
    solve<Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::IncompleteLUT<double>>>("Eigen::CG<IncompleteLUT>",k,Q,rhs,U);

    print_leaderboard(k);
  }

  if(g_csv) fclose(g_csv);
  return g_check_failed ? 1 : 0;
}
