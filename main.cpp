
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
#ifdef IGL_WITH_CUDSS
#include <cuda_runtime_api.h>
#include <cudss.h>
#include <cusolverSp.h>
#include <cusparse.h>
#endif
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
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
};

static std::vector<Result> g_results;
static FILE * g_csv = nullptr;
static bool g_check_mode = false;
// Looser tolerance for higher k: the k-harmonic systems get increasingly
// ill-conditioned/dense (see README), so residuals legitimately grow with k.
static double g_check_tol[4] = {0, 1e-4, 1e-1, 1e3};
static bool g_check_failed = false;

static void record(
  int k,
  const std::string & name,
  double t_factor,
  double t_solve,
  double residual,
  bool skipped = false,
  const std::string & skip_reason = "")
{
  g_results.push_back({k, name, t_factor, t_solve, residual, skipped, skip_reason});
  if(g_csv)
  {
    fprintf(g_csv,"%d,%s,%.9g,%.9g,%.9g,%d\n",
      k, name.c_str(), t_factor, t_solve, residual, skipped?1:0);
  }
  if(g_check_mode && !skipped && !(residual <= g_check_tol[k]))
  {
    fprintf(stderr,"CHECK FAILED: k=%d %s residual=%.6g exceeds tolerance %.6g\n",
      k, name.c_str(), residual, g_check_tol[k]);
    g_check_failed = true;
  }
}

// Iterative solvers (BiCGSTAB/ConjugateGradient) default to Eigen's built-in
// cap of twice the system size, which is effectively unbounded for a 720K-row
// mesh: on the k=3 triharmonic system (documented above as badly scaled) they
// can burn tens of minutes grinding through non-convergent iterations without
// ever tripping any other limit. Cap iterations so a hard system shows up as
// an honestly-reported large residual within a bounded time, not a hang.
// SFINAE dispatch: direct solvers (LLT/LDLT/LU/...) have no setMaxIterations,
// so the fallback (long) overload is selected for them and does nothing.
static const int kMaxIterativeIterations = 200;
template <typename Factor>
auto cap_iterations(Factor & factor, int) -> decltype(factor.setMaxIterations(0), void())
{
  factor.setMaxIterations(kMaxIterativeIterations);
}
template <typename Factor>
void cap_iterations(Factor &, long) {}

template <typename Factor>
void solve(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  Timer timer;
  Factor factor;
  cap_iterations(factor, 0);
  factor.compute(Q);
  const double t_factor = timer.toc();
  U = factor.solve(rhs);
  const double t_solve = timer.toc();
  record(k, name, t_factor, t_solve, (rhs-Q*U).array().abs().maxCoeff());
}

template <>
void solve<catamari::SparseLDL<double>>(
  const std::string & name,
  int k,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
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
  record(k, name, t_factor, t_solve, (rhs-Q*U).array().abs().maxCoeff());
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
  Eigen::MatrixXd & U)
{
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
  CUDSS_CHECK(cudssDataCreate(handle,&data));

  cudssMatrix_t A, B, X;
  CUDSS_CHECK(cudssMatrixCreateCsr(&A,n,n,A_csr.nnz,A_csr.row_ptr,nullptr,A_csr.col_idx,A_csr.val,
    CUDSS_R_32I,CUDSS_R_32I,CUDSS_R_64F,CUDSS_MTYPE_SPD,CUDSS_MVIEW_FULL,CUDSS_BASE_ZERO));
  CUDSS_CHECK(cudssMatrixCreateDn(&B,n,nrhs,n,d_b,CUDSS_R_64F,CUDSS_LAYOUT_COL_MAJOR));
  CUDSS_CHECK(cudssMatrixCreateDn(&X,n,nrhs,n,d_x,CUDSS_R_64F,CUDSS_LAYOUT_COL_MAJOR));

  Timer timer;
  CUDSS_CHECK(cudssExecute(handle,CUDSS_PHASE_REORDERING,config,data,A,X,B));
  CUDSS_CHECK(cudssExecute(handle,CUDSS_PHASE_SYMBOLIC_FACTORIZATION,config,data,A,X,B));
  CUDSS_CHECK(cudssExecute(handle,CUDSS_PHASE_FACTORIZATION,config,data,A,X,B));
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_factor = timer.toc();

  CUDSS_CHECK(cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,X,B));
  CUDA_CHECK(cudaDeviceSynchronize());
  const double t_solve = timer.toc();

  U.resize(n,nrhs);
  CUDA_CHECK(cudaMemcpy(U.data(),d_x,n*nrhs*sizeof(double),cudaMemcpyDeviceToHost));

  CUDSS_CHECK(cudssMatrixDestroy(A));
  CUDSS_CHECK(cudssMatrixDestroy(B));
  CUDSS_CHECK(cudssMatrixDestroy(X));
  CUDSS_CHECK(cudssDataDestroy(handle,data));
  CUDSS_CHECK(cudssConfigDestroy(config));
  CUDSS_CHECK(cudssDestroy(handle));
  cudaFree(d_b);
  cudaFree(d_x);

  record(k, name, t_factor, t_solve, (rhs-Q*U).array().abs().maxCoeff());
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

  record(k, name, 0, t_solve, (rhs-Q*U).array().abs().maxCoeff());
}

#endif

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

static void print_leaderboard(int k)
{
  std::vector<Result> rows;
  for(const auto & r : g_results) if(r.k==k) rows.push_back(r);
  std::stable_sort(rows.begin(),rows.end(),[](const Result & a,const Result & b)
  {
    if(a.skipped != b.skipped) return !a.skipped;
    return (a.t_factor+a.t_solve) < (b.t_factor+b.t_solve);
  });

  printf("\n");
  printf("| Rank |                          Method |      Factor |       Solve |     L∞ norm |\n");
  printf("|-----:|--------------------------------:|------------:|------------:|------------:|\n");
  int rank = 0;
  for(const auto & r : rows)
  {
    if(r.skipped)
    {
      printf("|    - | %32s |           - |           - | skipped: %s |\n",r.name.c_str(),r.skip_reason.c_str());
      continue;
    }
    rank++;
    const char * medal = rank==1 ? "\U0001F947" : rank==2 ? "\U0001F948" : rank==3 ? "\U0001F949" : "  ";
    printf("| %s%2d | %32s | %8.2g secs | %8.2g secs | %11.6g |\n",
      medal,rank,r.name.c_str(),r.t_factor,r.t_solve,r.residual);
  }
  printf("\n");
}

int main(int argc, char * argv[])
{
  setbuf(stdout, NULL);
  std::string mesh_path;
  std::string csv_path;
  int grid_n = 0;
  for(int i=1;i<argc;i++)
  {
    const std::string arg = argv[i];
    if(arg == "--csv" && i+1<argc) { csv_path = argv[++i]; }
    else if(arg == "--check") { g_check_mode = true; }
    else if(arg == "--grid" && i+1<argc) { grid_n = std::atoi(argv[++i]); }
    else { mesh_path = arg; }
  }
  if(mesh_path.empty() && grid_n<=0)
  {
    fprintf(stderr,"usage: %s [--csv results.csv] [--check] (<mesh> | --grid N)\n",argv[0]);
    return 1;
  }
  if(!csv_path.empty())
  {
    g_csv = fopen(csv_path.c_str(),"w");
    fprintf(g_csv,"k,method,factor_secs,solve_secs,linf_residual,skipped\n");
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

  for(int k = 1;k<=3;k++)
  {
    switch(k)
    {
      case 1: printf("# Harmonic\n"); break;
      case 2: printf("# Biharmonic\n"); break;
      case 3: printf("# Triharmonic\n"); break;
    }

    Eigen::SparseMatrix<double> W;
    igl::harmonic(L,M,k,W);
    Eigen::SparseMatrix<double> Q;
    Q = M+W;
    Eigen::MatrixXd rhs = M*V;

    Eigen::MatrixXd U;
#ifdef IGL_WITH_CHOLMOD
    solve<Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>>("Eigen::CholmodSupernodalLLT",k,Q,rhs,U);
    solve<Eigen::UmfPackLU<Eigen::SparseMatrix<double>>>("Eigen::UmfPackLU",k,Q,rhs,U);
#endif
    solve<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> >("Eigen::SimplicialLLT",k,Q,rhs,U);
    solve<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>("Eigen::SimplicialLDLT",k,Q,rhs,U);
    solve<catamari::SparseLDL<double>>("catamari::SparseLDL",k,Q,rhs,U);
#ifdef IGL_WITH_MKL
    solve<Eigen::PardisoLLT<Eigen::SparseMatrix<double>>>("Eigen::PardisoLLT",k,Q,rhs,U);
#endif
#ifdef IGL_WITH_CUDSS
    solve_cudss("NVIDIA cuDSS",k,Q,rhs,U);
    solve_cusolver("NVIDIA cuSOLVER (Sp Chol)",k,Q,rhs,U);
#endif
    solve<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>("Eigen::SparseLU",k,Q,rhs,U);
    solve<Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>>>("Eigen::BiCGSTAB<IncompleteLUT>",k,Q,rhs,U);
    solve<Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::IncompleteLUT<double>>>("Eigen::CG<IncompleteLUT>",k,Q,rhs,U);

    print_leaderboard(k);
  }

  if(g_csv) fclose(g_csv);
  return g_check_failed ? 1 : 0;
}
