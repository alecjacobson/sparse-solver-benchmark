
#include "catamari.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/harmonic.h>
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/get_seconds.h>
#include <igl/matlab_format.h>
#include <Eigen/PardisoSupport>
#include <Eigen/CholmodSupport>
#include <tuple>

template <typename Factor>
void solve(
  const std::string & name,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  printf("| %30s | ",name.c_str());
  tictoc();
  const Factor factor(Q);
  const double t_factor = tictoc();
  printf("%6.2g secs | ",t_factor);
  tictoc();
  U = factor.solve(rhs);
  const double t_solve = tictoc();
  printf("%6.2g secs | ",t_solve);
  printf("%6.6g |\n",(rhs-Q*U).array().abs().maxCoeff());
}

template <>
void solve<catamari::SparseLDL<double>>(
  const std::string & name,
  const Eigen::SparseMatrix<double> & Q,
  const Eigen::MatrixXd & rhs,
  Eigen::MatrixXd & U)
{
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  printf("| %30s | ",name.c_str());
  catamari::CoordinateMatrix<double> matrix;
  matrix.Resize(Q.rows(), Q.cols());
  matrix.ReserveEntryAdditions(Q.nonZeros());
  // Queue updates of entries in the sparse matrix using commands of the form:
  for(int k=0; k<Q.outerSize(); ++k)
  {
    for(Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it)
    {
      matrix.QueueEntryAddition(it.row(), it.col(), it.value());
    }
  }
  matrix.FlushEntryQueues();
    
  tictoc();
  // Fill the options for the factorization.
  catamari::SparseLDLControl<double> ldl_control;
  ldl_control.SetFactorizationType(catamari::kCholeskyFactorization);
    
  // Factor the matrix.
  catamari::SparseLDL<double> ldl;
  const catamari::SparseLDLResult<double> result = ldl.Factor(matrix, ldl_control);
  const double t_factor = tictoc();
  printf("%6.2g secs | ",t_factor);
    
  // copy rhs
  catamari::BlasMatrix<double> right_hand_sides;
  right_hand_sides.Resize(rhs.rows(), rhs.cols());
  // The (i, j) entry of the right-hand side can easily be read or modified, e.g.:
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      right_hand_sides(i, j) = rhs(i,j);
    }
  }

  // Solve a linear system using the factorization.
  tictoc();
  ldl.Solve(&right_hand_sides.view);
  const double t_solve = tictoc();
  printf("%6.2g secs | ",t_solve);

  // copy solution
  U.resize(rhs.rows(),rhs.cols());
  for(int i = 0;i<rhs.rows();i++)
  {
    for(int j = 0;j<rhs.cols();j++)
    {
      U(i,j) = right_hand_sides(i, j);
    }
  }
  printf("%6.6g |\n",(rhs-Q*U).array().abs().maxCoeff());
}


int main(int argc, char * argv[])
{
  setbuf(stdout, NULL);
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //igl::triangulated_grid(2,2,V,F);
  igl::read_triangle_mesh(argv[1],V,F);
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
    printf("\n");
    printf("|                         Method |      Factor |       Solve |     L∞ norm |\n");
    printf("|-------------------------------:|------------:|------------:|------------:|\n");
    solve<Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>>("Eigen::CholmodSupernodalLLT",Q,rhs,U);
    solve<Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> >("Eigen::SimplicialLLT",Q,rhs,U);
    solve<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>("Eigen::SimplicialLDLT",Q,rhs,U);
    solve<catamari::SparseLDL<double>>("catamari::SparseLDL",Q,rhs,U);
    solve<Eigen::PardisoLLT<Eigen::SparseMatrix<double>>>("Eigen::PardisoLLT",Q,rhs,U);
    solve<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>("Eigen::SparseLU",Q,rhs,U);
    printf("\n");
  }
}
