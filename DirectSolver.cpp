/*!
*  @file DirectSolver.cpp
*  @internal manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Fr�d�ric Magoul�s, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages

// project packages
#include "DirectSolver.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace Factor
namespace Factor {

________________________________________________________________________________


// -----------------------------------------------------------------------------
// -- Factorization: LU
// -----------------------------------------------------------------------------

//! @internal LU factorization m = L.U
int LU (
        MatrixDense<double,int>& LU,
        const MatrixDense<double,int>& A ) {

  int n = A.GetNumbRows();
  LU.Allocate(n,n);

  for(int i = 0;i < n;++i){
      for(int j = 0;j < i;++j){
        double sum = 0;

        for(int k = 0;k < j;++k){
          sum += LU(i,k) * LU(k,j);
        }

        LU(i,j) = (A(i,j) - sum) / LU(j,j);
      }

      for(int j = i;j < n;++j){
        double sum = 0;

        for(int k = 0;k < i;++k){
          sum += LU(i,k) * LU(k,j);
        }

        LU(i,j) = A(i,j) - sum;
      }
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Factorization: LDLt
// -----------------------------------------------------------------------------

//! @internal LDLt factorization A = L.D.Lt
int LDLt (
        MatrixDense<double,int>& LDLt,
        const MatrixDense<double,int>& A ) {

  int n = A.GetNumbRows();
  LDLt.Allocate(n,n);

  for(int i = 0;i < n;++i){
    for(int j = 0;j <= i;++j){
      LDLt(i,j) = A(i,j);
    }
  }

  for(int j = 0;j < n;++j){
    for(int k = j + 1;k < n;++k){
      LDLt(k,j) = LDLt(k,j) / LDLt(j,j);
    }

    for(int i = j + 1;i < n;++i){
      LDLt(i,i) = LDLt(i,i) - LDLt(i,j) * LDLt(i,j) * LDLt(j,j);

      for(int k = i + 1;k < n;++k){
        LDLt(k,i) = LDLt(k,i) - LDLt(k,j) * LDLt(i,j) * LDLt(j,j);
      }
    }
  }

  return 0;
}

________________________________________________________________________________


// -----------------------------------------------------------------------------
// -- Factorization: Cholesky
// -----------------------------------------------------------------------------

________________________________________________________________________________

//! @internal Cholesky factorization A = L.Lt
int Cholesky (
        MatrixDense<double,int>& LLt,
        const MatrixDense<double,int>& A,
        const char* where_L  ) {

  return 0;
}

________________________________________________________________________________


} // namespace Factor {

________________________________________________________________________________


//! @namespace DirectSolver
namespace DirectSolver {

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Solver: Forward and Backward
// -----------------------------------------------------------------------------

//! @internal forward substitution for Crout factorized matrix (L/U)
int Forward (
        Vector<double,int>& x,
        const MatrixDense<double,int>& L,
        const Vector<double,int>& rhs,
        const bool diag_one ) {
  int n = L.GetNumbRows();
  x.Allocate(n);

  for(int i = 0;i < n;++i){
    double b = rhs(i);

    for(int j = 0;j < i;++j)
      b -= L(i,j) * x(j);

    x(i) = b;// / L(i,i);
  }

  return 0;
}

________________________________________________________________________________

//! @internal backward substitution for Crout factorized matrix (L/U)
int Backward (
        Vector<double,int>& x,
        const MatrixDense<double,int>& U,
        const Vector<double,int>& rhs,
        const bool diag_one ) {
  int n = U.GetNumbRows();
  x.Allocate(n);

  for(int i = n - 1;i >= 0;--i){
    double b = rhs(i);

    for(int j = i + 1;j < n;++j)
      b -= U(i,j) * x(j);

    x(i) = b / U(i,i);
  }

  return 0;
}

________________________________________________________________________________

//! @internal solve with crout (LU) factorization
int SolveLU (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) {
  MatrixDense<double,int> lu;
  Factor::LU(lu, A);

  Vector<double,int> y;
  Forward(y, lu, rhs, true);

  Backward(x, lu, y, false);

  return 0;
}

________________________________________________________________________________

//! @internal solve with crout (LDLt) factorization
int SolveLDLt (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) {


  return 0;
}

________________________________________________________________________________

//! @internal solve with crout (Cholesky) factorization
int SolveCholesky (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) {

  return 0;
}

________________________________________________________________________________

} // namespace DirectSolver {
