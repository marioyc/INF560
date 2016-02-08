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

  return 0;
}

________________________________________________________________________________

//! @internal backward substitution for Crout factorized matrix (L/U)
int Backward (
        Vector<double,int>& x,
        const MatrixDense<double,int>& U,
        const Vector<double,int>& rhs,
        const bool diag_one ) {

  return 0;
}

________________________________________________________________________________

//! @internal solve with crout (LU) factorization
int SolveLU (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) {

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
