/*!
*  @file DirectSolver.hpp
*  @brief manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_DIRECTSOLVER_HPP_
#define GUARD_DIRECTSOLVER_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages

//! @namespace Factor
//! @brief box of factorization
namespace Factor {


// -----------------------------------------------------------------------------
// -- Factorization: LU
// -----------------------------------------------------------------------------

//! @brief LU factorization A = L.U
//! @param [in,out] LU = LU factorization
//!                 where L=lower(LU) with diag(L)=1, U=upper(U)
//! @param [in] A = matrix
//! @return error code
int LU (
        MatrixDense<double,int>& LU,
        const MatrixDense<double,int>& A ) ;

// -----------------------------------------------------------------------------
// -- Factorization: LDLt
// -----------------------------------------------------------------------------

//! @brief crout (LDLt) factorization
//! @param [in,out] LDLt = LDL' factorization
//!                 where L=lower(LDLt), D=diag(LDLt), Lt=upper(LDLt)
//! @param [in] A = matrix
//! @return error code
int LDLt (
        MatrixDense<double,int>& LDLt,
        const MatrixDense<double,int>& A ) ;

// -----------------------------------------------------------------------------
// -- Factorization: Cholesky
// -----------------------------------------------------------------------------

//! @brief Cholesky factorization A = L.Lt
//! @param [in,out] LLt = Cholesky factorization
//!                 where L=lower(C), Lt=upper(C)
//! @param [in] A = matrix
//! @param [in] where_Lt = lower (lower triangular matrix),
//!             upper (triangular matrix) and full all data
//! @remarks L = lower(LLt) and L' = upper(LLt) => Llt = L + L' - diag(L)
//! @return error code
int Cholesky (
        MatrixDense<double,int>& LLt,
        const MatrixDense<double,int>& A,
        const char* where_L = "full" ) ;

} // namespace Factor {


//! @namespace DirectSolver
namespace DirectSolver {


// -----------------------------------------------------------------------------
// -- Solver: Forward and Backward
// -----------------------------------------------------------------------------

//! @brief forward substitution for Crout factorized matrix (L/U)
//! @param[in,out] x = solution
//! @param [in] L = matrix
//! @param [in] rhs = right hand side
//! @param [in] diag_one = if true consider diag = 1, else diag = L(i,i)
//! @return error code
int Forward (
        Vector<double,int>& x,
        const MatrixDense<double,int>& L,
        const Vector<double,int>& rhs,
        const bool diag_one = false ) ;


//! @brief backward substitution for Crout factorized matrix (L/U)
//! @param[in,out] x = solution
//! @param [in] U = matrix
//! @param [in] rhs = right hand side
//! @param [in] diag_one = if true consider diag = 1, else diag = U(i,i)
//! @return error code
int Backward (
        Vector<double,int>& x,
        const MatrixDense<double,int>& U,
        const Vector<double,int>& rhs,
        const bool diag_one = false ) ;


//! @brief solve with crout (LU) factorization
//! @param [out] x = solution
//! @param [in] A = matrix
//! @param [in] rhs = right hand side
//! @param [in] is_free_memory = is free memory ?
//! @return error code
int SolveLU (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) ;


//! @brief solve with crout (Cholesky) factorization
//! @param [out] x = solution
//! @param [in] A = matrix
//! @param [in] rhs = right hand side
//! @param [in] is_free_memory = is free memory ?
//! @return error code
int SolveCholesky (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) ;

//! @brief solve with crout (LDLt) factorization
//! @param [out] x = solution
//! @param [in] A = matrix
//! @param [in] rhs = right hand side
//! @param [in] is_free_memory = is free memory ?
//! @return error code
int SolveLDLt (
        Vector<double,int>& x,
        MatrixDense<double,int>& A,
        const Vector<double,int>& rhs ) ;

} // namespace DirectSolver {


#endif // GUARD_DIRECTSOLVER_HPP_
