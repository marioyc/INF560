/*!
*  @file Jacobi.hpp
*  @brief header of class Jacobi
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Wed Dec 09 16:36:50 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_SOLVER_HPP_
#define GUARD_SOLVER_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace Solver
namespace Solver {

//! @brief solve A * x = b by Jacobi sequential algorithm
//! @param [in,out] x = solution
//! @param [in,out] numb_iter = number of iteration
//! @param [in,out] residual = history of residual
//! @param [in,out] temporal = temporal of residual
//! @param [in] A = matrix
//! @param [in] b = right hand side
//! @param [in] max_numb_iter = maximum number of iteration
//! @param [in] residual_threshold = residual threshold
//! @param [in] mpi_comm = MPI communicator
//! @note A_local is supposed to be a square matrix
//! @return error code
int SequentialJacobi (
        Vector<double,int>& x,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A,
        const Vector<double,int>& b,
        const int max_numb_iter,
        const double residual_threshold ) ;

//! @brief solve Ax = b_local by Jacobi synchronous iteration, synchronous comm algorithm
//! @param [in,out] x_local = local solution
//! @param [in,out] numb_iter = iteration number
//! @param [in,out] residual = history of residual
//! @param [in,out] temporal = temporal of residual
//! @param [in] A_local = local matrix
//! @param [in] b_local = local right hand side
//! @param [in] max_numb_iter = maximum iteration number
//! @param [in] residu_threshold = residual threshold
//! @param [in] mpi_comm = MPI communicator
//! @note A_local is a square matrix
//! remarks synchronous communication, synchronous iteration
//! @return error code
int JacobiBandRowSyncISyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residu_threshold,
        MPI_Comm& mpi_comm ) ;

// NOT VALIDATED
//! @brief solve Ax = b_local by Jacobi synchronous iteration, asynchronous comm algorithm
//! @param [in,out] x_local = local solution
//! @param [in,out] numb_iter = iteration number
//! @param [in,out] residual = history of residual
//! @param [in,out] temporal = temporal of residual
//! @param [in] A_local = local matrix
//! @param [in] b_local = local right hand side
//! @param [in] max_numb_iter = maximum iteration number
//! @param [in] residual_threshold = residual threshold
//! @param [in] mpi_comm = MPI communicator
//! @note A_local is a square matrix
//! remarks synchronous communication, synchronous iteration
//! @return error code
int JacobiBandRowSyncIAsyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residual_threshold,
        MPI_Comm& mpi_comm ) ;

// NOT VALIDATED
//! @brief solve Ax = b_local by Jacobi asynchronous iteration, asynchronous comm algorithm
//! @param [in,out] x_local = local solution
//! @param [in,out] numb_iter = iteration number
//! @param [in,out] residual = history of residual
//! @param [in,out] temporal = temporal of residual
//! @param [in] A_local = local matrix
//! @param [in] b_local = local right hand side
//! @param [in] max_numb_iter = maximum iteration number
//! @param [in] residual_threshold = residual threshold
//! @param [in] mpi_comm = MPI communicator
//! @note A_local is a square matrix
//! remarks synchronous communication, synchronous iteration
//! @return error code
int JacobiBandRowAsyncIAsyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residual_threshold,
        MPI_Comm& mpi_comm ) ;

} // namespace Solver {

#endif // GUARD_SOLVER_HPP_
