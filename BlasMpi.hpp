/*!
*  @file BlasMpi.hpp
*  @brief header of class BlasMpi
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_BLASMPI_HPP_
#define GUARD_BLASMPI_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace BlasMpi
namespace BlasMpi {

//! @brief compute matrix-vector product y:= A *x
//! @param [out] y_local = local result vector
//! @param [in] A_local = local matrix
//! @param [in] x_local = local vector
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int MatrixVectorProductBandRow (
        Vector<double,int>& y_local,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        MPI_Comm& mpi_comm ) ;

//! @brief compute matrix-vector product y:= A *x
//! @param [out] y_global = y global result vector
//! @param [in] A_local = local matrix
//! @param [in] x_local = local vector
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int MatrixVectorProductBandColumn (
        Vector<double,int>& y_global,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm ) ;

// NOT VALIDATED
//! @brief compute matrix-vector product y:= A *x
//! @param [out] y_local = local result vector
//! @param [in] A_local = local matrix
//! @param [in] x_local = local vector
//! @param [in] root = root processor
//! @param [in] mpi_comm_rows = grid rows communicator
//! @param [in] mpi_comm_columns = grid columns communicator
//! @return error code
int MatrixVectorProductBlock (
        Vector<double,int>& y_local,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) ;

// NOT VALIDATED
//! @brief compute matrix-matrix product C := A * B
//! @param [out] A = local result vector
//! @param [in] B = local matrix
//! @param [in] C = local vector
//! @param [in] root = root processor
//! @param [in] mpi_comm_rows = grid rows communicator
//! @param [in] mpi_comm_columns = grid columns communicator
//! @return error code
int MatrixMatrixProductBlock (
        MatrixDense<double,int>& C_local,
        const MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& B_local,
        const int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) ;

} // namespace BlasMpi {

#endif // GUARD_BLASMPI_HPP_
