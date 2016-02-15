/*!
*  @file Schur.hpp
*  @brief manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_SCHUR_HPP_
#define GUARD_SCHUR_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace Schur
namespace Schur {


// -----------------------------------------------------------------------------
// -- topology distribution: processors
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// -- split local matrix
// -----------------------------------------------------------------------------


//! @brief block splitting of a matrix
//! @param [in] A = matrix of size (ni+np, ni+np)
//! @param [in,out] Aii = block matrix associated to the internal nodes
//! @param [in,out] Aip = block matrix associated to the internal nodes
//! @param [in,out] Api =
//! @param [in,out] App =
//! @param [in] list_node_i = list of internal nodes
//! @param [in] numb_node_i = number of internal nodes
//! @param [in] list_node_p = list of interface (boundary) nodes
//! @param [in] numb_node_p = number of interface (boundary) nodes
int SplitMatrixToBlock (
        MatrixDense<double,int>& Aii,
        MatrixDense<double,int>& Aip,
        MatrixDense<double,int>& Api,
        MatrixDense<double,int>& App,
        const MatrixDense<double,int>& A,
        const int* list_node_i,
        const int numb_node_i,
        const int* list_node_p,
        const int numb_node_p ) ;

void CheckTranspose(
        const MatrixDense<double,int>& Kip,
        const MatrixDense<double,int>& Kpi );

void LocalMatrixToGlobalPositions(
        const MatrixDense<double,int>& A_local,
        MatrixDense<double,int>& A_global,
        int* l2g );

void ReconstructK(
        const MatrixDense<double, int>& K_local,
        const int* l2g,
        const int numb_global_node,
        const std::string gmatrix_filename,
        const MPI_Comm& mpi_comm );

void SolveSystem (
		Vector<double, int>& x_local,
		const MatrixDense<double, int>& K_local,
		const Vector<double, int>& b_local,
        const int numb_node_i,
        const int* list_node_i,
        const int* l2i,
        int numb_node_p,
        const int* list_node_p,
        const int* l2p,
        const int numb_global_node,
        const int numb_l2g,
        const int* l2g,
        const MPI_Comm& mpi_comm );


} // namespace Schur {


#endif // GUARD_SCHUR_HPP_
