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




} // namespace Schur {


#endif // GUARD_SCHUR_HPP_
