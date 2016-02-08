/*!
*  @file Schur.cpp
*  @internal manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages

// project packages
#include "Schur.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace Schur
namespace Schur {

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- topology distribution: processors
// -----------------------------------------------------------------------------


________________________________________________________________________________

// -----------------------------------------------------------------------------
// --
// -----------------------------------------------------------------------------

//! @internal block splitting of a matrix
int SplitMatrixToBlock (
        MatrixDense<double,int>& Aii,
        MatrixDense<double,int>& Aip,
        MatrixDense<double,int>& Api,
        MatrixDense<double,int>& App,
        const MatrixDense<double,int>& A,
        const int* list_node_i,
        const int numb_node_i,
        const int* list_node_p,
        const int numb_node_p ) {

  return 0;
}

________________________________________________________________________________


} // namespace Schur {
