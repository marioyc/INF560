/*!
*  @file Schur.cpp
*  @internal manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Fr�d�ric Magoul�s, Sonia Toubaline
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
  Aii.Allocate(numb_node_i,numb_node_i);

  for(int i = 0;i < numb_node_i;++i)
    for(int j = 0;j < numb_node_i;++j)
      Aii(i,j) = A(list_node_i[i],list_node_i[j]);

  App.Allocate(numb_node_p,numb_node_p);

  for(int i = 0;i < numb_node_p;++i)
    for(int j = 0;j < numb_node_[];++j)
      App(i,j) = A(list_node_p[i],list_node_p[j]);

  Aip.Allocate(numb_node_i,numb_node_p);

  for(int i = 0;i < numb_node_i;++i)
    for(int j = 0;j < numb_node_p;++j)
      Aip(i,j) = A(list_node_i[i],list_node_p[j]);

  Aip.Allocate(numb_node_p,numb_node_i);

  for(int i = 0;i < numb_node_p;++i)
    for(int j = 0;j < numb_node_i;++j)
      Api(i,j) = A(list_node_p[i],list_node_i[j]);

  return 0;
}

________________________________________________________________________________


} // namespace Schur {
