// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "Schur.hpp"

// third-party packages

int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------


  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------


  // -- matrix
  MatrixDense<double,int> A(5,5);
  A(0,0)= 1; A(0,1)=-1; A(0,2)=-1; A(0,3)=-1; A(0,4)=-1;
  A(1,0)=-1; A(1,1)= 2; A(1,2)= 0; A(1,3)= 0; A(1,4)= 0;
  A(2,0)=-1; A(2,1)= 0; A(2,2)= 3; A(2,3)= 1; A(2,4)= 1;
  A(3,0)=-1; A(3,1)= 0; A(3,2)= 1; A(3,3)= 4; A(3,4)= 2;
  A(4,0)=-1; A(4,1)= 0; A(4,2)= 1; A(4,3)= 2; A(1,4)= 5;

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  MatrixDense<double,int> Aii;
  MatrixDense<double,int> Aip;
  MatrixDense<double,int> Api;
  MatrixDense<double,int> App;

  const int numb_node_i= 3;
  const int numb_node_p= 2;

  int list_node_i[numb_node_i] = {4,1,2};

  int list_node_p[numb_node_p] = {3,0};


  Schur::SplitMatrixToBlock( Aii, Aip, Api, App, A,
                             list_node_i, numb_node_i,
                             list_node_p, numb_node_p );


  A.WriteToStdout('\n');
  printf("\n\n");
  Aii.WriteToStdout('\n');
  printf("\n\n");
  Aip.WriteToStdout('\n');
  printf("\n\n");
  Api.WriteToStdout('\n');
  printf("\n\n");
  App.WriteToStdout('\n');
  printf("\n\n");



  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------



  return 0;
}
