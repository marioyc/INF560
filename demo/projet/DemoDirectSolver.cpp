// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DirectSolver.hpp"

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
  A(4,0)=-1; A(4,1)= 0; A(4,2)= 1; A(4,3)= 2; A(4,4)= 5;

  // -- rhs
  Vector<double,int> rhs(5);
  rhs(0) = 1;
  rhs(1) = 1;
  rhs(2) = 1;
  rhs(3) = 1;
  rhs(4) = 1;

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------


  MatrixDense<double,int> LU;

  // -- LU decomposition A = L.U
  Factor::LU( LU, A );
  printf("-------------------->>>> LU \n");
  LU.WriteToStdout('\n');


//  Vector<double,int> lu;
//  DirectSolver::SolveLU( lu, A, rhs );

//  printf("\n lu \n");
//  lu.WriteToStdout('\n');

//  Vector<double,int> chol;
//  DirectSolver::SolveCholesky( chol, A, rhs );

//  printf("\n cholesky\n");
//  chol.WriteToStdout('\n');

  MatrixDense<double,int> ldlt;
  Factor::LDLt( ldlt, A );
  printf("-------------------->>>> LDLt \n");
  ldlt.WriteToStdout('\n');
  

//  printf("\n ldlt \n");
//  ldlt.WriteToStdout('\n');


  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------



  return 0;
}
