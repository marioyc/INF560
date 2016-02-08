// basic packages
#include <stdio.h>
#include <stdlib.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages

int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  // -- size of problem
  const int size = (argv[1]!=NULL) ? atoi(argv[1]) : 5;
  iomrg::printf("-- problem size: %d\n\n", size );

  // -- allocate and initialize Matrix and Vector
  MatrixDense<double,int> A;
  Vector<double,int> x;

  // A.ReadImageMatrixFromFileCsv("toto");

  // -- allocate and fill A
  A.Allocate( size, size );
  for( int i = 0; i < size; i++ ) {
    for( int j = 0; j < size; j++ ) {
      A(i,j) = i * size + j;
    }
  }
  // -- allocate and fill x
  x.Allocate( A.GetNumbColumns( ) );
  for( int i = 0; i < x.GetSize( ); i++ ) {
    x(i) = 1;
  }

 // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- compute y := A * x
  Vector<double,int> y( A.GetNumbRows( ) );
  A.MatrixVectorProduct( y, x );

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  // -- print
  if ( A.GetNumbRows( ) < 20 ) {
    iomrg::printf( ">>> print A \n" );
    A.WriteToStdout( );
    iomrg::printf( ">>> print x \n" );
    x.WriteToStdout( );
    iomrg::printf( ">>> print y \n" );
    y.WriteToStdout( );
  }
  // -- write to csv
  y.WriteToFileCsv("mvp_sequential.csv");

  return 0;
}
