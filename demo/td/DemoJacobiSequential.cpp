// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "Jacobi.hpp"

// third-party packages

int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------

  // -- argv command line arguments management
  if ( argc < 3 ) {
    printf( "Usage: %s max_numb_iter residual_threshold \
                           file_name file_type output_name", argv[0] );
    printf( " [image=false]\n");
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  int arg_num = 1;
  const char* argv_max_numb_iter = argc>arg_num?argv[arg_num++]:"100";
  const char* argv_residual_threshold = argc>arg_num?argv[arg_num++]:"1e-6";
  const char* argv_file_name = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_type = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_output_name = argc>arg_num?argv[arg_num++]:NULL;
  const bool argv_opt_csv_matrix = argc>arg_num?atoi(argv[arg_num++]):0;

  // sample: inn_matrix_name= "mattest.csv"
  std::string in_matrix_name = argv_file_name + std::string(".")
                                + argv_file_type;
  // sample: inn_vector_name = "mattest_b.csv"
  std::string in_vector_name = argv_file_name + std::string("_b.")
                                + argv_file_type;
  // sample: outt_vector_name = "mattest_x.csv"
  std::string out_name;
  if ( argv_output_name == NULL ) {
    out_name = argv_file_name;
  } else {
    out_name = argv_output_name;
  }

  std::string out_vector_name   = out_name + std::string("_x.") + argv_file_type;
  std::string out_residual_name = out_name +  std::string("_res.") + argv_file_type;
  std::string out_temporal_name = out_name +  std::string("_time.") + argv_file_type;

  // -- read (in) matrix
  MatrixDense<double,int> A;
  if ( argv_opt_csv_matrix == 0 ) {
    // -- read image matrix from file csv
    A.ReadFromFileCsv( in_matrix_name.c_str( ) );
  } else {
    // -- read csv matrix from file csv
    A.ReadImageMatrixFromFileCsv( in_matrix_name.c_str( ) );
  }

  // -- read (in) vector
  Vector<double,int> b;
  b.ReadFromFileCsv( in_vector_name.c_str( ) );

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- solve Ax = b, by Jacobi synchornous algorithmm
  // synchronous communication, synchronous iteration

  // -- conjugate gradient parameters
  const int max_numb_iter = atoi(argv_max_numb_iter);
  const double residual_threshold = atof(argv_residual_threshold);
  int numb_iter = 0;
  Vector<double,int> residual( max_numb_iter );
  Vector<double,int> temporal( max_numb_iter );

  printf("-- Paramètre du porblème: \n");
  printf("   .. max_numb_iter: %d \n", max_numb_iter);
  printf("   .. residual_threshold: %6e \n", residual_threshold);

  Vector<double,int> x( b.GetSize( ) );
  Solver::SequentialJacobi( x, numb_iter, residual, temporal, A, b,
                            max_numb_iter, residual_threshold );

  printf("-- Résultat: \n");
  printf("   .. numb_iter: %d \n", numb_iter);
  printf("   .. residual: %3e \n\n", residual[numb_iter-1]);


  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  // -- vérification:
  int norm_type = mathmrg::norm::c_LINF;
  // norm_b = ||b||
  double norm_b = b.NormLp( norm_type );

  // -- compute y := A * x - b
  Vector<double,int> y( A.GetNumbRows( ) );
  // y = A * x
  A.MatrixVectorProduct( y, x );
  // y = y - b
  y.Saxpy( -1, b );
  // norm_Axmb = ||y||
  double norm_Axmb = y.NormLp( norm_type );
  // norm_Axmb_s_b = ||Ax-b||/||b||
  double norm_Axmb_s_b = norm_Axmb/norm_b;

  printf("-- Norm ||Ax-b|| = %lf, ||b|| = %lf, ||Ax-b||/||b|| = %lf\n",
          norm_Axmb, norm_b, norm_Axmb_s_b );

  // -- print
  if ( x.GetSize( ) < 20 ) {
    printf( ">>> print A \n" );
    A.WriteToStdout( );
    printf( ">>> print x \n" );
    x.WriteToStdout( );
    printf( ">>> print b \n" );
    b.WriteToStdout( );
  }

  // -- write to csv
  x.WriteToFileCsv( out_vector_name.c_str( ) );
  // -- save residual
  residual.WriteToFileCsv(out_residual_name.c_str( ),'\n');
  // -- save temporal
  temporal.WriteToFileCsv(out_temporal_name.c_str( ),'\n');

  return 0;
}
