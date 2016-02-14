// basic packages
#include <stdio.h>
#include <stdlib.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DirectSolver.hpp"

int main(int argc, char** argv){
  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------

  // -- argv command line arguments management
  if(argc < 3){
    printf("Usage: %s file_name output_name\n", argv[0]);
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  int arg_num = 1;
  const char* argv_file_name = argv[arg_num++];
  const char* argv_output_name = argv[arg_num++];

  // sample: inn_matrix_name= "mattest.csv"
  std::string in_matrix_name = argv_file_name + std::string(".csv");
  // sample: inn_vector_name = "mattest_b.csv"
  std::string in_vector_name = argv_file_name + std::string("_b.csv");

  std::string out_vector_name = argv_output_name + std::string("_x.csv");

  MatrixDense<double,int> A;
  A.ReadFromFileCsv(in_matrix_name.c_str());

  Vector<double,int> b;
  b.ReadFromFileCsv(in_vector_name.c_str());

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------


  MatrixDense<double,int> LU;

  // -- LU decomposition A = L.U
  Factor::LU(LU, A);


  Vector<double,int> x;
  DirectSolver::SolveLU(x, A, b);

  x.WriteToFileCsv(out_vector_name.c_str(), '\n');

  return 0;
}