// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <string>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"

// third-party packages
#include "dllmrg.hpp"


int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  // -- argv command line arguments management
  if ( argc < 5 ) {
    printf( "Usage: %s type_mode file_group_name file_group_type numb_subdom",
            argv[0] );
    printf( " [csv=true]\n");
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------
  int arg_num = 1;
  const char* argv_type_mode = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_group_name = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_group_type = argc>arg_num?argv[arg_num++]:NULL;
  const int argv_numb_subdom = argc>arg_num?atoi(argv[arg_num++]):0;
  const bool argv_opt_csv_matrix = argc>arg_num?atoi(argv[arg_num++]):0;

  // -- type mode
  int type_mode = -1;
  // -- type_mode == 0 (all)
  if ( strcmp( argv_type_mode, "all" ) == 0 ) type_mode = 0;
  // -- type_mode == 1 (band-row)
  if ( strcmp( argv_type_mode, "band-row" ) == 0 ) type_mode = 1;
  // -- type_mode == 2 (band-column)
  if ( strcmp( argv_type_mode, "band-column" ) == 0 ) type_mode = 2;
  // -- type_mode == 3 (block)
  if ( strcmp( argv_type_mode, "block" ) == 0 ) type_mode = 3;

  // -- matrix file name
  std::string matrix_file_name = argv_file_group_name
                                 + std::string(".") + argv_file_group_type;

  // -- allocate and initialize Matrix
  MatrixDense<double,int> A;
  if ( argv_opt_csv_matrix == 0 ) {
    // -- read image matrix from file csv
    A.ReadFromFileCsv( matrix_file_name.c_str( ) );
  } else {
    // -- read csv matrix from file csv
    A.ReadImageMatrixFromFileCsv( matrix_file_name.c_str( ) );
  }

  std::string c="";
  const char separator_column = (argv_opt_csv_matrix==0) ? ' ' : '\0';
  const char separator_row = '\n';

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- for band-row
  if ( type_mode == 0 || type_mode == 1 ) {
    for ( int s = 0; s < argv_numb_subdom; s++ ) {
      MatrixDense<double,int> A_local;
      // -- band-row
      DataTopology::BuildMatrixBandRow( A_local, A, argv_numb_subdom, s );

      // -- format file name
      std::string proc_file_name = argv_file_group_name + std::string("-r");
      const char* output_file_name = NULL;
      output_file_name = DataTopology::GetProcFilename( proc_file_name.c_str( ),
                                                        argv_file_group_type,
                                                        s, argv_numb_subdom );
      A_local.WriteToFileCsv( output_file_name,
                              separator_column, separator_row );
    }
  }

  // -- for band-column
  if ( type_mode == 0 || type_mode == 2 ) {
    for ( int s = 0; s < argv_numb_subdom; s++ ) {
      MatrixDense<double,int> A_local;
      // -- band-column
      DataTopology::BuildMatrixBandColumn( A_local, A, argv_numb_subdom, s );

      // -- format file name
      std::string proc_file_name = argv_file_group_name + std::string("-c");
      const char* output_file_name = NULL;
      output_file_name = DataTopology::GetProcFilename( proc_file_name.c_str( ),
                                                        argv_file_group_type,
                                                        s, argv_numb_subdom );
      A_local.WriteToFileCsv( output_file_name,
                              separator_column, separator_row );
    }
  }

  // for matrix block
  if ( type_mode == 0 || type_mode == 3 ) {
    int numb_subdom_i = 0;
    int numb_subdom_j = 0;
    mathmrg::CreateDims( numb_subdom_i, numb_subdom_j , argv_numb_subdom );
    printf(" %d %d %d\n", numb_subdom_i, numb_subdom_j , argv_numb_subdom);
    for ( int s_i = 0; s_i < numb_subdom_i; s_i++ ) {
      for ( int s_j = 0; s_j < numb_subdom_j; s_j++ ) {
        MatrixDense<double,int> A_local;
        DataTopology::BuildMatrixMatrixBlock( A_local, A,
                                              numb_subdom_i,
                                              numb_subdom_j,
                                              s_i, s_j );

        // -- processor number
        int s = s_i * numb_subdom_j + s_j;
        // -- format file name
        std::string proc_file_name = argv_file_group_name + std::string("-k");
        const char* output_file_name = NULL;
        output_file_name = DataTopology::GetProcFilename(
                             proc_file_name.c_str( ), argv_file_group_type,
                             s, numb_subdom_i, numb_subdom_j );
        A_local.WriteToFileCsv( output_file_name,
                                separator_column, separator_row );
      }
    }
  }

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  return 0;
}
