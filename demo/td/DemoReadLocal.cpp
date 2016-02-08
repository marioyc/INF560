// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"
#include "BlasMpi.hpp"

// third-party packages

int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- initialize MPI
  // ---------------------------------------------------------------------------

  // -- number of processors
  int numb_procs;
  // -- process number (process rank)
  int proc_numb;
  // -- starts MPI
  MPI_Init( &argc, &argv );
  // -- get the communicator
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process rank
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- help for io printing
  iomrg::g_log_numb_procs = numb_procs;
  iomrg::g_log_proc_numb = proc_numb;

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  // -- argv command line arguments management
  if ( argc < 4 ) {
    iomrg::printf( "Usage: %s type_mode file_group_name file_group_type ",
            argv[0] );
    iomrg::printf( " [image=false]\n");
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------
  int arg_num = 1;
  const char* argv_type_mode = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_group_name = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_group_type = argc>arg_num?argv[arg_num++]:NULL;
  const bool argv_opt_csv_matrix = argc>arg_num?atoi(argv[arg_num++]):0;

  // -- type mode
  int type_mode = -1;
  // -- type_mode == 1 (band-row)
  if ( strcmp( argv_type_mode, "band-row" ) == 0 ) type_mode = 1;
  // -- type_mode == 2 (band-column)
  if ( strcmp( argv_type_mode, "band-column" ) == 0 ) type_mode = 2;
  // -- type_mode == 3 (block)
  if ( strcmp( argv_type_mode, "block" ) == 0 ) type_mode = 3;

  // -- sample
  // argv_file_group_name = "mattest", argv_file_group_type = "csv"

  // -- format file name
  const char* prefix_mode = "";
  int numb_procs_i = numb_procs;
  int numb_procs_j = 0;
  if ( strcmp( argv_type_mode, "band-row" ) == 0 ) {
    prefix_mode = "-r";
  } else if ( strcmp( argv_type_mode, "band-column" ) == 0 ) {
    prefix_mode = "-c";
  } else if ( strcmp( argv_type_mode, "block" ) == 0 ) {
    prefix_mode = "-k";
    mathmrg::CreateDims( numb_procs_i, numb_procs_j , numb_procs );
  }
  // sample: in_file_name = "mattest-r"
  std::string inn_matrix_name = argv_file_group_name + std::string(prefix_mode);
  // sample: in_file_name = "mattest-r_x"
  std::string inn_vector_name = inn_matrix_name + "-x";
  // sample: in_file_name = "mattest-r_out"
  std::string outt_matrix_name = argv_file_group_name +
                                 std::string(prefix_mode) + "-out";
  // sample: in_file_name = "mattest-r-out_x"
  std::string outt_vector_name = outt_matrix_name + "-x";

  // -- format "matrix" file name upon number of processors and proc. number
  // sample: in_matrix_name = "mattest-r_4_000.csv"
  const char* in_matrix_name = NULL;
  in_matrix_name = DataTopology::GetProcFilename( inn_matrix_name.c_str( ),
                                                  argv_file_group_type,
                                                  proc_numb,
                                                  numb_procs_i, numb_procs_j );
  // -- format "in-vector" file name upon number of processors and proc. number
  // sample: in_vector_name = "mattest-r-x_4_000.csv"
  const char* in_vector_name = NULL;
  in_vector_name = DataTopology::GetProcFilename( inn_vector_name.c_str( ),
                                                  argv_file_group_type,
                                                  proc_numb,
                                                  numb_procs_i, numb_procs_j );
  // -- format "matrix" file name upon number of processors and proc. number
  // sample: out_matrix_name = "mattest-r-out_4_000.csv"
  const char* out_matrix_name = NULL;
  out_matrix_name = DataTopology::GetProcFilename( outt_matrix_name.c_str( ),
                                                   argv_file_group_type,
                                                   proc_numb,
                                                  numb_procs_i, numb_procs_j );
  // -- format "out-vector" file name upon number of processors and proc. number
  // sample: out_vector_name = "mattest-r-out-x_4_000.csv"
  const char* out_vector_name = NULL;
  out_vector_name = DataTopology::GetProcFilename( outt_vector_name.c_str( ),
                                                   argv_file_group_type,
                                                   proc_numb,
                                                   numb_procs_i, numb_procs_j );

  iomrg::printf("-------- in_matrix_name: %s\n", in_matrix_name);
  iomrg::printf("-------- in_vector_name: %s\n", in_vector_name);
  iomrg::printf("-------- out_matrix_name: %s\n", out_matrix_name);
  iomrg::printf("-------- out_vector_name: %s\n", out_vector_name);
//  // -- read (in) matrix
//  MatrixDense<double,int> A_local;
//  if ( argv_opt_csv_matrix == 0 ) {
//    // -- read image matrix from file csv
//    A_local.ReadFromFileCsv( in_matrix_name );
//  } else {
//    // -- read csv matrix from file csv
//    A_local.ReadImageMatrixFromFileCsv( in_matrix_name );
//  }

//  // -- read (in) vector
//  Vector<double,int> x_local;
//  x_local.ReadFromFileCsv( in_vector_name );

//  // -- try to wait all processors
//  MPI_Barrier( mpi_comm );

//  // ---------------------------------------------------------------------------
//  // -- processing
//  // ---------------------------------------------------------------------------

//  // -- root processor (default: 0)
//  const int p_root = 0;

//  Vector<double,int> x_global;
//  DataTopology::AssembleVectorBand( x_global, x_local, p_root, mpi_comm );


//  MatrixDense<double,int> A_global;
//  if ( strcmp( argv_type_mode, "band-row" ) == 0 ) {
//    DataTopology::AssembleMatrixBandRow( A_global, A_local, p_root, mpi_comm );
//  } else if ( strcmp( argv_type_mode, "band-column" ) == 0 ) {
//    DataTopology::AssembleMatrixBandColumn( A_global, A_local, p_root, mpi_comm );
//  } else if ( strcmp( argv_type_mode, "block" ) == 0 ) {
//    // -- creation of two-dimensional grid communicator and communicators
//    //     for each row and each column of the grid
//    MPI_Comm mpi_comm_rows;
//    MPI_Comm mpi_comm_columns;
//    DataTopology::GridCartesianComm(mpi_comm_rows, mpi_comm_columns, mpi_comm);
//    DataTopology::AssembleMatrixBlock( A_global, A_local, p_root,
//                                     mpi_comm_rows, mpi_comm_columns );
//  }


//  // ---------------------------------------------------------------------------
//  // -- post-processing
//  // ---------------------------------------------------------------------------


//  // -- print
//  if ( proc_numb == p_root && A_global.GetNumbRows( ) < 20 ) {
//    iomrg::printf( ">>> print A \n" );
//    A_global.WriteToStdout( );
//    iomrg::printf( ">>> print x \n" );
//    x_global.WriteToStdout( );
//  }
//  // -- write to csv
//  if ( proc_numb == p_root ) {
//    A_global.WriteToFileCsv( out_matrix_name );
//    x_global.WriteToFileCsv( out_matrix_name );
//  }

  // ---------------------------------------------------------------------------
  // -- finalize MPI
  // ---------------------------------------------------------------------------

  // -- finalizes MPI
  MPI_Finalize( );

  return 0;
}
