// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"
#include "BlasMpi.hpp"
#include "BlasMpi1.hpp"
#include "Jacobi.hpp"

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
  // -- read command line arguments
  // ---------------------------------------------------------------------------

  // -- argv command line arguments management
  if ( argc < 3 ) {
    iomrg::printf( "Usage: %s max_numb_iter residual_threshold \
                           file_group_name file_group_type \
                           output_group_name", argv[0] );
    iomrg::printf( " [image=false]\n");
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  int arg_num = 1;
  const char* argv_max_numb_iter = argc>arg_num?argv[arg_num++]:"100";
  const char* argv_residual_threshold = argc>arg_num?argv[arg_num++]:"1e-6";
  const char* argv_file_group_name = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_file_group_type = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_output_group_name = argc>arg_num?argv[arg_num++]:NULL;
  const bool argv_opt_csv_matrix = argc>arg_num?atoi(argv[arg_num++]):0;
  const int argv_proc_root = argc>arg_num?atoi(argv[arg_num++]):0;

  // -- format file name
  const char* type_mode = "band-row";
  const char* prefix_mode = "-r";

  // sample: inn_matrix_name= "mattest-r"
  std::string inn_matrix_name = argv_file_group_name + std::string(prefix_mode);
  // sample: inn_vector_name = "mattest-r_x_local"
  std::string inn_vector_name = argv_file_group_name + std::string("_b")
                                + prefix_mode;
  // sample: outt_vector_name = "mattest-r-out_x_local"
  std::string outt_name;
  if ( argv_output_group_name == NULL ) {
    outt_name = argv_file_group_name;
  } else {
    outt_name = argv_output_group_name;
  }

  std::string outt_vector_name   = outt_name + std::string("_x") + prefix_mode;
  std::string outt_residual_name = outt_name +  std::string("_res") + prefix_mode;
  std::string outt_temporal_name = outt_name +  std::string("_time") + prefix_mode;

  // sample: outt_vector_name_global = "mattest-r-out_x_local.csv"
  std::string outt_vector_name_global = outt_vector_name + "." +
                                        argv_file_group_type;

  // -- format "matrix" file name upon number of processors and proc. number
  // sample: in_matrix_name = "mattest-r_4_000.csv"
  const char* in_matrix_name = NULL;
  in_matrix_name = DataTopology::GetProcFilename( inn_matrix_name.c_str( ),
                                                  argv_file_group_type,
                                                  proc_numb, numb_procs );
  // -- format "in-vector" file name upon number of processors and proc. number
  // sample: rhs_name = "mattest-r-x_local_4_000.csv"
  const char* rhs_name = NULL;
  rhs_name = DataTopology::GetProcFilename( inn_vector_name.c_str( ),
                                            argv_file_group_type,
                                            proc_numb, numb_procs );
  // -- format "out-vector" file name upon number of processors and proc. number
  // sample: solution_name = "mattest-r-out-x_local_4_000.csv"
  const char* solution_name = NULL;
  solution_name = DataTopology::GetProcFilename( outt_vector_name.c_str( ),
                                                 argv_file_group_type,
                                                 proc_numb, numb_procs );
  const char* residual_name = NULL;
  residual_name = DataTopology::GetProcFilename( outt_residual_name.c_str( ),
                                                 argv_file_group_type,

                                                 proc_numb, numb_procs );

  const char* temporal_name = NULL;
  temporal_name = DataTopology::GetProcFilename( outt_temporal_name.c_str( ),
                                                 argv_file_group_type,

                                                 proc_numb, numb_procs );

  // -- read (in) matrix: local
  MatrixDense<double,int> A_local;
  if ( argv_opt_csv_matrix == 0 ) {
    // -- read image matrix from file csv
    A_local.ReadFromFileCsv( in_matrix_name );
  } else {
    // -- read csv matrix from file csv
    A_local.ReadImageMatrixFromFileCsv( in_matrix_name );
  }

  // -- read (in) vector: local
  Vector<double,int> b_local;
  b_local.ReadFromFileCsv( rhs_name );

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- solve Ax = b_local, by Jacobi synchornous algorithmm
  // synchronous communication, synchronous iteration

  // -- conjugate gradient parameters
  const int max_numb_iter = atoi(argv_max_numb_iter);
  const double residual_threshold = atof(argv_residual_threshold);
  int numb_iter = 0;
  Vector<double,int> residual( max_numb_iter );
  Vector<double,int> temporal( max_numb_iter );

  if ( proc_numb == argv_proc_root ) {
    iomrg::printf("-- Paramètre du porblème: \n");
    iomrg::printf("   .. max_numb_iter: %d \n", max_numb_iter);
    iomrg::printf("   .. residual_threshold: %6e \n", residual_threshold);
  }
  Vector<double,int> x_local( b_local.GetSize( ) );
  Solver::JacobiBandRowAsyncIAsyncC( x_local,
                                     numb_iter, residual, temporal,
                                     A_local, b_local,
                                     max_numb_iter, residual_threshold,
                                     mpi_comm );

  if ( proc_numb == argv_proc_root ) {
    iomrg::printf("-- Résultat: \n");
    iomrg::printf("   .. numb_iter: %d \n", numb_iter);
    iomrg::printf("   .. residual: %3e \n\n", residual[numb_iter-1]);
  }

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );


  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  // -- vérification:
  int norm_type = mathmrg::norm::c_LINF;
  // norm_b_local = ||b_local||
  double norm_b_local = BlasMpi::NormLp( b_local, norm_type, mpi_comm );

  // -- compute y := A_local * x_local - b_local
  Vector<double,int> y( A_local.GetNumbRows( ) );
  // y = A_local * x_local
  BlasMpi::MatrixVectorProductBandRow( y, A_local, x_local, mpi_comm );
  // y = y - b_local
  y.Saxpy( -1, b_local );
  // norm_Axmb = ||y||
  double norm_Axmb = BlasMpi::NormLp( y, norm_type, mpi_comm );
  // norm_Axmb_s_b_local = ||Ax-b_local||/||b_local||
  double norm_Axmb_s_b_local = norm_Axmb/norm_b_local;

  if ( proc_numb == argv_proc_root ) {
    iomrg::printf("-- Norm ||Ax-b_local|| = %lf, ||b_local|| = %lf, ||Ax-b_local||/||b_local|| = %lf\n",
                  norm_Axmb, norm_b_local, norm_Axmb_s_b_local );
  }

  // -- assemble local to global
  Vector<double,int> x_global;
  DataTopology::AssembleVectorBand( x_global, x_local, argv_proc_root, mpi_comm );

  Vector<double,int> b_global;
  DataTopology::AssembleVectorBand( b_global, b_local, argv_proc_root, mpi_comm );

  MatrixDense<double,int> A_global;
  DataTopology::AssembleMatrixBandRow( A_global, A_local, argv_proc_root, mpi_comm );

  // -- print
  if ( proc_numb == argv_proc_root && x_global.GetSize( ) < 20 ) {
    iomrg::printf( ">>> print A_global \n" );
    A_global.WriteToStdout( );
    iomrg::printf( ">>> print x_global \n" );
    x_global.WriteToStdout( );
    iomrg::printf( ">>> print b_global \n" );
    b_global.WriteToStdout( );
  }
  // -- write to csv
  if ( proc_numb == argv_proc_root ) {
    x_global.WriteToFileCsv( outt_vector_name_global.c_str( ) );
  }
  // -- save residual
  residual.WriteToFileCsv(residual_name,'\n');
  // -- save temporal
  temporal.WriteToFileCsv(temporal_name,'\n');

  // ---------------------------------------------------------------------------
  // -- finalize MPI
  // ---------------------------------------------------------------------------

  // -- finalizes MPI
  MPI_Finalize( );

  return 0;
}
