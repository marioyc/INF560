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

  // -- size of problem
  const int size = (argv[1]!=NULL) ? atoi(argv[1]) : 5;
  // -- root processor (default: 0)
  const int proc_root = argv[2]!=NULL ? atoi(argv[2]) : 0;
  if( proc_numb == proc_root ) {
    iomrg::printf("-- problem size: %d [proc_root: %d]\n\n", size, proc_root );
  }

  // -- allocate and initialize Matrix and Vector
  MatrixDense<double,int> A_global;
  Vector<double,int> x_global;

  if ( proc_numb == proc_root ) {
    // -- allocate and fill A
    A_global.Allocate( size, size );
    for( int i = 0; i < size; i++ ) {
      for( int j = 0; j < size; j++ ) {
        A_global(i,j) = i * size + j;
      }
    }
    // -- allocate and fill x
    x_global.Allocate( size );
    for( int i = 0; i < size; i++ ) {
      x_global(i) = 1;
    }
  }

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // -- distribute matrix band-row
  MatrixDense<double,int> A_local;
  DataTopology::DistributeMatrixBandRow( A_local, A_global, proc_root, mpi_comm );

  // -- distribute vector
  Vector<double,int> x_local;
  DataTopology::DistributeVectorBand( x_local, x_global, proc_root, mpi_comm );

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- compute y := A * x
  Vector<double,int> y_local( A_local.GetNumbRows( ) );
  BlasMpi::MatrixVectorProductBandRow( y_local, A_local, x_local, mpi_comm );

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  Vector<double,int> y_global;
  DataTopology::AssembleVectorBand( y_global, y_local, proc_root, mpi_comm );

  // -- print
  if ( proc_numb == proc_root && size < 20 ) {
    iomrg::printf( ">>> print A \n" );
    A_global.WriteToStdout( );
    iomrg::printf( ">>> print x \n" );
    x_global.WriteToStdout( );
    iomrg::printf( ">>> print y \n" );
    y_global.WriteToStdout( );
  }
  // -- write to csv
  if ( proc_numb == proc_root ) {
    y_global.WriteToFileCsv("mvp_band_row.csv");
  }

  // ---------------------------------------------------------------------------
  // -- finalize MPI
  // ---------------------------------------------------------------------------

  // -- finalizes MPI
  MPI_Finalize( );

  return 0;
}
