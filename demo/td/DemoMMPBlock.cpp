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
  MatrixDense<double,int> B_global;

  if ( proc_numb == proc_root ) {
    // -- allocate and fill A
    A_global.Allocate( size, size );
    for( int i = 0; i < size; i++ ) {
      for( int j = 0; j < size; j++ ) {
        A_global(i,j) = i * size + j;
      }
    }
    // -- allocate and fill B
    B_global.Allocate( size, size );
    for( int i = 0; i < size; i++ ) {
      for( int j = 0; j < size; j++ ) {
        B_global(i,j) = i * size + j;
      }
    }

  }

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // -- consider a square matrix (first)

  // -- creation of two-dimensional grid communicator and communicators
  //     for each row and each column of the grid
  MPI_Comm mpi_comm_rows;
  MPI_Comm mpi_comm_columns;
  DataTopology::GridCartesianComm(mpi_comm_rows, mpi_comm_columns, mpi_comm);

  // -- distribute matrix band-row
  MatrixDense<double,int> A_local;
  DataTopology::DistributeMatrixBlock( A_local, A_global, proc_root,
                                       mpi_comm_rows, mpi_comm_columns );


  // -- distribute matrix band-row
  MatrixDense<double,int> B_local;
  DataTopology::DistributeMatrixBlock( B_local, B_global, proc_root,
                                       mpi_comm_rows, mpi_comm_columns );


  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- compute C := A * B
  MatrixDense<double,int> C_local( A_local.GetNumbRows( ),
                                   B_local.GetNumbRows( ) );
  BlasMpi::MatrixMatrixProductBlock( C_local, A_local, B_local,
                                     proc_root, mpi_comm_rows, mpi_comm_columns );

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  MatrixDense<double,int> C_global;
  // ...

  DataTopology::AssembleMatrixBlock( C_global, C_local, 0, mpi_comm_rows, mpi_comm_columns );

  // -- print
  if ( proc_numb == proc_root && size < 20 ) {
    iomrg::printf( ">>> print A \n" );
    A_global.WriteToStdout( );
    iomrg::printf( ">>> print B \n" );
    B_global.WriteToStdout( );
    iomrg::printf( ">>> print C \n" );
    C_global.WriteToStdout( );
  }
  // -- write to csv
  if ( proc_numb == proc_root ) {
    C_global.WriteToFileCsv("mmp_block.csv");
  }

  // ---------------------------------------------------------------------------
  // -- finalize MPI
  // ---------------------------------------------------------------------------

  // -- finalizes MPI
  MPI_Finalize( );

  return 0;
}
