/*!
*  @file DataTopology.cpp
*  @internal manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages
#include <cstdlib>

// project packages
#include "DataTopology.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace DataTopology
namespace DataTopology {

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- band distribution - topology
// -----------------------------------------------------------------------------

//! @internal index position of a given band
int BandPos (
        const int proc_numb,
        const int numb_procs,
        const int size ) {

  return (proc_numb * size) / numb_procs;
}

________________________________________________________________________________

//! @internal size of a given band
int BandSize (
        const int proc_numb,
        const int numb_procs,
        const int size ) {

  int band_size = BandPos( proc_numb+1, numb_procs, size ) -
                  BandPos( proc_numb, numb_procs, size );

  return band_size;
}

________________________________________________________________________________

//! @internal band list start and position
int BandTopology (
        int*& band_list_start,
        int*& band_list_size,
        const int size,
        const int numb_procs,
        const int shift ) {

  // -- band starts (displacement)
  band_list_start = new int[numb_procs];
  // -- band sizes (count)
  band_list_size = new int[numb_procs];

  for ( int k = 0; k < numb_procs; k++ ) {
    // -- start position
    band_list_start[k] = BandPos( k, numb_procs, size ) * shift;
    // -- size of the block
    band_list_size[k] = BandSize( k, numb_procs, size ) * shift;
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- processors distribution - topology
// -----------------------------------------------------------------------------

//! @internal creation of two-dimensional grid communicator and communicators
//!         for each row and each column of the grid
int GridCartesianComm (
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );

  // -- number of cartesian dimensions (integer)
  const int numb_dims = 2;

  // MPI_Dims_create - MPICH
  // MPI_Dims_create
  // Creates a division of processors in a cartesian grid
  //
  // Synopsis
  // int MPI_Dims_create(int nnodes, int ndims, int dims[])
  //
  // Input Parameters
  // nnodes
  //    number of nodes in a grid (integer)
  // ndims
  //    number of cartesian dimensions (integer)
  //
  // Input/Output Parameters
  // dims
  //    integer array of size ndims specifying the number of nodes in each dimension.
  //    A value of 0 indicates that MPI_Dims_create should fill in a suitable value.
  int grid_dims[numb_dims];
  grid_dims[0] = 0;
  grid_dims[1] = 0;
  MPI_Dims_create( numb_procs, numb_dims, grid_dims );

  // MPI_Cart_create - MPICH
  // MPI_Cart_create
  //   Makes a new communicator to which topology information has been attached
  //
  // Synopsis
  // int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[],
  //                    const int periods[], int reorder, MPI_Comm *comm_cart)
  //
  // Input Parameters
  // comm_old
  //    input communicator (handle)
  // ndims
  //    number of dimensions of cartesian grid (integer)
  // dims
  //    integer array of size ndims specifying the number of processes in each dimension
  // periods
  //    logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
  // reorder
  //    ranking may be reordered (true) or not (false) (logical)
  //
  // Output Parameters
  // comm_cart
  //    communicator with new cartesian topology (handle)
  MPI_Comm mpi_comm_cart;
  int periods[numb_dims];
  periods[0] = 0;
  periods[1] = 0;
  MPI_Cart_create( mpi_comm, numb_dims, grid_dims, periods, 1, &mpi_comm_cart);

  // -- process number (process proc_numb)
  int grid_proc_numb;
  // -- get current process proc_numb
  MPI_Comm_rank( mpi_comm_cart, &grid_proc_numb );

  // MPI_Cart_coords - MPICH
  // MPI_Cart_coords
  //    Determines process coords in cartesian topology given rank in group
  //
  // Synopsis
  // int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[])
  //
  // Input Parameters
  // comm
  //    communicator with cartesian structure (handle)
  // rank
  //    rank of a process within group of comm (integer)
  // maxdims
  //    length of vector coords in the calling program (integer)
  //
  // Output Parameters
  // coords
  //    integer array (of size ndims) containing the Cartesian coordinates of specified process (integer)
  // -- get local coordinates on grid; cartesian coordinates of specified process
  int grid_coords[numb_dims];
  grid_coords[0] = 0;
  grid_coords[1] = 0;
  MPI_Cart_coords( mpi_comm_cart, grid_proc_numb, numb_dims, grid_coords );

  // MPI_Comm_split - MPICH
  // MPI_Comm_split
  //    Creates new communicators based on colors and keys
  //
  // Synopsis
  // int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
  //
  // Input Parameters
  // comm
  //    communicator (handle)
  // color
  //    control of subset assignment (nonnegative integer). Processes with the same color are in the same new communicator
  // key
  //    control of rank assignment (integer)
  //
  // Output Parameters
  // newcomm
  //    new communicator (handle)
  // -- create communicators that work across entire rows or columns
  // Example: 3 x 3 grid
  // 0(0,0) 1(0,1) 2(0,2)
  // 3(1,0) 4(1,1) 5(1,2)
  // 6(2,0) 7(2,1) 8(2,2)
  MPI_Comm_split( mpi_comm_cart, grid_coords[0], grid_coords[1],
                  &mpi_comm_rows );
  MPI_Comm_split( mpi_comm_cart, grid_coords[1], grid_coords[0],
                  &mpi_comm_columns );

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Vector: BAND-ROW or BAND-COLUMN
// -----------------------------------------------------------------------------

//! @internal distribute vector upon processors
int DistributeVectorBand (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- size of global vector
  int size_global = x_global.GetSize( );

  // MPI_Bcast - MPICH
  // MPI_Bcast
  //    Broadcasts a message from the process with rank "root" to all other processes of the communicator
  //
  // Synopsis
  // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,
  //               MPI_Comm comm )
  //
  // Input/Output Parameters
  // buffer
  //    starting address of buffer (choice)
  //
  // Input Parameters
  // count
  //    number of entries in buffer (integer)
  // datatype
  //    data type of buffer (handle)
  // root
  //    rank of broadcast root (integer)
  // comm
  //    communicator (handle)
  MPI_Bcast( &size_global, 1, MPI_INT, root, mpi_comm );

  // -- compute local size (band-row splitting)
  int size_local = BandSize( proc_numb, numb_procs, size_global );

  // -- local vector
  x_local.Allocate( size_local );

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size, size_global, numb_procs );

  // MPI_Scatterv - MPICH
  // MPI_Scatterv
  //    Scatters a buffer in parts to all processes in a communicator
  //
  // Synopsis
  // int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
  //                  MPI_Datatype sendtype, void *recvbuf, int recvcount,
  //                  MPI_Datatype recvtype,
  //                  int root, MPI_Comm comm)
  //
  // Input Parameters
  // sendbuf
  //    address of send buffer (choice, significant only at root)
  // sendcounts
  //    integer array (of length group size) specifying the number of elements to send to each processor
  // displs
  //    integer array (of length group size). Entry i specifies the displacement (relative to sendbuf from which to take the outgoing data to process i
  // sendtype
  //    data type of send buffer elements (handle)
  // recvcount
  //    number of elements in receive buffer (integer)
  // recvtype
  //    data type of receive buffer elements (handle)
  // root
  //    rank of sending process (integer)
  // comm
  //    communicator (handle)
  //
  // Output Parameters
  // recvbuf
  //    address of receive buffer (choice)
  // -- distribute from 'root', the global vector upon processors
  MPI_Scatterv( x_global.GetCoef( ), band_list_size, band_list_start, MPI_DOUBLE,
                x_local.GetCoef( ), size_local, MPI_DOUBLE,
                root, mpi_comm );

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // MPI_Barrier - MPICH
  // MPI_Barrier
  //    Blocks until all processes in the communicator have reached this routine.
  //
  // Synopsis
  // int MPI_Barrier( MPI_Comm comm )
  //
  // Input Parameters
  // comm
  //    communicator (handle)
  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal assemble vector upon processors
int AssembleVectorBand (
        Vector<double,int>& x_global,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- size of local vector
  int size_local = x_local.GetSize( );

  // -- get global size (band-row splitting)
  int size_global = 0;

  // MPI_Reduce - MPICH
  // MPI_Reduce
  // Reduces values on all processes to a single value
  //
  // Synopsis
  // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
  //                MPI_Op op, int root, MPI_Comm comm)
  //
  // Input Parameters
  // sendbuf
  //     address of send buffer (choice)
  // count
  //     number of elements in send buffer (integer)
  // datatype
  //     data type of elements of send buffer (handle)
  // op
  //     reduce operation (handle)
  // root
  //     rank of root process (integer)
  // comm
  //     communicator (handle)
  //
  // Output Parameters
  // recvbuf
  //     address of receive buffer (choice, significant only at root)
  MPI_Reduce( &size_local, &size_global, 1, MPI_INT, MPI_SUM, root, mpi_comm );

  // -- global vector
  if ( proc_numb == root ) {
    x_global.Allocate( size_global );
  }

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size, size_global, numb_procs );

  // MPI_Gatherv - MPICH
  // MPI_Gatherv
  // Gathers into specified locations from all processes in a group
  // Synopsis
  // int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
  //                 void *recvbuf, const int *recvcounts, const int *displs,
  //                 MPI_Datatype recvtype, int root, MPI_Comm comm)
  //
  // Input Parameters
  // sendbuf
  //     starting address of send buffer (choice)
  // sendcount
  //     number of elements in send buffer (integer)
  // sendtype
  //     data type of send buffer elements (handle)
  // recvcounts
  //     integer array (of length group size) containing the number of elements that are received from each process (significant only at root)
  // displs
  //     integer array (of length group size). Entry i specifies the displacement relative to recvbuf at which to place the incoming data from process i (significant only at root)
  // recvtype
  //     data type of recv buffer elements (significant only at root) (handle)
  // root
  //     rank of receiving process (integer)
  // comm
  //     communicator (handle)
  //
  // Output Parameters
  // recvbuf
  //     address of receive buffer (choice, significant only at root)
  // -- distribute from 'root', the global vector upon processors
  MPI_Gatherv( x_local.GetCoef( ), size_local, MPI_DOUBLE,
               x_global.GetCoef( ), band_list_size, band_list_start, MPI_DOUBLE,
               root, mpi_comm );

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal build 'band_numb'-th band
int BuildVectorBand (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int numb_procs,
        const int proc_numb ) {

  // -- size of the global vector
  int size_global = x_global.GetSize( );

  // -- compute the number of rows per band
  int size_local = BandSize( proc_numb, numb_procs, size_global );
  // -- allocate and fill out the local vector
  x_local.Allocate( size_local );
  for( int i = 0; i < size_local; i++ ) {
    // local (i) to global indice (l2g_i)
    int l2g_i = BandPos( proc_numb, numb_procs, size_global ) + i;
    x_local(i) = x_global(l2g_i);
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Matrix: BAND-ROW
// -----------------------------------------------------------------------------

//! @internal distribute matrix upon processors (band row)
int DistributeMatrixBandRow (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process number proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- global matrix: number of rows
  // all processors except root get 0
  int numb_rows_global    = A_global.GetNumbRows( );
  // -- global matrix: number of columns
  // all processors except root get 0
  int numb_columns_global = A_global.GetNumbColumns( );

  // -- broadcast numb_rows_global
  MPI_Bcast( &numb_rows_global, 1, MPI_INT, root, mpi_comm );
  // -- broadcast numb_columns_global
  MPI_Bcast( &numb_columns_global, 1, MPI_INT, root, mpi_comm );

  // -- compute local size (numb_rows_local) (band-row splitting)
  int numb_rows_local = BandSize( proc_numb, numb_procs, numb_rows_global );

  // -- local matrix
  A_local.Allocate( numb_rows_local, numb_columns_global );

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size, numb_rows_global, numb_procs, numb_columns_global );

  // -- distribute from 'root', the global vector upon processors
  MPI_Scatterv( A_global.GetCoef( ), band_list_size, band_list_start, MPI_DOUBLE,
                A_local.GetCoef( ), numb_rows_local * numb_columns_global, MPI_DOUBLE,
                root, mpi_comm );

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal assemble matrix upon processors (band row)
int AssembleMatrixBandRow (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process number proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- local matrix: number of rows
  int numb_rows_local    = A_local.GetNumbRows( );
  // -- local matrix: number of columns
  int numb_columns_local = A_local.GetNumbColumns( );

  // -- global matrix: number of rows
  int numb_rows_global = 0;
  MPI_Reduce( &numb_rows_local, &numb_rows_global, 1, MPI_INT, MPI_SUM,
              root, mpi_comm );

  // -- global matrix
  if ( proc_numb == root ) {
    A_global.Allocate( numb_rows_global, numb_columns_local );
  }

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size,
                numb_rows_global, numb_procs, numb_columns_local );

  // -- assemble local matrices to 'root', upon processors
  MPI_Gatherv( A_local.GetCoef( ), numb_rows_local * numb_columns_local, MPI_DOUBLE,
               A_global.GetCoef( ), band_list_size, band_list_start, MPI_DOUBLE,
               root, mpi_comm );

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal build 'band_numb'-th band (row) of a matrix
int BuildMatrixBandRow (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs,
        const int proc_numb ) {

  // -- global matrix: number of rows
  int numb_rows_global    = A_global.GetNumbRows( );
  // -- global matrix: number of columns
  int numb_columns_global = A_global.GetNumbColumns( );

  // -- compute the number of rows per band
  int numb_rows_local = BandSize( proc_numb, numb_procs, numb_rows_global );
  // -- allocate and fill out the local matrix
  A_local.Allocate( numb_rows_local, numb_columns_global );
  // -- compute band matrix from global matrix
  for( int i = 0; i < numb_rows_local; i++ ) {
    // local (i) to global indice (l2g_i)
    int l2g_i = BandPos( proc_numb, numb_procs, numb_rows_global ) + i;
    for( int j = 0; j < numb_columns_global; j++ ) {
      A_local(i,j) = A_global(l2g_i, j);
    }
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Matrix: BAND-COLUMN
// -----------------------------------------------------------------------------

//! @internal distribute matrix upon processors (band column)
int DistributeMatrixBandColumn (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process number proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- global matrix: number of rows
  // all processors except root get 0
  int numb_rows_global    = A_global.GetNumbRows( );
  // -- global matrix: number of columns
  // all processors except root get 0
  int numb_columns_global = A_global.GetNumbColumns( );

  // -- broadcast numb_rows
  MPI_Bcast( &numb_rows_global, 1, MPI_INT, root, mpi_comm );
  // -- broadcast numb_columns
  MPI_Bcast( &numb_columns_global, 1, MPI_INT, root, mpi_comm );

  // -- compute local size (numb_columns_local) (band-column splitting)
  int numb_columns_local = BandSize( proc_numb, numb_procs, numb_columns_global );

  // -- local matrix
  A_local.Allocate( numb_rows_global, numb_columns_local );

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size, numb_columns_global, numb_procs );

  // -- for each row, do a distribute
  for ( int i = 0; i < numb_rows_global; i++ ) {
    // -- distribute from 'root', the global vector upon processors
    MPI_Scatterv( A_global.GetCoef( i ), band_list_size, band_list_start, MPI_DOUBLE,
                  A_local.GetCoef( i ), numb_columns_local, MPI_DOUBLE,
                  root, mpi_comm );

  }

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal assemble matrix upon processors (band column)
int AssembleMatrixBandColumn (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process number proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- local matrix: number of rows
  int numb_rows_local    = A_local.GetNumbRows( );
  // -- local matrix: number of columns
  int numb_columns_local = A_local.GetNumbColumns( );

  // -- global matrix: number of columns
  int numb_columns_global = 0;
  MPI_Reduce( &numb_columns_local, &numb_columns_global, 1, MPI_INT, MPI_SUM,
              root, mpi_comm );

  // -- global matrix
  if ( proc_numb == root ) {
    A_global.Allocate( numb_rows_local, numb_columns_global );
  }

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  BandTopology( band_list_start, band_list_size,
                numb_columns_global, numb_procs );

  // -- for each line, do a distribute
  for ( int i = 0; i < numb_rows_local; i++ ) {
    // -- assemble local matrices to 'root', upon processors
    MPI_Gatherv( A_local.GetCoef( i ), numb_columns_local, MPI_DOUBLE,
                 A_global.GetCoef( i ), band_list_size, band_list_start, MPI_DOUBLE,
                 root, mpi_comm );
  }

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal build 'band_numb'-th band (column) of a matrix
int BuildMatrixBandColumn (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs,
        const int proc_numb ) {

  // -- global matrix: number of rows
  int numb_rows_global    = A_global.GetNumbRows( );
  // -- global matrix: number of columns
  int numb_columns_global = A_global.GetNumbColumns( );

  // -- compute the number of columns per band
  int numb_columns_local = BandSize( proc_numb, numb_procs, numb_columns_global );
  // -- allocate and fill out the local matrix
  A_local.Allocate( numb_rows_global, numb_columns_local );
  // -- compute band matrix from global matrix
  for( int i = 0; i < numb_rows_global; i++ ) {
    for( int j = 0; j < numb_columns_local; j++ ) {
      // local (j) to global indice (l2g_j)
      int l2g_j = BandPos( proc_numb, numb_procs, numb_columns_global ) + j;
      A_local(i,j) = A_global(i,l2g_j);
    }
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Vector : BLOCK
// -----------------------------------------------------------------------------

//! @internal distribute vector upon processors (block)
int DistributeVectorBlock (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) {

  // -- number of processors (i-)
  int numb_procs_i;
  // -- process number (process proc_numb_i)
  int proc_numb_i;
  // -- get number of processes
  MPI_Comm_size( mpi_comm_rows, &numb_procs_i );
  // -- get current process proc_numb_i
  MPI_Comm_rank( mpi_comm_rows, &proc_numb_i );

  // -- number of processors (j-)
  int numb_procs_j;
  // -- process number (process proc_numb_j)
  int proc_numb_j;
  // -- get number of processes
  MPI_Comm_size( mpi_comm_columns, &numb_procs_j );
  // -- get current process proc_numb_j
  MPI_Comm_rank( mpi_comm_columns, &proc_numb_j );

  // -- global: size of the vector
  int size = x_global.GetSize( );

  // -- distribute will be done on matrix vector product
  // allocate x_local for processor in mpi_comm_rows
  DataTopology::DistributeVectorBand( x_local, x_global, root, mpi_comm_rows );

  // -- try to wait all row processors
  MPI_Barrier( mpi_comm_rows );

  // broadcast
  MPI_Bcast( &size, 1, MPI_INT, root, mpi_comm_rows );
  MPI_Bcast( &size, 1, MPI_INT, root, mpi_comm_columns );

  // -- compute local size (band-row splitting)
  int size_local = BandSize( proc_numb_i, numb_procs_i, size );

  // -- every process in the first row shares its share of x
  //    with every process in its column
  // -- only for processors in mpi_comm_rows
  if ( !x_local.Status( ) ) {
    x_local.Allocate( size_local );
  }
  MPI_Bcast( x_local.GetCoef( ), x_local.GetSize( ), MPI_DOUBLE,
             root, mpi_comm_columns );

  // -- try to wait all columns processors
  MPI_Barrier( mpi_comm_columns );

  return 0;
}

________________________________________________________________________________

//! @internal build the 'block_numb'-th band of a vector (block)
int BuildVectorBlock (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int numb_procs_i,
        const int numb_procs_j,
        const int proc_numb_i,
        const int proc_numb_j ) {

  // -- size of the global vector
  int size = x_global.GetSize( );

  // -- compute the number of block in j-
  int size_local = BandSize( proc_numb_j, numb_procs_j, size );
  // -- allocate and fill out the local vector
  x_local.Allocate( size_local );
  for( int j = 0; j < size_local; j++ ) {
    // local (j) to global indice (l2g_j)
    int l2g_j = BandPos( proc_numb_j, numb_procs_j, size ) + j;
    x_local(j) = x_global(l2g_j);
  }

  return 0;
}

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- Matrix: BLOCK
// -----------------------------------------------------------------------------

//! @internal distribute matrix upon processors (block)
//! @remarks number of processors must be perfect square
int DistributeMatrixBlock (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) {

  // -- local band row block matrix
  MatrixDense<double,int> A_row;
  int root_row = 0;
  DataTopology::DistributeMatrixBandRow( A_row, A_global, root_row,
                                         mpi_comm_columns );

  // -- distribute A_row to A_local
  DataTopology::DistributeMatrixBandColumn( A_local, A_row, root,
                                            mpi_comm_rows );

  // -- try to wait all row processors
  MPI_Barrier( mpi_comm_rows );

  // -- try to wait all columns processors
  MPI_Barrier( mpi_comm_columns );

  return 0;
}

________________________________________________________________________________

//! @internal assemble matrix upon processors (block)
//! @remarks number of processors must be perfect square
int AssembleMatrixBlock (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) {

  // -- gloabl band row block matrix
  MatrixDense<double,int> A_row;
  DataTopology::AssembleMatrixBandColumn( A_row, A_local, root, mpi_comm_rows );

  // -- assemble A to A_global
  int root_row = 0;
  DataTopology::AssembleMatrixBandRow( A_global, A_row, root_row,
                                       mpi_comm_columns );

  return 0;
}

________________________________________________________________________________

//! @internal build the 'band_numb'-th block of a matrix
//! @remarks proc_numb follows row major order
//! @remarks proc_numb = proc_numb_i * numb_procs_j + proc_numb_j
//! @note number of processors must be perfect square
int BuildMatrixMatrixBlock (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs_i,
        const int numb_procs_j,
        const int proc_numb_i,
        const int proc_numb_j ) {

  // -- global matrix: number of rows
  int numb_rows    = A_global.GetNumbRows( );
  // -- global matrix: number of columns
  int numb_columns = A_global.GetNumbColumns( );

  // -- compute the number of rows per block
  int numb_rows_local = BandSize( proc_numb_i, numb_procs_i, numb_rows );
  // -- compute the number of columns per block
  int numb_columns_local = BandSize( proc_numb_j, numb_procs_j, numb_columns );

  // -- allocate and fill out the local matrix
  A_local.Allocate( numb_rows_local, numb_columns_local );
  // -- compute band matrix from global matrix
  for( int i = 0; i < numb_rows_local; i++ ) {
    // local (i) to global indice (l2g_i)
    int l2g_i = BandPos( proc_numb_i, numb_procs_i, numb_rows ) + i;
    for( int j = 0; j < numb_columns_local; j++ ) {
      // local (j) to global indice (l2g_j)
      int l2g_j = BandPos( proc_numb_j, numb_procs_j, numb_columns ) + j;
      A_local(i,j) = A_global(l2g_i,l2g_j);
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
// -- Read Local
// -----------------------------------------------------------------------------

________________________________________________________________________________

//! @internal get filename of given processor among others
const char* GetProcFilename (
        const char* file_group_name,
        const char* file_group_type,
        const int proc_numb,
        const int numb_procs_i,
        const int numb_procs_j ) {

  // -- format file name
  std::string numb_procs_str = "000";

  int numb_procs = numb_procs_i;
  if ( numb_procs_j != 0 ) {
    numb_procs = numb_procs_i * numb_procs_j;
  }

  // sample: numb_procs = 9
  std::stringstream ss_numb_procs;
  // sample: ss_numb_procs = "9"
  ss_numb_procs << numb_procs;

  // sample: numb_procs_i = 3
  std::stringstream ss_numb_procs_i;
  // sample: ss_numb_procs_i = "3"
  ss_numb_procs_i << numb_procs_i;

  // sample: numb_procs_j = 3
  std::stringstream ss_numb_procs_j;
  // sample: ss_numb_procs_j = "3"
  ss_numb_procs_j << numb_procs_j;

  // -- format file name
  // sample: proc_numb = 1
  std::stringstream ss_proc_numb;
  // sample: ss_proc_numb = "1"
  ss_proc_numb << proc_numb;
  // sample: l_ns = 1
  int l_ns = std::string( ss_proc_numb.str() ).size();
  // sample: prefix_str = "3_001"
  std::string prefix_str =  ss_numb_procs.str() + "_";

  if ( numb_procs_j != 0 ) {
    prefix_str += ss_numb_procs_i.str() + "x" + ss_numb_procs_j.str() + "_";
  }
  prefix_str += numb_procs_str.substr(l_ns) + ss_proc_numb.str();

    // output filename
  std::string output_file_name = std::string(file_group_name)
                                       + "_" + prefix_str + "."
                                       + std::string(file_group_type);

  return strdup( output_file_name.c_str( ) );
}

________________________________________________________________________________

//! @internal read local vector
int ReadLocalVectorFromFile (
        Vector<double,int>& x_local,
        const char* file_group_name,
        const char* file_group_type,
        const int numb_procs,
        const int proc_numb ) {

  // -- format file name
  const char* output_file_name = NULL;
  output_file_name = GetProcFilename( file_group_name, file_group_type,
                                      numb_procs, proc_numb );

  // -- read vector
  x_local.ReadFromFileCsv( output_file_name );

  return 0;
}

________________________________________________________________________________

//! @internal read local matrix
int ReadLocalMatrixFromFile (
        MatrixDense<double,int>& A_local,
        const char* file_group_name,
        const char* file_group_type,
        const int numb_procs,
        const int proc_numb,
        const bool opt_image_matrix ) {

  // -- format file name
  const char* output_file_name = NULL;
  output_file_name = GetProcFilename( file_group_name, file_group_type,
                                      numb_procs, proc_numb );

  // -- read matrix
  if ( opt_image_matrix ) {
    // -- read image matrix from file csv
    A_local.ReadImageMatrixFromFileCsv( output_file_name );
  } else {
    // -- read image matrix from file csv
    A_local.ReadFromFileCsv( output_file_name );
  }

  return 0;
}

________________________________________________________________________________


//! @internal read l2g from file
int ReadL2gFromFile (
        int& numb_global,
        int& numb_l2g,
        int*& l2g,
        const char* file_name ) {

  iomrg::printf(".r. Reading L2g file: %s\n", file_name );

  // open file
  std::ifstream f_l2g( file_name, std::ios::in  );

  // read mrg information
  std::string current_line;
  int real_base = 0;
  const int standard_base = 0;
  iomrg::ReadMrgHeader( f_l2g, current_line, real_base, standard_base, '%' );

  // -- read new line
  std::getline( f_l2g, current_line );
  std::stringstream ss_current_line(current_line);
  // first line
  ss_current_line >> numb_l2g >> numb_global;

  int offset = real_base;

  l2g = new int[numb_l2g];
  for( int i = 0; i < numb_l2g; i++ ) {
    int value;
    f_l2g >> value;
    l2g[i] = value - offset;
  }

  // close file
  f_l2g.close( );

  return 0;
}

________________________________________________________________________________


//! @internal read interfaces from file
int ReadNeighb2InterfaceNodeFromFile (
        int& subdom_numb,
        int& numb_subdom,
        int& numb_neighb_subdom,
        int*& list_neighb_subdom,
        int*& p_neighb2interfnode,
        int*& neighb2interfnode,
        int*& neighb2interfnode_multiplicity,
        const char* file_name ) {

  iomrg::printf(".r. Reading interface file: %s\n", file_name );

  // -- open file
  std::ifstream f_interf( file_name, std::ios::in  );

  // -- field header variable: std::string field_header[3]
  // >> FIELD NAME_GROintP_FIELD #FIELD
  // >> field_header[0] := "FIELD"
  // >> field_header[1] := "NAME_GROintP_FIELD"
  // >> field_header[2] := "#FILED"
  std::string field_header[3];
  // -- field info variable: std::string field_info[4]
  // >> name_field dimension #line int
  // >> field_info[0] := "name_field"
  // >> field_info[1] := "dimension"
  // >> field_info[2] := "#line"
  // >> field_info[3] := "int"
  std::string field_info[4];

  // ---------------------------------------------------------------------------
  // -- first line
  // ---------------------------------------------------------------------------

  // read mrg information
  std::string current_line;
  int real_base = 0;
  const int standard_base = 0;
  iomrg::ReadMrgHeader( f_interf, current_line, real_base, standard_base, '%' );

  // -- read new line
  std::getline( f_interf, current_line );

  int offset = real_base;

  // ---------------------------------------------------------------------------
  // -- FIELD SintBDOMAIN 5
  // ---------------------------------------------------------------------------

  // -- FIELD NAME_GROintP_FIELD #FIELD
  std::stringstream ss_current_line(current_line);
  ss_current_line >> field_header[0] >> field_header[1] >> field_header[2];

  int numb_subdom_field = atoi( field_header[2].c_str( ) );

  // -- name_field dimension #line int
  // -- field #1: subdom_numb
  f_interf >> field_info[0] >> field_info[1] >> field_info[2] >> field_info[3];
  f_interf >> subdom_numb;

  // -- field #2: numb_subdom
  f_interf >> field_info[0] >> field_info[1] >> field_info[2] >> field_info[3];

  f_interf >> numb_subdom;

  // -- field #3: numb_neighb_subdom
  f_interf >> field_info[0] >> field_info[1] >> field_info[2] >> field_info[3];

  f_interf >> numb_neighb_subdom;

  // -- field #4: list_neighb_subdom
  f_interf >> field_info[0] >> field_info[1] >> field_info[2] >> field_info[3];

  int tmp_numb_neighb_subdom = atoi( field_info[2].c_str() );

  // read list of neighboring subdomain
  list_neighb_subdom = new int[numb_neighb_subdom];
  for( int i = 0; i < numb_neighb_subdom; i++ ) {
    int neighb_subdom_node;
    f_interf >> neighb_subdom_node;
    list_neighb_subdom[i] = neighb_subdom_node - offset;
  }

  // -- field #5: numb_total_neighb_node
  f_interf >> field_info[0] >> field_info[1] >> field_info[2] >> field_info[3];

  int numb_total_neighb_node = 0;
  f_interf >> numb_total_neighb_node;

  // ---------------------------------------------------------------------------
  // -- FIELD INTERFACES #INTERFACE
  // ---------------------------------------------------------------------------

  // -- FIELD INTERFACES #INTERFACE
  f_interf >> field_header[0] >> field_header[1] >> field_header[2];

  int numb_interface_field = atoi( field_header[2].c_str( ) );

  // -- build p_neighb2interfnode
  p_neighb2interfnode = new int[numb_neighb_subdom+1];
  // -- build for each neighbor the list of nodes
  neighb2interfnode = new int[numb_total_neighb_node];
  // -- build for each neighbor the list of node multiplicities
  neighb2interfnode_multiplicity = new int[numb_total_neighb_node];

  p_neighb2interfnode[0] = 0;
  for( int i = 0; i < numb_interface_field; i++ ) {
    int ii = p_neighb2interfnode[i];
    // -- field #i: interface#list_neighb_subdom[i]
    f_interf >> field_info[0] >> field_info[1]
             >> field_info[2] >> field_info[3];
    // convert neighboring subdomain to string
    std::stringstream ss_neighb_subdom_numb;
    ss_neighb_subdom_numb << list_neighb_subdom[i];
    std::string interface_name = "interface" + ss_neighb_subdom_numb.str();

    int numb_neighb_node = atoi( field_info[2].c_str( ) );

    // read node and multiplicity: node# multiplicity
    for( int j = 0; j < numb_neighb_node; j++ ) {
      int interfnode_node;
      f_interf >> interfnode_node;
      neighb2interfnode[j + ii] = interfnode_node - offset;
      f_interf >> neighb2interfnode_multiplicity[j + ii];
    }
    // increment the number total of neighboring node
    p_neighb2interfnode[i+1] = ii + numb_neighb_node;
  }

  // close file
  f_interf.close( );

  return 0;
}

________________________________________________________________________________

//! @internal write interfaces to stdout
int WriteNeighb2InterfaceNodeToStdout (
        const int subdom_numb,
        const int numb_subdom,
        const int numb_neighb_subdom,
        const int* list_neighb_subdom,
        const int* p_neighb2interfnode,
        const int* neighb2interfnode,
        const int* neighb2interfnode_multiplicity ) {

  iomrg::printf("-- subdom_numb        : %d\n", subdom_numb);
  iomrg::printf("-- numb_subdom        : %d\n", numb_subdom);
  iomrg::printf("-- numb_neighb_subdom : %d\n", numb_neighb_subdom);
  iomrg::printf("-- numb_neighb_nodes  : %d\n", p_neighb2interfnode[numb_neighb_subdom]);
  iomrg::printf("-- list_neighb_subdom : \n");
  for ( int s = 0; s < numb_neighb_subdom; s++ ) {
    printf("  .. neighb %d : ", list_neighb_subdom[s]);
    for ( int i = p_neighb2interfnode[s]; i < p_neighb2interfnode[s+1]; i++ ) {
      printf("(%d x %d) ",
             neighb2interfnode[i], neighb2interfnode_multiplicity[i]);
    }
    printf("\n");
  }
  printf("\n");

  return 0;
}

void LocalMatrixToGlobalPositions(
        const MatrixDense<double,int>& A_local,
        MatrixDense<double,int>& A_global,
        int* l2g){
  int n = A_local.GetNumbRows();

  for(int i = 0;i < n;++i){
    for(int j = 0;j < n;++j){
      A_global(l2g[i], l2g[j]) = A_local(i,j);
    }
  }
}

________________________________________________________________________________

} // namespace DataTopology {
