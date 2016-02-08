/*!
*  @file BlasMpi.cpp
*  @internal source of class BlasMpi
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages

// project packages
#include "BlasMpi.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"

// third-party packages


//! @namespace BlasMpi
namespace BlasMpi {

________________________________________________________________________________

//! @internal compute matrix-vector product
int MatrixVectorProductBandRow (
        Vector<double,int>& y_local,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        MPI_Comm& mpi_comm ) {

  // -- number of processors
  int numb_procs;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );

  // -- size of the vector
  int size_local = x_local.GetSize( );

  // -- get global size (band-row splitting)
  int size_global = 0;
  MPI_Allreduce( &size_local, &size_global, 1, MPI_INT, MPI_SUM, mpi_comm );

  // -- global vector
  Vector<double,int> x_global( size_global );

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  DataTopology::BandTopology( band_list_start, band_list_size,
                              size_global, numb_procs );

  // MPI_Allgatherv - MPICH
  // MPI_Allgatherv
  // Gathers data from all tasks and deliver the combined data to all tasks
  // Synopsis
  // int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
  //                    void *recvbuf, const int *recvcounts, const int *displs,
  //                    MPI_Datatype recvtype, MPI_Comm comm)
  // Input Parameters
  // sendbuf
  //     starting address of send buffer (choice)
  // sendcount
  //     number of elements in send buffer (integer)
  // sendtype
  //     data type of send buffer elements (handle)
  // recvcounts
  //     integer array (of length group size) containing the number of elements that are to be received from each process
  // displs
  //     integer array (of length group size). Entry i specifies the displacement (relative to recvbuf ) at which to place the incoming data from process i
  // recvtype
  //     data type of receive buffer elements (handle)
  // comm
  //     communicator (handle)
  //
  // Output Parameters
  // recvbuf
  //     address of receive buffer (choice)
  // -- assemble local vector to 'root', upon processors
  MPI_Allgatherv( x_local.GetCoef( ), size_local, MPI_DOUBLE,
                  x_global.GetCoef( ), band_list_size, band_list_start,
                  MPI_DOUBLE, mpi_comm );

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  // -- perform local matrix vector product: A_local * x_global
  A_local.MatrixVectorProduct( y_local, x_global );

  return 0;
}

________________________________________________________________________________

//! @internal compute matrix-vector product
int MatrixVectorProductBandColumn (
        Vector<double,int>& y_global,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm ) {

  // -- local matrix: number of rows
  int numb_rows_global    = A_local.GetNumbRows( );

  // -- global vector
  Vector<double,int> y_global_tmp( numb_rows_global );

  // -- perform local matrix vector product: A_local * x_local
  A_local.MatrixVectorProduct( y_global_tmp, x_local );

  // -- allocate global result
  y_global.Allocate( numb_rows_global );

  // reduce to root proc
  MPI_Reduce( y_global_tmp.GetCoef( ), y_global.GetCoef( ),
              numb_rows_global, MPI_DOUBLE, MPI_SUM, root, mpi_comm );

  return 0;
}

________________________________________________________________________________

//! @internal compute matrix-vector product y_local:= A_local *x_local
int MatrixVectorProductBlock (
        Vector<double,int>& y_local,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) {

  // -- local matrix: number of rows
  int numb_rows = A_local.GetNumbRows( );

  // -- tmp vector
  Vector<double,int> y_local_tmp( numb_rows );

  // -- perform local subdomain matrix vector product: A_local * x_local_global
  A_local.MatrixVectorProduct( y_local_tmp, x_local );

  // -- every row reduces its share of y_local
  MPI_Reduce( y_local_tmp.GetCoef( ), y_local.GetCoef( ),
              numb_rows, MPI_DOUBLE, MPI_SUM, root, mpi_comm_rows );

  return 0;
}

________________________________________________________________________________

//! @internal compute matrix-matrix product C := A_local * B
//! @note number of processor:  P = P_I * P_J
//!  page 60, livre "Calcul Scientifique Parallèle"
//! for K = 1 to P_j
//!  if J = K then
//!    B_Temp = B
//!  end if
//!  MPI_Bcast(B_temp, longueur, type, K-1, mpi_comm_rows)
//!  if I = K then
//!    C_Temp = C
//!  end if
//!  MPI_Bcast(C_temp, longueur, type, K-1, mpi_comm_columns)
//!  A_local = A_local + B_temp * C_temp
//!
int MatrixMatrixProductBlock (
        MatrixDense<double,int>& C,
        const MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& B,
        const int root,
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

  // -- local matrix: number of rows
  int A_local_numb_rows    = A_local.GetNumbRows( );
  // -- local matrix: number of columns
  int A_local_numb_columns = A_local.GetNumbColumns( );

  // -- local matrix: number of rows
  int B_numb_rows    = B.GetNumbRows( );
  // -- local matrix: number of columns
  int B_numb_columns = B.GetNumbColumns( );

  // -- result local matrix: C
  C.Allocate( A_local_numb_rows, B_numb_columns );

  // **************** optimize, declaration A_local et B TODO
  // -- temporary local matrix: A_local
  MatrixDense<double,int> A_local_temp;
  A_local_temp.Allocate( A_local_numb_rows, A_local_numb_columns );

  // -- temporary local matrix: B
  MatrixDense<double,int> B_temp;
  B_temp.Allocate( B_numb_rows, B_numb_columns );

  // -- initialize C to 0
  C.Initialize( 0 );

  // --
  for ( int k = 0; k < numb_procs_j; k++ ) {

    // -- A_local block communication
    if ( proc_numb_i == k ) {
      A_local_temp = A_local;
    }
    MPI_Bcast( A_local_temp.GetCoef( ), A_local_numb_rows*A_local_numb_columns,
               MPI_DOUBLE, k, mpi_comm_rows );

    // -- B block communication
    if ( proc_numb_j == k ) {
      B_temp = B;
    }
    MPI_Bcast( B_temp.GetCoef( ), B_numb_rows*B_numb_columns,
               MPI_DOUBLE, k, mpi_comm_columns );

    // -- local block multiplication: C = C + A_local_temp * B_temp
    A_local_temp.MatrixMatrixProductAdd( C, B_temp );

  }

  return 0;
}

________________________________________________________________________________

} // namespace BlasMpi {
