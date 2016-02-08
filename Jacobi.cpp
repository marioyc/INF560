/*!
*  @file Jacobi.cpp
*  @internal source of class Jacobi
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Wed Dec 09 16:36:50 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages
#include <time.h>

// project packages
#include "Jacobi.hpp"
#include "BlasMpi.hpp"
#include "BlasMpi1.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"

// third-party packages


#define GET_CLOCK( ) double(clock() / double(CLOCKS_PER_SEC))

//! @namespace Solver
namespace Solver {

________________________________________________________________________________

//! @brief solve A_local * x_local = b_local by Jacobi sequential algorithm
int SequentialJacobi (
        Vector<double,int>& x,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A,
        const Vector<double,int>& b,
        const int max_numb_iter,
        const double residual_threshold ) {

  // local size
  int size = A.GetNumbRows( );

  // old solution
  Vector<double,int> x_old( size );

  // x_local_old(:) = 0;
  x_old.Assign( 0, size-1, 0. );

  // LINF norm
  double error_x = 0.;
  double error_b = 0.;

  // -- compute LINF norm of b
  error_b = 0.;
  for( int i = 0; i < size; i++ ) {
    error_b = std::max( error_b, std::abs( b(i) ) );
  }

  // iteration parameters
  int iter_numb = 0;

  // -- loop until convergence
  while( iter_numb < max_numb_iter ) {

    // -- Jacobi algorithm
    // $x^{(k+1)}=D^{-1}(E+F) x^{(k)}+D^{-1}b$
    // $x^{(k+1)}_i= -\frac{1}{a_{ii}} \sum_{j=1 \atop j \ne i}^n a_{ij}x^{(k)}_j + \frac{b_i}{a_{ii}}$
    for( int i = 0; i < size; i++ ) {
      double temp = 0.;
      for( int j = 0; j < size; j++ ) {
        if ( j != i ) {
          temp = temp + A(i,j) * x_old(j);
        }
      }
      x(i) = (1.0 / A(i, i)) * (b(i) - temp);
    }

    // -- compute LINF norm of x
    error_x = 0.;
    for( int i = 0; i < size; i++ ) {
      error_x = std::max( error_x, std::abs( x(i) - x_old(i) ) );
   }

    // -- save current solution
    for( int i = 0; i < size; i++ ) {
      x_old(i) = x(i);
    }

    // -- save residual
    residual[iter_numb] = error_x; // / error_b;
    temporal[iter_numb] = GET_CLOCK();

    // -- check convergence
    if( (iter_numb >= max_numb_iter) /* ||
        (residual[iter_numb] <= residual_threshold) */ ) {
      iter_numb = iter_numb + 1;
      break;
    }

    // -- increment iterations number
    iter_numb = iter_numb + 1;
  }

  // -- set final number of iterations
  numb_iter = iter_numb;

  return 0;
}

________________________________________________________________________________

//! @internal solve Ax = b by Jacobi synchronous iteration, synchronous comm
int JacobiBandRowSyncISyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residual_threshold,
        MPI_Comm& mpi_comm ) {

  // number of processors
  int numb_procs;
  // process number
  int proc_numb;

  // get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // get current process number
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // local size
  int size_local = A_local.GetNumbRows( );
  // global size
  int size_global = A_local.GetNumbColumns( );

  // delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  DataTopology::BandTopology( band_list_start, band_list_size,
                              size_global, numb_procs );
  // band position
  int band_pos = band_list_start[proc_numb];
  // tag exchange value
  const int TAG_EXCHANGE_VALUE = 301;

  // global solution
  Vector<double,int> x_global_old( size_global );
  // x_global_old = 0;
  x_global_old.Assign( 0, size_global-1, 0. );

  // LINF norm
  double error_x_global = 0.;
  double error_b_global = 0.;
  double error_x_local = 0.;
  double error_b_local = 0.;

  // -- compute LINF norm of local b
  error_b_local = 0.;
  for( int i = 0; i < size_local; i++ ) {
    error_b_local = std::max( error_b_local, std::abs( b_local(i) ) );
  }

  // MPI all reduce
  MPI_Allreduce( &error_b_local, &error_b_global, 1, MPI_DOUBLE, MPI_MAX, mpi_comm );

  // iteration parameters
  int iter_numb = 0;

  // -- loop until convergence
  while( iter_numb < max_numb_iter ) {

    // -- Jacobi algorithm
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      double temp = 0.;
      for( int j = 0; j < size_global; j++ ) {
        if ( j != l2g_i ) {
          temp = temp + A_local(i,j) * x_global_old(j);
        }
      }
      x_local(i) = (1.0 / A_local(i, l2g_i)) * (b_local(i) - temp);
    }

    // -- compute LINF norm of local x
    error_x_local = 0.;
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      error_x_local = std::max( error_x_local, std::abs( x_local(i) - x_global_old(l2g_i) ) );
    }

    // -- save solution
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      x_global_old(l2g_i) = x_local(i);
    }

    // -- MPI send
    for( int k = 0; k < numb_procs; k++ ) {
      if ( k != proc_numb ) {
        MPI_Send( x_local.GetCoef( ), size_local,
                  MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm );
      }
    }

    // MPI receive
    for( int k = 0; k < numb_procs; k++ ) {
      if ( k != proc_numb ) {
         MPI_Status recv_status;
         MPI_Recv( x_global_old.GetCoef( ) + band_list_start[k], band_list_size[k],
                   MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &recv_status);
      }
    }

    // MPI barrier
    MPI_Barrier( mpi_comm );

    // MPI all reduce
    MPI_Allreduce( &error_x_local, &error_x_global, 1, MPI_DOUBLE, MPI_MAX, mpi_comm );

    // -- save residual
    residual[iter_numb] = error_x_global / error_b_global;
    temporal[iter_numb] = GET_CLOCK();

    // -- check convergence
    if( (iter_numb >= max_numb_iter) /*||
        (residual[iter_numb] <= residual_threshold) */ ) {
      iter_numb = iter_numb + 1;
      break;
    }

    // -- increment the number of iterations
    iter_numb = iter_numb + 1;
  }

  // MPI barrier
  MPI_Barrier( mpi_comm );

  // -- set final number of iteration
  numb_iter = iter_numb;

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  return 0;
}

________________________________________________________________________________

//! @internal solve Ax = b_local by Jacobi asynchronous algorithm
int JacobiBandRowSyncIAsyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residual_threshold,
        MPI_Comm& mpi_comm ) {


  // -- number of processors
  int numb_procs;
  // -- process number (process proc_numb)
  int proc_numb;
  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process proc_numb
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // -- local size
  int size = A_local.GetNumbRows( );
  // -- local matrix: number of columns
  int size_global = A_local.GetNumbColumns( );

  // -- delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  DataTopology::BandTopology( band_list_start, band_list_size,
                              size_global, numb_procs );

  // -- band position
  int band_pos = band_list_start[proc_numb];
  // -- tag exchange value
  const int TAG_EXCHANGE_VALUE = 301;

  // -- statuses
  MPI_Status* status = NULL;
  status = new MPI_Status[numb_procs-1];
  // -- receive request
  MPI_Request* recv_request = NULL;
  recv_request = new MPI_Request[numb_procs-1];
  // -- send request
  MPI_Request* send_request = NULL;
  send_request = new MPI_Request[numb_procs-1];

  // -- global solution
  Vector<double,int> x_global_old( size_global );
  // x_global_old = 0;
  x_global_old.Assign( 0, size_global-1, 0. );

  // -- local error
  double error = 0.;
  // -- LINF norm
  double error_max = 0.;

  // -- iteration parameters
  int iter_numb = 0;
  // -- loop until convergence
  while( iter_numb < max_numb_iter ) {

    // -- Jacobi algorithm
    for( int i = 0; i < size; i++ ) {
      int i_global = i + band_pos;
      double temp = 0.;
      for( int j = 0; j < size_global; j++ ) {
        if ( j != i_global ) {
          temp = temp + A_local(i,j) * x_global_old(j);
        }
      }
      x_local(i) = (1.0 / A_local(i, i_global)) * (b_local(i) - temp);
    }

    // -- compute local error and save solution
    error = 0.;
    for( int i = 0; i < size; i++ ) {
      int i_global = i + band_pos;
      error = std::max( error, std::abs( x_local(i) - x_global_old(i_global) ) );
      x_global_old(i_global) = x_local(i);
    }

    // -- sending
    int s = 0;
    for( int k = 0; k < numb_procs; k++ ) {
    int sent = 1;
      if ( k != proc_numb ) {
        // -- send data value
        MPI_Isend( x_local.GetCoef( ), size,
                  MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &send_request[s] );
        s++;
      }
    }

    // -- receiving
    int r = 0;
    for( int k = 0; k < numb_procs; k++ ) {
      if ( k != proc_numb ) {
         // -- receive data value
         MPI_Request rr;
         MPI_Irecv( x_global_old.GetCoef( ) + k * size, band_list_size[k],
                    MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &recv_request[r]);
        r++;
      }
    }

    // -- wait until all messages have been received
    MPI_Waitall(numb_procs-1, &recv_request[0], &status[0]);

    // all reduce
    MPI_Allreduce( &error, &error_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm );

    // -- save residual
    residual[iter_numb] = error_max;
    temporal[iter_numb] = GET_CLOCK();
    // -- check convergence
    if( (iter_numb >= max_numb_iter) ||
        (residual[iter_numb] <= residual_threshold) ) {
      iter_numb = iter_numb + 1;
      break;
    }

    if( (iter_numb >= max_numb_iter) ) {
      break;
    }

    // -- increment the number of iterations
    iter_numb = iter_numb + 1;
  }

  // -- set final number of iteration
  numb_iter = iter_numb;

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // -- delete request
  delete [] status;
  delete [] recv_request;
  delete [] send_request;

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  return 0;
}

________________________________________________________________________________

//! @internal solve Ax = b by Jacobi asynchronous iteration, asynchronous comm
int JacobiBandRowAsyncIAsyncC (
        Vector<double,int>& x_local,
        int& numb_iter,
        Vector<double,int>& residual,
        Vector<double,int>& temporal,
        const MatrixDense<double,int>& A_local,
        const Vector<double,int>& b_local,
        const int max_numb_iter,
        const double residual_threshold,
        MPI_Comm& mpi_comm ) {

  // number of processors
  int numb_procs;
  // process number
  int proc_numb;
  // get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // get current process number
  MPI_Comm_rank( mpi_comm, &proc_numb );

  // local size
  int size_local = A_local.GetNumbRows( );
  // global size
  int size_global = A_local.GetNumbColumns( );

  // delete band list position and size
  int* band_list_start = NULL;
  int* band_list_size = NULL;
  DataTopology::BandTopology( band_list_start, band_list_size,
                              size_global, numb_procs );
  // band position
  int band_pos = band_list_start[proc_numb];
  // tag exchange value
  const int TAG_EXCHANGE_VALUE = 301;

  // statuses
  MPI_Status* status = NULL;
  status = new MPI_Status[numb_procs];
  // receive request
  MPI_Request* recv_request = NULL;
  recv_request = new MPI_Request[numb_procs];
  // send request
  MPI_Request* send_request = NULL;
  send_request = new MPI_Request[numb_procs];

  // global solution
  Vector<double,int> x_global_old( size_global );
  // x_global_old = 0;
  x_global_old.Assign( 0, size_global-1, 0. );


  // -- LINF norm
  double error = 0.;
  double error_x_local = 0.;

   // -- first communications
  // MPI send
  for( int k = 0; k < numb_procs; k++ ) {
    if ( k != proc_numb ) {
      MPI_Isend( x_local.GetCoef( ), size_local,
                MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &send_request[k] );
    }
  }

  // MPI receive
  for( int k = 0; k < numb_procs; k++ ) {
    if ( k != proc_numb ) {
      MPI_Irecv( x_global_old.GetCoef( ) + band_list_start[k], band_list_size[k],
                 MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &recv_request[k]);
    }
  }

  // iteration parameters
  int iter_numb = 0;

  // -- loop until convergence
  while( iter_numb < max_numb_iter ) {

    // MPI receive
    for( int k = 0; k < numb_procs; k++ ) {
      if ( k != proc_numb ) {
        int is_received = 0;
        MPI_Status r_status;
        MPI_Test( &recv_request[k], &is_received, &r_status);
        if ( is_received == 1 ) {
          MPI_Irecv( x_global_old.GetCoef( ) + band_list_start[k], band_list_size[k],
                     MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &recv_request[k]);
        }
      }
    }

    // -- Jacobi algorithm
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      double temp = 0.;
      for( int j = 0; j < size_global; j++ ) {
        if ( j != l2g_i ) {
          temp = temp + A_local(i,j) * x_global_old(j);
        }
      }
      x_local(i) = (1.0 / A_local(i, l2g_i)) * (b_local(i) - temp);
    }

    // -- compute LINF norm of local x
    error_x_local = 0.;
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      error_x_local = std::max( error_x_local, std::abs( x_local(i) - x_global_old(l2g_i) ) );
    }

    // -- computation of error_x_global is missing but too complicated
    // -- computation of error_x_global done periodically

    // -- save solution
    for( int i = 0; i < size_local; i++ ) {
      int l2g_i = i + band_pos;
      x_global_old(l2g_i) = x_local(i);
    }

    // MPI send
    for( int k = 0; k < numb_procs; k++ ) {
      if ( k != proc_numb ) {
        int is_sent = 0;
        MPI_Status s_status;
        MPI_Test( &send_request[k], &is_sent, &s_status);
        if ( is_sent == 1 ) {
          MPI_Isend( x_local.GetCoef( ), size_local,
                    MPI_DOUBLE, k, TAG_EXCHANGE_VALUE, mpi_comm, &send_request[k] );
        }
      }
    }
    // -- save residual
    residual[iter_numb] = error_x_local;
    temporal[iter_numb] = GET_CLOCK();

    // -- check convergence
    if( (iter_numb >= max_numb_iter) ) {
      iter_numb = iter_numb + 1;
      break;
    }

    // -- increment the number of iterations
    iter_numb = iter_numb + 1;
  }

  // -- set final number of iteration
  numb_iter = iter_numb;

  // MPI barrier
  MPI_Barrier( mpi_comm );

  // -- delete request
  delete [] status;
  delete [] recv_request;
  delete [] send_request;

  // -- delete band list position and size
  delete [] band_list_start;
  delete [] band_list_size;

  return 0;
}

________________________________________________________________________________


} // namespace Solver {
