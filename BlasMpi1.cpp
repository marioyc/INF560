/*!
*  @file BlasMpi.cpp
*  @internal source of class BlasMpi part 2
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Wed Dec 09 16:36:50 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages

// project packages
#include "BlasMpi1.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"

// third-party packages


//! @namespace BlasMpi
namespace BlasMpi {

________________________________________________________________________________

//! @internal compute dot product
double Dot (
        const Vector<double,int>& x_local,
        const Vector<double,int>& y_local,
        MPI_Comm& mpi_comm ) {

  // -- compute local dot product
  double local_dot = x_local.Dot( y_local );

  // -- compute global dot product
  double global_dot = 0.;
  // all reduce
  MPI_Allreduce( &local_dot, &global_dot, 1, MPI_DOUBLE, MPI_SUM, mpi_comm );

  return global_dot;
}

________________________________________________________________________________

//! @internal compute xx-norm (xx = 2 | inf | ...)
double NormLp (
        const Vector<double,int>& x_local,
        const int norm_type,
        MPI_Comm& mpi_comm ) {

  // -- norm, only sum, without (.)^(1/p)
  double local_normp = x_local.NormpLp( norm_type );

  // -- compute global norm
  double global_normp = 0.;

  // all reduce
  if ( norm_type == mathmrg::norm::c_LINF ) { // -- norm-Linfty
    MPI_Allreduce( &local_normp, &global_normp, 1,
                   MPI_DOUBLE, MPI_MAX, mpi_comm );
  } else { // -- norm-Lp
    MPI_Allreduce( &local_normp, &global_normp, 1,
                   MPI_DOUBLE, MPI_SUM, mpi_comm );
    if ( norm_type == mathmrg::norm::c_L1 || norm_type == 1 ) {
      // global_normp is already the L1 norm
    } else if ( norm_type == mathmrg::norm::c_L2 || norm_type == 2 ) {
      global_normp = std::sqrt( global_normp );
    } else if ( norm_type > 2 ) {
      global_normp = std::pow( global_normp, 1. / double(norm_type) );
    }
  }

  return global_normp;
}

________________________________________________________________________________

//! @internal compute xx-norm, without (.)^(1/p), (xx = 2 | inf | ...)
double NormpLp (
        const Vector<double,int>& x_local,
        const int norm_type,
        MPI_Comm& mpi_comm ) {

  // -- norm, only sum, without (.)^(1/p)
  double local_normp = x_local.NormpLp( norm_type );

  // -- compute global norm
  double global_normp = 0.;

  // all reduce
  if ( norm_type == mathmrg::norm::c_LINF ) { // -- norm-Linfty
    MPI_Allreduce( &local_normp, &global_normp, 1,
                   MPI_DOUBLE, MPI_MAX, mpi_comm );
  } else { // -- norm-Lp
    MPI_Allreduce( &local_normp, &global_normp, 1,
                   MPI_DOUBLE, MPI_SUM, mpi_comm );
  }

  return global_normp;
}

________________________________________________________________________________


} // namespace BlasMpi {
