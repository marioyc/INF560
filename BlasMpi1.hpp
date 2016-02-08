/*!
*  @file BlasMpi.hpp
*  @brief header of BlasMpi part 2
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Wed Dec 09 16:36:50 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_BLASMPI1_HPP_
#define GUARD_BLASMPI1_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace BlasMpi
namespace BlasMpi {

//! @brief compute dot product
//! @param [in] x_local = local first vector
//! @param [in] y_local = local second vector
//! @param [in] mpi_comm = MPI communicator
//! @return dot product
double Dot (
        const Vector<double,int>& x_local,
        const Vector<double,int>& y_local,
        MPI_Comm& mpi_comm ) ;

//! @brief compute xx-norm (xx = 2 | inf | ...)
//! @param [in] x_local = local vector
//! @param [in] norm_type = type of the norm (x)
//! @param [in] mpi_comm = MPI communicator
//! @return norm Lp
double NormLp (
        const Vector<double,int>& x_local,
        const int norm_type,
        MPI_Comm& mpi_comm ) ;

//! @brief compute xx-norm, without (.)^(1/p), (xx = 2 | inf | ...)
//! @param [in] x_local = local vector
//! @param [in] norm_type = type of the norm (x)
//! @param [in] mpi_comm = MPI communicator
//! @return norm pLp, without (.)^(1/p)
double NormpLp (
        const Vector<double,int>& x_local,
        const int norm_type,
        MPI_Comm& mpi_comm ) ;

} // namespace BlasMpi {

#endif // GUARD_BLASMPI1_HPP_
