/*!
*  @file DataTopology.hpp
*  @brief manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Frédéric Magoulès, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

#ifndef GUARD_DATATOPOLOGY_HPP_
#define GUARD_DATATOPOLOGY_HPP_

// basic packages
#include <mpi.h>

// project packages
#include "dllmrg.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"

// third-party packages


//! @namespace DataTopology
namespace DataTopology {

// -----------------------------------------------------------------------------
// -- topology distribution: BAND-ROW or BAND-COLUMN
// -----------------------------------------------------------------------------

//! @brief index position of a given band
//! @param [in] proc_numb = processor number
//! @param [in] numb_procs = number of processors
//! @param [in] size = global size (= numb_columns or numb_rows)
//! @return index position of a given band
int BandPos (
        const int proc_numb,
        const int numb_procs,
        const int size ) ;

//! @brief size of a given band
//! @param [in] proc_numb = processor number
//! @param [in] numb_procs = number of processors
//! @param [in] size = global size (= numb_columns or numb_rows)
//! @return size of a given band
int BandSize (
        const int proc_numb,
        const int numb_procs,
        const int size ) ;

//! @brief band list start and position
//! @param [in,out] band_list_start = band list start position
//! @param [in,out] band_list_size = band list size
//! @param [in] size = global size (= numb_columns or numb_rows)
//! @param [in] numb_procs = number of processors
//! @param [in] shift = shift start and shift size (= numb_rows or numb_columns)
//! @return error code
int BandTopology (
        int*& band_list_start,
        int*& band_list_size,
        const int size,
        const int numb_procs,
        const int shift = 1 ) ;

// -----------------------------------------------------------------------------
// -- topology distribution: processors
// -----------------------------------------------------------------------------

//! @brief creation of two-dimensional grid communicator and communicators
//!         for each row and each column of the grid
//! @param [in,out] mpi_comm_rows = grid rows communicator
//! @param [in,out] mpi_comm_columns = grid columns communicator
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int GridCartesianComm (
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns,
        MPI_Comm& mpi_comm ) ;

// -----------------------------------------------------------------------------
// -- Vector: BAND-ROW or BAND-COLUMN
// -----------------------------------------------------------------------------

//! @brief distribute vector upon processors (band)
//! @param [in,out] x_local = local vector
//! @param [in] x_global = global vector
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int DistributeVectorBand (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief assemble vector upon processors (band)
//! @param [in,out] x_global = global vector
//! @param [in] x_local = local vector
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int AssembleVectorBand (
        Vector<double,int>& x_global,
        const Vector<double,int>& x_local,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief build 'band_numb'-th band of a vector (band)
//! @param [in,out] x_local = local vector
//! @param [in] x_global = global vector
//! @param [in] numb_procs = number of processors
//! @param [in] proc_numb = processor number
//! @return error code
int BuildVectorBand (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int numb_procs,
        const int proc_numb ) ;

// -----------------------------------------------------------------------------
// -- Matrix: BAND-ROW
// -----------------------------------------------------------------------------

//! @brief distribute matrix upon processors (band row)
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int DistributeMatrixBandRow (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief assemble matrix upon processors (band row)
//! @param [in,out] A_global = global matrix
//! @param [in] A_local = local matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int AssembleMatrixBandRow (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief build 'band_numb'-th band (row) of a matrix
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] numb_procs = number of processors
//! @param [in] proc_numb = processor number
//! @return error code
int BuildMatrixBandRow (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs,
        const int proc_numb ) ;

// -----------------------------------------------------------------------------
// -- Matrix: BAND-COLUMN
// -----------------------------------------------------------------------------

//! @brief distribute matrix upon processors (band column)
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int DistributeMatrixBandColumn (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief assemble matrix upon processors (band column)
//! @param [in,out] A_global = global matrix
//! @param [in] A_local = local matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm = MPI communicator
//! @return error code
int AssembleMatrixBandColumn (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        const int root,
        MPI_Comm& mpi_comm ) ;

//! @brief build 'band_numb'-th band (column) of a matrix
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] numb_procs = number of processors
//! @param [in] proc_numb = processor number
int BuildMatrixBandColumn (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs,
        const int proc_numb ) ;

// -----------------------------------------------------------------------------
// -- Vector: BLOCK
// -----------------------------------------------------------------------------

// NOT VALIDATED
//! @brief distribute vector upon processors (block)
//! @param [in,out] x_local = local vector
//! @param [in] x_global = global vector
//! @param [in] root = root processor
//! @param [in] root = root processor
//! @param [in] mpi_comm_rows = grid rows communicator
//! @param [in] mpi_comm_columns = grid columns communicator
//! @return error code
int DistributeVectorBlock (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) ;

// NOT VALIDATED
//! @brief build the 'block_numb'-th band of a vector (block)
//! @param [in,out] x_local = local vector
//! @param [in] x_global = global vector
//! @param [in] numb_procs_i = number of processors (i-)
//! @param [in] numb_procs_j = number of processors (j-)
//! @param [in] proc_numb_i = processor number (i-)
//! @param [in] proc_numb_j = processor number (j-)
//! @return error code
int BuildVectorBlock (
        Vector<double,int>& x_local,
        const Vector<double,int>& x_global,
        const int numb_procs_i,
        const int numb_procs_j,
        const int proc_numb_i,
        const int proc_numb_j ) ;

// -----------------------------------------------------------------------------
// -- Matrix: BLOCK
// -----------------------------------------------------------------------------

// NOT VALIDATED
//! @brief distribute matrix upon processors (block)
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm_rows = grid rows communicator
//! @param [in] mpi_comm_columns = grid columns communicator
//! @remarks checkerboard matrix decomposition
//! @return error code
int DistributeMatrixBlock (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) ;

// NOT VALIDATED
//! @brief assemble matrix upon processors (block)
//! @param [in,out] A_global = global matrix
//! @param [in] A_local = local matrix
//! @param [in] root = root processor
//! @param [in] mpi_comm_rows = grid rows communicator
//! @param [in] mpi_comm_columns = grid columns communicator
//! @remarks checkerboard matrix decomposition
//! @return error code
int AssembleMatrixBlock (
        MatrixDense<double,int>& A_global,
        const MatrixDense<double,int>& A_local,
        int root,
        MPI_Comm& mpi_comm_rows,
        MPI_Comm& mpi_comm_columns ) ;

// NOT VALIDATED
//! @brief build the 'band_numb'-th block of a matrix
//! @param [in,out] A_local = local matrix
//! @param [in] A_global = global matrix
//! @param [in] numb_procs_i = number of processors (i-)
//! @param [in] numb_procs_j = number of processors (j-)
//! @param [in] proc_numb_i = processor number (i-)
//! @param [in] proc_numb_j = processor number (j-)
//! @remarks proc_numb follows row major order
//! @remarks proc_numb = proc_numb_i * numb_procs_j + proc_numb_j
//! @note number of processors must be perfect square
//! @return error code
int BuildMatrixMatrixBlock (
        MatrixDense<double,int>& A_local,
        const MatrixDense<double,int>& A_global,
        const int numb_procs_i,
        const int numb_procs_j,
        const int proc_numb_i,
        const int proc_numb_j ) ;

// -----------------------------------------------------------------------------
// -- Input/Output
// -----------------------------------------------------------------------------

//! @brief get filename of given processor among others
//! @param [in] file_group_name = file groupname
//! @param [in] file_group_type = file grouptype
//! @param [in] proc_numb = processor number
//! @param [in] numb_procs_i = number of processors (i-)
//! @param [in] numb_procs_j = number of processors (j-)
//! @return error code
const char* GetProcFilename (
        const char* file_group_name,
        const char* file_group_type,
        const int proc_numb,
        const int numb_procs_i,
        const int numb_procs_j = 0 ) ;

//! @brief read local vector
//! @param [in,out] x_local = local matrix
//! @param [in] file_group_name = file groupname
//! @param [in] file_group_type = file grouptype
//! @param [in] numb_procs = number of processors
//! @param [in] proc_numb = processor number
//! @return error code
int ReadLocalVectorFromFile (
        Vector<double,int>& A_local,
        const char* file_group_name,
        const char* file_group_type,
        const int numb_procs,
        const int proc_numb ) ;

//! @brief read local matrix
//! @param [in,out] A_local = local matrix
//! @param [in] file_group_name = file groupname
//! @param [in] file_group_type = file grouptype
//! @param [in] numb_procs = number of processors
//! @param [in] proc_numb = processor number
//! @param [in] opt_image_matrix = matrix csv or matrix image
//! @return error code
int ReadLocalMatrixFromFile (
        MatrixDense<double,int>& A_local,
        const char* file_group_name,
        const char* file_group_type,
        const int numb_procs,
        const int proc_numb,
        const bool opt_image_matrix = false ) ;

//! @brief read l2g from file
//! @param [in] numb_global = size of the problem
//! @param [in] numb_l2g = size of l2g
//! @param [in] l2g = local to global pointer
//! @param [in] file_name = file name
//! @return error code
int ReadL2gFromFile (
        int& numb_global,
        int& numb_l2g,
        int*& l2g,
        const char* file_name ) ;

//! @brief read interfaces from file
//! @param [in] subdom_numb = subdom number
//! @param [in] numb_subdom = number of subdomain
//! @param [in] numb_neighb_subdom = number of neighboring subdomain
//! @param [in] list_neighb_subdom = list of neighboring subdomain
//! @param [in] p_neighb2interfnode = pointer to neighbor to interface node
//! @param [in] neighb2interfnode = neighbor to interface node
//! @param [in] neighb2interfnode_weight = neighbor to interface node weight
//! @param [in] file_name = file name
//! @return error code
int ReadNeighb2InterfaceNodeFromFile (
        int& subdom_numb,
        int& numb_subdom,
        int& numb_neighb_subdom,
        int*& list_neighb_subdom,
        int*& p_neighb2interfnode,
        int*& neighb2interfnode,
        int*& neighb2interfnode_multiplicity,
        const char* file_name ) ;

//! @brief write interfaces to stdout
//! @param [in] subdom_numb = subdom number
//! @param [in] numb_subdom = number of subdomain
//! @param [in] numb_neighb_subdom = number of neighboring subdomain
//! @param [in] list_neighb_subdom = list of neighboring subdomain
//! @param [in] p_neighb2interfnode = pointer to neighbor to interface node
//! @param [in] neighb2interfnode = neighbor to interface node
//! @param [in] neighb2interfnode_weight = neighbor to interface node weight
//! @return error code
int WriteNeighb2InterfaceNodeToStdout (
        const int subdom_numb,
        const int numb_subdom,
        const int numb_neighb_subdom,
        const int* list_neighb_subdom,
        const int* p_neighb2interfnode,
        const int* neighb2interfnode,
        const int* neighb2interfnode_multiplicity ) ;

void LocalMatrixToGlobalPositions(
        const MatrixDense<double,int>& A_local,
        MatrixDense<double,int>& A_global,
        int* l2g);

} // namespace DataTopology {


#endif // GUARD_DATATOPOLOGY_HPP_
