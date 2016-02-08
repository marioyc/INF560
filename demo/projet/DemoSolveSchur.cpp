// basic packages
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// project packages
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DataTopology.hpp"
#include "Schur.hpp"

// third-party packages

int main (
        int argc,
        char** argv ) {

  // ---------------------------------------------------------------------------
  // -- read command line arguments
  // ---------------------------------------------------------------------------


  // -- argv command line arguments management
  if( argc < 2 ) {
    printf("Usage: %s", argv[0]);
    printf(" matrix_group_name matrix_group_type ");
    printf("[rhs_group_name rhs_group_type output_group_name");
    printf(" output_group_type image=false root=0]\n");
    exit(-1);
  }

  // ---------------------------------------------------------------------------
  // -- pre-processing
  // ---------------------------------------------------------------------------

  // -- read command line arguments
  int arg_num = 1;
  const char* argv_matrix_group_name = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_matrix_group_type = argc>arg_num?argv[arg_num++]:NULL;
  const char* argv_rhs_group_name = argc>arg_num?argv[arg_num++]
                                    :argv_matrix_group_name;
  const char* argv_rhs_group_type = argc>arg_num?argv[arg_num++]
                                    :argv_matrix_group_type;
  const char* argv_output_group_name = argc>arg_num?argv[arg_num++]
                                       :argv_matrix_group_name;
  const char* argv_output_group_type = argc>arg_num?argv[arg_num++]
                                       :argv_matrix_group_type;
  const bool argv_opt_csv_matrix = argc>arg_num?atoi(argv[arg_num++]):0;
  const int argv_proc_root = argc>arg_num?atoi(argv[arg_num++]):0;

  // 0: not otherwise yes
  bool is_read_rhs = (argv_rhs_group_name != NULL) &&
                      (argv_rhs_group_type != NULL);


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


  // sample: proc_numb = 3
  std::string proc_numb_str = "000000";
  // sample: ss_proc_numb = "3"
  std::stringstream ss_proc_numb;
  ss_proc_numb << proc_numb;
  // sample: numb_procs = 5
  std::stringstream ss_numb_procs;
  // sample: ss_numb_procs = "5"
  ss_numb_procs << numb_procs;
  // sample: l_ns = 1
  int l_ns = std::string(ss_proc_numb.str()).size();
  // sample: prefix_str = "5_000003"
  std::string prefix_str = ss_numb_procs.str() + "_" +
                           proc_numb_str.substr(l_ns) + ss_proc_numb.str();
  // sample: matrix_group_name = "brain_5_000003"
  std::string matrix_group_name = std::string(argv_matrix_group_name)
                                         + "_" + prefix_str;

  // sample: matrix_group_name = "brain_5_000003"
  std::string rhs_group_name = std::string(argv_rhs_group_name)
                                      + "_" + prefix_str;

  // sample: mat_filename = "brain_5_000003.csv"
  std::string mat_filename = matrix_group_name
                              + "." + std::string(argv_matrix_group_type);
  // sample: rhs_filename = "brain_5_000003_b.csv"
  std::string rhs_filename = rhs_group_name
                             + "." + std::string(argv_rhs_group_type);
  // sample: l2g_filename = "brain_5_000003.l2g"
  std::string l2g_filename = matrix_group_name + ".l2g";
  // sample: interface_filename = "brain_5_000003.itf"
  std::string interface_filename = matrix_group_name + ".itf";

  // sample: residual_filename = "brain_5_000003_"
  std::string output_matrix_group_name = std::string(argv_output_group_name)
                                 + "_" + prefix_str;
  // sample: solution_filename = "brain_5_000003_x.csv"
  std::string solution_filename  = std::string(argv_output_group_name)
                                   + "_x_" + prefix_str
                                   + "." + argv_output_group_type;
  // sample: solution_filename = "brain_5_000003_xg.csv"
  std::string gsolution_filename  = std::string(argv_output_group_name)
                                   + "_xg_" + prefix_str
                                   + "." + argv_output_group_type;

  // -- read (in) matrix: local
  MatrixDense<double,int> K_local;
  if ( argv_opt_csv_matrix == 0 ) {
    // -- read image matrix from file csv
    K_local.ReadFromFileCsv( mat_filename.c_str( ) );
  } else {
    // -- read csv matrix from file csv
    K_local.ReadImageMatrixFromFileCsv( mat_filename.c_str( ) );
  }


  // -- read (in) vector: local
  Vector<double,int> b_local;
  if( is_read_rhs ) {
     b_local.ReadFromFileCsv( rhs_filename.c_str( ) );
  } else {
    b_local.Allocate( K_local.GetNumbColumns() );
    b_local.Assign( 0, b_local.GetSize( )-1, 1.0 );
  }

  // -- read neighb2interfnode
  int subdom_numb = 0;
  int numb_subdom = 0;
  int numb_neighb_subdom = 0;
  int* list_neighb_subdom = NULL;
  int* p_neighb2interfnode = NULL;
  int* neighb2interfnode = NULL;
  int* neighb2interfnode_multiplicity = NULL;
  DataTopology::ReadNeighb2InterfaceNodeFromFile( subdom_numb, numb_subdom,
                                       numb_neighb_subdom, list_neighb_subdom,
                                       p_neighb2interfnode, neighb2interfnode,
                                       neighb2interfnode_multiplicity,
                                       interface_filename.c_str() );

  // -- write interfaces to stdout
  DataTopology::WriteNeighb2InterfaceNodeToStdout( subdom_numb, numb_subdom,
                                     numb_neighb_subdom, list_neighb_subdom,
                                     p_neighb2interfnode, neighb2interfnode,
                                     neighb2interfnode_multiplicity );

  // -- read l2g
  int numb_global_node = 0;
  int numb_l2g = 0;
  int* l2g = NULL;
  DataTopology::ReadL2gFromFile( numb_global_node, numb_l2g, l2g,
                                 l2g_filename.c_str() );
  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // ---------------------------------------------------------------------------
  // -- processing
  // ---------------------------------------------------------------------------

  // -- solve Ax = b_local, by Schur Complement method
  // synchronous communication, synchronous iteration

  // -- local solution
  Vector<double,int> x_local( b_local.GetSize( ) );

  // -- solve with parallel schur complement method
//  Schur::SolveSystem( x_local, K_local, b_local,
//                      numb_neighb_subdom, list_neighb_subdom,
//                      p_neighb2interfnode, neighb2interfnode, mpi_comm );

  // -- try to wait all processors
  MPI_Barrier( mpi_comm );

  // ---------------------------------------------------------------------------
  // -- post-processing
  // ---------------------------------------------------------------------------

  // -- plot local solution
  x_local.WriteToFileCsv( solution_filename.c_str(), '\n' );

  // -- right only for 2 subdomains
  Vector<double, int> x_global;
  x_global.Allocate( numb_global_node );
  x_global.Assign( 0, numb_global_node-1, 0.0 );
  for ( int n = 0; n < numb_l2g; n++ ) {
    x_global( l2g[n] ) = x_local( n );
  }
  x_global.WriteToFileCsv( gsolution_filename.c_str(), '\n' );


  delete [] l2g;
  delete [] list_neighb_subdom;
  delete [] p_neighb2interfnode;
  delete [] neighb2interfnode;
  delete [] neighb2interfnode_multiplicity;

  // ---------------------------------------------------------------------------
  // -- finalize MPI
  // ---------------------------------------------------------------------------

  // -- finalizes MPI
  MPI_Finalize( );

  return 0;
}