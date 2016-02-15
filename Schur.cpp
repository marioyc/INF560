/*!
*  @file Schur.cpp
*  @internal manage data distribution
*  @author Abal-Kassim Cheik Ahamed, Fr�d�ric Magoul�s, Sonia Toubaline
*  @date Tue Nov 24 16:16:48 CET 2015
*  @version 1.0
*  @remarks
*/

// basic packages
#include <algorithm>
#include <vector>

// project packages
#include "Schur.hpp"
#include "Vector.hpp"
#include "MatrixDense.hpp"
#include "DirectSolver.hpp"

// third-party packages


//! @namespace Schur
namespace Schur {

________________________________________________________________________________

// -----------------------------------------------------------------------------
// -- topology distribution: processors
// -----------------------------------------------------------------------------


________________________________________________________________________________

// -----------------------------------------------------------------------------
// --
// -----------------------------------------------------------------------------

//! @internal block splitting of a matrix
int SplitMatrixToBlock (
        MatrixDense<double,int>& Aii,
        MatrixDense<double,int>& Aip,
        MatrixDense<double,int>& Api,
        MatrixDense<double,int>& App,
        const MatrixDense<double,int>& A,
        const int* list_node_i,
        const int numb_node_i,
        const int* list_node_p,
        const int numb_node_p ) {
  Aii.Allocate(numb_node_i,numb_node_i);

  for(int i = 0;i < numb_node_i;++i)
    for(int j = 0;j < numb_node_i;++j)
      Aii(i,j) = A(list_node_i[i],list_node_i[j]);

  App.Allocate(numb_node_p,numb_node_p);

  for(int i = 0;i < numb_node_p;++i)
    for(int j = 0;j < numb_node_p;++j)
      App(i,j) = A(list_node_p[i],list_node_p[j]);

  Aip.Allocate(numb_node_i,numb_node_p);

  for(int i = 0;i < numb_node_i;++i)
    for(int j = 0;j < numb_node_p;++j)
      Aip(i,j) = A(list_node_i[i],list_node_p[j]);

  Api.Allocate(numb_node_p,numb_node_i);

  for(int i = 0;i < numb_node_p;++i)
    for(int j = 0;j < numb_node_i;++j)
      Api(i,j) = A(list_node_p[i],list_node_i[j]);

  return 0;
}

void SolveSystem (
        Vector<double, int>& x_local,
        const MatrixDense<double, int>& K_local,
        const Vector<double, int>& b_local,
        const int numb_neighb_subdom,
        const int* list_neighb_subdom,
        const int* p_neighb2interfnode,
        const int* neighb2interfnode,
        const int numb_node_i,
        const int* list_node_i,
        const int* l2i,
        int numb_node_p,
        const int* list_node_p,
        const int* l2p,
        const int numb_global_node,
        const int numb_l2g,
        const int* l2g,
        const MPI_Comm& mpi_comm ) {
  // -- number of processors
  int numb_procs;
  // -- process number (process rank)
  int proc_numb;

  // -- get number of processes
  MPI_Comm_size( mpi_comm, &numb_procs );
  // -- get current process rank
  MPI_Comm_rank( mpi_comm, &proc_numb );


  MatrixDense<double,int> Kii;
  MatrixDense<double,int> Kip;
  MatrixDense<double,int> Kpi;
  MatrixDense<double,int> Kpp;

  Schur::SplitMatrixToBlock(Kii, Kip, Kpi, Kpp, K_local,
                            list_node_i, numb_node_i,
                            list_node_p, numb_node_p);

  // check transpose
  /*double maxd = 0;

  for(int i = 0;i < numb_node_i;++i){
    for(int j = 0;j < numb_node_p;++j){
      double d = Kip(i,j) - Kpi(j,i);
      if(d < 0) d = -d;
      maxd = std::max(maxd,d);
    }
  }

  iomrg::printf("max difference = %.10f\n", maxd);*/

  // Reconstruct K
  /*MatrixDense<double, int> K_global_partial,K_global_total;
  K_global_partial.Allocate(numb_global_node, numb_global_node);
  K_global_total.Allocate(numb_global_node, numb_global_node);
  K_global_partial.Initialize(0);

  DataTopology::LocalMatrixToGlobalPositions(K_local, K_global_partial, l2g);

  MPI_Reduce(K_global_partial.GetCoef(), K_global_total.GetCoef(), numb_global_node * numb_global_node,
            MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  if(proc_numb == 0){
    K_global_total.WriteToFileCsv(gmatrix_filename.c_str());
  }*/

  MatrixDense<double, int> Kii_lu;
  Factor::LU(Kii_lu, Kii);
  //Kii_lu.WriteToFileCsv( solution_filename.c_str());

  MatrixDense<double, int> Uii_inv,Lii_inv;
  Lii_inv.Allocate(numb_node_i, numb_node_i);
  Uii_inv.Allocate(numb_node_i, numb_node_i);

  Vector<double, int> x,rhs;
  rhs.Allocate(numb_node_i);

  for(int i = 0;i < numb_node_i;++i)
    rhs(i) = 0;

  // invert Lii and Uii
  for(int i = 0;i < numb_node_i;++i){
    if(i > 0) rhs(i - 1) = 0;
    rhs(i) = 1;
    DirectSolver::Forward(x, Kii_lu, rhs);

    for(int j = 0;j < numb_node_i;++j)
      Lii_inv(j,i) = x(j);

    DirectSolver::Backward(x,Kii_lu, rhs);

    for(int j = 0;j < numb_node_i;++j)
      Uii_inv(j,i) = x(j);
  }

  // calculate S_local
  MatrixDense<double, int> Lpi,Uip,prod;
  //Uii_inv.WriteToFileCsv( solution_filename.c_str());
  Kpi.MatrixMatrixProduct(Lpi, Uii_inv);
  Lii_inv.MatrixMatrixProduct(Uip, Kip);
  //Lii_inv.WriteToFileCsv( solution_filename.c_str());
  Lpi.MatrixMatrixProduct(prod, Uip);

  MatrixDense<double, int> S_local;
  Kpp.MatrixMatrixSubstraction(S_local, prod);

  // merge all numb_node_p in one list
  int length_list_p[numb_procs];

  MPI_Allgather(&numb_node_p, 1, MPI_INT, length_list_p, 1, MPI_INT, mpi_comm);

  std::stringstream out_length_list;

  for(int i = 0;i < numb_procs;++i){
    out_length_list << length_list_p[i] << " ";
  }
  out_length_list << "\n";

  //int list_global_p[numb_node_p];

  //for(int i = 0;i < numb_node_p;++i){
  //  list_global_p[i] = l2g[ list_node_p[i] ];
  //}

  int all_list_global_p[numb_procs][numb_global_node];

  for(int i = 0;i < numb_procs;++i){
    if(proc_numb == i){
      for(int j = 0;j < numb_node_p;++j){
        all_list_global_p[i][j] = l2g[ list_node_p[j] ];//list_global_p[j];
      }
    }

    MPI_Bcast(all_list_global_p[i], length_list_p[i], MPI_INT, i, mpi_comm);
  }

  for(int i = 0;i < numb_procs;++i){
    for(int j = 0;j < length_list_p[i];++j){
      out_length_list << all_list_global_p[i][j] << " ";
    }

    out_length_list << "\n";
  }

  //iomrg::printf("%s\n", out_length_list.str().c_str());

  // array from global positions to positions in S
  int pos_S[numb_global_node];

  for(int i = 0;i < numb_global_node;++i){
    pos_S[i] = -1;
  }

  std::vector<int> S_indices;

  for(int i = 0;i < numb_procs;++i){
    for(int j = 0;j < length_list_p[i];++j){
      S_indices.push_back(all_list_global_p[i][j]);
    }
  }

  std::sort(S_indices.begin(), S_indices.end());
  S_indices.erase(std::unique(S_indices.begin(), S_indices.end()), S_indices.end());
  int size_S = S_indices.size();

  for(int i = 0;i < size_S;++i){
    pos_S[ S_indices[i] ] = i;
  }

  // Assemble S
  MatrixDense<double, int> S;
  S.Allocate(size_S, size_S);
  S.Initialize(0);

  double aux_coef[size_S];

  iomrg::printf("size_S = %d\n",size_S);

  for(int i = 0;i < numb_procs;++i){
    int size_local_S = length_list_p[i];

    for(int j = 0;j < size_local_S;++j){
      if(proc_numb == i){
        for(int k = 0;k < size_local_S;++k){
          aux_coef[k] = S_local(j,k);
        }
      }

      MPI_Bcast(aux_coef, size_local_S, MPI_DOUBLE, i, mpi_comm);

      int r = pos_S[ all_list_global_p[i][j] ];

      for(int k = 0;k < size_local_S;++k){
        S(r, pos_S[ all_list_global_p[i][k] ]) += aux_coef[k];
      }
    }
  }

  //S.WriteToFileCsv(solution_filename.c_str());

  //MatrixDense<double, int> S_lu;
  //Factor::LU(S_lu, S);

  Vector<double, int> b_i,b_p;
  b_i.Allocate(numb_node_i);
  b_p.Allocate(numb_node_p);

  for(int i = 0;i < numb_l2g;++i){
    if(l2p[i] == -1){
      b_i( l2i[i] ) = b_local(i);
    }else{
      b_p( l2p[i] ) = b_local(i);
    }
  }
  //b_i.WriteToFileCsv( solution_filename.c_str(), '\n' );

  Vector<double, int> z;
  DirectSolver::Forward(z, Kii_lu, b_i);
  //z.WriteToFileCsv( solution_filename.c_str(), '\n' );
  //Kii_lu.WriteToFileCsv(solution_filename.c_str());

  Vector<double, int> aux_prod;
  //iomrg::printf("Lpi: %d, z: %d\n", Lpi.GetNumbColumns(), z.GetSize());
  Lpi.MatrixVectorProduct(aux_prod, z);
  //Lpi.WriteToFileCsv( solution_filename.c_str());
  //aux_prod.WriteToFileCsv( solution_filename.c_str(), '\n' );
  Vector<double, int> y_p_local;
  y_p_local.Allocate(numb_node_p);

  for(int i = 0;i < numb_node_p;++i){
    y_p_local(i) = b_p(i) - aux_prod(i);
  }

  Vector<double, int> y_p(size_S);

  for(int i = 0;i < size_S;++i){
    aux_coef[i] = 0;
  }

  for(int i = 0;i < numb_node_p;++i){
    aux_coef[ pos_S[ l2g[ list_node_p[i] ] ] ] += y_p_local(i);
  }

  MPI_Allreduce(aux_coef, y_p.GetCoef(), size_S, MPI_DOUBLE,  MPI_SUM,mpi_comm);

  Vector<double, int> x_p;
  DirectSolver::SolveLU(x_p, S, y_p);

  //iomrg::printf("%d %d %d\n",Uip.GetNumbRows(),Uip.GetNumbColumns(),x_p.GetSize());

  Vector<double, int> aux_prod_2;
  aux_prod_2.Allocate(numb_node_i);
  aux_prod_2.Assign(0, numb_node_i - 1, 0);
  iomrg::printf("Uip: %d, x_p: %d\n", Uip.GetNumbColumns(), x_p.GetSize());
  //Uip.MatrixVectorProduct(aux_prod_2, x_p);

  for(int i = 0;i < numb_node_i;++i){
    for(int j = 0;j < numb_node_p;++j){
      int id_S = pos_S[ l2g[ list_node_p[j] ] ];
      aux_prod_2(i) += Uip(i,j) * x_p(id_S);
    }
  }

  Vector<double, int> y_i;
  y_i.Allocate(numb_node_i);

  for(int i = 0;i < numb_node_i;++i){
    y_i(i) = z(i) - aux_prod_2(i);
  }
  //y_i.WriteToFileCsv( solution_filename.c_str(), '\n' );

  Vector<double, int> x_i;
  DirectSolver::Backward(x_i, Kii_lu, y_i);

  /*for(int i = 0;i < numb_procs;++i){
    for(int j = 0;j < numb_node_p;++j){
      aux_coef[j] = y_p_local(j);
    }

    MPI_Bcast(aux_coef, numb_node_p, MPI_DOUBLE, i, mpi_comm);
  }*/

  iomrg::printf("x_i.size = %d, x_p.size = %d\n",x_i.GetSize(),x_p.GetSize());

  //x_i.WriteToFileCsv( solution_filename.c_str(), '\n' );

  for(int i = 0;i < numb_node_i;++i){
    x_local(list_node_i[i]) = x_i(i);
  }

  for(int i = 0;i < numb_node_p;++i){
    x_local(list_node_p[i]) = x_p(i);
  }

}

________________________________________________________________________________


} // namespace Schur {
