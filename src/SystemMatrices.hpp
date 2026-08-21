#pragma once
#include "Parameters.hpp"
#include "MeshManager.hpp"
#include "PhysicsKernels.hpp"
#include <mpi.h>
#include <petscts.h>

extern Mat A, B_mat;

// This function saves the vector W to a file specified by Filename1. It opens the file in write mode, iterates over the elements of W, and writes each element to the file in scientific notation. Finally, it closes the file.
inline void save_W(string Filename1, double* W, ULLInt num_node)
{
    FILE * pFile1;
    pFile1 = fopen (Filename1.c_str(),"w");
    for (ULLInt i=0; i<num_node;i++){fprintf (pFile1, "%e\n",W[i]);}
    fclose(pFile1);
}

// This function saves the vector W to a file specified by Filename1. It opens the file in write mode, iterates over the elements of W, and writes each element to the file in scientific notation. Finally, it closes the file.
class system_matrix_B {
public:
    ULLInt my_start, my_end;

    // Constructor for the system_matrix_B class. It initializes the start and end indices for the local portion of the mesh that this process will handle, based on the total number of triangles and the rank of the process in the MPI communicator.
    system_matrix_B() {
        // Get the rank of the current process and the total number of processes in the MPI communicator
        int my_rank, num_cores;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

        ULLInt local_tri = mesh.num_tri / num_cores;
        my_start = my_rank * local_tri;
        my_end = (my_rank == num_cores - 1) ? mesh.num_tri : my_start + local_tri;
    }

    // This function fills the mass matrix B_mat by computing the integrals of the shape functions over the triangles in the local portion of the mesh. It uses Gaussian quadrature to evaluate the integrals and assembles the results into the global matrix B_mat. The function also synchronizes the assembly across all MPI processes.
    void fill_B_matrix() {
        // Loop over the local triangles assigned to this process
        for(ULLInt k = my_start; k < my_end; k++) {
            // Get the vertex indices of the current triangle
            ULLInt t_k_ull[] = {mesh.t_matrix[0][k], mesh.t_matrix[1][k], mesh.t_matrix[2][k]};
            PetscInt t_k[] = {(PetscInt)t_k_ull[0], (PetscInt)t_k_ull[1], (PetscInt)t_k_ull[2]};
            
            // Compute the edge vectors of the triangle in 3D space
            double Xu[] = {mesh.p_matrix[0][t_k[2]] - mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[2]] - mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[2]] - mesh.p_matrix[2][t_k[0]]};
            double Xv[] = {mesh.p_matrix[0][t_k[1]] - mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[1]] - mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[1]] - mesh.p_matrix[2][t_k[0]]};
            
            // Compute the metric tensor components for the triangle, which are used in the integral calculations. The metric tensor captures the geometric properties of the triangle in 3D space.
            double gA = dot(Xu, Xu); double gB = dot(Xu, Xv); double gC = dot(Xv, Xv);
            double g_ABC[] = {gA, gB, gC};
            
            PetscScalar local_B[9];
            int idx = 0;

            // Loop over the vertices of the triangle to compute the integrals of the shape functions and assemble them into the local mass matrix. The integrals are computed using the Integral_B function, which evaluates the contributions of each pair of vertices to the mass matrix.
            for (int ith = 0; ith < 3; ith++) {
                for (int jth = 0; jth < 3; jth++) {
                    local_B[idx++] = Integral_B(t_k_ull[ith], t_k_ull[jth], t_k_ull, g_ABC);
                }
            }
            PetscCallAbort(PETSC_COMM_WORLD, MatSetValues(B_mat, 3, t_k, 3, t_k, local_B, ADD_VALUES));
        }
        PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyBegin(B_mat, MAT_FINAL_ASSEMBLY));
        PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyEnd(B_mat, MAT_FINAL_ASSEMBLY));
        
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "\nMass matrix B computation and network sync complete.";
    }
};

// This class represents the system matrix A, which is used in the finite element method for solving partial differential equations. The constructor initializes the local portion of the mesh that this process will handle, and precomputes geometric quantities needed for assembling the matrix. The fill_A_matrix function computes the entries of the matrix based on the current state of the system and assembles them into the global matrix A.
class system_matrix_A {
public:

// Start and end indices for the local portion of the mesh that this process will handle
    ULLInt my_start , my_end;
    double *C1_table , *C2_table, *C3_table;
    double *P1_x, *P1_y, *P1_z;
    double *Xu_x, *Xu_y, *Xu_z;
    double *Xv_x, *Xv_y, *Xv_z;

    // Constructor for the system_matrix_A class. It initializes the start and end indices for the local portion of the mesh that this process will handle, based on the total number of triangles and the rank of the process in the MPI communicator. It also precomputes geometric quantities needed for assembling the matrix.
    system_matrix_A() {
        int my_rank, num_cores;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_cores);

        // Calculate the number of triangles assigned to this process and determine the start and end indices for the local portion of the mesh
        ULLInt local_tri = mesh.num_tri / num_cores;
        my_start = my_rank * local_tri;
        my_end = (my_rank == num_cores - 1) ? mesh.num_tri : my_start + local_tri;
        ULLInt my_num_tri = my_end - my_start;

        // Allocate memory for the geometric quantities needed for assembling the matrix. These include the coefficients C1, C2, and C3, as well as the coordinates of the triangle vertices and edge vectors.
        C1_table = new double[my_num_tri*9];
        C2_table = new double[my_num_tri*9];
        C3_table = new double[my_num_tri*9];
        P1_x = new double[my_num_tri]; P1_y = new double[my_num_tri]; P1_z = new double[my_num_tri];
        Xu_x = new double[my_num_tri]; Xu_y = new double[my_num_tri]; Xu_z = new double[my_num_tri];
        Xv_x = new double[my_num_tri]; Xv_y = new double[my_num_tri]; Xv_z = new double[my_num_tri];

        // Precompute the geometric quantities for each triangle in the local portion of the mesh. This involves calculating the coordinates of the triangle vertices, the edge vectors, and the metric tensor components. The results are stored in the corresponding arrays for later use in assembling the system matrix.
        for(ULLInt k = my_start; k < my_end; k++) {
           ULLInt local_k = k - my_start; 
           ULLInt t_k[]={mesh.t_matrix[0][k], mesh.t_matrix[1][k], mesh.t_matrix[2][k]}; 

           P1_x[local_k] = (mesh.p_matrix[0][t_k[0]]+mesh.p_matrix[0][t_k[1]]+mesh.p_matrix[0][t_k[2]])/3.0; 
           P1_y[local_k] = (mesh.p_matrix[1][t_k[0]]+mesh.p_matrix[1][t_k[1]]+mesh.p_matrix[1][t_k[2]])/3.0;
           P1_z[local_k] = (mesh.p_matrix[2][t_k[0]]+mesh.p_matrix[2][t_k[1]]+mesh.p_matrix[2][t_k[2]])/3.0;

           Xu_x[local_k] = mesh.p_matrix[0][t_k[2]] -mesh.p_matrix[0][t_k[0]];
           Xu_y[local_k] = mesh.p_matrix[1][t_k[2]] -mesh.p_matrix[1][t_k[0]];
           Xu_z[local_k] = mesh.p_matrix[2][t_k[2]] -mesh.p_matrix[2][t_k[0]];
           Xv_x[local_k] = mesh.p_matrix[0][t_k[1]] -mesh.p_matrix[0][t_k[0]];
           Xv_y[local_k] = mesh.p_matrix[1][t_k[1]] -mesh.p_matrix[1][t_k[0]];
           Xv_z[local_k] = mesh.p_matrix[2][t_k[1]] -mesh.p_matrix[2][t_k[0]];

           double tXu[] = {Xu_x[local_k], Xu_y[local_k], Xu_z[local_k]};
           double tXv[] = {Xv_x[local_k], Xv_y[local_k], Xv_z[local_k]};

           double gA = dot(tXu, tXu); double gB = dot(tXu, tXv); double gC = dot(tXv, tXv);
           double g_ABC[] = {gA, gB, gC};
           
           int idx=0;

           // Directly compute the three coefficients of Integral_A.
           // The previous MPI version obtained C1 and C2 by subtracting
           // two full Integral_A evaluations, which can cause subtractive
           // cancellation when the diffusion term is much larger than the
           // drift term. The direct coefficient calculation is algebraically
           // equivalent but avoids that cancellation.
           for(int ith=0;ith<3;ith++) {
             for(int jth=0;jth<3;jth++) {
                double C1, C2, C3;
                Integral_A_coefficients(t_k[ith], t_k[jth], t_k, g_ABC,
                                        C1, C2, C3);

                ULLInt global_idx = local_k * 9 + idx; 
                C1_table[global_idx] = C1;
                C2_table[global_idx] = C2;
                C3_table[global_idx] = C3;
                idx++;
             }
           }
        }
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "\nMatrix A Geometry precomputed natively.";
    }

    // This function fills the system matrix A by computing the contributions from each triangle in the local portion of the mesh. It uses the precomputed geometric quantities and the current state of the system to evaluate the entries of the matrix. The results are assembled into the global matrix A, and the assembly is synchronized across all MPI processes.
    void fill_A_matrix(double pulse_multiplier, double current_time_sec) {
        PetscCallAbort(PETSC_COMM_WORLD, MatZeroEntries(A));
        
        for(ULLInt k = my_start; k < my_end; k++) {
            ULLInt local_k = k - my_start; 
            ULLInt t_k_ull[] = {mesh.t_matrix[0][k], mesh.t_matrix[1][k], mesh.t_matrix[2][k]};
            PetscInt t_k[] = {(PetscInt)t_k_ull[0], (PetscInt)t_k_ull[1], (PetscInt)t_k_ull[2]};

            double P1[] = {P1_x[local_k], P1_y[local_k], P1_z[local_k]};
            double Xu[] = {Xu_x[local_k], Xu_y[local_k], Xu_z[local_k]};
            double Xv[] = {Xv_x[local_k], Xv_y[local_k], Xv_z[local_k]};
             
            auto temp = Torque(P1, Xu, Xv, pulse_multiplier, current_time_sec);

            PetscScalar local_A[9];
            int idx=0;
            for (int ith = 0; ith < 3; ith++) {
                for (int jth = 0; jth < 3; jth++) {
                    ULLInt global_idx = local_k * 9 + idx;
                    local_A[idx] = temp[0] * C1_table[global_idx] + 
                                   temp[1] * C2_table[global_idx] + 
                                   C3_table[global_idx];
                    idx++;
                }
            }
            PetscCallAbort(PETSC_COMM_WORLD, MatSetValues(A, 3, t_k, 3, t_k, local_A, ADD_VALUES));
        }
        PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    }

    ~system_matrix_A(){
        delete[] C1_table; delete[] C2_table; delete[] C3_table;
        delete[] P1_x; delete[] P1_y; delete[] P1_z;
        delete[] Xu_x; delete[] Xu_y; delete[] Xu_z;
        delete[] Xv_x; delete[] Xv_y; delete[] Xv_z;
    }
};

extern system_matrix_B* sysB;
extern system_matrix_A* sysA;

// W math functions REMAIN EXACTLY THE SAME BELOW
inline void InitializeW( double *W )
{
    double* Energy_kT = new double[mesh.num_node];
    double Emin_ref = 0.0;
    for (ULLInt i=0;i<mesh.num_node; i++)
    {
        double mx = mesh.p_matrix[0][i]; double my = mesh.p_matrix[1][i]; double mz = mesh.p_matrix[2][i];
        double mmag=sqrt(mx*mx+my*my+mz*mz); mx=mx/mmag; my=my/mmag; mz=mz/mmag;
        
        double Hext[] = {0.0,0.0,0.0};
        Hext[0]=FL.field_DC_Mag_SI*sin(FL.field_DC_theta_deg*M_PI/180.0)*cos(FL.field_DC_phi_deg*M_PI/180.0);
        Hext[1]=FL.field_DC_Mag_SI*sin(FL.field_DC_theta_deg*M_PI/180.0)*sin(FL.field_DC_phi_deg*M_PI/180.0);
        Hext[2]=FL.field_DC_Mag_SI*cos(FL.field_DC_theta_deg*M_PI/180.0);
        
        double E = (0.5*FL.mu0*FL.Ms*FL.Ms*FL.Vol/FL.kT) *(FL.Ndxx*mx*mx + FL.Ndyy*my*my + FL.Ndzz*mz*mz) + (FL.Ku*FL.Vol/FL.kT)*(1-mz*mz) - (FL.mu0*FL.Ms*FL.Vol/FL.kT)*(mx*Hext[0]+my*Hext[1]+mz*Hext[2]);
        Energy_kT[i] = E;
        if(i==1) {Emin_ref = E;}
        if (E < Emin_ref) { Emin_ref = E;}
        
        if (FL.PMA == 1) {
                if(mz>=0) { W[i] = exp(-E); } else  { W[i]= 0; }
        } else {
            if(mx>=0) { W[i] = exp(-E); } else {W[i]= 0;}
        }
    }
 
    if (Emin_ref !=0) {
        for (ULLInt i=0;i<mesh.num_node; i++) {
            double mx = mesh.p_matrix[0][i]; double my = mesh.p_matrix[1][i]; double mz = mesh.p_matrix[2][i];
            double mmag=sqrt(mx*mx+my*my+mz*mz); mx=mx/mmag; my=my/mmag; mz=mz/mmag;
            double E = Energy_kT[i]-Emin_ref;
            Energy_kT[i]=E;
            if (FL.PMA == 1) {
                    if(mz>=0) { W[i] = exp(-E); } else  { W[i]= 0; }
            } else {
                if(mx>=0) { W[i] = exp(-E); } else {W[i]= 0;}
            }
        }
    }
    save_W("Initial_E.txt", Energy_kT,mesh.num_node);
    delete[] Energy_kT;
    
    double k1[] = { 1.0/3.0 , 0.059715871789770, 0.470142064105115, 0.470142064105115  , 0.797426985353087 , 0.101286507323456 , 0.101286507323456 };
    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 ,  0.125939180544827 ,  0.125939180544827 };
    double sum_W_surface = 0.0;
    
    for (ULLInt i=0;i<mesh.num_tri;i++) {
        ULLInt t_k[]={mesh.t_matrix[0][i],mesh.t_matrix[1][i],mesh.t_matrix[2][i]};
        double Xu[] = {mesh.p_matrix[0][t_k[2]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[2]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[2]] -mesh.p_matrix[2][t_k[0]]};
        double Xv[] = {mesh.p_matrix[0][t_k[1]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[1]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[1]] -mesh.p_matrix[2][t_k[0]]};
        double Sum = 0.0 ;
        double gA = dot(Xu, Xu); double gB = dot(Xu, Xv); double gC = dot(Xv, Xv);
        for (int pp = 0; pp<7; pp++) {
             Sum = Sum + sqrt(gA*gC-gB*gB)*weight[pp] *( k1[pp]*W[t_k[0]] + k2[pp]*W[t_k[1]] +  k3[pp]*W[t_k[2]]);
        }
        sum_W_surface = sum_W_surface + Sum;
    }   
    
    for (ULLInt i=0;i<mesh.num_node; i++) {W[i]=W[i]/sum_W_surface;}
    sum_W_surface = 0.0;
    for (ULLInt i=0;i<mesh.num_tri;i++) {
        ULLInt t_k[]={mesh.t_matrix[0][i],mesh.t_matrix[1][i],mesh.t_matrix[2][i]};
        double Xu[] = {mesh.p_matrix[0][t_k[2]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[2]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[2]] -mesh.p_matrix[2][t_k[0]]};
        double Xv[] = {mesh.p_matrix[0][t_k[1]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[1]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[1]] -mesh.p_matrix[2][t_k[0]]};
        double Sum = 0.0 ;
        double gA = dot(Xu, Xu); double gB = dot(Xu, Xv); double gC = dot(Xv, Xv);
        for (int pp = 0; pp<7; pp++) {
            Sum = Sum + sqrt(gA*gC-gB*gB)*weight[pp] *( k1[pp]*W[t_k[0]] + k2[pp]*W[t_k[1]] +  k3[pp]*W[t_k[2]]);
        }
        sum_W_surface = sum_W_surface + Sum;
    }
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) cout<<"\nProbability Initialized. Norm = "<<sum_W_surface;
}

inline double checkWnormalization( const double* W)
{
    double k1[] = { 1.0/3.0 , 0.059715871789770, 0.470142064105115, 0.470142064105115  , 0.797426985353087 , 0.101286507323456 , 0.101286507323456 };
    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 ,  0.125939180544827 ,  0.125939180544827 };
    double sum_W_surface = 0.0;
    for (ULLInt i=0;i<mesh.num_tri;i++) {
        ULLInt t_k[]={mesh.t_matrix[0][i],mesh.t_matrix[1][i],mesh.t_matrix[2][i]};
        double Xu[] = {mesh.p_matrix[0][t_k[2]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[2]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[2]] -mesh.p_matrix[2][t_k[0]]};
        double Xv[] = {mesh.p_matrix[0][t_k[1]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[1]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[1]] -mesh.p_matrix[2][t_k[0]]};
        double Sum = 0.0 ;
        double gA = dot(Xu, Xu); double gB = dot(Xu, Xv); double gC = dot(Xv, Xv);
        for (int pp = 0; pp<7; pp++) {
            Sum = Sum + sqrt(gA*gC-gB*gB)*weight[pp] *( k1[pp]*W[t_k[0]] + k2[pp]*W[t_k[1]] +  k3[pp]*W[t_k[2]]);
        }
        sum_W_surface = sum_W_surface + Sum;
    }
    return sum_W_surface;
}

inline double calculate_WER(const double* W)
{
    double k1[] = { 1.0/3.0 , 0.059715871789770, 0.470142064105115, 0.470142064105115  , 0.797426985353087 , 0.101286507323456 , 0.101286507323456 };
    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 ,  0.125939180544827 ,  0.125939180544827 };
    double sum_W_surface = 0.0;
    for (ULLInt i=0;i<mesh.num_tri;i++) {
        ULLInt t_k[]={mesh.t_matrix[0][i],mesh.t_matrix[1][i],mesh.t_matrix[2][i]};
        double P1[] = {(1.0/3.0) *(mesh.p_matrix[0][t_k[0]]+mesh.p_matrix[0][t_k[1]]+mesh.p_matrix[0][t_k[2]]), (1.0/3.0) *(mesh.p_matrix[1][t_k[0]]+mesh.p_matrix[1][t_k[1]]+mesh.p_matrix[1][t_k[2]]),(1.0/3.0) *(mesh.p_matrix[2][t_k[0]]+mesh.p_matrix[2][t_k[1]]+mesh.p_matrix[2][t_k[2]])};
        
        if( ((FL.PMA==0) && (P1[0]>0.0)) || ((FL.PMA==1) && (P1[2]>FL.NS_threshold)) ){
                double Xu[] = {mesh.p_matrix[0][t_k[2]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[2]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[2]] -mesh.p_matrix[2][t_k[0]]};
                double Xv[] = {mesh.p_matrix[0][t_k[1]] -mesh.p_matrix[0][t_k[0]], mesh.p_matrix[1][t_k[1]] -mesh.p_matrix[1][t_k[0]], mesh.p_matrix[2][t_k[1]] -mesh.p_matrix[2][t_k[0]]};
                double Sum = 0.0 ;
                double gA = dot(Xu, Xu); double gB = dot(Xu, Xv); double gC = dot(Xv, Xv);
                for (int pp = 0; pp<7; pp++) {Sum = Sum + sqrt(gA*gC-gB*gB)*weight[pp] *( k1[pp]*W[t_k[0]] + k2[pp]*W[t_k[1]] +  k3[pp]*W[t_k[2]]);}
                sum_W_surface = sum_W_surface + Sum;
        }
    }
    return sum_W_surface;
}