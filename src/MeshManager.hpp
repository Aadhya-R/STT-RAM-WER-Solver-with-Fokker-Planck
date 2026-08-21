#pragma once
#include "Parameters.hpp"
#include <vector>
#include <petscmat.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

//function to compare two rows of the row-col pair for qsort
inline int cmpfunc (const void * a, const void * b)
{
    const ULLInt *rowA = *(const ULLInt **)a;
    const ULLInt *rowB = *(const ULLInt **)b;
    if (rowA[0] > rowB[0]) return 1;
    if (rowA[0] < rowB[0]) return -1;
    if (rowA[1] > rowB[1]) return 1;
    if (rowA[1] < rowB[1]) return -1;
    return 0;
}

// Class to manage the mesh data for the spherical magnet
class spherical_mesh {
public:
    int node_share_factor = 10; // This is a guess based on the maximum number of triangles that can share a node. Adjust as necessary.
    
    // These get overwritten by setup.sh automatically
    ULLInt num_node = 35512;
    ULLInt num_tri = 71020;
    std::string mesh_p_filename = "meshes/Coarse/p_opt.txt";
    std::string mesh_t_filename = "meshes/Coarse/t_opt.txt";

    // Mesh data structures
    double** p_matrix = new double*[3];
    ULLInt** t_matrix = new ULLInt*[3];
    ULLInt** sys_row_col =  new ULLInt*[2];
    ULLInt system_nnz = 0;
    
    // Vector to store local triangles for each node
    std::vector<ULLInt> my_local_triangles;
    
    spherical_mesh(){
        // Allocate memory for p_matrix and t_matrix
        for(ULLInt i = 0; i < 3; i++)
            p_matrix[i] = new double[num_node];
            
        std::ifstream infile; // Input file stream for reading mesh files
        infile.open(mesh_p_filename.c_str(), std::ifstream::in); // Open the mesh points (nodes) file
        
        if (!infile.is_open()) {
            std::cout << "\n[FATAL ERROR] Could not find Mesh P file at: " << mesh_p_filename << "\n";
            exit(1);
        }
        
        // Read the mesh points (nodes) from the file
        ULLInt c_count = 0;
        ULLInt r_count = 0;
        double p_val;
        
        //Parsing the mesh points file robustly for both Linux and Windows
        while (infile >> p_val) {
            p_matrix[c_count][r_count] = p_val;
            c_count++;
            if (c_count == 3) { r_count++; c_count = 0; }
            if (infile.peek() == ',') infile.ignore();
        }
        infile.close();
    
        // Allocate memory for t_matrix
        for(ULLInt i = 0; i < 3; i++)
            t_matrix[i] = new ULLInt[num_tri];
        
        infile.open(mesh_t_filename.c_str(), std::ifstream::in);
        
        if (!infile.is_open()) {
            std::cout << "\n[FATAL ERROR] Could not find Mesh T file at: " << mesh_t_filename << "\n";
            exit(1);
        }
        
        c_count = 0;
        r_count = 0;
        ULLInt t_val;
        
        //Parsing the mesh triangles file robustly for both Linux and Windows
        while (infile >> t_val) {
            t_matrix[c_count][r_count] = t_val - 1;
            c_count++;
            if (c_count == 3) { r_count++; c_count = 0; }
            if (infile.peek() == ',') infile.ignore(); // Safely skip commas
        }
        infile.close();
        
        std::cout << "Mesh Loaded Successfully.";
        
    //-----------------------------------------------------------------------------
    //-------------------Constructing the row-col pair for the system matrix-------------------
    //-----------------------------------------------------------------------------

        // Allocate memory for t_conn and t_conn_unique
        ULLInt** t_conn = new ULLInt*[num_tri*20];
        ULLInt** t_conn_unique = new ULLInt*[num_node*node_share_factor];
        // Fill t_conn with all the connections from triangles to nodes
        for(ULLInt i = 0; i < num_tri*9; i++)
            t_conn[i] = new ULLInt[2];
        ULLInt i_count = 0;
        // Fill t_conn with all the connections from triangles to nodes
        for(ULLInt k = 0; k < num_tri; k++)
        {
            // Each triangle contributes 9 connections (3 nodes, each connected to 3 nodes)
            t_conn[i_count][0]= t_matrix[0][k]; t_conn[i_count][1]= t_matrix[0][k]; i_count++;
            t_conn[i_count][0]= t_matrix[1][k]; t_conn[i_count][1]= t_matrix[1][k]; i_count++;
            t_conn[i_count][0]= t_matrix[2][k]; t_conn[i_count][1]= t_matrix[2][k]; i_count++;
            t_conn[i_count][0]= t_matrix[0][k]; t_conn[i_count][1]= t_matrix[1][k]; i_count++;
            t_conn[i_count][0]= t_matrix[0][k]; t_conn[i_count][1]= t_matrix[2][k]; i_count++;
            t_conn[i_count][0]= t_matrix[1][k]; t_conn[i_count][1]= t_matrix[2][k]; i_count++;
            t_conn[i_count][1]= t_matrix[0][k]; t_conn[i_count][0]= t_matrix[1][k]; i_count++;
            t_conn[i_count][1]= t_matrix[0][k]; t_conn[i_count][0]= t_matrix[2][k]; i_count++;
            t_conn[i_count][1]= t_matrix[1][k]; t_conn[i_count][0]= t_matrix[2][k]; i_count++;
        }
        std::cout << "\nTotal number of triangles: " << num_tri;

        // Sort the t_conn array to prepare for unique filtering 
        qsort(t_conn, (num_tri*9), sizeof t_conn[0], cmpfunc);

        // Filter out unique connections to create t_conn_unique
        for(ULLInt i = 0; i < num_node*node_share_factor; i++) {t_conn_unique[i] = new ULLInt[2];}
        
        // Initialize the first unique connection
        i_count = 1;    
        t_conn_unique[0][0]=t_conn[0][0]; t_conn_unique[0][1]=t_conn[0][1];

        // Loop through the sorted connections and keep only unique pairs
        for(ULLInt k = 1; k < num_tri*9; k++)
        {   
            // Check if the current connection is different from the last unique connection
            if ((t_conn[k][0]!=t_conn[k-1][0]) || ( t_conn[k][1]!=t_conn[k-1][1]))
            {
                t_conn_unique[i_count][0]=t_conn[k][0]; t_conn_unique[i_count][1]=t_conn[k][1];i_count++;
            }
        }

        // Store the number of unique connections and allocate memory for sys_row_col
        system_nnz = i_count;
        for(ULLInt i = 0; i < 2; i++) {sys_row_col[i] = new ULLInt[system_nnz];}
        // Fill sys_row_col with the unique connections
        for(ULLInt i = 0; i < system_nnz; i++)
        {
            sys_row_col[0][i]=t_conn_unique[i][0];
            sys_row_col[1][i]=t_conn_unique[i][1];
        }
        std::cout << "\nNumber of Nodes = Number of Equations = " << num_node << "\nNumber of Non-zero elements in system: " << system_nnz << "\n";
        
        // Clean up the temporary connection arrays to free memory
        for(ULLInt i = 0; i < num_tri*9; i++) {delete t_conn[i];}
        delete[] t_conn;
        for(ULLInt i = 0; i < num_node*node_share_factor; i++) delete t_conn_unique[i];
        delete[] t_conn_unique;
    }
    ~spherical_mesh() {}
};

extern spherical_mesh mesh;