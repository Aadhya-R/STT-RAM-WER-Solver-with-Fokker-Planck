#pragma once
#include "Parameters.hpp"
#include "MeshManager.hpp"

//-----------------------------------------------------------------------------
//-------------------Physics Kernels-------------------
//-----------------------------------------------------------------------------

// This function computes the dot product of two 3D vectors a and b. It takes two arrays of size 3 as input and returns a double representing the dot product.
inline double dot (const double a [3], const double b [3])
{
    double c;
    c=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];        
    return c;
}


// This function computes the cross product of two 3D vectors a and b. It takes two arrays of size 3 as input and returns a vector representing the cross product.
inline vector<double> cross (const double a [3], const double b [3])
{   
    std::vector<double> cproduct(3,0);
    cproduct[0]=a[1]*b[2]-a[2]*b[1];
    cproduct[1]=a[2]*b[0]-a[0]*b[2];
    cproduct[2]=a[0]*b[1]-a[1]*b[0];
    return cproduct;
}

// This function computes the cross product of two 3D vectors a and b. It takes two vectors of size 3 as input and returns a vector representing the cross product.
inline vector<double> Torque (double P1[3], double Xu[3], double Xv[3], double pulse_multiplier, double current_time_sec)
{
    // Normalize the magnetization vector P1 to get the unit vector mvec
    double R = sqrt( P1[0]*P1[0] + P1[1]*P1[1] + P1[2]*P1[2]) ;
    double mvec[] = {P1[0] / R, P1[1]/R, P1[2]/R};
    
    // Calculate the effective magnetic field Hext based on the applied DC and AC fields, and the current time
    double Hext[] = {0.0,0.0,0.0};
    Hext[0]=FL.field_DC_Mag_SI*sin(FL.field_DC_theta_deg*M_PI/180.0)*cos(FL.field_DC_phi_deg*M_PI/180.0);
    Hext[1]=FL.field_DC_Mag_SI*sin(FL.field_DC_theta_deg*M_PI/180.0)*sin(FL.field_DC_phi_deg*M_PI/180.0);
    Hext[2]=FL.field_DC_Mag_SI*cos(FL.field_DC_theta_deg*M_PI/180.0);
    // Add the AC field contribution to Hext based on the current time and the AC field parameters
    Hext[0]=Hext[0] + FL.field_AC_Mag_SI*sin(FL.field_AC_theta_deg*M_PI/180.0)*cos(FL.field_AC_phi_deg*M_PI/180.0)*sin(2*M_PI*FL.field_AC_freq_Hz*current_time_sec);
    Hext[1]=Hext[1] + FL.field_AC_Mag_SI*sin(FL.field_AC_theta_deg*M_PI/180.0)*sin(FL.field_AC_phi_deg*M_PI/180.0)*sin(2*M_PI*FL.field_AC_freq_Hz*current_time_sec);
    Hext[2]=Hext[2] + FL.field_AC_Mag_SI*cos(FL.field_AC_theta_deg*M_PI/180.0)*sin(2*M_PI*FL.field_AC_freq_Hz*current_time_sec);
    
    // Calculate the effective magnetic field Heff by combining the anisotropy field and the applied field Hext
    double Heff[] = {FL.H_coeff[0]*P1[0]+Hext[0],FL.H_coeff[1]*P1[1]+Hext[1],FL.H_coeff[2]*P1[2]+Hext[2]};
    double L_inplane[]={0,0,0};
    std::vector<double> temp(3,0);

    // Calculate the torque components based on the magnetization vector mvec, the effective field Heff, and the spin-transfer torque contributions from the applied current. The torque is computed using the Landau-Lifshitz-Gilbert equation with spin-transfer torque terms.
    temp=cross(mvec, Heff); double mxH[] = {temp[0],temp[1],temp[2]};
    temp=cross(mvec,mxH ); double mxmxH[]={temp[0],temp[1],temp[2]};
    temp=cross(mvec, FL.mp);double mxmp[] ={temp[0],temp[1],temp[2]};
    temp=cross(mvec, mxmp);double mxmxmp[] = {temp[0],temp[1],temp[2]};
    
    // Calculate the torque components in the plane of the magnetization vector, taking into account the damping and spin-transfer torque contributions. The torque is scaled by the LHS_scaling_factor to account for time scaling in the simulation.
    L_inplane[0] =( 1/(1+FL.alpha*FL.alpha) )*FL.LHS_scaling_factor* (  -FL.alpha*FL.gamma* mxmxH[0] - FL.gamma*mxH[0]  - ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps       + FL.alpha*FL.eps_prime) )*mxmxmp[0] -  ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps_prime - FL.alpha*FL.eps      ) )*mxmp[0]);
    L_inplane[1] =( 1/(1+FL.alpha*FL.alpha) )*FL.LHS_scaling_factor*(  -FL.alpha*FL.gamma* mxmxH[1] - FL.gamma*mxH[1]  - ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps       + FL.alpha*FL.eps_prime) )*mxmxmp[1] -  ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps_prime - FL.alpha*FL.eps      ) )*mxmp[1]);
    L_inplane[2] =( 1/(1+FL.alpha*FL.alpha) )*FL.LHS_scaling_factor*(  -FL.alpha*FL.gamma* mxmxH[2]- FL.gamma*mxH[2]  - ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps       + FL.alpha*FL.eps_prime) )*mxmxmp[2] -  ( FL.gamma*FL.BetaSTT_wo_J*pulse_multiplier*FL.Jdensity*(FL.eps_prime - FL.alpha*FL.eps      ) )*mxmp[2]);
    
    // Calculate the metric tensor components gA, gB, and gC based on the triangle edges Xu and Xv. These components are used to compute the inverse metric tensor and project the torque onto the local coordinate system defined by the triangle.
    double gA = dot(Xu, Xu); //parallel 
    double gB = dot(Xu, Xv); //parallel
    double gC = dot(Xv, Xv); //parallel
    double detg = (gA*gC - gB*gB);

    // Calculate the components of the torque in the local coordinate system defined by the triangle edges Xu and Xv. The torque is projected onto the local basis vectors using the inverse metric tensor.
    double fu[] ={ (gC / detg )*Xu[0] - (gB/detg) * Xv[0],  (gC / detg )*Xu[1] - (gB/detg) * Xv[1],  (gC / detg )*Xu[2] - (gB/detg) * Xv[2]};
    double fv[] ={ (-gB / detg )*Xu[0] + (gA/detg) * Xv[0],  (-gB / detg )*Xu[1] + (gA/detg) * Xv[1],  (-gB / detg )*Xu[2] + (gA/detg) * Xv[2]};
    std::vector<double> L_2comp(2,0);   
    L_2comp[0] = dot(L_inplane, fu);
    L_2comp[1] = dot(L_inplane, fv);
    return L_2comp;
}

// This function calculates the shape functions for a given vertex index i and triangle defined by t_k, at the local coordinates (u,v) within the triangle. The shape functions are used in the finite element method to interpolate values across the triangle based on the values at the vertices.
inline vector<double> fun_phi(ULLInt i,ULLInt t_k[3], double u, double v)
{
    // Determine which vertex of the triangle corresponds to the index i and assign the appropriate shape function values based on the local coordinates (u,v). The shape functions are linear and sum to 1 across the triangle.
    int Point_No = 0;
    for (int n = 0; n<3; n++)
    {
         if (i == t_k[n]) { Point_No = n; };
    }

    // Initialize a vector to hold the shape function values for the three vertices of the triangle
    std::vector<double> Phi_Arr(3,0);
    // Assign shape function values based on the vertex index (Point_No) and the local coordinates (u,v). The shape functions are defined such that they are 1 at the corresponding vertex and 0 at the other vertices.
    if (Point_No == 0) { Phi_Arr [ 0 ] = -1.; Phi_Arr [ 1 ] = -1.; Phi_Arr [ 2 ] = 1.0 - u - v; }
    if (Point_No == 1) { Phi_Arr[0] = 0.; Phi_Arr[1] = 1.; Phi_Arr[2] = v ; }
    if (Point_No == 2) { Phi_Arr[0] = 1.; Phi_Arr[1] = 0. ; Phi_Arr[2] = u ; }
    return Phi_Arr;
}

// This function computes the integral of the shape functions and their derivatives over a triangle defined by vertices i and j, using Gaussian quadrature. The integral is weighted by the metric tensor components g_ABC and the area of the triangle defined by L. The result is used in assembling the finite element system matrix.
inline double Integral_A(ULLInt i,ULLInt j, ULLInt t_k[3], double L[2], double g_ABC[3])
{
    // Extract the metric tensor components from the input array g_ABC
    double gA = g_ABC[0]; double gB = g_ABC[1]; double gC = g_ABC[2];
    double G_ij = sqrt(gA*gC - gB*gB);
    double g_IJ_11 = gC/(G_ij*G_ij);
    double g_IJ_12 = -gB/(G_ij*G_ij);
    double g_IJ_21 = -gB/(G_ij*G_ij);
    double g_IJ_22 = gA/(G_ij*G_ij);
    
    // Initialize the sum for the integral result
    double Sum_A = 0.;
    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 ,  0.125939180544827 ,  0.125939180544827 };
    
    // Loop over the Gaussian quadrature points to compute the weighted sum of the integrand evaluated at each point. The integrand involves the shape functions and their derivatives, as well as the metric tensor components.
    for (int pp=0; pp<7; pp++)
    {
        std::vector<double> temp(3,0); temp = fun_phi(i,t_k,k3[pp],k2[pp]);
        double Phi_i[]={temp[0], temp[1], temp[2]};
        temp = fun_phi(j,t_k,k3[pp],k2[pp]);
        double Phi_j[]={temp[0], temp[1], temp[2]};
        
        Sum_A = Sum_A + weight[pp] * G_ij * ( Phi_i[0]*L[0]*Phi_j[2] + Phi_i[1]*L[1]*Phi_j[2]  - FL.LHS_scaling_factor*FL.diffusion* (   Phi_i[0] * (g_IJ_11*Phi_j[0] + g_IJ_12*Phi_j[1]) + Phi_i[1] * (g_IJ_21*Phi_j[0] + g_IJ_22*Phi_j[1]) ) ) ;
    }
    return Sum_A;
}

// This function computes the integral of the shape functions over a triangle defined by vertices i and j, using Gaussian quadrature. The integral is weighted by the metric tensor components g_ABC. The result is used in assembling the finite element system matrix.
// Directly computes the three linear coefficients of Integral_A without
// obtaining C1 and C2 by subtracting two nearly equal integrals.
//
// Integral_A(i,j,L,g) = L[0]*C1 + L[1]*C2 + C3.
// C1 and C2 are the drift terms; C3 is the diffusion term.
// This is algebraically equivalent to Integral_A but avoids the
// numerical cancellation introduced by
//     Integral_A(L={1,0}) - Integral_A(L={0,0})
// and
//     Integral_A(L={0,1}) - Integral_A(L={0,0}).
inline void Integral_A_coefficients(ULLInt i, ULLInt j, ULLInt t_k[3],
                                     double g_ABC[3],
                                     double &C1, double &C2, double &C3)
{
    double gA = g_ABC[0];
    double gB = g_ABC[1];
    double gC = g_ABC[2];
    double G_ij = sqrt(gA*gC - gB*gB);
    double g_IJ_11 = gC/(G_ij*G_ij);
    double g_IJ_12 = -gB/(G_ij*G_ij);
    double g_IJ_21 = -gB/(G_ij*G_ij);
    double g_IJ_22 = gA/(G_ij*G_ij);

    C1 = 0.0;
    C2 = 0.0;
    C3 = 0.0;

    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 , 0.125939180544827 , 0.125939180544827 };

    for (int pp=0; pp<7; pp++)
    {
        std::vector<double> temp(3,0);
        temp = fun_phi(i,t_k,k3[pp],k2[pp]);
        double Phi_i[]={temp[0], temp[1], temp[2]};

        temp = fun_phi(j,t_k,k3[pp],k2[pp]);
        double Phi_j[]={temp[0], temp[1], temp[2]};

        double prefactor = weight[pp] * G_ij;

        C1 += prefactor * Phi_i[0] * Phi_j[2];
        C2 += prefactor * Phi_i[1] * Phi_j[2];
        C3 += prefactor * (
            -FL.LHS_scaling_factor*FL.diffusion *
            (Phi_i[0] * (g_IJ_11*Phi_j[0] + g_IJ_12*Phi_j[1]) +
             Phi_i[1] * (g_IJ_21*Phi_j[0] + g_IJ_22*Phi_j[1]))
        );
    }
}

inline double Integral_B(ULLInt i,ULLInt j, ULLInt t_k[3], double g_ABC[3])
{
    // Extract the metric tensor components from the input array g_ABC
    double gA = g_ABC[0]; double gB = g_ABC[1]; double gC = g_ABC[2];
    double G_ij = sqrt(gA*gC - gB*gB);
    double Sum_B = 0.;
    
    // Define the Gaussian quadrature points and weights for numerical integration over the triangle. The points are defined in barycentric coordinates, and the weights are used to compute the weighted sum of the integrand evaluated at each point.
    double k2[] = {1.0/3.0 , 0.470142064105115 , 0.059715871789770 , 0.470142064105115 , 0.101286507323456 , 0.797426985353087 , 0.101286507323456 };
    double k3[] = {1.0/3.0 , 0.470142064105115 , 0.470142064105115 , 0.059715871789770 , 0.101286507323456 , 0.101286507323456 , 0.797426985353087 };
    double weight[] = { 0.225 , 0.132394152788506 , 0.132394152788506 , 0.132394152788506,  0.125939180544827 ,  0.125939180544827 ,  0.125939180544827 };
    
    // Loop over the Gaussian quadrature points to compute the weighted sum of the integrand evaluated at each point. The integrand involves the shape functions and their derivatives, as well as the metric tensor components.
    for (int pp=0; pp<7; pp++)
    {
        std::vector<double> temp(3,0); temp = fun_phi(i,t_k,k3[pp],k2[pp]);
        double Phi_i[]={temp[0], temp[1], temp[2]};
        temp = fun_phi(j,t_k,k3[pp],k2[pp]);
        double Phi_j[]={temp[0], temp[1], temp[2]};
        
        Sum_B = Sum_B + weight[pp] * G_ij * Phi_i[2] * Phi_j[2];
    }
    return Sum_B;
}