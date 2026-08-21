#pragma once
#include <iostream>
#include <new>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;
typedef int ULLInt; 

//Global parameters for the magnet
class magnet_parameter {
public:
//-----------------------------------------------------------------------------
//-------------------Physical Constants and Simulation Parameters-------------------
//-----------------------------------------------------------------------------
    double C=1.602e-19; //elementary charge
    double mu0=4e-7*M_PI ; //permeability of free space
    double hbar=1.054571628e-34; //reduced Planck constant
    double gamma=2.21e5 ; //gyromagnetic ratio
    double mu_B=9.274e-24; //Bohr magneton
    double kT=1.3806504e-23*300.0; //Boltzmann constant times temperature

//-----------------------------------------------------------------------------
//-------------------Magnet Parameters-------------------
//-----------------------------------------------------------------------------
    double Lx = 56e-9; // length in x-direction
    double Ly = 56e-9; // length in y-direction
    double Lz = 1.0e-9; // length in z-direction
    double Vol = (M_PI/4.0)*Lx*Ly*Lz; // volume of the magnet
    double area  = (M_PI/4.0)*Lx*Ly; // area of the magnet
    double Ms = 1.1e6; // saturation magnetization
    double Ku=8.0e5; //uniaxial anisotropy constant
    double alpha = 0.2; //damping constant
    double pol = -0.32; //spin polarization
    double diffusion=alpha*kT*gamma/((1.0+alpha*alpha)*mu0*Ms*Vol); //diffusion constant for thermal noise
    double Ndxx = 0.0272; //changes if dimension changes - mention to use demag 
    double Ndyy = 0.0272; //changes if dimension changes - mention to use demag
    double Ndzz = 0.9456; //changes if dimension changes - mention to use demag
    double H_coeff[3]={-Ndxx*Ms,-Ms*Ndyy,-Ms*Ndzz+2.0*Ku/(mu0*Ms)}; //effective field coefficients for the magnet
    
    double BetaSTT_wo_J = hbar/(mu0*C*Lz*Ms);
    double eps = pol/2.0;
    double eta = -0.30;
    double eps_prime =eta*eps;

//-----------------------------------------------------------------------------
//-------------------Simulation Parameters-------------------
//-----------------------------------------------------------------------------

   double DRIVE_CURRENT = 170;
    double Jdensity = (DRIVE_CURRENT)*1e10;//120e10;//-(DRIVE_CURRENT)*38e-6/(M_PI*56*56*1e-18*0.25); // calculated applied current 

// DC Field Parameters
    double Hx_applied_Oe = 700;
    double field_DC_theta_deg = 90.0; //theta angle of the applied DC field in degrees
    double field_DC_phi_deg = 0.0; // phi angle of the applied DC field in degrees
    double field_DC_Mag_SI = Hx_applied_Oe*1e-4/mu0; //Applied DC field in SI units (A/m)

// AC Field Parameters    
    double Hext_AC_Oe = 0.00; //Applied AC field in Oersted
    double field_AC_Mag_SI = Hext_AC_Oe*1e-4/mu0; //Applied AC field in SI units (A/m)
    double field_AC_freq_Hz = 1e7; //AC field frequency in Hz
    double field_AC_theta_deg = 0; //theta angle of the applied AC field in degrees
    double field_AC_phi_deg = 0.0; //phi angle of the applied AC field in degrees

//Other Parameters    
    double NS_threshold=-0.0; //threshold for determining switching event based on the z-component of magnetization
    double theta_z = 90*3.14159/180.0; //angle of the applied field in radians
    double mp[3]={0,1,0}; //unit vector of the applied field
    int PMA = 1; //1 for PMA, 0 for IMA 

    double Hscale = 1.0; //scaling factor for the applied field, can be changed to scale the applied field
    double t_scale_factor_in_s = 1.0; //scaling factor for the time, can be changed to scale the time
    double LHS_scaling_factor=t_scale_factor_in_s; //scaling factor for the LHS of the LLG equation, can be changed to scale the LHS of the LLG equation
    
//-----------------------------------------------------------------------------
//-------------------Time and Solver Parameters-------------------
//-----------------------------------------------------------------------------

//Step Pulse Parameters
    bool step_pulse = false; //step pulse or ramp pulse, true for step pulse, false for ramp pulse 
    double solve_wo_current = 0e-9; //time to solve the LLG equation without current before applying the current pulse
    double pulse_width = 5e-9; //pulse width of the current pulse, can be changed
    
//Ramp Pulse Parameters
    bool ramp_up_down = true; // changable 
    double ramp_up_time = 1e-9; // changable 
    double ramp_down_time = 1e-9; //changable 
    double ON_time = 3e-9; //on time can be changed 
    
//Solver Parameters
    double min_step = 1e-13; //minimum step size for the solver 
    double max_step = 5e-12; //maximum step size for the solver
    double rtol = 1e-8; //relative tolerance for the solver
    double atol = 1e-10; //absolute tolerance for the solver
    double total_time = 15e-9; //total time for the simulation, can be changed
    
    int log_and_save_every_x_ps = 20; //log and save the results every x picoseconds
};

/*
//Global parameters for the magnet
class magnet_parameter {
public:
//-----------------------------------------------------------------------------
//-------------------Physical Constants and Simulation Parameters-------------------
//-----------------------------------------------------------------------------
    double C=1.602e-19; //elementary charge
    double mu0=4e-7*3.14159 ; //permeability of free space
    double hbar=1.054571628e-34; //reduced Planck constant
    double gamma=2.21e5 ; //gyromagnetic ratio
    double mu_B=9.274e-24; //Bohr magneton
    double kT=1.3806504e-23*300.0; //Boltzmann constant times temperature

//-----------------------------------------------------------------------------
//-------------------Magnet Parameters-------------------
//-----------------------------------------------------------------------------
    double Lx = 56e-9; // length in x-direction
    double Ly = 56e-9; // length in y-direction
    double Lz = 1.0e-9; // length in z-direction
    double Vol = (M_PI/4.0)*Lx*Ly*Lz; // volume of the magnet
    double area  = (M_PI/4.0)*Lx*Ly; // area of the magnet
    double Ms = 1.1e6; // saturation magnetization
    double Ku=8.0e5; //uniaxial anisotropy constant
    double alpha = 0.01; //damping constant
    double pol = 0.4; //spin polarization
    double diffusion=alpha*kT*gamma/((1.0+alpha*alpha)*mu0*Ms*Vol); //diffusion constant for thermal noise
    double Ndxx = 0.0279; //changes if dimension changes - mention to use demag 
    double Ndyy = 0.0279; //changes if dimension changes - mention to use demag
    double Ndzz = 0.9442; //changes if dimension changes - mention to use demag
    double H_coeff[3]={-Ndxx*Ms,-Ms*Ndyy,-Ms*Ndzz+2.0*Ku/(mu0*Ms)}; //effective field coefficients for the magnet
    
    double BetaSTT_wo_J = hbar/(mu0*C*Lz*Ms);
    double eps = pol/2;
    double eta = 0;
    double eps_prime =eta*eps;

//-----------------------------------------------------------------------------
//-------------------Simulation Parameters-------------------
//-----------------------------------------------------------------------------

    double DRIVE_CURRENT = 170;
    double Jdensity = -(DRIVE_CURRENT)*38e-6/(M_PI*56*56*1e-18*0.25); // calculated applied current 

// DC Field Parameters
    double Hx_applied_Oe = 700;
    double field_DC_theta_deg = 0; //theta angle of the applied DC field in degrees
    double field_DC_phi_deg = 0.0; // phi angle of the applied DC field in degrees
    double field_DC_Mag_SI = Hx_applied_Oe*1e-4/mu0; //Applied DC field in SI units (A/m)

// AC Field Parameters    
    double Hext_AC_Oe = 0.00; //Applied AC field in Oersted
    double field_AC_Mag_SI = Hext_AC_Oe*1e-4/mu0; //Applied AC field in SI units (A/m)
    double field_AC_freq_Hz = 1e7; //AC field frequency in Hz
    double field_AC_theta_deg = 0; //theta angle of the applied AC field in degrees
    double field_AC_phi_deg = 0.0; //phi angle of the applied AC field in degrees

//Other Parameters    
    double NS_threshold=-0.0; //threshold for determining switching event based on the z-component of magnetization
    double theta_z = 90*3.14159/180.0; //angle of the applied field in radians
    double mp[3]={0,0,1}; //unit vector of the applied field
    int PMA = 1; //1 for PMA, 0 for IMA 

    double Hscale = 1.0; //scaling factor for the applied field, can be changed to scale the applied field
    double t_scale_factor_in_s = 1.0; //scaling factor for the time, can be changed to scale the time
    double LHS_scaling_factor=t_scale_factor_in_s; //scaling factor for the LHS of the LLG equation, can be changed to scale the LHS of the LLG equation
    
//-----------------------------------------------------------------------------
//-------------------Time and Solver Parameters-------------------
//-----------------------------------------------------------------------------

//Step Pulse Parameters
    bool step_pulse = true; //step pulse or ramp pulse, true for step pulse, false for ramp pulse 
    double solve_wo_current = 0e-9; //time to solve the LLG equation without current before applying the current pulse
    double pulse_width = 50e-9; //pulse width of the current pulse, can be changed
    
//Ramp Pulse Parameters
    bool ramp_up_down = false; // changable 
    double ramp_up_time = 0.2e-9; // changable 
    double ramp_down_time = 0.2e-9; //changable 
    double ON_time = 1e-9; //on time can be changed 
    
//Solver Parameters
    double min_step = 1e-13; //minimum step size for the solver 
    double max_step = 5e-12; //maximum step size for the solver
    double rtol = 1e-8; //relative tolerance for the solver
    double atol = 1e-10; //absolute tolerance for the solver
    double total_time = 20e-9; //total time for the simulation, can be changed
    
    int log_and_save_every_x_ps = 20; //log and save the results every x picoseconds
};
*/



extern const magnet_parameter FL;
