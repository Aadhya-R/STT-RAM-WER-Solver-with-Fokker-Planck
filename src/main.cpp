#include <petscts.h>
#include "Parameters.hpp"
#include "MeshManager.hpp"
#include "PhysicsKernels.hpp"
#include "SystemMatrices.hpp"
#include <filesystem>
#include <string>
#include <fstream>
#include <sys/stat.h>
using namespace std;

// Native Global Matrices and Workspaces!
Mat A, B_mat;
Vec temp_Ax, temp_Bxdot;
VecScatter scatter_ctx;
Vec x_seq;

system_matrix_B* sysB = nullptr;
system_matrix_A* sysA = nullptr;

const magnet_parameter FL;
spherical_mesh mesh;

double pulse_previous = 0.0;
double last_log = 0.0;
bool create_log = true;

// NEW: Global string to hold the dynamic folder path
std::string actual_output_file = "Pns_vs_Time.txt";

inline double calculate_multiplier(double t_now) {
    double pulse_multiplier=0.0;
    if (FL.step_pulse) {
        if (t_now <= FL.solve_wo_current) {pulse_multiplier=0.0; }
        else if ( (t_now > FL.solve_wo_current) && (t_now<=(FL.solve_wo_current+FL.pulse_width)) ) {pulse_multiplier=1.0; }
        else if (t_now > (FL.solve_wo_current+FL.pulse_width)) {pulse_multiplier=0.0; }
        else {pulse_multiplier=0.0;}
    }
    if (FL.ramp_up_down) {
        if (t_now <= FL.ramp_up_time) {pulse_multiplier=t_now/FL.ramp_up_time; }
        else if ((t_now>=FL.ramp_up_time) && (t_now<=(FL.ON_time+FL.ramp_up_time))) {pulse_multiplier=1.0;}
        else if ((t_now > (FL.ON_time+FL.ramp_up_time)) && (t_now<=(FL.ON_time+FL.ramp_up_time+FL.ramp_down_time))) {pulse_multiplier=1.0 - ((t_now-(FL.ON_time+FL.ramp_up_time))/FL.ramp_down_time); }
        else if (t_now > (FL.ON_time+FL.ramp_up_time+FL.ramp_down_time)) {pulse_multiplier=0.0; }
        else {pulse_multiplier=0.0;}
    }
    return pulse_multiplier;
}

inline void decide_and_update_sysA(double pulse_multiplier, double current_time_sec) {
    if ((pulse_previous!=pulse_multiplier) || (FL.field_AC_Mag_SI!=0) ) {
        pulse_previous = pulse_multiplier;
        sysA->fill_A_matrix(pulse_multiplier, current_time_sec);
    } 
}

static PetscErrorCode FokkerPlanckDAE(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void *ctx) {
    PetscFunctionBeginUser;
    
    double mult_now = calculate_multiplier(t);
    decide_and_update_sysA(mult_now, t);
    
    //matrix multiplication
    PetscCall(MatMult(A, X, temp_Ax)); // A*X
    PetscCall(MatMult(B_mat, Xdot, temp_Bxdot)); // B*Xdot
    
    // F = (B * Xdot) - (A * X)
    PetscCall(VecWAXPY(F, -1.0, temp_Ax, temp_Bxdot)); //shared memory, so no need to free temp_Ax and temp_Bxdot
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FokkerPlanckJacobian(TS ts, PetscReal t, Vec X, Vec Xdot, PetscReal a, Mat Jac, Mat JacPre, void *ctx) {
    PetscFunctionBeginUser;
    
    double mult_now = calculate_multiplier(t);
    decide_and_update_sysA(mult_now, t);
    
    //JACOBIAN: Jac = a*B - A
    PetscCall(MatCopy(B_mat, Jac, SAME_NONZERO_PATTERN)); // Jac = B
    PetscCall(MatScale(Jac, a)); // Jac = a*B
    PetscCall(MatAXPY(Jac, -1.0, A, SAME_NONZERO_PATTERN)); // Jac = a*B - A
    
    if (Jac != JacPre) 
    { 
        PetscCall(MatAssemblyBegin(JacPre, MAT_FINAL_ASSEMBLY)); 
        PetscCall(MatAssemblyEnd(JacPre, MAT_FINAL_ASSEMBLY));
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FPSetSolution(PetscReal t, Vec X, void *ctx) {
    double* W_vec = new double[mesh.num_node];
    InitializeW( W_vec );
    PetscFunctionBeginUser;
    PetscCheck(t == 0, PETSC_COMM_WORLD, PETSC_ERR_SUP, "not implemented");
    
    PetscInt rstart, rend;
    PetscCall(VecGetOwnershipRange(X, &rstart, &rend));
    for(PetscInt i = rstart; i < rend; i++) {
        PetscCall(VecSetValue(X, i, W_vec[i], INSERT_VALUES));
    }
    
    PetscCall(VecAssemblyBegin(X));
    PetscCall(VecAssemblyEnd(X));
    
    delete[] W_vec;
    PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MonitorError(TS ts, PetscInt step, PetscReal t, Vec x, void *ctx) {
    double t_now=t; double time_elapsed = (t_now -last_log)*1e12;
    if ( (time_elapsed > 1.0*FL.log_and_save_every_x_ps ) || (create_log==true) ) {
        PetscReal h;
        PetscFunctionBeginUser;
        PetscCall(TSGetTimeStep(ts, &h));

        // We only scatter locally once every picosecond!
        PetscCall(VecScatterBegin(scatter_ctx, x, x_seq, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(scatter_ctx, x, x_seq, INSERT_VALUES, SCATTER_FORWARD));

        const PetscScalar *x_arr;
        PetscCall(VecGetArrayRead(x_seq, &x_arr));

        double min_W=0.0;
        double *x_val = new double[mesh.num_node];
        double *x_val_neg = new double[mesh.num_node];
        
        for (ULLInt iterm=0; iterm<mesh.num_node; iterm++) { 
            if (min_W > x_arr[iterm]) { min_W = x_arr[iterm]; }
            if (x_arr[iterm] < 0.0) { x_val_neg[iterm] = x_arr[iterm]; } else { x_val_neg[iterm] = 0.0; }
            x_val[iterm] = x_arr[iterm];
        }
        PetscCall(VecRestoreArrayRead(x_seq, &x_arr));
        
        double checknorm = checkWnormalization(x_val); 
        double checksum = checkWnormalization(x_val_neg);
        double Pns = calculate_WER(x_val); 
        double PnsN = calculate_WER(x_val_neg);
        
        int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            FILE * Pns_file;
            if (create_log)  { Pns_file = fopen(actual_output_file.c_str(),"w"); } 
            else {  Pns_file = fopen(actual_output_file.c_str(),"a");}
            fprintf(Pns_file, "%E\t%E\t%0.3f\t%E\t%E\t%0.3f\n",(double)t,Pns,checknorm,min_W,(double)h,calculate_multiplier(t_now));   
            double m=calculate_multiplier(t_now);
            printf("\nstep t=%12.8e h=% 8.2e min(W)=%e checkNsum=%e Sum(W)=%e  Pns=%e PnsN=%e Pulse=%e", (double)t, (double)h, min_W,checksum, checknorm, Pns,PnsN,m);
            fclose(Pns_file);
        }
        
        delete[] x_val;
        delete[] x_val_neg;
        
        create_log = false;
        last_log = t_now;
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EventFunction(TS ts, PetscReal t, Vec U, PetscScalar *fvalue, void *ctx) {
  PetscFunctionBeginUser;
  if (FL.step_pulse) {  fvalue[0] = (t-FL.pulse_width-FL.max_step)*1e12; }
  if (FL.ramp_up_down) { fvalue[0] = (t-FL.ON_time-FL.ramp_up_time-FL.max_step)*1e12;}
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PostEventFunction(TS ts, PetscInt nevents, PetscInt event_list[], PetscReal t, Vec U, PetscBool forwardsolve, void *ctx) {
  PetscFunctionBeginUser;
  PetscCall(TSSetTimeStep(ts, FL.min_step));
  PetscCall(TSSetFromOptions(ts));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) cout<<"\nPulse Event Detected.\nTime step set to minimum value at pulse event";
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
    TS ts; TSAdapt adapt; Vec x, r; Mat Jac;
    PetscReal ftime; void *dummy_data=NULL;
    
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    
    int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) cout<<"\nCurrent Density (x1e10 MA/cm2):"<<FL.Jdensity/1e10;

    // --- NEW: FOLDER AUTOMATION ---
    char out_dir[256] = "Simulation_Results"; 
    PetscOptionsGetString(NULL, NULL, "-out_dir", out_dir, sizeof(out_dir), NULL);
    std::string output_folder(out_dir);

    // Rank 0 creates the folder to prevent MPI crashes
    if (rank == 0) {
       mkdir(output_folder.c_str(), 0777);
    }
    // Force all other 35 cores to wait until Rank 0 finishes building the folder
    MPI_Barrier(PETSC_COMM_WORLD); 

    // Update the global file path to point inside the new folder
    actual_output_file = output_folder + "/Pns_vs_Time.txt";
    // ------------------------------

    // PREALLOCATION FOR MATRIX A
    PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
    PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mesh.num_node, mesh.num_node));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatMPIAIJSetPreallocation(A, 100, NULL, 100, NULL)); 
    PetscCall(MatSeqAIJSetPreallocation(A, 100, NULL));
    PetscCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCall(MatSetUp(A));

    // PREALLOCATION FOR MATRIX B
    PetscCall(MatCreate(PETSC_COMM_WORLD, &B_mat));
    PetscCall(MatSetSizes(B_mat, PETSC_DECIDE, PETSC_DECIDE, mesh.num_node, mesh.num_node));
    PetscCall(MatSetFromOptions(B_mat));
    PetscCall(MatMPIAIJSetPreallocation(B_mat, 100, NULL, 100, NULL)); 
    PetscCall(MatSeqAIJSetPreallocation(B_mat, 100, NULL));
    PetscCall(MatSetOption(B_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscCall(MatSetUp(B_mat));

    sysB = new system_matrix_B();
    sysB->fill_B_matrix();
    
    sysA = new system_matrix_A();
    sysA->fill_A_matrix(0.0, 0.0);

    // ALIGN GLOBAL VECTORS WITH MATRICES
    PetscCall(MatCreateVecs(A, &x, NULL));
    PetscCall(VecDuplicate(x, &r));
    PetscCall(VecDuplicate(x, &temp_Ax));
    PetscCall(VecDuplicate(x, &temp_Bxdot));
    
    // SINGLE SCATTER CONTEXT FOR THE LOGGING (Rank 0 gets full vector)
    PetscCall(VecScatterCreateToAll(x, &scatter_ctx, &x_seq));
    
    // NATIVE JACOBIAN ALLOCATION
    PetscCall(MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &Jac));

    PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
    PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
    PetscCall(TSSetType(ts, TSCN)); 
    PetscCall(TSSetIFunction(ts, r, &FokkerPlanckDAE, dummy_data));
    PetscCall(TSSetIJacobian(ts, Jac, Jac, &FokkerPlanckJacobian, dummy_data));
    PetscCall(TSSetMaxTime(ts, FL.total_time));
    PetscCall(TSGetAdapt(ts, &adapt));
    PetscCall(TSAdaptSetStepLimits(adapt, FL.min_step, FL.max_step));
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER));
    PetscCall(TSSetMaxStepRejections(ts, 8));
    PetscCall(TSSetMaxSNESFailures(ts, 3)); 
    PetscCall(TSMonitorSet(ts, &MonitorError, dummy_data, NULL));
    PetscCall(TSSetTolerances(ts, FL.atol, NULL, FL.rtol, NULL));
  
    PetscCall((FPSetSolution)(0, x, dummy_data));
    PetscCall(TSSetTimeStep(ts, FL.min_step));
    PetscCall(TSSetSolution(ts, x));
    PetscCall(TSSetFromOptions(ts));

    if (rank == 0) cout<<"\nStart solving now.\t"<<mesh.num_node;
    PetscCall(TSSolve(ts, x));
    
    PetscCall(TSGetSolveTime(ts, &ftime));

    delete sysA;
    delete sysB;
    PetscCall(VecScatterDestroy(&scatter_ctx));
    PetscCall(VecDestroy(&x_seq));
    PetscCall(VecDestroy(&temp_Ax));
    PetscCall(VecDestroy(&temp_Bxdot));
    PetscCall(MatDestroy(&A));
    PetscCall(MatDestroy(&B_mat));
    PetscCall(MatDestroy(&Jac));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&r));
    PetscCall(TSDestroy(&ts));
    PetscCall(PetscFinalize());
    return 0;
}
