#!/bin/bash

# 1. Generate a unique ID based on the exact date and time
RUN_ID=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="STT_Run_${RUN_ID}"

echo "Starting local 12-core simulation..."
echo "All data will be saved to: ${RUN_DIR}"

# --- Create the folder so Bash has a place to put the logs ---
mkdir -p $RUN_DIR

# --- RECOMPILE STEP ---
echo "Recompiling C++ solver with new Parameters.hpp..."
mpicxx -O3 main.cpp -o FP_mpi.o $(pkg-config --cflags --libs petsc)

if [ $? -ne 0 ] || [ ! -f "FP_mpi.o" ]; then
    echo "FATAL ERROR: Compilation failed! Please check the C++ errors above."
    exit 1
fi

echo "Compilation successful! Launching solver..."

# --- EXECUTION STEP ---
# Added -ts_max_steps to prevent PETSc from timing out after 10,000 steps!
# --- EXECUTION STEP ---
# Added -ts_max_steps to prevent PETSc from timing out after 10,000 steps!
time mpirun --bind-to core -n 12 ./FP_mpi.o \
    -out_dir $RUN_DIR \
    -ts_max_steps 100000000 \
    -ts_max_time 50e-9 \
    -ksp_type gmres -pc_type asm -pc_asm_overlap 0 -sub_pc_type ilu -sub_pc_factor_levels 0 \
    -ksp_rtol 1e-5 -ts_rtol 1e-5 -ts_atol 1e-5 \
    -ts_max_snes_failures -1 -ts_adapt_type basic -ts_dt_min 1e-12 \
    -log_view ascii:${RUN_DIR}/profiling_log.txt \
    -field 0.1 > ${RUN_DIR}/terminal_output.log 2>&1

echo "Simulation finished! Check the ${RUN_DIR} folder for your files."