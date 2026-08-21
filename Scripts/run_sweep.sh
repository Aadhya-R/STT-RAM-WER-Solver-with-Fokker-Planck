#!/bin/bash
# run_sweep.sh - Focused Angle Sweep (DC 400 Oe, AC 100 Oe @ 10 MHz)

export TMPDIR=$HOME/mpi_tmp
mkdir -p $TMPDIR

# ==========================================
# 1. SWEEP CONFIGURATION (LOCKED FIELDS, VARYING ANGLES)
# ==========================================
# Locking the current (Change this to whatever baseline current you are using, e.g., 2.0)
TEST_CURRENTS=(10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170)

# Locking the applied fields
TEST_DC_FIELDS=(300 400 500 600 700)
#TEST_AC_FIELDS=(100.0)

# The Sweeping Variable: Angles (Applied to BOTH AC and DC)
#TEST_DC_ANGLES=(0 15 30 45 60 75 90)
#TEST_AC_ANGLES=(0 15 30 45 60 75 90)

#TE


CORE_COUNT=4
# ==========================================

echo "============================================================"
echo " Starting Focused Angle Sweep (Synchronized AC & DC Angles)"
echo "============================================================"

SWEEP_MASTER_DIR="D56_ETAp0p3_3ns_1ns_$(date +"%Y%m%d_%H%M%S")"
mkdir -p $SWEEP_MASTER_DIR

for CURRENT in "${TEST_CURRENTS[@]}"; do
  for DC_FIELD in "${TEST_DC_FIELDS[@]}"; do
    #for AC_FIELD in "${TEST_AC_FIELDS[@]}"; do
    #  for DC_ANGLE in "${TEST_DC_ANGLES[@]}"; do
	#      for AC_ANGLE in "${TEST_AC_ANGLES[@]}";do
      
        
        #echo -e "\n➜ Testing: Angle = ${ANGLE}° (Both AC & DC) | DC = ${DC_FIELD} Oe | AC = ${AC_FIELD} Oe"
        
        # 1. Robotic text replacement (Locking frequency to 10 MHz / 1e7)
        sed -i "s/double DRIVE_CURRENT = .*/double DRIVE_CURRENT = ${CURRENT};/g" src/Parameters.hpp
        sed -i "s/double Hx_applied_Oe = .*/double Hx_applied_Oe = ${DC_FIELD};/g" src/Parameters.hpp
        #sed -i "s/double Hext_AC_Oe = .*/double Hext_AC_Oe = ${AC_FIELD};/g" Parameters.hpp
        #sed -i "s/double field_AC_freq_Hz = .*/double field_AC_freq_Hz = 1e7;/g" Parameters.hpp
        
        # Applying the angle to BOTH DC and AC fields simultaneously!
        #sed -i "s/double field_DC_theta_deg = .*/double field_DC_theta_deg = ${DC_ANGLE};/g" Parameters.hpp
       # sed -i "s/double field_AC_theta_deg = .*/double field_AC_theta_deg = ${AC_ANGLE};/g" Parameters.hpp
        
        # 2. Compile
        #mpicxx -O3 src/main.cpp -o FP_mpi.o $(pkg-config --cflags --libs petsc)
        mpicxx -O3 src/main.cpp -o FP_mpi.o -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -lpetsc

        if [ $? -ne 0 ]; then
            echo "✖ Compilation Failed. Aborting sweep."
            exit 1
        fi
        
        # 3. Setup directory
        RUN_DIR="${SWEEP_MASTER_DIR}"/"JSOT_${CURRENT}_Hx_${DC_FIELD}"
        mkdir -p $RUN_DIR
        cp src/Parameters.hpp $RUN_DIR/
        cp run_sweep.sh $RUN_DIR/
        # 4. Execute
        mpirun --bind-to core -n $CORE_COUNT ./FP_mpi.o \
            -out_dir $RUN_DIR \
            -ksp_type gmres -pc_type asm -sub_pc_type ilu \
            -ksp_rtol 1e-5 -ts_rtol 1e-5 -ts_atol 1e-5 \
            -ts_adapt_type basic \
            -log_view ascii:${RUN_DIR}/profiling_log.txt > ${RUN_DIR}/terminal_output.log 2>&1
            
        echo "✔ Run Completed"
        done
      done
   # done
 # done
#done

echo -e "\n✨ Sweep Complete! Data sets organized in: $SWEEP_MASTER_DIR"
