#!/bin/bash
# setup.sh - Interactive Installer & Runner for STT-RAM FP Solver

# --- ANSI Color Codes ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# --- Header ---
clear
echo -e "${CYAN}${BOLD}"
echo "  ____ _____ _____        ____   ___  _     __     _______ ____  "
echo " / ___|_   _|_   _|      / ___| / _ \| |    \ \   / / ____|  _ \ "
echo " \___ \ | |   | |  _____ \___ \| | | | |     \ \ / /|  _| | |_) |"
echo "  ___) || |   | | |_____| ___) | |_| | |___   \ V / | |___|  _ < "
echo " |____/ |_|   |_|        |____/ \___/|_____|   \_/  |_____|_| \_\\"
echo "================================================================="
echo "               Fokker-Planck Solver Setup Wizard                 "
echo "================================================================="
echo -e "${NC}"

# 1. Ask for cores
echo -e "${YELLOW}➜ Enter number of physical CPU cores to use ${NC}[Default: ${BOLD}12${NC}]: \c"
read CORE_COUNT
CORE_COUNT=${CORE_COUNT:-12}

# 2. Ask for mesh
echo -e "\n${YELLOW}➜ Select Mesh Density:${NC}"
echo -e "  ${MAGENTA}1)${NC} Coarse Mesh (Fast)"
echo -e "  ${MAGENTA}2)${NC} Fine Mesh (High Accuracy)"
echo -e "${YELLOW}Selection (1/2) ${NC}[Default: ${BOLD}1${NC}]: \c"
read MESH_SEL
MESH_SEL=${MESH_SEL:-1}

# Configure variables based on selection
if [ "$MESH_SEL" == "2" ]; then
    MESH_TYPE="Fine"
    MESH_P_PATH="meshes/Fine/p_opt.txt"
    MESH_T_PATH="meshes/Fine/t_opt.txt"
    NUM_NODE=878179
    NUM_TRI=1756354
else
    MESH_TYPE="Coarse"
    MESH_P_PATH="meshes/Coarse/p_opt.txt"
    MESH_T_PATH="meshes/Coarse/t_opt.txt"
    NUM_NODE=35512
    NUM_TRI=71020
fi

echo -e "\n${BLUE}[INFO]${NC} Selected ${BOLD}$MESH_TYPE${NC} mesh (${GREEN}$NUM_NODE${NC} Nodes, ${GREEN}$NUM_TRI${NC} Triangles)."

# 3. Read parameters dynamically from Parameters.hpp
# 3. Read parameters dynamically from Parameters.hpp
CURRENT=$(grep "double DRIVE_CURRENT" src/Parameters.hpp | grep -Eo '[+-]?[0-9]*\.?[0-9]+' | head -1)
DC_MAG=$(grep "double Hx_applied_Oe" src/Parameters.hpp | grep -Eo '[+-]?[0-9]*\.?[0-9]+' | head -1)
AC_MAG=$(grep "double Hext_AC_Oe" src/Parameters.hpp | grep -Eo '[+-]?[0-9]*\.?[0-9]+' | head -1)
DC_ANG=$(grep "double field_DC_theta_deg" src/Parameters.hpp | grep -Eo '[+-]?[0-9]*\.?[0-9]+' | head -1)
AC_ANG=$(grep "double field_AC_theta_deg" src/Parameters.hpp | grep -Eo '[+-]?[0-9]*\.?[0-9]+' | head -1)

# 4. Determine highly specific field string
FIELD_STR=$(awk -v dc="$DC_MAG" -v ac="$AC_MAG" -v dc_a="$DC_ANG" -v ac_a="$AC_ANG" '
BEGIN {
    if (dc == 0 && ac == 0) {
        printf "NoField"
    } else if (dc != 0 && ac == 0) {
        printf "DC_%sOe_Ang%s", dc, dc_a
    } else if (dc == 0 && ac != 0) {
        printf "AC_%sOe_Ang%s", ac, ac_a
    } else {
        printf "Both_DC%s_AC%s_Ang%s", dc, ac, dc_a
    }
}')

# 5. Generate Dynamic Output Names
RUN_ID=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="STT_Run_I${CURRENT}_${FIELD_STR}_${MESH_TYPE}_${RUN_ID}"
FINAL_PNS_NAME="Pns_I${CURRENT}_${FIELD_STR}.txt"

echo -e "${BLUE}[INFO]${NC} Output will be saved to: ${BOLD}$RUN_DIR${NC}\n"
mkdir -p $RUN_DIR

# 6. Inject Mesh settings into MeshManager.hpp automatically
# 6. Inject Mesh settings into MeshManager.hpp automatically
sed -i "s/ULLInt num_node = .*/ULLInt num_node = ${NUM_NODE};/g" src/MeshManager.hpp
sed -i "s/ULLInt num_tri = .*/ULLInt num_tri = ${NUM_TRI};/g" src/MeshManager.hpp
sed -i "s|std::string mesh_p_filename = .*|std::string mesh_p_filename = \"${MESH_P_PATH}\";|g" src/MeshManager.hpp
sed -i "s|std::string mesh_t_filename = .*|std::string mesh_t_filename = \"${MESH_T_PATH}\";|g" src/MeshManager.hpp

# 7. Compile
echo -e "${CYAN}⚙️  Compiling C++ solver...${NC}"
mpicxx -O3 src/main.cpp -o FP_mpi.o -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -lpetsc

if [ $? -ne 0 ] || [ ! -f "FP_mpi.o" ]; then
    echo -e "${RED}${BOLD}✖ FATAL ERROR: Compilation failed!${NC} Please check the C++ errors above."
    exit 1
fi
echo -e "${GREEN}✔ Compilation successful.${NC}\n"

# 8. Run
echo -e "${CYAN}🚀 Launching MPI Simulation on ${BOLD}${CORE_COUNT}${NC}${CYAN} cores...${NC}"
echo -e "${YELLOW}   (This may take a while depending on mesh size and total time. Tailing terminal_output.log...)${NC}"

export HYDRA_IFACE=lo
# Run PETSc
time mpirun --bind-to core -n $CORE_COUNT ./FP_mpi.o \
    -out_dir $RUN_DIR \
    -ksp_type gmres -pc_type asm -pc_asm_overlap 0 -sub_pc_type ilu -sub_pc_factor_levels 0 \
    -ksp_rtol 1e-5 -ts_rtol 1e-5 -ts_atol 1e-5 \
    -ts_adapt_type basic \
    -log_view ascii:${RUN_DIR}/profiling_log.txt > ${RUN_DIR}/terminal_output.log 2>&1

# 9. Rename the Pns file to reflect parameters
if [ -f "${RUN_DIR}/Pns_vs_Time.txt" ]; then
    mv "${RUN_DIR}/Pns_vs_Time.txt" "${RUN_DIR}/Pns.txt"
fi

cp src/Parameters.hpp $RUN_DIR/

# 10. Final Success Message
echo -e "\n${GREEN}${BOLD}✨ Sweep Complete!${NC}"
echo -e "------------------------------------------------------------"
echo -e "📄 ${BOLD}Results:${NC}          ${RUN_DIR}/${FINAL_PNS_NAME}"
echo -e "⏱️  ${BOLD}Profiling:${NC}        ${RUN_DIR}/profiling_log.txt"
echo -e "📋 ${BOLD}Terminal Output:${NC}  ${RUN_DIR}/terminal_output.log"
echo -e "------------------------------------------------------------\n"
