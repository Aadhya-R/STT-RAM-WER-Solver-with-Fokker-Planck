# A Highly Parallel Fokker–Planck Solver for STT-RAM

A high-performance, distributed-memory C++ solver for the Fokker–Planck (FP) equation, specifically engineered for estimating **Write Error Rates (WER)** in **Spin-Transfer Torque Random Access Memory (STT-RAM)**.

Unlike traditional **Stochastic Landau–Lifshitz–Gilbert (SLLG)** Monte Carlo simulations, this solver directly evolves the probability density function using a **Galerkin Finite Element Method (FEM)**, enabling significantly faster WER estimation while maintaining high physical fidelity.

---

##  Features

- **Distributed MPI Framework**
  - Built using **OpenMPI** and **PETSc**
  - Scales across multiple compute nodes
- **Automatic Domain Decomposition**
  - Uses **METIS** for mesh partitioning
- **Interactive Command Line Interface**
  - NPM-style guided setup
  - Automatic compilation
  - Easy local deployment
- **Automated Parameter Sweeping**
  - Batch execution across multiple operating conditions
  - Ideal for generating WER datasets and scalability benchmarks
---

# Prerequisites

The project uses **Conda** to manage all dependencies including PETSc, MPI, and C++ toolchains.

---

> **Where should these commands be run?**
>
> Execute the following commands in a Linux terminal on the machine where you intend to run the solver. This can be:
>
> - Your local Linux workstation
> - A Linux virtual machine (VM)
> - Windows using **WSL2 (Windows Subsystem for Linux)**
> - A remote Linux server (if Conda installations are permitted)
>
> The commands install Miniconda into your home directory (`$HOME/miniconda3`) and **do not require root (`sudo`) privileges**.
>
> **Note:** If you are using a managed HPC cluster, you typically **do not need to install Miniconda**. Instead, follow the instructions in **`README_HPC.md`**, which use the cluster's pre-installed software modules.

# Installation

## 1. Install Miniconda

Download and install Miniconda:

```bash
# Download installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
```

```bash
# Install
bash miniconda.sh -b -u -p $HOME/miniconda3
```
```bash
# Remove installer
rm miniconda.sh
```
```bash
# Activate Conda
source $HOME/miniconda3/bin/activate
```
```bash
# Initialize shell
conda init bash
```

> **Note**
>
> After running `conda init`, restart your terminal before continuing.

---

## 2. Clone the Repository

```bash
git clone https://github.com/Aadhya-R/STT-RAM-FP-Solver.git

cd STT-RAM-FP-Solver
```

---

## 3. Create the Conda Environment
Make sure all Terms of Service (TOS) are accepted 
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```
Then create an environment :
```bash
conda create \
    -n fp_env \
    -c conda-forge \
    petsc \
    mpi4py \
    cxx-compiler \
    -y
```
OR
```bash
conda create -n fp_env -c conda-forge petsc mpi4py cxx-compiler -y --solver classic
```

Activate the environment:

```bash
conda activate fp_env
```

---

# Usage

The solver supports two execution modes.

---

# Determining the Number of Physical CPU Cores

For best MPI performance, the solver should be launched using the number of **physical CPU cores**, **not** the number of logical processors (virtual threads).

## Number of actual physical CPU Cores

Run:

```bash
lscpu
```

Look for the following entries:

```text
Socket(s):              1
Core(s) per socket):    8
Thread(s) per core):    2
```

The number of **physical CPU cores** is calculated as:

```text
Physical Cores =
Socket(s) × Core(s) per socket
```

For the example above:

```text
1 × 8 = 8 Physical Cores
```

Use this value when selecting the number of MPI processes.

---

## Alternative Command to find number of Physical CPU Cores

The following command directly reports the number of physical cores:

```bash
lscpu | awk '/^Core\(s\) per socket:/ {cores=$4}
/^Socket\(s\):/ {sockets=$2}
END {print cores*sockets}'
```

Example output:

```text
8
```
---

# Mode 1 — Interactive Single Run

Ideal for testing a single device configuration.

### Step 1

Run 
```bash
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev
```

To Edit the Parameters , Open
```text
Parameters.hpp
```
and edit the required physical parameters.

Example:

```cpp
DRIVE_CURRENT
MAGNETIC_ANGLE
DC_FIELD
```

---

### Step 2

Run the setup script.

```bash
chmod +x setup.sh

./setup.sh
```

The interactive CLI will allow you to select:

- Mesh density
  - Coarse (35512 nodes and 71020 triangles)
  - Fine (878179 nodes and 1756354 triangles)

- Number of CPU cores

The script will then:

- Compile the solver
- Launch the MPI simulation
- Save results automatically

---

# Mode 2 — Automated Parameter Sweep

Ideal for generating datasets across multiple operating conditions.

Open:

```text
run_sweep.sh
```

Modify the configuration block:

```bash
TEST_CURRENTS=(
    1.5
    2.0
    3.0
    4.0
    6.0
)

TEST_ANGLES=(
    0
    45
    90
)
```

Execute:

```bash
chmod +x run_sweep.sh

./run_sweep.sh
```

The script automatically:

- Updates simulation parameters using `sed`
- Rebuilds the solver
- Executes every simulation
- Organizes outputs into timestamped folders

---

# Output

Each execution creates a directory similar to:

```text
STT_Run_20260704_143000/
```

Typical contents include:

```
STT_Run_20260704_143000/
│
├── profiling_log.txt
├── Pns_vs_Time.txt
└── terminal_output.out
```

Pns_vs_time.txt has the following :
(time,Pns,CheckNSum,PnsN,time_step pulse)

---

# Performance Profiling

The solver leverages PETSc's built-in profiling through:

```bash
-log_view
```

The generated `profiling_log.txt` contains detailed information on:

- Memory consumption
- MPI communication statistics
- Matrix assembly time
- Solver timings
- MFLOPs achieved
- Parallel efficiency

---

# Repository Structure

```
STT-RAM-FP-Solver/

├── src/
├── meshes/
├── Scripts/
├── setup.sh
├── run_sweep.sh
├── README_HPC.md
└── README.md
```

---

# Acknowledgments

This work was developed as part of the **SPARK Research Program** at the **Indian Institute of Technology (IIT) Roorkee**.

Special thanks to:

- **Prof. Tanmoy Pramanik** (https://scholar.google.com/citations?user=wl8G7iMAAAAJ&hl=en)
- **Sonalie Ahirwar** (https://scholar.google.com/citations?user=ff2E9KcAAAAJ&hl=en)
- **Susheel Kumar Arya** (https://scholar.google.com/citations?user=pSpyQYsAc1UC&hl=en)

for their invaluable mentorship, technical guidance, and continuous support throughout this research.

---

# Disclaimer — Numerical Artifacts

## Negative Probability Densities

The solver employs a **Galerkin Finite Element Method (FEM)** to solve the Fokker–Planck equation.

Under certain conditions—particularly:

- aggressive time stepping,
- highly advection-dominated regions, or
- extremely steep probability gradients,

small negative probability values may appear at isolated mesh nodes.

This is a well-known numerical artifact of Galerkin FEM and **does not affect** the overall Write Error Rate (WER) estimation or the physical validity of the simulation.

For visualization or post-processing, users may safely clamp the solution:

```cpp
W = std::max(W, 0.0);
```

before exporting VTK files or computing derived quantities.

---
# License

This project is released under the **MIT License**.
