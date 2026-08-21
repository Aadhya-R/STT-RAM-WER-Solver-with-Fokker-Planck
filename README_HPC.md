# Running on an HPC Cluster

This document describes how to compile and execute the STT-RAM Fokker–Planck Solver on a managed HPC cluster.

Unlike the local installation, all required software is provided through the cluster's module system.

---

# Prerequisites

Ensure that your HPC account has access to the required software modules.

The solver has been validated using:

- PETSc 3.20
- OpenMPI
- GCC/G++

---

# Loading Required Modules

Load the required software stack before compiling or submitting jobs.

```bash
module load petsc/3.20
module load openmpi
module load gcc
```

To verify the loaded modules:

```bash
module list
```

---

# Project Layout

Navigate to the project directory.

```bash
cd STT-RAM-FP-Solver
```

The HPC execution scripts are located in:

```text
hpc/
```

---

# Submitting a Job

Move into the HPC directory.

```bash
cd hpc
```

Submit the SLURM job using:

```bash
sbatch submit_HPC.sh
```

---

# Monitoring Jobs

Check the queue:

```bash
squeue -u $USER
```

View detailed job information:

```bash
scontrol show job <JOB_ID>
```

Cancel a running job if necessary:

```bash
scancel <JOB_ID>
```

---

# Output

SLURM automatically generates output files similar to:

```text
slurm-123456.out
```

Simulation results, VTK files, and solver logs are written to the output directory specified in the submission script.

---

# Notes

- Do not execute the solver directly using `mpirun`.
- Always submit jobs through SLURM using `sbatch`.
- The requested number of MPI ranks, CPUs, memory, and wall time can be modified inside `submit_HPC.sh` to match your allocation and workload.

---

# Example Workflow

```bash
git clone https://github.com/YOUR_USERNAME/STT-RAM-FP-Solver.git

cd STT-RAM-FP-Solver

module load petsc/3.20
module load openmpi
module load gcc

cd hpc

sbatch submit_HPC.sh
```
