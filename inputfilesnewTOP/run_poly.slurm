#!/bin/bash
#SBATCH --job-name=poly_top
#SBATCH --output=poly_top_%j.out
#SBATCH --error=poly_top_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
# TODO: add your usual '#SBATCH --partition=...' (and '--account=...' if you
#       use one) from your existing macro job script.

# Build/runtime environment (Lichtenberg, same as used to compile babbler).
# NOTE: mpich, NOT openmpi — the cluster PETSc/MOOSE was configured with
# MPICH and mixing them fails with "PETSc was configured with MPICH but..."
module purge
module load gcc/11.5.0-z7mc
module load mpich/4.3.2-hca2
module load cmake/3.31.9-btkq
module load python/3.11.14-xvq4
module load python-venv/1.0-py_3.11-qcz3
source ~/moose_env/bin/activate

# The polynomial model has no MultiApp, so a plain MPI run is fine.
# Recompile babbler on the cluster first:
#   cd ~/babbler && git pull && make -j8

BABBLER=$HOME/babbler/babbler-opt

# --------------------------------------------------------------------------
# The input's default currents (I_top = 0.03, ~TRUE 0.2C) are validated and
# converge robustly. Above I_top ~ 0.045 the parabolic closure is infeasible
# (see the input header) — do NOT raise to the old fast-cell 1.108.
# For a C-rate sweep, override all three I values together and keep
# I_catcc = 0.4249 * I_top, e.g.:
#   BCs/flux_c/I=0.02 BCs/flux_phi1/I=0.008498 BCs/flux_phi2/I=0.02
# --------------------------------------------------------------------------
srun $BABBLER -i macro_poly_top_1c.i \
  Outputs/file_base=macro_poly_out
