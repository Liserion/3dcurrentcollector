#!/bin/bash
#SBATCH -J poly02c                         # Job name
#SBATCH --mail-type=END                    # Email notification when job ends
#SBATCH -e babbler/%x.err.%j
#SBATCH -o babbler/%x.out.%j
#SBATCH --nodes=1
#SBATCH -n 8                               # Number of tasks (cores)
#SBATCH --mem-per-cpu=16G                  # Memory per CPU
#SBATCH -t 24:00:00                        # Max runtime
#SBATCH -A project02600

# Environment for the NEW moose stack (~/projects/moose): gcc 11.5 + MPICH.
# NOT openmpi — that belongs to the old ~/MOOSE stack; the new PETSc/libmesh
# were configured with MPICH and mixing MPIs fails at compile AND runtime.
module purge
module load gcc/11.5.0-z7mc
module load mpich/4.3.2-hca2
module load python/3.11.14-xvq4
module load python-venv/1.0-py_3.11-qcz3
source ~/moose_env/bin/activate

# Ensure gcc-11.5's libstdc++ wins over the old /shared/apps/gcc/10.2.0 one,
# otherwise: "version GLIBCXX_3.4.29 not found"
GCCLIB=$(dirname $(dirname $(which gcc)))/lib64
export LD_LIBRARY_PATH=$GCCLIB:$LD_LIBRARY_PATH

cd /work/scratch/jk98dago

# --------------------------------------------------------------------------
# The input's default currents (I_top = 0.03, ~TRUE 0.2C) are validated and
# converge robustly. Above I_top ~ 0.045 the parabolic closure is infeasible
# (see the input header) — do NOT raise to the old fast-cell 1.108.
# For a C-rate sweep, override all three I values together and keep
# I_catcc = 0.4249 * I_top, e.g.:
#   BCs/flux_c/I=0.02 BCs/flux_phi1/I=0.008498 BCs/flux_phi2/I=0.02
# --------------------------------------------------------------------------
srun $HOME/projects/babbler/babbler-opt \
  -i ./3dcurrentcollector/inputfilesnewTOP/macro_poly_top_1c.i \
  --disable-perf-graph-live --no-timing
EXITCODE=$?

exit $EXITCODE
