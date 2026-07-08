#!/bin/bash
#SBATCH --job-name=poly_top
#SBATCH --output=poly_top_%j.out
#SBATCH --error=poly_top_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
# TODO: add your usual '#SBATCH --partition=...' (and '--account=...' if you
#       use one) from your existing macro job script.

# TODO: copy the 'module load ...' / conda activation lines from your
#       existing macro job script here (whatever environment you compile
#       babbler with on the cluster).

# The polynomial model has no MultiApp, so a plain MPI run is fine.
# Recompile babbler on the cluster first:
#   cd ~/babbler && git pull && make -j8

BABBLER=$HOME/babbler/babbler-opt

# --------------------------------------------------------------------------
# IMPORTANT: at the full fast-cell current (I=1.108) the parabolic model is
# infeasible (demanded j_n ~ 25x the closure cap) and the solve stalls by
# construction — see the header of macro_poly_top_1c.i. The overrides below
# run the validated operating point (I_top = 0.0222, ~1/50 of fast-cell 1C).
# For a C-rate sweep, scale all three I values together and keep
# I_catcc = 0.4249 * I_top.
# --------------------------------------------------------------------------
srun $BABBLER -i macro_poly_top_1c.i \
  BCs/flux_c/I=0.0222 \
  BCs/flux_phi1/I=0.0094 \
  BCs/flux_phi2/I=0.0222 \
  Outputs/file_base=macro_poly_lowI_out
