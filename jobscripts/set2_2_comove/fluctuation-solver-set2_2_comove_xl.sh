#!/bin/bash

# set2_2_comove -- XL tier.  Companions: the other three tiers of
# fluctuation-solver-set2_2_comove_{lo,mid,hi,xl}.sh.
#
# Nint = 150, 151 -- n = 67500 / 68403, ~617-634 GiB and ~19 h per task.
#
# 2000G, not 1600G: 1600G is recorded in set2_2's _boxscan.sh as having been too
# tight for grids this size, and with only 8 tasks here there is no reason to
# economise.  3-0 rather than 2-0 for the same reason -- the n^3 estimate is a
# scaling from n = 11163, not a measurement at this size.
#
# THESE GRIDS CARRY THE ANSWER.  Both resolution criteria are met only here: the
# defect core ell_d = 0.8 reaches 4.1 h at box = 15, and the cutoff lNoise = 1.73
# reaches 8.8 h.  The other 24 tasks are the convergence curve.  Covariance files
# are ~35 GiB each, ~280 GiB for this tier.
#
# 8 tasks = 2 grids x 4 boxes (box = fbox = 15, 25, 30, 40), lNoise = 1.73,
# zeta = 8, comoving steady state and Lyapunov solve both on the reservoir
# (Dirichlet) wall.  See fluctuations.wls in this directory for what changed
# relative to the superseded zeta = 10 runs -- above all that the drift now
# CARRIES the comoving advection Div(dQ w_ss) that the old module dropped.
#
# Resource classes are split across four scripts (_lo, _mid, _hi, _xl) so that
# the cheap grids do not sit in the queue behind a 2000G reservation.  Peak RSS
# depends on n = 3 Nint^2 ALONE, not on the box, and follows this set's measured
# law peak GiB ~ 1.354e-7 n^2; the eigensolve is O(n^3), anchored at 199 s for
# n = 11163 with the rest of the pipeline adding ~50%.
#
# LICENSE THROTTLE: the %2 suffix caps this array at 2 concurrent tasks.  A
# Nint = 131 task once died with "The product exited because of a license error"
# when 18 kernels ran at once.  NB the cap is PER JOBSCRIPT, so submitting
# several of these together multiplies it.

# #SBATCH --exclusive

#SBATCH --array=0-7%2
#SBATCH --job-name="s22c_xl"

#SBATCH --mem=2000G
#SBATCH --partition=long
#SBATCH --time=3-0:00:00 # expected maximum runtime of job

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72

#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out

set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2_comove/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2_comove/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id

# One array id indexes the (Nint, box) product: Nint varies slowest, box fastest.
# Keep the id 0-BASED -- 1-N would silently skip the first Nint.
NintList=(150 151)
BoxList=(15 25 30 40)

i=${SLURM_ARRAY_TASK_ID}
nBox=${#BoxList[@]}
Nint=${NintList[$((i / nBox))]}
box=${BoxList[$((i % nBox))]}

if [ -z "$Nint" ] || [ -z "$box" ]; then
  echo "ERROR: array task $i is out of range; this script has ${#NintList[@]} x $nBox = $(( ${#NintList[@]} * nBox )) tasks"
  exit 1
fi

# Crop the Dirichlet boundary layer out of the figures.
plotBox=$(( box - 3 ))

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

# box is passed TWICE and must agree: solverParams box is the Real the comoving
# steady solver takes, lyapunovSolverParams fbox the Integer the Lyapunov module
# takes (it tests IntegerQ).  The driver aborts if they differ.
wolframscript -file fluctuations.wls \
  "modelParams={lNoise->1.73}; \
   solverParams={box->${box}.}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->${box}}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->${plotBox}}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = fbox = ${box}, lNoise = 1.73, zeta = 8, comoving, Dirichlet)"


# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
