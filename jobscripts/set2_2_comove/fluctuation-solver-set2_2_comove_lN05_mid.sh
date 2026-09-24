#!/bin/bash

# set2_2_comove, lNoise = 0.5 -- MID tier.  Companions: the other four tiers of
# fluctuation-solver-set2_2_comove_lN05_{lo,mid,hi,h2,xl}.sh.
#
# Nint = 100, 101 -- n = 30000 / 30603, ~122-127 GiB and ~1.7 h per task.
# 200G leaves ~60% headroom on the measured law.  Covariance files are ~7 GiB
# each, so this tier alone writes ~70 GiB.
#
# 10 tasks = 2 grids x 5 boxes (box = fbox = 10, 15, 25, 30, 40).  This is the
# coarsest tier that clears the lNoise = 0.5 calibration floor at every box,
# box 40 included (lNoise/h = 0.63 there, just above the 0.6 ::marginal line).
# ----------------------------------------------------------------------------
# WHAT THIS FAMILY IS
# ----------------------------------------------------------------------------
# lNoise = 0.5, the lNoise = 1.73 companion of the five
# fluctuation-solver-set2_2_comove_{lo,mid,hi,xl}.sh scripts in this same
# directory.  Same fluctuations.wls, same zeta = 8, same ell_d = 0.8, same
# comoving steady state and Lyapunov solve on the reservoir (Dirichlet) wall.
# ONLY the noise correlation length changes, 1.73 -> 0.5.
#
# Two things differ from the lNoise = 1.73 family beyond the noise length:
#   - box = 10 is added, so BoxList is (10 15 25 30 40) rather than
#     (15 25 30 40).  box 10 is where BOTH resolution criteria are most
#     comfortably met -- it is also the most wall-contaminated domain, since
#     the Dirichlet boundary sits at 10 rather than 40.  Read it as the
#     fine-grid anchor, not as the physical answer on its own.
#   - Nint = 140, 141 is added as a fifth tier (_h2), between _hi and _xl.
#
# Output goes to the SAME directory as the lNoise = 1.73 runs.  This cannot
# collide: the run tag carries the physical lNoise, so files land as
#   sigma_rho_box10_N141_Dir_lN0.5_z8_ld0.8.m
# alongside the _lN1.73_ ones.  Filter analysis on the tag, never on a glob.
#
# ----------------------------------------------------------------------------
# THE CALIBRATION FLOOR -- WHY box = 40 IS ABSENT FROM THE _lo TIER
# ----------------------------------------------------------------------------
# lyapunov_solver.m calibrates the cutoff by solving F(l/h) = 1/(Sqrt[2 Pi] t)
# with t = lNoise/h (dean.tex eq:lattice-bz).  That has a solution only for
#
#     lNoise >= h / Sqrt[2 Pi],        h = 2 box / (Nint + 1)
#
# (eq:ln-eff, recovered as the solvability condition).  Below it cutoffLength
# returns $Failed, SolveFluctuationsLyapunov issues ::uncalibratable and the
# task dies without writing anything.  At lNoise = 0.5 the floor is
# Nint >= Ceiling[2 box/(Sqrt[2 Pi] 0.5)] - 1, i.e.
#
#     box   10    15    25    30    40
#     Nint  15    23    39    47    63
#
# so Nint = 50, 51 clears every box EXCEPT 40, which needs 63.  Rather than
# ship two tasks that are certain to abort, the _lo tier drops box 40 and runs
# 2 x 4 = 8 tasks; the other four tiers run the full 2 x 5 = 10.  That is the
# only asymmetry in the set, and it is why nBox is read from BoxList rather
# than hard-coded.
#
# The ::marginal warning (lNoise < 0.6 h) still fires on _lo at box 25 and 30
# -- lNoise/h = 0.51 and 0.42.  Those tasks run and the pedestal is still
# calibrated, but they belong to the convergence curve, not to a measurement.
#
# lNoise/h over the whole set, for reference:
#
#          Nint:   50    51   100   101   130   131   140   141   150   151
#     box  10:   1.27  1.30  2.52  2.55  3.27  3.30  3.52  3.55  3.77  3.80
#     box  15:   0.85  0.87  1.68  1.70  2.18  2.20  2.35  2.37  2.52  2.53
#     box  25:   0.51* 0.52* 1.01  1.02  1.31  1.32  1.41  1.42  1.51  1.52
#     box  30:   0.42* 0.43* 0.84  0.85  1.09  1.10  1.18  1.18  1.26  1.27
#     box  40:    --    --   0.63  0.64  0.82  0.82  0.88  0.89  0.94  0.95
#                                                        (* = ::marginal)
#
# ----------------------------------------------------------------------------
# RESOURCES
# ----------------------------------------------------------------------------
# Peak RSS depends on n = 3 Nint^2 ALONE, not on the box, and follows this
# set's measured law peak GiB ~ 1.354e-7 n^2; the eigensolve is O(n^3),
# anchored at 199 s for n = 11163 with the rest of the pipeline adding ~50%.
# The five tiers are split so the cheap grids do not queue behind a 2000G
# reservation.  Adding box 10 does not change peak RSS, only wall time --
# the steady solve is cheaper on a smaller domain.
#
# LICENSE THROTTLE: the % suffix on --array caps concurrency.  A Nint = 131
# task once died with "The product exited because of a license error" when 18
# kernels ran at once.  NB the cap is PER JOBSCRIPT, so submitting several of
# these together multiplies it.

# #SBATCH --exclusive

#SBATCH --array=0-9%4
#SBATCH --job-name="s22c05_mid"

#SBATCH --mem=200G
#SBATCH --partition=medium,long
#SBATCH --time=0-8:00:00 # expected maximum runtime of job

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
NintList=(100 101)
BoxList=(10 15 25 30 40)

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
  "modelParams={lNoise->0.5}; \
   solverParams={box->${box}.}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->${box}}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->${plotBox}}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = fbox = ${box}, lNoise = 0.5, zeta = 8, comoving, Dirichlet)"


# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
