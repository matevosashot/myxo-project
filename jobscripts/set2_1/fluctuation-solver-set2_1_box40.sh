#!/bin/bash

# set2_1 at BOX 40 -- LOW-Nint half (Nint <= 131).  Companion: _box40_hi.sh.
#
# Same physics and the same CONSISTENT stencils (R == 1) as the box-20 sweep, in a
# box twice as wide, to separate genuine physics from finite-box artefacts.
#
# WHY solverParams={box->40} IS MANDATORY, not cosmetic.  fbox is the Lyapunov
# grid half-width; solverParams box is the FEM domain half-width of the steady
# state (Rectangle[{-box,-box},{box,box}] in iterative_solver_module.m).  The
# Lyapunov module samples rho_ss/Q1_ss/Q2_ss AT ITS OWN GRID POINTS, i.e. out to
# +-fbox.  Setting fbox->40 while leaving box at its default 20 evaluates the FEM
# InterpolatingFunction far outside its ElementMesh, which extrapolates.  ALWAYS
# set both, equal.
#
# GRID NOTE.  h = 2 fbox/(Nint+1) here, so at box 40 the SAME Nint gives TWICE the
# h of the box-20 run -- the two sweeps are not the same measurement at matched
# Nint.  They are at matched h: box40 Nint = 2 x box20 Nint, i.e. the pairs
# 40<->20, 60<->30, 80<->40, 160<->80 are already in the list and isolate the pure
# box-size effect with no discretisation confound.
#
# THE lo/hi SPLIT, from measured peak RSS (sacct on jobscripts/logs/56583960*):
#     Nint = 100 -> 128 GB      Nint = 130 -> 365 GB      Nint = 131 -> 376 GB
# and from this set's own failures at 600G (56558539*): 150, 151, 160, 161 ALL
# died at the "[0] densify A,Q" line, i.e. out of memory inside the eigen-solve.
# So Nint <= 131 fits comfortably in 600G -- and 600G is schedulable on every one
# of the 267 medium-partition nodes with >= 72 cpus -- while only the four largest
# tasks need to queue for the 1600G of _box40_hi.sh.
#
# cpus-per-task 96 -> 72: only the ten amun nodes have 96 cpus and they sit in
# "long" alone, whereas 267 nodes in "medium" have 72+.
#
# Output goes to the SAME directory as the box-20 sweep: runTag is
# "box<fbox>_N<Nint>_R1", so every data file is distinct, and the steady-state
# PNGs are now tagged too (visualizeSteadyState's 4th argument).

# #SBATCH --exclusive

# LICENSE THROTTLE: the %4 suffix caps this array at 4 concurrent tasks.  set2_1's
# Nint = 131 task died with "The product exited because of a license error"
# (jobscripts/logs/56558539_13) when 18 kernels ran at once.  NB the cap is
# PER JOBSCRIPT, so submitting several of these together multiplies it.
# Array indices stay 0-BASED: task id indexes NintValues directly, so 1-N would
# silently skip Nint = 20.
#SBATCH --array=0-13%4
#SBATCH --job-name="s21b40"

#SBATCH --mem=600G
#SBATCH --partition=medium,long
#SBATCH --time=2-0:00:00 # expected maximum runtime of job

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72

#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out


set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set2_1/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_1/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id


# Nint <= 131 only; 150 151 160 161 live in _box40_hi.sh.
NintValues=(20 21 30 31 40 41 60 61 80 81 100 101 130 131)

Nint=${NintValues[$SLURM_ARRAY_TASK_ID]}

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

# plotBox stays 15, as in the box-20 sweep, so the figures are directly comparable.
wolframscript -file fluctuations.wls \
  "solverParams={box->40}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->40}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->15}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = 40)"



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
