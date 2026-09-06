#!/bin/bash

# set3_1 at BOX 40 -- LOW-Nint half (Nint <= 131).  Companion: _box40_hi.sh.
#
# The SEALED-WALL set in a box twice as wide.  This is the set where box size
# matters most: the wall IS the observable, so how far the wall sits from the
# defect core is a physical parameter, not a convergence nuisance.
#
#   set2_1  reservoir : rho_ss pinned,     psi|_bdry = 0,      noise flux FREE
#   set3_1  sealed    : n.Grad rho_ss = 0, n.Grad d(rho) = 0,  noise flux BLOCKED
#
# Runs jobscripts/set3_1/fluctuations.wls, which loads
# lyapunov_solver_module_neumann.m and calls the steady state with "Neumann".
#
# WHY solverParams={box->40} IS MANDATORY, not cosmetic.  fbox is the Lyapunov
# grid half-width; solverParams box is the FEM domain half-width of the steady
# state (Rectangle[{-box,-box},{box,box}] in iterative_solver_module.m).  The
# Lyapunov module samples rho_ss/Q1_ss/Q2_ss AT ITS OWN GRID POINTS, i.e. out to
# +-fbox.  Setting fbox->40 while leaving box at its default 20 evaluates the FEM
# InterpolatingFunction far outside its ElementMesh, which extrapolates -- and here
# it would put the sealed wall of the fluctuation problem 20 units OUTSIDE the wall
# of the steady state it is linearized about.  ALWAYS set both, equal.
#
# GRID NOTE.  The sealed scheme is cell-centred with h = 2 fbox/Nint, so Nint here
# has EXACTLY the same h as Nint-1 in set2_1 -- the shared NintValues list supplies
# an h-matched Dirichlet partner for every odd run: 21<->20, 31<->30, and so on.
# Keep the list identical for that reason.  Separately, at box 40 the same Nint
# gives TWICE the h of the box-20 run, so the two sweeps are h-matched at
# box40 Nint = 2 x box20 Nint: 40<->20, 60<->30, 80<->40, 160<->80.
#
# THE lo/hi SPLIT, from measured peak RSS (sacct on this set's own box-20 run,
# jobscripts/logs/56583960*):
#     Nint = 100 -> 128 GB      Nint = 130 -> 365 GB      Nint = 131 -> 376 GB
# at 6 h 21 m for Nint = 131.  Only 150-161 need more than 600G, and 600G is
# schedulable on every one of the 267 medium-partition nodes with >= 72 cpus,
# whereas 1600G excludes the 91-node 1484G group in "long".

# #SBATCH --exclusive

# LICENSE THROTTLE: the %4 suffix caps this array at 4 concurrent tasks.  set2_1's
# Nint = 131 task died with "The product exited because of a license error"
# (jobscripts/logs/56558539_13) when 18 kernels ran at once.  NB the cap is
# PER JOBSCRIPT, so submitting several of these together multiplies it.
# Array indices stay 0-BASED: task id indexes NintValues directly, so 1-N would
# silently skip Nint = 20.
#SBATCH --array=0-13%4
#SBATCH --job-name="s31b40"

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

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set3_1/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set3_1/"

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

# plotBox = fbox: with sealed walls the boundary IS the point of the run, so
# nothing is cropped out of the figures.
wolframscript -file fluctuations.wls \
  "solverParams={box->40}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->40}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->40}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = 40)"



# icarus  96 2.3 TB

# amun 24 4TB   <- actually 96 cpus / 4.1 TB, partition "long"

# brahe 36 4TB  <- actually 72 cpus / 4.1 TB, partition "medium"

# eos01 72 cpus / 2 TB, feature "interactive" -- used for the validation runs



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
