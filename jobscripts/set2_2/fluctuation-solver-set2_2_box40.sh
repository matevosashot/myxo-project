#!/bin/bash

# set2_2 at BOX 40 -- LOW-Nint half (Nint <= 131).  Companion: _box40_hi.sh.
#
# Same physics as the box-20 sweep -- consistent stencils (R == 1) plus a PHYSICAL
# UV cutoff, the noise carrying a Gaussian spatial correlation length lNoise -- in
# a box twice as wide.  Runs jobscripts/set2_2/fluctuations.wls, which loads
# lyapunov_solver_module_cutoff.m.
#
# WHY solverParams={box->40} IS MANDATORY, not cosmetic.  fbox is the Lyapunov
# grid half-width; solverParams box is the FEM domain half-width of the steady
# state (Rectangle[{-box,-box},{box,box}] in iterative_solver_module.m).  The
# Lyapunov module samples rho_ss/Q1_ss/Q2_ss AT ITS OWN GRID POINTS, i.e. out to
# +-fbox.  Setting fbox->40 while leaving box at its default 20 evaluates the FEM
# InterpolatingFunction far outside its ElementMesh, which extrapolates.  ALWAYS
# set both, equal.
#
# THREE lNoise VALUES, NOT FIVE.  The box-20 script swept (0.1 0.5 1 Sqrt[3] 5);
# this one keeps (1 Sqrt[3] 5), matching set3_2.  42 tasks = 14 grids x 3 lNoise:
#     iN = task / 3,  iL = task % 3
# so the three lNoise values of one grid are consecutive tasks.  NB saveData moves
# from iL == 3 (the 5-value layout) to iL == 1.
#
# RESOLUTION, read this before interpreting the sweep.  The kernel is built on the
# DISCRETE Laplacian, whose mode weight Exp[-(l/h)^2 Sin[phi/2]^2] decays more
# slowly than the continuum Exp[-(l/h)^2 (phi/2)^2] because Sin[phi/2]^2 <
# (phi/2)^2.  It therefore keeps slightly too much short-wavelength noise, with a
# one-sided error of +(1/4)(h/l)^2 that depends on l/h ALONE: 1.0% at l = 5h, which
# is where the module's "lNoise >= 5 h" advice and its ::subgrid message come from.
# It is an accuracy tolerance, not a validity threshold.  At box 40, h = 80/(Nint+1)
# is twice the box-20 value, so the error is 4x -- at Nint = 131:
#     lNoise = 1        1.65 h  ->  9%    grid-limited, shows the crossover only
#     lNoise = Sqrt[3]  2.86 h  ->  3%    usable, with a known one-sided bias
#     lNoise = 5        8.25 h  -> 0.4%   converged
# Reaching 1% for Sqrt[3] would need Nint >= 231, which the dense eigen-solve
# cannot afford (n = 3 Nint^2, already ~860 GB at Nint = 161).  The clean
# comparison instead uses the h-MATCHED pairs box40 Nint = 2 x box20 Nint --
# 40<->20, 60<->30, 80<->40, 160<->80 -- where h and l/h are identical and the
# difference is pure box size.  The module prints lNoiseOverH for every run.
#
# THE lo/hi SPLIT, from measured peak RSS (sacct on jobscripts/logs/56583960*):
#     Nint = 100 -> 128 GB      Nint = 130 -> 365 GB      Nint = 131 -> 376 GB
# and from set2_1's failures at 600G (56558539*): 150, 151, 160, 161 ALL died at
# the "[0] densify A,Q" line.  So Nint <= 131 fits comfortably in 600G -- and 600G
# is schedulable on every one of the 267 medium-partition nodes with >= 72 cpus --
# while only the largest grids need the 1600G of _box40_hi.sh.
#
# cpus-per-task 96 -> 72: only the ten amun nodes have 96 cpus and they sit in
# "long" alone, whereas 267 nodes in "medium" have 72+.
#
# Output goes to the SAME directory as the box-20 sweep: runTag is
# "box<fbox>_N<Nint>_lN<value>", so every data file is distinct.  Disk: the full
# covariance is ~200 GiB for ONE lNoise over a full grid list, so saveData is True
# only for lNoise = Sqrt[3] (iL = 1), the physical value; the other two write the
# small sigma_rho_*.m only.

# #SBATCH --exclusive

# LICENSE THROTTLE: the %4 suffix caps this array at 4 concurrent tasks.  set2_1's
# Nint = 131 task died with "The product exited because of a license error"
# (jobscripts/logs/56558539_13) when 18 kernels ran at once.  NB the cap is
# PER JOBSCRIPT, so submitting several of these together multiplies it.
# Array indices stay 0-BASED: task id indexes NintValues directly, so 1-N would
# silently skip Nint = 20.
#SBATCH --array=0-41%4
#SBATCH --job-name="s22b40"

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

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id


# Nint <= 131 only; 150 151 160 161 live in _box40_hi.sh.
NintValues=(20 21 30 31 40 41 60 61 80 81 100 101 130 131)
lNoiseValues=(1 "Sqrt[3.]" 5)

iN=$(( SLURM_ARRAY_TASK_ID / 3 ))
iL=$(( SLURM_ARRAY_TASK_ID % 3 ))

Nint=${NintValues[$iN]}
lNoise=${lNoiseValues[$iL]}

# the full covariance only for the physical lNoise = Sqrt[3]; all three always
# write sigma_rho_*.m
if [ "$iL" -eq 1 ]; then saveData=True; else saveData=False; fi

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

# plotBox stays 15, as in the box-20 sweep, so the figures are directly comparable.
wolframscript -file fluctuations.wls \
  "modelParams={lNoise->${lNoise}}; \
   solverParams={box->40}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->40}; \
   otherParams={outputDir->\"${outputdir}\", saveData->${saveData}, plotBox->15}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = 40, lNoise = ${lNoise}, saveData = ${saveData})"



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
