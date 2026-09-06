#!/bin/bash

# set3_2 -- set3_1 (sealed, anchoring wall) plus a PHYSICAL UV cutoff: the noise
# gets a Gaussian spatial correlation length lNoise.  Runs
# jobscripts/set3_2/fluctuations.wls, which loads
# lyapunov_solver_module_neumann_cutoff.m.
#
# The kernel W = MatrixExp[(lNoise^2/4) Lap] is built on the NEUMANN Laplacian, so
# it REFLECTS off the sealed wall -- a rod cannot straddle a wall it cannot cross.
# Its diagonal is therefore ENHANCED near the wall, toward a factor 2, where
# set2_2's Dirichlet kernel is suppressed to zero.  Verified against the closed
# form Sigma = ped [W^2 - 11^T/N^2] to 1.4e-14 in
# claude_experiments/fluctuations_neumann/cutoff_analytic_probe.wls.
#
# 54 tasks = 18 grids x 3 lNoise values.  Task -> (Nint, lNoise) is
#     iN = task / 3,  iL = task % 3
# so the three lNoise values of one grid are consecutive tasks.
#
# Output files carry BOTH the grid and a "_Neu_lN<value>" tag, so nothing
# collides.  Disk: the full covariance is ~219 GiB for ONE lNoise over all 18
# grids, so saveData is True only for lNoise = Sqrt[3] (iL = 1), the physical
# value; the other two write the small sigma_rho_*.m only.
#
# RESOLUTION WARNING, read this before interpreting the sweep.  The module needs
# lNoise >= 5 h for the continuum value to within 1%, and h = 40/Nint here, so
# 5 h = 200/Nint:
#     lNoise = 1        needs Nint >= 200  -> GRID-LIMITED AT EVERY GRID IN THE LIST
#     lNoise = Sqrt[3]  needs Nint >= 116  -> converged only for Nint = 130..161
#     lNoise = 5        needs Nint >= 40   -> converged for Nint = 40..161
# The lNoise = 1 column still shows the crossover into the grid-limited regime,
# which is worth having, but it is not a measurement of the model.  The module
# prints lNoiseOverH and issues ::subgrid whenever lNoise < 5 h.
#
# RESOURCES: as set3_1 -- 1600G, 72 cpus, array throttled to 10 concurrent.  See
# the header of fluctuation-solver-set3_1.sh for why those differ from set2_1's.

# #SBATCH --exclusive

#SBATCH --array=0-53%10
#SBATCH --job-name="set3_2"

#SBATCH --mem=1600G
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

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set3_2/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set3_2/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id


NintValues=(20 21 30 31 40 41 60 61 80 81 100 101 130 131 150 151 160 161)
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

wolframscript -file fluctuations.wls \
  "modelParams={lNoise->${lNoise}}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->20}; \
   otherParams={outputDir->\"${outputdir}\", saveData->${saveData}, plotBox->20}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, lNoise = ${lNoise}, saveData = ${saveData})"



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
