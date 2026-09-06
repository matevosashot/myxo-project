#!/bin/bash

# set2_2 -- same physics as set-2/set2_1, consistent stencils (R == 1), plus a
# PHYSICAL UV cutoff: the noise gets a Gaussian spatial correlation length lNoise.
# Runs jobscripts/set2_2/fluctuations.wls, which loads
# lyapunov_solver_module_cutoff.m.
#
# 90 tasks = 18 grids x 5 lNoise values.  Task -> (Nint, lNoise) is
#     iN = task / 5,  iL = task % 5
# so the five lNoise values of one grid are consecutive tasks.
#
# Output files carry BOTH the grid and an "_lN<value>" tag, so nothing collides.
# Disk: the full covariance is 218 GiB for ONE lNoise over all 18 grids, so
# saveData is True only for lNoise = Sqrt[3] (iL = 3), the physical value; the
# other four write the small sigma_rho_*.m only.

# #SBATCH --exclusive

#SBATCH --array=0-89
#SBATCH --job-name="set2_2"

#SBATCH --mem=600G
#SBATCH --partition=medium,long
#SBATCH --time=2-0:00:00 # expected maximum runtime of job

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96


#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out


set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id


NintValues=(20 21 30 31 40 41 60 61 80 81 100 101 130 131 150 151 160 161)
lNoiseValues=(0.1 0.5 1 "Sqrt[3.]" 5)

iN=$(( SLURM_ARRAY_TASK_ID / 5 ))
iL=$(( SLURM_ARRAY_TASK_ID % 5 ))

Nint=${NintValues[$iN]}
lNoise=${lNoiseValues[$iL]}

# the full covariance only for the physical lNoise = Sqrt[3]; all five always
# write sigma_rho_*.m
if [ "$iL" -eq 3 ]; then saveData=True; else saveData=False; fi

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

wolframscript -file fluctuations.wls \
  "modelParams={lNoise->${lNoise}}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->20}; \
   otherParams={outputDir->\"${outputdir}\", saveData->${saveData}, plotBox->15}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, lNoise = ${lNoise}, saveData = ${saveData})"



# icarus  96 2.3 TB

# seth 8 6TB

# amun 24 4TB

# brahe 36 4TB

# flora 64 4TB

# hermes 36 2 TB



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# graphic	    2-00:00:00
# extra_long	28-00:00:00




#### ## # SBATCH --nodelist=snowy03


# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END # get notification to mail
