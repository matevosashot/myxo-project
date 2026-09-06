#!/bin/bash

# set3_2 at BOX 40 -- HIGH-Nint half (Nint > 131).  Companion: _box40.sh.
#
# Identical to _box40.sh except for the grid list and --mem.  See that file's
# header for why solverParams={box->40} must accompany fbox->40, for the l/h
# resolution table, and for the RSS measurements behind the lo/hi split.
#
# 12 tasks = 4 grids x 3 lNoise:  iN = task / 3,  iL = task % 3.
#
# 1600G, as in the box-20 run of this set: at 600G set2_1's four largest tasks all
# died at the "[0] densify A,Q" line (jobscripts/logs/56558539_14,15,16,17), while
# at 1600G this set's large tasks ran through.  Extrapolating the measured 376 GB
# at Nint = 131 by n^2 (n = 3 Nint^2) puts Nint = 161 near 860 GB.
#
# These are the best-resolved grids of the sweep: at Nint = 161, lNoise = Sqrt[3]
# reaches 3.5 h (~2% bias) and lNoise = 5 reaches 10 h (0.25%).
#
# TIME.  Nint = 131 took ~7 h at 72 cpus; the eigen-solve is ~n^3, so Nint = 161
# should land near 21 h -- inside the 2-day limit, but not by a wide margin.  If
# the box-20 tasks turn out to have hit the wall, raise this to --time=4-0:00:00
# and drop "medium" from --partition (medium caps at 2 days; "long" allows 14 and
# still has 120 nodes with >= 72 cpus and >= 1600G).

# #SBATCH --exclusive

# LICENSE THROTTLE: the %4 suffix caps this array at 4 concurrent tasks.  set2_1's
# Nint = 131 task died with "The product exited because of a license error"
# (jobscripts/logs/56558539_13) when 18 kernels ran at once.  NB the cap is
# PER JOBSCRIPT, so submitting several of these together multiplies it.
# Array indices stay 0-BASED: task id indexes NintValues directly, so 1-N would
# silently skip Nint = 20.
#SBATCH --array=0-11%4
#SBATCH --job-name="s32b40hi"

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


NintValues=(150 151 160 161)
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

# plotBox = fbox: with sealed walls the boundary IS the point of the run.
wolframscript -file fluctuations.wls \
  "modelParams={lNoise->${lNoise}}; \
   solverParams={box->40}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->40}; \
   otherParams={outputDir->\"${outputdir}\", saveData->${saveData}, plotBox->40}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint}, box = 40, lNoise = ${lNoise}, saveData = ${saveData})"



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
