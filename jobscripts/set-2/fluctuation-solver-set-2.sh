#!/bin/bash

# #SBATCH --exclusive

#SBATCH --array=0-17
#SBATCH --job-name="set 2"

#SBATCH --mem=600G
#SBATCH --partition=medium,long
#SBATCH --time=2-0:00:00 # expected maximum runtime of job

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96


#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out


mkdir -p /home/ashmat/Projects/myxo-project/jobscripts/slurm_logs

set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

source=/home/ashmat/Projects/myxo-project-dean
setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set-2/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set-2/"

export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id


NintValues=(20 21 30 31 40 41 60 61 80 81 100 101 130 131 150 151 160 161)

Nint=${NintValues[$SLURM_ARRAY_TASK_ID]}

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

wolframscript -file fluctuations.wls \
  "lyapunovSolverParams={Nint->${Nint}, fbox->20}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->15}"



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