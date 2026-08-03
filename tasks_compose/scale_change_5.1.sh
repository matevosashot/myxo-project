#!/bin/bash

#SBATCH --exclusive

#SBATCH --array=0-0
#SBATCH --job-name="scale_change_5.1"

#SBATCH --mem=1600G
#SBATCH --partition=medium
#SBATCH --time=0-12:00:00 # expected maximum runtime of job
 
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --output=/home/ashmat/Projects/myxo-project/tasks/slurm_logs/%A_%a.out
# #SBATCH --error=/home/ashmat/Projects/myxo-project/tasks/slurm_logs/%A_%a.err

mkdir -p /home/ashmat/Projects/myxo-project/tasks/slurm_logs

set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

data="/data/biophys/ashmat/data/myxo-project"
source=/home/ashmat/Projects/myxo-project
export SCRATCH_PATH="/scratch/$USER/$job_id"
export JOB_ID=$job_id

mkdir -p $SCRATCH_PATH



# setup environment
module load python/3.14
cd /home/ashmat/Projects/myxo-project
source .venv/bin/activate

python tasks_compose/scale_change_5.1.py





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
