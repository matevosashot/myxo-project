#!/bin/bash


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

#SBATCH --array=0-0

#SBATCH --job-name="ust"
#SBATCH --partition=medium
#SBATCH --time=0-03:00:00 # expected maximum runtime of job
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)

#### ## # SBATCH --nodelist=snowy03
# #SBATCH --cpus-per-task=64
# #SBATCH --mem=100000 # 100GB
#SBATCH --exclusive
#SBATCH --mem=1600G

#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END # get notification to mail

#SBATCH --output=/home/ashmat/Projects/myxo-project/tasks/slurm_logs/%j.out
#SBATCH --error=/home/ashmat/Projects/myxo-project/tasks/slurm_logs/%j.err

set -e
set -x

START=$(date +%s.%N)

job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

data="/data/biophys/ashmat/data/myxo-project"
source=/home/ashmat/Projects/myxo-project
scratch="/scratch/$USER/$job_id"


mkdir -p $scratch



# setup environment
module load python/3.14
cd /home/ashmat/Projects/myxo-project
source .venv/bin/activate

cd $scratch

toolbox worker \
    --task-base-path $source/tasks/ \
    --worker-name ust$job_id \
    --loop \
    --no-random 
