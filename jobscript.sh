#!/bin/bash


# icarus  96 2.3 TB

# seth 8 6TB

# amun 24 4TB

# brahe 36 4TB

# flora 64 4TB

# hermes 36 2 TB





#SBATCH --time=1-00:00:00 # expected maximum runtime of job
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --partition=long

#### ## # SBATCH --nodelist=snowy03
#SBATCH --cpus-per-task=64 # every task gets 2 CPUs
#SBATCH --mem=1000G

#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END # get notification to mail
#SBATCH --job-name="reserve"
#SBATCH --output=/home/ashmat/Projects/myxo-project/reserve%j.out
#SBATCH --error=/home/ashmat/Projects/myxo-project/reserve%j.err

# just infinitely running script to reserve the node for testing
# the node will be reserved for 14 days, but we can cancel the job when done

while true; do
    sleep 60
done
