#!/bin/bash

# set3_1 -- the SEALED-WALL counterpart of set2_1.
#
# Same physics and the same consistent stencils (R == 1), but the other consistent
# boundary condition of dean.tex sec-bc: the wall blocks mass instead of holding
# the density at a far-field value.
#
#   set2_1  reservoir : rho_ss pinned,     psi|_bdry = 0,      noise flux FREE
#   set3_1  sealed    : n.Grad rho_ss = 0, n.Grad d(rho) = 0,  noise flux BLOCKED
#
# The wall also ANCHORS the director (dQ|_bdry = 0), exactly as the steady state
# does.  Sealing constrains the rho sector only, because only rho is conserved;
# giving dQ Neumann as well leaves the +1/2 defect unconfined and the
# linearization carries three unstable modes that survive refinement.  See point
# (5) of the lyapunov_solver_module_neumann.m header.
#
# Runs jobscripts/set3_1/fluctuations.wls, which loads
# lyapunov_solver_module_neumann.m and calls the steady state with "Neumann".
# Output files carry a "_Neu" tag and go to their own directory, so nothing can
# collide with fluctuations-set-2, -set2_1 or -set2_2.
#
# GRID NOTE.  The sealed scheme is cell-centred with h = 2 fbox/Nint, so Nint here
# has EXACTLY the same h as Nint-1 in set2_1.  The shared NintValues list therefore
# supplies an h-matched Dirichlet partner for every odd run: 21<->20, 31<->30, and
# so on.  Keep the list identical for that reason.
#
# RESOURCES, changed from set2_1 on the evidence of its own logs.  set2_1 asked for
# 600G on 96 cpus; its four largest tasks (Nint = 150, 151, 160, 161) all died --
# three of them at the "[0] densify A,Q  ~90 GiB" line, i.e. out of memory inside
# Eigensystem, and one on a license error.  At Nint = 161 (n = 77763) the solve
# holds A and Q dense (2 x 48 GB), the complex eigenvector matrix and its transpose
# (2 x 97 GB), lam_i + lam_j (97 GB) and several transient n x n complex temporaries
# -- roughly 700 GB peak.  Hence:
#   --mem raised 600G -> 1600G;
#   --cpus-per-task lowered 96 -> 72, because ONLY the ten amun nodes have 96 cpus
#     and they sit in "long" alone, whereas 267 nodes in "medium" have 72+;
#   --array throttled to 10 concurrent, since 18 simultaneous kernels is what
#     produced set2_1's license failure.

# #SBATCH --exclusive

#SBATCH --array=0-17%10
#SBATCH --job-name="set3_1"

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

setdir="/home/ashmat/Projects/myxo-project-dean/jobscripts/set3_1/"
outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set3_1/"

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

# plotBox = fbox: with sealed walls the boundary IS the point of the run, so
# nothing is cropped out of the figures.
wolframscript -file fluctuations.wls \
  "lyapunovSolverParams={Nint->${Nint}, fbox->20}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, plotBox->20}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (Nint = ${Nint})"



# icarus  96 2.3 TB

# seth 8 6TB

# amun 24 4TB   <- actually 96 cpus / 4.1 TB, partition "long"

# brahe 36 4TB  <- actually 72 cpus / 4.1 TB, partition "medium"

# flora 64 4TB

# hermes 36 2 TB

# eos01 72 cpus / 2 TB, feature "interactive" -- used for the validation runs



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# graphic	    2-00:00:00
# extra_long	28-00:00:00
