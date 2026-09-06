#!/bin/bash

# set3_2 BOX SCAN -- one (fbox, Nint) pair per array task, lNoise = Sqrt[3] throughout.
#
# Unlike the _box40 scripts, which fix the box and sweep the grid, this one varies
# BOTH so that box size and mesh spacing can be separated.  Sealed, anchoring wall
# (Neumann); runs jobscripts/set3_2/fluctuations.wls ->
# lyapunov_solver_module_neumann_cutoff.m.
#
# THE ROW TABLE, ORDERED BY Nint so the cheap rows land first.  Cost is ~n^3 with
# n = 3 Nint^2 and does not depend on fbox, so tasks 0-3 finish in a few hours
# while tasks 6-7 take most of a day; with the %4 throttle the whole Nint = 120/121
# block clears before anything expensive starts.
#
# h = 2 fbox/Nint for this cell-centred scheme, and the noise correlation length is
# fixed at lNoise = Sqrt[3] = 1.7321, so l/h is what decides how well the Gaussian
# cutoff is resolved.  The module's accuracy rule is +(1/4)(h/l)^2, i.e. 1% at
# l = 5h (see the _box40.sh header for the derivation):
#
#   task  fbox  Nint     h       l/h    kernel err   n = 3 Nint^2   peak RSS
#   ----  ----  ----  -------  ------  -----------  -------------  ---------
#     0     60   120   1.0000   1.732      8.3%         43200        253 GiB
#     1    120   120   2.0000   0.866     33.3%         43200        253 GiB
#     2     60   121   0.9917   1.746      8.2%         43923        261 GiB
#     3    120   121   1.9835   0.873     32.8%         43923        261 GiB
#     4     40   140   0.5714   3.031      2.7%         58800        468 GiB
#     5     40   141   0.5674   3.053      2.7%         59643        482 GiB
#     6     75   150   1.0000   1.732      8.3%         67500        617 GiB
#     7     75   151   0.9934   1.744      8.2%         68403        634 GiB
#
# (fbox 40, Nint 160) and (fbox 40, Nint 161) are NOT here on purpose: they are
# already tasks iN=2,iL=1 and iN=3,iL=1 of fluctuation-solver-set3_2_box40_hi.sh,
# at the same box, the same grid and the same lNoise = Sqrt[3] with saveData=True,
# so runTag would collide exactly.  That script also gives them their lNoise = 1
# and lNoise = 5 columns, which this one does not, so it is the better home for
# them.  NB it requests only 1600G against a predicted 819 GiB at Nint = 161
# (1.8x headroom); consider raising it to 2000G to match this file.
#
# EVEN/ODD PAIRS.  Every box appears with both an even and an odd Nint.  In this
# cell-centred scheme x_i = -fbox + (i - 1/2) h, so ODD Nint puts the defect core
# (S = 0) exactly on a grid point and EVEN Nint straddles it -- the same reason the
# standard sweeps always pair 20/21, 30/31, and so on.  Ordering by Nint splits
# each pair by one task (0<->2, 1<->3, 4<->5, 6<->7), which costs nothing: the
# pairing is a property of the numbers, not of the submission order.
#
# TWO KINDS OF h-MATCHING, both usable.
#   * Within this file: tasks 0 and 6 share h = 1.0000 exactly (and their odd
#     partners 2 and 7 agree to 0.17%), so boxes 60 and 75 differ by box size
#     alone, with no discretisation confound.
#   * Against the Dirichlet partner: h here is 2 fbox/Nint, while set2_2's
#     node-centred h is 2 fbox/(Nint+1), so an ODD Nint here matches EVEN Nint-1
#     there at the same box, exactly -- (60,121)<->(60,120), (75,151)<->(75,150),
#     (120,121)<->(120,120), (40,141)<->(40,140).
#     fluctuation-solver-set2_2_boxscan.sh carries the identical row table in the
#     identical order, so the same task id is the same (fbox, Nint) in both files.
#
# READ THIS ABOUT TASKS 1 AND 3.  At l/h ~ 0.87 the Gaussian is narrower than one
# cell: it is essentially invisible to the grid, and the run reverts to the
# grid-limited no-cutoff (set3_1) answer.  Those rows measure the mesh, not the
# model.  They are worth having as the far end of the crossover, but they are not
# a measurement of lNoise = Sqrt[3].  The module prints lNoiseOverH and issues
# ::subgrid for them.
#
# WHY solverParams box TRACKS fbox PER ROW.  fbox is the Lyapunov grid half-width;
# solverParams box is the FEM domain half-width of the steady state
# (Rectangle[{-box,-box},{box,box}] in iterative_solver_module.m).  The Lyapunov
# module samples rho_ss/Q1_ss/Q2_ss AT ITS OWN GRID POINTS, i.e. out to +-fbox, so
# a box smaller than fbox evaluates the FEM InterpolatingFunction outside its
# ElementMesh, which extrapolates -- and here it would put the sealed wall of the
# fluctuation problem outside the wall of the steady state it is linearized about.
# Because fbox varies per row, box is set FROM the row, never hard-coded.
#
# MEMORY, measured not guessed.  sacct MaxRSS from the box-20 sealed run
# (jobscripts/logs/56583960*) is a clean n^2 law:
#     Nint =  60 ->  13.2 GiB     Nint = 100 -> 122.5 GiB
#     Nint =  80 ->  50.8 GiB     Nint = 130 -> 348.0 GiB
#                                 Nint = 131 -> 358.8 GiB
# The coefficient is 1.354e-7 GiB per n^2, i.e. 145 bytes per n^2 -- about NINE
# dense complex n x n matrices (16 B each), which is what the eigen-solve holds.
# So  peak GiB ~ 1.354e-7 * (3 Nint^2)^2.  It does NOT depend on fbox.
# The largest row here is Nint = 151 at 634 GiB, so --mem=2000G leaves 3.2x
# headroom.  2000G costs nothing in scheduling breadth: 424 node-partition entries
# with >= 72 cpus accept it, exactly as many as 1600G, and 267 of them are in
# "medium" (only past 2200G does the pool drop, to 256).
#
# TIME.  Nint = 131 took 6 h 21 m at 72 cpus and the eigen-solve is ~n^3, so the
# Nint = 150/151 rows extrapolate to ~15 h, comfortably inside the 2-day limit.
#
# DISK.  saveData is True for every row, since every row is the physical lNoise.
# The full covariance is n^2 x 8 bytes: ~15 GB per Nint = 120/121 row, ~28 GB per
# Nint = 140/141 row, ~37 GB per Nint = 150/151 row -- about 191 GB for the eight.
#
# NO OUTPUT COLLISIONS.  runTag is "box<fbox>_N<Nint>_Neu_lN<value>".  Boxes 60, 75
# and 120 appear in no other sweep, and the box-40 rows use Nint = 140/141, which
# the _box40 scripts (20..131 and 150..161) do not cover.

# #SBATCH --exclusive

# LICENSE THROTTLE: the %8 suffix caps this array at 8 concurrent tasks, i.e. all
# eight rows at once.  set2_1's Nint = 131 task died with "The product exited
# because of a license error" (jobscripts/logs/56558539_13) when 18 kernels ran at
# once, so 8 is comfortably clear -- but the cap is PER JOBSCRIPT, so running this
# together with set2_2's boxscan is 16 kernels, and adding any _box40 array on top
# will exceed what killed that task.
# Array indices stay 0-BASED: task id indexes rows directly, so 1-N would silently
# skip the first row.
#SBATCH --array=0-7%8
#SBATCH --job-name="s32bscan"

#SBATCH --mem=2000G
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


# ---- the row table: one "fbox Nint" pair per array task, 0-based -------------
# ORDERED BY Nint (cost ~ n^3, n = 3 Nint^2), so the cheap rows run first.
# Identical table, identical order, in fluctuation-solver-set2_2_boxscan.sh.
# (40 160) and (40 161) are deliberately absent -- see the header.
rows=(
  "60 120"
  "120 120"
  "60 121"
  "120 121"
  "40 140"
  "40 141"
  "75 150"
  "75 151"
)

row=${rows[$SLURM_ARRAY_TASK_ID]}

# Fail loudly rather than passing an empty Nint to wolframscript: an out-of-range
# task id would otherwise expand to "" and abort deep inside the parameter parse.
if [ -z "$row" ]; then
  echo "ERROR: no row defined for array task ${SLURM_ARRAY_TASK_ID}; rows has ${#rows[@]} entries"
  exit 1
fi

read -r fbox Nint <<< "$row"

# box tracks fbox -- see the header.  plotBox = fbox because with sealed walls the
# boundary IS the point of the run, so nothing is cropped out of the figures.
box=$fbox
plotBox=$fbox

# every row is the physical lNoise, so every row keeps its full covariance
lNoise="Sqrt[3.]"
saveData=True

mkdir -p $SCRATCH_PATH
mkdir -p $outputdir

module load mathematica/14

unset DISPLAY
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $setdir

wolframscript -file fluctuations.wls \
  "modelParams={lNoise->${lNoise}}; \
   solverParams={box->${box}}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->${fbox}}; \
   otherParams={outputDir->\"${outputdir}\", saveData->${saveData}, plotBox->${plotBox}}"

END=$(date +%s.%N)
echo "elapsed: $(echo "$END - $START" | bc) s   (fbox = ${fbox}, Nint = ${Nint}, lNoise = ${lNoise}, saveData = ${saveData})"



# The following partitions are available
# debug	        10:00
# short	        2:00:00
# medium	    2-00:00:00
# long	        14-00:00:00
# extra_long	28-00:00:00
