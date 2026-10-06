#!/bin/bash
#SBATCH --job-name=s22co2D_n50
#SBATCH --partition=short,medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=0-02:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out

# ===========================================================================
#  set2_2_comove_optimized_2, ssDir variant -- Nint = 50   (n = 3 Nint^2 = 7500)
#
#  2 tasks = 2 B at fbox = steady box = 20.  Parameters: params.md.
#  Steady state: free Q wall + DIRICHLET density wall (fluctuations_ssDir.wls).
#  Lyapunov: Dirichlet.  Own outputdir: the runTag does not encode the steady wall.
#
#  Pipeline check.  Measured ~3 min per task, ~3 GiB peak.
#  (sacct of set2_2_comove_optimized; the Lyapunov cost depends on Nint only.)
#
#  --cpus-per-task=24 is MEASURED: the Schur reduction peaks at 24 threads and
#  gets slower above it.  See lyapunov_solver/resource_measurements.md.
#
#  The array throttle (%5) is PER JOBSCRIPT.  Too many concurrent
#  Mathematica kernels have died with a license error before.  Stage them.
# ===========================================================================

set -u
cd /home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2_comove_optimized_2
unset DISPLAY
module load mathematica/14 2>/dev/null

# Node-local NVMe, NOT $TMPDIR -- that points at shared GPFS here.
# The backend stages A, Q and C through this: 3 * 8n^2 = 1 GiB.
scratchdir="/scratch/lyap_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$scratchdir" || { echo "FATAL: cannot create $scratchdir"; exit 1; }

# Always remove the scratch directory when the task ends, however it ends.
trap 'rm -rf "$scratchdir"' EXIT

outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2_comove_optimized_2_ssDir/"
mkdir -p "$outputdir"

Nint=50
# Both run at once (%2).
BList=(100000 90000)
box=20
ssBox=20
plotBox=$(( box - 3 ))

i=${SLURM_ARRAY_TASK_ID}
B=${BList[$i]:-}

if [ -z "$B" ]; then
  echo "FATAL: array id $i did not resolve (B='$B')"
  exit 1
fi

echo "=== task $i : Nint=$Nint  B=$B  fbox=$box  steady box=$ssBox  plotBox=$plotBox"
echo "=== host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}  scratch=$scratchdir"
df -h /scratch | tail -1

# Stagger kernel start-up.  Two tasks launched in the same second on one node
# (2006623, gaspra01) saw one wolframscript segfault before the driver printed
# anything (exit 139).
sleep $(( i * 30 ))

wolframscript -file fluctuations_ssDir.wls \
  "modelParams={lNoise->0.5, B->${B}.}; \
   solverParams={box->${ssBox}.}; \
   lyapunovSolverParams={Nint->${Nint}, fbox->${box}, \
                         scratchDir->\"${scratchdir}\"}; \
   otherParams={outputDir->\"${outputdir}\", saveData->True, \
                plotBox->${plotBox}}"
rc=$?
echo "=== wolframscript exit $rc"

# wolframscript returns 0 even when the script called Abort[], so its exit code
# cannot be trusted.  Check the artefacts.
shopt -s nullglob
produced=( "${outputdir}"covariance_box${box}_N${Nint}_Dir_*_B${B}_ss${ssBox}.bin "${outputdir}"covariance_box${box}_N${Nint}_Dir_*_B${B}_ss${ssBox}_unstable.bin )
if [ ${#produced[@]} -eq 0 ]; then
  echo "FATAL: no covariance_box${box}_N${Nint}_Dir_*_B${B}_ss${ssBox}[_unstable].bin was written."
  echo "       The run aborted despite exit $rc -- see above for the cause."
  exit 1
fi
expected=$(( 8 * (3 * Nint * Nint) * (3 * Nint * Nint) ))
actual=$(stat -c %s "${produced[0]}")
echo "=== wrote $(basename "${produced[0]}")  ${actual} bytes (expected ${expected})"
case "${produced[0]}" in
  *_unstable.bin) echo "=== NOTE: tagged _unstable -- max Re(lambda) >= 0, so this"
                  echo "===       result is NOT a covariance (Sigma_Q may be negative)."
                  echo "===       It is a convergence-curve point; exclude it from averages." ;;
esac
if [ "$actual" -ne "$expected" ]; then
  echo "FATAL: covariance file is ${actual} bytes, expected 8*n^2 = ${expected}."
  exit 1
fi
exit $rc
