#!/bin/bash
#SBATCH --job-name=s22co2_n230
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1700G
#SBATCH --time=8-00:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out

# ===========================================================================
#  set2_2_comove_optimized_2 -- Nint = 230   (n = 3 Nint^2 = 158700)
#
#  2 tasks = 2 B at fbox = steady box = 20.  Parameters: params.md.
#  Steady state: free Q wall + RhoNeumann density wall.  Lyapunov: Dirichlet.
#
#  Measured 3.5-4.6 days per task, 1506 GiB peak (94% of the old 1600G).
#  ~188 GiB per covariance file.
#  (sacct of set2_2_comove_optimized n230; the Lyapunov cost depends on Nint only.)
#
#  --cpus-per-task=24 is MEASURED: the Schur reduction peaks at 24 threads and
#  gets slower above it.  See lyapunov_solver/resource_measurements.md.
#
#  The array throttle (%2) is PER JOBSCRIPT.  Too many concurrent
#  Mathematica kernels have died with a license error before.  Stage them.
# ===========================================================================

set -u
cd /home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2_comove_optimized_2
unset DISPLAY
module load mathematica/14 2>/dev/null

# Node-local NVMe, NOT $TMPDIR -- that points at shared GPFS here.
# The backend stages A, Q and C through this: 3 * 8n^2 = 563 GiB.
scratchdir="/scratch/lyap_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$scratchdir" || { echo "FATAL: cannot create $scratchdir"; exit 1; }

# Always remove the scratch directory when the task ends, however it ends.
trap 'rm -rf "$scratchdir"' EXIT

outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2_comove_optimized_2/"
mkdir -p "$outputdir"

Nint=230
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

wolframscript -file fluctuations.wls \
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
