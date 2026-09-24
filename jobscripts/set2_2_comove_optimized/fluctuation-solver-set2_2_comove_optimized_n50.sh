#!/bin/bash
#SBATCH --job-name=s22co_n50
#SBATCH --partition=short,medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=0-02:00:00
#SBATCH --array=0-19%5
#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/jobscripts/logs/%A_%a.out

# ===========================================================================
#  set2_2_comove_optimized -- Nint = 50   (n = 3 Nint^2 = 7500)
#
#  20 tasks = 2 B  x  5 box  x  2 steady-box variants.
#
#  Estimated ~3 min per task, ~3 GiB peak.  See
#  lyapunov_solver/resource_measurements.md for where those come from.
#
#  --cpus-per-task=24 is MEASURED, not guessed.  The Schur reduction that
#  dominates the solve peaks at 24 threads and gets SLOWER above it; numpy's
#  OpenBLAS is capped at MAX_THREADS=64 so anything beyond that is silently
#  ignored.  At 96 threads, 63 cores spin for 2349 CPU-seconds to deliver the
#  same 37 s of wall that 8 threads deliver with 278.  Asking for more here
#  would only lengthen the queue wait.
#
#  The array throttle (%5) is PER JOBSCRIPT.  Submitting several of these
#  at once multiplies the number of concurrent Mathematica kernels, and a
#  Nint=131 task has previously died with "The product exited because of a
#  license error" when 18 ran together.  Stage them.
# ===========================================================================

set -u
cd /home/ashmat/Projects/myxo-project-dean/jobscripts/set2_2_comove_optimized
unset DISPLAY
module load mathematica/14 2>/dev/null

# Node-local NVMe (28 TB), NOT $TMPDIR -- that points at shared GPFS here.
# The backend stages A, Q and C through this: 3 * 8n^2 = 1 GiB.
scratchdir="/scratch/lyap_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$scratchdir" || { echo "FATAL: cannot create $scratchdir"; exit 1; }

# Always remove the scratch directory when the task ends, however it ends.
trap 'rm -rf "$scratchdir"' EXIT

outputdir="/data/biophys/ashmat/data/myxo-project/fluctuations-set2_2_comove_optimized/"
mkdir -p "$outputdir"

Nint=50
BList=(10000 200000)
BoxList=(8 10 15 20 25)
SSExtraList=(0 5)

# One array index -> (B, box, steady-box offset).  20 = 2 x 5 x 2, enumerated
# so that consecutive ids differ in the steady box, then the box, then B.
i=${SLURM_ARRAY_TASK_ID}
nBox=${#BoxList[@]}      # 5
nSS=${#SSExtraList[@]}   # 2
B=${BList[$(( i / (nBox * nSS) ))]}
box=${BoxList[$(( (i % (nBox * nSS)) / nSS ))]}
ssExtra=${SSExtraList[$(( i % nSS ))]}
ssBox=$(( box + ssExtra ))
plotBox=$(( box - 3 ))

if [ -z "$B" ] || [ -z "$box" ] || [ -z "$ssExtra" ]; then
  echo "FATAL: array id $i did not resolve (B='$B' box='$box' ssExtra='$ssExtra')"
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
# cannot be trusted -- a parameter error would otherwise be recorded as success
# and the missing run only noticed at analysis time.  Check the artefacts.
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
