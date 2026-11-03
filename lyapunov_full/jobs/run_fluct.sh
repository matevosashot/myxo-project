#!/bin/bash
#SBATCH --job-name=lfull
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/home/ashmat/Projects/myxo-project-dean/lyapunov_full/logs/%x_%j.out

# ===========================================================================
#  Full-model fluctuations, one run of run_fluct.wls.
#
#  sbatch [--job-name=..] [--partition=long --time=4-00:00:00] [--mem=..] \
#         jobs/run_fluct.sh <Nint> <milestone> [runTag] [saveC]
#
#  Memory: the backend peaks at ~5 x 8 n^2 (n = 3 Nint^2), plus the sparse
#  assembly: N100 ~ 40G, N150 ~ 190G, N180 ~ 400G.  Scratch: node-local
#  /scratch/$SLURM_JOB_ID (A, the five N_k as CSR, and C.bin if saveC).
#  24 cpus: the Schur step does not scale further with this OpenBLAS.
# ===========================================================================

set -u
cd /home/ashmat/Projects/myxo-project-dean/lyapunov_full
unset DISPLAY
WS=/usr/local/math/math1410/Executables/wolframscript
SCR=/scratch/${SLURM_JOB_ID:-$$}
mkdir -p "$SCR" && export TMPDIR="$SCR"
trap 'rm -rf "$SCR"' EXIT
touch "$SCR/start.stamp"

echo "=== host=$(hostname)  cpus=${SLURM_CPUS_PER_TASK:-?}  job=${SLURM_JOB_ID:-?}  args=$*"
date
# stagger start-up: concurrent kernels starting in the same second have segfaulted
sleep $(( ${SLURM_ARRAY_TASK_ID:-0} * 30 + RANDOM % 20 ))

$WS -file run_fluct.wls "$@" | tee "$SCR/run.log"
rc=${PIPESTATUS[0]}
echo "=== wolframscript exit $rc"
date

# wolframscript returns 0 even after Abort[]: check the artefacts instead.  The log
# can be truncated after the backend (seen at Nint = 180), so find the output by
# its time stamp, not by the "wrote" line.
res=$(find ../data/fluctuations-full -name result.json -newer "$SCR/start.stamp" | head -1)
out=$(dirname "${res:-/nonexistent/x}")
if [ -z "$res" ] || [ ! -s "$out/out.npz" ] || ! grep -q '"ok": true' "$res"; then
  echo "FATAL: no output written."
  exit 1
fi
ls -l "$out"
exit $rc
