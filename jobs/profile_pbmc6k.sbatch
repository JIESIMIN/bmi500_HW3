#!/bin/bash
#SBATCH --job-name=prof-pbmc6k
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=logs/profile_pbmc6k_%j.out
#SBATCH --error=logs/profile_pbmc6k_%j.err

# Profiling run for pbmc6k  (5,419 cells)

set -euo pipefail

DS=pbmc6k
HW=/home/rzhou73/BMI500/HW3
cd "$HW"
mkdir -p logs profiles

source /home/rzhou73/bmi500/bin/activate

# Pin every threading layer to 1 so timings are comparable across datasets
# rather than dependent on how busy the node is.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

COARSE="profiles/${DS}_coarse_time.txt"
INSTR="profiles/${DS}_instrumented_time.txt"
FINE="profiles/${DS}_rank_genes_groups_profile.txt"

echo "================================================================"
echo " dataset : $DS"
echo " node    : $(hostname)"
echo " job     : ${SLURM_JOB_ID}"
echo " cpus    : ${SLURM_CPUS_PER_TASK}    mem: 24G"
echo " started : $(date)"
echo "================================================================"


if [[ ! -d "data/${DS}/filtered_gene_bc_matrices" ]]; then
    echo "[setup] extracting data/${DS}.tgz"
    tar xzf "data/${DS}.tgz" -C data/
fi

# THE MEASURED RUN 
echo
echo "########## RUN: scanpy_pbmc.py under /usr/bin/time -v ##########"
/usr/bin/time -v -o "$COARSE" python scanpy_pbmc.py --data-set "$DS"

# ALL THREE EXPERIMENTS
echo
echo "########## 1. COARSE-GRAINED ($COARSE) ##########"
cat "$COARSE"

echo
echo "########## 2. INSTRUMENTED ($INSTR) ##########"
cat "$INSTR"

echo
echo "########## 3. FINE-GRAINED ($FINE) ##########"
head -20 "$FINE"

echo
echo "finished : $(date)"
ls -lh profiles/${DS}_* 2>/dev/null || true
