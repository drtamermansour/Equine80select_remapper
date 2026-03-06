#!/usr/bin/env bash
# submit_slurm.sh — SLURM wrapper for run_pipeline.sh.
#
# Submits the remapping pipeline as a SLURM job.
# Edit the SBATCH directives below to match your cluster's configuration,
# then pass the same arguments you would give to run_pipeline.sh.
#
# Usage:
#   bash submit_slurm.sh -i <manifest.csv> -r <reference.fa> [run_pipeline.sh options...]
#
# Example:
#   bash submit_slurm.sh \
#       -i backup_original/Equine80select_24_20067593_B1.csv \
#       -r equCab3/equCab3_genome.fa \
#       -a equCab3 \
#       -o results/ \
#       -t 64
# ──────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=remap_pipeline
#SBATCH --output=remap_pipeline_%j.out
#SBATCH --error=remap_pipeline_%j.err
#SBATCH --time=01-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G

conda activate remap
bash "$SCRIPT_DIR/run_pipeline.sh" "$@"
EOF

echo "SLURM job submitted. Monitor with: squeue -u \$USER"
