#!/bin/bash
#SBATCH --job-name=merge_reads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=merge_reads_%j.out
#SBATCH --error=merge_reads_%j.err

set -euo pipefail

# Input and output directories
working_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
source_dir="trimmed"  # Change to "fixed" if you ran fix_trimmed_fastq.sh
input_dir="$working_dir/output/$source_dir"
output_dir="$working_dir/output/merged"

# PEAR must be installed and available in PATH before submitting this job.
if ! command -v pear >/dev/null 2>&1; then
    echo "PEAR was not found. Load the PEAR module or install it first." >&2
    exit 1
fi

if [[ ! -d "$input_dir" ]]; then
    echo "Input directory not found: $input_dir" >&2
    exit 1
fi

mkdir -p "$output_dir"
shopt -s nullglob
r1_files=("$input_dir"/*_R1_001_trim.fastq)
if (( ${#r1_files[@]} == 0 )); then
    echo "No R1 trimmed FASTQ files found in $input_dir" >&2
    exit 1
fi

# Make sure every R1 and R2 file has a matching mate before merging.
for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1_001_trim.fastq)
    r2="$input_dir/${sample}_R2_001_trim.fastq"
    if [[ ! -f "$r2" ]]; then
        echo "Missing R2 file for sample: $sample" >&2
        exit 1
    fi
done
for r2 in "$input_dir"/*_R2_001_trim.fastq; do
    sample=$(basename "$r2" _R2_001_trim.fastq)
    if [[ ! -f "$input_dir/${sample}_R1_001_trim.fastq" ]]; then
        echo "Missing R1 file for sample: $sample" >&2
        exit 1
    fi
done

# PEAR writes .assembled.fastq and unassembled read files using this prefix.
for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1_001_trim.fastq)
    r2="$input_dir/${sample}_R2_001_trim.fastq"
    echo "Merging $sample"
    pear -f "$r1" -r "$r2" -o "$output_dir/$sample" -j "${SLURM_CPUS_PER_TASK:-2}"
done

echo "Finished merging ${#r1_files[@]} sample(s). Output: $output_dir"
