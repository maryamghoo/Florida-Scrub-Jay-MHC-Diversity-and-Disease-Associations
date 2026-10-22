#!/bin/bash
#SBATCH --job-name=fix_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=fix_fastq_%j.out
#SBATCH --error=fix_fastq_%j.err

set -euo pipefail

# Read trimmed files and write corrected files to a separate folder.
working_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"
input_dir="$working_dir/output/trimmed"
output_dir="$working_dir/output/fixed"

if [[ ! -d "$input_dir" ]]; then
    echo "Trimmed FASTQ folder not found: $input_dir" >&2
    exit 1
fi

mkdir -p "$output_dir"

shopt -s nullglob
files=("$input_dir"/*_R1_001_trim.fastq "$input_dir"/*_R2_001_trim.fastq)

if (( ${#files[@]} == 0 )); then
    echo "No trimmed FASTQ files found in $input_dir" >&2
    exit 1
fi

for file in "${files[@]}"; do
    # Keep the original filename so PEAR can use it without renaming.
    filename=$(basename "$file")
    sed 's/\r$//' "$file" > "$output_dir/$filename"
done

echo "Finished processing ${#files[@]} FASTQ files. Output: $output_dir"
