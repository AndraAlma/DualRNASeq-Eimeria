#!/bin/bash -l

#SBATCH -A snic2022-22-1221
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J htseq_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user alma.hansen.3759@student.uu.se

module load bioinfo-tools
module load htseq/0.9.1

# Your commands
bash scripts/htseq_run.sh results/star/mapped_reads/230119