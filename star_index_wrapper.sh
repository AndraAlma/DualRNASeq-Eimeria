#!/bin/bash -l

#SBATCH -A snic2022-22-1221
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J star_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user alma.hansen.3759@student.uu.se

module load bioinfo-tools
module load star/2.7.2b

# Your commands
bash scripts/star_index.sh data/reference/merged_reference.fna data/reference/fixed_copy.gtf