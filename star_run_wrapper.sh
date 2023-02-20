#!/bin/bash -l

#SBATCH -A snic2022-22-1221
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J star_run
#SBATCH --mail-type=ALL
#SBATCH --mail-user alma.hansen.3759@student.uu.se

module load bioinfo-tools
module load star/2.7.2b

# Your commands
bash scripts/star_run.sh results/trimmomatic/VJ3440/221111