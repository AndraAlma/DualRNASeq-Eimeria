# Script for merging two reference genomes for dual RNA-seq analysis.  Each reference should have a fna and gff file. also creates a gtf file from gff.
# Usage: bash reference_merger.sh /Path/to/reference1.fna /Path/to/reference1.gff /Path/to/reference2.fna /Path/to/reference2.gff /Path/to/merged.fna /Path/to/merged.gff /Path/to/merged.gtf

module load bioinfo-tools
module load gffread

cat ${1} ${3} > ${5}
touch ${6}
head -n -1 ${2} >> ${6}
tail -n +18 ${4} >> ${6}
gffread -T ${6} -o ${7}