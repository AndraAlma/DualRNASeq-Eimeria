# Script for producing an index for STAR
# Usage:  bash star_index.sh /Path/to/reference/genome.fna /Path/to/reference/genome.gff

mkdir -p results/star2/index

STAR --runMode genomeGenerate --genomeDir results/star2/index --genomeFastaFiles ${1} --sjdbGTFtagExonParentTranscript Parent \
	--sjdbGTFtagExonParentGene gene --runThreadN 8 --sjdbGTFfile ${2} --sjdbOverhang 150 --genomeChrBinNbits 18