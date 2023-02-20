# Script for running the STAR aligner for a set of read files and reference genomes
# Usage: bash star_run.sh /Path/to/read/files/folder

GROUP_NAME=$(echo ${1} | cut -d '/' -f 4)
mkdir -p results/star/mapped_reads/${GROUP_NAME}

for f in $(ls ${1}/*_paired.fastq.gz | sed 's/_R._paired.fastq.gz//' | sort -u)
do
    OUTPUT_NAME=$(echo ${f} | rev | cut -d '/' -f 1 | rev)
    STAR --genomeDir /proj/snic2022-23-635/private/results/star/index --readFilesIn ${f}_R1_paired.fastq.gz ${f}_R2_paired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    STAR --genomeDir /proj/snic2022-23-635/private/results/star/index --readFilesIn ${f}_R1_unpaired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_R1_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
    STAR --genomeDir /proj/snic2022-23-635/private/results/star/index --readFilesIn ${f}_R2_unpaired.fastq.gz --readFilesCommand gunzip -c \
		--runThreadN 8 --outFileNamePrefix results/star/mapped_reads/${GROUP_NAME}/${OUTPUT_NAME}_R2_unpaired_ --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate
done
