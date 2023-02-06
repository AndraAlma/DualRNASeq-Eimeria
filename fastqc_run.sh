# Script for running FastQC on a set of Trimmomatic trimmed read files and outputting both FastQC and MultiQC reports for them
# Usage: bash fastqc_run.sh /Path/to/folder/containing/trimmed/read/files1 /Path/to/folder/containing/trimmed/read/files2 ...

for FILE_PATH in "$@"
do
    GROUP_NAME=$(echo ${FILE_PATH} | cut -d '/' -f 3)
    mkdir -p results/fastqc/${GROUP_NAME} 
    fastqc -o results/fastqc/${GROUP_NAME} ${FILE_PATH}/*.fastq.gz
    multiqc -o results/fastqc/${GROUP_NAME} results/fastqc/${GROUP_NAME}
done
