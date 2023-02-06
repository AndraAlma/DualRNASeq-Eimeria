# Short bash script for using trimmomatic with user input for the location of the read files
# Usage: bash trimmomatic_run.sh /Path/to/folder/containing/read/files1 /Path/to/folder/containing/read/files2 ...

for FILE_PATH in "$@"
do
    GROUP_NAME=$(echo ${FILE_PATH} | cut -d '/' -f 2)
    mkdir -p results/trimmomatic/${GROUP_NAME}
    for f in $(ls ${FILE_PATH} | cut -d '/' -f 3 | sed 's/_.._001.fastq.gz//' | sort -u)
    do
        trimmomatic PE -threads 8 \
        ${FILE_PATH}/${f}_R1_001.fastq.gz \
        ${FILE_PATH}/${f}_R2_001.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R1_paired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R1_unpaired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R2_paired.fastq.gz \
        results/trimmomatic/${GROUP_NAME}/${f}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        SLIDINGWINDOW:4:20 \
        MINLEN:50
    done
done