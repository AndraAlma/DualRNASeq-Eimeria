# DualRNASeq-Eimeria
The folder **scripts** contains all scripts used for pre-processing, mapping, counting and differential expression analysis.

## Preprocessing
These scripts are used to perform quality control on the data, trimming the reads and merging the reference genome.

**fastqc_run** and **trimmomatic_run** defines the run settings.

**fastqc_wrapper** and **trimmomatic_wrapper** calls the respective run scripts, defines input files and softwares.

**reference_merger** concatenates two reference genomes for the reads to be mapped to.

## Read mapping & counting

These scripts are used to index the merged genome, map all reads to the genome and counting the reads.

**star_index** and **star_index_wrapper** indexes the merged reference genome.

**star_run** and **star_run_wrapper** maps the trimmed reads to the genome.

**htseq_run** and **htseq_wrapper** counts all mapped reads.

## Differential Expression



