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

**read_count_process** splits up read counts into two datasets (parasite and chicken) and summarizes the paired and unpaired counts for each sample.

**read_count_process_old** functions the same but is used for the old dataset with different sample names.

**read_count_concat** Concatenates the processed reads into one read count matrix including all samples.

**DEseq** Carries out read count normalization, differential expression analysis, GO and KEGG pathway enrichment, plots expression of specific identified gene groups and a WGCNA gene network analysis on the immune chicken dataset. Also saves some data in variables to plot old and new data side by side.

**DE_old_new** Carried out the same analysis, but on the naïve chicken data, and plots, for example, volcano plots side by side with the immune chicken data.

