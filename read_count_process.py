# Script that processes the read count data from htseq-count and splits the read counts by source as well
# as computing statistics about them, such as how large a fraction of the mapped reads map to the parasite
# Usage: python read_count_process.py /Path/to/metadata/file.tsv /Path/to/read/count/folder /Path/to/output/folder
# python scripts/read_count_process.py results/htseq/reverse/pilot/ results/processed_reads
# metadata: data/metadata/pilot_metadata.tsv

#!/usr/bin/python
import pandas as pd
from os import listdir, mkdir
import sys
import re
import numpy as np

# Read in the input and output paths from stdin and define folders for the output
#metadata_path = sys.argv[1]
read_count_path = sys.argv[1]
output_path = sys.argv[2]
mkdir(output_path+"/chicken")
mkdir(output_path+"/eimeria")
mkdir(output_path+"/stats")

# Read in the data
read_count_filenames = listdir(read_count_path)
read_count_filenames.sort()
#metadata_table = pd.read_csv(metadata_path, sep="\t")
#metadata_table = metadata_table.sort_values("File_name")
result_list = []
old_name = ""

i = 0
j = 0
while i < len(read_count_filenames):
    # The loop goes through each read count file, separates the chicken and E. tenella reads and sums the 
    # read counts from paired and unpaired reads for a given sample.  It also extracts statistics from the
    # bottom of the HTSeq output file
    if re.search("L00", read_count_filenames[i]):
        # Make sure only read count files are read in
        l_pos = re.search("L00", read_count_filenames[i]).start()
        curr_name = read_count_filenames[i][0:l_pos-1]

        if curr_name != old_name:
            sum_table = pd.read_csv(read_count_path+"/"+read_count_filenames[i], sep="\t", header=None, names=["gene_name","gene_count"], index_col=0)
        else:
            temp_table = pd.read_csv(read_count_path+"/"+read_count_filenames[i], sep="\t", header=None, names=["gene_name","gene_count"], index_col=0)
            sum_table += temp_table

            if i+1 == len(read_count_filenames) or not re.search(curr_name, read_count_filenames[i+1]):
                stat_table = sum_table.iloc[-5:]
                eimeria_table = sum_table.iloc[24389:-6]
                chicken_table = sum_table.iloc[np.r_[0:24389, -6]]
                

                total_read_num_sum = sum(sum_table.gene_count.iloc[:-5])
                eimeria_read_num = sum(eimeria_table.gene_count)
                chicken_read_num = sum(chicken_table.gene_count)

                eimeria_read_perc = eimeria_read_num/total_read_num_sum * 100
                no_feature_num = stat_table.gene_count.iloc[0]
                ambiguous_num = stat_table.gene_count.iloc[1]

                result_list.append([curr_name, total_read_num_sum, chicken_read_num, eimeria_read_num, eimeria_read_perc, no_feature_num, ambiguous_num])
                
                stat_table.to_csv(output_path+"/stats/"+curr_name+"_stat_table.csv", sep = ",")
                eimeria_table.to_csv(output_path+"/eimeria/"+curr_name+"_eimeria_table.csv", sep = ",")
                chicken_table.to_csv(output_path+"/chicken/"+curr_name+"_chicken_table.csv", sep = ",")
                j += 1
        old_name = curr_name
    i += 1

header_list = ["File_name","Total_number_of_mapped_reads", "Number_of_reads_mapped_to_chicken", "Number_of_reads_mapped_to_Eimeria", \
    "Percentage_of_Eimeria_reads", "Number_of_read_not_mapped_to_feature", "Number_of_reads_mapped_to_multiple_features"]
data_table = pd.DataFrame(result_list, columns=header_list)
data_table.to_csv(output_path+"/metadata_table_short.csv", sep = ",", index = False)
#out_table = pd.merge(metadata_table, data_table, on = "File_name")
#out_table.to_csv(output_path+"/metadata_table.csv", sep = ",", index = False)