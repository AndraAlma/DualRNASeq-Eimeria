# Concatenates the count files for each organism into a single file and fixes the gene names for both.  Also fixes the
# gene names of an Eimeria annotation file.
# Usage: python scripts/read_count_concat.py results/processed_reads/all_reads/chicken results/processed_reads/all_reads/eimeria data/reference/chicken_ref/chicken_genes.tsv data/reference/eimeria_ref/eimeria_genes.tsv data/reference/eimeria_ref/eimeria_gene_products.tsv results/processed_reads/concatenated

import pandas as pd
from os import listdir, mkdir
import sys
import re

chicken_path = sys.argv[1]
eimeria_path = sys.argv[2]
chicken_gene_path = sys.argv[3]
eimeria_gene_path = sys.argv[4]
eimeria_product_path = sys.argv[5]
output_path = sys.argv[6]

chicken_file_names = listdir(chicken_path)
eimeria_file_names = listdir(eimeria_path)
first_table = True

# Concatenate the chicken files, remove the gene- prefix in the gene symbols and replace the gene identifiers with Entrez Gene IDs
for file_name in chicken_file_names:
    sample_name = file_name[0:re.search("_chicken", file_name).start()]
    current_table = pd.read_csv(chicken_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

chicken_gene_table = pd.read_csv(chicken_gene_path, sep = "\t", header = None, names = ["gene_name","entrez_gene_id"])
old_table = pd.merge(old_table, chicken_gene_table, on = "gene_name")
temp_table = pd.concat([old_table.iloc[:,-1], old_table.iloc[:,1:-1]], axis = 1)
old_table = pd.concat([old_table.iloc[:,0], temp_table], axis = 1)

i = 0
while i < len(old_table):
    old_table.iloc[i,0] = old_table.iloc[i,0][0:]
    i += 1

old_table.to_csv(output_path+"/chicken_counts.csv", sep = ",", index = False)

# Concatenate the Eimeria files and replace the gene identifiers with Entrez Gene IDs and locus tags
first_table = True

for file_name in eimeria_file_names:
    sample_name = file_name[0:re.search("_eimeria", file_name).start()]
    current_table = pd.read_csv(eimeria_path+"/"+file_name, header = 0, names = ["gene_name",sample_name])
    if first_table:
        old_table = current_table
        first_table = False
    else:
        old_table = pd.merge(old_table, current_table, on = "gene_name")

eimeria_gene_table = pd.read_csv(eimeria_gene_path, sep = "\t", header = None, names = ["gene_name","entrez_gene_id","locus_tag"])
old_table = pd.merge(old_table, eimeria_gene_table, on = "gene_name")
old_table = pd.concat([old_table.iloc[:,-1], old_table.iloc[:,-2], old_table.iloc[:,1:-2]], axis = 1)
old_table.to_csv(output_path+"/eimeria_counts.csv", sep = ",", index = False)

# Replace Eimeria gene identifiers with Entrez Gene IDs and locus tags for an annotation file
eimeria_product_table = pd.read_csv(eimeria_product_path, sep = "\t", header = None, names = ["gene_name","product"])
eimeria_product_table = pd.merge(eimeria_product_table, eimeria_gene_table, on = "gene_name")
eimeria_product_table = pd.concat([eimeria_product_table.iloc[:,-1], eimeria_product_table.iloc[:,-2], eimeria_product_table.iloc[:,1]], axis = 1)
eimeria_product_table.to_csv(output_path+"/eimeria_gene_products.csv", sep = "\t", index = False)