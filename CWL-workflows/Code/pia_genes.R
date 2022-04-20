library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
file <- args[2]  

input_path <- paste(path, "proteins/", file, ".csv", sep = "")
output_path <- paste(path, "genes/", file, ".csv", sep = "")
fasta_headers <- paste(path, "fasta/fasta_headers.csv", sep = "")
annotation_path <- paste(path, "annotation.txt", sep = "")

#	read protein file
tbl_proteins = fread(input_path, sep = ',') %>% as_tibble

#   add gene_names
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = ("GN=")) %>% drop_na(gene) %>%
    separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
    select(V2, gene_name)
pia_gene <- merge(pia_proteins, gene_names, by.x = "accession", by.y = "V2") %>% as_tibble

#   add annotation and filter out removed files
annotation <- fread(annotation_path, sep= "\t")
pia_gene_annotation <- merge(pia_gene, annotation, by.x = "bait", by.y = "fileID") %>% as_tibble %>%
    filter(Comments=="")

fwrite(pia_gene_annotation, output_path, append = FALSE, col.names = TRUE)

##########################################################################################################################################

# #   only keep the one of multiple gene names per run with highest psms 
# pia_genes <- pia_gene %>% group_by(gene_name) %>% 
#     select(accession, description, distinct_peptide_sequences, bait, gene_name, number_of_psms) %>% 
#     top_n(1)

# #   save file
# fwrite(pia_genes, output_path, append = FALSE, col.names = TRUE)