library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# csv <- args[1]
# bait <- args[2]  

# # path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/proteasoom subnetwerk/pia/"
# # bait <- "b3712_FOXR2"

# # input_path <- paste(path, "proteins/", bait, "_proteins.csv", sep = "")
# output_path <- paste(bait, ".csv", sep = "")

# #	read protein file
# pia_proteins = fread(csv, sep = ',') %>% as_tibble

# #   add gene_names
# fasta_headers <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/proteasoom subnetwerk/fasta_headers.csv"
# unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
# gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
#     separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
#     select(V2, gene_name)
# pia_gene <- merge(pia_proteins, gene_names, by.x = "accession", by.y = "V2") %>% as_tibble

# fwrite(pia_gene, output_path, append = FALSE, col.names = TRUE)

#   only keep the one of multiple gene names per run with highest psms 
# pia_genes <- pia_gene %>% group_by(gene_name) %>% 
#     select(accession, description, distinct_peptide_sequences, bait, gene_name, number_of_psms) %>% 
#     top_n(1)

#   save file
# fwrite(pia_genes, output_path, append = FALSE, col.names = TRUE)

######################################################################################################################################################
path <- args[1]
file <- args[2]  

# path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/"
# bait <- "qx017077"

input_path <- paste(path, "proteins/", file, ".csv", sep = "")
output_path <- paste(path,"genes/", file, ".csv", sep = "")

# input_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/proteins/qx017077.csv"
# input_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/genes/qx017077.csv"

#	read protein file
pia_proteins = fread(input_path, sep = ',') %>% as_tibble

#   add gene_names
fasta_headers <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble #%>% filter(V2=="SARS_CoV_2_nsp10")
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = ("GN=")) %>% drop_na(gene) %>%
    separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
    select(V2, gene_name)
pia_gene <- merge(pia_proteins, gene_names, by.x = "accession", by.y = "V2") %>% as_tibble %>%
    select(accession, gene_name, bait, number_of_psms, distinct_peptide_sequences)
# pia_gene %>% filter(accession == "SARS_CoV_2_nsp10")
fwrite(pia_gene, output_path, append = FALSE, col.names = TRUE)

#   only keep the one of multiple gene names per run with highest psms 
# pia_genes <- pia_gene %>% group_by(gene_name) %>% 
#     select(accession, description, distinct_peptide_sequences, bait, gene_name, number_of_psms) %>% 
#     top_n(1)

#   save file
# fwrite(pia_genes, output_path, append = FALSE, col.names = TRUE)