library(protti)
library(tidyverse)
library(data.table)
library(ampir)
args <- commandArgs(trailingOnly = TRUE)

bait <- args[1]

# bait <- "qx017134"
path_psm <- paste("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl/" ,bait ,".csv", sep = "")
path_genes <- paste("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/genes/", bait, ".csv", sep="")
output_path <- paste("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/protein_coverage/", bait, ".csv", sep="")

#	Load data
# Load peptide-protein matches
tbl_psm_csv <- fread(path_psm) %>% as_tibble

# Load identified proteins
tbl_genes <- fread(path_genes) %>% as_tibble

# Load fasta file
path_fasta <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2_coverage.fasta"
df_fasta <- read_faa(file = path_fasta)

tbl_fasta <- df_fasta %>% as_tibble %>% separate(seq_name, into=c("sp", "accession_match", "gene"))

tbl_psm_accession <- tbl_psm_csv %>% 
	mutate(accession_match = str_replace(accession, "_CoV_2_", "")) %>%
	separate(accession_match, into=c("sp", "accession_match", "gene"))

tbl_genes_accession <- tbl_genes %>%
	mutate(accession_match = str_replace(accession, "_CoV_2_", ""))

tbl_psm_fasta <- merge(tbl_psm_accession, tbl_fasta, by="accession_match") %>% as_tibble

tbl_identified_fasta <- tbl_psm_fasta %>% filter(accession_match %in% tbl_genes_accession$accession_match)

tbl_coverage <- calculate_sequence_coverage(
	tbl_identified_fasta,
	protein_sequence = seq_aa,
	peptides = sequence
)

tbl_coverage_match <- tbl_coverage %>% select(accession_match, coverage) %>% unique()
tbl_genes_coverage <- merge(tbl_genes_accession, tbl_coverage_match, by="accession_match") %>% as_tibble
tbl_genes_selection <- tbl_genes_coverage %>% select(bait, accession, gene_name, distinct_peptide_sequences, coverage, number_of_psms)

tbl_genes_selection %>% filter(str_detect(accession, "SARS"))

fwrite(tbl_genes_selection, output_path, sep= "\t")