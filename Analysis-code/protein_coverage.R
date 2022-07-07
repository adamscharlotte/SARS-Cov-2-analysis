library(protti)
library(tidyverse)
library(data.table)
library(ampir)

#	Load data
# Load peptide-protein matches
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl"
setwd(path)
files <- dir(pattern = "*.csv")
tbl_psm_csv <- files %>%
	map(read_csv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind)			# reduce with rbind into one dataframe

# Load identified proteins
path_genes <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/genes"
files <- dir(path_genes, pattern = "*.csv")
files_dir <- paste(path_genes, "/", files, sep = "")
tbl_genes <- files_dir %>%
	map(read_csv) %>%			# read in all the files individually, using
								# the function read_tsv() from the readr package    
	purrr::reduce(rbind)		# reduce with rbind into one dataframe

# Load fasta file
path_fasta <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2_coverage.fasta"
df_fasta <- read_faa(file = path_fasta)

tbl_fasta <- df_fasta %>% as_tibble %>% separate(seq_name, into=c("sp", "accession", "gene"))

tbl_psm_accession <- psm_csv %>% 
	mutate(accession = str_replace(accession, "_CoV_2_", "")) %>%
	separate(accession, into=c("sp", "accession", "gene"))

tbl_genes_accession <- tbl_genes %>%
	mutate(accession = str_replace(accession, "_CoV_2_", ""))

tbl_psm_fasta <- merge(tbl_psm_accession, tbl_fasta, by="accession") %>% as_tibble

tbl_identified_fasta <- tbl_psm_fasta %>% filter(accession %in% tbl_genes_accession$accession)

tbl_coverage <- calculate_sequence_coverage(
	tbl_identified_fasta,
	protein_sequence = seq_aa,
	peptides = sequence
)
