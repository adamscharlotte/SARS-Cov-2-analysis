# install.packages("protti", dependencies = TRUE)
library(protti)
library(tidyverse)
library(data.table)
# library(seqinr)
library(ampir)
library("xlsx")

#	Load data
# Load peptide-protein matches
# path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl"
# setwd(path)
# files <- dir(pattern = "*.csv")
# psm_csv <- files %>%
# 	map(read_csv) %>%       # read in all the files individually, using
# 							# the function read_csv() from the readr package
# 	reduce(rbind)			# reduce with rbind into one dataframe

path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl_semi"
setwd(path)
files <- dir(pattern = "*.csv")
tbl_psm_csv <- files %>%
	map(read_csv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind)			# reduce with rbind into one dataframe

# # Load identified proteins
# path_proteins <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/proteins"
# files <- dir(path_proteins, pattern = "*.csv")
# files_dir <- paste(path_proteins, "/", files, sep = "")
# tbl_proteins <- files_dir %>%
#     map(read_csv) %>%       # read in all the files individually, using
#                             # the function read_tsv() from the readr package    
#     purrr::reduce(rbind)           # reduce with rbind into one dataframe

# Load identified genes
path_genes <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/genes"
files <- dir(path_genes, pattern = "*.csv")
files_dir <- paste(path_genes, "/", files, sep = "")
tbl_genes <- files_dir %>%
	map(read_csv) %>%			# read in all the files individually, using
								# the function read_tsv() from the readr package    
	purrr::reduce(rbind)		# reduce with rbind into one dataframe

# Load fasta file
path_fasta <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2_coverage.fasta"
# path_fasta <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/test.fasta"
# df_fasta <- read.fasta(path_fasta, seqtype="AA", as.string=TRUE)
df_fasta <- read_faa(file = path_fasta)
tbl_fasta <- df_fasta %>% as_tibble %>% separate(seq_name, into=c("sp", "accession", "gene"))

tbl_psm_accession <- tbl_psm_csv %>% 
	mutate(accession = str_replace(accession, "_CoV_2_", "")) %>%
	separate(accession, into=c("sp", "accession", "gene"))

tbl_genes_accession <- tbl_genes %>%
	mutate(accession = str_replace(accession, "_CoV_2_", ""))

tbl_psm_fasta <- merge(tbl_psm_accession, tbl_fasta, by="accession") %>% as_tibble

tbl_identified_fasta <- tbl_psm_fasta %>% filter(accession %in% tbl_genes_accession$accession)
tbl_identified_fasta %>% select(accession) %>% unique

tbl_coverage <- calculate_sequence_coverage(
	tbl_identified_fasta,
	protein_sequence = seq_aa,
	peptides = sequence
)

tbl_coverage %>% select(accession, coverage) %>% unique() %>% arrange(desc(coverage))
tbl_coverage %>% select(accession, coverage) %>% unique() %>% filter(str_detect(accession, "SARS")) %>% arrange(desc(coverage))
tbl_proteins %>% select(accession) %>% unique()
tbl_genes %>% select(accession) %>% unique()

tbl_identified_fasta %>% filter(str_detect(accession, "SARS"))
tbl_fasta %>% filter(str_detect(accession, "SARS"))
tbl_psm_accession %>% filter(str_detect(accession, "SARS"))
tbl_genes %>% filter(str_detect(accession, "SARS"))

identified <- tbl_coverage %>% select(accession, coverage) %>% unique()
write.xlsx(tbl_coverage, file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/supplementary.xlsx",
      sheetName = "proteins", append = FALSE)

# tbl_identified_fasta <- tbl_psm_fasta %>% filter(accession %in% tbl_proteins$accession)
# tbl_identified_fasta

# tbl_coverage_proteins <- calculate_sequence_coverage(
# 	tbl_identified_fasta,
# 	protein_sequence = seq_aa,
# 	peptides = sequence
# )

# tbl_coverage_proteins %>% select(accession, coverage) %>% unique()

tbl_psm_accession <- tbl_psm_csv %>% 
	separate(accession, into=c("sp", "accession", "gene"))

tbl_psm_accession %>% filter(accession %in% tbl_genes$accession)


tbl_genes %>% filter(accession %in% tbl_psm_accession$accession)



data <- data.frame(
	protein_sequence = c("abcdefghijklmnop", "abcdefghijklmnop", "abcdefghijklmnop"),
	pep_stripped_sequence = c("abc", "jkl", "bcde")
)

calculate_sequence_coverage(
	data,
	protein_sequence = protein_sequence,
	peptides = pep_stripped_sequence
)