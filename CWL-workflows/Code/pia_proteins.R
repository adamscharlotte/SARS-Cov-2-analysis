library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

csv <- args[1]
filename <- args[2]

# path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/ANN-SoLo runs/26.11.2019/pia/"
# filename <- "b10000_ZNF230"

# input_path <- paste(path, "csv/", filename, ".csv", sep = "")
output_path <- paste(filename, ".csv", sep = "")

#	read pia csv file
mztab = fread(csv, sep = '\t')

#	determine mass_tolerance boundaries
mztab_mtol <- mztab %>% as_tibble %>%
	select(accession, description,opt_global_number_of_psms, `opt_global_cv_MS:1001097_distinct_peptide_sequences`) %>%
    mutate(bait = filename) %>%
    rename(number_of_psms = opt_global_number_of_psms) %>%
    rename(distinct_peptide_sequences = `opt_global_cv_MS:1001097_distinct_peptide_sequences`)

#   save file
fwrite(mztab_mtol, output_path, append = FALSE, col.names = TRUE)