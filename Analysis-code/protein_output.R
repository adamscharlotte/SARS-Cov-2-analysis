library("xlsx")
library(tidyverse)
library(data.table)

#	Load data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/protein_coverage"
setwd(path)
files <- dir(pattern = "*.csv")
tbl_protein_coverage <- files %>%
	map(read_tsv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind)			# reduce with rbind into one dataframe

tbl_protein_coverage

#	Write excell
write.xlsx(tbl_protein_coverage, file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/excell/Protein identifications.xlsx",
      append = FALSE)