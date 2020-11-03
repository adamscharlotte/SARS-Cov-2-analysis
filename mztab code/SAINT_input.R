library(tidyverse)
library(data.table)

#   create bait input file
genes <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/genes"
bait_output <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/bait_file.txt"

files <- dir(genes, pattern = "*.csv")
files_dir <- paste(genes, "/", files, sep = "")
proteins <- files_dir %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package    
    purrr::reduce(rbind)           # reduce with rbind into one dataframe

baitFile <- proteins %>% as_tibble %>% mutate(C_or_T = str_replace(C_or_T, "TRUE", "T")) %>%
    select(SampleID, Condition, C_or_T) %>% unique

fwrite(baitFile, bait_output, sep= "\t", col.names=FALSE)

#   create prey input file
length <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/human_sars_length.txt", sep= "\t")
prey_output <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/prey_file.txt"

accession_length <- length %>% as_tibble %>%
    separate(V1, into=c("pre", "accession", "post"), sep=">>") %>% 
    rename(proteinlength = V2) %>% select(accession,proteinlength)

proteins_length <- merge(proteins, accession_length, by="accession") %>% as_tibble
preyFile <- proteins_length %>% select(accession, proteinlength, gene_name) %>% unique()

fwrite(preyFile, prey_output, sep= "\t", col.names=FALSE)

#   create interaction input file
interaction_output <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/interaction_file.txt"

interactionFile <- proteins %>% as_tibble %>%
    select(SampleID, Condition, accession, number_of_psms)

fwrite(interactionFile, interaction_output, sep= "\t", col.names=FALSE)