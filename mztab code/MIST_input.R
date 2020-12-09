library(tidyverse)
library(data.table)

genes <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/genes"
length <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/human_sars_length.txt", sep= "\t")
input_output <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/MiST/input_file.txt"

files <- dir(genes, pattern = "*.csv")
files_dir <- paste(genes, "/", files, sep = "")
proteins <- files_dir %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package    
    purrr::reduce(rbind)           # reduce with rbind into one dataframe

accession_length <- length %>% as_tibble %>%
    separate(V1, into=c("pre", "accession", "post"), sep=">>") %>% 
    rename(proteinlength = V2) %>% select(accession,proteinlength)

proteins_length <- merge(proteins, accession_length, by="accession") %>% as_tibble

inputFile <- proteins_length %>% 
    mutate(hash = "hash") %>%
    mutate(BaitSims = "hash") %>%
    select(accession, hash, proteinlength, BaitSims, Bioreplicate, number_of_psms) %>%
    spread(Bioreplicate, number_of_psms, fill=0)

fwrite(inputFile, input_output, sep= "\t", col.names=TRUE)

##########################################################################################################################################

# After this two extra lines need to be added and all the "hash" need to be replaced by "#"

# https://modbase.compbio.ucsf.edu/mist/

# HIV trained and no singleton filtering