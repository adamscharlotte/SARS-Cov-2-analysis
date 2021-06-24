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
    mutate_all(~gsub("SARS_CoV_2_|_low|NoRE", "", .)) %>%
    mutate_all(~gsub("protein14", "orf9c(protein14)", .)) %>%
    select(accession, hash, proteinlength, BaitSims, Bioreplicate, number_of_psms) %>%
    spread(Bioreplicate, number_of_psms, fill=0)

fwrite(inputFile, input_output, sep= "\t", col.names=TRUE)

##########################################################################################################################################

# After this two extra lines need to be added and all the "hash" need to be replaced by "#"

# #	#	#	Exps	E-1	E-2	E-3	EGFP-1	EGFP-2	EGFP-4	EGFP-5	EV-1	EV-2	EV-3	EV-4	EV-5	EV-6	EV-7	M-1	M-2	M-3	N-1	N-2	N-3	nsp1-1	nsp1-2	nsp1-3	nsp10-1	nsp10-2	nsp10-3	nsp11-1	nsp11-2	nsp11-3	nsp12-1	nsp12-2	nsp12-3	nsp13-1	nsp13-2	nsp13-3	nsp14-1	nsp14-2	nsp14-3	nsp15-1	nsp15-2	nsp15-3	nsp2-1	nsp2-2	nsp2-3	nsp4-1	nsp4-2	nsp4-3	nsp5_C145A-1	nsp5_C145A-2	nsp5_C145A-3	nsp5-1	nsp5-2	nsp5-3	nsp6-1	nsp6-2	nsp6-3	nsp7-1	nsp7-2	nsp7-3	nsp8-1	nsp8-2	nsp8-3	nsp9-1	nsp9-2	nsp9-3	orf10-1	orf10-2	orf10-3	orf3a-1	orf3a-2	orf3a-3	orf3b-1	orf3b-2	orf3b-3	orf6-1	orf6-2	orf6-3	orf7a-1	orf7a-2	orf7a-3	orf8-1	orf8-2	orf8-3	orf9b_low-1	orf9b_low-2	orf9b_low-3	orf9c(protein14)-1	orf9c(protein14)-2	orf9c(protein14)-3	S-1	S-2	S-3
# #	#	#	Baits	E	E	E	EGFP	EGFP	EGFP	EGFP	EV	EV	EV	EV	EV	EV	EV	M	M	M	N	N	N	nsp1	nsp1	nsp1	nsp10	nsp10	nsp10	nsp11	nsp11	nsp11	nsp12	nsp12	nsp12	nsp13	nsp13	nsp13	nsp14	nsp14	nsp14	nsp15	nsp15	nsp15	nsp2	nsp2	nsp2	nsp4	nsp4	nsp4	nsp5_C145A	nsp5_C145A	nsp5_C145A	nsp5	nsp5	nsp5	nsp6	nsp6	nsp6	nsp7	nsp7	nsp7	nsp8	nsp8	nsp8	nsp9	nsp9	nsp9	orf10	orf10	orf10	orf3a	orf3a	orf3a	orf3b	orf3b	orf3b	orf6	orf6	orf6	orf7a	orf7a	orf7a	orf8	orf8	orf8	orf9b	orf9b	orf9b	orf9c(protein14)	orf9c(protein14)	orf9c(protein14)	S	S	S
# accession	#	Length	BaitSims	E	E	E	EGFP	EGFP	EGFP	EGFP	EV	EV	EV	EV	EV	EV	EV	M	M	M	N	N	N	nsp1	nsp1	nsp1	nsp10	nsp10	nsp10	nsp11	nsp11	nsp11	nsp12	nsp12	nsp12	nsp13	nsp13	nsp13	nsp14	nsp14	nsp14	nsp15	nsp15	nsp15	nsp2	nsp2	nsp2	nsp4	nsp4	nsp4	nsp5_C145A	nsp5_C145A	nsp5_C145A	nsp5	nsp5	nsp5	nsp6	nsp6	nsp6	nsp7	nsp7	nsp7	nsp8	nsp8	nsp8	nsp9	nsp9	nsp9	orf10	orf10	orf10	orf3a	orf3a	orf3a	orf3b	orf3b	orf3b	orf6	orf6	orf6	orf7a	orf7a	orf7a	orf8	orf8	orf8	orf9b	orf9b	orf9b	orf9c(protein14)	orf9c(protein14)	orf9c(protein14)	S	S	S

# https://modbase.compbio.ucsf.edu/mist/

# HIV trained and no singleton filtering