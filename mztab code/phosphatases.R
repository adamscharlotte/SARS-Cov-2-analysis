library(tidyverse)
library(data.table)

#   Load PSM data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl"
setwd(path)
files <- dir(pattern = "*.csv")
psm_csv <- files %>%
	map(read_csv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind)			# reduce with rbind into one dataframe

path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/csv/psm"
setwd(path)
files <- dir(pattern = "*.csv")
csv <- files %>%
	map(read_tsv) %>%       # read in all the files individually, using
							# the function read_tsv() from the readr package
	reduce(rbind)           # reduce with rbind into one dataframe

#   Load modification data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/modifications_10ppm_expmass"
setwd(path)
files <- dir(pattern = "*.csv")
mod_data <- files %>%
	map(read_csv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind) %>%       # reduce with rbind into one dataframe
	select(bait, mod, mod_mass, mass_diff, everything())

csv_map <- csv %>% select(PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)
psm_map <- psm_csv %>% select(bait, sequence, accession, retention_time, exp_mass_to_charge)
psm_id <- merge(psm_map, csv_map) %>% as_tibble
psm_mod <- merge(psm_id, mod_data, by="PSM_ID") %>% as_tibble
psm_id %>% filter(!PSM_ID %in% psm_mod$PSM_ID) %>% select(PSM_ID) %>% unique

#   Add annotation
annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

psm_mod_annotation <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") #%>%      #Filter out samples that were not included in the hcip analysis.
	#filter(!accession=="null")
psm_annotations <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null")

#	Load HCIP results
HCIP_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv"
hcip <- fread(HCIP_path) %>% as_tibble


colnames(hcip)
colnames(psm_annotations)

#   Create column with BP_gene
psm_mod_annotation %>% select(accession) %>% unique

BP_psm_mod <- psm_mod_annotation %>% 
	separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
	filter(!name=="SARS") %>%
	separate(accession, into = c("pre", "Prey", "gene")) %>%
	separate(gene, into = c("PreyGene", "human")) %>%
	select(Condition, PreyGene, everything()) %>%
	unite(BP_gene, Condition:PreyGene ,sep ="_", remove=FALSE) %>% unique()
BP_psm <-  psm_annotations %>%
	separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
	filter(!name=="SARS") %>%
	separate(accession, into = c("pre", "Prey", "gene")) %>%
	separate(gene, into = c("PreyGene", "human")) %>%
	select(Condition, PreyGene, everything()) %>%
	unite(BP_gene, Condition:PreyGene ,sep ="_", remove=FALSE) %>% unique()

Orf8	RNGTT
Orf8	DUSP1
Nsp8	PPM1G
Nsp13	PTPN13
Nsp7	PTPMT1
