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

csv_map <- csv %>% select(PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)
psm_map <- psm_csv %>% select(bait, sequence, accession, retention_time, exp_mass_to_charge)
psm_id <- merge(psm_map, csv_map) %>% as_tibble

#	Load Gordon et al. PSMs
msms_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/msms.txt") %>% as_tibble
proteins_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/proteinGroups.txt") %>% as_tibble

#	Increase in spectral identifications
psm_id %>% select(PSM_ID, sequence) %>% unique
msms_gordon %>% select(RawFile, ScanNumber, Sequence, Reverse) %>% 
	filter(Reverse=="") %>% unique
830743/387535

#	Load modification data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/modifications_10ppm_expmass"
setwd(path)
files <- dir(pattern = "*.csv")
mod_data <- files %>%
	map(read_csv) %>%       # read in all the files individually, using
							# the function read_csv() from the readr package
	reduce(rbind) %>%       # reduce with rbind into one dataframe
	select(bait, mod, mod_mass, mass_diff, everything())

psm_mod <- merge(psm_id, mod_data, by="PSM_ID") %>% as_tibble

psm_mod %>% select(PSM_ID) %>% unique
psm_id %>% filter(!PSM_ID %in% mod_data$PSM_ID) %>% select(PSM_ID) %>% unique

psm_mod %>% count(mass_diff) %>% arrange(desc(n))
psm_mod %>% count(mod) %>% arrange(desc(n))

#   Add annotation
annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

annsolo_annotation <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%      # Filter out samples that were not included in the original analysis
	filter(!accession=="null")

gordon_annotation <- merge(msms_gordon, annotation, by.x = "RawFile", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%      # Filter out samples that were not included in the original analysis
	filter(!Proteins=="")

psm_annsolo <- annsolo_annotation %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait, scan, sequence, retention_time, scanID, accession) %>% 
	unite(SPEC, bait:scan, sep="|", remove=FALSE) %>%
	select(bait, scan, sequence, everything()) %>%
	unite(PSM_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(RawFile = bait) %>% unique()

psm_gordon <- gordon_annotation %>% select(RawFile, ScanNumber, Sequence, Retentiontime, Proteins, Reverse) %>%
	unite(SPEC, RawFile:ScanNumber, sep="|", remove=FALSE) %>% 
	select(RawFile, ScanNumber, Sequence, everything()) %>%
	unite(PSM_ID, RawFile:Sequence, sep="|", remove=FALSE) %>% 
	filter(!Reverse=="+") %>% unique()

scan_overlap <- psm_annsolo %>% filter(SPEC %in% psm_gordon$SPEC) %>% select(RawFile, PSM_ID, SPEC) %>% unique
scan_overlap_g <- psm_gordon %>% filter(SPEC %in% psm_annsolo$SPEC) %>% select(RawFile, PSM_ID, SPEC, Reverse) %>% unique
scan_overlap_g %>% count(Reverse)

ann <- psm_annsolo %>% filter(!PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "ANN-SoLo") %>% mutate(Level = "PSM") %>% unique %>% rename(ID=PSM_ID)
gordon <- psm_gordon %>% filter(!PSM_ID %in% psm_annsolo$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "Gordon et al.") %>% mutate(Level = "PSM") %>% unique %>% rename(ID=PSM_ID)
both <- psm_annsolo %>% filter(PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "Overlap") %>% mutate(Level = "PSM") %>% unique %>% rename(ID=PSM_ID)

comb_psm <- rbind(ann, gordon, both) 

psm <- ggplot(comb_psm, aes(Level, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="PSM overlap", x = "Level", y = "Number of PSMs") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563"), name="") +
	scale_y_continuous(labels = scales::comma) +
	theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank())

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PSM overlap baplot.png", psm, width = 7.5, height = 8, units = "cm")

pep_annsolo <- annsolo_annotation %>% separate(PSM_ID, into=c("a", "scan", "b", "c")) %>%
	select(bait, sequence, retention_time) %>% unite(PEP_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(RawFile = bait)

pep_gordon <- gordon_annotation %>% select(RawFile, Sequence, Retentiontime, Reverse) %>%
	filter(!Reverse=="+") %>%
	unite(PEP_ID, RawFile:Sequence, sep="|", remove=FALSE)

ann <- pep_annsolo %>% filter(!PEP_ID %in% pep_gordon$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "ANN-SoLo") %>% mutate(Level = "Peptide") %>% unique %>% rename(ID=PEP_ID)
gordon <- pep_gordon %>% filter(!PEP_ID %in% pep_annsolo$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "Gordon et al.") %>% mutate(Level = "Peptide") %>% unique %>% rename(ID=PEP_ID)
both <- pep_annsolo %>% filter(PEP_ID %in% pep_gordon$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "Overlap") %>% mutate(Level = "Peptide") %>% unique %>% rename(ID=PEP_ID)

comb_pep <- rbind(ann, gordon, both)
comb_overlap <- rbind(comb_psm, comb_pep)

bar <- ggplot(comb_overlap, aes(Level, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="Overlap", x = "Level", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563"), name="") +
	theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank())

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/Overlap baplot.png", bar, width = 10, height = 19, units = "cm")

prot_annsolo <- annsolo_annotation %>% separate(accession, into=c("pre", "Protein", "post")) %>%
	select(bait, Protein, sequence) %>% unite(PROT_ID, bait:Protein, sep="|", remove=FALSE) %>%
	rename(RawFile = bait) %>% rename(Sequence = sequence) %>% unique

prot_gordon <- gordon_annotation %>% select(RawFile, Proteins, Sequence, Reverse) %>%
	filter(!Reverse=="+") %>%
	unite(PEP_ID, RawFile:Proteins, sep="|", remove=FALSE) %>% unique

#   Load HCIP result file
HCIP_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv"
hcip <- fread(HCIP_path) %>% as_tibble

##############################################################################################################################

#	Load PPI data
SAINT_MIST_filtered <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv") %>% as_tibble
Gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/intact_gordon.txt", sep= "\t") %>% as_tibble
Gordon_map <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results//HCIP/AP-MS/Gordon/name_map_sars_cov_2.txt", sep= "\t")
Gordon_annot <- merge(Gordon, Gordon_map, by.x = "IDinteractorA", by.y = "accession") %>% as_tibble
#   remove space in protein14 (bait)
Gordon_annot$idBait <- gsub(" ", "", Gordon_annot$idBait)
Gordon_gene <- merge(Gordon_annot, gene_names, by.x = "IDinteractorB", by.y = "V2") %>% as_tibble

Gordon_result <- Gordon_gene %>% mutate(Prey = IDinteractorB) %>% mutate(PreyGene = gene_name) %>%
	select(idBait, Prey, everything()) %>% unite(BP, idBait:Prey, sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene, sep ="_", remove=FALSE)

ann <- SAINT_MIST_filtered %>% filter(!BP_gene %in% Gordon_result$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "ANN-SoLo") %>% mutate(Level = "PSM") %>% unique
gordon <- Gordon_result %>% filter(!BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "Gordon et al.") %>% mutate(Level = "PSM") %>% unique
both <- SAINT_MIST_filtered %>% filter(BP_gene %in% Gordon_result$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "Overlap") %>% mutate(Level = "PSM") %>% unique

comb_ppi <- rbind(ann, gordon, both) 

ppi <- ggplot(comb_ppi, aes(Level, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="PPI overlap", x = "Level", y = "Number of PPIs") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563"), name="") +
	scale_y_continuous(labels = scales::comma) +
	theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank())

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap baplot.png", ppi, width = 7, height = 8, units = "cm")
