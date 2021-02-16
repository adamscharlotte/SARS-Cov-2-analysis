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
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/modifications_10ppm"
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

#   Add annotation
annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

psm_mod_annotation <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%      #Filter out samples that were not included in the hcip analysis.
	filter(!accession=="null")
psm_annotations <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null")

#   Create an overview of the modifications present in the PSMs that could be matched to proteins


################################################################################################################################################

#   Only select modifications present in reported HCIPs
#   Load HCIP result file
HCIP_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv"
hcip <- fread(HCIP_path) %>% as_tibble

#   Create column with BPAccession
BP_hcip <- hcip %>% mutate(BPAccession = BP)
BP_psm_mod <- psm_mod_annotation %>% 
	separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
	filter(!name=="SARS") %>%
	separate(accession, into = c("pre", "Prey", "gene")) %>%
	select(Condition, Prey, everything()) %>%
	unite(BPAccession, Condition:Prey ,sep ="_", remove=FALSE) %>% unique
BP_psm <-  psm_annotations %>%
	separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
	filter(!name=="SARS") %>%
	separate(accession, into = c("pre", "Prey", "gene")) %>%
	select(Condition, Prey, everything()) %>%
	unite(BPAccession, Condition:Prey ,sep ="_", remove=FALSE) %>% unique

#   Link the modification information to the hcip information
hcip_mod <- merge(BP_psm_mod, BP_hcip, by="BPAccession") %>% as_tibble

#   How many PSMs belong to a HCIP
hcip_psm <- merge(BP_psm, BP_hcip, by="BPAccession") %>% as_tibble

hcip_psm %>% filter(!BP_gene %in% hcip_mod$BP_gene) %>% select(BP_gene) %>% unique
hcip_psm %>% select(BP_gene, idBait) %>% filter(idBait == "nsp5"| idBait == "nsp5_C145A") %>% unique
hcip_psm %>% filter(BP_gene == "nsp5_C145A_SORD") %>% select(Prey.x) %>% unique

hcip_psm %>% filter(BP_gene == "orf8_CNNM3") %>% select(Spec, sequence) %>% unique
hcip_mod %>% filter(BP_gene == "orf8_PLOD1") %>% select(mass_diff, sequence.x, PSM_ID, mod) %>% arrange(mass_diff) %>% print(n=200)
hcip_mod %>% filter(BP == "nsp5_C145A_Q6YN16") %>% count(mod) %>% arrange(desc(n)) %>% print(n=50)

hcip_mod %>% filter(BP_gene == "orf3a_MGAT1") %>% count(bait.x)
hcip_psm %>% filter(BP_gene == "orf3a_MGAT1") %>% count(bait, BPAccession)
hcip_psm %>% filter(BP_gene == "orf3a_MGAT1") %>% count(bait)
select(sequence, bait)
%>% select(Prey.x, BP_gene) %>% unique
hcip_mod %>% filter(BP_gene == "orf3a_MGAT1") %>% select(mass_diff, mod, bait.x, sequence.x)
BP_psm_mod %>% filter(BPAccession == "orf3a_P26572") %>% select(mass_diff, mod, bait.x, sequence.x)
BP_psm %>% filter(BPAccession == "orf3a_P26572") %>% select(bait, sequence)

#	Get non-modified peptides
hcip_nonmod <- hcip_psm %>% filter(!PSM_ID %in% hcip_mod$PSM_ID)
hcip_nonmod %>% filter(BP == "nsp5_C145A_Q9NXH9") %>% select(sequence, PSM_ID) %>% print(n=30)
hcip_psm %>% filter(BP == "nsp5_C145A_Q9NXH9") %>% count(sequence)

#   There are 50 BPs that do not contain any PSMs with a modifications
BP_nomod <- BP_hcip %>% filter(!BP_gene %in% hcip_mod$BP_gene) %>% select(BP_gene) %>% unique

################################################################################################################################################

#   Look at modifications in the SARS-CoV-2 proteins
SARS_psms <- psm_annotations %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14")) %>%
	filter(str_detect(Condition, prot))     #So nsp5 is matched to nsp5 and nsp5_C145A
SARS_mod <- psm_mod_annotation %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14")) %>%
	filter(str_detect(Condition, prot))

#   Look at modifications per SARS-CoV-2 protein
SARS_mod %>% count(Condition) %>% arrange(desc(n)) %>% print(n=30)
SARS_mod %>% filter(mod == "Sulfation / Phospho") %>% select(Condition, sequence.x, PSM_ID) %>% unique
SARS_mod %>% filter(str_detect(mod, "SNO")) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=100)

SARS_mod %>% filter(str_detect(mod, "SNO") & Condition=="N") %>% select(mod_mass, mod) %>% unique %>% print(n=100)
SARS_mod %>% filter(mod == "SNO / Val->Gln") %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=100)

SARS_mod %>% filter(str_detect(mod, "SNO") & Condition=="N") %>% pull(mass_diff)

SARS_mod %>% filter(str_detect(mod, "Gluratylation")) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30)
SARS_mod %>% filter(str_detect(mod, "GGQ")) %>% select(Condition, sequence.x, PSM_ID, mass_diff) %>% unique %>% print(n=30)

#	Get a unmodified PSM
SARS_psms %>% filter(!PSM_ID %in% SARS_mod$PSM_ID) %>% 
	filter(prot == "nsp15", sequence == "KPTETICAPLTVFFDGR", bait == "qx017215") %>% select(Condition, PSM_ID) 

#	General overview
SARS_mod %>% count(mod) %>% 
	arrange(desc(n)) %>%
	filter(n>500) %>% pull(mod)

SARS_mod %>% select(sequence.x, mod, everything()) %>%
	unite(mod_seq, sequence.x:mod, remove=FALSE) %>%
	filter(!mod=="No direct match found in Unimod") %>%
	filter(!str_detect(mod, "Deamidation")) %>%
	count(mod_seq) %>% arrange(desc(n)) %>% 
	# filter(n>=3) %>%
	# pull(mod_seq)  
	print(n=50) 

SARS_mod %>% filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% mutate(round_mass=round(mass_diff, 2)) %>%
	count(round_mass) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "E") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
	filter(sequence.x == "FNGIGVTQNVLYENQK") %>% count(mod) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "nsp9") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
	filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% count(mod) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "E") %>% select(sequence.y, mass_tol, mass_diff, mod) %>% pull(mod)

SARS_mod %>% as_tibble %>%
	# filter(str_detect(mod, "Gluratylation")) %>%
	select(sequence.y, mass_tol, mod, prot) %>% #count(prot)
	filter(prot == "orf3a") %>% #pull(sequence.y)
	pull(mod)
	unique %>% print(n=30)
