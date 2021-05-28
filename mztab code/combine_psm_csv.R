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

#	Get non-modified peptides
hcip_nonmod <- hcip_psm %>% filter(!PSM_ID %in% hcip_mod$PSM_ID)

#	Look at specific preys
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	select(V2, gene_name)

psm_genes <- merge(BP_psm, gene_names, by.x = "Prey", by.y = "V2") %>% as_tibble
psm_genes_mod <- merge(BP_psm_mod, gene_names, by.x = "Prey", by.y = "V2") %>% as_tibble


#   There are 6 BPs that do not contain any PSMs with a modifications
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
SARS_mod %>% filter(str_detect(mod, "Val->Ala")) %>% select(Condition, sequence.x, PSM_ID, mass_diff) %>% unique %>% print(n=100)
sno <- SARS_mod %>% filter(str_detect(mod, "Val->Ala")) %>% select(mass_diff, PSM_ID) %>% arrange(mass_diff)
sno %>% pull(mass_diff)
SARS_mod %>% filter(str_detect(mod, "Val->Ala")) %>% select(PSM_ID, mod) %>% print(n=100)
SARS_mod %>% filter(str_detect(mod, "Val->Ala")) %>% pull(mod) %>% unique

SARS_mod %>% filter(str_detect(mod, "Asn->Gly")) %>% select(Condition, sequence.x, PSM_ID, mass_diff) %>% unique %>% print(n=100)
SARS_mod %>% filter(str_detect(mod, "Asn->Gly")) %>% pull(sequence.x) %>% unique
SARS_psms %>% select(Condition, sequence) %>% 
	filter(Condition == "nsp5") %>%
	unique %>% print(n=100)

SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30) 
SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% select(mod, PSM_ID) %>% unique %>% print(n=30) 
SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% pull(mass_diff)

SARS_mod %>% filter(mass_diff > -28.0320 & mass_diff < -28.0304) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30) 

SARS_mod %>% filter(str_detect(mod, "Gluratylation")) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30)
SARS_mod %>% filter(str_detect(mod, "GGQ")) %>% select(Condition, sequence.x, PSM_ID, mass_diff) %>% unique %>% print(n=30)

#	Get a unmodified PSM
SARS_psms %>% filter(!PSM_ID %in% SARS_mod$PSM_ID) %>% 
	filter(prot == "orf9b", sequence == "KTLNSLEDK", bait == "qx017206") %>% select(Condition, PSM_ID) 

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

################################################################################################################################################

SARS_Shepherd_input <- SARS_psms %>% mutate(Spectrum = PSM_ID) %>% mutate(SpectrumFile = bait) %>%
	mutate(Peptide = sequence) %>% mutate(ModifiedPeptide = "NA") %>% mutate(PeptideLength = str_length(Peptide)) %>%
	
	select(Spectrum, SpectrumFile, Peptide, ModifiedPeptide, PeptideLength)


################################################################################################################################################

#	How many PSMs are there?				671,191
psm_annotations %>% select(PSM_ID) %>% unique
#	How many PSMs contain a Cys?			106,462
psm_annotations %>% filter(str_detect(sequence, "C")) %>% select(PSM_ID) %>% unique

#	How many modified PSMs are there?		307,641
psm_mod_annotation %>% select(PSM_ID) %>% unique
#	How many with Cys?						58,821
psm_mod_annotation %>% filter(str_detect(sequence.x, "C")) %>% select(PSM_ID) %>% unique

psm_cys <- psm_annotations %>% filter(str_detect(sequence, "C")) %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
sum(psm_cys_count) 

Cys_1_CAM <- psm_mod_annotation %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -57.021464 & mass_tol_neg < -57.021464) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_1_CAM_count <- str_count(Cys_1_CAM, pattern = "C")
table(Cys_1_CAM_count)
sum(Cys_1_CAM_count)

Cys_2_CAM <- psm_mod_annotation %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -114.0429 & mass_tol_neg < -114.0429) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_2_CAM_count <- str_count(Cys_2_CAM, pattern = "C")
table(Cys_2_CAM_count)


#	How many SARS PSMs?						40,828
SARS_psms %>% select(PSM_ID) %>% unique
#	How many PSMs contain a Cys?			17,858
SARS_psms %>% filter(str_detect(sequence, "C")) %>% select(PSM_ID) %>% unique

#	How many modified PSMs are there?		26,399
SARS_mod %>% select(PSM_ID) %>% unique
#	How many with Cys?						13,448
SARS_mod %>% filter(str_detect(sequence.x, "C")) %>% select(PSM_ID) %>% unique

psm_cys <- SARS_psms %>% filter(str_detect(sequence, "C")) %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
sum(psm_cys_count) 

Cys_1_CAM <- SARS_mod %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -57.021464 & mass_tol_neg < -57.021464) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_1_CAM_count <- str_count(Cys_1_CAM, pattern = "C")
table(Cys_1_CAM_count)
sum(Cys_1_CAM_count)

Cys_2_CAM <- SARS_mod %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -114.0429 & mass_tol_neg < -114.0429) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_2_CAM_count <- str_count(Cys_2_CAM, pattern = "C")
table(Cys_2_CAM_count)

Cys_3_CAM <- SARS_mod %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -171.0644 & mass_tol_neg < -171.0644) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_3_CAM_count <- str_count(Cys_3_CAM, pattern = "C")
table(Cys_3_CAM_count)

#	Count CAM Cys

csv_c_map <- csv %>% select(sequence, PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)
psm_c_map <- psm_csv %>% select(bait, accession, retention_time, exp_mass_to_charge)
psm_c_id <- merge(psm_c_map, csv_c_map) %>% as_tibble
psm_c_mod <- merge(psm_c_id, mod_data, by="PSM_ID") %>% as_tibble
psm_c_mod_annotation <- merge(psm_c_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%      #Filter out samples that were not included in the hcip analysis.
	filter(!accession=="null")
psm_c_annotations <- merge(psm_c_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null")

psm_c_annotations %>% select(sequence, PSM_ID) %>% unique %>% filter(str_detect(sequence, "C")) %>%
	filter(str_detect(sequence, "[160]"))

#	106,462 cys -> 106,462 [160]

psm_cys <- psm_c_annotations %>% filter(str_detect(sequence, "C")) %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
psm_160_count <- str_count(psm_cys, fixed("[160]"))
sum(psm_cys_count) 
sum(psm_160_count)


