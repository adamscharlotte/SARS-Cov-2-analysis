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

################################################################################################################################################

#	Mass difference histogram

psm_mod_his <- psm_mod %>% select(PSM_ID, mod, mod_mass, mass_diff) %>% unique %>%
	mutate(mass_diff_rounded = round(mass_diff, digits = 0)) %>% filter(!mass_diff_rounded==0)
	# filter(!mod=="No direct match found in Unimod")
	# filter(!(mass_diff > -0.05 & mass_diff < 0.05))
psm_mod_his %>% select(PSM_ID) %>% unique

ylabel <- paste("Number of PSMs")

mass_diff_histogram <- ggplot(psm_mod_his, aes(mass_diff)) + geom_histogram(binwidth = 0.2) + 
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank()
		) +
	scale_x_continuous(limits=c(-50,120), n.breaks = 8) +
	scale_y_continuous(n.breaks = 7, labels = scales::comma) +
	xlab('Precursor mass difference (Da)') + 
	geom_hline(yintercept = 0, size = 0.1) +
	# geom_vline(xintercept = 114) +
	ylab(ylabel)

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Precursor mass difference histogram.pdf", mass_diff_histogram,  width = 22, height = 10, units = "cm")

psm_mod_his %>% count(mass_diff) %>% arrange(desc(n))
psm_mod_his %>% filter(mass_diff> 31 & mass_diff< 33) %>% # count(mod) %>% arrange(desc(n))
	select(mod, mod_mass) %>% unique %>% print(n=50)

#   Create an overview of the modifications present in the PSMs that could be matched to proteins


################################################################################################################################################

#   Look at modifications in the SARS-CoV-2 proteins
SARS_psms <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14"))

SARS_mod <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14"))

#	Mass difference histogram
psm_mod_his <- SARS_mod %>% select(PSM_ID, mod, mod_mass, mass_diff) %>% unique %>%
	mutate(mass_diff_rounded = round(mass_diff, digits = 0)) %>% filter(!mass_diff_rounded==0)
	# filter(!mod=="No direct match found in Unimod")
	# filter(!(mass_diff > -0.05 & mass_diff < 0.05))
psm_mod_his %>% select(PSM_ID) %>% unique

ylabel <- paste("Number of PSMs")

mass_diff_histogram <- ggplot(psm_mod_his, aes(mass_diff)) + geom_histogram(binwidth = 0.2) + 
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank()
		) +
	scale_x_continuous(limits=c(-50,120), n.breaks = 8) +
	scale_y_continuous(n.breaks = 7, labels = scales::comma) +
	xlab('Precursor mass difference (Da)') + 
	geom_hline(yintercept = 0, size = 0.1) +
	# geom_vline(xintercept = 114) +
	ylab(ylabel)

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Precursor mass difference histogram SARS.pdf", mass_diff_histogram,  width = 22, height = 8, units = "cm")

#	Phosphorylation on SARS-CoV-2 proteins
phospho_sars <- SARS_mod %>% filter(mod=="Phospho"|mod=="Sulfation / Phospho") %>% 
	select(Condition, prot, sequence.x,bait.x, PSM_ID, mass_diff) %>%
	unite(seqbait, sequence.x:bait.x, remove=FALSE) %>%
	filter(str_detect(Condition, prot)) %>% unique %>% print(n=40)
# phospho_sars <- SARS_mod %>% filter(mod=="Phospho") %>% 
# 	select(Condition, prot, sequence.x,bait.x, PSM_ID) %>%
# 	unite(seqbait, sequence.x:bait.x, remove=FALSE) %>%
# 	filter(str_detect(Condition, prot)) %>% unique %>% pull(PSM_ID)

#	Ubiquitination on SARS-CoV-2 proteins
ub_sars <- SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% 
	select(Condition, prot, sequence.x, bait.x, PSM_ID, mod) %>% 
	unite(seqbait, sequence.x:bait.x, remove=FALSE) %>%
	unique %>% arrange(Condition) %>% print(n=40)
ub_sars %>% filter(Condition=="nsp5_C145A") %>%
	pull(mod)

#	S-nitrosylation on SARS-CoV-2 proteins
sno_sars <- SARS_mod %>% filter(str_detect(mod, "Val->Ala")) %>% 
	select(Condition, prot, sequence.x, bait.x, PSM_ID, mass_diff) %>% 
	unite(seqbait, sequence.x:bait.x, remove=FALSE) %>%
	unique %>% arrange(prot) 

#	Count unmodified PSMs for each 
tbl_unmodified_psms <- SARS_psms %>% filter(!PSM_ID %in% SARS_mod$PSM_ID) %>%
	select(Condition, sequence, everything()) %>%
	unite(protseq, Condition:sequence) %>% unique()
tbl_phospho_psms <- phospho_sars %>% select(Condition, sequence.x, everything()) %>%
	unite(protseq, Condition:sequence.x) %>% unique()
tbl_ub_psms <- ub_sars %>% select(Condition, sequence.x, everything()) %>%
	unite(protseq, Condition:sequence.x) %>% unique()
tbl_sno_psms <- sno_sars %>% select(Condition, sequence.x, everything()) %>%
	unite(protseq, Condition:sequence.x) %>% unique()

tbl_unmodified_psms %>% filter(protseq %in% tbl_phospho_psms$protseq) %>%
	count(protseq)
tbl_unmodified_psms %>% filter(protseq %in% tbl_ub_psms$protseq) %>%
	count(protseq)
tbl_unmodified_psms %>% filter(protseq %in% tbl_sno_psms$protseq) %>%
	count(protseq)

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
SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% dplyr::select(mod, PSM_ID) %>% unique %>% print(n=30) 
SARS_mod %>% filter(str_detect(mod, "GlyGly")) %>% pull(mass_diff)

SARS_mod %>% filter(mass_diff > -28.0320 & mass_diff < -28.0304) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30) 

SARS_mod %>% filter(str_detect(mod, "Gluratylation")) %>% select(Condition, sequence.x, PSM_ID) %>% unique %>% print(n=30)
SARS_mod %>% filter(str_detect(mod, "GGQ")) %>% select(Condition, sequence.x, PSM_ID, mass_diff) %>% unique %>% print(n=30)

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

#	How many PSMs? 							830,743
psm_id %>% select(PSM_ID) %>% unique
#	How many modified PSMs? 				402,586
psm_mod %>% select(PSM_ID) %>% unique
#	How many PSMs contain a Cys?			125,727
psm_id %>% filter(str_detect(sequence, "C")) %>% select(PSM_ID) %>% unique
#	How many modified PSMs with Cys?		58,821
psm_mod %>% filter(str_detect(sequence.x, "C")) %>% select(PSM_ID) %>% unique

#	How many cysteines are there?			152,766
psm_cys <- psm_id %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
sum(psm_cys_count)

Cys_1_CAM <- psm_mod %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -57.021464 & mass_tol_neg < -57.021464) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_1_CAM_count <- str_count(Cys_1_CAM, pattern = "C")
table(Cys_1_CAM_count)
sum(Cys_1_CAM_count)

Cys_2_CAM <- psm_mod %>% filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -114.0429 & mass_tol_neg < -114.0429) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_2_CAM_count <- str_count(Cys_2_CAM, pattern = "C")
table(Cys_2_CAM_count)

#	How many SARS PSMs?						43,781
psm_id %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% select(PSM_ID) %>% unique
#	How many modified SARS PSMs?			27,983
psm_mod %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% select(PSM_ID) %>% unique
#	How many PSMs contain a Cys?			18,868
psm_id %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% 
	filter(str_detect(sequence, "C")) %>% select(PSM_ID) %>% unique
#	How many modified PSMs with Cys?		14,102
psm_mod %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% 
	filter(str_detect(sequence.x, "C")) %>% select(PSM_ID) %>% unique

#	How many cysteines are there?			25,506
psm_cys <- psm_id %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
sum(psm_cys_count)

Cys_1_CAM <- psm_mod %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% 
	filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -57.021464 & mass_tol_neg < -57.021464) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_1_CAM_count <- str_count(Cys_1_CAM, pattern = "C")
table(Cys_1_CAM_count)
sum(Cys_1_CAM_count)

Cys_2_CAM <- psm_mod %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% 
	filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -114.0429 & mass_tol_neg < -114.0429) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_2_CAM_count <- str_count(Cys_2_CAM, pattern = "C")
table(Cys_2_CAM_count)

Cys_3_CAM <- psm_mod %>% filter(str_detect(accession,"SARS_CoV_2_")) %>% 
	filter(str_detect(sequence.x, "C")) %>% filter(mass_tol_pos > -171.0644 & mass_tol_neg < -171.0644) %>% 
	select(sequence.x, PSM_ID) %>% unique %>% pull(sequence.x) 
Cys_3_CAM_count <- str_count(Cys_3_CAM, pattern = "C")
table(Cys_3_CAM_count)

#	Count CAM Cys
csv_c_map <- csv %>% select(sequence, PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)
psm_c_map <- psm_csv %>% select(bait, accession, retention_time, exp_mass_to_charge)
psm_c_id <- merge(psm_c_map, csv_c_map) %>% as_tibble

psm_cys <- psm_c_id %>% select(sequence, PSM_ID) %>% unique %>% pull(sequence)
psm_cys_count <- str_count(psm_cys, pattern = "C")
psm_160_count <- str_count(psm_cys, fixed("[160]"))
sum(psm_cys_count) 
sum(psm_160_count)

################################################################################################################################################

annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

SARS_psms_SNO <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14")) %>%
	filter(str_detect(Condition, prot))     #So nsp5 is matched to nsp5 and nsp5_C145A
SARS_mod_SNO <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14")) %>%
	filter(str_detect(Condition, prot))

SARS_mod_SNO %>% as_tibble %>% filter(Condition=="nsp5") %>% 
	filter(prot=="nsp5") %>% count(sequence.x) %>% arrange(desc(n))
SARS_mod_SNO %>% as_tibble %>% filter(Condition=="nsp5_C145A") %>% 
	filter(prot=="nsp5") %>% count(sequence.x) %>% arrange(desc(n))

SARS_mod_SNO %>% as_tibble %>% filter(Condition=="nsp5") %>% 
	filter(prot=="nsp5") %>% filter(sequence.x=="HVICTSEDMLNPNYEDLLIR") %>%
	arrange(mass_diff) %>% pull(mass_diff)

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


psm_genes_mod %>% filter(str_detect(mod, "Val->Ala")) %>% count(Condition) %>% arrange(desc(n))
psm_genes_mod %>% filter(str_detect(mod, "SNO")) %>% count(Condition) %>% arrange(desc(n))

psm_genes_mod %>% filter(str_detect(mod, "Val->Ala")) %>% 
	#filter(str_detect(sequence.x, "C")) %>% count(Condition) %>% arrange(desc(n))
	filter(Condition=="nsp15") %>% 
	#filter(gene_name=="CDHR5") %>%
	#count(gene_name) %>% arrange(desc(n))
	count(sequence.x) %>% arrange(desc(n))

psm_genes_mod %>% filter(str_detect(mod, "SNO")) %>% 
	filter(Condition=="nsp5") %>% 
	#filter(gene_name=="CDHR5") %>%
	#count(gene_name) %>% arrange(desc(n))
	count(sequence.x) %>% arrange(desc(n))

hcip_mod %>% filter(str_detect(mod, "Val->Ala")) %>% count(Condition)

#	Look at phosphatases
psm_genes_mod %>% filter(gene_name == "RNGTT") %>% 
	filter(Condition == "orf8") %>% 
	select(BPAccession, mass_diff, mod) %>% unique()

psm_genes_mod %>% filter(gene_name == "RNGTT") %>% 
	filter(Condition == "orf8") %>% 
	select(BPAccession, mass_diff, mod) %>% unique()

psm_genes_mod %>% filter(gene_name == "RNGTT") %>% 
	filter(Condition == "orf8") %>% 
	select(BPAccession, mass_diff, mod) %>% unique()

################################################################################################################################################
#	Look at ubiquitination in different samples
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	select(V2, gene_name)

#	Load MiST output
# 	Before you can load the MistOutput the first couple of lines must be removed.
MIST <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/MiST/MistOutput.txt", sep= "\t") %>% as_tibble
MIST_output <- merge(MIST, gene_names, by.x = "Prey", by.y = "V2") %>% 		# Add gene names
	as_tibble %>% mutate(PreyGene = gene_name) %>% mutate(idBait=Bait) %>%
	filter(!MiST==0) %>%
	select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene ,sep ="_", remove=FALSE)

#	Load SAINT output
SAINT_output <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/list.txt", sep= "\t") %>% as_tibble %>%
	mutate(idBait=Bait) %>% select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene ,sep ="_", remove=FALSE)
 
#	Combine the MIST and SAINT output
MIST_filter <- MIST_output %>% filter(BP %in% SAINT_output$BP)
SAINT_filter <- SAINT_output %>% filter(BP %in% MIST_output$BP)
SAINT_MIST <- merge(MIST_filter, SAINT_filter) %>% as_tibble

#	Add modifications
BP_psm_mod <- psm_mod_annotation %>% 
	separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
	filter(!name=="SARS") %>%
	separate(accession, into = c("pre", "Prey", "gene")) %>%
	select(Condition, Prey, everything()) %>%
	unite(BP, Condition:Prey ,sep ="_", remove=FALSE) %>% unique
,#   Link the modification information to the hcip information
SAINT_MIST_mod <- merge(BP_psm_mod, SAINT_MIST, by="BP") %>% as_tibble

#   How many PSMs belong to a HCIP
SAINT_MIST_psm <- merge(BP_psm, BP_hcip, by="BP") %>% as_tibble

#	Get non-modified peptides
SAINT_MIST_nonmod <- SAINT_MIST_psm %>% filter(!PSM_ID %in% hcip_mod$PSM_ID)

SAINT_MIST_mod %>% filter(gene == "PTPMT1")

SAINT_MIST_mod %>% filter(str_detect(mod, "GlyGly")) %>% count(Condition) %>% arrange(desc(n)) %>% print(n=30)

SAINT_MIST_mod %>% filter(str_detect(mod, "GlyGly")) %>% 
	filter(MiST >= 0.7 & BFDR <= 0.05 & AvgSpec >= 2) %>% 
	dplyr::select(BP_gene, exp_mass_to_charge.y, calc_mass_to_charge, mass_diff, mod)
	
	dplyr::select(BP_gene, PSM_ID, sequence.x, mass_diff, mod)

SAINT_MIST_mod %>% filter(str_detect(mod, "GlyGly")) %>%
	select(BP_gene)
	# count(sequence.x)
	filter(Prey.x=="P62979") %>% count(bait.x) %>% arrange(desc(n))

SAINT_MIST_nonmod %>% filter(Prey.x=="Q9HB71") %>% filter(sequence=="ISNYGWDQSDKFVK") %>% count(Condition)

SAINT_MIST_mod %>% filter(mod=="Phospho"|mod=="Sulfation / Phospho") %>% 
	filter(MiST >= 0.7 & BFDR <= 0.05 & AvgSpec >= 2) %>% 
	select(BP_gene, PSM_ID, sequence.x) %>% print(n=40)

SAINT_MIST_mod %>% 
	select(Condition, BP_gene) %>%
	unique() %>%
	count(Condition) %>% 
	arrange(n)

################################################################################################################################################

psm_mod %>% select(PSM_ID) %>% unique

mass_difference_table <- psm_mod %>% mutate(mass_diff_rounded = round(mass_diff, 2)) %>% 
	select(PSM_ID, mass_diff_rounded, mod) %>% unique() %>%
	filter(!(mass_diff_rounded < 0.5 & mass_diff_rounded > -0.5)) %>%
	add_count(mass_diff_rounded) %>%
	select(mass_diff_rounded, n, mod) %>% unique() %>%
	arrange(desc(mass_diff_rounded))

fwrite(mass_difference_table, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/mass_differences.csv", append = FALSE, col.names = TRUE)
