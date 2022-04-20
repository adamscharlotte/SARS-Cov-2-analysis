library(tidyverse)
library(data.table)
chooseCRANmirror(4)
install.packages("VennDiagram")
library(VennDiagram)

#	How many MSMSspectra are there?			 2,503,010 of which 376,009 are identified
msms_scan_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/msmsScans.txt") %>% as_tibble
msms_scan_gordon %>% count(Identified) %>% unique

#   Load ANN-SoLo PSM data
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

#	How many identified MS/MS spectra?		830,743 PSMs, 91,553 are not mapped to a protein
#	830743 - 91553 = 739190
psm_id %>% select(PSM_ID) %>% unique
psm_id %>% count(accession) %>% arrange(desc(n))

#	Load Gordon et al. PSMs
msms_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/msms.txt") %>% as_tibble

#	How many identified MS/MS spectra?		388417 PSMs, 51830 are not mapped to a protein
#	388417 - 51830 = 336587
msms_gordon %>% filter(str_detect(Proteins, "CON__")) %>% select(GeneNames, Proteins) %>% unique
msms_gordon %>% filter(str_detect(Proteins, "CON__")) %>% count(GeneNames) %>% arrange(desc(n))
msms_gordon %>% filter(GeneNames=="") %>% count(Proteins) %>% arrange(desc(n))

#	PSM overlap
psm_annsolo <- psm_id %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait, scan, sequence, retention_time, scanID, accession) %>% 
	unite(SPEC, bait:scan, sep="|", remove=FALSE) %>%
	select(bait, scan, sequence, everything()) %>%
	unite(PSM_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait) %>% unique()
psm_gordon <- msms_gordon %>% select(RawFile, ScanNumber, Sequence, Retentiontime, Proteins) %>%
	unite(SPEC, RawFile:ScanNumber, sep="|", remove=FALSE) %>% 
	select(RawFile, ScanNumber, Sequence, everything()) %>%
	unite(PSM_ID, RawFile:Sequence, sep="|", remove=FALSE) %>% 
	rename(Rawfile = RawFile) %>% unique()

ann <- psm_annsolo %>% filter(!PSM_ID %in% psm_gordon$PSM_ID) %>% select(Rawfile, PSM_ID, SPEC) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- psm_gordon %>% filter(!PSM_ID %in% psm_annsolo$PSM_ID) %>% select(Rawfile, PSM_ID, SPEC) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- psm_annsolo %>% filter(PSM_ID %in% psm_gordon$PSM_ID) %>% select(Rawfile, PSM_ID, SPEC) %>% mutate(Identification = "Overlap") %>% unique
psm_annsolo %>% select(PSM_ID, SPEC) %>% unique

comb_psm <- rbind(ann, gordon, both) 

psm <- ggplot(comb_psm, aes(RawFile, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="PSM overlap", x = "Raw file", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=9, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PSM overlap.png", psm, width = 31, height = 9, units = "cm")

################################################################################################################################################################

#	Mapping of spectra to multiple sequences
msms_gordon_SPEC <- msms_gordon %>% select(RawFile, ScanNumber, everything()) %>% 
	unite(SPEC, RawFile:ScanNumber, remove=FALSE)
	
SPEC_3 <- msms_gordon_SPEC %>% count(SPEC) %>% filter(n==3)
SPEC_2 <- msms_gordon_SPEC %>% count(SPEC) %>% filter(n==2)

msms_gordon_SPEC %>% filter(SPEC %in% SPEC_3$SPEC) %>% select(SPEC, Sequence, Proteins) %>% arrange(SPEC)

psm_gordon %>% filter(SPEC %in% SPEC_2$SPEC)

psm_gordon %>% add_count(SPEC) %>% filter(n>1) %>% count(Proteins) %>% arrange(desc(n))
gordon_unique_psm <- psm_gordon %>% add_count(SPEC) %>% filter(n==1) 

gordon_unique_psm %>% filter(Proteins=="")
gordon_unique_psm %>% filter((Proteins, ))
gordon_unique_psm %>% filter(Proteins=="")


psm_ann_SPEC <- psm_id %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait, scan, sequence, retention_time, scanID, accession) %>% unite(SPEC, bait:scan, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait) %>% unique()

# psm_ann_SPEC %>% select(SPEC, sequence) %>% unique %>% count(SPEC) %>% count(n)

#	GeneName, Proteins == ""
msms_gordon %>% filter(GeneNames=="")						# 51,830
msms_gordon %>% filter(GeneNames=="") %>% 
	filter(str_detect(Proteins, "CON__")) %>%
	select(RawFile, ScanNumber, Sequence) %>% unique		# 17,484
msms_gordon %>% filter(GeneNames=="") %>% 
	filter(str_detect(Proteins, "SARS_")) %>%
	select(RawFile, ScanNumber, Sequence) %>% unique		# 33,159

CON_overlap <- psm_gordon %>% filter(str_detect(Proteins, "CON__")) %>% filter(PSM_ID %in% psm_annsolo$PSM_ID)
psm_annsolo %>% filter(PSM_ID %in% CON_overlap$PSM_ID) %>% count(accession) %>% arrange(desc(n))

SARS_CoV_overlap <- psm_gordon %>% filter(str_detect(Proteins, "SARS_CoV")) %>% filter(PSM_ID %in% psm_annsolo$PSM_ID)
psm_annsolo %>% filter(PSM_ID %in% SARS_CoV_overlap$PSM_ID) %>% count(accession) %>% arrange(desc(n))
SARS_CoV_overlap %>% select(PSM_ID) %>% unique

psm_annsolo %>% filter(str_detect(accession, "SARS_CoV")) %>% select(PSM_ID) %>% unique
43781 - 18250
33126 - 18250

msms_gordon %>% filter(GeneNames=="") %>% count(Proteins) %>% arrange(desc(n)) %>%
	filter(!str_detect(Proteins, "CON__")) %>% filter(!str_detect(Proteins, "SARS_")) 

msms_gordon %>% filter(Proteins=="")


#	Create the venn diagram
venn <- draw.pairwise.venn(388342, 830743, 306131,
	category = rep(c("Gordon et al.", "ANN-SoLo")), cex = 2, cat.cex = 2,
	fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
	rscaled = TRUE, cat.dist = c(0.04, 0.04),
	cat.pos = c(-26, 10))
# grid.draw(venn)

ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/PSM overlap.pdf", width = 17, height = 17, units = "cm")

gordon_proteins <- psm_gordon %>% pull(Proteins) %>% unique()
ann_proteins <- psm_annsolo %>% separate(accession, into=c("prefix", "protein", "suff")) %>% select(protein) %>% unique
psm_annsolo %>% select(accession) %>% unique()
psm_annsolo %>% filter(sequence == "HAVSEGTK") %>% select(accession) %>% unique
psm_gordon %>% filter(Sequence == "HAVSEGTK") %>% select(Proteins) %>% unique
psm_gordon %>% filter(Sequence == "HAVSEGTK") %>% filter(str_detect(Proteins, ann_proteins)) %>% select(Proteins) %>% unique
psm_annsolo %>% select(accession)

psm_annsolo %>% select(PSM_ID) %>% unique

################################################################################################################################################################

# evidence_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/evidence.txt") %>% as_tibble
# evidence_gordon %>% select('MS/MScount') %>% unique
# msms_gordon %>% select(id) %>% unique()
# mscount <- evidence_gordon %>% pull('MS/MScount')
# sum(mscount)
# Peptide overlap

#	Peptide overlap
pep_annsolo <- psm_id %>% separate(PSM_ID, into=c("a", "scan", "b", "c")) %>%
	select(bait, sequence, retention_time) %>% unite(PEP_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(RawFile = bait)

pep_gordon <- msms_gordon %>% select(RawFile, Sequence, Retentiontime) %>%
	unite(PEP_ID, RawFile:Sequence, sep="|", remove=FALSE)

ann <- pep_annsolo %>% filter(!PEP_ID %in% pep_gordon$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- pep_gordon %>% filter(!PEP_ID %in% pep_annsolo$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- pep_annsolo %>% filter(PEP_ID %in% pep_gordon$PEP_ID) %>% select(RawFile, PEP_ID) %>% mutate(Identification = "Overlap") %>% unique

comb <- rbind(ann, gordon, both) 

s <- ggplot(comb, aes(RawFile, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="Peptide overlap", x = "Raw file", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=9, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/Peptide overlap.png", s,width = 31, height = 9, units = "cm")

##################################################################################################################################################################

#	identified spectra overlap
psm_annsolo <- psm_id %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait, scan, sequence, retention_time, scanID) %>% unite(PSM_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(RawFile = bait) %>% unique()
psm_gordon <- msms_gordon %>% select(RawFile, ScanNumber, Sequence, Retentiontime) %>%
	filter(!Sequence=="") %>%
	unite(PSM_ID, RawFile:Sequence, sep="|", remove=FALSE) %>% unique()

ann <- psm_annsolo %>% filter(!PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- psm_gordon %>% filter(!PSM_ID %in% psm_annsolo$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- psm_annsolo %>% filter(PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "Overlap") %>% unique

comb_psm <- rbind(ann, gordon, both) 

psm <- ggplot(comb_psm, aes(RawFile, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="PSM overlap", x = "Raw file", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=9, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PSM overlap.png", psm, width = 31, height = 9, units = "cm")

##################################################################################################################################################################

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

#	How many of the ovelapping PSMs are modified?
psm_mod_annsolo <- psm_mod %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait.x, scan, sequence.x, retention_time, scanID, mod) %>% unite(PSM_ID, bait.x:sequence.x, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait.x) %>% unique()

psm_mod_annsolo %>% select(PSM_ID) %>% unique

ann %>% filter(PSM_ID %in% psm_mod_annsolo$PSM_ID) %>% select(PSM_ID) %>% unique

#	SARS
SARS_psms <- psm_id %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14"))

SARS_mod <- psm_mod %>% 
	separate(accession, into=c("pre","name","post","twee","prot")) %>% 
	filter(name=="SARS") %>%
	mutate(prot=str_replace(prot,"Protein14", "protein14"))

SARS_mod_annsolo <- SARS_mod %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait.x, scan, sequence.x, retention_time, scanID, mod) %>% unite(PSM_ID, bait.x:sequence.x, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait.x) %>% unique()

both %>% filter(PSM_ID %in% SARS_mod_annsolo$PSM_ID) %>% unique

### filter(str_detect(Condition, prot))



#	How many of the ovelapping PSMs are modified?
psm_mod_annotations <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null") %>% unique

psm_mod_annsolo <- psm_mod_annotations %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait.x, scan, sequence.x, retention_time, scanID, mod) %>% unite(PSM_ID, bait.x:sequence.x, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait.x) %>% unique()

psm_mod_annsolo %>% select(PSM_ID) %>% unique

ann %>% filter(PSM_ID %in% psm_mod_annsolo$PSM_ID) %>% select(PSM_ID) %>% unique


both %>% filter(PSM_ID %in% psm_mod$PSM_ID) %>% select(PSM_ID) %>% unique
psm_mod_annsolo %>% filter(PSM_ID %in% both$PSM_ID) %>% select(PSM_ID, mod) %>% unique
psm_mod_annsolo %>% filter(PSM_ID %in% ann$PSM_ID) %>% select(PSM_ID, mod) %>% unique

ann <- psm_annsolo %>% filter(!PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID, scanID) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- psm_gordon %>% filter(!PSM_ID %in% psm_annsolo$PSM_ID) %>% select(RawFile, PSM_ID) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- psm_annsolo %>% filter(PSM_ID %in% psm_gordon$PSM_ID) %>% select(RawFile, PSM_ID, scanID) %>% mutate(Identification = "Overlap") %>% unique

psm_mod_both <-  psm_mod %>% filter(PSM_ID %in% both$scanID)
psm_mod_ann <-  psm_mod %>% filter(PSM_ID %in% ann$scanID)

psm_mod_ann %>% select(PSM_ID, sequence.x) %>% unique
psm_mod %>% select(PSM_ID, sequence.x) %>% unique
# count(mod) %>% arrange(desc(n))

annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

psm_mod_annotation <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%      #Filter out samples that were not included in the hcip analysis.
	filter(!accession=="null")
psm_annotations <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null")

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

sars_mod_ann <-  SARS_mod %>% filter(PSM_ID %in% ann$scanID)
sars_mod_ann %>% select(PSM_ID, sequence.x) %>% unique

sars_mod_both <-  SARS_mod %>% filter(PSM_ID %in% both$scanID)
sars_mod_both %>% select(PSM_ID, sequence.x) %>% unique

SARS_mod %>% select(PSM_ID, sequence.x) %>% unique

SARS_psm_ann <- SARS_psms %>% filter(PSM_ID %in% ann$scanID)
SARS_psm_ann %>% select(PSM_ID, sequence) %>% unique

SARS_psm_both <- SARS_psms %>% filter(PSM_ID %in% both$scanID)
SARS_psm_both %>% select(PSM_ID, sequence) %>% unique




psm_mod_both <-  psm_mod %>% filter(PSM_ID %in% both$scanID)
psm_mod_ann <-  ann %>% filter(scanID %in% psm_mod$PSM_ID)

psm_mod_ann %>% select(PSM_ID, sequence.x) %>% unique
psm_mod %>% select(PSM_ID, sequence.x) %>% unique

psm_mod %>% select(PSM_ID, sequence.x) %>% unique
psm_mod_ann %>% select(PSM_ID, sequence.x) %>% unique

psm_mod_ann %>% filter(PSM_ID %in% psm_gordon$PSM_ID)

##################################################################################################################################################################

#	Add annotation

annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")
psm_annotations <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!accession=="null") %>% unique
psm_annotations %>% select(PSM_ID) %>% unique

gordon_annotations <- merge(msms_gordon, annotation, by.x = "RawFile", by.y = "fileID") %>% 
	as_tibble %>% filter(Comments=="") %>%
	filter(!GeneNames=="")

psm_annsolo_annotation <- psm_annotations %>% separate(PSM_ID, into=c("a", "scan", "b", "c"), remove=FALSE) %>%
	rename(scanID = PSM_ID) %>%
	select(bait, scan, sequence, retention_time, scanID, accession) %>% unite(PSM_ID, bait:sequence, sep="|", remove=FALSE) %>%
	rename(Rawfile = bait) %>% unique()
psm_gordon_annotation <- evidence_g_annotations %>% select(Rawfile, MSMSscannumber, Sequence, Retentiontime, Genenames) %>%
	filter(!Sequence=="") %>%
	unite(PSM_ID, Rawfile:Sequence, sep="|", remove=FALSE) %>% unique()

ann <- psm_annsolo_annotation %>% filter(!PSM_ID %in% psm_gordon_annotation$PSM_ID) %>% select(Rawfile, PSM_ID) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- psm_gordon_annotation %>% filter(!PSM_ID %in% psm_annsolo_annotation$PSM_ID) %>% select(Rawfile, PSM_ID) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- psm_annsolo_annotation %>% filter(PSM_ID %in% psm_gordon_annotation$PSM_ID) %>% select(Rawfile, PSM_ID) %>% mutate(Identification = "Overlap") %>% unique

comb_psm <- rbind(ann, gordon, both) 

psm <- ggplot(comb_psm, aes(RawFile, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="PSM overlap", x = "Raw file", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=9, angle=90))



##################################################################################################################################################################

#	Evidence Gordon

#	Load Gordon et al. evidence
evidence_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/20200318-qx-final/evidence.txt") %>% as_tibble
#	How many identified MS/MS spectra?		359232 PSMs, 38511 are not mapped to a protein
#	359232 - 38511 = 320721
evidence_gordon %>% filter(str_detect(Proteins, "CON__")) %>% select(Genenames, Proteins) %>% unique

