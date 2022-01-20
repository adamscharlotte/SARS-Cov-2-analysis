library(tidyverse)
library(data.table)
library(UpSetR)

#	Load PPI results -----------------------------------------------------------------------------------------------------------

gordon_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/gordon_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Original analysis") %>%
	unique()
annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Reanalysis") %>%
	unique()
Li_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/Li_bp.txt", header=FALSE) %>% 
	as_tibble %>%
	rename(PPI=V1) %>%
	mutate(Study="Li et al.") %>%
	unique()
Chen_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Chen/Chen_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Chen et al.") %>%
	unique()
Stukalov_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/Stukalov_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Stukalov et al.") %>%
	unique()
Germain_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Germain/Germain_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="St-Germain et al.") %>%
	unique()
Laurent_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/Laurent_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Laurent et al.") %>%
	unique()
Samavarchi_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Samavarchi/Samavarchi_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Samavarchi-Tehrani et al.") %>%
	unique()

#
#	Create input for upsetplot --------------------------------------------------------------------------------------------------
ppis_studies <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp)
ppis_studies_annsolo <- ppis_studies %>% filter(PPI %in% annsolo_bp$PPI)

#	Upset plot ------------------------------------------------------------------------------------------------------------------
ppi_upset <- ppis_studies_annsolo %>% unique %>%
	unnest(cols = Study) %>%
	mutate(StudyMember=1) %>%
	pivot_wider(names_from = Study, values_from = StudyMember, values_fill = list(StudyMember = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset(
		order.by = "freq",
		nsets=8)

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap upsetplot.pdf", width=7, height=5)
ppi_upset
dev.off()

#	Phosphorylation sites ------------------------------------------------------------------------------------------------------------
p_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/p_sites.txt", header=TRUE)

p_studies <- p_sites %>% as_tibble %>% 
	select(-proteinAccession) %>%
	unite(P_site, viralProtein:phosphosite) %>%
	gather(Study, Member, -P_site) %>%
	filter(Member) %>%
	select(-Member) %>% unique()

p_upset <- p_studies %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset()

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/phosphorylation site overlap upsetplot.pdf", width=5, height=5)
p_upset
dev.off()

#	Ubiquitination sites ------------------------------------------------------------------------------------------------------------
ub_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/ub_sites.txt", header=TRUE)
ub_sites %>% as_tibble
ub_studies <- ub_sites %>% as_tibble %>% 
	select(-proteinAccession) %>%
	unite(UB_site, viralProtein:ubisite) %>%
	gather(Study, Member, -UB_site) %>%
	filter(Member) %>%
	select(-Member) %>% unique()

ub_upset <- ub_studies %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset()

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/ubiquitination site overlap upsetplot.pdf", width=5, height=5)
ub_upset
dev.off()