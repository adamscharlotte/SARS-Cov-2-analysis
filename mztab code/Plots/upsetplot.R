# devtools::install_github("const-ae/ggupset")

library(ggplot2)
library(tidyverse)
library(data.table)
library(ggupset)

tidy_movies %>% pull(Genres)
tidy_movies %>%
	distinct(title, year, length, .keep_all=TRUE) %>%
	ggplot(aes(x=Genres)) +
		geom_bar() +
		scale_x_upset(n_intersections = 20)

data("gene_pathway_membership")
tidy_pathway_member <- gene_pathway_membership %>%
	as_tibble(rownames = "Pathway") %>%
	gather(Gene, Member, -Pathway) %>%
	filter(Member) %>%
	select(- Member)

tidy_pathway_member %>%
	group_by(Gene) %>%
	summarize(Pathways = list(Pathway)) %>%
	ggplot(aes(x = Pathways)) +
		geom_bar() +
		scale_x_upset()

#	Load BP lists
gordon_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/gordon_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Gordon et al.")
annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="ANN-SoLo")
Li_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/Li_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Li et al.")
Chen_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Chen/Chen_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Chen et al.")
Stukalov_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/Stukalov_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Stukalov et al.")
Germain_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Germain/Germain_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="St-Germain et al.")
Laurent_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/Laurent_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Laurent et al.")
Samavarchi_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Samavarchi/Samavarchi_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Samavarchi-Tehrani et al.")

ppis_studies <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp) %>%
	group_by(PPI) %>%
	summarize(Studies = list(Study))
ppis_studies_annsolo <- ppi_studies %>% filter(PPI %in% annsolo_bp$PPI)

ppi_upset <- ggplot(ppis_studies_annsolo, aes(x = Studies)) +
	geom_bar() +
	scale_x_upset()

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap upsetplot.png", ppi_upset, width = 31, height = 7, units = "cm")

#	Load phosphorylation sites
p_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/p_sites.txt", header=TRUE)

p_studies <- p_sites %>% as_tibble %>% 
	select(-proteinAccession) %>%
	unite(P_site, viralProtein:phosphosite) %>%
	gather(Study, Member, -P_site) %>%
	filter(Member) %>%
	select(-Member) %>% unique()

p_upset <- p_studies %>% roup_by(P_site) %>%
	summarize(Studies = list(Study)) %>%
	ggplot(aes(x = Studies)) +
	geom_bar() +
	scale_x_upset()

p_upset <- p_studies %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset()

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/phosphorylation site overlap upsetplot.pdf")
p_upset
dev.off()