# devtools::install_github("const-ae/ggupset")

library(ggplot2)
library(tidyverse)
library(data.table)
library(ggupset)
library(UpSetR)

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
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Gordon")
annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="ANN")
Li_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/Li_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Li")
Chen_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Chen/Chen_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Chen")
Stukalov_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/Stukalov_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Stukalov")
Germain_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Germain/Germain_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Germain")
Laurent_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/Laurent_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Laurent")
Samavarchi_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Samavarchi/Samavarchi_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Samavarchi")

ppis_studies <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp) %>%
	group_by(PPI) %>%
	summarize(Studies = list(Study))
ppis_studies_annsolo <- ppis_studies %>% filter(PPI %in% annsolo_bp$PPI)

ppi_upset <- ggplot(ppis_studies_annsolo, aes(x = Studies), show.numbers="yes") +
	geom_bar() +
	scale_x_upset()

ppis <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp) 

ppis_studies_annsolo %>% unnest(cols = Studies) %>% filter(Studies=="Chen")

ppi_plotter <- ppis_studies_annsolo %>%
	unnest(cols = Studies) %>%
	# mutate(StudyMember=as.integer(1)) %>%
	mutate(StudyMember=1) %>%
	pivot_wider(names_from = Studies, values_from = StudyMember, values_fn = list(count=list), 
	# values_fill = list(as.integer(0))) %>%
	values_fill = list(0)) %>%
	as.data.frame() #%>% head()

ppi_int <- ppis_studies_annsolo %>%
	unnest(cols = Studies) %>%
	mutate(StudyMember=as.integer(1)) %>%
	pivot_wider(names_from = Studies, values_from = StudyMember, values_fn = list(count=list), 
	values_fill = list(as.integer(0))) %>%
	mutate(Gordon=unlist(Gordon), ANN=unlist(ANN), Stukalov=unlist(Stukalov), Laurent=unlist(Laurent),
		Samavarchi=unlist(Samavarchi), Li=unlist(Li), Germain=unlist(Germain)) %>%
	mutate(Chen=unlist(Chen)) 
ppi_int %>% mutate(Chen=unlist(Chen)) 


for(i in 2:ncol(ppi_plotter)){ ppi_plotter[ , i] <- as.numeric(ppi_plotter[ , i]) }
for(i in 2:ncol(ppi_plotter)){ ppi_plotter[ , i] <- as.integer(ppi_plotter[ , i]) 
	for(j in 1:nrow(ppi_plotter)){ ppi_plotter[ j, i] <- ppi_plotter[j, i]-1 } }

head(ppi_plotter)

upset(ppi_plotter, sets=c(""))

?upset()


ppis_studies_annsolo %>% upset()

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap upsetplot.png", ppi_upset, width = 31, height = 7, units = "cm")

####################################################################################################################################################################################

#	Load PPI results

gordon_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/gordon_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Original analysis")
annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% rename(PPI=V1) %>% mutate(Study="Reanalysis")
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

#	Get set of studies with an overlap <5

ppis_studies_list <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp) %>%
	group_by(PPI) %>%
	summarize(Studies = list(Study))
ppis_studies_annsolo_list <- ppis_studies_list %>% filter(PPI %in% annsolo_bp$PPI)
Other_leq5 <- ppis_studies_annsolo_list %>% add_count(Studies) %>% filter(n<=5)

#	Create input for upsetplot
ppis_studies <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp)
ppis_studies_annsolo <- ppis_studies %>% filter(PPI %in% annsolo_bp$PPI)

ppis_studies_annsolo %>% add_count(PPI) %>% filter(n<5)
ppis_studies_annsolo_high <- ppis_studies_annsolo %>% filter(!PPI %in% Other_leq5$PPI)
ppis_studies_annsolo_low <- ppis_studies_annsolo %>% filter(PPI %in% Other_leq5$PPI)

#	Upset plot

ppi_int <- ppis_studies_annsolo_high %>% unique %>%
	unnest(cols = Study) %>%
	mutate(StudyMember=1) %>%
	pivot_wider(names_from = Study, values_from = StudyMember, values_fill = list(StudyMember = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset(order.by = "freq", cutoff = 10)

ppi_upset <- ppis_studies_annsolo %>% unique %>%
	unnest(cols = Study) %>%
	mutate(StudyMember=1) %>%
	pivot_wider(names_from = Study, values_from = StudyMember, values_fill = list(StudyMember = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset(order.by = "freq", nsets=8)

ppis_studies %>% select(Study) %>% unique
ppis_studies_annsolo_high %>% select(Study) %>% unique

ppi_upset <- ppis_studies_annsolo_low %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member) %>% unlist
	, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset(nintersects = 70, group.by = "sets", cutoff = 7)

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap upsetplot.pdf")
ppi_upset
dev.off()

####################################################################################################################################################################################

#	Load phosphorylation sites
p_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/p_sites.txt", header=TRUE)

p_studies <- p_sites %>% as_tibble %>% 
	select(-proteinAccession) %>%
	unite(P_site, viralProtein:phosphosite) %>%
	gather(Study, Member, -P_site) %>%
	filter(Member) %>%
	select(-Member) %>% unique()

p_upset <- p_studies %>% group_by(P_site) %>%
	summarize(Studies = list(Study)) %>%
	ggplot(aes(x = Studies)) +
	geom_bar() +
	scale_x_upset()

p_studies%>% select(Study) %>% unique

p_upset <- p_studies %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset()

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/phosphorylation site overlap upsetplot.pdf")
p_upset
dev.off()

####################################################################################################################################################################################

#	Load ubiquitination sites
ub_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/ub_sites.txt", header=TRUE)

ub_studies <- ub_sites %>% as_tibble %>% 
	select(-proteinAccession) %>%
	unite(UB_site, viralProtein:ubisite) %>%
	gather(Study, Member, -UB_site) %>%
	filter(Member) %>%
	select(-Member) %>% unique()

ub_upset <- ub_studies %>% group_by(UB_site) %>%
	summarize(Studies = list(Study)) %>%
	ggplot(aes(x = Studies)) +
	geom_bar() +
	scale_x_upset()

ub_studies %>% select(Study) %>% unique
ub_upset <- ub_studies %>% mutate(Member=1) %>%
	pivot_wider(names_from = Study, values_from = Member, values_fill = list(Member = 0)) %>%
	as.data.frame() %>%
	UpSetR::upset()

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/ubiquitination site overlap upsetplot.pdf")
ub_upset
dev.off()