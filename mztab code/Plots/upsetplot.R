library(tidyverse)
library(data.table)
library(UpSetR)

#	Load PPI results -----------------------------------------------------------------------------------------------------------
gordon_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Gordon/gordon_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Original analysis") %>%
	unique()
annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Reanalysis") %>%
	unique()
Li_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Li/Li_bp.txt", header=FALSE) %>% 
	as_tibble %>%
	rename(PPI=V1) %>%
	mutate(Study="Li et al.") %>%
	unique()
Chen_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Chen/Chen_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Chen et al.") %>%
	unique()
Stukalov_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Stukalov/Stukalov_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Stukalov et al.") %>%
	unique()
Germain_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Germain/Germain_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="St-Germain et al.") %>%
	unique()
Laurent_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Laurent/Laurent_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Laurent et al.") %>%
	unique()
Samavarchi_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Samavarchi/Samavarchi_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	rename(PPI=V1) %>% 
	mutate(Study="Samavarchi-Tehrani et al.") %>%
	unique()

#
#	Create input for upsetplot --------------------------------------------------------------------------------------------------
ppis_studies <- rbind(gordon_bp, annsolo_bp, Li_bp, Chen_bp, Stukalov_bp, Germain_bp, Laurent_bp, Samavarchi_bp)
ppis_studies_annsolo <- ppis_studies %>% filter(PPI %in% annsolo_bp$PPI)

#	Upset plot ------------------------------------------------------------------------------------------------------------------
tbl_ppi_upset <- ppis_studies_annsolo %>% unique %>%
	unnest(cols = Study) %>%
	mutate(StudyMember=1) %>%
	pivot_wider(names_from = Study, values_from = StudyMember, values_fill = list(StudyMember = 0)) 
	
plot_ppi_upset <- tbl_ppi_upset %>%
	as.data.frame() %>%
	UpSetR::upset(
		order.by = "freq",
		nsets=8)

pdf(file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Overlap/PPI overlap upsetplot.pdf", width=7, height=5)
plot_ppi_upset
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

########################################################################################################################################

## UPSETPLOT WITH COMPLEX HEATMAP

library(ComplexHeatmap)

#	Load PPI results -----------------------------------------------------------------------------------------------------------------
list_annsolo <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/ANN-SoLo/annsolo_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	pull(V1) %>% 
	unique()
list_gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Gordon/gordon_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Li <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Li/Li_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Chen <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Chen/Chen_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Stukalov <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Stukalov/Stukalov_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Germain <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Germain/Germain_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Laurent <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Laurent/Laurent_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
list_Samavarchi <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/BioID/Samavarchi/Samavarchi_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()

#
#	Create a binary matrix -----------------------------------------------------------------------------------------------------------
list_ppi <- list(
	"Original analysis" = list_gordon,
	"Reanalysis" = list_annsolo,
	"Li et al." = list_Li, 
	"Chen et al." = list_Chen,
	"Stukalov et al." = list_Stukalov,
	"St-Germain et al." = list_Germain,
	"Laurent et al." = list_Laurent,
	"Samavarchi-Tehrani et al." = list_Samavarchi
)
matrix_ppi <- list_to_matrix(list_ppi)

#	Create a combinarion matrix ------------------------------------------------------------------------------------------------------
matrix_ppi_combination <- make_comb_mat(matrix_ppi, top_n_sets = 8) # mode = "intersect")
set_size(matrix_ppi_combination)
set_name(matrix_ppi_combination)
matrix_ppi_combination <- matrix_ppi_combination[comb_degree(matrix_ppi_combination) > 0]

pdf(file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Upset/PPI other studies.pdf", height = 3.2, width = 10)
ss = set_size(matrix_ppi_combination)
cs = comb_size(matrix_ppi_combination)
ht = UpSet(matrix_ppi_combination, 
    set_order = order(ss),
    comb_order = order(comb_degree(matrix_ppi_combination), -cs),
    top_annotation = HeatmapAnnotation(
        "PPI intersection size" = anno_barplot(cs, 
            ylim = c(0, max(cs)*1.1),
            border = FALSE, 
            gp = gpar(fill = "black"), 
            height = unit(4, "cm")
        ), 
		annotation_name_gp= gpar(fontsize = 9),
        annotation_name_side = "left", 
        annotation_name_rot = 90),
    left_annotation = rowAnnotation(
        "Set size" = anno_barplot(-ss, 
            baseline = 0,
            axis_param = list(
                at = c(0, -100, -200, -300),
                labels = c(0, 100, 200, 300),
                labels_rot = 0),
            border = FALSE, 
            gp = gpar(fill = "black"), 
            width = unit(4, "cm")),
		set_name = anno_text(set_name(matrix_ppi_combination), 
            location = 0.5, 
            just = "center",
            width = max_text_width(set_name(matrix_ppi_combination)) + unit(0.01, "mm"),
			gp= gpar(fontsize = 8)
			),
			annotation_name_gp= gpar(fontsize = 9)
			),
    right_annotation = NULL,
    show_row_names = FALSE)

ht = draw(ht)
od = column_order(ht)
decorate_annotation("PPI intersection size", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
        default.units = "native", just = c("left", "bottom"), 
        gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})
dev.off()