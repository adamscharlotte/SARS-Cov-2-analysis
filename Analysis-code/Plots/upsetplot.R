library(tidyverse)
library(data.table)
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
list_Liu <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/Liu/Liu_bp.txt", header=FALSE) %>% 
	as_tibble %>% 
	filter(V1 %in% list_annsolo) %>%
	pull(V1) %>% 
	unique()
# list_Liu_apms <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/Liu/Liu_bp_apms.txt", header=FALSE) %>% 
# 	as_tibble %>% 
# 	# filter(V1 %in% list_annsolo) %>%
# 	pull(V1) %>% 
# 	unique()
# list_Chen_apms <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Other studies/AP-MS/Chen/Chen_bp_apms.txt", header=FALSE) %>% 
# 	as_tibble %>% 
# 	# filter(V1 %in% list_annsolo) %>%
# 	pull(V1) %>% 
# 	unique()

# #

#	Create a binary matrix -----------------------------------------------------------------------------------------------------------
list_ppi <- list(
	"Original analysis" = list_gordon,
	"Reanalysis" = list_annsolo,
	"Li et al." = list_Li, 
	"Chen et al." = list_Chen,
	"Stukalov et al." = list_Stukalov,
	"St-Germain et al." = list_Germain,
	"Laurent et al." = list_Laurent,
	"Samavarchi-Tehrani et al." = list_Samavarchi,
	"Liu et al." = list_Liu
)
matrix_ppi <- list_to_matrix(list_ppi)

list_other <- unique(c(list_Li, list_Chen, list_Stukalov, list_Germain, list_Laurent, list_Samavarchi, list_Liu))
list_ppi <- list(
	"Original analysis" = list_gordon,
	"Reanalysis" = list_annsolo,
	"Other studies" = list_other 
)
matrix_ppi <- list_to_matrix(list_ppi)

list_ppi <- list(
	"Original analysis" = list_gordon,
	"Reanalysis" = list_annsolo,
	"Li et al." = list_Li, 
	"Chen et al. AP-MS" = list_Chen_apms,
	"Stukalov et al." = list_Stukalov,
	"Liu et al. AP-MS" = list_Liu_apms
)
matrix_ppi <- list_to_matrix(list_ppi)

#	Create a combinarion matrix ------------------------------------------------------------------------------------------------------
matrix_ppi_combination <- make_comb_mat(matrix_ppi, top_n_sets = 16) # mode = "intersect")
set_size(matrix_ppi_combination)
set_name(matrix_ppi_combination)
matrix_ppi_combination <- matrix_ppi_combination[comb_degree(matrix_ppi_combination) > 0]

pdf(file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Upset/PPI overlap final.pdf", height = 3.4, width = 10)
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
				# at = c(0, -100, -200, -300),
				# labels = c(0, 100, 200, 300),
				at = c(0, -7000, -3500),
				labels = c(0, 7000, 3500),
				labels_rot = 0),
			border = FALSE, 
			gp = gpar(fill = "black"), 
			width = unit(3, "cm")),
		set_name = anno_text(set_name(matrix_ppi_combination), 
			location = 0.5, 
			just = "center",
			width = max_text_width(set_name(matrix_ppi_combination)) + unit(-11, "mm"),
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
		default.units = "native", just = c("center", "bottom"), 
		gp = gpar(fontsize = 6, col = "#404040"))
})
dev.off()

pdf(file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Upset/PPI other studies - small.pdf", height = 3.4, width = 4)
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
			width = unit(3, "cm")),
		set_name = anno_text(set_name(matrix_ppi_combination), 
			location = 0.5, 
			just = "center",
			width = max_text_width(set_name(matrix_ppi_combination)) + unit(-1, "mm"),
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
		default.units = "native", just = c("center", "bottom"), 
		gp = gpar(fontsize = 6, col = "#404040"))
})
dev.off()

#	Load phosphorylation sites -------------------------------------------------------------------------------------------------
df_p_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/p_sites.txt", header=TRUE)
tbl_p_sites <- df_p_sites %>% as_tibble %>%
	select(-proteinAccession) %>%
	unite(P_site, viralProtein:phosphosite)

#	Create a binary matrix
colnames(p_sites)
list_bouhaddou <- tbl_p_sites %>% 
	filter(`Bouhaddou et al.` == TRUE) %>%
	pull(P_site)
list_davidson <- tbl_p_sites %>% 
	filter(`Davidson et al.` == TRUE) %>%
	pull(P_site)
list_klan <- tbl_p_sites %>% 
	filter(`Klann et al.` == TRUE) %>%
	pull(P_site)
list_reanalysis <- tbl_p_sites %>% 
	filter(Reanalysis == TRUE) %>%
	pull(P_site)
#
list_p_site <- list(
	"Reanalysis" = list_reanalysis,
	"Bouhaddou et al." = list_bouhaddou, 
	"Klann et al." = list_klan,
	"Davidson et al." = list_davidson
)
matrix_p_site <- list_to_matrix(list_p_site)

#	Create a combination matrix
matrix_p_site_combination <- make_comb_mat(matrix_p_site, top_n_sets = 4) # mode = "intersect")
set_size(matrix_p_site_combination)
set_name(matrix_p_site_combination)
matrix_p_site_combination <- matrix_p_site_combination[comb_degree(matrix_p_site_combination) > 0]

#	Plot upset
pdf(file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Upset/P-sites.pdf", height = 3.2, width = 4.5)
ss = set_size(matrix_p_site_combination)
cs = comb_size(matrix_p_site_combination)
ht = UpSet(matrix_p_site_combination, 
	set_order = order(ss),
	comb_order = order(comb_degree(matrix_p_site_combination), -cs),
	top_annotation = HeatmapAnnotation(
		"P-site intersection size" = anno_barplot(cs, 
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
				at = c(0, -10, -20, -30, -40),
				labels = c(0, 10, 20, 30, 40),
				labels_rot = 0),
			border = FALSE, 
			gp = gpar(fill = "black"), 
			width = unit(3, "cm")),
		set_name = anno_text(set_name(matrix_p_site_combination), 
			location = 0.5, 
			just = "center",
			width = max_text_width(set_name(matrix_p_site_combination)) + unit(-8, "mm"),
			gp= gpar(fontsize = 8)
			),
			annotation_name_gp= gpar(fontsize = 9)
			),
	right_annotation = NULL,
	show_row_names = FALSE)

ht = draw(ht)
od = column_order(ht)
decorate_annotation("P-site intersection size", {
	grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
		default.units = "native", just = c("center", "bottom"), 
		gp = gpar(fontsize = 6, col = "#404040"))
})
dev.off()

#	Load ubiquitination sites --------------------------------------------------------------------------------------------------
df_ub_sites <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/PTM/ub_sites.txt", header=TRUE)
tbl_ub_sites <- df_ub_sites %>% as_tibble %>%
	select(-proteinAccession) %>%
	unite(UB_site, viralProtein:ubisite)

#	Create a binary matrix
colnames(tbl_ub_sites)
list_zang <- tbl_ub_sites %>% 
	filter(`Zang et al.` == TRUE) %>%
	pull(UB_site)
list_stukalov <- tbl_ub_sites %>% 
	filter(`Stukalov et al.` == TRUE) %>%
	pull(UB_site)
list_reanalysis <- tbl_ub_sites %>% 
	filter(Reanalysis == TRUE) %>%
	pull(UB_site)
#
list_ub_site <- list(
	"Reanalysis" = list_reanalysis,
	"Zang et al." = list_zang, 
	"Stukalov et al." = list_stukalov
)
matrix_ub_site <- list_to_matrix(list_ub_site)

#	Create a combination matrix
matrix_ub_site_combination <- make_comb_mat(matrix_ub_site, top_n_sets = 3) # mode = "intersect")
set_size(matrix_ub_site_combination)
set_name(matrix_ub_site_combination)
matrix_ub_site_combination <- matrix_ub_site_combination[comb_degree(matrix_ub_site_combination) > 0]

#	Plot upset
pdf(file = "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Upset/Ub-sites.pdf", height = 3, width = 3.5)
ss = set_size(matrix_ub_site_combination)
cs = comb_size(matrix_ub_site_combination)
ht = UpSet(matrix_ub_site_combination, 
	set_order = order(ss),
	comb_order = order(comb_degree(matrix_ub_site_combination), -cs),
	top_annotation = HeatmapAnnotation(
		"Intersection size" = anno_barplot(cs, 
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
				at = c(0, -50, -100),
				labels = c(0, 50, 100),
				labels_rot = 0),
			border = FALSE, 
			gp = gpar(fill = "black"), 
			width = unit(3, "cm")),
		set_name = anno_text(set_name(matrix_ub_site_combination), 
			location = 0.5, 
			just = "center",
			width = max_text_width(set_name(matrix_ub_site_combination)) + unit(-8, "mm"),
			gp= gpar(fontsize = 8)
			),
			annotation_name_gp= gpar(fontsize = 9)
			),
	right_annotation = NULL,
	show_row_names = FALSE)

ht = draw(ht)
od = column_order(ht)
decorate_annotation("Intersection size", {
	grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
		default.units = "native", just = c("center", "bottom"), 
		gp = gpar(fontsize = 6, col = "#404040"))
})
dev.off()

