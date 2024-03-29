library(data.table)
library(tidyverse)
library(annotables)			# Mapping
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(rrvgo)
library(ComplexHeatmap)

# Load data ------------------------------------------------------------------------------------------------------------------------

#	Load protein-protein interactions
tbl_ppi <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv") %>% as_tibble

#	Load genes background
path_fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
tbl_fasta_headers <- fread(path_fasta_headers, header = FALSE) %>% as_tibble
tbl_background <- tbl_fasta_headers %>% 
	separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% 
	drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	dplyr::rename(gene_symbol = gene_name) %>%
	dplyr::select(gene_symbol) %>% 
	unique()

# Map entrez gene names ------------------------------------------------------------------------------------------------------------
tbl_gene_map <- grch38 %>% 
	mutate(entrez_gene = as.character(entrez)) %>% 
	dplyr::rename(gene_symbol = symbol) %>% 
	dplyr::select(gene_symbol, entrez_gene) %>% 
	unique()

tbl_bait_prey <- tbl_ppi %>% 
	dplyr::rename(gene_symbol = PreyGene) %>% 
	dplyr::select(idBait, gene_symbol)

tbl_bait_entrez <- merge(tbl_bait_prey, tbl_gene_map) %>% as_tibble %>% 
	dplyr::select(idBait, entrez_gene) %>% 
	unique()

tbl_background_entrez <- merge(tbl_background, tbl_gene_map) %>% as_tibble %>% 
	dplyr::select(gene_symbol, entrez_gene) %>% 
	unique()

# GO enrichment analysis for each bait --------------------------------------------------------------------------------------------
Baits <- tbl_bait_entrez %>% dplyr::pull(idBait) %>% unique()
tbl_go_result <- data.frame()
for(bait in Baits) {
	tbl_selection <- tbl_bait_entrez %>% 
		filter(idBait == bait) %>% 
		unique()
	df_enrichgo_result <- enrichGO(
		gene = tbl_selection$entrez_gene, 
		org.Hs.eg.db, 
		ont="BP",
		universe = tbl_background_entrez$entrez_gene, 
		pAdjustMethod="BH",
		pvalueCutoff=0.01
		)
	tbl_enrichgo_result <- df_enrichgo_result %>% as_tibble %>% mutate(ViralProtein = bait)
	tbl_go_result <- rbind(tbl_go_result, tbl_enrichgo_result)
}

tbl_go_result_05 <- tbl_go_result
tbl_go_result <- tbl_go_result_01

# Simplify the GO enrichment results --------------------------- (rrvgo) ----------------------------------------------------------

#	Reduce GO terms 
matrix_similarity <- calculateSimMatrix(
	tbl_go_result$ID,
	orgdb="org.Hs.eg.db",
	ont="BP",
	method="Rel"
	)
scores <- setNames(-log10(tbl_go_result$p.adjust), tbl_go_result$ID)
df_reduced_terms <- reduceSimMatrix(
	matrix_similarity,
	scores,
	threshold=0.7,
	orgdb="org.Hs.eg.db"
	)
plot_scatter <- scatterPlot(matrix_similarity, df_reduced_terms)
df_reduced_terms %>% as_tibble %>% dplyr::select(term, parentTerm, score, size) %>% 
	arrange(desc(parentTerm)) %>%
	print(n=80) 

#	Reduce the GO terms for each viral protein ---------------------------------------------------------------------------------------
tbl_reduced_terms <- data.frame() %>% as_tibble
Baits <- tbl_go_result %>% filter(!ViralProtein == "orf7a") %>% dplyr::pull(ViralProtein) %>% unique()
for (bait in Baits) {
	tbl_go_bait <- tbl_go_result %>% filter(ViralProtein==bait)
	matrix_similarity_bait <- calculateSimMatrix(
		tbl_go_bait$ID,
		orgdb="org.Hs.eg.db",
		ont="BP",
		method="Rel"
		)
	scores <- setNames(-log10(tbl_go_bait$p.adjust), tbl_go_bait$ID)
	df_reduced_terms_bait <- reduceSimMatrix(
		matrix_similarity_bait,
		scores,
		threshold=0.7,
		orgdb="org.Hs.eg.db"
		)
	tbl_reduced_terms_bait <- df_reduced_terms_bait %>% as_tibble %>% 
		mutate(ID = go) %>% 
		mutate(parentID = parent) %>%
		dplyr::select(ID, parentID, parentTerm, score, size) %>%
		mutate(ViralProtein = bait)
	tbl_reduced_terms <- rbind(tbl_reduced_terms, tbl_reduced_terms_bait)
}

tbl_reduced_terms %>% dplyr::pull(parentTerm) %>% unique()

tbl_reduced_terms_orf7a <- data.frame(
	ID = "GO:0042273", 
	parentID = "GO:0042273",
	parentTerm = "ribosomal large subunit biogenesis",
	score = NA,
	size = NA,
	ViralProtein = "orf7a"
	) %>% as_tibble

tbl_reduced_go_result <- rbind(tbl_reduced_terms, tbl_reduced_terms_orf7a) %>%
	merge(tbl_go_result) %>% 
	as_tibble %>%
	group_by(parentTerm) %>% 
	top_n(2, score) %>%
	mutate(score_max = max(score)) %>%
	mutate(score_diff = score_max - score) %>%
	filter(score_diff == 0) %>%
	top_n(1, size) %>%
	ungroup() 

rbind(tbl_reduced_terms, tbl_reduced_terms_orf7a) %>%
	merge(tbl_go_result) %>% 
	as_tibble %>%
	# filter(Description == "rRNA processing") %>%
	# filter(parentTerm == "ncRNA metabolic process") %>%
	filter(parentTerm == "tRNA-containing ribonucleoprotein complex export from nucleus") %>%
	# filter(parentTerm == "granulocyte activation") %>%
	# filter(parentTerm == "oxidative phosphorylation") %>%
	# filter(parentTerm == "tRNA transport") %>%
	dplyr::select(score, size, Description) %>%
	unique()

# Get the GO term order --------------------------------------------------------------------------------------------------------------
#	Create similarity matrix
matrix_go_similarity <- calculateSimMatrix(
	tbl_reduced_go_result$ID,
	orgdb="org.Hs.eg.db",
	ont="BP",
	method="Rel"
	)

#	Cluster the GO terms based on similarity
plot_go_terms <- Heatmap(
	matrix_go_similarity,
	clustering_method_rows = "ward.D",
	clustering_method_columns = "ward.D"
	)

#	Get the GO terms from the rownames
df_go_terms <- matrix_go_similarity[, 0] %>% as.data.frame
go_terms <- rownames(df_go_terms)
rownames(df_go_terms) <- NULL
tbl_go_terms <- cbind(go_terms, df_go_terms) %>% 
	as_tibble %>% 
	mutate(index = 1:n())

#	Get the order of the GO terms
ht_list = draw(plot_go_terms)
index_order <- row_order(ht_list)
tbl_index_order <- data.frame(index_order = index_order) %>% 
	as_tibble %>%
	mutate(index = as.integer(tbl_go_terms$index))

#	Map the GO terms to the order
tbl_go_order <- merge(tbl_index_order, tbl_go_terms) %>% 
	as_tibble %>% rename(go_terms = "ID")


# Plot the GO terms -----------------------------------------------------------------------------------------------------------------
tbl_ggplot <- data.frame()
Baits <- tbl_reduced_go_result %>% dplyr::pull(ViralProtein) %>% unique()
tbl_all_descriptions <- tbl_reduced_go_result %>% dplyr::select(Description, ID) %>% unique
for (bait in Baits) {
	tbl_simple_go_bait <- tbl_reduced_go_result %>% filter(ViralProtein==bait) %>% dplyr::select(-ViralProtein)
	tbl_description_bait <- tbl_all_descriptions %>% left_join(tbl_simple_go_bait) %>% mutate(ViralProtein = bait)
	tbl_ggplot <- rbind(tbl_ggplot, tbl_description_bait)
}

tbl_ggplot_input <- tbl_ggplot %>% 
	merge(tbl_go_order) %>%
	as_tibble %>%
	arrange(index_order) %>%
	mutate(logpadjust = -log10(p.adjust)) %>% 
	dplyr::select(Description, logpadjust, ViralProtein, ID, index_order) %>% 
	unique() %>% 
	replace(is.na(.), 0)

first_to_upper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tbl_ggplot_upper <- tbl_ggplot_input %>% 
	mutate(ViralProtein = first_to_upper(ViralProtein)) %>%
	mutate(Description = first_to_upper(Description)) %>%
	mutate(Description = str_replace(Description, "TRNA", "tRNA")) %>%
	mutate(Description = str_replace(Description, "NcRNA", "ncRNA"))

tbl_ggplot_input %>% pull(ViralProtein) %>% unique()
tbl_ggplot_upper$ViralProtein <- factor(tbl_ggplot_upper$ViralProtein,
	levels = c("M", "N", "S", "Nsp1", "Nsp4", "Nsp6", "Nsp7", "Nsp8", "Nsp10", "Nsp11",
	"Nsp13", "Nsp14", "Orf3a", "Orf3b", "Orf6", "Orf7a", "Orf8", "Orf9c", "Orf10"))

plot <- ggplot(tbl_ggplot_upper, aes(y=Description, x=ViralProtein, fill=logpadjust)) + 
	geom_tile(colour = "black", size = 0.40) +
	scale_fill_viridis_c() +
	labs(title = "", x = "", y = "", fill=expression("-log"["10"]*"(adjusted p-value)")) +
	scale_fill_gradientn(
		colors = colorRampPalette(brewer.pal(8, "Blues"))(25),
		limits = c(0,10),
		na.value = "#08306C"
		) +
	theme(
		panel.background = element_blank(), 
		axis.ticks= element_blank(),
		# legend.key.size = unit(1, 'lines'),
		legend.key.size = unit(0.4, "cm"),
		# legend.spacing = unit(5, 'cm'),
		legend.title = element_text(size = 7, face="bold", colour="black"), 
		legend.text = element_text(size = 7),
		axis.text.x = element_text(size = 8, vjust=0.7, angle=90, colour="black"),
		axis.text.y = element_text(size = 8, colour="black"),
		panel.border = element_rect(colour = "black", fill=NA, size=1.2)
	) 

ggsave(
	plot, 
	file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Heatmap/GO analysis clustered terms.png", 
	width = 18.5, 
	height = 25, 
	units = "cm"
	)

# Select the top GO terms ----------------------------------------------------------------------------------------------------
tbl_top_go_result <- tbl_reduced_go_result %>%
	merge(tbl_go_order) %>%
	arrange(index_order) %>%
	group_by(ViralProtein) %>% 
	top_n(-1, p.adjust) %>%
	ungroup()

# Heatmap for the top GO terms -----------------------------------------------------------------------------------------------
tbl_top_ggplot <- data.frame()
Baits <- tbl_top_go_result %>% dplyr::pull(ViralProtein) %>% unique()
tbl_all_descriptions <- tbl_top_go_result %>% dplyr::select(Description, ID, index_order) %>% unique()
for (bait in Baits) {
	tbl_simple_go_bait <- tbl_top_go_result %>% 
		filter(ViralProtein == bait) %>% 
		dplyr::select(-ViralProtein)
	tbl_description_bait <- tbl_all_descriptions %>% 
		left_join(tbl_simple_go_bait) %>% 
		mutate(ViralProtein = bait)
	tbl_top_ggplot <- rbind(tbl_top_ggplot, tbl_description_bait)
}

tbl_top_ggplot_input <- tbl_top_ggplot %>% 
	mutate(logpadjust = -log10(p.adjust)) %>% 
	dplyr::select(Description, logpadjust, ViralProtein, ID, index_order) %>% 
	unique() %>% 
	replace(is.na(.), 0)

tbl_top_ggplot_upper <- tbl_top_ggplot_input %>% 
	mutate(ViralProtein = first_to_upper(ViralProtein)) %>%
	mutate(Description = first_to_upper(Description)) %>%
	mutate(Description = str_replace(Description, "TRNA", "tRNA")) %>%
	mutate(Description = str_replace(Description, "NcRNA", "ncRNA"))

tbl_top_ggplot_upper$ViralProtein <- factor(tbl_top_ggplot_upper$ViralProtein,
	levels = c("M", "Nsp1", "Nsp4", "Nsp6", "Nsp7", "Nsp8", "Nsp11",
	"Nsp13", "Nsp14", "Orf3b", "Orf6", "Orf7a", "Orf8", "Orf9c", "Orf10"))

plot <- ggplot(tbl_top_ggplot_upper, aes(y=Description, x=ViralProtein, fill=logpadjust)) + 
	geom_tile(colour = "black", size = 0.40) +
	labs(title = "", x = "", y = "", fill=expression("-log"["10"]*"(adjusted p-value)")) +
	scale_fill_gradientn(
		colors = colorRampPalette(brewer.pal(8, "Blues"))(25),
		limits = c(0,10),
		na.value = "#08306C"
		) +
	theme(
		panel.background = element_blank(), 
		axis.ticks= element_blank(),
		legend.key.size = unit(0.4, "cm"),
		legend.title = element_text(size = 7, face="bold", colour="black"), 
		legend.text = element_text(size = 7),
		axis.text.x = element_text(size = 8, vjust=0.7, angle=90, colour="black"),
		axis.text.y = element_text(size = 8, colour="black"),
		panel.border = element_rect(colour = "black", fill=NA, size=1.2)
	) 

ggsave(
	plot, 
	file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Heatmap/GO analysis top terms.png", 
	width = 18.5, 
	height = 9, 
	units = "cm"
	)

###################################################################################################################################################


heatmap(matrix_ggplot_input)

plot_heatmap <- heatmap(tbl_ggplot_input, Colv = NA, Rowv = NA, scale="logpadjust", ylab="Description", xlab="ViralProtein") + 
	# geom_raster() +
	scale_fill_viridis_c() +
	labs(title="", x ="", y = "") +
	scale_fill_gradientn(colors = coul) +
	# scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(25)) +
	theme(panel.background = element_blank(), axis.ticks= element_blank())

###################################################################################################################################################
#	Top parent terms
tbl_top_go_result <- merge(tbl_go_result, tbl_reduced_terms) %>% as_tibble %>%
	filter(parentTerm == Description) %>%
	group_by(ViralProtein) %>% top_n(-1, pvalue) %>%
	ungroup()

tbl_top_ggplot <- data.frame()
Baits <- tbl_top_go_result %>% dplyr::pull(ViralProtein) %>% unique()
tbl_all_descriptions <- tbl_top_go_result %>% dplyr::select(parentTerm) %>% unique
for(bait in Baits) {
	tbl_simple_go_bait <- tbl_reduced_go_result %>% filter(ViralProtein==bait) %>% dplyr::select(-ViralProtein)
	tbl_description_bait <- tbl_all_descriptions %>% left_join(tbl_simple_go_bait) %>% mutate(ViralProtein = bait)
	tbl_top_ggplot <- rbind(tbl_top_ggplot, tbl_description_bait)
}

tbl_top_ggplot_input <- tbl_top_ggplot %>% mutate(logpadjust = -log10(p.adjust)) %>% 
	dplyr::select(parentTerm, logpadjust, ViralProtein) %>% unique %>% replace(is.na(.), 0) 

tbl_top_ggplot_input$ViralProtein <- factor(tbl_top_ggplot_input$ViralProtein,levels = c("M", "nsp1", "nsp4", "nsp6", "nsp7", "nsp8", 
	"nsp11", "nsp13", "nsp14", "orf3b", "orf6", "orf7a", "orf8", "orf9c", "orf10"))

plot <- ggplot(tbl_top_ggplot_input, aes(y=parentTerm, x=ViralProtein, fill=logpadjust)) + 
	geom_tile(colour="grey",size=0.40) +
	scale_fill_viridis_c() +
	labs(title="", x ="", y = "") +
	# scale_fill_gradientn(colors = my_colors, na.value="white") +
	scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "Blues"))(25)) +
	# scale_fill_gradient(low = "white", high = "darkblue")
	theme(panel.background = element_blank(), axis.ticks= element_blank(), 
	axis.text.x = element_text(vjust=0.7, angle=90)) +
	labs(fill='-log10(adjusted p-value)')

plot <- plot + theme(legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 7, face="bold"), 
               legend.text = element_text(size = 7))
ggsave(plot, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Heatmap/GO analysis.png", width = 22.5, height = 11, units = "cm")

###################################################################################################################################################

library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="MF")
###################################################################################################################################################

#	Simplify the GO enrichment results
Baits <- tbl_go_result %>% dplyr::pull(ViralProtein) %>% unique()
tbl_simple_go_result <- data.frame()
for(bait in Baits) {
	tbl_selection <- tbl_bait_entrez %>% filter(idBait == bait) %>% unique()
	df_enrichgo_result <- enrichGO(gene = tbl_selection$entrez_gene, org.Hs.eg.db, ont="BP",
		universe = tbl_background_entrez$entrez_gene, pAdjustMethod="BH", pvalueCutoff=0.01)
	tbl_simple_enrichgo_result <- simplify(df_enrichgo_result, cutoff=0.1, by="p.adjust", select_fun=min) %>% as_tibble %>% 
		dplyr::select(Description, qvalue, p.adjust, pvalue) %>% mutate(ViralProtein = bait)
	tbl_simple_go_result <- rbind(tbl_simple_go_result, tbl_simple_enrichgo_result)
}

tbl_simple_go_result %>% dplyr::select(Description) %>% unique() 	# 205 at 0.8 cutoff less at 0.6 -> 309 without simplify
																	# 42 at 0.1 cutoff

###################################################################################################################################################

HM_matrix <- data.frame() %>% mutate(Description = "") %>% as_tibble
ALL_descriptions <- data.frame() %>% mutate(Description = "") %>% as_tibble

for(bait in Baits) {
	selection <- PPIs %>% filter(idBait == bait) %>% unique()
	enrichgo_result <- enrichGO(gene = selection$entrez_gene, org.Hs.eg.db, ont="BP",
		universe = background_genes$entrez_gene, pvalueCutoff=0.01)
	simple_enrichgo_result <- simplify(enrichgo_result) %>% as_tibble %>% 
		mutate(logQ = -log10(qvalue)) %>% dplyr::select(Description, logQ)
	simple_enrichgo_result[, bait] <- simple_enrichgo_result$logQ
	simple_enrichgo_result <- simple_enrichgo_result %>% dplyr::select(-c(logQ))
	Descriptions <- simple_enrichgo_result %>% dplyr::select(Description)
	ALL_descriptions <- rbind(ALL_descriptions,Descriptions) %>% as_tibble %>% unique()
	HM_matrix <- ALL_descriptions %>% left_join(HM_matrix) %>% left_join(simple_enrichgo_result)
}

HM_matrix_col <- HM_matrix %>% dplyr::select(-Description) %>% colnames()
HM_matrix_input <- HM_matrix %>% replace(is.na(.), 0) %>% mutate_at(HM_matrix_col, as.numeric) %>% 
	dplyr::select(where(~ is.numeric(.x) && sum(.x) != 0)) %>% as.matrix()
rownames(HM_matrix_input) = ALL_descriptions$Description

png("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/test.png",width=67,height=54,units="cm",res=1200)
Heatmap(HM_matrix_input,
	name = 'GO Enrichment', 
	show_row_dend = FALSE, show_column_dend = FALSE, 
	col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
	heatmap_legend_param = list(title="-log10(Q)")
)
dev.off()


