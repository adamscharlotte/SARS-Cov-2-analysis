library(data.table)
library(tidyverse)
# library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
# library(msigdbr)
# library(ComplexHeatmap)
library(RColorBrewer)
# library(GOxploreR)
library(rrvgo)

# Load data ---------------------------

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

#	Map entrez gene names
tbl_gene_map <- grch38 %>% 
	mutate(entrez_gene=as.character(entrez)) %>% 
	dplyr::rename(gene_symbol=symbol) %>% 
	dplyr::select(gene_symbol, entrez_gene) %>% 
	unique()
tbl_bait_prey <- tbl_ppi %>% dplyr::rename(gene_symbol = PreyGene) %>% dplyr::select(idBait, gene_symbol)
tbl_bait_entrez <- merge(tbl_bait_prey, tbl_gene_map) %>% as_tibble %>% dplyr::select(idBait, entrez_gene) %>% unique 
tbl_background_entrez <- merge(tbl_background, tbl_gene_map) %>% as_tibble %>% dplyr::select(gene_symbol, entrez_gene) %>% unique

#	Perform GO enrichment analysis for each bait
Baits <- tbl_bait_entrez %>% dplyr::pull(idBait) %>% unique()
tbl_go_result <- data.frame()
for(bait in Baits) {
	tbl_selection <- tbl_bait_entrez %>% filter(idBait == bait) %>% unique()
	df_enrichgo_result <- enrichGO(gene = tbl_selection$entrez_gene, org.Hs.eg.db, ont="BP",
		universe = tbl_background_entrez$entrez_gene, pAdjustMethod="BH", pvalueCutoff=0.01)
	tbl_enrichgo_result <- df_enrichgo_result %>% as_tibble %>% mutate(ViralProtein = bait)
	tbl_go_result <- rbind(tbl_go_result, tbl_enrichgo_result)
}

#	Simplify the GO enrichment results
matrix_similarity <- calculateSimMatrix(tbl_go_result$ID,
								orgdb="org.Hs.eg.db",
								ont="BP",
								method="Rel")
scores <- setNames(-log10(tbl_go_result$p.adjust), tbl_go_result$ID)
df_reduced_terms <- reduceSimMatrix(matrix_similarity,
								scores,
								threshold=0.7,
								orgdb="org.Hs.eg.db")
scatterPlot(matrix_similarity, df_reduced_terms)
treemapPlot(df_reduced_terms)
heatmapPlot(matrix_similarity,
            df_reduced_terms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

tbl_reduced_terms <- df_reduced_terms %>% as_tibble %>% mutate(ID = go) %>% dplyr::select(ID, parentTerm)

tbl_reduced_go_result <- merge(tbl_go_result, tbl_reduced_terms) %>% as_tibble %>%
	filter(parentTerm == Description)
	# dplyr::select(parentTerm, ViralProtein, everything()) %>% 
	# unite(VP, parentTerm:ViralProtein, sep ="_", remove=FALSE) %>% 
	# group_by(VP) %>% mutate(avg_value = mean(p.adjust)) %>%
	# ungroup() %>% dplyr::select(-VP)

tbl_ggplot <- data.frame()
Baits <- tbl_reduced_go_result %>% dplyr::pull(ViralProtein) %>% unique()
tbl_all_descriptions <- tbl_reduced_go_result %>% dplyr::select(parentTerm) %>% unique
for (bait in Baits) {
	tbl_simple_go_bait <- tbl_reduced_go_result %>% filter(ViralProtein==bait) %>% dplyr::select(-ViralProtein)
	tbl_description_bait <- tbl_all_descriptions %>% left_join(tbl_simple_go_bait) %>% mutate(ViralProtein = bait)
	tbl_ggplot <- rbind(tbl_ggplot, tbl_description_bait)
}

tbl_ggplot_input <- tbl_ggplot %>% mutate(logpadjust = -log10(p.adjust)) %>% 
	dplyr::select(parentTerm, logpadjust, ViralProtein) %>% unique %>% replace(is.na(.), 0) 

tbl_ggplot %>% filter(parentTerm == "cristae formation") %>% dplyr::pull(Description) %>% unique
# matrix_ggplot_input <- tbl_ggplot_input %>% as.matrix()

# GO_list <- c(enrichgo_results$ID)
# list_test <- c("GO:0010918", "GO:2000105", "GO:0006264")

# test <- prioritizedGOTerms(lst=GO_list,organism = "Human", sp = FALSE, domain="BP")

# Terms <-c("GO:0000278", "GO:0006414","GO:0022403","GO:0006415","GO:0006614",
# 	"GO:0045047","GO:0022411","GO:0001775","GO:0046649","GO:045321")
# test <- prioritizedGOTerms(Terms, organism = "Human", sp = FALSE, domain = "BP")

fun_color_range <- colorRampPalette(c("#ffffff", "#0e1c36"))   # Apply colorRampPalette
my_colors <- fun_color_range(25)                          # Extract 100 color codes

plot <- ggplot(tbl_ggplot_input, aes(y=parentTerm, x=ViralProtein, fill=logpadjust)) + 
	geom_tile(colour="grey",size=0.40) +
	scale_fill_viridis_c() +
	labs(title="", x ="", y = "") +
	# scale_fill_gradientn(colors = my_colors, na.value="white") +
	scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "Blues"))(25)) +
	# scale_fill_gradient(low = "white", high = "darkblue")
	labs(fill='-log10(adjusted p-value)') +
	theme(
		axis.ticks= element_blank(),
		axis.text.x = element_text(vjust=0.7, angle=90),
		legend.key.size = unit(0.5, "cm"),
		legend.title = element_text(size = 7, face="bold"),
		legend.text = element_text(size = 7),
		legend.position = c(-1.7, 0.67),
		panel.background = element_blank()
	)

# tbl_ggplot_input_2 <- tbl_ggplot_input %>% filter(!ViralProtein=="nsp8")
# plot_2

ggsave(plot, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Heatmap/GO analysis.png", width = 14.5, height = 17, units = "cm")


heatmap(matrix_ggplot_input)

plot_heatmap <- heatmap(tbl_ggplot_input, Colv = NA, Rowv = NA, scale="logpadjust", ylab="Description", xlab="ViralProtein") + 
	# geom_raster() +
	scale_fill_viridis_c() +
	labs(title="", x ="", y = "") +
	scale_fill_gradientn(colors = coul) +
	# scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(25)) +
	theme(panel.background = element_blank(), axis.ticks= element_blank())

##################################################################################################################################################################################################################################

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

##################################################################################################################################################################################################################################

library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="MF")
##################################################################################################################################################################################################################################

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

##################################################################################################################################################################################################################################

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


