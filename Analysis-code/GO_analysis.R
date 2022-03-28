library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ViSEAGO)
library(Rgraphviz)
library(RColorBrewer)

xx <- compareCluster(gcSample, fun="enrichKEGG",
				organism="hsa", pvalueCutoff=0.05)
xx %>% as_tibble

plot(xx, type="dot", caption="KEGG Enrichment Comparison")

as.data.frame(xx)

mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
                            '100127206', '100128071'),
                   group = c('A', 'A', 'A', 'B', 'B', 'B'),
                   othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))

xx.formula <- compareCluster(Entrez~group, data=mydf, fun='groupGO', OrgDb='org.Hs.eg.db')

as.data.frame(xx.formula)

xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf,
										fun='groupGO', OrgDb='org.Hs.eg.db')

as.data.frame(xx.formula.twogroups)


############################################################################################################################################
##### ENRICHER
#	Load gene selection
SAINT_MIST_filtered <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv") %>% as_tibble
selection <- SAINT_MIST_filtered %>% rename(gene_symbol=PreyGene) %>% select(gene_symbol)
selection_nsp1 <- SAINT_MIST_filtered %>% filter(idBait=="nsp1") %>% rename(gene_symbol=PreyGene) %>% select(gene_symbol)
selection_nsp7 <- SAINT_MIST_filtered %>% filter(idBait=="nsp7") %>% rename(gene_symbol=PreyGene) %>% select(gene_symbol)
selection_nsp8 <- SAINT_MIST_filtered %>% filter(idBait=="nsp8") %>% rename(gene_symbol=PreyGene) %>% select(gene_symbol)
selection_S <- SAINT_MIST_filtered %>% filter(idBait=="S") %>% rename(gene_symbol=PreyGene) %>% select(gene_symbol)

#	Load genes background
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
background <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	rename(gene_symbol = gene_name) %>%
	select(gene_symbol) %>% unique

#	Load term to gene
all_gene_sets = msigdbr(species = "Homo sapiens")
bobp_gene_sets = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")

gene_map <- all_gene_sets %>% mutate(entrez_gene = as.character(entrez_gene)) %>% select(gene_symbol, entrez_gene) %>% unique
selection_genes <- merge(selection, gene_map) %>% as_tibble %>% unique
selection_genes_nsp1 <- merge(selection_nsp1, gene_map) %>% as_tibble %>% unique
selection_genes_nsp7 <- merge(selection_nsp7, gene_map) %>% as_tibble %>% unique
selection_genes_nsp8 <- merge(selection_nsp8, gene_map) %>% as_tibble %>% unique
selection_genes_S <- merge(selection_S, gene_map) %>% as_tibble %>% unique

background_genes <- merge(background, gene_map) %>% as_tibble %>% unique
# enricher(gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
# 		minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE, TERM2NAME = NA)

# enrichGO(gene, bobp_gene_sets, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH",
# 	universe, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)

# msigdbr_t2g = bobp_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# go_result <- enricher(gene = SAINT_MIST_filtered$PreyGene, universe = background, TERM2GENE = msigdbr_t2g)
enrichgo_result <- enrichGO(gene = selection_genes$entrez_gene, org.Hs.eg.db, ont="BP",
	universe = background_genes$entrez_gene)

enrichgo_result_nsp1 <- enrichGO(gene = selection_genes_nsp1$entrez_gene, org.Hs.eg.db, ont="BP",
	universe = background_genes$entrez_gene, pvalueCutoff=0.01)
enrichgo_result_nsp7 <- enrichGO(gene = selection_genes_nsp7$entrez_gene, org.Hs.eg.db, ont="BP",
	universe = background_genes$entrez_gene, pvalueCutoff=0.01)
enrichgo_result_nsp8 <- enrichGO(gene = selection_genes_nsp8$entrez_gene, org.Hs.eg.db, ont="BP",
	universe = background_genes$entrez_gene, pvalueCutoff=0.01)
enrichgo_result_S <- enrichGO(gene = selection_genes_S$entrez_gene, org.Hs.eg.db, ont="BP",
	universe = background_genes$entrez_gene, pvalueCutoff=0.01)

head(enrichgo_result_nsp8)
plot(enrichgo_result)
?enrichGO
fwrite(go_result, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/GO/enricher/go_result.txt")
write.csv(go_result, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/GO/enricher/go_result.txt")

simple_nsp1 <- simplify(enrichgo_result_nsp1)
simple_nsp7 <- simplify(enrichgo_result_nsp7)
simple_nsp8 <- simplify(enrichgo_result_nsp8)

##### VISUALISATION
# Generate matrix
sNsp1 <- simple_nsp1 %>% as_tibble %>% mutate(Nsp1 = -log10(qvalue)) %>% select(Description, Nsp1)
sNsp7 <- simple_nsp7 %>% as_tibble %>% mutate(Nsp7 = -log10(qvalue)) %>% select(Description, Nsp7)
sNsp8 <- simple_nsp8 %>% as_tibble %>% mutate(Nsp8 = -log10(qvalue)) %>% select(Description, Nsp8)

sNsp1_D <- sNsp1 %>% select(Description)
sNsp7_D <- sNsp7 %>% select(Description)
sNsp8_D <- sNsp8 %>% select(Description)

ALL_descriptions <- rbind(sNsp1_D, sNsp7_D, sNsp8_D) %>% as_tibble %>% unique

HM_matrix <- left_join(ALL_descriptions, sNsp1) %>% left_join(sNsp7) %>% left_join(sNsp8) %>% 
	replace(is.na(.), 0) %>% select(-Description) %>% as.matrix()
rownames(HM_matrix) = ALL_descriptions$Description

ssNsp1 <- simple_nsp1 %>% as_tibble %>% mutate(Q = -log10(qvalue)) %>% select(Description, Q) %>% mutate(ViralProtein = "Nsp1")
ssNsp7 <- simple_nsp7 %>% as_tibble %>% mutate(Q = -log10(qvalue)) %>% select(Description, Q) %>% mutate(ViralProtein = "Nsp7")
ssNsp8 <- simple_nsp8 %>% as_tibble %>% mutate(Q = -log10(qvalue)) %>% select(Description, Q) %>% mutate(ViralProtein = "Nsp8")

HM_input <- rbind(ssNsp1, ssNsp7, ssNsp8) %>% as_tibble %>% unique %>% replace(is.na(.), 0)

col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red")

Heatmap(HM_matrix,
	name = 'GO Enrichment', 
	show_row_dend = FALSE, show_column_dend = FALSE, 
	col= colorRampPalette(brewer.pal(8, "Blues"))(25), 
	heatmap_legend_param = list(title="-log10(Q)"))

heatmap(HM_matrix,
	name = 'GO Enrichment', 
	Colv = NA, Rowv = NA, 
	col= colorRampPalette(brewer.pal(8, "Blues"))(25))

ggplot(HM_input, aes(y=Description, x=ViralProtein, fill=Q)) + 
	geom_raster() +
	scale_fill_viridis_c() +
	labs(title="GO enrichment analysis", x ="", y = "") +
	scale_fill_gradientn(colors = colorRampPalette(brewer.pal(8, "Blues"))(25)) +
	theme(panel.background = element_blank(), axis.ticks= element_blank())


enrichgo_list <- list(enrichgo_result_nsp1, enrichgo_result_nsp7, enrichgo_result_nsp8, inherit.name=TRUE)
merge_result(enrichgo_list)
plotGOgraph(enrichgo_result_nsp1)

?simplify()
ViSEAGO::show_table(go_result)

############################################################################################################################################
library(ViSEAGO)
library(msigdbr)

# load genes background
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
background <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	pull(gene_name) %>% unique

# load gene selection
HCIP_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv"
hcip <- fread(HCIP_path) %>% as_tibble
selection <- hcip %>% pull(PreyGene) %>% unique

# load GO annotations
myGENE2GO <- msigdbr(species = "Homo sapiens", subcategory = "GO:BP")

# EntrezGene<-ViSEAGO::EntrezGene2GO()
# myGENE2GO2<-ViSEAGO::annotate(
# 	"9606",
# 	EntrezGene
# )

Bioconductor<-ViSEAGO::Bioconductor2GO()
ViSEAGO::available_organisms(Bioconductor)
myGENE2GO<-ViSEAGO::annotate(
	"org.Hs.eg.db",
	Bioconductor
)

# create topGOdata for BP
BP <- ViSEAGO::create_topGOdata(
	geneSel=selection,
	allGenes=background,
	gene2GO=myGENE2GO, 
	ont="BP",
	nodeSize=5
)

head(myGENE2GO)

# TopGO
library(topGO)

# names(background) = background
sampleGOdata <- new("topGOdata",
	description = "Simple session", ontology = "BP",
	allGenes = background, geneSel = selection,
	nodeSize = 10,
	annot = annFUN.db, affyLib = affyLib)

background <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	pull(gene_name) %>% unique
