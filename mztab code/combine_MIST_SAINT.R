library(tidyverse)
library(data.table)
library(VennDiagram)

#   Gene names
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	select(V2, gene_name)

#   output MiST
# Before you can load the MistOutput the first couple of lines must be removed.
MIST <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/MiST/MistOutput.txt", sep= "\t") %>% as_tibble

MIST_result <- MIST %>% filter(MiST >= 0.7) %>%
	mutate(idBait=Bait) %>%
	select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE)

mist_gene <- merge(MIST, gene_names, by.x = "Prey", by.y = "V2") %>% as_tibble

MIST_gene <- mist_gene %>% filter(MiST >= 0.7) %>%
	mutate(idBait=Bait) %>%
	mutate(PreyGene = gene_name) %>%
	select(idBait, PreyGene, everything()) %>%
	unite(BP, idBait:PreyGene ,sep ="_", remove=FALSE)

#   output SAINT
SAINT <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/list.txt", sep= "\t") %>% as_tibble

SAINT_result <- SAINT %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>%
	mutate(idBait=Bait) %>%
	select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE)

SAINT_gene <- SAINT %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>%
	mutate(idBait=Bait) %>%
	select(idBait, PreyGene, everything()) %>%
	unite(BP, idBait:PreyGene ,sep ="_", remove=FALSE)

#   combine MIST and SAINT output
combined_result <- SAINT_result %>% filter(BP %in% MIST_result$BP)      #in both directions same amount
combined_gene <- SAINT_gene %>% filter(BP %in% MIST_gene$BP)
combined_gene_mist <- MIST_gene %>% filter(BP %in% SAINT_gene$BP)

output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/SAINT<=0.05_MiST>=0.7_AvgSpec>=2.csv"
fwrite(combined_gene, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/SAINT<=0.05_MiST>=0.7_AvgSpec>=2_MistColumns.csv"
fwrite(combined_gene_mist, output_path, append = FALSE, col.names = TRUE) 

################################################################################################################################################

#   load original Gordon et al. results
intact <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/intact_gordon.txt", sep= "\t") %>% as_tibble
intact_map <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results//HCIP/AP-MS/Gordon/name_map_sars_cov_2.txt", sep= "\t")
intact_annot <- merge(intact, intact_map, by.x = "IDinteractorA", by.y = "accession") %>% as_tibble

#   remove space in protein14 (bait)
intact_annot$idBait <- gsub(" ", "", intact_annot$idBait)

intact_result <- intact_annot %>%
	select(idBait, IDinteractorB, everything()) %>%
	unite(BP, idBait:IDinteractorB,sep ="_", remove=FALSE) #%>% select(BP) %>% unique()

intact_G <- merge(intact_annot, gene_names, by.x = "IDinteractorB", by.y = "V2") %>% as_tibble
intact_gene <- intact_G %>%
	select(idBait, gene_name, everything()) %>%
	unite(BP, idBait:gene_name,sep ="_", remove=FALSE)

#   Create plot with overlap
ann <- combined_gene %>% filter(!BP %in% intact_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "ANN-SoLo")
gordon <- intact_gene %>% filter(!BP %in% combined_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "Gordon et al.")
both <- combined_gene %>% filter(BP %in% intact_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "Overlap")

comb <- rbind(ann, gordon, both)

s <- ggplot(comb, aes(idBait, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="HCIP overlap", x = "Baits", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=10, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/HCIP/HCIP overlap.png", s,width = 24, height = 9, units = "cm")

################################################################################################################################################
combined_nsp5 <- combined_gene %>% filter(idBait == "nsp5") 
intact_nsp5 <- intact_gene %>% filter(idBait == "nsp5") %>%
	rename(PreyGene = gene_name) %>%
	select(idBait, BP, PreyGene)
all_nsp5 <- rbind(combined_nsp5, intact_nsp5) %>% unique

combined_nsp5_C145A <- combined_gene %>% filter(idBait == "nsp5_C145A")
intact_nsp5_C145A <- intact_gene %>% filter(idBait == "nsp5_C145A") %>%
	rename(PreyGene = gene_name) %>%
	select(idBait, BP, PreyGene)
all_nsp5_C145A <- rbind(combined_nsp5_C145A, intact_nsp5_C145A) %>% unique

combined_N <- combined_gene %>% filter(idBait == "N") 
intact_N <- intact_gene %>% filter(idBait == "N") %>%
	rename(PreyGene = gene_name) %>%
	select(idBait, BP, PreyGene)
all_N <- rbind(combined_N, intact_N) %>% unique

combined_nsp14 <- combined_gene %>% filter(idBait == "nsp14") 
intact_nsp14 <- intact_gene %>% filter(idBait == "nsp14") %>%
	rename(PreyGene = gene_name) %>%
	select(idBait, BP, PreyGene)
all_nsp14 <- rbind(combined_nsp14,intact_nsp14) %>% unique

output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp5.csv"
fwrite(all_nsp5, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp5_C145A.csv"
fwrite(all_nsp5_C145A, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/N.csv"
fwrite(all_N, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp14.csv"
fwrite(all_nsp14, output_path, append = FALSE, col.names = TRUE) 

combined_gene %>% filter(idBait == "nsp14") 
both %>% filter(idBait == "nsp14") %>% select(BP)
gordon %>% filter(idBait == "nsp14") %>% select(BP)

################################################################################################################################################

#   reproduce
nrow(both)/nrow(intact_result)
nrow(both)/nrow(intact_gene)

nrow(combined_result)/nrow(intact_result) -1
nrow(combined_gene)/nrow(intact_gene) -1

################################################################################################################################################
#   Write gene names to files, to color the nodes in network
combined_genes <- combined_gene %>% select(idBait, BP, PreyGene) %>% unique()

ann <- combined_genes %>% filter(!BP %in% intact_gene$BP) %>% mutate(identification = "ANN-SoLo") %>% unique()
gordon <- intact_gene %>% filter(!BP %in% combined_genes$BP) %>% mutate(identification = "Gordon et al.") %>% 
	mutate(PreyGene = gene_name) %>%
	select(idBait, BP, PreyGene) %>% unique()
both <- combined_genes %>% filter(BP %in% intact_gene$BP) %>% mutate(identification = "Overlap") %>% unique()

both_gene <- both %>% select(PreyGene) %>% unique()
ann_gene <- ann %>% select(PreyGene) %>% unique()
gordon_gene <- gordon %>% select(PreyGene) %>% unique()
bait_gene <- combined_genes %>% select(idBait) %>% unique()

output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/both_gene.csv"
fwrite(both_gene, output_path, append = FALSE, col.names = TRUE)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/ann_gene.csv"
fwrite(ann_gene, output_path, append = FALSE, col.names = TRUE)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/gordon_gene.csv"
fwrite(gordon_gene, output_path, append = FALSE, col.names = TRUE)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/baits.csv"
fwrite(bait_gene, output_path, append = FALSE, col.names = TRUE)

#   Save file that contains both the ANN-SoLo as the Gordon et al results to put into a network
gordon_combined_genes <- rbind(combined_genes, gordon)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/ann_gordon_gene.csv"
fwrite(gordon_combined_genes, output_path, append = FALSE, col.names = TRUE)

################################################################################################################################################

#   Save BP and gene name to create a venn diagram
gordon_bp <- intact_gene %>% select(BP)
gordon_gn <- intact_gene %>% select(gene_name) %>% unique
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/"
fwrite(gordon_bp, paste(output_path, "gordon_bp.txt", sep=""))
fwrite(gordon_gn, paste(output_path, "gordon_gn.txt", sep=""))

annsolo_bp <- combined_gene %>% select(BP)
annsolo_gn <- combined_gene %>% select(PreyGene) %>% unique
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/"
fwrite(annsolo_bp, paste(output_path, "annsolo_bp.txt", sep=""))
fwrite(annsolo_gn, paste(output_path, "annsolo_gn.txt", sep=""))

#   Load the Li et al. results
Li <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/AHCIP/P-MS/Li/SARS-CoV-2-host interactome list.tsv", sep= "\t") %>% as_tibble
Li_bp <- Li %>% select(Bait, GeneName) %>% 
	unite(BP, Bait:GeneName ,sep ="_", remove=TRUE)
Li_gn <- Li %>% select(GeneName) %>% unique()
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/"
fwrite(Li_bp, paste(output_path, "Li_bp.txt", sep=""))
fwrite(Li_gn, paste(output_path, "Li_gn.txt", sep=""))

#   Load the Stukalov et al. results
Stukalov <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/SARS-CoV-2_interactions.tsv", sep= "\t") %>% as_tibble
Stukalov_bp <- Stukalov %>% select(bait_name, gene_name) %>% 
	unite(BP, bait_name:gene_name ,sep ="_", remove=TRUE)
Stukalov_gn <- Stukalov %>% select(gene_name)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/"
fwrite(Stukalov_bp, paste(output_path, "Stukalov_bp.txt", sep=""))
fwrite(Stukalov_gn, paste(output_path, "Stukalov_gn.txt", sep=""))

#   Load the Laurent et al. results
Laurent <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/all_high_conf_interactions.tsv", sep= "\t") %>% as_tibble
Laurent_bp <- Laurent %>% select(Bait, Gene_name) %>% 
	unite(BP, Bait:Gene_name ,sep ="_", remove=TRUE)
Laurent_gn <- Laurent %>% select(Gene_name)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/"
fwrite(Laurent_bp, paste(output_path, "Laurent_bp.txt", sep=""))
fwrite(Laurent_gn, paste(output_path, "Laurent_gn.txt", sep=""))

#   Load the Gingras group results
Gingras <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Gingras/saint-virus.txt", sep= "\t") %>% as_tibble
Gingras_bp <- Gingras %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>%
	select(Bait, PreyGene) %>% 
	unite(BP, Bait:PreyGene ,sep ="_", remove=TRUE)
Gingras_gn <- Gingras %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>% select(PreyGene)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Gingras/"
fwrite(Gingras_bp, paste(output_path, "Gingras_bp.txt", sep=""))
fwrite(Gingras_gn, paste(output_path, "Gingras_gn.txt", sep=""))

#   Venn Diagram
# venn.plot <- draw.pairwise.venn(322, 612, 210, category = c("Gordon et al.", "ANN-SoLo"), 
# 	scaled = TRUE, inverted = TRUE, ext.text = TRUE, ext.percent = rep(0.05, 3),
# 	lwd = rep(2, 2), lty = rep("solid", 2), col = rep("black",2),
# 	fill = NULL, alpha = rep(0.5, 2), label.col = rep("black", 3),
# 	cex = rep(1, 3), fontface = rep("plain", 3), fontfamily = rep("Helvetica", 3),
# 	cat.pos = c(-50, 50), cat.dist = rep(0.025, 2), cat.cex = rep(1, 2), cat.col = rep("black", 2),
# 	cat.fontface = rep("plain", 2), cat.fontfamily = rep("Helvetica", 2), cat.just = rep(list(c(0.5, 0.5)), 2),
# 	cat.default.pos = "outer", cat.prompts = FALSE, ext.pos = rep(0, 2), ext.dist = rep(0, 2),
# 	ext.line.lty = "solid", ext.length = rep(0.95, 2), ext.line.lwd = 1, rotation.degree = 0,
# 	rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop = NULL,
# 	print.mode = "raw", sigdigs = 3)
# grid.draw(venn.plot);
# grid.newpage();

venn <- draw.pairwise.venn(322, 612, 210,
	category = rep(c("Gordon et al.", "ANN-SoLo")), cex = 2, cat.cex = 2,
	fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
	rscaled = TRUE, cat.dist = c(0.04, 0.04),
	cat.pos = c(-14, 10))
grid.draw(venn)

ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/double_Venn.pdf", width = 17, height = 17, units = "cm")

