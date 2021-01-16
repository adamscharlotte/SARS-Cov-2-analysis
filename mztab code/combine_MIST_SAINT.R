library(tidyverse)
library(data.table)
library(VennDiagram)

#   Load gene names
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	select(V2, gene_name)

#	Load MiST output
# 	Before you can load the MistOutput the first couple of lines must be removed.
MIST <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/MiST/MistOutput.txt", sep= "\t") %>% as_tibble
MIST_output <- merge(MIST, gene_names, by.x = "Prey", by.y = "V2") %>% 		# Add gene names
	as_tibble %>% mutate(PreyGene = gene_name) %>% mutate(idBait=Bait) %>%
	filter(!MiST == 0) %>% select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene ,sep ="_", remove=FALSE)

#	Load SAINT output
SAINT_output <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/list.txt", sep= "\t") %>% as_tibble %>%
	mutate(idBait=Bait) %>% select(idBait, Prey, everything()) %>%
	unite(BP, idBait:Prey ,sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene ,sep ="_", remove=FALSE)

#	Combine the MIST and SAINT output
MIST_filter <- MIST_output %>% filter(BP %in% SAINT_output$BP)
SAINT_filter <- SAINT_output %>% filter(BP %in% MIST_output$BP)
SAINT_MIST <- merge(MIST_filter, SAINT_filter) %>% as_tibble

#	Filter the results
SAINT_MIST_filtered <- SAINT_MIST %>% filter(MiST >= 0.7 & BFDR <= 0.05 & AvgSpec >= 2)

fwrite(SAINT_MIST, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/Scoring results.csv", append = FALSE, col.names = TRUE)
fwrite(SAINT_MIST_filtered, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/HCIPs.csv", append = FALSE, col.names = TRUE)

################################################################################################################################################

#   load original Gordon et al. results
Gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/intact_gordon.txt", sep= "\t") %>% as_tibble
Gordon_map <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results//HCIP/AP-MS/Gordon/name_map_sars_cov_2.txt", sep= "\t")
Gordon_annot <- merge(Gordon, Gordon_map, by.x = "IDinteractorA", by.y = "accession") %>% as_tibble

#   remove space in protein14 (bait)
Gordon_annot$idBait <- gsub(" ", "", Gordon_annot$idBait)
Gordon_gene <- merge(Gordon_annot, gene_names, by.x = "IDinteractorB", by.y = "V2") %>% as_tibble

Gordon_result <- Gordon_gene %>% mutate(Prey = IDinteractorB) %>% mutate(PreyGene = gene_name) %>%
	select(idBait, Prey, everything()) %>% unite(BP, idBait:Prey, sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene, sep ="_", remove=FALSE)

#   Create plot with overlap
ann <- SAINT_MIST_filtered %>% filter(!BP_gene %in% Gordon_result$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "ANN-SoLo") %>% unique
gordon <- Gordon_result %>% filter(!BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "Gordon et al.") %>% unique
both <- SAINT_MIST_filtered %>% filter(BP_gene %in% Gordon_result$BP_gene) %>% select(idBait, BP_gene) %>% mutate(Identification = "Overlap") %>% unique

comb <- rbind(ann, gordon, both) 

s <- ggplot(comb, aes(idBait, fill = Identification)) +
	geom_bar(position = "stack") +
	theme_minimal() +
	labs(title ="HCIP overlap", x = "Baits", y = "Number of identifications") + 
	scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
	theme(axis.text.x = element_text(size=10, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/HCIP/Overlap/HCIP overlap.png", s,width = 24, height = 9, units = "cm")

################################################################################################################################################

#   Write gene names to files, to color the nodes in network
both_gene <- SAINT_MIST_filtered %>% filter(BP_gene %in% Gordon_result$BP_gene) %>% select(PreyGene)
ann_gene <- SAINT_MIST_filtered %>% filter(!BP_gene %in% Gordon_result$BP_gene) %>% select(PreyGene)
gordon_gene <- Gordon_result %>% filter(!BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% select(PreyGene)
bait_gene <- SAINT_MIST_filtered %>% select(idBait) %>% unique()

ann_gene %>% filter(PreyGene %in% both_gene$PreyGene)
Gordon_result %>% filter(PreyGene==c("MOGS", "NARS2"))

output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/"
fwrite(both_gene, paste(output_path, "both_gene.txt", sep=""))
fwrite(ann_gene, paste(output_path, "ann_gene.txt", sep=""))
fwrite(gordon_gene, paste(output_path, "gordon_gene.txt", sep=""))
fwrite(bait_gene, paste(output_path, "bait_gene.txt", sep=""))

#   Save file that contains both the ANN-SoLo as the Gordon et al results to put into a network
ANNSoLo <- SAINT_MIST_filtered %>% select(BP_gene, idBait, PreyGene)
Gordonetal <- Gordon_result %>% select(BP_gene, idBait, PreyGene)
gordon_combined_genes <- rbind(ANNSoLo, Gordonetal) %>% unique
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/Network/ann_gordon_gene.csv"
fwrite(gordon_combined_genes, output_path, append = FALSE, col.names = TRUE)

################################################################################################################################################

#   Save BP and gene name to create a venn diagram
gordon_bp <- Gordon_result %>% select(BP_gene)
gordon_gn <- Gordon_result %>% select(PreyGene) %>% unique
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/"
fwrite(gordon_bp, paste(output_path, "gordon_bp.txt", sep=""))
fwrite(gordon_gn, paste(output_path, "gordon_gn.txt", sep=""))

annsolo_bp <- SAINT_MIST_filtered %>% select(BP_gene)
annsolo_gn <- SAINT_MIST_filtered %>% select(PreyGene) %>% unique
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

#	Create the venn diagram
venn <- draw.pairwise.venn(332, 553, 200,
	category = rep(c("Gordon et al.", "ANN-SoLo")), cex = 2, cat.cex = 2,
	fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
	rscaled = TRUE, cat.dist = c(0.04, 0.04),
	cat.pos = c(-14, 10))
grid.draw(venn)

ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Overlap.pdf", width = 17, height = 17, units = "cm")

################################################################################################################################################

#   Load original Gordon et al. results
Gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/Scoring results.csv", header = TRUE) %>% as_tibble
Gordon_output <- Gordon %>% separate(Bait, into=c("x", "idBait"), sep=" ") %>%
	mutate(idBait = str_replace(idBait,"Spike", "S")) %>%
	mutate(idBait = str_replace(idBait,"orf9c", "orf9c(protein14)")) %>%
	mutate(BFDR = Saint_BFDR) %>% mutate(MiST = MIST) %>% mutate(Prey = Preys) %>%
	select(idBait, Prey, everything()) %>% unite(BP, idBait:Prey, sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene, sep ="_", remove=FALSE)
Gordon_filtered <- Gordon_output %>% filter(BP_gene %in% Gordon_result$BP_gene)

#	Plot the MiST scoring against the BFDR
Gordon_input <- Gordon_output %>% mutate(MiST_Gordon = MiST) %>% mutate(BFDR_Gordon = BFDR) %>%
	filter(BP %in% Gordon_result$BP) %>% select(BP_gene, BP, idBait, PreyGene, Prey, MiST_Gordon, BFDR_Gordon)
ANNSoLo_input <- SAINT_MIST %>% mutate(MiST_ANNSoLo = MiST) %>% mutate(BFDR_ANNSoLo = BFDR) %>%
	filter(BP %in% SAINT_MIST_filtered$BP) %>% select(BP_gene, BP, idBait, PreyGene, Prey, MiST_ANNSoLo, BFDR_ANNSoLo)
Gordon_add_input <- Gordon_output %>% mutate(MiST_Gordon = MiST) %>% mutate(BFDR_Gordon = BFDR) %>%
	filter(BP %in% Gordon_result$BP | BP %in% SAINT_MIST_filtered$BP) %>% select(BP_gene, BP, idBait, PreyGene, Prey, MiST_Gordon, BFDR_Gordon)
ANNSoLo_add_input <- SAINT_MIST %>% mutate(MiST_ANNSoLo = MiST) %>% mutate(BFDR_ANNSoLo = BFDR) %>%
	filter(BP %in% Gordon_result$BP | BP %in% SAINT_MIST_filtered$BP) %>% select(BP_gene, BP, idBait, PreyGene, Prey, MiST_ANNSoLo, BFDR_ANNSoLo)
Plot_input <- merge(Gordon_add_input, ANNSoLo_add_input) %>% as_tibble
Plot_input_label <- merge(Plot_input, comb) %>% as_tibble

s <- ggplot(Plot_input_label, aes(MiST_Gordon, MiST_ANNSoLo, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="MiST scores of Gordon et al. plotted against those of ANN-SoLo", x = "MiST score Gordon et al.", y = "MiST score ANN-SoLo") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.7, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#071E22", "#849B96", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/MiST scores of Gordon et al. plotted against those of ANN-SoLo.pdf", width = 17, height = 17, units = "cm")

s <- ggplot(Plot_input_filtered_label, aes(BFDR_Gordon, BFDR_ANNSoLo, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="SAINT BFDR scores of Gordon et al. plotted against those of ANN-SoLo", x = "BFDR Gordon et al.", y = "BFDR ANN-SoLo") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#071E22", "#849B96", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINT BFDR scores of Gordon et al. plotted against those of ANN-SoLo.pdf", width = 17, height = 17, units = "cm")

s <- ggplot(Gordon_input, aes(BFDR_Gordon, MiST_Gordon)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="SAINT BFDR and MiST scores of Gordon et al. plotted against each other", x = "BFDR", y = "MiST score") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red")
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINT BFDR and MiST scores of Gordon et al. plotted against each other.pdf", width = 17, height = 17, units = "cm")

s <- ggplot(ANNSoLo_input, aes(BFDR_ANNSoLo, MiST_ANNSoLo)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="SAINT BFDR and MiST scores of ANN-SoLo plotted against each other", x = "BFDR", y = "MiST score") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red")
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINT BFDR and MiST scores of ANN-SoLo plotted against each other.pdf", width = 17, height = 17, units = "cm")


Gordon_output %>% summarise(a = mean(MiST), b = mean(BFDR))
SAINT_MIST %>% summarise(a = mean(MiST), b = mean(BFDR))