# devtools::install_github('andreacirilloac/updateR')
# library(updateR)
# chooseCRANmirror(4)
# install.packages('devtools')
# library(devtools)
# devtools::install_github("tidyverse/tidyverse")
# devtools::install_github("yanlinlin82/ggvenn")
# install.packages("data.table")
library(tidyverse)
library(data.table)
library("ggvenn")

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
Gordon_result %>% filter(PreyGene==c("NPTX1", "NARS2"))

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
fwrite(gordon_bp, paste(output_path, "gordon_bp.txt", sep=""), col.names = FALSE)
fwrite(gordon_gn, paste(output_path, "gordon_gn.txt", sep=""), col.names = FALSE)

annsolo_bp <- SAINT_MIST_filtered %>% select(BP_gene) %>% unique
annsolo_gn <- SAINT_MIST_filtered %>% select(PreyGene) %>% unique
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/"
fwrite(annsolo_bp, paste(output_path, "annsolo_bp.txt", sep=""), col.names = FALSE)
fwrite(annsolo_gn, paste(output_path, "annsolo_gn.txt", sep=""), col.names = FALSE)

#   Load the Li et al. results
Li <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/Li.txt", sep= "\t") %>% as_tibble
# Li %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Li_sars %>% pull(OfficialSymbolInteractorA) %>% unique
Li_sars <- Li %>% filter(OrganismNameInteractorB == "Homo sapiens")
Li_bp <- Li_sars %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>% 
	mutate(OfficialSymbolInteractorA = str_replace(OfficialSymbolInteractorA,"orf14","orf9c(protein14)")) %>%
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Li_gn <- Li_sars %>% select(OfficialSymbolInteractorB) %>% unique()
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Li/"
fwrite(Li_bp, paste(output_path, "Li_bp.txt", sep=""), col.names = FALSE)
fwrite(Li_gn, paste(output_path, "Li_gn.txt", sep=""), col.names = FALSE)

#   Load the Chen et al. results
Chen <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Chen/Chen.txt", sep= "\t") %>% as_tibble
# Chen %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Chen %>% pull(OfficialSymbolInteractorA) %>% unique
Chen_bp <- Chen %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>%
	mutate(OfficialSymbolInteractorA = str_replace(OfficialSymbolInteractorA,"orf14","orf9c(protein14)")) %>%
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Chen_gn <- Chen %>% select(OfficialSymbolInteractorB) %>% unique()
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Chen/"
fwrite(Chen_bp, paste(output_path, "Chen_bp.txt", sep=""), col.names = FALSE)
fwrite(Chen_gn, paste(output_path, "Chen_gn.txt", sep=""), col.names = FALSE)

#   Load the Stukalov et al. results
Stukalov <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/Stukalov.txt", sep= "\t") %>% as_tibble
# Stukalov %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Stukalov_sars %>% pull(OfficialSymbolInteractorA) %>% unique
Stukalov_sars <- Stukalov %>% filter(OrganismNameInteractorB=="Homo sapiens" & OrganismNameInteractorA=="Severe acute respiratory syndrome coronavirus 2")
Stukalov_bp <- Stukalov_sars %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>% 
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Stukalov_gn <- Stukalov_sars %>% select(OfficialSymbolInteractorB)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Stukalov/"
fwrite(Stukalov_bp, paste(output_path, "Stukalov_bp.txt", sep=""), col.names = FALSE)
fwrite(Stukalov_gn, paste(output_path, "Stukalov_gn.txt", sep=""), col.names = FALSE)

#   Load the Stukalov et al. results
Germain <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Germain/St-Germain.txt", sep= "\t") %>% as_tibble
# Germain %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Germain %>% pull(OfficialSymbolInteractorA) %>% unique
Germain_bp <- Germain %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>% 
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Germain_gn <- Germain %>% select(OfficialSymbolInteractorB)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Germain/"
fwrite(Germain_bp, paste(output_path, "Germain_bp.txt", sep=""), col.names = FALSE)
fwrite(Germain_gn, paste(output_path, "Germain_gn.txt", sep=""), col.names = FALSE)

#   Load the Laurent et al. results
Laurent <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/Laurent.txt", sep= "\t") %>% as_tibble
# Laurent %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
Laurent_bp <- Laurent %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>% 
	mutate(OfficialSymbolInteractorA = str_replace(OfficialSymbolInteractorA,"orf14","orf9c(protein14)")) %>%
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Laurent_gn <- Laurent %>% select(OfficialSymbolInteractorB)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Laurent/"
fwrite(Laurent_bp, paste(output_path, "Laurent_bp.txt", sep=""), col.names = FALSE)
fwrite(Laurent_gn, paste(output_path, "Laurent_gn.txt", sep=""), col.names = FALSE)

#   Load the Samavarchi-Tehrani et al. results
Samavarchi <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Samavarchi/Samavarchi.txt", sep= "\t") %>% as_tibble
# Samavarchi %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Samavarchi %>% pull(OfficialSymbolInteractorA) %>% unique
Samavarchi_bp <- Samavarchi %>% select(OfficialSymbolInteractorA, OfficialSymbolInteractorB) %>% 
	mutate(OfficialSymbolInteractorA = str_replace(OfficialSymbolInteractorA,"orf14","orf9c(protein14)")) %>%
	unite(BP, OfficialSymbolInteractorA:OfficialSymbolInteractorB ,sep ="_", remove=TRUE)
Samavarchi_gn <- Samavarchi %>% select(OfficialSymbolInteractorB)
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/BioID/Samavarchi/"
fwrite(Samavarchi_bp, paste(output_path, "Samavarchi_bp.txt", sep=""), col.names = FALSE)
fwrite(Samavarchi_gn, paste(output_path, "Samavarchi_gn.txt", sep=""), col.names = FALSE)

#	Create the venn diagram
venn <- draw.pairwise.venn(332, 375, 164,
	category = rep(c("Gordon et al.", "ANN-SoLo")), cex = 2, cat.cex = 2,
	fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
	rscaled = TRUE, cat.dist = c(0.04, 0.04),
	cat.pos = c(-14, 10))
grid.draw(venn)

ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Overlap.pdf", width = 17, height = 17, units = "cm")

#	Create a piechart
ann_data <- data.frame(Legend=c("Overlap with the original Gordon et al. results", "Overlap with other studies", "Interactions uniquely found with ANN-SoLo"),
	value=c(164,47,211))
ggplot(ann_data, aes(x="", y=value, fill=Legend)) +
	geom_bar(stat="identity", width=1, color="white") +
	coord_polar("y", start=0) +
	scale_fill_manual(values=c("#071E22", "#FF1E1E", "#FF8989")) +
	theme_void()

################################################################################################################################################

#   Load original Gordon et al. results
Gordon <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/Gordon/Scoring results.csv", header = TRUE) %>% as_tibble
Gordon_output <- Gordon %>% separate(Bait, into=c("x", "idBait"), sep=" ") %>%
	mutate(idBait = str_replace(idBait,"Spike", "S")) %>%
	mutate(idBait = str_replace(idBait,"orf9c", "orf9c(protein14)")) %>%
	mutate(BFDR = Saint_BFDR) %>% mutate(MiST = MIST) %>% mutate(Prey = Preys)
Gordon_output_gene <- merge(Gordon_output, gene_names, by.x = "Preys", by.y = "V2") %>% as_tibble %>%
	mutate(PreyGene = gene_name) %>%
	select(idBait, Prey, everything()) %>% unite(BP, idBait:Prey, sep ="_", remove=FALSE) %>%
	select(idBait, PreyGene, everything()) %>% unite(BP_gene, idBait:PreyGene, sep ="_", remove=FALSE) 
Gordon_filtered <- Gordon_output_gene %>% filter(BP_gene %in% Gordon_result$BP_gene)
	
#	Plot the MiST scoring against the BFDR
Gordon_add_input <- Gordon_output_gene %>% mutate(MiST_Gordon = MiST) %>% 
	mutate(BFDR_Gordon = BFDR) %>% mutate(AvgSpec_Gordon = AvgSpec) %>%
	filter(BP_gene %in% Gordon_result$BP_gene | BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% 
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_Gordon, BFDR_Gordon, AvgSpec_Gordon)
ANNSoLo_add_input <- SAINT_MIST %>% mutate(MiST_ANNSoLo = MiST) %>% 
	mutate(BFDR_ANNSoLo = BFDR) %>% mutate(AvgSpec_ANNSoLo = AvgSpec) %>%
	filter(BP_gene %in% Gordon_result$BP_gene | BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% 
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_ANNSoLo, BFDR_ANNSoLo, AvgSpec_ANNSoLo)
Plot_input <- merge(Gordon_add_input, ANNSoLo_add_input) %>% as_tibble
Plot_input_label <- merge(Plot_input, comb) %>% as_tibble
Plot_input_label_ann <- Plot_input_label %>% filter(Identification == "ANN-SoLo")

missed_int <- comb %>% filter(!BP_gene %in%Plot_input_label$BP_gene) %>% select(BP_gene)
ANNSoLo_add_input %>% filter(BP_gene %in%missed_int$BP_gene) %>% pull(BP_gene)

Plot_input_label_ann %>% summarize(avrg=mean(AvgSpec_ANNSoLo) ,max = max(AvgSpec_ANNSoLo), min = min(AvgSpec_ANNSoLo))

Plot_input_label_ann %>% select(AvgSpec_ANNSoLo, AvgSpec_Gordon) %>% 
	mutate(AvgSpecDif=AvgSpec_ANNSoLo-AvgSpec_Gordon) %>%
	summarize(avrg=mean(AvgSpecDif), max= max(AvgSpecDif), min=min(AvgSpecDif))

Venn_ann <- Plot_input_label_ann %>% select(BP_gene, AvgSpec_Gordon, BFDR_Gordon, MiST_Gordon) %>%
	mutate(MiST=if_else(!MiST_Gordon>=0.7, TRUE, FALSE)) %>%
	mutate(BFDR=if_else(!BFDR_Gordon<=0.05, TRUE, FALSE)) %>%
	mutate(AvgSpec=if_else(!AvgSpec_Gordon>=2, TRUE, FALSE))
	# mutate(label=case_when(!MiST_Gordon>=0.7 ~ "MiST", 
	# !BFDR_Gordon<=0.05 ~ "BFDR",
	# !AvgSpec_Gordon>=2 ~ "AvgSpec"))

venn <- ggvenn(Venn_ann, columns = c("MiST", "BFDR", "AvgSpec"),
	show_percentage=FALSE,
	fill_color=c("white","white","white"))

ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Filters.pdf", width = 17, height = 17, units = "cm")


	

s <- ggplot(Plot_input_label_ann, aes(MiST_Gordon, MiST_ANNSoLo, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="MiST scores of Gordon et al. plotted against those of ANN-SoLo", x = "MiST score Gordon et al.", y = "MiST score ANN-SoLo") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.7, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#071E22", "#849B96", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/MiST scores of Gordon et al. plotted against those of ANN-SoLo.pdf", width = 17, height = 17, units = "cm")

s <- ggplot(Plot_input_label_ann, aes(BFDR_Gordon, BFDR_ANNSoLo, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="SAINT BFDR scores of Gordon et al. plotted against those of ANN-SoLo", x = "BFDR Gordon et al.", y = "BFDR ANN-SoLo") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#071E22", "#849B96", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINT BFDR scores of Gordon et al. plotted against those of ANN-SoLo.pdf", width = 17, height = 17, units = "cm")


Gordon_output %>% filter(BP_gene == "nsp6_ATP5O") %>% pull(MiST)
SAINT_MIST_filtered %>% filter(BP_gene == "nsp6_ATP5O") %>% pull(MiST)

Gordon_output %>% filter(BP_gene == "orf8_GLG1") %>% pull(BFDR)
SAINT_MIST_filtered %>% filter(BP_gene == "orf8_GLG1") %>% pull(BFDR)

Gordon_output %>% filter(BP_gene == "nsp6_ATP5O") %>% select(CtrlCounts, Spec)
SAINT_MIST_filtered %>% filter(BP_gene == "nsp6_ATP5O") %>% select(ctrlCounts, Spec)

Plot_input %>% filter(BP_gene == "orf3a_ARMCX3") %>% pull(MiST_Gordon)
Plot_input %>% filter(BP_gene == "orf3a_ARMCX3") %>% pull(MiST_ANNSoLo)
Plot_input %>% filter(BP_gene == "orf3a_ARMCX3") %>% pull(BFDR_Gordon)
Plot_input %>% filter(BP_gene == "orf3a_ARMCX3") %>% pull(BFDR_ANNSoLo)

Gordon_output %>% filter(PreyGene == "SCAP") %>% select(idBait,Prey,  MIST, Spec)
SAINT_MIST %>% filter(PreyGene == "SCAP") %>% select(idBait, Prey, MiST, Spec)