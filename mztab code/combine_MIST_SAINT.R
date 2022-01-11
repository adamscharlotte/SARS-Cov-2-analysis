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
library("ggpubr")

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
	filter(!MiST==0) %>%
	select(idBait, Prey, everything()) %>%
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

#	"/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/GO/Viral interaction partners"

SAINT_MIST_filtered %>% select(idBait) %>% unique %>% print(n=40)
nsp14 <- SAINT_MIST_filtered %>% filter(idBait=="nsp14") %>% select(PreyGene)

fwrite(nsp14, "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/GO/Viral interaction partners/nsp14.txt", append = FALSE, col.names = FALSE)


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

annsolo_bp <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_bp.txt", sep= "\t", header=FALSE) %>% as_tibble
annsolo_gn <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/AP-MS/ANN-SoLo/annsolo_gn.txt", sep= "\t", header=FALSE) %>% as_tibble

annsolo_bp %>% filter (V1 %in% Kim_bp$BP)
annsolo_gn %>% filter (V1 %in% Kim_gn$HostProtein)
annsolo_bp %>% filter(str_detect(V1, "MKRN3"))
#   Load the Li et al. results
Kim <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/Y2H/Kim/Kim.txt", sep= "\t") %>% as_tibble
# Li %>% select(OrganismNameInteractorA, OrganismNameInteractorB) %>% unique
# Li_sars %>% pull(OfficialSymbolInteractorA) %>% unique
Kim_bp <- Kim %>% select(ViralProtein, HostProtein) %>% 
	mutate(ViralProtein = str_replace(ViralProtein,"ORF3B","orf3b")) %>%
	mutate(ViralProtein = str_replace(ViralProtein,"ORF6","orf6")) %>%
	mutate(ViralProtein = str_replace(ViralProtein,"ORF7B","orf7b")) %>%
	mutate(ViralProtein = str_replace(ViralProtein,"NSP","nsp")) %>%
	mutate(ViralProtein = str_replace(ViralProtein,"ORF9C","orf9c(protein14)")) %>%
	unite(BP, ViralProtein:HostProtein ,sep ="_", remove=TRUE)
Kim_gn <- Kim %>% select(HostProtein) %>% unique()
output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/HCIP/Y2H/Kim/"
fwrite(Kim_bp, paste(output_path, "Kim_bp.txt", sep=""), col.names = FALSE)
fwrite(Kim_gn, paste(output_path, "Kim_gn.txt", sep=""), col.names = FALSE)
Kim %>% count(ViralProtein) %>% arrange(desc(n))
annsolo_bp %>% print(n=400)

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

SAINT_MIST %>% select(BP_gene) %>% unique
#	Look at overlap in unfiltered data	
Gordon_missed <- Gordon_output_gene %>% filter(BP_gene %in% SAINT_MIST_filtered$BP_gene)
ANN_missed <- SAINT_MIST %>% filter(BP_gene %in% Gordon_filtered$BP_gene)
count(AvgSpec) %>% arrange(desc(n))
SAINT_MIST %>% filter(!BP_gene %in% Gordon_output_gene$BP_gene)

ANN_missed %>% filter(BP_gene %in% Gordon_our_filter$BP_gene)
Gordon_our_filter <- Gordon_output_gene %>% filter(MIST >= 0.7 & Saint_BFDR <= 0.05 & AvgSpec >= 2) 
Gordon_our_filter %>% filter(BP_gene %in% SAINT_MIST_filtered$BP_gene)

#	Plot the MiST scoring against the BFDR
Gordon_add_input <- Gordon_output_gene %>% mutate(MiST_Gordon = MiST) %>% 
	mutate(BFDR_Gordon = BFDR) %>% mutate(AvgSpec_Gordon = AvgSpec) %>%
	mutate(Spec_Gordon = Spec) %>% mutate(CtrlCounts_Gordon = CtrlCounts) %>%
	filter(BP_gene %in% Gordon_result$BP_gene | BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% 
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_Gordon, BFDR_Gordon, AvgSpec_Gordon, Spec_Gordon, CtrlCounts_Gordon)
ANNSoLo_add_input <- SAINT_MIST %>% mutate(MiST_ANNSoLo = MiST) %>% 
	mutate(BFDR_ANNSoLo = BFDR) %>% mutate(AvgSpec_ANNSoLo = AvgSpec) %>%
	mutate(Spec_ANNSoLo = Spec) %>% mutate(CtrlCounts_ANNSoLo = ctrlCounts) %>%
	filter(BP_gene %in% Gordon_result$BP_gene | BP_gene %in% SAINT_MIST_filtered$BP_gene) %>% 
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_ANNSoLo, BFDR_ANNSoLo, AvgSpec_ANNSoLo, Spec_ANNSoLo, CtrlCounts_ANNSoLo)
Plot_input <- merge(Gordon_add_input, ANNSoLo_add_input) %>% as_tibble
Plot_input_label <- merge(Plot_input, comb) %>% as_tibble
Plot_input_label_ann <- Plot_input_label %>% filter(Identification == "ANN-SoLo")
Plot_input_label_gorboth <- Plot_input_label %>% filter(Identification == "Gordon et al." | Identification == "Overlap")
Plot_input_label_annboth <- Plot_input_label %>% filter(Identification == "ANN-SoLo" | Identification == "Overlap")

#	Plot original scores for the identified PPIs

s <- ggplot(Plot_input_label_annboth, aes(BFDR_Gordon, MiST_Gordon, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="Original BFDR plotted against the original MiST scores", x = "BFDR Gordon et al.", y = "MiST score Gordon et al.") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#071E22", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/Original BFDR plotted against the original MiST scores.pdf", width = 17, height = 17, units = "cm")

Plot_input_label_gorboth
s <- ggplot(Plot_input_label_gorboth, aes(BFDR_ANNSoLo, MiST_ANNSoLo, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(title ="Our BFDR and MiST scores for the originally identified PPIs", x = "BFDR", y = "MiST score") + 
	theme(axis.text.x = element_text(size=10, angle=90)) +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	scale_color_manual(values=c("#849B96", "#377563"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/Our BFDR and MiST scores for the originally identified PPIs.pdf", width = 17, height = 10, units = "cm")

s <- ggplot(Plot_input_label, aes(BFDR_ANNSoLo, BFDR_Gordon, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(x = "reanalysis BFDR", y = "original BFDR") + 
	theme(axis.text.x = element_text(size=10, angle=90), legend.position="none") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
	xlim(0, 0.7) +
	ylim(0, 0.7) +
	scale_color_manual(values=c("#C0D2F7", "#F33B16", "#0E1C36"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINTexpress BFDR of Gordon et al. plotted against those from the reanalysis (0.7).pdf", width = 10, height = 10, units = "cm")

s <- ggplot(Plot_input_label, aes(MiST_ANNSoLo, MiST_Gordon, color=Identification)) +
	geom_point(alpha = 0.7) +
	theme_minimal() +
	labs(x = "reanalysis MiST", y = "original MiST") + 
	theme(axis.text.x = element_text(size=10, angle=90), legend.position="none") +
	geom_vline(xintercept=0.7, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	xlim(0.1, 1) +
	ylim(0.1, 1) +
	scale_color_manual(values=c("#C0D2F7", "#F33B16", "#0E1C36"))
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/MiST score of Gordon et al. plotted against those from the reanalysis.pdf", width = 10, height = 10, units = "cm")

################################################################################################################################################################################################################################################

Gordon_all <- Gordon_output_gene %>% filter(BP %in% SAINT_MIST$BP) %>%
	mutate(MiST_Gordon = MiST) %>% mutate(BFDR_Gordon = BFDR) %>% 
	mutate(AvgSpec_Gordon = AvgSpec) %>% mutate(Spec_Gordon = Spec) %>% 
	mutate(CtrlCounts_Gordon = CtrlCounts) %>%
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_Gordon, BFDR_Gordon, AvgSpec_Gordon, Spec_Gordon, CtrlCounts_Gordon)

ANNSoLo_all <- SAINT_MIST %>% filter(BP %in% Gordon_output_gene$BP) %>%
	mutate(MiST_ANNSoLo = MiST) %>% mutate(BFDR_ANNSoLo = BFDR) %>% 
	mutate(AvgSpec_ANNSoLo = AvgSpec) %>%
	mutate(Spec_ANNSoLo = Spec) %>% mutate(CtrlCounts_ANNSoLo = ctrlCounts) %>%
	select(BP_gene, BP, idBait, PreyGene, Prey, MiST_ANNSoLo, BFDR_ANNSoLo, AvgSpec_ANNSoLo, Spec_ANNSoLo, CtrlCounts_ANNSoLo)

All_unfiltered <- merge(Gordon_all, ANNSoLo_all) %>% as_tibble

#	Plot
s <- ggplot(All_unfiltered, aes(BFDR_ANNSoLo, BFDR_Gordon)) +
	geom_point(alpha = 0.2) +
	theme_minimal() +
	labs(x = "reanalysis BFDR", y = "original BFDR") + 
	theme(axis.text.x = element_text(size=10, angle=90), legend.position="none") +
	geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.05, linetype="dashed", color = "red") +
	xlim(0, 0.7) +
	ylim(0, 0.7) 
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/SAINTexpress BFDR of all unfiltered PPIs.pdf", width = 10, height = 10, units = "cm")

s <- ggplot(All_unfiltered, aes(MiST_ANNSoLo, MiST_Gordon)) +
	geom_point(alpha = 0.2) +
	theme_minimal() +
	labs(x = "reanalysis MiST", y = "original MiST") + 
	theme(axis.text.x = element_text(size=10, angle=90), legend.position="none") +
	geom_vline(xintercept=0.7, linetype="dashed", color = "red") +
	geom_hline(yintercept=0.7, linetype="dashed", color = "red") +
	xlim(0, 1) +
	ylim(0, 1)
ggsave(s, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Explorative plots/MiST score of all unfiltered PPIs.pdf", width = 10, height = 10, units = "cm")

################################################################################################################################################################################################################################################

#	Boxplot CNTR counts

control_Gordon <- All_unfiltered %>% select(CtrlCounts_Gordon) %>%
	separate(CtrlCounts_Gordon, into=c("a_g", "b", "c", "d_g", "e_g", "f_g", "g_g", "h_g", "i_g", "j_g", "k_g")) %>%
	mutate(a_g=as.numeric(a_g), b=as.numeric(b), c=as.numeric(c), d_g=as.numeric(d_g), e_g=as.numeric(e_g),
	f_g=as.numeric(f_g), g_g=as.numeric(g_g), h_g=as.numeric(h_g), i_g=as.numeric(i_g), j_g=as.numeric(j_g), k_g=as.numeric(k_g)) %>%
	mutate(sumofrow = rowSums(across(where(is.numeric)))) %>% pull(sumofrow)
sum(control_Gordon)

control_ANNSoLo <- All_unfiltered %>% select(CtrlCounts_ANNSoLo) %>%
	separate(CtrlCounts_ANNSoLo, into=c("a_a", "b_a", "c_a", "d_a", "e_a", "f_a", "g_a", "h_a", "i_a", "j_a", "k_a")) %>%
	mutate(a_a=as.numeric(a_a), b_a=as.numeric(b_a), c_a=as.numeric(c_a), d_a=as.numeric(d_a), e_a=as.numeric(e_a),
	f_a=as.numeric(f_a), g_a=as.numeric(g_a), h_a=as.numeric(h_a), i_a=as.numeric(i_a), j_a=as.numeric(j_a), k_a=as.numeric(k_a)) %>%
	mutate(sumofrow = rowSums(across(a_a:k_a))) %>% pull(sumofrow)
sum(control_ANNSoLo)

Count_ann <- All_unfiltered %>% 
	separate(CtrlCounts_ANNSoLo, into=c("a_a", "b_a", "c_a", "d_a", "e_a", "f_a", "g_a", "h_a", "i_a", "j_a", "k_a")) %>%
	mutate(a_a=as.numeric(a_a), b_a=as.numeric(b_a), c_a=as.numeric(c_a), d_a=as.numeric(d_a), e_a=as.numeric(e_a),
	f_a=as.numeric(f_a), g_a=as.numeric(g_a), h_a=as.numeric(h_a), i_a=as.numeric(i_a), j_a=as.numeric(j_a), k_a=as.numeric(k_a)) %>%
	mutate(AvgCntr = rowSums(across(a_a:k_a))/11) %>%
	rename(AvgSpec = AvgSpec_ANNSoLo) %>%
	rename(BFDR = BFDR_ANNSoLo) %>%
	rename(MiST = MiST_ANNSoLo) %>%
	separate(BP, into=c("bait", "prey"), remove=FALSE) %>% 
	select(BP, bait, AvgCntr, AvgSpec, BFDR, MiST) %>% mutate(Study="Reanalysis with ANN-SoLo")

Count_gor <- All_unfiltered %>%
	separate(CtrlCounts_Gordon, into=c("a_g", "b", "c", "d_g", "e_g", "f_g", "g_g", "h_g", "i_g", "j_g", "k_g")) %>%
	mutate(a_g=as.numeric(a_g), b=as.numeric(b), c=as.numeric(c), d_g=as.numeric(d_g), e_g=as.numeric(e_g),
	f_g=as.numeric(f_g), g_g=as.numeric(g_g), h_g=as.numeric(h_g), i_g=as.numeric(i_g), j_g=as.numeric(j_g), k_g=as.numeric(k_g)) %>%
	mutate(AvgCntr = rowSums(across(a_g:k_g))/11) %>%
	rename(AvgSpec = AvgSpec_Gordon) %>% 
	rename(BFDR = BFDR_Gordon) %>%
	rename(MiST = MiST_Gordon) %>%
	separate(BP, into=c("bait", "prey"), remove=FALSE) %>% 
	select(BP, bait, AvgCntr, AvgSpec, BFDR, MiST) %>% mutate(Study="Gordon et al.")

Counts_unfiltered <- rbind(Count_ann, Count_gor)

p <- ggplot(Counts_unfiltered, aes(x=Study, y=BFDR)) + 
	geom_boxplot()
p + geom_jitter(shape=16, position=position_jitter(0.2))

ggplot(Counts_unfiltered, aes(x=bait, y=BFDR, fill=Study)) +
	geom_boxplot()

#	T-tests

wilcox.test(BFDR ~ Study, data=Counts_unfiltered) 
wilcox.test(MiST ~ Study, data=Counts_unfiltered) 
wilcox.test(AvgCntr ~ Study, data=Counts_unfiltered) 
wilcox.test(AvgSpec ~ Study, data=Counts_unfiltered) 

ggplot(Counts_unfiltered, aes(x=AvgCntr, colour=Study)) +
	geom_density()

#	Correlation test

shapiro.test(All_unfiltered$BFDR_Gordon)
cor.test(All_unfiltered$MiST_Gordon, All_unfiltered$MiST_ANNSoLo, method="spearman")

################################################################################################################################################################################################################################################

#	create a venn diagram showing why HCIPs were missed by Gordon et al.
Venn_ann <- Plot_input_label_ann %>% select(BP_gene, AvgSpec_Gordon, BFDR_Gordon, MiST_Gordon) %>%
	mutate(MiST=if_else(!MiST_Gordon>=0.7, TRUE, FALSE)) %>%
	mutate(BFDR=if_else(!BFDR_Gordon<=0.05, TRUE, FALSE)) %>%
	mutate(AvgSpec=if_else(!AvgSpec_Gordon>=2, TRUE, FALSE))

venn <- ggvenn(Venn_ann, columns = c("MiST", "BFDR", "AvgSpec"),
	show_percentage=FALSE,
	fill_color=c("white","white","white"))
ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Filters.pdf", width = 17, height = 17, units = "cm")

#	create a venn diagram showing why HCIPs were missed by Gordon et al.
Venn_ann <- Plot_input_label_ann %>% select(BP_gene, AvgSpec_Gordon, BFDR_Gordon, MiST_Gordon) %>%
	mutate(MiST=if_else(!MiST_Gordon>=0.7, TRUE, FALSE)) %>%
	mutate(BFDR=if_else(!BFDR_Gordon<=0.05, TRUE, FALSE)) %>%
	mutate(AvgSpec=if_else(!AvgSpec_Gordon>=2, TRUE, FALSE))

venn_a <- ggvenn(Venn_ann, columns = c("MiST", "BFDR", "AvgSpec"),
	show_percentage=FALSE,
	fill_color=c("white","white","white"))
ggsave(venn, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Filters.pdf", width = 17, height = 17, units = "cm")

Venn_gor <- Plot_input_label_gor %>% select(BP_gene, AvgSpec_ANNSoLo, BFDR_ANNSoLo, MiST_ANNSoLo) %>%
	mutate(MiST=if_else(!MiST_ANNSoLo>=0.7, TRUE, FALSE)) %>%
	mutate(BFDR=if_else(!BFDR_ANNSoLo<=0.05, TRUE, FALSE)) %>%
	mutate(AvgSpec=if_else(!AvgSpec_ANNSoLo>=2, TRUE, FALSE))
Venn_gor %>% filter(MiST==TRUE) %>% filter(BP_gene %in% mis_gor$BP_gene)
Venn_gor2 <- Plot_input_label_gor %>% select(BP_gene, AvgSpec_Gordon, BFDR_Gordon, MiST_Gordon) %>%
	mutate(MiST=if_else(!MiST_Gordon>=0.7, TRUE, FALSE)) %>%
	mutate(BFDR=if_else(!BFDR_Gordon<=0.05, TRUE, FALSE)) %>%
	mutate(AvgSpec=if_else(!AvgSpec_Gordon>=2, TRUE, FALSE))
mis_gor <- Venn_gor2 %>% filter(MiST==TRUE)
Plot_input_label_gor %>% select(BP_gene, CtrlCounts_ANNSoLo, CtrlCounts_Gordon) %>% print(n=40)
mis_bfdr <- Venn_gor %>% filter(MiST==FALSE & BFDR==TRUE)
Plot_input_label_gor %>% select(BP_gene, AvgSpec_ANNSoLo, AvgSpec_Gordon, BFDR_ANNSoLo, BFDR_Gordon) %>% filter(BP_gene %in% mis_bfdr$BP_gene)



venn_g <- ggvenn(Venn_gor, columns = c("MiST", "BFDR", "AvgSpec"),
	show_percentage=FALSE,
	fill_color=c("white","white","white"))
ggsave(venn_g, file="/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Venn/Filters_Gordon_scores.pdf", width = 17, height = 17, units = "cm")

#	Interactions that were missed because the preys were not identified
missed_int <- comb %>% filter(!BP_gene %in%Plot_input_label$BP_gene) %>% select(BP_gene)
ANNSoLo_add_input %>% filter(BP_gene %in%missed_int$BP_gene) %>% pull(BP_gene)

#	Look at specific missed interactions
Plot_input_label_ann %>% filter(!MiST_Gordon>=0.7 & !BFDR_Gordon<=0.05 & AvgSpec_Gordon>=2) %>%
	select(BP_gene, CtrlCounts_ANNSoLo, CtrlCounts_Gordon, Spec_ANNSoLo, Spec_Gordon)


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


################################################################################################################################################

SAINT_MIST_filtered %>% filter(str_detect(PreyGene, "NDUF"))


