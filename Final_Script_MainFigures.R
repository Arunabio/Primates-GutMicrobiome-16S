library(vegan)
library (ape)
library (ade4)
library (factoextra)
library (cluster)

setwd("~/Work/AmpliconData_Primates/Analysis/mSphereRevision/Revised-Ana")
otu_table <- read.csv (file = "Genus_Rel_abun.txt", sep = "\t", row.names = 1, header = TRUE)
Groups <- otu_table$Group
Groups <- data.frame (Groups)

#-- Figure 1
bray_pcoa <-pcoa(vegdist(otu_table[,-1], "bray"))
bray_pcoa$values[1:2,]
mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, Groups)

Colors <- c("forestgreen","darkgoldenrod3","chartreuse1","purple3","darkmagenta","blue","darkcyan","deeppink3","black","darkorchid1","darkred","burlywood4","Brown1")

jpeg("Bray_PCoA.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, shape=Groups, colour=Groups)) + geom_point(size=1.5) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCo1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCo2 - ", mds.var.per[2], "%", sep="")) + ggtitle("Bray-Curtis Distances") + theme(axis.text.x = element_text(size = 16, colour = "black", face = "bold"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black", face = "bold"), legend.title = element_text(size = 18, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

pcoa1 <- Bray_PCoA_MATRIX_New[,c(1,3)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Groups")
jpeg("Bray_PCoA1_Distances.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(data = pcoa1_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo1") + theme_bw() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 7, colour = "black")) + coord_flip() + theme(legend.position='none')
dev.off ()

pcoa2 <- Bray_PCoA_MATRIX_New[,c(2,3)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Groups")
jpeg("Bray_PCoA2_Distances.jpg", height = 5, width = 4, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_bw() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 7, colour = "black"), axis.text.y = element_text(size = 10, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

#--- Significance Check
Stats_pcoa1 <- kruskalmc(pcoa1$Axis.1~pcoa1$Groups, probs=0.05)
Stats_pcoa1_p0.05 <- Stats_pcoa1$dif.com
write.table(Stats_pcoa1_p0.05, file = "Stats_pcoa1_p0.05.txt", sep= "\t")

Stats_pcoa2 <- kruskalmc(pcoa2$Axis.2~pcoa2$Groups, probs=0.05)
Stats_pcoa2_p0.05 <- Stats_pcoa2$dif.com
write.table(Stats_pcoa2_p0.05, file = "Stats_pcoa2_p0.05.txt", sep= "\t")

dist_bray <- vegdist(otu_table[,-1], "bray")
adonis(dist_bray ~ Group, data=otu_table, permutations=999)

#--Figure 2
dist_bray <- as.matrix (dist_bray)
jpeg("NumberOfOptimumClusters.jpg", height = 4, width = 5, units = 'in', res = 600)
fviz_nbclust(dist_bray, pam, method = "wss") + geom_vline(xintercept = 3, linetype = 2)
dev.off ()

pam <- pam(dist_bray, k=3)
clusters <- pam$clustering
clusters <- data.frame(clusters)
Bray_PCoA_MATRIX_Cluster <- cbind(Bray_PCoA_MATRIX, clusters)
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "1"] <- "Cluster1"
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "2"] <- "Cluster2"
Bray_PCoA_MATRIX_Cluster$clusters[Bray_PCoA_MATRIX_Cluster$clusters == "3"] <- "Cluster3"

jpeg("Bray_Pam_MetabotypesSplot.jpg", height = 4, width = 5, units = 'in', res = 600)
s.class(bray_pcoa$vectors, fac= as.factor(Bray_PCoA_MATRIX_Cluster$clusters), col = c("forestgreen", "darkgray", "red"), label = c("Cluster1", "Cluster2", "Cluster3"))
dev.off ()

write.table (Bray_PCoA_MATRIX_Cluster, file ="Bray_PCoA_MATRIX_Cluster_Group.txt", sep = "\t")

#-  Identified cluster membering of each primate species and graph plotted in Excel

#--- Figure 3 ---- Discriminating taxa were Identified using three approaches
#-- RandomForest
Bray_PCoA_MATRIX_Cluster <- read.csv(file = "Bray_PCoA_MATRIX_Cluster_Group.txt", sep = "\t", row.names = 1, header = T)
library (randomForest)
Clusters <- Bray_PCoA_MATRIX_Cluster$clusters
Clusters <- data.frame(Clusters)
Otu_table_withOutGroup <- otu_table[,-1]
Otu_table_Clusters <- cbind (Clusters, Otu_table_withOutGroup)

RF1 <- randomForest(Clusters ~., data=Otu_table_Clusters, ntree = 1000, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
jpeg("RF_MDA.jpg", height = 12, width = 15, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
MDA <- importance(RF1, type=1)
write.table (MDA, file = "MDA.txt", sep = "\t")

#---- Labdsv
library (labdsv)
Genus2 <- Otu_table_Clusters[ , which(!apply(Otu_table_Clusters==0,2,all))]
Genus2[is.na(Genus2)] <- 0
iva <- indval(Genus2[,2:525], Genus2$Clusters)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(Genus2[,2:525]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="indvalsummary.txt", sep = "\t")

#-- Krukshal-Wallis p-values

#--- Selected bacteria were picked from the otu table directly for plotting
#--Cluster1
grep("Collinsella", colnames(Otu_table_Clusters) )
Collinsella <- Otu_table_Clusters[,102]
Collinsella <- data.frame(Collinsella)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX, Collinsella)
grep("YS2", colnames(Otu_table_Clusters) )
YS2 <- Otu_table_Clusters[,189]
YS2 <- data.frame(YS2)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, YS2)
grep("Faecalibacterium", colnames(Otu_table_Clusters) )
Faecalibacterium <- Otu_table_Clusters[,304]
Faecalibacterium <- data.frame(Faecalibacterium)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX, Faecalibacterium)
grep("Coprococcus", colnames(Otu_table_Clusters) )
Coprococcus <- Otu_table_Clusters[,278]
Coprococcus <- data.frame(Coprococcus)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Coprococcus)
grep("Lachnospira", colnames(Otu_table_Clusters))
grep("g__Lachnospira", colnames(Otu_table_Clusters))
Lachnospira <- Otu_table_Clusters[,282]
Lachnospira <- data.frame(Lachnospira)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Lachnospira)
grep("g_Clostridium", colnames(Otu_table_Clusters) )
Clostridium <- Otu_table_Clusters[,260]
Clostridium <- data.frame(Clostridium)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Clostridium)
grep("RF32", colnames(Otu_table_Clusters) )
RF32 <- Otu_table_Clusters[,381]
RF32 <- data.frame(RF32)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, RF32)
grep("Victivallaceae", colnames(Otu_table_Clusters) )
Victivallaceae <- Otu_table_Clusters[,361]
Victivallaceae <- data.frame(Victivallaceae)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Victivallaceae)
grep("Prevotella", colnames(Otu_table_Clusters) )
Prevotella <- Otu_table_Clusters[,c(132,145)]
Prevotella <- rowSums(Prevotella)
Prevotella <- data.frame(Prevotella)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Prevotella)
Cluster1 <- rowSums(Bray_PCoA_MATRIX_Abundances[,3:9])
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Cluster1)

jpeg("Cluster1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Cluster1)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("Taxa_Cluster1_Sqrt_Per.jpg", height = 6, width = 6, units = 'in', res = 600)
boxplot(sqrt(Bray_PCoA_MATRIX_Abundances[c(8,7,5,3,6,4,9)]),data=Bray_PCoA_MATRIX_Abundances, main="Discriminating Taxa for Cluster 1", xlab="", ylab="sqrt (Total Abundance)", col = "lightgray", las=3, cex.axis=0.6)
dev.off ()

#--- Cluster2
grep("Sphaerochaeta", colnames(Otu_table_Clusters) )
Sphaerochaeta <- Otu_table_Clusters[,576]
Sphaerochaeta <- data.frame(Sphaerochaeta)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Sphaerochaeta)
grep("RFN20", colnames(Otu_table_Clusters) )
RFN20 <- Otu_table_Clusters[,350]
RFN20 <- data.frame(RFN20)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, RFN20)
grep("SHD.231", colnames(Otu_table_Clusters) )
SHD.231 <- Otu_table_Clusters[,178]
SHD.231 <- data.frame(SHD.231)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, SHD.231)
grep("Adlercreutzia", colnames(Otu_table_Clusters) )
Adlercreutzia <- Otu_table_Clusters[,100]
Adlercreutzia <- data.frame(Adlercreutzia)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Adlercreutzia)
grep("Mogibacterium", colnames(Otu_table_Clusters) )
Mogibacterium <- Otu_table_Clusters[,325]
Mogibacterium <- data.frame(Mogibacterium)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Mogibacterium)
grep("F16", colnames(Otu_table_Clusters) )
F16 <- Otu_table_Clusters[c(133,590)]
View(F16)
F16 <- Otu_table_Clusters[,590]
F16 <- data.frame(F16)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, F16)
grep("Treponema", colnames(Otu_table_Clusters) )
Treponema <- Otu_table_Clusters[,578]
Treponema <- data.frame(Treponema)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Treponema)
grep("Coriobacteriaceae.g__", colnames(Otu_table_Clusters) )
Coriobacteriaceae <- Otu_table_Clusters[c(99,100,101,102,103,104)]
View(Coriobacteriaceae)
Coriobacteriaceae <- Otu_table_Clusters[,99]
Coriobacteriaceae <- data.frame(Coriobacteriaceae)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Coriobacteriaceae)
grep("Unassigned.Other.Other.Other.Other.Other", colnames(Otu_table_Clusters) )
Unassigned <- Otu_table_Clusters[,2]
Unassigned <- data.frame(Unassigned)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Unassigned)
grep("Bulleidia", colnames(Otu_table_Clusters) )
Bulleidia <- Otu_table_Clusters[,342]
Bulleidia <- data.frame (Bulleidia)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Bulleidia)
grep("g__p.75.a5", colnames(Otu_table_Clusters) )
grep("p.75.a5", colnames(Otu_table_Clusters) )
p.75.a5 <- Otu_table_Clusters[,353]
p.75.a5 <- data.frame(p.75.a5)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, p.75.a5)
grep("Fibrobacter", colnames(Otu_table_Clusters) )
grep("g__Fibrobacter", colnames(Otu_table_Clusters) )
Fibrobacter <- Otu_table_Clusters[c(208,209)]
View(Fibrobacter)
Fibrobacter <- Otu_table_Clusters[,208]
Fibrobacter <- data.frame (Fibrobacter)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Fibrobacter)
grep("g__Butyrivibrio", colnames(Otu_table_Clusters) )
Butyrivibrio <- Otu_table_Clusters[,277]
Butyrivibrio <- data.frame (Butyrivibrio)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Butyrivibrio)

Cluster2 <- rowSums(Bray_PCoA_MATRIX_Abundances[,11:23])
Cluster2 <- data.frame(Cluster2)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Cluster2)

jpeg("Cluster2.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Cluster2)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("Taxa_Cluster2_Sqrt_Per.jpg", height = 6, width = 6, units = 'in', res = 600)
boxplot(sqrt(Bray_PCoA_MATRIX_Abundances[c(22,16,14,15,23,20,21,12,17,11,13,18,19)]),data=Bray_PCoA_MATRIX_Abundances, main="Discriminating Taxa for Cluster 2", xlab="", ylab="sqrt (Total Abundance)", col = "lightgray", las=3, cex.axis=0.6)
dev.off ()

#--Cluster3
grep("Bacteroides", colnames(Otu_table_Clusters) )
Bacteroides <- Otu_table_Clusters[,123]
Bacteroides <- data.frame(Bacteroides)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Bacteroides)
jpeg("Bacteroides.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Bacteroides)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
grep("Parabacteroides", colnames(Otu_table_Clusters) )
Parabacteroides <- Otu_table_Clusters[,128]
Parabacteroides <- data.frame(Parabacteroides)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Parabacteroides)
jpeg("Paracteroides.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Parabacteroides)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off()
grep("Sutterella", colnames(Otu_table_Clusters) )
Sutterella <- Otu_table_Clusters[,461]
Sutterella <- data.frame(Sutterella)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Sutterella)
jpeg("3_Sutterella.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Sutterella)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
grep("Streptococcus", colnames(Otu_table_Clusters) )
Streptococcus <- Otu_table_Clusters[,248]
Streptococcus <- data.frame(Streptococcus)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Streptococcus)
grep("Odoribacter", colnames(Otu_table_Clusters) )
grep("Odoribacter", colnames(Otu_table_Clusters) )
Odoribacter <- Otu_table_Clusters[,139]
Odoribacter <- data.frame(Odoribacter)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Odoribacter)

Cluster3 <- rowSums(Bray_PCoA_MATRIX_Abundances[,25:29])
Cluster3 <- data.frame(Cluster3)
Bray_PCoA_MATRIX_Abundances <- cbind(Bray_PCoA_MATRIX_Abundances, Cluster3)
jpeg("Cluster3.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_Abundances, aes(x=Axis.1, y=Axis.2, colour=Cluster3)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("Taxa_Cluster3_Sqrt_Per.jpg", height = 6, width = 6, units = 'in', res = 600)
boxplot(sqrt(Bray_PCoA_MATRIX_Abundances[c(27,25,26,28)]),data=Bray_PCoA_MATRIX_Abundances, main="Discriminating Taxa for Cluster 3", xlab="", ylab="sqrt (Total Abundance)", col = "lightgray", las=3, cex.axis=0.6)
dev.off ()

#---- Figure 4
library(phyloseq)
biom_file = "otu_table_mc2_w_tax_no_pynast_failures.biom" #File needed
biomot = import_biom(biom_file, parseFunction = parse_taxonomy_greengenes)
map_file = "Mapping.txt" 
bmsd = import_qiime_sample_data(map_file)
Bushman = merge_phyloseq(biomot, bmsd)
tree_file <- "rep_set.tre"
tree1 <- read_tree(tree_file)
Bushman1 <-merge_phyloseq(Bushman, tree1)
colSums(Bushman1@otu_table)

Bushman1.prune = prune_samples(sample_sums(Bushman1) > 1000, Bushman1)
set.seed(1234)
Bushman1.rare = rarefy_even_depth(Bushman1.prune, sample.size = 1000)
Bushman1.even.1000 = prune_taxa(taxa_sums(Bushman1.rare) > 0, Bushman1.rare)
#--eight samples having less than 1k depth were removed
#12.13, 12.169, 12.63A, HMP34, Vervet47	HMP11, Vervet55, HMP13

p <- plot_richness(Bushman1.even.1000, "Group", measures = c("Observed", "Chao1", "Shannon"))
p <- p + geom_boxplot(data = p$data, aes(x= Group, y = value, color = Group), alpha = 0.1) + scale_color_manual(values=Colors) + labs(x="",y="Alpha Diversity Measure") + theme(axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave("phyloseq_analysis-richness_estimates_Rarefied1k.pdf", p, width = 14, height = 7)
jpeg("Phyloseq_Richness_Rarefied1k.jpg", height = 7, width = 14, units = 'in', res = 600)
dev.off ()

erich_at1kdepth <- estimate_richness(Bushman1.even.1000, measures = c("Observed", "Chao1", "Shannon"))
write.table (erich_at1kdepth, file = "Richness_at1kdepth.txt", sep = "\t")
#-- Formatted by adding metadata information

richness_at1kdepth <- read.csv(file ="Richness_at1kdepth.txt.txt", sep = "\t", row.names = 1, header =T)
Observed_at1k <- richness_at1kdepth[,c(1,5)]
Observed_at1k_Melted <- melt(Observed_at1k, id.vars = "Group")
jpeg("Observed_at1k.png", height = 6, width = 6, units = 'in', res = 600)
#ggplot(Observed_at1k, aes(x=reorder(Group, Observed, FUN=median), y=Observed, color=Group, fill = Group), alpha = 0.1) + geom_boxplot() + ggtitle("Observed") + labs(x="",y="Alpha Diversity Measure") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
ggplot(data = Observed_at1k_Melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + ggtitle("Observed") + labs(x="",y="Alpha Diversity Measure") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
Stats_observed <- kruskalmc(richness_at1kdepth$Observed~richness_at1kdepth$Group, probs=0.05)
Stats_observed_p0.05 <- Stats_observed$dif.com

#--- Custer abundance and Observed OTUs was combined (file name: alpha_cluster) and used to do correlations 
jpeg("ObservedOtu_cluster1.jpg", height = 6, width = 6, units = 'in', res = 600)
plot (alpha_cluster$cluster1*100, alpha_cluster$Observed, col=c("forestgreen","gray", "red")[ alpha_cluster$Group], pch=19, xlab="%Abundance (Cluster1)",ylab="Observed_OTU")
cluster1_per <- alpha_cluster$cluster1*100
model1<-lm (alpha_cluster$Observed~ cluster1_per)
abline(model1)
dev.off ()
jpeg("ObservedOtu_cluster2.jpg", height = 6, width = 6, units = 'in', res = 600)
plot (alpha_cluster$cluster2*100, alpha_cluster$Observed, col=c("forestgreen","gray", "red")[ alpha_cluster$Group], pch=19, xlab="%Abundance (Cluster2)",ylab="Observed_OTU")
cluster2_per <- alpha_cluster$cluster2*100
model2<-lm (alpha_cluster$Observed~ cluster2_per)
abline(model2)
dev.off ()
jpeg("ObservedOtu_cluster3.jpg", height = 6, width = 6, units = 'in', res = 600)
plot (alpha_cluster$cluster3*100, alpha_cluster$Observed, col=c("forestgreen","gray", "red")[ alpha_cluster$Group], pch=19, xlab="%Abundance (Cluster3)",ylab="Observed_OTU")
cluster3_per <- alpha_cluster$cluster3*100
model3<-lm (alpha_cluster$Observed~ cluster3_per)
abline(model3)
dev.off ()