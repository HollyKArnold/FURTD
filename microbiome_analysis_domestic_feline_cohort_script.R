#AUTHOR: HOLLY ARNOLD
#DAY: September 25, 2022
#DATE: 20220925
#DESCRIPTION: 
# Provides script for analyses performed in paper chronic clinical signs of upper respiratory 
# microbiomes in a cohort of domestic felines

##############################################################################################################
#LOCATIONS
## 1. LOCAL
## /Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/originalData

## 2. SERVER
## /nfs3/Sharpton_Lab/prod/projects/arnoldhk/furtd2019/version2/2.dada2

#WRITE UP
#See 20201104PaperDraftVersion1.tex for further details on this section.
##############################################################################################################

#LIBRARIES
library(vegan)
library(ape)
library(xtable)
library(ggplot2); packageVersion("ggplot2")
library(gridExtra)
library(phyloseq); packageVersion("phyloseq")
library(ggtree)
library(randomForest)
library(stringr)
library(pROC)
library(reshape)
library(ggsignif)

#library(dada2)
#library(reshape2)
#library(pracma)
#library(devtools)
#library(ggbiplot)
#library(ggfortify)
#library(BiodiversityR)
#library(ggrepel)
#library(GUniFrac)
#library(dplyr)
#library(tidyr)
#library(adespatial)
#library(pROC)
#library(caret)

# DIRECTORIES
directory.data.original = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/originalData/"
directory.data.rarify = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/rarifiedData/"
directory.data.metadata = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/metaData/"
directory.data = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/"
directory.data.ctu.gut = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/cladalAnalysis/gut/"
directory.data.ctu.nasal = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/cladalAnalysis/nasal/"
directory.figures = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/figures/"
directory.scripts = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/scripts/"
directory.out = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/out/"

# LOAD FUNCTIONS
setwd(directory.scripts)
source("catMicrobiomeAnalysisVersion2Functions.R")

# SET COLOR PALLET
pal = RColorBrewer::display.brewer.all()
pal = RColorBrewer::brewer.pal(n = 7, "Dark2")

# LOAD DATA TABLES
## 1. Read in ASV tables
setwd(directory.data.rarify)
asvG = read.table("asvGutRarify.txt", header = T)
colnames(asvG) = sub(pattern = "X", replacement = "", x = colnames(asvG))
head(asvG)
dim(asvG)

asvN = read.table("asvNasalRarify.txt", header = T)
colnames(asvN) = sub(pattern = "X", replacement = "", x = colnames(asvN))
head(asvN)
dim(asvN)

## 2. Read in Tax Tables
taxG = read.table("asvGutTaxRarify.txt")
head(taxG)
dim(taxG)
taxN = read.table("asvNasalTaxRarify.txt")
head(taxN)
dim(taxN)

taxN.String = read.table("asvNasalTaxString.txt", sep = "\t", row.names = 1)
head(taxN.String)
taxG.String = read.table("asvGutTaxString.txt", sep = "\t", row.names = 1)
head(taxG.String)

## 3. Read in metadata tables
setwd(directory.out)
meta = read.table("metaData.txt", sep = "\t")
head(meta)
metaG = meta[colnames(asvG),]
head(metaG)  
rownames(meta) = meta$NasalID.1
metaN = meta[colnames(asvN),]
head(metaN)

# 4. Read in the cladal analysis files
setwd(directory.data.ctu.gut)
ctuG = t(read.table("ctuGut.txt", sep = "\t", header = T))
head(ctuG)
treeG = read.tree("new_prepped_tree.tre")
treeG

aG = getCladalAttributes(nodes2tax = "gutCladeStat_nodes2tax.txt", 
                         size = "gutCladeStat_clade_size.txt",
                         nodeTax = "taxCladesGut.txt",
                         groupTest = "gutGroupTestSignsVetReport.txt_stats.txt",
                         pTest = "gutPTest1000_stats.txt")
head(aG)

setwd(directory.data.ctu.nasal)
ctuN = t(read.table("ctuNasal.txt", sep = "\t", header = T))
head(ctuN)
treeN = read.tree("new_prepped_tree.tre")
treeN
aN = getCladalAttributes(nodes2tax = "nasalCladeStat_nodes2tax.txt",
                         size = "nasalCladeStat_clade_size.txt", 
                         nodeTax = "taxCladesNasal.txt",
                         groupTest = "nasalGroupTestSignsVetReport.txt_stats.txt",
                         pTest = "nasalPTest1000_stats.txt")

head(aN)

# 5. Read in PCR Results
setwd(directory.data)
pcr = read.table("metaData/PCRResults.txt", sep = "\t", header = T)
head(pcr)

##############################################################################################################
# 1. Get tables for owner reported signs
setwd(directory.data.metadata)
oCx = read.table("ownerReportedSigns.txt", sep = "\t", header = T)
head(oCx)

oCxP = oCx[which(oCx$OwnerReportedStatus == "P"),3:17]
oCxSumP = apply(oCxP, 2, sum)
oCxSumP = oCxSumP[order(oCxSumP, decreasing = T)]

oCxN = oCx[which(oCx$OwnerReportedStatus == "N"), 3:17]
oCxSumN = apply(oCxN, 2, sum)
oCxSumN = oCxSumN[names(oCxSumP)]

xtable(t(rbind(oCxSumP, oCxSumN)))

# 2. Get table for vet reported signs
vCx = read.table("vetReportedSigns.txt", sep = "\t", header = T)
head(vCx)

vCxP = vCx[which(vCx$OwnerReportedStatus == "P"), 3:15]
vCxP
vCxSumP = apply(vCxP, 2, sum)
vCxSumP = vCxSumP[order(vCxSumP, decreasing = T)]
vCxSumP

vCxN = vCx[which(vCx$OwnerReportedStatus == "N"), 3:15]
head(vCxN)
vCxSumN = apply(vCxN, 2, sum)
vCxSumN = vCxSumN[order(vCxSumN, decreasing = T)]
vCxSumN
xtable(t(rbind(vCxSumP, vCxSumN)))

rm(list = c("vCxN", "vCxP", "vCxSumN", "vCx", "vCxSumP", "oCxP", "oCx", "oCxSumP", "oCxN", "oCxSumN"))

##############################################################################################################

# Plot PCoAs

## 2A. Make a quantitative signs vector based on vet PE
metaG$sumSigns = metaG$Eyes + metaG$Ears + metaG$Nose + metaG$Oral + metaG$Ln + metaG$CV + metaG$LungsRespiratory + metaG$AbdominalGastrointestinal + metaG$UrogenitalPerineal + metaG$Musculoskeletal + metaG$Integument + metaG$Neurologic
metaG$sumSigns
metaN$sumSigns = metaN$Eyes + metaN$Ears + metaN$Nose + metaN$Oral + metaN$Ln + metaN$CV + metaN$LungsRespiratory + metaN$AbdominalGastrointestinal + metaN$UrogenitalPerineal + metaN$Musculoskeletal + metaN$Integument + metaN$Neurologic
metaN$sumSigns

## 2B. Distance matrix
bc.dist.gut.ab = vegdist(t(asvG), method = "bray", binary = F)
bc.dist.gut.pa = vegdist(t(asvG), method = "bray", binary = T)
bc.dist.nasal.ab = vegdist(t(asvN), method = "bray", binary = F)
bc.dist.nasal.pa = vegdist(t(asvN), method = "bray", binary = T)

pcoa.gut.ab = cmdscale(bc.dist.gut.ab, eig = T, k = 2)
pcoa.gut.pa = cmdscale(bc.dist.gut.pa, eig = T, k = 2)
pcoa.nasal.ab = cmdscale(bc.dist.nasal.ab, eig = T, k = 2)
pcoa.nasal.pa = cmdscale(bc.dist.nasal.pa, eig = T, k = 2)

pcoa.gut.ab.d = as.data.frame(pcoa.gut.ab$points)
colnames(pcoa.gut.ab.d) = c("PC1", "PC2")
pcoa.gut.ab.d$Household = metaG[rownames(pcoa.gut.ab.d), "Household"]
pcoa.gut.ab.d$SignsVetReported = metaG[rownames(pcoa.gut.ab.d), "SignsVetReported"]
pcoa.gut.ab.d$sumSigns = metaG[rownames(pcoa.gut.ab.d), "sumSigns"]

pcoa.gut.pa.d = as.data.frame(pcoa.gut.pa$points)
colnames(pcoa.gut.pa.d) = c("PC1", "PC2")
pcoa.gut.pa.d$Household = metaG[rownames(pcoa.gut.pa.d), "Household"]
pcoa.gut.pa.d$SignsVetReported = metaG[rownames(pcoa.gut.pa.d), "SignsVetReported"]
pcoa.gut.pa.d$sumSigns = metaG[rownames(pcoa.gut.pa.d), "sumSigns"]

pcoa.nasal.ab.d = as.data.frame(pcoa.nasal.ab$points)
colnames(pcoa.nasal.ab.d) = c("PC1", "PC2")
pcoa.nasal.ab.d$Household = metaN[rownames(pcoa.nasal.ab.d), "Household"]
pcoa.nasal.ab.d$SignsVetReported = metaN[rownames(pcoa.nasal.ab.d), "SignsVetReported"]
pcoa.nasal.ab.d$sumSigns = metaN[rownames(pcoa.nasal.ab.d), "sumSigns"]

pcoa.nasal.pa.d = as.data.frame(pcoa.nasal.pa$points)
colnames(pcoa.nasal.pa.d) = c("PC1", "PC2")
pcoa.nasal.pa.d$Household = metaN[rownames(pcoa.nasal.pa.d), "Household"]
pcoa.nasal.pa.d$SignsVetReported = metaN[rownames(pcoa.nasal.ab.d), "SignsVetReported"]
pcoa.nasal.pa.d$sumSigns = metaN[rownames(pcoa.nasal.pa.d), "sumSigns"]


## 2C. PCoAs by household and disease status
setwd(directory.figures)
#pdf("gutAbundancePCOAColoredByHouseholdAndStatus1.pdf", width = 14, height = 10)
ggplot(data = pcoa.gut.ab.d, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.gut.ab.d, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "Gut Abun BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("gutPAPCOAColoredByHouseholdAndStatus1.pdf", width = 14, height = 10)
ggplot(data = pcoa.gut.pa.d, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.gut.pa.d, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "Gut PA BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("nasalAbundancePCOAColoredByHouseholdAndStatus1.pdf", width = 14, height = 10)
ggplot(data = pcoa.nasal.ab.d, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.nasal.ab.d, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "Nasopharyngeal Abun BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("nasalPAPCOAColoredByHouseholdAndStatus1.pdf", width = 14, height = 10)
ggplot(data = pcoa.nasal.pa.d, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.nasal.pa.d, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "Nasopharyngeal PA BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status") 
#dev.off()

## 2C. Check that there aren't significant beta dispersion differences
### i. Gut Abundances by status
beta.bc.dist.gut.ab = betadisper(d = bc.dist.gut.ab, group = metaG[, "SignsVetReported"])
anova(beta.bc.dist.gut.ab) # Not significant 0.36
permutest(beta.bc.dist.gut.ab) #Not significant 0.33
labs = paste("Dimension", 1:4, "(", round(100*beta.bc.dist.gut.ab$eig/sum(beta.bc.dist.gut.ab$eig), 2), "% )")
setwd(directory.figures)

#pdf("betaDispersionGutAbundances68Conf.pdf")
par(mar = c(5, 5, 4, 2))
plot(beta.bc.dist.gut.ab, cex = 2, cex.lab = 2, cex.main = 2, pch = 15:16, xlab = labs[1], ylab = labs[2], hull = FALSE, ellipse = TRUE, conf = 0.68, lwd = 3, main = "Beta Dispersion \nBC Gut Abundance", col = c(pal[1], pal[3]) )
ordilabel(scores(beta.bc.dist.gut.ab, "centroids"), labels = c("Control", "FURTD"), col = c(pal[1], pal[3]), cex = 1.2)
#dev.off()

#pdf("betaDispersionGutAbundancesBoxplot.pdf")
par(mar = c(5, 5, 4, 2))
boxplot(beta.bc.dist.gut.ab, xlab = "Status", notch = F, cex.lab = 2, col = c(pal[1], pal[3]), main = "Beta Dispersion \nBC Gut Abundance", cex.main = 2, cex.axis = 2, names = c("Control", "FURTD"))
#dev.off()

### ii. Gut Presence absence by status
beta.bc.dist.gut.pa = betadisper(d = bc.dist.gut.pa, group = metaG[, "SignsVetReported"])
anova(beta.bc.dist.gut.pa) # Not significant 0.41
permutest(beta.bc.dist.gut.pa) # Not significant 0.5
labs = paste("Dimension", 1:4, "(", round(100*beta.bc.dist.gut.pa$eig/sum(beta.bc.dist.gut.pa$eig), 2), "% )")
setwd(directory.figures)

#pdf("betaDispersionGutPA68Conf.pdf")
par(mar = c(5, 5, 4, 2))
plot(beta.bc.dist.gut.pa, cex = 2, cex.main = 2, cex.lab = 2, pch = 15:16, xlab = labs[1], ylab = labs[2], hull = FALSE, ellipse = TRUE, conf = 0.68, lwd = 2, main = "Beta Dispersion\n BC Gut Presence Absence", col = c(pal[1], pal[3]))
ordilabel(scores(beta.bc.dist.gut.pa, "centroids"), labels = c("Control", "FURTD"), col = c(pal[1], pal[3]), cex = 1.2)
#dev.off()

#pdf("betaDispersionGutPABoxplot.pdf")
par(mar = c(5, 5, 4, 2))
boxplot(beta.bc.dist.gut.pa, xlab = "Status", notch = F, cex.lab = 2, col = c(pal[1], pal[3]), main = "Beta Dispersion \nBC Gut Presnece Absence", cex.main = 2, cex.axis = 2, names = c("Control", "FURTD"))
#dev.off()

### iii. Nasal Abundances
beta.bc.dist.nasal.ab = betadisper(d = bc.dist.nasal.ab, group = metaN[, "SignsVetReported"])
anova(beta.bc.dist.nasal.ab) # Not significant 0.19
permutest(beta.bc.dist.nasal.ab) # Not significant 0.2
labs = paste("Dimension", 1:4, "(", round(100*beta.bc.dist.nasal.ab$eig/sum(beta.bc.dist.nasal.ab$eig), 2), "% )")
setwd(directory.figures)

#pdf("betaDispersionNasalAbundances68Conf.pdf")
par(mar = c(5, 5, 4, 2))
plot(beta.bc.dist.nasal.ab, cex = 2, cex.lab = 2, cex.main = 2, pch = 15:16, xlab = labs[1], ylab = labs[2], hull = FALSE, ellipse = TRUE, conf = 0.68, lwd = 3, main = "Beta Dispersion \nBC Nasopharyngeal Abundance", col = c(pal[1], pal[3]))
ordilabel(scores(beta.bc.dist.nasal.ab, "centroids"), labels = c("Control", "FURTD"), col = c(pal[1], pal[3]), cex = 1.2)
#dev.off()

#pdf("betaDispersionNasalAbundancesBoxplot.pdf")
par(mar = c(5, 5, 4, 2))
boxplot(beta.bc.dist.nasal.ab, xlab = "Status", notch = F, cex.lab = 2, col = c(pal[1], pal[3]), main = "Beta Dispersion \nBC Nasopharyngeal Abundance", cex.main = 2, cex.axis = 2, names = c("Control", "FURTD"))
#dev.off()

### iv. Nasal Presence absence
beta.bc.dist.nasal.pa = betadisper(d = bc.dist.nasal.pa, group = metaN[, "SignsVetReported"])
anova(beta.bc.dist.nasal.pa) # Not significant 0.36
permutest(beta.bc.dist.nasal.pa) # Not significant 0.43
labs = paste("Dimension", 1:4, "(", round(100*beta.bc.dist.nasal.pa$eig/sum(beta.bc.dist.nasal.pa$eig), 2), "% )")

setwd(directory.figures)
#pdf("betaDispersionNasalPA68Conf.pdf")
par(mar = c(5, 5, 4, 2))
plot(beta.bc.dist.nasal.pa, cex = 2, cex.lab = 2, cex.main = 2, pch = 15:16, xlab = labs[1], ylab = labs[2], hull = FALSE, ellipse = TRUE, conf = 0.68, lwd = 2, main = "Beta Dispersion \n BC Nasopharyngeal Presence Absence", col = c(pal[1], pal[3]))
ordilabel(scores(beta.bc.dist.nasal.pa, "centroids"), labels = c("Control", "FURTD"), col = c(pal[1], pal[3]), cex = 1.2)
#dev.off()

#pdf("betaDispersionNasalPABoxplot.pdf")
par(mar = c(5, 5, 4, 2))
boxplot(beta.bc.dist.nasal.pa, xlab = "Status", notch = F, cex.lab = 2, col = c(pal[1], pal[3]), main = "Beta Dispersion BC \n Nasopharyngeal Presence Absence", cex.main = 2, cex.axis = 2, names = c("Control", "FURTD"))
#dev.off()

### v. Household
beta.bc.dist.nasal.ab = betadisper(d = bc.dist.nasal.ab, group = metaN[, "Household"])
anova(beta.bc.dist.nasal.ab) # Not significant 0.60
permutest(beta.bc.dist.nasal.ab) # Not significant 0.2


## 2D. Get significance between groups for bray curtis
dim(metaG)
adonis(formula = bc.dist.gut.ab ~ SignsVetReported, data = metaG) #p = 0.54
adonis(formula = bc.dist.gut.pa ~ SignsVetReported, data = metaG) #p = 0.39
adonis(formula = bc.dist.nasal.ab ~ SignsVetReported, data = metaN) #p = 0.92
adonis(formula = bc.dist.nasal.pa ~ SignsVetReported, data = metaN) # p = 0.99

adonis(formula = bc.dist.gut.ab ~ sumSigns, data = metaG) #p = 0.26
adonis(formula = bc.dist.gut.pa ~ sumSigns, data = metaG) #p = 0.15
adonis(formula = bc.dist.nasal.ab ~ sumSigns, data = metaN) #p = 0.6
adonis(formula = bc.dist.nasal.pa ~ sumSigns, data = metaN) # p = 0.95

adonis(formula = bc.dist.gut.ab ~ Household, data = metaG) #p = 0.14
adonis(formula = bc.dist.gut.pa ~ Household, data = metaG) #p = 0.008 **
adonis(formula = bc.dist.nasal.ab ~ Household, data = metaN) #p = 0.17
adonis(formula = bc.dist.nasal.pa ~ Household, data = metaN) #p = 0.22

adonis(bc.dist.gut.ab ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.03 *
adonis(bc.dist.gut.ab ~ SignsVetReported, data = metaG, strata = factor(metaG$Household)) #0.03 *


adonis(bc.dist.gut.pa ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.15
adonis(bc.dist.gut.pa ~ SignsVetReported, data = metaG, strata = factor(metaG$Household)) #0.15


adonis(bc.dist.nasal.ab ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household))  #0.81
adonis(bc.dist.nasal.ab ~ SignsVetReported, data = metaN, strata = factor(metaN$Household))  #0.81

adonis(bc.dist.nasal.pa ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household)) #0.88
adonis(bc.dist.nasal.pa ~ SignsVetReported , data = metaN, strata = factor(metaN$Household)) #0.88


adonis(bc.dist.gut.ab ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.02 ***
adonis(bc.dist.gut.ab ~ sumSigns, data = metaG, strata = factor(metaG$Household)) #0.02 ***

adonis(bc.dist.gut.pa ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.11
adonis(bc.dist.nasal.ab ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household))  #0.4
adonis(bc.dist.nasal.pa ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household)) #0.82


adonis(formula = bc.dist.gut.ab ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = bc.dist.gut.pa ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = bc.dist.nasal.ab ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)
adonis(formula = bc.dist.nasal.pa ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)

## 2E. Now, determine if Unifrac holds the same patterns.

## A. Make phyloseq objects
metaG$color = ifelse(metaG$SignsVetReported == "P", pal[3], pal[1])
otu.gut.physeq = otu_table(asvG, taxa_are_rows = T)
tax.gut.physeq = tax_table(as.matrix(taxG))
meta.gut.physeq = sample_data(metaG)
tree.gut.physeq = phy_tree(treeG)
psG = phyloseq(otu.gut.physeq, tax.gut.physeq, meta.gut.physeq, tree.gut.physeq)

metaN$color = ifelse(metaN$SignsVetReported == "P", pal[3], pal[1])
otu.nasal.physeq = otu_table(asvN, taxa_are_rows = T)
tax.nasal.physeq = tax_table(as.matrix(taxN))
meta.nasal.physeq = sample_data(metaN)
tree.nasal.physeq = phy_tree(treeN)
psN = phyloseq(otu.nasal.physeq, tax.nasal.physeq, meta.nasal.physeq, tree.nasal.physeq)

uni.dist.gut.ab = UniFrac(physeq = psG, weighted = TRUE, normalized = T, parallel = FALSE)
uni.dist.gut.pa = UniFrac(physeq = psG, weighted = FALSE, normalized = T, parallel = FALSE)
uni.dist.nasal.ab = UniFrac(physeq = psN, weighted = TRUE, normalized = T, parallel = FALSE)
uni.dist.nasal.pa = UniFrac(physeq = psN, weighted = FALSE, normalized= T, parallel = FALSE)

adonis(formula = uni.dist.gut.ab ~ SignsVetReported, data = metaG) #p = 0.52
adonis(formula = uni.dist.gut.pa ~ SignsVetReported, data = metaG) #p = 0.17
adonis(formula = uni.dist.nasal.ab ~ SignsVetReported, data = metaN) #p = 0.65
adonis(formula = uni.dist.nasal.pa ~ SignsVetReported, data = metaN) # p = 0.98

adonis(formula = uni.dist.gut.ab ~ sumSigns, data = metaG) #p = 0.41
adonis(formula = uni.dist.gut.pa ~ sumSigns, data = metaG) #p = 0.0.9
adonis(formula = uni.dist.nasal.ab ~ sumSigns, data = metaN) #p = 0.57
adonis(formula = uni.dist.nasal.pa ~ sumSigns, data = metaN) # p = 0.90

adonis(formula = uni.dist.gut.ab ~ Household, data = metaG) #p = 0.32
adonis(formula = uni.dist.gut.pa ~ Household, data = metaG) #p = 0.17
adonis(formula = uni.dist.nasal.ab ~ Household, data = metaN) #p = 0.37
adonis(formula = uni.dist.nasal.pa ~ Household, data = metaN) #p = 0.713

adonis(uni.dist.gut.ab ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.50
adonis(uni.dist.gut.ab ~ SignsVetReported, data = metaG, strata = factor(metaG$Household)) #0.5

adonis(uni.dist.gut.pa ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.15
adonis(uni.dist.gut.pa ~ SignsVetReported, data = metaG, strata = factor(metaG$Household)) #0.16

adonis(uni.dist.nasal.ab ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household))  #0.38
adonis(uni.dist.nasal.ab ~ SignsVetReported, data = metaN, strata = factor(metaN$Household))  #0.38

adonis(uni.dist.nasal.pa ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household)) #0.93
adonis(uni.dist.nasal.pa ~ SignsVetReported , data = metaN, strata = factor(metaN$Household)) #0.93


adonis(uni.dist.gut.ab ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.43
adonis(uni.dist.gut.ab ~ sumSigns, data = metaG, strata = factor(metaG$Household)) #0.40

adonis(uni.dist.gut.pa ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.06
adonis(uni.dist.nasal.ab ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household))  #0.46
adonis(uni.dist.nasal.pa ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household)) #0.66


adonis(formula = uni.dist.gut.ab ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = uni.dist.gut.pa ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = uni.dist.nasal.ab ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)
adonis(formula = uni.dist.nasal.pa ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)


### B. Test RDA of unifrac
colVector = factor((metaG$SignsVetReported))
col_palette = palette(value = c(pal[1], pal[3]))[colVector]
col_palette

dbrda.uni.dist.gut.ab = capscale(uni.dist.gut.ab ~ sumSigns, data = metaG)
plot(dbrda.uni.dist.gut.ab, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), tck = 0)
groups = metaG$SignsVetReported
household = metaG$Household
points(dbrda.uni.dist.gut.ab, col = as.numeric(factor(groups)), pch = 19, cex = 3)
ordispider(dbrda.uni.dist.gut.ab, groups, lty = 2, col = c(pal[1], pal[3]), label = F)
ordiellipse(dbrda.uni.dist.gut.ab, groups, lty = 2, col = c(pal[1], pal[3]), label = F, conf = 0.95)
text(y = dbrda.uni.dist.gut.ab$CCA$Xbar[,1], x = dbrda.uni.dist.gut.ab$CA$u[,1], labels = as.vector(metaG$Household))
legend('topright', legend=c("Control", "FURTD"), col=c(pal[1],pal[3]), pch = 16, pt.cex = 2.3, bty = "n")

anova.cca(dbrda.uni.dist.gut.ab, by = "terms") #0.39
anova.cca(dbrda.uni.dist.gut.ab, by = "axis") #0.4

dbrda.uni.dist.gut.pa = capscale(uni.dist.gut.pa ~ sumSigns + Condition(Household), data = metaG)
plot(dbrda.uni.dist.gut.pa)
plot(dbrda.uni.dist.gut.pa, type = "n")
groups = metaG$SignsVetReported
points(dbrda.uni.dist.gut.ab, col = as.numeric(as.factor(groups)), pch = as.numeric(as.factor(groups)))
ordispider(dbrda.uni.dist.gut.ab, groups, lty = 2, col = "grey", label = F)
ordiellipse(dbrda.uni.dist.gut.ab, groups, lty = 2, col = "grey", label = F)
anova.cca(dbrda.uni.dist.gut.ab, by = "terms")



plot(gut.ab.cap, type = "n")

points(gut.ab.cap, col = as.numeric(as.factor(groups)),
       pch = as.numeric(as.factor(groups)))

ordispider(gut.ab.cap, groups, lty=2, col="grey", label=F)
ordiellipse(gut.ab.cap, groups, lty=2, col="grey", label=F)



gut.ab = decostand(t(asvG), 'hell')
cca.gut.ab = cca(gut.ab ~ sumSigns, data = metaG)
cca.gut.ab
anova(cca.gut.ab)

anova(cca.gut.ab, by = 'axis', step = 1000)
anova(cca.gut.ab, by = 'term', step = 1000)
ordistep(cca(gut.ab ~1, data = metaG), scope = formula(cca.gut.ab), direction = "forward", pstep = 1000)
anova(cca(gut.ab[rownames(metaG),]~ sumSigns + Household + NumberNeutrophils, data = metaG, na.action = na.omit))

cca.gut.uni.ab = cca(uni.dist.gut.ab~sumSigns, data = metaG)

gut.ab.cap = capscale(t(asvG) ~ SignsVetReported + Household + OwnerReportedStatus, data = metaG, dist = "bray")

modNull.gut.ab.cap = capscale(t(asvG) ~ 1, data = metaG, distance = "bray")
anova(gut.ab.cap)
anova(gut.ab.cap, by = "axis", step = 1000)
anova(gut.ab.cap, by = "terms", step = 1000)


vasc <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/vasc_plants.txt', row.names = 1)
chem <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/chemistry.txt', row.names = 1)

# transform data using Hellinger transformation to prepare them for tb-RDA:
library (vegan)
vasc.hell <- decostand (vasc, 'hell')

# the last variable in the 'chem' dataset is 'slope', which is not a chemical variable - remove it:
chem1 <- chem[,-15]
head(vasc.hell)
head(chem1)
         
step = ordistep(gut.ab.cap, scope = modNull.gut.ab.cap, direction = "forward", permutations = 999)
vegan::anova.cca(step, by = "terms")
anova(step)
summary(step)
anova(step)
ordistep(modNull.gut.ab.cap, scope = formula(gut.ab.cap), perm.max = 200)
ordistep(gut.ab.cap, perm.max = 200)

asv.gut = t(asvG)
gut.ab.cap = capscale(t(asvG) ~ sumSigns, data = metaG, dist = "bray")
anova(gut.ab.cap)
anova(gut.ab.cap, by = "axis", step = 1000)
anova(gut.ab.cap, by = "terms", step = 1000)

nasal.ab.cap = capscale(t(asvN) ~ sumSigns, data = metaN, dist = "bray")
anova(nasal.ab.cap)
anova(nasal.ab.cap, by = "axis", step = 1000)
anova(nasal.ab.cap, by = "terms", step = 1000)

groups = metaG$SignsVetReported
plot(gut.ab.cap, type = "n")
points(gut.ab.cap, col = as.numeric(as.factor(groups)),
       pch = as.numeric(as.factor(groups)))

ordispider(gut.ab.cap, groups, lty=2, col="grey", label=F)
ordiellipse(gut.ab.cap, groups, lty=2, col="grey", label=F)

anova(gut.ab.cap, step = 1000, by = "terms")

###################################
#3. Redo analysis but with clades

# Plot Cladal PCoAs 
## A. Distance matrix clades
bc.dist.gut.ab.ctu = vegdist(t(ctuG), method = "bray", binary = F)
bc.dist.gut.pa.ctu = vegdist(t(ctuG), method = "bray", binary = T)
bc.dist.nasal.ab.ctu = vegdist(t(ctuN), method = "bray", binary = F)
bc.dist.nasal.pa.ctu = vegdist(t(ctuN), method = "bray", binary = T)

pcoa.gut.ab.ctu = cmdscale(bc.dist.gut.ab.ctu, eig = T, k = 2)
pcoa.gut.pa.ctu = cmdscale(bc.dist.gut.pa.ctu, eig = T, k = 2)
pcoa.nasal.ab.ctu = cmdscale(bc.dist.nasal.ab.ctu, eig = T, k = 2)
pcoa.nasal.pa.ctu = cmdscale(bc.dist.nasal.pa.ctu, eig = T, k = 2)

pcoa.gut.ab.d.ctu = as.data.frame(pcoa.gut.ab.ctu$points)
colnames(pcoa.gut.ab.d.ctu) = c("PC1", "PC2")
pcoa.gut.ab.d.ctu$Household = metaG[rownames(pcoa.gut.ab.d.ctu), "Household"]
pcoa.gut.ab.d.ctu$SignsVetReported = metaG[rownames(pcoa.gut.ab.d.ctu), "SignsVetReported"]
pcoa.gut.pa.d.ctu = as.data.frame(pcoa.gut.pa.ctu$points)
colnames(pcoa.gut.pa.d.ctu) = c("PC1", "PC2")
pcoa.gut.pa.d.ctu$Household = metaG[rownames(pcoa.gut.pa.d.ctu), "Household"]
pcoa.gut.pa.d.ctu$SignsVetReported = metaG[rownames(pcoa.gut.pa.d.ctu), "SignsVetReported"]
pcoa.nasal.ab.d.ctu = as.data.frame(pcoa.nasal.ab.ctu$points)
colnames(pcoa.nasal.ab.d.ctu) = c("PC1", "PC2")
pcoa.nasal.ab.d.ctu$Household = metaN[rownames(pcoa.nasal.ab.d.ctu), "Household"]
pcoa.nasal.ab.d.ctu$SignsVetReported = metaN[rownames(pcoa.nasal.ab.d.ctu), "SignsVetReported"]
pcoa.nasal.pa.d.ctu = as.data.frame(pcoa.nasal.pa.ctu$points)
colnames(pcoa.nasal.pa.d.ctu) = c("PC1", "PC2")
pcoa.nasal.pa.d.ctu$Household = metaN[rownames(pcoa.nasal.pa.d.ctu), "Household"]
pcoa.nasal.pa.d.ctu$SignsVetReported = metaN[rownames(pcoa.nasal.ab.d.ctu), "SignsVetReported"]


## 2B. PCoAs by household and disease status
setwd(directory.figures)
#pdf("gutAbundancePCOAColoredByHouseholdAndStatusCTU1.pdf", width = 14, height = 10)
ggplot(data = pcoa.gut.ab.d.ctu, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.gut.ab.d.ctu, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "CTU Gut Abun BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("gutPAPCOAColoredByHouseholdAndStatusCTU1.pdf", width = 14, height = 10)
ggplot(data = pcoa.gut.pa.d.ctu, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.gut.pa.d.ctu, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "CTU Gut PA BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("nasalAbundancePCOAColoredByHouseholdAndStatusCTU1.pdf", width = 14, height = 10)
ggplot(data = pcoa.nasal.ab.d.ctu, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.nasal.ab.d.ctu, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "CTU Nasopharyngeal Abun BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()

#pdf("nasalPAPCOAColoredByHouseholdAndStatusCTU1.pdf", width = 14, height = 10)
ggplot(data = pcoa.nasal.pa.d.ctu, aes(x = PC1, y = PC2)) +
  geom_point(data = pcoa.nasal.pa.d.ctu, aes(color = as.factor(SignsVetReported)), size = 30, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 20) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  labs(title = "CTU Nasopharyngeal PA BC", x = "PC1", y = "PC2") +
  theme(axis.text = element_text(size = 50), plot.title = element_text(size = 50), axis.title = element_text(size = 50), legend.text=element_text(size=50), legend.title = element_text(size = 50)) +
  labs(color = "Status")
#dev.off()


## C. Get significance between groups for bray curtis

adonis(formula = bc.dist.gut.ab.ctu ~ SignsVetReported, data = metaG) #p = 0.80
adonis(formula = bc.dist.gut.pa.ctu ~ SignsVetReported, data = metaG) #p = 0.40
adonis(formula = bc.dist.nasal.ab.ctu ~ SignsVetReported, data = metaN) #p = 0.37
adonis(formula = bc.dist.nasal.pa.ctu ~ SignsVetReported, data = metaN) # p = 0.08

adonis(formula = bc.dist.gut.ab.ctu ~ sumSigns, data = metaG) #p = 0.48
adonis(formula = bc.dist.gut.pa.ctu ~ sumSigns, data = metaG) #p = 0.63
adonis(formula = bc.dist.nasal.ab.ctu ~ sumSigns, data = metaN) #p = 0.34
adonis(formula = bc.dist.nasal.pa.ctu ~ sumSigns, data = metaN) # p = 0.014 *


adonis(formula = bc.dist.gut.ab.ctu ~ Household, data = metaG) #p = 0.82
adonis(formula = bc.dist.gut.pa.ctu ~ Household, data = metaG) #p = 0.79
adonis(formula = bc.dist.nasal.ab.ctu ~ Household, data = metaN) #p = 0.48
adonis(formula = bc.dist.nasal.pa.ctu ~ Household, data = metaN) #p = 0.86


adonis(bc.dist.gut.ab.ctu ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.94
adonis(bc.dist.gut.pa.ctu ~ SignsVetReported / Household, data = metaG, strata = factor(metaG$Household)) #0.47
adonis(bc.dist.nasal.ab.ctu ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household))  #0.48
adonis(bc.dist.nasal.pa.ctu ~ SignsVetReported / Household, data = metaN, strata = factor(metaN$Household)) #0.21

adonis(bc.dist.gut.ab.ctu ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.9
adonis(bc.dist.gut.pa.ctu ~ sumSigns / Household, data = metaG, strata = factor(metaG$Household)) #0.80
adonis(bc.dist.nasal.ab.ctu ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household))  #0.52
adonis(bc.dist.nasal.pa.ctu ~ sumSigns / Household, data = metaN, strata = factor(metaN$Household)) #0.01

adonis(formula = bc.dist.gut.ab.ctu ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = bc.dist.gut.pa.ctu ~ Household + SignsVetReported + Household*SignsVetReported, data = metaG)
adonis(formula = bc.dist.nasal.ab.ctu ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)
adonis(formula = bc.dist.nasal.pa.ctu ~ Household + SignsVetReported + Household*SignsVetReported, data = metaN)


# Make box plots of metadata

wilcoxTestVariables = c("WBCTotal", "RBC", "Hemoglobin", "Hematocrit", "PCV", "MCV", "MCH", "MCHC",
                        "PlateletCount", "PlasmaProtein", "NumberNeutrophils", "NumberLymphocytes",
                        "NumberMonocytes", "NumberEosinophils", "NumberBands", "BUN", "Creatinine",
                        "Glucose", "Cholesterol", "TotalProtein", "Albumin", "BilirubinTotal",
                        "CK", "AlkalinePhosphatase", "GGT", "ALT", "Sodium", "Potassium",
                        "Chloride", "Calcium", "Phosphorus", "tCO2", "AnionGap") 
d.wilcox.vet = data.frame(matrix(nrow = length(wilcoxTestVariables), ncol = 2))
rownames(d.wilcox.vet) = wilcoxTestVariables
colnames(d.wilcox.vet) = c("statistic", "pValue")

for(i in 1:length(wilcoxTestVariables)){
  curTest = wilcoxTestVariables[i]
 
  x = as.numeric(meta[which(meta$SignsVetReported == "P"), curTest])
  y = as.numeric(meta[which(meta$SignsVetReported == "N"), curTest])
  print(paste0(c("Testing variable for significance vet: ", curTest, ":"), sep = "", collapse = ""))
  w = wilcox.test(x, y)
  d.wilcox.vet[curTest, "statistic"] = w$statistic
  d.wilcox.vet[curTest, "pValue"] = w$p.value
}
d.wilcox.vet
d.wilcox.vet$qValue = p.adjust(d.wilcox.vet$pValue, method = "fdr")
d.wilcox.vet[which(d.wilcox.vet$pValue < 0.05),]
which(d.wilcox.vet$qValue < 0.05)

cont = rownames(meta)[which(meta$SignsVetReported == "P")]
furtd = rownames(meta)[which(meta$SignsVetReported == "N")]
delta = ave(as.vector(na.omit(meta[cont,"NumberNeutrophils"])))[1] - ave(as.vector(na.omit(meta[furtd,"NumberNeutrophils"])))[1]
sd = sd(meta$NumberNeutrophils, na.rm = TRUE)
sd
power.t.test(power = 0.95, delta = delta, sd = sd) # Power law for neutrophils

delta = abs(ave(as.vector(na.omit(meta[cont,"Albumin"])))[1] - ave(as.vector(na.omit(meta[furtd,"Albumin"])))[1])
delta
sd = sd(meta$Albumin, na.rm = TRUE)
sd
power.t.test(power = 0.95, delta = delta, sd = sd)


p = ggplot(meta, aes(x = SignsVetReported, y = NumberNeutrophils, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 2500, linetype = "dashed", color = pal[2])+
  geom_hline(yintercept = 12500, linetype = "dashed", color = pal[2])+
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Neutrophils", x = "Status", y = "#Neutrophils (cells / uL)") +
  theme(axis.text = element_text(size = 40), plot.title = element_text(size = 40), axis.title = element_text(size = 40))
p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5) + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p


p1 = ggplot(meta, aes(x = SignsVetReported, y = Albumin, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 2.6, linetype = "dashed", color = pal[2])+
  geom_hline(yintercept = 4, linetype = "dashed", color = pal[2])+
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Albumin", x = "Status", y = "Albumin (g/dL)") +
  theme(axis.text = element_text(size = 40), plot.title = element_text(size = 40), axis.title = element_text(size = 40))
p1 = p1 + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5) + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p1

p2 = ggplot(meta, aes(x = SignsVetReported, y = ALT, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 65, linetype = "dashed", color = pal[2])+
  geom_hline(yintercept = 5, linetype = "dashed", color = pal[2])+
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "ALT", x = "Status", y = "ALT (U/L)") +
  theme(axis.text = element_text(size = 40), plot.title = element_text(size = 40), axis.title = element_text(size = 40))
p2 = p2 + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5) + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p2

p3 = ggplot(meta, aes(x = SignsVetReported, y = AnionGap, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 13, linetype = "dashed", color = pal[2])+
  geom_hline(yintercept = 25, linetype = "dashed", color = pal[2])+
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Anion Gap", x = "Status", y = "Anion Gap (mEq/L)") +
  theme(axis.text = element_text(size = 40), plot.title = element_text(size = 40), axis.title = element_text(size = 40))
p3 = p3 + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5) + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p3

setwd(directory.figures)
#pdf("significantlyDifferentCBCChemParameters1.pdf", height = 14, width = 20)
grid.arrange(grobs = list(p, p1, p2, p3), layout_matrix = rbind(c(1, 2), c(4, 5)))
#dev.off()

# 6. Measures of alpha diversity
## A. Calculate alpha diversity for gut and nasal communities
alpha.G.shannon = diversity(t(asvG), index = "shannon", MARGIN = 1, base = exp(1))
alpha.G.simpson = diversity(t(asvG), index = "simpson", MARGIN = 1, base = exp(1))
alpha.N.shannon = diversity(t(asvN), index = "shannon", MARGIN = 1, base = exp(1))
alpha.N.simpson = diversity(t(asvN), index = "simpson", MARGIN = 1, base = exp(1))

# 3B. Make into a data frame
N = 2*length(alpha.G.shannon) + 2*length(alpha.N.shannon)
d = data.frame(matrix(ncol = 5, nrow = N ))
colnames(d) = c("sampleID", "diseaseState", "sampleSite", "diversityMetric", "diversityScore")
rownames(d)
head(d)
for(i in 1:length(alpha.G.shannon)){
  d[i, 1] = names(alpha.G.shannon)[i]
  d[i, 2] = as.vector(metaG[d[i,1], "SignsVetReported"])
  d[i, 3] = "G"
  d[i, 4] = "shannon"
  d[i, 5] = as.vector(alpha.G.shannon)[i]
}

start = length(alpha.G.shannon) + 1
stop = length(alpha.G.shannon) + length(alpha.G.simpson)
for(i in start:stop){
  j = i - length(alpha.G.shannon)
  d[i, 1] = names(alpha.G.simpson)[j]
  d[i, 2] = as.vector(metaG[d[i,1], "SignsVetReported"])
  d[i, 3] = "G"
  d[i, 4] = "simpson"
  d[i, 5] = as.vector(alpha.G.simpson)[j]
}
start = stop + 1
stop = start + length(alpha.N.simpson) - 1
for(i in start:stop){
  
  j = i - start + 1
  print(j)
  d[i, 1] = names(alpha.N.simpson)[j]
  d[i, 2] = as.vector(metaN[d[i,1], "SignsVetReported"])
  d[i, 3] = "N"
  d[i, 4] = "simpson"
  d[i, 5] = as.vector(alpha.N.simpson)[j]
}

start = stop + 1
stop = start + length(alpha.N.shannon) - 1
for(i in start:stop){
  j = i - start + 1
  d[i, 1] = names(alpha.N.shannon)[j]
  d[i, 2] = as.vector(metaN[d[i,1], "SignsVetReported"])
  d[i, 3] = "N"
  d[i, 4] = "shannon"
  d[i, 5] = as.vector(alpha.N.shannon)[j]
}

## 6B. Plot boxplots
### 6Bi. Plot GI samples Shannon diversity index
d.g.shannon = d[which(d$sampleSite == "G"), ]
d.g.shannon = d.g.shannon[which(d.g.shannon$diversityMetric == "shannon"), ]
d.g.shannon.p = d.g.shannon[which(d.g.shannon$diseaseState == "P"), ]
d.g.shannon.n = d.g.shannon[which(d.g.shannon$diseaseState == "N"), ]
wilcox.test(d.g.shannon.p$diversityScore, d.g.shannon.n$diversityScore) # welch two sample t-test. p = 0.22

p = ggplot(d.g.shannon, aes(x = diseaseState, y = diversityScore, color = diseaseState)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3])) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Gut Shannon Diversity", x = "Disease Status", y = "Shannon Diversity") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p = p + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))

p

setwd(directory.figures)
#pdf("shannonGut1.pdf")
p
#dev.off()


### 6Bii. Plot GI samples Simpson diversity index
d.g.simpson = d[which(d$sampleSite == "G"), ]
d.g.simpson = d.g.simpson[which(d.g.simpson$diversityMetric == "simpson"), ]
d.g.simpson.p = d.g.simpson[which(d.g.simpson$diseaseState == "P"), ]
d.g.simpson.n = d.g.simpson[which(d.g.simpson$diseaseState == "N"), ]
wilcox.test(d.g.simpson.p$diversityScore, d.g.simpson.n$diversityScore) # welch two sample t-test. p = 0.39

p = ggplot(d.g.simpson, aes(x = diseaseState, y = diversityScore, color = diseaseState)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3])) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Gut Simpson Diversity", x = "Disease Status", y = "Simpson Diversity") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p = p + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p

setwd(directory.figures)
#pdf("simpsonGut1.pdf")
p
#dev.off()

### 6Biii. Nasal samples Shannon diversity index
d.n.shannon = d[which(d$sampleSite == "N"), ]
d.n.shannon = d.n.shannon[which(d.n.shannon$diversityMetric == "shannon"), ]
d.n.shannon.p = d.n.shannon[which(d.n.shannon$diseaseState == "P"), ]
d.n.shannon.n = d.n.shannon[which(d.n.shannon$diseaseState == "N"), ]
wilcox.test(d.n.shannon.p$diversityScore, d.n.shannon.n$diversityScore) # welch two sample t-test. p = 1

p = ggplot(d.n.shannon, aes(x = diseaseState, y = diversityScore, color = diseaseState)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3])) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Nasopharyngeal Shannon Diversity", x = "Disease Status", y = "Shannon Diversity") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p = p + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))

p 

setwd(directory.figures)
#pdf("shannonNasal1.pdf")
p
#dev.off()

### 6Biv. Nasal samples Simpson diversity index
d.n.simpson = d[which(d$sampleSite == "N"), ]
d.n.simpson = d.n.simpson[which(d.n.simpson$diversityMetric == "simpson"), ]
d.n.simpson.p = d.n.simpson[which(d.n.simpson$diseaseState == "P"), ]
d.n.simpson.n = d.n.simpson[which(d.n.simpson$diseaseState == "N"), ]
wilcox.test(d.n.simpson.p$diversityScore, d.n.simpson.n$diversityScore) # welch two sample t-test. p = 0.42

p = ggplot(d.n.simpson, aes(x = diseaseState, y = diversityScore, color = diseaseState)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3])) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Nasopharyngeal Simpson Diversity", x = "Disease Status", y = "Simpson Diversity") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 22), axis.title = element_text(size = 22))
p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p = p + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
p 

setwd(directory.figures)
#pdf("simpsonNasal1.pdf")
p
#dev.off()


# 7. Describe taxonomic composition of gut and nasal communities.



## B. Make plot bars of all samples Phyla breakdown


### Bv. Make bar plots of all gut samples
setwd(directory.figures)
#pdf("PhylumPlotBarGut1.pdf", width = 12, height = 12)
p = plot_bar(psG, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")+
  labs(title = "Gut Microbiome Phylum Composition", plot.title = element_text(size = 22), axis.title = element_text(size = 32))
p$data$Sample = factor(p$data$Sample, levels = c(rownames(sample_data(psG))[which(sample_data(psG)[, "SignsVetReported"]$SignsVetReported == "P")], rownames(sample_data(psG))[which(sample_data(psG)[, "SignsVetReported"]$SignsVetReported == "N")]))
color = ifelse(metaG[levels(p$data$Sample),"SignsVetReported"] == "P", pal[3], pal[1])
p + theme(plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text.x = element_text(color = color, size = 30, vjust = 0.5),
          axis.text.y = element_text(size = 30), legend.text = element_text(size = 22), legend.title = element_text(size = 30))
p 
#dev.off()


### Bvi. Make bar plot of all nasal samples
setwd(directory.figures)
#pdf("PhylumPlotBarNasal1.pdf", width = 12, height = 12)
p = plot_bar(psN, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(title = "Nasopharyngeal Microbiome Phylum Composition", plot.title = element_text(size = 22), axis.title = element_text(size = 32))
p$data$Sample = factor(p$data$Sample, levels = c(rownames(sample_data(psN))[which(sample_data(psN)[, "SignsVetReported"]$SignsVetReported == "P")], rownames(sample_data(psN))[which(sample_data(psN)[, "SignsVetReported"]$SignsVetReported == "N")]))
color = ifelse(metaN[levels(p$data$Sample),"SignsVetReported"] == "P", pal[3], pal[1])
p = p + theme(plot.title = element_text(size = 30), axis.title = element_text(size = 30), axis.text.x = element_text(color = color, size = 30, vjust = 0.5),
          axis.text.y = element_text(size = 30), legend.text = element_text(size = 22), legend.title = element_text(size = 30))
p 
#dev.off()

## C. Look at breakdowns between positive and negatives across phyla

### Ci. Normalize samples
gut.scale = transform_sample_counts(psG, function(x) 100*x/sum(x))
glom.gut = tax_glom(gut.scale, taxrank = "Phylum")
dat.gut = psmelt(glom.gut)

p = ggplot(dat.gut, aes(x = Phylum, y = Abundance, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  labs(title = "Gut Relative Abundances", color = "Status") +
  theme(plot.title = element_text(size = 50), axis.title = element_text(size = 50), axis.text.y = element_text(size = 50), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 50), legend.title = element_text(size = 50), legend.text = element_text(size = 44), legend.key.size = unit(10, "line")) +
  geom_point(shape = 16, position = position_jitterdodge(), alpha = .6, size = 5)
p

setwd(directory.figures)
#pdf("PhylumBreakdownGut1.pdf", width = 16, height = 14)
p
#dev.off()

nasal.scale = transform_sample_counts(psN, function(x) 100*x/sum(x))
glom.nasal = tax_glom(nasal.scale, taxrank = "Phylum")
dat.nasal = psmelt(glom.nasal)


p = ggplot(dat.nasal, aes(x = Phylum, y = Abundance, color = SignsVetReported)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) + theme_classic() +
  labs(title = "Nasopharyngeal Relative Abundances", color = "Status") +
  theme(plot.title = element_text(size = 50), axis.title = element_text(size = 50), axis.text.y = element_text(size = 50), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 50), legend.title = element_text(size = 50), legend.text = element_text(size = 44), legend.key.size = unit(10, "line")) +
  geom_point(shape = 16, position = position_jitterdodge(), alpha = .6, size = 5)
p + scale_x_discrete(labels = c("N" = "Control", "P" = "FURTD"))
setwd(directory.figures)
#pdf("PhylumBreakdownNasal1.pdf", width = 16, height = 14)
p
#dev.off()

# 7. Look at significant clades accross all samples

## 7A. Get significant clades

allSigG = rownames(aG)[which(aG$pValueCladalConservationAll < 0.05)]
length(allSigG)
NG = length(treeG$tip.label) + length(treeG$node.label)
NG

allSigN = rownames(aN)[which(aN$pValueCladalConservationAll < 0.05)]
length(allSigN)
NN = length(treeN$tip.label) + length(treeN$node.label)
NN

## 7B. Make a dataframe with significant clades labeled gut

d = data.frame(node = 1:NG, nodeName = c(treeG$tip.label, treeG$node.label), pValue = rep(NA, NG), Phylum = rep(NA, NG))
rownames(d) = d$nodeName
d[1:length(treeG$tip.label),"Phylum"] = as.character(taxG[treeG$tip.label,2])

head(d)
d[allSigG, 'pValue'] = pal[1]

setwd(directory.figures)
#pdf("allSigTreeG1.pdf")
p = ggtree(treeG, layout = "circular", branch.length = "none", size = 0.5, aes(color = Phylum)) %<+% d
p + geom_point2(aes(subset = (pValue == pal[1])), color = "black", size = 2, alpha = 0.5) + theme(legend.position = "right")
#dev.off()

x = table(as.vector(aG[allSigG,]$nodeTax))
xtable(x[order(x, decreasing = T)])

## 7C. Make a dataframe with significant clades labeled nasal
d = data.frame(node = 1:NN, nodeName = c(treeN$tip.label, treeN$node.label), pValue = rep(NA, NN), Phylum = rep(NA, NN))
head(d)
rownames(d) = d$nodeName
d[1:length(treeN$tip.label),"Phylum"] = as.character(taxN[treeN$tip.label,2])
head(d)
d[allSigN, 'pValue'] = pal[1]

setwd(directory.figures)
#pdf("allSigTreeN1.pdf")
p = ggtree(treeN, layout = "circular", branch.length = "none", size = .5, aes(color = Phylum)) %<+% d
p + geom_point2(aes(subset = (pValue == pal[1])), color = "black", size = 2, alpha = 0.5) + theme(legend.position = "right")
#dev.off()

x = table(as.vector(aN[allSigN,]$nodeTax))
xtable(x[order(x, decreasing = T)])


# 8. Look at the distributions of clades conserved in FURTD and control samples.
## 8A. Gut samples
nSigG = rownames(aG)[which(aG$groupConservationPValue_N < 0.05)]
pSigG = rownames(aG)[which(aG$groupConservationPValue_P < 0.05)]
length(nSigG) # there are 12 nodes that are significant
length(pSigG) # 146 significantly conserved clades
length(intersect(nSigG, pSigG))

both = intersect(nSigG, pSigG)
nSigGOnly = setdiff(nSigG, pSigG)
pSigGOnly = setdiff(pSigG, nSigG)
d = data.frame(node = 1:NG, nodeName = c(treeG$tip.label, treeG$node.label), Both = rep(NA, NG), Control = rep(NA, NG), FURTD = rep(NA, NG), Phylum = rep(NA, NG))
rownames(d) = d$nodeName

d[1:length(treeG$tip.label),"Phylum"] = as.character(taxG[treeG$tip.label,2])
d[both, "Both"] = "Both"
d[nSigGOnly, "Control"] = "Control"
d[pSigGOnly, "FURTD"] = "FURTD"
head(d)
pal = RColorBrewer::brewer.pal(4, "Set1")
d[which(is.na(d$Phylum[1:length(treeG$tip.label)])), "Phylum"] = "Unlabeled"

l = as.data.frame(cbind(c("FURTD", "Control", "Both"), c(pal[1], pal[2], pal[4])))
l = data.frame(x = runif(n = 3), y = runif(n = 3), Status = c("Both", "Control", "FURTD"))
l = ggplot(l, aes(color = Status, x = x, y = y)) + 
  geom_point(size = 20) + 
  scale_color_manual(values = c(pal[4], pal[2], pal[1])) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30 ), legend.position = "bottom")
l


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
l = g_legend(l)


p = ggtree(treeG, layout = "circular", branch.length = "none", size = .2) %<+% d
p = p + geom_point2(aes(subset = (Both == "Both")), color = pal[4], size = 2, alpha = 0.7) 
p = p + geom_point2(aes(subset = (Control == "Control")), color = pal[2], size = 2, alpha = 0.7) 
p = p + geom_point2(aes(subset = (FURTD == "FURTD")), color = pal[1], size = 2, alpha = 0.7) 
p = gheatmap(p, d[,"Phylum", drop = F], width = 0.1, offset = -1, colnames = F, color = NA)
p = p + scale_color_manual(values = c("FURTD" = pal[1], "Control" = pal[2], "Both" = pal[4])) + 
  labs(title = "Gut Microbiome", subtitle = "Conserved Clades By Status")+
  theme(legend.position = "right", legend.text = element_text(size = 20), plot.title = element_text(size = 50), 
        plot.subtitle = element_text(size = 40)) 
  

lay = rbind(c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1), 
            c(2, 2, 2, 2, 2, 2))
pp = arrangeGrob(p, l, layout_matrix = lay)
grid.arrange(pp)

setwd(directory.figures)
#pdf("gutMicrobiomeConservedCladesByStatus.pdf", width = 10, height = 10)
grid.arrange(pp)
#dev.off()

## 8B. Get clade table of FURTD, Control, and Both. 
x = table(as.vector(aG[pSigGOnly,]$nodeTax))
xtable(x[order(x, decreasing = T)])

x = table(as.vector(aG[nSigGOnly,]$nodeTax))
xtable(x[order(x, decreasing = T)])

x = table(as.vector(aG[both,]$nodeTax))
xtable(x[order(x, decreasing = T)])

## 8C. Nasal samples
nSigN = rownames(aN)[which(aN$groupConservationPValue_N < 0.05)]
pSigN = rownames(aN)[which(aN$groupConservationPValue_P < 0.05)]
length(nSigN) # there are 139 nodes that are significant
length(pSigN) # 182 significantly conserved clades
length(intersect(nSigN, pSigN))

both = intersect(nSigN, pSigN)
nSigNOnly = setdiff(nSigN, pSigN)
length(nSigNOnly)
pSigNOnly = setdiff(pSigN, nSigN)
length(pSigNOnly)
d = data.frame(node = 1:NN, nodeName = c(treeN$tip.label, treeN$node.label), Both = rep(NA, NN), Control = rep(NA, NN), FURTD = rep(NA, NN), Phylum = rep(NA, NN))
rownames(d) = d$nodeName

d[1:length(treeN$tip.label),"Phylum"] = as.character(taxN[treeN$tip.label,2])
d[both, "Both"] = "Both"
d[nSigNOnly, "Control"] = "Control"
d[pSigNOnly, "FURTD"] = "FURTD"
head(d)
pal = RColorBrewer::brewer.pal(4, "Set1")
d[which(is.na(d$Phylum[1:length(treeN$tip.label)])), "Phylum"] = "Unlabeled"

l = as.data.frame(cbind(c("FURTD", "Control", "Both"), c(pal[1], pal[2], pal[4])))
l = data.frame(x = runif(n = 3), y = runif(n = 3), Status = c("Both", "Control", "FURTD"))
l = ggplot(l, aes(color = Status, x = x, y = y)) + 
  geom_point(size = 20) + 
  scale_color_manual(values = c(pal[4], pal[2], pal[1])) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30 ), legend.position = "bottom")
l


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
l = g_legend(l)


p = ggtree(treeN, layout = "circular", branch.length = "none", size = .2) %<+% d
p = p + geom_point2(aes(subset = (Both == "Both")), color = pal[4], size = 2, alpha = 0.7) 
p = p + geom_point2(aes(subset = (Control == "Control")), color = pal[2], size = 2, alpha = 0.7) 
p = p + geom_point2(aes(subset = (FURTD == "FURTD")), color = pal[1], size = 2, alpha = 0.7) 
p = gheatmap(p, d[,"Phylum", drop = F], width = 0.1, offset = -1, colnames = F, color = NA)
p = p + scale_color_manual(values = c("FURTD" = pal[1], "Control" = pal[2], "Both" = pal[4])) + 
  labs(title = "Nasopharyngeal Microbiome", subtitle = "Conserved Clades By Status")+
  theme(legend.position = "right", legend.text = element_text(size = 20), plot.title = element_text(size = 50), 
        plot.subtitle = element_text(size = 40)) 


lay = rbind(c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1), 
            c(2, 2, 2, 2, 2, 2))
pp = arrangeGrob(p, l, layout_matrix = lay)
grid.arrange(pp)

setwd(directory.figures)
#pdf("nasalMicrobiomeConservedCladesByStatus.pdf",  width = 10, height = 10)
grid.arrange(pp)
#dev.off()

## 8D. Get clade table of FURTD, Control, and Both for nasopharyngeal communities
x = table(as.vector(aN[pSigNOnly,]$nodeTax))
xtable(x[order(x, decreasing = T)])

x = table(as.vector(aN[nSigNOnly,]$nodeTax))
xtable(x[order(x, decreasing = T)])

x = table(as.vector(aN[both,]$nodeTax))
xtable(x[order(x, decreasing = T)])

# 9. Random Forest
d.gut.ctu = as.data.frame(t(ctuG))
d.gut.ctu = cbind(d.gut.ctu, as.vector(metaG[rownames(d.gut.ctu), "SignsVetReported"]))
colnames(d.gut.ctu)[ncol(d.gut.ctu)] = "status"
head(d.gut.ctu)
rfGut.ctu = randomForest(status~., data = d.gut.ctu)
rfGut.ctu
rfGut.ctu.roc = roc(rfGut.ctu$y, rfGut.ctu$votes[,1], levels = levels(d.gut.ctu$status))
rfGut.ctu.roc

setwd(directory.figures)
#pdf("GutCTUHist.pdf")
ggplot(data = as.data.frame(rfGut.ctu$importance), aes(rfGut.ctu$importance)) + geom_histogram(col = "black", fill = pal[1], alpha = 0.6) +
  labs(title = "CTU Gut RF Importance", x = "Importance Score", y = "Count") +
  theme(axis.title = element_text(size = 22), title = element_text(size = 22), axis.text = element_text(size = 18))
#dev.off()

d.gut.asv = as.data.frame(t(asvG))
d.gut.asv = cbind(d.gut.asv, as.vector(metaG[rownames(d.gut.asv), "SignsVetReported"]))
colnames(d.gut.asv)[ncol(d.gut.asv)] = "status"
rfGut.asv = randomForest(status~., data = d.gut.asv)
rfGut.asv
rfGut.asv.roc = roc(rfGut.asv$y, rfGut.asv$votes[,1], levels = levels(d.gut.asv$status))
rfGut.asv.roc

setwd(directory.figures)
#pdf("GutASVHist.pdf")
ggplot(data = as.data.frame(rfGut.asv$importance), aes(rfGut.asv$importance)) + geom_histogram(col = "black", fill = pal[2], alpha = 0.6) +
  labs(title = "ASV Gut RF Importance", x = "Importance Score", y = "Count") +
  theme(axis.title = element_text(size = 22), title = element_text(size = 22), axis.text = element_text(size = 18))
#dev.off()

d.nasal.ctu = as.data.frame(t(ctuN))
d.nasal.ctu =  cbind(d.nasal.ctu, as.vector(metaN[rownames(d.nasal.ctu), "SignsVetReported"]))
colnames(d.nasal.ctu)[ncol(d.nasal.ctu)] = "status"
head(d.nasal.ctu)
rfNasal.ctu = randomForest(status~., data = d.nasal.ctu)
rfNasal.ctu
rfNasal.ctu.roc = roc(rfNasal.ctu$y, rfNasal.ctu$votes[,1], levels = levels(d.nasal.ctu$status))
rfNasal.ctu.roc

setwd(directory.figures)
#pdf("NasalCTUHist.pdf")
ggplot(data = as.data.frame(rfNasal.ctu$importance), aes(rfNasal.ctu$importance)) + geom_histogram(col = "black", fill = pal[3], alpha = 0.6) +
  labs(title = "CTU Nasal Importance", x = "Importance Score", y = "Count") +
  theme(axis.title = element_text(size = 22), title = element_text(size = 22), axis.text = element_text(size = 18))
#dev.off()

d.nasal.asv = as.data.frame(t(asvN))
d.nasal.asv = cbind(d.nasal.asv, as.vector(metaN[rownames(d.nasal.asv), "SignsVetReported"]))
colnames(d.nasal.asv)[ncol(d.nasal.asv)] = "status"
rfNasal.asv = randomForest(status~., data = d.nasal.asv)
rfNasal.asv
rfNasal.asv.roc = roc(rfNasal.asv$y, rfNasal.asv$votes[,1], levels = levels(d.nasal.asv$status))
rfNasal.asv.roc

setwd(directory.figures)
#pdf("NasalASVHist.pdf")
ggplot(data = as.data.frame(rfNasal.asv$importance), aes(rfNasal.asv$importance)) + geom_histogram(col = "black", fill = pal[4], alpha = 0.6) +
  labs(title = "CTU Nasal Importance", x = "Importance Score", y = "Count") +
  theme(axis.title = element_text(size = 22), title = element_text(size = 22), axis.text = element_text(size = 18))
#dev.off()

#Name the gut PCR samples

PCR = pcr
rownames(PCR) = PCR$GutID
PCR = PCR[,-c(6,7)]
gutPCR = gutPCR[rownames(rfGut.asv$votes),]
gutPCR$status = metaG[rownames(gutPCR), "SignsVetReported"]
gutPCR.rf = randomForest(status~., data = gutPCR)
gutPCR.roc = roc(gutPCR.rf$y, gutPCR.rf$votes[,1], levels = levels(gutPCR$status))
gutPCR.roc


# Determine if incorporation of phylogenetic infomration includes classifier accuracy
rfGutASV.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.gut.asv, predictor = "status", ntree = 100, name = "gutASV")
rfGutCTU.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.gut.ctu, predictor = "status", ntree = 100, name = "gutCTU")
rfNasalASV.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.nasal.asv, predictor = "status", ntree = 100, name = "nasalASV")
rfNasalCTU.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.nasal.ctu, predictor = "status", ntree = 100, name = "nasalCTU")
rfPCR.100 = getROCCurveAndConfidenceIntervals(M = 100, data = gutPCR, predictor = "status", ntree = 100, name = "PCR")
d.nasal = as.data.frame(t(ctuN))
d.nasal = cbind(d.nasal[rownames(d.nasal),], t(asvN)[rownames(d.nasal),])
d.nasal =  cbind(d.nasal[rownames(d.nasal),], as.vector(metaN[rownames(d.nasal), "SignsVetReported"]))
colnames(d.nasal)[ncol(d.nasal)] = "status"
rfNasalBoth.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.nasal, predictor = "status", name = "nasalBoth", ntree = 100)

d.gut = as.data.frame(t(ctuG))
d.gut = cbind(d.gut[rownames(d.gut),], t(asvG)[rownames(d.gut),])
d.gut =  cbind(d.gut[rownames(d.gut),], as.vector(metaG[rownames(d.gut), "SignsVetReported"]))
colnames(d.gut)[ncol(d.gut)] = "status"
rfGutBoth.100 = getROCCurveAndConfidenceIntervals(M = 100, data = d.gut, predictor = "status", name = "gutBoth", ntree = 100)


gutASVAUC = melt(data.frame(rfGutASV.100$AUC, rep("GutASV", length(rfGutASV.100$AUC))))
colnames(gutASVAUC) = c("group", "variable", "value")
nasalASVAUC = melt(data.frame(rfNasalASV.100$AUC, rep("OropharyngealASV", length(rfNasalASV.100$AUC))))
colnames(nasalASVAUC) = c("group", "variable", "value")
gutCTUAUC = melt(data.frame(rfGutCTU.100$AUC, rep("GutCTU", length(rfGutCTU.100$AUC))))
colnames(gutCTUAUC) = c("group", "variable", "value")
nasalCTUAUC = melt(data.frame(rfNasalCTU.100$AUC, rep("OropharyngealCTU", length(rfNasalCTU.100$AUC))))
colnames(nasalCTUAUC) = c("group", "variable", "value")
gutAUC = melt(data.frame(rfGutBoth.100$AUC, rep("gutBoth", length(rfGutBoth.100$AUC))))
colnames(gutAUC) = c("group", "variable", "value")
nasalAUC = melt(data.frame(rfNasalBoth.100$AUC, rep("nasalBoth", length(rfNasalBoth.100$AUC))))
colnames(nasalAUC) = c("group", "variable", "value")
pcrAUC = melt(data.frame(rfPCR.100$AUC, rep("pcr", length(rfPCR.100$AUC))))
colnames(pcrAUC) = c("group", "variable", "value")

d = rbind(gutASVAUC, gutCTUAUC)
d = rbind(d, gutAUC)
d = rbind(d, nasalASVAUC)
d = rbind(d, nasalCTUAUC)
d = rbind(d, nasalAUC)
d = rbind(d, pcrAUC)
head(d)

shapiro.test(gutAUC$value) # p = 0.0002
shapiro.test(gutASVAUC$value) # p = 0.0002
shapiro.test(gutCTUAUC$value) # p = 0.4
shapiro.test(nasalASVAUC$value) # p = 0.0003
shapiro.test(nasalCTUAUC$value)# p = 0.03
shapiro.test(nasalAUC$value) # p < 2.04e-05

wilcox.test(gutAUC$value, gutASVAUC$value) # p = 3.523e-07
wilcox.test(gutAUC$value, gutCTUAUC$value) # p = 0.02921
wilcox.test(gutAUC$value, nasalAUC$value) # p = 0.135
wilcox.test(gutAUC$value, nasalASVAUC$value) # p = 0.0001138
wilcox.test(gutAUC$value, nasalCTUAUC$value) # p = 0.016
wilcox.test(gutAUC$value, pcrAUC$value) # p < 2.2e-16

wilcox.test(gutASVAUC$value, gutCTUAUC$value) # p = 6.483e-12
wilcox.test(gutASVAUC$value, nasalAUC$value) # p = 1.476e-10
wilcox.test(gutASVAUC$value, nasalASVAUC$value) # p = 4.35e-16
wilcox.test(gutASVAUC$value, nasalCTUAUC$value) # p = 1.126e-12
wilcox.test(gutASVAUC$value, pcrAUC$value) # p < 2.2e-16

wilcox.test(gutCTUAUC$value, nasalAUC$value) # p = 0.44
wilcox.test(gutCTUAUC$value, nasalASVAUC$value) # p = 0.042
wilcox.test(gutCTUAUC$value, nasalCTUAUC$value) # p = 0.77
wilcox.test(gutCTUAUC$value, pcrAUC$value) # p < 2.2e-16

wilcox.test(nasalAUC$value, nasalASVAUC$value) # p = 0.008
wilcox.test(nasalAUC$value, nasalCTUAUC$value) # p = 0.29
wilcox.test(gutCTUAUC$value, pcrAUC$value) # p < 2.2e-16

wilcox.test(nasalASVAUC$value, nasalCTUAUC$value) # p = 0.08
wilcox.test(nasalASVAUC$value, pcrAUC$value) # p < 2.2e-16

wilcox.test(nasalCTUAUC$value, pcrAUC$value) # p < 2.2e-16

pal = RColorBrewer::brewer.pal(7, "Dark2")
p = ggplot(data = d, aes(x = group, y = value, fill = group), ylim =c(0, 1.5)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, outlier.alpha = 0.6, notch = F) +
  geom_jitter(shape = 16, size = 0.2, alpha = 0.4) +
  labs(x = "Feature", title = "Phylogenetic Information", y = "AUC") + 
  theme(legend.position = "none", legend.text = element_text(size = 22), legend.title = element_blank(), axis.title = element_text(size = 22), title = element_text(size = 22), axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 16)) +
  scale_fill_manual(values = pal) +
  expand_limits(y = c(0.5, 1.4)) +
  scale_x_discrete(labels = c("Gut ASV", "Gut CTU", "Gut", "Oropharyngeal ASV", "Oropharyngeal CTU", "Oropharyngeal", "PCR")) 
p  


scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
geom_signif(comparisons = list(c("Gut", "Oropharyngeal")), annotations = "*** \n p = 4.076e-13", y_position = 1.10) +
  geom_signif(comparisons = list(c("Gut", "PCR")), annotations = "*** \n p<2.2e-16", y_position = 1.2) +
  geom_signif(comparisons = list(c("Oropharyngeal", "PCR")), annotations = "*** \n p<2.2e-16", y_position = 1.02) 

p

impAve.gut = apply(rfResults.gut$Imp, 2, ave)
impAve.gut = impAve.gut[1,]
impAve.gut = impAve.gut[order(impAve.gut, decreasing = T)]
impScoresOredered.gut = rfResults.gut$Imp[,names(impAve.gut)]
impScoresOredered.gut = melt(impScoresOredered.gut)
head(impScoresOredered.gut)
getFirst = function(x){
  return(substr(x, 1, 3))
}

##EDIT BELOW


# Now select out important variables for further modeling
rfNasalCTUImp = names(rfNasal.ctu$importance[which(rfNasal.ctu$importance >0.2),])
length(rfNasalCTUImp)
x = table(as.vector(aN[rfNasalCTUImp,]$nodeTax))
xtable(x[order(x, decreasing = T)])

rfNasalASVImp = names(rfNasal.asv$importance[which(rfNasal.asv$importance > 0.2),])
length(rfNasalASVImp)
x = table(as.vector(taxN.String[rfNasalASVImp,]))
xtable(x[order(x, decreasing = T)])

rfGutCTUImp = names(rfGut.ctu$importance[which(rfGut.ctu$importance > 0.2),])
length(rfGutCTUImp)
x = table(as.vector(aG[rfGutCTUImp,]$nodeTax))
xtable(x[order(x, decreasing = T)])

rfGutASVImp = names(rfGut.asv$importance[which(rfGut.asv$importance > 0.2),])
length(rfGutASVImp)
length(rfGutASVImp) + length(rfGutCTUImp) + length(rfNasalCTUImp) + length(rfNasalASVImp)
x = table(as.vector(taxG.String[rfGutASVImp,]))
xtable(x[order(x, decreasing = T)])

gutImp.asv = d.gut.asv[,c(rfGutASVImp, "status")]
gutImp.ctu = d.gut.ctu[,c(rfGutCTUImp)]
gutImp = cbind(gutImp.asv[rownames(gutImp.asv),], gutImp.ctu[rownames(gutImp.asv),])
head(gutImp)


nasalImp.asv = d.nasal.asv[,c(rfNasalASVImp, "status")]
nasalImp.ctu = d.nasal.ctu[,c(rfNasalCTUImp)]
nasalImp = cbind(nasalImp.asv[rownames(nasalImp.asv),], nasalImp.ctu[rownames(nasalImp.asv),])
head(nasalImp)

rfNasal.select = randomForest(status~., data = nasalImp)
rfNasal.select
rfNasal.select$importance[order(rfNasal.select$importance[,1], decreasing = T),]
rfNasal.select.roc = roc(rfNasal.select$y, rfNasal.select$votes[,1], levels = levels(nasalImp$status))
rfNasal.select.roc

rfGut.select = randomForest(status~., data = gutImp)
rfGut.select
rfGut.select$importance[order(rfGut.select$importance[,1], decreasing = T),]
rfGut.select.roc = roc(rfGut.select$y, rfGut.select$votes[,1], levels = levels(gutImp$status))
rfGut.select.roc

rownames(gutImp) = str_replace_all(rownames(gutImp), pattern = "f", replace = "")
rownames(nasalImp) = str_replace_all(rownames(nasalImp), pattern = "o", replace = "")  
commonSamps = intersect(rownames(gutImp), rownames(nasalImp))

intersect(colnames(gutImp), colnames(nasalImp)) #no common ASVs / CTUs
all = cbind(gutImp[commonSamps,], nasalImp[commonSamps,])
which(colnames(all) == "status")
all = all[,-24]
head(all)

all.select = randomForest(status~.,data = all)
all.select
all.select.roc = roc(all.select$y, all.select$votes[,1], levels = levels(all$status))
all.select.roc

setwd(directory.data.metadata)
pcr = read.table("PCRResults.txt", sep = "\t",  header = T)


rfResults.nasal = getROCCurveAndConfidenceIntervals(M = 100, data = nasalImp, predictor = "status", ntree = 10, name = "Oropharyngeal")
rfResults.gut = getROCCurveAndConfidenceIntervals(M = 100, data = gutImp, predictor = "status", ntree = 100, name = "Gut")
rfResults.pcr = getROCCurveAndConfidenceIntervals(M = 100, data = gutPCR, predictor = "status", ntree = 100, name = "PCR")

impAve.gut = apply(rfResults.gut$Imp, 2, ave)
impAve.gut = impAve.gut[1,]
impAve.gut = impAve.gut[order(impAve.gut, decreasing = T)]
impScoresOredered.gut = rfResults.gut$Imp[,names(impAve.gut)]
impScoresOredered.gut = melt(impScoresOredered.gut)
head(impScoresOredered.gut)
getFirst = function(x){
  return(substr(x, 1, 3))
}


#Name the gut PCR samples
gutPCR = pcr
rownames(gutPCR) = gutPCR$GutID
gutPCR = gutPCR[,-c(6,7)]
gutPCR = gutPCR[rownames(rfGut.asv$votes),]
gutPCR$status = metaG[rownames(gutPCR), "SignsVetReported"]
gutPCR.rf = randomForest(status~., data = gutPCR)
gutPCR.roc = roc(gutPCR.rf$y, gutPCR.rf$votes[,1], levels = levels(gutPCR$status))
gutPCR.roc


# Get the most important CTUs / ASVs
## Gut features
group = sapply(impScoresOredered.gut$variable, getFirst)
impScoresOredered.gut$group = group

p = ggplot(impScoresOredered.gut, aes(x = variable, y = value, fill = group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.2, outlier.alpha = 0.6, notch = F) +
  geom_jitter(shape = 16, size = 0.2, alpha = 0.4) +
  labs(x = "CTU or ASV", title = "Gut Importance Scores", y = "Mean Decrease in Accuracy") 
p = p + theme(legend.position = "bottom", legend.text = element_text(size = 22), legend.title = element_text(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 22), title = element_text(size = 22)) 
p = p +  scale_fill_discrete(labels = c("ASV", "CTU"), name = "") + scale_color_discrete(name = "group", guide = F) 
p = p + coord_flip() + scale_fill_manual(values = c(pal[1], pal[2]), labels = c("CTU", "ASV"), name = "")
p

setwd(directory.figures)
pdf("gutImportanceScores.pdf")
p
dev.off()

## Nasal features
impAve.nasal = apply(rfResults.nasal$Imp, 2, ave)
impAve.nasal = impAve.nasal[1,]
impAve.nasal = impAve.nasal[order(impAve.nasal, decreasing = T)]
impScoresOredered.nasal = rfResults.nasal$Imp[,names(impAve.nasal)]
impScoresOredered.nasal = melt(impScoresOredered.nasal)
head(impScoresOredered.nasal)
getFirst = function(x){
  return(substr(x, 1, 3))
}
group = sapply(impScoresOredered.nasal$variable, getFirst)
impScoresOredered.nasal$group = group

p = ggplot(impScoresOredered.nasal, aes(x = variable, y = value, fill = group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.2, outlier.alpha = 0.6, notch = F) +
  geom_jitter(shape = 16, size = 0.2, alpha = 0.4) +
  labs(x = "CTU or ASV", title = "Oropharygeal Importance Scores", y = "Mean Decrease in Accuracy") 
p = p + theme(legend.position = "bottom", legend.text = element_text(size = 22), legend.title = element_text(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 22), title = element_text(size = 22)) 
p = p +  scale_fill_discrete(labels = c("ASV", "CTU"), name = "") + scale_color_discrete(name = "group", guide = F) 
p = p + coord_flip() + scale_fill_manual(values = c(pal[1], pal[2]), labels = c("CTU", "ASV"), name = "")
p

setwd(directory.figures)
pdf("nasalImportanceScores.pdf")
p
dev.off()

# Now plot the ROC graph based on the sensitivities matrix
d = rbind(rfResults.gut$rocGraph, rfResults.nasal$rocGraph)
d = as.data.frame(rbind(d, rfResults.pcr$rocGraph))
d$aveTPR = as.double(d$aveTPR)
d$aveFPR = as.double(d$aveFPR)
head(d)
p = ggplot() + 
  geom_path(data = d, aes(x = aveFPR, y = aveTPR, group = name, color = name), size = 2) +
  scale_color_manual(values = pal[1:3]) + 
  theme_classic() + 
  theme(legend.position = "bottom", legend.text = element_text(size = 22), legend.title = element_blank(), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22), axis.title = element_text(size = 22), title = element_text(size = 22)) +
  guides(linetype = guide_legend(override.aes = list(size = 5))) + 
  geom_abline(intercept = 0, slope = 1, size = 2, linetype = "dotted") +
  xlab("FPR") + ylab("TPR")
p

setwd(directory.figures)
pdf("ROCGutvsNasalvsPCR.pdf")
p
dev.off()

p = ggplot() +
  geom_path(data = rfResults.gut$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveTPR)), color = pal[1], size = 2) +
  geom_ribbon(data = rfResults.gut$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveFPR), ymin = as.double(li), ymax = as.double(ui), fill = pal[1], alpha = 0.1)) +
  theme_classic() + 
  theme(legend.position = "none", legend.text = element_text(size = 22), legend.title = element_blank(), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22), axis.title = element_text(size = 22), title = element_text(size = 22)) +
  guides(linetype = guide_legend(override.aes = list(size = 5))) + 
  geom_abline(intercept = 0, slope = 1, size = 2, linetype = "dotted")+
  xlab("FPR") + ylab("TPR") + ggtitle(label = "Gut ROC") + 
  scale_fill_manual(values = pal[1])
p

setwd(directory.figures)
pdf("ROCGutCI.pdf")
p
dev.off()

p = ggplot() +
  geom_path(data = rfResults.nasal$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveTPR)), color = pal[2], size = 2) +
  geom_ribbon(data = rfResults.nasal$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveFPR), ymin = as.double(li), ymax = as.double(ui), fill = pal[2], alpha = 0.1)) +
  theme_classic() + 
  theme(legend.position = "none", legend.text = element_text(size = 22), legend.title = element_blank(), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22), axis.title = element_text(size = 22), title = element_text(size = 22)) +
  guides(linetype = guide_legend(override.aes = list(size = 5))) + 
  geom_abline(intercept = 0, slope = 1, size = 2, linetype = "dotted")+
  xlab("FPR") + ylab("TPR") + ggtitle(label = "Nasal ROC") +
  scale_fill_manual(values = c(pal[2]))
p

setwd(directory.figures)
pdf("ROCNasalCI.pdf")
p
dev.off()


p = ggplot() +
  geom_path(data = rfResults.pcr$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveTPR)), color = pal[3], size = 2) +
  geom_ribbon(data = rfResults.pcr$rocGraph, aes(x = as.double(aveFPR), y = as.double(aveFPR), ymin = as.double(li), ymax = as.double(ui), fill = pal[2], alpha = 0.1)) +
  theme_classic() + 
  theme(legend.position = "none", legend.text = element_text(size = 22), legend.title = element_blank(), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = 22), axis.title = element_text(size = 22), title = element_text(size = 22)) +
  guides(linetype = guide_legend(override.aes = list(size = 5))) + 
  geom_abline(intercept = 0, slope = 1, size = 2, linetype = "dotted")+
  xlab("FPR") + ylab("TPR") + ggtitle(label = "PCR ROC") +
  scale_fill_manual(values = c(pal[3]))
p

setwd(directory.figures)
pdf("ROCPCRCI.pdf")
p
dev.off()

gutAUC = melt(data.frame(rfResults.gut$AUC, rep("Gut", length(rfResults.gut$AUC))))
colnames(gutAUC) = c("group", "variable", "value")
nasalAUC = melt(data.frame(rfResults.nasal$AUC, rep("Oropharyngeal", length(rfResults.nasal$AUC))))
colnames(nasalAUC) = c("group", "variable", "value")
pcrAUC = melt(data.frame(rfResults.pcr$AUC, rep("PCR", length(rfResults.pcr$AUC))))
colnames(pcrAUC) = c("group", "variable", "value")
d = rbind(gutAUC, nasalAUC)
d = rbind(d, pcrAUC)


p = ggplot(data = d, aes(x = group, y = value, fill = group), ylim =c(0, 1.5)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, outlier.alpha = 0.6, notch = F) +
  geom_jitter(shape = 16, size = 0.2, alpha = 0.4) +
  labs(x = "Feature", title = "Microbiome vs PCR Predictor Performance", y = "AUC") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 22), legend.title = element_blank(), axis.text.x = element_text(size = 20), axis.title = element_text(size = 22), title = element_text(size = 22), axis.text.y = element_text(size = 20))+
  scale_fill_manual(values = pal[1:3]) +
  expand_limits(y = c(0.5, 1.4)) +
  geom_signif(comparisons = list(c("Gut", "Oropharyngeal")), annotations = "*** \n p = 4.076e-13", y_position = 1.10) +
  geom_signif(comparisons = list(c("Gut", "PCR")), annotations = "*** \n p<2.2e-16", y_position = 1.2) +
  geom_signif(comparisons = list(c("Oropharyngeal", "PCR")), annotations = "*** \n p<2.2e-16", y_position = 1.02) 

p

setwd(directory.figures)
pdf("AUC.pdf")
p
dev.off()

shapiro.test(gutAUC$value) #0.02
shapiro.test(nasalAUC$value) #0.2
shapiro.test(pcrAUC$value) # < 0.05

wilcox.test(gutAUC$value, nasalAUC$value) #3 *10 ^-13
wilcox.test(gutAUC$value, pcrAUC$value) #<2.2e-16
wilcox.test(nasalAUC$value, pcrAUC$value) #<2.2e-16

## ***EDIT


head(meta)


d.gut.glm.ctu = as.data.frame(d.gut.glm.ctu)
d.ctu = as.data.frame(matrix(nrow = ncol(d.gut.glm.ctu) - 1, ncol = 2))
colnames(d.ctu) = c("node", "pValue")

for(i in 1:948){
  cur = colnames(d.gut.glm.ctu)[i]
  formula = as.formula(paste0("status~", cur))
  g = glm(formula = formula, family = binomial, data = d.gut.glm.ctu)
  d.ctu[i, "node"] = colnames(d.gut.glm.ctu)[i]
  d.ctu[i, "pValue"] = as.vector((summary(g)$coefficients[,4])[2])
  
}

d.ctu$qValue = p.adjust(p = as.vector(d.ctu$pValue), method = "fdr")
d.ctu$pValue[152] = 1
min(d.ctu$pValue)
min(d.ctu$qValue)
length(which(d.ctu$pValue < 0.1))
sigCTUs = d.ctu[which(d.ctu$pValue < 0.1),]
sigCTUs = sigCTUs$node
xtable(aG[sigCTUs,c("cladeSize", "nodeTax")])
head(d.ctu)
hist(d.ctu$qValue)


d.gut.glm.asv = as.data.frame(t(asvG))
d.gut.glm.asv = cbind(d.gut.glm.asv, as.vector(metaG[rownames(d.gut.glm.asv), "SignsVetReported"]))
colnames(d.gut.glm.asv)[949] = "status"
head(d.gut.glm.asv)

d.gut.glm.asv = as.data.frame(d.gut.glm.asv)
d.asv = as.data.frame(matrix(nrow = ncol(asvG), ncol = 2))
colnames(d.asv) = c("asv", "pValue")

for(i in 1:948){
  cur = colnames(d.gut.glm.asv)[i]
  formula = as.formula(paste0("status~", cur))
  g = glm(formula = formula, family = binomial, data = d.gut.glm.asv)
  d.asv[i, "asv"] = colnames(d.gut.glm.asv)[i]
  d.asv[i, "pValue"] = as.vector((summary(g)$coefficients[,4])[2])
  
}

d.asv$qValue = p.adjust(p = as.vector(d.asv$pValue), method = "fdr")
min(d.asv$pValue)
min(d.asv$qValue)
which(d.asv$pValue < 0.05)
hist(d.asv$pValue)
head(d.asv)
hist(d.asv$qValue)

ctu.gut.RF = randomForest(status~., data = d.gut.glm.ctu[,c(sigCTUs, "status")])
ctu.gut.RF
roc(ctu.gut.RF$y, ctu.gut.RF$votes[,1])
ctu.gut.RF

# Now do the same thing but for oropharyngeal communities
d.nasal.glm.ctu = as.data.frame(t(ctuN))
d.nasal.glm.ctu = cbind(d.nasal.glm.ctu, as.vector(metaN[rownames(d.nasal.glm.ctu), "SignsVetReported"]))
N = ncol(d.nasal.glm.ctu)
colnames(d.nasal.glm.ctu)[N] = "status"
head(d.nasal.glm.ctu)

d.nasal.glm.ctu = as.data.frame(d.nasal.glm.ctu)
d.ctu = as.data.frame(matrix(nrow = ncol(d.nasal.glm.ctu) - 1, ncol = 2))
colnames(d.ctu) = c("node", "pValue")

for(i in 1:948){
  cur = colnames(d.nasal.glm.ctu)[i]
  formula = as.formula(paste0("status~", cur))
  g = glm(formula = formula, family = binomial, data = d.nasal.glm.ctu)
  d.ctu[i, "node"] = colnames(d.nasal.glm.ctu)[i]
  d.ctu[i, "pValue"] = as.vector((summary(g)$coefficients[,4])[2])
  
}

d.ctu[which(is.na(d.ctu$pValue)),"pValue"] = 1
d.ctu$qValue = p.adjust(p = as.vector(d.ctu$pValue), method = "fdr")
min(d.ctu$pValue)
max(d.ctu$pValue)
min(d.ctu$qValue)
max(d.ctu$qValue)
length(which(d.ctu$pValue < 0.1))
sigCTUs = d.ctu[which(d.ctu$pValue < 0.1),]
sigCTUs = sigCTUs$node
head(d.ctu)
hist(d.ctu$qValue)


d.nasal.glm.asv = as.data.frame(t(asvN))
d.nasal.glm.asv = cbind(d.nasal.glm.asv, as.vector(metaN[rownames(d.nasal.glm.asv), "SignsVetReported"]))
colnames(d.nasal.glm.asv)[ncol(d.nasal.glm.asv)] = "status"
head(d.gut.glm.asv)

d.nasal.glm.asv = as.data.frame(d.nasal.glm.asv)
d.asv = as.data.frame(matrix(nrow = ncol(asvN), ncol = 2))
colnames(d.asv) = c("asv", "pValue")

for(i in 1:948){
  cur = colnames(d.nasal.glm.asv)[i]
  formula = as.formula(paste0("status~", cur))
  g = glm(formula = formula, family = binomial, data = d.nasal.glm.asv)
  d.asv[i, "asv"] = colnames(d.nasal.glm.asv)[i]
  d.asv[i, "pValue"] = as.vector((summary(g)$coefficients[,4])[2])
  
}

d.asv$qValue = p.adjust(p = as.vector(d.asv$pValue), method = "fdr")
sigASVs = d.asv[which(d.asv$pValue < 0.1),]
sigASVs = sigASVs$asv
length(sigASVs)
min(d.asv$pValue)
max(d.asv$pValue)
min(d.asv$qValue)
max(d.asv$qValue)


ctu.nasal.RF = randomForest(status~., data = d.nasal.glm.ctu[,c(sigCTUs, "status")])
asv.nasal.RF = randomForest(status~., data = d.nasal.glm.asv[,c(sigASVs, "status")])
ctu.nasal.RF
asv.nasal.RF
roc(ctu.nasal.RF$y, ctu.nasal.RF$votes[,1])
roc(asv.nasal.RF$y, asv.nasal.RF$votes[,1])

data.asv.ctu = cbind(d.nasal.glm.ctu[rownames(d.nasal.glm.ctu),], d.nasal.glm.asv[rownames(d.nasal.glm.ctu),])
ctu.asv.nasal.RF = randomForest(status~., data = data.asv.ctu[,c(sigCTUs, sigASVs, "status")])
ctu.asv.nasal.RF
roc(ctu.asv.nasal.RF$y, ctu.asv.nasal.RF$votes[,1])

meta$PCRResultOverview
meta$SignsVetReported
meta[,c("PCRResultOverview", "SignsVetReported")]

confusionMatrix(data = meta$SignsVetReported, reference = meta$PCRResultOverview)

setwd(directory.out)
ss = read.table("sensitivitySpecificty.txt", sep = "\t", col.names = 1)
ss
xtable(ss)
