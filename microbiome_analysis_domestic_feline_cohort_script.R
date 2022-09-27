###############################################################################
#                                                                             #
###############################################################################

#AUTHOR: HOLLY ARNOLD
#DAY: September 25, 2022
#DATE: 20220925
#DESCRIPTION: 
# Provides script for analyses performed in paper chronic clinical signs of 
# upper respiratory microbiomes in a cohort of domestic felines

###############################################################################
#  LOAD LIBRARIES                                                             #
###############################################################################

library(dplyr)
library(phyloseq)
library(ggplot2)
#library(ggtree)

#library(vegan)
#library(ape)
#library(xtable)
#library(gridExtra)
#library(randomForest)
#library(stringr)
#library(pROC)
#library(reshape)
#library(ggsignif)

###############################################################################
#  SET DIRECTORIES                                                            #
###############################################################################

directory.data = file.path(getwd(), "data/")
directory.figures = file.path(getwd(), "figures/")
directory.scripts = getwd()
directory.out = file.path(getwd(), "out/")


###############################################################################
#  LOAD SCRIPTS                                                               #
###############################################################################

setwd(directory.scripts)
source("microbiome_analysis_domestic_feline_cohort_script_functions.R")


###############################################################################
#  SET COLOR PALLETE                                                          #
###############################################################################

pal = RColorBrewer::display.brewer.all()
pal = RColorBrewer::brewer.pal(n = 7, "Dark2")

###############################################################################
#  INPUT DATA                                                                 #
###############################################################################
setwd(directory.data)

# Metadata
meta = readxl::read_excel(path = "Supplementary_Data.xlsx", 
                          sheet = "Animal Metadata",
                          col_types = c("text", "text", "text", "text", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "text", "date", "text", 
                                        "numeric", "numeric", "text", "text", 
                                        "text", "text","text", "text", "text", 
                                        "text", "text", "text", "text", "text", 
                                        "text", "text", "text", "text", "text", 
                                        "date", "text", "text", "text", "text", 
                                        "text", "date", "date", "text", "text", 
                                        "text", "text", "text", "text", "text",
                                        "text", "text", "text", "text", "text", 
                                        "text", "text", "text", "text", "text"),
                          na = "N/A")
meta$GutID = sub(pattern = "$", replacement = "f", x = meta$Animal_ID)
meta$NasalID = sub(pattern = "$", replacement = "o", x = meta$Animal_ID)




# ASV tables
asvG = utils::read.table("asvGutRarify.txt", header = T)
colnames(asvG) = sub(pattern = "X", replacement = "", x = colnames(asvG))
head(asvG)
dim(asvG)

asvN = utils::read.table("asvNasalRarify.txt", header = T)
colnames(asvN) = sub(pattern = "X", replacement = "", x = colnames(asvN))
head(asvN)
dim(asvN)

## ASV Tax tables
taxG = utils::read.table("asvGutTaxRarify.txt")
head(taxG)
dim(taxG)
taxN = utils::read.table("asvNasalTaxRarify.txt")
head(taxN)
dim(taxN)

taxN.String = utils::read.table("asvNasalTaxString.txt", sep = "\t", row.names = 1)
head(taxN.String)
taxG.String = utils::read.table("asvGutTaxString.txt", sep = "\t", row.names = 1)
head(taxG.String)

# CTU Tables
ctuG = t(utils::read.table("ctuGut.txt", sep = "\t", header = T))
head(ctuG)
treeG = ape::read.tree("treeGut.tre")
treeG

ctuN = t(utils::read.table("ctuNasal.txt", sep = "\t", header = T))
head(ctuN)
treeN = ape::read.tree("treeNasal.tre")
treeN

# Convert relevant bloodwork values comparable units across laboratories
bloodwork = readxl::read_excel(path = "Supplementary_Data.xlsx",
                               sheet = "Animal Laboratory Database",
                               range = "A1:Y60",
                               col_names = c("Description", "Units_OSU_VDL",
                                             "Reference_range_OSU_VDL_Low",
                                             "Reference_range_OSU_VDL_High",
                                             "Units_IDEXX",
                                             "Reference_range_IDEXX_Low",
                                             "Reference_range_IDEXX_High",
                                             "1", "2", "6", "5", "11", "12",
                                             "10", "7", "8", "9", "3", "4",
                                             "13", "14", "15", "16", "17",
                                             "18"))

meta = 
  bloodwork %>% 
  filter(Description %in% 
           c("Animal_ID", "Laboratory", "# Neutrophils","Albumin"))  %>% 
  select(!c("Units_OSU_VDL", "Units_IDEXX", "Reference_range_OSU_VDL_Low", 
            "Reference_range_OSU_VDL_High", "Reference_range_IDEXX_Low", 
            "Reference_range_IDEXX_High", "Description")) %>%
  t(.) %>% 
  as.data.frame(.) %>%
  magrittr::set_colnames(c("Animal_ID", "Laboratory", "Neutrophils", "Albumin")) %>%
  tibble(.) %>%
  mutate_at(c("Neutrophils", "Albumin"), as.numeric) %>%
  mutate(Neutrophils = case_when(Neutrophils < 10 ~ Neutrophils*1000,
                                 TRUE ~ Neutrophils)) %>%
  left_join(meta, by = "Animal_ID") 


# Phyloseq objects
metaG = as.data.frame(meta)
rownames(metaG) = metaG$GutID
metaG = metaG[colnames(asvG),]
asvG_ps = phyloseq::phyloseq(phyloseq::otu_table(object = asvG, taxa_are_rows = TRUE), 
                        phyloseq::sample_data(metaG), 
                        phyloseq::tax_table(as.matrix(taxG)),
                        phyloseq::phy_tree(treeG))
phyloseq::sample_data(asvG_ps)$Status = factor(phyloseq::sample_data(asvG_ps)$Status, levels = c("Control", "Case"))
phyloseq::sample_data(asvG_ps)$Household = factor(phyloseq::sample_data(asvG_ps)$Household)

asvG_ps %>% sample_data %>% count(Status)


metaN = as.data.frame(meta)
rownames(metaN) = metaN$NasalID
metaN = metaN[colnames(asvN),]
asvN_ps = phyloseq::phyloseq(phyloseq::otu_table(object = asvN, taxa_are_rows = TRUE), 
                                    phyloseq::sample_data(metaN), 
                                    phyloseq::tax_table(as.matrix(taxN)),
                                    phyloseq::phy_tree(treeN))
phyloseq::sample_data(asvN_ps)$Status = factor(phyloseq::sample_data(asvN_ps)$Status, levels = c("Control", "Case"))
phyloseq::sample_data(asvN_ps)$Household = factor(phyloseq::sample_data(asvN_ps)$Household)

asvN_ps %>% sample_data %>% count(Status)

ctuG_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuG, taxa_are_rows = TRUE), 
                             phyloseq::sample_data(metaG))
phyloseq::sample_data(ctuG_ps)$Status = factor(phyloseq::sample_data(ctuG_ps)$Status, levels = c("Control", "Case"))

ctuN_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuN, taxa_are_rows = TRUE), 
                             phyloseq::sample_data(metaN))
phyloseq::sample_data(ctuN_ps)$Status = factor(phyloseq::sample_data(ctuN_ps)$Status, levels = c("Control", "Case"))

###############################################################################
#  Beta Diversity Gut Microbiome                                              #
###############################################################################

# Center log transform and ordinate
asvG_ps_clr = microbiome::transform(asvG_ps, "clr")
ordG_clr_unifrac = phyloseq::ordinate(asvG_ps_clr, "RDA",  formula = asvG_ps_clr~Status, distance = "unifrac")

# Scree plot
screeG = phyloseq::plot_scree(ordG_clr_unifrac) + 
  geom_bar(stat = "identity", fill = "black") + 
  labs (x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle(label = "Unifrac")

# RDA plot
clr1G = ordG_clr_unifrac$CA$eig[1] / sum(ordG_clr_unifrac$CA$eig)
clr2G = ordG_clr_unifrac$CA$eig[2] / sum(ordG_clr_unifrac$CA$eig)
rdaG = phyloseq::plot_ordination(asvG_ps_clr, ordG_clr_unifrac, type="samples", color="Status", axes = c(1, 2)) + 
  geom_point(size = 6, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 2, color = "black") +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  coord_fixed(clr2G / clr1G) +
  stat_ellipse(aes(group = Status), linetype = 2) +
  ggtitle("Gut Microbiome Unifrac")

beta_diversity_gut = gridExtra::grid.arrange(grobs = list(rdaG, screeG), 
                                             layout_matrix = cbind(c(1, 1, 1, 2, 2)))
beta_diversity_gut
setwd(directory.figures)
ggsave(file = "gut_beta_diversity.pdf", beta_diversity_gut)

# PERMANOVA
clr_dist_matrix_gut = phyloseq::distance(asvG_ps_clr, method = "unifrac")
clr_adonis = vegan::adonis(clr_dist_matrix_gut ~ phyloseq::sample_data(asvG_ps_clr)$Status)
clr_adonis
clr_adonis = vegan::adonis(clr_dist_matrix_gut ~ phyloseq::sample_data(asvG_ps_clr)$Status, strata = phyloseq::sample_data(asvG_ps_clr)$Household)
clr_adonis

# Beta dispersion
dispr_G = vegan::betadisper(clr_dist_matrix_gut, phyloseq::sample_data(asvG_ps_clr)$Status)

# PCOA beta dispersion plot
setwd(directory.figures)
pdf(file = "Beta_Dispersion_PCOA_Gut.pdf")
plot(dispr_G, main = "Gut Microbiome Beta Dispersion", sub = "")
dev.off()

# Box plot beta dispersion
beta_dispersion_gut_boxplot = 
  data.frame("Distance" = dispr_G$distances, 
           "Status" = as.vector(dispr_G$group)) %>%                            
  ggplot(aes(x = Status, y = Distance)) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(aes(color = Status), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c(pal[1], pal[3])) +
  ylab("Distance to Centroid") +
  theme(legend.position = "none") +
  ggtitle("Beta Dispersion Gut Microbiome")

setwd(directory.figures)
ggsave(file = "beta_dispersion_boxplot_Gut.pdf", beta_dispersion_gut_boxplot)


# Beta dispersion permutation
vegan::permutest(dispr_G)

###############################################################################
#  Beta Diversity Nasal Microbiome                                            #
###############################################################################

# Center log transform and ordinate
asvN_ps_clr = microbiome::transform(asvN_ps, "clr")
ordN_clr_unifrac = phyloseq::ordinate(asvN_ps_clr, "RDA",  formula = asvN_ps_clr~Status, distance = "unifrac")

# Scree plot
screeN = phyloseq::plot_scree(ordN_clr_unifrac) + 
  geom_bar(stat = "identity", fill = "black") + 
  labs (x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle(label = "Unifrac")

# RDA plot 
clr1G = ordN_clr_unifrac$CA$eig[1] / sum(ordN_clr_unifrac$CA$eig)
clr2G = ordN_clr_unifrac$CA$eig[2] / sum(ordN_clr_unifrac$CA$eig)
rdaN = phyloseq::plot_ordination(asvN_ps_clr, ordN_clr_unifrac, type="samples", color="Status", axes = c(1, 2)) + 
  geom_point(size = 6, alpha = 0.8) +
  geom_text(aes(label = Household), hjust = 0.5, vjust = 0.5, size = 2, color = "black") +
  scale_color_manual(values = c(pal[1], pal[3]), labels = c("Control", "FURTD")) +
  coord_fixed(clr2G / clr1G) +
  stat_ellipse(aes(group = Status), linetype = 2) +
  ggtitle("Nasal Microbiome Unifrac")

beta_diversity_nasal = gridExtra::grid.arrange(grobs = list(rdaN, screeN), 
                                               layout_matrix = cbind(c(1, 1, 1, 2, 2)))
setwd(directory.figures)
ggsave(file = "nasal_beta_diversity.pdf", beta_diversity_nasal)

# PERMANOVA
clr_dist_matrix_nasal = phyloseq::distance(asvN_ps_clr, method = "unifrac")
clr_adonis = vegan::adonis(clr_dist_matrix_nasal ~ phyloseq::sample_data(asvN_ps_clr)$Status)
clr_adonis
clr_adonis = vegan::adonis(clr_dist_matrix_nasal ~ phyloseq::sample_data(asvN_ps_clr)$Status, strata = phyloseq::sample_data(asvN_ps_clr)$Household)
clr_adonis

# Beta dispersion
dispr_N = vegan::betadisper(clr_dist_matrix_nasal, phyloseq::sample_data(asvN_ps_clr)$Status)

# PCOA beta dispersion plot
setwd(directory.figures)
pdf(file = "Beta_Dispersion_PCOA_Nasal.pdf")
plot(dispr_N, main = "Nasal Microbiome Beta Dispersion", sub = "")
dev.off()

# Box plot beta dispersion
beta_dispersion_nasal_boxplot = 
  data.frame("Distance" = dispr_N$distances, 
             "Status" = as.vector(dispr_N$group)) %>%                            
  ggplot(aes(x = Status, y = Distance)) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(aes(color = Status), size = 2, alpha = 0.8) + 
  scale_color_manual(values = c(pal[1], pal[3])) +
  ylab("Distance to Centroid") +
  theme(legend.position = "none") +
  ggtitle("Beta Dispersion Nasal Microbiome")

setwd(directory.figures)
ggsave(file = "beta_dispersion_boxplot_Nasal.pdf", beta_dispersion_nasal_boxplot)


# Beta dispersion permutation
vegan::permutest(dispr_N)

###############################################################################
#  Alpha Diversity Gut Microbiome                                             #
###############################################################################

# Calculate alpha diversity for three different measures
alpha_diversity_gut = data.frame(
  "Observed" = phyloseq::estimate_richness(asvG_ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(asvG_ps, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(asvG_ps)))), tree = phyloseq::phy_tree(asvG_ps))[, 1],
  "Status" = phyloseq::sample_data(asvG_ps)$Status)

# Plot alpha diversity
alpha_diveristy_gut_plot = 
  alpha_diversity_gut %>%
  tidyr::gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2, size = 2, alpha = 0.8) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none", plot.title = element_text(size = 20)) +
  ggtitle("Alpha Diversity Gut Microbiome") +
  scale_color_manual(values = c(pal[1], pal[3])) 
setwd(directory.figures)
ggsave(filename = "alpha_diversity_gut.pdf", plot = alpha_diveristy_gut_plot)

# Summarize
alpha_diversity_gut %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))

# Test for significantly different alpha diversity measures based on status
wilcox.test(Observed ~ Status, data = alpha_diversity_gut, alternative = "two.sided")
wilcox.test(Shannon ~ Status, data = alpha_diversity_gut, alternative = "two.sided")              
wilcox.test(PD ~ Status, data = alpha_diversity_gut, alternative = "two.sided")



###############################################################################
#  Alpha Diversity Nasal Microbiome                                           #
###############################################################################

# Calculate alpha diversity for three different measures
alpha_diversity_nasal = data.frame(
  "Observed" = phyloseq::estimate_richness(asvN_ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(asvN_ps, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(asvN_ps)))), tree = phyloseq::phy_tree(asvN_ps))[, 1],
  "Status" = phyloseq::sample_data(asvN_ps)$Status)

# Plot alpha diversity
alpha_diveristy_nasal_plot = 
  alpha_diversity_nasal %>%
  tidyr::gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2, size = 2, alpha = 0.8) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none", plot.title = element_text(size = 20)) +
  ggtitle("Alpha Diversity Nasal Microbiome") +
  scale_color_manual(values = c(pal[1], pal[3])) 
setwd(directory.figures)
ggsave(filename = "alpha_diversity_nasal.pdf", plot = alpha_diveristy_nasal_plot)

# Summarize
alpha_diversity_nasal %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))

# Test for significantly different alpha diversity measures based on status
wilcox.test(Observed ~ Status, data = alpha_diversity_nasal, alternative = "two.sided")
wilcox.test(Shannon ~ Status, data = alpha_diversity_nasal, alternative = "two.sided")              
wilcox.test(PD ~ Status, data = alpha_diversity_nasal, alternative = "two.sided")

###############################################################################
#  Phyla Differential Abundance Testing Gut Microbiome                        #
###############################################################################

#Agglomerate to phylum-level and rename
ps_phylum_gut <- phyloseq::tax_glom(asvG_ps, "Phylum")
phyloseq::taxa_names(ps_phylum_gut) <- phyloseq::tax_table(ps_phylum_gut)[, "Phylum"]

# Get top 6 most abundant phyla for differential 
taxa_counts_gut = apply(phyloseq::otu_table(ps_phylum_gut), 1, sum)[order(apply(phyloseq::otu_table(ps_phylum_gut), 1, sum), decreasing = TRUE)]
top_taxa_gut = names(taxa_counts_gut[1:6])
ps_phylum_gut_top = phyloseq::subset_taxa(ps_phylum_gut, Phylum %in% top_taxa_gut)

# Plot the top 6 most abundant phyla
top_phyla_gut = 
  phyloseq::psmelt(ps_phylum_gut_top)
top_phyla_gut$Phylum = factor(top_phyla_gut$Phylum, levels = top_taxa_gut, ordered = FALSE)
top_phyla_gut$Phylum  

top_phyla_gut_plot = 
  top_phyla_gut %>%
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free", nrow = 1) +
  theme(legend.position = "none")
top_phyla_gut_plot
setwd(directory.figures)
ggsave(file = "phylum_abundances_gut.pdf", top_phyla_gut_plot)

# Test for significance between abundance between top 6 phyla
differential_abundance_gut = data.frame(matrix(nrow = 6, ncol = 4))
colnames(differential_abundance_gut) = c("Phylum", "W", "p", "q")
for(i in 1:length(top_taxa_gut)){
  print(paste0(c("Testing Phylum: ", top_taxa_gut[i]), sep = "", collapse = ""))
  
  case_abundance = top_phyla_gut %>% filter(Phylum == top_taxa_gut[i], Status == "Case") %>% select(Abundance) %>% pull(.)
  control_abundance = top_phyla_gut %>% filter(Phylum == top_taxa_gut[i], Status == "Control") %>% select(Abundance) %>% pull(.)
  
  test = wilcox.test(x = case_abundance, y = control_abundance, alternative = "two.sided", exact = FALSE)
  differential_abundance_gut[i, 1] = top_taxa_gut[i]
  differential_abundance_gut[i, 2] = as.vector(test$statistic)
  differential_abundance_gut[i, 3] = as.vector(test$p.value)
  
  print(test)
  
}
differential_abundance_gut$q = p.adjust(p = differential_abundance_gut$p, method = "fdr")
differential_abundance_gut



###############################################################################
#  Phyla Differential Abundance Testing Nasal Microbiome                     #
###############################################################################

#Agglomerate to phylum-level and rename
ps_phylum_nasal <- phyloseq::tax_glom(asvN_ps, "Phylum")
phyloseq::taxa_names(ps_phylum_nasal) <- phyloseq::tax_table(ps_phylum_nasal)[, "Phylum"]

# Get top 6 most abundant phyla for differential 
taxa_counts_nasal = apply(phyloseq::otu_table(ps_phylum_nasal), 1, sum)[order(apply(phyloseq::otu_table(ps_phylum_nasal), 1, sum), decreasing = TRUE)]
top_taxa_nasal = names(taxa_counts_nasal[1:6])
ps_phylum_nasal_top = phyloseq::subset_taxa(ps_phylum_nasal, Phylum %in% top_taxa_nasal)
ps_phylum_nasal_top

# Plot the top 6 most abundant phyla
top_phyla_nasal = 
  phyloseq::psmelt(ps_phylum_nasal_top)
top_phyla_nasal$Phylum = factor(top_phyla_nasal$Phylum, levels = top_taxa_nasal, ordered = FALSE)
top_phyla_nasal$Phylum  

top_phyla_nasal_plot = 
  top_phyla_nasal %>%
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free", nrow = 1) +
  theme(legend.position = "none")
top_phyla_nasal_plot
setwd(directory.figures)
ggsave(file = "phylum_abundances_nasal.pdf", top_phyla_nasal_plot)

# Test for significance between abundance between top 6 phyla
differential_abundance_nasal = data.frame(matrix(nrow = 6, ncol = 4))
colnames(differential_abundance_nasal) = c("Phylum", "W", "p", "q")
for(i in 1:length(top_taxa_nasal)){
  print(paste0(c("Testing Phylum: ", top_taxa_nasal[i]), sep = "", collapse = ""))
  
  case_abundance = top_phyla_nasal %>% filter(Phylum == top_taxa_nasal[i], Status == "Case") %>% select(Abundance) %>% pull(.)
  control_abundance = top_phyla_nasal %>% filter(Phylum == top_taxa_nasal[i], Status == "Control") %>% select(Abundance) %>% pull(.)
  
  test = wilcox.test(x = case_abundance, y = control_abundance, alternative = "two.sided", exact = FALSE)
  differential_abundance_nasal[i, 1] = top_taxa_nasal[i]
  differential_abundance_nasal[i, 2] = as.vector(test$statistic)
  differential_abundance_nasal[i, 3] = as.vector(test$p.value)
  
  print(test)
  
}
differential_abundance_nasal$q = p.adjust(p = differential_abundance_nasal$p, method = "fdr")
differential_abundance_nasal

###############################################################################
#  Diversity Summary Figure                                                   #
###############################################################################


lay = rbind(c(1, 1, 1, 1, 2, 2, 4, 4, 4, 4), 
            c(1, 1, 1, 1, 2, 2, 4, 4, 4, 4), 
            c(1, 1, 1, 1, 3, 3, 4, 4, 4, 4),
            c(1, 1, 1, 1, 3, 3, 4, 4, 4, 4), 
            c(5, 5, 5, 5, 6, 6, 8, 8, 8, 8), 
            c(5, 5, 5, 5, 6, 6, 8, 8, 8, 8), 
            c(5, 5, 5, 5, 7, 7, 8, 8, 8, 8),
            c(5, 5, 5, 5, 7, 7, 8, 8, 8, 8))
            
           
pp = arrangeGrob(layout_matrix = lay)
grid.arrange(pp)


## TO DO

# Alpha Diversity

# Figure 1

# Figure 2

# Figure 3

# Figure 4

# Write out session info

######################################################################



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
