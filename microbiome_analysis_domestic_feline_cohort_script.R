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

taxN.String = utils::read.table("asvNasalTaxString.txt", 
                                sep = "\t", 
                                row.names = 1)
head(taxN.String)
taxG.String = utils::read.table("asvGutTaxString.txt", 
                                sep = "\t", 
                                row.names = 1)
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
asvG_ps = phyloseq::phyloseq(phyloseq::otu_table(object = asvG, 
                                                 taxa_are_rows = TRUE), 
                        phyloseq::sample_data(metaG), 
                        phyloseq::tax_table(as.matrix(taxG)),
                        phyloseq::phy_tree(treeG))
phyloseq::sample_data(asvG_ps)$Status = 
  factor(phyloseq::sample_data(asvG_ps)$Status, levels = c("Control", "Case"))
phyloseq::sample_data(asvG_ps)$Household = 
  factor(phyloseq::sample_data(asvG_ps)$Household)

asvG_ps %>% sample_data %>% count(Status)

metaN = as.data.frame(meta)
rownames(metaN) = metaN$NasalID
metaN = metaN[colnames(asvN),]
asvN_ps = phyloseq::phyloseq(phyloseq::otu_table(object = asvN, 
                                                 taxa_are_rows = TRUE), 
                                    phyloseq::sample_data(metaN), 
                                    phyloseq::tax_table(as.matrix(taxN)),
                                    phyloseq::phy_tree(treeN))
phyloseq::sample_data(asvN_ps)$Status = 
  factor(phyloseq::sample_data(asvN_ps)$Status, levels = c("Control", "Case"))
phyloseq::sample_data(asvN_ps)$Household = 
  factor(phyloseq::sample_data(asvN_ps)$Household)

asvN_ps %>% sample_data %>% count(Status)

ctuG_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuG, 
                                                 taxa_are_rows = TRUE), 
                             phyloseq::sample_data(metaG))
phyloseq::sample_data(ctuG_ps)$Status = 
  factor(phyloseq::sample_data(ctuG_ps)$Status, levels = c("Control", "Case"))

ctuN_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuN, 
                                                 taxa_are_rows = TRUE), 
                             phyloseq::sample_data(metaN))
phyloseq::sample_data(ctuN_ps)$Status = 
  factor(phyloseq::sample_data(ctuN_ps)$Status, levels = c("Control", "Case"))

###############################################################################
#  Beta Diversity Gut Microbiome                                              #
###############################################################################

# Center log transform and ordinate
asvG_ps_clr = microbiome::transform(asvG_ps, "clr")
ordG_clr_unifrac = phyloseq::ordinate(asvG_ps_clr, 
                                      "RDA",  
                                      formula = asvG_ps_clr~Status, 
                                      distance = "unifrac")

# Scree plot
pc_variance_gut = 
  vegan::eigenvals(ordG_clr_unifrac) / sum(vegan::eigenvals(ordG_clr_unifrac))
pc_variance_gut = data.frame("Axis" =  names(pc_variance_gut), 
                             "Variance" = as.vector(pc_variance_gut)) 
pc_variance_gut$Variance = round(x = pc_variance_gut$Variance*100, digits = 1)

pc_variance_gut$Axis = factor(pc_variance_gut$Axis, 
                              levels = rev(pc_variance_gut$Axis))
screeG = pc_variance_gut %>%
  ggplot(aes(y = Axis, x = Variance, fill = Axis)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs (x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle(label = "Unifrac") +
  theme(legend.position = "none")
screeG

# RDA plot
clr1G = ordG_clr_unifrac$CA$eig[1] / sum(ordG_clr_unifrac$CA$eig)
clr2G = ordG_clr_unifrac$CA$eig[2] / sum(ordG_clr_unifrac$CA$eig)
rdaG = phyloseq::plot_ordination(asvG_ps_clr, 
                                 ordG_clr_unifrac, 
                                 type="samples", 
                                 color="Status", 
                                 axes = c(1, 2)) + 
  geom_point(size = 8, alpha = 0.8) +
  geom_text(aes(label = Household), 
            hjust = 0.5, 
            vjust = 0.5, 
            size = 3, 
            color = "black") +
  scale_color_manual(values = c(pal[1], pal[3]), 
                     labels = c("Control", "FURTD")) +
  coord_fixed(clr2G / clr1G) +
  stat_ellipse(aes(group = Status), linetype = 2) +
  ggtitle("Gut Microbiome Unifrac")

beta_diversity_gut = 
  gridExtra::grid.arrange(grobs = list(rdaG, screeG),
                          layout_matrix = cbind(c(1, 1, 1, 2, 2)))
beta_diversity_gut
setwd(directory.figures)
ggsave(file = "gut_beta_diversity.pdf", beta_diversity_gut)

# PERMANOVA
clr_dist_matrix_gut = phyloseq::distance(asvG_ps_clr, method = "unifrac")
clr_adonis = 
  vegan::adonis(clr_dist_matrix_gut ~ phyloseq::sample_data(asvG_ps_clr)$Status)
clr_adonis
clr_adonis = 
  vegan::adonis(clr_dist_matrix_gut ~ phyloseq::sample_data(asvG_ps_clr)$Status, 
                strata = phyloseq::sample_data(asvG_ps_clr)$Household)
clr_adonis

# Beta dispersion
dispr_G = vegan::betadisper(clr_dist_matrix_gut, 
                            phyloseq::sample_data(asvG_ps_clr)$Status)

# PCOA beta dispersion plot
setwd(directory.figures)
pdf(file = "Beta_Dispersion_PCOA_Gut.pdf")
plot(dispr_G, main = "Gut Microbiome Beta Dispersion", sub = "")
dev.off()

# Box plot beta dispersion
beta_dispersion_gut_boxplot = 
  data.frame("Distance" = dispr_G$distances, 
           "Status" = as.vector(dispr_G$group)) %>%                            
  ggplot(aes(x = Status, y = Distance, color = Status)) +
  geom_boxplot(outlier.size = 0, size = 3) +
  geom_jitter(color = "black", size = 1, alpha = 1) + 
  scale_color_manual(values = c(pal[1], pal[3])) +
  ylab("Distance to Centroid") +
  theme(legend.position = "none") +
  ggtitle("Beta Dispersion Gut Microbiome")
beta_dispersion_gut_boxplot

setwd(directory.figures)
ggsave(file = "beta_dispersion_boxplot_Gut.pdf", beta_dispersion_gut_boxplot)


# Beta dispersion permutation
vegan::permutest(dispr_G)

###############################################################################
#  Beta Diversity Nasal Microbiome                                            #
###############################################################################

# Center log transform and ordinate
asvN_ps_clr = microbiome::transform(asvN_ps, "clr")
ordN_clr_unifrac = phyloseq::ordinate(asvN_ps_clr, 
                                      "RDA",  
                                      formula = asvN_ps_clr~Status, 
                                      distance = "unifrac")

# Scree plot
pc_variance_nasal = 
  vegan::eigenvals(ordN_clr_unifrac) / sum(vegan::eigenvals(ordN_clr_unifrac))
pc_variance_nasal = data.frame("Axis" =  names(pc_variance_nasal), 
                               "Variance" = as.vector(pc_variance_nasal)) 
pc_variance_nasal$Variance = round(x = pc_variance_nasal$Variance*100, 
                                   digits = 1)

pc_variance_nasal$Axis = factor(pc_variance_nasal$Axis, 
                                levels = rev(pc_variance_nasal$Axis))
screeN = pc_variance_nasal %>%
  ggplot(aes(y = Axis, x = Variance, fill = Axis)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs (x = "\nAxis", y = "Proportion of Variance\n") +
  ggtitle(label = "Unifrac") +
  theme(legend.position = "none")
screeN

# RDA plot 
clr1G = ordN_clr_unifrac$CA$eig[1] / sum(ordN_clr_unifrac$CA$eig)
clr2G = ordN_clr_unifrac$CA$eig[2] / sum(ordN_clr_unifrac$CA$eig)
rdaN = phyloseq::plot_ordination(asvN_ps_clr, 
                                 ordN_clr_unifrac, 
                                 type="samples", 
                                 color="Status", 
                                 axes = c(1, 2)) + 
  geom_point(size = 8, alpha = 0.8) +
  geom_text(aes(label = Household), 
            hjust = 0.5, 
            vjust = 0.5, 
            size = 3, 
            color = "black") +
  scale_color_manual(values = c(pal[1], pal[3]), 
                     labels = c("Control", "FURTD")) +
  coord_fixed(clr2G / clr1G) +
  stat_ellipse(aes(group = Status), linetype = 2) +
  ggtitle("Nasal Microbiome Unifrac")

beta_diversity_nasal = 
  gridExtra::grid.arrange(grobs = list(rdaN, screeN), 
                          layout_matrix = cbind(c(1, 1, 1, 2, 2)))
setwd(directory.figures)
ggsave(file = "nasal_beta_diversity.pdf", beta_diversity_nasal)

# PERMANOVA
clr_dist_matrix_nasal = phyloseq::distance(asvN_ps_clr, method = "unifrac")
clr_adonis = 
  vegan::adonis(clr_dist_matrix_nasal~phyloseq::sample_data(asvN_ps_clr)$Status)
clr_adonis
clr_adonis = 
  vegan::adonis(clr_dist_matrix_nasal~phyloseq::sample_data(asvN_ps_clr)$Status, 
                strata = phyloseq::sample_data(asvN_ps_clr)$Household)
clr_adonis

# Beta dispersion
dispr_N = vegan::betadisper(clr_dist_matrix_nasal, 
                            phyloseq::sample_data(asvN_ps_clr)$Status)

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
  geom_boxplot(outlier.size = 0, aes(color = Status), size = 3) +
  geom_jitter(color = "black", size = 1, alpha = 0.8) + 
  scale_color_manual(values = c(pal[1], pal[3])) +
  ylab("Distance to Centroid") +
  theme(legend.position = "none") +
  ggtitle("Beta Dispersion Nasal Microbiome")

setwd(directory.figures)
ggsave(file = "beta_dispersion_boxplot_Nasal.pdf", 
       beta_dispersion_nasal_boxplot)


# Beta dispersion permutation
vegan::permutest(dispr_N)

###############################################################################
#  Alpha Diversity Gut Microbiome                                             #
###############################################################################

# Calculate alpha diversity for three different measures
alpha_diversity_gut = data.frame(
  "Observed" = phyloseq::estimate_richness(asvG_ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(asvG_ps, measures = "Shannon"),
  "PD" = 
    picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(asvG_ps)))), 
                tree = phyloseq::phy_tree(asvG_ps))[, 1],
  "Status" = phyloseq::sample_data(asvG_ps)$Status)

# Plot alpha diversity
alpha_diveristy_gut_plot = 
  alpha_diversity_gut %>%
  tidyr::gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value, color = Status)) +
  geom_boxplot(outlier.color = NA, size = 3) +
  geom_jitter(color = "black", height = 0, width = .2, size = 1, alpha = 0.8) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none", plot.title = element_text(size = 20)) +
  ggtitle("Alpha Diversity Gut Microbiome") +
  scale_color_manual(values = c(pal[1], pal[3])) 
alpha_diveristy_gut_plot
setwd(directory.figures)
ggsave(filename = "alpha_diversity_gut.pdf", plot = alpha_diveristy_gut_plot)

# Summarize
alpha_diversity_gut %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))

# Test for significantly different alpha diversity measures based on status
wilcox.test(Observed ~ Status, data = alpha_diversity_gut, 
            alternative = "two.sided")
wilcox.test(Shannon ~ Status, data = alpha_diversity_gut, 
            alternative = "two.sided")              
wilcox.test(PD ~ Status, data = alpha_diversity_gut, 
            alternative = "two.sided")



###############################################################################
#  Alpha Diversity Nasal Microbiome                                           #
###############################################################################

# Calculate alpha diversity for three different measures
alpha_diversity_nasal = data.frame(
  "Observed" = phyloseq::estimate_richness(asvN_ps, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(asvN_ps, measures = "Shannon"),
  "PD" = 
    picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(asvN_ps)))), 
                tree = phyloseq::phy_tree(asvN_ps))[, 1],
  "Status" = phyloseq::sample_data(asvN_ps)$Status)

# Plot alpha diversity
alpha_diveristy_nasal_plot = 
  alpha_diversity_nasal %>%
  tidyr::gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value, color = Status)) +
  geom_boxplot(outlier.color = NA, size = 3) +
  geom_jitter(color = "black", height = 0, width = .2, size = 1, alpha = 0.8) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none", plot.title = element_text(size = 20)) +
  ggtitle("Alpha Diversity Nasal Microbiome") +
  scale_color_manual(values = c(pal[1], pal[3])) 
alpha_diveristy_nasal_plot
setwd(directory.figures)
ggsave(filename = "alpha_diversity_nasal.pdf", 
       plot = alpha_diveristy_nasal_plot)

# Summarize
alpha_diversity_nasal %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))

# Test for significantly different alpha diversity measures based on status
wilcox.test(Observed ~ Status, 
            data = alpha_diversity_nasal, 
            alternative = "two.sided")
wilcox.test(Shannon ~ Status, 
            data = alpha_diversity_nasal, 
            alternative = "two.sided")              
wilcox.test(PD ~ Status, 
            data = alpha_diversity_nasal, 
            alternative = "two.sided")

###############################################################################
#  Phyla Differential Abundance Testing Gut Microbiome                        #
###############################################################################

#Agglomerate to phylum-level and rename
ps_phylum_gut <- phyloseq::tax_glom(asvG_ps, "Phylum")
phyloseq::taxa_names(ps_phylum_gut) <- 
  phyloseq::tax_table(ps_phylum_gut)[, "Phylum"]

# Get top 6 most abundant phyla for differential 
taxa_counts_gut = 
  apply(phyloseq::otu_table(ps_phylum_gut), 
        1, 
        sum)[order(apply(phyloseq::otu_table(ps_phylum_gut), 
                         1, 
                         sum), 
                   decreasing = TRUE)]
top_taxa_gut = names(taxa_counts_gut[1:6])
ps_phylum_gut_top = phyloseq::subset_taxa(ps_phylum_gut, 
                                          Phylum %in% top_taxa_gut)

# Plot the top 6 most abundant phyla
top_phyla_gut = 
  phyloseq::psmelt(ps_phylum_gut_top)
top_phyla_gut$Phylum = factor(top_phyla_gut$Phylum, 
                              levels = top_taxa_gut, ordered = FALSE)
top_phyla_gut$Phylum  

top_phyla_gut_plot = 
  top_phyla_gut %>%
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), 
              height = 0, 
              size = 3, 
              alpha = 0.8, 
              width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free", nrow = 1) 
top_phyla_gut_plot
setwd(directory.figures)
ggsave(file = "phylum_abundances_gut.pdf", top_phyla_gut_plot)

# Test for significance between abundance between top 6 phyla
differential_abundance_gut = data.frame(matrix(nrow = 6, ncol = 4))
colnames(differential_abundance_gut) = c("Phylum", "W", "p", "q")
for(i in 1:length(top_taxa_gut)){
  print(paste0(c("Testing Phylum: ", top_taxa_gut[i]), sep = "", collapse = ""))
  
  case_abundance = top_phyla_gut %>% 
    filter(Phylum == top_taxa_gut[i], Status == "Case") %>% 
    select(Abundance) %>% pull(.)
  control_abundance = top_phyla_gut %>% 
    filter(Phylum == top_taxa_gut[i], Status == "Control") %>% 
    select(Abundance) %>% pull(.)
  
  test = wilcox.test(x = case_abundance, 
                     y = control_abundance, 
                     alternative = "two.sided", 
                     exact = FALSE)
  differential_abundance_gut[i, 1] = top_taxa_gut[i]
  differential_abundance_gut[i, 2] = as.vector(test$statistic)
  differential_abundance_gut[i, 3] = as.vector(test$p.value)
  
  print(test)
  
}
differential_abundance_gut$q = p.adjust(p = differential_abundance_gut$p, 
                                        method = "fdr")
differential_abundance_gut



###############################################################################
#  Phyla Differential Abundance Testing Nasal Microbiome                     #
###############################################################################

#Agglomerate to phylum-level and rename
ps_phylum_nasal <- phyloseq::tax_glom(asvN_ps, "Phylum")
phyloseq::taxa_names(ps_phylum_nasal) <- 
  phyloseq::tax_table(ps_phylum_nasal)[, "Phylum"]

# Get top 6 most abundant phyla for differential 
taxa_counts_nasal = apply(phyloseq::otu_table(ps_phylum_nasal), 
                          1, 
                          sum)[order(apply(phyloseq::otu_table(ps_phylum_nasal), 
                                           1, 
                                           sum), 
                                     decreasing = TRUE)]
top_taxa_nasal = names(taxa_counts_nasal[1:6])
ps_phylum_nasal_top = phyloseq::subset_taxa(ps_phylum_nasal, 
                                            Phylum %in% top_taxa_nasal)
ps_phylum_nasal_top

# Plot the top 6 most abundant phyla
top_phyla_nasal = 
  phyloseq::psmelt(ps_phylum_nasal_top)
top_phyla_nasal$Phylum = factor(top_phyla_nasal$Phylum, 
                                levels = top_taxa_nasal, ordered = FALSE)
top_phyla_nasal$Phylum  

top_phyla_nasal_plot = 
  top_phyla_nasal %>%
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), 
              height = 0, 
              size = 3, 
              alpha = 0.8, 
              width = .2) +
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
  print(paste0(c("Testing Phylum: ", top_taxa_nasal[i]), 
               sep = "", 
               collapse = ""))
  
  case_abundance = top_phyla_nasal %>% 
    filter(Phylum == top_taxa_nasal[i], Status == "Case") %>% 
    select(Abundance) %>% pull(.)
  control_abundance = top_phyla_nasal %>% 
    filter(Phylum == top_taxa_nasal[i], Status == "Control") %>% 
    select(Abundance) %>% pull(.)
  
  test = wilcox.test(x = case_abundance, 
                     y = control_abundance, 
                     alternative = "two.sided", 
                     exact = FALSE)
  differential_abundance_nasal[i, 1] = top_taxa_nasal[i]
  differential_abundance_nasal[i, 2] = as.vector(test$statistic)
  differential_abundance_nasal[i, 3] = as.vector(test$p.value)
  
  print(test)
  
}
differential_abundance_nasal$q = p.adjust(p = differential_abundance_nasal$p, 
                                          method = "fdr")
differential_abundance_nasal

###############################################################################
#  Diversity Summary Figure                                                   #
###############################################################################

# Figure layout
lay = rbind(c(15, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12),
            c(15, 1, 1, 1, 2, 2, 3, 4, 4, 4, 4), 
            c(15, 1, 1, 1, 2, 2, 3, 4, 4, 4, 4), 
            c(15, 1, 1, 1, 2, 2, 3, 4, 4, 4, 4),
            c(15, 1, 1, 1, 2, 2, 3, 4, 4, 4, 4),
            c(15, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17),
            c(15, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9),
            c(15, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9),
            c(15, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9),
            
            c(16, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14),
            c(16, 5, 5, 5, 6, 6, 7, 8, 8, 8, 8), 
            c(16, 5, 5, 5, 6, 6, 7, 8, 8, 8, 8), 
            c(16, 5, 5, 5, 6, 6, 7, 8, 8, 8, 8),
            c(16, 5, 5, 5, 6, 6, 7, 8, 8, 8, 8),
            c(15, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18),
            c(16, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10),
            c(16, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10),
            c(16, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10))

# Make some minimalist figures for the summary figure         
rdaG_diversity = rdaG + 
  ggtitle("Unifrac") +
  theme(legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 18),
        axis.title.x = element_text(color = pal[2]),
        title = element_text(size = 18), 
        legend.position = "bottom") 
screeG_diversity = 
  screeG + 
  ggtitle("PC % of Variance") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 18)) +
  scale_fill_manual(values = c(rep("grey", times = 13), pal[2])) +
  geom_text(aes(label = Variance), hjust = 1, size = 7)
  
beta_dispersion_gut_boxplot_diversity =
  beta_dispersion_gut_boxplot +
  ggtitle("") +
  ylab("Beta Dispersion") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        title = element_text(size = 18)) 
beta_dispersion_gut_boxplot_diversity       

alpha_diveristy_gut_plot_diversity = 
  alpha_diveristy_gut_plot +
  ggtitle("") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        title = element_text(size = 18), 
        strip.text = element_text(size = 20))

alpha_diveristy_gut_plot_diversity
top_phyla_gut_plot_diversity = 
  top_phyla_gut_plot +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 50),
        axis.text.y = element_text(angle = 50, size = 8),
        strip.text = element_text(size = 20),
        legend.position = "none")
top_phyla_gut_plot_diversity

rdaN_diversity = rdaN + 
  ggtitle("Unifrac") +
  theme(legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 18),
        axis.title.x = element_text(color = pal[2]),
        title = element_text(size = 18), 
        legend.position = "bottom") 
rdaN_diversity

screeN_diversity = 
  screeN + 
  ggtitle("PC % of Variance") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 18)) +
  scale_fill_manual(values = c(rep("grey", times = 15), pal[2])) +
  geom_text(aes(label = Variance), hjust = 1, size = 7)
screeN_diversity 

beta_dispersion_nasal_boxplot_diversity =
  beta_dispersion_nasal_boxplot +
  ggtitle("") +
  ylab("Beta Dispersion") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        title = element_text(size = 18)) 
beta_dispersion_nasal_boxplot_diversity       

alpha_diveristy_nasal_plot_diversity = 
  alpha_diveristy_nasal_plot +
  ggtitle("") +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        title = element_text(size = 18), 
        strip.text = element_text(size = 20))
alpha_diveristy_nasal_plot_diversity

top_phyla_nasal_plot_diversity = 
  top_phyla_nasal_plot +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 50),
        axis.text.y = element_text(angle = 50, size = 8),
        strip.text = element_text(size = 20),
        legend.position = "none")
top_phyla_nasal_plot_diversity

sub_title_size = 30
title_size = 40 

alpha_grob = ggpubr::as_ggplot(ggpubr::text_grob("Alpha Diversity", 
                                                 face = "italic", 
                                                 color = "black", 
                                                 hjust = 1.5, 
                                                 size = sub_title_size))
beta_grob = ggpubr::as_ggplot(ggpubr::text_grob("Beta Diversity", 
                                                face = "italic", 
                                                color = "black", 
                                                hjust = 2.5, 
                                                size = sub_title_size))
nasal_grob = ggpubr::as_ggplot(ggpubr::text_grob("Nasal Microbiome", 
                                                 face = "italic", 
                                                 rot = 90, 
                                                 color = "black", 
                                                 hjust = 0.5, 
                                                 size = title_size))
gut_grob = ggpubr::as_ggplot(ggpubr::text_grob("Gut Microbiome", 
                                               face = "italic", 
                                               rot = 90, 
                                               color = "black", 
                                               hjust = -.4, 
                                               size = title_size))
top_taxa_grob = ggpubr::as_ggplot(ggpubr::text_grob("Taxonomic Composition", 
                                                    face = "italic", 
                                                    color = "black", 
                                                    hjust = 2.1, 
                                                    size = sub_title_size))



pp = gridExtra::arrangeGrob(rdaG_diversity, 
                            screeG_diversity, 
                            beta_dispersion_gut_boxplot_diversity, 
                            alpha_diveristy_gut_plot_diversity, 
                            rdaN_diversity, 
                            screeN_diversity, 
                            beta_dispersion_nasal_boxplot_diversity, 
                            alpha_diveristy_nasal_plot_diversity, 
                            top_phyla_gut_plot_diversity, 
                            top_phyla_nasal_plot_diversity, 
                            beta_grob, 
                            alpha_grob, 
                            beta_grob, 
                            alpha_grob, 
                            gut_grob, 
                            nasal_grob, 
                            top_taxa_grob, 
                            top_taxa_grob, 
                            layout_matrix = lay)

setwd(directory.figures)
ggsave(pp, filename = "Diversity_Summary.pdf", height = 20, width = 20)
 




## TO DO


# Figure 1

# Figure 2

# Figure 3

# Figure 4

# Write out session info

######################################################################

#AUTHOR: ARNOLD
#DAY: August 4th, 2021
#DATE: 20210804
#PURPOSE: 
# 1. Metadata analysis
# 2. PCOA correlate with signs
# 3. ASV analysis
# 4. CTU analysis

##############################################################################################################
#LOCATIONS
## 1. LOCAL
## /Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/originalData

## 2. SERVER
## /nfs3/Sharpton_Lab/prod/projects/arnoldhk/furtd2019/version2/2.dada2

#WRITE UP
#See 20200909.tex for further details on this section.
##############################################################################################################

#LIBRARIES
#LIBRARIES
#library(dada2)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(reshape2)
library(pracma)
library(ape)
library(vegan)
#library(devtools)
#library(ggbiplot)
#library(ggfortify)
#library(BiodiversityR)
#library(ggrepel)
library(xtable)
#library(GUniFrac)
library(dplyr)
#library(tidyr)
#library(adespatial)
library(gridExtra)
library(ggtree)
library(randomForest)
library(pROC)
#library(caret)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

# DIRECTORIES
directory.data.original = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/originalData/"
directory.data.rarify = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/rarifiedData/"
directory.data.metadata = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version3/data/"
directory.data = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/"
directory.data.ctu.gut = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/cladalAnalysis/gut/"
directory.data.ctu.nasal = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version2/data/cladalAnalysis/nasal/"
directory.figures = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version3/figures/"
directory.scripts = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version3/scripts/"
directory.out = "/Users/arnoldhk/Desktop/Research/2019CatMicrobiomes/version3/out/"

# LOAD FUNCTIONS
setwd(directory.scripts)
source("catMicrobiomeAnalysisVersion3Functions.R")

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

## 3. Read in metadata tables
setwd(directory.out)
meta = read.table("../../version2/out/metaData.txt", sep = "\t")
head(meta)

# Clean up the metadata table
metaAll = as_tibble(meta)
colnames(metaAll)

library(dplyr)
metaAll = metaAll %>% 
  mutate(
    NasalID = NasalID.x,
    case_control = case_when(
      SignsOwnerReported == "P" ~ "case",
      SignsOwnerReported == "N" ~ "control"
    ), 
    study = case_when( # 
      Household == 1 ~ "A",
      Household == 2 ~ "B",
      Household == 3 ~ "C",
      Household == 4 ~ "D",
      Household == 5 ~ "E",
      Household == 6 ~ "F",
      Household == 7 ~ "G",
      Household == 8 ~ "H"
    ),
    FCV = FelineCalcivirusPCR,
    Bordetella = BordetellaBronchisepticaPCR,
    Chlamydophila = ChlamydophilaFelisPsittaciPCR,
    FHV = FelineHerpesvirus1PCR,
    Mycoplasma = MycoplasmaFelisPCR,
    PCR = PCRResultOverview,
    
    
  ) %>%
  dplyr::select(-PatientName.x, 
                -NasalID.y, 
                -GutID.1, 
                -NasalID.1, 
                -NasalIDOwnerSigns, 
                -GutIDOwnerSigns, 
                -RecordNumber, 
                -PatientName.y, 
                -PatientName.x,
                -NasalID.x, 
  ) %>%
  mutate(
    OcularHx = OcularDischarge + OcularCrusts + Conjunctivitis +OcularSwelling + OcularClouding + CornealUlceration + Uveitis,
    OralHx = Stomatitis + Glossitis + Faucitis,
    RespiratoryHx = NasalDischarge + NasalCrusts + Rhinitis + LungsRespiratory,
    SignsHx = (1/7)*OcularHx + (1/3)*OralHx + (1/4)*RespiratoryHx
    
  ) %>% 
  
  print(width = Inf)

signsCovariates = c("Eyes", "Nose", "Oral", "Respiratory")
bloodworkCovariates = c("WBCTotal", "PCV", "PlasmaProtein", "NumberNeutrophils", "NumberLymphocytes", "NumberMonocytes", "NumberEosinophils", "TotalProtein", "Albumin")

metaG = as.data.frame(metaAll)
rownames(metaG) = metaG$GutID
metaG = metaG[colnames(asvG),]
metaG$study = factor(metaG$study)
metaG$case_control = factor(metaG$case_control)
metaG$ID = rownames(metaG)


metaN = as.data.frame(metaAll)
rownames(metaN) = metaN$NasalID
metaN = metaN[colnames(asvN),]
metaN$study = factor(metaN$study)
metaN$case_control = factor(metaN$case_control)
metaN$ID = rownames(metaN)

# 4. Read in the cladal analysis files
setwd(directory.data.ctu.gut)
ctuG = t(read.table("ctuGut.txt", sep = "\t", header = T))
head(ctuG)
treeG = read.tree("new_prepped_tree.tre")
treeG

#aG = getCladalAttributes(nodes2tax = "gutCladeStat_nodes2tax.txt", size = "gutCladeStat_clade_size.txt",nodeTax = "taxCladesGut.txt", groupTest = "gutGroupTestSignsVetReport.txt_stats.txt", pTest = "gutPTest1000_stats.txt")
#head(aG)

setwd(directory.data.ctu.nasal)
ctuN = t(read.table("ctuNasal.txt", sep = "\t", header = T))
head(ctuN)
treeN = read.tree("new_prepped_tree.tre")
treeN
#aN = getCladalAttributes(nodes2tax = "nasalCladeStat_nodes2tax.txt", size = "nasalCladeStat_clade_size.txt", nodeTax = "taxCladesNasal.txt",groupTest = "nasalGroupTestSignsVetReport.txt_stats.txt", pTest = "nasalPTest1000_stats.txt")
#head(aN)

# Make phyloseq objects
library(stringr)

psNASV = phyloseq(otu_table(asvN, taxa_are_rows = T), tax_table(as.matrix(taxN)), sample_data(metaN), phy_tree(treeN))
taxNCTU = cladifierGetRefTreeTaxonomy(tree = treeN, ref = taxN)
taxNCTU = taxNCTU[rownames(ctuN),]
psNCTU = phyloseq(otu_table(ctuN, taxa_are_rows = T), tax_table(as.matrix(taxNCTU)), sample_data(metaN))
psNCTU
psGASV = phyloseq(otu_table(asvG, taxa_are_rows = T), tax_table(as.matrix(taxG)), sample_data(metaG), phy_tree(treeG))
psGASV
taxGCTU = cladifierGetRefTreeTaxonomy(tree = treeG, ref = taxG)
taxGCTU = taxGCTU[rownames(ctuG),]
psGCTU = phyloseq(otu_table(ctuG, taxa_are_rows = T), tax_table(as.matrix(taxGCTU)), sample_data(metaG))


##############################################################################################################
# 1. SIGNS
##############################################################################################################

setwd(directory.data.metadata)
oCx = read.table("ownerReportedSigns.txt", sep = "\t", header = T)
oCx = oCx[1:18,1:20]

oCxP = oCx[which(oCx$OwnerReportedStatus == "P"), 4:17]
oCxN = oCx[which(oCx$OwnerReportedStatus == "N"), 4:17]
oCxP = apply(oCxP, 2, sum)
oCxN = apply(oCxN, 2, sum)
hx = data.frame("Sign" = names(oCxP), "Case" = oCxP, "Control" = oCxN)
print(xtable(hx),  include.rownames = FALSE)

# 2. Get table for vet reported signs
vCx = read.table("vetReportedSigns.txt", sep = "\t", header = T)
head(vCx)

vCxP = vCx[which(vCx$OwnerReportedStatus == "P"), c("Eyes", "Nose", "Oral", "LungsRespiratory")]
vCxSumP = apply(vCxP, 2, sum)

vCxN = vCx[which(vCx$OwnerReportedStatus == "N"), c("Eyes", "Nose", "Oral", "LungsRespiratory")]
vCxSumN = apply(vCxN, 2, sum)
vCxSumN

vx = data.frame("Sign" = names(vCxSumN), "Case" = vCxSumP, "Control" = vCxSumN)
print(xtable(vx),  include.rownames = FALSE)
rm(list = c("vCxN", "vCxP", "vCxSumN", "vCx", "vCxSumP", "oCxP", "oCx", "oCxN"))

##############################################################################################################
# 2. BETA DIVERSITY PLOTS
##############################################################################################################


metaG %>% group_by(study, case_control) %>% select(study, case_control) %>% arrange(study)
metaN %>% group_by(study, case_control) %>% select(study, case_control) %>% arrange(study)

bc.dist.n = vegdist(t(otu_table(psNASV)), method = "bray", binary = FALSE)
bc.dist.g = vegdist(t(otu_table(psGASV)), method = "bray", binary = FALSE)
jac.dist.g = vegdist(t(otu_table(psGASV)), method = "jaccard")
jac.dist.n = vegdist(t(otu_table(psNASV)), method = "jaccard")
wUni.dist.n = UniFrac(physeq = psNASV, weighted = TRUE)
wUni.dist.g = UniFrac(physeq = psGASV, weighted = TRUE)
uUni.dist.n = UniFrac(physeq = psNASV, weighted = FALSE)
uUni.dist.g = UniFrac(physeq = psGASV, weighted = FALSE)

### Gut
#### Nasal Cx
adonis(formula = bc.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = jac.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = wUni.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*

adonis(formula = bc.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = jac.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = wUni.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) 
adonis(formula = uUni.dist.g ~ study*Nose, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*

#### Ocular signs
adonis(formula = bc.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = jac.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = wUni.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*

adonis(formula = bc.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = jac.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = wUni.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) 
adonis(formula = uUni.dist.g ~ study*Eyes, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*

#### Oral Cx
adonis(formula = bc.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = jac.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*
adonis(formula = wUni.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum") #*

adonis(formula = bc.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = jac.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*
adonis(formula = wUni.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) 
adonis(formula = uUni.dist.g ~ study*Oral, data = sample.data.frame(psGASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psGASV)$study)) #*

### Nasal Microbiomes
#### Ocular Signs
adonis(formula = bc.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum" )
adonis(formula = jac.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = wUni.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 

adonis(formula = bc.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study))
adonis(formula = jac.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = wUni.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = uUni.dist.n ~ study*Eyes, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 

#### Nasal Signs
adonis(formula = bc.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum" )
adonis(formula = jac.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = wUni.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 

adonis(formula = bc.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum",  strata = factor(sample.data.frame(psNASV)$study))
adonis(formula = jac.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum",  strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = wUni.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum",  strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = uUni.dist.n ~ study*Nose, data = sample.data.frame(psNASV), contr.unordered = "contr.sum",  strata = factor(sample.data.frame(psNASV)$study)) 

#### Oral Signs
adonis(formula = bc.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum" )
adonis(formula = jac.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = wUni.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 
adonis(formula = uUni.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum") 

adonis(formula = bc.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study))
adonis(formula = jac.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = wUni.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 
adonis(formula = uUni.dist.n ~ study*Oral, data = sample.data.frame(psNASV), contr.unordered = "contr.sum", strata = factor(sample.data.frame(psNASV)$study)) 

## Beta Diversity Plots
### Gut
#### Bray Curtis
bc.cap.gut = capscale(bc.dist.g ~ Nose + Condition(study), data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
p1 = plotDBRDAGut1Constraint(cap = bc.cap.gut, 
                             meta = sample.data.frame(psGASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Gut Taxonomic Abundance (Bray Curtis)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p1

setwd(directory.figures)
png(file = "GutBCConditionStudy.png")
p1
dev.off()

#### Jaccard
jac.cap.gut = capscale(jac.dist.g ~ Nose + Condition(study), data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
p2 = plotDBRDAGut1Constraint(cap = jac.cap.gut, 
                             meta = sample.data.frame(psGASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Gut Taxonomic Presence / Absence (Jaccard)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p2

setwd(directory.figures)
png(file = "GutJacConditionStudy.png")
p2
dev.off()

#### Weighted Unifrac
wUni.cap.gut = capscale(wUni.dist.g ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
p3 = plotDBRDAGut2Constraint(cap = wUni.cap.gut, 
                             meta = sample.data.frame(psGASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Gut Phylogenetic Abundance\n(Weighted UNIFRAC)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p3

setwd(directory.figures)
png(file = "GutWUnifracConditionStudy.png")
p3
dev.off()

uUni.cap.gut = capscale(uUni.dist.g ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
p4 = plotDBRDAGut2Constraint(cap = uUni.cap.gut, 
                             meta = sample.data.frame(psGASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Gut Phylogenetic Presence / Absence\n(Unweighted UNIFRAC)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p4

setwd(directory.figures)
png(file = "GutUUnifracConditionStudy.png")
p4
dev.off()

### Nasal Microbiomes
#### Bray Curtis
bc.cap.nasal = capscale(bc.dist.n ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p5 = plotDBRDAGut2Constraint(cap = bc.cap.nasal, 
                             meta = sample.data.frame(psNASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Nasal Taxonomic Abundance (Bray Curtis)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p5

setwd(directory.figures)
png(file = "NasalBCStudy.png")
p5
dev.off()

#### Jaccard
jac.cap.nasal = capscale(jac.dist.n ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p6 = plotDBRDAGut2Constraint(cap = jac.cap.nasal, 
                             meta = sample.data.frame(psNASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Nasal Taxonomic Presence / Absence (Jaccard)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p6

setwd(directory.figures)
png(file = "NasalJacStudy.png")
p6
dev.off()

bc.cap.nasal = capscale(bc.dist.n ~ Nose + Condition(study), data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p7 = plotDBRDAGut1Constraint(cap = bc.cap.nasal, 
                             meta = sample.data.frame(psNASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Nasal Taxonomic Abundance (Bray Curtis)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p7

setwd(directory.figures)
png(file = "NasalBCConditionStudy.png")
p7
dev.off()

#### Jaccard
jac.cap.nasal = capscale(jac.dist.n ~ Nose + Condition(study), data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p8 = plotDBRDAGut1Constraint(cap = jac.cap.nasal, 
                             meta = sample.data.frame(psNASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Nasal Taxonomic Presence / Absence (Jaccard)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p8

setwd(directory.figures)
png(file = "NasalJacConditionStudy.png")
p8
dev.off()

#### Weighted Unifrac
wUni.cap.nasal = capscale(wUni.dist.n ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p9 = plotDBRDAGut2Constraint(cap = wUni.cap.nasal, 
                             meta = sample.data.frame(psNASV), 
                             ID = "ID", 
                             colorGroup = "case_control", 
                             colorLab = "Case|Control", 
                             sizeGroup = "Nose",
                             sizeLab = "Nasal Signs",
                             shapeGroup = "study", 
                             shapeLab = "Household", 
                             title = "Nasal Phylogenetic Abundance\n(Weighted UNIFRAC)",
                             sizeValues = c(6, 12),
                             shapeValues = c(7:14))
p9

setwd(directory.figures)
png(file = "NasalWUnifracConditionStudy.png")
p9
dev.off()

uUni.cap.nasal = capscale(uUni.dist.n ~ Nose + study, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
p10 = plotDBRDAGut2Constraint(cap = uUni.cap.nasal, 
                              meta = sample.data.frame(psNASV), 
                              ID = "ID", 
                              colorGroup = "case_control", 
                              colorLab = "Case|Control", 
                              sizeGroup = "Nose",
                              sizeLab = "Nasal Signs",
                              shapeGroup = "study", 
                              shapeLab = "Household", 
                              title = "Nasal Phylogenetic Presence / Absence\n(Unweighted UNIFRAC)",
                              sizeValues = c(6, 12),
                              shapeValues = c(7:14))
p10

setwd(directory.figures)
png(file = "NasalUUnifracConditionStudy.png")
p10
dev.off()

## LDA of microbial communities
setwd(directory.out)
phyloseqCompanion::phyloseq2lefse(psNASV, file.name = "lefseNasal.txt", covars = c("Nose", "Household"))
phyloseqCompanion::phyloseq2lefse(psGASV, file.name = "lefseGut.txt", covars = c("Nose", "Household"))

bc.cap.gut = capscale(bc.dist.g ~ Nose + Oral + Eyes + Condition(study), data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
bc.cap.gut.0 = capscale(bc.dist.g ~ 1, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
plot(bc.cap.gut)
anova(bc.cap.gut)
bc.cap.selected.gut = ordistep(bc.cap.gut.0, Pin = 0.1, scope = formula(bc.cap.gut), direction = 'both', permutations = 999)
bc.cap.selected.gut$anova
bc.cap.gut = capscale(bc.dist.g ~ Nose + Condition(study) , data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
anova(bc.cap.gut)

jac.cap.gut = capscale(jac.dist.g ~ Nose + Oral + Eyes + Condition(study), data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
jac.cap.gut.0 = capscale(jac.dist.g ~ 1, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
plot(jac.cap.gut)
anova(jac.cap.gut)
jac.cap.selected.gut = ordistep(jac.cap.gut.0, Pin = 0.1, scope = formula(jac.cap.gut), direction = 'both', permutations = 999)
jac.cap.selected.gut$anova
jac.cap.gut = capscale(jac.dist.g ~ Nose + Condition(study) , data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
anova(jac.cap.gut)

wUni.cap.gut = capscale(wUni.dist.g ~ Nose + Oral + Eyes + Condition(study) , data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
wUni.cap.gut.0 = capscale(wUni.dist.g ~ 1, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
plot(wUni.cap.gut)
anova(wUni.cap.gut)
wUni.cap.selected.gut = ordistep(wUni.cap.gut.0, Pin = 0.1, scope = formula(wUni.cap.gut), direction = 'both', permutations = 999)
wUni.cap.selected.gut$anova

uUni.cap.gut = capscale(uUni.dist.g ~ Nose + Oral + Eyes + Condition(study) , data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
uUni.cap.gut.0 = capscale(uUni.dist.g ~ 1, data = phyloseqCompanion::sample.data.frame(psGASV), dist = "bray", na.action = "na.omit")
plot(uUni.cap.gut)
anova(uUni.cap.gut)
uUni.cap.selected.gut = ordistep(uUni.cap.gut.0, Pin = 0.1, scope = formula(uUni.cap.gut), direction = 'both', permutations = 999)
uUni.cap.selected.gut$anova


bc.cap.nasal = capscale(bc.dist.n ~ Nose + Oral + Eyes + Condition(study), data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
bc.cap.nasal.0 = capscale(bc.dist.n ~ 1, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
plot(bc.cap.nasal)
anova(bc.cap.nasal)
bc.cap.selected.nasal = ordistep(bc.cap.nasal.0, Pin = 0.1, scope = formula(bc.cap.nasal), direction = 'both', permutations = 999)
bc.cap.selected.nasal$anova


jac.cap.nasal = capscale(jac.dist.n ~ Nose + Oral + Eyes + Condition(study) , data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
jac.cap.nasal.0 = capscale(jac.dist.n ~ 1, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
plot(jac.cap.nasal)
anova(jac.cap.nasal)
jac.cap.selected.nasal = ordistep(bc.cap.nasal.0, Pin = 0.1, scope = formula(jac.cap.nasal), direction = 'both', permutations = 999)
jac.cap.selected.nasal$anova

wUni.cap.nasal = capscale(wUni.dist.n ~ Nose + Oral + Eyes + Condition(study), data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
wUni.cap.nasal.0 = capscale(wUni.dist.n ~ 1, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
wUni.cap.selected.nasal = ordistep(wUni.cap.nasal.0, Pin = 0.1, scope = formula(bc.cap.nasal), direction = 'both', permutations = 999)
wUni.cap.selected.nasal$anova

uUni.cap.nasal = capscale(uUni.dist.n ~ Nose + Oral + Eyes + Condition(study), data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
uUni.cap.nasal.0 = capscale(uUni.dist.n ~ 1, data = phyloseqCompanion::sample.data.frame(psNASV), dist = "bray", na.action = "na.omit")
uUni.cap.selected.nasal = ordistep(uUni.cap.nasal.0, Pin = 0.1, scope = formula(bc.cap.nasal), direction = 'both', permutations = 999)
uUni.cap.selected.nasal$anova


##############################################################################################################
# 3. ALPHA DIVERSITY PLOTS
##############################################################################################################
## Nasal Signs
### Gut
library(phyloseqCompanion)
alphaColors = RColorBrewer::brewer.pal(n = 4, name = "Paired")
rust = "#B7410E"
gutAlpha = estimate_richness(psGASV, measures = c("Shannon", "Simpson", "Fisher"), split = TRUE) 
gutAlpha$ID = rownames(gutAlpha)
gutAlpha = gutAlpha %>% left_join(sample.data.frame(psGASV), by = c("ID" = "GutID")) %>%
  mutate(NasalCx = case_when(
    Nose == 1 ~ "Present",
    Nose == 0 ~ "Absent"
  ),
  OcularCx = case_when(
    Eyes == 1 ~ "Present",
    Eyes == 0 ~ "Absent"
  ))
nasalAlpha = estimate_richness(psNASV, measures = c("Shannon", "Simpson", "Fisher"), split = TRUE) 
nasalAlpha$ID = rownames(nasalAlpha)
nasalAlpha = nasalAlpha %>% left_join(sample.data.frame(psNASV), by = c("ID" = "NasalID")) %>%
  mutate(
    NasalCx = case_when(
      Nose == 1 ~ "Present",
      Nose == 0 ~ "Absent"
    ),
    OcularCx = case_when(
      Eyes == 1 ~ "Present",
      Eyes == 0 ~ "Absent"
    )
  )
library(ggpubr)
gutShannonNasal =  ggplot(aes(x = NasalCx, y = Shannon), data = gutAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Shannon, shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Shannon") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutShannonNasal

gutSimpsonNasal =  ggplot(aes(x = NasalCx, y = Simpson), data = gutAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Simpson, shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Simpson") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutSimpsonNasal

gutFischerNasal =  ggplot(aes(x = NasalCx, y = Fisher), data = gutAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Fisher, shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Fisher") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutFischerNasal

### Nasal
nasalShannonNasal =  ggplot(aes(x = NasalCx, y = Shannon), data = nasalAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Shannon, shape = case_control), color = "black", alpha = 0.8, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal Microbiome Shannon") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalShannonNasal

nasalSimpsonNasal =  ggplot(aes(x = NasalCx, y = Simpson), data = nasalAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Simpson, shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal  Microbiome Simpson") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalSimpsonNasal

nasalFischerNasal =  ggplot(aes(x = NasalCx, y = Fisher), data = nasalAlpha)+
  geom_boxplot(aes(fill = NasalCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  #geom_point(aes(x = Nose + 1, y = Fisher, shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Nasal Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal Microbiome Fisher") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalFischerNasal

lay = rbind(c(1, 1, 2, 2, 3, 3),
            c(1, 1, 2, 2, 3, 3),
            c(4, 4, 5, 5, 6, 6),
            c(4, 4, 5, 5, 6, 6))
pp = arrangeGrob(gutShannonNasal, gutSimpsonNasal, gutFischerNasal, nasalShannonNasal, nasalSimpsonNasal, nasalFischerNasal, layout_matrix = lay)
plot(pp)
setwd(directory.figures)
#ggsave(pp, file = "alphaDiversityNasalSigns.png", height = 10, width = 20)

## Ocular Signs
### Gut
gutShannonOcular =  ggplot(aes(x = OcularCx, y = Shannon), data = gutAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Shannon") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutShannonOcular

gutSimpsonOcular =  ggplot(aes(x = OcularCx, y = Simpson), data = gutAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Simpson") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutSimpsonOcular

gutFischerOcular =  ggplot(aes(x = OcularCx, y = Fisher), data = gutAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[1:2]) +
  ggtitle("Gut Microbiome Fisher") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
gutFischerOcular

### Nasal
nasalShannonOcular =  ggplot(aes(x = OcularCx, y = Shannon), data = nasalAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal Microbiome Shannon") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalShannonOcular

nasalSimpsonOcular =  ggplot(aes(x = OcularCx, y = Simpson), data = nasalAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal  Microbiome Simpson") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalSimpsonOcular

nasalFischerOcular =  ggplot(aes(x = OcularCx, y = Fisher), data = nasalAlpha)+
  geom_boxplot(aes(fill = OcularCx)) +
  geom_jitter(width = 0.2, aes(shape = case_control), color = "black", alpha = 0.6, size = 3) +
  labs(shape = "Case|Control", fill = "Ocular Cx") +
  scale_fill_manual(values = alphaColors[3:4]) +
  ggtitle("Nasal Microbiome Fisher") +
  theme(axis.title = element_text(size = 16), title = element_text(size = 18), axis.text.x = element_text(size = 15), legend.text = element_text(size = 15)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Absent", "Present")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
nasalFischerOcular

lay = rbind(c(1, 1, 2, 2, 3, 3),
            c(1, 1, 2, 2, 3, 3),
            c(4, 4, 5, 5, 6, 6),
            c(4, 4, 5, 5, 6, 6))
pp = arrangeGrob(gutShannonOcular, gutSimpsonOcular, gutFischerOcular, nasalShannonOcular, nasalSimpsonOcular, nasalFischerOcular, layout_matrix = lay)
plot(pp)
setwd(directory.figures)
#ggsave(pp, file = "alphaDiversityOccularSigns.png", height = 10, width = 20)


##############################################################################################################
# 4. INFLAMMATORY MARKERS
##############################################################################################################

wilcoxTestVariables = c("WBCTotal", "RBC", "Hemoglobin", "Hematocrit", "PCV", "MCV", "MCH", "MCHC",
                        "PlateletCount", "PlasmaProtein", "NumberNeutrophils", "NumberLymphocytes",
                        "NumberMonocytes", "NumberEosinophils", "NumberBands", "BUN", "Creatinine",
                        "Glucose", "Cholesterol", "TotalProtein", "Albumin", "BilirubinTotal",
                        "CK", "AlkalinePhosphatase", "ALT", "Sodium", "Potassium",
                        "Chloride", "Calcium") 

d.wilcox.vet = data.frame(matrix(nrow = length(wilcoxTestVariables), ncol = 2))
rownames(d.wilcox.vet) = wilcoxTestVariables
colnames(d.wilcox.vet) = c("statistic", "pValue")
metaG = metaG %>% mutate(
  status = case_when(
    Nose + Eyes + Oral + LungsRespiratory == 0 ~ "Control",
    Nose + Eyes + Oral + LungsRespiratory > 0 ~ "FURTD"
  )
)
metaN = metaN%>% mutate(
  status = case_when(
    Nose + Eyes + Oral + LungsRespiratory == 0 ~ "Control",
    Nose + Eyes + Oral + LungsRespiratory > 0 ~ "FURTD"
  )
)

for(i in 1:length(wilcoxTestVariables)){
  curTest = wilcoxTestVariables[i]
  
  x = metaG[which(metaG$status == "FURTD"), curTest]
  y = metaG[which(metaG$status == "Control"), curTest]
  print(paste0(c("Testing variable for significance vet: ", curTest, ":"), sep = "", collapse = ""))
  w = wilcox.test(x, y)
  d.wilcox.vet[curTest, "statistic"] = w$statistic
  d.wilcox.vet[curTest, "pValue"] = w$p.value
}

d.wilcox.vet$q = p.adjust(d.wilcox.vet$pValue, method = "fdr")
d.wilcox.vet[which(d.wilcox.vet$q < 0.5),]

p = ggplot(metaG, aes(x = status, y = NumberNeutrophils, color = status)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 12500, linetype = "dashed", color = "grey")+
  scale_color_manual(values = c("steelblue", rust)) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Neutrophils", x = "", y = "Number Neutrophils") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 30), axis.title = element_text(size = 22)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Control", "FURTD")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 

p = p + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p
setwd(directory.figures)
#png(file  = "AbsoluteNeutrophilNumber.png")
p
#dev.off()

p1 = ggplot(metaG, aes(x = status, y = Albumin, color = status)) +
  geom_boxplot(outlier.color = NA, outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 2.6, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey")+
  scale_color_manual(values = c("steelblue", rust)) + theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Albumin", x = "", y = "Albumin") +
  theme(axis.text = element_text(size = 18), plot.title = element_text(size = 30), axis.title = element_text(size = 22)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Control", "FURTD")), paired = FALSE, label.x.npc = "top", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) 
p1 = p1 + geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
p1

#png(file  = "Albumin.png")
p1
#dev.off()


##############################################################################################################
# 5. DETERMINE CLADES LINIKED TO CLINICAL SIGNS
##############################################################################################################
## GUT MICROBIOME

### Make a status variable
metaG = metaG %>% mutate(
  status = case_when(
    Nose + Eyes + Oral + LungsRespiratory == 0 ~ "Control",
    Nose + Eyes + Oral + LungsRespiratory > 0 ~ "FURTD"
  )
)
metaG = metaG %>% mutate(
  Status = case_when(status == "FURTD" ~ 1, status == "Control"~0)
)
metaN = metaN%>% mutate(
  status = case_when(
    Nose + Eyes + Oral + LungsRespiratory == 0 ~ "Control",
    Nose + Eyes + Oral + LungsRespiratory > 0 ~ "FURTD"
  )
)
metaN = metaN %>% mutate(
  Status = case_when(status == "FURTD" ~ 1, status == "Control"~0)
)

### Determine clades / ASVs that associate with Nasal / Ocular 
library(pscl)
qThresh = 0.05
prevGASV = apply(X = otu_table(psGASV), MARGIN = 1, FUN = function(x){sum(x > 0)})
psGASV.nb = phyloseq(otu_table(psGASV), tax_table(psGASV), sample_data(metaG), phy_tree(psGASV))
psGASV.nb = prune_taxa(x = psGASV.nb, taxa = names(prevGASV)[which(prevGASV > 6)]) # filter so that we only consider ASVs present in 6/15 of samples
nb.g.asv = getNBGLMS(covariates = c("Eyes", "Nose", "Status"), phyloseq = psGASV.nb)
dim(nb.g.asv[which(nb.g.asv$q < qThresh),])
sigGutASV = nb.g.asv[which(nb.g.asv$q < qThresh & nb.g.asv$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigGutASV

prevGCTU = apply(X = otu_table(psGCTU), MARGIN = 1, FUN = function(x){sum(x > 0)})
psGCTU.nb = phyloseq(otu_table(psGCTU), tax_table(psGCTU), sample_data(metaG))
psGCTU.nb = prune_taxa(x = psGCTU.nb, taxa = names(prevGCTU)[which(prevGCTU > 6)]) # filter so that we only consider ASVs present in 6/15 of samples (>40%)
toPrune = vector() # now get rid of taxa who have any sample sums equal to root.
for(i in 1:length(taxa_names(psGCTU.nb))){
  root = max(otu_table(psGCTU.nb)["root",])
  if(max(otu_table(psGCTU.nb)[i,]) == root){
    toPrune = c(taxa_names(psGCTU.nb)[i], toPrune)
  }
}
toPrune # gets rid of "root", "node2" "node3" and "node1" which cause GLMS to fail
psGCTU.nb = prune_taxa(x = psGCTU.nb, taxa = taxa_names(psGCTU.nb)[!taxa_names(psGCTU.nb) %in% toPrune]) # filter so that we only consider ASVs present in 6/15 of samples
nb.g.ctu = getNBGLMS(covariates = c("Nose", "Eyes", "Status"), phyloseq = psGCTU.nb)
dim(nb.g.ctu[which(nb.g.ctu$q < qThresh &nb.g.ctu$model %in%  c("glm.nb", "glm.nb.plus.01")),])
sigGutCTU = nb.g.ctu[which(nb.g.ctu$q < qThresh & nb.g.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigGutCTU


prevNASV = apply(X = otu_table(psNASV), MARGIN = 1, FUN = function(x){sum(x > 0)})
psNASV.nb = phyloseq(otu_table(psNASV), tax_table(psNASV), phy_tree(psNASV), sample_data(metaN))
psNASV.nb = prune_taxa(x = psNASV.nb, taxa = names(prevNASV)[which(prevNASV > 7)]) # filter so that we only consider ASVs present in 7/15 of samples
nb.n.asv = getNBGLMS(covariates = c("Nose", "Eyes", "Status"), phyloseq = psNASV.nb)
dim(nb.n.asv[which(nb.n.asv$q < qThresh & nb.n.asv$model == "NB"),])
sigNasalASV = nb.n.asv[which(nb.n.asv$q < qThresh & nb.n.asv$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigNasalASV

prevNCTU = apply(X = otu_table(psNCTU), MARGIN = 1, FUN = function(x){sum(x > 0)})
psNCTU.nb = phyloseq(otu_table(psNCTU), tax_table(psNCTU), sample_data(metaN))
psNCTU.nb = prune_taxa(x = psNCTU.nb, taxa = names(prevNCTU)[which(prevNCTU > 7)]) # filter so that we only consider ASVs present in 7/15 of samples
toPrune = vector() # now get rid of taxa who have any sample sums equal to root.
for(i in 1:length(taxa_names(psNCTU.nb))){
  root = max(otu_table(psNCTU.nb)["root",])
  if(max(otu_table(psNCTU.nb)[i,]) == root){
    toPrune = c(taxa_names(psNCTU.nb)[i], toPrune)
  }
}
toPrune # gets rid of "root", "node4", "node5", and "node3". Nodes close to root are failing GLMS, but we expect their distribution to be close to uniform anyway
psNCTU.nb = prune_taxa(x = psNCTU.nb, taxa = taxa_names(psNCTU.nb)[!taxa_names(psNCTU.nb) %in% toPrune]) # filter so that we only consider ASVs present in 6/15 of samples
nb.n.ctu = getNBGLMS(covariates = c("Nose", "Eyes", "Status"), phyloseq = psNCTU.nb)
dim(nb.n.ctu[which(nb.n.ctu$q < qThresh & nb.n.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),])
sigNasalCTU = nb.n.ctu[which(nb.n.ctu$q < qThresh & nb.n.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigNasalCTU


# Visulaize clades linked to clincial signs 
## Combine data into one frame
sigGutASV$type = "ASV"
sigGutCTU$type = "CTU"
sigGut = rbind(sigGutASV, sigGutCTU)
sigGut$microbiome = "Gut"

sigNasalASV$type = "ASV"
sigNasalCTU$type = "CTU"
sigNasal = rbind(sigNasalASV, sigNasalCTU)
sigNasal$microbiome = "Nasal"

sig = rbind(sigGut, sigNasal)
sig = as_tibble(sig)
sig
library(phangorn)

treeGLevel = getRefTreeCladeLevel(tree = treeG)
rownames(treeGLevel) = treeGLevel$node
rootsG = RemoveNestedClades(cladeList = treeGLevel[unique(sigGutCTU %>% filter(covariate %in% c("Status"))) %>% pull(clade) %>% unique(),], tree = treeG)
rootsG
taxGAll = as_tibble(rbind(taxa.data.table(ps = psGASV), taxa.data.table(ps = (psGCTU))))
taxGAll %>% filter(Taxon %in% rootsG) %>% print(n = 35)
as_tibble(sigGut) %>% filter(covariate == "Status") %>% left_join(taxGAll, by = c("clade" = "Taxon")) %>% arrange(Family) %>% print(n = 107)

dG = treeDataFrame("Status", tree = treeG)
head(dG)
dG$nodeName = rownames(dG)
dG = dG %>% left_join(taxGAll, by = c("nodeName" = "Taxon"))
rownames(dG) = dG$nodeName
head(dG)
dG[sigGut[which(sigGut$covariate == "Status"), "clade"], "Status"] = sigGut[which(sigGut$covariate == "Status"), "Estimate"]
dG = dG %>% 
  mutate(
    Status = case_when(
      Status < 0 ~ "-",
      Status >= 0 ~ "+"
    ))

dG
dG %>% filter(nodeName %in% rootsG) 

palEarth = RColorBrewer::brewer.pal(n = 11, "BrBG")
idxsG = c(treeG$tip.label, treeG$node.label)
paste(c(
  cladeLabels(dataFrame = dG %>% filter(nodeName %in% rootsG) %>% filter(Status == "+"), idxs = "idxsG", nodeName = "nodeName", color = "rust", fontsize = "10"),
  cladeLabels(dataFrame = dG %>% filter(!is.na(Status)) %>% filter(nodeName %in% rootsG) %>%  filter(Status == "-"), idxs = "idxsG",  nodeName = "nodeName", color = "steelblue", fontsize = "10")
  
), sep = "", collapse = "")

cxGutTree = ggtree(treeG, layout = "rectangular", size = 0.15) %<+% dG
fsize = 4
cxGutTree = cxGutTree +  
  geom_point2(aes(color = Phylum), size = 2, alpha = .5) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 18), legend.position = "left") +
  geom_point2(aes(subset = (!is.na(Status)), fill = Status), alpha = 1, size = 5, shape = 21) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_fill_manual(values = c(steelblue, rust)) +  
  geom_cladelabel(node = which(idxsG=='node5'), label = 'Bacteria', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node84'), label = 'Fusobacterium', align = F, geom = 'label', fill = rust, fontsize = fsize) +
  geom_cladelabel(node = which(idxsG=='node186'), label = 'Bacteroidales', align = F, geom = 'label', fill = rust, fontsize = fsize) +
  geom_cladelabel(node = which(idxsG=='node241'), label = 'Marinifilaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) +
  geom_cladelabel(node = which(idxsG=='node264'), label = 'Bacteroides', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node276'), label = 'Bacteroides', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node284'), label = 'Bacteroides', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node293'), label = 'Parabacteroides', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node336'), label = 'Alistipes', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node462'), label = 'Enterobacteriaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node508'), label = 'Betaproteobacteriales', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node536'), label = 'Clostridium_sensu_stricto_1', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node556'), label = 'Peptostreptococcus', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node562'), label = 'Peptostreptococcaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node600'), label = 'Peptoniphilus', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node617'), label = 'Finegoldia', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node625'), label = 'Bacteria', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node678'), label = 'Anaerostipes', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node684'), label = 'Lachnospiraceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node791'), label = 'Lachnospiraceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node873'), label = 'Subdoligranulum', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node889'), label = 'Ruminococcaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  #geom_cladelabel(node = which(idxsG=='node927'), label = 'Streptococcus', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node934'), label = 'Lactobacillaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node216'), label = 'Prevotella_9', align = F, geom = 'label', fill = steelblue, fontsize = fsize) +
  #geom_cladelabel(node = which(idxsG=='node218'), label = 'Prevotella_9', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node281'), label = 'Bacteroides', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node294'), label = 'Porphyromonas', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node366'), label = 'Desulfovibrionales', align = F, geom = 'label', fill = steelblue, fontsize = fsize) +
  geom_cladelabel(node = which(idxsG=='node387'), label = 'Bacteria', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node444'), label = 'Gammaproteobacteria', align = F, geom = 'label', fill = steelblue, fontsize = fsize) +
  geom_cladelabel(node = which(idxsG=='node652'), label = 'Peptococcus', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node708'), label = 'Clostridiales', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node782'), label = 'Lachnospiraceae', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node804'), label = 'Ruminococcaceae_UCG-014', align = F, geom = 'label', fill = steelblue, fontsize = fsize) +
  #geom_cladelabel(node = which(idxsG=='node830'), label = 'Ruminococcaceae', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + 
  geom_cladelabel(node = which(idxsG=='node836'), label = 'Ruminococcaceae_UCG-004', align = F, geom = 'label', fill = steelblue, fontsize = fsize)


setwd(directory.figures)
pdf(file = "statusGutMap.pdf", height = 10, width = 10)
cxGutTree
dev.off()



tipsG = vector()
for(i in 1:length(rootsG)){
  curNode = idxsG[which(idxsG == rootsG[i])]
  curClade = extract.clade(treeG, curNode)
  tipsG = c(tipsG, curClade$tip.label)
  tipsG = unique(tipsG)
}
tipsG = unique(c(tipsG, sigGut %>% filter(type == "ASV") %>% filter(covariate %in% c("Status")) %>% pull(clade)))

treeGSub = ape::drop.tip(phy = treeG, tip = setdiff(treeG$tip.label, tipsG))  
dGSub = dG[c(treeGSub$tip.label, treeGSub$node.label),]
dGSub$node = seq(from = 1, to = length(treeGSub$tip.label) + length(treeGSub$node.label), by = 1)
idxsG = seq(from = 1, to = length(treeGSub$tip.label) + length(treeGSub$node.label), by = 1)
names(idxsG) = c(treeGSub$tip.label, treeGSub$node.label)

nodeids = idxsG[rootsG]
nodedf <- data.frame(node=nodeids,
                     clade = names(idxsG)[nodeids],
                     pos = c(rep(0.2, times = length(nodeids))))

nodedf = nodedf %>% left_join(as_tibble(taxGAll), by = c("clade" = "Taxon"))
nodedf = nodedf %>% 
  mutate(label = case_when(
    is.na(Phylum) ~ Kingdom,
    is.na(Class) ~ Phylum,
    is.na(Order)~ Class,
    is.na(Family)~ Order,
    is.na(Genus) ~ Family,
    is.na(Species) ~ Genus
    
  )
  )

# Move labels outward from smaller clades
nodedf$size = NA
for(i in 1:nrow(nodedf)){
  curClade = extract.clade(treeGSub, node = nodedf[i, "node"])
  nodedf[i, "size"] = length(curClade$tip.label)
}
nodedf = nodedf %>% 
  mutate(
    pos = case_when(
      size < 5 ~ .7,
      size < 15 ~ .6,
      size >= 15 ~ .2
    )
  )

lac = extract.clade(phy = treeGSub, node = idxsG["node934"] )
as_tibble(taxGAll) %>% filter(Taxon %in% lac$tip.label) %>% print(n = 100)


cxGutTree = ggtree(treeGSub, layout = "fan", size = 0.2, open.angle = 5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=2.2, alpha=0.2, fill="grey", color="grey50",size=0.05) + 
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=3.3, alpha=0.3, fill="white", color="grey50",size=0.05) + 
  geom_cladelab(data=nodedf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=TRUE, 
                fontsize=1.7,
                fontface="italic")
cxGutTree

colnames(dGSub)[which(colnames(dGSub) == "Status" )] = "Association"
cxGutTree = cxGutTree %<+% dGSub + 
  geom_point2(aes(color = Family), alpha = 0.5, size = 1) +
  geom_star(mapping=aes(fill=Association, starshape=Association), position="identity", starstroke=0.1, size = 2, alpha = 1) +
  scale_fill_manual(values = c(steelblue, rust)) +
  guides(starshape = guide_legend(override.aes = list(size = 8))) 
cxGutTree


tmpPhyloseq = phyloseq(otu_table(psGASV), phy_tree(psGASV), tax_table(psGASV), sample_data(metaG))
tmpPhyloseq = prune_taxa(taxa = treeGSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = prune_samples(samples = rownames(sample_data(tmpPhyloseq)[which(sample_data(tmpPhyloseq)$status == "FURTD"),]), x = tmpPhyloseq)
furtdAve = apply(otu_table(tmpPhyloseq)[treeGSub$tip.label,], 1, mean)


furtdAve = data.frame("Average" = log(furtdAve + 0.001),
                      "Cx" = rep("FURTD", times = length(furtdAve)))
furtdAve$ID = rownames(furtdAve)

tmpPhyloseq = phyloseq(otu_table(psGASV), phy_tree(psGASV), tax_table(psGASV), sample_data(metaG))
tmpPhyloseq = prune_taxa(taxa = treeGSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = prune_samples(samples = rownames(sample_data(tmpPhyloseq)[which(sample_data(tmpPhyloseq)$status == "Control"),]), x = tmpPhyloseq)
controlAve = apply(otu_table(tmpPhyloseq)[treeGSub$tip.label,], 1, mean)


controlAve = data.frame("Average" = log(controlAve + 0.001),
                        "Cx" = rep("Control", times = length(controlAve)))
controlAve$ID = rownames(controlAve)
datAve = rbind(furtdAve, controlAve)


cxGutTree = cxGutTree + 
  new_scale_fill() +
  geom_fruit(data=datAve, geom=geom_tile,
             mapping=aes(y=ID, x=Cx, fill=Cx, alpha = Average),
             color = "grey50", offset = 0.45,size = 0.01) +
  scale_fill_manual(values = c(steelblue, rust)) 
cxGutTree

# Make a higher abundance table
colnames(furtdAve) = c("AveF", "FStatus", "ID")
colnames(controlAve) = c("AveC", "CStatus", "ID")
scaleBarsToZero = min(min(furtdAve$AveF), min(controlAve$AveC))

dat3 = as_tibble(furtdAve) %>% 
  left_join(controlAve, by = "ID") %>% 
  mutate(
    HigherAve = case_when(
      AveF > AveC ~ AveF - scaleBarsToZero,
      AveF < AveC ~ AveC - scaleBarsToZero,
      AveF == AveC ~ 0),
    HigherClass = case_when(
      AveF > AveC ~ "FURTD",
      AveF < AveC ~ "Control",
      AveF == AveC ~ "Tie"
    )
    
  )
dat3 = as.data.frame(dat3)
as.data.frame(dat3)
rownames(dat3) = dat3$ID
dat3$entropy = NA

for(i in 1:nrow(nodedf)){
  curClade = extract.clade(treeGSub, node = nodedf[i, "node"])
  curDat =  dat3 %>% filter(ID %in% curClade$tip.label)
  curRatio = table(curDat$HigherClass)
  if(length(curRatio) == 1){
    dat3[rownames(curDat),"entropy"] = 1
  }else if(as.vector(curRatio["Control"]) == as.vector(curRatio["FURTD"])){
    dat3[rownames(curDat),"entropy"] = 0.5
  }else if(as.vector(curRatio["Control"]) > as.vector(curRatio["FURTD"])){
    dat3[rownames(curDat), "entropy"] = as.vector(curRatio["Control"]) / sum(as.vector(curRatio["Control"]) + as.vector(curRatio["FURTD"]))
  }else if(as.vector(curRatio["Control"]) < as.vector(curRatio["FURTD"])){
    dat3[rownames(curDat), "entropy"] = as.vector(curRatio["FURTD"]) / sum(as.vector(curRatio["Control"]) + as.vector(curRatio["FURTD"]))
  }
  else{
    print("unexpected output found")
  }
}

dat3[is.na(dat3$entropy), "entropy"] = 1
dat3 = as_tibble(dat3)
dat3 = dat3 %>% mutate(
  entropy = case_when(entropy <= 0.5 ~ 50,
                      entropy < 0.7 ~ 70,
                      entropy < 0.9 ~ 90,
                      entropy <=1 ~ 100)
)
cxGutTree = cxGutTree + 
  geom_fruit(data=dat3, geom=geom_bar, mapping=aes(y=ID, x=HigherAve, fill = HigherClass) , pwidth=0.38, offset = .1,orientation="y",  stat="identity") 
cxGutTree = cxGutTree + 
  guides(color = guide_legend(override.aes = list(size=8))) +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size = 12), legend.position = "right") +
  guides(alpha=guide_legend(title="Log Average\nAbundance"))

cxGutTreeLegend = get_legend(cxGutTree)

lay = rbind(c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1))


# Arrange plots
pp = arrangeGrob(cxGutTree + theme(legend.position = "none"), layout_matrix = lay)
plot(pp)
ggsave(pp, file = "cladesClinicalSignsReducedGut.pdf", height = 10, width = 10, units = "in")

lay = rbind(c( 1, 1, 1),
            c(1, 1, 1),
            c(1, 1, 1),
            c(2, 3, 4))

pp = arrangeGrob(cxGutTreeLegend$grobs[[4]], cxGutTreeLegend$grobs[[1]], cxGutTreeLegend$grobs[[2]], cxGutTreeLegend$grobs[[3]], layout_matrix = lay)
plot(pp)
ggsave(pp, file = "cladesClinicalSignsReducedGutKey.pdf", height = 9, width = 5)

### Nasal Tree
treeNLevel = getRefTreeCladeLevel(tree = treeN)
rownames(treeNLevel) = treeNLevel$node
rootsN = RemoveNestedClades(cladeList = treeGLevel[unique(sigNasalCTU %>% filter(covariate %in% c("Status"))) %>% pull(clade) %>% unique(),], tree = treeN)
rootsN
taxNAll = as_tibble(rbind(taxa.data.table(ps = psNASV), taxa.data.table(ps = (psNCTU))))
taxNAll %>% filter(Taxon %in% rootsG) %>% print(n = 35)
as_tibble(sigNasal) %>% filter(covariate == "Status") %>% left_join(taxNAll, by = c("clade" = "Taxon")) %>% arrange(Family) %>% print(n = 107)

dN = treeDataFrame("Status", tree = treeN)
head(dN)
dN$nodeName = rownames(dN)
dN = dN %>% left_join(taxNAll, by = c("nodeName" = "Taxon"))
rownames(dN) = dN$nodeName
head(dN)
dN[sigNasal[which(sigNasal$covariate == "Status"), "clade"], "Status"] = sigNasal[which(sigNasal$covariate == "Status"), "Estimate"]
dN = dN %>% 
  mutate(
    Status = case_when(
      Status < 0 ~ "-",
      Status >= 0 ~ "+"
    ))

dN
dN %>% filter(nodeName %in% rootsN) 
idxsN = c(treeN$tip.label, treeN$node.label)
paste(c(
  cladeLabels(dataFrame = dN %>% filter(nodeName %in% rootsN) %>% filter(Status == "+"), idxs = "idxsN", nodeName = "nodeName", color = "rust", fontsize = "10"),
  cladeLabels(dataFrame = dN %>% filter(!is.na(Status)) %>% filter(nodeName %in% rootsN) %>%  filter(Status == "-"), idxs = "idxsN",  nodeName = "nodeName", color = "steelblue", fontsize = "10")
  
), sep = "", collapse = "")

cxNasalTree = ggtree(treeN, layout = "rectangular", size = 0.15) %<+% dN
fsize = 4
cxNasalTree = cxNasalTree +  
  geom_point2(aes(color = Phylum), size = 2, alpha = .5) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 18), legend.position = "left") +
  geom_point2(aes(subset = (!is.na(Status)), fill = Status), alpha = 1, size = 5, shape = 21) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_fill_manual(values = c(steelblue, rust)) +geom_cladelabel(node = which(idxsN=='node31'), label = 'Actinomyces', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node150'), label = 'Catonella', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node168'), label = 'Ezakiella', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node281'), label = 'Streptococcus', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node385'), label = 'Leptotrichiaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node396'), label = 'Fusobacteriaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node452'), label = 'Treponema_2', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node489'), label = 'Porphyromonas', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node572'), label = 'Capnocytophaga', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node664'), label = 'Alloprevotella', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node757'), label = 'Bacteria', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node864'), label = 'Gammaproteobacteria', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node937'), label = 'Pasteurellaceae', align = F, geom = 'label', fill = rust, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node126'), label = 'Johnsonella', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node200'), label = 'Fusibacter', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node218'), label = 'Proteocatella', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node231'), label = 'Clostridiales', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node293'), label = 'MVP-15', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node373'), label = 'JGI_0000069-P22', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node428'), label = 'Treponema_2', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node553'), label = 'Bergeyella', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node575'), label = 'Flavobacteriaceae', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node634'), label = 'Porphyromonas', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node677'), label = 'Alloprevotella', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node730'), label = 'Absconditabacteriales_(SR1)', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node884'), label = 'Moraxella', align = F, geom = 'label', fill = steelblue, fontsize = fsize) + geom_cladelabel(node = which(idxsN=='node931'), label = 'Bibersteinia', align = F, geom = 'label', fill = steelblue, fontsize = fsize)

cxNasalTree


setwd(directory.figures)
pdf(file = "statusNasalMap.pdf", height = 18, width = 10)
cxNasalTree
dev.off()

tipsN = vector()
for(i in 1:length(rootsN)){
  curNode = idxsN[which(idxsN == rootsN[i])]
  curClade = extract.clade(treeN, curNode)
  tipsN = c(tipsN, curClade$tip.label)
  tipsN = unique(tipsN)
}
tipsN = unique(c(tipsN, sigNasal %>% filter(type == "ASV") %>% filter(covariate %in% c("Status")) %>% pull(clade)))

treeNSub = ape::drop.tip(phy = treeN, tip = setdiff(treeN$tip.label, tipsN))  
dNSub = dN[c(treeNSub$tip.label, treeNSub$node.label),]
dNSub$node = seq(from = 1, to = length(treeNSub$tip.label) + length(treeNSub$node.label), by = 1)
idxsN = seq(from = 1, to = length(treeNSub$tip.label) + length(treeNSub$node.label), by = 1)
names(idxsN) = c(treeNSub$tip.label, treeNSub$node.label)

nodeids = idxsN[rootsN]
nodedf <- data.frame(node=nodeids,
                     clade = names(idxsN)[nodeids],
                     pos = c(rep(0.2, times = length(nodeids))))

nodedf = nodedf %>% left_join(as_tibble(taxNAll), by = c("clade" = "Taxon"))
nodedf = nodedf %>% 
  mutate(label = case_when(
    is.na(Phylum) ~ Kingdom,
    is.na(Class) ~ Phylum,
    is.na(Order)~ Class,
    is.na(Family)~ Order,
    is.na(Genus) ~ Family,
    is.na(Species) ~ Genus
    
  )
  )

# Move labels outward from smaller clades
nodedf$size = NA
for(i in 1:nrow(nodedf)){
  curClade = extract.clade(treeNSub, node = nodedf[i, "node"])
  nodedf[i, "size"] = length(curClade$tip.label)
}
nodedf = nodedf %>% 
  mutate(
    pos = case_when(
      size < 5 ~ .1,
      size < 15 ~ .1,
      size >= 15 ~ .1
    )
  )


cxNasalTree = ggtree(treeNSub, layout = "fan", size = 0.2, open.angle = 5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=2.3, alpha=0.2, fill="grey", color="grey50",size=0.05) + 
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=3.5, alpha=0.3, fill="white", color="grey50",size=0.05) + 
  geom_cladelab(data=nodedf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=TRUE, 
                fontsize=1.7,
                fontface="italic")
cxNasalTree

colnames(dNSub)[which(colnames(dNSub) == "Status" )] = "Association"
cxNasalTree = cxNasalTree %<+% dNSub + 
  geom_point2(aes(color = Family), alpha = 0.5, size = 1) +
  geom_star(mapping=aes(fill=Association, starshape=Association), position="identity", starstroke=0.1, size = 2, alpha = 1) +
  scale_fill_manual(values = c(steelblue, rust)) +
  guides(starshape = guide_legend(override.aes = list(size = 8))) 
cxNasalTree


tmpPhyloseq = phyloseq(otu_table(psNASV), phy_tree(psNASV), tax_table(psNASV), sample_data(metaN))
tmpPhyloseq = prune_taxa(taxa = treeNSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = prune_samples(samples = rownames(sample_data(tmpPhyloseq)[which(sample_data(tmpPhyloseq)$status == "FURTD"),]), x = tmpPhyloseq)
furtdAve = apply(otu_table(tmpPhyloseq)[treeNSub$tip.label,], 1, mean)


furtdAve = data.frame("Average" = log(furtdAve + 0.001),
                      "Cx" = rep("FURTD", times = length(furtdAve)))
furtdAve$ID = rownames(furtdAve)

tmpPhyloseq = phyloseq(otu_table(psNASV), phy_tree(psNASV), tax_table(psNASV), sample_data(metaN))
tmpPhyloseq = prune_taxa(taxa = treeNSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = prune_samples(samples = rownames(sample_data(tmpPhyloseq)[which(sample_data(tmpPhyloseq)$status == "Control"),]), x = tmpPhyloseq)
controlAve = apply(otu_table(tmpPhyloseq)[treeNSub$tip.label,], 1, mean)


controlAve = data.frame("Average" = log(controlAve + 0.001),
                        "Cx" = rep("Control", times = length(controlAve)))
controlAve$ID = rownames(controlAve)
datAve = rbind(furtdAve, controlAve)


cxNasalTree = cxNasalTree + 
  new_scale_fill() +
  geom_fruit(data=datAve, geom=geom_tile,
             mapping=aes(y=ID, x=Cx, fill=Cx, alpha = Average),
             color = "grey50", offset = 0.31,size = 0.01) +
  scale_fill_manual(values = c(steelblue, rust, "grey")) 
cxNasalTree

# Make a higher abundance table
colnames(furtdAve) = c("AveF", "FStatus", "ID")
colnames(controlAve) = c("AveC", "CStatus", "ID")
scaleBarsToZero = min(min(furtdAve$AveF), min(controlAve$AveC))

dat3 = as_tibble(furtdAve) %>% 
  left_join(controlAve, by = "ID") %>% 
  mutate(
    HigherAve = case_when(
      AveF > AveC ~ AveF - scaleBarsToZero,
      AveF < AveC ~ AveC - scaleBarsToZero,
      AveF == AveC ~ 0),
    HigherClass = case_when(
      AveF > AveC ~ "FURTD",
      AveF < AveC ~ "Control",
      AveF == AveC ~ "Tie"
    )
    
  )



cxNasalTree = cxNasalTree + 
  geom_fruit(data=dat3, geom=geom_bar, mapping=aes(y=ID, x=HigherAve, fill = HigherClass) , pwidth=0.38, offset = .1,orientation="y",  stat="identity") 
cxNasalTree

cxNasalTree = cxNasalTree + 
  guides(color = guide_legend(override.aes = list(size=8))) +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size = 12), legend.position = "right") +
  guides(alpha=guide_legend(title="Log Average\nAbundance"))

cxNasalTreeLegend = get_legend(cxNasalTree)

lay = rbind(c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1), 
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1),
            c(1, 1, 1, 1, 1, 1))


# Arrange plots
pp = arrangeGrob(cxNasalTree + theme(legend.position = "none"), layout_matrix = lay)
plot(pp)
ggsave(pp, file = "cladesClinicalSignsReducedNasal.pdf", height = 10, width = 10, units = "in")

lay = rbind(c(6, 6,  1, 1),
            c(2, 2, 1, 1),
            c(3, 3, 1, 1),
            c(4, 4, 1, 1), 
            c(5, 5, 1, 1))

pp = arrangeGrob(cxNasalTreeLegend$grobs[[4]], cxNasalTreeLegend$grobs[[1]],  cxNasalTreeLegend$grobs[[3]], cxNasalTreeLegend$grobs[[2]],  layout_matrix = lay)
plot(pp)
ggsave(pp, file = "cladesClinicalSignsReducedNasalKey.pdf", height = 9, width = 5)



# See if the significant clades are clustered.
sigCladesIdxG = seq(from = 1, to = length(c(treeG$tip.label, treeG$node.label)), by = 1)
names(sigCladesIdxG) = c(treeG$tip.label, treeG$node.label)
sigCladesIdxN = seq(from = 1, to = length(c(treeN$tip.label, treeN$node.label)), by = 1)
names(sigCladesIdxN) = c(treeN$tip.label, treeN$node.label)


sigCladesPositiveG = sigCladesIdxG[sigGut %>% filter(covariate %in% c("Nose", "Eyes")) %>% filter(Estimate > 0) %>% pull(clade)]
sigCladesNegativeG = sigCladesIdxG[sigGut %>% filter(covariate %in% c("Nose", "Eyes")) %>%filter(Estimate < 0) %>% pull(clade)]

sigCladesPositiveN = sigCladesIdxN[sigNasal %>% filter(covariate %in% c("Nose", "Eyes")) %>% filter(Estimate > 0) %>% pull(clade)]
sigCladesNegativeN = sigCladesIdxN[sigNasal %>% filter(covariate %in% c("Nose", "Eyes")) %>% filter(Estimate < 0) %>% pull(clade)]

combCladesPositiveG = combn(x = sigCladesPositiveG, m = 2)
combCladesNegativeG = combn(x = sigCladesNegativeG, m = 2)

combCladesPositiveN = combn(x = sigCladesPositiveN, m = 2)
combCladesNegativeN = combn(x = sigCladesNegativeN, m = 2)

distG = dist.nodes(treeG)
distN = dist.nodes(treeN)

distPositiveG = getDistanceCombos(combos = combCladesPositiveG, distMatrix = distG)
distNegativeG = getDistanceCombos(combos = combCladesNegativeG, distMatrix = distG)
distRandG = distRandom(dist = distG, n = length(distPositiveG))

distPositiveN = getDistanceCombos(combos = combCladesPositiveN, distMatrix = distN)
distNegativeN = getDistanceCombos(combos = combCladesNegativeN, distMatrix = distN)
distRandN = distRandom(dist = distN, n = length(distPositiveN))

ks.test(x = distPositiveG, y = distRandG, alternative = "two.sided")
ks.test(x = distNegativeG, y = distRandG, alternative = "two.sided")
ks.test(x = distNegativeG, y = distPositiveG, alternative = "two.sided")

ks.test(x = distPositiveN, y = distRandN, alternative = "two.sided")
ks.test(x = distNegativeN, y = distRandN, alternative = "two.sided")
ks.test(x = distNegativeN, y = distPositiveN, alternative = "two.sided")

wilcox.test(distRandG, distNegativeG, method = "two.sided")
wilcox.test(distRandG, distPositiveG, method = "two.sided")
wilcox.test(x = distNegativeG, y = distPositiveG, method = "two.sided")
kruskal.test(formula = dist~group, 
             data = data.frame("dist" = c(distRandG, distNegativeG, distPositiveG), 
                               "group" = c(rep("rand", times = length(distRandG)), 
                                           rep("negative", times = length(distNegativeG)), 
                                           rep("positive", times = length(distPositiveG)))))

wilcox.test(distRandN, distNegativeN, method = "two.sided")
wilcox.test(distRandN, distPositiveN, method = "two.sided")
wilcox.test(x = distNegativeN, y = distPositiveN, method = "two.sided")
kruskal.test(formula = dist~group, 
             data = data.frame("dist" = c(distRandN, distNegativeN, distPositiveN), 
                               "group" = c(rep("rand", times = length(distRandN)), 
                                           rep("negative", times = length(distNegativeN)), 
                                           rep("positive", times = length(distPositiveN)))))


# Stats on average pairwise phylogenetic distances.
# Stats on average pairwise phylogenetic distances.
mean(distRandG)
sd(distRandG)
mean(distPositiveG)
sd(distPositiveG)
mean(distNegativeG)
sd(distNegativeG)

mean(distRandN)
sd(distRandN)
mean(distPositiveN)
sd(distPositiveN)
mean(distNegativeN)
sd(distNegativeN)

as_tibble(sigGut) %>% filter(covariate == "Status")
as_tibble(sigGut) %>% filter(covariate == "Status") %>%left_join(taxGAll, by = c("clade" = "Taxon")) %>% filter(Estimate <0) %>% arrange(Family) %>% print(n = 123) 
as_tibble(sigGut) %>% filter(covariate == "Status") %>% left_join(taxGAll, by = c("clade" = "Taxon")) %>% filter(Estimate <0) %>% arrange(Family) %>% group_by(Family) %>% tally() %>% print(n = 123) 

as_tibble(sigGut) %>% filter(covariate == "Status") %>%left_join(taxGAll, by = c("clade" = "Taxon")) %>% filter(Estimate >0) %>% arrange(Family) %>% print(n = 123) 
as_tibble(sigGut) %>% filter(covariate == "Status") %>% left_join(taxGAll, by = c("clade" = "Taxon")) %>% filter(Estimate >0) %>% arrange(Family) %>% group_by(Family) %>% tally() %>% print(n = 123) 
taxGAll %>% filter(Taxon %in% rootsG) %>% arrange(Class) %>% print( n = 37) 


as_tibble(sigNasal) %>% filter(covariate == "Status")
as_tibble(sigNasal) %>% filter(covariate == "Status") %>%left_join(taxNAll, by = c("clade" = "Taxon")) %>% filter(Estimate <0) %>% arrange(Kingdom, Phylum, Class, Order, Family)%>% print(n = 123) 
as_tibble(sigNasal) %>% filter(covariate == "Status") %>% left_join(taxNAll, by = c("clade" = "Taxon")) %>% filter(Estimate <0) %>% arrange(Family) %>% group_by(Kingdom, Phylum, Class, Order, Family) %>% tally() %>% arrange(desc(n)) %>% print(n = 123) 

as_tibble(sigNasal) %>% filter(covariate == "Status") %>%left_join(taxNAll, by = c("clade" = "Taxon")) %>% filter(Estimate >0) %>% arrange(Kingdom, Phylum, Class, Order, Family)%>% print(n = 123) 
as_tibble(sigNasal) %>% filter(covariate == "Status") %>% left_join(taxNAll, by = c("clade" = "Taxon")) %>% filter(Estimate >0) %>% arrange(Family) %>% group_by(Kingdom, Phylum, Class, Order, Family) %>% tally() %>% arrange(desc(n)) %>% print(n = 123) 
taxNAll %>% filter(Taxon %in% rootsN) %>% arrange(Family) %>% print( n = 37) 


gamma = extract.clade(treeG, which(treeG$node.label == "node444") + length(treeG$tip.label)) 
taxGAll %>% filter(Taxon %in% gamma$tip.label) %>% print(n = 103)

morax = extract.clade(treeN, which(treeN$node.label == "node884") + length(treeN$tip.label)) 
taxNAll %>% filter(Taxon %in% morax$tip.label)
apply(otu_table(psNASV)[morax$tip.label,rownames(metaN)[which(metaN$status == "Control")]], 1, mean)
apply(otu_table(psNASV)[morax$tip.label,rownames(metaN)[which(metaN$status == "FURTD")]], 1, mean)

gammaNose = extract.clade(treeN, which(treeN$node.label == "node864") + length(treeN$tip.label)) 
taxNAll %>% filter(Taxon %in% gammaNose$tip.label)

fusoNose = extract.clade(treeN, which(treeN$node.label == "node396") + length(treeN$tip.label)) 
taxNAll %>% filter(Taxon %in% fusoNose$tip.label)


pal = RColorBrewer::brewer.pal(n = 4, name = "Dark2")
dG = rbind(data.frame(cbind("Dist" = distPositiveG, 
                            "Node" = rep("+", times = length(distPositiveG)))), 
           data.frame(cbind("Dist" = distNegativeG, 
                            "Node" = rep("-", times = length(distNegativeG)))), 
           data.frame(cbind("Dist" = distRandG, "Node" = rep("Random", times = length(distRandG)))))
dG$Dist = as.numeric(dG$Dist)

setwd(directory.figures)
pdf(file = "densityPlotPhylogeneticClustersGut.pdf")
densityGutPerm = ggplot(data = dG, aes(x = Dist, group = Node, fill = Node)) +
  scale_fill_manual(values=c("blue", "red", "grey")) +
  geom_density(adjust = 1.5, alpha = 0.5) +
  labs(x = "Phylogenetic Distance", y = "Phylogenetic Distance") +
  geom_boxplot(width=0.1, outlier.size = .05) + 
  ggtitle(label = "Pairwise Phylogenetic Distance")
densityGutPerm
dev.off()

library(ggpubr)
setwd(directory.figures)
pdf(file = "violinPlotPhylogeneticClustersGut.pdf")
violinGutPerm = ggplot(dG, aes(x=Node, y=Dist, fill = Node)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(aes(fill = Node), width=0.1, fill = "white", outlier.alpha = .1, outlier.size = .5) + 
  scale_fill_manual(values=c("blue", "red", "grey")) +
  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  labs(x = element_blank(), y = "Phylogenetic Distance") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("+", "-"), c("+", "Random"), c("-", "Random"))) +
  stat_compare_means(method = "kruskal.test", label.y = 4)
violinGutPerm
dev.off()



dN = rbind(data.frame(cbind("Dist" = distPositiveN, 
                            "Node" = rep("+", times = length(distPositiveN)))), 
           data.frame(cbind("Dist" = distNegativeN, 
                            "Node" = rep("-", times = length(distNegativeN)))), 
           data.frame(cbind("Dist" = distRandN, 
                            "Node" = rep("Random", times = length(distRandN)))))
dN$Dist = as.numeric(dN$Dist)

setwd(directory.figures)
pdf(file = "densityPlotPhylogeneticClustersNasal.pdf")
densityNasalPerm = ggplot(data = dN, aes(x = Dist, group = Node, fill = Node)) +
  scale_fill_manual(values=c("blue", "red", "grey")) +
  geom_density(adjust = 1.5, alpha = 0.5) +
  labs(x = "Phylogenetic Distance", y = "Phylogenetic Distance") +
  geom_boxplot(width=0.1, outlier.size = .05) + 
  ggtitle(label = "Pairwise Phylogenetic Distance")
densityNasalPerm
dev.off()

library(ggpubr)
setwd(directory.figures)
pdf(file = "violinPlotPhylogeneticClustersNasal.pdf")
violinNasalPerm = ggplot(dN, aes(x=Node, y=Dist, fill = Node)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(aes(fill = Node), width=0.1, fill = "white", outlier.alpha = .1, outlier.size = .5) + 
  scale_fill_manual(values=c("blue", "red", "grey")) +
  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  labs(x = element_blank(), y = "Phylogenetic Distance") +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("+", "-"), c("+", "Random"), c("-", "Random"))) +
  stat_compare_means(method = "kruskal.test", label.y = 4)
violinNasalPerm
dev.off()

###
# Random Forest
library(randomForest)
library(Boruta)
library(caret)

psNASV = phyloseq(otu_table(psNASV), tax_table(psNASV), sample_data(metaN))
psNCTU = phyloseq(otu_table(psNCTU), tax_table(psNCTU), sample_data(metaN))
psGASV = phyloseq(otu_table(psGASV), tax_table(psGASV), sample_data(metaG))
psGCTU = phyloseq(otu_table(psGCTU), tax_table(psGCTU), sample_data(metaG))

data.asv.n = getModelingTable(psNASV)
data.asv.n = data.asv.n[,c(taxa_names(psNASV), "Status")]
data.asv.n$Status = factor(data.asv.n$Status)
data.ctu.n = getModelingTable(psNCTU)
data.ctu.n$Status = factor(data.ctu.n$Status)
data.ctu.n= data.ctu.n[,c(taxa_names(psNCTU), "Status")]

data.asv.g = getModelingTable(psGASV)
data.asv.g = data.asv.g[,c(taxa_names(psGASV), "Status")]
data.asv.g$Status = factor(data.asv.g$Status)
data.ctu.g = getModelingTable(psGCTU)
data.ctu.g = data.ctu.g[,c(taxa_names(psGCTU), "Status")]
data.ctu.g$Status = factor(data.ctu.g$Status)

asv.n.selected = Boruta(Status ~., data = data.asv.n, ntree = 100)
asv.n.selected = names(asv.n.selected$finalDecision)[which(asv.n.selected$finalDecision %in% c("Tentative", "Confirmed"))]
ctu.n.selected = Boruta(Status~., data = data.ctu.n, ntree = 100)
ctu.n.selected = names(ctu.n.selected$finalDecision)[which(ctu.n.selected$finalDecision %in% c("Tentative", "Confirmed"))]
asv.g.selected = Boruta(Status~., data = data.asv.g, ntree = 100)
asv.g.selected = names(asv.g.selected$finalDecision)[which(asv.g.selected$finalDecision %in% c("Tentative", "Confirmed"))]
ctu.g.selected = Boruta(Status~., data = data.ctu.g, ntree = 100)
ctu.g.selected = names(ctu.g.selected$finalDecision)[which(ctu.g.selected$finalDecision %in% c("Tentative", "Confirmed"))]

# Get combined data table
microbiomeDataAlone = as_tibble(metaAll[,c("GutID", "NasalID", "FelineCalcivirusPCR")]) %>%
  filter(GutID %in% metaG$ID) %>%
  filter(NasalID %in% metaN$ID) %>%
  dplyr::filter(!is.na(FelineCalcivirusPCR)) %>%
  select(-FelineCalcivirusPCR) %>%
  left_join(getModelingTable(psNASV) %>% select(ID, all_of(asv.n.selected)), by = c("NasalID" = "ID")) %>%
  left_join(getModelingTable(psNCTU) %>% select(ID, all_of(ctu.n.selected)), by = c("NasalID" = "ID")) %>%
  left_join(getModelingTable(psGASV) %>% select(ID, all_of(asv.g.selected)), by = c("GutID" = "ID")) %>%
  left_join(getModelingTable(psGCTU) %>% select(ID, all_of(ctu.g.selected)), by = c("GutID" = "ID")) %>%
  left_join(metaN %>% select(ID, Status), by = c("NasalID" = "ID"))
microbiomeDataAlone$Status = factor(microbiomeDataAlone$Status, levels = c("1", "0"))
microbiomeDataAlone

pcrNames = c( "PCR")
microbiomePCRData = microbiomeDataAlone %>% 
  left_join(metaAll %>% select(GutID, all_of(pcrNames))) 

microbiomeDataAlone = microbiomeDataAlone %>% select(-GutID, -NasalID)
microbiomePCRData = microbiomePCRData %>% select(-GutID, -NasalID)
microbiomeDataAlone
microbiomePCRData

mb.rf = randomForest(Status~., data = microbiomeDataAlone, ntree = 1000, importance = T)
mb.pcr.rf = randomForest(Status~., data = microbiomePCRData, ntree = 1000, importance = T)
confusionMatrix(data = mb.rf$predicted, reference = mb.rf$y)
confusionMatrix(data = mb.pcr.rf$predicted, reference = mb.pcr.rf$y)

pcr.gut = metaG[,c("Status", "FCV", "Bordetella", "Chlamydophila", "FHV", "Mycoplasma", "PCR")]
pcr.gut = pcr.gut[complete.cases(pcr.gut),]
pcr.gut = pcr %>% mutate(
  pcrResult = case_when(
    PCR == "P"~ 1,
    PCR == "N"~ 0
  )
)
pcr.gut$Status = factor(pcr.gut$Status, levels = c("1", "0"))
pcr.gut$pcrResult = factor(pcr.gut$pcrResult, levels = c("1", "0"))
confusionMatrix(reference = pcr.gut$Status, pcr.gut$pcrResult)

taxGAll %>% filter(Taxon %in% c(asv.g.selected, ctu.g.selected))
taxNAll %>% filter(Taxon %in% c(asv.n.selected, ctu.n.selected))

impScores = as_tibble(mb.pcr.rf$importance, rownames = "ID") %>%
  left_join(taxGAll %>% filter(Taxon %in% c(asv.g.selected, ctu.g.selected)), by = c("ID" = "Taxon")) %>%
  left_join(taxNAll %>% filter(Taxon %in% c(asv.n.selected, ctu.n.selected)), by = c("ID" = "Taxon")) %>%
  mutate(
    Kingdom = case_when(
      is.na(Kingdom.y) ~ Kingdom.x,
      is.na(Kingdom.x) ~ Kingdom.y
    ),
    Phylum = case_when(
      is.na(Phylum.y) ~ Phylum.x,
      is.na(Phylum.x) ~ Phylum.y
    ),
    Class = case_when(
      is.na(Class.x) ~ Class.y,
      is.na(Class.y) ~Class.x
    ),
    Order = case_when(
      is.na(Order.x)~ Order.y,
      is.na(Order.y) ~ Order.x
    ),
    Family = case_when(
      is.na(Family.x) ~Family.y,
      is.na(Family.y) ~ Family.x
    ),
    Genus = case_when(
      is.na(Genus.x) ~ Genus.y,
      is.na(Genus.y) ~ Genus.x
    ),
    Species = case_when(
      is.na(Species.x) ~ Species.y,
      is.na(Species.y) ~Species.x
    )
    
  ) %>% select(-Kingdom.x, -Kingdom.y, -Phylum.x, -Phylum.y, -Class.x, -Class.y, -Order.x, -Order.y, -Family.x, -Family.y, -Genus.x, -Genus.y, -Species.x, -Species.y) %>%
  arrange (desc(MeanDecreaseAccuracy)) %>%
  mutate(
    label = case_when(
      is.na(Phylum) ~Kingdom,
      is.na(Class) ~Phylum,
      is.na(Order) ~ Class,
      is.na(Family) ~ Order,
      is.na(Genus) ~ Family,
      is.na(Species) ~ Genus,
      !is.na(Species) & !is.na(Genus) ~ Genus
    )
  )
impScores = as.data.frame(impScores)
rownames(impScores) = impScores$ID
impScores$System = NA  
impScores$Feature = NA
impScores[c(asv.n.selected, ctu.n.selected), "System"] = "Nasal"
impScores[c(asv.n.selected, asv.g.selected), "Feature"] = "ASV"
impScores[c(ctu.n.selected, ctu.g.selected), "Feature"] = "CTU"
impScores[c(asv.g.selected, ctu.g.selected), "System"] = "Gut"
impScores[which(impScores$ID == "PCR"),"label"] = "PCR"
impScores[which(impScores$label == "Family_XIII"),"label"] =  c("Family_XIII A", "Family_XIII B")
impScores[which(impScores$label == "Bacteria"),"label"] =  c("Bacteria A", "Bacteria B", "Bacteria C")
impScores[which(impScores$label == "Bacteroides"),"label"] = c("Bacteroides A", "Bacteroides B")

impScores[which(impScores$label == "Porphyromonas"),"label"] =  c("Porphyromonas A", "Porphyromonas B")
impScores[which(impScores$ID == "PCR"),"System"] = "PCR"
impScores

impScores = as_tibble(impScores) %>% arrange(desc(MeanDecreaseGini))
impScores$label = factor(impScores$label, levels = impScores$label)

setwd(directory.figures)
pdf(file = "ImporatanceScores.pdf")
ggplot(data = impScores, aes(x = reorder(label, -MeanDecreaseGini), y = MeanDecreaseGini)) + 
  geom_bar(stat = "identity", aes(fill = System)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(3, 1,2)]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 12)) +
  xlab("Taxonomic Label") + 
  ylab("Importance Score\n(Mean Decrease Accuracy)") +
  geom_text(aes(label=Feature), position=position_dodge(width=0.9), vjust=-0.25)
dev.off()


setwd(directory.figures)
save.image(file = "MondayAug16ThesisChapterEnv.RData", compress = TRUE)

