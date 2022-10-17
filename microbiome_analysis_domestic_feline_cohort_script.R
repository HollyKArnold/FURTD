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
library(ggtree)

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

# Remove animals based on exclusion criteria
meta = meta %>% filter(Animal_ID != "7") 


meta = meta %>% mutate(Status = case_when(Status == "Case" ~ "FURTD",
                                   Status == "Control" ~ "Control"))
meta$Status = factor(meta$Status)

# Create identifiers for gut and nasal sample accessions. 
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
taxGCTU = cladifierGetRefTreeTaxonomy(tree = treeG, ref = taxG)
taxGCTU = taxGCTU[rownames(ctuG),]

ctuN = t(utils::read.table("ctuNasal.txt", sep = "\t", header = T))
head(ctuN)
treeN = ape::read.tree("treeNasal.tre")
treeN
taxNCTU = cladifierGetRefTreeTaxonomy(tree = treeN, ref = taxN)
taxNCTU = taxNCTU[rownames(ctuN),]

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
# Remove individual based on exclusion criteria
bloodwork = bloodwork %>% select(!`7`)



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
phyloseq::sample_data(asvN_ps)$Household = 
  factor(phyloseq::sample_data(asvN_ps)$Household)
asvN_ps %>% sample_data %>% count(Status)

ctuG_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuG, 
                                                 taxa_are_rows = TRUE), 
                             phyloseq::tax_table(object = as.matrix(taxGCTU)),
                             phyloseq::sample_data(metaG))
phyloseq::sample_data(ctuG_ps)$Status = 
  factor(phyloseq::sample_data(ctuG_ps)$Status, levels = c("Control", "Case"))

ctuN_ps = phyloseq::phyloseq(phyloseq::otu_table(object = ctuN, 
                                                 taxa_are_rows = TRUE), 
                             phyloseq::tax_table(object = as.matrix(taxNCTU)),
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
                                      formula = asvG_ps_clr ~ Status, 
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
rdaG

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
    filter(Phylum == top_taxa_gut[i], Status == "FURTD") %>% 
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
    filter(Phylum == top_taxa_nasal[i], Status == "FURTD") %>% 
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
 

###############################################################################
#  Inflammatory Markers                                                       #
###############################################################################

# Test if difference between groups
wilcox.test(x = meta %>% 
              filter(Status == "FURTD") %>% 
              select(Albumin) %>% 
              pull(.), 
            y = meta %>% 
              filter(Status == "Control") %>% 
              select(Albumin) %>% 
              pull(.),
            exact = FALSE)

wilcox.test(x = meta %>% 
              filter(Status == "FURTD") %>% 
              select(Neutrophils) %>% pull(.), 
            y = meta %>% 
              filter(Status == "Control") %>% 
              select(Neutrophils) %>% 
              pull(.),
            exact = FALSE)


# Plot inflammatory markers
neutrophils = ggplot(meta, aes(x = Status, y = Neutrophils, color = Status)) +
  geom_boxplot(outlier.color = NA, 
               outlier.shape = 16, 
               outlier.size = 2, 
               notch = FALSE) +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 12500, linetype = "dashed", color = "grey")+
  scale_color_manual(values = c(pal[1], pal[3])) + 
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Neutrophils", x = "", y = "cells/uL") +
  theme(axis.text = element_text(size = 18), 
        plot.title = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c("Control", "FURTD")), 
                             paired = FALSE, 
                             label = "p.signif",
                             label.y = 12500,
                             #label.x.npc = "top", 
                             symnum.args = list(
                               cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                               symbols = c("****", "***", "**", "*", "ns")),
                             vjust = 0.6, 
                             size = 10)  + 
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
neutrophils
albumin = ggplot(meta, aes(x = Status, y = Albumin, color = Status)) +
  geom_boxplot(outlier.color = NA, 
               outlier.shape = 16, 
               outlier.size = 2, 
               notch = FALSE) +
  geom_hline(yintercept = 2.6, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c(pal[1], pal[3])) + 
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Albumin", x = "", y = "g/dL") +
  theme(axis.text = element_text(size = 18), 
        plot.title = element_text(size = 20), 
        axis.title = element_text(size = 20)) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             comparisons = list(c("Control", "FURTD")), 
                             paired = FALSE, 
                             label.x.npc = "top", 
                             symnum.args = list(
                               cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                               symbols = c("****", "***", "**", "*", "ns")),
                             vjust = 0.6,
                             size = 10) + 
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = .6, size = 5)
albumin


inflammatory_markers = cowplot::plot_grid(neutrophils, 
                                          albumin, 
                                          labels = c("A", "B"))
setwd(directory.figures)
ggsave(filename = "inflammatory_markers.pdf", 
       inflammatory_markers, 
       height = 5, 
       width = 8)

###############################################################################
#  Clades associating with clinical signs                                     #
###############################################################################

# Factor Status
metaG = 
  metaG %>% 
  mutate(Status = case_when(Status == "FURTD" ~ 1, 
                            Status == "Control" ~ 0)
)
metaG = metaG %>% mutate("ID" = GutID)

metaN = 
  metaN %>% 
  mutate(Status = case_when(Status == "FURTD" ~ 1, 
                            Status == "Control" ~ 0)
  )
metaN = metaN %>% mutate("ID" = NasalID)

# Set filtering thresholds
qThresh = 0.05 # fdr threshold
prevThresGut = 6 # prevalence threshold 40% of individuals sampled
prevThresNasal = 7 # prevalence threshold 40% of individuals sampled

# Calculate significant gut ASVs
prevGASV = apply(X = otu_table(asvG_ps), 
                 MARGIN = 1, 
                 FUN = function(x){sum(x > 0)})
psGASV.nb = phyloseq(otu_table(asvG_ps), 
                     tax_table(asvG_ps), 
                     sample_data(metaG), 
                     phy_tree(asvG_ps))
psGASV.nb = prune_taxa(x = psGASV.nb, 
                       taxa = names(prevGASV)[which(prevGASV > prevThresGut)]) 
nb.g.asv = getNBGLMS(covariates = c("Status"), phyloseq = psGASV.nb)
dim(nb.g.asv[which(nb.g.asv$q < qThresh),])
sigGutASV = 
  nb.g.asv[which(nb.g.asv$q < qThresh & 
                   nb.g.asv$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigGutASV

# Calculate significant gut CTUs
prevGCTU = apply(X = otu_table(ctuG_ps), 
                 MARGIN = 1, 
                 FUN = function(x){sum(x > 0)})
psGCTU.nb = phyloseq(otu_table(ctuG_ps), 
                     tax_table(ctuG_ps), 
                     sample_data(metaG))
psGCTU.nb = prune_taxa(x = psGCTU.nb, 
                       taxa = names(prevGCTU)[which(prevGCTU > prevThresGut)]) 

# Prune CTUs that are too close to root to model with NB. Nodes close
# to root are expected to have close to uniform distribution 
toPrune = vector() 
for(i in 1:length(taxa_names(psGCTU.nb))){
  root = max(otu_table(psGCTU.nb)["root",])
  if(max(otu_table(psGCTU.nb)[i,]) == root){
    toPrune = c(taxa_names(psGCTU.nb)[i], toPrune)
  }
}
toPrune # removes "root", "node2" "node3" and "node1" 
psGCTU.nb = 
  prune_taxa(x = psGCTU.nb, 
             taxa = taxa_names(psGCTU.nb)[!taxa_names(psGCTU.nb) %in% toPrune]) 
nb.g.ctu = getNBGLMS(covariates = c("Status"), phyloseq = psGCTU.nb)
dim(nb.g.ctu[which(nb.g.ctu$q < qThresh & 
                     nb.g.ctu$model %in%  c("glm.nb", "glm.nb.plus.01")),])
sigGutCTU = 
  nb.g.ctu[which(nb.g.ctu$q < qThresh & 
                   nb.g.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigGutCTU

# Calculate significant ASVs nasal microbiome
prevNASV = apply(X = otu_table(asvN_ps), 
                 MARGIN = 1, 
                 FUN = function(x){sum(x > 0)})
psNASV.nb = phyloseq(otu_table(asvN_ps), 
                     tax_table(asvN_ps), 
                     phy_tree(asvN_ps), 
                     sample_data(metaN))
psNASV.nb = 
  phyloseq::prune_taxa(x = psNASV.nb, 
                       taxa = names(prevNASV)[which(prevNASV > prevThresNasal)])
nb.n.asv = getNBGLMS(covariates = c("Status"), phyloseq = psNASV.nb)
dim(nb.n.asv[which(nb.n.asv$q < qThresh & 
                     nb.n.asv$model %in%  c("glm.nb", "glm.nb.plus.01")),])
sigNasalASV = 
  nb.n.asv[which(nb.n.asv$q < qThresh & 
                   nb.n.asv$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigNasalASV

# Calculate significant CTUs in nasal microbiome
prevNCTU = apply(X = otu_table(ctuN_ps), 
                 MARGIN = 1, 
                 FUN = function(x){sum(x > 0)})
psNCTU.nb = phyloseq(otu_table(ctuN_ps), 
                     tax_table(ctuN_ps), 
                     sample_data(metaN))
psNCTU.nb = 
  prune_taxa(x = psNCTU.nb, 
             taxa = names(prevNCTU)[which(prevNCTU > prevThresNasal)]) 

# Prune CTUs that are too close to root to model with NB. Nodes close
# to root are expected to have close to uniform distribution 
toPrune = vector() 
for(i in 1:length(taxa_names(psNCTU.nb))){
  root = max(otu_table(psNCTU.nb)["root",])
  if(max(otu_table(psNCTU.nb)[i,]) == root){
    toPrune = c(taxa_names(psNCTU.nb)[i], toPrune)
  }
}
toPrune # Removes "root", "node4", "node5", and "node3". 
psNCTU.nb = 
  prune_taxa(x = psNCTU.nb, 
             taxa = taxa_names(psNCTU.nb)[!taxa_names(psNCTU.nb) %in% toPrune]) 
nb.n.ctu = getNBGLMS(covariates = c("Status"), phyloseq = psNCTU.nb)
dim(nb.n.ctu[which(nb.n.ctu$q < qThresh & 
                     nb.n.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),])
sigNasalCTU = 
  nb.n.ctu[which(nb.n.ctu$q < qThresh & 
                   nb.n.ctu$model %in% c("glm.nb", "glm.nb.plus.01")),]
sigNasalCTU

# Combine significant features into one frame
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

sig %>% filter(microbiome == "Gut") %>% nrow() #136
sig %>% filter(microbiome == "Nasal")%>% nrow() #89

###############################################################################
#  Visualize significant microbial features gut                               #
###############################################################################


# Get significant clade roots gut
treeGLevel = getRefTreeCladeLevel(tree = treeG)
rownames(treeGLevel) = treeGLevel$node
rootsG = RemoveNestedClades(
  cladeList = treeGLevel[unique(sigGutCTU %>% 
                                  filter(covariate %in% c("Status"))) %>% 
                           pull(clade) %>% unique(),], 
  tree = treeG)
rootsG
taxGAll = tidyr::as_tibble(rbind(
  phyloseqCompanion::taxa.data.table(ps = asvG_ps), 
  phyloseqCompanion::taxa.data.table(ps = (ctuG_ps))))
taxGAll %>% filter(Taxon %in% rootsG) %>% print(n = 35)
as_tibble(sigGut) %>% 
  filter(covariate == "Status") %>% 
  left_join(taxGAll, by = c("clade" = "Taxon")) %>% 
  arrange(Family) %>% 
  print(n = 107)

dG = treeDataFrame("Status", tree = treeG)
head(dG)
dG$nodeName = rownames(dG)
dG = dG %>% left_join(taxGAll, by = c("nodeName" = "Taxon"))
rownames(dG) = dG$nodeName
head(dG)
dG[sigGut[which(sigGut$covariate == "Status"), "clade"], "Status"] = 
  sigGut[which(sigGut$covariate == "Status"), "Estimate"]
dG = dG %>% 
  mutate(
    Status = case_when(
      Status < 0 ~ "-",
      Status >= 0 ~ "+"
    ))

dG
dG %>% filter(nodeName %in% rootsG) 

idxsG = c(treeG$tip.label, treeG$node.label)


fsize = 3
bsize = 1.5
offset.bar = 1

# Tree with labels for reference
cxGutTreeLabel = ggtree::ggtree(treeG, layout = "circular", 
                                size = 0.1, branch.length = "none") %<+% dG +
  geom_point2(aes(color = Phylum), size = 1, alpha = .3, shape = 19) +
  geom_tiplab(size = 0.8) +
  geom_nodelab(size = 0.8) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
cxGutTreeLabel

setwd(directory.figures)
ggsave(filename = "labeled_tree_gut.pdf",
       cxGutTreeLabel, 
       height = 20, 
       width = 20)


# Color palette for labels
palPhylumLab = RColorBrewer::brewer.pal(n = 6, "Dark2")
palPhylumLab = palPhylumLab[c(2, 4, 5, 6)]
names(palPhylumLab) = c("orange", "pink", "green", "yellow")

# Create a column to use for color coding figure
dG %>% filter(nodeName %in% rootsG) %>% pull(Phylum) %>% unique()

dG = dG %>%
  mutate(phylumBarKeyColor = 
           case_when(nodeName %in% rootsG & Phylum == "Fusobacteria" ~ 
                       paste(c("'", palPhylumLab[4], "'"), 
                             sep = "", 
                             collapse = ""),
                                       
                     nodeName %in% rootsG & Phylum == "Bacteroidetes" 
                     ~ paste(c("'", palPhylumLab[1], "'"), 
                             sep = "", 
                             collapse = ""),
                                       
                     nodeName %in% rootsG & Phylum == "Proteobacteria" ~ 
                       paste(c("'", palPhylumLab[3], "'"), 
                             sep = "", 
                             collapse = ""),
                                       
                     nodeName %in% rootsG & Phylum == "Firmicutes" ~ 
                       paste(c("'", palPhylumLab[2], "'"), 
                             sep = "", 
                             collapse = ""),
                                       
                     nodeName %in% rootsG & is.na(Phylum) ~ "'grey'"))

# Get labels 
paste(c(
  cladeLabels(dataFrame = dG %>% 
                filter(nodeName %in% rootsG) %>% 
                filter(Status == "+"), 
              idxs = "idxsG", 
              nodeName = "nodeName",  
              fontsize = "10", 
              bsize = 5, 
              offset.bar = 0.8),
  cladeLabels(dataFrame = dG %>% filter(!is.na(Status)) %>% 
                filter(nodeName %in% rootsG) %>%  
                filter(Status == "-"), 
              idxs = "idxsG",  
              nodeName = "nodeName", 
              fontsize = "10", 
              bsize = 5, 
              offset.bar = 0.8)
  
), sep = "", collapse = "")

# For reference tree with labeled clades
cxGutTree = ggtree::ggtree(treeG, 
                           layout = "circular", 
                           size = 0.3, 
                           branch.length = "none") %<+% dG
cxGutTree = cxGutTree + 
  geom_point2(aes(subset = (!is.na(Status)), 
                  color = Status, 
                  fill = Status, shape = Status), alpha = 0.8, size = 14) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = c(pal[1], pal[3])) +
  theme(legend.position = "left", legend.text = element_text(size = 14)) + 
  geom_cladelabel(node = which(idxsG=='node5'), label = 'Bacteria', align = T,
                  geom = 'text', angle = 0, color = 'grey', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node84'), label = 'Fusobacterium', 
                  align = T, geom = 'text', angle = 0, color = '#E6AB02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node186'), label = 'Bacteroidales',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node241'), label = 'Marinifilaceae',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node264'), label = 'Bacteroides', 
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsG=='node276'), label = 'Bacteroides',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node284'), label = 'Bacteroides', 
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node293'), label = 'Parabacteroides',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node336'), label = 'Alistipes', 
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node434'), 
                  label = 'Succinivibrionaceae', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node462'), 
                  label = 'Enterobacteriaceae', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node508'), 
                  label = 'Betaproteobacteriales', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node536'), 
                  label = 'Clostridium_sensu_stricto_1', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', 
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node556'), 
                  label = 'Peptostreptococcus', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', 
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node562'), 
                  label = 'Peptostreptococcaceae', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', 
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node600'), label = 'Peptoniphilus',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node617'), label = 'Finegoldia', 
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsG=='node625'), label = 'Bacteria', 
                  align = T, geom = 'text', angle = 0, color = 'grey',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node678'), label = 'Anaerostipes',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node684'), label = 'Lachnospiraceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node791'), label = 'Lachnospiraceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node865'), label = 'Ruminococcaceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node873'), label = 'Subdoligranulum',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node889'), label = 'Ruminococcaceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node908'), label = 'Megamonas', 
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node927'), label = 'Streptococcus',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node934'), label = 'Lactobacillaceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node215'), label = 'Prevotella_9',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node281'), label = 'Bacteroides', 
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node294'), label = 'Porphyromonas',
                  align = T, geom = 'text', angle = 0, color = '#D95F02',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node366'), 
                  label = 'Desulfovibrionales', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node387'), label = 'Bacteria', 
                  align = T, geom = 'text', angle = 0, color = 'grey',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node444'), 
                  label = 'Gammaproteobacteria', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node650'), label = 'Peptococcus', 
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node708'), label = 'Clostridiales',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node782'), label = 'Lachnospiraceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node804'), 
                  label = 'Ruminococcaceae_UCG-014', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node830'), label = 'Ruminococcaceae',
                  align = T, geom = 'text', angle = 0, color = '#E7298A',
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node836'), 
                  label = 'Ruminococcaceae_UCG-004', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', 
                  fontsize = 10, barsize = 5, offset = 0.8)
cxGutTree

# Strip legend
status_legend = ggpubr::get_legend(cxGutTree)
ggpubr::as_ggplot(status_legend)
setwd(directory.figures)
ggsave(filename = "labeled_Gut_Tree.pdf", cxGutTree, height = 40, width = 40)
ggsave(filename = "legend_status_gut.png", ggpubr::as_ggplot(status_legend))

# Get phylum legend
dG %>% filter(nodeName %in% rootsG) 
legend_phylum_gut = ggplot(data = dG %>% filter(nodeName %in% rootsG), ) + 
  geom_point(aes(color = phylumBarKeyColor, x = node, y = node)) +
  scale_color_manual(values = c(palPhylumLab, "grey"), 
                     labels = c("Bacteroidetes", 
                                "Firmicutes", 
                                "Proteobacteria", 
                                "Fusobacteria", 
                                ">Phylum")) +
  labs(color = "Phylum") +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
legend_phylum_gut = ggpubr::get_legend(legend_phylum_gut)
ggpubr::as_ggplot(legend_phylum_gut)
setwd(directory.figures)
ggsave(filename = "legend_phylum_gut.png", 
       ggpubr::as_ggplot(legend_phylum_gut))

  
cxGutTree = ggtree::ggtree(treeG, layout = "circular", size = 0.3, 
                           branch.length = "none") %<+% dG
cxGutTree = cxGutTree + 
  geom_point2(aes(subset = (!is.na(Status)), color = Status, fill = Status, 
                  shape = Status), alpha = 0.8, size = 14) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  scale_color_manual(values = c(pal[1], pal[3])) +
  theme(legend.position = "none") + 
  geom_cladelabel(node = which(idxsG=='node5'), label = '', align = T, 
                  geom = 'text', angle = 0, color = 'grey', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node84'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E6AB02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node186'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node241'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node264'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsG=='node276'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node284'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node293'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node336'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node434'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#66A61E', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node462'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#66A61E', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node508'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#66A61E', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node536'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node556'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node562'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node600'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node617'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsG=='node625'), label = '', align = T, 
                  geom = 'text', angle = 0, color = 'grey', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node678'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node684'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node791'), label = '', align = T,
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node865'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node873'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node889'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node908'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node927'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node934'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node215'), label = '', align = T,
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node281'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node294'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#D95F02', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node366'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#66A61E', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node387'), label = '', align = T, 
                  geom = 'text', angle = 0, color = 'grey', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node444'), label = '', align = T,
                  geom = 'text', angle = 0, color = '#66A61E', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node650'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node708'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node782'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node804'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node830'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10, 
                  barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsG=='node836'), label = '', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', fontsize = 10,
                  barsize = 5, offset = 0.8)
cxGutTree
setwd(directory.figures)
ggsave(filename = "full_gut_tree_labeless.png", 
       cxGutTree, 
       height = 40, 
       width = 40)


###############################################################################
#  Get reduced tree for visualization gut                                     #
###############################################################################

# Get vector of tips in all significant clades
tipsG = vector()
for(i in 1:length(rootsG)){
  curNode = idxsG[which(idxsG == rootsG[i])]
  curClade = ape::extract.clade(treeG, curNode)
  tipsG = c(tipsG, curClade$tip.label)
  tipsG = unique(tipsG)
}
tipsG = unique(c(tipsG, sigGut %>% 
                   filter(type == "ASV") %>% 
                   filter(covariate %in% c("Status")) %>% 
                   pull(clade)))

# Drop tips which are not descendants of significant ancestors
treeGSub = ape::drop.tip(phy = treeG, tip = setdiff(treeG$tip.label, tipsG))  
dGSub = dG[c(treeGSub$tip.label, treeGSub$node.label),]
dGSub$node = seq(from = 1, to = length(treeGSub$tip.label) + 
                   length(treeGSub$node.label), by = 1)
idxsG = seq(from = 1, to = length(treeGSub$tip.label) + 
              length(treeGSub$node.label), by = 1)
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
    
  ))

# Move labels outward from smaller clades
nodedf$size = NA
for(i in 1:nrow(nodedf)){
  curClade = ape::extract.clade(treeGSub, node = nodedf[i, "node"])
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

lac = ape::extract.clade(phy = treeGSub, node = idxsG["node934"] )
as_tibble(taxGAll) %>% filter(Taxon %in% lac$tip.label) %>% print(n = 100)


cxGutTree = ggtree(treeGSub, layout = "fan", size = 0.2, open.angle = 5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=2.2, alpha=0.2, 
               fill="grey", color="grey50",size=0.05) + 
  geom_hilight(data=nodedf, mapping=aes(node=node),extendto=3.3, alpha=0.3, 
               fill="white", color="grey50",size=0.05) + 
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

# Output a labeled version for checking tips / labels match graph
setwd(directory.figures)
pdf(file = "tree_reduced_gut_tip_node_labels.pdf")
cxGutTree + geom_tiplab(size = 0.5) + geom_nodelab(size = 0.5)
dev.off()

colnames(dGSub)[which(colnames(dGSub) == "Status" )] = "Association"
cxGutTree = cxGutTree %<+% dGSub + 
  geom_point2(aes(subset = (!is.na(Association)), color = Association), 
              alpha = 0.6, size = 2, shape = 16) +
  scale_color_manual(values = c(pal[1], pal[3])) +
  guides(starshape = guide_legend(override.aes = list(size = 8))) 
cxGutTree


tmpPhyloseq = phyloseq(otu_table(asvG_ps), 
                       phy_tree(asvG_ps), 
                       tax_table(asvG_ps), 
                       sample_data(asvG_ps))
tmpPhyloseq = prune_taxa(taxa = treeGSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = 
  prune_samples(samples = 
                  rownames(sample_data(tmpPhyloseq)
                           [which(sample_data(tmpPhyloseq)$Status == 
                                    "FURTD"),]), x = tmpPhyloseq)
furtdAve = apply(otu_table(tmpPhyloseq)[treeGSub$tip.label,], 1, mean)


furtdAve = data.frame("Average" = log(furtdAve + 0.001),
                      "Status" = rep("FURTD", times = length(furtdAve)))
furtdAve$ID = rownames(furtdAve)

tmpPhyloseq = phyloseq(otu_table(asvG_ps), 
                       phy_tree(asvG_ps), 
                       tax_table(asvG_ps), 
                       sample_data(asvG_ps))
tmpPhyloseq = prune_taxa(taxa = treeGSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = 
  prune_samples(samples = 
                              rownames(sample_data(tmpPhyloseq)[
                                which(sample_data(tmpPhyloseq)$Status == 
                                        "Control"),]), x = tmpPhyloseq)
controlAve = apply(otu_table(tmpPhyloseq)[treeGSub$tip.label,], 1, mean)


controlAve = data.frame("Average" = log(controlAve + 0.001),
                        "Status" = rep("Control", times = length(controlAve)))
controlAve$ID = rownames(controlAve)
datAve = rbind(furtdAve, controlAve)


cxGutTree = cxGutTree + 
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(data=datAve, geom=geom_tile,
             mapping=aes(y=ID, x=Status, fill=Status, alpha = Average),
             color = "grey50", offset = 0.45,size = 0.01) +
  scale_fill_manual(values = c(pal[1], pal[3])) 
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

cxGutTree = cxGutTree + 
  ggtreeExtra::geom_fruit(data=dat3, 
                          geom=geom_bar, 
                          mapping=aes(y=ID, x=HigherAve, fill = HigherClass), 
                          pwidth=0.38, 
                          offset = .1,
                          orientation="y",  
                          stat="identity") 
cxGutTree = cxGutTree + 
  guides(color = guide_legend(override.aes = list(size=8))) +
  theme(legend.text = element_text(size = 9), 
        legend.title = element_text(size = 12), legend.position = "bottom") +
  guides(alpha=guide_legend(title="Log Average\nAbundance"))
cxGutTree
cxGutTreeLegend = ggpubr::get_legend(cxGutTree)

# Plot reduced tree
setwd(directory.figures)
ggsave(cxGutTree, file = "cladesClinicalSignsReducedGut.pdf", 
       height = 25, width = 23, units = "cm")
ggsave(cxGutTree, file = "cladesClinicalSignsReducedGut.png", 
       height = 25, width = 23, units = "cm")
lay = rbind(c(1, 1),
            c(1, 1),
            c(2, 2),
            c(2, 2))

pp = gridExtra::arrangeGrob(cxGutTreeLegend$grobs[[2]], 
                            cxGutTreeLegend$grobs[[3]],  
                            layout_matrix = lay)
plot(pp)
ggsave(pp, 
       file = "cladesClinicalSignsReducedGutKey.pdf", 
       height = 9, 
       width = 5)


###############################################################################
#  Visualize significant microbial features nasal                             #
###############################################################################

treeNLevel = getRefTreeCladeLevel(tree = treeN)
rownames(treeNLevel) = treeNLevel$node
rootsN = 
  RemoveNestedClades(cladeList = 
                       treeNLevel[unique(sigNasalCTU %>% 
                                           filter(covariate %in% c("Status")))
                                  %>% pull(clade) %>% unique(),], 
                     tree = treeN)
rootsN
taxNAll = as_tibble(rbind(phyloseqCompanion::taxa.data.table(ps = asvN_ps), 
                          phyloseqCompanion::taxa.data.table(ps = (ctuN_ps))))
taxNAll %>% filter(Taxon %in% rootsN) %>% print(n = 35)
as_tibble(sigNasal) %>% 
  filter(covariate == "Status") %>% 
  left_join(taxNAll, by = c("clade" = "Taxon")) %>% 
  arrange(Family) %>% print(n = 107)

dN = treeDataFrame("Status", tree = treeN)
head(dN)
dN$nodeName = rownames(dN)
dN = dN %>% left_join(taxNAll, by = c("nodeName" = "Taxon"))
rownames(dN) = dN$nodeName
head(dN)
dN[sigNasal[which(sigNasal$covariate == "Status"), "clade"], "Status"] = 
  sigNasal[which(sigNasal$covariate == "Status"), "Estimate"]
dN = dN %>% 
  mutate(
    Status = case_when(
      Status < 0 ~ "-",
      Status >= 0 ~ "+"
    ))

dN
dN %>% filter(nodeName %in% rootsN) 
idxsN = c(treeN$tip.label, treeN$node.label)


# Color palette for labels
palPhylumLab = c(RColorBrewer::brewer.pal(n = 8, "Dark2"), 
                 RColorBrewer::brewer.pal(n = 11, "Spectral"))
palPhylumLab = palPhylumLab[c(2, 4, 5, 6, 8, 9, 18, 19)]
scales::show_col(palPhylumLab)
palPhylumLab
names(palPhylumLab) = 
  c("orange", "pink", "green", "yellow", "grey", "red", "blue", "purple")

# Create a column to use for color coding figure
dN %>% filter(nodeName %in% rootsN) %>% pull(Phylum) %>% unique()

dN = dN %>%
  mutate(phylumBarKeyColor = 
           case_when(nodeName %in% rootsN & Phylum == "Fusobacteria" ~ 
                       paste(c("'", 
                               palPhylumLab[which(names(palPhylumLab) == 
                                                    "yellow")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Bacteroidetes" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "orange")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Proteobacteria" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "green")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Firmicutes" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "pink")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Spirochaetes" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "red")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Patescibacteria" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "blue")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     nodeName %in% rootsN & Phylum == "Actinobacteria" ~ 
                       paste(c("'", palPhylumLab[which(names(palPhylumLab) ==
                                                         "purple")], "'"), 
                             sep = "", 
                             collapse = ""),
                     
                     
                     nodeName %in% rootsN & is.na(Phylum) ~ "'grey'")) 
  


paste(c(
  cladeLabels(dataFrame = dN %>% 
              filter(nodeName %in% rootsN) %>% 
              filter(Status == "+"), 
              idxs = "idxsN", 
              nodeName = "nodeName",  
              fontsize = "10", 
              bsize = 5, 
              offset.bar = 0.8),
  cladeLabels(dataFrame = dN %>% filter(!is.na(Status)) %>% 
                filter(nodeName %in% rootsN) %>%  
                filter(Status == "-"), 
              idxs = "idxsN",  
              nodeName = "nodeName", 
              fontsize = "10", 
              bsize = 5, 
              offset.bar = 0.8)
  
), sep = "", collapse = "")




cxNasalTree = ggtree(treeN, 
                     layout = "circular", 
                     branch.length = 'none', 
                     size = 0.15) %<+% dN
cxNasalTree = cxNasalTree +  
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 18), legend.position = "right") +
  geom_point2(aes(subset = (!is.na(Status)), color = Status, fill = Status, shape = Status), alpha = 0.8, size = 14) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  scale_color_manual(values = c(pal[1], pal[3])) + 
  geom_cladelabel(node = which(idxsN=='node31'), 
                  label = 'Actinomyces', align = T, geom = 'text', 
                  angle = 0, color = '#5E4FA2', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node100'), 
                  label = 'Defluviitaleaceae_UCG-011', align = T, 
                  geom = 'text', angle = 0, color = '#E7298A', 
                  fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node150'), 
                  label = 'Catonella', align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node168'), 
                  label = 'Ezakiella', align = T, geom = 'text', angle = 0,
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node217'), 
                  label = 'Proteocatella', align = T, geom = 'text', 
                  angle = 0, color = '#E7298A', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node281'), 
                  label = 'Streptococcus', align = T, geom = 'text',
                  angle = 0, color = '#E7298A', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node385'), 
                  label = 'Leptotrichiaceae', align = T, geom = 'text', 
                  angle = 0, color = '#E6AB02', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node396'), 
                  label = 'Fusobacteriaceae', align = T, geom = 'text', 
                  angle = 0, color = '#E6AB02', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node414'), 
                  label = 'Sediminispirochaeta', align = T, geom = 'text', 
                  angle = 0, color = '#9E0142', fontsize = 10, barsize = 5,
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node446'), 
                  label = 'Treponema_2', align = T, geom = 'text', 
                  angle = 0, color = '#9E0142', fontsize = 10, 
                  barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node452'), 
                  label = 'Treponema_2', align = T, geom = 'text', 
                  angle = 0, color = '#9E0142', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node489'), 
                  label = 'Porphyromonas', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5,
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node572'), 
                  label = 'Capnocytophaga', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node664'), 
                  label = 'Alloprevotella', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node757'), 
                  label = 'Bacteria', align = T, geom = 'text', angle = 0, 
                  color = 'grey', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node806'), 
                  label = 'Burkholderiaceae', align = T, geom = 'text', 
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node864'), 
                  label = 'Gammaproteobacteria', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node937'), 
                  label = 'Pasteurellaceae', align = T, geom = 'text', 
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node126'), 
                  label = 'Johnsonella', align = T, geom = 'text', angle = 0,
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node200'), 
                  label = 'Fusibacter', align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node218'), 
                  label = 'Proteocatella', align = T, geom = 'text', 
                  angle = 0, color = '#E7298A', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node231'), 
                  label = 'Clostridiales', align = T, geom = 'text', 
                  angle = 0, color = '#E7298A', fontsize = 10, barsize = 5, 
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node293'),
                  label = 'MVP-15', align = T, geom = 'text', angle = 0,
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node373'), 
                  label = 'JGI_0000069-P22', align = T, geom = 'text', 
                  angle = 0, color = '#3288BD', fontsize = 10, barsize = 5,
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node428'), 
                  label = 'Treponema_2', align = T, geom = 'text', 
                  angle = 0, color = '#9E0142', fontsize = 10, barsize = 5, 
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node553'), 
                  label = 'Bergeyella', align = T, geom = 'text',
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node575'), 
                  label = 'Flavobacteriaceae', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node614'), 
                  label = 'Porphyromonas', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, 
                  barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node634'), 
                  label = 'Porphyromonas', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5, 
                  offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node672'), 
                  label = 'Bacteroides', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node677'), 
                  label = 'Alloprevotella', align = T, geom = 'text', 
                  angle = 0, color = '#D95F02', fontsize = 10, barsize = 5, 
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node730'), 
                  label = 'Absconditabacteriales_(SR1)', align = T, 
                  geom = 'text', angle = 0, color = '#3288BD', 
                  fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node884'), 
                  label = 'Moraxella', align = T, geom = 'text', 
                  angle = 0, color = '#66A61E', fontsize = 10, barsize = 5,
                  offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node931'), 
                  label = 'Bibersteinia', align = T, geom = 'text',
                  angle = 0, color = '#66A61E', fontsize = 10, 
                  barsize = 5, offset = 0.8)
cxNasalTree

# Strip legend
status_legend = ggpubr::get_legend(cxNasalTree)
cxNasalTree = cxNasalTree + theme(legend.position = "none")
ggpubr::as_ggplot(status_legend)

setwd(directory.figures)
ggsave(filename = "labeled_Nasal_Tree.pdf", 
       cxNasalTree, 
       height = 40, 
       width = 40)
ggsave(filename = "legend_status_nasal.png", ggpubr::as_ggplot(status_legend))


# Get phylum legend
dN %>% filter(nodeName %in% rootsN) 
legend_phylum_nasal = ggplot(data = dN %>% filter(nodeName %in% rootsN)) + 
  geom_point(aes(color = phylumBarKeyColor, x = node, y = node)) +
  scale_color_manual(values = 
                       c(palPhylumLab[which(names(palPhylumLab) == "orange")],
                         palPhylumLab[which(names(palPhylumLab) == "pink")],
                         palPhylumLab[which(names(palPhylumLab) == "green")],
                         palPhylumLab[which(names(palPhylumLab) == "yellow")],
                         palPhylumLab[which(names(palPhylumLab) == "blue")],
                         palPhylumLab[which(names(palPhylumLab) == "red")],
                         palPhylumLab[which(names(palPhylumLab) == "purple")],
                         "grey"),
                     labels = c("Bacteroidetes", 
                                "Firmicutes", 
                                "Proteobacteria", 
                                "Fusobacteria", 
                                "Patescibacteria",
                                "Spirochaetes",
                                "Actinobacteria",
                                ">Phylum")) +
  labs(color = "Phylum") +
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 5))) 
legend_phylum_nasal
legend_phylum_nasal = ggpubr::get_legend(legend_phylum_nasal)
ggpubr::as_ggplot(legend_phylum_nasal)
setwd(directory.figures)
ggsave(filename = "legend_phylum_nasal.png", 
       ggpubr::as_ggplot(legend_phylum_nasal))


# Unlabeled tree
cxNasalTree = ggtree(treeN, 
                     layout = "circular", 
                     branch.length = 'none', 
                     size = 0.3) %<+% dN
cxNasalTree = cxNasalTree +  
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 18), legend.position = "none") +
  geom_point2(aes(subset = (!is.na(Status)), color = Status, fill = Status, shape = Status), alpha = 0.8, size = 14) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  scale_color_manual(values = c(pal[1], pal[3])) + 
  geom_cladelabel(node = which(idxsN=='node31'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#5E4FA2', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node100'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node150'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node168'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node217'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node281'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node385'), label = '',
                  align = T, geom = 'text', angle = 0, 
                  color = '#E6AB02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node396'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E6AB02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node414'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node446'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node452'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node489'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node572'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node664'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node757'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = 'grey', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node806'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#66A61E', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node864'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#66A61E', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node937'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#66A61E', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node126'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node200'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node218'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node231'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#E7298A', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node293'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node373'), label = '-', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#3288BD', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node428'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#9E0142', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node553'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node575'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node614'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node634'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node672'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node677'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#D95F02', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node730'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#3288BD', fontsize = 10, barsize = 5, offset = 0.8) +
  geom_cladelabel(node = which(idxsN=='node884'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#66A61E', fontsize = 10, barsize = 5, offset = 0.8) + 
  geom_cladelabel(node = which(idxsN=='node931'), label = '', 
                  align = T, geom = 'text', angle = 0, 
                  color = '#66A61E', fontsize = 10, barsize = 5, offset = 0.8)
cxNasalTree

setwd(directory.figures)
ggsave(filename = "full_nasal_tree_labeless.png", 
       cxNasalTree, 
       height = 40, 
       width = 40)

###############################################################################
#  Visualize significant microbial features nasal reduced                     #
###############################################################################

# Get a vector of descendant tips of significant clades
tipsN = vector()

for(i in 1:length(rootsN)){
  curNode = idxsN[which(idxsN == rootsN[i])]
  curClade = ape::extract.clade(treeN, curNode)
  tipsN = c(tipsN, curClade$tip.label)
  tipsN = unique(tipsN)
}
tipsN = unique(c(tipsN, sigNasal %>% 
                   filter(type == "ASV") %>% 
                   filter(covariate %in% c("Status")) %>%
                   pull(clade)))

# Drop unwanted tips
treeNSub = ape::drop.tip(phy = treeN, tip = setdiff(treeN$tip.label, tipsN))  
dNSub = dN[c(treeNSub$tip.label, treeNSub$node.label),]
dNSub$node = seq(from = 1, 
                 to = length(treeNSub$tip.label) + length(treeNSub$node.label),
                 by = 1)
idxsN = seq(from = 1, 
            to = length(treeNSub$tip.label) + length(treeNSub$node.label), 
            by = 1)
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
    
  ))

# Move labels outward from smaller clades
nodedf$size = NA
for(i in 1:nrow(nodedf)){
  curClade = ape::extract.clade(treeNSub, node = nodedf[i, "node"])
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
  geom_hilight(data=nodedf, 
               mapping=aes(node=node),
               extendto=2.3, 
               alpha=0.2, 
               fill="grey", 
               color="grey50",
               size=0.05) + 
  geom_hilight(data=nodedf, 
               mapping=aes(node=node),
               extendto=3.5, 
               alpha=0.3, 
               fill="white", 
               color="grey50",
               size=0.05) + 
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


# Output a labeled version for checking tips / labels match graph
setwd(directory.figures)
pdf(file = "tree_reduced_nasal_tip_node_labels.pdf")
cxNasalTree + geom_tiplab(size = 0.5) + geom_nodelab(size = 0.5)
dev.off()

# Decorate tree with the significant associations 
colnames(dNSub)[which(colnames(dNSub) == "Status" )] = "Association"
cxNasalTree = cxNasalTree %<+% dNSub + 
  geom_point2(aes(subset = (!is.na(Association)), color = Association), 
              alpha = 0.6, size = 2, shape = 16) +
  scale_color_manual(values = c(pal[1], pal[3])) 
cxNasalTree

# Create tiles 
tmpPhyloseq = phyloseq(otu_table(asvN_ps), 
                       phy_tree(asvN_ps), 
                       tax_table(asvN_ps), 
                       sample_data(asvN_ps))
tmpPhyloseq = prune_taxa(taxa = treeNSub$tip.label, 
                         x = tmpPhyloseq)
tmpPhyloseq = 
  prune_samples(samples = 
                  rownames(
                    sample_data(tmpPhyloseq)[which(sample_data(
                      tmpPhyloseq)$Status == "FURTD"),]), x = tmpPhyloseq)
furtdAve = apply(otu_table(tmpPhyloseq)[treeNSub$tip.label,], 1, mean)


furtdAve = data.frame("Average" = log(furtdAve + 0.001),
                      "Status" = rep("FURTD", times = length(furtdAve)))
furtdAve$ID = rownames(furtdAve)

tmpPhyloseq = phyloseq(otu_table(asvN_ps), 
                       phy_tree(asvN_ps), 
                       tax_table(asvN_ps), 
                       sample_data(asvN_ps))
tmpPhyloseq = prune_taxa(taxa = treeNSub$tip.label, x = tmpPhyloseq)
tmpPhyloseq = 
  prune_samples(samples = 
                  rownames(sample_data(
                    tmpPhyloseq)[which(sample_data(
                      tmpPhyloseq)$Status == "Control"),]), x = tmpPhyloseq)
controlAve = apply(otu_table(tmpPhyloseq)[treeNSub$tip.label,], 1, mean)


controlAve = data.frame("Average" = log(controlAve + 0.001),
                        "Status" = rep("Control", times = length(controlAve)))
controlAve$ID = rownames(controlAve)
datAve = rbind(furtdAve, controlAve)

# Plot tiles
cxNasalTree = cxNasalTree + 
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(data=datAve, geom=geom_tile,
             mapping=aes(y=ID, x=Status, fill=Status, alpha = Average),
             color = "grey50", offset = 0.31,size = 0.01) +
  scale_fill_manual(values = c(pal[1], pal[3], "grey")) 
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

dat3 = as.data.frame(dat3)
as.data.frame(dat3)
rownames(dat3) = dat3$ID

cxNasalTree = cxNasalTree + 
  ggtreeExtra::geom_fruit(data = dat3, 
                          geom = geom_bar, 
                          mapping = aes(y = ID, 
                                      x = HigherAve, 
                                      fill = HigherClass) , 
                          pwidth=0.38, offset = .1,orientation="y",  
                          stat="identity") 
cxNasalTree

cxNasalTree = cxNasalTree + 
  guides(color = guide_legend(override.aes = list(size=8))) +
  theme(legend.text = element_text(size = 9), 
        legend.title = element_text(size = 12), legend.position = "bottom") +
  guides(alpha=guide_legend(title="Log Average\nAbundance"))
cxNasalTree
cxNasalTreeLegend = ggpubr::get_legend(cxNasalTree)

# Plot reduced tree
setwd(directory.figures)
ggsave(cxNasalTree, file = "cladesClinicalSignsReducedNasal.pdf", 
       height = 25, 
       width = 23, 
       units = "cm")
ggsave(cxNasalTree, file = "cladesClinicalSignsReducedNasal.png", 
       height = 25, 
       width = 23, units = "cm")

###############################################################################
#  Calculate mean pairwise distances between sig clades                       #
###############################################################################

# Subset out nodes that are significantly associated with clinical signs 
sigCladesIdxG = seq(from = 1, 
                    to = length(c(treeG$tip.label, treeG$node.label)),
                    by = 1)
names(sigCladesIdxG) = c(treeG$tip.label, treeG$node.label)
sigCladesIdxN = seq(from = 1, 
                    to = length(c(treeN$tip.label, treeN$node.label)), 
                    by = 1)
names(sigCladesIdxN) = c(treeN$tip.label, treeN$node.label)

sigCladesNegativeG = sigCladesIdxG[sigGut %>% 
                                     filter(covariate %in% c("Status")) %>% 
                                     filter(Estimate < 0) %>% 
                                     pull(clade)]
sigCladesPositiveG = sigCladesIdxG[sigGut %>% 
                                     filter(covariate %in% c("Status")) %>% 
                                     filter(Estimate > 0) %>% pull(clade)]

sigCladesPositiveN = sigCladesIdxN[sigNasal %>% 
                                     filter(covariate %in% c("Status")) %>% 
                                     filter(Estimate > 0) %>% pull(clade)]
sigCladesNegativeN = sigCladesIdxN[sigNasal %>% 
                                     filter(covariate %in% c("Status")) %>% 
                                     filter(Estimate < 0) %>% pull(clade)]

# All pairwise distances between significant nodes 
combCladesPositiveG = combn(x = sigCladesPositiveG, m = 2)
combCladesNegativeG = combn(x = sigCladesNegativeG, m = 2)
combCladesPositiveN = combn(x = sigCladesPositiveN, m = 2)
combCladesNegativeN = combn(x = sigCladesNegativeN, m = 2)

# Distance between all nodes / tips on trees
distG = ape::dist.nodes(treeG)
distN = ape::dist.nodes(treeN)

# Distances between all significant nodes
distPositiveG = getDistanceCombos(combos = combCladesPositiveG, 
                                  distMatrix = distG)
distNegativeG = getDistanceCombos(combos = combCladesNegativeG, 
                                  distMatrix = distG)
distPositiveN = getDistanceCombos(combos = combCladesPositiveN, 
                                  distMatrix = distN)
distNegativeN = getDistanceCombos(combos = combCladesNegativeN, 
                                  distMatrix = distN)

# Calculate z score of mean distance between positively correlated nodes and random bootstrap
mean_positive_gut = mean(distPositiveG)
mean_positive_gut
sd(distPositiveG)
mean_negative_gut = mean(distNegativeG)
mean_negative_gut
sd(distNegativeG)
mean_positive_nasal = mean(distPositiveN)
mean_positive_nasal
sd(distPositiveN)
mean_negative_nasal = mean(distNegativeN)
mean_negative_nasal
sd(distNegativeN)

###############################################################################
#  Determine z-score mean pairwise distance vs. null dist                     #
###############################################################################

# Gut microbiome - nodes POSITIVELY correlated with clinical signs. 
N = 1000
random_mean_dists = vector()
for(i in 1:N){
  print(i)
  dist_cur = distRandom(dist = distG, n = length(distPositiveG))
  random_mean_dists = c(random_mean_dists, mean(dist_cur)) 
}
zscore_gut_positive = 
  (mean_positive_gut - mean(random_mean_dists)) / sd(random_mean_dists)
zscore_gut_positive
pnorm(q = zscore_gut_positive, mean = 0, sd = 1, lower.tail = TRUE)

random_mean_dists = data.frame("Pairwise_Distance" = random_mean_dists)
hist_gut_positive =  ggplot(random_mean_dists, aes(x = Pairwise_Distance)) + 
  geom_histogram(aes(y = ..density..))  +
  geom_density() +
  xlab("Average Pairwise Distance") +
  ylab("Density") +
  ggtitle("Gut Microbiome +") +
  xlim(c(1.1, 1.4)) +
  geom_vline(aes(xintercept = mean(Pairwise_Distance)), linetype = "dotted") +
  geom_vline(aes(xintercept = mean_positive_gut), 
             col = pal[3], 
             linetype = "solid", 
             size = 1, 
             alpha = 1) 
hist_gut_positive


# Gut microbiome - nodes NEGATIVELY correlated with clinical signs. 
N = 1000
random_mean_dists = vector()
for(i in 1:N){
  print(i)
  dist_cur = distRandom(dist = distG, n = length(distNegativeG))
  random_mean_dists = c(random_mean_dists, mean(dist_cur)) 
}
zscore_gut_negative = 
  (mean_negative_gut - mean(random_mean_dists)) / sd(random_mean_dists)
zscore_gut_negative
pnorm(q = zscore_gut_negative, mean = 0, sd = 1, lower.tail = TRUE)

random_mean_dists = data.frame("Pairwise_Distance" = random_mean_dists)
hist_gut_negative =  ggplot(random_mean_dists, aes(x = Pairwise_Distance)) + 
  geom_histogram(aes(y = ..density..))  +
  geom_density() +
  xlab("Average Pairwise Distance") +
  ggtitle("Gut Microbiome -") +
  ylab("Density") +
  geom_vline(aes(xintercept = mean(Pairwise_Distance)), linetype = "dotted") +
  geom_vline(aes(xintercept = mean_negative_gut), 
             col = pal[1], 
             linetype = "solid", 
             size = 1, 
             alpha = 1) +
  xlim(c(1.1, 1.4)) 
hist_gut_negative

# Nasal microbiome - nodes POSITIVELY correlated with clinical signs. 
random_mean_dists = vector()
for(i in 1:N){
  print(i)
  dist_cur = distRandom(dist = distN, n = length(distPositiveN))
  random_mean_dists = c(random_mean_dists, mean(dist_cur)) 
}
zscore_nasal_positive = 
  (mean_positive_nasal - mean(random_mean_dists)) / sd(random_mean_dists)
zscore_nasal_positive
pnorm(q = zscore_nasal_positive, mean = 0, sd = 1, lower.tail = FALSE)

random_mean_dists = data.frame("Pairwise_Distance" = random_mean_dists)
hist_nasal_positive = ggplot(random_mean_dists, aes(x = Pairwise_Distance)) + 
  geom_histogram(aes(y = ..density..))  +
  geom_density() +
  xlab("Average Pairwise Distance") +
  ggtitle("Nasal Microbiome +") +
  ylab("Density") +
  geom_vline(aes(xintercept = mean(Pairwise_Distance)), linetype = "dotted") +
  geom_vline(aes(xintercept = mean_positive_nasal), 
             col = pal[3], 
             linetype = "solid", 
             size = 1, 
             alpha = 1) +
  xlim(c(1.1, 1.4)) 
hist_nasal_positive

# Nasal microbiome - nodes NEGATIVELY correlated with clinical signs. 
random_mean_dists = vector()
for(i in 1:N){
  print(i)
  dist_cur = distRandom(dist = distN, n = length(distNegativeN))
  random_mean_dists = c(random_mean_dists, mean(dist_cur)) 
}
zscore_nasal_negative = 
  (mean_negative_nasal - mean(random_mean_dists)) / sd(random_mean_dists)
zscore_nasal_negative
pnorm(q = zscore_nasal_negative, mean = 0, sd = 1, lower.tail = FALSE)

random_mean_dists = data.frame("Pairwise_Distance" = random_mean_dists)
hist_nasal_negative = ggplot(random_mean_dists, aes(x = Pairwise_Distance)) + 
  geom_histogram(aes(y = ..density..))  +
  geom_density() +
  xlab("Average Pairwise Distance") +
  ggtitle("Nasal Microbiome -") +
  ylab("Density") +
  geom_vline(aes(xintercept = mean(Pairwise_Distance)), linetype = "dotted") +
  geom_vline(aes(xintercept = mean_negative_nasal), 
             col = pal[1], 
             linetype = "solid", 
             size = 1, 
             alpha = 1) +
  xlim(c(1.1, 1.4)) 
hist_nasal_negative

# Plot combined figure
hist_means = cowplot::plot_grid(hist_gut_positive, 
                                hist_gut_negative, 
                                hist_nasal_positive, 
                                hist_nasal_negative, 
                                labels = c("A", "B", "C", "D"))
setwd(directory.figures)
ggsave(filename = "average_pairwise_distance.pdf", hist_means)


###############################################################################
#  Pairwise distance null vs. significant clades                              #
###############################################################################


# DIST GUT microbiome - nodes POSITIVELY correlated with clinical signs. 
N = 1000
test_statistics_rand_positive = vector()
test_statistics_random_df = data.frame(matrix (nrow = length(distPositiveG), 
                                               ncol = N))

for(i in 1:N){
  
  # Generate a random set of pairwise distances
  rand_cur = distRandom(dist = distG, n = length(distPositiveG))
  test_statistics_random_df[,i] = rand_cur
  
  # Get test statistic between positive nodes and the current random stat
  ks_test_cur = ks.test(rand_cur, distPositiveG, exact = FALSE)
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_rand_positive = c(test_statistics_rand_positive, ks_test_cur)
  
}

# Randomly sample random distributions for computational tractability
test_statistics_random = vector()
test_statistics_random_combos = combn(x = N, m = 2)
test_statistics_random_combos_idx = 
  sample.int(n = ncol(test_statistics_random_combos), 
             size = length(test_statistics_rand_positive), 
             replace = FALSE)
test_statistics_random_combos = 
  test_statistics_random_combos[,test_statistics_random_combos_idx]
for(i in 1:ncol(test_statistics_random_combos)){
  ks_test_cur = ks.test(
    as.vector(test_statistics_random_df[,test_statistics_random_combos[1, i]]), 
    as.vector(test_statistics_random_df[,test_statistics_random_combos[2, i]]))
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_random = c(test_statistics_random, ks_test_cur)
}

ks.test(test_statistics_random, test_statistics_rand_positive)
ks.result = ks.test(test_statistics_random, test_statistics_rand_positive)
ks.result

# Plot
d_positive_random_df = 
  rbind(data.frame("D" = test_statistics_random, 
                   "Class" = rep("Random", 
                                 times = length(test_statistics_random))),
        data.frame("D" = test_statistics_rand_positive, 
                   "Class" = rep("Positive", 
                                 times = length(test_statistics_rand_positive))))
ks_test_D_statistics_histogram_gut_positive = 
  ggplot(d_positive_random_df, aes(x = `D`, fill = Class, color = Class)) +
  geom_histogram(position = "identity", bins = 50) +
  scale_fill_manual(
    values = c(colorspace::darken(col = pal[3], amount = 0.4), "grey")) +
  scale_color_manual(
    values = c(colorspace::darken(col = pal[3], amount = 0.6), 
               colorspace::darken("grey", amount = 0.2))) +
  xlab("D Statistic") +
  ylab("Frequency") 
ks_test_D_statistics_histogram_gut_positive
setwd(directory.figures)
ggplot(filename = "ks_test_D_statistics_histogram_gut_positive.pdf", 
       ks_test_D_statistics_histogram_gut_positive)


# DIST GUT microbiome - nodes NEGATIVELY correlated with clinical signs. 
N = 1000
test_statistics_rand_negative = vector()
test_statistics_random_df = data.frame(matrix(nrow = length(distNegativeG), 
                                              ncol = N))

for(i in 1:N){
  
  # Generate a random set of pairwise distances
  rand_cur = distRandom(dist = distG, n = length(distNegativeG))
  test_statistics_random_df[,i] = rand_cur
  
  # Get test statistic between positive nodes and the current random stat
  ks_test_cur = ks.test(rand_cur, distNegativeG, exact = FALSE)
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_rand_negative = c(test_statistics_rand_negative, ks_test_cur)
  
}

test_statistics_random = vector()
test_statistics_random_combos = combn(x = N, m = 2)
test_statistics_random_combos_idx = 
  sample.int(n = ncol(test_statistics_random_combos), 
             size = length(test_statistics_rand_negative), replace = FALSE)
test_statistics_random_combos = 
  test_statistics_random_combos[,test_statistics_random_combos_idx]
for(i in 1:ncol(test_statistics_random_combos)){
  ks_test_cur = 
    ks.test(as.vector(
      test_statistics_random_df[,test_statistics_random_combos[1, i]]),
      as.vector(
        test_statistics_random_df[,test_statistics_random_combos[2, i]]))
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_random = c(test_statistics_random, ks_test_cur)
}

ks.test(test_statistics_random, test_statistics_rand_negative)
ks.result = ks.test(test_statistics_random, test_statistics_rand_negative)
ks.result

# Plot
d_negative_random_df = 
  rbind(data.frame("D" = test_statistics_random, 
                   "Class" = rep("Random", 
                                 times = length(test_statistics_random))),
        data.frame("D" = test_statistics_rand_negative, 
                   "Class" = rep ("Negative", 
                                  times = length(test_statistics_rand_negative))))
ks_test_D_statistics_histogram_gut_negative = 
  ggplot(d_negative_random_df, aes(x = `D`, fill = Class, color = Class)) +
  geom_histogram(position = "identity", bins = 50) +
  scale_fill_manual(values = 
                      c(colorspace::darken(col = pal[1], amount = 0.4), "grey")) +
  scale_color_manual(values = 
                       c(colorspace::darken(col = pal[1], amount = 0.6), 
                         colorspace::darken("grey", amount = 0.2))) +
  xlab("D Statistic") +
  ylab("Frequency") 
ks_test_D_statistics_histogram_gut_negative
setwd(directory.figures)
ggplot(filename = "ks_test_D_statistics_histogram_gut_negative.pdf", 
       ks_test_D_statistics_histogram_gut_negative)



# DIST NASAL microbiome - nodes POSITIVELY correlated with clinical signs. 
N = 1000
test_statistics_rand_positive = vector()
test_statistics_random_df = data.frame(matrix(nrow = length(distPositiveN), 
                                              ncol = N))

for(i in 1:N){
  
  # Generate a random set of pairwise distances
  rand_cur = distRandom(dist = distN, n = length(distPositiveN))
  test_statistics_random_df[,i] = rand_cur
  
  # Get test statistic between positive nodes and the current random stat
  ks_test_cur = ks.test(rand_cur, distPositiveN, exact = FALSE)
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_rand_positive = c(test_statistics_rand_positive, ks_test_cur)
  
}

test_statistics_random = vector()
test_statistics_random_combos = combn(x = N, m = 2)
test_statistics_random_combos_idx = 
  sample.int(n = ncol(test_statistics_random_combos), 
             size = length(test_statistics_rand_positive), 
             replace = FALSE)
test_statistics_random_combos = 
  test_statistics_random_combos[,test_statistics_random_combos_idx]
for(i in 1:ncol(test_statistics_random_combos)){
  ks_test_cur = 
    ks.test(as.vector(
      test_statistics_random_df[,test_statistics_random_combos[1, i]]), 
      as.vector(
        test_statistics_random_df[,test_statistics_random_combos[2, i]]))
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_random = c(test_statistics_random, ks_test_cur)
}

ks.result = ks.test(test_statistics_random, test_statistics_rand_positive)

# Plot
d_positive_random_df = 
  rbind(data.frame("D" = test_statistics_random, 
                   "Class" = rep("Random", 
                                 times = length(test_statistics_random))),
        data.frame("D" = test_statistics_rand_positive, 
                   "Class" = rep ("Positive", 
                                  times = length(test_statistics_rand_positive))))
ks_test_D_statistics_histogram_nasal_positive = 
  ggplot(d_positive_random_df, aes(x = `D`, fill = Class, color = Class)) +
  geom_histogram(position = "identity", bins = 50) +
  scale_fill_manual(values = 
                      c(colorspace::lighten(col = pal[3], amount = 0.2), "grey")) +
  scale_color_manual(values = 
                       c(colorspace::lighten(col = pal[3], amount = 0.4), 
                         colorspace::darken(col = "grey", amount = 0.2))) +
  xlab("D Statistic") +
  ylab("Frequency") 
ks_test_D_statistics_histogram_nasal_positive
setwd(directory.figures)
ggplot(filename = "ks_test_D_statistics_histogram_nasal_positive.pdf", 
       ks_test_D_statistics_histogram_nasal_positive)


# Calculate if distributions are significantly different - nasal - negative
N = 1000
test_statistics_rand_negative = vector()
test_statistics_random_df = data.frame(matrix(nrow = length(distNegativeN), 
                                              ncol = N))

for(i in 1:N){
  
  # Generate a random set of pairwise distances
  rand_cur = distRandom(dist = distN, n = length(distNegativeN))
  test_statistics_random_df[,i] = rand_cur
  
  # Get test statistic between positive nodes and the current random stat
  ks_test_cur = ks.test(rand_cur, distNegativeN, exact = FALSE)
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_rand_negative = c(test_statistics_rand_negative, ks_test_cur)
  
}

test_statistics_random = vector()
test_statistics_random_combos = combn(x = N, m = 2)
test_statistics_random_combos_idx = 
  sample.int(n = ncol(test_statistics_random_combos), 
             size = length(test_statistics_rand_negative), replace = FALSE)
test_statistics_random_combos = 
  test_statistics_random_combos[,test_statistics_random_combos_idx]
for(i in 1:ncol(test_statistics_random_combos)){
  ks_test_cur = 
    ks.test(as.vector(
      test_statistics_random_df[,test_statistics_random_combos[1, i]]), 
      as.vector(
        test_statistics_random_df[,test_statistics_random_combos[2, i]]))
  ks_test_cur = as.vector(ks_test_cur$statistic)
  test_statistics_random = c(test_statistics_random, ks_test_cur)
}


ks.result = ks.test(test_statistics_random, test_statistics_rand_negative)

# Plot
d_negative_random_df = 
  rbind(data.frame("D" = test_statistics_random, 
                   "Class" = rep("Random", times = length(test_statistics_random))),
        data.frame("D" = test_statistics_rand_negative, 
                   "Class" = rep ("Negative", 
                                  times = length(test_statistics_rand_negative))))
ks_test_D_statistics_histogram_nasal_negative = 
  ggplot(d_negative_random_df, aes(x = `D`, fill = Class, color = Class)) +
  geom_histogram(position = "identity", bins = 50) +
  scale_fill_manual(values = 
                      c(colorspace::lighten(col = pal[1], amount = 0.2), "grey")) +
  scale_color_manual(values = 
                       c(colorspace::lighten(col = pal[1], amount = 0.4), 
                         colorspace::darken(col = "grey", amount = 0.2))) +
  xlab("D Statistic") +
  ylab("Frequency")

ks_test_D_statistics_histogram_nasal_negative
setwd(directory.figures)
ggplot(filename = "ks_test_D_statistics_histogram_nasal_negative.pdf", 
       ks_test_D_statistics_histogram_nasal_negative)

ks_test_D_statistics_histogram_gut_positive = 
  ks_test_D_statistics_histogram_gut_positive + 
  theme(legend.position = "none")
ks_test_D_statistics_histogram_gut_negative = 
  ks_test_D_statistics_histogram_gut_negative + 
  theme(legend.position = "none")
ks_test_D_statistics_histogram_nasal_positive = 
  ks_test_D_statistics_histogram_nasal_positive + 
  theme(legend.position = "none")
ks_test_D_statistics_histogram_nasal_negative = 
  ks_test_D_statistics_histogram_nasal_negative + 
  theme(legend.position = "none")

d_stat_legend = data.frame("Status" = c("Gut +", 
                                        "Gut -", 
                                        "Nasal +", 
                                        "Nasal -", 
                                        "Random"),
                           "Values" = c(1, 2, 3, 4, 5))
d_stat_legend$Status = 
  factor(d_stat_legend$Status, 
         levels = c("Gut +", "Gut -", "Nasal +", "Nasal -", "Random"))
d_stat_legend = ggplot(d_stat_legend, 
                       aes(x = Values, fill = Status, color = Status)) +
  geom_histogram(position = "identity") +
  scale_fill_manual(values = 
                      c(colorspace::darken(col = pal[3], amount = 0.4),
                        colorspace::darken(col = pal[1], amount = 0.4), 
                        colorspace::lighten(col = pal[3], amount = 0.2),
                        colorspace::lighten(col = pal[1], amount = 0.2), 
                        "grey")) + 
  scale_color_manual(values = 
                       c(colorspace::darken(col = pal[3], amount = 0.6),
                         colorspace::darken(col = pal[1], amount = 0.6),
                         colorspace::lighten(col = pal[3], amount = 0.4),
                         colorspace::lighten(col = pal[1], amount = 0.4),
                         colorspace::darken(col = "grey", amount = 0.2))) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
d_stat_legend = ggpubr::get_legend(d_stat_legend)
ggpubr::as_ggplot(d_stat_legend)


p_gut = cowplot::plot_grid(ks_test_D_statistics_histogram_gut_positive, 
                           ks_test_D_statistics_histogram_gut_negative, 
                           labels = c ("A", "B"))
p_nasal = cowplot::plot_grid(ks_test_D_statistics_histogram_nasal_positive, 
                             ks_test_D_statistics_histogram_nasal_negative, 
                             labels = c ("C", "D"))
p_legend = ggpubr::as_ggplot(d_stat_legend)

d_stat_plot = 
  cowplot::plot_grid(p_gut, p_legend, p_nasal, rel_widths = c(4, 1, 4))                   
d_stat_plot
setwd(directory.figures)
ggsave(filename = "d_stat_plot.pdf", d_stat_plot)

###############################################################################
#  Save Session Info                                                          #
###############################################################################

setwd(directory.out)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")




