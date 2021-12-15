# load libraries ----
library(tidyverse)
library(reshape2)
library(clusterProfiler)
library(cowplot)
library(jtools)
library(lme4)

# ggplot theme ----
my_theme <- theme(panel.grid = element_blank(),
                  panel.border = element_rect(fill = NA, color = 'black', size = 1),
                  panel.background = element_blank(),
                  axis.title.x.top = element_blank(),
                  axis.text.x.top = element_blank(),
                  axis.title.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.background = element_blank(),
                  axis.ticks.length=unit(2, "mm"))

# load data ----
proteins <- read_rds('data/breast_cells_proteins.rds')
genes <- read_rds('data/breast_cells_genes.rds')
growth <- read_rds('data/breast_cells_growth.rds')
hallmarks <- read_rds('data/hallmarks_sets.rds')

common_cells <- intersect(colnames(genes), colnames(proteins))
common_symbols <- intersect(rownames(genes), rownames(proteins))

genes <- genes[common_symbols, common_cells]
proteins <- proteins[common_symbols, common_cells]

dim(genes); dim(proteins)

all(colnames(genes) == colnames(genes))
all(rownames(genes) == rownames(proteins))

# gene protein correlation ----
# calculate the correlation between genes and proteins in all cells
gene_protein_corr <- cor(genes, proteins, use = 'complete') %>%
  melt() %>%
  as_tibble() %>%
  mutate(group = ifelse(Var1 == Var2, 'Same cell', 'Other cells'))

# test the correlation within cell type and in between
ll <- split(gene_protein_corr$value, gene_protein_corr$group)
ks.test(ll$`Other cells`, ll$`Same cell`)

# plot the correlations within cell type and in between
gene_protein_corr_plot <- gene_protein_corr %>%
  ggplot(aes(x = value, fill = group)) +
  geom_bar(stat="bin", aes(y=..density..), alpha = .5) +
  scale_x_continuous(sec.axis = dup_axis()) +
  scale_y_continuous(sec.axis = dup_axis()) +
  labs(x = 'Gene-protein Correlation',
       y = 'Percent of cell lines',
       fill = '') +
  theme(legend.position = c(.3, .8),
        legend.background = element_blank()) +
  my_theme
gene_protein_corr_plot

gene_protein_test <- ks.test(ll$`Other cells`, ll$`Same cell`)
cdf1 <- ecdf(ll$`Same cell`) 
cdf2 <- ecdf(ll$`Other cells`)
minMax <- seq(min(ll$`Same cell`, ll$`Other cells`), max(ll$`Same cell`, ll$`Other cells`), length.out=length(ll$`Same cell`)) 
x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
y0 <- cdf1(x0) 
y1 <- cdf2(x0) 

tt1 <- paste0(
  'D^+ = ',
  round(gene_protein_test$statistic, 2),
  '; P ',
  ifelse(gene_protein_test$p.value < .0001,
         '< 0.0001',
         paste0('= ', round(gene_protein_test$p.value, 3)))
)

gene_protein_test_plot <- gene_protein_corr %>%
  ggplot() +
  stat_ecdf(aes(x = value, color = group)) +
  geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
               linetype = "dashed" ) +
  geom_point(aes(x = x0[1] , y= y0[1]), size=3) +
  geom_point(aes(x = x0[1] , y= y1[1]), size=3) +
  ggtitle(tt1) + 
  labs(x = 'Gene-protein Correlation',
       y = 'ECDF',
       color = '') +
  theme(legend.position = c(.25,.85),
        legend.background = element_blank())

# calculate the correlation for the gene products in each hallmark
gene_protein_hallmark <- map_df(
  hallmarks,
  function(x) {
    as_tibble(
      melt(
        cor(
          genes[intersect(rownames(genes), x),],
          proteins[intersect(rownames(proteins), x),],
          use = 'complete'
        )
      )
    )
  }, .id = 'hallmark') %>%
  filter(Var1 == Var2)

gene_protein_hallmark_spearman <- map_df(
  hallmarks,
  function(x) {
    as_tibble(
      melt(
        cor(
          genes[intersect(rownames(genes), x),],
          proteins[intersect(rownames(proteins), x),],
          use = 'complete',
          method = 'spearman'
        )
      )
    )
  }, .id = 'hallmark') %>%
  filter(Var1 == Var2)

# generate a random set of genes
set.seed(123)
random_set <- sample(rownames(genes), 500)

# calculate the correlation for the random set
gene_protein_random <- cor(
  genes[intersect(rownames(genes), random_set),],
  proteins[intersect(rownames(proteins), random_set),],
  use = 'complete'
) %>%
  melt() %>%
  as_tibble() %>%
  filter(Var1 == Var2) %>%
  mutate(hallmark = 'random')

# test the correlation in the random set
ll <- split(gene_protein_hallmark$value, gene_protein_hallmark$hallmark)
ll$random <- gene_protein_random$value

# test the difference between the random set and hallmarks sets
df <- melt(ll) %>% as_tibble()
df$L1  <- relevel(factor(df$L1), ref = 'random')

plot(lm(df$value ~ df$L1))

cor_lm <- lm(df$value ~ df$L1)
summary(cor_lm)

corr_anova <- anova(cor_lm)

tt2 <- paste0(
  'F = ',
  round(corr_anova$`F value`[1], 2),
  '; P ',
  ifelse(corr_anova$`Pr(>F)`[1] < .0001,
         '< 0.0001',
         paste0('= ', round(corr_anova$`Pr(>F)`[1], 3)))
)

corr_lm_plot <- confint(cor_lm) %>%
  as.data.frame() %>%
  set_names(c('lower', 'upper')) %>%
  rownames_to_column('term') %>%
  left_join(broom::tidy(cor_lm)) %>%
  filter(term != '(Intercept)') %>%
  mutate(hallmark = as.numeric(fct_reorder(term, (estimate)))) %>%
  ggplot(aes(x = hallmark, y = estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = 1:10, labels = paste0('#', 1:10)) +
  coord_flip() +
  labs(x = 'Cancer Hallmarks',
       y = 'Gene-protein Correlation (Estimate)') +
  ggtitle(tt2)
  
# enrichment of hallmarks in breast cancer cells (genes) ----
# term_gene <- melt(hallmarks) %>%
#   as_tibble() %>%
#   dplyr::select(term = L1, gene = value)
# 
# stats <- lapply(seq_len(ncol(genes)), function(i) genes[,i])
# stats <- map(stats, sort, decreasing = TRUE)
# names(stats) <- colnames(genes)
# 
# enrich_hallmark_genes <- map_df(stats, function(x) {
#   enrich <- GSEA(x,
#                  pAdjustMethod = 'fdr',
#                  TERM2GENE = term_gene,
#                  pvalueCutoff = 1,
#                  maxGSSize = Inf)
#   as_tibble(enrich)
# }, .id = 'cell_id')
# 
# write_rds(enrich_hallmark_genes, 'data/enrich_hallmark_genes.rds')

enrich_hallmark_genes <- read_rds('data/enrich_hallmark_genes.rds')

correlation_ernichment <- gene_protein_hallmark %>%
  full_join(enrich_hallmark_genes,
            by = c('hallmark' = 'ID',
                   'Var1' = 'cell_id')) %>%
  dplyr::select(hallmark, Var1, value, NES) %>%
  mutate(hallmark = as.numeric(fct_reorder(hallmark, (value))))

hallmark_enrichment_plot <- correlation_ernichment %>%
  ggplot(aes(x = hallmark, y = NES, group = hallmark)) +
  geom_boxplot() +
  labs(x = 'Cancer Hallmarks',
       y = 'Enrichment Score') +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = 1:10, labels = paste0('#', 1:10)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme

hallmark_correlation_plot <- correlation_ernichment %>%
  ggplot(aes(x = hallmark, y = value, group = hallmark)) +
  geom_boxplot() +
  labs(x = 'Cancer Hallmarks', y = 'Gene-protein Correlation') +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = 1:10, labels = paste0('#', 1:10)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme


hallmark_enrichment_points <- enrich_hallmark_genes %>%
  mutate(ID = as.numeric(as.factor(ID))) %>%
  ggplot(aes(x = ID, y = NES, color = p.adjust < .2)) +
  geom_point() +
  scale_x_continuous(breaks = 1:10, labels = paste0('#', 1:10)) +
  labs(x = 'Cancer Hallmarks',
       y = 'Enrichment Score',
       color = 'FDR < 0.2: ') +
  theme(legend.position = 'top')

# enrichment of hallmarks in breast cancer cells (proteins) ----
# stats <- lapply(seq_len(ncol(proteins)), function(i) proteins[,i])
# stats <- map(stats, sort, decreasing = TRUE)
# names(stats) <- colnames(proteins)
# 
# enrich_hallmark_proteins <- map_df(stats, function(x) {
#   enrich <- GSEA(x,
#                  pAdjustMethod = 'fdr',
#                  TERM2GENE = term_gene,
#                  pvalueCutoff = 1,
#                  maxGSSize = Inf)
#   as_tibble(enrich)
# }, .id = 'cell_id')
# 
# write_rds(enrich_hallmark_proteins, 'data/enrich_hallmark_proteins.rds')

enrich_hallmark_proteins <- read_rds('data/enrich_hallmark_proteins.rds')

summary(enrich_hallmark_genes$NES)
table(enrich_hallmark_genes$p.adjust < .2)
table(enrich_hallmark_genes$ID, enrich_hallmark_genes$p.adjust < .2)
table(enrich_hallmark_genes$cell_id, enrich_hallmark_genes$p.adjust < .2)

# difference in enrichment ----
enrichment_difference <- full_join(enrich_hallmark_genes,
                                   enrich_hallmark_proteins,
                                   by = c('cell_id', 'ID', 'Description')) %>%
  mutate(ID = fct_reorder(ID, NES.x),
         difference = NES.y - NES.x) 

enrichment_difference_plot <- enrichment_difference %>%
  ggplot(aes(x = difference, y = as.numeric(ID), group = ID)) +
  geom_point() +
  labs(x = 'Difference in enrichment', y = 'Cancer Hallmarks') +
  scale_y_continuous(sec.axis = dup_axis(),
                     breaks = 1:10, labels = paste0('#', 1:10)) +
  scale_x_continuous(sec.axis = dup_axis()) +
  my_theme

plot_grid(hallmark_enrichment_plot,
          enrichment_difference_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/hallmark_enrichment.png',
         width = 6, height = 3)

plot_grid(gene_protein_corr_plot,
          hallmark_correlation_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/gene_protein_correlation.png',
         width = 6, height = 3)

plot_grid(gene_protein_test_plot,
          corr_lm_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/gene_protein_correlation_tests.png',
         width = 6, height = 3)

# replication ----
difference_division  <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  group_by(cell_id, ID) %>%
  summarise(difference = unique(difference),
            ndr = mean(Nominal.Division.Rate))

# model the relation between division and difference in enrichment
df <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name'))

division_lm <- lme(Nominal.Division.Rate ~ difference , random = ~1|ID, data = df)
summary(division_lm)

division_lm_list <- lmList(Nominal.Division.Rate ~ difference | ID, data = df)
division_conf <- confint(division_lm_list)
summary(division_lm_list) 

tt3 <- paste0(
  'FE = ',
  round(unclass(summary(division_lm))$tTable[2,'Value'], 2),
  '; P ',
  ifelse(corr_anova$`Pr(>F)`[1] < .0001,
         '< 0.0001',
         paste0('= ', round(unclass(summary(division_lm))$tTable[2,'p-value'], 3)))
)

division_test_plot <- bind_cols(
  as.data.frame(unclass(summary(division_lm_list))$coefficients),
  as.data.frame(unclass((division_conf)))
) %>%
  rownames_to_column('hallmark') %>%
  mutate(hallmark = as.numeric(fct_reorder(hallmark, (Estimate.difference)))) %>%
  ggplot(aes(x = hallmark, y = Estimate.difference)) +
  geom_point() +
  geom_linerange(aes(ymin = `2.5 %.difference`, ymax = `97.5 %.difference`)) +
  scale_x_continuous(breaks = 1:10, labels = paste0('#', 1:10)) +
  coord_flip() +
  labs(x = 'Cancer Hallmarks',
       y = 'Division Rate (Estimate)') +
  ggtitle(tt3)
  
difference_division_plot <- difference_division %>%
  ggplot(aes(x = difference, y = ndr)) +
  geom_density2d_filled()+
  labs(x = 'Difference in enrichment',
       y = 'Division Rate') +
  theme(legend.position = 'none') +
  scale_x_continuous(sec.axis = dup_axis(),
                     expand = c(0, 0)) +
  scale_y_continuous(sec.axis = dup_axis(), expand = c(0, 0)) +
  my_theme

# resilience ----
difference_inhibition <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  group_by(cell_id, ID, Small.Molecule.Name) %>%
  summarise(difference = unique(difference),
            grmax = mean(Normalized.Growth.Rate.Inhibition.Value))

difference_inhibition_plot <- difference_inhibition %>%
  ggplot(aes(x = difference, y = grmax)) +
  geom_density2d_filled() +
  labs(x = 'Difference in enrichment',
       y = 'Growth Inhibition') +
  theme(legend.position = 'none') +
  scale_x_continuous(sec.axis = dup_axis(),
                     expand = c(0, 0)) +
  scale_y_continuous(sec.axis = dup_axis(),expand = c(0, 0)) +
  my_theme

plot_grid(difference_division_plot,
          difference_inhibition_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/division_inhibition.png',
         width = 6, height = 3)

cor(difference_division$difference, difference_division$ndr)
cor.test(difference_division$difference, difference_division$ndr)

cor(difference_inhibition$difference, difference_inhibition$grmax)
cor.test(difference_inhibition$difference, difference_inhibition$grmax)

# cost and benefit ----
# the higher the difference the higher the cost, the benefit is being robust and resilient
difference_growth <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  group_by(cell_id, ID, Small.Molecule.Name) %>%
  summarise(difference = unique(difference),
            grmax = min(Normalized.Growth.Rate.Inhibition.Value),
            ndr = mean(Nominal.Division.Rate)) %>%
  ungroup()

division_growth <- difference_growth %>%
  group_by(ID) %>%
  summarise('Division Rate' = cor(difference, ndr),
            'Growth Inhibition' = cor(difference, grmax)) %>%
  mutate(ID = as.numeric(fct_reorder(ID, `Division Rate`)),
         ID = paste0('#', as.factor(as.numeric(ID))),
         ID = factor(ID, levels = paste0('#', 1:10))) %>%
  gather(cor, value, `Division Rate`, `Growth Inhibition`) %>%
  ggplot(aes(x = value, y = ID)) +
  geom_col() +
  facet_wrap(~cor) +
  scale_x_continuous(breaks = -c(0, .25, .5), labels = c('0', '-0.25', '-0.5')) +
  labs(x = "Pearson's Correlation",
       y = 'Cancer Hallmarks') +
  my_theme +
  theme(panel.spacing = unit(-.2, 'mm'),
        strip.background = element_blank())

division_growth
# ggsave(division_growth,
#        filename = 'output/manuscript/figures/division_growth.png',
#        width = 6, height = 3)

difference_growth %>%
  group_by(ID) %>%
  summarise(dr = cor(difference, ndr),
            dr_pval = cor.test(difference, ndr)$p.value,
            gr = cor(difference, grmax),
            gr_pval = cor.test(difference, grmax)$p.value)

# prediction
enrichment_range <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  ungroup() %>%
  group_by(cell_id, Small.Molecule.Name) %>%
  summarise(
    gene_min = min(NES.x),
    gene_max = max(NES.x),
    gene_diff = gene_max - gene_min,
    protein_min = min(NES.y),
    protein_max = max(NES.y),
    protein_diff = protein_max - protein_min,
    ndr_max = max(Relative.Cell.Count),
    ndr_med = median(Nominal.Division.Rate),
    grmax_med = median(Normalized.Growth.Rate.Inhibition.Value)
  ) %>%
  unique()

# range of enrichment correlation with grmax_med
range_ndr_correaltion <- enrichment_range %>%
  group_by(Small.Molecule.Name) %>%
  summarise(Genes = cor(gene_diff, ndr_med),
            Proteins = cor(protein_diff, ndr_med))

enrichment_range_plot <- range_ndr_correaltion %>%
  gather(key, value, -Small.Molecule.Name) %>%
  ggplot(aes(x = value, fill = key)) +
  geom_bar(stat="bin", aes(y=..density..), alpha = .5) +
  labs(x = "Pearson's Correlation",
       y = 'Percent of drugs',
       fill = '') +
  theme(legend.position = c(.5,.85),
        legend.background = element_blank()) +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = c(-.6,-.3, 0, .3)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme

range_ndr_correaltion %>%
  gather(key, value, -Small.Molecule.Name) %>%
  ggplot(aes(x = value, color = key)) +
  stat_ecdf()

plot_grid(division_growth,
          enrichment_range_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/difference_range.png',
         width = 6, height = 3)

set.seed(123)
ks.test(range_ndr_correaltion$Genes, 'rnorm', alternative = 'less')
ks.test(range_ndr_correaltion$Proteins, 'rnorm', alternative = 'greater')

range_ndr_test <- ks.test(range_ndr_correaltion$Genes, range_ndr_correaltion$Proteins)
cdf1 <- ecdf(range_ndr_correaltion$Genes) 
cdf2 <- ecdf(range_ndr_correaltion$Proteins)

minMax <- seq(min(range_ndr_correaltion$Genes, range_ndr_correaltion$Proteins, na.rm = TRUE),
              max(range_ndr_correaltion$Genes, range_ndr_correaltion$Proteins, na.rm = TRUE),
              length.out=length(range_ndr_correaltion$Genes)) 

x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
y0 <- cdf1(x0) 
y1 <- cdf2(x0) 

tt4 <- paste0(
  'D^+ = ',
  round(range_ndr_test$statistic, 2),
  '; P ',
  ifelse(range_ndr_test$p.value < .0001,
         '< 0.0001',
         paste0('= ', round(range_ndr_test$p.value, 3)))
)

range_ndr_test_plot <- range_ndr_correaltion %>%
  gather(key, value, -Small.Molecule.Name) %>%
  ggplot() +
  stat_ecdf(aes(x = value, color = key)) +
  geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
               linetype = "dashed" ) +
  geom_point(aes(x = x0[1] , y= y0[1]), size=3) +
  geom_point(aes(x = x0[1] , y= y1[1]), size=3) +
  ggtitle(tt4) + 
  labs(x = 'Gene-protein Correlation',
       y = 'ECDF',
       color = '') +
  theme(legend.position = c(.47,.5),
        legend.background = element_blank())

division_test_plot
range_ndr_test_plot

plot_grid(division_test_plot,
          range_ndr_test_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/difference_range_tests.png',
         width = 6, height = 3)

# dose_groups <- enrichment_difference %>%
#   na.omit() %>%
#   inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
#   filter(ID == 'Resisting Cell Death') %>%
#   ungroup() %>%
#   mutate(dose = cut(scale(Perturbagen.Concentration),
#                     3,
#                     labels = c('L', 'M', 'H')),
#          gene = cut(NES.x, 2, labels = c('Lower', 'Higher')),
#          protein = cut(NES.y, 2, labels = c('Lower', 'Higher'))) %>%
#   group_by(cell_id, Small.Molecule.Name, gene, protein, dose) %>%
#   summarise(grmax = min(Normalized.Growth.Rate.Inhibition.Value))
# 
# dose_groups_plot <- dose_groups %>%
#   group_by(gene, protein, dose) %>%
#   summarise(ave = mean(grmax),
#             # sd = sd(grmax)/sqrt(n()),
#             sd = sd(grmax),
#             upper = ave + sd,
#             lower = ave - sd) %>%
#   ggplot(aes(x = dose, y = ave, ymax = upper, ymin = lower)) +
#   geom_point() +
#   geom_errorbar() +
#   facet_grid(protein~gene, ) +
#   labs(x = 'Dose', y = 'Growth Inhinition') +
#   theme(legend.position = 'none')
# 
# dose_groups <- enrichment_difference %>%
#   na.omit() %>%
#   inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
#   ungroup() %>%
#   group_by(Small.Molecule.Name) %>%
#   mutate(dose = cut((Perturbagen.Concentration),
#                     3,
#                     labels = paste0(1:3))) %>%
#   ungroup() %>%
#   mutate(gene = NES.x,
#          protein = NES.y,
#          # gene = cut(NES.x, 2, labels = c('L', 'H')),
#          # protein = cut(NES.y, 2, labels = c('L', 'H')),
#          grmax = Normalized.Growth.Rate.Inhibition.Value,
#          drug = Small.Molecule.Name) %>%
#   group_by(cell_id, drug, gene, protein, dose) %>%
#   summarise(grmax = mean(Normalized.Growth.Rate.Inhibition.Value))

# dose_groups_plot <- dose_groups %>%
#   ggplot() +
#   geom_smooth(aes(x = scale(protein), y= grmax, color = 'Proteins'), method = 'lm') +
#   geom_smooth(aes(x = scale(gene), y= grmax, color = 'Genes'), method = 'lm') +
#   labs(x = 'Enrichment score (predictor)',
#        y = 'Growth Inhibition (response)',
#        color = '') +
#   theme(legend.position = c(.5, .85),
#         legend.background = element_blank())

# plot_grid(enrichment_range_plot,
#           dose_groups_plot,
#           nrow = 1,
#           scale = .9,
#           labels = 'AUTO',
#           label_fontface = 'plain') %>%
#   ggsave(plot = .,
#          filename = 'output/manuscript/figures/enrichment_range.png',
#          width = 6, height = 3)

# # range of enrichment and response per cell
# ind <- enrichment_range %>%
#   ungroup() %>%
#   dplyr::select(Small.Molecule.Name, grmax_med) %>%
#   unique() %>%
#   arrange(grmax_med) %>%
#   pull(Small.Molecule.Name)
# 
# p1 <- enrichment_range %>%
#   group_by(Small.Molecule.Name) %>%
#   mutate(Small.Molecule.Name = factor(Small.Molecule.Name, ind)) %>%
#   ggplot(aes(x = grmax_med, y = Small.Molecule.Name)) +
#   geom_boxplot() +
#   labs(x = 'Growth Inhibition (Median)', y = '')
# 
# p2 <- enrichment_range %>%
#   dplyr::select(Small.Molecule.Name, protein_diff) %>%
#   unique() %>%
#   mutate(cell_id = factor(cell_id, ind)) %>%
#   ggplot(aes(y = Small.Molecule.Name, x = protein_diff)) +
#   geom_col() +
#   labs(x = 'Range (Protein)', y = '') +
#   scale_y_discrete(labels = rep('', 11))
# 
# p3 <- enrichment_range %>%
#   dplyr::select(Small.Molecule.Name, gene_diff) %>%
#   unique() %>%
#   mutate(Small.Molecule.Name = factor(Small.Molecule.Name, ind)) %>%
#   ggplot(aes(y = Small.Molecule.Name, x = gene_diff)) +
#   geom_col() +
#   labs(x = 'Range (Gene)', y = '') +
#   scale_y_discrete(labels = rep('', 11))
# 
# plot_grid(p1, p3, p2, nrow = 1, rel_widths = c(2, 1, 1)) %>%
#   ggsave(plot = .,
#          filename = 'output/manuscript/figures/ranges.png',
#          width = 8, height = 3)
