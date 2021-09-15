# load libraries ----
library(tidyverse)
library(reshape2)
library(clusterProfiler)
library(cowplot)

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
gene_protein_corr <- melt(cor((genes), (proteins), use = 'complete')) %>%
  as_tibble() %>%
  mutate(group = ifelse(Var1 == Var2, 'Same cell', 'Other cells'))

ll <- split(gene_protein_corr$value, gene_protein_corr$group)
ks.test(ll$`Other cells`, ll$`Same cell`)

gene_protein_corr_plot <- gene_protein_corr %>%
  ggplot(aes(x = value, fill = group)) +
  geom_bar(stat="bin", aes(y=..density..), alpha = .5) +
  labs(x = 'Transcription-translation Correlation',
       y = 'Percent of cell lines',
       fill = '') +
  theme(legend.position = c(.3, .8),
        legend.background = element_blank())

gene_protein_hallmark <- map_df(hallmarks, function(x) {
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
lengths(hallmarks)

set.seed(123)
random_set <- sample(rownames(genes), 500)
gene_protein_random <- cor(genes[intersect(rownames(genes), random_set),],
                           proteins[intersect(rownames(proteins), random_set),],
                           use = 'complete') %>%
  melt() %>%
  as_tibble() %>%
  filter(Var1 == Var2) %>%
  mutate(hallmark = 'random')

ll <- split(gene_protein_hallmark$value, gene_protein_hallmark$hallmark)
ll$random <- gene_protein_random$value

wt <- wilcox.test(ll$random, ll$`Activating Invasion and Metastasis`,
            paired = TRUE, alternative = 'less')
df <- melt(ll)
df$L1 <- relevel(factor(df$L1), ref = 'random')
cor_lm <- lm(df$value ~ df$L1)
broom::tidy(cor_lm)
summary(cor_lm)

anova(cor_lm)

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
  scale_x_continuous(breaks = 1:10, labels = paste0('#', 1:10))

hallmark_correlation_plot <- correlation_ernichment %>%
  ggplot(aes(x = hallmark, y = value, group = hallmark)) +
  geom_boxplot() +
  labs(x = 'Cancer Hallmarks', y = 'Transcription-translation Correlation') +
  scale_x_continuous(breaks = 1:10, labels = paste0('#', 1:10))


plot_grid(hallmark_enrichment_plot,
          enrichment_difference_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/hallmark_enrichment.png',
         width = 6, height = 3)

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
enrichment_difference <- full_join(enrich_hallmark_genes, enrich_hallmark_proteins,
                                   by = c('cell_id', 'ID', 'Description')) %>%
  mutate(ID = fct_reorder(ID, NES.x),
         difference = NES.x - NES.y) 

enrichment_difference_plot <- enrichment_difference %>%
  ggplot(aes(x = difference, y = as.numeric(ID), group = ID)) +
  geom_point() +
  labs(x = 'Difference in enrichment', y = 'Cancer Hallmarks') +
  scale_y_continuous(breaks = 1:10, labels = paste0('#', 1:10))

plot_grid(gene_protein_corr_plot,
          hallmark_correlation_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/gene_protein_correlation.png',
         width = 6, height = 3)

# replication ----
difference_division  <- enrichment_difference %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  group_by(cell_id, ID) %>%
  summarise(difference = unique(difference),
            ndr = mean(Nominal.Division.Rate))

difference_division_plot <- difference_division %>%
  ggplot(aes(x = difference, y = ndr)) +
  geom_density2d_filled()+
  labs(x = 'Difference in enrichment',
       y = 'Division Rate') +
  theme(legend.position = 'none')

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
  theme(legend.position = 'none')

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
  mutate(ID = fct_reorder(ID, `Division Rate`)) %>%
  gather(cor, value, `Division Rate`, `Growth Inhibition`) %>%
  ggplot(aes(x = value, y = ID)) +
  geom_col() +
  facet_wrap(~cor) +
  labs(x = "Pearson's Correlation", y = '')

ggsave(division_growth,
       filename = 'output/manuscript/figures/division_growth.png',
       width = 6, height = 3)

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
        legend.background = element_blank())

range_ndr_correaltion %>%
  gather(key, value, -Small.Molecule.Name) %>%
  ggplot(aes(x = value, color = key)) +
  stat_ecdf()

set.seed(123)
ks.test(range_ndr_correaltion$Genes, 'rnorm', alternative = 'less')
ks.test(range_ndr_correaltion$Proteins, 'rnorm', alternative = 'greater')

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

dose_groups <- enrichment_difference %>%
  na.omit() %>%
  inner_join(growth, by = c('cell_id' = 'Cell.Name')) %>%
  ungroup() %>%
  group_by(Small.Molecule.Name) %>%
  mutate(dose = cut((Perturbagen.Concentration),
                    3,
                    labels = paste0(1:3))) %>%
  ungroup() %>%
  mutate(gene = NES.x,
         protein = NES.y,
         # gene = cut(NES.x, 2, labels = c('L', 'H')),
         # protein = cut(NES.y, 2, labels = c('L', 'H')),
         grmax = Normalized.Growth.Rate.Inhibition.Value,
         drug = Small.Molecule.Name) %>%
  group_by(cell_id, drug, gene, protein, dose) %>%
  summarise(grmax = mean(Normalized.Growth.Rate.Inhibition.Value))

dose_groups_plot <- dose_groups %>%
  ggplot() +
  geom_smooth(aes(x = scale(protein), y= grmax, color = 'Proteins'), method = 'lm') +
  geom_smooth(aes(x = scale(gene), y= grmax, color = 'Genes'), method = 'lm') +
  labs(x = 'Enrichment score (predictor)',
       y = 'Growth Inhibition (response)',
       color = '') +
  theme(legend.position = c(.5, .85),
        legend.background = element_blank())

plot_grid(enrichment_range_plot,
          dose_groups_plot,
          nrow = 1,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain') %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/enrichment_range.png',
         width = 6, height = 3)

summary(lm(grmax ~ drug + gene + protein + dose, data = dose_groups))

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
