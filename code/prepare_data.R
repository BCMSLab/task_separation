# load libraries
library(tidyverse)
library(reshape2)
library(slinky)
library(org.Hs.eg.db)

# # Breast Cancer Profiling Project, Gene Expression 1: Baseline mRNA sequencing on 35 breast cell lines - Dataset (ID:20348)
# url <- 'http://lincs.hms.harvard.edu/data/HMS_Dataset_20348.zip'
# download.file(url, destfile = 'data/HMS_Dataset_20348.zip')
# unzip('data/HMS_Dataset_20348.zip', exdir = 'data')
# 
# # Breast Cancer Profiling Project – Proteomics 1: 1 total proteome dataset for a 35-cell line breast cancer panel under basal conditions - Dataset (ID:20352)
# url <- 'http://lincs.hms.harvard.edu/data/HMS_Dataset_20352.zip'
# download.file(url, destfile = 'data/HMS_Dataset_20352.zip')
# unzip('data/HMS_Dataset_20352.zip', exdir = 'data')
# 
# # Growth rate-corrected (GR) dose-response metrics across a panel of 71 breast cancer cell lines treated with a library of small molecule and antibody perturbagens. Dataset 1 of 4: Relative cell counts and normalized growth rate inhibition values across technical replicates. - Dataset (ID:20268)
# url <- 'https://lincs.hms.harvard.edu/db/datasets/20268/results?search=&output_type=.csv'
# download.file(url, destfile = 'data/HMS_Dataset_20268.csv')

# gene expression
breast_cells_genes_df <- read_csv('data/HMS_Dataset_20348/HMS_Dataset_20348_DataFile.csv')
breast_cells_genes <- as.matrix(breast_cells_genes_df[, -1])
rownames(breast_cells_genes) <- breast_cells_genes_df$id
breast_cells_genes <- breast_cells_genes[!duplicated(rownames(breast_cells_genes)),]

write_rds(breast_cells_genes, 'data_clean/breast_cells_genes.rds')

rm(breast_cells_genes_df)

# proteoms
breast_cells_proteins_df <- read_csv('data/HMS_Dataset_20352/HMS_Dataset_20352_DataFile.csv')
breast_cells_proteins <- as.matrix(breast_cells_proteins_df[, -c(1,2)])
rownames(breast_cells_proteins) <- breast_cells_proteins_df$Gene_symbol
breast_cells_proteins <- breast_cells_proteins[!duplicated(rownames(breast_cells_proteins)),]

write_rds(breast_cells_proteins, 'data_clean/breast_cells_proteins.rds')

rm(breast_cells_proteins_df)

# growth rate
breast_cells_growth <- read_csv('data/HMS_Dataset_20268.csv', name_repair = 'universal')

write_rds(breast_cells_growth, 'data_clean/breast_cells_growth.rds')

# hallmarks
# source: CHG: A Systematically Integrated Database of Cancer Hallmark Genes
# source: http://bio-bigdata.hrbmu.edu.cn/CHG/index.html
url <- 'http://bio-bigdata.hrbmu.edu.cn/CHG/download/Supplementary%20Table%202%20Genesets%20of%2010%20hallmarker.xlsx'

download.file(url, 'data/hallmarks_gene_sets.xlsx')
symbols <- keys(org.Hs.eg.db, 'SYMBOL')

hallmarks <- readxl::read_excel('data/hallmarks_gene_sets.xlsx',
                   skip = 2) %>%
  as.list() %>%
  map(function(x) {
    l <- x[!is.na(x)]
    l <- toupper(l)
    intersect(l, symbols)
  })

write_rds(hallmarks, 'data_clean/hallmarks_sets.rds')

