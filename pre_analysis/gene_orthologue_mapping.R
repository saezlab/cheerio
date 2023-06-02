## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-06-02
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## create ortholog mapping for all contrasts






library(tidyverse)
library(biomaRt)


hs <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
mm<- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
rn <- useMart('ensembl', dataset = 'rnorvegicus_gene_ensembl')

##mm - rt

mm_rt <- getLDS(
  mart = mm,
  attributes = c('ensembl_gene_id',"mgi_symbol", 'external_gene_name'),
  martL = rn,
  attributesL = c('ensembl_gene_id','external_gene_name'),
  filters = 'external_gene_name',
  values = rat_genes)

## ms -hs 

mm_hs <- getLDS(
  mart = mm,
  attributes = c('ensembl_gene_id',"mgi_symbol",'external_gene_name'),
  martL = hs,
  attributesL = c('ensembl_gene_id','external_gene_name'))

saveRDS(list("mm_rt"= mm_rt,
             "mm_hs"= mm_hs))
