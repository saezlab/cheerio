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


## ensembl way
library(tidyverse)
library(biomaRt)

human.R95 <- useMart(host = "jan2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
hs <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
rn <- useMart('ensembl', dataset = 'rnorvegicus_gene_ensembl')
mm<- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

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


##liana

library(tidyverse)
library(OmnipathR)
library(liana)
library(magrittr)

liana_path <- system.file(package = "liana")
testdata <-
  readRDS(file.path(liana_path, "testdata", "input", "testdata.rds"))

# Convert testdata to putative mouse symbols
# Note that we explicitly provide this tutorial to avoid any such conversions when working with real data
# We simply do this here to provide an example
rownames(testdata@assays$RNA@counts) <- stringr::str_to_title(rownames(testdata@assays$RNA@counts))
rownames(testdata@assays$RNA@data) <- stringr::str_to_title(rownames(testdata@assays$RNA@data))


## orthogene
BiocManager::install("orthogene")
