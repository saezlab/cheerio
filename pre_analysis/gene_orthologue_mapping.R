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
hs= useEnsembl("ensembl", dataset = 'hsapiens_gene_ensembl', mirror= "useast" )
mm= useEnsembl("ensembl", dataset = 'mmusculus_gene_ensembl', mirror= "useast" )

hs <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', mirror= "useast")
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
show_homologene()


symbols_dict <- liana::get_homologene_dict(entities = entities, 
                                    target_organism = target_organism)

op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse

# Run LIANA with the orthologous resource
liana_res <- liana_wrap(testdata,
                        resource = 'custom', # resource has to be set to 'custom' to work with external resources
                        external_resource = ortholog_resource, # provide orthologous resource
                        method=c('sca', 'natmi') # run only with sca and natmi for comp. time
)


## orthogene
BiocManager::install("liana")


# rough way go through decoupler py. ---------------------------------------------------------------

contrast_df = readRDS("data/contrasts_query_df.rds")

gene.tbl = contrast_df %>% 
  ungroup()%>%
  mutate(gene = ifelse(grepl("Rn|Mm",gene), str_to_title(gene), gene),
         spec= substr(contrast_id,1,2))%>% 
  #select(-contrast_id)%>%
  distinct(gene, spec)

gene.tbl%>% write_csv("data/gene_names.csv")
# 
# rat_genes= contrast_df %>% 
#   filter(grepl("Rn", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()
# 
# mouse_genes= contrast_df %>% 
#   filter(grepl("Mm", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()
# 
# human_genes= contrast_df %>% 
#   filter(grepl("Hs", contrast_id))%>% pull(gene) %>% unique()
# 


