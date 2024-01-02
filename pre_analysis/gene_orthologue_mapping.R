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

human.R95 <- useMart(host = "https://jan2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
hs= useEnsembl("ensembl", dataset = 'hsapiens_gene_ensembl',host="www.ensembl.org")
mm= useEnsembl("ensembl", dataset = 'mmusculus_gene_ensembl',host="www.ensembl.org")

hs <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl',host="https://www.ensembl.org")
rn <- useMart('ensembl', dataset = 'rnorvegicus_gene_ensembl',host="https://www.ensembl.org")
mm<- useMart('ensembl', dataset = 'mmusculus_gene_ensembl',host="https://www.ensembl.org")

##mm - rt

mm_rt <- getLDS(
  #mart = mm,
  #attributes = c('ensembl_gene_id',"mgi_symbol", 'external_gene_name'),
  mart = rn,
  attributes = c('ensembl_gene_id','external_gene_name'),
  martL = hs,
  attributesL = c('ensembl_gene_id','hgnc_symbol')
  #filters = 'external_gene_name',
  #values = rat_genes)
)
## ms -hs 

mm_hs <- getLDS(
  mart = mm,
  attributes = c('ensembl_gene_id',"mgi_symbol",'external_gene_name'),
  martL = hs,
  attributesL = c('ensembl_gene_id','hgnc_symbol'))

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

rat_genes= contrast_df %>%
  filter(grepl("Rn", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()

mouse_genes= contrast_df %>%
  filter(grepl("Mm", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()

human_genes= contrast_df %>%
  filter(grepl("Hs", contrast_id))%>% pull(gene) %>% unique()



BiocManager::install(c("vitkl/orthologsBioMART"), dependencies=T)
library(orthologsBioMART)
findOrthologsHsMm(from_filters = "hgnc_symbol",
                  from_values = c("TP53","TERT"), 
                  to_attributes = "external_gene_name")
mm.genes= gene.tbl %>% filter(spec == "Mm")%>% pull(gene)

hs.genes= findOrthologsMmHs(from_filters = "ensembl_gene_id",
                  from_values = c("Postn"),
                  to_attributes = "hgnc_symbol")


chicken <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
chimpz<- useMart('ensembl', dataset = 'ptroglodytes_gene_ensembl')

nnot_table <- getLDS(
  mart = chicken,
  attributes = c('ensembl_gene_id','external_gene_name','chromosome_name'),
  martL = chimpz,
  attributesL = c('ensembl_gene_id','external_gene_name','chromosome_name','gene_biotype'))



# JAX fix -----------------------------------------------------------------


contrast_df = readRDS("data/contrasts_query_df.rds")

gene.tbl = contrast_df %>% 
  ungroup()%>%
  mutate(gene = ifelse(grepl("Rn|Mm",gene), str_to_title(gene), gene),
         spec= substr(contrast_id,1,2))%>% 
  #select(-contrast_id)%>%
  distinct(gene, spec)

gene.tbl%>% write_csv("data/gene_names.csv")

rat_genes= contrast_df %>%
  filter(grepl("Rn", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()

mouse_genes= contrast_df %>%
  filter(grepl("Mm", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

# we want to map mouse to human and rat to human, then display only the translated genes
#we have 1 to 1 maps, no action
# 1 to many, 
# many to 1 ,
# 1 to none,
# quantify those categories for our DF
mouse_human_genes <- as_tibble(mouse_human_genes)
class_keys= mouse_human_genes %>%
  filter(Symbol %in% mouse_genes)%>% 
  pull(DB.Class.Key)

mappings.df= as_tibble(mouse_human_genes) %>% 
  group_by(DB.Class.Key)%>%
  filter(DB.Class.Key %in% class_keys)%>%
  distinct(DB.Class.Key, 
           Common.Organism.Name, 
           Symbol)%>%
  summarise(n= n())


mappings.df%>%
  count(n)%>% 
  ggplot(.,aes(x= n, y= nn))+
  geom_col()

big_map= mappings.df%>% filter(n== 17)%>% pull(DB.Class.Key)
mouse_human_genes%>% filter(DB.Class.Key %in% big_map)

one_to_one= mappings.df%>% filter(n== 2)%>%
  pull(DB.Class.Key)%>% unique()
one_to_none= mappings.df%>% filter(n== 1)%>%
  pull(DB.Class.Key)%>% unique()
ambiguous = mappings.df%>% filter(n>2)%>%
  pull(DB.Class.Key)%>% unique()

length(mouse_genes)
length(class_keys)
lapply(list(one_to_none, one_to_one, ambiguous), length)#%>% unlist() %>% sum()

# add the one_to_one mappings now to contrast_df

translate_mh <- mouse_human_genes%>%
  filter(DB.Class.Key %in% one_to_one)%>%
  group_by(DB.Class.Key)%>%
  mutate(Common.Organism.Name= ifelse(grepl("mouse", Common.Organism.Name), "mouse", Common.Organism.Name))%>%
  distinct(Common.Organism.Name, Symbol)%>%
pivot_wider(names_from = Common.Organism.Name, 
              values_from= Symbol, id_expand = F)%>%
  distinct(mouse, human)%>% ungroup()

contrast_df2= contrast_df2%>% 
  left_join(translate_mh%>%
              distinct(mouse, human)%>%
              rename(gene_orig= mouse))%>%
  mutate(gene= human)%>%
  dplyr::select(-human)


library(msigdbr)
genes= rat_genes
translate_species_to_hs = function(species= "Rattus norvegicus",
                                   genes){
  require(msigdbr)
  
  m_df = msigdbr(species = "Rattus norvegicus")
  
  m_df_2=m_df %>% filter(gene_symbol %in% genes)%>%
    distinct(gene_symbol,human_gene_symbol)
  
  one_to_none <- genes[!genes %in% m_df_2$gene_symbol]
  length(one_to_none)
  
  m_df_count= m_df_2 %>% 
    group_by(gene_symbol)%>% count()
  
  one_to_one = m_df_count %>% filter(n==1)%>% pull(gene_symbol)
  length(one_to_one)
  
  ambiguous = m_df_count %>% filter(n!=1)%>% pull(gene_symbol)
  length(ambiguous)
  
  print( paste("one to one:", length(one_to_one)))
  print( paste("one to none:", length(one_to_none)))
  print( paste("ambiguous maps:", length(ambiguous)))
  
  return(list("df"= m_df_2,
              "genes"= list("oto"= one_to_one,
                            "otn"=one_to_none,
                            "amb"= ambiguous)
              )
         )
  }


m_df = msigdbr(species = "Rattus norvegicus")
rat_genes
m_df_2=m_df %>% filter(gene_symbol %in% rat_genes)%>%
  distinct(gene_symbol,human_gene_symbol)

length(rat_genes)
length(unique(m_df_2$gene_symbol))
length(unique(m_df_2$human_gene_symbol))
m_df_sub = m_df[,c(5,8)] # select only the gene ids necessary
m_df_sub_dis = distinct(m_df_sub) # remove non unqique rows

df = merge(x = df, y = m_df_sub_dis, 'gene_symbol', all.x = TRUE)

m_df = msigdbr(species = "Mus musculus")
m_df_2=m_df %>% filter(gene_symbol %in% mouse_genes)%>%
  distinct(gene_symbol,human_gene_symbol)
length(unique(mouse_genes))
