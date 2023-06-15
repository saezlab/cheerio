## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-06-15
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## read and process IPMC phenotype data



library(tidyverse)


dat= read_tsv("raw_data/pheno/IMPC/IMPC_Cardiovascular_System.tsv")
dat

unique(dat$Phenotype)
View(dat)
unique(dat$Parameter)


impc_data = dat%>% 
  mutate(Cardiac_Hypertrophy= ifelse(Phenotype %in% c("enlarged heart", "small heart", "decreased heart weight", "increased heart weight"), Phenotype, "-"), 
         Other_Phenotypes= ifelse(Cardiac_Hypertrophy == "-", Phenotype, "-"))%>%
  mutate(IPMC_link = paste0("<a href='", Data,"' target='_blank'>", "link","</a>"),
       Allele= str_replace_all(Allele, "<", "<sup"), 
       Allele= str_replace_all(Allele, ">", "</sup>"), 
       Allele= str_replace_all(Allele, "<sup", "<sup>"))

impc_data %>% saveRDS("data/ipmc_data.rds")


impc_data %>% distinct(Gene, Phenotype)
