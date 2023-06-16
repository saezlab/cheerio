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
  mutate(Cardiac_Hypertrophy= factor(ifelse(Phenotype %in% c("decreased heart weight", "increased heart weight"), Phenotype, "-")), 
         Other_Phenotypes= factor(ifelse(Cardiac_Hypertrophy == "-", Phenotype, "-")))%>%
  rename(P_val= `P Value`)%>%
  mutate(IPMC_link= ifelse(P_val!= 0,  modified_string <- str_extract(Data, ".*(?=&[^&]*$)"), Data),  
         # the links lead to the phenotyping center, we can remove the last part to link to thedata overview. the pval=0 condition is since exactly 0 p-vals carry different links
         IPMC_link = paste0("<a href='", IPMC_link,"' target='_blank'>", "link","</a>"), #add linkable html
          Allele= str_replace_all(Allele, "<", "<sup"), #make the allele superscript
          Allele= str_replace_all(Allele, ">", "</sup>"), 
          Allele= str_replace_all(Allele, "<sup", "<sup>"))

## to reduce amount of info , we take the median p-value if multiple exist dependning on sex or zygositiy or phenotype center
ipmc_data2= impc_data %>% 
  select(Gene, Allele, Cardiac_Hypertrophy, P_val,Other_Phenotypes,IPMC_link)%>% 
  group_by(Gene, Allele, Cardiac_Hypertrophy, Other_Phenotypes, IPMC_link)%>%
  summarise(median_p_val = median(P_val))%>%
  mutate(median_p_val= scientific(median_p_val))


ipmc_data2%>% filter(Gene== "Mybpc3")

ipmc_data2 %>% 
  ungroup()%>%
  saveRDS("data/ipmc_data.rds")


impc_data %>% distinct(Gene, Phenotype)
