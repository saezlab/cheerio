

library(tidyverse)

hcm <- read_csv("raw_data/sc_chaffin/HCMvsNF.csv")
colnames(hcm)<- str_replace_all(colnames(hcm), '\\n', "")
hcm <- hcm %>%
  mutate(Comparison ="HCMvsNF")%>%
  rename(logFC = "CellBender:logFC", 
         pval=`CellRanger:P-Value`,       
         FDR= "CellBender:Adjusted P-Value",
         Cell_type = "Cell Type")%>%
  filter(BackgroundContaminationFlag==0)%>%
  select(Comparison, 
         Gene, 
         Cell_type,
         logFC, 
         pval, 
         FDR )


dcm <- read_csv("raw_data/sc_chaffin/DCMvsHCM.csv")
colnames(dcm)<- str_replace_all(colnames(dcm), '\\n', "")
dcm <- dcm %>%
  rename(logFC = "CellBender:logFC", 
         pval=`CellRanger:P-Value`,       
         FDR= "CellBender:Adjusted P-Value",
         Cell_type = "Cell Type")%>%
  mutate(Comparison ="HCMvsDCM")%>%
  filter(BackgroundContaminationFlag==0)%>%
  select(Comparison, 
         Gene, 
         Cell_type,
         logFC, 
         pval, 
         FDR, )%>%
  mutate(logFC= -1 * logFC)

range(hcm$pval, na.rm = T)

rbind(dcm, hcm)%>%
  write_csv("data/chaffinetal_degs.csv")


