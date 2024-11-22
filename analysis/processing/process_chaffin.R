

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

df <- rbind(dcm, hcm)
df%>%
  write_csv("data/chaffinetal_degs.csv")


hcm%>% filter(Cell_type =="Endothelial I")%>% pull()
unique(hcm$Cell_type)

df %>% 
  ggplot(aes(x= Cell_type, y= logFC, color=FDR<0.05))+
  geom_jitter(width= 0.3)+
  facet_grid(~Comparison)
