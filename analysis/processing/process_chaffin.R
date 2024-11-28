

library(tidyverse)

source_data <- read_csv("../cheerio_data/raw_data/sc_chaffin/41586_2022_4817_MOESM5_ESM.csv")

sc.gex2= source_data %>% 
  mutate(pval= 10^(-`log10(P)`), 
         logFC = ifelse(grepl("DCMvsHCM", Comparison), logFC*-1, logFC))   %>% 
  group_by(Comparison, CellType)%>%
  mutate(FDR= p.adjust(pval))%>%
  rename(Cell_type= CellType)%>%
  select("Comparison" ,"Gene"     ,  "Cell_type"  ,"logFC"   ,  "pval"  ,     "FDR"    )
hist(sc.gex2$FDR)
hist(sc.gex2$pval)

#flip direction of contrast comp (logFC flipped above)
sc.gex2$Comparison= str_replace_all(pattern = "DCMvsHCM",replacement =  "HCMvsDCM",sc.gex2$Comparison)
sc.gex2%>%
  write_csv("data/chaffinetal_degs.csv")


# 
# 
# hcm <- read_csv("../cheerio_data/raw_data/sc_chaffin/HCMvsNF.csv")
# colnames(hcm)<- str_replace_all(colnames(hcm), '\\n', "")
# 
# 
# hcm <- hcm %>%
#   mutate(Comparison ="HCMvsNF")%>%
#   rename(logFC = "CellBender:logFC", 
#          pval=`CellRanger:P-Value`,       
#          FDR= "CellBender:Adjusted P-Value",
#          Cell_type = "Cell Type")%>%
#   filter(BackgroundContaminationFlag==0)%>%
#   select(Comparison, 
#          Gene, 
#          Cell_type,
#          logFC, 
#          pval, 
#          FDR )
# 
# dcm <- read_csv("../cheerio_data/raw_data/sc_chaffin/DCMvsHCM.csv")
# colnames(dcm)<- str_replace_all(colnames(dcm), '\\n', "")
# dcm <- dcm %>%
#   rename(logFC = "CellBender:logFC", 
#          pval=`CellRanger:P-Value`,       
#          FDR= "CellBender:Adjusted P-Value",
#          Cell_type = "Cell Type")%>%
#   mutate(Comparison ="HCMvsDCM")%>%
#   filter(BackgroundContaminationFlag==0)%>%
#   select(Comparison, 
#          Gene, 
#          Cell_type,
#          logFC, 
#          pval, 
#          FDR, )%>%
#   mutate(logFC= -1 * logFC)
# 
# range(hcm$pval, na.rm = T)
# 
# df <- rbind(dcm, hcm)
# df%>%
#   write_csv("data/chaffinetal_degs.csv")
# 
# 
# hcm%>% filter(Cell_type =="Endothelial I")%>% pull()
# unique(hcm$Cell_type)
# 
# df %>% 
#   ggplot(aes(x= Cell_type, y= logFC, color=FDR<0.05))+
#   geom_jitter(width= 0.3)+
#   facet_grid(~Comparison)
# 
# df%>% filter(grepl("COL", Gene))
# # 