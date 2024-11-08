library(BiocManager)
options(repos = BiocManager::repositories())
library(rclipboard)
library(shiny)
library(shinyWidgets)
library(shinyhelper)
library(shinyjs)
library(DT)
# remotes::install_github("christianholland/AachenColorPalette")
library(forcats)
library(scales)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(tidyr)
library(fgsea)
#library(decoupleR)
library(readr)
library(stringr)
library(shinycssloaders)
#plotting: 
library(RColorBrewer)
library(eulerr)
library(plotly)
library(AachenColorPalette) 
library(cowplot)
library(ggpmisc)
library(broom)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
#install.packages(c("shinytest", "showimage", "webdriver"))
theme_set(theme_cowplot())


# load static data

#gex.obj= readRDS("data/GEX.list.hypertrophy.rds")


##hypertrophy mice
#contrasts = readRDS("data/contrasts.hypertrophy.rds")
#tf_hypertrophy = readRDS("data/dorothea_hypertrophyMM.rds")


#sc chaffin dcm, hcm
# sc.gex= read.csv("data/sc_gex_chaffin.csv")%>%
#   as_tibble()%>%
#   mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)),
#          Comparison= str_replace(Comparison, "vs", "vs\n"),
#          Comparison = factor(Comparison, levels= c("DCMvs\nNF", "HCMvs\nNF", "DCMvs\nHCM")))

##reheat
contrasts_HF = readRDS("app_data/study_contrasts.rds")
ranks = readRDS("app_data/study_ranks.rds")
directed_signature = readRDS("app_data/signature.rds")

#example
example_geneset = read_csv("app_data/multiple_geneset.csv", show_col_types = FALSE)

#progeny:
prog.res= readRDS("app_data/progeny_results_all.rds")

#dorothea:
df_tf= readRDS("app_data/dorothea_results_all.rds")


#contrast query:
joint_contrast_df= readRDS( "app_data/contrasts_query_df.rds")


TFs= sapply(df_tf, function(x){
  unique(toupper(x$source))
})%>% unlist(use.names = F)%>% unique()

ipmc_data= readRDS("app_data/ipmc_data.rds")
HW_DF= readRDS("app_data/heart_weight_gex.rds")

error_text = "Queried gene(s) were not captured\nin data"
