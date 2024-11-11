library(tidyverse)
prosser<- read.csv("raw_data/Proteom Analyse fu╠Иr Jan/Human Prosser Paper/Limma_results_V1.csv")

proteome_mouse <- read.csv("raw_data/Proteom Analyse fu╠Иr Jan/Maus Swim/Limma_results_V1.csv")
unique(proteome_mouse$comparison.label)

unique(prosser$comparison.label)
prosser %>% group_by(comparison.label)%>% table()
