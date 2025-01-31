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
## prepare phenotype data

library(tidyverse)
library(ggpmisc)
library(broom)
library(ggrepel)

obj.list= readRDS("data/GEX.list.hypertrophy.rds")

source("analysis/utils_pre_analysis.R")

##extract contrast data into single tidy format
single.sample.exp= lapply(obj.list, function(mod){
  
  lapply(mod, function(dat){
    #print(colnames(dat))
    
    lapply(dat, function(tp){
      print(colnames(tp))
      tp %>% select(contains("CPM"),MgiSymbol, tp, modal, model)%>% select(-logCPM)%>% 
        pivot_longer(!c(MgiSymbol, tp, modal, model),names_to = "sid", values_to= "exp",)
    })%>%
      do.call(rbind,.)
  })%>%
    do.call(rbind,.)
})%>%
  do.call(rbind,.)

single.sample.exp = single.sample.exp %>% 
  mutate(sid= str_replace_all(sid, "_CPM", ""))%>%
  mutate(id_clean= substr(sid, nchar(sid)-2, nchar(sid)))

unique(single.sample.exp$sid)

pheno.data= read_csv("raw_data/pheno/Updated batches_TAC_Swim_revision sorted_mod.csv")

pheno.data= pheno.data %>% select(ID, 'HW/BW')%>%
  rename(HW_BW='HW/BW' )
pheno.data2= pheno.data %>%
  rowwise()%>%
  mutate(id_clean= substr(ID, nchar(ID)-2, nchar(ID)))

HW_DF= single.sample.exp %>% left_join(pheno.data2, by= "id_clean")%>%
  mutate(exp= as.numeric(exp),
         exp.group= ifelse(grepl("sham|sedent", sid), "ct", "exp"),
         exp.group= paste0(exp.group, "_", tp),
         modal= factor(paste0(toupper(modal), "seq"), levels= c("RNAseq", "RIBOseq"))
         )

HW_DF%>%
  ggplot(aes(y= log10(exp), x=sid))+
  geom_boxplot()+
  coord_flip()+
  facet_grid(~model+ modal)

mm_hs= translate_species_to_hs("Mus musculus", unique(HW_DF$MgiSymbol))

HW_DF2= HW_DF%>%
  filter(MgiSymbol %in% mm_hs$genes$oto)%>%
  rename(gene_orig = MgiSymbol)%>%
  #dplyr::select(-gene)%>%
  left_join(mm_hs$df%>%
              distinct(gene_symbol, human_gene_symbol)%>%
              rename(gene_orig= gene_symbol,
                     gene= human_gene_symbol))


saveRDS(HW_DF2, "data/heart_weight_gex.rds")
HW_DF2<- readRDS("data/heart_weight_gex.rds")

# pre-calculate for each group --------------------------------------------


results <- HW_DF2 %>% 
  filter(modal =="RNAseq")%>% # we focus only on RNA for the korrelation
  mutate(logcpm = log10(exp)) %>%
  group_by(gene_orig, tp, model) %>%
  do({
    fit <- lm(HW_BW ~ logcpm, data = .)
    summary_fit <- summary(fit)
    
    # Extracting R2, coefficient for logcpm, and its p-value
    tibble(
      gene_orig = unique(.$gene_orig),
      tp = unique(.$tp),
      model = unique(.$model),
      R2 = summary_fit$r.squared,
      logcpm_coef = coef(fit)["logcpm"],
      logcpm_p_value = summary_fit$coefficients["logcpm", "Pr(>|t|)"]
    )
  }) %>%
  ungroup()

results <- results %>% 
  group_by(model, tp)%>%
  mutate(FDR = p.adjust(logcpm_p_value, method="BH"))%>%
  ungroup()

# use performed translation: 
c.df<- readRDS("data/contrasts_query_df_translated3.rds")

x<- c.df %>% 
  filter(grepl("mm", contrast_id))%>%
  filter(gene_orig %in% results$gene_orig)%>%
  distinct(gene, gene_orig)

results<-results %>% left_join(x, by= "gene_orig")

results %>% saveRDS("app_data/heart_weight_precalc.rds")

genes= c("Mylk4", "Nppa")
x= "Nppa"

results %>%
  filter(gene_orig %in% genes)%>%
  mutate(group= paste(model, tp,  sep = "_"))%>%
  ggplot(aes(x= group, y= logcpm_coef, 
             #size= -log10(logcpm_p_value),
             color= R2, shape= logcpm_p_value>.05))+
  geom_point(aes(shape= logcpm_p_value>.05), show.legend = T,
             size= 3)+
  facet_grid(~gene_orig)+
  theme(axis.text.x= element_text(angle= 90, hjust= 1, vjust = 0.5),
        panel.grid.major = element_line(color = "gray80", size = 0.5), # Major grid lines
        panel.grid.minor = element_line(color = "gray90", size = 0.25)  # Minor grid lines
        )+
  geom_hline(yintercept = 0)+
  scale_shape_manual(values = rev(c("square", "circle")))+
  scale_size_continuous(range = c(1, 8))+
  scale_color_gradient(low= "grey", high= "red")+
  labs(x= "", y= "Coefficient")
  



# fit model on demand -----------------------------------------------------


genes= c("Mylk4", "Nppa")
x= "Nppa"
my.formula <- y~x
map(genes , function(x){
  #map(c("rna", "ribo"), function(y){
  HW_DF %>%
    filter(gene_orig == x )%>%
    ggplot(., aes(x= HW_BW, y= exp, color= model))+
    geom_point(aes(shape= exp.group), size = 3, alpha= 0.6)+
    facet_grid(rows= vars(modal), scales="free_y")+
    stat_smooth(fullrange = T, method = "lm", formula = my.formula, se = F, linewidth= 0.4) +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))), 
                label.x = "left", label.y = "top",
               formula = my.formula, parse = TRUE, size = 4)+
     stat_fit_glance(method = 'lm',
                     method.args = list(formula =my.formula),
                    geom = 'label_repel', 
                     aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                     label.x = 'right', label.y = "top", size = 4, alpha= 0.7)+
    # stat_fit_glance(method = 'lm',
    #                 method.args = list(formula = formula),
    #                 geom = 'text',
    #                 aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
    #                 label.x.npc = 'right', label.y.npc = 0.35, size = 3)+
    #geom_smooth(method = "lm")+
    ggtitle(x)+
    scale_color_manual(values = c("swim" = "darkblue",
                                 "tac"="darkred", 
                                 drop= FALSE))+
    labs(y= "Normalized gene expression",
         x= "Normalized heart weight", 
         shape= "Experimental group\n(treatment_timepoint)")+
    theme(panel.grid.major = element_line(color = "grey",
                                          linewidth = 0.1,
                                          linetype = 1),
          panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) 
    })


library(ggrepel)
my.formula <- y~x
library(scales)

plot_dataframe <- function(data, facet_variable, color_variable) {
  ggplot(datas, aes(x = HW_BW, y = exp, color = model)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_grid(rows = vars(modal)) +
    #xlab(x_column) +
    #ylab(y_column) +
    annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = 0,
             label = paste0("R^2 = ", percent(summary(lm(exp ~ HW_BW, data = datas))$r.squared)),
             parse = TRUE) +
    annotate("text", x = Inf, y = -Inf, hjust = 1, vjust = -1,
             label = paste0("p-value = ", signif(summary(lm(exp ~ HW_BW, data = datas))$coefficients[2, "Pr(>|t|)"], digits = 3)))
}
    
summary(lm(exp ~ HW_BW, data = datas))$r.squared





str_replace_all(pheno.data$ID, pattern = "-", replacement = "")


s.s[grepl(pattern = "CPM", s.s)]


s.s = colnames(gex.obj$TAC$RNA$'2d')
x= lapply(c(colnames(gex.obj$TAC$RNA$'2d'),colnames(gex.obj$TAC$RNA$'2wk'),colnames(gex.obj$Swim$RNA$'2d'),colnames(gex.obj$Swim$RNA$'2wk')),
       function(x){
      y=  x[grepl(pattern = "CPM", x)]
      y= y[!grepl(pattern = "log", x)]
      str_replace_all(y, "_CPM", "")
       })

x= unlist(x)
sample.ids= x[!is.na(x)]
sample.ids
s.id2= substr(sample.ids, nchar(sample.ids)-2, nchar(sample.ids))


table(s.id2 %in% pheno.data2$id_clean)
table(pheno.data2$id_clean %in% s.id2)

map(pheno.data$ID, function(x){
  
  substr(x, length(x)-4, length(x))
  
})
y= do.call(rbind, gex.obj)

