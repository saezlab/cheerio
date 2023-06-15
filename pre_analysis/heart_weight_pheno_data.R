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



obj.list= readRDS("data/GEX.list.hypertrophy.rds")

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
single.sample.exp
unique(single.sample.exp$sid)

pheno.data= read_csv("raw_data/pheno/Updated batches_TAC_Swim_revision sorted_mod.csv")

pheno.data= pheno.data %>% select(ID, 'HW/BW')%>%
  rename(HW_BW='HW/BW' )
pheno.data2= pheno.data %>%
  rowwise()%>%
  mutate(id_clean= substr(ID, nchar(ID)-2, nchar(ID)))


HW_DF= single.sample.exp %>% left_join(pheno.data2, by= "id_clean")%>%
  mutate(exp= as.numeric(exp))

genes= c("Mylk4", "Nppa")
library(ggpmisc)
library(broom)
map(genes , function(x){
  #map(c("rna", "ribo"), function(y){
  HW_DF %>%
    filter(MgiSymbol == x
           )%>%
    #plot_dataframe(., "HW_BW", "exp", "modal", "model")
    ggplot(., aes(x= HW_BW, y= exp, color= model))+
    geom_point(aes(shape= tp), size = 4)+
    facet_grid(rows= vars(modal), scales="free_y")+
    stat_smooth(method = "lm", formula = my.formula, se = T) +
    stat_poly_eq(aes(label = paste(..rr.label..)), 
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
    ggtitle(x)
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


x= "Nppa"
map(genes, function(x){
  HW_DF2= HW_DF %>% filter(MgiSymbol==x)
  lm(HW_BW~ )
})
  



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

