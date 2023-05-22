

library(decoupleR)
library(progeny)
library(dorothea)

data("dorothea_hs")
data("dorothea_mm")

contrasts = readRDS("data/contrasts.hypertrophy.rds")


# prepare contrast data in matrix format ------------------------------------------------------

#mouse 

x= contrasts %>%
  filter(model != "fetal")%>%
  mutate(exp.group= paste(tp, modal, model, sep= "_"))%>%
  select(logFC, exp.group, MgiSymbol)%>% 
  pivot_wider(names_from = exp.group, values_from= logFC,values_fn= mean) 

# we run rna and ribo separately, they have very different gene coverage 
rna.mat= x %>% select( MgiSymbol, grep("rna", colnames(x))) %>% drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()
ribo.mat=x %>% select( MgiSymbol, grep("ribo", colnames(x)))%>% drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()

#human- reheat


reheat_mat= contrasts_HF %>%
  select(study, t, gene)%>% 
  pivot_wider(names_from = study, values_from= t,values_fn= mean)

reheat_mat= reheat_mat %>% column_to_rownames("gene")%>%as.matrix()

#human- sc 
sc.gex%>% mutate(effect= logFC*log10.P.)
sc.df= sc.gex%>%
  mutate(logFC= logFC*log10.P.)%>% 
  select(logFC, Gene, Comparison, CellType)
sc.df.l= split(sc.df, sc.df$CellType) 
sc.df.m= lapply(sc.df.l, function(x) {
  x %>% 
    select(-CellType) %>% 
    pivot_wider(names_from = Comparison, values_from= logFC ,values_fn= mean)%>% 
    column_to_rownames("Gene")%>%
    as.matrix()
})

# human- fetal

fetal.m= contrasts %>%
  filter(model == "fetal")%>%
  select(logFC, tp, MgiSymbol)%>% 
  pivot_wider(names_from = tp, values_from= logFC,values_fn= mean) %>% 
  drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()




# run DoRoTHEA ------------------------------------------------------------

net= dorothea_mm %>% rename(source= tf) %>% filter(confidence %in% c("A", "B", "C"))

ulm.rna= decoupleR::run_ulm(mat= rna.mat, network = net ) 
ulm.ribo= decoupleR::run_ulm(mat= ribo.mat, network = net )

tf_hypertrophy= rbind(ulm.rna %>% mutate("modal"= "rna"),
                      ulm.ribo %>% mutate("modal"= "ribo"))

tf_hypertrophy= 
  tf_hypertrophy%>% mutate(model= ifelse(grepl("swim", condition), "swim", "tac"),
                         tp = ifelse(grepl("2d", condition), "2d", "2wk"),
                         sig= p_value<0.05)%>%
  select(-condition)

saveRDS(tf_hypertrophy, "data/dorothea_hypertrophyMM.rds")
tf_hypertrophy= readRDS("data/dorothea_hypertrophyMM.rds")


## reheat
net= dorothea_hs %>% rename(source= tf) %>% filter(confidence %in% c("A", "B", "C"))

reheat_mat
x.df= apply(reheat_mat, 2, FUN = function(x){
  
  x= as.matrix(x)
  rownames(x)= rownames(reheat_mat)
  x <- x[!rowSums(is.na(x)),]
  x= as.matrix(x)
  ulm.reheat= decoupleR::run_ulm(mat= x, network = net )
  
})
reheat.tfs= lapply(names(x.df), function(x){
  x.df[[x]] %>% mutate(study = x)
}) %>% do.call(rbind, .)

##single cell
sc.df.tf= lapply(sc.df.m, function(x){
  x <- x[!rowSums(is.na(x)),]
  decoupleR::run_ulm(mat= x, network = net )
})

sc.df.tf2= lapply(names(sc.df.tf), function(x){
  sc.df.tf[[x]] %>% mutate(celltype = x)
}) %>% do.call(rbind, .)

## fetal

rownames(fetal.m) = toupper(rownames(fetal.m))
x <- fetal.m[!rowSums(is.na(fetal.m)),]
fetal.tf= decoupleR::run_ulm(mat= x, network = net )



saveRDS(list("mm"= tf_hypertrophy, 
             "hs_reheat"= reheat.tfs%>% mutate(sig= ifelse(p_value<0.05, T, F)), 
             "hs_sc"= sc.df.tf2%>% mutate(sig= ifelse(p_value<0.05, T, F)), 
             "hs_fetal"= fetal.tf%>% mutate(sig= ifelse(p_value<0.05, T, F))),
             "data/dorothea_results_all.rds")
df_tf= readRDS("data/dorothea_results_all.rds")


tf_hypertrophy= tf_hypertrophy%>% pivot_wider(-p_value, names_from = "condition", 
                              values_from = "score")

ggplot(tf_hypertrophy, aes(x= `2wk_rna_tac`, y= `2wk_ribo_tac`))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot(tf_hypertrophy, aes(x= `2d_rna_tac`, y= `2d_ribo_tac`))+
  geom_point()+
  geom_smooth(method = "lm")


## plot: 

tfs= c("Ar", "Gata4")
pls= map(tfs, function(x){
  tf_hypertrophy %>% 
    filter(source ==x)%>%
    mutate(condition= paste(tp, modal,sep =  "_"))%>%
    ggplot(aes(x = condition, y =  score, fill = sig)) +
    facet_grid(rows= vars(model))+
    #geom_boxplot()+
    geom_hline(yintercept = 0, color= "black")+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot() +
    scale_fill_manual(values = c("TRUE" = "darkgreen",
                                 "FALSE"="orange"))+
    labs(x = "experimental group", y = "TF activity", fill = "p<0.05") +
    theme(panel.grid.major = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 1),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(x)
  
})
p1= cowplot::plot_grid(plotlist =  pls)

p1= tf_hypertrophy %>%
  filter(model!= "fetal")%>%
  select(-PValue, -FDR)%>%
  pivot_wider(names_from= modal, values_from = logFC, values_fn= mean)%>%
  mutate(labels= ifelse(MgiSymbol %in% input$select_gene, MgiSymbol, ""),
         labls= factor(labels, levels= c("", input$select_gene)),
         alphas= factor(ifelse(labels=="", "bg","normal"))
  )%>%
  ggplot(aes(x= rna, y= ribo, color = labels,size= alphas, alpha= alphas))+
  facet_grid(rows= vars(model), 
             cols= vars(tp))+
  geom_point(show.legend = T)+
  scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1))+
  scale_size_manual(values=c("bg"= 0.5, "normal"= 2))+
  #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
  theme(panel.grid.major = element_line(color = "grey",
                                        size = 0.1,
                                        linetype = 1))+
  labs(alpha= "")


pls= map(tfs, function(x){
  contrasts %>% 
    filter(MgiSymbol==x)%>%
    mutate(exp.group= paste(modal, tp, sep= "_"),
           sig= FDR<0.05)%>%
    ggplot(aes(x = exp.group, y = logFC, fill = sig)) +
    facet_grid(rows= vars(model))+
    #geom_boxplot()+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot() +
    scale_fill_manual(values = c("TRUE" = "darkgreen",
                                 "FALSE"="orange"))+
    labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size= 11), 
      axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(x)
  
})
pls  
cowplot::plot_grid(plotlist = pls)



# run progeny -------------------------------------------------------------------------------------

library(ComplexHeatmap)
#pathway:

##mouse=
p1.rna= progeny::progeny(rna.mat, organism = "Mouse", z_scores= T, perm= 500)
p2.ribo= progeny::progeny(ribo.mat, organism = "Mouse", z_scores= T, perm= 500)


rownames(p1.rna)= str_replace_all(rownames(p1.rna), "X", "")
rownames(p2.ribo)= str_replace_all(rownames(p2.ribo), "X", "")

## human

# apply(reheat_mat, 2, function(x){
#   print(as.matrix(x[1:10]))
#   progeny::progeny(as.matrix(x), organism = "Human", z_scores= T, perm= 500)
#   
# })
p3.reheat= progeny::progeny(reheat_mat, organism = "Human", z_scores= T, perm= 500)

##single cell
sc.df.prg= lapply(sc.df.m, function(x){
  progeny::progeny(x, organism = "Human", z_scores= T,top=500,  perm= 500)
})

sc.df.prg

# max(sapply(sc.df.prg, max))
# min(sapply(sc.df.prg, min))

sc.hmaps= lapply(names(sc.df.prg), function(x){
  plot_hmap(sc.df.prg[[x]], x,
            max.ps = c(min(sapply(sc.df.prg, min)),
                       max(sapply(sc.df.prg, max)))
            )
  })
names(sc.hmaps)= names(sc.df.prg)
#sc.hmaps= sc.hmaps[1:4]

eval(parse(text= paste(paste0("sc.hmaps$",paste0("`", names(sc.hmaps), "`")),  collapse = " + ")))
sc.hmaps$`Adipocyte`+ sc.hmaps$Cardiomyocyte+sc.hmaps$`Endothelial I`

## fetal

p4.fetal= progeny::progeny(fetal.m, organism = "Mouse", z_scores= T, perm= 500)

#plot full matrix

plot_hmap(p1.rna)
plot_hmap(p2.ribo)
plot_hmap(p3.reheat)
plot_hmap(p4.fetal)

saveRDS(list("mmRNA"= p1.rna,
              "mmRibo"=p2.ribo, 
              "hsReheat"=p3.reheat, 
              "hsfetal"= p4.fetal, 
              "hsSC"= sc.df.prg), "data/progeny_results_all.rds")

model.m =progeny::getModel("Mouse")
pway= "TGFb"

df= merge(model.m, rna.mat, by= "row.names" )%>% as_tibble()%>%
  rename(gene = Row.names)%>%
  pivot_longer(cols= all_of(colnames(model.m)), names_to = "pathway", values_to= "weight")

df%>% 
  filter(pathway== pway,
         weight !=0)%>% 
  ggplot(., aes(x= weight, y= `2d_rna_swim`))+
  geom_point()




