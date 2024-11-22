
library(decoupleR)
library(progeny)
library(dorothea)

data("dorothea_hs")
data("dorothea_mm")

contrasts = readRDS("data/contrasts.hypertrophy.rds")
joint_contrast_df= readRDS("data/contrasts_query_df.rds")

# prepare contrast data in matrix format ------------------------------------------------------

#mouse 

x= contrasts %>%
  filter(model != "fetal")%>%
  mutate(exp.group= paste(tp, modal, model, sep= "_"))%>%
  dplyr::select(logFC, exp.group, MgiSymbol)%>% 
  pivot_wider(names_from = exp.group, values_from= logFC,values_fn= mean) 

# we run rna and ribo separately, they have very different gene coverage 
rna.mat= x %>% dplyr::select( MgiSymbol, grep("rna", colnames(x))) %>% drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()
ribo.mat=x %>% dplyr::select( MgiSymbol, grep("ribo", colnames(x)))%>% drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()

#human- reheat


reheat_mat= contrasts_HF %>%
  dplyr::select(study, t, gene)%>% 
  pivot_wider(names_from = study, values_from= t,values_fn= mean)

reheat_mat= reheat_mat %>% column_to_rownames("gene")%>%as.matrix()

#human- sc 

sc.df= sc.gex%>%
  mutate(logFC= logFC*log10.P.)%>% 
  dplyr::select(logFC, Gene, Comparison, CellType)
sc.df.l= split(sc.df, sc.df$CellType) 
sc.df.m= lapply(sc.df.l, function(x) {
  x %>% 
    dplyr::select(-CellType) %>% 
    pivot_wider(names_from = Comparison, values_from= logFC ,values_fn= mean)%>% 
    column_to_rownames("Gene")%>%
    as.matrix()
})

#human - magnet

magnet.df= joint_contrast_df%>% 
  mutate(effect= logFC*-log10(FDR))%>%
  filter(grepl("bulk", contrast_id))%>%
  dplyr::select(gene, contrast_id, effect)%>%
  pivot_wider(names_from = contrast_id, values_from= effect ,values_fn= mean)%>% 
  column_to_rownames("gene")

# rat in vitro: 
## deal with infinity

rat.df= joint_contrast_df%>% 
  filter(grepl("invitro", contrast_id))%>%
  mutate(FDR_min= min(FDR[FDR>0]), 
    effect= logFC*-log10(FDR_min))%>%
  dplyr::select(gene, contrast_id, effect)%>%
  pivot_wider(names_from = contrast_id, values_from= effect ,values_fn= mean)%>% 
  filter(!is.na(gene))%>%drop_na()%>%mutate(gene= str_to_title(gene))%>%
  column_to_rownames("gene")
  

rat.df  
#sc 


# human- fetal

fetal.m= contrasts %>%
  filter(model == "fetal")%>%
  dplyr::select(logFC, tp, MgiSymbol)%>% 
  pivot_wider(names_from = tp, values_from= logFC,values_fn= mean) %>% 
  drop_na()%>% column_to_rownames("MgiSymbol")%>%as.matrix()




# run DoRoTHEA ------------------------------------------------------------

net= dorothea_mm %>% rename(source= tf) %>% filter(confidence %in% c("A", "B", "C"))

ulm.rna= decoupleR::run_ulm(mat= rna.mat, network = net ) 
ulm.ribo= decoupleR::run_ulm(mat= ribo.mat, network = net )
ulm.rt= decoupleR::run_ulm(mat= as.matrix(rat.df), network = net )

tf_hypertrophy= rbind(ulm.rna %>% mutate("modal"= "rna"),
                      ulm.ribo %>% mutate("modal"= "ribo"), 
                      ulm.rt%>% mutate("modal"= "ribo"))

tf_hypertrophy= 
  tf_hypertrophy%>%
  mutate(model= ifelse(grepl("swim", condition), "swim",
                       ifelse(grepl("tac", condition), "tac", "invitro")),
         modality= ifelse(grepl("rna", condition), "rna", "ribo"),
         tp = ifelse(grepl("2d", condition), "2d", "2wk"),
         tp = ifelse(grepl("invitro", condition), "", tp),
         sig= p_value<0.05)%>%
  dplyr::select(-condition)
View(tf_hypertrophy)
tf_hypertrophy%>% arrange(source)

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

##magnet 
magnet.tf = decoupleR::run_ulm(mat= magnet.df, network = net ) %>% 
  mutate(sig= ifelse(p_value<0.05, T, F))



## fetal

rownames(fetal.m) = toupper(rownames(fetal.m))
x <- fetal.m[!rowSums(is.na(fetal.m)),]
fetal.tf= decoupleR::run_ulm(mat= x, network = net )



saveRDS(list("mm"= tf_hypertrophy, 
             "hs_reheat"= reheat.tfs%>% mutate(sig= ifelse(p_value<0.05, T, F)), 
             "hs_sc"= sc.df.tf2%>% mutate(sig= ifelse(p_value<0.05, T, F)), 
             "hs_magnet"= magnet.tf, 
             "hs_fetal"= fetal.tf%>% mutate(sig= ifelse(p_value<0.05, T, F))),
             "data/dorothea_results_all.rds")

df_tf= readRDS("data/dorothea_results_all.rds")

# df_tf$hs_magnet
# 
# tf_hypertrophy= tf_hypertrophy%>% pivot_wider(-p_value, names_from = "condition", 
#                               values_from = "score")
# 
# ggplot(tf_hypertrophy, aes(x= `2wk_rna_tac`, y= `2wk_ribo_tac`))+
#   geom_point()+
#   geom_smooth(method = "lm")
# 
# ggplot(tf_hypertrophy, aes(x= `2d_rna_tac`, y= `2d_ribo_tac`))+
#   geom_point()+
#   geom_smooth(method = "lm")
# 

## plot: 


# run progeny -------------------------------------------------------------------------------------

library(ComplexHeatmap)
#pathway:

##mouse=
p1.rna= progeny::progeny(rna.mat, organism = "Mouse", z_scores= T, perm= 500)
p2.ribo= progeny::progeny(ribo.mat, organism = "Mouse", z_scores= T, perm= 500)

p3.rat= progeny::progeny(as.matrix(rat.df), organism = "Mouse", z_scores= T, perm= 500)

rownames(p1.rna)= str_replace_all(rownames(p1.rna), "X", "")
rownames(p2.ribo)= str_replace_all(rownames(p2.ribo), "X", "")
#rownames(p1.rna)= str_replace_all(rownames(p1.rna), "X", "")

## human

# apply(reheat_mat, 2, function(x){
#   print(as.matrix(x[1:10]))
#   progeny::progeny(as.matrix(x), organism = "Human", z_scores= T, perm= 500)
#   
# })
p4.reheat= progeny::progeny(reheat_mat, organism = "Human", z_scores= T, perm= 500)

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
#sc.hmaps$`Adipocyte`+ sc.hmaps$Cardiomyocyte+sc.hmaps$`Endothelial I`


##magenet 
p5.magnet = progeny::progeny(as.matrix(magnet.df), organism = "Human", z_scores= T, perm= 500)

## fetal
rownames(fetal.m)= toupper(rownames(fetal.m))
p6.fetal= progeny::progeny(fetal.m, organism = "Human", z_scores= T, perm= 500)

#plot full matrix

plot_hmap(p1.rna)
plot_hmap(p2.ribo)
plot_hmap(p3.rat)
plot_hmap(p4.reheat)
plot_hmap(p6.fetal)
plot_hmap(p5.magnet)

saveRDS(list("mmRNA"= p1.rna,
              "mmRibo"=p2.ribo, 
             "rn"= p3.rat, 
              "hsReheat"=p4.reheat, 
              "hsfetal"= p6.fetal, 
             "hsmagnet"= p5.magnet, 
              "hsSC"= sc.df.prg), "data/progeny_results_all.rds")

