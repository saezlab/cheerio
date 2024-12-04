## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-05-26
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## 
# ---------------------------------------------------------------------------------------------



library(tidyverse)
library(biomaRt)
#library(limma)
#library(edgeR)
#library(ggrepel)

pw<- "~/R-projects/Collaborations/cheerio/"
gc= read.csv(paste0(pw, "../cheerio_data/raw_data/magnet/gene_count_matrix.csv"))
#tc= read.csv("data/raw_data/magnet/transcript_count_matrix.csv")
meta= read.csv(paste0(pw,"../cheerio_data/raw_data/magnet/MAGE_metadata.txt"))

meta$Sample.Name
colnames(gc)
dim(gc)
rownames(gc)
head(gc)
gc$gene_id

#translate gene ID to symbol
hs <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#listFilters(ensembl)
foo <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name'),
                          #'mgi_symbol')
             filters = 'ensembl_gene_id',
             values = gc$gene_id,
             mart = hs)



gex.df= as_tibble(gc)%>% 
  left_join(foo %>% rename(gene_id = ensembl_gene_id,
                           gene= external_gene_name))#%>% 
summed.gex= gex.df %>%
  dplyr::select(gene, everything(), -gene_id)%>%
  filter(gene != "")%>%
  group_by(gene)%>%
  summarise_all(sum)
  
df= column_to_rownames(summed.gex, "gene")

colnames(df)= substr(colnames(df), 1,11)

table(colnames(df) %in% meta$Run)

meta= meta %>% arrange(Sample.Name)

colnames(df)[match( colnames(df), meta$Run,)]

meta= meta %>% filter(Run %in% colnames(df), 
                      etiology %in% c("Non-Failing Donor",
                                      "Hypertrophic cardiomyopathy (HCM)" ,
                                      "Dilated cardiomyopathy (DCM)")
                      )%>%
  as_tibble()

dim(df)
df= df[,meta$Run]
dim(df)

unique(meta$etiology)
dge<- DGEList(counts=df, group= meta$etiology)#, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.count	= 15, min.total.count= 20, min.prop = 0.75)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)



# OCA -----------------------------------------------------------------------------------------


PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = F)

plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Run") %>%
  as_tibble()%>%
  left_join(meta, by= "Run" )

p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2, col= etiology))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))
  geom_text_repel(aes(label= sample),show.legend = FALSE)

p.pca



#Adjust a linear model to the expression data
group= meta$group
meta  = meta %>% mutate(group= ifelse(grepl("Non", etiology), "CT", 
                                      ifelse(grepl("DCM", etiology), "DCM", 
                                             "HCM")))

f = factor(meta$group, levels= c("CT", "DCM", "HCM"))
f2= factor(meta$race)
f3= (meta$Age)
f4= meta$sex
design = model.matrix(~0+f+f2+f3+f4)
#colnames(design) = c("Ct","HFpEF")

fit = lmFit(v$E, design)

#Define contrasts
cont.matrix = makeContrasts(HCMvsNF = fHCM-fCT,
                            DCMvsNF = fDCM-fCT,
                            HCMvsDCM = fHCM-fDCM,
                            levels=design)

#Empirical Bayes
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

tble.list= 
map(colnames(fit2$coefficients), function(x){
  DE_results = as.data.frame(topTable(fit2,coef = x,
                                      adjust.method = "BH",
                                      number = Inf)) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(t))) %>%
    as_tibble()
})
names(tble.list)= colnames(fit2$coefficients)

map(tble.list, function(x){
  x%>%mutate(sig= ifelse(adj.P.Val<0.05,T, F))%>%
    ggplot(.,aes(x= logFC, y= -log10(P.Value), col =sig))+
    geom_point()
    
})

map(tble.list, function(x){
  x%>%mutate(sig= ifelse(adj.P.Val<0.05,T, F))%>%
    ggplot(.,aes(x= (P.Value)))+
    geom_histogram()
  
})

saveRDS(tble.list, "raw_data/magnet/contrast_list.rds")


tble.list <- readRDS("../cheerio_data/raw_data/magnet/contrast_list.rds")

x <-readRDS("~/Downloads/CPMS_SVA_corrected.RDS")

boxplot(log10(x[, 1:10]))



unique(meta$etiology)
dge<- DGEList(counts=log2(x+1), group= meta$etiology)#, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.count	= 15, min.total.count= 20, min.prop = 0.75)
table(keep)

dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

