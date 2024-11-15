
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggdendro)

c.df= readRDS("data/contrasts_query_df_translated3.rds")

contrasts_oi <- c.df %>%
  #filter(!grepl("DCMvsNF", contrast_id)) %>%
  filter(!grepl("HCMvsDCM", contrast_id)) %>%
  pull(contrast_id)%>% unique()


# plot histogram of contrasts and genes -----------------------------------

length(unique(c.df$gene))
contrast_comparisons<- c.df%>% 
  group_by(gene)%>%
  summarize(contrast_count = n_distinct(contrast_id))

sum(contrast_comparisons$contrast_count >9)

p.hist.counts <- contrast_comparisons%>% 
  ggplot(aes(x= contrast_count))+
  geom_histogram(bins = 30)+
  scale_y_log10()+
  labs(x= "number of contrasts", 
       y= "number of genes")+
  theme_cowplot()

p.hist.counts

# plot contrast overview --------------------------------------------------

#create a filtered version of the contrast data frame
c.df2<- c.df  %>%
  filter(contrast_id %in% contrasts_oi)

# Filter and count genes with FDR < 0.05 for each contrast
count_df <- c.df2 %>%
  filter(contrast_id %in% contrasts_oi) %>%
  filter(FDR < 0.05) %>%
  group_by(contrast_id) %>%
  summarise(num_genes = n())%>%
  complete(contrast_id = unique(contrasts_oi), fill = list(num_genes = 0))#%>%
  #mutate(contrast_id = factor(contrast_id, levels = ordered_contrast_ids)) # Set order based on clustering

# Bar plot displaying the number of genes with FDR < 0.05
bar_plot <- ggplot(count_df, aes(x = contrast_id, y = num_genes)) +
  theme_cowplot() +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, 
                                                   suffix = " × 10³"),
                     breaks = seq(0, max(count_df$num_genes), by = 4000)) +
  geom_hline(yintercept = seq(2000, max(count_df$num_genes), by = 2000), 
             color = "black", linetype = "dashed") +  # Horizontal lines
  geom_bar(stat = "identity", fill = "grey") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10)) +
  labs(y = "DEG\ncount")
bar_plot
# Main plot
main_plot <- c.df %>%
  filter(contrast_id %in% contrasts_oi) %>%
  mutate(color = factor(ifelse(FDR < 0.05 & logFC > 0, "up",
                               ifelse(FDR < 0.05 & logFC < 0, "dn", "unreg")),
                        levels = c("unreg", "up", "dn"))) %>%  # Set order of plotting
#  mutate(contrast_id = factor(contrast_id, levels = ordered_contrast_ids))%>% # Set order based on clustering
  ggplot(aes(x = contrast_id, y = logFC, color = color)) +
  geom_jitter(size = 0.2, width = 0.3) +
  theme_cowplot() +
  scale_color_manual(values = c("grey", "red", "blue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5),
        legend.position = "none") +
  labs(x = "")

# Display the combined plot



combined_plot <- plot_grid(
  #dendrogram_plot,
  bar_plot,
  main_plot,
  ncol = 1,
  align = "v",
  rel_heights = c(0.2, 1)
)
combined_plot

print(combined_plot)

pdf(file = "pre_analysis/manhatten_plot.pdf", 
    height= 7, 
    width=7)
print(combined_plot)
dev.off()
#



# plot n# genes -----------------------------------------------------------

p.coverag <- c.df%>% 
  distinct(contrast_id, n_genes)%>%
  ggplot(., aes(x= reorder(contrast_id,max.ranks), y= n_genes))+
  geom_col()+
  coord_flip()+
  theme_cowplot()+
  labs(y= "gene coverage", 
       x= "")

p.deg <- c.df %>%
  count(sig) %>%
  filter(sig)%>%
  ggplot(., aes(x= reorder(contrast_id,n), y= n))+
  geom_col()+
  coord_flip()+
  theme_minimal()+
  labs(y= "number of DEGs", 
       x= "")

pdf("plots/plot_coverage.pdf",
    width = 15, 
    height= 10)
plot_grid(p.coverag, p.deg)
dev.off()


