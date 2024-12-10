
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggdendro)

c.df= readRDS("app_data/contrasts_query_df_translated3.rds")

contrasts_oi <- c.df %>%
  filter(!grepl("snRNA_", contrast_id) | grepl("snRNA_CM|snRNA_FB", contrast_id)) %>%
  #filter(!grepl("HCMvsDCM", contrast_id)) %>%
  pull(contrast_id)%>% unique()

contrasts_oi

c.df2 <- c.df %>% 
  filter(contrast_id %in% contrasts_oi)
# plot histogram of contrasts and genes -----------------------------------

length(unique(c.df$gene))
contrast_comparisons<- c.df%>% 
  group_by(gene)%>%
  summarize(contrast_count = n_distinct(contrast_id))

sum(contrast_comparisons$contrast_count >9)

p.hist.counts <- contrast_comparisons%>% 
  ggplot(aes(x= contrast_count))+
  geom_histogram(bins = 30)+
  #scale_y_log10()+
  labs(x= "number of contrasts", 
       y= "number of genes")+
  theme_cowplot()

pdf("figures/histogram_number_of_contrasts.pdf", 
    width= 3, height= 3)
p.hist.counts
dev.off()
# plot contrast overview --------------------------------------------------

#create a filtered version of the contrast data frame
c.df2<- c.df  %>%
  filter(contrast_id %in% contrasts_oi)

cc_colors <- c("A" = "#3A6EA5","B" = "#E29F91",  "C" = "#AB82A4", "D" = "#6F9F8A")

label_colors <- c.df2 %>%
  distinct(contrast_id, cc) %>%
  arrange(cc, contrast_id)%>%
  mutate(color = cc_colors[as.character(cc)])

unname(label_colors$color)
label_colors$contrast_id

# Filter and count genes with FDR < 0.05 for each contrast
count_df <- c.df2 %>%
  filter(contrast_id %in% contrasts_oi) %>%
  filter(FDR < 0.05) %>%
  group_by(contrast_id) %>%
  summarise(num_genes = n())%>%
  complete(contrast_id = unique(contrasts_oi), fill = list(num_genes = 0))%>%
  left_join(c.df2 %>%distinct(contrast_id, cc), by= "contrast_id")%>%
  mutate(contrast_id = factor(contrast_id, levels= label_colors$contrast_id)) %>%
  arrange(cc, contrast_id)
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
        axis.title.y = element_text(size = 10)
        ) +
  labs(y = "DEG\ncount")
bar_plot

# Main plot
main_plot <- c.df %>%
  filter(contrast_id %in% contrasts_oi) %>%
  mutate(color = factor(ifelse(FDR < 0.05 & logFC > 0, "up",
                               ifelse(FDR < 0.05 & logFC < 0, "dn", "unreg")),
                        levels = c("unreg", "up", "dn"))) %>% 
  arrange(color) %>% # Set order of plotting
  mutate(contrast_id = factor(contrast_id, levels = label_colors$contrast_id))%>% # Set order based on clustering
  ggplot(aes(x = contrast_id, y = logFC, color = color)) +
  geom_jitter(size = 0.2, width = 0.3, alpha = 0.3) +
  theme_cowplot() +
  scale_color_manual(values = c("grey", "red", "blue")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust= 0.5,
                                   color = label_colors$color),
        legend.position = "none") +
  labs(x = "")
  
main_plot

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
#c.df2 %>% filter(grepl("Endo", contrast_id))
print(combined_plot)

pdf(file = "figures/lfc_jitter_plot.pdf", 
    height= 6, 
    width=6)
print(combined_plot)
dev.off()
#



# plot n# genes -----------------------------------------------------------

p.coverag <- c.df%>% 
  distinct(contrast_id, n_genes)%>%
  #count(contrast_id)%>%
  ggplot(., aes(x= reorder(contrast_id,n_genes), y= n_genes))+
  geom_col(color="black", fill="grey", width= 0.5)+
  coord_flip()+
  theme_cowplot()+
  labs(y= "number of genes", 
       x= "")

p.deg <- c.df %>%
  group_by(contrast_id)%>%
  count(sig) %>%
  filter(sig)%>%
  ggplot(., aes(x= reorder(contrast_id,n), y= n))+
  geom_col()+
  coord_flip()+
  theme_cowplot()+
  labs(y= "number of DEGs", 
       x= "")

pdf("figures/plot_coverage.pdf",
    width = 8, 
    height= 9)
p.coverag
dev.off()


# plots for UCK2 ----------------------------------------------------------

genes_oi<- c("NPPA", "UCK2")
sub_ranks = ranks %>%
  filter(gene %in% genes_oi)

max_rank = max(ranks$rank)
library(ggplot2)

# Improved plot
dens <- ranks %>%
  ggplot(aes(x = mean_t)) +
  # Add density curve with a more distinct line style
  stat_density(geom = "line", color = "steelblue", size = 1.2) +
  # Highlight specific genes with a rug plot
  geom_rug(data = sub_ranks, aes(x = mean_t), color = "darkred", size = 1.5) +
  # Annotate the highlighted genes
  geom_text(data = sub_ranks, aes(x = mean_t, y = 0.03, label = gene), 
            color = "darkred", fontface = "bold", hjust = 0.5, vjust= 1, size = 4) +
  # Classic theme with subtle gridlines
  #theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  # Adjust labels and title for clarity
  labs(
    #title = "Distribution of Mean T-Values Across Genes",
    x = "Mean t-Value",
    y = "Density"
  )
dens
uck2_hw <- plot_hw_association(HW_DF, genes_oi)


##save pdfs
pdf("figures/uck2_hw_asso.pdf",
    height= 3.2, width = 5)
uck2_hw
dev.off()

##save pdfs
pdf("figures/uck2_t_val_dit_hf.pdf",
    height= 3.2, width = 3.5)
dens
dev.off()
