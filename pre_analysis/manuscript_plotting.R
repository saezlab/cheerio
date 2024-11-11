
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggdendro)

c.df= readRDS("data/contrasts_query_df_translated.rds")

contrasts_oi <- c.df %>%
  filter(!grepl("DCMvsNF", contrast_id)) %>%
  filter(!grepl("HCMvsDCM", contrast_id)) %>%
  pull(contrast_id)%>% unique()

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
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(x = "")

# Display the combined plot



combined_plot <- plot_grid(
  #dendrogram_plot,
  bar_plot,
  main_plot,
  ncol = 1,
  align = "v",
  rel_heights = c(0.1, 1)
)
combined_plot

print(combined_plot)

pdf(file = "pre_analysis/manhatten_plot.pdf", 
    height=10, 
    width=10)
print(combined_plot)
dev.off()
#


# contrast_meta -----------------------------------------------------------

meta= read_csv("~/Downloads/contrast and data overview - Sheet1.csv")
str_split(meta$contrast_id,"_")[][4]
meta %>% mutate(contrast_id2 = paste(species, 
                                     tissue,
                                     modality,
                                     resolution,
                                     `disease/model`, 
                                     sep= "_"))%>% pull(contrast_id2)

