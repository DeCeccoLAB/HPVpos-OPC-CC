library(tibble)
library(immunedeconv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
##################### XCELL CC488 ################################
install.packages("devtools")
devtools::install_github('dviraran/xCell')
library(xCell)

install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)

install.packages("devtools")
devtools::install_github('dviraran/xCell',force = TRUE)
library(xCell)

install.packages("EPIC")
library(EPIC)

#--------------------------------------------------------------------------
ACC488<-read.delim("ACC488.txt")
ACC488_clinics<-read.csv("ACC488_clinics.csv")
ACC488_clusters <- ACC488_clinics %>% dplyr::select(c(ID, Cluster)) %>% arrange(Cluster)

immunepopulation<- ACC488 %>% dplyr::select(-c(EntrezID, gene_id)) %>% column_to_rownames("gene_symbol")


xCell_data <- suppressWarnings(immunedeconv::deconvolute_xcell(immunepopulation, arrays = TRUE))
write.csv(xCell_data, file='xCell_CC.csv')

xCell_data <- xCell_data %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "score") %>%
  left_join(ACC488_clusters, by = c("sample" = "ID"))
xCell_data$Cluster <- gsub('Cl1','CC1',xCell_data$Cluster)
xCell_data$Cluster <- gsub('Cl2','CC2',xCell_data$Cluster)
xCell_data$Cluster <- gsub('Cl3','CC3',xCell_data$Cluster)

##------------------------plot--------------------------------
kw_results <- xCell_data %>%
  filter(cell_type %in% c("ImmuneScore", 
                          "StromaScore", 
                          "MicroenvironmentScore")) %>%
  group_by(cell_type) %>%
  summarise(
    kw = list(kruskal.test(score ~ Cluster)),
    n = sum(!is.na(score) & !is.na(Cluster)),
    k = n_distinct(Cluster),
    .groups = "drop"
  ) %>%
  mutate(
    H = sapply(kw, function(x) as.numeric(x$statistic)),
    p = sapply(kw, function(x) x$p.value),
    eta2 = pmax((H - k + 1) / (n - k), 0),
    effect_size = round(eta2, 4),
    p_label = formatC(p, format = "e", digits = 1),
    label = paste0(
      "Kruskal-Wallis:~p==", p_label,
      "~'('~eta[H]^2~'='~", effect_size, "~')'"
    )
  )
kw_results

#Total
#Kruskal test & Wilcoxon test graph
my_comparisons <- list(c("CC1","CC2"),c("CC1","CC3"),c("CC2","CC3")) #--> to define the pair-wise comparisons
xCell_data %>% dplyr::filter(cell_type %in% c("ImmuneScore", "StromaScore","MicroenvironmentScore")) %>%  ggplot(aes(x = Cluster, y = score)) +
  geom_boxplot(aes(color = Cluster), outlier.shape = NA, width = 0.8, lwd =0.8) +  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_color_manual(values = c("red", "blue", "green")) +
  facet_wrap(vars(cell_type), scales = "free", ncol = 3) +
  stat_compare_means(comparisons = my_comparisons, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test',size=10) + #--> add wilcoxon test to the boxplot
  #stat_compare_means(method = 'kruskal.test',label.x.npc = 'centre', label.y.npc = 'top') +     
  theme_bw() + ggtitle("xCell scores_CC") +
  geom_text(data = kw_results,aes(x = 2,y = Inf,label = label),inherit.aes = FALSE,parse = TRUE,vjust = 5,
            size = 5
  ) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16),
        strip.text = element_text(size = 16))
ggsave('plots/boxplot_xCellScore.pdf', height = 5, width = 15)


#Testing single immune population in all the clusters
#eliminating the ImmuneScore etc.. from the table
m_woScore <- xCell_data[-grep('Score',xCell_data$cell_type),]

p <- unique(m_woScore$cell_type)

# Kruskal-wallis #preparing input for the test
pval <- c()
for(i in p){
  print(i)
  subm <- m_woScore[which(m_woScore$cell_type == i),]
  k <- kruskal.test(list(subm$score[which(subm$Cluster == 'Cl1')],
                         subm$score[which(subm$Cluster== 'Cl2')],
                         subm$score[which(subm$Cluster== 'Cl3')]))
  pval <- c(pval,k$p.value)
}
kw_table <- data.frame(p,pval)
colnames(kw_table) <- c('cell_type','pval')
length(which(kw_table$pval < 0.05))
kw_table$adjPval <- p.adjust(kw_table$pval, method = 'fdr')
length(which(kw_table$adjPval < 0.05))
write.csv(kw_table, file = 'Kruskal-Wallis_xCell_population_table.csv',row.names = F)
# retrieving population with KW adjusted p-value < 0.05
kw_table <- kw_table[which(kw_table$adjPval < 0.05),]


#Kruskal test & Wilcoxon test graph
my_comparisons <- list(c("Cl1","Cl2"),c("Cl1","Cl3"),c("Cl2","Cl3")) #--> to define the pair-wise comparisons
xCell_data %>% dplyr::filter(cell_type %in% (kw_table$cell_type)) %>%  ggplot(aes(x = Cluster, y = score)) +
  geom_boxplot(aes(color = Cluster), outlier.shape = NA, width = 0.8, lwd =0.8) +  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_color_manual(values = c("red", "blue", "green")) +
  facet_wrap(vars(cell_type), scales = "free", ncol = 3) +
  stat_compare_means(comparisons = my_comparisons, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test') + #--> add wilcoxon test to the boxplot
  stat_compare_means(method = 'kruskal.test',label.x.npc = 'centre', label.y.npc = 'top') +     
  theme_bw() + ggtitle("xCell_significant_immune population")
theme(axis.title.x = element_blank(),legend.position = 'none', axis.text = element_text(size=11),axis.title.y = element_text(size=11))
ggsave('boxplot_xcell_signif_impop.pdf', height = 48, width = 15)


################################################################################
#Boxplot of single immune population separately
library(ggplot2)
library(dplyr)
library(ggpubr)

# Getting significant immune populations
signif_celltypes <- kw_table$cell_type

# Loop through each significant immune cell type
for (ct in signif_celltypes) { #filtering data to one cell type, and loop over these cell types, drawing one boxplot per loop
  sub_data <- xCell_data %>% filter(cell_type == ct)
  
  p <- ggplot(sub_data, aes(x = Cluster, y = score)) +
    geom_boxplot(aes(color = Cluster), outlier.shape = NA, width = 0.8, lwd = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.4) +
    scale_color_manual(values = c("red", "blue", "green")) +
    stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test') +
    stat_compare_means(method = 'kruskal.test') +
    theme_bw() +
    labs(title = ct, x = "Cluster", y = "Score") +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none")
  # Save each plot as a separate PDF file
  ggsave(filename = paste0("boxplot_", gsub(" ", "_", ct), ".pdf"), plot = p, width = 6, height = 5)
}

