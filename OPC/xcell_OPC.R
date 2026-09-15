library(tibble)
library(immunedeconv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

########################### XCELL OPC553 ###################################
opc_553<- read.delim("HPV553.txt")
cl_opc_553<-read.delim("Clusters_HPV553.txt")
cl_opc_553_arranged <- cl_opc_553 %>% dplyr::select(c(ID, HPVpos_Clusters)) %>% arrange(HPVpos_Clusters)

# Removing rows with missing GeneSymbol
opc_553 <- opc_553 %>% filter(!is.na(GeneSymbol) & GeneSymbol != "#N/A") 

immunepopulation <- opc_553 %>% 
  dplyr::select(-ProbeID) %>%  
  column_to_rownames("GeneSymbol")

xCell_data <- suppressWarnings(immunedeconv::deconvolute_xcell(immunepopulation, arrays = TRUE))
write.csv(xCell_data, file='xCell_OPC.csv')

xCell_data <- xCell_data %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "cell_type", values_to = "score") %>%
  left_join(cl_opc_553_arranged, by = c("sample" = "ID"))
xCell_data$HPVpos_Clusters <- gsub('Cl1','OPC1',xCell_data$HPVpos_Clusters)
xCell_data$HPVpos_Clusters <- gsub('Cl2','OPC2',xCell_data$HPVpos_Clusters)
xCell_data$HPVpos_Clusters <- gsub('Cl3','OPC3',xCell_data$HPVpos_Clusters)

###--------------------------plot------------------------------------

kw_results <- xCell_data %>%
  filter(cell_type %in% c("ImmuneScore", 
                          "StromaScore", 
                          "MicroenvironmentScore")) %>%
  group_by(cell_type) %>%
  summarise(
    kw = list(kruskal.test(score ~ HPVpos_Clusters)),
    n = sum(!is.na(score) & !is.na(HPVpos_Clusters)),
    k = n_distinct(HPVpos_Clusters),
    .groups = "drop"
  ) %>%
  mutate(
    H = sapply(kw, function(x) as.numeric(x$statistic)),
    p = sapply(kw, function(x) x$p.value),
    eta2 = (H - k + 1) / (n - k),
    effect_size = round(eta2, 4),
    p_label = formatC(p, format = "e", digits = 1),
    label = paste0(
      "Kruskal-Wallis:~p==", p_label,
      "~'('~eta[H]^2~'='~", effect_size, "~')'"
    )
  )

#Total
#Kruskal test & Wilcoxon test graph
my_comparisons <- list(c("OPC2","OPC1"),c("OPC2","OPC3"),c("OPC3","OPC1")) #--> to define the pair-wise comparisons
xCell_data %>% dplyr::filter(cell_type %in% c("ImmuneScore", "StromaScore","MicroenvironmentScore")) %>%  ggplot(aes(x = HPVpos_Clusters, y = score)) +
  geom_boxplot(aes(color = HPVpos_Clusters), outlier.shape = NA, width = 0.8, lwd =0.8) +  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_color_manual(values = c("purple", "red", "blue")) +
  facet_wrap(vars(cell_type), scales = "free", ncol = 3) +
  stat_compare_means(comparisons = my_comparisons, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test', size = 10) + #--> add wilcoxon test to the boxplot
  #stat_compare_means(method = 'kruskal.test',label.x.npc = 'centre', label.y.npc = 'top') +     
  theme_bw() + ggtitle("xCell scores_OPC553") +
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
ggsave('plots/boxplot_xCellScore_OPC553.pdf', height = 5, width = 15)



#Testing single immune population in all the clusters
#eliminating the ImmuneScore etc.. from the table
m_woScore <- xCell_data[-grep('Score',xCell_data$cell_type),]

p <- unique(m_woScore$cell_type)

# Kruskal-wallis 
pval <- c()
for(i in p){
  print(i)
  subm <- m_woScore[which(m_woScore$cell_type == i),]
  k <- kruskal.test(list(subm$score[which(subm$HPVpos_Clusters == 'OPC1')],
                         subm$score[which(subm$HPVpos_Clusters== 'OPC2')],
                         subm$score[which(subm$HPVpos_Clusters== 'OPC3')]))
  pval <- c(pval,k$p.value)
}
kw_table <- data.frame(p,pval)
colnames(kw_table) <- c('cell_type','pval')
length(which(kw_table$pval < 0.05))
kw_table$adjPval <- p.adjust(kw_table$pval, method = 'fdr')
length(which(kw_table$adjPval < 0.05))
write.csv(kw_table, file = 'Kruskal-Wallis_xCell_population_table_opc553.csv',row.names = F)
# retrieving population with KW adjusted p-value < 0.05
kw_table <- kw_table[which(kw_table$adjPval < 0.05),]


#Kruskal test & Wilcoxon test graph
my_comparisons <- list(c("OPC2","OPC1"),c("OPC2","OPC3"),c("OPC3","OPC1")) #--> to define the pair-wise comparisons
xCell_data %>% dplyr::filter(cell_type %in% (kw_table$cell_type)) %>%  ggplot(aes(x = HPVpos_Clusters, y = score)) +
  geom_boxplot(aes(color = HPVpos_Clusters), outlier.shape = NA, width = 0.8, lwd =0.8) +  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_color_manual(values = c("purple", "red", "blue")) +
  facet_wrap(vars(cell_type), scales = "free", ncol = 3) +
  stat_compare_means(comparisons = my_comparisons, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test') + #--> add wilcoxon test to the boxplot
  stat_compare_means(method = 'kruskal.test',label.x.npc = 'centre', label.y.npc = 'top') +     
  theme_bw() + ggtitle("xCell_significant_immune population_OPC553") +
theme(axis.title.x = element_blank(),legend.position = 'none', axis.text = element_text(size=11),axis.title.y = element_text(size=11))
ggsave('plots/boxplot_xcell_signif_impop_OPC553.pdf', height = 65, width = 15,limitsize = FALSE)


