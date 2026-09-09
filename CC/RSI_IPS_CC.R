
########################### RSI in CC ###################################
install.packages("BiocManager")
BiocManager::install("hacksig")
library(hacksig)
library(dplyr)
library(ggplot2)
library(ggpubr)

ACC488<-read.delim("ACC488.txt")
ACC488<-ACC488[,-c(1,2)]
# read cluster info
ACC488_clinics<-read.csv("ACC488_clinics.csv")
ACC488_clusters <- ACC488_clinics %>% dplyr::select(c(ID, Cluster)) %>% arrange(Cluster)

# Removing rows where GeneSymbol is NA, #N/A, or empty and duplicated names
ACC488<- ACC488 %>%
  filter(gene_symbol != "" & gene_symbol != "#N/A" & !is.na(gene_symbol))
ACC488<- ACC488[!duplicated(ACC488$gene_symbol), ]
rownames(ACC488) <- ACC488$gene_symbol
ACC488 <- ACC488[,-1]


rsi_score_ACC<-hack_sig(
  ACC488,
  signatures = "all",
  method = "original",
  direction = "none",
  sample_norm = "raw",
  rank_norm = "none",
  alpha = 0.25)

#immunephenoscore
himp_ACC <- hack_immunophenoscore(ACC488, extract = "ips")


#Extracting the rsi score and combine it with immunephenoscore
rsi_only_ACC <- rsi_score_ACC %>%
  select(sample_id, eschrich2009_rsi) 

rsimps_matrix_ACC <- himp_ACC %>%
  left_join(rsi_only_ACC, by = c("sample_id" = "sample_id"))


# joining the rsi and imps with cluster info
final_data_ACC <- rsimps_matrix_ACC %>%
  left_join(ACC488_clusters, by = c("sample_id" = "ID"))


#Replacing cluster with cc
final_data_ACC$Cluster <- sub("Cl", "CC", final_data_ACC$Cluster)
colnames(final_data_ACC)[colnames(final_data_ACC) == "eschrich2009_rsi"] <- "RSI"
#write_xlsx(as.data.frame(final_data_ACC), "CC488_RSIdata.xlsx")
write.csv(final_data_ACC, "CC488_rsi&ips_scores.csv", row.names = TRUE)


# Defining pairwise comparisons for the clusters
my_comparisons_CC <- list(c("CC1", "CC2"), c("CC2", "CC3"), c("CC1", "CC3"))

# calculate KW effect size
kw <- kruskal.test(RSI ~ Cluster,data = final_data_ACC)
H <- as.numeric(kw$statistic)
k <- length(unique(final_data_ACC$Cluster))
n <- sum(!is.na(final_data_ACC$RSI) & !is.na(final_data_ACC$Cluster))
eta2 <- (H - k + 1) / (n - k)
eta2

effect_size <- round(eta2, 4)
p_label <- formatC(kw$p.value, format = "e", digits = 1)

label <- paste0("Kruskal-Wallis:~p==",p_label,
  "~'('~eta[H]^2~'='~",
  effect_size,
  "~')'"
)

# Boxplot for RSI across clusters with Kruskal-Wallis & Wilcoxon tests
final_data_ACC %>%
  ggplot(aes(x = Cluster, y = RSI, fill = Cluster, color = Cluster)) +  # Fill and color based on clusters
  geom_boxplot(outlier.shape = NA, width = 0.8, lwd = 1, aes(color = Cluster)) +  # Color for the box outline, white inside
  geom_jitter(width = 0.15, alpha = 0.4, color = "black") +  # Jitter for individual points with black color
  scale_fill_manual(values = c("white", "white", "white")) +  # Set the fill color of the boxes to white
  scale_color_manual(values = c("red", "blue", "green")) +  # Set the outline colors for the boxes
  stat_compare_means(comparisons = my_comparisons_CC, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test',size =5) +  # Pairwise comparisons (Wilcoxon)
  #stat_compare_means(method = 'kruskal.test', label.x.npc = 'centre', label.y.npc = 'top') +  # Kruskal-Wallis test for overall comparison
  theme_bw() +  # Clean white background theme
  annotate("text",x = 2,y = Inf,label = label,parse=T,
           vjust = 5,size = 4) +
  ggtitle("Radiosensitivity Index (RSI) Across CC Clusters") +
  theme(
    axis.title.x = element_blank(),
    legend.position = 'none',
    axis.text = element_text(size = 16),
    axis.text.y = element_text(size = 16), 
    axis.title.y = element_text(size = 16))
ggsave('plots/boxplot_RSI_CC488_Clusters_UPdated.pdf', height = 5, width = 7)



############################ IMMUNEPHENOGRAM ###################################

# 
# 
# ggplot(final_data_ACC, aes(x = reorder(sample_id, ips_score), y = ips_score)) +
#   geom_segment(aes(xend = sample_id, y = 0, yend = ips_score), color = "skyblue", linewidth = 0.6) +
#   geom_point(aes(color = IPS_status), size = 3) +
#   geom_hline(yintercept = 8, linetype = "dashed", color = "black") +
#   facet_wrap(~ Cluster, scales = "free_x") +
#   scale_color_manual(values = c("High" = "#1a9850", "Low" = "#fdae61")) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 20),
#     axis.title.y = element_text(size = 20),
#     panel.grid.major.y = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.y = element_blank(),
#     strip.text = element_text(size = 15, face = "bold"),
#     strip.background = element_rect(fill = "white")) +
#   labs(title = "ImmunePhenoScore Distribution by Cluster_CC", x = "", y = "IPS Score")
# 
# ggsave('Lolipop_IPS_CC488_updated.pdf', height = 10, width = 15)


########################### IPS with Statistics ################################
################################################################################
########################## CHI_square_test_CC ##################################
#Defining IPS status
final_data_ACC <- final_data_ACC %>%
  mutate( ips_score = as.numeric(as.character(ips_score)),  # ensure numeric
    IPS_status = ifelse(ips_score >=8 , "High", "Low"))
table(final_data_ACC$IPS_status, final_data_ACC$ips_score >= 8) #High 147_Low 341

table_ACC <- table(final_data_ACC$Cluster, final_data_ACC$IPS_status)
table_ACC

set.seed(1)
#Chi-square test
chi_ACC<-chisq.test(table_ACC,simulate.p.value = TRUE, B = 10000)
#simulate.p.value=TRUE ---> if some cells have low counts the normal approximation becomes unreliable.
#It tells R to ignore the theoretical chi-square distribution, and instead estimate the p value by simulation
#R randomly permutes the table(rows and columns) many times
#It calculate a chi-square statistics for each simulated table
#This gives a Monte Carlo p-value
#B = 10000 sets the number of random permutations (or simulations) to run.
#It means the Chi-square p-value is based on 10,000 simulated datasets, making the result more reliable when cell counts are low.

chisq.test(table_ACC)$expected #Expected count in each cell should be >5
#Expected counts are the numbers we expect under the null hypothesis of independence between variables

#chisq.test(table_ACC)
#view(table_ACC)
table_ACC<-as.data.frame(table_ACC)

# calculate KW effect size
kw <- kruskal.test(ips_score ~ Cluster,data = final_data_ACC)
H <- as.numeric(kw$statistic)
k <- length(unique(final_data_ACC$Cluster))
n <- sum(!is.na(final_data_ACC$ips_score) & !is.na(final_data_ACC$Cluster))
eta2 <- (H - k + 1) / (n - k)
eta2

effect_size <- round(eta2, 4)
p_label <- formatC(kw$p.value, format = "e", digits = 1)

label <- paste0("Kruskal-Wallis:~p==",p_label,
                "~'('~eta[H]^2~'='~",
                effect_size,
                "~')'")

ggplot(final_data_ACC, aes(x = Cluster, y = ips_score)) +
  geom_boxplot(aes(color = Cluster), fill = "white", outlier.shape = NA, alpha = 1, width = 0.8) +
  geom_jitter(aes(color = IPS_status), width = 0.2,height = 0, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black") +  # <-- threshold line
  #stat_compare_means(method = "kruskal.test", label.y = max(final_data_ACC$ips_score, na.rm = TRUE) + 1.5) +  # Global
  stat_compare_means(comparisons = my_comparisons_CC, method = "wilcox.test", label = "p.format", hide.ns = TRUE,size=10) +
  scale_color_manual(values = c(
    "CC1" = "red",
    "CC2" = "blue",
    "CC3" = "green",
    "High" = "#EE30A7",
    "Low" = "#999999" )) +
  labs(title = "IPS Score Distribution Across CC Clusters",
       x = "Cluster", y = "IPS Score",
       caption = paste("Chi-square test (IPS >= 8): p =", signif(chi_ACC$p.value, 3))) +
  annotate("text",x = 2,y = Inf,label = label,parse=T,vjust = 5,size = 6) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.y = element_text(size = 16), 
    axis.title.y = element_text(size = 16),
    plot.caption = element_text(size = 14),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 16))

ggsave("plots/Boxplot_IPS&Chitest_CC_updatedcolor.pdf",width = 8,height = 7)

#"Chi-square test (IPS > 8): p =" tells that the p-value shown is from a Chi-square test comparing IPS > 8 across groups.
#signif(chi_ACC$p.value):pulls the p-value result from the object chi_res, which is assumed to be the output of a previous chisq.test().
#signif(chi_ACC$p.value, 3) rounds the p-value to 3 significant digits. This keeps the number concise and readable (e.g., 0.08329 becomes 0.0833).

summary(final_data_ACC)
summary(rsimps_matrix_ACC)

table(final_data_ACC$ips_score)
# Count of IPS scores per cluster
table(final_data_ACC$ips_score, final_data_ACC$Cluster)


