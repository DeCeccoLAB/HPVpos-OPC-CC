

####################### Radiosensitivity for OPC553 ############################

install.packages("BiocManager")
BiocManager::install("hacksig")
library(hacksig)

opc_553<- read.delim("HPV553.txt")
cl_opc_553<-read.delim("Clusters_HPV553.txt")
cl_opc_553_arranged <- cl_opc_553 %>% dplyr::select(c(ID, HPVpos_Clusters)) %>% arrange(HPVpos_Clusters)

    
opc_553<-opc_553[,-1]

# Remove rows where GeneSymbol is NA, or empty and duplicated names
opc_553 <- opc_553 %>%
  filter(GeneSymbol != "" & GeneSymbol != "#N/A" & !is.na(GeneSymbol))
#opc_553 <- opc_553[!duplicated(opc_553$GeneSymbol), ]
rownames(opc_553) <- opc_553$GeneSymbol
opc_553<-opc_553[,-1]

rsi_score<-hack_sig(
  opc_553,
  signatures = "all",
  method = "original",
  direction = "none",
  sample_norm = "raw",
  rank_norm = "none",
  alpha = 0.25)


#immunephenoscore
himp <- hack_immunophenoscore(opc_553, extract = "ips")


#Extracting the rsi score and combine it with immunephenoscore
rsi_only <- rsi_score %>%
  select(sample_id, eschrich2009_rsi) 

rsimps_matrix <- himp %>%
  left_join(rsi_only, by = c("sample_id" = "sample_id"))


# joining the rsi and imps with cluster info
final_data <- rsimps_matrix %>%
  left_join(cl_opc_553, by = c("sample_id" = "ID"))

colnames(final_data)[colnames(final_data) == "HPVpos_Clusters"] <- "Cluster"


#Replacing cluster with opc
final_data$Cluster <- sub("Cl", "OPC", final_data$Cluster)
colnames(final_data)[colnames(final_data) == "eschrich2009_rsi"] <- "RSI"
write.csv(final_data, "OPC553_ips&rsi_scores.csv", row.names = TRUE)

my_comparisons_OPC<- list(c("OPC2", "OPC1"), c("OPC2", "OPC3"), c("OPC3", "OPC1"))

# calculate KW effect size
kw <- kruskal.test(RSI ~ Cluster,data = final_data)
H <- as.numeric(kw$statistic)
k <- length(unique(final_data$Cluster))
n <- sum(!is.na(final_data$RSI) & !is.na(final_data$Cluster))
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
final_data %>%
  ggplot(aes(x = Cluster, y = RSI, fill = Cluster, color = Cluster)) +  # Fill and color based on clusters
  geom_boxplot(outlier.shape = NA, width = 0.8, lwd = 1, aes(color = Cluster)) +  # Color for the box outline, white inside
  geom_jitter(width = 0.15, alpha = 0.4, color = "black") +  # Jitter for individual points with black color
  scale_fill_manual(values = c("white", "white", "white")) +  # Set the fill color of the boxes to white
  scale_color_manual(values = c("purple", "red", "blue")) +  # Set the outline colors for the boxes
  stat_compare_means(comparisons = my_comparisons_OPC, size = 4, label.x.npc = 0.3, label.y.npc = 0.5, method = 'wilcox.test',size =10) +  # Pairwise comparisons (Wilcoxon)
  #stat_compare_means(method = 'kruskal.test', label.x.npc = 'centre', label.y.npc = 'top') +  # Kruskal-Wallis test for overall comparison
  theme_bw() +  # Clean white background theme
  annotate("text",x = 2,y = Inf,label = label,parse=T,vjust = 5,size = 6) +
  ggtitle("Radiosensitivity Index (RSI) Across OPC Clusters") +  # Title
  theme(
    axis.title.x = element_blank(),
    legend.position = 'none',
    axis.text = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16))

ggsave('plots/boxplot_RSI_OPC553_Clusters_updated_colors.pdf', height = 5, width = 7)



############################ IMMUNEPHENOGRAM ###################################

# colnames(final_data)
# final_data<-final_data[,-5]
# 
# 
# ggplot(final_data, aes(x = reorder(sample_id, ips_score), y = ips_score)) +
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
#   labs(title = "ImmunePhenoScore Distribution by Cluster_OPC",x = "", y = "IPS Score")
# 
# ggsave("Lolipop_IPS_OPC_updated.pdf",width = 10,height = 10)


######################### CHI_square_test_OPC ##################################

#Defining IPS status
final_data <- final_data %>%
  mutate(
    ips_score = as.numeric(as.character(ips_score)),
    IPS_status = ifelse(ips_score >= 8 , "High", "Low"),
    Cluster = as.factor(Cluster))


table_OPC <- table(final_data$Cluster,final_data$IPS_status)

set.seed(1)
chi_OPC<-chisq.test(table_OPC,simulate.p.value = TRUE, B = 10000)
table_OPC

table(final_data$ips_score)
# Count of IPS scores per cluster
table(final_data$ips_score, final_data$Cluster)


#With chi_square test we examine observed values and expected values(what you would see if there is no association between the variables)
#Expected frequency is >5 ,so that each cell in table has enough data+the expected counts in each cell are not too small

#If p < 0.05: There is a significant association between cluster and IPS level.
#If p ≥ 0.05: No significant association.

#There is no evidence that the distribution of HIGH vs. LOW IPS differs significantly across the clusters


table_OPC<-as.data.frame(table_OPC)
#Stackbar plot
ggplot(table_OPC, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Cluster", y = "Proportion", fill = "IPS Level") +
  ggtitle("Proportion of HIGH/LOW IPS across OPC Clusters") +
  theme_minimal()


#Balloon plot 
ggballoonplot(table_OPC, x = "Var1", y = "Var2", size = "Freq",
              fill = "Freq", color = "black") +
  scale_size_area(max_size = 15) +
  ggtitle("Balloon Plot of IPS vs OPC Cluster") +
  theme_minimal()


########################## IPS with Statistics #################################
################################################################################
# calculate KW effect size
kw <- kruskal.test(ips_score ~ Cluster,data = final_data)
H <- as.numeric(kw$statistic)
k <- length(unique(final_data$Cluster))
n <- sum(!is.na(final_data$ips_score) & !is.na(final_data$Cluster))
eta2 <- (H - k + 1) / (n - k)
eta2

effect_size <- round(eta2, 4)
p_label <- formatC(kw$p.value, format = "e", digits = 1)

label <- paste0("Kruskal-Wallis:~p==",p_label,
                "~'('~eta[H]^2~'='~",
                effect_size,
                "~')'"
)
my_comparisons_OPC <- combn(levels(final_data$Cluster), 2, simplify = FALSE)
#simplify=FALSE ensures results are list of pairs not matrix, what we need for stat_compare_means
#2 for pairwise combination

# Plot
ggplot(final_data, aes(x = Cluster, y = ips_score)) +
  geom_boxplot(aes(color = Cluster), fill = "white", outlier.shape = NA, alpha = 1, width = 0.8) +
  geom_jitter(aes(color = IPS_status), width = 0.2,height = 0, size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "black") + 
  #stat_compare_means(method = "kruskal.test", label.y = max(final_data$ips_score, na.rm = TRUE) + 1.5) +
  stat_compare_means(comparisons = my_comparisons_OPC, method = "wilcox.test", label = "p.format", hide.ns = TRUE,size=5) +
  scale_color_manual(values = c(
    "OPC1" = "purple",
    "OPC2" = "red",
    "OPC3" = "blue",
    "High" = "#EE30A7",
    "Low" = "#999999")) +
  labs(title = "IPS Score Distribution Across OPC Clusters",
       x = "Clusters", y = "IPS Score",
       caption = paste("Chi-square test (IPS >=8): p =", signif(chi_OPC$p.value, 3))) +
  theme_bw() +
  annotate("text",x = 2,y = Inf,label = label,parse=T,vjust = 5,size = 5) +
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
ggsave("plots/boxplot_IPS&CHI_OPC_updatedcolors.pdf",width = 8,height = 7)


table(final_data$IPS_status)
summary(final_data$ips_score[final_data$IPS_status == "Low"])

unique(final_data$ips_score)
summary(final_data)

