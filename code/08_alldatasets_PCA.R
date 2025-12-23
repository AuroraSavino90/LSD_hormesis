library(dplyr)
library(ggplot2)
library(ggrepel)
date        <- "20251124"
results_dir <- file.path("results", date)

load(file.path(results_dir, "NES_mat_C5.RData"))
load(file.path(results_dir, "NES_dataset_type.RData"))

# hormesis = fasting + preconditioning + exercise
hormesis_levels <- c("fasting", "preconditioning", "exercise")

dataset_type_pca <- dplyr::case_when(
  dataset_type %in% hormesis_levels ~ "Hormesis",
  dataset_type == "stress"          ~ "Chronic stress",
  dataset_type == "LSD"             ~ "LSD",
  TRUE                              ~ "Other"
)

#########################################################
## Preprocessing: remove pathways with many NAs
#########################################################

NES_filt <- NES_mat %>%
  filter(rowSums(is.na(.)) < ncol(.) * 0.5)

# NA -> 0
NES_filt[is.na(NES_filt)] <- 0

#########################################################
## PCA
#########################################################



pca_res <- prcomp(t(NES_filt),   
                  center = TRUE,
                  scale. = TRUE)

var_expl <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)


pca_df <- data.frame(
  dataset = rownames(pca_res$x),
  PC1     = pca_res$x[, 1],
  PC2     = pca_res$x[, 2],
  type    = dataset_type_pca
)

pca_df$type <- factor(
  pca_df$type,
  levels = c("Hormesis", "Chronic stress", "LSD", "Other")
)

#########################################################
## Plot PCA
#########################################################

cols <- c(
  "Hormesis"        = "#1b9e77",
  "Chronic stress"  = "#d95f02",
  "LSD"             = "#7570b3",
  "Other"           = "grey70"
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = cols) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of NES_mat (C5 pathways)",
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)"),
    color = "Group"
  )

print(p_pca)

ggsave(
  file.path(results_dir, "figures", "PCA_NES_C5.png"),
  p_pca, width = 8, height = 6, dpi = 300
)


ggsave(
  file.path(results_dir, "figures", "PCA_NES_C5.pdf"),
  p_pca, width = 8, height = 6
)
