load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")
source("code/Helpers.R")

library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
  library(stringr)
  library(mixOmics)
  library(ggrepel)
  library(patchwork)

results_dir <- "results/20251124"
dir.create(file.path(results_dir, "figures"),
           recursive = TRUE, showWarnings = FALSE)



##############################################
## GSEA for all DE/DEl objects
##############################################

species <- "Mus musculus"  

# Gene sets (C5 biol. process)
msig <- msigdbr(species = species, category = "C5")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# build ranked list (signed –log10 p)
get_ranked <- function(df) {
  x <- get_fc_p(df)
  
  p <- x$pval
  p[!is.finite(p) | p <= 0] <- .Machine$double.xmin
  
  lfc <- x$logFC
  lfc[!is.finite(lfc)] <- 0
  
  stat <- sign(lfc) * -log10(p)
  keep <- is.finite(stat) & !is.na(stat)
  stat  <- stat[keep]
  names(stat) <- rownames(df)[keep]
  
  stat <- sort(stat, decreasing = TRUE)
  stat
}

all_names <- ls(pattern = "^(DE|DEl)_")

# ranked lists for each DE
ranked_lists <- lapply(all_names, function(nm){
  get_ranked(get(nm))
})
names(ranked_lists) <- all_names

# GSEA/fgsea
fgsea_res_list <- lapply(names(ranked_lists), function(nm){
  rnk <- ranked_lists[[nm]]
  fgseaRes <- fgsea(pathways = gene_sets,
                    stats     = rnk,
                    nperm     = 5000,
                    minSize   = 10,
                    maxSize   = 500)
  fgseaRes <- as.data.frame(fgseaRes)
  fgseaRes$dataset <- nm
  fgseaRes
})
fgsea_res_list <- bind_rows(fgsea_res_list)
fgsea_res_list_tosave <- apply(fgsea_res_list,2,as.character)
write.csv(
  fgsea_res_list_tosave,
  file.path(results_dir, "GSEA_all_datasets_C5.csv"),
  row.names = FALSE
)

##############################################
## Build a matrix with all GSEA results
##############################################

NES_mat <- fgsea_res_list %>%
  dplyr::select(dataset, pathway = pathway, NES) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = NES) %>%
  as.data.frame()

row.names(NES_mat) <- NES_mat$pathway
NES_mat$pathway <- NULL

# dataset annotation
dataset_type <- case_when(
  str_detect(colnames(NES_mat), regex("DE_90min|DE_24h", ignore_case=TRUE)) ~ "LSD",
  str_detect(colnames(NES_mat), regex("exer", ignore_case=TRUE))       ~ "exercise",
  str_detect(colnames(NES_mat), regex("prec|hyp", ignore_case=TRUE))   ~ "preconditioning",
  str_detect(colnames(NES_mat), regex("fasting", ignore_case=TRUE))    ~ "fasting",
  str_detect(colnames(NES_mat), regex("stress",  ignore_case=TRUE))    ~ "stress",
  TRUE ~ "other"
)
names(dataset_type) <- colnames(NES_mat)

save(NES_mat,      file=file.path(results_dir, "NES_mat_C5.RData"))
save(dataset_type, file=file.path(results_dir, "NES_dataset_type.RData"))

##############################################
## sPLS-DA: hormesis vs chronic stress + LSD projection
##############################################

# hormesis = fasting + preconditioning + exercise
hormesis_levels <- c("fasting", "preconditioning", "exercise")

dataset_type_simpl <- dataset_type
dataset_type_simpl[dataset_type_simpl %in% hormesis_levels] <- "hormesis"

keep_train <- dataset_type_simpl %in% c("hormesis", "stress")
X_train    <- t(NES_mat[, keep_train, drop = FALSE])  # samples x pathways
Y_train    <- factor(dataset_type_simpl[keep_train],
                     levels = c("hormesis", "stress"))

set.seed(576906830)
splsda_res <- splsda(X_train, Y_train,
                     ncomp = 2,
                     keepX = c(20, 20))

# used pathways
feat_used <- rownames(splsda_res$loadings$X)

# projection of all datasets
keep_types <- c("fasting", "preconditioning", "exercise",
                "other", "stress", "LSD")
keep_cols  <- colnames(NES_mat)[dataset_type %in% keep_types]

X_all    <- t(NES_mat[feat_used, keep_cols, drop = FALSE])
type_all <- dataset_type[keep_cols]

group3 <- dplyr::case_when(
  type_all %in% hormesis_levels ~ "Hormesis",
  type_all == "stress"          ~ "Chronic stress",
  type_all == "LSD"             ~ "LSD",
  TRUE                          ~ "Other"
)
group3 <- factor(group3, levels = c("Hormesis", "Chronic stress", "LSD"))

##############################################
## Projection of all datasets 
##############################################

center_tr <- colMeans(X_train, na.rm = TRUE)
sd_tr     <- apply(X_train, 2, sd, na.rm = TRUE)
sd_tr[!is.finite(sd_tr) | sd_tr == 0] <- 1

## Align features between X_all and training loadings
loadX <- splsda_res$loadings$X[, 1:2, drop = FALSE]  # pathways x components (only 2 comps here)
needed_feats <- rownames(loadX)
missing_feats <- setdiff(needed_feats, colnames(X_all))
if (length(missing_feats) > 0) {
  stop("X_all is missing ", length(missing_feats), " features used by the model. Example: ",
       paste(head(missing_feats, 5), collapse = ", "))
}
X_all2 <- X_all[, needed_feats, drop = FALSE]

## Replace NA/Inf with 0
X_all2[!is.finite(X_all2)] <- NA_real_
X_all2[is.na(X_all2)] <- 0

## Scale X_all using TRAINING params
X_all_scaled <- sweep(sweep(X_all2, 2, center_tr[needed_feats], "-"), 2, sd_tr[needed_feats], "/")

## Manual projection: scores = X_scaled %*% loadings
scores_mat <- X_all_scaled %*% loadX
scores_mat[!is.finite(scores_mat)] <- 0
colnames(scores_mat) <- c("comp1", "comp2")

scores_all <- as.data.frame(scores_mat)
scores_all$dataset  <- keep_cols
scores_all$type_raw <- type_all
scores_all$group3   <- group3


##############################################
## LSD distance from centroids
##############################################

centroids <- scores_all %>%
  filter(group3 %in% c("Hormesis", "Chronic stress")) %>%
  group_by(group3) %>%
  summarise(
    comp1 = mean(comp1),
    comp2 = mean(comp2),
    .groups = "drop"
  )

##############################################
## Pairwise distances: each point = (LSD, hormesis/stress dataset)
##############################################

# LSD scores
lsd_scores <- scores_all %>%
  filter(group3 == "LSD") %>%
  dplyr::select(lsd_dataset = dataset,
                lsd_comp1  = comp1,
                lsd_comp2  = comp2)

# Hormesis + Chronic stress scores
hs_scores <- scores_all %>%
  filter(group3 %in% c("Hormesis", "Chronic stress")) %>%
  dplyr::select(hs_dataset = dataset,
                hs_comp1   = comp1,
                hs_comp2   = comp2,
                target     = group3)

pair_dist <- tidyr::expand_grid(
  lsd_dataset = lsd_scores$lsd_dataset,
  hs_dataset  = hs_scores$hs_dataset
) %>%
  left_join(lsd_scores, by = "lsd_dataset") %>%
  left_join(hs_scores, by = "hs_dataset") %>%
  mutate(
    distance = sqrt(
      (lsd_comp1 - hs_comp1)^2 +
        (lsd_comp2 - hs_comp2)^2
    ),
    target = factor(target,
                    levels = c("Hormesis", "Chronic stress"))
  )

wilcox_dist <- wilcox.test(distance ~ target, data = pair_dist, exact = FALSE)
wilcox_p    <- wilcox_dist$p.value

y_max <- max(pair_dist$distance, na.rm = TRUE)

dist_summary <- pair_dist %>%
  group_by(target) %>%
  summarise(
    mean_distance = mean(distance),
    sd_distance   = sd(distance),
    n_pairs       = n(),
    .groups = "drop"
  )

write.csv(
  pair_dist,
  file.path(results_dir, "sPLSDA_LSD_pairwise_distances.csv"),
  row.names = FALSE
)
write.csv(
  dist_summary,
  file.path(results_dir, "sPLSDA_LSD_pairwise_distances_summary.csv"),
  row.names = FALSE
)

##############################################
## Plot
##############################################

cols <- c("Hormesis" = "#1b9e77",
          "Chronic stress" = "#d95f02",
          "LSD" = "#7570b3")

shapes <- c("Hormesis" = 16,
            "Chronic stress" = 16,
            "LSD" = 17)

p_spls <- ggplot(scores_all, aes(x = comp1, y = comp2)) +
  stat_ellipse(
    data = subset(scores_all, group3 %in% c("Hormesis", "Chronic stress")),
    aes(color = group3),
    type = "norm",
    linetype = "dashed",
    size = 0.5,
    level = 0.68
  ) +
  geom_point(aes(color = group3, shape = group3),
             size = 2.8, alpha = 0.9) +
  geom_point(
    data = centroids,
    aes(x = comp1, y = comp2, color = group3),
    size = 4, shape = 4, stroke = 1.1
  ) +
  ggplot2::annotate("segment",
                    x = centroids$comp1[centroids$group3=="Hormesis"],
                    y = centroids$comp2[centroids$group3=="Hormesis"],
                    xend = centroids$comp1[centroids$group3=="Chronic stress"],
                    yend = centroids$comp2[centroids$group3=="Chronic stress"],
                    linetype = "dotted", size = 0.4) +
  geom_text_repel(
    data = subset(scores_all, group3 == "LSD"),
    aes(label = dataset),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.size = 0.2
  ) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(
    x = "sPLS-DA component 1",
    y = "sPLS-DA component 2",
    color = "Group",
    shape = "Group",
    title = "A  |  Hormesis vs chronic stress in sPLS-DA space",
    subtitle = "LSD conditions projected onto the hormesis–stress discriminant space"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9)
  )


p_dist <- ggplot(pair_dist,
                 aes(x = target, y = distance, fill = target)) +
  geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 2) +
  stat_summary(fun = median, geom = "point",
               size = 3, color = "black") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y))),
    vjust = -0.8,
    size = 3
  ) +
 
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_max * 1.05, yend = y_max * 1.05
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_max * 1.05, yend = y_max * 1.02
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_max * 1.05, yend = y_max * 1.02
  ) +
  annotate(
    "text",
    x = 1.5, y = y_max * 1.08,
    label = paste0("Wilcoxon p = ", signif(wilcox_p, 2)),
    size = 3
  ) +
  scale_fill_manual(values = cols[c("Hormesis", "Chronic stress")],
                    guide = "none") +
  labs(
    x = NULL,
    y = "Distance (sPLS-DA space)",
    title = "B  |  Proximity of LSD to hormesis vs chronic stress"
  ) +
  coord_cartesian(ylim = c(0, y_max * 1.2)) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



p_combined <- p_spls + p_dist + plot_layout(widths = c(2.2, 1))

p_combined

ggsave(
  filename = file.path(results_dir, "figures",
                       "sPLSDA_LSD_hormesis_vs_chronic_stress.png"),
  plot     = p_combined,
  width    = 9,
  height   = 4.5,
  dpi      = 300
)

ggsave(
  filename = file.path(results_dir, "figures",
                       "sPLSDA_LSD_hormesis_vs_chronic_stress.pdf"),
  plot     = p_combined,
  width    = 9,
  height   = 4.5
)



cent_stress <- centroids %>% 
  filter(group3 == "Chronic stress")

dist_axis <- scores_all %>%
  filter(group3 == "LSD") %>%
  mutate(
    dist_comp1 = abs(comp1 - cent_stress$comp1),
    dist_comp2 = abs(comp2 - cent_stress$comp2)
  )

dist_axis

##LSD distance to stress is higher for comp1

##############################################
## Pathways most relevant for comp1 and comp2
##############################################

loadings <- splsda_res$loadings$X  # matrix: pathways x components

load_df <- as.data.frame(loadings)
load_df$pathway <- rownames(load_df)

TOP_N <- 20

# comp1
top_comp1 <- load_df %>%
  filter(comp1 != 0) %>%                 
  arrange(desc(abs(comp1))) %>%          
  mutate(component = "comp1")

# comp2
top_comp2 <- load_df %>%
  filter(comp2 != 0) %>%                 
  arrange(desc(abs(comp2))) %>%
  mutate(component = "comp2")

top_pathways <- bind_rows(top_comp1, top_comp2)


shorten_pathway <- function(x, max_words = 4) {
  x <- gsub("^GO(BP|CC|MF)?_", "", x)
  words <- unlist(strsplit(x, "_"))
  if (length(words) > max_words) words <- c(words[1:max_words], "...")
  paste(words, collapse = " ")
}

top_pathways$pathway_short <- sapply(top_pathways$pathway, shorten_pathway)

top_pathways <- top_pathways %>%
  mutate(value = ifelse(component == "comp1", comp1, comp2))

write.xlsx(top_pathways, file="top_pathways_sPLSD.xlsx")


df_comp1 <- top_pathways %>% 
  filter(component == "comp1") %>%
  mutate(pathway = factor(pathway, levels = pathway[order(value)])) 

p_comp1 <- ggplot(df_comp1,
                  aes(x = value,              
                      y = pathway,
                      color = value > 0)) +    
  geom_point(size = 3) +
  scale_color_manual(
    values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02"),
    labels = c("FALSE" = "Negative loading", "TRUE" = "Positive loading")
  ) +
  labs(
    x = "Loading",
    y = "Pathway",
    color = "Direction",
    title = "Top pathways contributing to sPLS-DA — Component 1"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9)
  )

p_comp1

ggsave(
  filename = file.path(results_dir, "figures", "sPLSDA_loadings_comp1.png"),
  plot     = p_comp1,
  width    = 10,
  height   = 6,
  dpi      = 300
)

ggsave(
  filename = file.path(results_dir, "figures", "sPLSDA_loadings_comp1.pdf"),
  plot     = p_comp1,
  width    = 10,
  height   = 6
)


##############################################
## Pathway-level summary: NES + mean log2FC
##############################################

pathway_selected <- "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"
#"GOBP_POSITIVE_REGULATION_OF_HORMONE_SECRETION"

df_nes <- fgsea_res_list %>%
  dplyr::filter(pathway == pathway_selected) %>%
  dplyr::mutate(
    dataset_type       = dataset_type[dataset],
    dataset_type_simpl = dataset_type_simpl[dataset]
  )


library(ggpubr)

df_nes$dataset_type_simpl <- factor(
  df_nes$dataset_type_simpl,
  levels = c("LSD", "hormesis", "stress")
)

p_nes_box <- ggplot(df_nes,
                    aes(x = dataset_type_simpl,
                        y = NES,
                        fill = dataset_type_simpl)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("hormesis", "stress")),
    label = "p.format",  
    label.y = max(df_nes$NES, na.rm = TRUE) * 1.05
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "NES (GSEA)",
    title = paste0("NES of ", pathway_selected),
    subtitle = "Per dataset type"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

p_nes_box


ggsave(
  filename = file.path(
    results_dir, "figures",
    paste0("Pathway_", pathway_selected, "_NES_boxplot.png")
  ),
  plot = p_nes_box,
  width = 4,
  height = 4,
  dpi = 300
)


ggsave(
  filename = file.path(
    results_dir, "figures",
    paste0("Pathway_", pathway_selected, "_NES_boxplot.pdf")
  ),
  plot = p_nes_box,
  width = 4,
  height = 4
)

