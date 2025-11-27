load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")

######################
## remove non-cortex datasets
##########################

rm(DEl_exer_GSE39697_4d, DEl_exer_GSE39697_30d, DEl_exer_GSE39697_30d_gc, DEl_exer_GSE39697_4d_gc,
   DEl_stress_GSE151807_hipp, DEl_stress_GSE151807_amy, DEl_stress_GSE109315_vHPC_C57BL,
   DEl_stress_GSE109315_vHPC_DBA, DEl_stress_GSE109315_BNST_C57BL,
   DEl_stress_GSE109315_BNST_DBA, DEl_stress_GSE84183_DG, DEl_stress_GSE84183_CG,
   DEl_stress_GSE8262_hippx_CCH, DEl_stress_GSE8262_hippx_CIH,
   DE_stress_GSE198725_BLA, DE_stress_GSE198725_VTA )

## packages

suppressPackageStartupMessages({
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
})

## folder

results_dir <- "results/20251124"
dir.create(file.path(results_dir, "figures"),
           recursive = TRUE, showWarnings = FALSE)



##############################################
## GSEA for all DE/DEl objects
##############################################

species <- "Mus musculus"  # cambia in "Homo sapiens" se necessario

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
if (length(all_names) == 0) stop("Nessun oggetto DE/DEl_ trovato in ambiente.")

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

save(NES_mat,      file.path(results_dir, "NES_mat_C5.RData"))
save(dataset_type, file.path(results_dir, "NES_dataset_type.RData"))

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

pred_all <- predict(splsda_res, newdata = X_all)

# extract scores
if (is.list(pred_all) && "variates" %in% names(pred_all)) {
  var_comp <- pred_all[["variates"]]
  if (is.list(var_comp) && "X" %in% names(var_comp)) {
    scores_arr <- var_comp[["X"]]
  } else {
    scores_arr <- var_comp
  }
} else if (is.list(pred_all) && "x" %in% names(pred_all)) {
  scores_arr <- pred_all[["x"]]
} else {
  scores_arr <- pred_all
}

if (length(dim(scores_arr)) == 3L) {
  scores_mat <- scores_arr[, 1:2, 1, drop = FALSE]
} else {
  scores_mat <- scores_arr[, 1:2, drop = FALSE]
}

scores_all <- as.data.frame(scores_mat)
colnames(scores_all) <- c("comp1", "comp2")
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

distances <- scores_all %>%
  filter(group3 == "LSD") %>%
  rowwise() %>%
  mutate(
    dist_hormesis = sqrt(
      (comp1 - centroids$comp1[centroids$group3 == "Hormesis"])^2 +
        (comp2 - centroids$comp2[centroids$group3 == "Hormesis"])^2
    ),
    dist_stress = sqrt(
      (comp1 - centroids$comp1[centroids$group3 == "Chronic stress"])^2 +
        (comp2 - centroids$comp2[centroids$group3 == "Chronic stress"])^2
    )
  ) %>%
  ungroup()

dist_summary <- distances %>%
  summarise(
    mean_dist_hormesis = mean(dist_hormesis),
    mean_dist_stress   = mean(dist_stress)
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "target",
    values_to = "mean_distance"
  ) %>%
  mutate(
    target = recode(target,
                    mean_dist_hormesis = "Hormesis",
                    mean_dist_stress   = "Chronic stress")
  )

write.csv(
  distances,
  file.path(results_dir, "sPLSDA_LSD_distances_per_sample.csv"),
  row.names = FALSE
)
write.csv(
  dist_summary,
  file.path(results_dir, "sPLSDA_LSD_distances_summary.csv"),
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

p_dist <- ggplot(dist_summary,
                 aes(x = target, y = mean_distance, fill = target)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", mean_distance)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = cols[c("Hormesis", "Chronic stress")],
                    guide = "none") +
  labs(
    x = NULL,
    y = "Mean distance (sPLS-DA space)",
    title = "B  |  Proximity of LSD to hormesis vs chronic stress"
  ) +
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
