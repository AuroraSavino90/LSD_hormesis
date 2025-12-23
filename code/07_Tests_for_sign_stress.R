library(dplyr)
library(purrr)
library(tidyr)
library(fgsea)
library(ggplot2)

load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")

## General params & paths
date        <- "20251124"
results_dir <- file.path("results", date)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)


all_obj <- ls(pattern = "^(DE|DEl)_")

refs <- c("DE_90min", "DE_24h")

non_ref <- setdiff(all_obj, refs)

fasting_targets <- non_ref[grepl("fasting", non_ref, ignore.case = TRUE)]
exercise_targets <- non_ref[grepl("exer", non_ref, ignore.case = TRUE)]
precond_targets <- non_ref[grepl("prec", non_ref, ignore.case = TRUE) |
                             grepl("hyp",  non_ref, ignore.case = TRUE)]
stress_targets <- non_ref[grepl("stress", non_ref, ignore.case = TRUE)]
hormesis_targets <- unique(c(fasting_targets, exercise_targets, precond_targets))


############################################
## GSEA of -log10(p) * sign(log2FC)
############################################


.col_by2 <- function(df, pat) {
  i <- grep(pat, tolower(colnames(df)))
  if (length(i) == 0) return(NA_character_)
  colnames(df)[i[1]]
}

get_fc_p2 <- function(df) {
  stopifnot(!is.null(rownames(df)))
  fc_col <- .col_by2(df, "log2foldchange|logfc")
  p_col  <- .col_by2(df, "^(pvalue|p\\.val|p\\.value|p)$|[^a-z]pvalue[^a-z]|[^a-z]p\\.value[^a-z]")
  if (is.na(fc_col)) stop("No fold-change column (log2FoldChange/logFC) in one object.")
  if (is.na(p_col))  stop("No p-value column (pvalue/P.Value/p) in one object.")
  out <- df[, c(fc_col, p_col), drop = FALSE]
  colnames(out) <- c("logFC", "pval")
  out
}

get_lsd_sets <- function(ref_name,
                         p_thresh = 0.05,
                         fc_thresh = 0,
                         env = .GlobalEnv) {
  stopifnot(exists(ref_name, envir = env))
  ref_df <- get(ref_name, envir = env)
  stopifnot(!is.null(rownames(ref_df)))
  
  
  fc_col <- .col_by2(ref_df, "log2foldchange|logfc")
  padj_col <- .col_by2(ref_df, "padj|adj\\.p|fdr|qvalue")
  
  if (is.na(fc_col))  stop("No fold-change column (log2FoldChange/logFC) in LSD reference object.")
  if (is.na(padj_col)) stop("No adjusted p-value column (padj/adj.P/FDR/qvalue) in LSD reference object.")
  
  x <- ref_df[, c(fc_col, padj_col), drop = FALSE]
  colnames(x) <- c("logFC", "padj")
  
  sig_idx <- which(abs(x$logFC) > fc_thresh & x$padj < p_thresh)
  
  genes_sig <- rownames(ref_df)[sig_idx]
  up_idx    <- sig_idx[x$logFC[sig_idx] > 0]
  down_idx  <- sig_idx[x$logFC[sig_idx] < 0]
  
  list(
    LSD_up   = rownames(ref_df)[up_idx],
    LSD_down = rownames(ref_df)[down_idx]
  )
}



make_rank_vector <- function(df) {
  x <- get_fc_p2(df)
  stat <- -log10(x$pval + 1e-300) * sign(x$logFC)
  names(stat) <- rownames(df)
  sort(stat, decreasing = TRUE)
}

run_gsea_block <- function(ref_name,
                           target_names,
                           p_thresh = 0.05,
                           fc_thresh = 0,
                           nperm = 5000,
                           env = .GlobalEnv) {
  lsd_sets <- get_lsd_sets(ref_name, p_thresh = p_thresh,
                           fc_thresh = fc_thresh, env = env)
  gene_sets <- list(
    LSD_up   = lsd_sets$LSD_up,
    LSD_down = lsd_sets$LSD_down
  )
  
  rows <- list()
  
  for (tgt_name in target_names) {
    
    tgt_df <- get(tgt_name, envir = env)
    stats_vec <- make_rank_vector(tgt_df)
    
    gsea_res <- fgsea(
      pathways = gene_sets,
      stats    = stats_vec,
      nperm    = nperm
    )
    
    rows[[paste(ref_name, tgt_name, sep = "_")]] <-
      gsea_res %>%
      as_tibble() %>%
      mutate(
        lsd_ref = ref_name,
        dataset = tgt_name
      ) %>%
      dplyr::select(lsd_ref, dataset, pathway, NES, pval, padj, size, leadingEdge)
  }
  
  bind_rows(rows)
}

############################################
## Run
############################################

stress_targets <- stress_targets[stress_targets %in% ls()]


gsea_24 <- run_gsea_block(
  ref_name     = "DE_24h",
  target_names = stress_targets,
  p_thresh     = 0.05,
  fc_thresh    = 0,
  nperm        = 5000
)

gsea_90 <- run_gsea_block(
  ref_name     = "DE_90min",
  target_names = stress_targets,
  p_thresh     = 0.05,
  fc_thresh    = 0,
  nperm        = 5000
)



gsea_24$lsd_ref <- "24h"
gsea_90$lsd_ref <- "90min"

gsea_all <- bind_rows(gsea_24, gsea_90) %>%
  mutate(
    lsd_ref   = factor(lsd_ref, levels = c("90min", "24h")),
    pathway   = factor(pathway, levels = c("LSD_up", "LSD_down")),
    dataset_clean = gsub("^DEl?_+", "", dataset) %>% gsub("_+", " ", .),
    minus_log10_p = -log10(pval + 1e-300),
    minus_log10_p_capped = pmin(minus_log10_p, 10)
  )

# (a) Heatmap NES per dataset × gene-set, separando 90 vs 24
ordering_gsea <- gsea_all %>%
  group_by(lsd_ref, dataset_clean) %>%
  summarise(
    rank_key = max(abs(NES) * (padj < 0.25), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(lsd_ref, desc(rank_key))

gsea_all <- gsea_all %>%
  left_join(ordering_gsea, by = c("lsd_ref", "dataset_clean")) %>%
  group_by(lsd_ref) %>%
  mutate(
    dataset_clean = factor(
      dataset_clean,
      levels = unique(dataset_clean[order(-rank_key, dataset_clean)])
    )
  ) %>%
  ungroup()

p_gsea_heat <- ggplot(
  gsea_all,
  aes(x = pathway, y = dataset_clean)
) +
  geom_tile(aes(fill = NES), width = 0.95, height = 0.95, color = "grey80") +
  geom_point(aes(size = minus_log10_p_capped),
             shape = 21, fill = "white", color = "black", stroke = 0.2) +
  scale_fill_gradient2(
    name   = "NES",
    low    = "#3B4CC0",
    mid    = "white",
    high   = "#B40426",
    limits = c(-3, 3),
    oob    = scales::squish
  ) +
  scale_size_continuous(
    name = expression(-log[10](p)),
    range = c(0.7, 4)
  ) +
  facet_wrap(~ lsd_ref, ncol = 2, scales = "free_y") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    legend.key.height = unit(0.4, "cm")
  ) 

ggsave(
  file.path(results_dir, "figures", "LSD_vs_stress_GSEA_signedStats_ALL.png"),
  p_gsea_heat,
  width = 7, height = max(4.5, 0.35 * length(unique(gsea_all$dataset_clean)) + 2),
  dpi = 300
)

ggsave(
  file.path(results_dir, "figures", "LSD_vs_stress_GSEA_signedStats_ALL.pdf"),
  p_gsea_heat,
  width = 7, height = max(4.5, 0.35 * length(unique(gsea_all$dataset_clean)) + 2),
  dpi = 300
)


############################################
## Proportion of NES>0 / NES<0
############################################

gsea_sig <- gsea_all %>%
  mutate(
    significant = pval < 0.05,
    NES_dir = case_when(
      NES >  0 ~ "NES > 0",
      NES <  0 ~ "NES < 0",
      TRUE      ~ "NES = 0"
    )
  ) %>%
  filter(significant, NES_dir != "NES = 0")

# Conta per tempo, gene set, direzione NES
gsea_counts <- gsea_sig %>%
  group_by(lsd_ref, pathway, NES_dir) %>%
  summarise(
    n_tests = n(),
    .groups = "drop"
  )

gsea_counts




gsea_sig <- gsea_all %>%
  mutate(
    significant = pval < 0.05,
    NES_dir = case_when(
      NES >  0 ~ "NES > 0",
      NES <  0 ~ "NES < 0",
      TRUE      ~ "NES = 0"
    )
  ) %>%
  filter(significant, NES_dir != "NES = 0")

gsea_counts <- gsea_sig %>%
  group_by(lsd_ref, pathway, NES_dir) %>%
  summarise(
    n_tests = n(),
    .groups = "drop"
  )

p_gsea_bar_prop <- ggplot(gsea_counts,
                          aes(x = pathway, y = n_tests, fill = NES_dir)) +
  geom_col(position = "fill", color = "black") +
  facet_wrap(~ lsd_ref) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Proportion of significant tests (p < 0.05)",
    fill = "NES direction",
    title = "Direction of sign (NES>0 / NES<0) in significant GSEA"
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_signedStats_sigCounts_barplot_100pct.png"),
  p_gsea_bar_prop,
  width = 6, height = 4.5, dpi = 300
)

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_signedStats_sigCounts_barplot_100pct.pdf"),
  p_gsea_bar_prop,
  width = 6, height = 4.5, dpi = 300
)



############################################
## Chisquare
############################################

chi_sq_tests <- gsea_counts %>%
  pivot_wider(
    names_from = NES_dir,
    values_from = n_tests,
    values_fill = 0
  ) %>%
  group_by(lsd_ref) %>%
  summarise(
    chi_sq_p = chisq.test(
      matrix(c(
        `NES > 0`[pathway == "LSD_up"],
        `NES < 0`[pathway == "LSD_up"],
        `NES > 0`[pathway == "LSD_down"],
        `NES < 0`[pathway == "LSD_down"]
      ), nrow = 2, byrow = TRUE)
    )$p.value,
    .groups = "drop"
  )

print(chi_sq_tests)

############################################
## Boxplot sign(NES) * -log10(pvalue)
############################################

gsea_signed <- gsea_all %>%
  mutate(
    signed_stat = sign(NES) * -log10(pval + 1e-300),
    pathway = factor(pathway, levels = c("LSD_up", "LSD_down"))
  )

p_gsea_box_signed <- gsea_signed %>%
  ggplot(aes(x = pathway, y = signed_stat, fill = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.6, outlier.size = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  facet_wrap(~ lsd_ref) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = expression(sign(NES) %*% -log[10](pvalue)),
    title = "Strength and direction of GSEA enrichment"
  )

gsea_signed <- gsea_signed %>%
  mutate(
    pathway = factor(pathway, levels = c("LSD_up", "LSD_down"))
  )

wilcox_gsea <- gsea_signed %>%
  group_by(lsd_ref) %>%
  summarise(
    p_wilcox = wilcox.test(
      signed_stat ~ pathway,
      data = cur_data(),
      alternative = "greater"   # H1: LSD_up > LSD_down
    )$p.value,
    y_max = max(signed_stat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Wilcoxon (LSD_up > LSD_down) p = ", signif(p_wilcox, 3)),
    x_pos = 1.5,
    y_pos = y_max * 1.05
  )

p_gsea_box_signed_w <- p_gsea_box_signed +
  geom_text(
    data = wilcox_gsea,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    vjust = 0
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_signedStat_boxplot_wilcoxon_greater.png"),
  p_gsea_box_signed_w,
  width = 6, height = 4.7, dpi = 300
)

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_signedStat_boxplot_wilcoxon_greater.pdf"),
  p_gsea_box_signed_w,
  width = 6, height = 4.7, dpi = 300
)

############################################
## Boxplot NES (p < 0.05)
############################################

gsea_sig_NES <- gsea_all %>%
  filter(pval < 0.05) %>%
  mutate(
    pathway = factor(pathway, levels = c("LSD_up", "LSD_down"))
  )

# Boxplot NES (solo significativi)
p_gsea_box_NES_sig <- gsea_sig_NES %>%
  ggplot(aes(x = pathway, y = NES, fill = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(alpha = 0.6, outlier.size = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  facet_wrap(~ lsd_ref) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "NES (only p < 0.05)",
    title = "NES distribution for LSD_up / LSD_down"
  )


wilcox_gsea_NES <- gsea_sig_NES %>%
  group_by(lsd_ref) %>%
  summarise(
    p_wilcox = wilcox.test(
      NES ~ pathway,
      data = cur_data(),
      alternative = "greater"   # H1: LSD_up > LSD_down
    )$p.value,
    y_max = max(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("Wilcoxon (LSD_up > LSD_down) p = ", signif(p_wilcox, 3)),
    x_pos = 1.5,
    y_pos = y_max * 1.05
  )

p_gsea_box_NES_sig_w <- p_gsea_box_NES_sig +
  geom_text(
    data = wilcox_gsea_NES,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    vjust = 0
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_NES_sigOnly_boxplot_wilcoxon_greater.png"),
  p_gsea_box_NES_sig_w,
  width = 6, height = 4.7, dpi = 300
)

ggsave(
  file.path(results_dir, "figures",
            "LSD_vs_stress_GSEA_NES_sigOnly_boxplot_wilcoxon_greater.pdf"),
  p_gsea_box_NES_sig_w,
  width = 6, height = 4.7, dpi = 300
)
