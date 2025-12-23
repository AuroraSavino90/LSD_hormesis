library(dplyr)
library(ggplot2)
library(metap)
library(tidyr)

load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")

## General params & paths
date <- "20251124"
results_dir <- file.path("results", date)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)


all_obj <- ls(pattern = "^(DE|DEl)_")

refs <- c("DE_90min", "DE_24h")
non_ref <- setdiff(all_obj, refs)

fasting_targets   <- non_ref[grepl("fasting", non_ref, ignore.case = TRUE)]
exercise_targets  <- non_ref[grepl("exer",    non_ref, ignore.case = TRUE)]
precond_targets   <- non_ref[grepl("prec", non_ref, ignore.case = TRUE) |
                               grepl("hyp",  non_ref, ignore.case = TRUE)]
stress_targets    <- non_ref[grepl("stress",  non_ref, ignore.case = TRUE)]
hormesis_targets  <- unique(c(fasting_targets, exercise_targets, precond_targets))

cat("N fasting: ", length(fasting_targets),  "\n")
cat("N exercise:", length(exercise_targets), "\n")
cat("N precond: ", length(precond_targets),  "\n")
cat("N stress:  ", length(stress_targets),   "\n")
cat("N hormesis:", length(hormesis_targets), "\n")

###############################
## Directional Fisher 
###############################

run_fisher_block <- function(ref_names,
                             target_pattern,
                             block_label,
                             results_prefix,
                             results_dir,
                             env = .GlobalEnv) {
  message("=== Running block: ", block_label, " ===")
  
   all_obj <- ls(envir = env, pattern = "^(DE|DEl)_")
  message("Found DE objects in env: ", paste(all_obj, collapse = ", "))
  
  refs_found <- intersect(ref_names, all_obj)
  message("Using references: ", paste(refs_found, collapse = ", "))
  
  non_ref <- setdiff(all_obj, refs_found)
  targets <- non_ref[grepl(target_pattern, non_ref, ignore.case = TRUE)]
   
   .col_by <- function(df, pat) {
    i <- grep(pat, tolower(colnames(df)))
    if (length(i) == 0) return(NA_character_)
    colnames(df)[i[1]]
  }
  
  get_fc_p <- function(df) {
    stopifnot(!is.null(rownames(df)))
    fc_col <- .col_by(df, "log2foldchange|logfc")
    p_col  <- .col_by(df, "^(pvalue|p\\.val|p\\.value|p)$|[^a-z]pvalue[^a-z]|[^a-z]p\\.value[^a-z]")
    if (is.na(fc_col)) stop("No fold-change column (log2FoldChange/logFC) in one object.")
    if (is.na(p_col))  stop("No NOMINAL p-value column (pvalue/P.Value/p) in one object.")
    out <- df[, c(fc_col, p_col), drop = FALSE]
    colnames(out) <- c("logFC","pval")
    out
  }
  
  get_sig_dir <- function(df, p_thresh = 0.05, fc_thresh = 0) {
    x <- get_fc_p(df)
    up   <- rownames(df)[which(x$logFC >  fc_thresh & x$pval < p_thresh)]
    down <- rownames(df)[which(x$logFC < -fc_thresh & x$pval < p_thresh)]
    list(
      up   = unique(up),
      down = unique(down)
    )
  }
  
  fisher_enrich <- function(A, B, bg) {
    A <- intersect(A, bg); B <- intersect(B, bg)
    a <- length(intersect(A, B))
    b <- length(setdiff(A, B))
    c <- length(setdiff(B, A))
    d <- length(bg) - a - b - c
    mat <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
    ft <- suppressWarnings(fisher.test(mat, alternative = "greater"))
    list(a = a, b = b, c = c, d = d,
         N = length(bg),
         OR = unname(ft$estimate),
         p  = ft$p.value)
  }
  
  fisher_rows <- list()
  
  for (ref_name in refs_found) {
    
    ref_df  <- get(ref_name, envir = env)
    ref_sig <- get_sig_dir(ref_df)
    
    for (tgt_name in targets) {
      tgt_df  <- get(tgt_name, envir = env)
      tgt_sig <- get_sig_dir(tgt_df)
      
      bg <- intersect(rownames(ref_df), rownames(tgt_df))
      if (length(bg) < 50) {
        message("Skip ", ref_name, " vs ", tgt_name,
                " (background < 50, N=", length(bg), ")")
        next
      }
      
       pairs_dir <- list(
        up_vs_up      = list(A = ref_sig$up,   B = tgt_sig$up),
        down_vs_down  = list(A = ref_sig$down, B = tgt_sig$down),
        up_vs_down    = list(A = ref_sig$up,   B = tgt_sig$down),
        down_vs_up    = list(A = ref_sig$down, B = tgt_sig$up)
      )
      
      for (direction in names(pairs_dir)) {
        stats_overlap <- fisher_enrich(pairs_dir[[direction]]$A,
                                       pairs_dir[[direction]]$B,
                                       bg)
        
        fisher_rows[[paste(ref_name, tgt_name, block_label, direction, sep = "_")]] <-
          data.frame(
            block      = block_label,
            reference  = ref_name,
            dataset    = tgt_name,
            direction  = direction,
            test       = "DE_overlap",
            overlap    = stats_overlap$a,
            only_ref   = stats_overlap$b,
            only_tgt   = stats_overlap$c,
            neither    = stats_overlap$d,
            background = stats_overlap$N,
            odds_ratio = stats_overlap$OR,
            p_value    = stats_overlap$p,
            stringsAsFactors = FALSE
          )
      }
    }
  }
  
  fisher_summary <- dplyr::bind_rows(fisher_rows)
  fisher_summary$padj_BH <- p.adjust(fisher_summary$p_value, method = "BH")
  
  out_file <- file.path(
    results_dir,
    paste0(results_prefix, "_", block_label, "_Fisher_DEonly_directional.csv")
  )
  write.csv(fisher_summary, out_file, row.names = FALSE)
  
  
  #############################################
  ## plots
  #############################################
  
  
   fisher_plot_df_DE <- fisher_summary %>%
    dplyr::mutate(
      dataset_clean = gsub("^DEl?_+", "", dataset) %>% gsub("_+", " ", .),
      log2OR        = log2(pmax(odds_ratio, .Machine$double.eps)),
      nlp           = -log10(pmax(padj_BH, .Machine$double.xmin)),
      nlp_capped    = pmin(nlp, 10),
      direction     = factor(direction,
                             levels = c("up_vs_up", "down_vs_down",
                                        "up_vs_down", "down_vs_up"))
    )
  
  ordering_DE <- fisher_plot_df_DE %>%
    dplyr::group_by(reference, dataset_clean) %>%
    dplyr::summarise(
      rank_key = max(abs(log2OR) * (padj_BH < 0.05), na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(reference, dplyr::desc(rank_key))
  
  fisher_plot_df_DE <- fisher_plot_df_DE %>%
    dplyr::left_join(ordering_DE, by = c("reference","dataset_clean")) %>%
    dplyr::group_by(reference) %>%
    dplyr::mutate(
      dataset_clean = factor(
        dataset_clean,
        levels = unique(dataset_clean[order(-rank_key, dataset_clean)])
      )
    ) %>%
    dplyr::ungroup()
  
  p_title <- paste0("LSD vs ", block_label, " – directional DE gene overlap")
  
  p_fisher_DE <- ggplot(fisher_plot_df_DE,
                        aes(x = direction, y = dataset_clean)) +
    geom_tile(aes(fill = log2OR), width = 0.98, height = 0.98) +
    geom_point(aes(size = nlp_capped),
               shape = 21, color = "black", fill = "white", stroke = 0.2) +
    scale_fill_gradient2(
      name   = "log2(OR)",
      low    = "#3B4CC0",
      mid    = "white",
      high   = "#B40426",
      limits = c(-4, 4),
      oob    = scales::squish
    ) +
    scale_size_continuous(
      name  = expression(-log[10]("FDR")),
      range = c(0.7, 5)
    ) +
    labs(
      x = NULL, y = NULL,
      title    = p_title,
      subtitle = "Tile: log2(OR) · Dot: –log10(FDR), directional overlaps"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(face = "bold")
    ) +
    facet_wrap(~ reference, ncol = 2, scales = "free_y")
  
  h_fisher_DE <- 0.35 * length(unique(fisher_plot_df_DE$dataset_clean)) + 2
  
  fig_png <- file.path(
    results_dir, "figures",
    paste0(results_prefix, "_", block_label, "_Fisher_DEoverlap_matrix_directional.png")
  )
  fig_pdf <- file.path(
    results_dir, "figures",
    paste0(results_prefix, "_", block_label, "_Fisher_DEoverlap_matrix_directional.pdf")
  )
  
  ggsave(fig_png, p_fisher_DE,
         width = 7.5, height = max(4.5, h_fisher_DE), dpi = 300)
  ggsave(fig_pdf, p_fisher_DE,
         width = 7.5, height = max(4.5, h_fisher_DE))
  
  
  ##########################################################################
  ##  Aggregated meta-Fisher per reference and direction
  ################################################################################
  
  
  agg_DE <- fisher_summary %>%
    dplyr::group_by(reference, direction) %>%
    dplyr::summarise(
      k          = sum(!is.na(p_value)),
      meta       = list(metap::sumlog(p_value)),  # H1: OR>1
      mean_logOR = mean(log(pmax(odds_ratio, .Machine$double.eps)), na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      chi2   = meta$chisq,
      df     = meta$df,
      p_meta = meta$p
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(reference, direction, k, p_meta, mean_logOR, chi2, df)
  
  agg_file <- file.path(
    results_dir,
    paste0(results_prefix, "_", block_label, "_Aggregated_Fisher_DEoverlap_by_reference_directional.csv")
  )
  write.csv(agg_DE, agg_file, row.names = FALSE)
 
  
  invisible(list(
    fisher_summary = fisher_summary,
    agg_DE         = agg_DE,
    plot           = p_fisher_DE
  ))
}

###############################
## Execute functions
###############################

refs <- c("DE_90min", "DE_24h")

# Fasting/CR
res_fast <- run_fisher_block(
  ref_names      = refs,
  target_pattern = "fasting",
  block_label    = "fasting",
  results_prefix = "LSD_vs_fasting",
  results_dir    = results_dir
)

# Exercise
res_exer <- run_fisher_block(
  ref_names      = refs,
  target_pattern = "exer",
  block_label    = "exercise",
  results_prefix = "LSD_vs_exercise",
  results_dir    = results_dir
)

# Preconditioning
res_prec <- run_fisher_block(
  ref_names      = refs,
  target_pattern = "prec|hyp",
  block_label    = "preconditioning",
  results_prefix = "LSD_vs_prec",
  results_dir    = results_dir
)

# Chronic stress
res_stress <- run_fisher_block(
  ref_names      = refs,
  target_pattern = "stress",
  block_label    = "chronic_stress",
  results_prefix = "LSD_vs_stress",
  results_dir    = results_dir
)

# Hormesis (fasting + exercise + preconditioning)
res_hormesis <- run_fisher_block(
  ref_names      = refs,
  target_pattern = "fasting|exer|prec|hyp",
  block_label    = "hormesis",
  results_prefix = "LSD_vs_hormesis",
  results_dir    = results_dir
)





###############################
## Combine results
###############################


fisher_all <- bind_rows(
  res_fast$fisher_summary,
  res_exer$fisher_summary,
  res_prec$fisher_summary,
  res_stress$fisher_summary,
  res_hormesis$fisher_summary
)


comp_all <- fisher_all %>%
  filter(block %in% c("hormesis", "chronic_stress")) %>%
  mutate(
    log2OR = log2(pmax(odds_ratio, .Machine$double.eps)),
    block  = factor(block, levels = c("hormesis", "chronic_stress")),
    direction = factor(
      direction,
      levels = c("up_vs_up", "down_vs_down", "up_vs_down", "down_vs_up")
    )
  )

comp_all_clean <- comp_all %>%
  group_by(reference, direction, block) %>%
  mutate(
    Q1  = quantile(log2OR, 0.25, na.rm = TRUE),
    Q3  = quantile(log2OR, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(log2OR >= lower_bound & log2OR <= upper_bound) %>%
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound) %>%
  ungroup()

comp_all_clean <- comp_all_clean %>%
  mutate(
    dir_group = case_when(
      direction %in% c("up_vs_up", "down_vs_down") ~ "same_direction",
      direction %in% c("up_vs_down", "down_vs_up") ~ "opposite_direction",
      TRUE                                         ~ NA_character_
    ),
    dir_group = factor(
      dir_group,
      levels = c("same_direction", "opposite_direction"),
      labels = c("same direction", "opposite direction")
    )
  ) %>%
  filter(!is.na(dir_group))

############################################
##  FIGURE 1
############################################

summary_block_group <- comp_all_clean %>%
  group_by(block, reference, dir_group) %>%
  summarise(
    median_log2OR = median(log2OR, na.rm = TRUE),
    q1            = quantile(log2OR, 0.25, na.rm = TRUE),
    q3            = quantile(log2OR, 0.75, na.rm = TRUE),
    n_datasets    = n(),
    prop_sig      = mean(padj_BH < 0.05 & odds_ratio > 1, na.rm = TRUE),
    .groups       = "drop"
  )

## Wilcoxon: hormesis > chronic_stress 
wilcox_group_ref <- comp_all_clean %>%
  group_by(reference, dir_group) %>%
  summarise(
    p_wilcox = wilcox.test(
      log2OR ~ block,
      data = cur_data(),
      alternative = "greater"  # H1: hormesis > chronic_stress
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p_wilcox, 3)))

y_annot <- summary_block_group %>%
  group_by(reference) %>%
  summarise(y = max(q3, na.rm = TRUE) + 0.4, .groups = "drop")

wilcox_labels <- wilcox_group_ref %>%
  left_join(y_annot, by = "reference")

position_dodge2 <- position_dodge(width = 0.6)


p_summary_groups <- ggplot() +
  geom_col(
    data = summary_block_group,
    aes(x = dir_group, y = median_log2OR, fill = block),
    position = position_dodge2,
    width = 0.5,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_block_group,
    aes(x = dir_group, ymin = q1, ymax = q3, group = block),
    position = position_dodge2,
    width = 0.2
  ) +
  geom_jitter(
    data = comp_all_clean,
    aes(x = dir_group, y = log2OR, color = block),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.6
    ),
    size = 1.4,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ reference, nrow = 1) +
  geom_text(
    data = wilcox_labels,
    aes(x = dir_group, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0,
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = "median log2(odds ratio) · DE overlap",
    title = "Directional enrichment of LSD DEGs in hormesis vs chronic stress DEGs",
    subtitle = paste(
      "Same direction = (up↔up + down↔down);",
      "Opposite direction = (up↔down + down↔up).",
      "Bars: median (Q1–Q3); points: datasets; Wilcoxon one-sided (hormesis > chronic_stress)."
    )
  ) +
  scale_fill_manual(values = c("hormesis" = "#E69F00", "chronic_stress" = "#56B4E9")) +
  scale_color_manual(values = c("hormesis" = "#E69F00", "chronic_stress" = "#56B4E9")) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "top",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 20, hjust = 0.9)
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_hormesis_vs_chronicStress_same_vs_opposite_direction_outliersFiltered.png"),
  p_summary_groups,
  width = 7.5, height = 6, dpi = 300
)

ggsave(
  file.path(results_dir, "figures",
            "LSD_hormesis_vs_chronicStress_same_vs_opposite_direction_outliersFiltered.pdf"),
  p_summary_groups,
  width = 7.5, height = 6
)

############################################
## FIGURE 2: p_summary_directions
##   
############################################



summary_block_dir <- comp_all_clean %>%
  group_by(block, reference, direction) %>%
  summarise(
    median_log2OR = median(log2OR, na.rm = TRUE),
    q1            = quantile(log2OR, 0.25, na.rm = TRUE),
    q3            = quantile(log2OR, 0.75, na.rm = TRUE),
    n_datasets    = n(),
    prop_sig      = mean(padj_BH < 0.05 & odds_ratio > 1, na.rm = TRUE),
    .groups       = "drop"
  )

## Wilcoxon: hormesis > chronic_stress 

wilcox_dir_ref <- comp_all_clean %>%
  group_by(reference, direction) %>%
  summarise(
    p_wilcox = wilcox.test(
      log2OR ~ block,
      data = cur_data(),
      alternative = "greater"   # H1: hormesis > chronic_stress
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p_wilcox, 3)))

y_annot_dir <- summary_block_dir %>%
  group_by(reference) %>%
  summarise(y = max(q3, na.rm = TRUE) + 0.4, .groups = "drop")

wilcox_labels_dir <- wilcox_dir_ref %>%
  left_join(y_annot_dir, by = "reference")

pos_dodge <- position_dodge(width = 0.6)


p_summary_directions <- ggplot() +
  geom_col(
    data = summary_block_dir,
    aes(x = direction, y = median_log2OR, fill = block),
    position = pos_dodge,
    width = 0.5,
    color = "black"
  ) +
  geom_errorbar(
    data = summary_block_dir,
    aes(x = direction, ymin = q1, ymax = q3, group = block),
    position = pos_dodge,
    width = 0.2
  ) +
  geom_jitter(
    data = comp_all_clean,
    aes(x = direction, y = log2OR, color = block),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.6
    ),
    size = 1.2,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ reference, nrow = 1) +
  geom_text(
    data = wilcox_labels_dir,
    aes(x = direction, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0,
    fontface = "bold"
  ) +
  labs(
    x = NULL,
    y = "median log2(odds ratio) · DE overlap",
    title = "Directional enrichment of LSD DEGs (4 directions)",
    subtitle = paste(
      "up_vs_up, down_vs_down, up_vs_down, down_vs_up;",
      "bar = median (Q1–Q3); punti = singoli dataset;",
      "Wilcoxon one-sided (hormesis > chronic_stress) per direzione."
    )
  ) +
  scale_fill_manual(values = c("hormesis" = "#E69F00", "chronic_stress" = "#56B4E9")) +
  scale_color_manual(values = c("hormesis" = "#E69F00", "chronic_stress" = "#56B4E9")) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "top",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 25, hjust = 0.9)
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_hormesis_vs_chronicStress_4directions_log2OR_outliersFiltered.png"),
  p_summary_directions,
  width = 10, height = 5, dpi = 300
)
ggsave(
  file.path(results_dir, "figures",
            "LSD_hormesis_vs_chronicStress_4directions_log2OR_outliersFiltered.pdf"),
  p_summary_directions,
  width = 10, height = 5
)

############################################
## FIGURE 3: p_log2OR_dir
##   
############################################


summary_dir <- comp_all_clean %>%
  group_by(block, reference, direction) %>%
  summarise(
    median_log2OR = median(log2OR, na.rm = TRUE),
    q1            = quantile(log2OR, 0.25, na.rm = TRUE),
    q3            = quantile(log2OR, 0.75, na.rm = TRUE),
    n_datasets    = n(),
    .groups       = "drop"
  )

## Test: up_vs_up / down_vs_down > (up_vs_down + down_vs_up)

wilcox_dir_same_vs_opp <- comp_all_clean %>%
  filter(direction %in% c("up_vs_up", "down_vs_down",
                          "up_vs_down", "down_vs_up")) %>%
  group_by(block, reference) %>%
  summarise(
    p_up_vs_opp = if (
      sum(direction == "up_vs_up", na.rm = TRUE) > 0 &&
      sum(direction %in% c("up_vs_down", "down_vs_up"), na.rm = TRUE) > 0
    ) {
      wilcox.test(
        log2OR[direction == "up_vs_up"],
        log2OR[direction %in% c("up_vs_down", "down_vs_up")],
        alternative = "greater"          # H1: up_vs_up > opposite
      )$p.value
    } else {
      NA_real_
    },
    p_down_vs_opp = if (
      sum(direction == "down_vs_down", na.rm = TRUE) > 0 &&
      sum(direction %in% c("up_vs_down", "down_vs_up"), na.rm = TRUE) > 0
    ) {
      wilcox.test(
        log2OR[direction == "down_vs_down"],
        log2OR[direction %in% c("up_vs_down", "down_vs_up")],
        alternative = "greater"          # H1: down_vs_down > opposite
      )$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

y_pos_dir <- summary_dir %>%
  filter(direction %in% c("up_vs_up", "down_vs_down")) %>%
  group_by(block, reference, direction) %>%
  summarise(
    y = max(q3, na.rm = TRUE) + 0.3,
    .groups = "drop"
  )

wilcox_labels_dir2 <- wilcox_dir_same_vs_opp %>%
  pivot_longer(
    cols      = c(p_up_vs_opp, p_down_vs_opp),
    names_to  = "comparison",
    values_to = "p_value"
  ) %>%
  mutate(
    direction = dplyr::case_when(
      comparison == "p_up_vs_opp"   ~ "up_vs_up",
      comparison == "p_down_vs_opp" ~ "down_vs_down",
      TRUE                          ~ NA_character_
    ),
    label = ifelse(
      is.na(p_value),
      "n/a",
      paste0("p = ", signif(p_value, 3))
    )
  ) %>%
  filter(!is.na(direction)) %>%
  left_join(y_pos_dir, by = c("block", "reference", "direction"))


p_log2OR_dir <- ggplot(summary_dir,
                       aes(x = direction,
                           y = median_log2OR,
                           fill = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = q1, ymax = q3),
                width = 0.2) +
  geom_text(
    data = wilcox_labels_dir2,
    aes(x = direction, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0,
    fontface = "bold"
  ) +
  facet_grid(block ~ reference) +
  labs(
    x = NULL,
    y = "median log2(odds ratio) · DE overlap",
    title = "Directional enrichment of LSD DEGs",
    subtitle = paste(
      "up_vs_up e down_vs_down testati vs (up_vs_down + down_vs_up)",
      "con Wilcoxon one-sided (same > opposite)."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 25, hjust = 0.9)
  )

ggsave(
  file.path(results_dir, "figures",
            "LSD_4directions_log2OR_by_block_reference_withWilcoxon.png"),
  p_log2OR_dir,
  width = 8, height = 7, dpi = 300
)

ggsave(
  file.path(results_dir, "figures",
            "LSD_4directions_log2OR_by_block_reference_withWilcoxon.pdf"),
  p_log2OR_dir,
  width = 8, height = 7
)
