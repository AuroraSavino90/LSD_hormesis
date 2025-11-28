library(dplyr)
library(ggplot2)

load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")

## General params & paths
date        <- "20251124"
results_dir <- file.path("results", date)

######################
## remove non-cortex datasets
##########################

rm(DEl_exer_GSE39697_4d, DEl_exer_GSE39697_30d, DEl_exer_GSE39697_30d_gc, DEl_exer_GSE39697_4d_gc,
   DEl_stress_GSE151807_hipp, DEl_stress_GSE151807_amy, DEl_stress_GSE109315_vHPC_C57BL,
   DEl_stress_GSE109315_vHPC_DBA, DEl_stress_GSE109315_BNST_C57BL,
   DEl_stress_GSE109315_BNST_DBA, DEl_stress_GSE84183_DG, DEl_stress_GSE84183_CG,
   DEl_stress_GSE8262_hippx_CCH, DEl_stress_GSE8262_hippx_CIH,
   DE_stress_GSE198725_BLA, DE_stress_GSE198725_VTA )


###############################
## define main objects
###############################

all_obj <- ls(pattern = "^(DE|DEl)_")

# Riferimenti LSD (devono esistere)
refs <- c("DE_90min", "DE_24h")
stopifnot(all(refs %in% all_obj))

# Escludi i riferimenti dal resto
non_ref <- setdiff(all_obj, refs)

# Gruppi di dataset in base al nome dell'oggetto
fasting_targets <- non_ref[grepl("fasting", non_ref, ignore.case = TRUE)]
exercise_targets <- non_ref[grepl("exer", non_ref, ignore.case = TRUE)]
# Preconditioning: ischemic precond + mild hypoxia "hormetic"
precond_targets <- non_ref[grepl("prec", non_ref, ignore.case = TRUE) |
                             grepl("hyp",  non_ref, ignore.case = TRUE)]
stress_targets <- non_ref[grepl("stress", non_ref, ignore.case = TRUE)]

# "hormesis" 
hormesis_targets <- unique(c(fasting_targets, exercise_targets, precond_targets))


cat("N fasting: ", length(fasting_targets), "\n")
cat("N exercise:", length(exercise_targets), "\n")
cat("N precond: ", length(precond_targets), "\n")
cat("N stress:  ", length(stress_targets), "\n")
cat("N hormesis:", length(hormesis_targets), "\n")


run_fisher_block <- function(ref_names,
                             target_pattern,
                             block_label,
                             results_prefix,
                             results_dir,
                             env = .GlobalEnv) {
  message("=== Running block: ", block_label, " ===")
  
  ## 1) Trova tutti gli oggetti DE/DEl nell'ambiente specificato
  all_obj <- ls(envir = env, pattern = "^(DE|DEl)_")
  message("Found DE objects in env: ", paste(all_obj, collapse = ", "))
  
  ## 2) Reference effettivamente presenti
  refs_found <- intersect(ref_names, all_obj)
  if (length(refs_found) == 0) {
    stop("Nessun oggetto di riferimento trovato tra: ",
         paste(ref_names, collapse = ", "),
         "\nOggetti DE disponibili: ",
         paste(all_obj, collapse = ", "))
  }
  if (length(refs_found) < length(ref_names)) {
    missing <- setdiff(ref_names, refs_found)
    warning("I seguenti reference non sono stati trovati e verranno ignorati: ",
            paste(missing, collapse = ", "))
  }
  message("Using references: ", paste(refs_found, collapse = ", "))
  
  ## 3) Target: tutti gli altri oggetti che matchano il pattern
  non_ref <- setdiff(all_obj, refs_found)
  targets <- non_ref[grepl(target_pattern, non_ref, ignore.case = TRUE)]
  if (length(targets) == 0) {
    warning("Nessun target trovato per pattern '", target_pattern,
            "' (block: ", block_label, ").")
    return(invisible(NULL))
  }
  message("Targets (", block_label, "): ", paste(targets, collapse = ", "))
  
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
  
  get_sig <- function(df, p_thresh = 0.05, fc_thresh = 0) {
    x <- get_fc_p(df)
    sig <- rownames(df)[which(abs(x$logFC) > fc_thresh & x$pval < p_thresh)]
    unique(sig)
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
    if (!exists(ref_name, envir = env)) {
      warning("Reference '", ref_name, "' non trovato in env; salto.")
      next
    }
    
    ref_df  <- get(ref_name, envir = env)
    ref_sig <- get_sig(ref_df)
    
    for (tgt_name in targets) {
      if (!exists(tgt_name, envir = env)) {
        warning("Target '", tgt_name, "' non trovato in env; salto.")
        next
      }
      
      tgt_df  <- get(tgt_name, envir = env)
      tgt_sig <- get_sig(tgt_df)
      
      bg <- intersect(rownames(ref_df), rownames(tgt_df))
      if (length(bg) < 50) {
        message("Skip ", ref_name, " vs ", tgt_name,
                " (background < 50, N=", length(bg), ")")
        next
      }
      
      stats_overlap <- fisher_enrich(ref_sig, tgt_sig, bg)
      
      fisher_rows[[paste(ref_name, tgt_name, block_label, sep = "_")]] <-
        data.frame(
          block      = block_label,
          reference  = ref_name,
          dataset    = tgt_name,
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
  
  if (length(fisher_rows) == 0) {
    warning("Nessun risultato Fisher per il block '", block_label, "'.")
    return(invisible(NULL))
  }
  
  fisher_summary <- dplyr::bind_rows(fisher_rows)
  fisher_summary$padj_BH <- p.adjust(fisher_summary$p_value, method = "BH")
 
   out_file <- file.path(
    results_dir,
    paste0(results_prefix, "_", block_label, "_Fisher_DEonly.csv")
  )
  write.csv(fisher_summary, out_file, row.names = FALSE)
  message("Salvato: ", out_file)
  
  fisher_plot_df_DE <- fisher_summary %>%
    dplyr::mutate(
      dataset_clean = gsub("^DEl?_+", "", dataset) %>% gsub("_+", " ", .),
      test_label    = "DE∩DE",
      log2OR        = log2(pmax(odds_ratio, .Machine$double.eps)),
      nlp           = -log10(pmax(padj_BH, .Machine$double.xmin)),
      nlp_capped    = pmin(nlp, 10)
    )
  
  # Global Ordering based on max |log2OR| on any reference
  ordering_DE <- fisher_plot_df_DE %>%
    dplyr::group_by(dataset_clean) %>%
    dplyr::summarise(
      rank_key = max(abs(log2OR) * (padj_BH < 0.05), na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(rank_key))
  
  fisher_plot_df_DE <- fisher_plot_df_DE %>%
    dplyr::left_join(ordering_DE, by = "dataset_clean") %>%
    dplyr::mutate(
      dataset_clean = factor(
        dataset_clean,
        levels = unique(dataset_clean[order(-rank_key, dataset_clean)])
      ),
      reference = factor(reference, levels = sort(unique(reference)))
    )
  
  p_title <- paste0("LSD vs ", block_label, " – DE gene overlap")
  
   p_fisher_DE <- ggplot(fisher_plot_df_DE,
                        aes(x = reference, y = dataset_clean)) +
    geom_tile(aes(fill = log2OR),
              width  = 0.99,
              height = 0.99) +
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
    scale_x_discrete(
      expand = c(0, 0),
      labels = c(
        "DE_90min" = "90min",
        "DE_24h"   = "24h"
      )
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = NULL, y = NULL,
      title    = p_title,
      subtitle = "Tile: log2(OR) for DE∩DE · Dot: –log10(FDR), H1: OR>1"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid  = element_blank(),
  
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      
      plot.title  = element_text(face = "bold", hjust = 0.5),
      plot.title.position = "plot"
    )
  
  
  n_rows <- length(unique(fisher_plot_df_DE$dataset_clean))
  h_fisher_DE <- 0.22 * n_rows + 1.5  
  
  fig_png <- file.path(
    results_dir, "figures",
    paste0(results_prefix, "_", block_label, "_Fisher_DEoverlap_matrix.png")
  )
  fig_pdf <- file.path(
    results_dir, "figures",
    paste0(results_prefix, "_", block_label, "_Fisher_DEoverlap_matrix.pdf")
  )
  
  ggsave(fig_png, p_fisher_DE,
         width = 3.8,                    
         height = max(4, h_fisher_DE),   
         dpi = 300)
  
  ggsave(fig_pdf, p_fisher_DE,
         width = 3.8,
         height = max(4, h_fisher_DE))
  
  
  ## Aggregated meta-Fisher per reference
  agg_DE <- fisher_summary %>%
    dplyr::group_by(reference) %>%
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
    dplyr::select(reference, k, p_meta, mean_logOR, chi2, df)
  
  agg_file <- file.path(
    results_dir,
    paste0(results_prefix, "_", block_label, "_Aggregated_Fisher_DEoverlap_by_reference.csv")
  )
  write.csv(agg_DE, agg_file, row.names = FALSE)
 
  
  invisible(list(
    fisher_summary = fisher_summary,
    agg_DE         = agg_DE,
    plot           = p_fisher_DE
  ))
}


###############################
## Run for each group
###############################

refs <- c("DE_90min", "DE_24h")

# Fasting/CR
res_fast <- run_fisher_block(
  ref_names     = refs,
  target_pattern = "fasting",
  block_label    = "fasting",
  results_prefix = "LSD_vs_fasting",
  results_dir    = results_dir
)

# Exercise
res_exer <- run_fisher_block(
  ref_names     = refs,
  target_pattern = "exer",
  block_label    = "exercise",
  results_prefix = "LSD_vs_exercise",
  results_dir    = results_dir
)

# Preconditioning
res_prec <- run_fisher_block(
  ref_names     = refs,
  target_pattern = "prec|hyp",
  block_label    = "preconditioning",
  results_prefix = "LSD_vs_prec",
  results_dir    = results_dir
)

# Chronic stress
res_stress <- run_fisher_block(
  ref_names     = refs,
  target_pattern = "stress",
  block_label    = "chronic_stress",
  results_prefix = "LSD_vs_stress",
  results_dir    = results_dir
)

# hormesis
res_hormesis <- run_fisher_block(
  ref_names     = refs,
  target_pattern = "fasting|exer|prec|hyp",
  block_label    = "hormesis",
  results_prefix = "LSD_vs_hormesis",
  results_dir    = results_dir
)



fisher_all <- bind_rows(
  res_fast$fisher_summary,
  res_exer$fisher_summary,
  res_prec$fisher_summary,
  res_stress$fisher_summary,
  res_hormesis$fisher_summary
)

agg_all <- bind_rows(
  res_fast$agg_DE      %>% mutate(block = "fasting"),
  res_exer$agg_DE      %>% mutate(block = "exercise"),
  res_prec$agg_DE      %>% mutate(block = "preconditioning"),
  res_stress$agg_DE    %>% mutate(block = "chronic_stress"),
  res_hormesis$agg_DE  %>% mutate(block = "hormesis")
)


comp_ds <- fisher_all %>%
  filter(block %in% c("hormesis", "chronic_stress")) %>%
  mutate(
    log2OR = log2(pmax(odds_ratio, .Machine$double.eps)),
    block  = factor(block, levels = c("hormesis", "chronic_stress"))
  )

summary_block <- comp_ds %>%
  filter(odds_ratio > 1) %>%  # solo arricchimenti "positivi"
  group_by(block, reference) %>%
  summarise(
    n_datasets  = n(),
    n_signif    = sum(padj_BH < 0.05, na.rm = TRUE),
    prop_signif = n_signif / n_datasets,
    median_log2OR = median(log2OR[padj_BH < 0.05], na.rm = TRUE),
    .groups = "drop"
  )


p_log2OR <- ggplot(comp_ds, aes(x = block, y = log2OR, color = block)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  facet_wrap(~ reference, ncol = 2) +
  labs(
    x = NULL,
    y = "log2(odds ratio) DE overlap",
    title = "Enrichment of LSD DEGs in hormesis vs chronic stress DEGs",
    subtitle = "Points: individual datasets, black dot: median"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )


wilcox_by_reference <- comp_ds %>%
  group_by(reference) %>%
  summarise(
    wilcox_p = wilcox.test(
      log2OR ~ block,
      data = cur_data(),
      alternative = "greater"   # perché ipotesi: hormesis > stress
    )$p.value
  )

wilcox_labels <- wilcox_by_reference %>%
  mutate(label = paste0("Wilcoxon p = ", signif(wilcox_p, 3)))

p_log2OR_labeled <- p_log2OR +
  geom_text(
    data = wilcox_labels,
    aes(
      x = 1.5,
      y = max(comp_ds$log2OR, na.rm = TRUE),
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1,
    size = 4,
    fontface = "bold"
  )



########################################
## Meta-Fisher: hormesis vs stress 
########################################

compare_hormesis_vs_stress <- function(fisher_all,
                                       block_hormesis = "hormesis",
                                       block_stress   = "chronic_stress",
                                       reference      = "DE_24h",
                                       B = 5000,
                                       seed = 123) {
  set.seed(seed)
  
  df_h <- fisher_all %>%
    filter(block == block_hormesis, reference == !!reference)
  
  df_s <- fisher_all %>%
    filter(block == block_stress, reference == !!reference)
  
  n_h <- n_distinct(df_h$dataset)
  n_s <- n_distinct(df_s$dataset)
  
  if (n_h < n_s) {
    warning("Fewer hormesis datasets (", n_h,
            ") than chronic-stress datasets (", n_s,
            "). Subsampling may be suboptimal.")
  }
  
  # meta-Fisher for chronic stress (all datasets)
  p_meta_stress <- metap::sumlog(df_s$p_value)$p
  
  # Bootstrap meta-Fisher for hormesis using n_s datasets each time
  boot_horm <- replicate(B, {
    ds <- sample(unique(df_h$dataset), n_s, replace = FALSE)
    ps <- df_h %>% filter(dataset %in% ds) %>% pull(p_value)
    metap::sumlog(ps)$p
  })
  
  # Permutation p-value: how often hormesis meta-p < stress meta-p
  n_less <- sum(boot_horm < p_meta_stress)
  p_perm <- n_less / length(boot_horm)
  
  list(
    boot_hormesis = boot_horm,
    p_meta_stress = p_meta_stress,
    n_hormesis    = n_h,
    n_stress      = n_s,
    perm_p        = p_perm,
    n_less        = n_less
  )
}

res_comp <- compare_hormesis_vs_stress(
  fisher_all,
  block_hormesis = "hormesis",
  block_stress   = "chronic_stress",
  reference      = "DE_24h",
  B = 5000
)

########################################
## Bootstrap density plot ##
########################################

df_boot <- tibble(
  block  = "hormesis",
  p_meta = res_comp$boot_hormesis
)

p_meta_stress <- res_comp$p_meta_stress
perm_p        <- res_comp$perm_p

p_boot <- ggplot(df_boot, aes(x = -log10(p_meta))) +
  geom_density(alpha = 0.4) +
  geom_vline(
    xintercept = -log10(p_meta_stress),
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(-log[10](p[meta])),
    y = "Density",
    title    = "Meta-analysis of LSD DEG overlap: hormesis vs chronic stress",
    subtitle = paste0(
      "Bootstrap distribution of hormesis meta-p (n = ", res_comp$n_stress,
      " datasets) vs observed chronic-stress meta-p (vertical line);\n",
      "permutation p = ", signif(perm_p, 3)
    )
  )

boot_png <- file.path(
  results_dir, "figures",
  "LSD_vs_hormesis_vs_chronic_stress_metaFisher_bootstrap.png"
)
boot_pdf <- file.path(
  results_dir, "figures",
  "LSD_vs_hormesis_vs_chronic_stress_metaFisher_bootstrap.pdf"
)

ggsave(boot_png, p_boot, width = 5.0, height = 4.0, dpi = 600)
ggsave(boot_pdf, p_boot, width = 5.0, height = 4.0)


########################################
## Running-threshold significance   ##
########################################

significance_curve <- function(df, n_points = 200) {
  thr <- seq(0, 0.2, length.out = n_points)  # extend to 1 if needed
  
  curves <- lapply(thr, function(t) {
    tibble(
      threshold = t,
      block     = df$block[1],
      reference = df$reference[1],
      prop      = mean(df$p_value < t, na.rm = TRUE)
    )
  })
  
  bind_rows(curves)
}

sig_hormesis <- fisher_all %>%
  filter(block == "hormesis", reference == "DE_24h") %>%
  significance_curve()

sig_stress <- fisher_all %>%
  filter(block == "chronic_stress", reference == "DE_24h") %>%
  significance_curve()

sig_curve_all <- bind_rows(sig_hormesis, sig_stress)

p_curve <- ggplot(sig_curve_all,
                  aes(x = threshold, y = prop, color = block)) +
  geom_line(linewidth = 1.1) +
  theme_classic(base_size = 14) +
  labs(
    x = "P-value threshold",
    y = "Fraction of datasets with p < threshold",
    title    = "Cumulative fraction of significant overlaps",
    subtitle = "LSD DEG enrichment in hormetic vs chronic stress signatures"
  ) +
  scale_color_manual(
    name   = "Block",
    values = c("hormesis" = "#1b9e77", "chronic_stress" = "#d95f02")
  )

curve_png <- file.path(
  results_dir, "figures",
  "LSD_vs_hormesis_vs_chronic_stress_running_pcurve.png"
)
curve_pdf <- file.path(
  results_dir, "figures",
  "LSD_vs_hormesis_vs_chronic_stress_running_pcurve.pdf"
)

ggsave(curve_png, p_curve, width = 5.0, height = 4.2, dpi = 600)
ggsave(curve_pdf, p_curve, width = 5.0, height = 4.2)

## permutation p-value
res_comp$perm_p
