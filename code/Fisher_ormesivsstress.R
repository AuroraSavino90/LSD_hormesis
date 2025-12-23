library(dplyr)
library(ggplot2)
library(metap)
library(stringr)

load("results/20251124/DE_24h.RData")
load("results/20251124/DE_90min.RData")
load("results/20251124/External_datasets.RData")

## General params & paths
date        <- "20251124"
results_dir <- file.path("results", date)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

get_fc_p <- function(df, nm = NULL) {
  cn     <- colnames(df)
  cn_low <- tolower(cn)
  
  fc_idx <- grep("log2foldchange|logfc", cn_low)
  fc_col <- cn[fc_idx[1]]
  
  p_idx <- grep("^p$|^pvalue$|^p\\.value$|^p\\.val$|^p_value$", cn_low)
   p_col <- cn[p_idx[1]]
  
  fc_raw <- df[[fc_col]]
  p_raw  <- df[[p_col]]
  
  fc <- suppressWarnings(as.numeric(as.character(fc_raw)))
  p  <- suppressWarnings(as.numeric(as.character(p_raw)))
  
  if (length(fc) != nrow(df) || length(p) != nrow(df)) {
    stop("Length of logFC or p-value is not consistent with nrow(df)",
         if (!is.null(nm)) paste0(" (", nm, ")") else "")
  }
  
  tibble(
    gene  = rownames(df),
    logFC = fc,
    pval  = p
  )
}

get_sig_genes <- function(df, nm = NULL, p_thresh = 0.05, fc_thresh = 0) {
  x <- get_fc_p(df, nm = nm)
  keep <- which(
    !is.na(x$logFC) & !is.na(x$pval) &
      x$pval < p_thresh &
      abs(x$logFC) > fc_thresh
  )
  x$gene[keep]
}

fisher_enrich <- function(A, B, bg) {
  A <- intersect(A, bg); B <- intersect(B, bg)
  a <- length(intersect(A, B))
  b <- length(setdiff(A, B))
  c <- length(setdiff(B, A))
  d <- length(bg) - a - b - c
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft  <- suppressWarnings(fisher.test(mat, alternative = "greater"))
  list(
    a = a, b = b, c = c, d = d,
    N = length(bg),
    OR = unname(ft$estimate),
    p  = ft$p.value
  )
}


all_obj <- ls(pattern = "^(DE|DEl)_")

hormesis_objs <- all_obj[grepl("fasting|exer|prec|hyp", all_obj, ignore.case = TRUE)]
stress_objs   <- all_obj[grepl("stress",               all_obj, ignore.case = TRUE)]


fisher_rows <- list()
err_list <- c()

for (h in hormesis_objs) {
  for (s in stress_objs) {
    
    ds1 <- get(h)
    ds2 <- get(s)
    
    bg <- intersect(rownames(ds1), rownames(ds2))
    if (length(bg) < 50) {
      message("Skip ", h, " vs ", s, " (background < 50, N = ", length(bg), ")")
      next
    }
    
    res_pair <- tryCatch({
      A <- get_sig_genes(ds1, nm = h)
      B <- get_sig_genes(ds2, nm = s)
      
      stats <- fisher_enrich(A, B, bg)
      
      data.frame(
        dataset_hormesis = h,
        dataset_stress   = s,
        overlap          = stats$a,
        only_hormesis    = stats$b,
        only_stress      = stats$c,
        neither          = stats$d,
        background       = stats$N,
        odds_ratio       = stats$OR,
        p_value          = stats$p,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Errore in ", h, " vs ", s, ": ", e$message)
      err_list <<- c(err_list, paste(h, s, ":", e$message))
      NULL
    })
    
    if (!is.null(res_pair)) {
      key <- paste(h, s, sep = "_vs_")
      fisher_rows[[key]] <- res_pair
    }
  }
}

fisher_hormesis_stress <- bind_rows(fisher_rows)
fisher_hormesis_stress$padj_BH <- p.adjust(fisher_hormesis_stress$p_value, method = "BH")
agg_hs <- fisher_hormesis_stress %>%
  summarise(
    k      = sum(!is.na(p_value)),
    meta   = list(sumlog(p_value)),  # H1: OR>1
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    chi2   = meta$chisq,
    df     = meta$df,
    p_meta = meta$p
  ) %>%
  ungroup() %>%
  select(k, chi2, df, p_meta)

out_fisher_file <- file.path(results_dir,
                             "Hormesis_vs_stress_Fisher_DEoverlap_per_pair.csv")
write.csv(fisher_hormesis_stress, out_fisher_file, row.names = FALSE)

out_agg_file <- file.path(results_dir,
                          "Hormesis_vs_stress_Fisher_DEoverlap_meta.csv")
write.csv(agg_hs, out_agg_file, row.names = FALSE)




########################################
##  Heatmap hormesis vs stress
########################################


if (!"padj_BH" %in% colnames(fisher_hormesis_stress)) {
  fisher_hormesis_stress <- fisher_hormesis_stress %>%
    mutate(padj_BH = p.adjust(p_value, method = "BH"))
}

fisher_plot_hs <- fisher_hormesis_stress %>%
  mutate(
    # Clean dataset names for plotting
    dataset_hormesis_clean = dataset_hormesis %>%
      gsub("^DEl?_+", "", .) %>%    # remove DE_ / DEl_
      gsub("_+", " ", .),           # underscores -> spaces
    dataset_stress_clean   = dataset_stress %>%
      gsub("^DEl?_+", "", .) %>%
      gsub("_+", " ", .),
    log2OR     = log2(pmax(odds_ratio, .Machine$double.eps)),
    nlp        = -log10(pmax(padj_BH, .Machine$double.xmin)),
    nlp_capped = pmin(nlp, 10)
  )

# Ordering: rank rows/cols by max |log2OR| among significant pairs
row_order <- fisher_plot_hs %>%
  group_by(dataset_hormesis_clean) %>%
  summarise(
    rank_row = max(abs(log2OR) * (padj_BH < 0.05), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(desc(rank_row))

col_order <- fisher_plot_hs %>%
  group_by(dataset_stress_clean) %>%
  summarise(
    rank_col = max(abs(log2OR) * (padj_BH < 0.05), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(desc(rank_col))

fisher_plot_hs <- fisher_plot_hs %>%
  mutate(
    dataset_hormesis_clean = factor(
      dataset_hormesis_clean,
      levels = row_order$dataset_hormesis_clean
    ),
    dataset_stress_clean = factor(
      dataset_stress_clean,
      levels = col_order$dataset_stress_clean
    )
  )

p_hs <- ggplot(
  fisher_plot_hs,
  aes(x = dataset_stress_clean, y = dataset_hormesis_clean)
) +
  geom_tile(aes(fill = log2OR), width = 0.95, height = 0.95) +
  geom_point(
    aes(size = nlp_capped),
    shape  = 21,
    color  = "black",
    fill   = "white",
    stroke = 0.25
  ) +
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
    range = c(0.8, 4.5)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Chronic stress datasets",
    y = "Hormesis datasets",
    title    = "Pairwise overlap between hormetic and chronic stress signatures",
    subtitle = "Tile: log2(odds ratio) of DE gene overlap · Dot: –log10(FDR), H1: OR > 1"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(hjust = 1),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.title     = element_text(face = "bold")
  )

# Height proportional to number of hormesis datasets
n_rows_h <- length(levels(fisher_plot_hs$dataset_hormesis_clean))
n_cols_s <- length(levels(fisher_plot_hs$dataset_stress_clean))

h_hs <- 0.28 * n_rows_h + 1.8
w_hs <- max(4.5, 0.35 * n_cols_s + 2)

fig_png_hs <- file.path(
  results_dir, "figures",
  "Hormesis_vs_chronic_stress_Fisher_DEoverlap_matrix.png"
)
fig_pdf_hs <- file.path(
  results_dir, "figures",
  "Hormesis_vs_chronic_stress_Fisher_DEoverlap_matrix.pdf"
)

ggsave(fig_png_hs, p_hs, width = w_hs, height = max(4.5, h_hs), dpi = 600)
ggsave(fig_pdf_hs, p_hs, width = w_hs, height = max(4.5, h_hs))


