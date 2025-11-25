#############################################
## LSD PFC – Basic DE & Enrichment Pipeline
#############################################

## Packages ----
library(DESeq2)
library(FactoMineR)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(enrichR)

theme_set(theme_classic())

## General params & paths

date        <- "20251124"
results_dir <- file.path("results", date)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

## Helper functions 

# Rename genes using annotation, merging duplicates by highest average expression
changenames <- function(data, anno){
  stopifnot(!is.null(rownames(data)))
  ids <- rownames(data)
  anno_sub <- anno[anno[,1] %in% ids, , drop = FALSE]
  sym_clean <- anno_sub[,2]; sym_clean[is.na(sym_clean)] <- ""
  anno_sub[,2] <- sym_clean
  ids_ok  <- as.character(anno_sub[, 1])
  syms_ok <- as.character(anno_sub[, 2])
  keep    <- !is.na(ids_ok) & nzchar(ids_ok)
  ids_ok  <- ids_ok[keep]
  syms_ok <- syms_ok[keep]
  
  sym_by_id <- tapply(
    syms_ok, ids_ok,
    function(x) {
      u <- unique(x[!is.na(x) & nzchar(x)])
      u
    },
    simplify = FALSE
  )
  
  # 2) Harden the resolver
  resolve_symbol <- function(id) {
    id <- as.character(id)
    if (is.na(id) || !nzchar(id)) return(NA_character_)  # invalid name
    syms <- sym_by_id[[id]]
    if (is.null(syms) || length(syms) == 0) return(NA_character_)
    if (length(syms) == 1) return(syms)
    non_ps <- syms[!grepl("-ps", syms)]
    if (length(non_ps) == 1) return(non_ps)
    NA_character_
  }
  
  # 3) Coerce row ids to character before vapply
  ids <- as.character(rownames(data))
  map_vec <- vapply(ids, resolve_symbol, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  keep <- !is.na(map_vec) & nzchar(map_vec)
  if (!any(keep)) stop("After filtering, no rows remain with a valid symbol.")
  data <- data[keep, , drop = FALSE]; map_vec <- map_vec[keep]
  rms <- rowMeans(data, na.rm = TRUE)
  groups <- split(seq_len(nrow(data)), map_vec)
  keep_idx <- unlist(lapply(groups, function(ix) ix[which.max(rms[ix])]), use.names = FALSE)
  data <- data[keep_idx, , drop = FALSE]
  rownames(data) <- map_vec[keep_idx]
  data
}

## Volcano & summary helpers ----

create_volcano <- function(res,
                           contrast_name,
                           results_dir,
                           lfc_thresh = 1,
                           alpha = 0.05) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj), ]
  
  df$negLog10Padj <- -log10(df$padj)
  df$status <- "NS"
  df$status[df$padj < alpha & df$log2FoldChange >  lfc_thresh] <- "Up"
  df$status[df$padj < alpha & df$log2FoldChange < -lfc_thresh] <- "Down"
  
  p <- ggplot(df, aes(x = log2FoldChange,
                      y = negLog10Padj,
                      colour = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = c("Down" = "blue",
                                   "NS"   = "grey70",
                                   "Up"   = "red")) +
    theme_bw() +
    labs(
      title = paste("Volcano –", contrast_name),
      x = "log2 fold change",
      y = "-log10(adj p-value)"
    ) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed")
  
  ggsave(
    file.path(results_dir, paste0("Volcano_", contrast_name, ".png")),
    p, width = 6, height = 5, dpi = 300
  )
}

summarize_de <- function(res,
                         contrast_name,
                         results_dir,
                         alpha = 0.05,
                         n = 50) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  
  # significant
  df_sig <- df[!is.na(df$padj) & df$padj < alpha, ]
  df_sig <- df_sig[order(df_sig$padj,
                         -abs(df_sig$log2FoldChange)), ]
  
  # save all significant
  write.csv(
    df_sig,
    file = file.path(results_dir,
                     paste0("DE_", contrast_name, "_all_sig.csv")),
    row.names = FALSE
  )
  
  # top up / down
  up <- df_sig[df_sig$log2FoldChange > 0, ]
  dn <- df_sig[df_sig$log2FoldChange < 0, ]
  
  write.csv(
    head(up, n),
    file = file.path(results_dir,
                     paste0("DE_", contrast_name, "_top", n, "_up.csv")),
    row.names = FALSE
  )
  write.csv(
    head(dn, n),
    file = file.path(results_dir,
                     paste0("DE_", contrast_name, "_top", n, "_down.csv")),
    row.names = FALSE
  )
}

# Simple ggplot boxplot helper for twitch counts or gene expression
save_boxplot <- function(df, x, y, fill = NULL, filename, title = NULL, ylab = NULL) {
  p <- ggplot(df, aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    geom_boxplot() +
    theme_bw() +
    labs(title = title, y = ylab, x = NULL)
  
  ggsave(filename, p, width = 5, height = 7, dpi = 300)
}


#############################################
## 1. Load & normalize counts
#############################################

counts <- read.csv("data/merged_gene_counts.txt", sep = "\t")
rownames(counts) <- counts[, 1]
anno <- counts[, 1:2]

# Rename genes and drop annotation columns
counts_name <- changenames(counts[, -c(1, 2)], anno = anno)

# RPM
RPM <- t(t(counts_name) / colSums(counts_name)) * 1e6

# Filter lowly expressed genes
keep_genes <- rowSums(counts_name >= 10) > 2
counts_name <- counts_name[keep_genes, ]
RPM         <- RPM[keep_genes, ]

# Log2 transform
RPMlog <- log2(RPM + 1)


#############################################
## 2. Metadata & PCA
#############################################

meta <- data.frame(
  treatment = c("LSD", "Saline", "Saline", "Saline", "Saline", "LSD",
                "LSD", "LSD", "LSD", "Saline", "Saline", "LSD"),
  time      = c("90min", "24h", "24h", "90min", "90min", "90min",
                "90min", "24h", "24h", "90min", "24h", "24h")
)
rownames(meta) <- colnames(counts_name)

pca <- PCA(t(RPMlog), graph = FALSE)

pca_df <- data.frame(
  PC1 = pca$ind$coord[, 1],
  PC2 = pca$ind$coord[, 2],
  PC3 = pca$ind$coord[, 3],
  PC4 = pca$ind$coord[, 4],
  PC5 = pca$ind$coord[, 5],
  meta
)

p_pca_treat <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(size = 3) +
  labs(title = "PCA – treatment")

p_pca_time <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = time)) +
  geom_point(size = 3) +
  labs(title = "PCA – time")

ggsave(file.path(results_dir, "PCA_treatment.pdf"), p_pca_treat, width = 5, height = 5)
ggsave(file.path(results_dir, "PCA_time.pdf"),      p_pca_time, width = 5, height = 5)
ggsave(file.path(results_dir, "PCA_treatment.png"), p_pca_treat, width = 5, height = 5, dpi = 300)
ggsave(file.path(results_dir, "PCA_time.png"),      p_pca_time,  width = 5, height = 5, dpi = 300)


#############################################
## 3. IEG QC at 90min 
#############################################

# Indices for 90min LSD vs Saline 
idx_90min <- c(1, 6, 7, 10, 4, 5) # first three LSD, last three Saline
group_90min <- rep(c("LSD", "Saline"), each = 3)

# Indices for 24h LSD vs Saline 
idx_24h      <- c(9, 12, 8, 11, 3, 2)
group_24h    <- rep(c("LSD", "Saline"), each = 3)


plot_ieg_box <- function(gene, filename_prefix, time) {
  png(file.path(results_dir, paste0(filename_prefix, ".png")), res = 300,
      width = 700, height = 1500)
  
  if(time=="90min"){
    boxplot(
      RPMlog[gene, idx_90min] ~ group_90min,
      main = paste(gene, time),
      ylab = "log2 RPM + 1")
  } else {
    boxplot(
      RPMlog[gene, idx_24h] ~ group_24h,
      main = paste(gene, time),
      ylab = "log2 RPM + 1")
  }
  dev.off()
}

plot_ieg_box("Fos",  "Fos90min", "90min")
plot_ieg_box("Arc",  "Arc90min", "90min")
plot_ieg_box("Bdnf", "Bdnf90min", "90min")
plot_ieg_box("Fos",  "Fos24h", "24h")
plot_ieg_box("Arc",  "Arc24h", "24h")
plot_ieg_box("Bdnf", "Bdnf24h", "24h")


#############################################
## 4. DESeq2 – Separate 90min and 24h analyses
#############################################

## 90min
dds_90min <- DESeqDataSetFromMatrix(
  countData = counts_name[,idx_90min],
  colData   = data.frame(treatment = group_90min,
                         row.names = colnames(counts_name)[idx_90min]),
  design    = ~ treatment
)
dds_90min <- DESeq(dds_90min)
DE_90min  <- results(dds_90min, contrast = c("treatment", "LSD", "Saline"))

up_90min <- rownames(DE_90min)[which(DE_90min$padj < 0.05 & DE_90min$log2FoldChange > 0)]
dn_90min <- rownames(DE_90min)[which(DE_90min$padj < 0.05 & DE_90min$log2FoldChange < 0)]

write.csv(up_90min,
          file = file.path(results_dir, "90min_up.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(dn_90min,
          file = file.path(results_dir, "90min_dn.csv"),
          quote = FALSE, row.names = FALSE)
save(DE_90min, file = file.path(results_dir, "DE_90min.RData"))

create_volcano(DE_90min, "90min", results_dir)
summarize_de(DE_90min, "90min", results_dir)

## 24h
dds_24h <- DESeqDataSetFromMatrix(
  countData = counts_name[, idx_24h],
  colData   = data.frame(treatment = group_24h,
                         row.names = colnames(counts_name)[idx_24h]),
  design    = ~ treatment
)
dds_24h <- DESeq(dds_24h)
DE_24h  <- results(dds_24h, contrast = c("treatment", "LSD", "Saline"))

up_24h <- rownames(DE_24h)[which(DE_24h$padj < 0.05 & DE_24h$log2FoldChange > 0)]
dn_24h <- rownames(DE_24h)[which(DE_24h$padj < 0.05 & DE_24h$log2FoldChange < 0)]

write.csv(up_24h,
          file = file.path(results_dir, "24h_up.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(dn_24h,
          file = file.path(results_dir, "24h_dn.csv"),
          quote = FALSE, row.names = FALSE)
save(DE_24h, file = file.path(results_dir, "DE_24h.RData"))

create_volcano(DE_24h, "24h", results_dir)
summarize_de(DE_24h, "24h", results_dir)


#############################################
## 5. GSEA with msigdbr (BP)
#############################################

m_df <- msigdbr(species = "Mus musculus",
                category = "C5",
                subcategory = "BP") %>%
  dplyr::select(gs_name, gene_symbol)

# 90min
forgesea_90 <- DE_90min$log2FoldChange
names(forgesea_90) <- rownames(DE_90min)
forgesea_90 <- forgesea_90[!is.na(forgesea_90)]

gsea_90 <- GSEA(
  geneList     = sort(forgesea_90, decreasing = TRUE),
  TERM2GENE    = m_df,
  pvalueCutoff = 1,
  maxGSSize    = 10000
)

# 24h  
forgesea_24 <- DE_24h$log2FoldChange
names(forgesea_24) <- rownames(DE_24h)
forgesea_24 <- forgesea_24[!is.na(forgesea_24)]

gsea_24 <- GSEA(
  geneList     = sort(forgesea_24, decreasing = TRUE),
  TERM2GENE    = m_df,
  pvalueCutoff = 1,
  maxGSSize    = 10000
)

# Dotplots
p_gsea_90 <- dotplot(gsea_90, showCategory = 20) + ggtitle("GSEA at 90min")
p_gsea_24 <- dotplot(gsea_24, showCategory = 20) + ggtitle("GSEA at 24h")

ggsave(file.path(results_dir, "GSEA_90min.png"), p_gsea_90, width = 7, height = 10, dpi = 300)
ggsave(file.path(results_dir, "GSEA_24h.png"), p_gsea_24, width = 7, height = 10, dpi = 300)


library(ggplot2)
library(dplyr)

custom_gsea_dotplot <- function(gsea_res, n = 20, title = "") {
  
  df <- as.data.frame(gsea_res)
  
  # Prendi i top n per |NES|
  df <- df[order(abs(df$NES), decreasing = TRUE), ] %>% 
    head(n)
  
  df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
  
  ggplot(df, aes(x = NES, y = Description)) +
    geom_point(aes(size = -log10(p.adjust), color = NES)) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = title,
      x = "NES (Normalized Enrichment Score)",
      y = "",
      color = "NES",
      size = "-log10(adj p)"
    )
}

# Plots
p_gsea_90_signed <- custom_gsea_dotplot(gsea_90, n = 20, 
                                        title = "GSEA 90min (signed NES)")
p_gsea_24_signed <- custom_gsea_dotplot(gsea_24, n = 20, 
                                        title = "GSEA 24h (signed NES)")

ggsave(file.path(results_dir, "GSEA_90min_signedNES.png"),
       p_gsea_90_signed, width = 12, height = 10, dpi = 300)

ggsave(file.path(results_dir, "GSEA_24h_signedNES.png"),
       p_gsea_24_signed, width = 12, height = 10, dpi = 300)



#############################################
## 6. GO over-representation (up/down)
#############################################

ego90minup <- enrichGO(
  gene          = up_90min,
  universe      = rownames(DE_90min),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

ego90mindn <- enrichGO(
  gene          = dn_90min,
  universe      = rownames(DE_90min),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

ego24hup <- enrichGO(
  gene          = up_24h,
  universe      = rownames(DE_24h),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

ego24hdn <- enrichGO(
  gene          = dn_24h,
  universe      = rownames(DE_24h),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

# Dotplots
ggsave(file.path(results_dir, "GO_ego90min_up.png"),
       dotplot(ego90minup, showCategory = 20) + ggtitle("GO – 90min up"),
       width = 7, height = 7, dpi = 300)

ggsave(file.path(results_dir, "GO_ego90min_dn.png"),
       dotplot(ego90mindn, showCategory = 20) + ggtitle("GO – 90min down"),
       width = 7, height = 7, dpi = 300)

ggsave(file.path(results_dir, "GO_ego24h_up.png"),
       dotplot(ego24hup, showCategory = 10) + ggtitle("GO – 24h up"),
       width = 7, height = 7, dpi = 300)

ggsave(file.path(results_dir, "GO_ego24h_dn.png"),
       dotplot(ego24hdn, showCategory = 20) + ggtitle("GO – 24h down"),
       width = 7, height = 7, dpi = 300)


#############################################
## 7. Common DEGs between 90min and 24h
#############################################

DE_common <- intersect(
  rownames(DE_90min)[DE_90min$padj < 0.05],
  rownames(DE_24h)[DE_24h$padj < 0.05]
)

write.csv(DE_common,
          file = file.path(results_dir, "DE_common.csv"),
          quote = FALSE, row.names = FALSE)

egocommon <- enrichGO(
  gene          = DE_common,
  universe      = rownames(DE_24h),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

ggsave(file.path(results_dir, "GO_common.png"),
       dotplot(egocommon, showCategory = 20) +
         ggtitle("GO – DEGs in both 90min & 24h"),
       width = 7, height = 12, dpi = 300)


#############################################
## 8. Expression boxplots (genes vs time/treatment)
#############################################

expr_df <- cbind.data.frame(t(RPMlog), meta)
expr_df$time <- factor(expr_df$time, levels = c("90min", "24h"))

genes_to_plot <- c("Arc", "Mbp", "Tcf4", "Tet3", "Per2")

for (g in genes_to_plot) {
  p <- ggplot(expr_df, aes(x = time, y = .data[[g]], fill = treatment)) +
    geom_boxplot() +
    theme_bw() +
    labs(title = paste(g, "expression"), y = "log2 RPM + 1", x = "Time")
  
  ggsave(file.path(results_dir, paste0("Expr_", g, ".png")),
         p, width = 5, height = 5, dpi = 300)
}


