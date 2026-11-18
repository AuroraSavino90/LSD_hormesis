#############################################
## LSD PFC – QC, Differential expression and enrichment
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
library(ggrepel)
library(openxlsx)

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

create_volcano <- function(res,
                           contrast_name,
                           results_dir,
                           alpha = 0.05,
                           n_label = 10) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj), ]
  
  df$negLog10Padj <- -log10(df$padj)
  
  
  df$status <- "NS"
  df$status[df$padj < alpha & df$log2FoldChange >  0] <- "Up"
  df$status[df$padj < alpha & df$log2FoldChange <  0] <- "Down"
  
  ## Select the most significant genes to label
  df_label <- df[order(df$padj), ]
  df_label <- head(df_label, n_label)
  
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
    
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed") +
    
    geom_text_repel(
      data = df_label,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    )
  
  ggsave(
    file.path(results_dir, paste0("Volcano_", contrast_name, ".png")),
    p, width = 6, height = 5, dpi = 300
  )
  ggsave(
    file.path(results_dir, paste0("Volcano_", contrast_name, ".pdf")),
    p, width = 6, height = 5
  )
}

###gsea plot with color corresponding to NES sign
custom_gsea_dotplot <- function(gsea_res, n = 20, title = "") {
  
  df <- as.data.frame(gsea_res)
  
  # Take top for |NES|
  df <- df[order(abs(df$NES), decreasing = TRUE), ] %>% 
    head(n)
  
  df$Description <- factor(df$Description, levels = df$Description[order(df$NES)])
  
  ggplot(df, aes(x = NES, y = Description)) +
    geom_point(aes(size = -log10(p.adjust), color = NES)) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0
    ) +
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


#############################################
## Load & normalize counts
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
## PCA
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

p_pca_both <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = treatment, shape=time)) +
  geom_point(size = 3) +
  labs(title = "PCA – time")

ggsave(file.path(results_dir, "PCA_treatment.pdf"), p_pca_treat, width = 5, height = 5)
ggsave(file.path(results_dir, "PCA_time.pdf"),      p_pca_time, width = 5, height = 5)
ggsave(file.path(results_dir, "PCA_treatment.png"), p_pca_treat, width = 5, height = 5, dpi = 300)
ggsave(file.path(results_dir, "PCA_time.png"),      p_pca_time,  width = 5, height = 5, dpi = 300)
ggsave(file.path(results_dir, "PCA_both.png"), p_pca_both, width = 5, height = 5, dpi = 300)
ggsave(file.path(results_dir, "PCA_both.pdf"),      p_pca_both,  width = 5, height = 5)


#############################################
## IEG expression 
#############################################

expr_df <- cbind.data.frame(t(RPMlog), meta)
expr_df$time <- factor(expr_df$time, levels = c("90min", "24h"))

genes_to_plot <- c("Arc", "Bdnf","Egr1", "Fos")

for (g in genes_to_plot) {
  p <- ggplot(expr_df, aes(x = time, y = .data[[g]], fill = treatment)) +
    geom_boxplot() +
    labs(title = paste(g, "expression"), y = "log2 RPM + 1", x = "Time")
  
  ggsave(file.path(results_dir, paste0("Expr_", g, ".png")),
         p, width = 4, height = 5, dpi = 300)
  ggsave(file.path(results_dir, paste0("Expr_", g, ".pdf")),      
         p, width = 4, height = 5)
  
}


#############################################
## DEGs
#############################################

## 90min
dds_90min <- DESeqDataSetFromMatrix(
  countData = counts_name[,meta$time=="90min"],
  colData   = data.frame(treatment = meta$treatment[meta$time=="90min"],
                         row.names = colnames(counts_name)[meta$time=="90min"]),
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

## 24h
dds_24h <- DESeqDataSetFromMatrix(
  countData = counts_name[, meta$time=="24h"],
  colData   = data.frame(treatment = meta$treatment[meta$time=="24h"],
                         row.names = colnames(counts_name)[meta$time=="24h"]),
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



df_DE_90min <- as.data.frame(DE_90min)
df_DE_90min$gene <- rownames(df_DE_90min)
df_DE_90min <- df_DE_90min[, c("gene", setdiff(colnames(df_DE_90min), "gene"))]

df_DE_24h <- as.data.frame(DE_24h)
df_DE_24h$gene <- rownames(df_DE_24h)
df_DE_24h <- df_DE_24h[, c("gene", setdiff(colnames(df_DE_24h), "gene"))]

wb <- createWorkbook()

addWorksheet(wb, "DESeq2_90min")
writeData(wb, sheet = "DESeq2_90min", x = df_DE_90min)

addWorksheet(wb, "DESeq2_24h")
writeData(wb, sheet = "DESeq2_24h", x = df_DE_24h)

saveWorkbook(
  wb,
  file = file.path(results_dir, "DESeq2_results_90min_24h.xlsx"),
  overwrite = TRUE
)


create_volcano(DE_24h, "24h", results_dir)


#############################################
## GSEA with msigdbr (BP)
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

# Plots
p_gsea_90_signed <- custom_gsea_dotplot(gsea_90, n = 20, 
                                        title = "GSEA 90min (signed NES)")
p_gsea_24_signed <- custom_gsea_dotplot(gsea_24, n = 20, 
                                        title = "GSEA 24h (signed NES)")

ggsave(file.path(results_dir, "GSEA_90min_signedNES.png"),
       p_gsea_90_signed, width = 9, height = 8, dpi = 300)

ggsave(file.path(results_dir, "GSEA_24h_signedNES.png"),
       p_gsea_24_signed, width = 11, height = 8, dpi = 300)

ggsave(file.path(results_dir, "GSEA_90min_signedNES.pdf"),
       p_gsea_90_signed, width = 9, height = 8)

ggsave(file.path(results_dir, "GSEA_24h_signedNES.pdf"),
       p_gsea_24_signed, width = 11, height = 8)



#############################################
## GO over-representation (up/down separately)
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


ggsave(file.path(results_dir, "GO_ego24h_up.png"),
       dotplot(ego24hup, showCategory = 20) + ggtitle("GO – 24h up"),
       width = 7, height = 7, dpi = 300)


ggsave(file.path(results_dir, "GO_ego90min_up.pdf"),
       dotplot(ego90minup, showCategory = 20) + ggtitle("GO – 90min up"),
       width = 7, height = 12)


ggsave(file.path(results_dir, "GO_ego24h_up.pdf"),
       dotplot(ego24hup, showCategory = 20) + ggtitle("GO – 24h up"),
       width = 7, height = 12)


#############################################
## Common DEGs between 90min and 24h
#############################################

DE_common <- intersect(
  rownames(DE_90min)[DE_90min$log2FoldChange>0 & DE_90min$padj < 0.05],
  rownames(DE_24h)[DE_24h$log2FoldChange>0 & DE_24h$padj < 0.05]
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
       width = 7, height = 7, dpi = 300)

ggsave(file.path(results_dir, "GO_common.pdf"),
       dotplot(egocommon, showCategory = 20) +
         ggtitle("GO – DEGs in both 90min & 24h"),
       width = 7, height = 12)


#############save xlsx

df_ego90minup  <- as.data.frame(ego90minup)
df_ego90mindn  <- as.data.frame(ego90mindn)
df_ego24hup    <- as.data.frame(ego24hup)
df_ego24hdn    <- as.data.frame(ego24hdn)
df_egocommon   <- as.data.frame(egocommon)

wb <- createWorkbook()

addWorksheet(wb, "90min_up")
writeData(wb, sheet = "90min_up",   x = df_ego90minup)

addWorksheet(wb, "90min_down")
writeData(wb, sheet = "90min_down", x = df_ego90mindn)

addWorksheet(wb, "24h_up")
writeData(wb, sheet = "24h_up",     x = df_ego24hup)

addWorksheet(wb, "24h_down")
writeData(wb, sheet = "24h_down",   x = df_ego24hdn)

addWorksheet(wb, "common_90min_24h")
writeData(wb, sheet = "common_90min_24h", x = df_egocommon)

saveWorkbook(
  wb,
  file = file.path(results_dir, "GO_enrichment_results.xlsx"),
  overwrite = TRUE
)


