library(openxlsx)
library(GEOquery)
library(DESeq2)
library(limma)
library(biomaRt)

## General params & paths
date        <- "20251124"
results_dir <- file.path("results", date)

# ID -> symbol 
# - anno must be a 2-column object: anno[,1] = ID, anno[,2] = symbol
# - returns a numeric matrix with unique HGNC symbol rownames
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
    if (is.na(id) || !nzchar(id)) return(NA_character_)
    syms <- sym_by_id[[id]]
    if (is.null(syms) || length(syms) == 0) return(NA_character_)
    if (length(syms) == 1) return(syms)
    non_ps <- syms[!grepl("-ps", syms)]
    if (length(non_ps) == 1) return(non_ps)
    NA_character_
  }
  
  ids     <- as.character(rownames(data))
  map_vec <- vapply(ids, resolve_symbol, FUN.VALUE = character(1), USE.NAMES = FALSE)
  
  keep <- !is.na(map_vec) & nzchar(map_vec)
  if (!any(keep)) stop("After filtering, no rows remain with a valid symbol.")
  
  data    <- data[keep, , drop = FALSE]
  map_vec <- map_vec[keep]
  
  rms    <- rowMeans(data, na.rm = TRUE)
  groups <- split(seq_len(nrow(data)), map_vec)
  keep_idx <- unlist(lapply(groups, function(ix) ix[which.max(rms[ix])]), use.names = FALSE)
  
  data <- data[keep_idx, , drop = FALSE]
  rownames(data) <- map_vec[keep_idx]
  data
}

##### VARIOUS FORMS OF HORMETIC STRESS, IMPACT ON BRAIN
# download data
# change rownames into gene symbols
# define conditions of interest based on metadata
# compute DEGs with DESeq2 or limma

##########################################################################
########## Fasting
##########################################################################

####################################
###### GSE236077
####################################

geo <- getGEO("GSE236077", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE236077      <- read.csv("data/Fasting/GSE236077_counts.csv.gz",
                           row.names = 1)
GSE236077_meta <- read.csv("data/Fasting/GSE236077_gsm_meta_long.csv.gz")

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE236077),
             mart       = ensembl.mouse)

GSE236077 <- changenames(data = as.matrix(GSE236077), anno = ids)

normid <- with(GSE236077_meta,
               which(region == "PFC" & diet == "chow"))
fastid <- with(GSE236077_meta,
               which(region == "PFC" & diet == "fast"))

normid <- intersect(colnames(GSE236077), GSE236077_meta[normid, "GSM"])
fastid <- intersect(colnames(GSE236077), GSE236077_meta[fastid, "GSM"])

## Differentially expressed genes (counts -> DESeq2)
dds_fasting <- DESeqDataSetFromMatrix(
  countData = GSE236077[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting <- DESeq(dds_fasting)
DE_fasting_GSE236077 <- results(dds_fasting,
                                contrast = c("treatment", "Fast", "Norm"))


####################################
###### GSE252873
####################################

geo <- getGEO("GSE252873", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")


GSE252873      <- read.csv("data/Fasting/GSE252873_RNAseq_fasting_sham.csv.gz",
                           row.names = 1)
GSE252873_meta <- pData(geo$GSE252873_series_matrix.txt.gz)
# no way of matching data and metadata 


####################################
###### GSE253612
####################################

geo <- getGEO("GSE253612", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE253612      <- read.csv("data/Fasting/GSE253612_Newman_2020_Mouse_RNAseq_counts_all.star.tsv.gz",
                           row.names = 1, sep = "\t")
GSE253612_meta <- pData(geo$GSE253612_series_matrix.txt.gz)

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE253612),
             mart       = ensembl.mouse)

GSE253612 <- changenames(data = as.matrix(GSE253612), anno = ids)


GSE253612_meta <- GSE253612_meta[
  match(colnames(GSE253612), GSE253612_meta$title),
]

normid <- which(GSE253612_meta$`tissue:ch1` == "brain" &
                  GSE253612_meta$`treatment:ch1` == "control diet - 12 months")
fastid <- which(GSE253612_meta$`tissue:ch1` == "brain" &
                  GSE253612_meta$`treatment:ch1` == "keto diet - 12 months")

normid <- intersect(colnames(GSE253612), GSE253612_meta[normid, "title"])
fastid <- intersect(colnames(GSE253612), GSE253612_meta[fastid, "title"])

## Differentially expressed genes
dds_fasting_GSE253612 <- DESeqDataSetFromMatrix(
  countData = GSE253612[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE253612 <- DESeq(dds_fasting_GSE253612)
DE_fasting_GSE253612_12m <- results(
  dds_fasting_GSE253612,
  contrast = c("treatment", "Fast", "Norm")
)

## 26 months
normid <- which(GSE253612_meta$`tissue:ch1` == "brain" &
                  GSE253612_meta$`treatment:ch1` == "control diet - 26 months")
fastid <- which(GSE253612_meta$`tissue:ch1` == "brain" &
                  GSE253612_meta$`treatment:ch1` == "keto diet - 26 months")

normid <- intersect(colnames(GSE253612), GSE253612_meta[normid, "title"])
fastid <- intersect(colnames(GSE253612), GSE253612_meta[fastid, "title"])

dds_fasting_GSE253612_26 <- DESeqDataSetFromMatrix(
  countData = GSE253612[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE253612_26 <- DESeq(dds_fasting_GSE253612_26)
DE_fasting_GSE253612_26m <- results(
  dds_fasting_GSE253612_26,
  contrast = c("treatment", "Fast", "Norm")
)


####################################
###### GSE139148 (microarray, limma)
####################################

geo <- getGEO("GSE139148", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE139148      <- exprs(geo$GSE139148_series_matrix.txt.gz)
GSE139148_meta <- pData(geo$GSE139148_series_matrix.txt.gz)
GSE139148_anno <- fData(geo$GSE139148_series_matrix.txt.gz)

GSE139148 <- changenames(
  data = GSE139148,
  anno = GSE139148_anno[, c("ID", "Gene symbol")]
)

normid <- which(GSE139148_meta$`genotype/variation:ch1` == "AgRP+/-" &
                  GSE139148_meta$`treatment:ch1` == "Ad libitum fed")
fastid <- which(GSE139148_meta$`genotype/variation:ch1` == "AgRP+/-" &
                  GSE139148_meta$`treatment:ch1` == "fasted for 30h")

normid <- intersect(colnames(GSE139148),
                    GSE139148_meta[normid, "geo_accession"])
fastid <- intersect(colnames(GSE139148),
                    GSE139148_meta[fastid, "geo_accession"])

data      <- GSE139148[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit         <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE139148 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                  adjust = "BH")

library(GEOquery)
library(limma)


####################################
###### GSE75574 – calorie restriction
####################################

# Scarica serie
geo <- getGEO("GSE75574", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/Fasting")

GSE75574_expr  <- exprs(geo$GSE75574_series_matrix.txt.gz)
GSE75574_meta  <- pData(geo$GSE75574_series_matrix.txt.gz)
GSE75574_anno  <- fData(geo$GSE75574_series_matrix.txt.gz)

# Map probe ID -> gene symbol
GSE75574_expr <- changenames(
  data = GSE75574_expr,
  anno = GSE75574_anno[, c("ID", "Gene symbol")]
)

GSE75574_expr <- log2(GSE75574_expr + 1)

# filter for "cerebral cortex" 
cortex_idx <- grepl("^cerebral cortex$", GSE75574_meta$`tissue:ch1`, ignore.case = TRUE)
GSE75574_expr <- GSE75574_expr[, cortex_idx, drop = FALSE]
GSE75574_meta <- GSE75574_meta[cortex_idx, , drop = FALSE]

# define strain and treatment
strain <- factor(GSE75574_meta$`strain:ch1`)

treatment <- ifelse(
  grepl("Calorie\\s*Restricted", GSE75574_meta$source_name_ch1, ignore.case = TRUE),
  "CR",
  ifelse(
    grepl("Control\\s*diet", GSE75574_meta$source_name_ch1, ignore.case = TRUE),
    "Control",
    NA
  )
)
treatment <- factor(treatment, levels = c("Control", "CR"))

# Contrast CR vs Control for each strain

out_list <- vector("list", length(levels(strain)))
names(out_list) <- levels(strain)

for (st in levels(strain)) {
  keep <- strain == st & !is.na(treatment)
 
  expr_s <- GSE75574_expr[, keep, drop = FALSE]
  trt_s  <- droplevels(treatment[keep])
  
  design_s <- model.matrix(~ 0 + trt_s)
  colnames(design_s) <- levels(trt_s)  # "Control","CR"
  
  fit_s  <- lmFit(expr_s, design_s)
  contr  <- makeContrasts(CRvsControl = CR - Control, levels = design_s)
  fit_s2 <- eBayes(contrasts.fit(fit_s, contr))
  
  tt <- topTable(fit_s2, coef = "CRvsControl", number = Inf, adjust.method = "BH")
  tt$feature <- rownames(tt)
  tt$strain  <- st
  
  out_list[[st]] <- tt
}


clean_strain_name <- function(x) gsub("[^A-Za-z0-9]+", "", x)

for (st in names(out_list)) {
  tt <- out_list[[st]]
  if (is.null(tt)) next
  
  obj_name <- paste0("DE_fasting_GSE75574_", clean_strain_name(st))
  assign(obj_name, tt, envir = .GlobalEnv)
  message("Creato oggetto: ", obj_name, " (n = ", nrow(tt), " righe)")
}


####################################
###### GSE57110 (rat, counts, DESeq2)
####################################

geo<-getGEO("GSE57110", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/Fasting")
GSE57110<-exprs(geo$GSE57110_series_matrix.txt.gz)
GSE57110_meta<-pData(geo$GSE57110_series_matrix.txt.gz)

files<-list.files("Data/Fasting/GSE57110_RAW")
GSE57110<-read.csv(paste("Data/Fasting/GSE57110_RAW/", files[1], sep=""), sep="\t", row.names = 1)
for(f in files[-1]){
  GSE57110<-cbind(GSE57110, read.csv(paste("Data/Fasting/GSE57110_RAW/", f, sep=""), sep="\t", row.names = 1))
}
colnames(GSE57110)<-gsm_ids <- sub("^(GSM[0-9]+).*", "\\1", files)
match(colnames(GSE57110), GSE57110_meta[,2])


ensembl.rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

ids<-getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
           filters = 'ensembl_gene_id',
           values = rownames(GSE57110),
           mart = ensembl.rat)
GSE57110<-changenames(data=GSE57110, anno=ids)

normid<-which(GSE57110_meta$`treatment:ch1`=="6 months ad libitum")
fastid<-which(GSE57110_meta$`treatment:ch1`=="6 months calorie restriction")
normid<-intersect(colnames(GSE57110), GSE57110_meta[normid,"geo_accession"])
fastid<-intersect(colnames(GSE57110), GSE57110_meta[fastid,"geo_accession"])

## Differentially expressed genes
dds_fasting_GSE57110_6<- DESeqDataSetFromMatrix(countData = GSE57110[,c(fastid, normid)],
                                                colData = data.frame(treatment=c(rep("Fast", length(fastid)), rep("Norm", length(normid)))),
                                                design= ~ treatment)
dds_fasting_GSE57110_6 <- DESeq(dds_fasting_GSE57110_6)
DE_fasting_GSE57110_6 <- results(dds_fasting_GSE57110_6, contrast = c("treatment","Fast", "Norm"))

normid<-which(GSE57110_meta$`treatment:ch1`=="12 months ad libitum")
fastid<-which(GSE57110_meta$`treatment:ch1`=="12 months calorie restriction")
normid<-intersect(colnames(GSE57110), GSE57110_meta[normid,"geo_accession"])
fastid<-intersect(colnames(GSE57110), GSE57110_meta[fastid,"geo_accession"])

## Differentially expressed genes
dds_fasting_GSE57110_12<- DESeqDataSetFromMatrix(countData = GSE57110[,c(fastid, normid)],
                                                colData = data.frame(treatment=c(rep("Fast", length(fastid)), rep("Norm", length(normid)))),
                                                design= ~ treatment)
dds_fasting_GSE57110_12 <- DESeq(dds_fasting_GSE57110_12)
DE_fasting_GSE57110_12 <- results(dds_fasting_GSE57110_12, contrast = c("treatment","Fast", "Norm"))

normid<-which(GSE57110_meta$`treatment:ch1`=="28 months ad libitum")
fastid<-which(GSE57110_meta$`treatment:ch1`=="28 months calorie restriction")
normid<-intersect(colnames(GSE57110), GSE57110_meta[normid,"geo_accession"])
fastid<-intersect(colnames(GSE57110), GSE57110_meta[fastid,"geo_accession"])

## Differentially expressed genes
dds_fasting_GSE57110_28<- DESeqDataSetFromMatrix(countData = GSE57110[,c(fastid, normid)],
                                                 colData = data.frame(treatment=c(rep("Fast", length(fastid)), rep("Norm", length(normid)))),
                                                 design= ~ treatment)
dds_fasting_GSE57110_28 <- DESeq(dds_fasting_GSE57110_28)
DE_fasting_GSE57110_28 <- results(dds_fasting_GSE57110_28, contrast = c("treatment","Fast", "Norm"))



####################################
###### GSE164461 (FPKM, log2, limma)
####################################

geo <- getGEO("GSE164461", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE164461_meta <- pData(geo$GSE164461_series_matrix.txt.gz)
GSE164461      <- read.xlsx("data/Fasting/GSE164461_FPKM.xlsx", rowNames = TRUE)
GSE164461 <- changenames(
  data = as.matrix(GSE164461[, -1]),
  anno = cbind(rownames(GSE164461), GSE164461[, 1])
)
GSE164461 <- log2(GSE164461 + 1)

GSE164461_meta <- GSE164461_meta[
  match(colnames(GSE164461), GSE164461_meta$title),
]

## Female (FWT vs FWF)
normid <- grep("FWT", GSE164461_meta$`treatment:ch1`)
fastid <- grep("FWF", GSE164461_meta$`treatment:ch1`)

normid <- intersect(colnames(GSE164461), GSE164461_meta[normid, "title"])
fastid <- intersect(colnames(GSE164461), GSE164461_meta[fastid, "title"])

data      <- GSE164461[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE164461_F <- topTable(fit_contrast, number = nrow(fit_contrast),
                                    adjust = "BH")

## Male (MWT vs MWF)
normid <- grep("MWT", GSE164461_meta$`treatment:ch1`)
fastid <- grep("MWF", GSE164461_meta$`treatment:ch1`)

normid <- intersect(colnames(GSE164461), GSE164461_meta[normid, "title"])
fastid <- intersect(colnames(GSE164461), GSE164461_meta[fastid, "title"])

data      <- GSE164461[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE164461_M <- topTable(fit_contrast, number = nrow(fit_contrast),
                                    adjust = "BH")



####################################
###### GSE135945 (counts, DESeq2)
####################################

geo <- getGEO("GSE135945", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE135945_meta <- pData(geo$GSE135945_series_matrix.txt.gz)

GSE135945_IF     <- read.csv("data/Fasting/GSE135945_readcount_genename_IF_RNA.txt.gz",
                             sep = "\t", row.names = 1)
GSE135945_Refeed <- read.csv("data/Fasting/GSE135945_readcount_genename_Refeed_RNA.txt.gz",
                             sep = "\t", row.names = 1)

GSE135945 <- cbind(GSE135945_IF[, 1:15], GSE135945_Refeed[, 1:15])

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE135945),
             mart       = ensembl.mouse)

GSE135945 <- changenames(data = as.matrix(GSE135945), anno = ids)

normid <- which(GSE135945_meta$`treatment:ch1` == "AL")
fastid <- which(GSE135945_meta$`treatment:ch1` == "IF16")

normid <- intersect(colnames(GSE135945), GSE135945_meta[normid, "description"])
fastid <- intersect(colnames(GSE135945), GSE135945_meta[fastid, "description"])

dds_fasting_GSE135945 <- DESeqDataSetFromMatrix(
  countData = GSE135945[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE135945 <- DESeq(dds_fasting_GSE135945)
DE_fasting_GSE135945 <- results(
  dds_fasting_GSE135945,
  contrast = c("treatment", "Fast", "Norm")
)


####################################
###### GSE116878 (counts, DESeq2)
####################################

geo <- getGEO("GSE116878", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE116878_meta <- pData(geo$GSE116878_series_matrix.txt.gz)
GSE116878      <- read.xlsx("data/Fasting/GSE116878_Readcount.xlsx",
                            rowNames = TRUE)

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE116878),
             mart       = ensembl.mouse)

GSE116878 <- changenames(data = as.matrix(GSE116878), anno = ids)

normid <- intersect(
  grep("Sham", GSE116878_meta$title),
  which(GSE116878_meta$`diet:ch1` == "Ad libitum")
)
fastid <- intersect(
  grep("Sham", GSE116878_meta$title),
  which(GSE116878_meta$`diet:ch1` == "every-other-day fasting")
)

normid <- intersect(colnames(GSE116878),
                    gsub("0", "", GSE116878_meta[normid, "description"]))
fastid <- intersect(colnames(GSE116878),
                    gsub("0", "", GSE116878_meta[fastid, "description"]))

dds_fasting_GSE116878 <- DESeqDataSetFromMatrix(
  countData = GSE116878[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE116878 <- DESeq(dds_fasting_GSE116878)
DE_fasting_GSE116878 <- results(
  dds_fasting_GSE116878,
  contrast = c("treatment", "Fast", "Norm")
)


####################################
###### GSE96644
####################################

geo <- getGEO("GSE96644", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE96644_meta <- pData(geo$GSE96644_series_matrix.txt.gz)

files <- list.files("data/Fasting/GSE96644_RAW")

first <- read.csv(
  file.path("data/Fasting/GSE96644_RAW", files[1]),
  sep = "\t", row.names = 1, header = FALSE
)[-c(1:2), 6, drop = FALSE]

GSE96644 <- first
for (f in files[-1]) {
  tmp <- read.csv(
    file.path("data/Fasting/GSE96644_RAW", f),
    sep = "\t", row.names = 1, header = FALSE
  )[-c(1:2), 6, drop = FALSE]
  GSE96644 <- cbind(GSE96644, tmp)
}
colnames(GSE96644) <- gsm_ids <- sub("^(GSM[0-9]+).*", "\\1", files)

raw_first <- read.csv(
  file.path("data/Fasting/GSE96644_RAW", files[1]),
  sep = "\t", row.names = 1, header = FALSE
)
rownames(GSE96644) <- rownames(raw_first)[-c(1:2)]

GSE96644 <- as.matrix(GSE96644)
storage.mode(GSE96644) <- "numeric"

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE96644),
             mart       = ensembl.mouse)

GSE96644 <- changenames(data = GSE96644, anno = ids)

stopifnot(all(colnames(GSE96644) %in% GSE96644_meta$geo_accession))

## Young
normid <- which(
  GSE96644_meta$description ==
    "Ad libitum fed control young male C57BL/6J mouse"
)
fastid <- which(
  GSE96644_meta$description ==
    "Intermittent fasted young male C57BL/6J mouse"
)
normid <- intersect(colnames(GSE96644), GSE96644_meta[normid, "geo_accession"])
fastid <- intersect(colnames(GSE96644), GSE96644_meta[fastid, "geo_accession"])

dds_fasting_GSE96644_Y <- DESeqDataSetFromMatrix(
  countData = GSE96644[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE96644_Y <- DESeq(dds_fasting_GSE96644_Y)
DE_fasting_GSE96644_Y <- results(
  dds_fasting_GSE96644_Y,
  contrast = c("treatment", "Fast", "Norm")
)

## Old
normid <- which(
  GSE96644_meta$description ==
    "Ad libitum fed control old male C57BL/6J mouse"
)
fastid <- which(
  GSE96644_meta$description ==
    "Intermittent fasted old male C57BL/6J mouse"
)
normid <- intersect(colnames(GSE96644), GSE96644_meta[normid, "geo_accession"])
fastid <- intersect(colnames(GSE96644), GSE96644_meta[fastid, "geo_accession"])

dds_fasting_GSE96644_O <- DESeqDataSetFromMatrix(
  countData = GSE96644[, c(fastid, normid)],
  colData   = data.frame(
    treatment = c(rep("Fast", length(fastid)), rep("Norm", length(normid)))
  ),
  design    = ~ treatment
)
dds_fasting_GSE96644_O <- DESeq(dds_fasting_GSE96644_O)
DE_fasting_GSE96644_O <- results(
  dds_fasting_GSE96644_O,
  contrast = c("treatment", "Fast", "Norm")
)


####################################
###### GSE236731 
####################################

geo <- getGEO("GSE236731", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE236731      <- read.csv("data/Fasting/GSE236731_deseq_normalized_counts.txt.gz",
                           sep = "\t", row.names = 1)
GSE236731      <- log2(as.matrix(GSE236731) + 1)
GSE236731_meta <- pData(geo$GSE236731_series_matrix.txt.gz)

## 0h
normid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "ALF")
fastid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "TRF")

normid <- intersect(colnames(GSE236731), GSE236731_meta[normid, "title"])
fastid <- intersect(colnames(GSE236731), GSE236731_meta[fastid, "title"])
normid <- normid[grep("_0_", normid)]
fastid <- fastid[grep("_0_", fastid)]

data      <- GSE236731[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE236731_0 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                    adjust = "BH")

## 12h
normid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "ALF")
fastid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "TRF")

normid <- intersect(colnames(GSE236731), GSE236731_meta[normid, "title"])
fastid <- intersect(colnames(GSE236731), GSE236731_meta[fastid, "title"])
normid <- normid[grep("_12_", normid)]
fastid <- fastid[grep("_12_", fastid)]

data      <- GSE236731[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE236731_12 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                     adjust = "BH")

## 18h
normid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "ALF")
fastid <- which(GSE236731_meta$`genotype:ch1` == "NTG" &
                  GSE236731_meta$`treatment:ch1` == "TRF")

normid <- intersect(colnames(GSE236731), GSE236731_meta[normid, "title"])
fastid <- intersect(colnames(GSE236731), GSE236731_meta[fastid, "title"])
normid <- normid[grep("_18_", normid)]
fastid <- fastid[grep("_18_", fastid)]

data      <- GSE236731[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE236731_18 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                     adjust = "BH")


####################################
###### GSE64226 (microarray, limma)
####################################

geo <- getGEO("GSE64226", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE64226      <- exprs(geo$GSE64226_series_matrix.txt.gz)
GSE64226_meta <- pData(geo$GSE64226_series_matrix.txt.gz)
GSE64226_anno <- fData(geo$GSE64226_series_matrix.txt.gz)

GSE64226 <- changenames(
  data = GSE64226,
  anno = GSE64226_anno[, c("ID", "Gene symbol")]
)
GSE64226 <- log2(GSE64226 + 1)

normid <- which(GSE64226_meta$`diet:ch1` == "control")
fastid <- which(GSE64226_meta$`diet:ch1` == "calorie restricted")

normid <- intersect(colnames(GSE64226), GSE64226_meta[normid, "geo_accession"])
fastid <- intersect(colnames(GSE64226), GSE64226_meta[fastid, "geo_accession"])

data      <- GSE64226[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE64226 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                 adjust = "BH")


####################################
###### GSE11291 (microarray, limma)
####################################

geo <- getGEO("GSE11291", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Fasting")

GSE11291      <- exprs(geo$GSE11291_series_matrix.txt.gz)
GSE11291_meta <- pData(geo$GSE11291_series_matrix.txt.gz)
GSE11291_anno <- fData(geo$GSE11291_series_matrix.txt.gz)

GSE11291 <- changenames(
  data = GSE11291,
  anno = GSE11291_anno[, c("ID", "Gene symbol")]
)
GSE11291 <- log2(GSE11291 + 1)

normid <- which(GSE11291_meta$source_name_ch1 ==
                  "neocortex, 30 months of age, control diet")
fastid <- which(GSE11291_meta$source_name_ch1 ==
                  "neocortex, 30 months of age, CR diet")

normid <- intersect(colnames(GSE11291), GSE11291_meta[normid, "geo_accession"])
fastid <- intersect(colnames(GSE11291), GSE11291_meta[fastid, "geo_accession"])

data      <- GSE11291[, c(normid, fastid)]
treatment <- factor(
  c(rep("Norm", length(normid)), rep("Fast", length(fastid))),
  levels = c("Norm", "Fast")
)

design <- model.matrix(~ 0 + treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast - Norm, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_fasting_GSE11291_O <- topTable(fit_contrast, number = nrow(fit_contrast),
                                   adjust = "BH")

####################################
###### GSE200820
####################################

geo<-getGEO("GSE200820", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/Fasting")
GSE200820<-exprs(geo$GSE200820_series_matrix.txt.gz)
GSE200820_meta<-pData(geo$GSE200820_series_matrix.txt.gz)
GSE200820<-read.csv("data/Fasting/GSE200820_CRvCD_normalised_counts_table.txt.gz", sep="\t")
GSE200820<-changenames(data=GSE200820[,-1], anno=cbind(rownames(GSE200820), GSE200820[,1]))

GSE200820<-log2(GSE200820+1)

normid<-grep("CT", colnames(GSE200820))
fastid<-grep("CR", colnames(GSE200820))
normid<-colnames(GSE200820)[normid]
fastid<-colnames(GSE200820)[fastid]

data<-GSE200820[,c(normid, fastid)]
treatment<-c(rep("Norm", length(normid)), rep("Fast", length(fastid)))
treatment<-factor(treatment, levels=c("Norm", "Fast"))

design <- model.matrix(~ 0+treatment)
colnames(design) <- c("Norm", "Fast")
cont_matrix <- makeContrasts(FastvsControl = Fast-Norm, levels=design)

# Fit the expression matrix to a linear model
fit <- lmFit(data, design)
# Compute contrast
fit_contrast <- contrasts.fit(fit, cont_matrix)
# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast <- eBayes(fit_contrast)
DEl_fasting_GSE200820 <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")

###################################################################################################
# HORMESIS: preconditioning, hypoxia, exercise
###################################################################################################

####################################
###### GSE122107 – Ischemic preconditioning (counts, DESeq2)
####################################

geo <- getGEO("GSE122107", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Hormesis")

GSE122107_meta <- pData(geo$GSE122107_series_matrix.txt.gz)
files <- list.files("data/Hormesis/GSE122107_RAW/")

GSE122107 <- read.csv(
  file.path("data/Hormesis/GSE122107_RAW", files[1]),
  sep = "\t", header = FALSE
)
for (f in files[-1]) {
  tmp <- read.csv(
    file.path("data/Hormesis/GSE122107_RAW", f),
    sep = "\t", header = FALSE
  )[, 2, drop = FALSE]
  GSE122107 <- cbind(GSE122107, tmp)
}

rownames(GSE122107) <- GSE122107[, 1]
GSE122107 <- GSE122107[, -1]
GSE122107 <- as.matrix(GSE122107)
storage.mode(GSE122107) <- "numeric"

ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
x_clean <- sub("\\..*", "", rownames(GSE122107))
rownames(GSE122107) <- x_clean

ids <- getBM(attributes = c('ensembl_gene_id', 'mgi_symbol'),
             filters    = 'ensembl_gene_id',
             values     = rownames(GSE122107),
             mart       = ensembl.mouse)

GSE122107 <- changenames(data = GSE122107, anno = ids)

dds_prec <- DESeqDataSetFromMatrix(
  countData = GSE122107,
  colData   = data.frame(
    treatment = rep(c("Preconditioning", "Ctrl"), each = 3)
  ),
  design    = ~ treatment
)
dds_prec <- DESeq(dds_prec)
DE_prec_GSE122107 <- results(
  dds_prec,
  contrast = c("treatment", "Preconditioning", "Ctrl")
)


####################################
###### GSE19709 – Hypoxia (microarray, limma)
####################################

geo <- getGEO("GSE19709", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Hormesis")

GSE19709      <- exprs(geo$GSE19709_series_matrix.txt.gz)
GSE19709_meta <- pData(geo$GSE19709_series_matrix.txt.gz)
GSE19709_anno <- fData(geo$GSE19709_series_matrix.txt.gz)

GSE19709 <- changenames(
  data = GSE19709,
  anno = GSE19709_anno[, c("ID", "Gene symbol")]
)
GSE19709 <- log2(GSE19709 + 1)

## Cortex: hypoxia vs control (Stim vs Ctrl)
normid <- which(GSE19709_meta$`brain region:ch1` == "cerebral cortex" &
                  GSE19709_meta$`treatment protocol:ch1` == "control")
stimid <- which(GSE19709_meta$`brain region:ch1` == "cerebral cortex" &
                  GSE19709_meta$`treatment protocol:ch1` ==
                  "hypoxia, followed by variable-length normoxia")

normid <- intersect(colnames(GSE19709), GSE19709_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE19709), GSE19709_meta[stimid, "geo_accession"])

data  <- GSE19709[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
# DEl_hyp_GSE19709 <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")

## Cortex: hypoxia only vs control (Stim vs Ctrl)
normid <- which(GSE19709_meta$`brain region:ch1` == "cerebral cortex" &
                  GSE19709_meta$`treatment protocol:ch1` == "control")
stimid <- which(GSE19709_meta$`brain region:ch1` == "cerebral cortex" &
                  GSE19709_meta$`treatment protocol:ch1` == "hypoxia")

normid <- intersect(colnames(GSE19709), GSE19709_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE19709), GSE19709_meta[stimid, "geo_accession"])

data  <- GSE19709[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
# DEl_hyp_GSE19709_2 <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")


####################################
###### GSE32529 – brief ischemia (microarray, limma)
####################################

geo <- getGEO("GSE32529", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Hormesis")

GSE32529      <- exprs(geo$GSE32529_series_matrix.txt.gz)
GSE32529_meta <- pData(geo$GSE32529_series_matrix.txt.gz)
GSE32529_anno <- fData(geo$GSE32529_series_matrix.txt.gz)

GSE32529 <- changenames(
  data = GSE32529,
  anno = GSE32529_anno[, c("ID", "Gene symbol")]
)

## 3h
normid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "sham 3h")
stimid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "brief ischemia 3h")

normid <- intersect(colnames(GSE32529), GSE32529_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE32529), GSE32529_meta[stimid, "geo_accession"])

data  <- GSE32529[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_hyp_GSE32529_3h <- topTable(fit_contrast, number = nrow(fit_contrast),
                                adjust = "BH")

## 24h
normid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "sham 24h")
stimid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "brief ischemia 24h")

normid <- intersect(colnames(GSE32529), GSE32529_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE32529), GSE32529_meta[stimid, "geo_accession"])

data  <- GSE32529[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_hyp_GSE32529_24h <- topTable(fit_contrast, number = nrow(fit_contrast),
                                 adjust = "BH")

## 72h
normid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "sham 72h")
stimid <- which(GSE32529_meta$`tissue:ch1` == "brain cortex" &
                  GSE32529_meta$`treatment:ch1` == "brief ischemia 72h")

normid <- intersect(colnames(GSE32529), GSE32529_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE32529), GSE32529_meta[stimid, "geo_accession"])

data  <- GSE32529[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_hyp_GSE32529_72h <- topTable(fit_contrast, number = nrow(fit_contrast),
                                 adjust = "BH")



####################################
###### GSE39697 – Physical hormesis (exercise, microarray)
####################################

geo <- getGEO("GSE39697", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Exercise")

GSE39697      <- exprs(geo$GSE39697_series_matrix.txt.gz)
GSE39697_meta <- pData(geo$GSE39697_series_matrix.txt.gz)
GSE39697_anno <- fData(geo$GSE39697_series_matrix.txt.gz)

GSE39697 <- changenames(
  data = GSE39697,
  anno = GSE39697_anno[, c("ID", "Gene symbol")]
)
GSE39697 <- log2(GSE39697 + 1)

## 4d – SGZ
normid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-subgranular zone" &
                  GSE39697_meta$`days_running:ch1` == "0")
stimid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-subgranular zone" &
                  GSE39697_meta$`days_running:ch1` == "4")

normid <- intersect(colnames(GSE39697), GSE39697_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE39697), GSE39697_meta[stimid, "geo_accession"])

data  <- GSE39697[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_exer_GSE39697_4d <- topTable(fit_contrast, number = nrow(fit_contrast),
                                 adjust = "BH")

## 30d – SGZ
normid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-subgranular zone" &
                  GSE39697_meta$`days_running:ch1` == "0")
stimid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-subgranular zone" &
                  GSE39697_meta$`days_running:ch1` == "30")

normid <- intersect(colnames(GSE39697), GSE39697_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE39697), GSE39697_meta[stimid, "geo_accession"])

data  <- GSE39697[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_exer_GSE39697_30d <- topTable(fit_contrast, number = nrow(fit_contrast),
                                  adjust = "BH")

## 4d – GCL
normid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-granule cell layer" &
                  GSE39697_meta$`days_running:ch1` == "0")
stimid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-granule cell layer" &
                  GSE39697_meta$`days_running:ch1` == "4")

normid <- intersect(colnames(GSE39697), GSE39697_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE39697), GSE39697_meta[stimid, "geo_accession"])

data  <- GSE39697[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_exer_GSE39697_4d_gc <- topTable(fit_contrast, number = nrow(fit_contrast),
                                    adjust = "BH")

## 30d – GCL
normid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-granule cell layer" &
                  GSE39697_meta$`days_running:ch1` == "0")
stimid <- which(GSE39697_meta$`region:ch1` == "dentate gyrus-granule cell layer" &
                  GSE39697_meta$`days_running:ch1` == "30")

normid <- intersect(colnames(GSE39697), GSE39697_meta[normid, "geo_accession"])
stimid <- intersect(colnames(GSE39697), GSE39697_meta[stimid, "geo_accession"])

data  <- GSE39697[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_exer_GSE39697_30d_gc <- topTable(fit_contrast, number = nrow(fit_contrast),
                                     adjust = "BH")


##############################
## GSE205907
################################

geo <- getGEO("GSE205907", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Exercise/")
GSE205907<-exprs(geo$GSE205907_series_matrix.txt.gz)
GSE205907 <- read.xlsx("data/Exercise/GSE205907_fpkm.xlsx")
GSE205907<-changenames(data=GSE205907[,-1], anno=cbind(rownames(GSE205907), GSE205907[,1]))

GSE205907_meta<-pData(geo$GSE205907_series_matrix.txt.gz)

ctrl<-which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP")
trt<-which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP+exercise")

normid <- which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP")
stimid <- which(GSE205907_meta$`tissue:ch1`=="cortex" & GSE205907_meta$`treatment:ch1`=="MPTP+exercise")

normid <- colnames(GSE205907)[normid]
stimid <- colnames(GSE205907)[stimid]

data  <- GSE205907[, c(normid, stimid)]
group <- factor(
  c(rep("Ctrl", length(normid)), rep("Stim", length(stimid))),
  levels = c("Ctrl", "Stim")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("Ctrl", "Stim")
cont_matrix <- makeContrasts(StimVsCtrl = Stim - Ctrl, levels = design)

fit          <- lmFit(data, design)
fit_contrast <- eBayes(contrasts.fit(fit, cont_matrix))
DEl_exer_GSE205907 <- topTable(fit_contrast, number = nrow(fit_contrast),
                                     adjust = "BH")



ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

####################################
###### GSE271865 – chronic stress (RNA-seq counts + TPM/FPKM)
####################################

geo <- getGEO("GSE271865", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE271865_meta <- pData(geo$GSE271865_series_matrix.txt.gz)

GSE271865 <- read.csv("data/Stress/GSE271865_counts_TPM_FPKM.txt.gz",
                      row.names = 1, sep = "\t")
GSE271865 <- GSE271865[, grep("count", colnames(GSE271865))]
GSE271865 <- as.matrix(GSE271865)
storage.mode(GSE271865) <- "numeric"

ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE271865),
             mart       = ensembl.mouse)
GSE271865 <- changenames(data = GSE271865, anno = ids)
GSE271865 <- as.matrix(GSE271865)
storage.mode(GSE271865) <- "numeric"

ctrl_ids   <- colnames(GSE271865)[1:3]
stress_ids <- colnames(GSE271865)[4:6]

dds_stress_GSE271865 <- DESeqDataSetFromMatrix(
  countData = GSE271865[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE271865 <- DESeq(dds_stress_GSE271865)
DE_stress_GSE271865 <- results(
  dds_stress_GSE271865,
  contrast = c("condition", "Stress", "Ctrl")
)


####################################
###### GSE292290 – female, RNA-seq counts
####################################

geo <- getGEO("GSE292290", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE292290_meta <- pData(geo$GSE292290_series_matrix.txt.gz)

GSE292290 <- read.csv("data/Stress/GSE292290_female_count.txt.gz",
                      row.names = 1, sep = "\t")
GSE292290 <- as.matrix(round(GSE292290))
storage.mode(GSE292290) <- "numeric"

ctrl_ids   <- colnames(GSE292290)[1:3]
stress_ids <- colnames(GSE292290)[4:6]

dds_stress_GSE292290 <- DESeqDataSetFromMatrix(
  countData = GSE292290[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE292290 <- DESeq(dds_stress_GSE292290)
DE_stress_GSE292290 <- results(
  dds_stress_GSE292290,
  contrast = c("condition", "Stress", "Ctrl")
)



####################################
###### GSE226576 – RNA-seq counts (Processed)
####################################

geo <- getGEO("GSE226576", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE226576_meta <- pData(geo$GSE226576_series_matrix.txt.gz)

GSE226576 <- read.csv("data/Stress/GSE226576_Processed.csv",
                      row.names = 3)
GSE226576 <- GSE226576[, grep("Count", colnames(GSE226576))]
GSE226576 <- as.matrix(round(GSE226576))
storage.mode(GSE226576) <- "numeric"


ctrl_ids   <- colnames(GSE226576)[1:3]
stress_ids <- colnames(GSE226576)[7:9]

dds_stress_GSE226576 <- DESeqDataSetFromMatrix(
  countData = GSE226576[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE226576 <- DESeq(dds_stress_GSE226576)
DE_stress_GSE226576 <- results(
  dds_stress_GSE226576,
  contrast = c("condition", "Stress", "Ctrl")
)


####################################
###### GSE233225 – CSS vs non-stressed (log2, limma)
####################################

geo <- getGEO("GSE233225", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE233225_meta <- pData(geo$GSE233225_series_matrix.txt.gz)

GSE233225 <- read.csv("data/Stress/GSE233225_GeneExpression_counts_CSS.csv.gz",
                      row.names = 1)
GSE233225 <- changenames(
  data = as.matrix(GSE233225[, -c(1:5)]),
  anno = cbind(rownames(GSE233225), GSE233225[, 1])
)
GSE233225 <- log2(GSE233225 + 1)

ctrl_ids   <- colnames(GSE233225)[9:16]
stress_ids <- colnames(GSE233225)[1:8]

expr <- GSE233225[, c(ctrl_ids, stress_ids)]
condition <- factor(
  c(rep("Ctrl", length(ctrl_ids)), rep("Stress", length(stress_ids))),
  levels = c("Ctrl", "Stress")
)

design <- model.matrix(~ 0 + condition)
colnames(design) <- c("Ctrl", "Stress")
cont_matrix <- makeContrasts(StressVsCtrl = Stress - Ctrl, levels = design)

fit          <- lmFit(expr, design)
fit_contrast <- contrasts.fit(fit, cont_matrix)
fit_contrast <- eBayes(fit_contrast)
DEl_stress_GSE233225 <- topTable(
  fit_contrast, number = nrow(fit_contrast), adjust = "BH"
)


####################################
###### GSE211276 – chronic stress
####################################

geo <- getGEO("GSE211276", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE211276_meta <- pData(geo$GSE211276_series_matrix.txt.gz)

GSE211276 <- read.xlsx("data/Stress/GSE211276_all.fpkm_anno.xlsx")
GSE211276 <- GSE211276[, -grep("FPKM", colnames(GSE211276)), drop = FALSE]
rownames(GSE211276) <- GSE211276[, 1]
GSE211276 <- GSE211276[, 3:18]
GSE211276 <- as.matrix(GSE211276)
storage.mode(GSE211276) <- "numeric"

ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE211276),
             mart       = ensembl.mouse)
GSE211276 <- changenames(data = GSE211276, anno = ids)
GSE211276 <- as.matrix(GSE211276)
storage.mode(GSE211276) <- "numeric"

## Female: non-stressed vs chronic-stressed
ctrl_idx   <- which(GSE211276_meta$`conditions:ch1` == "non-stressed" &
                      GSE211276_meta$`Sex:ch1` == "female")
stress_idx <- which(GSE211276_meta$`conditions:ch1` == "chronic-stressed" &
                      GSE211276_meta$`Sex:ch1` == "female")

ctrl_ids   <- intersect(colnames(GSE211276), GSE211276_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE211276), GSE211276_meta[stress_idx, "title"])

dds_stress_GSE211276_F <- DESeqDataSetFromMatrix(
  countData = GSE211276[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE211276_F <- DESeq(dds_stress_GSE211276_F)
DE_stress_GSE211276_F <- results(
  dds_stress_GSE211276_F,
  contrast = c("condition", "Stress", "Ctrl")
)

## Male
ctrl_idx   <- which(GSE211276_meta$`conditions:ch1` == "non-stressed" &
                      GSE211276_meta$`Sex:ch1` == "male")
stress_idx <- which(GSE211276_meta$`conditions:ch1` == "chronic-stressed" &
                      GSE211276_meta$`Sex:ch1` == "male")

ctrl_ids   <- intersect(colnames(GSE211276), GSE211276_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE211276), GSE211276_meta[stress_idx, "title"])

dds_stress_GSE211276_M <- DESeqDataSetFromMatrix(
  countData = GSE211276[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE211276_M <- DESeq(dds_stress_GSE211276_M)
DE_stress_GSE211276_M <- results(
  dds_stress_GSE211276_M,
  contrast = c("condition", "Stress", "Ctrl")
)



####################################
###### GSE148629 – chronic restraint stress (counts, DESeq2)
####################################

geo <- getGEO("GSE148629", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE148629_meta <- pData(geo$GSE148629_series_matrix.txt.gz)

GSE148629 <- read.csv("data/Stress/GSE148629_count_table.txt.gz",
                      sep = "\t", row.names = 1)
GSE148629 <- as.matrix(GSE148629)
storage.mode(GSE148629) <- "numeric"

ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE148629),
             mart       = ensembl.mouse)
GSE148629 <- changenames(data = GSE148629, anno = ids)
GSE148629 <- as.matrix(GSE148629)
storage.mode(GSE148629) <- "numeric"

## PFC
ctrl_idx   <- which(GSE148629_meta$`tissue:ch1` == "prefrontal cortex" &
                      GSE148629_meta$`treatment:ch1` == "control (no stress)")
stress_idx <- which(GSE148629_meta$`tissue:ch1` == "prefrontal cortex" &
                      GSE148629_meta$`treatment:ch1` == "multimodal chronic restraint stress")

ctrl_ids   <- intersect(colnames(GSE148629), GSE148629_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE148629), GSE148629_meta[stress_idx, "title"])

dds_stress_GSE148629_PFC <- DESeqDataSetFromMatrix(
  countData = GSE148629[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE148629_PFC <- DESeq(dds_stress_GSE148629_PFC)
DE_stress_GSE148629 <- results(
  dds_stress_GSE148629_PFC,
  contrast = c("condition", "Stress", "Ctrl")
)

## NAc
ctrl_idx   <- which(GSE148629_meta$`tissue:ch1` == "nucleus accumbens" &
                      GSE148629_meta$`treatment:ch1` == "control (no stress)")
stress_idx <- which(GSE148629_meta$`tissue:ch1` == "nucleus accumbens" &
                      GSE148629_meta$`treatment:ch1` == "multimodal chronic restraint stress")

ctrl_ids   <- intersect(colnames(GSE148629), GSE148629_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE148629), GSE148629_meta[stress_idx, "title"])

dds_stress_GSE148629_NAc <- DESeqDataSetFromMatrix(
  countData = GSE148629[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE148629_NAc <- DESeq(dds_stress_GSE148629_NAc)
DE_stress_GSE148629_NAc <- results(
  dds_stress_GSE148629_NAc,
  contrast = c("condition", "Stress", "Ctrl")
)


#####################################################
####### GSE151807 – chronic mild stress (microarray, limma)
#######################################################

geo <- getGEO("GSE151807", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE151807      <- exprs(geo$GSE151807_series_matrix.txt.gz)
GSE151807_meta <- pData(geo$GSE151807_series_matrix.txt.gz)
GSE151807_anno <- fData(geo$GSE151807_series_matrix.txt.gz)

GSE151807 <- changenames(
  data = GSE151807,
  anno = GSE151807_anno[, c("ID", "Gene symbol")]
)

run_limma_stress <- function(expr_mat, meta, ctrl_pattern, stress_pattern) {
  ctrl_idx   <- which(meta$source_name_ch1 == ctrl_pattern)
  stress_idx <- which(meta$source_name_ch1 == stress_pattern)
  ctrl_ids   <- intersect(colnames(expr_mat), meta[ctrl_idx, "geo_accession"])
  stress_ids <- intersect(colnames(expr_mat), meta[stress_idx, "geo_accession"])
  
  expr      <- expr_mat[, c(ctrl_ids, stress_ids), drop = FALSE]
  condition <- factor(
    c(rep("Ctrl", length(ctrl_ids)), rep("Stress", length(stress_ids))),
    levels = c("Ctrl", "Stress")
  )
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("Ctrl", "Stress")
  cont_matrix <- makeContrasts(StressVsCtrl = Stress - Ctrl, levels = design)
  
  fit          <- lmFit(expr, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
}

DEl_stress_GSE151807_PFC <- run_limma_stress(
  GSE151807, GSE151807_meta,
  "Prefrontal cortex from wild-type mice without chronic mild stress treatments.",
  "Prefrontal cortex from wild-type mice with chronic mild stress treatments."
)

DEl_stress_GSE151807_CC <- run_limma_stress(
  GSE151807, GSE151807_meta,
  "Cerebral cortex from wild-type mice without chronic mild stress treatments.",
  "Cerebral cortex from wild-type mice with chronic mild stress treatments."
)

DEl_stress_GSE151807_hipp <- run_limma_stress(
  GSE151807, GSE151807_meta,
  "Hippocampus from wild-type mice without chronic mild stress treatments.",
  "Hippocampus from wild-type mice with chronic mild stress treatments."
)

DEl_stress_GSE151807_amy <- run_limma_stress(
  GSE151807, GSE151807_meta,
  "Amygdala from wild-type mice without chronic mild stress treatments.",
  "Amygdala from wild-type mice with chronic mild stress treatments."
)



#####################################################
####### GSE109315 – CSDS (logCPM, limma)
#######################################################

geo <- getGEO("GSE109315", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")

GSE109315_PFC   <- read.csv("data/Stress/GSE109315_mPFC_counts_logCPM.txt.gz",
                            sep = "\t", row.names = 1)
GSE109315_meta1 <- pData(geo$`GSE109315-GPL13112_series_matrix.txt.gz`)
GSE109315_meta2 <- pData(geo$`GSE109315-GPL19057_series_matrix.txt.gz`)

ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE109315_PFC),
             mart       = ensembl.mouse)
GSE109315_PFC <- changenames(data = as.matrix(GSE109315_PFC), anno = ids)

GSE109315_meta2$title <- gsub("_defeated_mPFC|_control_mPFC", "", GSE109315_meta2$title)

run_limma_csds <- function(expr_mat, meta, strain_value) {
  ctrl_idx   <- which(meta$`strain:ch1` == strain_value &
                        meta$`csds phenotype:ch1` == "control")
  stress_idx <- which(meta$`strain:ch1` == strain_value &
                        meta$`csds phenotype:ch1` == "susceptible")
  ctrl_ids   <- intersect(colnames(expr_mat), meta[ctrl_idx, "title"])
  stress_ids <- intersect(colnames(expr_mat), meta[stress_idx, "title"])
  
  expr      <- expr_mat[, c(ctrl_ids, stress_ids), drop = FALSE]
  condition <- factor(
    c(rep("Ctrl", length(ctrl_ids)), rep("Stress", length(stress_ids))),
    levels = c("Ctrl", "Stress")
  )
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("Ctrl", "Stress")
  cont_matrix <- makeContrasts(StressVsCtrl = Stress - Ctrl, levels = design)
  
  fit          <- lmFit(expr, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
}

DEl_stress_GSE109315_PFC_C57BL <- run_limma_csds(
  GSE109315_PFC, GSE109315_meta2, "C57BL/6Crl"
)
DEl_stress_GSE109315_PFC_DBA <- run_limma_csds(
  GSE109315_PFC, GSE109315_meta2, "DBA/2Crl"
)

## vHPC
GSE109315_vHPC <- read.csv("data/Stress/GSE109315_vHPC_counts_logCPM.txt.gz",
                           sep = "\t", row.names = 1)
ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE109315_vHPC),
             mart       = ensembl.mouse)
GSE109315_vHPC <- changenames(data = as.matrix(GSE109315_vHPC), anno = ids)

GSE109315_meta1$title <- gsub("_defeated_vHPC|_control_vHPC", "", GSE109315_meta1$title)

DEl_stress_GSE109315_vHPC_C57BL <- run_limma_csds(
  GSE109315_vHPC, GSE109315_meta1, "C57BL/6Crl"
)
DEl_stress_GSE109315_vHPC_DBA <- run_limma_csds(
  GSE109315_vHPC, GSE109315_meta1, "DBA/2Crl"
)

## BNST
GSE109315_BNST <- read.csv("data/Stress/GSE109315_BNST_counts_logCPM.txt.gz",
                           sep = "\t", row.names = 1)
ids <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
             filters    = "ensembl_gene_id",
             values     = rownames(GSE109315_BNST),
             mart       = ensembl.mouse)
GSE109315_BNST <- changenames(data = as.matrix(GSE109315_BNST), anno = ids)

GSE109315_meta2$title <- gsub("_defeated_BNST|_control_BNST", "", GSE109315_meta2$title)

DEl_stress_GSE109315_BNST_C57BL <- run_limma_csds(
  GSE109315_BNST, GSE109315_meta2, "C57BL/6Crl"
)
DEl_stress_GSE109315_BNST_DBA <- run_limma_csds(
  GSE109315_BNST, GSE109315_meta2, "DBA/2Crl"
)


##############################
## GSE84183 – UCMS, microarray, limma
################################

geo <- getGEO("GSE84183", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE84183      <- exprs(geo$GSE84183_series_matrix.txt.gz)
GSE84183_meta <- pData(geo$GSE84183_series_matrix.txt.gz)
GSE84183_anno <- fData(geo$GSE84183_series_matrix.txt.gz)

GSE84183 <- changenames(
  data = GSE84183,
  anno = GSE84183_anno[, c("ID", "GENE_SYMBOL")]
)

run_limma_region <- function(expr_mat, meta, tissue_name) {
  ctrl_idx   <- which(meta$`tissue:ch1` == tissue_name &
                        meta$`ucms:ch1`   == "Unstressed")
  stress_idx <- which(meta$`tissue:ch1` == tissue_name &
                        meta$`ucms:ch1`   == "Stressed")
  ctrl_ids   <- intersect(colnames(expr_mat), meta[ctrl_idx, "geo_accession"])
  stress_ids <- intersect(colnames(expr_mat), meta[stress_idx, "geo_accession"])
  
  expr      <- expr_mat[, c(ctrl_ids, stress_ids), drop = FALSE]
  condition <- factor(
    c(rep("Ctrl", length(ctrl_ids)), rep("Stress", length(stress_ids))),
    levels = c("Ctrl", "Stress")
  )
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("Ctrl", "Stress")
  cont_matrix <- makeContrasts(StressVsCtrl = Stress - Ctrl, levels = design)
  
  fit          <- lmFit(expr, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
}

DEl_stress_GSE84183_DG <- run_limma_region(GSE84183, GSE84183_meta, "Dentate Gyrus")
DEl_stress_GSE84183_CG <- run_limma_region(GSE84183, GSE84183_meta, "Cingulate Gyrus")



##############################
## GSE8262 – intermittent/chronic hypoxia (microarray, limma)
################################

geo <- getGEO("GSE8262", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE8262      <- exprs(geo$GSE8262_series_matrix.txt.gz)
GSE8262_meta <- pData(geo$GSE8262_series_matrix.txt.gz)
GSE8262_anno <- fData(geo$GSE8262_series_matrix.txt.gz)
GSE8262 <- changenames(
  data = GSE8262,
  anno = GSE8262_anno[, c("ID", "Gene symbol")]
)

run_limma_indexed <- function(expr_mat, meta, ctrl_idx, stress_idx, label) {
  ctrl_ids   <- intersect(colnames(expr_mat), meta[ctrl_idx, "geo_accession"])
  stress_ids <- intersect(colnames(expr_mat), meta[stress_idx, "geo_accession"])
  
  expr      <- expr_mat[, c(ctrl_ids, stress_ids), drop = FALSE]
  condition <- factor(
    c(rep("Ctrl", length(ctrl_ids)), rep("Stress", length(stress_ids))),
    levels = c("Ctrl", "Stress")
  )
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- c("Ctrl", "Stress")
  cont_matrix <- makeContrasts(StressVsCtrl = Stress - Ctrl, levels = design)
  
  fit          <- lmFit(expr, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  res <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")
  res
}

## Nota: usi indici 1:3,4:6,... come nel tuo script originale
DEl_stress_GSE8262_cortex_CCH <- run_limma_indexed(
  GSE8262, GSE8262_meta, 1:3, 4:6, "cortex_CCH"
)
DEl_stress_GSE8262_cortex_CIH <- run_limma_indexed(
  GSE8262, GSE8262_meta, 1:3, 7:9, "cortex_CIH"
)
DEl_stress_GSE8262_hippx_CCH <- run_limma_indexed(
  GSE8262, GSE8262_meta, 10:12, 13:15, "hippocampus_CCH"
)
DEl_stress_GSE8262_hippx_CIH <- run_limma_indexed(
  GSE8262, GSE8262_meta, 10:12, 16:18, "hippocampus_CIH"
)


##############################
## GSE198725 – VTA/BLA counts (DESeq2)
################################

geo <- getGEO("GSE198725", GSEMatrix = TRUE, AnnotGPL = TRUE,
              destdir = "data/Stress")
GSE198725<-read.csv("data/Stress/GSE198725_HF5YFDRXX_Counts_VTA_BLA.txt.gz", sep="\t")
GSE198725_meta<-pData(geo$GSE198725_series_matrix.txt.gz)
GSE198725<-changenames(data=GSE198725[,-1], anno=cbind(rownames(GSE198725), GSE198725[,1]))


## BLA
ctrl_idx   <- 4:6
stress_idx <- 1:3
ctrl_ids   <- intersect(colnames(GSE198725), GSE198725_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE198725), GSE198725_meta[stress_idx, "title"])

dds_stress_GSE198725_BLA <- DESeqDataSetFromMatrix(
  countData = GSE198725[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE198725_BLA <- DESeq(dds_stress_GSE198725_BLA)
DE_stress_GSE198725_BLA <- results(
  dds_stress_GSE198725_BLA,
  contrast = c("condition", "Stress", "Ctrl")
)

## VTA
ctrl_idx   <- 10:12
stress_idx <- 7:9
ctrl_ids   <- intersect(colnames(GSE198725), GSE198725_meta[ctrl_idx, "title"])
stress_ids <- intersect(colnames(GSE198725), GSE198725_meta[stress_idx, "title"])

dds_stress_GSE198725_VTA <- DESeqDataSetFromMatrix(
  countData = GSE198725[, c(stress_ids, ctrl_ids)],
  colData   = data.frame(
    condition = factor(
      c(rep("Stress", length(stress_ids)), rep("Ctrl", length(ctrl_ids))),
      levels = c("Ctrl", "Stress")
    )
  ),
  design = ~ condition
)
dds_stress_GSE198725_VTA <- DESeq(dds_stress_GSE198725_VTA)
DE_stress_GSE198725_VTA <- results(
  dds_stress_GSE198725_VTA,
  contrast = c("condition", "Stress", "Ctrl")
)




save.image(file=paste(results_dir, "/External_datasets.RData", sep=""))
