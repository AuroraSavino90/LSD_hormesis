#Helpers

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

## ---------------------------------------------------------
## Shared parameters for pathway-level GSEA and sPLS-DA
## (used by 09_NES_sPLSDA.R and 10_LODO.R)
## ---------------------------------------------------------

SEED         <- 576906830
SPLSDA_NCOMP <- 2
SPLSDA_KEEPX <- c(20, 20)

## Fit sPLS-DA on training NES (samples x pathways; NAs allowed) and
## return a projection function that maps any NES matrix onto the
## learned space: NA/Inf -> 0, scaling with TRAINING mean/sd,
## scores = X_scaled %*% loadings. Training centroids are computed
## with the same projection, so training and projected samples share
## exactly the same coordinates system.
fit_project_splsda <- function(X_train, Y_train,
                               ncomp = SPLSDA_NCOMP,
                               keepX = SPLSDA_KEEPX) {
  X_train <- as.matrix(X_train)
  X_train[!is.finite(X_train)] <- NA_real_

  ## drop pathways with no NES in any training dataset: they carry no
  ## information and their NaN training mean would turn every projected
  ## score into NaN
  X_train <- X_train[, colSums(!is.na(X_train)) > 0, drop = FALSE]

  model <- mixOmics::splsda(X_train, Y_train, ncomp = ncomp, keepX = keepX)

  loadX <- model$loadings$X[, seq_len(ncomp), drop = FALSE]
  feats <- rownames(loadX)

  center_tr <- colMeans(X_train[, feats, drop = FALSE], na.rm = TRUE)
  sd_tr     <- apply(X_train[, feats, drop = FALSE], 2, sd, na.rm = TRUE)
  sd_tr[!is.finite(sd_tr) | sd_tr == 0] <- 1

  project <- function(X) {
    X <- as.matrix(X)
    missing_feats <- setdiff(feats, colnames(X))
    if (length(missing_feats) > 0)
      stop("Projection matrix is missing ", length(missing_feats),
           " features used by the model, e.g. ",
           paste(head(missing_feats, 5), collapse = ", "))
    X <- X[, feats, drop = FALSE]
    X[!is.finite(X)] <- 0
    S <- sweep(sweep(X, 2, center_tr, "-"), 2, sd_tr, "/") %*% loadX
    S[!is.finite(S)] <- 0
    colnames(S) <- paste0("comp", seq_len(ncomp))
    S
  }

  scores_train <- project(X_train)
  centroids <- aggregate(as.data.frame(scores_train),
                         by = list(group = Y_train), FUN = mean)

  list(model = model, project = project,
       scores_train = scores_train, centroids = centroids)
}