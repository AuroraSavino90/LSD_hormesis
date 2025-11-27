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