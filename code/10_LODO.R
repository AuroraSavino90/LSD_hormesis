library(mixOmics)
library(dplyr)
library(tibble)
library(purrr)
library(dplyr)
library(ggplot2)
source("code/NES_sPLSDA.R")

hormesis_label<-"hormesis"
stress_label<-"stress"

train_cols <- names(dataset_type_simpl)[dataset_type_simpl %in% c(hormesis_label, stress_label)]
LSD_cols <- names(dataset_type_simpl)[dataset_type_simpl == "LSD"]

X_train <- t(NES_mat[, train_cols, drop = FALSE])   
X_LSD   <- t(NES_mat[, LSD_cols,   drop = FALSE])

rownames(X_train) <- train_cols
rownames(X_LSD)   <- LSD_cols

Y_train <- factor(dataset_type_simpl[train_cols])


## =========================================================
## fit sPLS-DA and distances LSD from ormesi/stress
## =========================================================

compute_model_lsd_dist <- function(X_train_sub, Y_train_sub, X_LSD,
                                   hormesis_label, stress_label,
                                   ncomp = 2, keepX = rep(50, 2),
                                   scale_train = TRUE) {
  
  X_train_sub <- as.matrix(X_train_sub)
  X_LSD       <- as.matrix(X_LSD)
  
  X_train_sub[!is.finite(X_train_sub)] <- NA_real_
  X_LSD[!is.finite(X_LSD)]             <- NA_real_
  X_train_sub[is.na(X_train_sub)] <- 0
  X_LSD[is.na(X_LSD)]             <- 0
 
    common_feats <- intersect(colnames(X_train_sub), colnames(X_LSD))
  common_feats <- colnames(X_train_sub)[colnames(X_train_sub) %in% common_feats]
  if (length(common_feats) < 10) {
    stop("Too few common features between X_train_sub and X_LSD (n=", length(common_feats), ").")
  }
  X_train_sub <- X_train_sub[, common_feats, drop = FALSE]
  X_LSD       <- X_LSD[,       common_feats, drop = FALSE]
 
   spls_fit <- splsda(X_train_sub, Y_train_sub, ncomp = ncomp, keepX = keepX, scale = scale_train)
  
  
   scores_tr <- as.data.frame(spls_fit$variates$X, check.names = FALSE)
  comp_cols_local <- paste0("comp", seq_len(ncomp))
  colnames(scores_tr)[seq_len(ncomp)] <- comp_cols_local
  scores_tr$train_dataset <- rownames(X_train_sub)
  scores_tr$group         <- Y_train_sub
  
  
   centroids <- scores_tr %>%
    group_by(group) %>%
    summarise(across(all_of(comp_cols_local), mean, na.rm = TRUE),
              .groups = "drop")
  
   center_tr <- colMeans(X_train_sub, na.rm = TRUE)
  
  if (scale_train) {
    sd_tr <- apply(X_train_sub, 2, sd, na.rm = TRUE)
    sd_tr[!is.finite(sd_tr) | sd_tr == 0] <- 1
    X_LSD_scaled <- sweep(sweep(X_LSD, 2, center_tr, "-"), 2, sd_tr, "/")
  } else {
    X_LSD_scaled <- sweep(X_LSD, 2, center_tr, "-")
  }
  
  loadX <- spls_fit$loadings$X[, seq_len(ncomp), drop = FALSE]
  
  scores_LSD_mat <- X_LSD_scaled %*% loadX
  scores_LSD_mat[!is.finite(scores_LSD_mat)] <- 0
  colnames(scores_LSD_mat) <- comp_cols_local
  
  scores_LSD <- as.data.frame(scores_LSD_mat, check.names = FALSE)
  scores_LSD$LSD_dataset <- rownames(X_LSD)
  
  
  dist_to_centroid <- function(point, center) sqrt(sum((point - center)^2))
  
  dist_df <- scores_LSD %>%
    rowwise() %>%
    mutate(
      dist_to_hormesis = {
        c_h <- as.numeric(centroids[centroids$group == hormesis_label,
                                    comp_cols_local, drop = TRUE])
        dist_to_centroid(as.numeric(c_across(all_of(comp_cols_local))), c_h)
      },
      dist_to_stress = {
        c_s <- as.numeric(centroids[centroids$group == stress_label,
                                    comp_cols_local, drop = TRUE])
        dist_to_centroid(as.numeric(c_across(all_of(comp_cols_local))), c_s)
      }
    ) %>%
    ungroup() %>%
    dplyr::select(LSD_dataset, dist_to_hormesis, dist_to_stress)
  
  dist_df
}


ncomp <- 2
keepX <- rep(50, ncomp)

base_dist <- compute_model_lsd_dist(
  X_train_sub    = X_train,
  Y_train_sub    = Y_train,
  X_LSD          = X_LSD,
  hormesis_label = hormesis_label,
  stress_label   = stress_label,
  ncomp          = ncomp,
  keepX          = keepX
)

base_dist


## =========================================================
## Leave-one-out
## =========================================================

train_datasets <- rownames(X_train)

influence_list <- purrr::map(train_datasets, function(ds_rm) {
  keep_ds <- setdiff(train_datasets, ds_rm)
  
  X_tr_sub <- X_train[keep_ds, , drop = FALSE]
  Y_tr_sub <- droplevels(Y_train[keep_ds])
  
  dist_sub <- compute_model_lsd_dist(
    X_train_sub    = X_tr_sub,
    Y_train_sub    = Y_tr_sub,
    X_LSD          = X_LSD,
    hormesis_label = hormesis_label,
    stress_label   = stress_label,
    ncomp          = ncomp,
    keepX          = keepX
  )
  
  dist_sub %>%
    mutate(
      removed_dataset = ds_rm,
      removed_group   = as.character(Y_train[ds_rm])
    )
})

influence_df <- bind_rows(influence_list)

## =========================================================
## Compare with full model
## =========================================================

influence_delta <- influence_df %>%
  left_join(base_dist,
            by = "LSD_dataset",
            suffix = c("_sub", "_base")) %>%
  mutate(
    delta_dist_hormesis = dist_to_hormesis_sub - dist_to_hormesis_base,
    delta_dist_stress   = dist_to_stress_sub   - dist_to_stress_base
  )

influence_summary <- influence_delta %>%
  group_by(removed_dataset, removed_group, LSD_dataset) %>%
  summarise(
    delta_dist_hormesis_mean = mean(delta_dist_hormesis),
    delta_dist_stress_mean   = mean(delta_dist_stress),
    .groups = "drop"
  ) %>%
  arrange(desc(delta_dist_stress_mean))

influence_summary

influence_summary %>%
  filter(LSD_dataset == "DE_24h") %>%
  arrange(desc(delta_dist_stress_mean)) %>%
  head(10)

influence_summary %>%
  filter(LSD_dataset == "DE_90min") %>%
  arrange(desc(delta_dist_stress_mean)) %>%
  head(10)


base_dist_long <- base_dist %>%
  mutate(
    model_id        = "full",
    removed_dataset = NA_character_,
    removed_group   = NA_character_
  )

loo_dist_long <- influence_df %>%
  mutate(
    model_id = paste0("minus_", removed_dataset)
  )

all_models_dist <- bind_rows(base_dist_long, loo_dist_long)

all_models_dist <- all_models_dist %>%
  mutate(
    closer_group = case_when(
      dist_to_hormesis < dist_to_stress ~ "hormesis",
      dist_to_hormesis > dist_to_stress ~ "stress",
      TRUE                              ~ "tie"
    )
  )

all_models_dist


lsd_closeness_summary <- all_models_dist %>%
  group_by(LSD_dataset, closer_group) %>%
  summarise(n_models = n(), .groups = "drop") %>%
  group_by(LSD_dataset) %>%
  mutate(freq = n_models / sum(n_models)) %>%
  arrange(LSD_dataset, desc(freq))

lsd_closeness_summary






dist_wide <- all_models_dist %>%
  filter(model_id != "full") %>%
  select(LSD_dataset, removed_dataset, removed_group,
         dist_to_hormesis, dist_to_stress)

p<-ggplot(dist_wide, aes(x = dist_to_hormesis, y = dist_to_stress)) +
  geom_abline(linetype = "dashed", linewidth = 0.6, color = "grey50") +
  geom_point(alpha = 0.7, size = 2) +
  facet_wrap(~ LSD_dataset) +
  theme_classic(base_size = 13) +
  labs(
    x = "Distance to hormesis",
    y = "Distance to chronic stress",
    title = "Paired LSD distances across leave-one-dataset-out models",
    subtitle = "Each point is one LODO model; dashed line indicates equal distance"
  )


ggsave(
  filename = file.path(results_dir, "figures", "sPLSDA_LODO.png"),
  plot     = p,
  width    = 10,
  height   = 6,
  dpi      = 300
)

ggsave(
  filename = file.path(results_dir, "figures", "sPLSDA_LODO.pdf"),
  plot     = p,
  width    = 10,
  height   = 6
)
