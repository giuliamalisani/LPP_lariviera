library(terra)
library(ranger)

# Run the validated notebook workflow first so the experiment uses exactly the
# same Bregaglia samples, spatial folds and untouched assessment fold.
rmd_script <- knitr::purl(
  "data/practical_7.Rmd",
  output = tempfile(fileext = ".R"),
  documentation = 0,
  quiet = TRUE
)
pdf(tempfile(fileext = ".pdf"))
source(rmd_script)
dev.off()

# Extract additional predictors at the already accepted sample locations.
sample_points <- vect(
  samples[, c("x", "y")], geom = c("x", "y"), crs = crs(image_2025)
)

image_2018_bregaglia <- crop(image_2018, bregaglia_satellite, mask = TRUE)
names(image_2018_bregaglia) <- paste0(bands, "_2018")

ndvi_2018 <- ratio_index(image_2018_bregaglia$NIR_2018, image_2018_bregaglia$RED_2018)
ndmi_2018 <- ratio_index(image_2018_bregaglia$NIR_2018, image_2018_bregaglia$SWIR1_2018)
nbr_2018 <- ratio_index(image_2018_bregaglia$NIR_2018, image_2018_bregaglia$SWIR2_2018)
ndwi_2018 <- ratio_index(image_2018_bregaglia$GREEN_2018, image_2018_bregaglia$NIR_2018)
bsi_denominator_2018 <- image_2018_bregaglia$SWIR1_2018 +
  image_2018_bregaglia$RED_2018 + image_2018_bregaglia$NIR_2018 +
  image_2018_bregaglia$BLUE_2018
bsi_2018 <- ifel(
  abs(bsi_denominator_2018) < 1e-06, NA,
  ((image_2018_bregaglia$SWIR1_2018 + image_2018_bregaglia$RED_2018) -
     (image_2018_bregaglia$NIR_2018 + image_2018_bregaglia$BLUE_2018)) /
    bsi_denominator_2018
)
names(ndvi_2018) <- "NDVI_2018"
names(ndmi_2018) <- "NDMI_2018"
names(nbr_2018) <- "NBR_2018"
names(ndwi_2018) <- "NDWI_2018"
names(bsi_2018) <- "BSI_2018"

temporal_stack <- c(
  image_2018_bregaglia, ndvi_2018, ndmi_2018, nbr_2018, ndwi_2018, bsi_2018
)
temporal_values <- extract(temporal_stack, sample_points, ID = FALSE)
samples <- cbind(samples, temporal_values)

index_2025_names <- c("NDVI", "NDMI", "NBR", "NDWI", "BSI")
index_2018_names <- paste0(index_2025_names, "_2018")
delta_band_names <- paste0("D_", bands)
delta_index_names <- paste0("D_", index_2025_names)
for (i in seq_along(bands)) {
  samples[[delta_band_names[i]]] <- samples[[bands[i]]] -
    samples[[paste0(bands[i], "_2018")]]
}
for (i in seq_along(index_2025_names)) {
  samples[[delta_index_names[i]]] <- samples[[index_2025_names[i]]] -
    samples[[index_2018_names[i]]]
}

# Three-pixel neighbourhood summaries add landscape context without using the
# target map or scenario rasters as predictors.
context_source <- c(ndvi, ndwi, bsi, image_2025_bregaglia$NIR,
                    image_2025_bregaglia$SWIR1)
names(context_source) <- c("NDVI", "NDWI", "BSI", "NIR", "SWIR1")
context_mean <- focal(context_source, w = 3, fun = mean, na.rm = TRUE)
context_sd <- focal(context_source, w = 3, fun = sd, na.rm = TRUE)
names(context_mean) <- paste0(names(context_source), "_MEAN3")
names(context_sd) <- paste0(names(context_source), "_SD3")
context_names <- c(names(context_mean), names(context_sd))
samples <- cbind(
  samples,
  extract(c(context_mean, context_sd), sample_points, ID = FALSE)
)

# Refresh the two data frames after attaching the experimental predictors.
development <- samples[samples$spatial_fold != 5, ]
testing <- samples[samples$spatial_fold == 5, ]

temporal_names <- c(
  names(temporal_stack), delta_band_names, delta_index_names
)
feature_sets <- list(
  spectral_2025 = predictor_names,
  temporal = unique(c(predictor_names, temporal_names)),
  context = unique(c(predictor_names, context_names)),
  temporal_context = unique(c(predictor_names, temporal_names, context_names))
)

macro_f1 <- function(observed, predicted) {
  classes <- sort(unique(as.character(observed)))
  observed <- factor(observed, levels = classes)
  predicted <- factor(predicted, levels = classes)
  cm <- table(observed, predicted)
  recall <- diag(cm) / rowSums(cm)
  precision <- ifelse(colSums(cm) > 0, diag(cm) / colSums(cm), 0)
  f1 <- ifelse(
    recall + precision > 0,
    2 * recall * precision / (recall + precision), 0
  )
  mean(f1, na.rm = TRUE)
}

cv_score <- function(features, weighted, settings, trees = 250L) {
  fold_ids <- sort(unique(development$spatial_fold))
  scores <- rep(NA_real_, length(fold_ids))
  for (j in seq_along(fold_ids)) {
    train_index <- development$spatial_fold != fold_ids[j]
    valid_index <- !train_index
    train_complete <- train_index & complete.cases(development[, features, drop = FALSE])
    valid_complete <- valid_index & complete.cases(development[, features, drop = FALSE])
    y_train <- droplevels(development$SSP_CH[train_complete])
    weights <- if (weighted) make_class_weights(y_train) else NULL
    model <- ranger(
      x = development[train_complete, features, drop = FALSE],
      y = y_train,
      num.trees = trees,
      mtry = min(settings$mtry, length(features)),
      min.node.size = settings$min.node.size,
      sample.fraction = settings$sample.fraction,
      splitrule = settings$splitrule,
      class.weights = weights,
      num.threads = threads,
      seed = 700 + j
    )
    prediction <- predict(
      model, development[valid_complete, features, drop = FALSE],
      num.threads = threads
    )$predictions
    scores[j] <- macro_f1(development$SSP_CH[valid_complete], prediction)
  }
  mean(scores)
}

# Stage 1: choose the predictor family and whether class weighting helps,
# without consulting the assessment fold.
stage1 <- expand.grid(
  feature_set = names(feature_sets),
  weighted = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
stage1$cv_macro_f1 <- NA_real_
initial_settings <- list(
  mtry = 3L, min.node.size = 20L,
  sample.fraction = 0.632, splitrule = "extratrees"
)
for (i in seq_len(nrow(stage1))) {
  stage1$cv_macro_f1[i] <- cv_score(
    feature_sets[[stage1$feature_set[i]]], stage1$weighted[i],
    initial_settings
  )
}
stage1 <- stage1[order(stage1$cv_macro_f1, decreasing = TRUE), ]
print(stage1)

selected_features <- feature_sets[[stage1$feature_set[1]]]
selected_weighting <- stage1$weighted[1]

# Stage 2: tune the best predictor family on development folds only.
p_selected <- length(selected_features)
stage2 <- expand.grid(
  mtry = unique(c(2L, floor(sqrt(p_selected)), ceiling(p_selected / 4),
                  ceiling(p_selected / 2), p_selected)),
  min.node.size = c(1L, 5L, 10L, 20L, 30L),
  sample.fraction = c(0.5, 0.632, 0.8),
  splitrule = c("gini", "extratrees"),
  stringsAsFactors = FALSE
)
stage2$cv_macro_f1 <- NA_real_
for (i in seq_len(nrow(stage2))) {
  stage2$cv_macro_f1[i] <- cv_score(
    selected_features, selected_weighting,
    as.list(stage2[i, c("mtry", "min.node.size", "sample.fraction", "splitrule")]),
    trees = 300L
  )
}
stage2 <- stage2[order(stage2$cv_macro_f1, decreasing = TRUE), ]
print(head(stage2, 20))

# One assessment after feature-family and parameter selection.
chosen <- stage2[1, ]
development_complete <- complete.cases(development[, selected_features, drop = FALSE])
testing_complete <- complete.cases(testing[, selected_features, drop = FALSE])
y_train <- droplevels(development$SSP_CH[development_complete])
final_weights <- if (selected_weighting) make_class_weights(y_train) else NULL
experiment_model <- ranger(
  x = development[development_complete, selected_features, drop = FALSE],
  y = y_train,
  num.trees = 1500,
  mtry = chosen$mtry,
  min.node.size = chosen$min.node.size,
  sample.fraction = chosen$sample.fraction,
  splitrule = chosen$splitrule,
  class.weights = final_weights,
  importance = "permutation",
  num.threads = threads,
  seed = 909
)
experiment_prediction <- predict(
  experiment_model,
  testing[testing_complete, selected_features, drop = FALSE],
  num.threads = threads
)$predictions
experiment_observed <- testing$SSP_CH[testing_complete]
experiment_results <- data.frame(
  selected_feature_set = stage1$feature_set[1],
  weighted = selected_weighting,
  n_development = sum(development_complete),
  n_test = sum(testing_complete),
  overall_accuracy_pct = 100 * mean(experiment_prediction == experiment_observed),
  balanced_accuracy_pct = 100 * balanced_accuracy(
    experiment_observed, experiment_prediction
  ),
  macro_f1_pct = 100 * macro_f1(experiment_observed, experiment_prediction),
  mtry = chosen$mtry,
  min.node.size = chosen$min.node.size,
  sample.fraction = chosen$sample.fraction,
  splitrule = chosen$splitrule
)
print(experiment_results)
print(sort(experiment_model$variable.importance, decreasing = TRUE))

write.csv(stage1, "data/bregaglia_rf_feature_screen.csv", row.names = FALSE)
write.csv(stage2, "data/bregaglia_rf_tuning_screen.csv", row.names = FALSE)
write.csv(experiment_results, "data/bregaglia_rf_experiment_result.csv", row.names = FALSE)
