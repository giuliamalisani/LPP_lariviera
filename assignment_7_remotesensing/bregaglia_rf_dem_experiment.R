library(terra)
library(ranger)

# Reuse exactly the same samples, spatial blocks and assessment fold as the
# main Bregaglia workflow.
rmd_script <- knitr::purl(
  "data/practical_7.Rmd", output = tempfile(fileext = ".R"),
  documentation = 0, quiet = TRUE
)
pdf(tempfile(fileext = ".pdf"))
source(rmd_script)
dev.off()

# DHM25 is supplied without embedded CRS metadata. Its coordinate range is the
# legacy Swiss LV03 grid. The local PROJ installation lacks the optional Swiss
# datum grid needed for an LV03-to-LV95 transformation, so express the grid in
# LV95 by applying its 2,000,000/1,000,000 m false-origin shift. Residual LV03
# distortion is well below the 25 m DEM and 100 m modelling resolutions.
dem <- rast("data/Data_earlierass/dhm25_grid_raster.asc")
dem <- shift(dem, dx = 2000000, dy = 1000000)
crs(dem) <- "EPSG:2056"
bregaglia_lv95 <- project(bregaglia, crs(dem))
dem_local <- crop(dem, buffer(bregaglia_lv95, width = 1000), mask = TRUE)
names(dem_local) <- "elevation"

terrain_native <- terrain(
  dem_local,
  v = c("slope", "aspect", "TPI", "TRI", "roughness"),
  unit = "degrees",
  neighbors = 8
)
northness <- cos(terrain_native$aspect * pi / 180)
eastness <- sin(terrain_native$aspect * pi / 180)
names(northness) <- "northness"
names(eastness) <- "eastness"

terrain_native <- c(
  dem_local, terrain_native$slope, northness, eastness,
  terrain_native$TPI, terrain_native$TRI, terrain_native$roughness
)
names(terrain_native) <- c(
  "elevation", "slope", "northness", "eastness",
  "TPI", "TRI", "roughness"
)
terrain_bregaglia <- project(
  terrain_native, image_2025_bregaglia, method = "bilinear"
)

sample_points <- vect(
  samples[, c("x", "y")], geom = c("x", "y"), crs = crs(image_2025)
)
samples <- cbind(
  samples, extract(terrain_bregaglia, sample_points, ID = FALSE)
)
development <- samples[samples$spatial_fold != 5, ]
testing <- samples[samples$spatial_fold == 5, ]

terrain_names <- names(terrain_bregaglia)
feature_sets <- list(
  spectral_2025 = predictor_names,
  spectral_terrain = c(predictor_names, terrain_names)
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

cv_score <- function(features, weighted, settings, trees = 300L) {
  fold_ids <- sort(unique(development$spatial_fold))
  balanced_scores <- macro_scores <- numeric(length(fold_ids))
  for (j in seq_along(fold_ids)) {
    train_index <- development$spatial_fold != fold_ids[j]
    valid_index <- !train_index
    train_index <- train_index & complete.cases(development[, features, drop = FALSE])
    valid_index <- valid_index & complete.cases(development[, features, drop = FALSE])
    y_train <- droplevels(development$SSP_CH[train_index])
    weights <- if (weighted) make_class_weights(y_train) else NULL
    model <- ranger(
      x = development[train_index, features, drop = FALSE],
      y = y_train,
      num.trees = trees,
      mtry = min(settings$mtry, length(features)),
      min.node.size = settings$min.node.size,
      sample.fraction = settings$sample.fraction,
      splitrule = settings$splitrule,
      class.weights = weights,
      num.threads = threads,
      seed = 1200 + j
    )
    prediction <- predict(
      model, development[valid_index, features, drop = FALSE],
      num.threads = threads
    )$predictions
    balanced_scores[j] <- balanced_accuracy(
      development$SSP_CH[valid_index], prediction
    )
    macro_scores[j] <- macro_f1(development$SSP_CH[valid_index], prediction)
  }
  c(
    balanced_accuracy = mean(balanced_scores),
    macro_f1 = mean(macro_scores)
  )
}

# Screen whether terrain helps before the larger parameter search.
screen <- expand.grid(
  feature_set = names(feature_sets), weighted = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
screen$cv_balanced_accuracy <- screen$cv_macro_f1 <- NA_real_
initial_settings <- list(
  mtry = 3L, min.node.size = 20L,
  sample.fraction = 0.632, splitrule = "extratrees"
)
for (i in seq_len(nrow(screen))) {
  score <- cv_score(
    feature_sets[[screen$feature_set[i]]], screen$weighted[i],
    initial_settings
  )
  screen$cv_balanced_accuracy[i] <- score["balanced_accuracy"]
  screen$cv_macro_f1[i] <- score["macro_f1"]
}
screen <- screen[order(screen$cv_macro_f1, decreasing = TRUE), ]
print(screen)

selected_features <- feature_sets[[screen$feature_set[1]]]
selected_weighting <- screen$weighted[1]
p_selected <- length(selected_features)
tuning <- expand.grid(
  mtry = unique(c(2L, floor(sqrt(p_selected)), ceiling(p_selected / 4),
                  ceiling(p_selected / 2), p_selected)),
  min.node.size = c(1L, 5L, 10L, 20L, 30L),
  sample.fraction = c(0.5, 0.632, 0.8),
  splitrule = c("gini", "extratrees"),
  stringsAsFactors = FALSE
)
tuning$cv_balanced_accuracy <- tuning$cv_macro_f1 <- NA_real_
for (i in seq_len(nrow(tuning))) {
  settings <- as.list(tuning[i, c(
    "mtry", "min.node.size", "sample.fraction", "splitrule"
  )])
  score <- cv_score(
    selected_features, selected_weighting, settings, trees = 300L
  )
  tuning$cv_balanced_accuracy[i] <- score["balanced_accuracy"]
  tuning$cv_macro_f1[i] <- score["macro_f1"]
}
tuning <- tuning[order(tuning$cv_macro_f1, decreasing = TRUE), ]
print(head(tuning, 20))

chosen <- tuning[1, ]
development_complete <- complete.cases(development[, selected_features, drop = FALSE])
testing_complete <- complete.cases(testing[, selected_features, drop = FALSE])
y_train <- droplevels(development$SSP_CH[development_complete])
final_weights <- if (selected_weighting) make_class_weights(y_train) else NULL
dem_model <- ranger(
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
  seed = 1400
)
dem_prediction <- predict(
  dem_model, testing[testing_complete, selected_features, drop = FALSE],
  num.threads = threads
)$predictions
dem_observed <- testing$SSP_CH[testing_complete]
dem_result <- data.frame(
  selected_feature_set = screen$feature_set[1],
  weighted = selected_weighting,
  n_development = sum(development_complete),
  n_test = sum(testing_complete),
  overall_accuracy_pct = 100 * mean(dem_prediction == dem_observed),
  balanced_accuracy_pct = 100 * balanced_accuracy(dem_observed, dem_prediction),
  macro_f1_pct = 100 * macro_f1(dem_observed, dem_prediction),
  mtry = chosen$mtry,
  min.node.size = chosen$min.node.size,
  sample.fraction = chosen$sample.fraction,
  splitrule = chosen$splitrule
)
print(dem_result)
print(sort(dem_model$variable.importance, decreasing = TRUE))

write.csv(screen, "data/bregaglia_rf_dem_feature_screen.csv", row.names = FALSE)
write.csv(tuning, "data/bregaglia_rf_dem_tuning_screen.csv", row.names = FALSE)
write.csv(dem_result, "data/bregaglia_rf_dem_result.csv", row.names = FALSE)
