library(terra)
library(ranger)

# Establish the local samples, fixed spatial folds, spectral predictors and
# native-resolution Bregaglia terrain predictors from the main workflow.
rmd_script <- knitr::purl(
  "data/practical_7.Rmd", output = tempfile(fileext = ".R"),
  documentation = 0, quiet = TRUE
)
pdf(tempfile(fileext = ".pdf"))
source(rmd_script)
dev.off()

metric_bundle <- function(observed, predicted) {
  classes <- levels(samples$SSP_CH)
  observed <- factor(observed, levels = classes)
  predicted <- factor(predicted, levels = classes)
  cm <- table(observed, predicted)
  recall <- diag(cm) / rowSums(cm)
  precision <- ifelse(colSums(cm) > 0, diag(cm) / colSums(cm), 0)
  f1 <- ifelse(
    recall + precision > 0,
    2 * recall * precision / (recall + precision), 0
  )
  assessed <- rowSums(cm) > 0
  list(
    summary = c(
      overall_accuracy_pct = 100 * sum(diag(cm)) / sum(cm),
      balanced_accuracy_pct = 100 * mean(recall[assessed]),
      macro_f1_pct = 100 * mean(f1[assessed])
    ),
    per_class = data.frame(
      class = classes[assessed], n_test = rowSums(cm)[assessed],
      recall_pct = 100 * recall[assessed],
      precision_pct = 100 * precision[assessed],
      f1_pct = 100 * f1[assessed]
    )
  )
}

comparison_rows <- list()
per_class_rows <- list()
add_result <- function(model_name, training_scope, predictors, tuned,
                       n_training, observed, predicted) {
  metrics <- metric_bundle(observed, predicted)
  comparison_rows[[model_name]] <<- data.frame(
    model = model_name,
    training_scope = training_scope,
    predictors = predictors,
    tuned = tuned,
    n_training = n_training,
    n_test = length(observed),
    overall_accuracy_pct = metrics$summary["overall_accuracy_pct"],
    balanced_accuracy_pct = metrics$summary["balanced_accuracy_pct"],
    macro_f1_pct = metrics$summary["macro_f1_pct"],
    row.names = NULL
  )
  metrics$per_class$model <- model_name
  per_class_rows[[model_name]] <<- metrics$per_class
}

test_observed <- testing$SSP_CH

# 1. Original-style local RF: bands only, fixed defaults from the initial
# notebook, assessed on the same spatial holdout as every model below.
basic_local <- ranger(
  x = development[, bands, drop = FALSE],
  y = droplevels(development$SSP_CH),
  num.trees = 500, min.node.size = 10,
  num.threads = threads, seed = 2001
)
basic_local_prediction <- predict(
  basic_local, testing[, bands, drop = FALSE], num.threads = threads
)$predictions
add_result(
  "Basic Bregaglia", "Bregaglia development blocks", "10 bands", "No",
  nrow(development), test_observed, basic_local_prediction
)

# 2. Tuned local spectral model, without terrain, for an isolated estimate of
# the gain from indices and RF tuning.
spectral_names <- c(bands, "NDVI", "NDMI", "NBR", "NDWI", "BSI")
tuned_local_spectral <- ranger(
  x = development[, spectral_names, drop = FALSE],
  y = droplevels(development$SSP_CH),
  num.trees = 1000, mtry = 3, min.node.size = 20,
  sample.fraction = 0.632, splitrule = "extratrees",
  class.weights = make_class_weights(development$SSP_CH),
  num.threads = threads, seed = 2002
)
tuned_local_spectral_prediction <- predict(
  tuned_local_spectral, testing[, spectral_names, drop = FALSE],
  num.threads = threads
)$predictions
add_result(
  "Tuned Bregaglia spectral", "Bregaglia development blocks",
  "10 bands + 5 indices", "Yes", nrow(development), test_observed,
  tuned_local_spectral_prediction
)

# 3. Main tuned local terrain model already selected without using fold 5.
add_result(
  "Tuned Bregaglia terrain", "Bregaglia development blocks",
  "Bands + indices + DHM25 terrain", "Yes", nrow(development),
  test_observed, test_prediction
)

# Create and cache a national 100 m terrain stack. National terrain is derived
# at 100 m to keep the full-Switzerland computation tractable.
terrain_ch_path <- file.path(
  data_dir, "Data_earlierass", "dhm25_terrain_100m_EPSG3035.tif"
)
if (!file.exists(terrain_ch_path)) {
  dem_ch <- rast(file.path(
    data_dir, "Data_earlierass", "dhm25_grid_raster.asc"
  ))
  dem_ch <- shift(dem_ch, dx = 2000000, dy = 1000000)
  crs(dem_ch) <- "EPSG:2056"
  dem_ch_100 <- aggregate(dem_ch, fact = 4, fun = mean, na.rm = TRUE)
  names(dem_ch_100) <- "elevation"
  terrain_ch_native <- terrain(
    dem_ch_100,
    v = c("slope", "aspect", "TPI", "TRI", "roughness"),
    unit = "degrees", neighbors = 8
  )
  northness_ch <- cos(terrain_ch_native$aspect * pi / 180)
  eastness_ch <- sin(terrain_ch_native$aspect * pi / 180)
  terrain_ch_native <- c(
    dem_ch_100, terrain_ch_native$slope, northness_ch, eastness_ch,
    terrain_ch_native$TPI, terrain_ch_native$TRI,
    terrain_ch_native$roughness
  )
  names(terrain_ch_native) <- c(
    "elevation", "slope", "northness", "eastness",
    "TPI", "TRI", "roughness"
  )
  terrain_ch <- project(
    terrain_ch_native, image_2025, method = "bilinear",
    filename = terrain_ch_path, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  )
} else {
  terrain_ch <- rast(terrain_ch_path)
}

# Build the national predictor table with exactly the same class legend.
ndvi_ch <- ratio_index(image_2025$NIR, image_2025$RED)
ndmi_ch <- ratio_index(image_2025$NIR, image_2025$SWIR1)
nbr_ch <- ratio_index(image_2025$NIR, image_2025$SWIR2)
ndwi_ch <- ratio_index(image_2025$GREEN, image_2025$NIR)
bsi_denominator_ch <- image_2025$SWIR1 + image_2025$RED +
  image_2025$NIR + image_2025$BLUE
bsi_ch <- ifel(
  abs(bsi_denominator_ch) < 1e-06, NA,
  ((image_2025$SWIR1 + image_2025$RED) -
     (image_2025$NIR + image_2025$BLUE)) / bsi_denominator_ch
)
names(ndvi_ch) <- "NDVI"
names(ndmi_ch) <- "NDMI"
names(nbr_ch) <- "NBR"
names(ndwi_ch) <- "NDWI"
names(bsi_ch) <- "BSI"
national_predictors <- c(
  image_2025, ndvi_ch, ndmi_ch, nbr_ch, ndwi_ch, bsi_ch, terrain_ch
)
names(national_predictors) <- c(spectral_names, names(terrain_ch))

national_xy <- crds(reference)
national_key <- paste(round(national_xy[, 1], 3), round(national_xy[, 2], 3))
bregaglia_key <- paste(round(samples$x, 3), round(samples$y, 3))
national_class <- ssp_lookup[as.integer(as.character(reference$AS_27))]
national_data <- data.frame(
  SSP_CH = factor(national_class, levels = levels(samples$SSP_CH)),
  key = national_key,
  extract(national_predictors, reference, ID = FALSE)
)
national_data <- na.omit(national_data)
national_outside <- national_data[!national_data$key %in% bregaglia_key, ]

# A capped stratified sample preserves the whole Swiss spatial domain without
# allowing very common classes to dominate computation or voting.
set.seed(2100)
national_index <- unlist(lapply(
  split(seq_len(nrow(national_outside)), national_outside$SSP_CH),
  function(i) sample(i, min(length(i), 5000L))
))
national_train <- droplevels(national_outside[national_index, ])

# 4. Basic national transfer model: no Bregaglia labels in training.
basic_ch <- ranger(
  x = national_train[, bands, drop = FALSE],
  y = national_train$SSP_CH,
  num.trees = 500, min.node.size = 10,
  num.threads = threads, seed = 2101
)
basic_ch_prediction <- predict(
  basic_ch, testing[, bands, drop = FALSE], num.threads = threads
)$predictions
add_result(
  "Basic CH transfer", "CH-wide, Bregaglia excluded", "10 bands", "No",
  nrow(national_train), test_observed, basic_ch_prediction
)

# 5. Tune national terrain transfer using Bregaglia development observations
# only as validation. They are not part of national_train.
ch_grid <- expand.grid(
  mtry = c(4L, 11L, 22L), min.node.size = c(5L, 20L),
  sample.fraction = c(0.632, 0.8),
  splitrule = "extratrees", stringsAsFactors = FALSE
)
ch_grid$development_macro_f1 <- NA_real_
for (i in seq_len(nrow(ch_grid))) {
  ch_candidate <- ranger(
    x = national_train[, predictor_names, drop = FALSE],
    y = national_train$SSP_CH,
    num.trees = 300,
    mtry = ch_grid$mtry[i],
    min.node.size = ch_grid$min.node.size[i],
    sample.fraction = ch_grid$sample.fraction[i],
    splitrule = ch_grid$splitrule[i],
    class.weights = make_class_weights(national_train$SSP_CH),
    num.threads = threads, seed = 2200 + i
  )
  candidate_prediction <- predict(
    ch_candidate, development[, predictor_names, drop = FALSE],
    num.threads = threads
  )$predictions
  ch_grid$development_macro_f1[i] <-
    metric_bundle(development$SSP_CH, candidate_prediction)$summary["macro_f1_pct"]
}
ch_grid <- ch_grid[order(ch_grid$development_macro_f1, decreasing = TRUE), ]
ch_best <- ch_grid[1, ]

tuned_ch <- ranger(
  x = national_train[, predictor_names, drop = FALSE],
  y = national_train$SSP_CH,
  num.trees = 1000,
  mtry = ch_best$mtry,
  min.node.size = ch_best$min.node.size,
  sample.fraction = ch_best$sample.fraction,
  splitrule = ch_best$splitrule,
  class.weights = make_class_weights(national_train$SSP_CH),
  num.threads = threads, seed = 2300
)
tuned_ch_prediction <- predict(
  tuned_ch, testing[, predictor_names, drop = FALSE], num.threads = threads
)$predictions
add_result(
  "Tuned CH terrain transfer", "CH-wide, Bregaglia excluded",
  "Bands + indices + DHM25 terrain", "Yes", nrow(national_train),
  test_observed, tuned_ch_prediction
)

# 6. Hybrid model: national coverage plus all non-holdout Bregaglia blocks.
hybrid_train <- rbind(
  national_train[, c("SSP_CH", predictor_names)],
  development[, c("SSP_CH", predictor_names)]
)
hybrid_train$SSP_CH <- factor(
  hybrid_train$SSP_CH, levels = levels(samples$SSP_CH)
)
hybrid <- ranger(
  x = hybrid_train[, predictor_names, drop = FALSE],
  y = hybrid_train$SSP_CH,
  num.trees = 1000,
  mtry = ch_best$mtry,
  min.node.size = ch_best$min.node.size,
  sample.fraction = ch_best$sample.fraction,
  splitrule = ch_best$splitrule,
  class.weights = make_class_weights(hybrid_train$SSP_CH),
  case.weights = c(rep(1, nrow(national_train)), rep(5, nrow(development))),
  num.threads = threads, seed = 2400
)
hybrid_prediction <- predict(
  hybrid, testing[, predictor_names, drop = FALSE], num.threads = threads
)$predictions
add_result(
  "Hybrid CH + Bregaglia terrain",
  "CH-wide plus weighted Bregaglia development blocks",
  "Bands + indices + DHM25 terrain", "Yes", nrow(hybrid_train),
  test_observed, hybrid_prediction
)

comparison <- do.call(rbind, comparison_rows)
rownames(comparison) <- NULL
comparison[, c(
  "overall_accuracy_pct", "balanced_accuracy_pct", "macro_f1_pct"
)] <- round(comparison[, c(
  "overall_accuracy_pct", "balanced_accuracy_pct", "macro_f1_pct"
)], 1)
comparison <- comparison[order(comparison$macro_f1_pct, decreasing = TRUE), ]
per_class_comparison <- do.call(rbind, per_class_rows)
rownames(per_class_comparison) <- NULL
per_class_comparison[, c("recall_pct", "precision_pct", "f1_pct")] <- round(
  per_class_comparison[, c("recall_pct", "precision_pct", "f1_pct")], 1
)

print(comparison)
print(per_class_comparison)
print(ch_grid)
write.csv(comparison, "data/bregaglia_rf_model_comparison.csv", row.names = FALSE)
write.csv(
  per_class_comparison,
  "data/bregaglia_rf_model_comparison_by_class.csv", row.names = FALSE
)
write.csv(ch_grid, "data/bregaglia_rf_ch_tuning.csv", row.names = FALSE)
