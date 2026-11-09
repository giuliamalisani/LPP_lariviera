library(terra)
library(ranger)

# Reproduce the tuned nine-class Bregaglia terrain model and its fixed spatial
# holdout before comparing alternative class abstractions.
rmd_script <- knitr::purl(
  "data/practical_7.Rmd", output = tempfile(fileext = ".R"),
  documentation = 0, quiet = TRUE
)
pdf(tempfile(fileext = ".pdf"))
source(rmd_script)
dev.off()

# Six defensible higher-level groups. Artificial remains separate from the
# existing Static class, which already mixes bare ground with some transport
# and infrastructure surfaces in the source legend.
fine_to_coarse <- c(
  "10" = 1L,  # Artificial/amenities
  "11" = 2L,  # Static/bare
  "12" = 3L, "13" = 3L,  # Forest
  "16" = 4L, "17" = 4L,  # Grassland/pasture
  "19" = 5L,              # Glacier/snow
  "20" = 6L, "21" = 6L   # Water
)
coarse_names <- c(
  "Artificial/amenities", "Static/bare", "Forest",
  "Grassland/pasture", "Glacier/snow", "Water"
)
coarse_levels <- as.character(seq_along(coarse_names))

to_coarse <- function(x) {
  factor(
    unname(fine_to_coarse[as.character(x)]),
    levels = seq_along(coarse_names), labels = coarse_names
  )
}

classification_metrics <- function(observed, predicted, class_levels) {
  observed <- factor(observed, levels = class_levels)
  predicted <- factor(predicted, levels = class_levels)
  cm <- table(observed, predicted)
  recall <- diag(cm) / rowSums(cm)
  precision <- ifelse(colSums(cm) > 0, diag(cm) / colSums(cm), 0)
  f1 <- ifelse(
    recall + precision > 0,
    2 * recall * precision / (recall + precision), 0
  )
  assessed <- rowSums(cm) > 0
  list(
    confusion = cm,
    summary = c(
      overall_accuracy_pct = 100 * sum(diag(cm)) / sum(cm),
      balanced_accuracy_pct = 100 * mean(recall[assessed]),
      macro_f1_pct = 100 * mean(f1[assessed])
    ),
    per_class = data.frame(
      class = class_levels[assessed],
      n_test = rowSums(cm)[assessed],
      recall_pct = 100 * recall[assessed],
      precision_pct = 100 * precision[assessed],
      f1_pct = 100 * f1[assessed]
    )
  )
}

# Test A: retain the nine-class model, then aggregate both observations and
# predictions after classification. This measures how much fine-level error is
# confusion within the same broader group.
fine_observed <- testing$SSP_CH
fine_prediction <- test_prediction
fine_metrics <- classification_metrics(
  fine_observed, fine_prediction, levels(samples$SSP_CH)
)
posthoc_observed <- to_coarse(fine_observed)
posthoc_prediction <- to_coarse(fine_prediction)
posthoc_metrics <- classification_metrics(
  posthoc_observed, posthoc_prediction, coarse_names
)

# Test B: recode labels first and tune a new forest specifically for the six
# broader classes, using the same development folds and untouched fold 5.
development$COARSE <- to_coarse(development$SSP_CH)
testing$COARSE <- to_coarse(testing$SSP_CH)

coarse_weights <- function(y) {
  counts <- table(droplevels(y))
  weights <- sum(counts) / (length(counts) * counts)
  setNames(as.numeric(weights), names(counts))
}

coarse_macro_f1 <- function(observed, predicted) {
  classification_metrics(observed, predicted, coarse_names)$summary[
    "macro_f1_pct"
  ] / 100
}

p <- length(predictor_names)
coarse_grid <- expand.grid(
  mtry = unique(c(2L, floor(sqrt(p)), ceiling(p / 2), p)),
  min.node.size = c(1L, 5L, 10L, 20L),
  sample.fraction = c(0.632, 0.8),
  splitrule = c("gini", "extratrees"),
  weighted = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
coarse_grid$cv_macro_f1 <- NA_real_
tuning_folds <- sort(unique(development$spatial_fold))
for (g in seq_len(nrow(coarse_grid))) {
  fold_scores <- numeric(length(tuning_folds))
  for (j in seq_along(tuning_folds)) {
    validation_fold <- tuning_folds[j]
    train_index <- development$spatial_fold != validation_fold
    validation_index <- !train_index
    y_train <- droplevels(development$COARSE[train_index])
    weights <- if (coarse_grid$weighted[g]) coarse_weights(y_train) else NULL
    model <- ranger(
      x = development[train_index, predictor_names, drop = FALSE],
      y = y_train,
      num.trees = 250,
      mtry = coarse_grid$mtry[g],
      min.node.size = coarse_grid$min.node.size[g],
      sample.fraction = coarse_grid$sample.fraction[g],
      splitrule = coarse_grid$splitrule[g],
      class.weights = weights,
      num.threads = threads,
      seed = 3000 + g + j
    )
    prediction <- predict(
      model,
      development[validation_index, predictor_names, drop = FALSE],
      num.threads = threads
    )$predictions
    fold_scores[j] <- coarse_macro_f1(
      development$COARSE[validation_index], prediction
    )
  }
  coarse_grid$cv_macro_f1[g] <- mean(fold_scores)
}
coarse_grid <- coarse_grid[order(coarse_grid$cv_macro_f1, decreasing = TRUE), ]
coarse_best <- coarse_grid[1, ]

coarse_development_y <- droplevels(development$COARSE)
assessment_weights <- if (coarse_best$weighted) {
  coarse_weights(coarse_development_y)
} else {
  NULL
}
coarse_assessment_model <- ranger(
  x = development[, predictor_names, drop = FALSE],
  y = coarse_development_y,
  num.trees = 1000,
  mtry = coarse_best$mtry,
  min.node.size = coarse_best$min.node.size,
  sample.fraction = coarse_best$sample.fraction,
  splitrule = coarse_best$splitrule,
  class.weights = assessment_weights,
  importance = "permutation",
  num.threads = threads,
  seed = 3200
)
coarse_prediction <- predict(
  coarse_assessment_model,
  testing[, predictor_names, drop = FALSE],
  num.threads = threads
)$predictions
coarse_metrics <- classification_metrics(
  testing$COARSE, coarse_prediction, coarse_names
)

hierarchy_comparison <- rbind(
  data.frame(
    test = "Fine model, fine evaluation",
    model_classes = 9L, evaluation_classes = 9L,
    overall_accuracy_pct = fine_metrics$summary["overall_accuracy_pct"],
    balanced_accuracy_pct = fine_metrics$summary["balanced_accuracy_pct"],
    macro_f1_pct = fine_metrics$summary["macro_f1_pct"]
  ),
  data.frame(
    test = "Fine model, post-hoc coarse evaluation",
    model_classes = 9L, evaluation_classes = 6L,
    overall_accuracy_pct = posthoc_metrics$summary["overall_accuracy_pct"],
    balanced_accuracy_pct = posthoc_metrics$summary["balanced_accuracy_pct"],
    macro_f1_pct = posthoc_metrics$summary["macro_f1_pct"]
  ),
  data.frame(
    test = "Coarse model, coarse evaluation",
    model_classes = 6L, evaluation_classes = 6L,
    overall_accuracy_pct = coarse_metrics$summary["overall_accuracy_pct"],
    balanced_accuracy_pct = coarse_metrics$summary["balanced_accuracy_pct"],
    macro_f1_pct = coarse_metrics$summary["macro_f1_pct"]
  )
)
hierarchy_comparison[, c(
  "overall_accuracy_pct", "balanced_accuracy_pct", "macro_f1_pct"
)] <- round(hierarchy_comparison[, c(
  "overall_accuracy_pct", "balanced_accuracy_pct", "macro_f1_pct"
)], 1)

posthoc_class_metrics <- posthoc_metrics$per_class
posthoc_class_metrics$test <- "Post-hoc aggregation"
coarse_class_metrics <- coarse_metrics$per_class
coarse_class_metrics$test <- "Coarse retraining"
hierarchy_by_class <- rbind(posthoc_class_metrics, coarse_class_metrics)
hierarchy_by_class[, c("recall_pct", "precision_pct", "f1_pct")] <- round(
  hierarchy_by_class[, c("recall_pct", "precision_pct", "f1_pct")], 1
)

print(hierarchy_comparison)
print(hierarchy_by_class)
print(head(coarse_grid, 20))
print(fine_metrics$confusion)
print(posthoc_metrics$confusion)
print(coarse_metrics$confusion)

# Write a post-hoc six-class map by recomposing the existing nine-class output.
prediction_bregaglia_posthoc <- subst(
  prediction_bregaglia,
  from = as.integer(names(fine_to_coarse)),
  to = as.integer(fine_to_coarse)
)
coarse_palette <- data.frame(
  value = seq_along(coarse_names), class = coarse_names,
  colour = c("#d73027", "#bdbdbd", "#1b7837", "#a6d96a", "#80b1d3", "#2166ac")
)
levels(prediction_bregaglia_posthoc) <- coarse_palette[, c("value", "class")]
writeRaster(
  prediction_bregaglia_posthoc,
  file.path(data_dir, "landcover_RF_Bregaglia_2025_posthoc_6class.tif"),
  overwrite = TRUE, datatype = "INT1U"
)

# Fit the selected coarse model on every eligible Bregaglia sample and write a
# separate map. This does not replace the nine-class or post-hoc products.
samples$COARSE <- to_coarse(samples$SSP_CH)
coarse_all_y <- droplevels(samples$COARSE)
production_weights <- if (coarse_best$weighted) {
  coarse_weights(coarse_all_y)
} else {
  NULL
}
coarse_model <- ranger(
  x = samples[, predictor_names, drop = FALSE],
  y = coarse_all_y,
  num.trees = 1000,
  mtry = coarse_best$mtry,
  min.node.size = coarse_best$min.node.size,
  sample.fraction = coarse_best$sample.fraction,
  splitrule = coarse_best$splitrule,
  class.weights = production_weights,
  importance = "permutation",
  num.threads = threads,
  seed = 3300
)
predict_coarse <- function(model, data) {
  prediction <- predict(
    model, data = data, num.threads = threads
  )$predictions
  match(as.character(prediction), coarse_names)
}
prediction_bregaglia_coarse <- terra::predict(
  predictors_bregaglia, coarse_model, fun = predict_coarse, na.rm = TRUE,
  filename = file.path(
    data_dir, "landcover_RF_Bregaglia_2025_coarse_6class.tif"
  ),
  overwrite = TRUE, wopt = list(datatype = "INT1U")
)
levels(prediction_bregaglia_coarse) <- coarse_palette[, c("value", "class")]

write.csv(
  hierarchy_comparison,
  file.path(data_dir, "bregaglia_rf_hierarchy_comparison.csv"),
  row.names = FALSE
)
write.csv(
  hierarchy_by_class,
  file.path(data_dir, "bregaglia_rf_hierarchy_by_class.csv"),
  row.names = FALSE
)
write.csv(
  coarse_grid,
  file.path(data_dir, "bregaglia_rf_coarse_tuning.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame.matrix(fine_metrics$confusion),
  file.path(data_dir, "confusion_fine_9class.csv")
)
write.csv(
  as.data.frame.matrix(posthoc_metrics$confusion),
  file.path(data_dir, "confusion_posthoc_6class.csv")
)
write.csv(
  as.data.frame.matrix(coarse_metrics$confusion),
  file.path(data_dir, "confusion_retrained_6class.csv")
)
