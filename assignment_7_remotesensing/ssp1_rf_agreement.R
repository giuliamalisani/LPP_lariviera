library(terra)

data_dir <- "."
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

ssp1 <- rast(file.path(
  data_dir, "NCCS-SSP-results", "raster_data", "SSP1", "SSP1-2025.tif"
))
national_full <- rast(file.path(data_dir, "landcover_RF_SSPCH_2025.tif"))
local_rf <- rast(file.path(
  data_dir, "landcover_RF_SSPCH_Bregaglia_2025.tif"
))
boundary <- vect(file.path(
  data_dir, "communes", "square_20x20km_Bregaglia_commune.gpkg"
))
boundary_rf <- project(boundary, crs(local_rf))

# Use the local 100-m grid as the common geometry. The scenario map is a
# categorical raster, so nearest-neighbour projection preserves class codes.
ssp1_local <- project(ssp1, local_rf, method = "near")
ssp1_local <- mask(ssp1_local, boundary_rf)
national_local <- crop(national_full, local_rf)
national_local <- mask(national_local, boundary_rf)

class_names <- c(
  "10" = "Urban/amenities", "11" = "Static", "12" = "Open Forest",
  "13" = "Closed forest", "14" = "Overgrown/shrubland",
  "15" = "Agriculture", "16" = "Alpine pastures",
  "17" = "Grassland or meadows", "18" = "Permanent crops",
  "19" = "Glacier", "20" = "Lake", "21" = "River"
)
class_levels <- as.integer(names(class_names))

agreement_metrics <- function(reference_map, comparison_map, model_name,
                              supported_classes) {
  paired <- na.omit(cbind(
    ssp1 = values(reference_map, mat = FALSE),
    rf = values(comparison_map, mat = FALSE)
  ))
  cm <- table(
    SSP1 = factor(paired[, "ssp1"], levels = class_levels),
    RF = factor(paired[, "rf"], levels = class_levels)
  )

  true_positive <- diag(cm)
  reference_total <- rowSums(cm)
  comparison_total <- colSums(cm)
  recall <- ifelse(reference_total > 0, true_positive / reference_total, NA)
  precision <- ifelse(comparison_total > 0, true_positive / comparison_total, NA)
  iou_denominator <- reference_total + comparison_total - true_positive
  iou <- ifelse(iou_denominator > 0, true_positive / iou_denominator, NA)

  observed_agreement <- sum(true_positive) / sum(cm)
  expected_agreement <- sum(reference_total * comparison_total) / sum(cm)^2
  kappa <- (observed_agreement - expected_agreement) /
    (1 - expected_agreement)

  supported <- paired[, "ssp1"] %in% supported_classes
  summary <- data.frame(
    model = model_name,
    n_compared_cells = nrow(paired),
    overall_agreement_pct = 100 * observed_agreement,
    balanced_agreement_pct = 100 * mean(recall, na.rm = TRUE),
    macro_iou_pct = 100 * mean(iou, na.rm = TRUE),
    cohen_kappa = kappa,
    ssp1_cells_in_unsupported_classes_pct = 100 * mean(!supported),
    agreement_within_supported_ssp1_classes_pct =
      100 * mean(paired[supported, "ssp1"] == paired[supported, "rf"])
  )

  by_class <- data.frame(
    model = model_name,
    class_value = class_levels,
    class_name = unname(class_names[as.character(class_levels)]),
    n_ssp1 = as.integer(reference_total),
    n_rf = as.integer(comparison_total),
    agreement_recall_pct = 100 * recall,
    agreement_precision_pct = 100 * precision,
    iou_pct = 100 * iou
  )

  list(summary = summary, by_class = by_class, confusion = cm)
}

national_supported <- sort(unique(values(national_full, mat = FALSE)))
national_supported <- national_supported[!is.na(national_supported)]
local_supported <- sort(unique(values(local_rf, mat = FALSE)))
local_supported <- local_supported[!is.na(local_supported)]

national_result <- agreement_metrics(
  ssp1_local, national_local,
  "Original national RF output", national_supported
)
local_result <- agreement_metrics(
  ssp1_local, local_rf,
  "Local terrain RF output", local_supported
)

summary_results <- rbind(
  national_result$summary,
  local_result$summary
)
class_results <- rbind(
  national_result$by_class,
  local_result$by_class
)

write.csv(
  summary_results,
  file.path(output_dir, "ssp1_2025_rf_agreement_summary.csv"),
  row.names = FALSE
)
write.csv(
  class_results,
  file.path(output_dir, "ssp1_2025_rf_agreement_by_class.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame.matrix(national_result$confusion),
  file.path(output_dir, "ssp1_2025_national_rf_confusion.csv")
)
write.csv(
  as.data.frame.matrix(local_result$confusion),
  file.path(output_dir, "ssp1_2025_local_rf_confusion.csv")
)

print(summary_results)
print(class_results)

