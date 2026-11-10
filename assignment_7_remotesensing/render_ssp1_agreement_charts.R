summary_data <- read.csv("results/ssp1_2025_rf_agreement_summary.csv")
class_data <- read.csv("results/ssp1_2025_rf_agreement_by_class.csv")
library(jsonlite)

figure_dir <- "figures"
dir.create(figure_dir, showWarnings = FALSE)

model_labels <- c("National RF", "Local terrain RF")
model_colours <- c("#7f8c8d", "#1b7837")

# Overall comparison across complementary agreement measures.
metric_matrix <- rbind(
  summary_data$overall_agreement_pct,
  summary_data$balanced_agreement_pct,
  summary_data$macro_iou_pct
)
rownames(metric_matrix) <- c(
  "Overall agreement", "Balanced agreement", "Macro-IoU"
)

png(
  file.path(figure_dir, "ssp1_2025_rf_agreement_summary.png"),
  width = 2200, height = 1450, res = 240, bg = "white"
)
par(mar = c(5.5, 5, 2.2, 1))
barplot(
  t(metric_matrix),
  beside = TRUE,
  col = model_colours,
  border = NA,
  ylim = c(0, 65),
  names.arg = rownames(metric_matrix),
  ylab = "Agreement with SSP1-2025 (%)",
  las = 1,
  cex.names = 0.9
)
abline(h = seq(0, 60, by = 10), col = "#dddddd", lty = 1)
legend(
  "topright", legend = model_labels, fill = model_colours,
  border = NA, bty = "n", cex = 0.9
)
dev.off()

# Class-wise intersection over union. Zero values are retained because they
# expose SSP1 classes that a classifier cannot produce.
national_iou <- class_data$iou_pct[
  class_data$model == "Original national RF output"
]
local_iou <- class_data$iou_pct[
  class_data$model == "Local terrain RF output"
]
class_names <- class_data$class_name[
  class_data$model == "Original national RF output"
]
ssp_colours <- fromJSON("NCCS-SSP-results/LULC_pal.json")
class_colours <- ssp_colours$colour[
  match(
    class_data$class_value[class_data$model == "Original national RF output"],
    ssp_colours$class_value
  )
]

png(
  file.path(figure_dir, "ssp1_2025_rf_class_iou.png"),
  width = 2400, height = 1900, res = 240, bg = "white"
)
par(mfrow = c(1, 2), mar = c(4.5, 9.5, 2.5, 1))
for (i in seq_along(model_labels)) {
  values <- if (i == 1) national_iou else local_iou
  barplot(
    rev(values),
    horiz = TRUE,
    col = rev(class_colours),
    border = NA,
    xlim = c(0, 65),
    names.arg = rev(class_names),
    xlab = "IoU with SSP1-2025 (%)",
    main = model_labels[i],
    las = 1,
    cex.names = 0.72
  )
  abline(v = seq(0, 60, by = 10), col = "#dddddd", lty = 1)
}
dev.off()

