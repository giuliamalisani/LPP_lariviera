library(terra)
library(jsonlite)

data_dir <- "."
figure_dir <- "figures"
dir.create(figure_dir, showWarnings = FALSE)

national_cover <- rast(file.path(
  data_dir, "landcover_RF_SSPCH_2025.tif"
))
fine_cover <- rast(file.path(
  data_dir, "landcover_RF_SSPCH_Bregaglia_2025.tif"
))
coarse_cover <- rast(file.path(
  data_dir, "landcover_RF_Bregaglia_2025_posthoc_6class.tif"
))
boundary <- vect(file.path(
  data_dir, "communes", "square_20x20km_Bregaglia_commune.gpkg"
))
boundary <- project(boundary, crs(fine_cover))
national_cover <- crop(national_cover, boundary, mask = TRUE)

# Load the SSP-CH class names and colours and correct the reversed water labels so that
# it matches the exported raster values (20 = Lake; 21 = River).
ssp_palette <- fromJSON(file.path(
  data_dir, "NCCS-SSP-results", "LULC_pal.json"
))
ssp_palette$class_name[ssp_palette$class_value == 20] <- "Lake"
ssp_palette$class_name[ssp_palette$class_value == 21] <- "River"

fine_palette <- ssp_palette[, c("class_value", "class_name", "colour")]
levels(fine_cover) <- fine_palette[, c("class_value", "class_name")]
levels(national_cover) <- fine_palette[, c("class_value", "class_name")]

coarse_names <- levels(coarse_cover)[[1]]$class
coarse_colours <- c(
  "#d73027", "#bdbdbd", "#1b7837",
  "#a6d96a", "#80b1d3", "#2166ac"
)

draw_map <- function(x, colours, title) {
  par(mar = c(0.3, 0.3, 2.1, 0.3))
  plot(
    x,
    col = colours,
    type = "classes",
    axes = FALSE,
    main = title,
    cex.main = 1.15,
    all_levels = TRUE,
    legend = FALSE
  )
  plot(boundary, add = TRUE, border = "#252525", lwd = 1.1)
}

draw_legend <- function(labels, colours, columns = 1L) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = labels,
    fill = colours,
    border = "#252525",
    ncol = columns,
    bty = "n",
    cex = 0.82,
    x.intersp = 0.7,
    y.intersp = 0.9
  )
}

# Standalone detailed map, useful if the two panels are separated later.
png(
  file.path(figure_dir, "chapter2_rf_bregaglia_2025_fine.png"),
  width = 2200, height = 1650, res = 240, bg = "white"
)
layout(matrix(c(1, 2), nrow = 1), widths = c(4.4, 1.6))
draw_map(
  fine_cover, fine_palette$colour,
  "Detailed classification (9 predicted classes)"
)
draw_legend(fine_palette$class_name, fine_palette$colour)
dev.off()

# Report figure: both thematic resolutions from the same fine prediction.
png(
  file.path(figure_dir, "chapter2_rf_bregaglia_2025_both_levels.png"),
  width = 3200, height = 1900, res = 240, bg = "white"
)
layout(
  matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
  widths = c(1, 1), heights = c(4.4, 1.25)
)
draw_map(
  fine_cover, fine_palette$colour,
  "Detailed classification (9 predicted classes)"
)
draw_map(coarse_cover, coarse_colours, "Aggregated classification (6 classes)")
draw_legend(fine_palette$class_name, fine_palette$colour, columns = 2L)
draw_legend(coarse_names, coarse_colours, columns = 2L)
dev.off()

# Full workflow figure: original practical baseline followed by the two added
# local products. The first two panels share the complete SSP-CH class scheme.
png(
  file.path(
    figure_dir,
    "chapter2_rf_bregaglia_2025_baseline_and_additions.png"
  ),
  width = 4200, height = 1900, res = 240, bg = "white"
)
layout(
  matrix(c(1, 2, 3, 4, 4, 5), nrow = 2, byrow = TRUE),
  widths = c(1, 1, 1), heights = c(4.4, 1.25)
)
draw_map(
  national_cover, fine_palette$colour,
  "Original practical: national RF (11 classes)"
)
draw_map(
  fine_cover, fine_palette$colour,
  "Addition: local terrain RF (9 classes)"
)
draw_map(
  coarse_cover, coarse_colours,
  "Addition: post-hoc aggregation (6 classes)"
)
draw_legend(fine_palette$class_name, fine_palette$colour, columns = 3L)
draw_legend(coarse_names, coarse_colours, columns = 2L)
dev.off()

