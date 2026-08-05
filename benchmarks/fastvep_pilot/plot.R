#!/usr/bin/env Rscript

local_library <- Sys.getenv(
  "VCFCACHE_R_LIBRARY",
  unset = file.path(
    path.expand("~"), ".local", "lib", "R",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
)
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot.R SUMMARY_TSV OUTPUT_DIR")
}

input <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "panels"), recursive = TRUE, showWarnings = FALSE)

data <- read.delim(input, stringsAsFactors = FALSE, check.names = FALSE)
data$hit_label <- ifelse(
  data$condition == "direct",
  "Direct",
  paste0(round(100 * data$observed_hit_rate), "% hits")
)
data$hit_label <- factor(
  data$hit_label,
  levels = c("Direct", "80% hits", "90% hits", "100% hits")
)

runtime <- ggplot(data, aes(hit_label, wall_seconds, fill = hit_label)) +
  geom_col(width = 0.68, colour = "#263238", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.1f s", wall_seconds)),
    vjust = -0.4,
    size = 3.7
  ) +
  scale_fill_manual(values = c("#78909C", "#80CBC4", "#26A69A", "#00796B")) +
  labs(x = NULL, y = "End-to-end wall time (seconds)", title = "A  Runtime") +
  coord_cartesian(ylim = c(0, max(data$wall_seconds) * 1.15), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

speed <- subset(data, condition != "direct")
targets <- data.frame(
  observed_hit_rate = c(0.8, 0.9),
  speedup = c(2, 3),
  series = "Product target"
)
speed$series <- "Measured (exact output)"
speedup <- ggplot(speed, aes(observed_hit_rate, speedup)) +
  geom_hline(yintercept = 1, colour = "#78909C", linewidth = 0.5) +
  geom_line(colour = "#00796B", linewidth = 1) +
  geom_point(aes(shape = series, colour = series), size = 3.4) +
  geom_point(
    data = targets,
    aes(shape = series, colour = series),
    size = 3.4,
    stroke = 1.2
  ) +
  geom_text(
    aes(label = sprintf("%.2fx", speedup)),
    nudge_y = -0.12,
    size = 3.7
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.75, 1.01)) +
  scale_shape_manual(values = c("Measured (exact output)" = 16, "Product target" = 4)) +
  scale_colour_manual(values = c("Measured (exact output)" = "#00796B", "Product target" = "#C62828")) +
  labs(
    x = "Observed cache-hit rate",
    y = "Speedup over direct fastVEP",
    shape = NULL,
    colour = NULL,
    title = "B  Incremental benefit",
    subtitle = "Targets: at least 2x at 80% and 3x at 90%"
  ) +
  theme_minimal(base_size = 12)

save_plot <- function(plot, stem, width, height) {
  prefix <- file.path(output_dir, "panels", stem)
  grDevices::svg(paste0(prefix, ".svg"), width = width, height = height, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width, height = height, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::png(
    paste0(prefix, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 300,
    type = "cairo",
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

save_plot(runtime, "fastvep_pilot_runtime", 6.2, 4.3)
save_plot(speedup, "fastvep_pilot_speedup", 6.2, 4.3)

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- runtime + speedup + patchwork::plot_layout(widths = c(1, 1))
  grDevices::svg(file.path(output_dir, "fastvep_pilot.svg"), width = 12.4, height = 4.3, bg = "white")
  print(combined)
  grDevices::dev.off()
  grDevices::cairo_pdf(file.path(output_dir, "fastvep_pilot.pdf"), width = 12.4, height = 4.3, bg = "white")
  print(combined)
  grDevices::dev.off()
  grDevices::png(
    file.path(output_dir, "fastvep_pilot.png"),
    width = 12.4,
    height = 4.3,
    units = "in",
    res = 300,
    type = "cairo",
    bg = "white"
  )
  print(combined)
  grDevices::dev.off()
}
