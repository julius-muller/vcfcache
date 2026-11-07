#!/usr/bin/env Rscript
#### PATHS AND DATA ####
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, args[[2]])
if (dir.exists(out)) stop("Choose a new output folder; existing figures are preserved")
.libPaths(c(file.path(path.expand("~"), ".local/lib/R/4.5"), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(tibble)})
cells <- as_tibble(read.delim(file.path(root, "conditions.tsv"), na.strings = c("", "NA")))
data <- cells[cells$case != "reader", ]
stopifnot(nrow(data) > 0)
labels <- c(`direct-125` = "Direct\n1.25", `baseline-121` = "Adaptive\n1.21",
            `baseline-125` = "Adaptive\n1.25", `prepared-auto125` = "Prepared\nauto",
            `prepared-contigs2` = "Contigs\n2 shards", `prepared-windows2` = "Windows\n2 shards",
            `prepared-windows4` = "Windows\n4 shards")
data$condition <- factor(labels[data$mode], levels = unname(labels))
data$version <- ifelse(grepl("^baseline", data$mode), "Unmodified adaptive", "Prepared-input pilot")
dir.create(out, recursive = TRUE)
saveRDS(list(conditions = cells, annotation_conditions = data), file.path(out, "pilot_data.rds"))

#### ELAPSED TIME: INDIVIDUAL CONDITIONS ####
p <- ggplot(data, aes(condition, seconds, fill = version)) +
  geom_col(width = .7) + facet_wrap(~ case, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Unmodified adaptive" = "#777777", "Prepared-input pilot" = "#27658B")) +
  labs(x = NULL, y = "Elapsed time (s)", fill = NULL) + theme_classic(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        axis.text.x = element_text(size = 8))
save_plot <- function(plot, name, width, height) {
  path <- file.path(out, name)
  cairo_pdf(paste0(path, ".pdf"), width = width, height = height); print(plot); dev.off()
  svg(paste0(path, ".svg"), width = width, height = height); print(plot); dev.off()
  png(paste0(path, ".png"), width = width, height = height, units = "in", res = 600,
      type = "cairo", bg = "white"); print(plot); dev.off()
}
save_plot(p, "Figure1_elapsed_time", 12, 3 * ceiling(length(unique(data$case)) / 2))
dir.create(file.path(out, "panels"))
for (case in unique(data$case)) {
  save_plot(p + data[data$case == case, ], file.path("panels", case), 7, 4)
}
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "R_sessionInfo.txt"))
