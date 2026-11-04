#!/usr/bin/env Rscript
#### PATHS AND DATA ####
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, args[[2]])
if (dir.exists(out)) stop("Choose a new figure directory; existing plots are preserved")
.libPaths(c(file.path(path.expand("~"), ".local/lib/R/4.5"), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(tibble)})
read_data <- function(name) as_tibble(read.delim(file.path(root, name), na.strings = c("", "NA"), check.names = FALSE))
conditions <- read_data("conditions.tsv")
pairs <- read_data("paired_ratios.tsv")
summary <- read_data("summary_by_node.tsv")
attempts <- read_data("attempts.tsv")
stopifnot(nrow(conditions) > 0)
dir.create(out, recursive = TRUE)
saveRDS(list(conditions = conditions, pairs = pairs, summary_by_node = summary,
             attempts = attempts), file.path(out, "pilot_data.rds"))

#### FIXED-NODE PAIRED ELAPSED TIMES ####
labels <- c(direct = "Direct", `control-auto` = "Original\ndefault", control4 = "Original\n4 shards",
            fused4 = "Fused\n4 shards", fused2 = "Fused\n2 shards")
conditions$condition <- factor(labels[conditions$mode], levels = unname(labels))
conditions$panel <- paste(conditions$node, conditions$sample, sep = " / ")
conditions$replicate <- factor(conditions[["repeat"]])
p <- ggplot(conditions, aes(condition, wall_seconds, group = replicate,
                            colour = replicate, shape = replicate)) +
  geom_line(linewidth = .35, alpha = .6) + geom_point(size = 1.7) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("#306B8D", "#A36037")) +
  labs(x = NULL, y = "Elapsed time (s)", colour = "Technical repeat", shape = "Technical repeat") +
  theme_classic(base_size = 10) + theme(strip.background = element_blank(),
    strip.text = element_text(size = 9), legend.position = "bottom")
save_plot <- function(plot, name, width, height) {
  prefix <- file.path(out, name)
  cairo_pdf(paste0(prefix, ".pdf"), width = width, height = height); print(plot); dev.off()
  svg(paste0(prefix, ".svg"), width = width, height = height); print(plot); dev.off()
  png(paste0(prefix, ".png"), width = width, height = height, units = "in", res = 600,
      type = "cairo", bg = "white"); print(plot); dev.off()
}
save_plot(p, "Figure1_same_node_pass_fusion", 12, 3 * ceiling(length(unique(conditions$panel)) / 3))
dir.create(file.path(out, "panels"))
for (panel in unique(conditions$panel)) {
  one <- p %+% conditions[conditions$panel == panel, ]
  save_plot(one, file.path("panels", gsub("[^A-Za-z0-9_-]", "_", panel)), 6, 4)
}
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
