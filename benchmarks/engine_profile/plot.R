#!/usr/bin/env Rscript
# Engineering diagnostics only. Panels are not new cohort efficacy claims.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- args[[1]]
cells <- read.delim(file.path(root, "cells.tsv"), check.names = FALSE)
phases <- read.delim(file.path(root, "phases.tsv"), check.names = FALSE)
if (!nrow(phases)) quit(status = 0)
keep <- c("native.lookup", "native.construct_encode_output", "native.finalize_output",
          "native.prepare_decode_split_stream", "pipeline.preprocess", "pipeline.route",
          "pipeline.annotate-misses", "pipeline.index")
phases <- phases[phases$phase %in% keep, ]
phases <- merge(phases, cells[, c("cell", "purpose", "state")], by = "cell")
phases <- phases[phases$purpose == "primary" & grepl("^complete", phases$state), ]
if (nrow(phases)) {
  p <- ggplot(phases, aes(union_seconds, phase, color = effective_engine)) +
    geom_point(size = 2, position = position_dodge(width = 0.45)) +
    facet_grid(assay ~ recipe, scales = "free_x") +
    scale_color_manual(values = c(native = "#0072B2", adaptive = "#D55E00", fused = "#D55E00", legacy = "#666666")) +
    labs(x = "Inclusive stage wall time (s; stages may overlap)", y = NULL, color = "Executed backend") +
    theme_classic(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(root, "stage_profile.pdf"), p, width = 11, height = 7)
  ggsave(file.path(root, "stage_profile.png"), p, width = 11, height = 7, dpi = 400)
}
scaling <- cells[cells$purpose %in% c("primary", "cpu-scaling", "worker-scaling") &
                   !is.na(cells$wall_seconds) & cells$role == "measure", ]
if (nrow(scaling)) {
  p <- ggplot(scaling, aes(cpus, wall_seconds, color = effective_engine, shape = purpose)) +
    geom_point(size = 2) + facet_grid(assay ~ recipe, scales = "free_y") +
    labs(x = "Enforced CPU allocation", y = "End-to-end wall time (s)", color = "Executed backend", shape = "Condition") +
    theme_classic(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(root, "cpu_scaling.pdf"), p, width = 11, height = 6)
  ggsave(file.path(root, "cpu_scaling.png"), p, width = 11, height = 6, dpi = 400)
}
