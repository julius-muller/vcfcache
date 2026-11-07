#!/usr/bin/env Rscript
#### PATHS AND DATA ####
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, args[[2]])
if (dir.exists(out)) stop("Choose a new output folder; previous figures are preserved")
.libPaths(c(file.path(path.expand("~"), ".local/lib/R/4.5"), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(tibble)})
cells <- as_tibble(read.delim(file.path(root, "conditions.tsv"), na.strings = c("", "NA")))
stopifnot(nrow(cells) > 0)
labels <- c(`direct-prepared125`="Direct\nprepared", `adaptive-control121`="Adaptive\ncontrol",
            `adaptive-external125`="Adaptive\nbcftools 1.25", `adaptive-linked125`="Adaptive\nboth 1.25",
            `adaptive-prepared125`="Adaptive\nprepared", `resident-control121`="Resident\ncontrol",
            `resident-external125`="Resident\nbcftools 1.25", `resident-linked125`="Resident\nboth 1.25",
            `resident-prepared125`="Resident\nprepared")
cells$condition <- factor(labels[cells$mode], levels=unname(labels))
cells$engine <- ifelse(grepl("^resident",cells$mode),"Resident Rust",
                      ifelse(grepl("^adaptive",cells$mode),"Rust-adaptive","Direct"))
dir.create(out, recursive=TRUE)
saveRDS(list(conditions=cells),file.path(out,"pilot_data.rds"))

#### ELAPSED TIMES: ALL MEASURED CONDITIONS ####
p <- ggplot(cells,aes(condition,seconds,fill=engine)) + geom_col(width=.7) +
  facet_wrap(~case,scales="free",ncol=2) +
  scale_fill_manual(values=c("Direct"="#777777","Resident Rust"="#9D6339","Rust-adaptive"="#27658B")) +
  labs(x=NULL,y="Elapsed time (s)",fill=NULL) + theme_classic(base_size=10) +
  theme(legend.position="bottom",strip.background=element_blank(),axis.text.x=element_text(size=7))
save_plot <- function(plot,name,width,height) {
  path <- file.path(out,name)
  cairo_pdf(paste0(path,".pdf"),width=width,height=height); print(plot); dev.off()
  svg(paste0(path,".svg"),width=width,height=height); print(plot); dev.off()
  png(paste0(path,".png"),width=width,height=height,units="in",res=600,type="cairo",bg="white")
  print(plot); dev.off()
}
save_plot(p,"Figure1_rust_library_input_comparison",15,3.6*ceiling(length(unique(cells$case))/2))
dir.create(file.path(out,"panels"))
for (case in unique(cells$case)) {
  save_plot(p + cells[cells$case==case,],file.path("panels",case),10,4)
}
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"R_sessionInfo.txt"))
