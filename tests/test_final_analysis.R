# Standalone R helper/inference checks (R dependencies are not CI requirements).
.libPaths(c(file.path(path.expand("~"),".local/lib/R/4.5"),.libPaths()))
source("benchmarks/performance_sprint/final_analysis_helpers.R")
root <- "benchmarks/performance_sprint/final-analysis-20261029"
x <- enrich(read_table(file.path(root,"tables/cells.tsv")))
stopifnot(nrow(x)==600L,sum(x$validated)==572L)
# A singleton stratum must resample its index, not sample integers 1:index.
set.seed(1)
ci <- bootstrap_interval(c(1,1,100),c("A","A","B"),replicates=100)
stopifnot(identical(ci,c(1,1)))
stopifnot(all(is.na(bootstrap_interval(c(1,2),c("A","A")))))
p <- phase_data(subset(x,validated))
fractions <- aggregate(fraction~condition_id,p,sum)$fraction
stopifnot(max(abs(fractions-1))<1e-7)
duplicates <- x[c(1,1),]
stopifnot(inherits(try(summarise_conditions(duplicates,c("recipe","assay")),silent=TRUE),"try-error"))
cat("Final R analysis checks passed\n")
