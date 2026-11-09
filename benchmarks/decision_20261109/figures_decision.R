#!/usr/bin/env Rscript
# All plotting is R/ggplot2. Inputs and PNGs live together in this flat folder.
# No package installation, benchmark submission, or changes to earlier figures.

#### PATHS AND INPUT DATA ####
args <- commandArgs(TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
source_file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
root <- if (length(args)) args[[1]] else if(length(file_arg)) {
  dirname(sub("^--file=", "", file_arg[[1]]))
} else if(!is.null(source_file)) {
  dirname(source_file)
} else {
  # For selected-line execution in RStudio, edit this single path if needed.
  "/home/j380r/projects/vcfcache/benchmarks/decision_20261109"
}
root <- normalizePath(root, mustWork=TRUE)
.libPaths(c(Sys.getenv("VCFCACHE_R_LIBRARY", unset="/home/j380r/.local/lib/R/4.5"), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(tibble); library(scales)})
read_input <- function(name) as_tibble(read.delim(file.path(root,name),na.strings=c("","NA"),check.names=FALSE))
measurements <- read_input("measurements.tsv")
summaries <- read_input("summary.tsv")
configurations <- read_input("configurations.tsv")
sources <- read_input("sources.tsv")

# Add labels only; every numerical value was derived from the documented raw
# exports. Never combine complete-cohort, archival, and pilot evidence strata.
s <- merge(summaries,configurations,by="config_id",all.x=TRUE,sort=FALSE)
m <- merge(measurements,configurations,by="config_id",all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(s$label), !anyNA(m$label), all(m$cached_s>0))
assays <- c("Panel","WES","WGS")
recipe_order <- c("vep-rs core","fastVEP basic","fastVEP clinical","SnpEff core",
                  "vep-rs clinical","SnpEff clinical","VEP everything","VEP everything (4 forks)")
wrap <- function(x,width=48) vapply(x,function(z) paste(strwrap(z,width=width),collapse="\n"),character(1))
labels <- setNames(configurations$label,configurations$config_id)
palette <- setNames(hcl.colors(nrow(configurations),"Dark 3"),configurations$config_id)
# Keep the two decision-relevant cohort engines visually distinct everywhere.
palette[configurations$config_id[configurations$plot_group=="cohort-primary" & configurations$engine=="Python-adaptive"]] <- "#226F9D"
palette[configurations$config_id[configurations$plot_group=="cohort-primary" & configurations$engine=="Rust-adaptive"]] <- "#C56635"
matrix_ids <- configurations$config_id[configurations$plot_group=="matrix"]
palette[matrix_ids] <- c("#555555","#999999","#9467BD","#D77720","#18785E","#AB963B")[seq_along(matrix_ids)]
for (i in which(configurations$plot_group=="delay")) {
  is_python <- configurations$engine[[i]]=="Python-adaptive"
  af10 <- grepl("10%",configurations$label[[i]],fixed=TRUE)
  palette[configurations$config_id[[i]]] <- if(is_python) {
    if(af10) "#79A8C9" else "#226F9D"
  } else {
    if(af10) "#DBA263" else "#A95224"
  }
}
plot_data <- function(data) {
  data$mark <- ifelse(grepl("timing-only|failed|unknown",data$validation),
                      "Annotation differences / timing only","Accepted under recorded comparator")
  data$fallback <- grepl("fallback",data$engine)
  data
}
s <- plot_data(s); m <- plot_data(m)
theme_set(theme_classic(base_size=11) +
            theme(legend.position="bottom",legend.box="vertical",strip.background=element_blank(),
                  strip.text=element_text(face="bold"),axis.title=element_text(face="plain"),
                  panel.spacing=grid::unit(1.4,"lines")))
colour_key <- function(ids,ncol=2) {
  scale_colour_manual(values=palette[ids],breaks=ids,labels=wrap(labels[ids]),
                      name="Branch / benchmark configuration",
                      guide=guide_legend(ncol=ncol,byrow=TRUE,override.aes=list(alpha=1,shape=16)))
}
shape_key <- function() scale_shape_manual(values=c("Accepted under recorded comparator"=16,
                                                   "Annotation differences / timing only"=1),name=NULL)
fold_axis <- function() scale_y_continuous(trans="log2",breaks=c(.125,.25,.5,1,2,4,8,16,32),
                                          labels=function(x) paste0(format(x,trim=TRUE),"x"))
exports <- character()
save_png <- function(plot,name,width=13,height=7,dpi=240) {
  destination <- file.path(root,paste0(name,".png"))
  # This is a new decision folder. Existing publication figures are never used
  # as output destinations; rerunning replaces only these derived decision PNGs.
  ggsave(destination,plot,width=width,height=height,units="in",dpi=dpi,bg="white",limitsize=FALSE)
  exports <<- c(exports,basename(destination))
}

#### TYPE 1: COMPLETE-COHORT SPEEDUP BY ANNOTATION RECIPE ####
# One line per engine/configuration, 12 distinct GRCh38 inputs per recipe/assay.
# Points underneath expose the input-level spread, including slow samples.
# These are observed distributions, NOT confidence intervals or a randomized
# Python-vs-Rust trial; the direct denominators were retained across campaigns.
for (current_assay in assays) {
  d <- s[s$plot_group=="cohort-primary" & s$assay==current_assay,]
  raw <- m[m$plot_group=="cohort-primary" & m$assay==current_assay,]
  d$recipe <- factor(d$recipe,levels=recipe_order)
  raw$recipe <- factor(raw$recipe,levels=recipe_order)
  ids <- configurations$config_id[configurations$plot_group=="cohort-primary"]
  p <- ggplot(d,aes(recipe,speedup_median,colour=config_id,group=config_id)) +
    geom_hline(yintercept=1,colour="grey55",linetype="dashed") +
    geom_point(data=raw,aes(recipe,speedup,shape=mark),alpha=.27,size=1.5,
               position=position_jitter(width=.045,height=0,seed=14)) +
    geom_line(linewidth=.85) + geom_point(aes(shape=mark),size=3,stroke=1) +
    colour_key(ids) + shape_key() + fold_axis() +
    scale_x_discrete(labels=function(x) wrap(x,16)) +
    labs(x="Annotation recipe (different information content)",y="Speedup vs own direct baseline") +
    theme(legend.text=element_text(size=9))
  save_png(p,paste0("01_cohort_recipe_speedup_",current_assay),height=7)
}

#### TYPE 2: HISTORICAL BRANCH MATRIX, MATCHED 16-CPU CONTROLS ####
# This is engineering evidence from ONE SGDP input, not the cohort above.
# A native-request fallback is its OWN labelled series, never resident WGS.
for (current_assay in assays) {
  d <- s[s$plot_group=="matrix" & s$assay==current_assay & !is.na(s$speedup_median),]
  d$recipe <- factor(d$recipe,levels=recipe_order)
  ids <- configurations$config_id[configurations$config_id %in% d$config_id]
  p <- ggplot(d,aes(recipe,speedup_median,colour=config_id,group=config_id)) +
    geom_hline(yintercept=1,colour="grey55",linetype="dashed") +
    geom_line(aes(linetype=fallback),linewidth=.8) +
    geom_point(aes(shape=mark),size=3,stroke=1) +
    scale_linetype_manual(values=c(`FALSE`="solid",`TRUE`="dotted"),guide="none") +
    colour_key(ids) + shape_key() + fold_axis() +
    scale_x_discrete(labels=function(x) wrap(x,16)) +
    labs(x="Annotation recipe — archived one-input, 16-CPU comparison",y="Speedup vs own direct baseline") +
    theme(legend.text=element_text(size=8.5))
  save_png(p,paste0("02_historical_branch_speedup_",current_assay),height=8.6)
}

#### TYPE 3: EFFECTIVE ANNOTATION COST VERSUS CACHING BENEFIT ####
# Observed direct seconds/N includes startup, I/O and output work; it is NOT
# pure annotation cost, a biological-complexity score, or a causal model fit.
# Lines connect recipe medians for navigation only. Strata remain faceted.
for (current_assay in assays) {
  d <- s[s$plot_group %in% c("cohort-primary","matrix","delay") & s$assay==current_assay &
           !is.na(s$speedup_median) & !is.na(s$direct_us_per_variant_median),]
  d$evidence <- ifelse(d$plot_group=="cohort-primary","12-input cohort; 8 CPUs",
                       ifelse(d$plot_group=="matrix","One-input branch matrix; 16 CPUs",
                              "One WGS input; controlled delay; 8 CPUs"))
  d <- d[order(d$config_id,d$direct_us_per_variant_median),]
  ids <- configurations$config_id[configurations$config_id %in% d$config_id]
  p <- ggplot(d,aes(direct_us_per_variant_median,speedup_median,colour=config_id,group=config_id)) +
    geom_hline(yintercept=1,colour="grey55",linetype="dashed") +
    geom_line(aes(linetype=fallback),linewidth=.75) +
    geom_point(aes(shape=mark),size=2.7,stroke=.9) +
    facet_wrap(~evidence,scales="free_x",ncol=2) +
    scale_linetype_manual(values=c(`FALSE`="solid",`TRUE`="dotted"),guide="none") +
    scale_x_log10(labels=label_number(accuracy=.1,big.mark=",")) + fold_axis() +
    colour_key(ids) + shape_key() +
    labs(x="Observed direct runtime / variant (microseconds; includes startup and I/O)",
         y="Speedup vs own direct baseline") + theme(legend.text=element_text(size=8))
  save_png(p,paste0("03_observed_cost_speedup_",current_assay),width=15,height=if(current_assay=="WGS") 12 else 9)
}

#### OVERALL: ALL CONFIGURATIONS, INCLUDING BLANK/FAILED EVIDENCE ####
# This deliberately large overview is zoomable; the workbook is the readable
# detailed version. Asterisks mark annotation differences, NOT significance.
# Delay-only rows are represented in Type 3 and the workbook, not blanked here.
heat_configs <- configurations[configurations$plot_group!="delay",]
grid <- expand.grid(config_id=heat_configs$config_id,assay=assays,recipe=recipe_order,
                    stringsAsFactors=FALSE)
h <- merge(grid,s[,c("config_id","assay","recipe","speedup_median","validation","runtime_median_s")],
           by=c("config_id","assay","recipe"),all.x=TRUE,sort=FALSE)
h$y <- factor(h$config_id,levels=rev(heat_configs$config_id))
h$recipe <- factor(h$recipe,levels=recipe_order)
h$assay <- factor(h$assay,levels=assays)
h$txt <- ifelse(is.na(h$speedup_median),ifelse(is.na(h$runtime_median_s),"","s"),
                paste0(formatC(h$speedup_median,format="f",digits=2),
                       ifelse(grepl("timing-only|failed|unknown",h$validation),"*","")))
p <- ggplot(h,aes(recipe,y,fill=log2(speedup_median))) +
  geom_tile(colour="white",linewidth=.2) + geom_text(aes(label=txt),size=2.4,na.rm=TRUE) +
  facet_grid(~assay) +
  scale_fill_gradient2(low="#BC665B",mid="#F6F6F2",high="#267872",midpoint=0,
                       limits=c(-2,5),oob=squish,na.value="#E6E6E6",
                       breaks=c(-2,-1,0,1,2,3,4,5),labels=c("0.25x","0.5x","1x","2x","4x","8x","16x","32x"),
                       name="Speedup vs direct") +
  scale_y_discrete(labels=function(x) wrap(labels[x],68)) +
  scale_x_discrete(labels=function(x) wrap(x,14)) +
  labs(x=NULL,y=NULL) + theme(axis.text.y=element_text(size=6.5),axis.text.x=element_text(size=8),
                              axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")
save_png(p,"04_all_configurations_overview",width=25,height=max(18,nrow(heat_configs)*.34),dpi=200)

#### OVERALL: HISTORICAL PYTHON VERSUS RUST ON IDENTICAL CONDITIONS ####
# Condition IDs, denominators, inputs and cache payloads match; execution dates
# and workers do NOT always match. This is a historical comparison, not causal.
cohort <- subset(m,plot_group=="cohort-primary")
a <- subset(cohort,engine=="Python-adaptive")
b <- subset(cohort,engine=="Rust-adaptive")
paired <- merge(a,b,by=c("condition_id","assay","recipe","sample"),suffixes=c("_python","_rust"))
stopifnot(nrow(paired)==180,all(abs(paired$direct_s_python-paired$direct_s_rust)<1e-6))
paired$gain <- paired$cached_s_python/paired$cached_s_rust
paired$recipe <- factor(paired$recipe,levels=recipe_order)
paired$mark <- ifelse(grepl("timing-only|failed|unknown",paired$validation_python),
                       "Annotation differences / timing only","Accepted under recorded comparator")
p <- ggplot(paired,aes(recipe,gain)) + geom_hline(yintercept=1,colour="grey55",linetype="dashed") +
  geom_boxplot(width=.5,outlier.shape=NA,fill="#DCE8EC",linewidth=.5) +
  geom_point(aes(shape=mark),position=position_jitter(width=.09,height=0,seed=23),size=2,alpha=.65) +
  facet_wrap(~assay,ncol=3) + scale_y_log10(breaks=c(.5,.75,1,1.25,1.5,2,3),labels=label_number(accuracy=.01)) +
  shape_key() + scale_x_discrete(labels=function(x) wrap(x,12)) +
  labs(x="Annotation recipe",y="Python cached time / Rust cached time\nAbove 1: Rust faster") +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=8))
save_png(p,"05_python_vs_rust_same_conditions",width=16,height=6.5)

#### OVERALL: DOES THE RECENT RUST 1.25 CHANGE ACTUALLY HELP? ####
# Compare each new configuration to its OWN same-block pre-1.25 engine, never
# compare resident directly to adaptive and call it a library gain.
recent <- subset(m,plot_group=="rust125")
recent$mode <- sub(".*; ([^;]+)\\]$","\\1",recent$label)
recent$engine_class <- ifelse(grepl("^resident",recent$mode),"Resident Rust","Rust-adaptive")
control <- subset(recent,grepl("control121$",mode))
d <- merge(subset(recent,!grepl("control121$",mode)),
           control[,c("sample","assay","recipe","engine_class","cached_s")],
           by=c("sample","assay","recipe","engine_class"),suffixes=c("","_before"))
d$gain <- d$cached_s_before/d$cached_s
d$change <- ifelse(grepl("external125$",d$mode),"External bcftools 1.25",
                  ifelse(grepl("linked125$",d$mode),"External + linked HTSlib 1.25",
                         "Both 1.25 + input preserved"))
d$condition <- factor(d$config_id,levels=rev(configurations$config_id[configurations$plot_group=="rust125"]))
d$case <- paste(d$recipe,d$assay,d$sample,sep=" / ")
p <- ggplot(d,aes(gain,condition,colour=engine_class)) +
  geom_vline(xintercept=1,colour="grey50",linetype="dashed") + geom_point(size=2.8) +
  facet_wrap(~case,ncol=2,scales="free_y",labeller=label_wrap_gen(48)) +
  scale_y_discrete(labels=function(x) wrap(labels[x],58)) +
  scale_colour_manual(values=c("Resident Rust"="#18785E","Rust-adaptive"="#D77720"),name=NULL) +
  scale_x_continuous(breaks=c(.9,.95,1,1.05,1.1,1.15),labels=label_number(accuracy=.01)) +
  labs(x="Original cached runtime / updated cached runtime\nAbove 1: update faster",y=NULL) +
  theme(axis.text.y=element_text(size=7),strip.text=element_text(size=9))
save_png(p,"06_rust_125_incremental_gain",width=18,height=13)

#### TYPE 4: ABSOLUTE CACHED VERSUS DIRECT LATENCY ####
# Below the parity diagonal means faster with caching. Keeping absolute seconds
# visible prevents a fold-only graphic from concealing short-pipeline latency.
for (current_assay in assays) {
  latency <- s[s$plot_group %in% c("cohort-primary","matrix") &
                 s$assay==current_assay & !is.na(s$direct_median_s),]
  latency$evidence <- ifelse(latency$plot_group=="cohort-primary",
                             "12-input cohort; 8 CPUs","One-input matrix; 16 CPUs")
  latency <- latency[order(latency$config_id,latency$direct_median_s),]
  ids <- configurations$config_id[configurations$config_id %in% latency$config_id]
  p <- ggplot(latency,aes(direct_median_s,runtime_median_s,colour=config_id,group=config_id)) +
    geom_abline(data=data.frame(slope=1,intercept=0,reference="Equal runtime"),
                aes(slope=slope,intercept=intercept,linetype=reference),colour="grey55",show.legend=TRUE) +
    geom_line(linewidth=.75) + geom_point(aes(shape=mark),size=2.8,stroke=.9) +
    facet_wrap(~evidence,ncol=2) + scale_x_log10(labels=label_number(big.mark=",")) +
    scale_y_log10(labels=label_number(big.mark=",")) +
    scale_linetype_manual(values=c("Equal runtime"="dashed"),name=NULL) +
    colour_key(ids) + shape_key() +
    labs(x="Direct runtime (s; recipe median)",y="Cached runtime (s; recipe median)") +
    theme(legend.text=element_text(size=8))
  save_png(p,paste0("07_latency_vs_direct_",current_assay),width=15,height=8.8)
}

#### RSTUDIO OBJECTS AND REPRODUCIBILITY ####
# A named tibble list makes exploration possible without rerunning collection.
saveRDS(list(measurements=as_tibble(m),summaries=as_tibble(s),configurations=configurations,
             historical_python_rust_pairs=as_tibble(paired),recent_rust_upgrade_pairs=as_tibble(d),
             sources=sources),file.path(root,"decision_data.rds"))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(root,"R_sessionInfo.txt"))
writeLines(exports,file.path(root,"PNG_FILES.txt"))
cat(length(exports),"PNG figures written to",root,"\n")
