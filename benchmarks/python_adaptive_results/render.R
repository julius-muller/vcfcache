#!/usr/bin/env Rscript
# Plot only measured Python-adaptive results. No old figure is overwritten.

#### PATHS AND PACKAGES ####
args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Usage: render.R <snapshot directory> [new figure directory name]")
root <- normalizePath(args[[1]], mustWork = TRUE)
.libPaths(c(file.path(path.expand("~"), ".local/lib/R/4.5"), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(tibble)})
out <- file.path(root, if(length(args)==2L) args[[2]] else "figures")
if (dir.exists(out)) stop("Choose a new output directory; old figures are immutable")
# Carry the unchanged conceptual workflow alongside the regenerated empirical
# figures. This is a copy of a vector figure, not a newly timed measurement.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
workflow <- file.path(dirname(normalizePath(script_file)),"../figures/output/final/2026-08-27T225600Z-story-campaign-final-v4/manuscript")
dir.create(file.path(out,"manuscript/panels"),recursive=TRUE)
for (ext in c("pdf","svg","png")) {
  stopifnot(file.copy(file.path(workflow,paste0("figure_workflow_and_value.",ext)),
    file.path(out,"manuscript",paste0("Figure1_workflow.",ext)),overwrite=FALSE))
  for(panel in c("A","B")) stopifnot(file.copy(
    file.path(workflow,"panels",paste0("figure_workflow_",panel,".",ext)),
    file.path(out,"manuscript/panels",paste0("figure_workflow_",panel,".",ext)),overwrite=FALSE))
}
read_tsv <- function(p) as_tibble(read.delim(p, quote="", comment.char="", na.strings=c("","NA")))
write_tsv <- function(x, name) write.table(x, file.path(root,"tables",name), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
save_plot <- function(p, name, w=10, h=5) {
  prefix <- file.path(out,name); dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  cairo_pdf(paste0(prefix,".pdf"),width=w,height=h,family="sans");print(p);dev.off()
  svg(paste0(prefix,".svg"),width=w,height=h,family="sans");print(p);dev.off()
  png(paste0(prefix,".png"),width=w,height=h,units="in",res=600,type="cairo",bg="white");print(p);dev.off()
}
theme_journal <- function(size=9) theme_classic(base_size=size) + theme(
  axis.text=element_text(colour="#24343E"), axis.line=element_line(linewidth=.3),
  strip.background=element_blank(),strip.text=element_text(face="bold"),
  plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),
  plot.tag=element_text(face="bold",size=12),legend.position="bottom",
  plot.margin=margin(8,10,8,8))
cohort_colours <- c(KPGP="#306B8D",SGDP="#A36037",PGP="#78618E")
cohort_shapes <- c(KPGP=16,SGDP=17,PGP=15)
strategy_colours <- c("AF ≥ 10%"="#47718B","AF ≥ 1%"="#2F806A","Three genomes"="#B47837")
point_style <- function() list(scale_colour_manual(values=cohort_colours),
  scale_shape_manual(values=cohort_shapes),guides(colour=guide_legend(override.aes=list(alpha=1,size=2.5))))
ratio_axis <- function() scale_x_log10(labels=scales::label_number(suffix="×"))

#### DATA LOADING AND INTEGRITY ####
cells <- read_tsv(file.path(root,"tables/cells.tsv"))
views <- read_tsv(file.path(root,"tables/plot_views.tsv"))
enrich <- function(x) {
  x$validated <- tolower(x$validated)=="true"
  x$is_alias <- tolower(x$is_alias)=="true"
  x$recipe <- ifelse(x$tool=="vep","VEP everything",
    ifelse(x$tool=="fastvep",paste("fastVEP",ifelse(x$level=="clinical","clinical","basic")),
    ifelse(x$tool=="snpeff",paste("SnpEff",x$level),paste("vep-rs",x$level))))
  x$assay <- factor(toupper(x$assay),levels=c("PANEL","WES","WGS"),labels=c("Panel","WES","WGS"))
  x$cohort <- factor(toupper(x$cohort),levels=names(cohort_colours))
  x$strategy_label <- factor(ifelse(grepl("0.01$|1%",x$strategy)&!grepl("10%",x$strategy),"AF ≥ 1%",
    ifelse(grepl("0.1$|10%",x$strategy),"AF ≥ 10%","Three genomes")),levels=names(strategy_colours))
  x$saved_hours <- x$saved_seconds/3600
  x$cached_ratio <- x$cached_seconds/x$direct_seconds
  x
}
cells <- enrich(cells); views <- enrich(views)
stopifnot(nrow(cells)==600,nrow(views)==612,sum(cells$validated)==572,
  !anyDuplicated(cells$measurement_id),sum(views$is_alias)==12,!anyNA(cells$cohort))
assays <- subset(cells,family %in% c("native_clinical_assays","vep_matched_assays","veprs_extension"))
accepted <- subset(assays,validated)
primary <- subset(accepted,assembly=="GRCh38")
wgs <- subset(views,family=="full52_wgs_strategies" & validated)
stopifnot(nrow(wgs)==312,all(table(wgs$tool,wgs$strategy)==52))
synthetic <- subset(cells,family=="synthetic_pipeline_spectrum")
stopifnot(nrow(synthetic)==12,length(unique(synthetic$sample))==1,all(synthetic$validated))

#### SAMPLE-LEVEL SUMMARIES ####
# One sample/condition, not technical repeats. Preserve cohort allocation when
# bootstrapping pooled medians; separate summaries retain assembly and recipe.
set.seed(20261103)
interval <- function(x, cohort) {
  if(length(x)<3) return(c(NA_real_,NA_real_))
  strata <- split(seq_along(x),cohort,drop=TRUE)
  quantile(replicate(4000,median(x[unlist(lapply(strata,function(i) i[sample.int(length(i),length(i),replace=TRUE)]))])),c(.025,.975),names=FALSE)
}
summarise <- function(x, groups) {
  split_index <- split(seq_len(nrow(x)),interaction(x[,groups],drop=TRUE,lex.order=TRUE))
  as_tibble(do.call(rbind,lapply(split_index,function(i) {
    z<-x[i,];ci<-interval(z$speedup,z$cohort)
    cbind(z[1,groups,drop=FALSE],data.frame(n=nrow(z),median_speedup=median(z$speedup),
      ci_low=ci[1],ci_high=ci[2],min_speedup=min(z$speedup),max_speedup=max(z$speedup),
      fraction_faster=mean(z$speedup>1),median_saved_hours=median(z$saved_hours),
      median_direct_seconds=median(z$direct_seconds),median_cached_seconds=median(z$cached_seconds),
      median_hit_rate=median(z$cache_hit_rate)))
  })))
}
assay_summary <- summarise(accepted,c("recipe","assembly","assay"))
primary_summary <- summarise(primary,c("recipe","assembly","assay"))
wgs_summary <- summarise(wgs,c("recipe","cohort","assembly","strategy_label"))
write_tsv(assay_summary,"assay_summary.tsv");write_tsv(primary_summary,"primary_grch38_summary.tsv")
write_tsv(wgs_summary,"wgs_strategy_summary.tsv")
write_tsv(summarise(accepted,c("recipe","cohort","assembly","assay")),"assay_summary_by_cohort.tsv")

#### FIGURE 2: EMPIRICAL TIME AND THROUGHPUT ####
# The lead comparison uses GRCh38 only. PGP/GRCh37 remains visible separately
# in the supplementary figure rather than treating assemblies as matched recipes.
recipes <- c("VEP everything","fastVEP clinical","SnpEff clinical","vep-rs core")
main <- subset(primary,recipe %in% recipes); main$recipe <- factor(main$recipe,levels=rev(recipes))
ms <- subset(primary_summary,recipe %in% recipes);ms$recipe<-factor(ms$recipe,levels=rev(recipes))
p2a <- ggplot(main,aes(speedup,recipe)) +
  geom_vline(xintercept=1,linetype=2,colour="#A3A9AD",linewidth=.35) +
  geom_point(aes(colour=cohort,shape=cohort),position=position_jitter(height=.09,width=0,seed=3),size=1.7,alpha=.8) +
  geom_errorbar(data=ms,aes(x=median_speedup,xmin=ci_low,xmax=ci_high),orientation="y",width=.24,linewidth=.6) +
  geom_point(data=ms,aes(x=median_speedup),shape=18,size=2.6) +
  facet_wrap(~assay,nrow=1) + ratio_axis() + point_style() +
  labs(x="Speedup over the same recipe without caching",y=NULL,colour=NULL,shape=NULL) + theme_journal()
p2b <- ggplot(subset(main,assay=="WGS"),aes(saved_hours,recipe)) +
  geom_boxplot(outlier.shape=NA,width=.45,fill="#EDF1F3",linewidth=.4) +
  geom_point(aes(colour=cohort,shape=cohort),position=position_jitter(width=0,height=.1,seed=4),size=1.7) +
  point_style() + labs(x="Time saved per WGS sample (h)",y=NULL,colour=NULL,shape=NULL) + theme_journal()
p2c <- ggplot(subset(main,assay=="WGS"),aes(samples_per_day_direct,samples_per_day_cached,colour=cohort,shape=cohort)) +
  geom_abline(slope=1,intercept=0,linetype=2,colour="#A3A9AD",linewidth=.35) +
  geom_point(size=1.9) + facet_wrap(~recipe,nrow=2,scales="free") +
  scale_x_log10() + scale_y_log10() + point_style() +
  labs(x="Direct samples per 24 h",y="Cached samples per 24 h",colour=NULL,shape=NULL) + theme_journal()
save_plot(p2a+labs(tag="A"),"manuscript/panels/Figure2A_assay_speedup",10.5,3.7)
save_plot(p2b+labs(tag="B"),"manuscript/panels/Figure2B_WGS_time_saved",6,3.7)
save_plot(p2c+labs(tag="C"),"manuscript/panels/Figure2C_WGS_throughput",6,4.2)
save_plot((p2a/(p2b|p2c))+plot_layout(heights=c(1,1.2))+plot_annotation(tag_levels="A"),"manuscript/Figure2_empirical_performance",11,8.2)

#### FIGURE 3: WGS CACHE STRATEGIES ####
hit <- subset(wgs,tool=="vep")
p3a <- ggplot(hit,aes(strategy_label,cache_hit_rate,colour=strategy_label)) +
  geom_boxplot(width=.5,outlier.shape=NA,linewidth=.4) +
  geom_point(position=position_jitter(width=.12,height=0,seed=5),size=1.3,alpha=.7) +
  facet_wrap(~cohort,nrow=1) + scale_colour_manual(values=strategy_colours) +
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1)) +
  labs(x=NULL,y="Cache hit rate")+theme_journal()+theme(legend.position="none")
p3b <- ggplot(wgs,aes(speedup,colour=strategy_label,linetype=strategy_label))+
  stat_ecdf(linewidth=.65)+facet_grid(recipe~cohort)+ratio_axis()+
  scale_y_continuous(labels=scales::label_percent(),breaks=c(0,.5,1))+
  scale_colour_manual(values=strategy_colours)+labs(x="Per-sample speedup",y="Cumulative fraction of genomes",colour=NULL,linetype=NULL)+theme_journal()
save_plot(p3a+labs(tag="A"),"manuscript/panels/Figure3A_WGS_cache_hits",10.5,3.2)
save_plot(p3b+labs(tag="B"),"manuscript/panels/Figure3B_WGS_speedup_distributions",10.5,5)
save_plot((p3a/p3b)+plot_layout(heights=c(1,1.5))+plot_annotation(tag_levels="A"),"manuscript/Figure3_WGS_cache_strategies",10.5,8)

#### FIGURE S1: CONFIGURABLE PIPELINE COMPLEXITY ####
synthetic$delay_ms <- as.numeric(sub(".*delay_([0-9]+)us.*","\\1",synthetic$level))/1000
stopifnot(!anyNA(synthetic$delay_ms))
synthetic$delay <- factor(synthetic$delay_ms,levels=sort(unique(synthetic$delay_ms)))
synthetic$scenario <- paste(as.character(synthetic$strategy_label),sprintf("(%.1f%% hits)",100*synthetic$cache_hit_rate))
runtime <- rbind(transform(synthetic,method="Direct",hours=direct_seconds/3600),transform(synthetic,method="VCFcache",hours=cached_seconds/3600))
ps1a <- ggplot(runtime,aes(delay,hours,fill=method))+geom_col(position=position_dodge(.75),width=.65)+
  facet_wrap(~scenario,nrow=1)+scale_fill_manual(values=c(Direct="#BCC5CA",VCFcache="#347F69"))+
  labs(x="Configured wait per consequence (ms)",y="Elapsed time (h)",fill=NULL)+theme_journal()
ps1b <- ggplot(synthetic,aes(delay_ms,speedup,colour=strategy_label))+
  geom_line(linewidth=.7)+geom_point(size=2)+facet_wrap(~scenario,nrow=1)+
  scale_colour_manual(values=strategy_colours)+scale_y_continuous(labels=scales::label_number(suffix="×"))+
  labs(x="Configured wait per consequence (ms)",y="Observed speedup")+theme_journal()+theme(legend.position="none")
save_plot(ps1a+labs(tag="A"),"supplements/panels/FigureS1A_pipeline_runtime",10.5,3.8)
save_plot(ps1b+labs(tag="B"),"supplements/panels/FigureS1B_pipeline_speedup",10.5,3.4)
save_plot((ps1a/ps1b)+plot_annotation(tag_levels="A"),"supplements/FigureS1_pipeline_complexity",10.5,7.2)

#### SUPPLEMENTARY ASSAYS, RESOURCES AND DIAGNOSTICS ####
ps2 <- ggplot(accepted,aes(speedup,recipe,colour=cohort,shape=cohort))+
  geom_vline(xintercept=1,linetype=2,colour="#A3A9AD",linewidth=.35)+
  geom_point(position=position_jitter(width=0,height=.1,seed=7),size=1.6,alpha=.8)+
  facet_grid(assembly~assay)+ratio_axis()+point_style()+
  labs(x="Per-sample speedup over the same direct recipe",y=NULL,colour=NULL,shape=NULL)+theme_journal()
save_plot(ps2,"supplements/FigureS2_assays_by_assembly",12,7)
resources <- subset(accepted,direct_kind=="new_paired_reference")
resources$cpu_ratio <- resources$cached_cpu_seconds/resources$direct_cpu_seconds
resources$memory_ratio <- resources$cached_peak_gib/resources$direct_peak_gib
long <- rbind(transform(resources,metric="Aggregate CPU time",ratio=cpu_ratio),transform(resources,metric="Aggregate peak memory",ratio=memory_ratio))
ps3 <- ggplot(long,aes(recipe,ratio,colour=assay))+
  geom_hline(yintercept=1,linetype=2,colour="#A3A9AD",linewidth=.35)+
  geom_boxplot(outlier.shape=NA,linewidth=.4)+facet_grid(metric~assembly)+
  scale_y_log10()+labs(x=NULL,y="Cached / direct",colour="Assay")+theme_journal()+theme(axis.text.x=element_text(angle=30,hjust=1))
save_plot(ps3,"supplements/FigureS3_resources",12,7)
diag <- ggplot(subset(assays,tool=="veprs"),aes(speedup,level,colour=cohort,shape=evidence_class))+
  geom_vline(xintercept=1,linetype=2,colour="#A3A9AD",linewidth=.35)+
  geom_point(position=position_jitter(width=0,height=.08,seed=8),size=1.8)+facet_wrap(~assay,nrow=1)+
  ratio_axis()+scale_colour_manual(values=cohort_colours)+
  scale_shape_manual(values=c(validated_with_recorded_scope=16,veprs_clinical_timing_only=1),
    labels=c(validated_with_recorded_scope="Accepted",veprs_clinical_timing_only="Timing only: CSQ differs"))+
  labs(x="Timing ratio (hollow points are not validated speedups)",y=NULL,colour=NULL,shape=NULL)+theme_journal()
save_plot(diag,"diagnostics/Diagnostic_veprs_clinical_not_validated",10.5,3.7)

#### SIMPLE REPOSITORY HEADLINER ####
repo <- subset(ms,assay=="WGS" & recipe %in% c("VEP everything","fastVEP clinical","SnpEff clinical"))
repo$remaining <- 1/repo$median_speedup
headliner <- ggplot(repo,aes(y=recipe))+
  geom_segment(aes(x=0,xend=1,yend=recipe),colour="#E1E7E9",linewidth=11,lineend="round")+
  geom_segment(aes(x=0,xend=remaining,yend=recipe),colour="#347F69",linewidth=11,lineend="round")+
  geom_text(aes(x=1.08,label=sprintf("%.0f%% less time",100*(1-remaining))),hjust=0,size=4)+
  scale_x_continuous(labels=scales::label_percent(),limits=c(0,1.65),breaks=c(0,.5,1))+
  labs(x="Runtime relative to direct annotation",y=NULL)+theme_journal(12)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
save_plot(headliner,"repository/WGS_time_saved",10,3.5)

#### REUSABLE RSTUDIO TIBBLES AND SESSION RECORD ####
saveRDS(list(measurements=cells,plot_views=views,assays=assays,validated_assays=accepted,
  primary_grch38=primary,wgs_strategies=wgs,synthetic=synthetic,assay_summary=assay_summary,
  primary_summary=primary_summary,wgs_summary=wgs_summary,resources=resources),file.path(root,"tables/manuscript_data.rds"))
writeLines(capture.output(sessionInfo()),file.path(root,"source/R_sessionInfo.txt"))
print(primary_summary,width=Inf)
