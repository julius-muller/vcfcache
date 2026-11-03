#!/usr/bin/env Rscript
# Measured vep-rs evidence and explicitly hypothetical sensitivity analysis.
# All output is new. Nothing in the frozen final dataset is overwritten.

#### PATHS AND PACKAGES ####
args <- commandArgs(TRUE)
if(length(args)<2L || length(args)>3L) stop("Usage: extend_analysis.R <frozen snapshot> <output> [new figure directory]")
root <- normalizePath(args[[1]],mustWork=TRUE)
out <- args[[2]]
figure_dir<-if(length(args)==3L) args[[3]] else "figures"
.libPaths(c(file.path(path.expand("~"),".local/lib/R/4.5"),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(tibble)})
if(dir.exists(file.path(out,figure_dir))) stop("Use a new figure directory")
dir.create(file.path(out,"tables"),recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,name) write.table(x,file.path(out,"tables",name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
save_plot <- function(p,name,w=10,h=4) {
  prefix<-file.path(out,figure_dir,name);dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
  stopifnot(!any(file.exists(paste0(prefix,c(".pdf",".svg",".png")))) )
  cairo_pdf(paste0(prefix,".pdf"),width=w,height=h);print(p);dev.off()
  svg(paste0(prefix,".svg"),width=w,height=h);print(p);dev.off()
  png(paste0(prefix,".png"),width=w,height=h,units="in",res=600,type="cairo",bg="white");print(p);dev.off()
}
theme_journal<-function() theme_classic(base_size=10)+theme(
  strip.background=element_blank(),strip.text=element_text(face="bold"),
  axis.text=element_text(colour="#24343E"),plot.tag=element_text(face="bold",size=13),
  plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),
  legend.position="bottom",plot.margin=margin(8,10,8,8))
colours<-c(KPGP="#306B8D",SGDP="#A36037")
data<-readRDS(file.path(root,"tables/manuscript_data.rds"))

#### VEP-RS CORE AND ANNOTATION-RICH RECIPES: MEASURED TIMINGS ####
# Clinical WES/WGS records have unresolved CSQ differences. They are displayed
# explicitly as timing evidence, not promoted into the validated headline data.
x<-subset(data$assays,tool=="veprs")
stopifnot(nrow(x)==72,sum(!x$validated)==24)
x$recipe_label<-factor(ifelse(x$level=="core","Core","Clinical: CSQ differs in WES/WGS"),
  levels=c("Core","Clinical: CSQ differs in WES/WGS"))
x$equivalence<-factor(ifelse(x$validated,"Accepted comparison","Timing only: CSQ differs"),
  levels=c("Accepted comparison","Timing only: CSQ differs"))
summary<-do.call(rbind,lapply(split(seq_len(nrow(x)),interaction(x$level,x$assay,drop=TRUE)),function(i) {
  z<-x[i,];data.frame(recipe=z$level[1],assay=z$assay[1],n=nrow(z),
    median_direct_seconds=median(z$direct_seconds),median_cached_seconds=median(z$cached_seconds),
    median_timing_ratio=median(z$speedup),min_ratio=min(z$speedup),max_ratio=max(z$speedup),
    accepted=sum(z$validated),median_prior_rust_seconds=median(z$previous_rust_cached_seconds),
    median_prior_rust_ratio=median(z$direct_seconds/z$previous_rust_cached_seconds),
    median_cached_runtime_cut_for_parity=median(pmax(0,1-z$speedup)))
}))
write_tsv(summary,"veprs_core_and_clinical_summary.tsv")
write_tsv(x,"veprs_all_conditions.tsv")

pa<-ggplot(x,aes(direct_seconds,cached_seconds,colour=cohort,shape=equivalence))+
  geom_abline(slope=1,intercept=0,linetype=2,colour="#9AA4AA",linewidth=.35)+
  geom_point(size=2)+facet_grid(recipe_label~assay,scales="free")+
  scale_x_log10()+scale_y_log10()+scale_colour_manual(values=colours)+
  scale_shape_manual(values=c(16,1))+labs(x="Direct elapsed time (s)",y="Cached elapsed time (s)",colour=NULL,shape=NULL)+theme_journal()
pb<-ggplot(x,aes(speedup,recipe_label,colour=cohort,shape=equivalence))+
  geom_vline(xintercept=1,linetype=2,colour="#9AA4AA",linewidth=.35)+
  geom_point(position=position_jitter(width=0,height=.09,seed=20261103),size=2)+
  facet_wrap(~assay,nrow=1)+scale_x_log10(labels=scales::label_number(suffix="×"))+
  scale_colour_manual(values=colours)+scale_shape_manual(values=c(16,1))+
  labs(x="Direct / cached elapsed time (hollow points are timing evidence only)",y=NULL,colour=NULL,shape=NULL)+theme_journal()
save_plot(pa+labs(tag="A"),"provisional/panels/Figure4A_veprs_elapsed_times",11,6)
save_plot(pb+labs(tag="B"),"provisional/panels/Figure4B_veprs_timing_ratios",11,3.6)
save_plot((pa/pb)+plot_layout(heights=c(1.5,1))+plot_annotation(tag_levels="A"),"provisional/Figure4_veprs_core_and_annotation_rich",11,9)

#### EXPLICIT HYPOTHETICAL RUNTIME MODEL ####
# a and c are SERIAL costs per variant, not observed costs divided by N at
# eight CPUs. Efficiency is modeled by an Amdahl serial fraction sigma.
# Startup is paid once in each arm whenever misses exist, not multiplied by
# the miss fraction. Cache setup O and per-input-variant routing c remain.
predict_runtime<-function(N=5e6,f=.8,a_us=100,c_us=20,startup=5,O=1,p=8,
                         sigma_ann=0,sigma_cache=0) {
  ann_scale<-1/ (sigma_ann+(1-sigma_ann)/p)
  cache_scale<-1/ (sigma_cache+(1-sigma_cache)/p)
  direct<-startup+N*a_us*1e-6/ann_scale
  miss_work<-startup+(1-f)*N*a_us*1e-6/ann_scale
  miss<-ifelse(rep_len(f,length(miss_work))<1,miss_work,0)
  cached<-O+N*c_us*1e-6/cache_scale+miss
  data.frame(direct=direct,cached=cached,speedup=direct/cached,
    eta_ann=ann_scale/p,eta_cache=cache_scale/p)
}
# Algebraic sanity checks, including the discontinuous all-hit startup shortcut.
stopifnot(predict_runtime(N=100,f=0,O=0,c_us=0)$speedup==1,
          predict_runtime(N=100,f=1,O=1,c_us=0)$cached==1,
          predict_runtime(p=1,sigma_ann=.5)$eta_ann==1)
# Vector inputs must agree with separate scalar calls; a scalar ifelse test
# would otherwise silently recycle just the first miss-annotation runtime.
check<-predict_runtime(p=c(1,8,32))
stopifnot(isTRUE(all.equal(check$speedup,vapply(c(1,8,32),function(p) predict_runtime(p=p)$speedup,numeric(1)))))
check<-predict_runtime(a_us=c(1,100,1000))
stopifnot(isTRUE(all.equal(check$speedup,vapply(c(1,100,1000),function(a) predict_runtime(a_us=a)$speedup,numeric(1)))))
scenario<-function(d) {
  defaults<-list(N=5e6,f=.8,a_us=100,c_us=20,startup=5,O=1,p=8,sigma_ann=0,sigma_cache=0)
  for(k in names(d)) if(k %in% names(defaults)) defaults[[k]]<-d[[k]]
  cbind(d,do.call(predict_runtime,defaults))
}
coverage_colours<-c("50%"="#45738F","80%"="#32816A","93%"="#AB6C37")
ma<-scenario(expand.grid(a_us=10^seq(0,3,length.out=121),f=c(.5,.8,.93),N=c(82714,5e6)))
mb<-scenario(expand.grid(startup=seq(0,120,length.out=121),f=c(.5,.8,.93),N=c(82714,5e6)))
for(name in c("ma","mb")) {
  d<-get(name);d$coverage<-factor(paste0(100*d$f,"%"),levels=names(coverage_colours))
  d$input<-ifelse(d$N<1e6,"WES-sized: 82,714 variants","WGS-sized: 5 million variants")
  assign(name,d)
}
mA<-ggplot(ma,aes(a_us,speedup,colour=coverage))+
  geom_hline(yintercept=1,linetype=2,colour="#9AA4AA",linewidth=.35)+geom_line(linewidth=.7)+
  facet_wrap(~input,nrow=1)+scale_x_log10()+scale_y_log10(labels=scales::label_number(suffix="×"))+
  scale_colour_manual(values=coverage_colours)+labs(x="Serial annotation cost per variant (µs)",y="Modeled speedup",colour="Cache hits")+theme_journal()
mB<-ggplot(mb,aes(startup,speedup,colour=coverage))+
  geom_hline(yintercept=1,linetype=2,colour="#9AA4AA",linewidth=.35)+geom_line(linewidth=.7)+
  facet_wrap(~input,nrow=1)+scale_y_continuous(labels=scales::label_number(suffix="×"))+
  scale_colour_manual(values=coverage_colours)+labs(x="Annotator startup per invocation (s)",y="Modeled speedup",colour="Cache hits")+theme_journal()
mc<-expand.grid(p=1:32,case=c("Equal ideal scaling","Cache scales poorly","Annotator scales poorly"))
mc$sigma_ann<-ifelse(mc$case=="Annotator scales poorly",.2,0)
mc$sigma_cache<-ifelse(mc$case=="Cache scales poorly",.2,0)
mc<-scenario(mc)
mC<-ggplot(mc,aes(p,speedup,colour=case,linetype=case))+
  geom_hline(yintercept=1,linetype=2,colour="#9AA4AA",linewidth=.35)+geom_line(linewidth=.7)+
  scale_y_continuous(labels=scales::label_number(suffix="×"))+
  scale_x_continuous(breaks=c(1,8,16,32))+
  scale_colour_manual(values=c("Equal ideal scaling"="#306B8D","Cache scales poorly"="#32816A","Annotator scales poorly"="#AB6C37"))+
  labs(x="Configured CPUs (model, not a CPU-scaling experiment)",y="Modeled speedup",colour=NULL,linetype=NULL)+theme_journal()
md<-scenario(expand.grid(a_us=10^seq(0,3,length.out=101),c_us=10^seq(0,2.5,length.out=101)))
mD<-ggplot(md,aes(a_us,c_us))+
  geom_raster(aes(fill=log2(speedup)))+
  geom_contour(aes(z=speedup),breaks=1,colour="black",linewidth=.45)+
  scale_x_log10()+scale_y_log10()+scale_fill_gradient2(low="#A45443",mid="white",high="#2C7864",midpoint=0,
    breaks=c(-2,0,2),labels=c("0.25×","1×","4×"))+
  labs(x="Serial annotation cost per variant (µs)",y="Serial cache-path cost per variant (µs)",fill="Modeled speedup")+theme_journal()
save_plot(mA+labs(tag="A"),"model/panels/Figure5A_annotation_cost",10.5,3.6)
save_plot(mB+labs(tag="B"),"model/panels/Figure5B_startup_cost",10.5,3.6)
save_plot(mC+labs(tag="C"),"model/panels/Figure5C_parallel_scaling",6,4)
save_plot(mD+labs(tag="D"),"model/panels/Figure5D_break_even",6,4)
save_plot((mA/mB/(mC|mD))+plot_layout(heights=c(1,1,1.25))+plot_annotation(tag_levels="A"),"model/Figure5_runtime_model_sensitivity",11,11)
write_tsv(ma,"model_annotation_cost.tsv");write_tsv(mb,"model_startup_cost.tsv")
write_tsv(mc,"model_parallelization.tsv");write_tsv(md,"model_break_even.tsv")
saveRDS(list(measured_veprs=x,measured_summary=as_tibble(summary),model_annotation=as_tibble(ma),
  model_startup=as_tibble(mb),model_parallelization=as_tibble(mc),model_break_even=as_tibble(md)),file.path(out,"tables/interpretation_data.rds"))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"R_sessionInfo.txt"))
print(summary,row.names=FALSE)
