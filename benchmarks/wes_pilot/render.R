#!/usr/bin/env Rscript
#### PATHS AND PACKAGES ####
args<-commandArgs(TRUE);stopifnot(length(args)>=1,length(args)<=2)
root<-normalizePath(args[[1]],mustWork=TRUE)
.libPaths(c(file.path(path.expand("~"),".local/lib/R/4.5"),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(tibble)})
out<-file.path(root,if(length(args)==2)args[[2]] else "figures");if(dir.exists(out))stop("Choose a new output directory")
read_tsv<-function(name) as_tibble(read.delim(file.path(root,"tables",name),quote="",comment.char="",na.strings=c("","NA")))
write_tsv<-function(x,name)write.table(x,file.path(root,"tables",name),sep="\t",quote=FALSE,row.names=FALSE)
save_plot<-function(p,name,w=10,h=5){
 prefix<-file.path(out,name);dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
 cairo_pdf(paste0(prefix,".pdf"),width=w,height=h);print(p);dev.off()
 svg(paste0(prefix,".svg"),width=w,height=h);print(p);dev.off()
 png(paste0(prefix,".png"),width=w,height=h,units="in",res=600,type="cairo",bg="white");print(p);dev.off()
}
theme_journal<-function()theme_classic(base_size=10)+theme(strip.background=element_blank(),strip.text=element_text(face="bold"),legend.position="bottom",plot.tag=element_text(face="bold",size=13),plot.title=element_blank(),plot.caption=element_blank())

#### COLD APPLICATION-START CLI COMPARISONS ####
x<-read_tsv("cli_conditions.tsv")
if(file.exists(file.path(root,"tables/worker_tuning.tsv"))){
 tuning<-read_tsv("worker_tuning.tsv")
 w8<-subset(tuning,workers==8)
 # A single globally selected configuration, never the best arm per sample.
 stopifnot(nrow(w8)==24,!anyDuplicated(w8$id))
 fields<-union(names(x),names(w8))
 for(n in setdiff(fields,names(x)))x[[n]]<-NA
 for(n in setdiff(fields,names(w8)))w8[[n]]<-NA
 x<-rbind(x[,fields],w8[,fields])
}
names<-c(adaptive="Rust-adaptive",resident="Resident input: 4 workers",coalesced="Coalesced input: 4 workers",`resident-w8`="Resident input: 8 workers")
x<-subset(x,mode!="direct");x$backend<-factor(names[x$mode],levels=unname(names))
x$recipe<-ifelse(x$tool=="veprs","vep-rs core","fastVEP clinical")
x$cohort<-toupper(x$cohort)
summary<-do.call(rbind,lapply(split(x,interaction(x$tool,x$mode,drop=TRUE)),function(d)data.frame(tool=d$tool[1],backend=as.character(d$backend[1]),n=nrow(d),median_ratio=median(d$ratio),min_ratio=min(d$ratio),max_ratio=max(d$ratio),n_faster=sum(d$ratio>=1),n_clear_wins=sum(d$ratio>=1.05),median_seconds=median(d$cached_seconds))))
write_tsv(summary,"cli_summary.tsv")
p1<-ggplot(x,aes(ratio,backend,colour=cohort,shape=cohort))+
 geom_vline(xintercept=1,linetype=2,colour="#929DA5",linewidth=.35)+
 geom_point(position=position_jitter(width=0,height=.1,seed=19),size=2)+facet_wrap(~recipe,nrow=1)+
 scale_x_log10(labels=scales::label_number(suffix="×"))+scale_colour_manual(values=c(KPGP="#306B8D",SGDP="#A36037"))+
 labs(x="Direct / cached CLI elapsed time",y=NULL,colour=NULL,shape=NULL)+theme_journal()
save_plot(p1+labs(tag="A"),"panels/Figure1A_cold_CLI_WES",10.5,4)
if(all(c("annotation_seconds","lookup_seconds") %in% names(x))){
 z<-subset(x,mode %in% c("resident","coalesced","resident-w8"))
 z$other_seconds<-z$cached_seconds-z$annotation_seconds-z$lookup_seconds
 stopifnot(all(z$other_seconds>=0))
 stages<-rbind(transform(z,phase="Miss-recipe invocation",fraction=annotation_seconds/cached_seconds),
  transform(z,phase="Cache lookup",fraction=lookup_seconds/cached_seconds),
  transform(z,phase="Other CLI / I/O",fraction=other_seconds/cached_seconds))
 stages<-aggregate(fraction~recipe+backend+phase,stages,mean)
 p1b<-ggplot(stages,aes(fraction,backend,fill=phase))+geom_col(width=.65)+facet_wrap(~recipe,nrow=1)+
  scale_x_continuous(labels=scales::label_percent())+scale_fill_manual(values=c("Miss-recipe invocation"="#76608D","Cache lookup"="#306B8D","Other CLI / I/O"="#B08A48"))+
  labs(x="Mean fraction of cached elapsed time",y=NULL,fill=NULL)+theme_journal()
 save_plot(p1b+labs(tag="B"),"panels/Figure1B_measured_time_fractions",10.5,3.5)
 save_plot((p1/p1b)+plot_annotation(tag_levels="A"),"Figure1_cold_CLI_WES",10.5,7.5)
}else save_plot(p1,"Figure1_cold_CLI_WES",10.5,4)

#### GENERIC LOOKUP WORKER ALLOCATION: SAME EIGHT-CPU BUDGET ####
if(exists("tuning")){
 # Show every paired confirmation, including regressions. No sample-specific
 # or assay-specific dispatch policy is introduced by this engineering pilot.
 paired<-subset(x,mode %in% c("resident","resident-w8"))
 paired$configuration<-factor(ifelse(paired$mode=="resident","4 lookup + 4 codec","8 lookup + 0 codec"),levels=c("4 lookup + 4 codec","8 lookup + 0 codec"))
 p3<-ggplot(paired,aes(configuration,cached_seconds,group=sample,colour=cohort,shape=cohort))+
  geom_line(linewidth=.4,alpha=.65)+geom_point(size=2)+facet_wrap(~recipe,nrow=1)+
  scale_colour_manual(values=c(KPGP="#306B8D",SGDP="#A36037"))+
  labs(x=NULL,y="Cached CLI elapsed time (s)",colour=NULL,shape=NULL)+theme_journal()
 save_plot(p3,"Figure2_paired_worker_allocation",10.5,4)
 worker_summary<-do.call(rbind,lapply(split(subset(tuning,workers==8),subset(tuning,workers==8)$tool),function(d)data.frame(tool=d$tool[1],n=nrow(d),median_gain_over_resident4=median(d$gain_over_resident4),n_faster_than_resident4=sum(d$gain_over_resident4>1))))
 write_tsv(worker_summary,"worker_summary.tsv")
}

#### RESIDENT SESSION: WARM LATENCY IS NOT COLD PERFORMANCE ####
if(file.exists(file.path(root,"tables/resident_warm.tsv"))){
 w<-read_tsv("resident_warm.tsv");c<-read_tsv("resident_cumulative.tsv")
 w$recipe<-ifelse(w$tool=="veprs","vep-rs core","fastVEP clinical")
 c$recipe<-ifelse(c$tool=="veprs","vep-rs core","fastVEP clinical")
 p2a<-ggplot(w,aes(warm_ratio,recipe))+geom_vline(xintercept=1,linetype=2,colour="#929DA5",linewidth=.35)+
 geom_point(position=position_jitter(width=0,height=.07,seed=20),colour="#306B8D",size=2)+
 scale_x_log10(labels=scales::label_number(suffix="×"))+labs(x="Warm resident ratio versus equally batched direct API",y=NULL)+theme_journal()
 p2b<-ggplot(c,aes(samples,ratio,colour=recipe))+geom_hline(yintercept=1,linetype=2,colour="#929DA5",linewidth=.35)+geom_line(linewidth=.7)+geom_point(size=1.5)+
 scale_y_log10(labels=scales::label_number(suffix="×"))+scale_x_continuous(breaks=c(1,4,8,12))+
 labs(x="Successive WES samples",y="Cumulative speedup including loading",colour=NULL)+theme_journal()
 save_plot(p2a+labs(tag="A"),"panels/Figure2A_resident_warm",10,3.5)
 save_plot(p2b+labs(tag="B"),"panels/Figure2B_resident_amortization",10,3.5)
 save_plot((p2a/p2b)+plot_annotation(tag_levels="A"),"Figure2_resident_session",10,7)
 ws<-do.call(rbind,lapply(split(w,w$tool),function(d)data.frame(tool=d$tool[1],n=nrow(d),median_warm_ratio=median(d$warm_ratio),n_faster=sum(d$warm_ratio>=1),loader_seconds=d$loader_seconds[1],session_overhead_seconds=d$session_overhead_seconds[1],session_peak_gib=d$session_peak_gib[1])))
 write_tsv(ws,"resident_summary.tsv")
 saveRDS(list(cli=x,cli_summary=as_tibble(summary),warm=w,cumulative=c,warm_summary=as_tibble(ws)),file.path(root,"tables/pilot_data.rds"))
}else saveRDS(list(cli=x,cli_summary=as_tibble(summary)),file.path(root,"tables/pilot_data.rds"))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(root,"R_sessionInfo.txt"))
print(summary,row.names=FALSE)
