#!/usr/bin/env Rscript
#### PATHS AND PACKAGES ####
args<-commandArgs(TRUE);stopifnot(length(args)==1)
root<-normalizePath(args[[1]],mustWork=TRUE)
.libPaths(c(file.path(path.expand("~"),".local/lib/R/4.5"),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(tibble)})
out<-file.path(root,"figures");if(dir.exists(out))stop("Choose a new output directory")
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
names<-c(adaptive="Rust-adaptive",resident="Resident input",coalesced="Coalesced input")
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
save_plot(p1,"Figure1_cold_CLI_WES",10.5,4)

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
