#!/usr/bin/env Rscript
# Final publication-campaign analysis. Plotting is exclusively R / ggplot2.
# Source evidence is immutable; rerenders require a NEW figure directory.

#### PATHS AND PACKAGES ####

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: analyse_final.R <snapshot_directory> [new_figure_subdirectory]")
root <- normalizePath(args[[1]], mustWork = TRUE)
figure_name <- if (length(args) >= 2L) args[[2]] else "figures"
figure_root <- file.path(root, figure_name)
script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_argument[[1]])))
.libPaths(c(Sys.getenv("VCFCACHE_R_LIBRARY", file.path(path.expand("~"), ".local/lib/R/4.5")), .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
source(file.path(script_dir, "final_analysis_helpers.R"))
if (!file.exists(file.path(root, "source/SNAPSHOT.json"))) stop("Source snapshot has not completed")
if (dir.exists(figure_root)) stop("Use a new figure subdirectory; existing figures are immutable")

#### DATA LOADING AND VALIDATION ####

cells <- enrich(read_table(file.path(root, "tables/cells.tsv")))
views <- enrich(read_table(file.path(root, "tables/plot_views.tsv")))
attempts <- read_table(file.path(root, "tables/attempts.tsv"))
stopifnot(nrow(cells) == 600L, nrow(views) == 612L,
          !anyDuplicated(cells$measurement_id), sum(views$is_alias) == 12L,
          sum(cells$validated) == 572L)

# These are derived views of the SAME measurements, never extra samples.
assay_families <- c("vep_matched_assays", "native_clinical_assays", "veprs_extension")
assays <- subset(cells, family %in% assay_families)
validated_assays <- subset(assays, validated)
primary <- subset(validated_assays, assembly == "GRCh38")
wgs <- subset(views, family == "full52_wgs_strategies" & validated)
synthetic <- subset(cells, family == "synthetic_pipeline_spectrum")
stopifnot(nrow(wgs) == 312L, all(table(wgs$tool, wgs$strategy) == 52L))

# Identical recipe names do not imply identical biological information across
# annotators or assemblies. Every ratio uses its own sample/recipe baseline.
assay_summary <- summarise_conditions(validated_assays, c("recipe", "assay", "assembly"))
primary_summary <- summarise_conditions(primary, c("recipe", "assay", "assembly"))
wgs_summary <- summarise_conditions(wgs, c("recipe", "assembly", "cohort", "strategy_label"))
all_summary <- summarise_conditions(views, c("family", "recipe", "level", "assay", "assembly", "strategy_label", "evidence_class"))
cohort_summary <- summarise_conditions(validated_assays, c("recipe", "assay", "assembly", "cohort"))
write_table(assay_summary, file.path(root, "tables/assay_summary_validated.tsv"))
write_table(primary_summary, file.path(root, "tables/primary_grch38_summary.tsv"))
write_table(wgs_summary, file.path(root, "tables/wgs_strategy_summary.tsv"))
write_table(all_summary, file.path(root, "tables/all_summary_by_evidence.tsv"))
write_table(cohort_summary, file.path(root, "tables/assay_summary_by_cohort.tsv"))

#### OPTIMIZATION HEADROOM: MEASURED PHASES, NOT INFERRED STARTUP ####

phase_long <- phase_data(primary)
phase_summary <- aggregate(fraction ~ recipe + assay + phase, phase_long, mean)
stopifnot(max(abs(aggregate(fraction ~ recipe + assay, phase_summary, sum)$fraction - 1)) < 1e-6)
write_table(phase_summary, file.path(root, "tables/phase_fraction_summary.tsv"))
phase_seconds <- aggregate(cbind(cli_other_seconds, preprocess_seconds, route_seconds,
  assemble_seconds, engine_other_seconds, annotate_misses_seconds, cached_seconds) ~
  recipe + assay, primary, median)
write_table(phase_seconds, file.path(root, "tables/phase_seconds_median.tsv"))

# Optimistic counterfactuals hold miss annotation fixed and remove known costs.
# They are NOT measured speedups and cannot be used as publication benchmarks.
headroom <- do.call(rbind, lapply(c(0, .25, .5, .75, 1), function(q) {
  x <- primary
  x$removed_non_annotation_fraction <- q
  x$counterfactual_speedup <- x$direct_seconds /
    (x$annotate_misses_seconds + (1-q) * x$cache_side_residual_seconds)
  x
}))
headroom_summary <- aggregate(counterfactual_speedup ~ recipe + assay +
  removed_non_annotation_fraction, headroom, median)
write_table(headroom_summary, file.path(root, "tables/cache_side_what_if.tsv"))
one_stage <- primary
one_stage$zero_route_speedup <- one_stage$direct_seconds /
  (one_stage$cached_seconds - one_stage$route_seconds)
one_stage$zero_route_assemble_speedup <- one_stage$direct_seconds /
  (one_stage$cached_seconds - one_stage$route_seconds - one_stage$assemble_seconds)
one_stage_summary <- aggregate(cbind(zero_route_speedup, zero_route_assemble_speedup,
  required_cache_side_reduction_for_parity) ~ recipe + assay, one_stage, median)
write_table(one_stage_summary, file.path(root, "tables/individual_stage_what_if.tsv"))
goals <- do.call(rbind,lapply(c(1,2,3,4),function(target) {
  x <- primary
  x$target_speedup <- target
  x$required_side_reduction <- pmax(0,(x$cached_seconds-x$direct_seconds/target) /
    x$cache_side_residual_seconds)
  x$attainable_if_annotation_unchanged <- x$annotate_misses_seconds < x$direct_seconds/target
  x
}))
goal_summary <- aggregate(cbind(required_side_reduction,attainable_if_annotation_unchanged) ~
  recipe+assay+target_speedup,goals,median)
write_table(goal_summary,file.path(root,"tables/target_speedup_headroom.tsv"))

#### PAIRED HISTORICAL COMPARISON AND CACHE STRATEGY CONTRASTS ####

# Historical-to-current cached differences are not a fresh randomized A/B trial.
# The direct denominators, annotation settings and resource controls are frozen,
# but execution date/node/OS-cache state can still affect elapsed differences.
historical <- subset(cells, validated & family != "veprs_extension" &
                       family != "synthetic_pipeline_spectrum")
historical_summary <- summarise_conditions(historical,
  c("family", "recipe", "assay", "assembly", "strategy_label"))
write_table(historical_summary, file.path(root, "tables/historical_comparison.tsv"))
contrasts <- list()
for (recipe in unique(as.character(wgs$recipe))) {
  for (cohort in levels(wgs$cohort)) {
    x <- wgs[as.character(wgs$recipe) == recipe & as.character(wgs$cohort) == cohort, ]
    if (!nrow(x)) next
    base <- subset(x, strategy == "gnomad_af_0.1")
    for (candidate in c("gnomad_af_0.01", "cohort_3_genomes")) {
      b <- subset(x, strategy == candidate)
      pair <- merge(base, b, by = "sample", suffixes = c("_af10", "_candidate"))
      stopifnot(nrow(pair) == nrow(base), !anyDuplicated(pair$sample))
      ratio <- pair$cached_seconds_af10 / pair$cached_seconds_candidate
      ci <- bootstrap_interval(ratio, rep(cohort, length(ratio)))
      contrasts[[length(contrasts)+1L]] <- data.frame(recipe = recipe, cohort = cohort,
        candidate = candidate, n = length(ratio), median_cached_gain_vs_af10 = median(ratio),
        ci_low = ci[1], ci_high = ci[2], fraction_faster = mean(ratio > 1),
        median_hit_increase = median(pair$cache_hit_rate_candidate-pair$cache_hit_rate_af10))
    }
  }
}
contrasts <- tibble::as_tibble(do.call(rbind, contrasts))
write_table(contrasts, file.path(root, "tables/paired_strategy_contrasts.tsv"))

#### RESOURCE AND OUTLIER TABLES ####

resources <- subset(validated_assays, !is.na(direct_peak_gib) & direct_kind == "new_paired_reference")
resources$cpu_cost_ratio <- resources$cached_cpu_seconds / resources$direct_cpu_seconds
resources$memory_ratio <- resources$cached_peak_gib / resources$direct_peak_gib
resource_summary <- aggregate(cbind(cached_peak_gib, direct_peak_gib, cpu_cost_ratio,
  memory_ratio) ~ recipe + assay + assembly, resources, median)
write_table(resources, file.path(root, "tables/paired_resources.tsv"))
write_table(resource_summary, file.path(root, "tables/resource_summary.tsv"))
worst <- do.call(rbind, lapply(split_groups(validated_assays, c("recipe", "assay", "assembly")), function(i) {
  x <- validated_assays[i, ]; x[order(x$speedup)[seq_len(min(3, nrow(x)))], ]
}))
write_table(worst, file.path(root, "tables/lowest_speedup_samples.tsv"))
write_table(subset(cells, !validated), file.path(root, "tables/validation_exclusions.tsv"))

#### MANUSCRIPT FIGURE 1: ASSAY PERFORMANCE AND ABSOLUTE SAVINGS ####

primary_recipes <- c("VEP everything", "fastVEP clinical", "SnpEff clinical", "vep-rs core")
main <- subset(primary, recipe %in% primary_recipes)
main$recipe <- factor(as.character(main$recipe), levels = rev(primary_recipes))
ms <- subset(primary_summary, recipe %in% primary_recipes)
ms$recipe <- factor(as.character(ms$recipe), levels = rev(primary_recipes))
position <- position_jitter(width = 0, height = .09, seed = 10)
p1a <- ggplot(main, aes(speedup, recipe)) +
  geom_vline(xintercept = 1, colour = "#9A9FA3", linetype = 2, linewidth = .35) +
  geom_point(aes(colour = cohort, shape = cohort), position = position, alpha = .75, size = 1.6) +
  geom_errorbar(data = ms, aes(xmin = speedup_ci_low, xmax = speedup_ci_high,
    x = median_speedup), orientation = "y", width = .24, linewidth = .65) +
  geom_point(data = ms, aes(x = median_speedup), shape = 18, size = 2.8) +
  scale_x_log10(breaks = c(.25, .5, 1, 2, 4, 8, 16), labels = scales::label_number(suffix = "×")) +
  scale_colour_manual(values = cohort_colours) + scale_shape_manual(values = c(16, 17, 15)) +
  facet_wrap(~assay, nrow = 1) + labs(x = "Speedup over the same recipe's direct run", y = NULL,
    colour = NULL, shape = NULL) + journal_theme()
p1b <- ggplot(subset(main, assay == "WGS"), aes(saved_hours, recipe)) +
  geom_vline(xintercept = 0, colour = "#9A9FA3", linewidth = .35) +
  geom_boxplot(outlier.shape = NA, width = .42, colour = "#33414B", fill = "#EDF1F3", linewidth = .4) +
  geom_point(aes(colour = cohort, shape = cohort), position = position, size = 1.8, alpha = .85) +
  scale_colour_manual(values = cohort_colours) + scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = "Time saved per WGS sample (h)", y = NULL, colour = NULL, shape = NULL) + journal_theme()
save_figure(p1a + labs(tag = "A"), file.path(figure_root,"manuscript/panels/Figure1A_assay_speedups"), 10.5, 3.8)
save_figure(p1b + labs(tag = "B"), file.path(figure_root,"manuscript/panels/Figure1B_wgs_time_saved"), 7, 3.2)
save_figure((p1a / p1b) + plot_layout(heights = c(1.2,1)) + plot_annotation(tag_levels = "A"),
  file.path(figure_root,"manuscript/Figure1_assay_performance"), 10.5, 7)

#### MANUSCRIPT FIGURE 2: PAIRED WGS CACHE STRATEGIES ####

# Hit rates are displayed for one annotator rather than duplicated as if extra genomes.
hit <- subset(wgs, tool == "vep")
p2a <- ggplot(hit, aes(strategy_label, cache_hit_rate, colour = strategy_label)) +
  geom_boxplot(width = .5, outlier.shape = NA, linewidth = .4) +
  geom_point(position = position_jitter(width = .13, height = 0, seed = 11), size = 1.4, alpha = .75) +
  facet_wrap(~cohort, nrow = 1) + scale_colour_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0,1)) +
  scale_x_discrete(labels = c("AF ≥ 10%", "AF ≥ 1%", "Three genomes")) +
  labs(x = NULL, y = "Cache hit rate") + journal_theme() + theme(legend.position = "none")
p2b <- ggplot(wgs, aes(speedup, colour = strategy_label, linetype = strategy_label)) +
  stat_ecdf(linewidth = .65) + geom_vline(xintercept = 1, colour = "#A5AAAE", linetype = 2, linewidth = .3) +
  facet_grid(recipe ~ cohort) + scale_colour_manual(values = strategy_colours) +
  scale_x_log10(labels = scales::label_number(suffix = "×")) +
  scale_y_continuous(labels = scales::label_percent(), breaks = c(0,.5,1)) +
  labs(x = "Per-sample speedup over the same direct baseline", y = "Cumulative fraction of genomes",
    colour = NULL, linetype = NULL) + journal_theme()
save_figure(p2a + labs(tag = "A"), file.path(figure_root,"manuscript/panels/Figure2A_cache_hits"), 10.5, 3.2)
save_figure(p2b + labs(tag = "B"), file.path(figure_root,"manuscript/panels/Figure2B_speedup_ecdf"), 10.5, 5)
save_figure((p2a / p2b) + plot_layout(heights = c(1,1.6)) + plot_annotation(tag_levels = "A"),
  file.path(figure_root,"manuscript/Figure2_wgs_cache_strategies"), 10.5, 8.2)

#### MANUSCRIPT FIGURE 3: WHERE THE REMAINING TIME GOES ####

phase_colours <- c("Annotate misses"="#776094", "Lookup / routing"="#316E8E",
  "Preprocessing"="#A9B5BB", "Assembly"="#38816F", "Other engine"="#D8DFE2", "CLI / other"="#C98A43")
ps <- subset(phase_summary, recipe %in% primary_recipes)
ps$recipe <- factor(as.character(ps$recipe), levels = rev(primary_recipes))
p3a <- ggplot(ps, aes(fraction, recipe, fill = phase)) + geom_col(width = .7) +
  facet_wrap(~assay, nrow = 1) + scale_fill_manual(values = phase_colours) +
  scale_x_continuous(labels = scales::label_percent(), breaks = c(0,.5,1)) +
  labs(x = "Fraction of cached elapsed time", y = NULL, fill = NULL) + journal_theme()
wh <- subset(headroom_summary, recipe %in% primary_recipes & assay %in% c("WES","WGS"))
p3b <- ggplot(wh, aes(removed_non_annotation_fraction, counterfactual_speedup,
                      colour = recipe, linetype = recipe)) + geom_line(linewidth = .7) +
  geom_hline(yintercept = 1, colour = "#A5AAAE", linetype = 2, linewidth = .3) +
  facet_wrap(~assay, nrow = 1) + scale_x_continuous(labels = scales::label_percent()) +
  scale_y_log10(labels = scales::label_number(suffix = "×")) +
  scale_colour_manual(values = c("VEP everything"="#316E8E", "fastVEP clinical"="#BE7B37",
    "SnpEff clinical"="#766093", "vep-rs core"="#38816F")) +
  labs(x = "Non-annotation time hypothetically eliminated", y = "Counterfactual median speedup",
       colour = NULL, linetype = NULL) + journal_theme()
save_figure(p3a + labs(tag = "A"), file.path(figure_root,"manuscript/panels/Figure3A_measured_phases"), 10.5, 4)
save_figure(p3b + labs(tag = "B"), file.path(figure_root,"manuscript/panels/Figure3B_optimistic_headroom"), 10.5, 4)
save_figure((p3a / p3b) + plot_annotation(tag_levels = "A"),
  file.path(figure_root,"manuscript/Figure3_optimization_headroom"), 10.5, 8)

#### SUPPLEMENTARY FIGURES: SYNTHETIC LOAD, ASSEMBLY AND VALIDITY ####

synthetic$delay_label <- factor(synthetic$delay_us / 1000, levels = sort(unique(synthetic$delay_us/1000)))
synthetic$scenario <- factor(paste0(as.character(synthetic$strategy_label), "\n",
  sprintf("%.1f%% hits", 100*synthetic$cache_hit_rate)), levels = unique(paste0(
    as.character(synthetic$strategy_label), "\n", sprintf("%.1f%% hits",100*synthetic$cache_hit_rate))))
runtime <- rbind(transform(synthetic, route="Direct", hours=direct_seconds/3600),
                 transform(synthetic, route="VCFcache", hours=cached_seconds/3600))
p4a <- ggplot(runtime, aes(delay_label, hours, fill = route)) +
  geom_col(position = position_dodge(width = .75), width = .65) + coord_flip() +
  facet_wrap(~scenario, nrow = 1) + scale_fill_manual(values = c(Direct="#BBC3C8",VCFcache="#347E6A")) +
  labs(x = "Simulated wait per consequence (ms)", y = "Elapsed time (h)", fill = NULL) + journal_theme()
p4b <- ggplot(synthetic, aes(delay_us/1000, speedup, colour = strategy_label)) +
  geom_line(linewidth = .7) + geom_point(size = 2) + facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual(values = strategy_colours) +
  scale_y_continuous(labels = scales::label_number(suffix = "×")) +
  labs(x = "Simulated wait per consequence (ms)", y = "Observed speedup") + journal_theme() +
  theme(legend.position = "none")
save_figure(p4a + labs(tag="A"),file.path(figure_root,"supplements/panels/FigureS1A_pipeline_runtimes"),10.5,4.2)
save_figure(p4b + labs(tag="B"),file.path(figure_root,"supplements/panels/FigureS1B_pipeline_speedup"),10.5,3.5)
save_figure((p4a / p4b) + plot_annotation(tag_levels="A"),file.path(figure_root,"supplements/FigureS1_pipeline_complexity"),10.5,7.7)

all_assay_summary <- summarise_conditions(assays, c("recipe","assay","assembly","evidence_class"))
write_table(all_assay_summary,file.path(root,"tables/assay_all_evidence_summary.tsv"))
p5 <- ggplot(validated_assays, aes(speedup, recipe, colour = cohort, shape = cohort)) +
  geom_vline(xintercept=1, colour="#A5AAAE",linetype=2,linewidth=.3) +
  geom_point(position=position_jitter(width=0,height=.12,seed=15),size=1.5,alpha=.7) +
  facet_grid(assembly ~ assay, scales="free_y", space="free_y") +
  scale_x_log10(labels=scales::label_number(suffix="×")) +
  scale_colour_manual(values=cohort_colours) +
  scale_shape_manual(values=c(KPGP=16,SGDP=17,PGP=15)) +
  labs(x="Speedup over each recipe's own direct baseline",y=NULL,colour=NULL,shape=NULL) + journal_theme()
save_figure(p5,file.path(figure_root,"supplements/FigureS2_all_assays_by_assembly"),12,8.5)

# These excluded clinical timings are diagnostic decision evidence, not a
# validated manuscript-performance panel. Keep them in a separate directory.
vp <- subset(assays,tool=="veprs")
diagnostic <- ggplot(vp,aes(speedup,level,colour=cohort,shape=evidence_class)) +
  geom_vline(xintercept=1,colour="#A5AAAE",linetype=2,linewidth=.3) +
  geom_point(position=position_jitter(width=0,height=.1,seed=22),size=1.8) +
  facet_wrap(~assay,nrow=1) + scale_x_log10(labels=scales::label_number(suffix="×")) +
  scale_colour_manual(values=cohort_colours) +
  scale_shape_manual(values=c(validated=16,veprs_csq_timing_only=1),
    labels=c(validated="Accepted comparison",veprs_csq_timing_only="Timing only: CSQ differs")) +
  labs(x="Observed ratio; hollow points are not validated speedups",y=NULL,colour=NULL,shape=NULL) + journal_theme()
save_figure(diagnostic,file.path(figure_root,"diagnostics/DiagnosticD1_veprs_recipe_timings"),10.5,3.5)

resource_long <- rbind(transform(resources,metric="Peak aggregate memory",ratio=memory_ratio),
                       transform(resources,metric="Aggregate CPU cost",ratio=cpu_cost_ratio))
p6 <- ggplot(resource_long,aes(recipe,ratio,colour=assay)) +
  geom_hline(yintercept=1,colour="#A5AAAE",linetype=2,linewidth=.3) +
  geom_boxplot(outlier.shape=NA,position=position_dodge(width=.75),linewidth=.4) +
  facet_grid(metric ~ assembly,scales="free_x",space="free_x") +
  scale_y_log10(labels=scales::label_number(suffix="×")) +
  scale_colour_manual(values=c(Panel="#316E8E",WES="#BE7B37",WGS="#38816F")) +
  labs(x=NULL,y="Cached / direct",colour="Assay") + journal_theme() +
  theme(axis.text.x=element_text(angle=25,hjust=1))
save_figure(p6,file.path(figure_root,"supplements/FigureS3_cpu_and_memory"),11,7)

p7 <- ggplot(subset(historical,assay=="WGS"),aes(historical_to_new_speedup,recipe,colour=cohort)) +
  geom_vline(xintercept=1,colour="#A5AAAE",linetype=2,linewidth=.3) +
  geom_boxplot(outlier.shape=NA,linewidth=.4) +
  facet_wrap(~strategy_label,nrow=1,scales="free_y") +
  scale_x_log10(labels=scales::label_number(suffix="×")) + scale_colour_manual(values=cohort_colours) +
  labs(x="Historical cached time / current cached time",y=NULL,colour=NULL) + journal_theme()
save_figure(p7,file.path(figure_root,"supplements/FigureS4_historical_cached_comparison"),12,5)

#### SIMPLE REPOSITORY SUMMARY AND EXPLICITLY ASSUMED COHORT MODEL ####

repo <- subset(primary_summary,assay=="WGS" & recipe %in% c("VEP everything","fastVEP clinical","vep-rs core"))
repo$recipe <- factor(as.character(repo$recipe),levels=rev(c("VEP everything","fastVEP clinical","vep-rs core")))
repo$fraction <- 1/repo$median_speedup
p8 <- ggplot(repo,aes(y=recipe)) +
  geom_segment(aes(x=0,xend=1,yend=recipe),colour="#E0E6E8",linewidth=12,lineend="round") +
  geom_segment(aes(x=0,xend=fraction,yend=recipe),colour="#347E6A",linewidth=12,lineend="round") +
  geom_text(aes(x=1.09,label=sprintf("%.0f%% less time",100*(1-fraction))),hjust=0,size=4) +
  scale_x_continuous(labels=scales::label_percent(),breaks=c(0,.5,1),limits=c(0,1.65)) +
  labs(x="Runtime relative to the same pipeline without caching",y=NULL) +
  journal_theme(12) + theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
save_figure(p8,file.path(figure_root,"repository/Repository_WGS_time_reduction"),10,3.6)

model <- expand.grid(samples=c(5,10,20,50,100,200,500,1000),
  build_in_direct_samples=c(0,10,100),recipe=as.character(repo$recipe),stringsAsFactors=FALSE)
model$cached_ratio <- repo$fraction[match(model$recipe,as.character(repo$recipe))]
model$effective_ratio <- model$cached_ratio + model$build_in_direct_samples/model$samples
model$build_label <- factor(paste0("Build cost: ",model$build_in_direct_samples," direct-sample runtimes (assumed)"))
write_table(model,file.path(root,"tables/illustrative_cohort_amortization.tsv"))
p9 <- ggplot(model,aes(samples,effective_ratio,colour=recipe,linetype=recipe)) +
  geom_line(linewidth=.7) + geom_hline(yintercept=1,colour="#A5AAAE",linetype=2,linewidth=.3) +
  facet_wrap(~build_label,nrow=1) + scale_x_log10(breaks=c(5,10,100,1000)) +
  scale_y_log10(labels=scales::label_number(suffix="×")) +
  labs(x="Samples processed",y="Effective cost per sample / direct cost",colour=NULL,linetype=NULL) + journal_theme()
save_figure(p9,file.path(figure_root,"supplements/FigureS5_illustrative_amortization"),12,4)

#### PAYLOAD-BOUND HISTORICAL BUILD COSTS: A CONDITIONAL SCENARIO ####

# audit_final.py rejects build costs that are not bound to the executed payload.
# These historical builds are not hardware-matched new measurements. Downloads
# and blueprint preparation are not measured: no zero-setup claim is possible.
cache_inventory <- read_table(file.path(root,"tables/cache_inventory.tsv"))
comparator_audit <- read_table(file.path(root,"tables/comparator_audit.tsv"))
seed_audit <- read_table(file.path(root,"tables/cohort_seed_audit.tsv"))
recorded_builds <- read_table(file.path(root,"tables/recorded_build_amortization.tsv"))
build_model <- do.call(rbind,lapply(c(5,10,20,50,100,200,500,1000),function(samples) {
  x <- subset(recorded_builds,tool=="fastvep")
  x$samples <- samples
  # Mean over observed sample pairs is the appropriate additive total-time
  # estimator. Medians and repayment counts remain descriptive in the audit TSV.
  x$mean_direct_seconds <- mapply(function(cohort,strategy) {
    mean(wgs$direct_seconds[as.character(wgs$cohort)==toupper(cohort) &
      wgs$strategy==strategy & wgs$tool=="fastvep"])
  },x$cohort,x$strategy)
  x$mean_cached_seconds <- mapply(function(cohort,strategy) {
    mean(wgs$cached_seconds[as.character(wgs$cohort)==toupper(cohort) &
      wgs$strategy==strategy & wgs$tool=="fastvep"])
  },x$cohort,x$strategy)
  x$effective_speedup <- x$mean_direct_seconds / (x$mean_cached_seconds+
    x$payload_bound_historical_build_seconds/samples)
  x
}))
build_model$strategy_label <- factor(ifelse(build_model$strategy=="cohort_3_genomes",
  "Three-genome cohort cache",ifelse(build_model$strategy=="gnomad_af_0.01",
  "gnomAD AF ≥ 1%","gnomAD AF ≥ 10%")),levels=names(strategy_colours))
build_model$cohort <- factor(toupper(build_model$cohort),levels=c("KPGP","SGDP","PGP"))
write_table(build_model,file.path(root,"tables/recorded_fastvep_build_scenario.tsv"))
p10 <- ggplot(build_model,aes(samples,effective_speedup,colour=strategy_label,linetype=strategy_label)) +
  geom_line(linewidth=.7) + geom_hline(yintercept=1,colour="#A5AAAE",linetype=2,linewidth=.3) +
  facet_wrap(~cohort,nrow=1) + scale_x_log10(breaks=c(5,10,100,1000)) +
  scale_y_continuous(labels=scales::label_number(suffix="×")) +
  scale_colour_manual(values=strategy_colours) +
  labs(x="New genomes processed after building the cache",y="Conditional amortized fastVEP speedup",
    colour=NULL,linetype=NULL) + journal_theme()
save_figure(p10,file.path(figure_root,"supplements/FigureS6_recorded_build_scenario"),11,4)

#### REUSABLE TIBBLES AND REPRODUCIBILITY RECORD ####

saveRDS(list(measurements=cells,plot_views=views,attempts=attempts,assays=assays,
  validated_assays=validated_assays,primary_grch38=primary,wgs_strategies=wgs,
  synthetic=synthetic,assay_summary=assay_summary,primary_summary=primary_summary,
  wgs_summary=wgs_summary,cohort_summary=cohort_summary,phase_long=phase_long,
  phase_summary=phase_summary,headroom=headroom,resources=resources,
  strategy_contrasts=contrasts,illustrative_cohort_model=model,cache_inventory=cache_inventory,
  comparator_audit=comparator_audit,cohort_seed_audit=seed_audit,recorded_build_costs=recorded_builds,
  recorded_build_scenario=build_model,target_headroom=goals),
  file.path(root,"tables/manuscript_data.rds"))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(root,"source/R_sessionInfo.txt"))
writeLines(c("Figures are new outputs; no previous manuscript figure was overwritten.",
 "All plots were generated with R/ggplot2. PNG: 600 dpi; PDF and SVG: vector.",
 "Bootstrap intervals resample independent samples within cohort, not technical runs.",
 "One sample in synthetic scenarios: descriptive points, no inferential interval.",
 "Counterfactual headroom and cohort amortization panels are explicitly models, not new measurements."),
 file.path(figure_root,"RENDER_NOTES.txt"))
print(primary_summary)
