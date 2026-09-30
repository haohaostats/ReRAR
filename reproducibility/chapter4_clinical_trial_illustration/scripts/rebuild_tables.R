source("chapter4_clinical_trial_illustration/scripts/run_all_methods.R")
results <- readRDS(file.path(chapter4_root, "results", "replicate_results.rds"))
trace <- readRDS(file.path(chapter4_root, "results", "monthly_trace.rds"))
results <- results[setdiff(names(results), c("cr_expected_failures", "failures_avoided_vs_cr"))]
out <- summarize_chapter4(list(list(summary = results, trace = trace)), chapter4_config())
write.csv(out$main, file.path(chapter4_root, "tables", "table2_main_results.csv"), row.names = FALSE)
write.csv(out$trajectory, file.path(chapter4_root, "results", "monthly_trajectory.csv"), row.names = FALSE)
