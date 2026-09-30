


source("simulation/seven_method_comparison.R")

chapter4_root <- "chapter4_clinical_trial_illustration"
chapter4_methods <- SEVEN_METHODS

chapter4_config <- function() {


  p_y <- c(`0` = 55 / 111, `1` = 150 / 219)
  list(
    n = 330L,
    n_months = 15L,
    patients_per_month = 22L,
    batch_sizes = rep(22L, 15L),
    p_y = p_y,
    design_p_y = p_y,
    sensitivity = c(`0` = 36 / 55, `1` = 131 / 150),
    specificity = c(`0` = 48 / 56, `1` = 44 / 69),
    investigator_delay_months = 2,
    central_delay_median_months = 3.05,
    central_delay_q90_months = 7.50,
    alpha_one_sided = 0.025,
    final_critical_value = 2.00,
    target_power = 0.80,
    positivity_lower = 0.15,
    positivity_upper = 0.85,


    power_envelope_lower = 0.20,
    power_envelope_upper = 0.78
  )
}

generate_chapter4_stream <- function(config, seed) {
  set.seed(seed)
  n <- config$n
  entry <- rep(seq_len(config$n_months), config$batch_sizes)
  stopifnot(length(entry) == n)
  generate_arm <- function(a) {
    p <- config$p_y[as.character(a)]
    se <- config$sensitivity[as.character(a)]
    sp <- config$specificity[as.character(a)]
    y <- stats::rbinom(n, 1L, p)
    z <- stats::rbinom(n, 1L, ifelse(y == 1L, se, 1 - sp))
    delay <- stats::rlnorm(
      n, log(config$central_delay_median_months),
      central_delay_sdlog(config)
    )
    list(y = y, z = z, delay = delay)
  }
  arm0 <- generate_arm(0)
  arm1 <- generate_arm(1)
  list(
    entry_month = entry,
    z_available_month = entry + config$investigator_delay_months,
    y0 = arm0$y,
    z0 = arm0$z,
    y0_available_month = entry + config$investigator_delay_months + arm0$delay,
    y1 = arm1$y,
    z1 = arm1$z,
    y1_available_month = entry + config$investigator_delay_months + arm1$delay
  )
}

bind_rows_with_missing <- function(rows) {
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (label in missing) row[[label]] <- NA
    row[columns]
  })
  do.call(rbind, rows)
}

simulate_chapter4_replicate <- function(replicate, config,
                                        seed = 2026092801L) {
  stream_seed <- seed + 10000L * replicate
  stream <- generate_chapter4_stream(config, stream_seed)
  gate <- safe_gate_config(config, eta = 1)
  fits <- lapply(seq_along(chapter4_methods), function(j) {
    simulate_seven_method(
      chapter4_methods[j], stream, config, gate,
      assignment_seed = stream_seed + j,
      return_trace = TRUE
    )
  })
  summaries <- bind_rows_with_missing(lapply(fits, `[[`, "summary"))
  traces <- bind_rows_with_missing(lapply(fits, `[[`, "trace"))
  summaries$replicate <- replicate
  traces$replicate <- replicate
  list(summary = summaries, trace = traces)
}

summarize_chapter4 <- function(replications, config) {
  results <- bind_rows_with_missing(lapply(replications, `[[`, "summary"))
  trace <- bind_rows_with_missing(lapply(replications, `[[`, "trace"))
  cr <- results[results$method == "CR", c("replicate", "expected_failures")]
  names(cr)[2] <- "cr_expected_failures"
  results <- merge(results, cr, by = "replicate", all.x = TRUE, sort = FALSE)
  results$failures_avoided_vs_cr <-
    results$cr_expected_failures - results$expected_failures

  pieces <- split(results, factor(results$method, levels = chapter4_methods))
  main <- do.call(rbind, lapply(pieces, function(x) data.frame(
    method = x$method[1],
    allocation_experimental = mean(x$allocation_1),
    expected_progressions = mean(x$expected_failures),
    progressions_avoided_vs_cr = mean(x$failures_avoided_vs_cr),
    bias = mean(x$contrast_stabilized_estimate) - unname(diff(config$p_y)),
    empirical_se = stats::sd(x$contrast_stabilized_estimate),
    mean_estimated_se = mean(x$contrast_stabilized_se),
    power = mean(x$contrast_stabilized_reject),
    coverage = mean(x$contrast_stabilized_cover)
  )))
  rownames(main) <- NULL

  trajectory <- do.call(rbind, lapply(
    split(trace, interaction(factor(trace$method, levels = chapter4_methods),
                             trace$month, drop = TRUE)),
    function(x) data.frame(
      method = x$method[1], month = x$month[1],
      mean_probability = mean(x$randomization_probability),
      q10_probability = unname(stats::quantile(x$randomization_probability, .10)),
      q90_probability = unname(stats::quantile(x$randomization_probability, .90)),
      mean_cumulative_allocation = mean(x$cumulative_allocation_1),
      q10_cumulative_allocation = unname(stats::quantile(x$cumulative_allocation_1, .10)),
      q90_cumulative_allocation = unname(stats::quantile(x$cumulative_allocation_1, .90))
    )
  ))
  trajectory$method <- factor(trajectory$method, levels = chapter4_methods)
  trajectory <- trajectory[order(trajectory$method, trajectory$month), ]
  trajectory$method <- as.character(trajectory$method)
  rownames(trajectory) <- NULL
  list(results = results, trace = trace, main = main, trajectory = trajectory)
}

write_chapter4_outputs <- function(out, config, n_rep) {
  dirs <- file.path(chapter4_root, c("data", "results", "tables"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  saveRDS(out$results, file.path(chapter4_root, "results", "replicate_results.rds"))
  saveRDS(out$trace, file.path(chapter4_root, "results", "monthly_trace.rds"))
  utils::write.csv(out$trajectory,
                   file.path(chapter4_root, "results", "monthly_trajectory.csv"),
                   row.names = FALSE)
  utils::write.csv(out$main,
                   file.path(chapter4_root, "tables", "table2_main_results.csv"),
                   row.names = FALSE)
  meta <- data.frame(
    item = c("Number of methods", "Methods", "Replicates per method",
             "Total trial size", "Accrual", "Final nonprogression rates",
             "Provisional outcome availability", "Additional final-outcome delay",
             "Case-specific allocation envelope"),
    value = c(
      "7", paste(chapter4_methods, collapse = ", "), as.character(n_rep),
      as.character(config$n), "15 monthly batches of 22",
      sprintf("placebo %.4f; cabozantinib %.4f", config$p_y["0"], config$p_y["1"]),
      "2 months after enrollment",
      "log-normal; median 3.05 months; 90th percentile 7.50 months",
      sprintf("[%.2f, %.2f]", config$power_envelope_lower,
              config$power_envelope_upper)
    )
  )
  utils::write.csv(meta, file.path(chapter4_root, "tables", "computation_metadata.csv"),
                   row.names = FALSE)
}

run_chapter4 <- function(n_rep = 10000L, workers = 1L) {
  config <- chapter4_config()
  ids <- seq_len(n_rep)
  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv),
                            envir = .GlobalEnv)
    replications <- parallel::parLapply(
      cl, ids, simulate_chapter4_replicate, config = config
    )
  } else {
    replications <- lapply(ids, simulate_chapter4_replicate, config = config)
  }
  out <- summarize_chapter4(replications, config)
  write_chapter4_outputs(out, config, n_rep)
  print(out$main, row.names = FALSE)
  invisible(out)
}

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}
if (sys.nframe() == 0L) run_chapter4(
  n_rep = as.integer(arg_value("--n-rep=", "10000")),
  workers = as.integer(arg_value("--workers=", "1"))
)
